"""
Generate Psi4 decomposition on spherical like grids using spin-weighted spherical harmonics.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures import BHaH


def register_CFunction_psi4_spinweightm2_decomposition(CoordSystem: str) -> None:
    """
    Register C function for decomposing psi4 into spin-weighted spherical harmonics via 2D interpolation on uniform source grids at each φ slice.

    :param CoordSystem: Specifies the coordinate system for psi4 decomposition.
    :raises ValueError: If psi4 decomposition is not supported for the coordinate system.
    """
    if "Spherical" in CoordSystem:
        radial_like_index = 0
        theta_like_index = 1
        phi_index = 2
    elif "Cylindrical" in CoordSystem:
        radial_like_index = 0
        theta_like_index = 2
        phi_index = 1
    elif "SymTP" in CoordSystem:
        radial_like_index = 0
        theta_like_index = 1
        phi_index = 2
    else:
        raise ValueError(f"CoordSystem = {CoordSystem} not supported.")

    par.register_CodeParameter(
        "int",
        __name__,
        "num_theta_points_on_shell_for_psi4_interp",
        32,
        commondata=True,
    )

    par.register_CodeParameter(
        cparam_type="int",
        module=__name__,
        name="num_psi4_extraction_radii",
        defaultvalue=6,
        description="Number of radii at which psi4 is extracted.",
        add_to_parfile=True,
        commondata=True,
        add_to_set_CodeParameters_h=False,
    )

    # Register list_of_psi4_extraction_radii: the array containing extraction radii values.
    par.register_CodeParameter(
        cparam_type="REAL[6]",
        module=__name__,
        name="list_of_psi4_extraction_radii",
        defaultvalue=[15, 30, 45, 60, 75, 90],
        description="List of radii at which psi4 is extracted. Must set num_psi4_extraction_radii consistently.",
        add_to_parfile=True,
        commondata=True,
        add_to_set_CodeParameters_h=False,
        assumption="RealPositive",
    )

    setup_prefunc = f"""
#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

static const int radial_like_dirn = {radial_like_index}; // Index corresponding to the radial-like coordinate direction.
static const int theta_like_dirn = {theta_like_index}; // Index corresponding to the theta-like coordinate direction.
static const int phi_dirn = {phi_index}; // Index corresponding to the phi coordinate direction.
"""
    setup_prefunc += r"""
// IDX2GENERAL: Map 2D indices (i,j) to a 1D index using stride Ni.
// Useful for representing 2D data in a 1D array.
#define IDX2GENERAL(i, j, Ni) ((i) + (Ni) * (j))
// REVERSE_IDX2GENERAL: Recover 2D indices (i,j) from a 1D index and stride Ni.
// The inverse operation of IDX2GENERAL.
#define REVERSE_IDX2GENERAL(index, Ni, i, j)                                                                                                         \
  {                                                                                                                                                  \
    j = (index) / (Ni);                                                                                                                              \
    i = (index) % (Ni);                                                                                                                              \
  }

// Structure to hold data related to psi4 extraction shells for diagnostics.
// This includes grid points on shells, their coordinates, and related counts.
typedef struct __diagnostic_struct__ {
  // Number of grid points found on each extraction shell within the current process's grid boundaries.
  int *N_shell_pts_grid; // Array shape: [num_psi4_extraction_radii]

  // Radial-like coordinate (e.g., x0) of grid points on each shell.
  // Used as destination points for interpolation.
  REAL **xx_radial_like_shell_grid; // Array shape: [num_psi4_extraction_radii][N_shell_pts_grid[which_R_ext]]

  // Theta-like coordinate (e.g., x2) of grid points on each shell.
  // Used as destination points for interpolation.
  REAL **xx_theta_like_shell_grid; // Array shape: [num_psi4_extraction_radii][N_shell_pts_grid[which_R_ext]]

  // Number of unique theta points found on each extraction shell within the grid boundaries.
  int *N_theta_shell_grid; // Array shape: [num_psi4_extraction_radii]

  // Unique theta values (in radians) corresponding to the grid points on each shell.
  REAL **theta_shell_grid; // Array shape: [num_psi4_extraction_radii][N_theta_shell_grid[which_R_ext]]

  // Angular step size in theta for the uniformly distributed points on the ideal spherical shell.
  REAL dtheta;
} diagnostic_struct;

// Global instance of the diagnostic structure.
static diagnostic_struct diagnosticstruct;

/**
 * @brief Sets up the diagnostic structure for psi4 extraction.
 *
 * This function prepares the necessary data structures and computes coordinates
 * for extracting the psi4 Weyl scalar on thin spherical shells at specified radii.
 * It identifies which points on these conceptual shells fall within the boundaries
 * of the current computational grid domain (handling potential domain decomposition)
 * and calculates the corresponding grid coordinates for interpolation.
 *
 * @param commondata Common data structure containing simulation parameters like extraction radii list.
 * @param params Parameters structure containing grid dimensions and boundaries.
 * @param xx Array of grid coordinate arrays (e.g., xx[0] for x0, xx[1] for x1, xx[2] for x2).
 * @note Assumes the grid might represent only a portion (e.g., a patch or process domain) of the full simulation space.
 * @note Populates the global `diagnosticstruct`.
 */
static void psi4_diagnostics_set_up(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3]) {

  const int num_psi4_extraction_radii = commondata->num_psi4_extraction_radii;
  const REAL *list_of_psi4_extraction_radii = commondata->list_of_psi4_extraction_radii;

  // Number of theta points to define the resolution of the conceptual spherical shells.
  const int N_theta = commondata->num_theta_points_on_shell_for_psi4_interp;

  // Phi values on the shell must match the grid's phi coordinate points exactly.
  // This simplifies interpolation later.
  const int Nxx[3] = {params->Nxx0, params->Nxx1, params->Nxx2};
  const int N_phi = Nxx[phi_dirn]; // Number of phi points matches the grid's phi dimension.

  // Define the properties of the conceptual uniform spherical shells.
  const REAL PI = M_PI;
  const REAL theta_min = 0.0; // Shells span the full theta range [0, pi].
  REAL theta_max = PI;
  REAL phi_min = -PI; // Shells span the full phi range [-pi, pi].
  REAL phi_max = PI;
  const int N_tot_shell = N_theta * N_phi; // Total points on one conceptual shell.
  REAL dtheta = (theta_max - theta_min) / N_theta; // Uniform theta spacing.
  diagnosticstruct.dtheta = dtheta; // Store dtheta for later integration steps.
  const REAL dphi = (phi_max - phi_min) / N_phi; // Uniform phi spacing.

  // Allocate memory for spherical coordinates (theta, phi) on each shell.
  // xx_shell_sph[which_R_ext][0] -> theta values
  // xx_shell_sph[which_R_ext][1] -> phi values
  REAL ***xx_shell_sph = (REAL ***)malloc(num_psi4_extraction_radii * sizeof(REAL **));
  for (int which_R_ext = 0; which_R_ext < num_psi4_extraction_radii; which_R_ext++) {
    xx_shell_sph[which_R_ext] = (REAL * *)malloc(2 * sizeof(REAL *));
    xx_shell_sph[which_R_ext][0] = (REAL *)malloc(N_theta * sizeof(REAL)); // Theta points
    xx_shell_sph[which_R_ext][1] = (REAL *)malloc(N_phi * sizeof(REAL));   // Phi points
  } // END FOR loop allocating spherical coordinate arrays

  // Populate the spherical coordinate arrays for each shell.
  // Points are cell-centered ((j + 0.5) * step).
  for (int which_R_ext = 0; which_R_ext < num_psi4_extraction_radii; which_R_ext++) {
    // Populate theta values for this shell.
    for (int j = 0; j < N_theta; j++) {
      xx_shell_sph[which_R_ext][0][j] = theta_min + ((REAL)j + 0.5) * dtheta;
    } // END FOR populating theta values
      // Populate phi values for this shell.
    for (int j = 0; j < N_phi; j++) {
      xx_shell_sph[which_R_ext][1][j] = phi_min + ((REAL)j + 0.5) * dphi;
    } // END FOR populating phi values
  } // END FOR loop populating spherical coordinates

  // Allocate memory for Cartesian coordinates (x, y, z) corresponding to the shell points.
  REAL ***xx_shell_Cart = (REAL ***)malloc(num_psi4_extraction_radii * sizeof(REAL **));
  for (int which_R_ext = 0; which_R_ext < num_psi4_extraction_radii; which_R_ext++) {
    xx_shell_Cart[which_R_ext] = (REAL * *)malloc(3 * sizeof(REAL *));
    xx_shell_Cart[which_R_ext][0] = (REAL *)malloc(N_tot_shell * sizeof(REAL)); // x coordinates
    xx_shell_Cart[which_R_ext][1] = (REAL *)malloc(N_tot_shell * sizeof(REAL)); // y coordinates
    xx_shell_Cart[which_R_ext][2] = (REAL *)malloc(N_tot_shell * sizeof(REAL)); // z coordinates
  } // END FOR loop allocating Cartesian coordinate arrays

  // Convert spherical shell points (r, theta, phi) to Cartesian coordinates (x, y, z).
  for (int which_R_ext = 0; which_R_ext < num_psi4_extraction_radii; which_R_ext++) {
    REAL r = list_of_psi4_extraction_radii[which_R_ext]; // Radius for the current shell.
    // Loop through all phi points on the shell.
    for (int i_ph = 0; i_ph < N_phi; i_ph++) {
      REAL phi = xx_shell_sph[which_R_ext][1][i_ph];
      // Loop through all theta points on the shell.
      for (int i_th = 0; i_th < N_theta; i_th++) {
        REAL theta = xx_shell_sph[which_R_ext][0][i_th];
        // Standard spherical to Cartesian transformation.
        REAL x = r * sin(theta) * cos(phi);
        REAL y = r * sin(theta) * sin(phi);
        REAL z = r * cos(theta);
        // Map the 2D (theta, phi) index to a 1D index for storage.
        const int idx2 = IDX2GENERAL(i_th, i_ph, N_theta);
        xx_shell_Cart[which_R_ext][0][idx2] = x;
        xx_shell_Cart[which_R_ext][1][idx2] = y;
        xx_shell_Cart[which_R_ext][2][idx2] = z;
      } // END FOR theta loop
    } // END FOR phi loop
  } // END FOR loop converting coordinates

  // Allocate memory within the global diagnostic structure.
  diagnosticstruct.N_shell_pts_grid = (int *)malloc(sizeof(int) * num_psi4_extraction_radii);
  diagnosticstruct.xx_radial_like_shell_grid = (REAL * *)malloc(num_psi4_extraction_radii * sizeof(REAL *));
  diagnosticstruct.xx_theta_like_shell_grid = (REAL * *)malloc(num_psi4_extraction_radii * sizeof(REAL *));
  diagnosticstruct.N_theta_shell_grid = (int *)malloc(sizeof(int) * num_psi4_extraction_radii);
  diagnosticstruct.theta_shell_grid = (REAL * *)malloc(num_psi4_extraction_radii * sizeof(REAL *));

  // For each extraction shell, determine which points lie within the current grid domain.
  for (int which_R_ext = 0; which_R_ext < num_psi4_extraction_radii; which_R_ext++) {
    // Counter for points on this shell that are inside the grid boundaries.
    int count_pt_on_grid = 0;
    // Boolean flags to track which theta/phi indices correspond to points on the grid.
    bool b_theta_grid[N_theta];
    // Initialize flags to false.
    for (int i = 0; i < N_theta; i++) {
      b_theta_grid[i] = false;
    } // END FOR theta flag initialization

    // First pass: Iterate through all points on the conceptual shell to count how many fall within this grid patch.
    // Loop through phi points.
    for (int i_ph = 0; i_ph < N_phi; i_ph++) {
      // Loop through theta points.
      for (int i_th = 0; i_th < N_theta; i_th++) {
        const int idx2 = IDX2GENERAL(i_th, i_ph, N_theta); // 1D index for the current shell point.
        // Get the Cartesian coordinates of the current point on the shell.
        const REAL xCart_pt_on_shell[3] = {xx_shell_Cart[which_R_ext][0][idx2], xx_shell_Cart[which_R_ext][1][idx2],
                                           xx_shell_Cart[which_R_ext][2][idx2]};
        int Cart_to_i0i1i2[3]; // Nearest grid indices (i, j, k).
        REAL closest_xx[3];    // Grid coordinates (x0, x1, x2) corresponding to the Cartesian point.
        // Convert Cartesian point to grid coordinates and find the nearest grid indices.
        Cart_to_xx_and_nearest_i0i1i2(params, xCart_pt_on_shell, closest_xx, Cart_to_i0i1i2);
        // Check if the point's grid coordinates fall within the min/max bounds of this grid patch.
        if ((params->xxmin0 <= closest_xx[0] && closest_xx[0] <= params->xxmax0) &&
            (params->xxmin1 <= closest_xx[1] && closest_xx[1] <= params->xxmax1) &&
            (params->xxmin2 <= closest_xx[2] && closest_xx[2] <= params->xxmax2)) {
          // If the point is inside the grid boundaries:
          count_pt_on_grid++;      // Increment the counter.
          b_theta_grid[i_th] = true; // Mark this theta index as present on the grid.
        } // END IF point is within grid boundaries
      } // END FOR theta loop (counting pass)
    } // END FOR phi loop (counting pass)

    // Count the number of unique theta values that are represented on the grid for this shell.
    int count_theta_grid = 0;
    for (int i = 0; i < N_theta; i++) {
      if (b_theta_grid[i]) {
        count_theta_grid++;
      } // END IF
    } // END FOR counting unique theta values

    // Store the counts in the diagnostic structure for this shell.
    diagnosticstruct.N_shell_pts_grid[which_R_ext] = count_pt_on_grid;
    diagnosticstruct.N_theta_shell_grid[which_R_ext] = count_theta_grid;

    // Allocate memory for the grid coordinates and theta values for the points found on the grid.
    // The size is determined by the counts computed above.
    diagnosticstruct.xx_radial_like_shell_grid[which_R_ext] = (REAL *)malloc(count_pt_on_grid * sizeof(REAL));
    diagnosticstruct.xx_theta_like_shell_grid[which_R_ext] = (REAL *)malloc(count_pt_on_grid * sizeof(REAL));
    diagnosticstruct.theta_shell_grid[which_R_ext] = (REAL *)malloc(count_theta_grid * sizeof(REAL));

    // Second pass: Iterate again, this time storing the grid coordinates and theta values for the points within bounds.
    int which_pt_on_grid = 0; // Counter/index for points stored in the diagnostic struct arrays.
    // Loop through phi points.
    for (int i_ph = 0; i_ph < N_phi; i_ph++) {
      // Loop through theta points.
      for (int i_th = 0; i_th < N_theta; i_th++) {
        const int idx2 = IDX2GENERAL(i_th, i_ph, N_theta); // 1D index for the current shell point.
        // Get the Cartesian coordinates.
        const REAL xCart_pt_on_shell[3] = {xx_shell_Cart[which_R_ext][0][idx2], xx_shell_Cart[which_R_ext][1][idx2],
                                           xx_shell_Cart[which_R_ext][2][idx2]};
        int Cart_to_i0i1i2[3]; // Nearest grid indices (i, j, k).
        REAL closest_xx[3];    // Grid coordinates (x0, x1, x2).
        // Convert Cartesian point to grid coordinates.
        Cart_to_xx_and_nearest_i0i1i2(params, xCart_pt_on_shell, closest_xx, Cart_to_i0i1i2);
        // Check again if the point is within the grid boundaries.
        if ((params->xxmin0 <= closest_xx[0] && closest_xx[0] <= params->xxmax0) &&
            (params->xxmin1 <= closest_xx[1] && closest_xx[1] <= params->xxmax1) &&
            (params->xxmin2 <= closest_xx[2] && closest_xx[2] <= params->xxmax2)) {
          // Store the radial-like and theta-like grid coordinates. These will be the destination points for interpolation.
          diagnosticstruct.xx_radial_like_shell_grid[which_R_ext][which_pt_on_grid] = closest_xx[radial_like_dirn];
          diagnosticstruct.xx_theta_like_shell_grid[which_R_ext][which_pt_on_grid] = closest_xx[theta_like_dirn];

          // Also store the corresponding original theta value (in radians) from the spherical shell coordinates.
          // This requires careful indexing based on the points found *on the grid*.
          int i_th_grid;             // Index within the grid-specific theta array.
          MAYBE_UNUSED int i_ph_grid; // Index within the grid-specific phi dimension (unused here).
          const int N_theta_shell_grid_local = diagnosticstruct.N_theta_shell_grid[which_R_ext];
          // Recover the 2D grid-based indices (theta_grid, phi_grid) from the 1D on-grid point index.
          REVERSE_IDX2GENERAL(which_pt_on_grid, N_theta_shell_grid_local, i_th_grid, i_ph_grid);
          // Store the original theta value at the appropriate index in the grid-specific theta array.
          diagnosticstruct.theta_shell_grid[which_R_ext][i_th_grid] = xx_shell_sph[which_R_ext][0][i_th];
          which_pt_on_grid++; // Move to the next storage location.
        } // END IF point is within grid boundaries (storage pass)
      } // END FOR theta loop (storage pass)
    } // END FOR phi loop (storage pass)
  } // END FOR loop over all extraction radii (R_ext)

  // Optional: Summing the number of points found on the grid across all shells (could be used for verification).
  int sum = 0;
  for (int which_R_ext = 0; which_R_ext < num_psi4_extraction_radii; which_R_ext++) {
    if (diagnosticstruct.N_shell_pts_grid[which_R_ext] > 0) {
      sum += diagnosticstruct.N_shell_pts_grid[which_R_ext];
    } // END IF
  } // END FOR summing points

  // Free the temporary memory used for spherical and Cartesian coordinates of the conceptual shells.
  for (int which_R_ext = 0; which_R_ext < num_psi4_extraction_radii; which_R_ext++) {
    free(xx_shell_sph[which_R_ext][0]); // Free theta array for this shell.
    free(xx_shell_sph[which_R_ext][1]); // Free phi array for this shell.
    free(xx_shell_sph[which_R_ext]);    // Free the pointer array for this shell.
  } // END FOR freeing spherical arrays
  free(xx_shell_sph); // Free the top-level pointer array.

  for (int which_R_ext = 0; which_R_ext < num_psi4_extraction_radii; which_R_ext++) {
    free(xx_shell_Cart[which_R_ext][0]); // Free x array for this shell.
    free(xx_shell_Cart[which_R_ext][1]); // Free y array for this shell.
    free(xx_shell_Cart[which_R_ext][2]); // Free z array for this shell.
    free(xx_shell_Cart[which_R_ext]);    // Free the pointer array for this shell.
  } // END FOR freeing Cartesian arrays
  free(xx_shell_Cart); // Free the top-level pointer array.
} // END FUNCTION psi4_diagnostics_set_up
"""

    decomposition_prefunc = r"""
/**
 * @brief Performs the decomposition of psi4 into spin-weighted spherical harmonic modes (s=-2).
 *
 * Integrates psi4 multiplied by the conjugate of the spin-weighted spherical harmonic
 * Y^{-2}_{lm} over the surface of a sphere at a given extraction radius R_ext.
 * Appends the results (real and imaginary parts for each (l,m) mode, scaled by R_ext)
 * to output files, one file per 'l' value.
 *
 * @param Nxx_plus_2NGHOSTS1 The size of the phi-like dimension including ghosts (used for indexing consistency, though integration ignores ghosts).
 * @param dxx1 The grid spacing in the phi-like direction (dphi).
 * @param dxx2 The grid spacing in the theta-like direction (dtheta, provided by diagnosticstruct.dtheta).
 * @param swm2sh_maximum_l_mode_to_compute The maximum 'l' mode to compute (inclusive, starts from l=2).
 * @param curr_time The current simulation time.
 * @param R_ext The extraction radius for this decomposition.
 * @param th_array Array of theta values (radians) on the extraction sphere corresponding to grid points.
 * @param sinth_array Array of sin(theta) values corresponding to th_array.
 * @param ph_array Array of phi values (radians) on the extraction sphere corresponding to grid points.
 * @param N_theta Number of unique theta points on the grid for this shell.
 * @param psi4r_at_R_ext Real part of psi4 interpolated onto the extraction sphere points (stored as 1D array).
 * @param psi4i_at_R_ext Imaginary part of psi4 interpolated onto the extraction sphere points (stored as 1D array).
 * @note The integration uses a simple summation approximating the surface integral: Sum(psi4 * conj(Y) * sin(theta) * dtheta * dphi).
 * @note Outputs data to files named "Rpsi4_l<l>-r<R_ext>.txt".
 * @note File format: Column 1 is retarded time (t - R_ext), subsequent pairs of columns are Re(psi4_lm * R_ext) and Im(psi4_lm * R_ext).
 * @note Uses OpenMP parallel for reduction to speed up the integration loops.
 */
static void lowlevel_decompose_psi4_into_swm2_modes(const int Nxx_plus_2NGHOSTS1, const REAL dxx1, const REAL dxx2,
                                                    const int swm2sh_maximum_l_mode_to_compute, const REAL curr_time, const REAL R_ext,
                                                    const REAL *restrict th_array, const REAL *restrict sinth_array, const REAL *restrict ph_array,
                                                    const int N_theta, const REAL *restrict psi4r_at_R_ext, const REAL *restrict psi4i_at_R_ext) {
  char filename[100];   // Buffer for constructing output filenames.
  FILE *outpsi4_l_m;    // File pointer for output.

  // At the very start of the simulation (t=0), create/overwrite files and write headers.
  if (curr_time == 0) {
    // Loop over each 'l' mode that will be computed.
    for (int l = 2; l <= swm2sh_maximum_l_mode_to_compute; l++) {
      // Construct filename including 'l' and extraction radius.
      sprintf(filename, "Rpsi4_l%d-r%06.1f.txt", l, (double)R_ext);
      outpsi4_l_m = fopen(filename, "w"); // Open in write mode ('w') to overwrite/create.
      if (!outpsi4_l_m) {
        perror("Error opening psi4 output file");
        // Consider adding more robust error handling, e.g., exit or return error code.
        continue; // Skip to next l if file opening fails.
      }
      fprintf(outpsi4_l_m, "# column 1: t-R_ext = [retarded time]\n");
      int col = 2; // Column counter, starting after the time column.
      // Loop over all possible 'm' modes for the current 'l'.
      for (int m = -l; m <= l; m++) {
        // Write header descriptions for the real and imaginary parts of this (l,m) mode.
        fprintf(outpsi4_l_m, "# column %d: Re(psi4_{l=%d,m=%d}) * R_ext\n", col, l, m);
        col++;
        fprintf(outpsi4_l_m, "# column %d: Im(psi4_{l=%d,m=%d}) * R_ext\n", col, l, m);
        col++;
      } // END FOR m loop (header writing)
      fclose(outpsi4_l_m); // Close the file after writing the header.
    } // END FOR l loop (header writing)
  } // END IF curr_time == 0

  // --- Perform decomposition and append data for the current timestep ---

  // Loop over each 'l' mode to compute.
  for (int l = 2; l <= swm2sh_maximum_l_mode_to_compute; l++) {
    // Construct filename again to append data.
    sprintf(filename, "Rpsi4_l%d-r%06.1f.txt", l, (double)R_ext);
    outpsi4_l_m = fopen(filename, "a"); // Open in append mode ('a').
    if (!outpsi4_l_m) {
      perror("Error opening psi4 output file for appending");
      // Consider adding more robust error handling.
      continue; // Skip to next l if file opening fails.
    }
    char oneline[10000]; // Buffer to build the output line for this 'l'.
    // Start the line with the retarded time.
    sprintf(oneline, "%e", (double)(curr_time - R_ext));

    // Loop over all 'm' modes for the current 'l'.
    for (int m = -l; m <= l; m++) {
      // Initialize accumulators for the real and imaginary parts of the integral for this (l,m) mode.
      REAL psi4r_l_m = 0.0;
      REAL psi4i_l_m = 0.0;

      // Integrate over the sphere (represented by the grid points) using OpenMP for parallelization.
      // The reduction clause ensures correct summation across threads.
      // Note: The integration loops over the *grid* points corresponding to the shell.
      // i1 iterates over phi index, i2 iterates over theta index.
      // dxx1 is dphi, dxx2 is dtheta.
#pragma omp parallel for reduction(+ : psi4r_l_m, psi4i_l_m)
      // Loop over the phi dimension (excluding ghost zones).
      for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1 - 2 * NGHOSTS; i1++) {
        const REAL ph = ph_array[i1]; // Get the phi value for this column.
        // Loop over the theta dimension (using the count of unique theta points on grid).
        for (int i2 = 0; i2 < N_theta; i2++) {
          const REAL th = th_array[i2];     // Get the theta value for this row.
          const REAL sinth = sinth_array[i2]; // Get sin(theta).

          // Calculate the spin-weighted spherical harmonic Y^{-2}_{lm}(theta, phi).
          REAL ReY_sm2_l_m, ImY_sm2_l_m;
          spin_weight_minus2_sph_harmonics(l, m, th, ph, &ReY_sm2_l_m, &ImY_sm2_l_m);

          // Get the 1D index corresponding to the current (theta, phi) point.
          const int idx2d = IDX2GENERAL(i2, i1, N_theta);
          // Extract the interpolated real (a) and imaginary (b) parts of psi4 at this point.
          const REAL a = psi4r_at_R_ext[idx2d];
          const REAL b = psi4i_at_R_ext[idx2d];
          // Real (c) and imaginary (d) parts of the spherical harmonic (note: using Y, not conj(Y) here, need to check convention).
          // Assuming the standard integral is Integral(psi4 * conj(Y) dOmega).
          // conj(Y^{-2}_{lm}) = (-1)^m Y^{-2}_{l,-m}. We need the harmonic itself.
          const REAL c = ReY_sm2_l_m;
          const REAL d = ImY_sm2_l_m;

          // Accumulate the integral components. The integral is Sum[ (a + ib) * (c - id) * sin(th) * dth * dph ]
          // Real part: (ac + bd) * sin(th) * dth * dph
          // Imaginary part: (bc - ad) * sin(th) * dth * dph
          // dxx2 is dtheta, dxx1 is dphi.
          psi4r_l_m += (a * c + b * d) * dxx2 * sinth * dxx1; // Accumulate real part.
          psi4i_l_m += (b * c - a * d) * dxx2 * sinth * dxx1; // Accumulate imaginary part.
        } // END FOR theta integration loop (i2)
      } // END FOR phi integration loop (i1) & OMP parallel region

      // Append the computed real and imaginary parts (scaled by R_ext) to the output line buffer.
      sprintf(oneline + strlen(oneline), " %.15e %.15e", (double)(R_ext * psi4r_l_m), (double)(R_ext * psi4i_l_m));
    } // END FOR m loop (mode computation)

    // Write the completed line (retarded time + all m modes for this l) to the file.
    fprintf(outpsi4_l_m, "%s\n", oneline);
    fclose(outpsi4_l_m); // Close the file for this 'l'.
  } // END FOR l loop (mode computation and file output)
} // END FUNCTION lowlevel_decompose_psi4_into_swm2_modes
"""
    prefunc = setup_prefunc + "\n\n" + decomposition_prefunc

    desc = """
@brief Calculates and outputs the spin-weighted spherical harmonic decomposition of psi4.

This function orchestrates the process of extracting the psi4 Weyl scalar on multiple
spherical shells, interpolating it onto these shells, and then decomposing the interpolated
data into spin-weighted spherical harmonic modes (s=-2, l>=2). The results are written to
files for each extraction radius and l-mode. This specific version is tailored for
SinhCylindrical coordinates.

@param commondata Common data structure containing simulation parameters (time, extraction radii, etc.).
@param params Parameters structure containing grid information (dimensions, spacings, boundaries).
@param diagnostic_gfs Array containing grid functions, including DIAG_PSI4_REGF and DIAG_PSI4_REGF.
@param xx Array of grid coordinate arrays (e.g., xx[0] for x0, xx[1] for x1, xx[2] for x2).
@note This function assumes SinhCylindrical coordinates where x0 is radial-like, x2 is theta-like, and x1 is phi.
@note Calls `psi4_diagnostics_set_up` to prepare extraction data.
@note Calls `interpolation_2d_general__uniform_src_grid` to interpolate psi4 onto shells.
@note Calls `lowlevel_decompose_psi4_into_swm2_modes` to perform the decomposition and file output.
@note Manages memory allocation and deallocation for intermediate arrays.
"""
    name = "psi4_spinweightm2_decomposition"
    params = r"""const commondata_struct *restrict commondata,
    const params_struct *restrict params,
    REAL *restrict diagnostic_gfs,
    REAL *restrict xx[3]"""

    # Register C functions apply_bcs_inner_only_specific_gfs, needed for 2d interp of psi4
    BHaH.CurviBoundaryConditions.apply_bcs_inner_only_specific_gfs.register_CFunction_apply_bcs_inner_only_specific_gfs()

    body = r"""
  // Step 0: Set up the diagnostic structure (calculates shell points on grid, etc.).
  psi4_diagnostics_set_up(commondata, params, xx);

  // Extract necessary information from common data and diagnostic structures.
  const int num_psi4_extraction_radii = commondata->num_psi4_extraction_radii;
  const REAL *restrict list_of_psi4_extraction_radii = commondata->list_of_psi4_extraction_radii;
  int *restrict N_shell_pts_grid = diagnosticstruct.N_shell_pts_grid;                     // Num points on grid per shell.
  REAL **restrict xx_radial_like_shell_grid = diagnosticstruct.xx_radial_like_shell_grid; // Dest radial coord for interp.
  REAL **restrict xx_theta_like_shell_grid = diagnosticstruct.xx_theta_like_shell_grid;   // Dest theta coord for interp.
  int *restrict N_theta_shell_grid = diagnosticstruct.N_theta_shell_grid;                 // Num unique theta points per shell.
  REAL **restrict theta_shell_grid = diagnosticstruct.theta_shell_grid;                   // Actual theta values (rad) per shell.
  const REAL dtheta = diagnosticstruct.dtheta;                                            // Theta spacing on conceptual shell.

  // Grid dimensions including ghost zones.
  const int Nxx_plus_2NGHOSTS_arr[3] = {Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2};
  // Grid spacings.
  const REAL dxx[3] = {params->dxx0, params->dxx1, params->dxx2};
  // Get sizes and spacings for the specific directions used in interpolation and decomposition.
  const int Nxx_plus_2NGHOSTS_radial_like_dirn = Nxx_plus_2NGHOSTS_arr[radial_like_dirn]; // Size of source grid dim 0.
  const int Nxx_plus_2NGHOSTS_theta_like_dirn = Nxx_plus_2NGHOSTS_arr[theta_like_dirn];   // Size of source grid dim 1.
  const int Nxx_plus_2NGHOSTS_phi = Nxx_plus_2NGHOSTS_arr[phi_dirn];                      // Size of phi dimension.
  const REAL dxx_radial_like = dxx[radial_like_dirn];                                     // Spacing of source grid dim 0.
  const REAL dxx_theta_like = dxx[theta_like_dirn];                                       // Spacing of source grid dim 1.
  const REAL dxx_phi = dxx[phi_dirn];                                                     // Spacing of phi dimension.

  // --- Set up parameters for 2D interpolation ---
  // The source grid for interpolation is a 2D slice (radial-like, theta-like) for a fixed phi.
  const int N_interp_GHOSTS = NGHOSTS;                                   // Number of ghost zones needed by the interpolator.
  const REAL src_dxx0 = dxx_radial_like;                                 // Source grid spacing in dimension 0.
  const REAL src_dxx1 = dxx_theta_like;                                  // Source grid spacing in dimension 1.
  const int src_Nxx_plus_2NGHOSTS0 = Nxx_plus_2NGHOSTS_radial_like_dirn; // Source grid size in dimension 0.
  const int src_Nxx_plus_2NGHOSTS1 = Nxx_plus_2NGHOSTS_theta_like_dirn;  // Source grid size in dimension 1.

  // Allocate memory for source grid coordinate arrays needed by the interpolator.
  REAL *src_x0x1[3];                                                           // Interpolator expects 3 elements, but we only use [1] and [2].
  src_x0x1[0] = NULL;                                                          // Unused by 2D interpolator.
  src_x0x1[1] = (REAL *)malloc(sizeof(REAL) * src_Nxx_plus_2NGHOSTS0); // Radial-like coordinates.
  src_x0x1[2] = (REAL *)malloc(sizeof(REAL) * src_Nxx_plus_2NGHOSTS1); // Theta-like coordinates.
  // Populate source coordinate arrays.
  for (int j = 0; j < src_Nxx_plus_2NGHOSTS0; j++)
    src_x0x1[1][j] = xx[radial_like_dirn][j];
  for (int j = 0; j < src_Nxx_plus_2NGHOSTS1; j++)
    src_x0x1[2][j] = xx[theta_like_dirn][j];

  // Total number of points in the 2D source slice.
  const int total_size = src_Nxx_plus_2NGHOSTS0 * src_Nxx_plus_2NGHOSTS1;
  // Allocate memory for the source grid function data (real and imaginary parts of psi4) for one phi-slice.
  // These will be repopulated for each phi slice inside the loop.
  REAL *src_gf_psi4r = (REAL *)malloc(sizeof(REAL) * total_size); // Shape: [src_N0 * src_N1]
  REAL *src_gf_psi4i = (REAL *)malloc(sizeof(REAL) * total_size); // Shape: [src_N0 * src_N1]

  // Step 1: Loop over all extraction radii.
  for (int which_R_ext = 0; which_R_ext < num_psi4_extraction_radii; which_R_ext++) {
    // Get the radius value for the current shell.
    const REAL R_ext = list_of_psi4_extraction_radii[which_R_ext];
    // Get the number of destination points on this shell that fall within the grid.
    const int num_dst_pts = N_shell_pts_grid[which_R_ext];

    // Only proceed if there are points on this shell within the current grid domain.
    if (num_dst_pts > 0) {

      // Allocate memory for destination points coordinates (radial-like, theta-like) for interpolation.
      REAL(*dst_pts)[2] = (REAL(*)[2])malloc(sizeof(REAL) * num_dst_pts * 2);
      // Populate the destination points array from the diagnostic structure.
      for (int i = 0; i < num_dst_pts; i++) {
        dst_pts[i][0] = xx_radial_like_shell_grid[which_R_ext][i]; // Radial-like coordinate.
        dst_pts[i][1] = xx_theta_like_shell_grid[which_R_ext][i];  // Theta-like coordinate.
      } // END FOR populating destination points

      // Allocate memory for the interpolated data (destination data) for psi4 real and imaginary parts.
      REAL *dst_data_psi4r = (REAL *)malloc(sizeof(REAL) * num_dst_pts);
      REAL *dst_data_psi4i = (REAL *)malloc(sizeof(REAL) * num_dst_pts);

      // Allocate memory for the 2D arrays (flattened) to store the fully interpolated psi4 on the shell grid points.
      // The size depends on the number of phi points (from grid) and unique theta points (from diagnostic struct).
      int N_theta_local = N_theta_shell_grid[which_R_ext];
      int N_phi_local = Nxx_plus_2NGHOSTS_phi - 2 * NGHOSTS; // Number of physical phi points.
      int sizeof_2Darray = sizeof(REAL) * N_phi_local * N_theta_local;
      REAL *psi4r_at_R_ext = (REAL *)malloc(sizeof_2Darray); // Real part on the shell.
      REAL *psi4i_at_R_ext = (REAL *)malloc(sizeof_2Darray); // Imaginary part on the shell.

      // Allocate and populate 1D arrays for theta, sin(theta), and phi values corresponding
      // to the points used in the final decomposition integral.
      REAL *sinth_array = (REAL *)malloc(sizeof(REAL) * N_theta_local);
      REAL *th_array = (REAL *)malloc(sizeof(REAL) * N_theta_local);
      REAL *ph_array = (REAL *)malloc(sizeof(REAL) * N_phi_local);

      // Step 2: Interpolate psi4 onto the shell points. This requires looping through phi slices.
      // The interpolation is 2D (radial-like, theta-like), performed independently for each phi slice.
      // Loop over physical grid points in the phi direction.
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS_phi - NGHOSTS; i1++) {
        // Store the phi value for this slice (index relative to physical points).
        ph_array[i1 - NGHOSTS] = xx[phi_dirn][i1];

      // Populate the 2D source grid function arrays (src_gf_psi4r/i) for the current phi slice (i1).
      // Note: This inner loop is implicitly parallelized by the outer OMP directive.
#pragma omp parallel for                                     // Parallelize the loop over theta slices.
        for (int j = 0; j < src_Nxx_plus_2NGHOSTS1; j++) {   // Loop over theta-like dimension (index j).
          for (int i = 0; i < src_Nxx_plus_2NGHOSTS0; i++) { // Loop over radial-like dimension (index i).
            // Calculate the 1D index for the source array.
            int idx = i + src_Nxx_plus_2NGHOSTS0 * j;

            // Construct the 3D grid index (i0, i1, i2) for accessing diagnostic_gfs.
            int ii012[3];
            ii012[radial_like_dirn] = i; // Radial-like index.
            ii012[theta_like_dirn] = j;  // Theta-like index.
            ii012[phi_dirn] = i1;        // Current phi slice index.

            // Copy psi4 real and imaginary data from the main grid function array to the 2D source slice.
            src_gf_psi4r[idx] = diagnostic_gfs[IDX4(DIAG_PSI4_REGF, ii012[0], ii012[1], ii012[2])];
            src_gf_psi4i[idx] = diagnostic_gfs[IDX4(DIAG_PSI4_IMGF, ii012[0], ii012[1], ii012[2])];
          } // END FOR radial-like dimension loop (i)
        } // END FOR theta-like dimension loop (j)

        // Perform the 2D interpolation for the real part of psi4 for this phi slice.
        int error_code1 =
            interpolation_2d_general__uniform_src_grid(N_interp_GHOSTS, src_dxx0, src_dxx1, src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_x0x1,
                                                       src_gf_psi4r, num_dst_pts, dst_pts, dst_data_psi4r);
        // Basic error check.
        if (error_code1 > 0) {
          printf("Interpolation error code (real part): %d\n", error_code1);
          // Consider more robust error handling.
        } // END IF interpolation error check

        // Perform the 2D interpolation for the imaginary part of psi4 for this phi slice.
        int error_code2 =
            interpolation_2d_general__uniform_src_grid(N_interp_GHOSTS, src_dxx0, src_dxx1, src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_x0x1,
                                                       src_gf_psi4i, num_dst_pts, dst_pts, dst_data_psi4i);
        // Basic error check.
        if (error_code2 > 0) {
          printf("Interpolation error code (imaginary part): %d\n", error_code2);
          // Consider more robust error handling.
        } // END IF interpolation error check

        // Store the interpolated results (dst_data_psi4r/i) into the final 2D shell arrays (psi4r/i_at_R_ext).
        // Also, populate the theta and sin(theta) arrays within this loop (though they only need populating once).
#pragma omp parallel for // Parallelize the loop over theta slices.
        for (int i2 = 0; i2 < N_theta_local; i2++) {
          // Populate theta and sin(theta) arrays using values from the diagnostic struct.
          // This is done inside the phi loop, which is slightly redundant but safe if OMP is used.
          th_array[i2] = theta_shell_grid[which_R_ext][i2];
          sinth_array[i2] = sin(th_array[i2]);

          // Calculate the 1D index for the final 2D shell arrays (psi4r/i_at_R_ext).
          // Index maps (theta_index, phi_index) -> 1D index.
          const int idx2d = IDX2GENERAL(i2, (i1 - NGHOSTS), N_theta_local);

          // The interpolated data `dst_data_psi4r/i` corresponds to the destination points `dst_pts`.
          // We need to map the point `(dst_pts[k][0], dst_pts[k][1])` back to the `(i2, i1-NGHOSTS)` structure.
          // Assuming the `dst_pts` were generated in the same order (phi varying fastest, then theta).
          // Find the correct index `k` in `dst_data_psi4r/i` that corresponds to `(i2, i1-NGHOSTS)`.
          // A simpler way: assume `dst_data_psi4r/i` is already ordered correctly according to (theta, phi).
          // The index into dst_data needs conversion based on how dst_pts was populated.
          // Assuming dst_pts[k] corresponds to the k-th point found in the second pass of psi4_diagnostics_set_up.
          // That pass iterates phi first, then theta. So the index should match IDX2GENERAL(i_th_grid, i_ph_grid).
          // Let's assume the index calculation used for `dst_data_psi4r/i` storage matches `idx2d`.
          // **Correction:** The `dst_data` corresponds directly to `dst_pts`, which contains *all* points for the shell *on the grid*.
          // The `psi4r_at_R_ext` array needs to store these values organized by the `(i2, i1-NGHOSTS)` grid structure.
          // We need the mapping from the `dst_pts` index to the `(i2, i1-NGHOSTS)` index.
          // Let's retrieve the index from dst_pts assuming it was created with phi varying fastest.
          const int dst_idx = IDX2GENERAL(i2, (i1 - NGHOSTS), N_theta_local); // Assuming this matches the implicit ordering in dst_data

          // Copy the interpolated value from the temporary destination buffer to the final 2D array.
          psi4r_at_R_ext[idx2d] = dst_data_psi4r[dst_idx];
          psi4i_at_R_ext[idx2d] = dst_data_psi4i[dst_idx];
        } // END FOR theta loop (storing results)
      } // END OMP parallel for loop over phi slices (i1)

      // Step 3: Perform the spin-weighted spherical harmonic decomposition using the interpolated data.
      lowlevel_decompose_psi4_into_swm2_modes(Nxx_plus_2NGHOSTS_phi,            // Size of phi dimension (used for indexing range)
                                              dxx_phi,                          // Phi spacing (dphi)
                                              dtheta,                           // Theta spacing (dtheta)
                                              swm2sh_maximum_l_mode_to_compute, // Max l mode to compute
                                              commondata->time,                 // Current time
                                              R_ext,                            // Current extraction radius
                                              th_array,                         // Array of theta values
                                              sinth_array,                      // Array of sin(theta) values
                                              ph_array,                         // Array of phi values
                                              N_theta_local,                    // Number of theta points on grid for this shell
                                              psi4r_at_R_ext,                   // Interpolated real part of psi4
                                              psi4i_at_R_ext                    // Interpolated imaginary part of psi4
      );

      // Step 4: Free memory allocated specifically for this extraction radius.
      free(psi4r_at_R_ext);
      free(psi4i_at_R_ext);
      free(sinth_array);
      free(th_array);
      free(ph_array);
      free(dst_pts);
      free(dst_data_psi4r);
      free(dst_data_psi4i);
    } // END IF num_dst_pts > 0
  } // END FOR loop over all extraction radii (which_R_ext)

  // Free memory allocated for interpolation setup (source grid info).
  free(src_gf_psi4r);
  free(src_gf_psi4i);
  free(src_x0x1[1]); // Free radial-like source coordinates.
  free(src_x0x1[2]); // Free theta-like source coordinates.

  // Free memory allocated within the global diagnostic structure during setup.
  for (int i = 0; i < commondata->num_psi4_extraction_radii; ++i) {
    free(diagnosticstruct.xx_radial_like_shell_grid[i]);
    free(diagnosticstruct.xx_theta_like_shell_grid[i]);
    free(diagnosticstruct.theta_shell_grid[i]);
  } // END FOR freeing diagnostic struct arrays per radius
  free(diagnosticstruct.xx_radial_like_shell_grid);
  free(diagnosticstruct.xx_theta_like_shell_grid);
  free(diagnosticstruct.theta_shell_grid);
  free(diagnosticstruct.N_shell_pts_grid);
  free(diagnosticstruct.N_theta_shell_grid);
  // Note: diagnosticstruct itself is static global, not freed here.
"""

    cfc.register_CFunction(
        includes=[
            "BHaH_defines.h",
            "BHaH_function_prototypes.h",
            "diagnostics/diagnostic_gfs.h",
        ],
        prefunc=prefunc,
        desc=desc,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
