"""
Register C functions for 2D angular interpolation at arbitrary sets of points.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.BHaH.BHaHAHA import error_message


def register_CFunction_interpolation_2d_general__uniform_src_grid(
    enable_simd: bool,
    project_dir: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for general-purpose 2D Lagrange interpolation.
    :param enable_simd: Whether the rest of the code enables SIMD optimizations, as this code requires simd_intrinsics.h (which includes SIMD-disabled options).
    :param project_dir: Directory of the project, to set the destination for simd_instrinsics.h .

    >>> result = register_CFunction_interpolation_2d_general__uniform_src_grid()
    >>> result is None or isinstance(result, pcg.NRPyEnv_type)
    True

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    if not enable_simd:
        copy_files(
            package="nrpy.helpers",
            filenames_list=["simd_intrinsics.h"],
            project_dir=project_dir,
            subdirectory="intrinsics",
        )
    copy_files(
        package="nrpy.infrastructures.BHaH.interpolation",
        filenames_list=["interpolation_lagrange_uniform.h"],
        project_dir=project_dir,
        subdirectory="./",
    )

    includes = ["stdio.h", "stdlib.h", "math.h", "interpolation_lagrange_uniform.h"]

    prefunc = """
#ifndef REAL
#define REAL double
#endif
#define DEBUG

#ifdef STANDALONE
// Remove the bah_ prefix if compiling as a standalone code, as this function goes by other names in other codes.
#define bah_interpolation_2d_general__uniform_src_grid interpolation_2d_general__uniform_src_grid
#endif

// In case this code is compiled as C++:
#ifdef __cplusplus
#ifndef restrict
#define restrict __restrict__
#endif
#endif
"""
    prefunc += """//===============================================
// BHaHAHA Error Handling
//===============================================
// Error codes, set in error_message.py
typedef enum {
"""
    for item in error_message.error_code_msg_tuples_list:
        prefunc += f"  {item[0]},\n"
    prefunc += """} bhahaha_error_codes;
#pragma GCC optimize("unroll-loops")"""
    desc = r"""
Performs 2D Lagrange interpolation on a uniform source grid.

This function interpolates a scalar grid function from a uniform source grid defined in theta and phi
to a set of arbitrary destination points using Lagrange interpolation of a specified order.

@param N_interp_GHOSTS          Number of ghost zones around each source point. Determines the interpolation order as 2 * N_interp_GHOSTS + 1.
@param src_dxx1                 Grid spacing in the theta direction on the source grid.
@param src_dxx2                 Grid spacing in the phi direction on the source grid.
@param src_Nxx_plus_2NGHOSTS1  Total number of grid points in the theta direction, including ghost zones.
@param src_Nxx_plus_2NGHOSTS2  Total number of grid points in the phi direction, including ghost zones.
@param src_r_theta_phi          Arrays containing coordinate values for r, theta, and phi on the source grid.
@param src_gf                   Pointer to the source grid function data, organized as a flattened 2D array.
@param num_dst_pts              Number of destination points where interpolation is performed.
@param dst_pts                  Array of destination points' coordinates, each consisting of (theta, phi).
@param dst_data                 Output array to store interpolated values at each destination point.

@return                         BHAHAHA_SUCCESS on successful interpolation.
                                 Appropriate error code if an error is encountered.

@note
- Assumes that the source and destination grids are uniform in theta and phi directions.
- Ensures that destination points lie within the bounds of the source grid to prevent memory access violations.
"""
    cfunc_type = "int"
    name = "interpolation_2d_general__uniform_src_grid"
    params = """const int n_interp_ghosts, const REAL src_dxx1, const REAL src_dxx2,
                const int src_Nxx_plus_2NGHOSTS1, const int src_Nxx_plus_2NGHOSTS2,
                REAL *restrict src_r_theta_phi[3], const REAL *restrict src_gf, const int num_dst_pts,
                const REAL dst_pts[][2], REAL *restrict dst_data"""
    body = r"""
  // Define the interpolation order based on the number of ghost zones.
  const int INTERP_ORDER = (2 * n_interp_ghosts + 1); // Interpolation order corresponds to the number of points in the stencil per dimension.

  // Calculate inverse grid spacings for efficient index calculations.
  const REAL src_invdxx1 = 1.0 / src_dxx1;
  const REAL src_invdxx2 = 1.0 / src_dxx2;

  // Compute normalization factor to scale the interpolated result appropriately.
  const REAL src_invdxx12_INTERP_ORDERm1 = pow(src_dxx1 * src_dxx2, -(INTERP_ORDER - 1));

  // Validate input pointers to prevent segmentation faults.
  if (src_r_theta_phi[1] == NULL || src_r_theta_phi[2] == NULL || src_gf == NULL || dst_data == NULL)
    return INTERP2D_GENERAL_NULL_PTRS; // Exit if any required pointer is NULL.

  // Ensure that the interpolation order does not exceed the grid dimensions in either direction.
  if (INTERP_ORDER > src_Nxx_plus_2NGHOSTS1 || INTERP_ORDER > src_Nxx_plus_2NGHOSTS2)
    return INTERP2D_GENERAL_INTERP_ORDER_GT_NXX_PLUS_2NGHOSTS12; // Exit if interpolation order is too high.

  // Precompute inverse denominators for Lagrange interpolation coefficients to reduce redundant calculations.
  REAL inv_denom[INTERP_ORDER];
  compute_inv_denom(INTERP_ORDER, inv_denom);

  // Define the minimum coordinate values including ghost zones for both theta and phi.
  const REAL xxmin_incl_ghosts1 = src_r_theta_phi[1][0];
  const REAL xxmin_incl_ghosts2 = src_r_theta_phi[2][0];

  // Initialize the error flag to track any interpolation issues.
  int error_flag = BHAHAHA_SUCCESS;

#pragma omp parallel for
  for (int dst_pt = 0; dst_pt < num_dst_pts; dst_pt++) {
    const REAL theta_dst = dst_pts[dst_pt][0]; // Destination point's theta coordinate.
    const REAL phi_dst = dst_pts[dst_pt][1];   // Destination point's phi coordinate.

    // Determine the central grid indices in theta and phi for the interpolation stencil.
    const int idx_center_th = (int)((theta_dst - xxmin_incl_ghosts1) * src_invdxx1 + 0.5);
    const int idx_center_ph = (int)((phi_dst - xxmin_incl_ghosts2) * src_invdxx2 + 0.5);

    {
      // Verify that the interpolation stencil does not exceed grid boundaries.
      if ((idx_center_th - n_interp_ghosts < 0) || (idx_center_th + n_interp_ghosts >= src_Nxx_plus_2NGHOSTS1) ||
          (idx_center_ph - n_interp_ghosts < 0) || (idx_center_ph + n_interp_ghosts >= src_Nxx_plus_2NGHOSTS2)) {
#ifdef DEBUG
        // Provide detailed error messages in debug mode for easier troubleshooting.
        fprintf(stderr, "ERROR: Interpolation stencil exceeds grid boundaries. %d %d %d %d\n", (idx_center_th - n_interp_ghosts < 0),
                (idx_center_th + n_interp_ghosts >= src_Nxx_plus_2NGHOSTS1), (idx_center_ph - n_interp_ghosts < 0),
                (idx_center_ph + n_interp_ghosts >= src_Nxx_plus_2NGHOSTS2));
        fprintf(stderr, "(th_dst, ph_dst) = (%.6f, %.6f). (dst_pt) = %d.\n", theta_dst, phi_dst, dst_pt);
        fprintf(stderr, "Grid bounds along theta direction: [0, %d], stencil indices: [%d, %d]\n", src_Nxx_plus_2NGHOSTS1 - 1,
                idx_center_th - n_interp_ghosts, idx_center_th + n_interp_ghosts);
        fprintf(stderr, "Grid bounds along phi direction: [0, %d], stencil indices: [%d, %d]\n", src_Nxx_plus_2NGHOSTS2 - 1,
                idx_center_ph - n_interp_ghosts, idx_center_ph + n_interp_ghosts);
        fprintf(stderr, "Ensure that the destination point is within grid bounds or adjust the interpolation stencil.\n");
#endif // DEBUG
#pragma omp critical
        {
          error_flag = INTERP2D_GENERAL_HORIZON_OUT_OF_BOUNDS; // Set error flag if stencil is out of bounds.
        }
        continue; // Skip interpolation for this destination point to prevent invalid memory access.
      } // END IF: Check stencil boundaries.

      // Additional sanity checks to ensure central index is correctly positioned.
#ifdef DEBUG
      const REAL TOLERANCE = 1e-13; // Tolerance to account for floating-point precision.
      if (fabs(src_r_theta_phi[1][idx_center_th] - theta_dst) > src_dxx1 * (0.5 + TOLERANCE)) {
        fprintf(stderr, "ERROR: theta center index too far from destination point! %.15e > %.15e\n",
                fabs(src_r_theta_phi[1][idx_center_th] - theta_dst), src_dxx1 * (0.5 + TOLERANCE));
      }
      if (fabs(src_r_theta_phi[2][idx_center_ph] - phi_dst) > src_dxx2 * (0.5 + TOLERANCE)) {
        fprintf(stderr, "ERROR: phi center index too far from destination point! %.15e > %.15e\n", fabs(src_r_theta_phi[2][idx_center_ph] - phi_dst),
                src_dxx2 * (0.5 + TOLERANCE));
      } // END IF: Central index is properly centered
#endif // DEBUG
    } // END SANITY CHECKS: Ensure stencil is valid and central index is correct.

    // Calculate the starting indices for the interpolation stencil in theta and phi directions.
    const int base_idx_th = idx_center_th - n_interp_ghosts;
    const int base_idx_ph = idx_center_ph - n_interp_ghosts;

    // Step 1: Precompute differences between destination theta and source grid theta points within the stencil.
    REAL diffs_th[INTERP_ORDER], diffs_ph[INTERP_ORDER];
    compute_diffs_xi(INTERP_ORDER, theta_dst, &src_r_theta_phi[1][base_idx_th], diffs_th);
    compute_diffs_xi(INTERP_ORDER, phi_dst, &src_r_theta_phi[2][base_idx_ph], diffs_ph);

    // Step 2: Precompute combined Lagrange coefficients to reduce computations
    REAL lagrange_basis_coeffs_th[INTERP_ORDER], lagrange_basis_coeffs_ph[INTERP_ORDER];
    compute_lagrange_basis_coeffs_xi(INTERP_ORDER, inv_denom, diffs_th, lagrange_basis_coeffs_th);
    compute_lagrange_basis_coeffs_xi(INTERP_ORDER, inv_denom, diffs_ph, lagrange_basis_coeffs_ph);
    REAL coeff_2d[INTERP_ORDER][INTERP_ORDER];
    for (int iph = 0; iph < INTERP_ORDER; iph++) {
      const REAL coeff_ph_i = lagrange_basis_coeffs_ph[iph];
      for (int ith = 0; ith < INTERP_ORDER; ith++) {
        coeff_2d[iph][ith] = coeff_ph_i * lagrange_basis_coeffs_th[ith];
      } // END LOOP over theta
    } // END LOOP over phi

    // Define a macro to calculate the flattened index for accessing the source grid function.
#define SRC_IDX2(j, k) ((j) + src_Nxx_plus_2NGHOSTS1 * (k))
      // Step 3: Perform the 1D Lagrange interpolation along the radial direction.
    REAL sum = 0.0;

    for (int iph = 0; iph < INTERP_ORDER; iph++) {
      const int idx_ph = base_idx_ph + iph;
      const int base_offset = base_idx_th + src_Nxx_plus_2NGHOSTS1 * idx_ph;

      sum += sum_lagrange_x0_simd(INTERP_ORDER, &src_gf[base_offset], &coeff_2d[iph][0]);
    } // END LOOP phi direction

    // Store the interpolated value for this grid function and destination point.
    dst_data[dst_pt] = sum * src_invdxx12_INTERP_ORDERm1;

  } // END LOOP: Interpolate all destination points.

  return error_flag; // Return the status of the interpolation process.
"""
    postfunc = r"""
#pragma GCC reset_options // Reset compiler optimizations after the function.

#ifdef STANDALONE

#include <omp.h>

#define NUM_INTERP_GFS 2
#define NUM_RESOLUTIONS 3
#define NUM_DST_PTS 40000000

// Analytic functions.
static inline REAL analytic_function1(REAL x0, REAL x1) { return sin(x0) * cos(x1); }
static inline REAL analytic_function2(REAL x0, REAL x1) { return cos(x0) * sin(x1); }

/**
 * @brief Initializes the 1D coordinate arrays for a 2D uniform source grid.
 *
 * This function calculates the grid spacing (dx) for each dimension and allocates memory for
 * and populates the 1D coordinate arrays. The coordinate arrays include ghost zones. The first
 * coordinate pointer is unused and set to NULL to match the expected 3-pointer array format.
 *
 * @param n_interp_ghosts The number of ghost zones on each side for interpolation.
 * @param N_x0 The number of interior grid points in the x0-dimension.
 * @param N_x1 The number of interior grid points in the x1-dimension.
 * @param[out] src_x0x1 An array of 3 pointers. The function allocates memory for indices 1 and 2. Index 0 is unused.
 * @param[out] src_dxx0 Pointer to store the calculated grid spacing in the x0-dimension.
 * @param[out] src_dxx1 Pointer to store the calculated grid spacing in the x1-dimension.
 * @param src_Nxx_plus_2NGHOSTS0 The total number of points in the x0-dimension, including ghost zones.
 * @param src_Nxx_plus_2NGHOSTS1 The total number of points in the x1-dimension, including ghost zones.
 * @return 0 on success, -1 on memory allocation failure.
 */
int initialize_coordinates(const int n_interp_ghosts, const int N_x0, const int N_x1, REAL *src_x0x1[3], REAL *src_dxx0, REAL *src_dxx1,
                           const int src_Nxx_plus_2NGHOSTS0, const int src_Nxx_plus_2NGHOSTS1) {
  *src_dxx0 = (M_PI) / N_x0;
  *src_dxx1 = (2.0 * M_PI) / N_x1;
  // Index 0 is unused, mimicking the original 2D code's r-theta-phi structure.
  src_x0x1[0] = NULL;
  src_x0x1[1] = (REAL *)malloc(sizeof(REAL) * src_Nxx_plus_2NGHOSTS0);
  src_x0x1[2] = (REAL *)malloc(sizeof(REAL) * src_Nxx_plus_2NGHOSTS1);
  if (!src_x0x1[1] || !src_x0x1[2]) {
    free(src_x0x1[1]);
    free(src_x0x1[2]);
    src_x0x1[1] = src_x0x1[2] = NULL;
    return -1;
  } // END IF: check for allocation failure.
  for (int i = 0; i < src_Nxx_plus_2NGHOSTS0; i++)
    src_x0x1[1][i] = (i - n_interp_ghosts) * (*src_dxx0);
  // END LOOP: initialize x0 coordinates.
  for (int i = 0; i < src_Nxx_plus_2NGHOSTS1; i++)
    src_x0x1[2][i] = (i - n_interp_ghosts) * (*src_dxx1);
  // END LOOP: initialize x1 coordinates.
  return 0;
} // END FUNCTION initialize_coordinates()

/**
 * @brief Populates a 2D source grid function with values from a given analytic function.
 *
 * This function iterates through all points of a 2D grid (including ghost zones)
 * and computes the value of the grid function at each point using the provided
 * analytic function pointer.
 *
 * @param src_Nxx_plus_2NGHOSTS0 The total number of points in the x0-dimension.
 * @param src_Nxx_plus_2NGHOSTS1 The total number of points in the x1-dimension.
 * @param src_x0x1 The pre-initialized 1D coordinate arrays (as an array of 3 pointers).
 * @param[out] src_gf The 2D source grid function data array to be populated.
 * @param func A function pointer to the analytic function used to compute the values.
 */
void initialize_src_gf(const int src_Nxx_plus_2NGHOSTS0, const int src_Nxx_plus_2NGHOSTS1, REAL *src_x0x1[3], REAL *src_gf,
                       REAL (*func)(REAL, REAL)) {
  for (int i1 = 0; i1 < src_Nxx_plus_2NGHOSTS1; i1++) {
    for (int i0 = 0; i0 < src_Nxx_plus_2NGHOSTS0; i0++) {
      const int idx = i0 + src_Nxx_plus_2NGHOSTS0 * i1;
      src_gf[idx] = func(src_x0x1[1][i0], src_x0x1[2][i1]);
    } // END LOOP: over i0.
  } // END LOOP: over i1.
} // END FUNCTION initialize_src_gf()

/**
 * @brief Main driver for testing 2D Lagrange interpolation.
 *
 * This program tests a 2D Lagrange interpolation routine by performing the following steps:
 * 1. Sets up source grids at multiple resolutions.
 * 2. Populates the source grids with data from known analytic functions.
 * 3. Generates a set of random destination points within the source grid domain.
 * 4. Calculates the exact function values at these destination points.
 * 5. Calls the interpolation routine to compute interpolated values at the destination points.
 * 6. Measures the L2 norm of the error between the interpolated and exact values.
 * 7. Calculates and prints the observed order of convergence to verify the accuracy of the interpolator.
 * 8. Prints performance benchmarks.
 *
 * @return EXIT_SUCCESS on successful completion, EXIT_FAILURE otherwise.
 */
int main() {
  int return_code = EXIT_SUCCESS;
  REAL(*dst_pts)[2] = NULL;
  REAL *f_exact[NUM_INTERP_GFS] = {NULL};
  // Changed to an array of 3 pointers to match the interpolator's function signature.
  REAL *src_x0x1[3] = {NULL, NULL, NULL};
  REAL *src_gf[NUM_INTERP_GFS] = {NULL};
  REAL *dst_data[NUM_INTERP_GFS] = {NULL};

  const int n_interp_ghosts = 3;
  const int INTERP_ORDER = 2 * n_interp_ghosts + 1;

  int N_x0_arr[NUM_RESOLUTIONS] = {16, 32, 64};
  int N_x1_arr[NUM_RESOLUTIONS] = {16, 32, 64};
  REAL h_arr[NUM_RESOLUTIONS];
  REAL error_L2_norm[NUM_INTERP_GFS][NUM_RESOLUTIONS];

  dst_pts = (REAL(*)[2])malloc(sizeof(REAL) * NUM_DST_PTS * 2);
  if (!dst_pts) {
    fprintf(stderr, "malloc failed for dst_pts.\n");
    return_code = EXIT_FAILURE;
    goto cleanup;
  } // END IF: check malloc for dst_pts.

  for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
    f_exact[gf] = (REAL *)malloc(sizeof(REAL) * NUM_DST_PTS);
    if (!f_exact[gf]) {
      fprintf(stderr, "malloc failed for f_exact.\n");
      return_code = EXIT_FAILURE;
      goto cleanup;
    } // END IF: check malloc for f_exact.
  } // END LOOP: allocate exact function value arrays.

  for (int res = 0; res < NUM_RESOLUTIONS; res++) {
    // --- Main Loop Over Resolutions ---
    int N_x0 = N_x0_arr[res];
    int N_x1 = N_x1_arr[res];
    h_arr[res] = (N_x0 > 0) ? ((REAL)(M_PI) / N_x0) : 0.0;
    int src_Nxx_plus_2NGHOSTS0 = N_x0 + 2 * n_interp_ghosts;
    int src_Nxx_plus_2NGHOSTS1 = N_x1 + 2 * n_interp_ghosts;
    REAL src_dxx0_val, src_dxx1_val;
    if (initialize_coordinates(n_interp_ghosts, N_x0, N_x1, src_x0x1, &src_dxx0_val, &src_dxx1_val, src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1) !=
        0) {
      fprintf(stderr, "malloc failed for coordinates.\n");
      return_code = EXIT_FAILURE;
      goto cleanup;
    } // END IF: initialize coordinates.
    REAL x0_min_safe = src_x0x1[1][n_interp_ghosts] + 1e-6;
    REAL x0_max_safe = src_x0x1[1][src_Nxx_plus_2NGHOSTS0 - n_interp_ghosts - 1] - 1e-6;
    REAL x1_min_safe = src_x0x1[2][n_interp_ghosts] + 1e-6;
    REAL x1_max_safe = src_x0x1[2][src_Nxx_plus_2NGHOSTS1 - n_interp_ghosts - 1] - 1e-6;
    srand(42 + res);
    for (int i = 0; i < NUM_DST_PTS; i++) {
      dst_pts[i][0] = x0_min_safe + ((REAL)rand() / RAND_MAX) * (x0_max_safe - x0_min_safe);
      dst_pts[i][1] = x1_min_safe + ((REAL)rand() / RAND_MAX) * (x1_max_safe - x1_min_safe);
      f_exact[0][i] = analytic_function1(dst_pts[i][0], dst_pts[i][1]);
      f_exact[1][i] = analytic_function2(dst_pts[i][0], dst_pts[i][1]);
    } // END LOOP: initialize destination points and exact values.
    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      const size_t size = (size_t)src_Nxx_plus_2NGHOSTS0 * src_Nxx_plus_2NGHOSTS1;
      src_gf[gf] = (REAL *)malloc(sizeof(REAL) * size);
      if (!src_gf[gf]) {
        fprintf(stderr, "malloc failed for src_gf.\n");
        return_code = EXIT_FAILURE;
        goto cleanup;
      } // END IF: check malloc for src_gf.
      dst_data[gf] = (REAL *)malloc(sizeof(REAL) * NUM_DST_PTS);
      if (!dst_data[gf]) {
        fprintf(stderr, "malloc failed for dst_data.\n");
        return_code = EXIT_FAILURE;
        goto cleanup;
      } // END IF: check malloc for dst_data.
    } // END LOOP: allocate source grid and destination data arrays.
    initialize_src_gf(src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_x0x1, src_gf[0], analytic_function1);
    initialize_src_gf(src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_x0x1, src_gf[1], analytic_function2);

#ifdef _OPENMP
    double start_time = omp_get_wtime();
#endif
    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      int error_code =
          bah_interpolation_2d_general__uniform_src_grid(n_interp_ghosts, src_dxx0_val, src_dxx1_val, src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1,
                                                         src_x0x1, src_gf[gf], NUM_DST_PTS, dst_pts, dst_data[gf]);
      if (error_code != BHAHAHA_SUCCESS) {
        fprintf(stderr, "Interpolation error code: %d for GF %d\n", error_code, gf + 1);
        return_code = error_code;
        goto cleanup;
      } // END IF: check for interpolation error.
    } // END LOOP: over grid functions for interpolation.
#ifdef _OPENMP
    double elapsed_time = omp_get_wtime() - start_time;
#endif

    printf("\n--- Benchmarking for Resolution %d (%dx%d) ---\n", res, N_x0, N_x1);
#ifdef _OPENMP
    printf("Interpolated %d points for %d GFs in %.4f seconds.\n", NUM_DST_PTS, NUM_INTERP_GFS, elapsed_time);
    printf("Performance: %.4f million points per second.\n", (double)(NUM_DST_PTS * NUM_INTERP_GFS) / elapsed_time / 1e6);
#endif
    printf("--------------------------------------------------\n");

    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      REAL error_sum = 0.0;
      for (int i = 0; i < NUM_DST_PTS; i++) {
        REAL error = dst_data[gf][i] - f_exact[gf][i];
        error_sum += error * error;
      } // END LOOP: calculate squared error sum.
      error_L2_norm[gf][res] = sqrt(error_sum / NUM_DST_PTS);
      printf("Resolution %d: N_x0=%d, h=%.5e, GF %d, L2 error=%.5e\n", res, N_x0, h_arr[res], gf + 1, error_L2_norm[gf][res]);
    } // END LOOP: calculate L2 error norm for each grid function.
    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      free(src_gf[gf]);
      src_gf[gf] = NULL;
      free(dst_data[gf]);
      dst_data[gf] = NULL;
    } // END LOOP: free data for current resolution.
    // Changed loop to free coordinate arrays up to index 2
    for (int dim = 1; dim < 3; dim++) {
      free(src_x0x1[dim]);
      src_x0x1[dim] = NULL;
    } // END LOOP: free coordinate arrays for current resolution.
  } // END LOOP: over resolutions.

  printf("\n--- Convergence Results ---\n");
  for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
    for (int res = 1; res < NUM_RESOLUTIONS; res++) {
      REAL observed_order = log(error_L2_norm[gf][res - 1] / error_L2_norm[gf][res]) / log(h_arr[res - 1] / h_arr[res]);
      printf("Observed order of convergence for GF %d between res %d and %d: %.2f\n", gf + 1, res - 1, res, observed_order);
    } // END LOOP: calculate observed convergence order.
    printf("Expected order of convergence for GF %d: %d\n", gf + 1, INTERP_ORDER);
  } // END LOOP: print convergence results for each grid function.

cleanup:
  if (return_code == EXIT_FAILURE)
    printf("\nAn error occurred. Cleaning up...\n");
  else
    printf("\nProgram finished successfully. Cleaning up...\n");
  // END IF: print final status message.
  for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
    free(src_gf[gf]);
    free(dst_data[gf]);
  } // END LOOP: cleanup src_gf and dst_data.
  // Changed loop to clean up all three pointers
  for (int dim = 0; dim < 3; dim++) {
    free(src_x0x1[dim]);
  } // END LOOP: cleanup src_x0x1.
  for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
    free(f_exact[gf]);
  } // END LOOP: cleanup f_exact.
  free(dst_pts);
  return return_code;

} // END FUNCTION main()

#endif // STANDALONE
"""
    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        postfunc=postfunc,
    )
    return pcg.NRPyEnv()
