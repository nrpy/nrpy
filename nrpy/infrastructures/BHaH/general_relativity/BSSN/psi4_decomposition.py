"""
Generate Psi4 decomposition on spherical like grids using spin-weighted spherical harmonics.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

import nrpy.c_function as cfc
import nrpy.params as par
import numpy as np
from nrpy.infrastructures.BHaH import BHaH_defines_h, griddata_commondata
from nrpy.infrastructures.superB.CurviBoundaryConditions import (
    register_CFunction_apply_bcs_inner_only_specific_gfs,
    register_CFunction_apply_bcs_outerextrap_and_inner_specific_gfs,
)


par.register_CodeParameter(
    "int", __name__, "num_theta_points_on_shell_for_psi4_interp", 32, commondata=True
)

# Define governing parameters for extraction radii (R_exts)
R_ext_min, R_ext_max, num_R_exts = 36.923, 150.0, 8

# Calculate the list of R_exts values.
# The distribution is linear in 1/R_ext, ensuring denser sampling at smaller radii.
# Reciprocal values span 1/R_ext_max to 1/R_ext_min for R_exts to span R_ext_min to R_ext_max.
reciprocal_R_ext_values = np.linspace(1.0 / R_ext_max, 1.0 / R_ext_min, num_R_exts)
# Convert back to R_ext, sort, and round for consistent representation.
list_of_R_exts_values = sorted([round(1.0 / val, 3) for val in reciprocal_R_ext_values])

# Register NUM_OF_R_EXTS: the number of extraction radii.
par.register_CodeParameter(
    cparam_type="#define",
    module=__name__,
    name="NUM_OF_R_EXTS",
    defaultvalue=num_R_exts,
    description=f"Number of extraction radii ({num_R_exts})."
)

# Register list_of_R_exts: the array containing extraction radii values.
par.register_CodeParameter(
    cparam_type=f"REAL[{num_R_exts}]",
    module=__name__,
    name="list_of_R_exts",
    defaultvalue=list_of_R_exts_values,
    description=f"List of radii at which psi4 is extracted.",
    add_to_parfile=True,
    commondata=True,
    add_to_set_CodeParameters_h=False,
    assumption="RealPositive"
)

def register_griddata() -> None:
    """Register the diagnostic_struct's contribution to the griddata_struct."""
    griddata_commondata.register_griddata_commondata(
        __name__,
        "diagnostic_struct diagnosticstruct",
        "data needed to do psi4 decomposition in cylindrical-like coordinates",
    )

# Register the diagnostic_struct to the BHaH_defines.h
def register_BHaH_defines_h() -> None:
    """Register the diagnostic_struct's contribution to the BHaH_defines.h file."""
    BHd_str = r"""
    // for psi4 decomposition
    typedef struct __diagnostic_struct__ {

      int *restrict N_shell_pts_grid; // of shape int [NUM_OF_R_EXTS]
      REAL **restrict xx_radial_like_shell_grid; // of shape [NUM_OF_R_EXTS][N_shell_pts_grid]
      REAL **restrict xx_theta_like_shell_grid; // of shape [NUM_OF_R_EXTS][N_shell_pts_grid]
      int *restrict N_theta_shell_grid; // of shape int [NUM_OF_R_EXTS]
      REAL **restrict theta_shell_grid; // of shape [NUM_OF_R_EXTS][N_theta_shell_grid]
      REAL dtheta;
    } diagnostic_struct;
    """

    BHaH_defines_h.register_BHaH_defines(
        __name__,
        BHd_str,
    )


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

    setup_prefunc = rf"""

  // IDX2GENERAL: Map 2D indices (i,j) to a 1D index using stride Ni
  #define IDX2GENERAL(i, j, Ni) ((i) + (Ni) * (j))
  // REVERSE_IDX2GENERAL: Recover (i,j) from a 1D index and stride Ni
  #define REVERSE_IDX2GENERAL(index, Ni, i, j) \
  {{                                            \
      j = (index) / (Ni);                      \
      i = (index) % (Ni);                      \
  }}

  static void psi4_diagnostics_set_up(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],
                                                     diagnostic_struct *restrict diagnosticstruct) {{

  const REAL *restrict list_of_R_exts = commondata->list_of_R_exts;

  //Thin shells with N_theta and N_phi number of points are constructed at each extraction radius R_ext
  const int N_theta = commondata->num_theta_points_on_shell_for_psi4_interp;

  // phi values of shell need to be exactly the same as phi values of grid, all coordinate systems supported have a phi coordinate
  const int N_phi = params->Nxx{phi_index};

  // Set up uniform 2d grid in theta and phi at R_ext (2d shells at different R_ext)
  const REAL PI = M_PI;
  const REAL theta_min = 0.0;
  REAL theta_max = PI;
  REAL phi_min = -PI;
  REAL phi_max = PI;
  const int N_tot_shell = N_theta * N_phi;
  REAL dtheta = (theta_max - theta_min) / N_theta;
  diagnosticstruct->dtheta = dtheta;
  const REAL dphi = (phi_max - phi_min) / N_phi;
  REAL ***restrict xx_shell_sph = (REAL * **restrict)malloc(NUM_OF_R_EXTS * sizeof(REAL **));
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {{
    xx_shell_sph[which_R_ext] = (REAL * *restrict)malloc(2 * sizeof(REAL *));
    xx_shell_sph[which_R_ext][0] = (REAL *restrict)malloc(N_theta * sizeof(REAL));
    xx_shell_sph[which_R_ext][1] = (REAL *restrict)malloc(N_phi * sizeof(REAL));
  }}
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {{
    for (int j = 0; j < N_theta; j++) {{
      xx_shell_sph[which_R_ext][0][j] = theta_min + ((REAL)j + 0.5) * dtheta;
    }}
    for (int j = 0; j < N_phi; j++) {{
      xx_shell_sph[which_R_ext][1][j] = phi_min + ((REAL)j + 0.5) * dphi;
    }}
  }}

  // Convert points on shells to Cartesian coordinates
  REAL ***restrict xx_shell_Cart = (REAL * **restrict)malloc(NUM_OF_R_EXTS * sizeof(REAL **));
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {{
    xx_shell_Cart[which_R_ext] = (REAL * *restrict)malloc(3 * sizeof(REAL *));
    xx_shell_Cart[which_R_ext][0] = (REAL *restrict)malloc(N_tot_shell * sizeof(REAL));
    xx_shell_Cart[which_R_ext][1] = (REAL *restrict)malloc(N_tot_shell * sizeof(REAL));
    xx_shell_Cart[which_R_ext][2] = (REAL *restrict)malloc(N_tot_shell * sizeof(REAL));
  }}
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {{
    REAL r = list_of_R_exts[which_R_ext];
    for (int i_ph = 0; i_ph < N_phi; i_ph++) {{
      REAL phi = xx_shell_sph[which_R_ext][1][i_ph];
      for (int i_th = 0; i_th < N_theta; i_th++) {{
        REAL theta = xx_shell_sph[which_R_ext][0][i_th];
        REAL x = r * sin(theta) * cos(phi);
        REAL y = r * sin(theta) * sin(phi);
        REAL z = r * cos(theta);
        const int idx2 = IDX2GENERAL(i_th, i_ph, N_theta);
        xx_shell_Cart[which_R_ext][0][idx2] = x;
        xx_shell_Cart[which_R_ext][1][idx2] = y;
        xx_shell_Cart[which_R_ext][2][idx2] = z;
      }}
    }}
  }}

  // For each pt on each spherical shell at R_ext find if pt lies within the grid
  //(it should for bhah single grids, but not for multipatch...)
  // set N_shell_pts_grid and xx_shell_grid in diagnostic struct
  diagnosticstruct->N_shell_pts_grid = (int *restrict)malloc(sizeof(int) * NUM_OF_R_EXTS);
  diagnosticstruct->xx_radial_like_shell_grid = (REAL **restrict)malloc(NUM_OF_R_EXTS * sizeof(REAL *));
  diagnosticstruct->xx_theta_like_shell_grid = (REAL **restrict)malloc(NUM_OF_R_EXTS * sizeof(REAL *));
  diagnosticstruct->N_theta_shell_grid = (int *restrict)malloc(sizeof(int) * NUM_OF_R_EXTS);
  diagnosticstruct->theta_shell_grid = (REAL **restrict)malloc(NUM_OF_R_EXTS * sizeof(REAL *));

  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {{
    // Count number of pts that lie within the extent of this grid
    int count_pt_on_grid = 0;
    bool b_theta_grid[N_theta];
    bool b_phi_grid[N_phi];
    for (int i = 0; i < N_theta; i++) {{
        b_theta_grid[i] = false;
    }}
    for (int j = 0; j < N_phi; j++) {{
        b_phi_grid[j] = false;
    }}
    for (int i_ph = 0; i_ph < N_phi; i_ph++) {{
      for (int i_th = 0; i_th < N_theta; i_th++) {{
        const int idx2 = IDX2GENERAL(i_th, i_ph, N_theta);
        const REAL xCart_pt_on_shell[3] = {{xx_shell_Cart[which_R_ext][0][idx2], xx_shell_Cart[which_R_ext][1][idx2],
                                           xx_shell_Cart[which_R_ext][2][idx2]}};
        int Cart_to_i0i1i2[3];
        REAL closest_xx[3];
        Cart_to_xx_and_nearest_i0i1i2(params, xCart_pt_on_shell, closest_xx, Cart_to_i0i1i2); // convert from Cart to xx
        if ((params->xxmin0 <= closest_xx[0] && closest_xx[0] <= params->xxmax0) &&
            (params->xxmin1 <= closest_xx[1] && closest_xx[1] <= params->xxmax1) &&
            (params->xxmin2 <= closest_xx[2] && closest_xx[2] <= params->xxmax2)) {{
          count_pt_on_grid++;
          b_theta_grid[i_th] = true;
          b_phi_grid[i_ph] = true;
        }}
      }}
    }}
    // Count the number of unique theta values on grid
    int count_theta_grid = 0;
    for (int i = 0; i < N_theta; i++) {{
      if (b_theta_grid[i]) {{
        count_theta_grid++;
      }}
    }}
    // Count the number of unique phi values on grid (should be the same as Nxx1)
    int count_phi_grid = 0;
    for (int i = 0; i < N_phi; i++) {{
      if (b_phi_grid[i]) {{
        count_phi_grid++;
      }}
    }}
    // Set values in diagnosticstruct
    diagnosticstruct->N_shell_pts_grid[which_R_ext] = count_pt_on_grid;
    diagnosticstruct->N_theta_shell_grid[which_R_ext] = count_theta_grid;

    // Allocate memory after counting
    diagnosticstruct->xx_radial_like_shell_grid[which_R_ext] = (REAL *restrict)malloc(count_pt_on_grid * sizeof(REAL));
    diagnosticstruct->xx_theta_like_shell_grid[which_R_ext] = (REAL *restrict)malloc(count_pt_on_grid * sizeof(REAL));
    diagnosticstruct->theta_shell_grid[which_R_ext] = (REAL *restrict)malloc(count_theta_grid * sizeof(REAL));

    // Now set them
    int which_pt_on_grid = 0;
    for (int i_ph = 0; i_ph < N_phi; i_ph++) {{
      for (int i_th = 0; i_th < N_theta; i_th++) {{
        const int idx2 = IDX2GENERAL(i_th, i_ph, N_theta);
        const REAL xCart_pt_on_shell[3] = {{xx_shell_Cart[which_R_ext][0][idx2], xx_shell_Cart[which_R_ext][1][idx2],
                                           xx_shell_Cart[which_R_ext][2][idx2]}};
        int Cart_to_i0i1i2[3];
        REAL closest_xx[3];
        Cart_to_xx_and_nearest_i0i1i2(params, xCart_pt_on_shell, closest_xx, Cart_to_i0i1i2);
        if ((params->xxmin0 <= closest_xx[0] && closest_xx[0] <= params->xxmax0) &&
            (params->xxmin1 <= closest_xx[1] && closest_xx[1] <= params->xxmax1) &&
            (params->xxmin2 <= closest_xx[2] && closest_xx[2] <= params->xxmax2)) {{

          diagnosticstruct->xx_radial_like_shell_grid[which_R_ext][which_pt_on_grid] = closest_xx[{radial_like_index}];
          diagnosticstruct->xx_theta_like_shell_grid[which_R_ext][which_pt_on_grid] = closest_xx[{theta_like_index}];

          // also save theta values
          int i_th_grid;
          MAYBE_UNUSED int i_ph_grid;
          const int N_theta_shell_grid = diagnosticstruct->N_theta_shell_grid[which_R_ext];
          REVERSE_IDX2GENERAL(which_pt_on_grid, N_theta_shell_grid, i_th_grid, i_ph_grid);
          diagnosticstruct->theta_shell_grid[which_R_ext][i_th_grid] = xx_shell_sph[which_R_ext][0][i_th];
          which_pt_on_grid++;
        }}
      }}
    }}
  }} // end loop over all R_ext

  int sum = 0;
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {{
    if (diagnosticstruct->N_shell_pts_grid[which_R_ext] > 0) {{
      sum += diagnosticstruct->N_shell_pts_grid[which_R_ext];
    }}
  }}

  // Free memory
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {{
    free(xx_shell_sph[which_R_ext][0]);
    free(xx_shell_sph[which_R_ext][1]);
    free(xx_shell_sph[which_R_ext]);
  }}
  free(xx_shell_sph);
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {{
    free(xx_shell_Cart[which_R_ext][0]);
    free(xx_shell_Cart[which_R_ext][1]);
    free(xx_shell_Cart[which_R_ext][2]);
    free(xx_shell_Cart[which_R_ext]);
  }}
  free(xx_shell_Cart);
}}
"""

    decomposition_prefunc = r"""
static void lowlevel_decompose_psi4_into_swm2_modes(const int Nxx_plus_2NGHOSTS1, const REAL dxx1, const REAL dxx2,
                                                    const int swm2sh_maximum_l_mode_to_compute, const REAL curr_time, const REAL R_ext,
                                                    const REAL *restrict th_array, const REAL *restrict sinth_array, const REAL *restrict ph_array,
                                                    const int N_theta, const REAL *restrict psi4r_at_R_ext, const REAL *restrict psi4i_at_R_ext) {
  char filename[100];
  FILE *outpsi4_l_m;
  // Output header at t=0:
  if (curr_time == 0) {
    for (int l = 2; l <= swm2sh_maximum_l_mode_to_compute; l++) {
      sprintf(filename, "Rpsi4_l%d-r%06.1f.txt", l, (double)R_ext);
      outpsi4_l_m = fopen(filename, "w");
      fprintf(outpsi4_l_m, "# column 1: t-R_ext = [retarded time]\n");
      int col = 2;
      for (int m = -l; m <= l; m++) {
        fprintf(outpsi4_l_m, "# column %d: Re(psi4_{l=%d,m=%d}) * R_ext\n", col, l, m);
        col++;
        fprintf(outpsi4_l_m, "# column %d: Im(psi4_{l=%d,m=%d}) * R_ext\n", col, l, m);
        col++;
      }
      fclose(outpsi4_l_m);
    }
  }

  // Output one file per l mode; each column represents a unique complex component of l,m
  for (int l = 2; l <= swm2sh_maximum_l_mode_to_compute; l++) {
    sprintf(filename, "Rpsi4_l%d-r%06.1f.txt", l, (double)R_ext);
    outpsi4_l_m = fopen(filename, "a");
    char oneline[10000];
    sprintf(oneline, "%e", (double)(curr_time - R_ext));
    for (int m = -l; m <= l; m++) {
      // Parallelize the integration loop:
      REAL psi4r_l_m = 0.0;
      REAL psi4i_l_m = 0.0;
#pragma omp parallel for reduction(+ : psi4r_l_m, psi4i_l_m)
      for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1 - 2 * NGHOSTS; i1++) {
        const REAL ph = ph_array[i1];
        for (int i2 = 0; i2 < N_theta; i2++) {
          const REAL th = th_array[i2];
          const REAL sinth = sinth_array[i2];
          // Construct integrand for psi4 spin-weight s=-2 spherical harmonic
          REAL ReY_sm2_l_m, ImY_sm2_l_m;
          spin_weight_minus2_sph_harmonics(l, m, th, ph, &ReY_sm2_l_m, &ImY_sm2_l_m);

          const int idx2d = IDX2GENERAL(i2, i1, N_theta);
          const REAL a = psi4r_at_R_ext[idx2d];
          const REAL b = psi4i_at_R_ext[idx2d];
          const REAL c = ReY_sm2_l_m;
          const REAL d = ImY_sm2_l_m;
          psi4r_l_m += (a * c + b * d) * dxx2 * sinth * dxx1;
          psi4i_l_m += (b * c - a * d) * dxx2 * sinth * dxx1;
        }
      }
      sprintf(oneline + strlen(oneline), " %.15e %.15e", (double)(R_ext * psi4r_l_m), (double)(R_ext * psi4i_l_m));
    }
    fprintf(outpsi4_l_m, "%s\n", oneline);
    fclose(outpsi4_l_m);
  }
}
"""
    prefunc = setup_prefunc + "\n\n" + decomposition_prefunc

    desc = "Decompose psi4 across all l,m modes from l=2 up to and including L_MAX (global variable)."
    name = "psi4_spinweightm2_decomposition"
    params = r"""const commondata_struct *restrict commondata,
    const params_struct *restrict params,
    REAL *restrict diagnostic_output_gfs,
    REAL *restrict xx[3],
    diagnostic_struct *restrict diagnosticstruct"""


    # Register diagnostic_struct's contribution to BHaH_defines.h:
    register_BHaH_defines_h()

    # Register diagnostic_struct's contribution to griddata_struct:
    register_griddata()

    # Register C functions apply_bcs_inner_only_specific_gfs and apply_bcs_outerextrap_and_inner_specific_gfs, needed for 2d interp of psi4
    register_CFunction_apply_bcs_inner_only_specific_gfs()
    register_CFunction_apply_bcs_outerextrap_and_inner_specific_gfs()


    body = r"""
  static bool is_first_call = true;
  if (is_first_call) {
    psi4_diagnostics_set_up(commondata, params, xx, diagnosticstruct);
    is_first_call = false;
  }

  const REAL *restrict list_of_R_exts = commondata->list_of_R_exts;
  int *restrict N_shell_pts_grid = diagnosticstruct->N_shell_pts_grid;
  REAL **restrict xx_radial_like_shell_grid = diagnosticstruct->xx_radial_like_shell_grid;
  REAL **restrict xx_theta_like_shell_grid = diagnosticstruct->xx_theta_like_shell_grid;
  int *restrict N_theta_shell_grid = diagnosticstruct->N_theta_shell_grid;
  REAL **restrict theta_shell_grid = diagnosticstruct->theta_shell_grid;
  const REAL dtheta = diagnosticstruct->dtheta;
"""

    body += rf"""

  const int radial_like_index = {radial_like_index};
  const int Nxx_plus_2NGHOSTS_radial_like_index = Nxx_plus_2NGHOSTS{radial_like_index};
  const REAL dxx_radial_like_index = params->dxx{radial_like_index};

  const int theta_like_index = {theta_like_index};
  const int Nxx_plus_2NGHOSTS_theta_like_index = Nxx_plus_2NGHOSTS{theta_like_index};
  const REAL dxx_theta_like_index = params->dxx{theta_like_index};

  const int phi_index = {phi_index};
  const int Nxx_plus_2NGHOSTS_phi_index = Nxx_plus_2NGHOSTS{phi_index};
  const REAL dxx_phi_index = params->dxx{phi_index};
"""

    body += r"""

  // set parameters for interpolation
  // src grid is the whole grid in radial-like and theta-like coordinates
  const int N_interp_GHOSTS = NGHOSTS;
  const REAL src_dxx0 = dxx_radial_like_index;
  const REAL src_dxx1 = dxx_theta_like_index;
  const int src_Nxx_plus_2NGHOSTS0 = Nxx_plus_2NGHOSTS_radial_like_index;
  const int src_Nxx_plus_2NGHOSTS1 = Nxx_plus_2NGHOSTS_theta_like_index;
  REAL *restrict src_x0x1[3];
  src_x0x1[1] = (REAL *restrict)malloc(sizeof(REAL) * src_Nxx_plus_2NGHOSTS0);
  src_x0x1[2] = (REAL *restrict)malloc(sizeof(REAL) * src_Nxx_plus_2NGHOSTS1);
  for (int j = 0; j < src_Nxx_plus_2NGHOSTS0; j++)
    src_x0x1[1][j] = xx[radial_like_index][j];
  for (int j = 0; j < src_Nxx_plus_2NGHOSTS1; j++)
    src_x0x1[2][j] = xx[theta_like_index][j];
  const int total_size = src_Nxx_plus_2NGHOSTS0 * src_Nxx_plus_2NGHOSTS1;
  REAL *restrict src_gf_psi4r = (REAL *)malloc(sizeof(REAL) * total_size); // with shape [src_Nxx_plus_2NGHOSTS0 * src_Nxx_plus_2NGHOSTS1]
  REAL *restrict src_gf_psi4i = (REAL *)malloc(sizeof(REAL) * total_size); // with shape [src_Nxx_plus_2NGHOSTS0 * src_Nxx_plus_2NGHOSTS1]

  // Step 2: Loop over all extraction indices:
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {
    // Step 2.a: Set the extraction radius R_ext based on the radial index R_ext_idx
    const REAL R_ext = list_of_R_exts[which_R_ext];
    const int num_dst_pts = N_shell_pts_grid[which_R_ext];
    if (num_dst_pts > 0) {

      REAL(*dst_pts)[2] = (REAL(*)[2])malloc(sizeof(REAL) * num_dst_pts * 2);
      // Destination points
      for (int i = 0; i < num_dst_pts; i++) {
        dst_pts[i][0] = xx_radial_like_shell_grid[which_R_ext][i];
        dst_pts[i][1] = xx_theta_like_shell_grid[which_R_ext][i];
      }
      REAL *dst_data_psi4r = (REAL *)malloc(sizeof(REAL) * num_dst_pts);
      REAL *dst_data_psi4i = (REAL *)malloc(sizeof(REAL) * num_dst_pts);

      // Step 1: Allocate memory for 2D arrays used to store psi4, theta, sin(theta), and phi.
      int sizeof_2Darray = sizeof(REAL) * (Nxx_plus_2NGHOSTS_phi_index - 2 * NGHOSTS) * N_theta_shell_grid[which_R_ext];
      REAL *restrict psi4r_at_R_ext = (REAL *restrict)malloc(sizeof_2Darray);
      REAL *restrict psi4i_at_R_ext = (REAL *restrict)malloc(sizeof_2Darray);
      //         ... also store theta, sin(theta), and phi to corresponding 1D arrays.
      REAL *restrict sinth_array = (REAL *restrict)malloc(sizeof(REAL) * N_theta_shell_grid[which_R_ext]);
      REAL *restrict th_array = (REAL *restrict)malloc(sizeof(REAL) * N_theta_shell_grid[which_R_ext]);
      REAL *restrict ph_array = (REAL *restrict)malloc(sizeof(REAL) * (Nxx_plus_2NGHOSTS_phi_index - 2 * NGHOSTS));

// Step 2.b: Compute psi_4 at this extraction radius and store to a local 2D array.
      #pragma omp parallel for
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS_phi_index - NGHOSTS; i1++) {
        ph_array[i1 - NGHOSTS] = xx[phi_index][i1];
// Initialize src_gf
#pragma omp parallel for
        for (int j = 0; j < src_Nxx_plus_2NGHOSTS1; j++) {
          for (int i = 0; i < src_Nxx_plus_2NGHOSTS0; i++) {
            int idx = i + src_Nxx_plus_2NGHOSTS0 * j;

            int ii012[3];
            ii012[radial_like_index] = i;
            ii012[theta_like_index] = j;
            ii012[phi_index] = i1;

            src_gf_psi4r[idx] = diagnostic_output_gfs[IDX4(PSI4_REGF, ii012[0], ii012[1], ii012[2])];
            src_gf_psi4i[idx] = diagnostic_output_gfs[IDX4(PSI4_IMGF, ii012[0], ii012[1], ii012[2])];
          }
        }
        // Call the interpolation function
        int error_code1 =
            interpolation_2d_general__uniform_src_grid(N_interp_GHOSTS, src_dxx0, src_dxx1, src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_x0x1,
                                                       src_gf_psi4r, num_dst_pts, dst_pts, dst_data_psi4r);

        if (error_code1 > 0) {
          printf("Interpolation error code: %d\n", error_code1);
        }

        int error_code2 =
            interpolation_2d_general__uniform_src_grid(N_interp_GHOSTS, src_dxx0, src_dxx1, src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_x0x1,
                                                       src_gf_psi4i, num_dst_pts, dst_pts, dst_data_psi4i);
        if (error_code2 > 0) {
          printf("Interpolation error code: %d\n", error_code2);
        }

        for (int i2 = 0; i2 < N_theta_shell_grid[which_R_ext]; i2++) {
          th_array[i2] = theta_shell_grid[which_R_ext][i2];
          sinth_array[i2] = sin(th_array[i2]);
          // Store result to "2D" array (actually 1D array with 2D storage):
          const int idx2d = IDX2GENERAL(i2, (i1 - NGHOSTS), N_theta_shell_grid[which_R_ext]);
          psi4r_at_R_ext[idx2d] = dst_data_psi4r[IDX2GENERAL(i2, i1 - NGHOSTS, N_theta_shell_grid[which_R_ext])];
          psi4i_at_R_ext[idx2d] = dst_data_psi4i[IDX2GENERAL(i2, i1 - NGHOSTS, N_theta_shell_grid[which_R_ext])];
        }
      } // end loop over all phi values of grid

      // Step 3: Perform integrations across all l,m modes from l=2 up to and including L_MAX (global variable):
      lowlevel_decompose_psi4_into_swm2_modes(Nxx_plus_2NGHOSTS_phi_index, dxx_phi_index, dtheta, swm2sh_maximum_l_mode_to_compute, commondata->time,
                                              R_ext, th_array, sinth_array, ph_array, N_theta_shell_grid[which_R_ext], psi4r_at_R_ext,
                                              psi4i_at_R_ext);

      // Step 4: Free all allocated memory:
      free(psi4r_at_R_ext);
      free(psi4i_at_R_ext);
      free(sinth_array);
      free(th_array);
      free(ph_array);
      free(dst_pts);
      free(dst_data_psi4r);
      free(dst_data_psi4i);
    } // end if num dst pts > 0
  } // end for loop over all R_ext
  free(src_gf_psi4r);
  free(src_gf_psi4i);
  free(src_x0x1[1]);
  free(src_x0x1[2]);
"""

    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        prefunc=prefunc,
        desc=desc,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
