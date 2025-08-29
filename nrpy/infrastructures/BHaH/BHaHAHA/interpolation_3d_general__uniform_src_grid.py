"""
Register C function for 3D Lagrange interpolation at arbitrary sets of points.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
from nrpy.helpers.generic import copy_files


def register_CFunction_interpolation_3d_general__uniform_src_grid(
    enable_simd: bool,
    project_dir: str,
    use_cpp: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for general-purpose 3D Lagrange interpolation.

    Even if `enable_simd` is False, `intrinsics/simd_intrinsics.h` is still required.

    :param enable_simd: Whether the rest of the code enables SIMD optimizations, as this code requires simd_intrinsics.h (which includes SIMD-disabled options).
    :param project_dir: Directory of the project, to set the destination for simd_instrinsics.h .
    :param use_cpp: If True, emit a C++-compatible variant: map 'restrict'→'__restrict__' under C++, and switch sized array params to unsized pointer params (src_x0x1x2[], src_gf_ptrs[], dst_data[]).

    :return: None if in registration phase, else the updated NRPy environment.

    DocTests:
    >>> env = register_CFunction_interpolation_3d_general__uniform_src_grid(enable_simd=False, project_dir=".")
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

enum {
  INTERP_SUCCESS,
  INTERP3D_GENERAL_NULL_PTRS,
  INTERP3D_GENERAL_INTERP_ORDER_GT_NXX123,
  INTERP3D_GENERAL_HORIZON_OUT_OF_BOUNDS
} error_codes;

#pragma GCC optimize("unroll-loops")
"""
    if use_cpp:
        prefunc += r"""
#ifndef restrict
#define restrict __restrict__
#endif
"""

    desc = r"""Performs 3D Lagrange interpolation from a set of uniform grid points on the source grid to arbitrary destination points.

This function interpolates scalar grid functions from a source grid to a set of destination points in the x0, x1, and x2 directions,
using Lagrange interpolation of order INTERP_ORDER.

@param N_interp_GHOSTS - Number of ghost zones from the center of source point; interpolation order = 2 * N_interp_GHOSTS + 1.
@param src_dxx0 - Grid spacing in the x0 direction on the source grid.
@param src_dxx1 - Grid spacing in the x1 direction on the source grid.
@param src_dxx2 - Grid spacing in the x2 direction on the source grid.
@param src_Nxx_plus_2NGHOSTS0 - Dimension of the source grid in x0, including ghost zones.
@param src_Nxx_plus_2NGHOSTS1 - Dimension of the source grid in x1, including ghost zones.
@param src_Nxx_plus_2NGHOSTS2 - Dimension of the source grid in x2, including ghost zones.
@param NUM_INTERP_GFS - Number of grid functions to interpolate.
@param src_x0x1x2 - Arrays of coordinate values for x0, x1, and x2 on the source grid.
@param src_gf_ptrs - Array of pointers to source grid functions data.
@param num_dst_pts - Number of destination points.
@param dst_x0x1x2 - Destination points' coordinates (x0, x1, x2).
@param dst_data - Output interpolated data for each grid function at the destination points, of size [NUM_INTERP_GFS][num_dst_pts].

@return - Error code indicating success or type of error encountered.

@note - The function interpolates each grid function separately and stores the results independently.
The source and destination grids are assumed to be uniform in x0, x1, and x2 directions.
The function assumes that the destination grid points are within the range of the source grid.
"""

    cfunc_type = "int"
    name = "interpolation_3d_general__uniform_src_grid"

    # Parameter suffixes: C++ uses unsized pointers ([]); C uses sized arrays ([3], [NUM_INTERP_GFS])
    src_x0x1x2_suffix = "[]" if use_cpp else "[3]"
    interp_gfs_suffix = "[]" if use_cpp else "[NUM_INTERP_GFS]"
    params = f"""
        const int n_interp_ghosts, const REAL src_dxx0, const REAL src_dxx1, const REAL src_dxx2,
        const int src_Nxx_plus_2NGHOSTS0, const int src_Nxx_plus_2NGHOSTS1, const int src_Nxx_plus_2NGHOSTS2,
        const int NUM_INTERP_GFS, REAL *restrict src_x0x1x2{src_x0x1x2_suffix}, const REAL *restrict src_gf_ptrs{interp_gfs_suffix},
        const int num_dst_pts, const REAL dst_x0x1x2[][3], REAL *restrict dst_data{interp_gfs_suffix}"""

    body = r"""
  // Unpack parameters.
  const int NinterpGHOSTS = n_interp_ghosts;
  const int INTERP_ORDER = (2 * NinterpGHOSTS + 1); // Interpolation order (number of points in stencil in each dimension).

  const REAL src_invdxx0 = 1.0 / src_dxx0;
  const REAL src_invdxx1 = 1.0 / src_dxx1;
  const REAL src_invdxx2 = 1.0 / src_dxx2;

  // Compute normalization factor once to avoid repeated expensive pow() operations.
  const REAL src_invdxx012_INTERP_ORDERm1 = pow(src_dxx0 * src_dxx1 * src_dxx2, -(INTERP_ORDER - 1));

  // Check for null pointers in source coordinates and output data.
  if (src_x0x1x2[0] == NULL || src_x0x1x2[1] == NULL || src_x0x1x2[2] == NULL)
    return INTERP3D_GENERAL_NULL_PTRS;
  for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
    if (dst_data[gf] == NULL)
      return INTERP3D_GENERAL_NULL_PTRS;
  }

  // Ensure interpolation order does not exceed grid dimensions.
  if (INTERP_ORDER > src_Nxx_plus_2NGHOSTS0 || INTERP_ORDER > src_Nxx_plus_2NGHOSTS1 || INTERP_ORDER > src_Nxx_plus_2NGHOSTS2)
    return INTERP3D_GENERAL_INTERP_ORDER_GT_NXX123;

  // Precompute inverse denominators for Lagrange interpolation coefficients to optimize performance.
  REAL inv_denom[INTERP_ORDER];
  compute_inv_denom(INTERP_ORDER, inv_denom);

  // Perform interpolation for each destination point (x0, x1, x2).
  const REAL xxmin_incl_ghosts0 = src_x0x1x2[0][0];
  const REAL xxmin_incl_ghosts1 = src_x0x1x2[1][0];
  const REAL xxmin_incl_ghosts2 = src_x0x1x2[2][0];
  int error_flag = INTERP_SUCCESS;

#pragma omp parallel for
  for (int dst_pt = 0; dst_pt < num_dst_pts; dst_pt++) {
    // Extract destination point coordinates.
    const REAL x0_dst = dst_x0x1x2[dst_pt][0];
    const REAL x1_dst = dst_x0x1x2[dst_pt][1];
    const REAL x2_dst = dst_x0x1x2[dst_pt][2];

    // Compute the central index of the interpolation stencil in each dimension.
    int idx_center_x0 = (int)((x0_dst - xxmin_incl_ghosts0) * src_invdxx0 + 0.5);
    int idx_center_x1 = (int)((x1_dst - xxmin_incl_ghosts1) * src_invdxx1 + 0.5);
    int idx_center_x2 = (int)((x2_dst - xxmin_incl_ghosts2) * src_invdxx2 + 0.5);

    // Check if the interpolation stencil goes out of bounds, and adjust indices to prevent memory corruption.
    {
      if ((idx_center_x0 - NinterpGHOSTS < 0) || (idx_center_x0 + NinterpGHOSTS >= src_Nxx_plus_2NGHOSTS0) || (idx_center_x1 - NinterpGHOSTS < 0) ||
          (idx_center_x1 + NinterpGHOSTS >= src_Nxx_plus_2NGHOSTS1) || (idx_center_x2 - NinterpGHOSTS < 0) ||
          (idx_center_x2 + NinterpGHOSTS >= src_Nxx_plus_2NGHOSTS2)) {
#pragma omp critical
        {
          error_flag = INTERP3D_GENERAL_HORIZON_OUT_OF_BOUNDS;
          // If out of bounds, set indices to NinterpGHOSTS to avoid accessing invalid memory.
          idx_center_x0 = idx_center_x1 = idx_center_x2 = NinterpGHOSTS;
        }
        continue; // Skip further work for this iteration
      }
    } // END IF: Set ERROR code & adjust indices if stencil is out of bounds.

    // Compute base indices for interpolation stencil.
    const int base_idx_x0 = idx_center_x0 - NinterpGHOSTS;
    const int base_idx_x1 = idx_center_x1 - NinterpGHOSTS;
    const int base_idx_x2 = idx_center_x2 - NinterpGHOSTS;

    // Compute differences for Lagrange interpolation.
    REAL diffs_x0[INTERP_ORDER], diffs_x1[INTERP_ORDER], diffs_x2[INTERP_ORDER];
    compute_diffs_xi(INTERP_ORDER, x0_dst, &src_x0x1x2[0][base_idx_x0], diffs_x0);
    compute_diffs_xi(INTERP_ORDER, x1_dst, &src_x0x1x2[1][base_idx_x1], diffs_x1);
    compute_diffs_xi(INTERP_ORDER, x2_dst, &src_x0x1x2[2][base_idx_x2], diffs_x2);

    REAL coeff_x0[INTERP_ORDER], coeff_x1[INTERP_ORDER], coeff_x2[INTERP_ORDER];
    compute_lagrange_basis_coeffs_xi(INTERP_ORDER, inv_denom, diffs_x0, coeff_x0);
    compute_lagrange_basis_coeffs_xi(INTERP_ORDER, inv_denom, diffs_x1, coeff_x1);
    compute_lagrange_basis_coeffs_xi(INTERP_ORDER, inv_denom, diffs_x2, coeff_x2);

    // Compute the combined 3D Lagrange coefficients with reordered indices.
    REAL coeff_3d[INTERP_ORDER][INTERP_ORDER][INTERP_ORDER];
    for (int ix2 = 0; ix2 < INTERP_ORDER; ix2++) {
      const REAL coeff_x2_i = coeff_x2[ix2];
      for (int ix1 = 0; ix1 < INTERP_ORDER; ix1++) {
        const REAL coeff_x1_i = coeff_x1[ix1];
        for (int ix0 = 0; ix0 < INTERP_ORDER; ix0++) {
          coeff_3d[ix2][ix1][ix0] = coeff_x0[ix0] * coeff_x1_i * coeff_x2_i;
        } // END LOOP x0 direction
      } // END LOOP x1 direction
    } // END LOOP x2 direction

#define SRC_IDX3(i0, i1, i2) ((i0) + src_Nxx_plus_2NGHOSTS0 * ((i1) + src_Nxx_plus_2NGHOSTS1 * (i2)))

    // For each grid function, compute the interpolated value.
    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      REAL sum = 0.0;

      for (int ix2 = 0; ix2 < INTERP_ORDER; ix2++) {
        const int idx2 = base_idx_x2 + ix2;
        const int src_Nxx_plus_2NGHOSTS1_times_idx2 = src_Nxx_plus_2NGHOSTS1 * idx2;
        for (int ix1 = 0; ix1 < INTERP_ORDER; ix1++) {
          const int idx1 = base_idx_x1 + ix1;
          const int base_offset = base_idx_x0 + src_Nxx_plus_2NGHOSTS0 * (idx1 + src_Nxx_plus_2NGHOSTS1_times_idx2);

          sum += sum_lagrange_x0_simd(INTERP_ORDER, &src_gf_ptrs[gf][base_offset], &coeff_3d[ix2][ix1][0]);

        } // END LOOP x1 direction
      } // END LOOP x2 direction

      // Store the interpolated value for this grid function and destination point.
      dst_data[gf][dst_pt] = sum * src_invdxx012_INTERP_ORDERm1;
    } // END LOOP: Over grid functions.

  } // END PARALLEL FOR: Interpolate all destination points.

  return error_flag;
"""

    postfunc = r"""
#pragma GCC reset_options // Reset compiler optimizations after the function

#ifdef STANDALONE

#include <omp.h>

#define NUM_INTERP_GFS 4
#define NUM_RESOLUTIONS 3
#define NUM_DST_PTS 8000000

// Analytic functions.
static inline REAL analytic_function1(REAL x0, REAL x1, REAL x2) { return sin(x0) * cos(x1) * exp(-x2 * x2); }
static inline REAL analytic_function2(REAL x0, REAL x1, REAL x2) { return cos(x0) * sin(x1) * exp(-x2); }
static inline REAL analytic_function3(REAL x0, REAL x1, REAL x2) { return sin(x0 + x1 + x2); }
static inline REAL analytic_function4(REAL x0, REAL x1, REAL x2) { return cos(x1) * sin(x0) + x2 * x2 * x2; }

/**
 * @brief Initializes the 1D coordinate arrays for a 3D uniform source grid.
 *
 * This function calculates the grid spacing (dx) for each dimension and allocates memory for
 * and populates the 1D coordinate arrays. The coordinate arrays include ghost zones.
 *
 * @param n_interp_ghosts The number of ghost zones on each side for interpolation.
 * @param N_x0 The number of interior grid points in the x0-dimension.
 * @param N_x1 The number of interior grid points in the x1-dimension.
 * @param N_x2 The number of interior grid points in the x2-dimension.
 * @param[out] src_x0x1x2 An array of 3 pointers. The function allocates memory for each and fills it with coordinate values.
 * @param[out] src_dxx0 Pointer to store the calculated grid spacing in the x0-dimension.
 * @param[out] src_dxx1 Pointer to store the calculated grid spacing in the x1-dimension.
 * @param[out] src_dxx2 Pointer to store the calculated grid spacing in the x2-dimension.
 * @param src_Nxx_plus_2NGHOSTS0 The total number of points in the x0-dimension, including ghost zones.
 * @param src_Nxx_plus_2NGHOSTS1 The total number of points in the x1-dimension, including ghost zones.
 * @param src_Nxx_plus_2NGHOSTS2 The total number of points in the x2-dimension, including ghost zones.
 * @return 0 on success, -1 on memory allocation failure.
 */
int initialize_coordinates(const int n_interp_ghosts, const int N_x0, const int N_x1, const int N_x2, REAL *src_x0x1x2[3], REAL *src_dxx0,
                           REAL *src_dxx1, REAL *src_dxx2, const int src_Nxx_plus_2NGHOSTS0, const int src_Nxx_plus_2NGHOSTS1,
                           const int src_Nxx_plus_2NGHOSTS2) {
  *src_dxx0 = (2.0 * M_PI) / N_x0;
  *src_dxx1 = (2.0 * M_PI) / N_x1;
  *src_dxx2 = (2.0) / N_x2;
  src_x0x1x2[0] = (REAL *)malloc(sizeof(REAL) * src_Nxx_plus_2NGHOSTS0);
  src_x0x1x2[1] = (REAL *)malloc(sizeof(REAL) * src_Nxx_plus_2NGHOSTS1);
  src_x0x1x2[2] = (REAL *)malloc(sizeof(REAL) * src_Nxx_plus_2NGHOSTS2);
  if (!src_x0x1x2[0] || !src_x0x1x2[1] || !src_x0x1x2[2]) {
    free(src_x0x1x2[0]);
    free(src_x0x1x2[1]);
    free(src_x0x1x2[2]);
    src_x0x1x2[0] = src_x0x1x2[1] = src_x0x1x2[2] = NULL;
    return -1;
  } // END IF: check for allocation failure.
  for (int i = 0; i < src_Nxx_plus_2NGHOSTS0; i++)
    src_x0x1x2[0][i] = (i - n_interp_ghosts) * (*src_dxx0);
  // END LOOP: initialize x0 coordinates.
  for (int i = 0; i < src_Nxx_plus_2NGHOSTS1; i++)
    src_x0x1x2[1][i] = (i - n_interp_ghosts) * (*src_dxx1);
  // END LOOP: initialize x1 coordinates.
  for (int i = 0; i < src_Nxx_plus_2NGHOSTS2; i++)
    src_x0x1x2[2][i] = -1.0 + (i - n_interp_ghosts) * (*src_dxx2);
  // END LOOP: initialize x2 coordinates.
  return 0;
} // END FUNCTION initialize_coordinates()

/**
 * @brief Populates a 3D source grid function with values from a given analytic function.
 *
 * This function iterates through all points of a 3D grid (including ghost zones)
 * and computes the value of the grid function at each point using the provided
 * analytic function pointer.
 *
 * @param src_Nxx_plus_2NGHOSTS0 The total number of points in the x0-dimension.
 * @param src_Nxx_plus_2NGHOSTS1 The total number of points in the x1-dimension.
 * @param src_Nxx_plus_2NGHOSTS2 The total number of points in the x2-dimension.
 * @param src_x0x1x2 The pre-initialized 1D coordinate arrays for each dimension.
 * @param[out] src_gf The 3D source grid function data array to be populated.
 * @param func A function pointer to the analytic function used to compute the values.
 */
void initialize_src_gf(const int src_Nxx_plus_2NGHOSTS0, const int src_Nxx_plus_2NGHOSTS1, const int src_Nxx_plus_2NGHOSTS2, REAL *src_x0x1x2[3],
                       REAL *src_gf, REAL (*func)(REAL, REAL, REAL)) {
  for (int i2 = 0; i2 < src_Nxx_plus_2NGHOSTS2; i2++) {
    for (int i1 = 0; i1 < src_Nxx_plus_2NGHOSTS1; i1++) {
      for (int i0 = 0; i0 < src_Nxx_plus_2NGHOSTS0; i0++) {
        const int idx = i0 + src_Nxx_plus_2NGHOSTS0 * (i1 + src_Nxx_plus_2NGHOSTS1 * i2);
        src_gf[idx] = func(src_x0x1x2[0][i0], src_x0x1x2[1][i1], src_x0x1x2[2][i2]);
      } // END LOOP: over i0.
    } // END LOOP: over i1.
  } // END LOOP: over i2.
} // END FUNCTION initialize_src_gf()

/**
 * @brief Main driver for testing 3D Lagrange interpolation.
 *
 * This program tests a 3D Lagrange interpolation routine by performing the following steps:
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
  REAL(*dst_pts)[3] = NULL;
  REAL *f_exact[NUM_INTERP_GFS] = {NULL};
  REAL *src_x0x1x2[3] = {NULL, NULL, NULL};
  REAL *src_gf[NUM_INTERP_GFS] = {NULL};
  REAL *dst_data[NUM_INTERP_GFS] = {NULL};

  const int n_interp_ghosts = 3;
  const int INTERP_ORDER = 2 * n_interp_ghosts + 1;

  int N_x0_arr[NUM_RESOLUTIONS] = {16, 32, 64};
  int N_x1_arr[NUM_RESOLUTIONS] = {16, 32, 64};
  int N_x2_arr[NUM_RESOLUTIONS] = {16, 32, 64};
  REAL h_arr[NUM_RESOLUTIONS];
  REAL error_L2_norm[NUM_INTERP_GFS][NUM_RESOLUTIONS];

  dst_pts = (REAL(*)[3])malloc(sizeof(REAL) * NUM_DST_PTS * 3);
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
    int N_x2 = N_x2_arr[res];
    h_arr[res] = (N_x0 > 0) ? ((REAL)(2.0 * M_PI) / N_x0) : 0.0;
    int src_Nxx_plus_2NGHOSTS0 = N_x0 + 2 * n_interp_ghosts;
    int src_Nxx_plus_2NGHOSTS1 = N_x1 + 2 * n_interp_ghosts;
    int src_Nxx_plus_2NGHOSTS2 = N_x2 + 2 * n_interp_ghosts;
    REAL src_dxx0_val, src_dxx1_val, src_dxx2_val;
    if (initialize_coordinates(n_interp_ghosts, N_x0, N_x1, N_x2, src_x0x1x2, &src_dxx0_val, &src_dxx1_val, &src_dxx2_val, src_Nxx_plus_2NGHOSTS0,
                               src_Nxx_plus_2NGHOSTS1, src_Nxx_plus_2NGHOSTS2) != 0) {
      fprintf(stderr, "malloc failed for coordinates.\n");
      return_code = EXIT_FAILURE;
      goto cleanup;
    } // END IF: initialize coordinates.
    REAL x0_min_safe = src_x0x1x2[0][n_interp_ghosts] + 1e-6;
    REAL x0_max_safe = src_x0x1x2[0][src_Nxx_plus_2NGHOSTS0 - n_interp_ghosts - 1] - 1e-6;
    REAL x1_min_safe = src_x0x1x2[1][n_interp_ghosts] + 1e-6;
    REAL x1_max_safe = src_x0x1x2[1][src_Nxx_plus_2NGHOSTS1 - n_interp_ghosts - 1] - 1e-6;
    REAL x2_min_safe = src_x0x1x2[2][n_interp_ghosts] + 1e-6;
    REAL x2_max_safe = src_x0x1x2[2][src_Nxx_plus_2NGHOSTS2 - n_interp_ghosts - 1] - 1e-6;
    srand(42 + res);
    for (int i = 0; i < NUM_DST_PTS; i++) {
      dst_pts[i][0] = x0_min_safe + ((REAL)rand() / RAND_MAX) * (x0_max_safe - x0_min_safe);
      dst_pts[i][1] = x1_min_safe + ((REAL)rand() / RAND_MAX) * (x1_max_safe - x1_min_safe);
      dst_pts[i][2] = x2_min_safe + ((REAL)rand() / RAND_MAX) * (x2_max_safe - x2_min_safe);
      f_exact[0][i] = analytic_function1(dst_pts[i][0], dst_pts[i][1], dst_pts[i][2]);
      f_exact[1][i] = analytic_function2(dst_pts[i][0], dst_pts[i][1], dst_pts[i][2]);
      f_exact[2][i] = analytic_function3(dst_pts[i][0], dst_pts[i][1], dst_pts[i][2]);
      f_exact[3][i] = analytic_function4(dst_pts[i][0], dst_pts[i][1], dst_pts[i][2]);
    } // END LOOP: initialize destination points and exact values.
    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      const size_t size = (size_t)src_Nxx_plus_2NGHOSTS0 * src_Nxx_plus_2NGHOSTS1 * src_Nxx_plus_2NGHOSTS2;
      src_gf[gf] = (REAL *)malloc(sizeof(REAL) * size);
      if (!src_gf[gf]) {
        fprintf(stderr, "malloc failed for src_gf.\n");
        return_code = EXIT_FAILURE;
        goto cleanup;
      } // END IF: check malloc for src_gf.
    } // END LOOP: allocate source grid functions.
    initialize_src_gf(src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_Nxx_plus_2NGHOSTS2, src_x0x1x2, src_gf[0], analytic_function1);
    initialize_src_gf(src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_Nxx_plus_2NGHOSTS2, src_x0x1x2, src_gf[1], analytic_function2);
    initialize_src_gf(src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_Nxx_plus_2NGHOSTS2, src_x0x1x2, src_gf[2], analytic_function3);
    initialize_src_gf(src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_Nxx_plus_2NGHOSTS2, src_x0x1x2, src_gf[3], analytic_function4);
    const REAL *restrict src_gf_ptrs[NUM_INTERP_GFS] = {src_gf[0], src_gf[1], src_gf[2], src_gf[3]};
    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      dst_data[gf] = (REAL *)malloc(sizeof(REAL) * NUM_DST_PTS);
      if (!dst_data[gf]) {
        fprintf(stderr, "malloc failed for dst_data.\n");
        return_code = EXIT_FAILURE;
        goto cleanup;
      } // END IF: check malloc for dst_data.
    } // END LOOP: allocate destination data arrays.
#ifdef _OPENMP
    double start_time = omp_get_wtime();
#endif
    int error_code = interpolation_3d_general__uniform_src_grid(n_interp_ghosts, src_dxx0_val, src_dxx1_val, src_dxx2_val, src_Nxx_plus_2NGHOSTS0,
                                                                src_Nxx_plus_2NGHOSTS1, src_Nxx_plus_2NGHOSTS2, NUM_INTERP_GFS, src_x0x1x2,
                                                                src_gf_ptrs, NUM_DST_PTS, dst_pts, dst_data);
#ifdef _OPENMP
    double elapsed_time = omp_get_wtime() - start_time;
#endif
    printf("\n--- Benchmarking for Resolution %d (%dx%dx%d) ---\n", res, N_x0, N_x1, N_x2);
#ifdef _OPENMP
    printf("Interpolated %d points in %.4f seconds.\n", NUM_DST_PTS, elapsed_time);
    printf("Performance: %.4f million points per second.\n", (double)NUM_DST_PTS / elapsed_time / 1e6);
#endif
    printf("--------------------------------------------------\n");
    if (error_code != INTERP_SUCCESS) {
      fprintf(stderr, "Interpolation error code: %d\n", error_code);
      return_code = error_code;
      goto cleanup;
    } // END IF: check for interpolation error.
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
    for (int dim = 0; dim < 3; dim++) {
      free(src_x0x1x2[dim]);
      src_x0x1x2[dim] = NULL;
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
  for (int dim = 0; dim < 3; dim++) {
    free(src_x0x1x2[dim]);
  } // END LOOP: cleanup src_x0x1x2.
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
