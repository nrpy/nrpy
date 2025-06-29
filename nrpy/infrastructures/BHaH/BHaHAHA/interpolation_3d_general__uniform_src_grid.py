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
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for general-purpose 3D Lagrange interpolation.

    Even if `enable_simd` is False, `intrinsics/simd_intrinsics.h` is still required.

    :param enable_simd: Whether the rest of the code enables SIMD optimizations, as this code requires simd_intrinsics.h (which includes SIMD-disabled options).
    :param project_dir: Directory of the project, to set the destination for simd_instrinsics.h .
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

    includes = ["stdio.h", "stdlib.h", "math.h", "intrinsics/simd_intrinsics.h"]

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
    params = """
    const int N_interp_GHOSTS, const REAL src_dxx0, const REAL src_dxx1, const REAL src_dxx2,
    const int src_Nxx_plus_2NGHOSTS0, const int src_Nxx_plus_2NGHOSTS1, const int src_Nxx_plus_2NGHOSTS2,
    const int NUM_INTERP_GFS, REAL *restrict src_x0x1x2[3], const REAL *restrict src_gf_ptrs[NUM_INTERP_GFS],
    const int num_dst_pts, const REAL dst_x0x1x2[][3], REAL *restrict dst_data[NUM_INTERP_GFS]"""

    body = r"""
  // Unpack parameters.
  const int NinterpGHOSTS = N_interp_GHOSTS;
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
  for (int i = 0; i < INTERP_ORDER; i++) {
    REAL denom = 1.0;
    for (int j = 0; j < i; j++)
      denom *= (REAL)(i - j);
    for (int j = i + 1; j < INTERP_ORDER; j++)
      denom *= (REAL)(i - j);
    inv_denom[i] = 1.0 / denom; // Divisions are expensive, so we do them only once.
  } // END LOOP: Precompute inverse denominators.

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
    REAL coeff_x0[INTERP_ORDER], coeff_x1[INTERP_ORDER], coeff_x2[INTERP_ORDER];
    REAL diffs_x0[INTERP_ORDER], diffs_x1[INTERP_ORDER], diffs_x2[INTERP_ORDER];

#pragma omp simd
    for (int j = 0; j < INTERP_ORDER; j++) {
      diffs_x0[j] = x0_dst - src_x0x1x2[0][base_idx_x0 + j];
      diffs_x1[j] = x1_dst - src_x0x1x2[1][base_idx_x1 + j];
      diffs_x2[j] = x2_dst - src_x0x1x2[2][base_idx_x2 + j];
    } // END LOOP: Compute differences for Lagrange interpolation.

    // Compute the numerator of the Lagrange basis polynomials.
#pragma omp simd
    for (int i = 0; i < INTERP_ORDER; i++) {
      REAL numer_i_x0 = 1.0, numer_i_x1 = 1.0, numer_i_x2 = 1.0;
      for (int j = 0; j < i; j++) {
        numer_i_x0 *= diffs_x0[j];
        numer_i_x1 *= diffs_x1[j];
        numer_i_x2 *= diffs_x2[j];
      }
      for (int j = i + 1; j < INTERP_ORDER; j++) {
        numer_i_x0 *= diffs_x0[j];
        numer_i_x1 *= diffs_x1[j];
        numer_i_x2 *= diffs_x2[j];
      }
      coeff_x0[i] = numer_i_x0 * inv_denom[i];
      coeff_x1[i] = numer_i_x1 * inv_denom[i];
      coeff_x2[i] = numer_i_x2 * inv_denom[i];
    } // END LOOP: Compute Lagrange basis polynomials.

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
        for (int ix1 = 0; ix1 < INTERP_ORDER; ix1++) {
          const int idx1 = base_idx_x1 + ix1;

          int ix0 = 0;
          REAL_SIMD_ARRAY vec_sum = SetZeroSIMD; // Initialize vector sum to zero

          // Precompute the base offset for the current ix2 and ix1
          // This avoids recalculating the constant part inside the vector loop
          const int base_offset = base_idx_x0 + src_Nxx_plus_2NGHOSTS0 * (idx1 + src_Nxx_plus_2NGHOSTS1 * idx2);

          // Vectorized loop using SIMD with FMA, if available
          for (; ix0 <= INTERP_ORDER - simd_width; ix0 += simd_width) { // Process simd_width doubles at a time
            // Calculate the flat index for the current set of ix0
            // Ensure that ix0 is added correctly to the base_offset
            const int current_idx0 = base_offset + ix0;

            // Load simd_width elements from src_gf_ptrs and coeff_3d
            REAL_SIMD_ARRAY vec_src = ReadSIMD(&src_gf_ptrs[gf][current_idx0]);
            REAL_SIMD_ARRAY vec_coeff = ReadSIMD(&coeff_3d[ix2][ix1][ix0]);
            // Use FMA to multiply src and coeff and add to vec_sum
            vec_sum = FusedMulAddSIMD(vec_src, vec_coeff, vec_sum);
          } // END LOOP x0 direction up to integer number of SIMD widths

          sum += HorizAddSIMD(vec_sum);

          // Handle remaining elements that don't fit into a full AVX register
          for (; ix0 < INTERP_ORDER; ix0++) {
            const int current_idx0 = base_offset + ix0;
            sum += src_gf_ptrs[gf][current_idx0] * coeff_3d[ix2][ix1][ix0];
          } // END LOOP remainder of x0 direction
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

#include <omp.h> // Include OpenMP header for the timer function.

// Define the number of grid functions to interpolate as a macro for compile-time constant.
#define NUM_INTERP_GFS 4 // Number of grid functions to interpolate.

// Define multiple analytic functions.
static inline REAL analytic_function1(REAL x0, REAL x1, REAL x2) { return sin(x0) * cos(x1) * exp(-x2 * x2); }

static inline REAL analytic_function2(REAL x0, REAL x1, REAL x2) { return cos(x0) * sin(x1) * exp(-x2); }

static inline REAL analytic_function3(REAL x0, REAL x1, REAL x2) { return sin(x0 + x1 + x2); }

static inline REAL analytic_function4(REAL x0, REAL x1, REAL x2) { return cos(x1) * sin(x0) + x2 * x2 * x2; }

/**
 * Initializes the coordinates for the source grid.
 *
 * @param N_interp_GHOSTS - Number of ghost zones.
 * @param N_x0, N_x1, N_x2 - Number of grid points in x0, x1, x2 directions.
 * @param src_x0x1x2 - Arrays to store coordinate values.
 * @param src_dxx0, src_dxx1, src_dxx2 - Pointers to grid spacings.
 * @param src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_Nxx_plus_2NGHOSTS2 - Dimensions including ghost zones.
 */
void initialize_coordinates(const int N_interp_GHOSTS, const int N_x0, const int N_x1, const int N_x2, REAL *src_x0x1x2[3], REAL *src_dxx0,
                            REAL *src_dxx1, REAL *src_dxx2, const int src_Nxx_plus_2NGHOSTS0, const int src_Nxx_plus_2NGHOSTS1,
                            const int src_Nxx_plus_2NGHOSTS2) {
  // Initialize grid spacings.
  *src_dxx0 = (2.0 * M_PI) / N_x0;
  *src_dxx1 = (2.0 * M_PI) / N_x1;
  *src_dxx2 = (2.0) / N_x2;

  // Allocate memory for coordinate arrays.
  src_x0x1x2[0] = (REAL *)malloc(sizeof(REAL) * src_Nxx_plus_2NGHOSTS0);
  src_x0x1x2[1] = (REAL *)malloc(sizeof(REAL) * src_Nxx_plus_2NGHOSTS1);
  src_x0x1x2[2] = (REAL *)malloc(sizeof(REAL) * src_Nxx_plus_2NGHOSTS2);

  // Initialize coordinates in each dimension.
  for (int i = 0; i < src_Nxx_plus_2NGHOSTS0; i++) {
    src_x0x1x2[0][i] = (i - N_interp_GHOSTS) * (*src_dxx0);
  }
  for (int i = 0; i < src_Nxx_plus_2NGHOSTS1; i++) {
    src_x0x1x2[1][i] = (i - N_interp_GHOSTS) * (*src_dxx1);
  }
  for (int i = 0; i < src_Nxx_plus_2NGHOSTS2; i++) {
    src_x0x1x2[2][i] = -1.0 + (i - N_interp_GHOSTS) * (*src_dxx2);
  }
} // END FUNCTION: initialize_coordinates.

/**
 * Initializes the source grid function with an analytic function.
 *
 * @param src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_Nxx_plus_2NGHOSTS2 - Grid dimensions including ghost zones.
 * @param src_x0x1x2 - Arrays of coordinate values.
 * @param src_gf - Source grid function to initialize.
 * @param func - Analytic function to compute grid values.
 */
void initialize_src_gf(const int src_Nxx_plus_2NGHOSTS0, const int src_Nxx_plus_2NGHOSTS1, const int src_Nxx_plus_2NGHOSTS2, REAL *src_x0x1x2[3],
                       REAL *src_gf, REAL (*func)(REAL, REAL, REAL)) {
  // Initialize source grid function using the provided analytic function.
  for (int i2 = 0; i2 < src_Nxx_plus_2NGHOSTS2; i2++) {
    for (int i1 = 0; i1 < src_Nxx_plus_2NGHOSTS1; i1++) {
      for (int i0 = 0; i0 < src_Nxx_plus_2NGHOSTS0; i0++) {
        const int idx = i0 + src_Nxx_plus_2NGHOSTS0 * (i1 + src_Nxx_plus_2NGHOSTS1 * i2);
        src_gf[idx] = func(src_x0x1x2[0][i0], src_x0x1x2[1][i1], src_x0x1x2[2][i2]);
      } // END LOOP: Over i0.
    } // END LOOP: Over i1.
  } // END LOOP: Over i2.
} // END FUNCTION: initialize_src_gf.

int main() {
  const int N_interp_GHOSTS = 4;                    // For 9th order interpolation.
  const int INTERP_ORDER = 2 * N_interp_GHOSTS + 1; // 9th order.
  const int num_resolutions = 3;                    // Number of resolutions to test.
  const int num_dst_pts = 3000000;                  // Number of destination points.

  int N_x0_arr[num_resolutions];
  int N_x1_arr[num_resolutions];
  int N_x2_arr[num_resolutions];
  REAL h_arr[num_resolutions];
  REAL error_L2_norm[NUM_INTERP_GFS][num_resolutions];

  // Initialize the resolutions.
  N_x0_arr[0] = 16;
  N_x1_arr[0] = 16;
  N_x2_arr[0] = 16;

  N_x0_arr[1] = 32;
  N_x1_arr[1] = 32;
  N_x2_arr[1] = 32;

  N_x0_arr[2] = 64;
  N_x1_arr[2] = 64;
  N_x2_arr[2] = 64;

  // Allocate memory for destination points.
  REAL(*dst_pts)[3] = (REAL(*)[3])malloc(sizeof(REAL) * num_dst_pts * 3);
  if (dst_pts == NULL) {
    fprintf(stderr, "Memory allocation failed for destination points.\n");
    return EXIT_FAILURE;
  }

  // Allocate exact solution arrays for each grid function.
  REAL *f_exact[NUM_INTERP_GFS];
  for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
    f_exact[gf] = (REAL *)malloc(sizeof(REAL) * num_dst_pts);
    if (f_exact[gf] == NULL) {
      fprintf(stderr, "Memory allocation failed for f_exact[%d].\n", gf);
      // Free previously allocated memory before exiting.
      for (int g = 0; g < gf; g++) {
        free(f_exact[g]);
      }
      free(dst_pts);
      return EXIT_FAILURE;
    }
  } // END LOOP: Allocate exact solution arrays.

  // Loop over resolutions.
  for (int res = 0; res < num_resolutions; res++) {
    int N_x0 = N_x0_arr[res];
    int N_x1 = N_x1_arr[res];
    int N_x2 = N_x2_arr[res];

    h_arr[res] = (N_x0 > 0) ? ((REAL)(2.0 * M_PI) / N_x0) : 0.0; // Assuming src_dxx0 == src_dxx1 == src_dxx2.

    // Define source grid dimensions including ghost zones.
    int src_Nxx_plus_2NGHOSTS0 = N_x0 + 2 * N_interp_GHOSTS;
    int src_Nxx_plus_2NGHOSTS1 = N_x1 + 2 * N_interp_GHOSTS;
    int src_Nxx_plus_2NGHOSTS2 = N_x2 + 2 * N_interp_GHOSTS;

    // Allocate and initialize coordinate arrays.
    REAL *src_x0x1x2[3];
    REAL src_dxx0_val, src_dxx1_val, src_dxx2_val;

    initialize_coordinates(N_interp_GHOSTS, N_x0, N_x1, N_x2, src_x0x1x2, &src_dxx0_val, &src_dxx1_val, &src_dxx2_val, src_Nxx_plus_2NGHOSTS0,
                           src_Nxx_plus_2NGHOSTS1, src_Nxx_plus_2NGHOSTS2);

    // Compute safe domain for destination points
    REAL x0_min_safe = src_x0x1x2[0][N_interp_GHOSTS] + 1e-6;
    REAL x0_max_safe = src_x0x1x2[0][src_Nxx_plus_2NGHOSTS0 - N_interp_GHOSTS - 1] - 1e-6;
    REAL x1_min_safe = src_x0x1x2[1][N_interp_GHOSTS] + 1e-6;
    REAL x1_max_safe = src_x0x1x2[1][src_Nxx_plus_2NGHOSTS1 - N_interp_GHOSTS - 1] - 1e-6;
    REAL x2_min_safe = src_x0x1x2[2][N_interp_GHOSTS] + 1e-6;
    REAL x2_max_safe = src_x0x1x2[2][src_Nxx_plus_2NGHOSTS2 - N_interp_GHOSTS - 1] - 1e-6;

    // Seed the random number generator.
    srand(42 + res); // Use different seed for each resolution if desired

    // Generate random destination points and compute the exact function values for each grid function.
    for (int i = 0; i < num_dst_pts; i++) {
      REAL x0 = x0_min_safe + ((REAL)rand() / RAND_MAX) * (x0_max_safe - x0_min_safe);
      REAL x1 = x1_min_safe + ((REAL)rand() / RAND_MAX) * (x1_max_safe - x1_min_safe);
      REAL x2 = x2_min_safe + ((REAL)rand() / RAND_MAX) * (x2_max_safe - x2_min_safe);
      dst_pts[i][0] = x0;
      dst_pts[i][1] = x1;
      dst_pts[i][2] = x2;
      // Compute exact values for each grid function.
      f_exact[0][i] = analytic_function1(x0, x1, x2);
      f_exact[1][i] = analytic_function2(x0, x1, x2);
      f_exact[2][i] = analytic_function3(x0, x1, x2);
      f_exact[3][i] = analytic_function4(x0, x1, x2);
    }

    // Allocate and initialize grid functions.
    REAL *src_gf[NUM_INTERP_GFS];
    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      src_gf[gf] = (REAL *)malloc(sizeof(REAL) * src_Nxx_plus_2NGHOSTS0 * src_Nxx_plus_2NGHOSTS1 * src_Nxx_plus_2NGHOSTS2);
      if (src_gf[gf] == NULL) {
        fprintf(stderr, "Memory allocation failed for src_gf[%d].\n", gf);
        // Free previously allocated memory before exiting.
        for (int g = 0; g < gf; g++) {
          free(src_gf[g]);
        }
        for (int dim = 0; dim < 3; dim++) {
          free(src_x0x1x2[dim]);
        }
        free(dst_pts);
        for (int g = 0; g < NUM_INTERP_GFS; g++) {
          free(f_exact[g]);
        }
        return EXIT_FAILURE;
      }
    } // END LOOP: Allocate grid functions.

    // Initialize each grid function with its respective analytic function.
    initialize_src_gf(src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_Nxx_plus_2NGHOSTS2, src_x0x1x2, src_gf[0], analytic_function1);
    initialize_src_gf(src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_Nxx_plus_2NGHOSTS2, src_x0x1x2, src_gf[1], analytic_function2);
    initialize_src_gf(src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_Nxx_plus_2NGHOSTS2, src_x0x1x2, src_gf[2], analytic_function3);
    initialize_src_gf(src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_Nxx_plus_2NGHOSTS2, src_x0x1x2, src_gf[3], analytic_function4);

    // Create an array of pointers to src_gf.
    const REAL *restrict src_gf_ptrs[NUM_INTERP_GFS] = {src_gf[0], src_gf[1], src_gf[2], src_gf[3]};

    // Allocate memory for interpolated data for all grid functions.
    REAL *dst_data[NUM_INTERP_GFS];
    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      dst_data[gf] = (REAL *)malloc(sizeof(REAL) * num_dst_pts);
      if (dst_data[gf] == NULL) {
        fprintf(stderr, "Memory allocation failed for dst_data[%d].\n", gf);
        // Free previously allocated memory before exiting.
        for (int g = 0; g < gf; g++) {
          free(dst_data[g]);
        }
        for (int g = 0; g < NUM_INTERP_GFS; g++) {
          free(src_gf[g]);
        }
        for (int dim = 0; dim < 3; dim++) {
          free(src_x0x1x2[dim]);
        }
        free(dst_pts);
        for (int g = 0; g < NUM_INTERP_GFS; g++) {
          free(f_exact[g]);
        }
        return EXIT_FAILURE;
      }
    }

    // --- Start of benchmarking ---
    double start_time = omp_get_wtime();

    // Call the interpolation function.
    int error_code = interpolation_3d_general__uniform_src_grid(N_interp_GHOSTS, src_dxx0_val, src_dxx1_val, src_dxx2_val, src_Nxx_plus_2NGHOSTS0,
                                                                src_Nxx_plus_2NGHOSTS1, src_Nxx_plus_2NGHOSTS2, NUM_INTERP_GFS, src_x0x1x2,
                                                                src_gf_ptrs, num_dst_pts, dst_pts, dst_data);

    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;
    // --- End of benchmarking ---

    // Report performance metrics.
    printf("\n--- Benchmarking for Resolution %d (%dx%dx%d) ---\n", res, N_x0, N_x1, N_x2);
    printf("Interpolated %d points in %.4f seconds.\n", num_dst_pts, elapsed_time);
    printf("Performance: %.4f million points per second.\n", (double)num_dst_pts / elapsed_time / 1e6);
    printf("--------------------------------------------------\n");

    if (error_code != INTERP_SUCCESS) {
      fprintf(stderr, "Interpolation error code: %d\n", error_code);
      // Free allocated memory before exiting.
      for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
        free(src_gf[gf]);
        free(dst_data[gf]);
      }
      for (int dim = 0; dim < 3; dim++) {
        free(src_x0x1x2[dim]);
      }
      for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
        free(f_exact[gf]);
      }
      free(dst_pts);
      return error_code;
    }

    // Compute the L2 norm of the error for each grid function.
    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      REAL error_sum = 0.0;
      for (int i = 0; i < num_dst_pts; i++) {
        REAL error = dst_data[gf][i] - f_exact[gf][i];
        error_sum += error * error;
      }
      error_L2_norm[gf][res] = sqrt(error_sum / num_dst_pts);
    } // END LOOP: Compute L2 norms.

    // Output the error for each grid function.
    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      printf("Resolution %d: N_x0 = %d, N_x1 = %d, N_x2 = %d, h = %.5e, Grid Function %d, L2 error = %.5e\n", res, N_x0, N_x1, N_x2, h_arr[res],
             gf + 1, error_L2_norm[gf][res]);
    } // END LOOP: Output errors.

    // Free allocated memory for this resolution.
    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      free(src_gf[gf]);
      free(dst_data[gf]);
    }
    for (int dim = 0; dim < 3; dim++) {
      free(src_x0x1x2[dim]);
    }
  } // END LOOP: Over resolutions.

  printf("\n--- Convergence Results ---\n");
  // Compute the observed order of convergence for each grid function.
  for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
    for (int res = 1; res < num_resolutions; res++) {
      REAL observed_order = log(error_L2_norm[gf][res - 1] / error_L2_norm[gf][res]) / log(h_arr[res - 1] / h_arr[res]);
      printf("Observed order of convergence for Grid Function %d between resolutions %d and %d: %.2f\n", gf + 1, res - 1, res, observed_order);
    }
    // Expected order is INTERP_ORDER (since we are using 9th order interpolation).
    printf("Expected order of convergence for Grid Function %d: %d\n", gf + 1, INTERP_ORDER);
  } // END LOOP: Compute observed orders.

  // Clean up.
  for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
    free(f_exact[gf]);
  }
  free(dst_pts);

  return 0;
} // END FUNCTION: main.

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
