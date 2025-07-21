/**
 * @file interpolation_lagrange_uniform.h
 * @brief Lagrange interpolation helpers on a uniform source grid.
 */

#ifndef INTERPOLATION_LAGRANGE_UNIFORM_H_
#define INTERPOLATION_LAGRANGE_UNIFORM_H_

#include "intrinsics/simd_intrinsics.h"

/**
 * @brief Define REAL data type. Defaults to double.
 */
#ifndef REAL
#define REAL double
#endif

/**
 * @brief Precompute inverse denominators for Lagrange interpolation coefficients.
 *
 * Precomputes inverse denominators to optimize performance by avoiding repeated divisions during interpolation.  This significantly speeds up the
 * computation of Lagrange basis coefficients.
 *
 * @param INTERP_ORDER The order of the interpolation.
 * @param inv_denom Array to store the precomputed inverse denominators.  Must be of size INTERP_ORDER.
 */
static inline void compute_inv_denom(const int INTERP_ORDER, REAL inv_denom[INTERP_ORDER]) {
  for (int i = 0; i < INTERP_ORDER; i++) {
    REAL denom = 1.0;
    for (int j = 0; j < i; j++)
      denom *= (REAL)(i - j);
    for (int j = i + 1; j < INTERP_ORDER; j++)
      denom *= (REAL)(i - j);
    inv_denom[i] = 1.0 / denom;
  } // END LOOP: Precompute inverse denominators.
} // END FUNCTION compute_inv_denom()

/**
 * @brief Compute differences between a destination point and source stencil points.
 *
 * Calculates the differences between a destination coordinate and each coordinate in the source stencil.  These differences are crucial for computing
 * Lagrange basis coefficients.
 *
 * @param INTERP_ORDER The order of the interpolation.
 * @param dst_xi The destination coordinate.
 * @param src_xi_stencil Pointer to the array of source stencil coordinates.
 * @param diffs Array to store the computed differences. Must be of size INTERP_ORDER.
 */
static inline void compute_diffs_xi(const int INTERP_ORDER, const REAL dst_xi, const REAL *restrict src_xi_stencil, REAL diffs[INTERP_ORDER]) {
#pragma omp simd
  for (int j = 0; j < INTERP_ORDER; j++) {
    diffs[j] = dst_xi - src_xi_stencil[j];
  } // END LOOP over j.
} // END FUNCTION compute_diffs_xi()

/**
 * @brief Compute Lagrange basis coefficients.
 *
 * Computes the Lagrange basis coefficients using precomputed inverse denominators and the differences between the destination and source stencil
 * coordinates.
 *
 * @param INTERP_ORDER The order of the interpolation.
 * @param inv_denom Array of precomputed inverse denominators.
 * @param diffs Array of differences between destination and source stencil coordinates.
 * @param lagrange_basis_coeffs_xi Array to store the computed Lagrange basis coefficients. Must be of size INTERP_ORDER.
 */
static inline void compute_lagrange_basis_coeffs_xi(const int INTERP_ORDER, const REAL inv_denom[INTERP_ORDER], const REAL diffs[INTERP_ORDER],
                                                    REAL lagrange_basis_coeffs_xi[INTERP_ORDER]) {
#pragma omp simd
  for (int i = 0; i < INTERP_ORDER; i++) {
    REAL numer_i = 1.0;
    // Compute product for j < i (scalar loop).
    for (int j = 0; j < i; j++) {
      numer_i *= diffs[j];
    } // END LOOP over j < i.
    // Compute product for j > i (scalar loop).
    for (int j = i + 1; j < INTERP_ORDER; j++) {
      numer_i *= diffs[j];
    } // END LOOP over j > i.
    lagrange_basis_coeffs_xi[i] = numer_i * inv_denom[i];
  } // END LOOP over i.
} // END FUNCTION compute_lagrange_basis_coeffs_xi()

/**
 * @brief Compute the weighted sum of Lagrange basis coefficients and source data using SIMD instructions where possible.
 *
 * This function computes the weighted sum of the Lagrange basis coefficients and the corresponding source data values.  It leverages SIMD
 * instructions for improved performance when the interpolation order allows.  The code handles cases where the interpolation order is not a multiple
 * of the SIMD vector width.
 *
 * @param INTERP_ORDER The order of the interpolation.
 * @param src_gf_base_idx Pointer to the array of source data values.
 * @param lagrange_basis_coeffs_x0_base_idx Pointer to the array of Lagrange basis coefficients.
 * @return The weighted sum of the Lagrange interpolation.
 */
static inline REAL sum_lagrange_x0_simd(const int INTERP_ORDER, const REAL *restrict src_gf_base_idx,
                                        const REAL *restrict lagrange_basis_coeffs_x0_base_idx) {
  REAL sum = 0;

  // Vectorized loop over ix0 using SIMD with FMA, if available
  int ix0 = 0;
  // For AVX-256/AVX-512 CPUs, when 3 < INTERP_ORDER = INTERP_ORDER < 8, then use AVX-256 SIMD instructions.
#if INTERP_ORDER > 3 && INTERP_ORDER < 8 && ((simd_width == 8) || (simd_width == 4))
  // AVX-256 MulSIMD():
  REAL_SIMD_ARRAY vec_sum = _mm256_mul_pd(_mm256_loadu_pd(&src_gf_base_idx[ix0]), _mm256_loadu_pd(&lagrange_basis_coeffs_x0_base_idx[ix0]));
  // AVX-256 HorizAddSIMD():
  sum += ({
    const __m128d low_128 = _mm256_castpd256_pd128(vec_sum);                   /* Extract lower 128 bits */
    const __m128d high_128 = _mm256_extractf128_pd(vec_sum, 1);                /* Extract upper 128 bits */
    const __m128d sum_128 = _mm_add_pd(low_128, high_128);                     /* Add low and high parts */
    _mm_cvtsd_f64(sum_128) + _mm_cvtsd_f64(_mm_unpackhi_pd(sum_128, sum_128)); /* Final scalar sum */
  });
  ix0 = 4;
#else
  REAL_SIMD_ARRAY vec_sum = SetZeroSIMD;
  for (; ix0 <= INTERP_ORDER - simd_width; ix0 += simd_width) {
    vec_sum = FusedMulAddSIMD(ReadSIMD(&src_gf_base_idx[ix0]), ReadSIMD(&lagrange_basis_coeffs_x0_base_idx[ix0]), vec_sum);
  } // END LOOP x0 direction up to integer number of SIMD widths
  // Accumulate SIMD result
  sum += HorizAddSIMD(vec_sum);
#endif

  // Handle remaining elements that don't fit into a full SIMD register
  for (; ix0 < INTERP_ORDER; ix0++) {
    sum += src_gf_base_idx[ix0] * lagrange_basis_coeffs_x0_base_idx[ix0];
  } // END LOOP remainder of x0 direction
  return sum;
} // END FUNCTION lagrange_sum_x0()

#endif // INTERPOLATION_LAGRANGE_UNIFORM_H_
