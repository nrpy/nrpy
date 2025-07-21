#include "intrinsics/simd_intrinsics.h"

#ifndef REAL
#define REAL double
#endif

static inline void compute_inv_denom(const int n_interp_ghosts, REAL inv_denom[(2 * n_interp_ghosts + 1)]) {
  // Precompute inverse denominators for Lagrange interpolation coefficients to optimize performance.
  for (int i = 0; i < (2 * n_interp_ghosts + 1); i++) {
    REAL denom = 1.0;
    for (int j = 0; j < i; j++)
      denom *= (REAL)(i - j);
    for (int j = i + 1; j < (2 * n_interp_ghosts + 1); j++)
      denom *= (REAL)(i - j);
    inv_denom[i] = 1.0 / denom; // Store the inverse to avoid repeated division operations.
  } // END LOOP: Precompute inverse denominators.
} // END FUNCTION compute_inv_denom()

static inline void compute_diffs_xi(const int n_interp_ghosts, const REAL dst_xi, const REAL *restrict src_xi_stencil,
                                    REAL diffs[(2 * n_interp_ghosts + 1)]) {
#pragma omp simd
  for (int j = 0; j < (2 * n_interp_ghosts + 1); j++) {
    diffs[j] = dst_xi - src_xi_stencil[j];
  } // END LOOP over j.
} // END FUNCTION compute_diffs_xi()

static inline void compute_lagrange_basis_coeffs_xi(const int n_interp_ghosts, const REAL inv_denom[(2 * n_interp_ghosts + 1)],
                                                    const REAL diffs[(2 * n_interp_ghosts + 1)],
                                                    REAL lagrange_basis_coeffs_xi[(2 * n_interp_ghosts + 1)]) {
#pragma omp simd
  for (int i = 0; i < (2 * n_interp_ghosts + 1); i++) {
    REAL numer_i = 1.0;
    // Compute product for j < i (scalar loop).
    for (int j = 0; j < i; j++) {
      numer_i *= diffs[j];
    } // END LOOP over j < i.
    // Compute product for j > i (scalar loop).
    for (int j = i + 1; j < (2 * n_interp_ghosts + 1); j++) {
      numer_i *= diffs[j];
    } // END LOOP over j > i.

    lagrange_basis_coeffs_xi[i] = numer_i * inv_denom[i]; // Scale by the inverse denominator.
  } // END LOOP over i.
} // END FUNCTION compute_lagrange_basis_coeffs_xi()

static inline REAL sum_lagrange_x0_simd(const int n_interp_ghosts, const REAL *restrict src_gf_base_idx,
                                        const REAL *restrict lagrange_basis_coeffs_x0_base_idx) {
  REAL sum = 0;

  // Vectorized loop using SIMD with FMA, if available
  int ix0 = 0;

  // For AVX-256/AVX-512 CPUs, when 3 < INTERP_ORDER = (2 * n_interp_ghosts + 1) < 8, then use AVX-256 SIMD instructions.
#if (2 * n_interp_ghosts + 1) > 3 && (2 * n_interp_ghosts + 1) < 8 && ((simd_width == 8) || (simd_width == 4))
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
  for (; ix0 <= (2 * n_interp_ghosts + 1) - simd_width; ix0 += simd_width) {
    vec_sum = FusedMulAddSIMD(ReadSIMD(&src_gf_base_idx[ix0]), ReadSIMD(&lagrange_basis_coeffs_x0_base_idx[ix0]), vec_sum);
  } // END LOOP x0 direction up to integer number of SIMD widths
  // Accumulate SIMD result
  sum += HorizAddSIMD(vec_sum);
#endif

  // Handle remaining elements that don't fit into a full SIMD register
  for (; ix0 < (2 * n_interp_ghosts + 1); ix0++) {
    sum += src_gf_base_idx[ix0] * lagrange_basis_coeffs_x0_base_idx[ix0];
  } // END LOOP remainder of x0 direction
  return sum;
} // END FUNCTION lagrange_sum_x0()
