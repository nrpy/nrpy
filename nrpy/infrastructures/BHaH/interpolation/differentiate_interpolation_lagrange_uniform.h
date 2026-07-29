/**
 * @file differentiate_interpolation_lagrange_uniform.h
 * @brief Analytic derivatives of uniform-grid Lagrange basis coefficients.
 */

#ifndef DIFFERENTIATE_INTERPOLATION_LAGRANGE_UNIFORM_H_
#define DIFFERENTIATE_INTERPOLATION_LAGRANGE_UNIFORM_H_

#include "interpolation_lagrange_uniform.h"

/**
 * Compute derivatives of the unnormalized uniform-grid Lagrange basis.
 *
 * The companion value helper returns
 * `inv_denom[i] * product_{j != i}(x - x_j)`. This helper returns its
 * analytic derivative by summing products with one factor omitted. Callers
 * apply the same uniform-grid normalization used for the value basis.
 *
 * @param interp_order Number of stencil nodes.
 * @param[in] inv_denom Precomputed index-space inverse denominators.
 * @param[in] diffs Differences x - x_j for all stencil nodes.
 * @param[out] derivative_coeffs Analytic derivative coefficients.
 */
static inline void compute_lagrange_basis_derivative_coeffs_xi(
    const int interp_order,
    const REAL *RESTRICT inv_denom,
    const REAL *RESTRICT diffs,
    REAL *RESTRICT derivative_coeffs) {
  for (int i = 0; i < interp_order; i++) {
    REAL derivative_sum = 0.0;
    for (int omitted = 0; omitted < interp_order; omitted++) {
      if (omitted == i)
        continue;
      REAL product = 1.0;
      for (int j = 0; j < interp_order; j++) {
        if (j != i && j != omitted)
          product *= diffs[j];
      } // END LOOP: for j over retained Lagrange numerator factors
      derivative_sum += product;
    } // END LOOP: for omitted over differentiated Lagrange numerator factors
    derivative_coeffs[i] = inv_denom[i] * derivative_sum;
  } // END LOOP: for i over Lagrange derivative coefficients
} // END FUNCTION: compute_lagrange_basis_derivative_coeffs_xi

#endif // DIFFERENTIATE_INTERPOLATION_LAGRANGE_UNIFORM_H_
