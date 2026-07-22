#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "interpolation_lagrange_uniform.h"

/**
 * Interpolate the selected geometry bundle in physical time.
 *
 * The caller supplies flat per-slice spatial-interpolation outputs for the
 * configured temporal stencil, plus the corresponding physical `slice_times`.
 * This helper assumes those times are trusted, finite, and strictly increasing,
 * derives the actual number of time nodes from
 * `commondata->numerical_spacetime_temporal_interp_order`, and builds one shared
 * nonuniform barycentric Lagrange basis from the exact physical slice times.
 *
 * The 10 serialized `g4DD` components are interpolated independently to
 * `t_target`. For `g4DD`, analytic derivatives of that same temporal
 * interpolation basis obtain the 10 coordinate-time metric derivatives, while
 * the 30 Cartesian spatial metric derivatives are interpolated independently.
 * For `g4DD_d0` and `GammaUDD`, all forty components of the secondary geometry
 * bundle are interpolated independently.
 *
 * The metric bundle ordering is:
 * `g4DD00, g4DD01, g4DD02, g4DD03, g4DD11, g4DD12, g4DD13, g4DD22, g4DD23,
 * g4DD33`.
 *
 * For metric methods, the 40-component secondary bundle stores derivatives with
 * metric-pair outermost and derivative direction innermost:
 * `g4DD_dD000, g4DD_dD001, g4DD_dD002, g4DD_dD003, g4DD_dD010, ...,
 * g4DD_dD333`. The final index is `(t, x, y, z) = (0, 1, 2, 3)`. Thus bundle
 * slot `4*p` receives the temporal derivative of metric component `p`, while
 * slots `4*p+1`, `4*p+2`, and `4*p+3` receive the temporally interpolated
 * Cartesian spatial derivatives.
 *
 * @param[in] commondata Common runtime parameters.
 * @param[in] slice_times Trusted, finite, strictly increasing physical coordinate
 * times for the supplied slices.
 * @param[in] g4dd_slices Flat per-slice metric components in kernel order.
 * @param[in] geometry_slices Flat per-slice metric derivatives or Christoffel
 * components in the selected method's storage order.
 * @param[in] t_target Target physical coordinate time.
 * @param[out] g4dd_out Final interpolated metric components.
 * @param[out] rhs_geometry_out Final metric-derivative or Christoffel bundle.
 * @return `TEMPORAL_LAGRANGE_INTERP_SUCCESS` on success, or
 * `TEMPORAL_LAGRANGE_INTERP_INVALID_ORDER` if the configured interpolation
 * half-width is outside the supported range.
 *
 * @note Callers must provide trusted, finite, strictly increasing physical
 * `slice_times`, not abstract slot indices.
 */
int temporal_lagrange_interpolation(const commondata_struct *restrict commondata, const REAL *restrict slice_times, const REAL *restrict g4dd_slices,
                                    const REAL *restrict geometry_slices, const REAL t_target, REAL *restrict g4dd_out,
                                    REAL *restrict rhs_geometry_out) {
  // Step 1: Build the shared nonuniform barycentric Lagrange basis and its
  // analytic derivative in normalized physical time.
  const int temporal_half_width = commondata->numerical_spacetime_temporal_interp_order;
  if (temporal_half_width < 0 || temporal_half_width > TEMPORAL_LAGRANGE_INTERP_MAX_HALF_WIDTH)
    return TEMPORAL_LAGRANGE_INTERP_INVALID_ORDER;
  const int interp_order = 2 * temporal_half_width + 1;
  REAL barycentric_weights[interp_order];
  REAL normalized_slice_times[interp_order];
  REAL coeff_t[interp_order];
  REAL coeff_dt_normalized[interp_order];
  const REAL time_origin = slice_times[0];
  const REAL time_scale = (interp_order > 1) ? (slice_times[interp_order - 1] - slice_times[0]) : 1.0;
  const REAL normalized_t_target = (t_target - time_origin) / time_scale;

  // Step 1.a: Normalize temporal nodes to one O(1) interval.
  for (int i = 0; i < interp_order; i++)
    normalized_slice_times[i] = (slice_times[i] - time_origin) / time_scale;

  // Step 1.b: Build barycentric interpolation weights from the normalized
  // temporal nodes.
  for (int i = 0; i < interp_order; i++) {
    REAL weight_denom = 1.0;
    for (int j = 0; j < interp_order; j++) {
      if (j != i) {
        const REAL time_diff = normalized_slice_times[i] - normalized_slice_times[j];
        weight_denom *= time_diff;
      } // END IF: multiplying one nontrivial barycentric denominator factor
    } // END LOOP: for j over temporal nodes while building one barycentric weight
    barycentric_weights[i] = 1.0 / weight_denom;
  } // END LOOP: for i over temporal nodes while building barycentric weights

  // Step 1.c: Check whether the target time exactly matches one supplied node.
  int exact_node = -1;
  for (int i = 0; i < interp_order; i++) {
    const REAL target_diff = normalized_t_target - normalized_slice_times[i];
    if (target_diff == 0.0) {
      exact_node = i;
      break;
    } // END IF: target time exactly matched one supplied slice time
  } // END LOOP: for i over temporal nodes while checking for an exact target match

  if (exact_node >= 0) {
    // At an exact node the value basis is one-hot, but its derivative is not
    // generally zero. Use the barycentric differentiation-matrix row for the
    // matched node.
    REAL exact_diagonal_derivative = 0.0;
    for (int i = 0; i < interp_order; i++) {
      coeff_t[i] = (i == exact_node) ? 1.0 : 0.0;
      if (i == exact_node) {
        coeff_dt_normalized[i] = 0.0;
      } else {
        coeff_dt_normalized[i] =
            barycentric_weights[i] / (barycentric_weights[exact_node] * (normalized_slice_times[exact_node] - normalized_slice_times[i]));
        exact_diagonal_derivative -= coeff_dt_normalized[i];
      } // END ELSE: one off-diagonal exact-node derivative coefficient
    } // END LOOP: for i over exact-node value and derivative coefficients
    coeff_dt_normalized[exact_node] = exact_diagonal_derivative;
  } else {
    REAL barycentric_sum = 0.0;
    REAL barycentric_inverse_square_sum = 0.0;

    for (int i = 0; i < interp_order; i++) {
      const REAL target_diff = normalized_t_target - normalized_slice_times[i];
      const REAL weighted_term = barycentric_weights[i] / target_diff;
      coeff_t[i] = weighted_term;
      barycentric_sum += weighted_term;
      barycentric_inverse_square_sum += weighted_term / target_diff;
    } // END LOOP: for i over temporal nodes while summing barycentric terms

    for (int i = 0; i < interp_order; i++) {
      const REAL target_diff = normalized_t_target - normalized_slice_times[i];
      coeff_t[i] /= barycentric_sum;
      coeff_dt_normalized[i] = coeff_t[i] * (barycentric_inverse_square_sum / barycentric_sum - 1.0 / target_diff);
    } // END LOOP: for i over normalized value and derivative coefficients
  } // END ELSE: target time required a full barycentric basis evaluation

  // Step 2: Interpolate each metric component and differentiate the same
  // temporal reconstruction with respect to physical coordinate time.
  for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
    REAL component_series[interp_order];

    for (int s = 0; s < interp_order; s++) {
      component_series[s] = g4dd_slices[s * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT + comp];
    } // END LOOP: for s over trusted metric slices for one component
    g4dd_out[comp] = sum_lagrange_x0_simd(interp_order, component_series, coeff_t);

    const int temporal_derivative_slot =
        TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_DIRECTION_COUNT * comp + TEMPORAL_LAGRANGE_INTERP_TEMPORAL_DERIVATIVE_INDEX;
    rhs_geometry_out[temporal_derivative_slot] = sum_lagrange_x0_simd(interp_order, component_series, coeff_dt_normalized) / time_scale;
  } // END LOOP: for comp over serialized metric components

  // Step 3: Interpolate only the 30 Cartesian spatial metric derivatives.
  // The derivative direction is the fastest-changing bundle index, so each
  // metric component occupies four consecutive slots ordered as (t, x, y, z).
  for (int metric_comp = 0; metric_comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; metric_comp++) {
    for (int derivative_direction = 1; derivative_direction < TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_DIRECTION_COUNT; derivative_direction++) {
      REAL component_series[interp_order];
      const int derivative_slot = TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_DIRECTION_COUNT * metric_comp + derivative_direction;

      for (int s = 0; s < interp_order; s++) {
        component_series[s] = geometry_slices[s * TEMPORAL_LAGRANGE_INTERP_GEOMETRY_COMPONENT_COUNT + derivative_slot];
      } // END LOOP: for s over trusted spatial-derivative slices for one component
      rhs_geometry_out[derivative_slot] = sum_lagrange_x0_simd(interp_order, component_series, coeff_t);
    } // END LOOP: for derivative_direction over Cartesian spatial directions
  } // END LOOP: for metric_comp over serialized metric components

  return TEMPORAL_LAGRANGE_INTERP_SUCCESS;
} // END FUNCTION: temporal_lagrange_interpolation
