#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "interpolation_lagrange_uniform.h"

/**
 * Interpolate trusted Cartesian geodesic tensors in physical time.
 *
 * The caller supplies flat per-slice spatial-interpolation outputs for the
 * configured temporal stencil, plus the corresponding physical `slice_times`.
 * This helper assumes those times are trusted, finite, and strictly increasing,
 * derives the actual number of time nodes from
 * `commondata->numerical_spacetime_temporal_interp_order`, builds one shared
 * nonuniform barycentric Lagrange basis in time from the exact physical
 * `slice_times`, and interpolates each serialized `g4DD` and `Gamma4UDD`
 * component independently to `t_target`.
 *
 * The metric bundle ordering matches the geodesic interpolation-kernel contract:
 * `g4DD00, g4DD01, g4DD02, g4DD03, g4DD11, g4DD12, g4DD13, g4DD22, g4DD23,
 * g4DD33`.
 *
 * The Christoffel bundle ordering also matches the geodesic interpolation-kernel
 * contract: `Gamma4UDD<alpha><mu><nu>` with `alpha` outermost and `(mu, nu)` in
 * upper-triangular order, i.e.
 * `Gamma4UDD000, Gamma4UDD001, Gamma4UDD002, Gamma4UDD003, Gamma4UDD011,
 * Gamma4UDD012, ..., Gamma4UDD333`.
 *
 * @param[in] commondata Common runtime parameters.
 * @param[in] slice_times Trusted, finite, strictly increasing physical coordinate
 * times for the supplied slices.
 * @param[in] g4dd_slices Flat per-slice metric components in kernel order.
 * @param[in] gamma4udd_slices Flat per-slice Christoffel components in kernel order.
 * @param[in] t_target Target physical coordinate time.
 * @param[out] g4dd_out Final interpolated metric components.
 * @param[out] gamma4udd_out Final interpolated Christoffel components.
 * @return `TEMPORAL_LAGRANGE_INTERP_SUCCESS` on success, or
 * `TEMPORAL_LAGRANGE_INTERP_INVALID_ORDER` if the configured interpolation
 * half-width is outside the supported range.
 *
 * @note Callers must provide trusted, finite, strictly increasing physical
 * `slice_times`, not abstract slot indices.
 */
int temporal_lagrange_interpolation(const commondata_struct *restrict commondata, const REAL *restrict slice_times, const REAL *restrict g4dd_slices,
                                    const REAL *restrict gamma4udd_slices, const REAL t_target, REAL *restrict g4dd_out,
                                    REAL *restrict gamma4udd_out) {
  // Step 1: Build the shared nonuniform barycentric Lagrange basis in
  // normalized physical time.
  const int temporal_half_width = commondata->numerical_spacetime_temporal_interp_order;
  if (temporal_half_width < 0 || temporal_half_width > TEMPORAL_LAGRANGE_INTERP_MAX_HALF_WIDTH)
    return TEMPORAL_LAGRANGE_INTERP_INVALID_ORDER;
  const int interp_order = 2 * temporal_half_width + 1;
  REAL barycentric_weights[interp_order];
  REAL normalized_slice_times[interp_order];
  REAL coeff_t[interp_order];
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

  // Step 1.c: Evaluate the barycentric basis at the requested target time.
  int exact_node = -1;
  for (int i = 0; i < interp_order; i++) {
    const REAL target_diff = normalized_t_target - normalized_slice_times[i];
    if (target_diff == 0.0) {
      exact_node = i;
      break;
    } // END IF: target time exactly matched one supplied slice time
  } // END LOOP: for i over temporal nodes while checking for an exact target match

  if (exact_node >= 0) {
    for (int i = 0; i < interp_order; i++)
      coeff_t[i] = (i == exact_node) ? 1.0 : 0.0;
  } else {
    REAL barycentric_sum = 0.0;

    for (int i = 0; i < interp_order; i++) {
      const REAL target_diff = normalized_t_target - normalized_slice_times[i];
      const REAL weighted_term = barycentric_weights[i] / target_diff;
      coeff_t[i] = weighted_term;
      barycentric_sum += weighted_term;
    } // END LOOP: for i over temporal nodes while summing barycentric terms

    for (int i = 0; i < interp_order; i++)
      coeff_t[i] /= barycentric_sum;
  } // END ELSE: target time required a full barycentric basis evaluation

  // Step 2: Interpolate the serialized metric components independently in time.
  for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
    REAL component_series[interp_order];

    for (int s = 0; s < interp_order; s++) {
      component_series[s] = g4dd_slices[s * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT + comp];
    } // END LOOP: for s over trusted metric slices for one component
    g4dd_out[comp] = sum_lagrange_x0_simd(interp_order, component_series, coeff_t);
  } // END LOOP: for comp over serialized metric components

  // Step 3: Interpolate the serialized Christoffel components independently in time.
  for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
    REAL component_series[interp_order];

    for (int s = 0; s < interp_order; s++) {
      component_series[s] = gamma4udd_slices[s * TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT + comp];
    } // END LOOP: for s over trusted Christoffel slices for one component
    gamma4udd_out[comp] = sum_lagrange_x0_simd(interp_order, component_series, coeff_t);
  } // END LOOP: for comp over serialized Christoffel components

  return TEMPORAL_LAGRANGE_INTERP_SUCCESS;
} // END FUNCTION: temporal_lagrange_interpolation
