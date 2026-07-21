#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "interpolation/interpolation_lagrange_uniform.h"

/**
 * Interpolate trusted Cartesian geodesic tensors in physical time.
 *
 * The caller supplies flat per-slice spatial-interpolation outputs for the
 * configured temporal stencil, plus the corresponding physical `slice_times`.
 * This helper assumes those times are strictly increasing and uniformly spaced,
 * derives the actual number of time nodes from
 * `commondata->numerical_spacetime_temporal_interp_order`, builds one shared 1D
 * Lagrange basis in time, and interpolates each serialized `g4DD` and
 * `Gamma4UDD` component independently to `t_target`.
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
 * @param[in] slice_times Trusted physical coordinate times for the supplied slices.
 * @param[in] g4dd_slices Flat per-slice metric components in kernel order.
 * @param[in] gamma4udd_slices Flat per-slice Christoffel components in kernel order.
 * @param[in] t_target Target physical coordinate time.
 * @param[out] g4dd_out Final interpolated metric components.
 * @param[out] gamma4udd_out Final interpolated Christoffel components.
 * @return Status code indicating success or invalid runtime interpolation order.
 *
 * @note Callers must provide trusted, strictly increasing, uniformly spaced
 * physical `slice_times`, not abstract slot indices.
 */
int temporal_lagrange_interpolation(const commondata_struct *restrict commondata, const REAL *restrict slice_times, const REAL *restrict g4dd_slices,
                                    const REAL *restrict gamma4udd_slices, const REAL t_target, REAL *restrict g4dd_out,
                                    REAL *restrict gamma4udd_out) {
  // Step 1: Build the shared 1D Lagrange basis in physical time.
  const int temporal_half_width = commondata->numerical_spacetime_temporal_interp_order;
  if (temporal_half_width < 0 || temporal_half_width > TEMPORAL_LAGRANGE_INTERP_MAX_HALF_WIDTH)
    return TEMPORAL_LAGRANGE_INTERP_INVALID_ORDER;
  const int interp_order = 2 * temporal_half_width + 1;
  REAL inv_denom[interp_order];
  REAL diffs_t[interp_order];
  REAL coeff_t[interp_order];
  const REAL normalization_time = interp_order > 1 ? pow(slice_times[1] - slice_times[0], -(interp_order - 1)) : 1.0;

  compute_inv_denom(interp_order, inv_denom);
  compute_diffs_xi(interp_order, t_target, slice_times, diffs_t);
  compute_lagrange_basis_coeffs_xi(interp_order, inv_denom, diffs_t, coeff_t);

  // Step 2: Interpolate the serialized metric components independently in time.
  for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
    REAL component_series[interp_order];

    for (int s = 0; s < interp_order; s++) {
      component_series[s] = g4dd_slices[s * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT + comp];
    } // END LOOP: for s over trusted metric slices for one component
    g4dd_out[comp] = normalization_time * sum_lagrange_x0_simd(interp_order, component_series, coeff_t);
  } // END LOOP: for comp over serialized metric components

  // Step 3: Interpolate the serialized Christoffel components independently in time.
  for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
    REAL component_series[interp_order];

    for (int s = 0; s < interp_order; s++) {
      component_series[s] = gamma4udd_slices[s * TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT + comp];
    } // END LOOP: for s over trusted Christoffel slices for one component
    gamma4udd_out[comp] = normalization_time * sum_lagrange_x0_simd(interp_order, component_series, coeff_t);
  } // END LOOP: for comp over serialized Christoffel components

  return TEMPORAL_LAGRANGE_INTERP_SUCCESS;
} // END FUNCTION: temporal_lagrange_interpolation
