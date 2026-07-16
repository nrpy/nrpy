"""
Register C1-constrained temporal metric interpolation near the initial slice.

This module emits the lower-temporal-boundary stage of the numerical-spacetime
interpolation pipeline used by the geodesic integrators. The generated C helper
accepts one complete, trusted temporal stencil of spatially interpolated
Cartesian-basis metric values and spatial metric derivatives. It reconstructs
the metric while preserving every supplied metric-node value and enforcing that
the reconstructed metric time derivative vanishes at ``t_metric_0``.

The constraint is imposed by correcting the ordinary temporal Lagrange
polynomial with a node-vanishing polynomial whose derivative is one at the
initial-time node. Applying the same correction to the spatial metric
derivative series ensures that all returned quantities remain derivatives of
one reconstructed four-metric. This helper is intended only for lower-boundary
stencils that contain ``t_metric_0``; ordinary temporal interpolation remains
the appropriate choice away from that boundary.

The module relies on the component-count macros and temporal interpolation
order CodeParameter registered by
``temporal_lagrange_interpolation_metric_derivatives``. The code-generation
caller must therefore register that ordinary temporal helper before this one.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
from nrpy.helpers.generic import copy_files

TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_DEFINES = r"""
// Status codes returned by the C1 startup temporal interpolation helper.
typedef enum {
  TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_SUCCESS = 0,
  TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_INVALID_ORDER = 1,
  TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_INVALID_INITIAL_TIME_NODE = 2
} temporal_lagrange_interp_c1_startup_status; // END ENUM: temporal_lagrange_interp_c1_startup_status
"""


def register_CFunction_temporal_lagrange_interpolation_c1_startup(
    enable_simd: bool,
    project_dir: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C1-constrained startup temporal interpolation helper.

    The generated helper accepts the same flat metric and metric-derivative
    bundles as the ordinary temporal helper, but it additionally receives the
    initial coordinate time ``t_metric_0``. The supplied strictly increasing
    stencil times must contain this initial time exactly once. The helper
    first forms ordinary barycentric Lagrange value and derivative bases, then
    subtracts a polynomial correction that is zero at every supplied node and
    has unit physical-time derivative at ``t_metric_0``. This preserves every
    supplied node value while enforcing zero metric time derivatives at the
    initial-time node.

    The correction is applied consistently to the metric and its Cartesian
    spatial derivative series. Consequently, the output ``g4DD_dD`` bundle
    remains the derivative of the corrected metric reconstruction rather than
    a postprocessed collection of unrelated derivatives.

    :param enable_simd: Whether SIMD helper headers are already available.
    :param project_dir: Destination project directory for copied headers.
    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    Bdefines_h.register_BHaH_defines(
        "temporal_lagrange_interpolation_c1_startup",
        TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_DEFINES,
    )

    # Step 1: Ensure required helper headers are copied into the project.
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

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "interpolation_lagrange_uniform.h",
    ]

    desc = r"""Interpolate metric data with a C1 lower-time startup constraint.

The caller supplies one trusted, strictly increasing temporal stencil of flat
per-slice spatial-interpolation outputs, the target coordinate time, and the
initial numerical-spacetime time ``t_metric_0``. The stencil must contain
``t_metric_0`` exactly once. This helper preserves every supplied metric value
at its corresponding time node and enforces that the reconstructed metric has
zero coordinate-time derivative at ``t_metric_0``.

The helper first builds ordinary nonuniform barycentric Lagrange bases. It then
subtracts the metric's derivative at ``t_metric_0`` times a correction
polynomial that vanishes at every supplied node and has unit physical-time
derivative at ``t_metric_0``. The same correction is applied to each Cartesian
spatial metric-derivative series, preserving consistency between the metric
and all returned first derivatives.

The metric bundle ordering is:
`g4DD00, g4DD01, g4DD02, g4DD03, g4DD11, g4DD12, g4DD13, g4DD22, g4DD23,
g4DD33`.

The 40-component `gamma4udd`-named bundle stores metric derivatives with
metric-pair outermost and derivative direction innermost:
`g4DD_dD000, g4DD_dD001, g4DD_dD002, g4DD_dD003, g4DD_dD010, ...,
g4DD_dD333`. The final index is `(t, x, y, z) = (0, 1, 2, 3)`.

@param[in] commondata Common runtime parameters.
@param[in] slice_times Trusted, finite, strictly increasing physical coordinate
times for the supplied slices.
@param[in] g4dd_slices Flat per-slice metric components in kernel order.
@param[in] gamma4udd_slices Flat per-slice metric derivatives in temporary
connection-bundle storage order.
@param t_metric_0 Initial numerical-spacetime coordinate time in the stencil.
@param t_target Target physical coordinate time at or above `t_metric_0`.
@param[out] g4dd_out Corrected interpolated metric components.
@param[out] gamma4udd_out Corrected metric first derivatives in temporary
connection-bundle storage order.
@return `TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_SUCCESS` on success,
`TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_INVALID_ORDER` for an unsupported order,
or `TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_INVALID_INITIAL_TIME_NODE` when the
stencil does not contain exactly one initial-time node.

@pre The ordinary temporal interpolation registration has already defined the
shared component-count macros and temporal interpolation-order CodeParameter.
"""
    cfunc_type = "int"
    name = "temporal_lagrange_interpolation_c1_startup"
    params = """const commondata_struct *restrict commondata,
                const REAL *restrict slice_times,
                const REAL *restrict g4dd_slices,
                const REAL *restrict gamma4udd_slices,
                const REAL t_metric_0,
                const REAL t_target,
                REAL *restrict g4dd_out,
                REAL *restrict gamma4udd_out"""
    body = r"""
  // Step 1: Validate the configured stencil size and locate the initial-time
  // node used to impose the C1 startup constraint.
  const int temporal_half_width =
      commondata->numerical_spacetime_temporal_interp_order;
  if (temporal_half_width < 0 ||
      temporal_half_width > TEMPORAL_LAGRANGE_INTERP_MAX_HALF_WIDTH)
    return TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_INVALID_ORDER;
  const int interp_order = 2 * temporal_half_width + 1;
  const REAL time_origin = t_metric_0;
  const REAL time_scale =
      (interp_order > 1) ?
          (slice_times[interp_order - 1] - slice_times[0]) :
          1.0;
  if (!(time_scale > 0.0) || t_target < t_metric_0)
    return TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_INVALID_INITIAL_TIME_NODE;

  int initial_time_node = -1;
  for (int i = 0; i < interp_order; i++) {
    if (slice_times[i] == t_metric_0) {
      if (initial_time_node >= 0)
        return TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_INVALID_INITIAL_TIME_NODE;
      initial_time_node = i;
    } // END IF: one stencil time matched the initial numerical time
  } // END LOOP: for i over temporal stencil nodes while locating the initial time
  if (initial_time_node < 0)
    return TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_INVALID_INITIAL_TIME_NODE;

  // Step 2: Normalize temporal coordinates about t_metric_0 and build the
  // barycentric weights shared by all metric and derivative components.
  REAL barycentric_weights[interp_order];
  REAL normalized_slice_times[interp_order];
  REAL coeff_t[interp_order];
  REAL coeff_dt_target_normalized[interp_order];
  REAL coeff_dt_initial_normalized[interp_order];
  const REAL normalized_t_target = (t_target - time_origin) / time_scale;

  for (int i = 0; i < interp_order; i++)
    normalized_slice_times[i] =
        (slice_times[i] - time_origin) / time_scale;

  for (int i = 0; i < interp_order; i++) {
    REAL weight_denom = 1.0;
    for (int j = 0; j < interp_order; j++) {
      if (j != i) {
        const REAL time_difference =
            normalized_slice_times[i] - normalized_slice_times[j];
        if (time_difference == 0.0)
          return TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_INVALID_INITIAL_TIME_NODE;
        weight_denom *= time_difference;
      } // END IF: accumulating one barycentric weight denominator factor
    } // END LOOP: for j over temporal nodes while building one barycentric weight
    barycentric_weights[i] = 1.0 / weight_denom;
  } // END LOOP: for i over temporal nodes while building barycentric weights

  // Step 3: Build value and derivative bases at the target time.
  int target_node = -1;
  for (int i = 0; i < interp_order; i++) {
    if (normalized_t_target == normalized_slice_times[i]) {
      target_node = i;
      break;
    } // END IF: target time exactly matched one supplied node
  } // END LOOP: for i over temporal nodes while checking the target time

  if (target_node >= 0) {
    REAL target_diagonal_derivative = 0.0;
    for (int i = 0; i < interp_order; i++) {
      coeff_t[i] = (i == target_node) ? 1.0 : 0.0;
      if (i == target_node) {
        coeff_dt_target_normalized[i] = 0.0;
      } else {
        coeff_dt_target_normalized[i] =
            barycentric_weights[i] /
            (barycentric_weights[target_node] *
             (normalized_slice_times[target_node] - normalized_slice_times[i]));
        target_diagonal_derivative -= coeff_dt_target_normalized[i];
      } // END ELSE: one off-diagonal target-node derivative coefficient
    } // END LOOP: for i over exact-target value and derivative coefficients
    coeff_dt_target_normalized[target_node] = target_diagonal_derivative;
  } else {
    REAL barycentric_sum = 0.0;
    REAL barycentric_inverse_square_sum = 0.0;

    for (int i = 0; i < interp_order; i++) {
      const REAL target_difference =
          normalized_t_target - normalized_slice_times[i];
      const REAL weighted_term = barycentric_weights[i] / target_difference;
      coeff_t[i] = weighted_term;
      barycentric_sum += weighted_term;
      barycentric_inverse_square_sum += weighted_term / target_difference;
    } // END LOOP: for i over temporal nodes while summing barycentric terms

    for (int i = 0; i < interp_order; i++) {
      const REAL target_difference =
          normalized_t_target - normalized_slice_times[i];
      coeff_t[i] /= barycentric_sum;
      coeff_dt_target_normalized[i] =
          coeff_t[i] *
          (barycentric_inverse_square_sum / barycentric_sum -
           1.0 / target_difference);
    } // END LOOP: for i over target-time value and derivative coefficients
  } // END ELSE: target time required a full barycentric basis evaluation

  // Step 4: Build the ordinary interpolation derivative basis at t_metric_0.
  REAL initial_diagonal_derivative = 0.0;
  for (int i = 0; i < interp_order; i++) {
    if (i == initial_time_node) {
      coeff_dt_initial_normalized[i] = 0.0;
    } else {
      coeff_dt_initial_normalized[i] =
          barycentric_weights[i] /
          (barycentric_weights[initial_time_node] *
           (normalized_slice_times[initial_time_node] -
            normalized_slice_times[i]));
      initial_diagonal_derivative -= coeff_dt_initial_normalized[i];
    } // END ELSE: one off-diagonal initial-time derivative coefficient
  } // END LOOP: for i over initial-time derivative coefficients
  coeff_dt_initial_normalized[initial_time_node] = initial_diagonal_derivative;

  // Step 5: Evaluate the node-vanishing correction and its derivative. In
  // normalized time u, r_hat(u) vanishes at every stencil node and has
  // r_hat'(0) = 1. The physical correction is time_scale * r_hat.
  REAL correction_numerator = 1.0;
  REAL correction_denominator = 1.0;
  for (int i = 0; i < interp_order; i++) {
    if (i != initial_time_node) {
      correction_numerator *=
          normalized_t_target - normalized_slice_times[i];
      correction_denominator *= -normalized_slice_times[i];
    } // END IF: one non-initial node contributed to the correction polynomial
  } // END LOOP: for i over non-initial nodes while building the correction
  correction_numerator *= normalized_t_target;
  if (correction_denominator == 0.0)
    return TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_INVALID_INITIAL_TIME_NODE;
  const REAL correction_value = correction_numerator / correction_denominator;

  REAL correction_derivative_numerator = 0.0;
  for (int omitted_factor = -1; omitted_factor < interp_order;
       omitted_factor++) {
    if (omitted_factor == initial_time_node)
      continue;
    REAL product_except_one_factor = 1.0;
    if (omitted_factor != -1)
      product_except_one_factor *= normalized_t_target;
    for (int i = 0; i < interp_order; i++) {
      if (i != initial_time_node && i != omitted_factor)
        product_except_one_factor *=
            normalized_t_target - normalized_slice_times[i];
    } // END LOOP: for i over correction factors other than the omitted factor
    correction_derivative_numerator += product_except_one_factor;
  } // END LOOP: for omitted_factor over correction-polynomial derivative terms
  const REAL correction_derivative =
      correction_derivative_numerator / correction_denominator;

  // Step 6: Correct each metric series so the reconstructed time derivative
  // vanishes at t_metric_0 while retaining every supplied metric value.
  for (int metric_component = 0;
       metric_component < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT;
       metric_component++) {
    REAL metric_series[interp_order];
    for (int s = 0; s < interp_order; s++)
      metric_series[s] =
          g4dd_slices[
              s * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT + metric_component];

    const REAL ordinary_value =
        sum_lagrange_x0_simd(interp_order, metric_series, coeff_t);
    const REAL ordinary_target_derivative =
        sum_lagrange_x0_simd(
            interp_order, metric_series, coeff_dt_target_normalized) /
        time_scale;
    const REAL ordinary_initial_derivative =
        sum_lagrange_x0_simd(
            interp_order, metric_series, coeff_dt_initial_normalized) /
        time_scale;

    g4dd_out[metric_component] =
        ordinary_value - ordinary_initial_derivative * time_scale * correction_value;
    const int temporal_derivative_slot =
        TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_DIRECTION_COUNT * metric_component +
        TEMPORAL_LAGRANGE_INTERP_TEMPORAL_DERIVATIVE_INDEX;
    gamma4udd_out[temporal_derivative_slot] =
        ordinary_target_derivative -
        ordinary_initial_derivative * correction_derivative;
  } // END LOOP: for metric_component over metric temporal reconstructions

  // Step 7: Apply the same correction to the Cartesian spatial-derivative
  // series, preserving their relationship to the corrected metric.
  for (int metric_component = 0;
       metric_component < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT;
       metric_component++) {
    for (int derivative_direction = 1;
         derivative_direction <
             TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_DIRECTION_COUNT;
         derivative_direction++) {
      REAL spatial_derivative_series[interp_order];
      const int derivative_slot =
          TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_DIRECTION_COUNT *
              metric_component +
          derivative_direction;
      for (int s = 0; s < interp_order; s++)
        spatial_derivative_series[s] =
            gamma4udd_slices[
                s * TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT +
                derivative_slot];

      const REAL ordinary_value = sum_lagrange_x0_simd(
          interp_order, spatial_derivative_series, coeff_t);
      const REAL ordinary_initial_time_derivative =
          sum_lagrange_x0_simd(
              interp_order, spatial_derivative_series,
              coeff_dt_initial_normalized) /
          time_scale;
      gamma4udd_out[derivative_slot] =
          ordinary_value - ordinary_initial_time_derivative * time_scale *
              correction_value;
    } // END LOOP: for derivative_direction over Cartesian spatial directions
  } // END LOOP: for metric_component over spatial-derivative reconstructions

  return TEMPORAL_LAGRANGE_INTERP_C1_STARTUP_SUCCESS;
"""

    cfc.register_CFunction(
        subdirectory="interpolation",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
