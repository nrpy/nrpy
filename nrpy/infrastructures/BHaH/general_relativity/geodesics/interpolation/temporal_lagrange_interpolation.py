"""
Register trusted temporal Lagrange interpolation.

This module emits the temporal stage of the numerical-spacetime interpolation
pipeline used by the geodesic integrators. The generated C API consumes flat
per-slice spatial-interpolation outputs for Cartesian-basis `g4DD` and metric
first derivatives at one fixed spatial point, plus the corresponding physical
slice times. It interpolates the 10 metric values and 30 Cartesian spatial
metric derivatives to one target coordinate time, and analytically
differentiates the same temporal metric reconstruction to obtain the 10
coordinate-time derivatives.

The helper deliberately assumes trusted inputs. In particular, callers must
provide physical `slice_times` that are finite and strictly increasing in
time, and they must supply exactly the same number of per-slice tensor records
as interpolation nodes. The metric bundle follows the same 10-component
upper-triangular ordering consumed by the geodesic interpolation kernels:
`g4DD00, g4DD01, g4DD02, g4DD03, g4DD11, g4DD12, g4DD13, g4DD22, g4DD23,
g4DD33`.

The 40-component `g4DD_dD` bundle stores metric first derivatives. Its
ordering is metric-pair outermost and
derivative direction innermost:
`g4DD_dD000, g4DD_dD001, g4DD_dD002, g4DD_dD003, g4DD_dD010, ...,
g4DD_dD333`, where the final index is `(t, x, y, z) = (0, 1, 2, 3)`.
The spatial stage supplies zero in each temporal-derivative slot and supplies
the 30 Cartesian spatial derivatives. This temporal helper replaces the 10
temporal slots with derivatives of the metric's temporal Lagrange interpolant
and interpolates the 30 spatial derivative series normally in time.

This module performs no file I/O, no stencil selection, no coordinate mapping,
and no tensor rotation. Higher-level code is responsible for selecting the
temporal stencil and ensuring that the supplied slices already correspond to
one fixed spatial point.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par
from nrpy.helpers.generic import copy_files

TEMPORAL_LAGRANGE_INTERP_DEFINES = r"""
// Constants for temporal Lagrange interpolation.
#define TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT 10
#define TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT 40
#define TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_DIRECTION_COUNT 4
#define TEMPORAL_LAGRANGE_INTERP_TEMPORAL_DERIVATIVE_INDEX 0
#define TEMPORAL_LAGRANGE_INTERP_TENSOR_COMPONENT_COUNT \
  (TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT + \
   TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT)
#define TEMPORAL_LAGRANGE_INTERP_MAX_HALF_WIDTH 32

// Status codes returned by the temporal Lagrange interpolation helper.
typedef enum {
  TEMPORAL_LAGRANGE_INTERP_SUCCESS = 0,
  TEMPORAL_LAGRANGE_INTERP_INVALID_ORDER = 1
} temporal_lagrange_interp_status; // END ENUM: temporal_lagrange_interp_status
"""


def register_CFunction_temporal_lagrange_interpolation(
    enable_simd: bool,
    project_dir: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the trusted temporal Lagrange interpolation helper.

    The generated C helper assumes the caller already ran the spatial
    interpolation stage for the desired temporal stencil and now wants a
    lightweight 1D interpolation in physical coordinate time. In the full
    pipeline, this helper is called after the spatial helper has produced one
    metric and metric-derivative bundle per mapped numerical time slice. It
    assumes the supplied `slice_times` are trusted, finite, strictly
    increasing, and contain exactly `2*n+1` entries, where `n` is
    `commondata->numerical_spacetime_temporal_interp_order`. The flat bundles
    must contain one entry per supplied time node.

    The generated helper builds barycentric Lagrange coefficients directly
    from the exact physical slice times, so temporal nodes need not be
    uniformly spaced. It also builds analytic derivative coefficients for the
    same barycentric reconstruction. These coefficients are applied to the 10
    metric series to obtain their physical coordinate-time derivatives. The
    30 Cartesian spatial metric-derivative series are interpolated normally in
    time and retain the metric-pair-outermost, derivative-direction-innermost
    packing used by the spatial helper.

    :param enable_simd: Whether SIMD helper headers are already available.
    :param project_dir: Destination project directory for copied headers.
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import os
    >>> import tempfile
    >>> import nrpy.c_function as cfc
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> cfc.CFunction_dict.clear()
    >>> with tempfile.TemporaryDirectory(dir=os.getcwd()) as project_dir:
    ...     old_cache_home = os.environ.get("XDG_CACHE_HOME")
    ...     _ = os.environ.__setitem__("XDG_CACHE_HOME", project_dir)
    ...     with contextlib.redirect_stdout(io.StringIO()):
    ...         _ = register_CFunction_temporal_lagrange_interpolation(
    ...             enable_simd=True, project_dir=project_dir
    ...         )
    ...         generated = clang_format(
    ...             cfc.CFunction_dict["temporal_lagrange_interpolation"].full_function
    ...         )
    ...         _ = validate_strings(
    ...             generated, "temporal_lagrange_interpolation", file_ext="c"
    ...         )
    ...     if old_cache_home is None:
    ...         _ = os.environ.pop("XDG_CACHE_HOME", None)
    ...     else:
    ...         _ = os.environ.__setitem__("XDG_CACHE_HOME", old_cache_home)
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    _ = par.register_CodeParameter(
        "int",
        __name__,
        "numerical_spacetime_temporal_interp_order",
        2,
        commondata=True,
        add_to_parfile=True,
        description=(
            "Centered temporal Lagrange interpolation half-width n; the generated "
            "helper uses 2*n+1 mapped time slices."
        ),
    )
    Bdefines_h.register_BHaH_defines(
        "temporal_lagrange_interpolation", TEMPORAL_LAGRANGE_INTERP_DEFINES
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

    desc = r"""Interpolate metric values and first derivatives in physical time.

The caller supplies flat per-slice spatial-interpolation outputs for the
configured temporal stencil, plus the corresponding physical `slice_times`.
This helper assumes those times are trusted, finite, and strictly increasing,
derives the actual number of time nodes from
`commondata->numerical_spacetime_temporal_interp_order`, and builds one shared
nonuniform barycentric Lagrange basis from the exact physical slice times.

The 10 serialized `g4DD` components are interpolated independently to
`t_target`. Analytic derivatives of that same temporal interpolation basis are
applied to the metric series to obtain the 10 coordinate-time metric
derivatives. The 30 Cartesian spatial metric derivatives supplied by the
spatial stage are interpolated independently in time.

The metric bundle ordering is:
`g4DD00, g4DD01, g4DD02, g4DD03, g4DD11, g4DD12, g4DD13, g4DD22, g4DD23,
g4DD33`.

The 40-component `g4dd_dD` bundle stores metric derivatives with metric-pair
outermost and derivative direction innermost:
`g4DD_dD000, g4DD_dD001, g4DD_dD002, g4DD_dD003, g4DD_dD010, ...,
g4DD_dD333`. The final index is `(t, x, y, z) = (0, 1, 2, 3)`. Thus bundle
slot `4*p` receives the temporal derivative of metric component `p`, while
slots `4*p+1`, `4*p+2`, and `4*p+3` receive the temporally interpolated
Cartesian spatial derivatives.

@param[in] commondata Common runtime parameters.
@param[in] slice_times Trusted, finite, strictly increasing physical coordinate
times for the supplied slices.
@param[in] g4dd_slices Flat per-slice metric components in kernel order.
@param[in] g4dd_dD_slices Flat per-slice metric derivatives in metric-pair
outermost storage order.
@param[in] t_target Target physical coordinate time.
@param[out] g4dd_out Final interpolated metric components.
@param[out] g4dd_dD_out Final metric first derivatives in metric-pair
outermost storage order.
@return `TEMPORAL_LAGRANGE_INTERP_SUCCESS` on success, or
`TEMPORAL_LAGRANGE_INTERP_INVALID_ORDER` if the configured interpolation
half-width is outside the supported range.

@note Callers must provide trusted, finite, strictly increasing physical
`slice_times`, not abstract slot indices.
"""
    cfunc_type = "int"
    name = "temporal_lagrange_interpolation"
    params = """const commondata_struct *restrict commondata,
                const REAL *restrict slice_times,
                const REAL *restrict g4dd_slices,
                const REAL *restrict g4dd_dD_slices,
                const REAL t_target,
                REAL *restrict g4dd_out,
                REAL *restrict g4dd_dD_out"""
    body = r"""
  // Step 1: Build the shared nonuniform barycentric Lagrange basis and its
  // analytic derivative in normalized physical time.
  const int temporal_half_width =
      commondata->numerical_spacetime_temporal_interp_order;
  if (temporal_half_width < 0 ||
      temporal_half_width > TEMPORAL_LAGRANGE_INTERP_MAX_HALF_WIDTH)
    return TEMPORAL_LAGRANGE_INTERP_INVALID_ORDER;
  const int interp_order = 2 * temporal_half_width + 1;
  REAL barycentric_weights[interp_order];
  REAL normalized_slice_times[interp_order];
  REAL coeff_t[interp_order];
  REAL coeff_dt_normalized[interp_order];
  const REAL time_origin = slice_times[0];
  const REAL time_scale =
      (interp_order > 1) ? (slice_times[interp_order - 1] - slice_times[0]) : 1.0;
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
        const REAL time_diff =
            normalized_slice_times[i] - normalized_slice_times[j];
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
            barycentric_weights[i] /
            (barycentric_weights[exact_node] *
             (normalized_slice_times[exact_node] - normalized_slice_times[i]));
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
      coeff_dt_normalized[i] =
          coeff_t[i] *
          (barycentric_inverse_square_sum / barycentric_sum - 1.0 / target_diff);
    } // END LOOP: for i over normalized value and derivative coefficients
  } // END ELSE: target time required a full barycentric basis evaluation

  // Step 2: Interpolate each metric component and differentiate the same
  // temporal reconstruction with respect to physical coordinate time.
  for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
    REAL component_series[interp_order];

    for (int s = 0; s < interp_order; s++) {
      component_series[s] =
          g4dd_slices[s * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT + comp];
    } // END LOOP: for s over trusted metric slices for one component
    g4dd_out[comp] = sum_lagrange_x0_simd(interp_order, component_series, coeff_t);

    const int temporal_derivative_slot =
        TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_DIRECTION_COUNT * comp +
        TEMPORAL_LAGRANGE_INTERP_TEMPORAL_DERIVATIVE_INDEX;
    g4dd_dD_out[temporal_derivative_slot] =
        sum_lagrange_x0_simd(
            interp_order, component_series, coeff_dt_normalized) /
        time_scale;
  } // END LOOP: for comp over serialized metric components

  // Step 3: Interpolate only the 30 Cartesian spatial metric derivatives.
  // The derivative direction is the fastest-changing bundle index, so each
  // metric component occupies four consecutive slots ordered as (t, x, y, z).
  for (int metric_comp = 0;
       metric_comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT;
       metric_comp++) {
    for (int derivative_direction = 1;
         derivative_direction <
             TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_DIRECTION_COUNT;
         derivative_direction++) {
      REAL component_series[interp_order];
      const int derivative_slot =
          TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_DIRECTION_COUNT * metric_comp +
          derivative_direction;

      for (int s = 0; s < interp_order; s++) {
        component_series[s] =
            g4dd_dD_slices[
                s * TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT +
                derivative_slot];
      } // END LOOP: for s over trusted spatial-derivative slices for one component
      g4dd_dD_out[derivative_slot] =
          sum_lagrange_x0_simd(interp_order, component_series, coeff_t);
    } // END LOOP: for derivative_direction over Cartesian spatial directions
  } // END LOOP: for metric_comp over serialized metric components

  return TEMPORAL_LAGRANGE_INTERP_SUCCESS;
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
