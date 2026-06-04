"""
Register trusted temporal Lagrange interpolation.

This module emits the temporal stage of the numerical-spacetime interpolation
pipeline used by the geodesic integrators. The generated C API consumes flat
per-slice spatial-interpolation outputs for Cartesian-basis `g4DD` and
`Gamma4UDD` components at one fixed spatial point, plus the corresponding
physical slice times, and interpolates all 50 serialized tensor components to
one target coordinate time.

The helper deliberately assumes trusted inputs. In particular, callers must
provide physical `slice_times` that are finite and strictly increasing in
time, and they must supply exactly the same number of per-slice tensor
records as interpolation nodes. The metric bundle follows the same
10-component upper-triangular ordering consumed by the geodesic interpolation
kernels:
`g4DD00, g4DD01, g4DD02, g4DD03, g4DD11, g4DD12, g4DD13, g4DD22, g4DD23,
g4DD33`. The Christoffel bundle follows the same 40-component ordering used by
the geodesic kernels, with `alpha` outermost and `(mu, nu)` serialized in
upper-triangular order:
`Gamma4UDD000, Gamma4UDD001, ..., Gamma4UDD333`.
This module performs no file I/O, no stencil selection, no coordinate mapping,
no tensor rotation, and only a lightweight runtime check that the
requested interpolation half-width is within the supported range.
Higher-level code is responsible for selecting the temporal stencil and
for ensuring that the supplied slices already correspond to one fixed
spatial point.

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
#define TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT 40
#define TEMPORAL_LAGRANGE_INTERP_TENSOR_COMPONENT_COUNT \
  (TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT + \
   TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT)
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
    tensor bundle per mapped numerical time slice. It assumes the supplied
    `slice_times` are trusted, finite, strictly increasing, and contain
    exactly `2*n+1` entries, where `n` is
    `commondata->numerical_spacetime_temporal_interp_order`. The flat tensor
    bundles must contain one entry per supplied time node. The generated
    helper builds barycentric Lagrange coefficients directly from the exact
    physical slice times, so temporal nodes need not be uniformly spaced.

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

    desc = r"""Interpolate trusted Cartesian geodesic tensors in physical time.

The caller supplies flat per-slice spatial-interpolation outputs for the
configured temporal stencil, plus the corresponding physical `slice_times`.
This helper assumes those times are trusted and strictly increasing, derives
the actual number of time nodes from
`commondata->numerical_spacetime_temporal_interp_order`, builds one shared
nonuniform barycentric Lagrange basis in time from the exact physical
`slice_times`, and interpolates each serialized `g4DD` and `Gamma4UDD`
component independently to `t_target`.

The metric bundle ordering matches the geodesic interpolation-kernel contract:
`g4DD00, g4DD01, g4DD02, g4DD03, g4DD11, g4DD12, g4DD13, g4DD22, g4DD23,
g4DD33`.

The Christoffel bundle ordering also matches the geodesic interpolation-kernel
contract: `Gamma4UDD<alpha><mu><nu>` with `alpha` outermost and `(mu, nu)` in
upper-triangular order, i.e.
`Gamma4UDD000, Gamma4UDD001, Gamma4UDD002, Gamma4UDD003, Gamma4UDD011,
Gamma4UDD012, ..., Gamma4UDD333`.

@param[in] commondata Common runtime parameters.
@param[in] slice_times Trusted physical coordinate times for the supplied slices.
@param[in] g4dd_slices Flat per-slice metric components in kernel order.
@param[in] gamma4udd_slices Flat per-slice Christoffel components in kernel order.
@param[in] t_target Target physical coordinate time.
@param[out] g4dd_out Final interpolated metric components.
@param[out] gamma4udd_out Final interpolated Christoffel components.
@return Status code indicating success or invalid runtime interpolation order.

@note Callers must provide trusted, strictly increasing physical `slice_times`,
not abstract slot indices.
"""
    cfunc_type = "int"
    name = "temporal_lagrange_interpolation"
    params = """const commondata_struct *restrict commondata,
                const REAL *restrict slice_times,
                const REAL *restrict g4dd_slices,
                const REAL *restrict gamma4udd_slices,
                const REAL t_target,
                REAL *restrict g4dd_out,
                REAL *restrict gamma4udd_out"""
    body = r"""
  // Step 1: Build the shared 1D Lagrange basis in physical time.
  const int temporal_half_width =
      commondata->numerical_spacetime_temporal_interp_order;
  if (temporal_half_width < 0 ||
      temporal_half_width > TEMPORAL_LAGRANGE_INTERP_MAX_HALF_WIDTH)
    return TEMPORAL_LAGRANGE_INTERP_INVALID_ORDER;
  const int interp_order = 2 * temporal_half_width + 1;
  REAL barycentric_weights[interp_order];
  REAL normalized_slice_times[interp_order];
  REAL coeff_t[interp_order];
  const REAL time_origin = slice_times[0];
  const REAL time_scale =
      (interp_order > 1) ? (slice_times[interp_order - 1] - slice_times[0]) : 1.0;
  const REAL normalized_t_target = (t_target - time_origin) / time_scale;

  // Step 1.a: Normalize temporal nodes to one O(1) interval.
  for (int i = 0; i < interp_order; i++) {
    normalized_slice_times[i] = (slice_times[i] - time_origin) / time_scale;
  } // END LOOP: for i over temporal nodes while normalizing times

  // Step 1.b: Build barycentric interpolation weights after all nodes exist.
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
    for (int i = 0; i < interp_order; i++) {
      coeff_t[i] = (i == exact_node) ? 1.0 : 0.0;
    } // END LOOP: for i over temporal nodes while building an exact-match basis
  } else {
    REAL barycentric_sum = 0.0;

    for (int i = 0; i < interp_order; i++) {
      const REAL target_diff = normalized_t_target - normalized_slice_times[i];
      const REAL weighted_term = barycentric_weights[i] / target_diff;
      coeff_t[i] = weighted_term;
      barycentric_sum += weighted_term;
    } // END LOOP: for i over temporal nodes while summing barycentric terms

    for (int i = 0; i < interp_order; i++) {
      coeff_t[i] /= barycentric_sum;
    } // END LOOP: for i over temporal nodes while normalizing barycentric coefficients
  } // END ELSE: target time required a full barycentric basis evaluation

  // Step 2: Interpolate the serialized metric components independently in time.
  for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
    REAL component_series[interp_order];

    for (int s = 0; s < interp_order; s++) {
      component_series[s] =
          g4dd_slices[s * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT + comp];
    } // END LOOP: for s over trusted metric slices for one component
    g4dd_out[comp] = sum_lagrange_x0_simd(interp_order, component_series, coeff_t);
  } // END LOOP: for comp over serialized metric components

  // Step 3: Interpolate the serialized Christoffel components independently in time.
  for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
    REAL component_series[interp_order];

    for (int s = 0; s < interp_order; s++) {
      component_series[s] =
          gamma4udd_slices[s * TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT + comp];
    } // END LOOP: for s over trusted Christoffel slices for one component
    gamma4udd_out[comp] =
        sum_lagrange_x0_simd(interp_order, component_series, coeff_t);
  } // END LOOP: for comp over serialized Christoffel components

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
