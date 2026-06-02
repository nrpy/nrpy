"""
Register chunk-based numerical-spacetime interpolation.

This module emits the host-side numerical interpolation wrapper used by
geodesic integrators. The generated C function mirrors the analytic
interpolation-kernel bundle contract: it consumes one chunk of photon states,
parallelizes over rays on the CPU, writes the 10-component metric bundle, and
writes the 40-component Christoffel bundle only when requested.

The wrapper does not own file or mmap lifetime. A NumericalTimeWindowManager
must already have an active mapped time window for the slot being processed.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.azimuthal_symmetry_spatial_lagrange_interpolation import (
    register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.numerical_time_window_manager import (
    numerical_time_window_manager,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.temporal_lagrange_interpolation import (
    register_CFunction_temporal_lagrange_interpolation,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.time_slot_manager_helpers import (
    time_slot_manager_helpers,
)


def register_CFunction_numerical_interpolation(
    CoordSystem: str,
    enable_simd: bool,
    project_dir: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the CPU numerical-spacetime interpolation wrapper.

    :param CoordSystem: Coordinate system used by the mapped numerical dataset.
    :param enable_simd: Whether SIMD helper headers are already available.
    :param project_dir: Destination project directory for copied headers.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If `CoordSystem` is not supported.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import os
    >>> import tempfile
    >>> import nrpy.c_function as cfc
    >>> from nrpy.helpers.generic import validate_strings
    >>> cfc.CFunction_dict.clear()
    >>> with tempfile.TemporaryDirectory() as project_dir:
    ...     old_cache_home = os.environ.get("XDG_CACHE_HOME")
    ...     _ = os.environ.__setitem__("XDG_CACHE_HOME", project_dir)
    ...     with contextlib.redirect_stdout(io.StringIO()):
    ...         _ = register_CFunction_numerical_interpolation(
    ...             "Spherical", enable_simd=True, project_dir=project_dir
    ...         )
    ...         generated = cfc.CFunction_dict["numerical_interpolation"].full_function
    ...         _ = validate_strings(generated, "numerical_interpolation", file_ext="c")
    ...     if old_cache_home is None:
    ...         _ = os.environ.pop("XDG_CACHE_HOME", None)
    ...     else:
    ...         _ = os.environ.__setitem__("XDG_CACHE_HOME", old_cache_home)
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    if CoordSystem != "Spherical":
        raise ValueError(
            "numerical_interpolation currently supports only CoordSystem='Spherical'; "
            f"found '{CoordSystem}'."
        )

    if "time_slot_manager" not in par.glb_extras_dict.get("BHaH_defines", {}):
        time_slot_manager_helpers()
    if "numerical_time_window_manager" not in par.glb_extras_dict.get(
        "BHaH_defines", {}
    ):
        numerical_time_window_manager()

    spatial_name = "azimuthal_symmetry_spatial_lagrange_interpolation__rfm__Spherical"
    if spatial_name not in cfc.CFunction_dict:
        register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation(
            CoordSystem, enable_simd=enable_simd, project_dir=project_dir
        )
    if "temporal_lagrange_interpolation" not in cfc.CFunction_dict:
        register_CFunction_temporal_lagrange_interpolation(
            enable_simd=enable_simd, project_dir=project_dir
        )

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<stdint.h>",
    ]

    desc = r"""Interpolate numerical-spacetime tensors for one photon chunk.

The caller supplies an active numerical time window, a spatial interpolation
context, and one chunk of photon states in the same Structure-of-Arrays bundle
layout used by the analytic geodesic interpolation kernel. This CPU wrapper
parallelizes over rays, selects each photon's mapped temporal stencil, performs
spatial interpolation on every stencil slice, performs temporal interpolation
at the photon coordinate time, and writes the final metric and optional
Christoffel bundles.

@param[in] commondata Common runtime parameters.
@param[in] params Generated BHaH grid parameters for the mapped numerical data.
@param[in] spatial_context Trusted azimuthal-symmetry spatial interpolation context.
@param[in] numerical_window Active mapped numerical time-window manager.
@param[in] d_f_bundle Photon state bundle.
@param[out] d_metric_bundle Destination metric bundle.
@param[out] d_connection_bundle Destination Christoffel bundle, or NULL.
@param chunk_size Number of active rays in the chunk.
@param stream_idx Analytic-kernel compatibility argument; ignored on CPU.

@pre `numerical_window` must already map all slices needed by the active chunk.
"""
    cfunc_type = "void"
    name = "numerical_interpolation"
    params = """const commondata_struct *restrict commondata,
                const params_struct *restrict params,
                const azimuthal_symmetry_spatial_lagrange_context_struct *restrict spatial_context,
                const NumericalTimeWindowManager *restrict numerical_window,
                const double *restrict d_f_bundle,
                double *restrict d_metric_bundle,
                double *restrict d_connection_bundle,
                const long int chunk_size,
                const int stream_idx"""
    body = r"""
  (void)stream_idx;

  #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
  #define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
  #define IDX_CONN(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))

  const int temporal_half_width =
      commondata->numerical_spacetime_temporal_interp_order;
  const int temporal_num_points = 2 * temporal_half_width + 1;

  #pragma omp parallel for
  for (long int i = 0; i < chunk_size; i++) {
    const REAL t = (REAL)d_f_bundle[IDX_F(0, i)];
    const REAL x = (REAL)d_f_bundle[IDX_F(1, i)];
    const REAL y = (REAL)d_f_bundle[IDX_F(2, i)];
    const REAL z = (REAL)d_f_bundle[IDX_F(3, i)];
    uint64_t slice_indices[temporal_num_points];
    REAL slice_times[temporal_num_points];
    const double *slice_payloads[temporal_num_points];
    REAL g4dd_slices[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL gamma4udd_slices[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT];
    REAL g4dd_local[TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL gamma4udd_local[TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT];

    const int window_status = numerical_time_window_manager_stencil_for_time(
        numerical_window, (double)t, temporal_half_width, slice_indices,
        slice_times, slice_payloads);
    if (window_status != NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      } // END LOOP: for comp over metric failure outputs
      if (d_connection_bundle != NULL) {
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
        } // END LOOP: for comp over connection failure outputs
      } // END IF: connection output bundle was requested
      continue;
    } // END IF: temporal stencil was not inside the mapped numerical window

    const int spatial_status =
        azimuthal_symmetry_spatial_lagrange_interpolation__rfm__Spherical(
            spatial_context, commondata, params, x, y, z, temporal_num_points,
            slice_payloads, g4dd_slices, gamma4udd_slices);
    if (spatial_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      } // END LOOP: for comp over metric failure outputs
      if (d_connection_bundle != NULL) {
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
        } // END LOOP: for comp over connection failure outputs
      } // END IF: connection output bundle was requested
      continue;
    } // END IF: spatial interpolation failed for this ray

    temporal_lagrange_interpolation(
        commondata, slice_times, g4dd_slices, gamma4udd_slices, t, g4dd_local,
        gamma4udd_local);

    for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
      d_metric_bundle[IDX_METRIC(comp, i)] = (double)g4dd_local[comp];
    } // END LOOP: for comp over final metric components
    if (d_connection_bundle != NULL) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
        d_connection_bundle[IDX_CONN(comp, i)] = (double)gamma4udd_local[comp];
      } // END LOOP: for comp over final Christoffel components
    } // END IF: connection output bundle was requested
  } // END LOOP: for i over rays in chunk

  #undef IDX_F
  #undef IDX_METRIC
  #undef IDX_CONN
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
