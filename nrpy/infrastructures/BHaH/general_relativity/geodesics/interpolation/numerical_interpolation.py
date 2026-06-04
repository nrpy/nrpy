"""
Register chunk-based numerical-spacetime interpolation.

This module emits the host-side numerical interpolation wrapper used by
geodesic integrators. The generated C function mirrors the analytic
interpolation-kernel bundle contract: it consumes one chunk of photon states,
parallelizes over rays on the CPU, writes the 10-component metric bundle, and
writes the 40-component Christoffel bundle only when requested.

Operationally, this wrapper is the bridge between the data-management layer
and the two interpolation helpers. For each photon, it asks the numerical
time-window manager for the mapped temporal stencil, runs the spatial helper on
every slice in that stencil, and then runs the temporal helper once to recover
the final tensors at the photon coordinate time.

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
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.temporal_lagrange_interpolation import (
    register_CFunction_temporal_lagrange_interpolation,
)


def register_CFunction_numerical_interpolation(
    CoordSystem: str,
    enable_simd: bool,
    project_dir: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the CPU numerical-spacetime interpolation wrapper.

    This wrapper owns the per-photon orchestration of the numerical-spacetime
    interpolation pipeline. It does not select or map the active numerical
    window itself; instead, it assumes a caller already mapped a conservative
    slot-based window and then interpolates every ray in one chunk against that
    shared mapped data.

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
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> cfc.CFunction_dict.clear()
    >>> with tempfile.TemporaryDirectory(dir=os.getcwd()) as project_dir:
    ...     old_cache_home = os.environ.get("XDG_CACHE_HOME")
    ...     _ = os.environ.__setitem__("XDG_CACHE_HOME", project_dir)
    ...     with contextlib.redirect_stdout(io.StringIO()):
    ...         _ = register_CFunction_numerical_interpolation(
    ...             "Spherical", enable_simd=True, project_dir=project_dir
    ...         )
    ...         generated = clang_format(
    ...             cfc.CFunction_dict["numerical_interpolation"].full_function
    ...         )
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

    from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.time_window_manager_numerical import (  # pylint: disable=import-outside-toplevel
        time_window_manager_numerical,
    )
    from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.time_slot_manager_helpers import (  # pylint: disable=import-outside-toplevel
        time_slot_manager_helpers,
    )

    if "time_slot_manager" not in par.glb_extras_dict.get("BHaH_defines", {}):
        time_slot_manager_helpers()
    if "time_window_manager_numerical" not in par.glb_extras_dict.get(
        "BHaH_defines", {}
    ):
        time_window_manager_numerical()

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

The design goal is to let all photons in the chunk reuse the same mapped
numerical-spacetime payload window rather than loading numerical grids
independently ray-by-ray.

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
  // The mapped numerical window and the temporal helper must agree on the
  // centered temporal stencil width before any variable-length arrays are
  // sized from runtime data.
  if (temporal_half_width < 0 ||
      temporal_half_width > TEMPORAL_LAGRANGE_INTERP_MAX_HALF_WIDTH ||
      temporal_half_width != numerical_window->temporal_interp_half_width) {
    #pragma omp parallel for
    for (long int i = 0; i < chunk_size; i++) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      } // END LOOP: for comp over metric outputs after invalid temporal order
      if (d_connection_bundle != NULL) {
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
        } // END LOOP: for comp over connection outputs after invalid temporal order
      } // END IF: connection output bundle was requested
    } // END LOOP: for i over rays after invalid temporal order
    #undef IDX_F
    #undef IDX_METRIC
    #undef IDX_CONN
    return;
  } // END IF: runtime temporal interpolation half-width was invalid
  const int temporal_num_points = 2 * temporal_half_width + 1;
  if (temporal_num_points != numerical_window->temporal_interp_num_points) {
    #pragma omp parallel for
    for (long int i = 0; i < chunk_size; i++) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      } // END LOOP: for comp over metric outputs after inconsistent temporal stencil size
      if (d_connection_bundle != NULL) {
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
        } // END LOOP: for comp over connection outputs after inconsistent temporal stencil size
      } // END IF: connection output bundle was requested
    } // END LOOP: for i over rays after inconsistent temporal stencil size
    #undef IDX_F
    #undef IDX_METRIC
    #undef IDX_CONN
    return;
  } // END IF: runtime temporal stencil size did not match the mapped numerical window

  #pragma omp parallel for
  for (long int i = 0; i < chunk_size; i++) {
    const REAL t = (REAL)d_f_bundle[IDX_F(0, i)];
    const REAL x = (REAL)d_f_bundle[IDX_F(1, i)];
    const REAL y = (REAL)d_f_bundle[IDX_F(2, i)];
    const REAL z = (REAL)d_f_bundle[IDX_F(3, i)];
    int ray_failed = 0;
    uint64_t slice_indices[temporal_num_points];
    REAL slice_times[temporal_num_points];
    const double *slice_payloads[temporal_num_points];
    REAL g4dd_slices[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL gamma4udd_slices[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT];
    REAL g4dd_local[TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL gamma4udd_local[TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT];

    // Step 1: Recover the mapped temporal stencil for this photon from the
    // slot-level numerical window shared by the whole chunk.
    const int window_status = time_window_manager_numerical_stencil_for_time(
        numerical_window, (double)t, temporal_half_width, slice_indices,
        slice_times, slice_payloads);
    if (window_status != TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
      ray_failed = 1;
    } else {
      // Step 2: Interpolate each mapped time slice in space at the photon
      // position, producing one tensor bundle per temporal node.
      const int spatial_status =
          azimuthal_symmetry_spatial_lagrange_interpolation__rfm__Spherical(
              spatial_context, commondata, params, x, y, z, temporal_num_points,
              slice_payloads, g4dd_slices, gamma4udd_slices);
      if (spatial_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS) {
        ray_failed = 1;
      } else {
        // Step 3: Interpolate the per-slice tensor bundles in physical time to the
        // photon coordinate time.
        const int temporal_status = temporal_lagrange_interpolation(
            commondata, slice_times, g4dd_slices, gamma4udd_slices, t, g4dd_local,
            gamma4udd_local);
        if (temporal_status != TEMPORAL_LAGRANGE_INTERP_SUCCESS) {
          ray_failed = 1;
        } // END IF: temporal interpolation failed for this ray
      } // END ELSE: spatial interpolation succeeded for this ray
    } // END ELSE: mapped temporal stencil was available for this ray

    if (ray_failed) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      } // END LOOP: for comp over metric failure outputs
      if (d_connection_bundle != NULL) {
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
        } // END LOOP: for comp over connection failure outputs
      } // END IF: connection output bundle was requested
      continue;
    } // END IF: at least one interpolation stage failed for this ray

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
