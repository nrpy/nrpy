"""
Register chunk-based numerical-spacetime interpolation.

This module emits the host-side numerical interpolation wrapper used by
geodesic integrators. The generated C function mirrors the analytic
interpolation-kernel bundle contract: it consumes one chunk of photon states,
parallelizes over rays on the CPU, writes the 10-component metric bundle, and
writes the 40-component Christoffel bundle only when requested.

Operationally, this wrapper is the bridge between the data-management layer
and the analytic/numerical interpolation helpers. For each photon, it either
evaluates one of two static analytic metrics directly or asks the numerical
time-window manager for one adaptive centered temporal stencil. In the mixed
case, the wrapper spatially interpolates only the mapped numerical stencil
subset, fills the one missing stencil edge from the appropriate static
analytic metric, and then runs the temporal helper once on the reconstructed
full stencil to recover the final tensors at the photon coordinate time.

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
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.time_window_manager_numerical import (
    time_window_manager_numerical,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.time_slot_manager_helpers import (
    time_slot_manager_helpers,
)


def register_CFunction_numerical_interpolation(
    CoordSystem: str,
    enable_simd: bool,
    project_dir: str,
    analytical_metric_0: str,
    analytical_metric_1: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the CPU numerical-spacetime interpolation wrapper.

    This wrapper owns the per-photon orchestration of the numerical-spacetime
    interpolation pipeline. It does not select or map the active numerical
    window itself; instead, it assumes a caller already mapped a conservative
    slot-based window and then evaluates every ray in one chunk against that
    shared mapped data. Rays strictly below `t_metric_0` or strictly above
    `t_metric_1` are handled directly by the requested static analytic
    metrics. Rays inside the numerical interval use the adaptive
    `time_window_manager_numerical_stencil_for_time()` contract from
    `time_window_manager_numerical`, spatially interpolate only the mapped
    numerical stencil subset, fill the one missing stencil edge analytically,
    and then perform temporal interpolation on the reconstructed full stencil.

    :param CoordSystem: Coordinate system used by the mapped numerical dataset;
        must be `"Spherical"`, `"SinhSpherical"`, `"Cylindrical"`, or
        `"SinhCylindrical"`.
    :param enable_simd: Whether SIMD helper headers are already available.
    :param project_dir: Destination project directory for copied headers.
    :param analytical_metric_0: Name of the lower-time static analytic metric
        used for photon times strictly below `t_metric_0`.
    :param analytical_metric_1: Name of the upper-time static analytic metric
        used for photon times strictly above `t_metric_1`.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If `CoordSystem` is not supported or if either
        analytic metric name is empty.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import os
    >>> import tempfile
    >>> import nrpy.c_function as cfc
    >>> from nrpy.equations.general_relativity.geodesics import geodesics as geo
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> from nrpy.infrastructures.BHaH.general_relativity.geodesics import connections, g4DD_metric
    >>> cfc.CFunction_dict.clear()
    >>> with tempfile.TemporaryDirectory(dir=os.getcwd()) as project_dir:
    ...     old_cache_home = os.environ.get("XDG_CACHE_HOME")
    ...     _ = os.environ.__setitem__("XDG_CACHE_HOME", project_dir)
    ...     with contextlib.redirect_stdout(io.StringIO()):
    ...         metric_data = geo.Geodesic_Equations["KerrSchild_Cartesian_photon"]
    ...         g4DD_metric.g4DD_metric(
    ...             metric_data.g4DD, "KerrSchild_Cartesian", "photon"
    ...         )
    ...         connections.connections(
    ...             metric_data.Gamma4UDD, "KerrSchild_Cartesian", "photon"
    ...         )
    ...         _ = register_CFunction_numerical_interpolation(
    ...             "Spherical",
    ...             enable_simd=True,
    ...             project_dir=project_dir,
    ...             analytical_metric_0="KerrSchild_Cartesian",
    ...             analytical_metric_1="KerrSchild_Cartesian",
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

    if CoordSystem not in (
        "Spherical",
        "SinhSpherical",
        "Cylindrical",
        "SinhCylindrical",
    ):
        raise ValueError(
            "numerical_interpolation currently supports only CoordSystem in "
            "('Spherical', 'SinhSpherical', 'Cylindrical', 'SinhCylindrical'); "
            f"found '{CoordSystem}'."
        )
    if not analytical_metric_0 or not analytical_metric_1:
        raise ValueError("Analytic metric names must be non-empty strings.")

    if "time_slot_manager" not in par.glb_extras_dict.get("BHaH_defines", {}):
        time_slot_manager_helpers()
    if "time_window_manager_numerical" not in par.glb_extras_dict.get(
        "BHaH_defines", {}
    ):
        time_window_manager_numerical()
    _ = par.register_CodeParameter(
        "REAL",
        __name__,
        "t_metric_0",
        0.0,
        commondata=True,
        add_to_parfile=True,
        description=(
            "Lower coordinate-time transition from analytical_metric_0 to the "
            "numerical spacetime interpolation region."
        ),
    )
    _ = par.register_CodeParameter(
        "REAL",
        __name__,
        "t_metric_1",
        1.0,
        commondata=True,
        add_to_parfile=True,
        description=(
            "Upper coordinate-time transition from the numerical spacetime "
            "interpolation region to analytical_metric_1."
        ),
    )

    spatial_name = (
        "azimuthal_symmetry_spatial_lagrange_interpolation__rfm__" f"{CoordSystem}"
    )
    if spatial_name not in cfc.CFunction_dict:
        register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation(
            CoordSystem, enable_simd=enable_simd, project_dir=project_dir
        )
    if "temporal_lagrange_interpolation" not in cfc.CFunction_dict:
        register_CFunction_temporal_lagrange_interpolation(
            enable_simd=enable_simd, project_dir=project_dir
        )
    metric_worker_0 = f"g4DD_metric_{analytical_metric_0}"
    metric_worker_1 = f"g4DD_metric_{analytical_metric_1}"
    conn_worker_0 = f"connections_{analytical_metric_0}"
    conn_worker_1 = f"connections_{analytical_metric_1}"
    prefunc_lines = []
    for worker_name in [metric_worker_0, metric_worker_1, conn_worker_0, conn_worker_1]:
        if worker_name not in cfc.CFunction_dict:
            raise ValueError(
                f"Analytic worker '{worker_name}' must be registered before "
                "register_CFunction_numerical_interpolation() is called."
            )
        worker_code = cfc.CFunction_dict[worker_name].full_function
        if worker_code not in prefunc_lines:
            prefunc_lines.append(worker_code)
    prefunc = "\n\n".join(prefunc_lines)

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<stdint.h>",
    ]

    desc = rf"""Interpolate piecewise analytic/numerical spacetime tensors for one photon chunk.

The caller supplies an active numerical time window, a spatial interpolation
context, and one chunk of photon states in the same Structure-of-Arrays bundle
layout used by the analytic geodesic interpolation kernel. This CPU wrapper
parallelizes over rays, evaluates `{analytical_metric_0}` directly for times
below `t_metric_0`, evaluates `{analytical_metric_1}` directly for times above
`t_metric_1`, and otherwise reconstructs one mixed temporal stencil by
combining spatial interpolation on the mapped numerical slices with one static
analytic fill on the missing stencil edge before performing temporal
interpolation at the photon coordinate time.

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
    body = (
        r"""
  (void)stream_idx;

  #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
  #define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
  #define IDX_CONN(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
  #define G4_SLICE(base_ptr, slice_idx, comp_idx) \
    ((base_ptr)[(slice_idx) * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT + (comp_idx)])
  #define GAMMA_SLICE(base_ptr, slice_idx, comp_idx) \
    ((base_ptr)[(slice_idx) * TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT + (comp_idx)])

  const int temporal_half_width =
      commondata->numerical_spacetime_temporal_interp_order;
  const REAL t_metric_0 = (REAL)commondata->t_metric_0;
  const REAL t_metric_1 = (REAL)commondata->t_metric_1;
  // The mapped numerical window and the temporal helper must agree on the
  // centered temporal stencil width before any variable-length arrays are
  // sized from runtime data.
  if (temporal_half_width < 0 ||
      temporal_half_width > TEMPORAL_LAGRANGE_INTERP_MAX_HALF_WIDTH ||
      temporal_half_width != numerical_window->temporal_interp_half_width ||
      !isfinite((double)t_metric_0) || !isfinite((double)t_metric_1) ||
      t_metric_0 >= t_metric_1) {
    #pragma omp parallel for
    for (long int i = 0; i < chunk_size; i++) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      if (d_connection_bundle != NULL)
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++)
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
    } // END LOOP: for i over rays after invalid temporal order
    return;
  } // END IF: runtime temporal interpolation half-width was invalid
  const int temporal_num_points = 2 * temporal_half_width + 1;
  if (temporal_num_points != numerical_window->temporal_interp_num_points) {
    #pragma omp parallel for
    for (long int i = 0; i < chunk_size; i++) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      if (d_connection_bundle != NULL)
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++)
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
    } // END LOOP: for i over rays after inconsistent temporal stencil size
    return;
  } // END IF: runtime temporal stencil size did not match the mapped numerical window

  #pragma omp parallel for
  for (long int i = 0; i < chunk_size; i++) {
    double f_local[9];
    for (int comp = 0; comp < 9; comp++)
      f_local[comp] = d_f_bundle[IDX_F(comp, i)];
    const REAL t = (REAL)f_local[0];
    const REAL x = (REAL)f_local[1];
    const REAL y = (REAL)f_local[2];
    const REAL z = (REAL)f_local[3];
    int ray_failed = 0;
    uint64_t available_slice_indices[temporal_num_points];
    REAL available_slice_times[temporal_num_points];
    const double *available_slice_payloads[temporal_num_points];
    int num_available_slices = 0;
    REAL missing_slice_times[temporal_num_points];
    int num_missing_slices = 0;
    REAL g4dd_available[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL gamma4udd_available[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT];
    REAL full_slice_times[temporal_num_points];
    REAL g4dd_slices[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL gamma4udd_slices[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT];
    REAL g4dd_local[TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL gamma4udd_local[TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT];
    REAL g4dd_missing_local[TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT] = {0};
    REAL gamma4udd_missing_local[TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT] = {0};

    // Step 1: Dispatch directly to the static analytic metrics when this ray
    // is fully outside the numerical time interval.
    if (t < t_metric_0) {
      {metric_worker_0}(commondata, f_local, g4dd_local);
      if (d_connection_bundle != NULL)
        {conn_worker_0}(commondata, f_local, gamma4udd_local);
    } else if (t > t_metric_1) {
      {metric_worker_1}(commondata, f_local, g4dd_local);
      if (d_connection_bundle != NULL)
        {conn_worker_1}(commondata, f_local, gamma4udd_local);
    } else {
      // Step 2: Recover one adaptive numerical stencil from the slot-level
      // time window shared by the whole chunk.
      const int window_status = time_window_manager_numerical_stencil_for_time(
          numerical_window, (double)t, temporal_half_width,
          available_slice_indices, available_slice_times,
          available_slice_payloads, &num_available_slices, missing_slice_times,
          &num_missing_slices);
      if (window_status != TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS ||
          num_available_slices <= 0 ||
          num_available_slices + num_missing_slices != temporal_num_points) {
        ray_failed = 1;
      } else {
        // Step 3: Interpolate only the available mapped numerical slices in
        // space at the photon position.
        const int spatial_status =
            {spatial_name}(
                spatial_context, commondata, params, x, y, z,
                num_available_slices, available_slice_payloads, g4dd_available,
                gamma4udd_available);
        if (spatial_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS)
          ray_failed = 1;
        else {
          int missing_is_upper = 0;
          int missing_is_lower = 0;

          // Step 4: If the centered stencil is only partially mapped, fill
          // the missing edge with one static analytic metric evaluation.
          if (num_missing_slices > 0) {
            missing_is_upper =
                missing_slice_times[0] >
                available_slice_times[num_available_slices - 1];
            missing_is_lower =
                missing_slice_times[num_missing_slices - 1] <
                available_slice_times[0];

            if (missing_is_upper == missing_is_lower) {
              ray_failed = 1;
            } else {
              f_local[0] = (double)missing_slice_times[0];
              if (missing_is_upper) {
                {metric_worker_1}(commondata, f_local, g4dd_missing_local);
                {conn_worker_1}(commondata, f_local, gamma4udd_missing_local);
              } else {
                {metric_worker_0}(commondata, f_local, g4dd_missing_local);
                {conn_worker_0}(commondata, f_local, gamma4udd_missing_local);
              } // END ELSE: lower missing stencil edge uses analytical_metric_0
              f_local[0] = (double)t;
            } // END ELSE: missing stencil edge classification was usable
          } // END IF: one analytic stencil edge must be synthesized

          if (!ray_failed) {
            // Step 5: Reconstruct the full ordered temporal stencil expected
            // by the temporal interpolation helper.
            if (num_missing_slices > 0 && missing_is_lower) {
              for (int s = 0; s < num_missing_slices; s++) {
                full_slice_times[s] = missing_slice_times[s];
                for (int comp = 0;
                     comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                  G4_SLICE(g4dd_slices, s, comp) = g4dd_missing_local[comp];
                } // END LOOP: for comp over missing analytic metric components
                for (int comp = 0;
                     comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
                  GAMMA_SLICE(gamma4udd_slices, s, comp) =
                      gamma4udd_missing_local[comp];
                } // END LOOP: for comp over missing analytic Christoffel components
              } // END LOOP: for s over lower missing stencil nodes
              for (int s = 0; s < num_available_slices; s++) {
                const int full_slot = num_missing_slices + s;
                full_slice_times[full_slot] = available_slice_times[s];
                for (int comp = 0;
                     comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                  G4_SLICE(g4dd_slices, full_slot, comp) =
                      G4_SLICE(g4dd_available, s, comp);
                } // END LOOP: for comp over mapped numerical metric components after lower analytic padding
                for (int comp = 0;
                     comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
                  GAMMA_SLICE(gamma4udd_slices, full_slot, comp) =
                      GAMMA_SLICE(gamma4udd_available, s, comp);
                } // END LOOP: for comp over mapped numerical Christoffel components after lower analytic padding
              } // END LOOP: for s over available numerical stencil nodes after lower analytic padding
            } else {
              for (int s = 0; s < num_available_slices; s++) {
                full_slice_times[s] = available_slice_times[s];
                for (int comp = 0;
                     comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                  G4_SLICE(g4dd_slices, s, comp) =
                      G4_SLICE(g4dd_available, s, comp);
                } // END LOOP: for comp over mapped numerical metric components
                for (int comp = 0;
                     comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
                  GAMMA_SLICE(gamma4udd_slices, s, comp) =
                      GAMMA_SLICE(gamma4udd_available, s, comp);
                } // END LOOP: for comp over mapped numerical Christoffel components
              } // END LOOP: for s over available numerical stencil nodes
              for (int s = 0; s < num_missing_slices; s++) {
                const int full_slot = num_available_slices + s;
                full_slice_times[full_slot] = missing_slice_times[s];
                for (int comp = 0;
                     comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                  G4_SLICE(g4dd_slices, full_slot, comp) = g4dd_missing_local[comp];
                } // END LOOP: for comp over upper missing analytic metric components
                for (int comp = 0;
                     comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
                  GAMMA_SLICE(gamma4udd_slices, full_slot, comp) =
                      gamma4udd_missing_local[comp];
                } // END LOOP: for comp over upper missing analytic Christoffel components
              } // END LOOP: for s over upper missing stencil nodes
            } // END ELSE: upper missing edge or fully numerical centered stencil
          } // END IF: full ordered stencil reconstruction remained valid

          if (!ray_failed) {
            // Step 6: Interpolate the reconstructed per-slice tensor bundles
            // in physical time to the photon coordinate time.
            const int temporal_status = temporal_lagrange_interpolation(
                commondata, full_slice_times, g4dd_slices, gamma4udd_slices, t,
                g4dd_local, gamma4udd_local);
            if (temporal_status != TEMPORAL_LAGRANGE_INTERP_SUCCESS)
              ray_failed = 1;
          } // END IF: reconstructed stencil was ready for temporal interpolation
        } // END ELSE: spatial interpolation succeeded for the mapped numerical stencil subset
      } // END ELSE: adaptive stencil query succeeded for this mixed ray
    } // END ELSE: photon required numerical or mixed interpolation

    if (ray_failed) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      if (d_connection_bundle != NULL)
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++)
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
      continue;
    } // END IF: at least one interpolation stage failed for this ray

    for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
      d_metric_bundle[IDX_METRIC(comp, i)] = (double)g4dd_local[comp];
    if (d_connection_bundle != NULL)
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++)
        d_connection_bundle[IDX_CONN(comp, i)] = (double)gamma4udd_local[comp];
  } // END LOOP: for i over rays in chunk

  #undef IDX_F
  #undef IDX_METRIC
  #undef IDX_CONN
  #undef G4_SLICE
  #undef GAMMA_SLICE
""".replace("{spatial_name}", spatial_name)
        .replace("{metric_worker_0}", metric_worker_0)
        .replace("{metric_worker_1}", metric_worker_1)
        .replace("{conn_worker_0}", conn_worker_0)
        .replace("{conn_worker_1}", conn_worker_1)
    )

    cfc.register_CFunction(
        subdirectory="interpolation",
        prefunc=prefunc,
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
