"""
Register chunk-based numerical-spacetime interpolation.

This module emits the host-side numerical interpolation wrapper used by
geodesic integrators. The generated C function mirrors the analytic
interpolation-kernel bundle contract: it consumes one chunk of photon states,
parallelizes over rays on the CPU, writes the 10-component metric bundle, and
writes the 40-component Christoffel bundle only when requested.

Operationally, this wrapper is the bridge between the data-management layer
and the numerical interpolation helpers. For each photon, it either spatially
interpolates the first selected numerical slice directly or asks the numerical
time-window manager for one adaptive centered temporal stencil. In the mixed
case, the wrapper spatially interpolates only the mapped numerical stencil
subset, fills lower missing stencil nodes from the first selected numerical
slice, fills upper missing stencil nodes by freezing the final selected
numerical slice, and then runs the temporal helper once on the reconstructed
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
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.azimuthal_symmetry_spatial_lagrange_interpolation_metric_derivatives import (
    register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.temporal_lagrange_interpolation_metric_derivatives import (
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
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the CPU numerical-spacetime interpolation wrapper.

    This wrapper owns the per-photon orchestration of the numerical-spacetime
    interpolation pipeline. It does not select or map the active numerical
    window itself; instead, it assumes a caller already mapped a conservative
    slot-based window and then evaluates every ray in one chunk against that
    shared mapped data. Rays at or below `t_metric_0` spatially interpolate
    the first selected numerical slice directly, treating that slice as static.
    Rays at or above the final selected numerical slice time use that final
    numerical slice directly, treating the numerical spacetime as frozen in
    time thereafter. Rays between those bounds use the adaptive
    `time_window_manager_numerical_stencil_for_time()` contract from
    `time_window_manager_numerical`, spatially interpolate only the mapped
    numerical stencil subset, fill lower missing stencil nodes by spatially
    interpolating the first selected numerical slice, fill upper missing stencil
    nodes by reusing the final selected numerical slice, and then perform
    ordinary temporal interpolation on the reconstructed full stencil.

    :param CoordSystem: Coordinate system used by the mapped numerical dataset;
        must be `"Spherical"`, `"SinhSpherical"`, `"Cylindrical"`,
        `"SinhCylindrical"`, or `"SinhCylindricalv2n2"`.
    :param enable_simd: Whether SIMD helper headers are already available.
    :param project_dir: Destination project directory for copied headers.
    :param analytical_metric_0: Retained generator API argument. The current
        first-slice extension does not call this analytic metric.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If `CoordSystem` is not supported or if the retained
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
        "SinhCylindricalv2n2",
    ):
        raise ValueError(
            "numerical_interpolation currently supports only CoordSystem in "
            "('Spherical', 'SinhSpherical', 'Cylindrical', "
            "'SinhCylindrical', 'SinhCylindricalv2n2'); "
            f"found '{CoordSystem}'."
        )
    if not analytical_metric_0:
        raise ValueError("The lower-time analytic metric name must be non-empty.")

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
            "Coordinate time at and below which the first stored numerical "
            "slice is spatially interpolated as a static spacetime."
        ),
    )
    _ = par.register_CodeParameter(
        "REAL",
        __name__,
        "t_numerical_end",
        1.0,
        commondata=True,
        add_to_parfile=True,
        description=(
            "Coordinate time of the final selected numerical slice. The "
            "numerical interpolation wrapper validates this user-facing value "
            "against the selected stride endpoint within two selected time "
            "steps and freezes the numerical spacetime thereafter."
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
    conn_worker_0 = f"connections_{analytical_metric_0}"
    prefunc_lines = []
    for worker_name in [metric_worker_0, conn_worker_0]:
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
        "<stdio.h>",
        "<stdint.h>",
        "<stdlib.h>",
    ]

    desc = r"""Interpolate piecewise static/numerical spacetime tensors for one photon chunk.

The caller supplies an active numerical time window, a spatial interpolation
context, and one chunk of photon states in the same Structure-of-Arrays bundle
layout used by the analytic geodesic interpolation kernel. This CPU wrapper
parallelizes over rays, spatially interpolates the first selected numerical slice
for times at or below `t_metric_0`, freezes the numerical spacetime to the
final selected slice for times at or above that final slice time, and otherwise
reconstructs one mixed temporal stencil. Lower missing stencil nodes use a
fresh spatial interpolation of the first selected numerical slice, whereas upper
missing nodes reuse the final selected numerical slice. Every reconstructed
stencil uses ordinary temporal interpolation.

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
  const REAL t_numerical_end = (REAL)commondata->t_numerical_end;
  const REAL dt_numerical_spacetime_data =
      (REAL)commondata->dt_numerical_spacetime_data;
  // Step 0: Validate the mapped numerical window pointers before reading any
  // of their fields.
  if (numerical_window == NULL || numerical_window->slice_times == NULL ||
      numerical_window->num_time_slices < 1ULL ||
      numerical_window->time_slice_stride == 0ULL) {
    #pragma omp parallel for
    for (long int i = 0; i < chunk_size; i++) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      if (d_connection_bundle != NULL)
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++)
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
    } // END LOOP: for i over rays after invalid mapped numerical window metadata
    return;
  } // END IF: mapped numerical window metadata was unavailable
  const uint64_t first_slice_index = 0ULL;
  const uint64_t final_slice_index =
      time_window_manager_numerical_final_selected_slice(numerical_window);
  const REAL t_final_numerical_slice =
      (REAL)numerical_window->slice_times[final_slice_index];
  const REAL selected_slice_dt =
      dt_numerical_spacetime_data *
      (REAL)numerical_window->time_slice_stride;
  // The user-facing t_numerical_end must identify the final selected numerical
  // slice to within two selected time steps.
  const REAL t_numerical_end_tolerance =
      2.0 * selected_slice_dt;
  // The mapped numerical window and the temporal helper must agree on the
  // centered temporal stencil width before any variable-length arrays are
  // sized from runtime data.
  if (temporal_half_width < 0 ||
      temporal_half_width > TEMPORAL_LAGRANGE_INTERP_MAX_HALF_WIDTH ||
      temporal_half_width != numerical_window->temporal_interp_half_width ||
      !isfinite((double)t_metric_0) ||
      !isfinite((double)t_numerical_end) ||
      !isfinite((double)dt_numerical_spacetime_data) ||
      dt_numerical_spacetime_data <= 0.0 ||
      !isfinite((double)selected_slice_dt) ||
      selected_slice_dt <= 0.0 ||
      !isfinite((double)t_final_numerical_slice) ||
      t_metric_0 >= t_numerical_end) {
    #pragma omp parallel for
    for (long int i = 0; i < chunk_size; i++) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      if (d_connection_bundle != NULL)
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++)
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
    } // END LOOP: for i over rays after invalid temporal bounds or stencil metadata
    return;
  } // END IF: temporal bounds or mapped stencil metadata were invalid
  if (fabs((double)(t_numerical_end - t_final_numerical_slice)) >
      (double)t_numerical_end_tolerance) {
    fprintf(stderr,
            "ERROR: commondata->t_numerical_end=%e differs from the final selected numerical slice time=%e "
            "by more than 2*selected_slice_dt=%e.\n",
            (double)t_numerical_end, (double)t_final_numerical_slice,
            (double)t_numerical_end_tolerance);
    exit(1);
  } // END IF: t_numerical_end was inconsistent with the final selected numerical slice time
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

    // Step 1: Dispatch directly to a static numerical endpoint when this ray
    // is outside the mixed temporal-interpolation region. At and below
    // t_metric_0, interpolate the first stored slice in space only, so its
    // metric time derivatives remain zero. At or above the final selected
    // slice time, likewise use the final slice in space only.
    if (t <= t_metric_0) {
      const double *first_slice_payloads[1];
      first_slice_payloads[0] =
          time_window_manager_numerical_grid_ptr(numerical_window, first_slice_index);
      if (first_slice_payloads[0] == NULL) {
        ray_failed = 1;
      } else {
        const int spatial_status =
            {spatial_name}(
                spatial_context, commondata, params, x, y, z,
                1, first_slice_payloads, g4dd_local, gamma4udd_local);
        if (spatial_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS)
          ray_failed = 1;
      } // END ELSE: first selected numerical slice was available for direct spatial interpolation
    } else if (t >= t_final_numerical_slice) {
      const double *final_slice_payloads[1];
      final_slice_payloads[0] =
          time_window_manager_numerical_grid_ptr(numerical_window, final_slice_index);
      if (final_slice_payloads[0] == NULL) {
        ray_failed = 1;
      } else {
        const int spatial_status =
            {spatial_name}(
                spatial_context, commondata, params, x, y, z,
                1, final_slice_payloads, g4dd_local, gamma4udd_local);
        if (spatial_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS)
          ray_failed = 1;
      } // END ELSE: final selected numerical slice was available for direct spatial interpolation
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

          // Step 4: Centered lower-boundary stencils pad their missing negative
          // time nodes from the static first numerical slice. Freeze the final
          // selected numerical slice for ordinary upper missing nodes.
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
              if (missing_is_upper) {
                const int last_available_slot = num_available_slices - 1;
                if (available_slice_indices[last_available_slot] != final_slice_index) {
                  ray_failed = 1;
                } else {
                  for (int comp = 0;
                       comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                    g4dd_missing_local[comp] =
                        G4_SLICE(g4dd_available, last_available_slot, comp);
                  } // END LOOP: for comp over upper frozen metric components
                  for (int comp = 0;
                       comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
                    gamma4udd_missing_local[comp] =
                        GAMMA_SLICE(gamma4udd_available, last_available_slot, comp);
                  } // END LOOP: for comp over upper frozen Christoffel components
                } // END ELSE: upper missing stencil edge reached the final selected numerical slice
              } else {
                const double *first_slice_payloads[1];
                first_slice_payloads[0] =
                    time_window_manager_numerical_grid_ptr(
                        numerical_window, first_slice_index);
                if (first_slice_payloads[0] == NULL) {
                  ray_failed = 1;
                } else {
                  const int first_slice_spatial_status =
                      {spatial_name}(
                          spatial_context, commondata, params, x, y, z,
                          1, first_slice_payloads, g4dd_missing_local,
                          gamma4udd_missing_local);
                  if (first_slice_spatial_status !=
                      AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS)
                    ray_failed = 1;
                } // END ELSE: first selected numerical slice was available for lower stencil padding
              } // END ELSE: lower missing stencil edge uses the first selected numerical slice
            } // END ELSE: missing stencil edge classification was usable
          } // END IF: one stencil edge must be synthesized

          if (!ray_failed) {
            // Step 5: Reconstruct the full ordered temporal stencil expected
            // by the temporal interpolation helper.
            if (num_missing_slices > 0 && missing_is_lower) {
              for (int s = 0; s < num_missing_slices; s++) {
                full_slice_times[s] = missing_slice_times[s];
                for (int comp = 0;
                     comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                  G4_SLICE(g4dd_slices, s, comp) = g4dd_missing_local[comp];
                } // END LOOP: for comp over missing first-slice metric components
                for (int comp = 0;
                     comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
                  GAMMA_SLICE(gamma4udd_slices, s, comp) =
                      gamma4udd_missing_local[comp];
                } // END LOOP: for comp over missing first-slice metric derivatives
              } // END LOOP: for s over lower missing static first-slice nodes
              for (int s = 0; s < num_available_slices; s++) {
                const int full_slot = num_missing_slices + s;
                full_slice_times[full_slot] = available_slice_times[s];
                for (int comp = 0;
                     comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                  G4_SLICE(g4dd_slices, full_slot, comp) =
                      G4_SLICE(g4dd_available, s, comp);
                } // END LOOP: for comp over mapped numerical metric components after lower first-slice padding
                for (int comp = 0;
                     comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
                  GAMMA_SLICE(gamma4udd_slices, full_slot, comp) =
                      GAMMA_SLICE(gamma4udd_available, s, comp);
                } // END LOOP: for comp over mapped numerical metric derivatives after lower first-slice padding
              } // END LOOP: for s over available numerical stencil nodes after lower first-slice padding
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
                } // END LOOP: for comp over upper missing frozen metric components
                for (int comp = 0;
                     comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
                  GAMMA_SLICE(gamma4udd_slices, full_slot, comp) =
                      gamma4udd_missing_local[comp];
                } // END LOOP: for comp over upper missing frozen Christoffel components
              } // END LOOP: for s over upper missing stencil nodes
            } // END ELSE: upper missing edge or fully numerical centered stencil
          } // END IF: full ordered stencil reconstruction remained valid

          if (!ray_failed) {
            // Step 6: Interpolate the reconstructed stencil in physical time.
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
        .replace("{conn_worker_0}", conn_worker_0)
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
