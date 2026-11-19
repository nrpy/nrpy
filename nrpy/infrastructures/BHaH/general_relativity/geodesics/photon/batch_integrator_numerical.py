# infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_numerical.py
r"""
Emit the CPU/OpenMP numerical-spacetime photon batch integrator.

This module registers the C orchestrator for batched photon geodesics in a
numerical spacetime. The generated integrator consumes a validated combined
numerical raytracing ``.bin`` path stored on ``commondata``, initializes a
``NumericalTimeWindowManager``, maps numerical time windows per time slot, and
calls ``numerical_interpolation()`` for the metric and the 40 components used
by the photon equations of motion. The observer metric is interpolated once at the common
initial event; the numerical initializer constructs one validated observer
tetrad per tile-batch call and passes that same tetrad to every ray in the
batch.

Fast-light generation keeps the numerical window at ``t_start`` while photon
coordinate times and integration parameters advance through scheduling slots.
The numerical interpolation registration must use the same sampling option.

The numerical time-window logic assumes the generator also registers
``rkf45_finalize_and_control_kernel(enable_numerical_time_window_step_cap=True)``
so accepted RKF45 steps remain inside the mapped numerical window.

Each RK attempt selects one spatial mask across its locked temporal nodes.
All stages and direct-EOM endpoint metric evaluations reuse that mask.
Delayed norm and angular-momentum checks rebuild the mask from the accepted
stencil selection or the saved plane-crossing candidate.
Both EOM modes then pass the common ``ln|alpha p^0|`` scalar to event
detection; normalized mode reads it directly from ``u``.

The terminal display counts retryable rejections in the latest processed
photon chunk and separates RKF45 error control from spatial and temporal
interpolation failures. All rejection percentages use that chunk's attempt
count. Terminal failures and accepted-step time caps do not enter these counts.

The script keeps the broad photon orchestration, RKF45 stepping, and
``TimeSlotManager`` structure used by the analytic photon batch integrator, but
all geometry evaluations go through the numerical ``.bin`` data path. When the
generator declares axial symmetry about the Cartesian ``z`` axis, optional
initial, nonterminal plane, and terminal ``L_z`` evaluations measure
angular-momentum conservation.

An optional ``apparent_horizon_bin_path`` enables stopping inside saved apparent
horizons. The puncture output records the latest horizon ID, containment status,
and signed radial gap for each photon. The separate directory diagnostic
classifies final evolution-measure stops when binary horizon stopping is disabled.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.find_plane_crossings_brent import (
    register_plane_crossing_candidates,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.handle_non_terminal_plane_intersection import (
    register_non_terminal_plane_crossing_result,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.handle_terminal_plane_intersection import (
    register_terminal_plane_crossing_result,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.puncture_diagnostics import (
    puncture_diagnostics,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.time_slot_manager_helpers import (
    time_slot_manager_helpers,
)


def batch_integrator_numerical(
    spacetime_name: str,
    dataset_coord_system: str,
    interpolation_method: str = "g4DD",
    normalized_eom: bool = False,
    axisymmetric_about_z: bool = False,
    fast_light: bool = False,
    enable_metric_gradient_tolerances: bool = False,
) -> None:
    r"""
    Construct the CPU numerical-spacetime photon batch integrator.

    :param spacetime_name: Spacetime identifier used to validate registration context.
    :param dataset_coord_system: Coordinate system used by the numerical dataset.
    :param interpolation_method: Numerical metric or connection interpolation method.
    :param normalized_eom: Whether to evolve normalized coordinate-time photon equations.
    :param axisymmetric_about_z: Whether the numerical spacetime has a rotational
        Killing vector about the Cartesian ``z`` axis, permitting ``L_z`` checks.
    :param fast_light: Whether to map geometry at ``t_start`` for integration
        and diagnostics; requires matching numerical interpolation registration.
    :param enable_metric_gradient_tolerances: Whether to accumulate the maximum
        absolute metric derivative over each trial's stages and pass it to
        RKF45 finalization. Requires ``g4DD`` or ``g4DD_d0`` interpolation and
        matching finalization registration. A minimum scale of one or zero
        sensing length skips the derivative scan and its storage allocation.
    :raises ValueError: If the parallelization mode, interpolation method, or
        dataset coordinate system is unsupported.

    Doctests:
    >>> import os
    >>> import tempfile
    >>> from unittest.mock import patch
    >>> import nrpy.c_function as cfc
    >>> with tempfile.TemporaryDirectory() as cache_dir, patch.dict(os.environ, {"XDG_CACHE_HOME": cache_dir}):
    ...     cfc.CFunction_dict.clear()
    ...     batch_integrator_numerical(
    ...         "Schwarzschild", "SinhCylindricalv2n2", interpolation_method="g4DD"
    ...     )
    ...     generated = cfc.CFunction_dict["batch_integrator_numerical"].full_function
    >>> "@param[in,out] commondata" in generated
    True
    >>> "time_window_manager_numerical_mmap_for_slot" in generated
    True
    >>> "observer_metric[10]" in generated
    True
    >>> "normal_observer_log_energy(" in generated
    True
    >>> "d_spatial_stencil_center_i0[current]" in generated
    True
    >>> "d_log_energy_bundle[current]" in generated
    True
    >>> generated_modes = [generated]
    >>> with tempfile.TemporaryDirectory() as cache_dir, patch.dict(os.environ, {"XDG_CACHE_HOME": cache_dir}):
    ...     cfc.CFunction_dict.clear()
    ...     batch_integrator_numerical(
    ...         "Schwarzschild", "SinhCylindricalv2n2", interpolation_method="g4DD",
    ...         normalized_eom=True
    ...     )
    ...     generated = cfc.CFunction_dict["batch_integrator_numerical"].full_function
    >>> "fabs(exp(2.0 * d_f_bundle[0][4 * BUNDLE_CAPACITY + norm_i]) * (d_norm_bundle[norm_i].C - 1.0))" in generated
    True
    >>> "const double current_norm_err = fabs(d_norm_bundle[norm_i].C - 1.0);" in generated
    True
    >>> "const double sidecar_norm_err = fabs(exp(2.0 * d_f_bundle[0][4 * BUNDLE_CAPACITY + norm_i]) * (d_norm_bundle[norm_i].C - 1.0));" in generated
    True
    >>> "normalization_abs_by_ray[master_idx] = sidecar_norm_err;" in generated
    True
    >>> "non_terminal_plane_norm_recorded_by_ray" in generated
    True
    >>> "normalization_abs_non_terminal_by_ray" in generated
    True
    >>> "all_photons_host.non_terminal_norm_f" in generated
    False
    >>> "normal_observer_log_energy(" in generated
    False
    >>> "d_f_bundle[current][4 * BUNDLE_CAPACITY + log_energy_i]" in generated
    True
    >>> generated_modes.append(generated)
    >>> for code in generated_modes:
    ...     source = "".join(code.split())
    ...     for buffer in ("current", "next"):
    ...         assert f"min_steps_bridge[{buffer}][bridge_i]=all_photons_host.consecutive_min_steps[m_idx];" in source
    ...         assert f"d_retries[{buffer}],d_min_steps[{buffer}],active_chunks[{buffer}]" in source
    ...     assert "all_photons_host.consecutive_min_steps[m_idx]=min_steps_bridge[current][fin_i];" in source
    ...     limit_check = source.index("constintmin_step_limit=commondata->rkf45_max_consecutive_min_steps;")
    ...     assert source.rfind("event_detection_manager_kernel(", 0, limit_check) != -1
    ...     assert limit_check < source.index("slot_add_photon(", limit_check)
    ...     assert "status_bridge[current][fin_i]=FAILURE_RKF45_MIN_STEP_LIMIT;all_photons_host.status[m_idx]=FAILURE_RKF45_MIN_STEP_LIMIT;" in source
    >>> expected_counting = '''
    ... for (int sum_i = 0; sum_i < active_chunks[current]; ++sum_i) {
    ...     if (status_bridge[current][sum_i] == REJECTED) {
    ...         batch_rejections++;
    ...         if (trial_status[current][sum_i] == FAILURE_SPATIAL_INTERPOLATION) {
    ...             spatial_rejections++;
    ...         } else if (trial_status[current][sum_i] == FAILURE_TEMPORAL_INTERPOLATION) {
    ...             temporal_rejections++;
    ...         } else {
    ...             rkf45_rejections++;
    ...         }
    ...     }
    ... }
    ... '''
    >>> for generated in generated_modes:
    ...     source = "".join(line.split("//", 1)[0].strip() for line in generated.splitlines())
    ...     source = "".join(source.split())
    ...     for buffer in ("current", "next"):
    ...         snapshot = f"memcpy(trial_status[{buffer}],d_status[{buffer}],sizeof(termination_type_t)*active_chunks[{buffer}]);"
    ...         finalize = f"rkf45_finalize_and_control(commondata,d_f_bundle[{buffer}]"
    ...         assert source.count(snapshot) == 1
    ...         snapshot_index = source.index(snapshot)
    ...         finalize_index = source.index(finalize)
    ...         assert source.rfind("rkf45_stage_update(", 0, finalize_index) < snapshot_index < finalize_index
    ...     loop_start = source.index("for(intsum_i=0;sum_i<active_chunks[current];++sum_i)")
    ...     loop_end = source.index("constdoublereject_percent", loop_start)
    ...     assert source[loop_start:loop_end] == "".join(expected_counting.split())
    ...     for count in ("batch", "rkf45", "spatial", "temporal"):
    ...         percentage = f"(active_chunks[current]>0)?(100.0*(double){count}_rejections/(double)active_chunks[current]):0.0;"
    ...         assert percentage in source
    ...     assert source.index("BHAH_MALLOC(trial_status[s],sizeof(termination_type_t)*BUNDLE_CAPACITY);") < loop_start
    ...     assert source.index("BHAH_FREE(trial_status[s]);") > loop_end
    ...     assert source.count(r'printf("\n\n\n\n\n\n\n\n");') == 1
    ...     assert source.count(r'printf("\033[8A");') == 1
    >>> mapping_call = "time_window_manager_numerical_mmap_for_slot("
    >>> expected_managers = ["tsm", "tsm", "non_terminal_norm_tsm", "diagnostic_tsm"]
    >>> for generated in generated_modes:
    ...     source = "".join(generated.split())
    ...     mappings = [tuple(call.split(")", 1)[0].split(",")[1:])
    ...                 for call in source.split(mapping_call)[1:]]
    ...     assert sorted(mappings) == sorted(zip(["&" + manager for manager in expected_managers], [
    ...         "initial_slot_idx", "slot_idx", "norm_slot_idx", "diagnostic_slot_idx"
    ...     ]))
    >>> for method in ("g4DD", "g4DD_d0", "GammaUDD"):
    ...     for normalized in (False, True):
    ...         with tempfile.TemporaryDirectory() as cache_dir, patch.dict(os.environ, {"XDG_CACHE_HOME": cache_dir}):
    ...             cfc.CFunction_dict.clear()
    ...             batch_integrator_numerical(
    ...                 "Numerical", "SinhCylindricalv2n2", interpolation_method=method,
    ...                 normalized_eom=normalized, fast_light=True
    ...             )
    ...             generated = cfc.CFunction_dict["batch_integrator_numerical"].full_function
    ...         source = "".join(generated.split())
    ...         mappings = [tuple(call.split(")", 1)[0].split(",")[1:])
    ...                     for call in source.split(mapping_call)[1:]]
    ...         assert sorted(mappings) == sorted([
    ...             ("&" + manager, "initial_slot_idx") for manager in expected_managers
    ...         ])
    ...         coordinate_time = (
    ...             "all_photons_host.integration_param[m_idx]" if normalized
    ...             else "all_photons_host.f[0*num_rays+m_idx]"
    ...         )
    ...         assert f"slot_get_index(&tsm,{coordinate_time})" in source
    >>> with tempfile.TemporaryDirectory() as cache_dir, patch.dict(os.environ, {"XDG_CACHE_HOME": cache_dir}):
    ...     batch_integrator_numerical(
    ...         "Numerical", "SinhCylindricalv2n2", interpolation_method="GammaUDD",
    ...         enable_metric_gradient_tolerances=True
    ...     )
    Traceback (most recent call last):
        ...
    ValueError: Metric-gradient tolerances require g4DD or g4DD_d0.
    """
    if "time_slot_manager" not in par.glb_extras_dict.get("BHaH_defines", {}):
        time_slot_manager_helpers()
    register_plane_crossing_candidates()
    register_non_terminal_plane_crossing_result()
    register_terminal_plane_crossing_result()

    if not isinstance(axisymmetric_about_z, bool):
        raise ValueError(
            "axisymmetric_about_z must be a bool, got "
            f"{type(axisymmetric_about_z).__name__}."
        )
    if not spacetime_name:
        raise ValueError("spacetime_name must contain a valid identifier.")

    parallelization = par.parval_from_str("parallelization")
    if parallelization != "openmp":
        raise ValueError(
            "batch_integrator_numerical currently supports only "
            "parallelization='openmp'."
        )
    if interpolation_method not in ("g4DD", "g4DD_d0", "GammaUDD"):
        raise ValueError(
            "interpolation_method must be one of ('g4DD', 'g4DD_d0', 'GammaUDD'); "
            f"found '{interpolation_method}'."
        )
    if enable_metric_gradient_tolerances and interpolation_method == "GammaUDD":
        raise ValueError("Metric-gradient tolerances require g4DD or g4DD_d0.")
    if dataset_coord_system != "SinhCylindricalv2n2":
        raise ValueError(
            "batch_integrator_numerical currently supports only "
            "dataset_coord_system='SinhCylindricalv2n2'; "
            f"found '{dataset_coord_system}'."
        )
    phi_dim = 1
    geometry_slot_idx = "initial_slot_idx" if fast_light else "slot_idx"
    non_terminal_geometry_slot_idx = (
        "initial_slot_idx" if fast_light else "norm_slot_idx"
    )
    terminal_geometry_slot_idx = (
        "initial_slot_idx" if fast_light else "diagnostic_slot_idx"
    )

    # Supplied-name mapping: the checklist aliases numerical_time_window_manager_set_inert,
    # numerical_time_window_manager_init, numerical_time_window_manager_mmap_for_slot,
    # and numerical_time_window_manager_free correspond to the supplied helper names
    # time_window_manager_numerical_set_inert/init/mmap_for_slot/free.
    # Core physics and numerical simulation parameters for the global spacetime struct.
    par.register_CodeParameters(
        "REAL",
        __name__,
        ["r_escape"],
        [150.0],
        commondata=True,
        add_to_parfile=True,
    )
    par.register_CodeParameter(
        "char[4096]",
        __name__,
        "numerical_spacetime_bin_path",
        "",
        commondata=True,
        add_to_parfile=True,
        description=(
            "Path to the validated combined numerical raytracing .bin file used by "
            "the numerical photon batch integrator."
        ),
    )
    par.register_CodeParameters(
        "bool",
        __name__,
        ["perform_normalization_check"],
        [False],
        commondata=True,
        add_to_parfile=True,
    )
    par.register_CodeParameter(
        "bool",
        __name__,
        "perform_synthetic_slice_check",
        False,
        commondata=True,
        add_to_parfile=True,
        description=(
            "Record first RKF45 interpolation request times that use synthetic "
            "temporal-stencil nodes below t=0 or above t_numerical_end, then "
            "write one per-photon tile sidecar."
        ),
    )
    if axisymmetric_about_z:
        par.register_CodeParameter(
            "bool",
            __name__,
            "perform_Lz_check",
            False,
            commondata=True,
            add_to_parfile=True,
            description=(
                "Evaluate initial and final axial angular momentum L_z and report "
                "the largest absolute change. Valid only for a numerical spacetime "
                "with rotational symmetry about the Cartesian z axis."
            ),
        )

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<ctype.h>",
        "<limits.h>",
        "<stdio.h>",
        "<stdint.h>",
        "<stdlib.h>",
        "<string.h>",
        "<time.h>",
    ]

    angular_momentum_desc = (
        r"""
    @param[in] Lz_bin_path Optional output filename containing photon index,
                           initial $L_z$, and final $L_z$ for every ray."""
        if axisymmetric_about_z
        else ""
    )
    desc = rf"""CPU numerical-spacetime photon batch integrator.

    This function bins active rays by coordinate time using TimeSlotManager,
    maps combined numerical-spacetime .bin time windows through
    NumericalTimeWindowManager, interpolates the common observer metric once,
    initializes all rays from one metric-orthonormal observer tetrad,
    interpolates the metric and the components needed by the photon RHS,
    advances photons with the existing RKF45 kernels, and writes final
    blueprint results.

    @param[in,out] commondata Struct containing global numerical-spacetime,
                              observer, and integration parameters. The numerical
                              time-window manager fills its authoritative first
                              stored slice time during initialization.
    @param num_rays Number of rays in this tile batch. The caller derives this
                     from the width-side scan density and the selected
                     height-side scan density; this integrator does
                     not calculate image placement or tile metadata.
    @param scan_density_height Caller-selected rows in this tile batch.
    @param[out] results_buffer Host array storing the final physical
                               intersections.
    @param[in] norm_abs_bin_path Optional output filename for the raw float64
                                 absolute direct-scale-equivalent null-constraint
                                 error sidecar.
    @param[in] norm_abs_non_terminal_bin_path Optional output filename for the
                                              sparse nonterminal plane
                                              normalization sidecar.
    @param[in] synthetic_slice_usage_bin_path Optional output filename for the
                                               per-photon synthetic temporal-
                                               stencil usage sidecar.
    @param[in] non_terminal_crossings_path Per-tile nonterminal crossing file.
    @param[in] terminal_crossings_path Per-tile terminal crossing file.{angular_momentum_desc}"""
    desc += "\n    @param[in] puncture_bin_path Optional per-tile puncture-ID, horizon-status, and radial-gap file."

    cfunc_type = "void"

    name = "batch_integrator_numerical"

    angular_momentum_param = (
        ", const char *restrict Lz_bin_path" if axisymmetric_about_z else ""
    )
    params = (
        "commondata_struct *restrict commondata, "
        "long int num_rays, const int scan_density_height, "
        "blueprint_data_t *restrict results_buffer, "
        "const char *restrict norm_abs_bin_path, "
        "const char *restrict norm_abs_non_terminal_bin_path, "
        "const char *restrict synthetic_slice_usage_bin_path, "
        "const char *restrict non_terminal_crossings_path, "
        "const char *restrict terminal_crossings_path"
        ", const char *restrict puncture_bin_path"
        f"{angular_momentum_param}"
    )

    plane_crossing_writer = r"""
static void write_plane_crossings(
    const char *restrict path,
    const long int num_rays,
    const non_terminal_plane_crossing_result_t *restrict non_terminal_results,
    const terminal_plane_crossing_result_t *restrict terminal_results,
    const bool terminal_plane)
{
    FILE *restrict output = fopen(path, "wb");
    if (output == NULL) {
        fprintf(stderr, "ERROR: Could not open plane crossing file '%s'.\n", path);
        exit(EXIT_FAILURE);
    }
    for (long int photon_index = 0; photon_index < num_rays; ++photon_index) {
        const non_terminal_plane_crossing_result_t *non_terminal =
            &non_terminal_results[photon_index];
        const terminal_plane_crossing_result_t *terminal =
            &terminal_results[photon_index];
        if (!(terminal_plane ? terminal->found : non_terminal->found)) continue;
        plane_crossing_record_t record = {0};
        record.photon_index = (uint64_t)photon_index;
        record.brent_iterations = terminal_plane
            ? terminal->brent_iterations : non_terminal->brent_iterations;
        record.integration_param = terminal_plane
            ? terminal->integration_param : non_terminal->integration_param;
        record.y_local = terminal_plane
            ? terminal->y_local : non_terminal->y_local;
        record.z_local = terminal_plane
            ? terminal->z_local : non_terminal->z_local;
        for (int component = 0; component < 9; ++component) {
            record.state[component] = terminal_plane
                ? terminal->state[component] : non_terminal->state[component];
        }
        if (fwrite(&record, sizeof(record), 1, output) != 1) {
            fprintf(stderr, "ERROR: Could not write plane crossing to '%s'.\n", path);
            fclose(output);
            exit(EXIT_FAILURE);
        }
    }
    if (fclose(output) != 0) {
        fprintf(stderr, "ERROR: Could not close plane crossing file '%s'.\n", path);
        exit(EXIT_FAILURE);
    }
}
"""

    include_CodeParameters_h = True

    malloc_pinned = "BHAH_MALLOC"
    malloc_device = "BHAH_MALLOC"
    free_pinned = "BHAH_FREE"
    free_device = "BHAH_FREE"
    pin_comment = "Host memory allocation: CPU RAM mapped for"
    bridge_alloc_comment = (
        "Allocate photon state and integration arrays for this chunk."
    )
    scratch_alloc_comment = (
        "Allocate 1D Host scratchpad arrays for temporal data staging."
    )
    calc_blueprint = (
        "calculate_and_fill_blueprint_data_universal("
        "&all_photons_host, "
        "num_rays, "
        "results_buffer, "
        "NULL, "
        "NULL, "
        "NULL, "
        "NULL, "
        "NULL, "
        "0);"
    )
    calc_blueprint_with_norm_abs = (
        "calculate_and_fill_blueprint_data_universal("
        "&all_photons_host, "
        "num_rays, "
        "results_buffer, "
        "normalization_abs_by_ray, "
        "norm_abs_bin_path, "
        "normalization_abs_non_terminal_by_ray, "
        "non_terminal_plane_norm_recorded_by_ray, "
        "norm_abs_non_terminal_bin_path, "
        "0);"
    )
    set_initial_conditions_call = (
        " set_initial_conditions_kernel"
        "(commondata, num_rays, scan_density_height, &all_photons_host, observer_metric, observer_tetrad);"
    )

    normalized_momentum_conversion = (
        "photon_momentum_to_normalized_kernel("
        "d_f_bundle[0], d_metric_bundle[0], chunk_size);"
        if normalized_eom
        else ""
    )
    initial_integration_param = "commondata->t_start" if normalized_eom else "0.0"
    initial_coordinate_time = (
        "all_photons_host.integration_param[sync_i]"
        if normalized_eom
        else "all_photons_host.f[0 * num_rays + sync_i]"
    )
    updated_coordinate_time = (
        "all_photons_host.integration_param[m_idx]"
        if normalized_eom
        else "all_photons_host.f[0 * num_rays + m_idx]"
    )
    terminal_mask_sampling_time = "accepted_mask_selection[norm_ray].sampling_time"
    interpolation_synthetic_usage_args_current = (
        "commondata->perform_synthetic_slice_check ? chunk_buffer[current] : NULL, "
        "commondata->perform_synthetic_slice_check ? synthetic_slice_usage_by_photon : NULL, "
    )
    interpolation_synthetic_usage_args_next = (
        "commondata->perform_synthetic_slice_check ? chunk_buffer[next] : NULL, "
        "commondata->perform_synthetic_slice_check ? synthetic_slice_usage_by_photon : NULL, "
    )
    interpolation_integration_param_args = (
        "d_spatial_stencil_center_i0[current], "
        "d_spatial_stencil_center_i2[current], "
        "d_temporal_stencil_center[current], "
        "d_spatial_mask[current], d_has_excluded_nodes[current], "
        + (
            "d_integration_param_bundle[current], d_h[current], stage,"
            if normalized_eom
            else ""
        )
        + interpolation_synthetic_usage_args_current
    )
    interpolation_integration_param_args_next = (
        "d_spatial_stencil_center_i0[next], "
        "d_spatial_stencil_center_i2[next], "
        "d_temporal_stencil_center[next], "
        "d_spatial_mask[next], d_has_excluded_nodes[next], "
        + (
            "d_integration_param_bundle[next], d_h[next], stage,"
            if normalized_eom
            else ""
        )
        + interpolation_synthetic_usage_args_next
    )
    interpolation_initial_integration_param_args = (
        "NULL, NULL, NULL, NULL, NULL, "
        + ("d_integration_param_bundle[0], d_h[0], 1," if normalized_eom else "")
        + "NULL, NULL, "
    )
    interpolation_terminal_integration_param_args = (
        "d_spatial_stencil_center_i0[0], d_spatial_stencil_center_i2[0], "
        "d_temporal_stencil_center[0], d_spatial_mask[0], "
        "d_has_excluded_nodes[0], "
        + ("d_integration_param_bundle[0], d_h[0], 1," if normalized_eom else "")
        + "NULL, NULL, "
    )
    initial_integration_param_setup = (
        """
        for (long int init_i = 0; init_i < chunk_size; ++init_i) {
            d_integration_param_bundle[0][init_i] = commondata->t_start;
            d_h[0][init_i] = commondata->initial_h;
        } // END LOOP: for init_i over initial numerical
        """
        if normalized_eom
        else ""
    )
    observer_state_coordinate_time = "0.0" if normalized_eom else "commondata->t_start"
    observer_state_integration_setup = (
        """
    // The normalized interpolation wrapper reads coordinate time from this
    // integration-parameter array while f[0] stores the affine parameter.
    d_integration_param_bundle[0][0] = commondata->t_start;
    d_h[0][0] = commondata->initial_h;
"""
        if normalized_eom
        else ""
    )
    terminal_integration_param_setup = (
        """
                    for (long int norm_i = 0; norm_i < chunk_size; ++norm_i) {
                        const long int master_idx = chunk_buffer[0][norm_i];
                        integration_param_bridge[0][norm_i] =
                            terminal_crossing_results[master_idx].found
                                ? terminal_crossing_results[master_idx].integration_param
                                : all_photons_host.integration_param[master_idx];
                    } // END LOOP: for norm_i over terminal integration
                    memcpy(
                        d_integration_param_bundle[0],
                        integration_param_bridge[0],
                        sizeof(double) * chunk_size);
        """
        if normalized_eom
        else ""
    )
    rhs_integration_param_args = (
        "d_integration_param_bundle[current], d_h[current]," if normalized_eom else ""
    )
    rhs_integration_param_args_next = (
        "d_integration_param_bundle[next], d_h[next]," if normalized_eom else ""
    )
    gradient_validation = ""
    gradient_declarations = ""
    gradient_allocations = ""
    gradient_cleanup = ""
    gradient_update_current = ""
    gradient_argument_current = ""
    if enable_metric_gradient_tolerances:
        gradient_validation = r"""
    if (!isfinite(commondata->L_sense) || commondata->L_sense < 0.0 ||
        !isfinite(commondata->tol_scale_min) ||
        commondata->tol_scale_min <= 0.0 || commondata->tol_scale_min > 1.0) {
        fprintf(stderr, "ERROR: metric-gradient tolerances require finite L_sense >= 0 and 0 < tol_scale_min <= 1.\n");
        exit(EXIT_FAILURE);
    } // END IF: invalid runtime tolerance parameters
    const bool metric_gradient_tolerances_active =
        commondata->tol_scale_min < 1.0 && commondata->L_sense > 0.0;
"""
        gradient_declarations = "double *d_trial_gradient_max[2] = {NULL, NULL};"
        gradient_allocations = r"""
        if (metric_gradient_tolerances_active) {
            d_trial_gradient_max[s] = (double *)malloc(sizeof(double) * BUNDLE_CAPACITY);
            if (d_trial_gradient_max[s] == NULL) {
                fprintf(stderr, "ERROR: could not allocate photon metric-gradient maxima.\n");
                exit(EXIT_FAILURE);
            } // END IF: metric-gradient array allocation failed
        } // END IF: metric-gradient tolerance scaling enabled
"""
        gradient_cleanup = "free(d_trial_gradient_max[s]);"
        gradient_argument_current = "d_trial_gradient_max[current], "
        gradient_update_current = r"""
                if (metric_gradient_tolerances_active) {
                    #pragma omp parallel for
                    for (long int ray = 0; ray < active_chunks[current]; ++ray) {
                        if (d_status[current][ray] != ACTIVE) continue;
                        double stage_gradient_max = 0.0;
                        for (int comp = 0; comp < 40; ++comp) {
                            stage_gradient_max = fmax(stage_gradient_max,
                                fabs(d_rhs_geometry_bundle[current][comp * BUNDLE_CAPACITY + ray]));
                        } // END LOOP: maximum absolute Cartesian metric derivative
                        d_trial_gradient_max[current][ray] = stage == 1
                            ? stage_gradient_max
                            : fmax(d_trial_gradient_max[current][ray], stage_gradient_max);
                    } // END LOOP: accumulate per-photon derivative maxima
                } // END IF: metric-gradient tolerance scaling enabled
"""
    gradient_argument_next = gradient_argument_current.replace("[current]", "[next]")
    gradient_update_next = gradient_update_current.replace("[current]", "[next]")
    spatial_center_declarations = """
    // Trial-locked spatial stencil centers, one pair per double-buffered ray.
    int *d_spatial_stencil_center_i0[2];
    int *d_spatial_stencil_center_i2[2];
    uint64_t *d_temporal_stencil_center[2];
    unsigned char *d_spatial_mask[2];
    bool *d_has_excluded_nodes[2];
    batch_spatial_mask_selection_struct *d_mask_selection[2];
    batch_spatial_mask_selection_struct *accepted_mask_selection = NULL;
"""
    spatial_center_allocations = f"""
        {malloc_device}(d_spatial_stencil_center_i0[s], sizeof(int) * BUNDLE_CAPACITY); // Allocate native dimension-0 trial centers.
        {malloc_device}(d_spatial_stencil_center_i2[s], sizeof(int) * BUNDLE_CAPACITY); // Allocate native dimension-2 trial centers.
        {malloc_device}(d_temporal_stencil_center[s], sizeof(uint64_t) * BUNDLE_CAPACITY); // Allocate retained temporal centers.
        {malloc_device}(d_mask_selection[s], sizeof(batch_spatial_mask_selection_struct) * BUNDLE_CAPACITY);
        {malloc_device}(d_has_excluded_nodes[s], sizeof(bool) * BUNDLE_CAPACITY);
        if (d_mask_selection[s] == NULL || d_has_excluded_nodes[s] == NULL) {{
            fprintf(stderr, "ERROR: could not allocate photon spatial-mask selections.\\n");
            exit(EXIT_FAILURE);
        }} // END IF: mask selection allocation failed
"""
    spatial_center_frees = f"""
            {free_device}(d_spatial_stencil_center_i0[s]); // Purges native dimension-0 trial centers.
            {free_device}(d_spatial_stencil_center_i2[s]); // Purges native dimension-2 trial centers.
            {free_device}(d_temporal_stencil_center[s]); // Release retained temporal centers.
            {free_device}(d_mask_selection[s]);
            {free_device}(d_has_excluded_nodes[s]);
            {free_device}(d_spatial_mask[s]);
"""
    trial_spatial_center_setup_current = r"""
            for (long int stencil_i = 0; stencil_i < active_chunks[current]; ++stencil_i) {
                const REAL trial_cartesian[3] = {
                    (REAL)d_f_start_bundle[current][1 * BUNDLE_CAPACITY + stencil_i],
                    (REAL)d_f_start_bundle[current][2 * BUNDLE_CAPACITY + stencil_i],
                    (REAL)d_f_start_bundle[current][3 * BUNDLE_CAPACITY + stencil_i]};
                REAL trial_native[3];
                int automatic_center_idx[3];
                int selected_center_idx[3];
                d_has_excluded_nodes[current][stencil_i] = false;
                if (time_window_manager_numerical_resolve_spatial_target_and_stencil(
                        &numerical_params, trial_cartesian, NULL, trial_native,
                        automatic_center_idx, selected_center_idx) !=
                    TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
                    // Coordinate inversion belongs to spatial interpolation.
                    // Mark only this photon; the wrapper will emit NAN scratch
                    // outputs and the finalizer will preserve its accepted state.
                    d_status[current][stencil_i] =
                        FAILURE_SPATIAL_INTERPOLATION;
                    d_spatial_stencil_center_i0[current][stencil_i] = 0;
                    d_spatial_stencil_center_i2[current][stencil_i] = 0;
                } else {
                    d_spatial_stencil_center_i0[current][stencil_i] =
                        selected_center_idx[0];
                    d_spatial_stencil_center_i2[current][stencil_i] =
                        selected_center_idx[2];
                    d_mask_selection[current][stencil_i].phi_plane =
                        trial_native[1] < 0.0 ||
                        trial_native[1] >= 3.14159265358979323846 - 1.0e-12 ? 0 : 1;
                } // END ELSE: current trial center resolved
            } // END LOOP: for stencil_i over trial rays
"""
    trial_spatial_center_setup_next = r"""
                for (long int stencil_i = 0; stencil_i < active_chunks[next]; ++stencil_i) {
                    const REAL trial_cartesian[3] = {
                        (REAL)d_f_start_bundle[next][1 * BUNDLE_CAPACITY + stencil_i],
                        (REAL)d_f_start_bundle[next][2 * BUNDLE_CAPACITY + stencil_i],
                        (REAL)d_f_start_bundle[next][3 * BUNDLE_CAPACITY + stencil_i]};
                    REAL trial_native[3];
                    int automatic_center_idx[3];
                    int selected_center_idx[3];
                    d_has_excluded_nodes[next][stencil_i] = false;
                    if (time_window_manager_numerical_resolve_spatial_target_and_stencil(
                            &numerical_params, trial_cartesian, NULL, trial_native,
                            automatic_center_idx, selected_center_idx) !=
                        TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
                        // Coordinate inversion belongs to spatial interpolation.
                        // Mark only this photon; later stages remain lane-local.
                        d_status[next][stencil_i] =
                            FAILURE_SPATIAL_INTERPOLATION;
                        d_spatial_stencil_center_i0[next][stencil_i] = 0;
                        d_spatial_stencil_center_i2[next][stencil_i] = 0;
                    } else {
                        d_spatial_stencil_center_i0[next][stencil_i] =
                            selected_center_idx[0];
                        d_spatial_stencil_center_i2[next][stencil_i] =
                            selected_center_idx[2];
                        d_mask_selection[next][stencil_i].phi_plane =
                            trial_native[1] < 0.0 ||
                            trial_native[1] >= 3.14159265358979323846 - 1.0e-12 ? 0 : 1;
                    } // END ELSE: next trial center resolved
                } // END LOOP: for stencil_i over trial rays
"""

    temporal_midpoint_current = (
        "commondata->t_start"
        if fast_light
        else (
            "d_integration_param_bundle[current][stencil_i] + "
            "0.5 * d_h[current][stencil_i]"
            if normalized_eom
            else "d_f_start_bundle[current][stencil_i] + "
            "0.5 * d_h[current][stencil_i] * "
            "d_f_start_bundle[current][4 * BUNDLE_CAPACITY + stencil_i]"
        )
    )
    temporal_center_setup_current = r"""
            for (long int stencil_i = 0; stencil_i < active_chunks[current]; ++stencil_i) {
                if (d_status[current][stencil_i] != ACTIVE) continue;
                d_temporal_stencil_center[current][stencil_i] =
                    time_window_manager_numerical_nearest_selected_slice_for_time(
                        &numerical_window, @MIDPOINT@);
                batch_spatial_mask_selection_struct *selection =
                    &d_mask_selection[current][stencil_i];
                selection->spatial_center.i0 = d_spatial_stencil_center_i0[current][stencil_i];
                selection->spatial_center.i2 = d_spatial_stencil_center_i2[current][stencil_i];
                selection->temporal_center = d_temporal_stencil_center[current][stencil_i];
                selection->sampling_time = @SAMPLING_TIME@;
                if (time_window_manager_numerical_build_spatial_mask(
                        commondata, &numerical_params, &numerical_window,
                        &selection->spatial_center, selection->phi_plane,
                        selection->temporal_center, selection->sampling_time,
                        d_spatial_mask[current] + (size_t)stencil_i * spatial_mask_node_count,
                        &d_has_excluded_nodes[current][stencil_i]) !=
                    TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
                    d_status[current][stencil_i] = FAILURE_SPATIAL_INTERPOLATION;
                } // END IF: trial mask unavailable
            } // END LOOP: select trial stencils and masks
""".replace("@MIDPOINT@", temporal_midpoint_current).replace(
        "@SAMPLING_TIME@",
        (
            "commondata->t_start"
            if fast_light
            else (
                "d_integration_param_bundle[current][stencil_i]"
                if normalized_eom
                else "d_f_start_bundle[current][stencil_i]"
            )
        ),
    )
    trial_spatial_center_setup_current += temporal_center_setup_current
    trial_spatial_center_setup_next += temporal_center_setup_current.replace(
        "[current]", "[next]"
    )
    accepted_mask_save_current = r"""
            if (accepted_mask_selection != NULL) {
                for (long int mask_i = 0; mask_i < active_chunks[current]; ++mask_i) {
                    if (d_status[current][mask_i] == ACTIVE)
                        accepted_mask_selection[chunk_buffer[current][mask_i]] =
                            d_mask_selection[current][mask_i];
                } // END LOOP: retain accepted stencil selections
            } // END IF: delayed metric diagnostics enabled
"""
    accepted_mask_save_next = accepted_mask_save_current.replace("[current]", "[next]")
    temporal_cap_disable = (
        "commondata->numerical_spacetime_temporal_max_delta_t = 0.0;"
        if fast_light
        else ""
    )

    normalization_kernel_name = (
        "normalization_constraint_photon_normalized"
        if normalized_eom
        else "normalization_constraint_photon"
    )
    normalization_error_expr = (
        "fabs(d_norm_bundle[norm_i].C - 1.0)"
        if normalized_eom
        else "fabs(d_norm_bundle[norm_i].C)"
    )
    # Keep the terminal report in its existing convention. Only the optional
    # per-ray .bin sidecar uses the direct-scale-equivalent normalized error.
    normalization_sidecar_error_expr = (
        "fabs(exp(2.0 * d_f_bundle[0][4 * BUNDLE_CAPACITY + norm_i]) * "
        "(d_norm_bundle[norm_i].C - 1.0))"
        if normalized_eom
        else "fabs(d_norm_bundle[norm_i].C)"
    )
    initial_constraint_target = "1.0" if normalized_eom else "0.0"
    angular_momentum_metric_argument = "" if normalized_eom else "d_metric_bundle[0], "
    angular_momentum_allocation = (
        r"""
    // Axial-angular-momentum arrays exist only in projects generated for a
    // numerical spacetime with rotational symmetry about the Cartesian z axis.
    axial_angular_momentum_t *d_angular_momentum_bundle = NULL;
    axial_angular_momentum_t *initial_angular_momentum = NULL;
    axial_angular_momentum_t *final_angular_momentum = NULL;
    if (commondata->perform_Lz_check) {
        BHAH_MALLOC(
            d_angular_momentum_bundle,
            sizeof(axial_angular_momentum_t) * BUNDLE_CAPACITY);
        BHAH_MALLOC(
            initial_angular_momentum,
            sizeof(axial_angular_momentum_t) * num_rays);
        BHAH_MALLOC(
            final_angular_momentum,
            sizeof(axial_angular_momentum_t) * num_rays);
        for (long int Lz_init_i = 0; Lz_init_i < num_rays; ++Lz_init_i) {
            initial_angular_momentum[Lz_init_i].Lz = NAN;
            final_angular_momentum[Lz_init_i].Lz = NAN;
        }
    } // END IF: allocate axial-angular-momentum arrays
"""
        if axisymmetric_about_z
        else ""
    )
    initial_angular_momentum_evaluation = (
        f"""
        if (commondata->perform_Lz_check) {{
            axial_angular_momentum_z(
                d_f_bundle[0],
                {angular_momentum_metric_argument}d_angular_momentum_bundle,
                chunk_size,
                0);
            for (long int Lz_i = 0; Lz_i < chunk_size; ++Lz_i) {{
                const long int master_idx = start_idx + Lz_i;
                const double initial_Lz = d_angular_momentum_bundle[Lz_i].Lz;
                if (!isfinite(initial_Lz)) {{
                    fprintf(
                        stderr,
                        "ERROR: Initial L_z is non-finite for ray %ld.\\n",
                        master_idx);
                    time_window_manager_numerical_free(&numerical_window);
                    slot_manager_free(&tsm);
                    exit(1);
                }}
                initial_angular_momentum[master_idx].Lz = initial_Lz;
            }} // END LOOP: record initial photon L_z
        }} // END IF: evaluate initial L_z
"""
        if axisymmetric_about_z
        else ""
    )
    finalize_current = r"""
            // Preserve interpolation failure reasons before retryable trials
            // receive the common REJECTED status during finalization.
            memcpy(trial_status[current], d_status[current], sizeof(termination_type_t) * active_chunks[current]);
            // Finalize step: apply RKF45 error control and update the
            // integration-parameter baseline and step size $h$.
            rkf45_finalize_and_control(commondata, d_f_bundle[current], d_f_start_bundle[current], d_k_bundle[current], d_h[current], d_status[current], d_integration_param_bundle[current], d_retries[current], d_min_steps[current], @GRADIENT_ARGUMENT@active_chunks[current], current);""".replace(
        "@GRADIENT_ARGUMENT@", gradient_argument_current
    )
    finalize_next = r"""
                // Keep this chunk's reasons separate from the current buffer.
                memcpy(trial_status[next], d_status[next], sizeof(termination_type_t) * active_chunks[next]);
                // Finalize step: apply RKF45 error control and update the
                // integration-parameter baseline and step size $h$.
                rkf45_finalize_and_control(commondata, d_f_bundle[next], d_f_start_bundle[next], d_k_bundle[next], d_h[next], d_status[next], d_integration_param_bundle[next], d_retries[next], d_min_steps[next], @GRADIENT_ARGUMENT@active_chunks[next], next);""".replace(
        "@GRADIENT_ARGUMENT@", gradient_argument_next
    )

    accepted_metric_refresh_current = (
        """
            // Direct-EOM diagnostic refresh: evaluate g4DD at the finalized
            // state while reusing the exact spatial centers selected for this
            // RK trial and every substep.
            numerical_interpolation(
                commondata,
                &numerical_params,
                &spatial_context,
                &numerical_window,
                d_f_bundle[current],
                    d_status[current],
                    d_spatial_stencil_center_i0[current],
                    d_spatial_stencil_center_i2[current],
                    d_temporal_stencil_center[current],
                    d_spatial_mask[current],
                    d_has_excluded_nodes[current],
                    NULL,
                    NULL,
                    d_metric_bundle[current],
                NULL,
                active_chunks[current],
                current);
"""
        if not normalized_eom
        else ""
    )
    accepted_metric_refresh_next = (
        """
                // Direct-EOM diagnostic refresh: evaluate g4DD at the finalized
                // state while reusing the exact spatial centers selected for this
                // RK trial and every substep.
                numerical_interpolation(
                    commondata,
                    &numerical_params,
                    &spatial_context,
                    &numerical_window,
                    d_f_bundle[next],
                    d_status[next],
                    d_spatial_stencil_center_i0[next],
                    d_spatial_stencil_center_i2[next],
                    d_temporal_stencil_center[next],
                    d_spatial_mask[next],
                    d_has_excluded_nodes[next],
                    NULL,
                    NULL,
                    d_metric_bundle[next],
                    NULL,
                    active_chunks[next],
                    next);
"""
        if not normalized_eom
        else ""
    )
    log_energy_evaluation_current = (
        """
            for (long int log_energy_i = 0;
                 log_energy_i < active_chunks[current]; ++log_energy_i) {
                if (d_status[current][log_energy_i] != ACTIVE) {
                    continue;
                } // END IF: skip log-energy for non-active photons
                d_log_energy_bundle[current][log_energy_i] =
                    d_f_bundle[current][4 * BUNDLE_CAPACITY + log_energy_i];
            } // END LOOP: for i over normalized log-energy
"""
        if normalized_eom
        else """
            normal_observer_log_energy(
                d_f_bundle[current],
                d_metric_bundle[current],
                d_status[current],
                d_log_energy_bundle[current],
                active_chunks[current],
                0);
"""
    )
    log_energy_evaluation_next = (
        """
                    for (long int log_energy_i = 0;
                         log_energy_i < active_chunks[next]; ++log_energy_i) {
                        if (d_status[next][log_energy_i] != ACTIVE) {
                            continue;
                        } // END IF: skip log-energy for non-active photons
                        d_log_energy_bundle[next][log_energy_i] =
                            d_f_bundle[next][4 * BUNDLE_CAPACITY + log_energy_i];
                } // END LOOP: for i over normalized log-energy
"""
        if normalized_eom
        else """
                normal_observer_log_energy(
                    d_f_bundle[next],
                    d_metric_bundle[next],
                    d_status[next],
                    d_log_energy_bundle[next],
                    active_chunks[next],
                    0);
        """
    )

    non_terminal_norm_integration_param_setup = (
        """
                    for (long int norm_i = 0; norm_i < chunk_size; ++norm_i) {
                        const long int master_idx = chunk_buffer[0][norm_i];
                        integration_param_bridge[0][norm_i] =
                            non_terminal_crossing_results[master_idx].integration_param;
                    } // END LOOP: load nonterminal plane parameters
                    memcpy(
                        d_integration_param_bundle[0],
                        integration_param_bridge[0],
                        sizeof(double) * chunk_size);
        """
        if normalized_eom
        else ""
    )

    # Evaluate both diagnostics at the refined nonterminal plane state.
    non_terminal_diagnostic_condition = "commondata->perform_normalization_check"
    non_terminal_Lz_kernel_call = ""
    non_terminal_Lz_store = ""
    if axisymmetric_about_z:
        non_terminal_diagnostic_condition += " || commondata->perform_Lz_check"
        non_terminal_Lz_kernel_call = f"""
                    if (commondata->perform_Lz_check) {{
                        axial_angular_momentum_z(
                            d_f_bundle[0],
                            {angular_momentum_metric_argument}d_angular_momentum_bundle,
                            chunk_size,
                            0);
                    }}
"""
        non_terminal_Lz_store = r"""
                        if (commondata->perform_Lz_check) {
                            const double plane_Lz =
                                d_angular_momentum_bundle[norm_i].Lz;
                            if (isfinite(plane_Lz)) {
                                non_terminal_crossing_results[master_idx].Lz = plane_Lz;
                                non_terminal_crossing_results[master_idx].Lz_recorded = true;
                                results_buffer[master_idx].non_terminal_plane_Lz = plane_Lz;
                                results_buffer[master_idx].non_terminal_plane_Lz_time =
                                    non_terminal_crossing_results[master_idx].coordinate_time;
                            } else {
                                non_terminal_skipped_count++;
                            }
                        } // END IF: record L_z at nonterminal plane
"""

    non_terminal_normalization_block = r"""
        //==========================================
        // NONTERMINAL PLANE NORM AND ANGULAR MOMENTUM DIAGNOSTICS
        //==========================================
        if ({NON_TERMINAL_DIAGNOSTIC_CONDITION}) {
            TimeSlotManager non_terminal_norm_tsm;
            long int non_terminal_recorded_count = 0;
            long int non_terminal_skipped_count = 0;
            slot_manager_init(
                &non_terminal_norm_tsm,
                commondata->slot_manager_t_min,
                slot_manager_t_max,
                commondata->slot_manager_delta_t,
                num_rays);

            for (long int norm_ray = 0; norm_ray < num_rays; ++norm_ray) {
                if (!non_terminal_crossing_results[norm_ray].found) {
                    continue;
                }
                non_terminal_recorded_count++;
                const int norm_slot_idx = slot_get_index(
                    &non_terminal_norm_tsm,
                    @NON_TERMINAL_MASK_TIME@);
                if (norm_slot_idx < 0) {
                    non_terminal_skipped_count++;
                    continue;
                }
                slot_add_photon(&non_terminal_norm_tsm, norm_slot_idx, norm_ray);
            } // END LOOP: bin nonterminal plane states

            for (int norm_slot_idx = non_terminal_norm_tsm.num_slots - 1;
                 norm_slot_idx >= 0;
                 --norm_slot_idx) {
                if (non_terminal_norm_tsm.slot_counts[norm_slot_idx] <= 0) {
                    continue;
                }
                if (time_window_manager_numerical_mmap_for_slot(
                        &numerical_window,
                        &non_terminal_norm_tsm,
                        {NON_TERMINAL_GEOMETRY_SLOT}) != TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
                    non_terminal_skipped_count +=
                        non_terminal_norm_tsm.slot_counts[norm_slot_idx];
                    continue;
                }

                while (non_terminal_norm_tsm.slot_counts[norm_slot_idx] > 0) {
                    const long int chunk_size = NRPYMIN(
                        (long int)BUNDLE_CAPACITY,
                        non_terminal_norm_tsm.slot_counts[norm_slot_idx]);
                    slot_remove_chunk(
                        &non_terminal_norm_tsm,
                        norm_slot_idx,
                        chunk_buffer[0],
                        chunk_size);

                    for (int norm_k = 0; norm_k < 9; ++norm_k) {
                        for (long int norm_i = 0; norm_i < chunk_size; ++norm_i) {
                            const long int master_idx = chunk_buffer[0][norm_i];
                            f_bridge[0][norm_k * BUNDLE_CAPACITY + norm_i] =
                                non_terminal_crossing_results[master_idx].state[norm_k];
                        }
                    } // END LOOP: copy nonterminal plane states

                    for (int norm_k = 0; norm_k < 9; ++norm_k) {
                        memcpy(
                            d_f_bundle[0] + norm_k * BUNDLE_CAPACITY,
                            f_bridge[0] + norm_k * BUNDLE_CAPACITY,
                            sizeof(double) * chunk_size);
                    }
                    for (long int norm_i = 0; norm_i < chunk_size; ++norm_i) {
                        // This status buffer is diagnostic scratch space. An
                        // interpolation failure skips only this diagnostic and
                        // never overwrites the photon's physical termination.
                        d_status[0][norm_i] = ACTIVE;
                    }
{NON_TERMINAL_INTEGRATION_PARAM_SETUP}
{NON_TERMINAL_MASK_SETUP}

                    numerical_interpolation(
                        commondata,
                        &numerical_params,
                        &spatial_context,
                        &numerical_window,
                        d_f_bundle[0],
                        d_status[0],
                        {INTERPOLATION_TERMINAL_ARGS}
                        d_metric_bundle[0],
                        NULL,
                        chunk_size,
                        0);

                    if (commondata->perform_normalization_check) {
                        {NORMALIZATION_KERNEL_NAME}(
                            d_f_bundle[0],
                            d_metric_bundle[0],
                            d_norm_bundle,
                            chunk_size,
                            0);
                    }
{NON_TERMINAL_LZ_KERNEL_CALL}

                    for (long int norm_i = 0; norm_i < chunk_size; ++norm_i) {
                        const long int master_idx = chunk_buffer[0][norm_i];
                        if (d_status[0][norm_i] == FAILURE_SPATIAL_INTERPOLATION ||
                            d_status[0][norm_i] == FAILURE_TEMPORAL_INTERPOLATION) {
                            non_terminal_skipped_count++;
                            continue;
                        }
                        if (commondata->perform_normalization_check) {
                            const double current_norm_err = {NORMALIZATION_ERROR_EXPR};
                            const double sidecar_norm_err = {NORMALIZATION_SIDECAR_ERROR_EXPR};
                            if (!isfinite(current_norm_err) || !isfinite(sidecar_norm_err)) {
                                non_terminal_skipped_count++;
                            } else {
                                normalization_abs_non_terminal_by_ray[master_idx] =
                                    sidecar_norm_err;
                                non_terminal_plane_norm_recorded_by_ray[master_idx] = true;
                                non_terminal_crossing_results[master_idx].null_constraint_abs =
                                    sidecar_norm_err;
                                non_terminal_crossing_results[master_idx].null_constraint_recorded = true;
                            }
                        }
{NON_TERMINAL_LZ_STORE}
                    } // END LOOP: evaluate nonterminal plane diagnostics
                } // END WHILE: evaluate nonterminal normalization slot
            } // END LOOP: evaluate nonterminal normalization slots

            slot_manager_free(&non_terminal_norm_tsm);

            printf(
                "Nonterminal plane diagnostic records: %ld; skipped diagnostics: %ld\n",
                non_terminal_recorded_count,
                non_terminal_skipped_count);
        } // END IF: evaluate nonterminal plane diagnostics
"""
    non_terminal_normalization_block = (
        non_terminal_normalization_block.replace(
            "{NON_TERMINAL_GEOMETRY_SLOT}", non_terminal_geometry_slot_idx
        )
        .replace(
            "{NON_TERMINAL_INTEGRATION_PARAM_SETUP}",
            non_terminal_norm_integration_param_setup,
        )
        .replace(
            "{INTERPOLATION_TERMINAL_ARGS}",
            interpolation_terminal_integration_param_args,
        )
        .replace("{NORMALIZATION_KERNEL_NAME}", normalization_kernel_name)
        .replace("{NORMALIZATION_ERROR_EXPR}", normalization_error_expr)
        .replace("{NORMALIZATION_SIDECAR_ERROR_EXPR}", normalization_sidecar_error_expr)
        .replace(
            "{NON_TERMINAL_DIAGNOSTIC_CONDITION}", non_terminal_diagnostic_condition
        )
        .replace("{NON_TERMINAL_LZ_KERNEL_CALL}", non_terminal_Lz_kernel_call)
        .replace("{NON_TERMINAL_LZ_STORE}", non_terminal_Lz_store)
        .replace(
            "@NON_TERMINAL_MASK_TIME@",
            (
                "commondata->t_start"
                if fast_light
                else (
                    "plane_crossing_candidates[norm_ray].non_terminal_candidate.parameter_before"
                    if normalized_eom
                    else "plane_crossing_candidates[norm_ray].non_terminal_candidate.state_before[0]"
                )
            ),
        )
    )

    terminal_diagnostic_condition = "commondata->perform_normalization_check"
    terminal_diagnostic_declarations = ""
    terminal_skip_one = r"""
                    if (commondata->perform_normalization_check) {
                        normalization_skipped_count++;
                    }
"""
    terminal_skip_count = r"""
                    if (commondata->perform_normalization_check) {
                        normalization_skipped_count += diagnostic_count;
                    }
"""
    terminal_angular_momentum_kernel_call = ""
    terminal_angular_momentum_evaluation = ""
    terminal_angular_momentum_report = ""
    if axisymmetric_about_z:
        terminal_diagnostic_condition += " || commondata->perform_Lz_check"
        terminal_diagnostic_declarations = r"""
            long int angular_momentum_skipped_count = 0;
            double max_abs_delta_Lz = NAN;
            long int worst_ray_Lz = -1;
            double worst_initial_Lz = NAN;
            double worst_final_Lz = NAN;
"""
        terminal_skip_one += r"""
                    if (commondata->perform_Lz_check) {
                        angular_momentum_skipped_count++;
                    }
"""
        terminal_skip_count += r"""
                    if (commondata->perform_Lz_check) {
                        angular_momentum_skipped_count += diagnostic_count;
                    }
"""
        terminal_angular_momentum_kernel_call = f"""
                    if (commondata->perform_Lz_check) {{
                        axial_angular_momentum_z(
                            d_f_bundle[0],
                            {angular_momentum_metric_argument}d_angular_momentum_bundle,
                            chunk_size,
                            0);
                    }}
"""
        terminal_angular_momentum_evaluation = r"""
                        if (commondata->perform_Lz_check) {
                            const double initial_Lz =
                                initial_angular_momentum[master_idx].Lz;
                            const double final_Lz =
                                d_angular_momentum_bundle[diagnostic_i].Lz;
                            const double abs_delta_Lz = fabs(final_Lz - initial_Lz);
                            if (!isfinite(initial_Lz) || !isfinite(final_Lz) ||
                                !isfinite(abs_delta_Lz)) {
                                angular_momentum_skipped_count++;
                            } else {
                                final_angular_momentum[master_idx].Lz = final_Lz;
                                if (terminal_crossing_results[master_idx].found) {
                                    terminal_crossing_results[master_idx].Lz = final_Lz;
                                    terminal_crossing_results[master_idx].Lz_recorded = true;
                                }
                                if (worst_ray_Lz < 0 ||
                                    abs_delta_Lz > max_abs_delta_Lz) {
                                    max_abs_delta_Lz = abs_delta_Lz;
                                    worst_ray_Lz = master_idx;
                                    worst_initial_Lz = initial_Lz;
                                    worst_final_Lz = final_Lz;
                                }
                            }
                        } // END IF: evaluate terminal L_z
"""
        terminal_angular_momentum_report = r"""
            if (commondata->perform_Lz_check) {
                printf("\n=================================================\n");
                printf(" AXIAL ANGULAR MOMENTUM DIAGNOSTIC REPORT\n");
                printf("=================================================\n");
                printf(
                    "  Max Absolute Delta L_z: %e (Ray %ld)\n",
                    max_abs_delta_Lz,
                    worst_ray_Lz);
                printf("  Initial L_z for worst ray: %e\n", worst_initial_Lz);
                printf("  Final L_z for worst ray: %e\n", worst_final_Lz);
                printf(
                    "  Skipped terminal L_z diagnostics: %ld\n",
                    angular_momentum_skipped_count);
            }
"""

    terminal_diagnostics_block = r"""
        //==========================================
        // TERMINAL PHOTON DIAGNOSTICS
        //==========================================
        if ({TERMINAL_DIAGNOSTIC_CONDITION}) {
            TimeSlotManager diagnostic_tsm;
            long int normalization_skipped_count = 0;
            double max_err_norm = 0.0;
            long int worst_ray_norm = -1;
            double max_err_norm_excluding_failures = 0.0;
            long int worst_ray_norm_excluding_failures = -1;
{TERMINAL_DIAGNOSTIC_DECLARATIONS}
            slot_manager_init(
                &diagnostic_tsm,
                commondata->slot_manager_t_min,
                slot_manager_t_max,
                commondata->slot_manager_delta_t,
                num_rays);

            for (long int diagnostic_ray = 0;
                 diagnostic_ray < num_rays;
                 ++diagnostic_ray) {
                if (terminal_crossing_results[diagnostic_ray].found) {
                    const int candidate_number =
                        terminal_crossing_results[diagnostic_ray].selected_candidate_number;
                    const photon_plane_crossing_candidates_t *candidates =
                        &plane_crossing_candidates[diagnostic_ray];
                    const photon_plane_crossing_candidate_t *candidate =
                        candidate_number == 1 ? &candidates->terminal_candidate_1 :
                        candidate_number == 2 ? &candidates->terminal_candidate_2 :
                                                &candidates->terminal_candidate_3;
                    accepted_mask_selection[diagnostic_ray] =
                        batch_numerical_mask_selection(
                            commondata, &numerical_params, &numerical_window, candidate);
                } // END IF: use refined crossing selection
                const int diagnostic_slot_idx = slot_get_index(
                    &diagnostic_tsm, {TERMINAL_MASK_SAMPLING_TIME});
                if (diagnostic_slot_idx < 0) {
{TERMINAL_SKIP_ONE}
                    continue;
                }
                slot_add_photon(
                    &diagnostic_tsm, diagnostic_slot_idx, diagnostic_ray);
            } // END LOOP: bin terminal photons by time

            for (int diagnostic_slot_idx = diagnostic_tsm.num_slots - 1;
                 diagnostic_slot_idx >= 0;
                 --diagnostic_slot_idx) {
                if (diagnostic_tsm.slot_counts[diagnostic_slot_idx] <= 0) {
                    continue;
                }

                if (time_window_manager_numerical_mmap_for_slot(
                        &numerical_window,
                        &diagnostic_tsm,
                        {TERMINAL_GEOMETRY_SLOT}) !=
                    TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
                    const long int diagnostic_count =
                        diagnostic_tsm.slot_counts[diagnostic_slot_idx];
{TERMINAL_SKIP_COUNT}
                    continue;
                }

                while (diagnostic_tsm.slot_counts[diagnostic_slot_idx] > 0) {
                    const long int chunk_size = NRPYMIN(
                        (long int)BUNDLE_CAPACITY,
                        diagnostic_tsm.slot_counts[diagnostic_slot_idx]);
                    slot_remove_chunk(
                        &diagnostic_tsm,
                        diagnostic_slot_idx,
                        chunk_buffer[0],
                        chunk_size);

                    for (int diagnostic_k = 0; diagnostic_k < 9; ++diagnostic_k) {
                        for (long int diagnostic_i = 0;
                             diagnostic_i < chunk_size;
                             ++diagnostic_i) {
                            const long int master_idx =
                                chunk_buffer[0][diagnostic_i];
                            f_bridge[0][
                                diagnostic_k * BUNDLE_CAPACITY + diagnostic_i] =
                                terminal_crossing_results[master_idx].found
                                    ? terminal_crossing_results[master_idx].state[diagnostic_k]
                                    : all_photons_host.f[
                                        diagnostic_k * num_rays + master_idx];
                        }
                    }

                    for (int diagnostic_k = 0; diagnostic_k < 9; ++diagnostic_k) {
                        memcpy(
                            d_f_bundle[0] + diagnostic_k * BUNDLE_CAPACITY,
                            f_bridge[0] + diagnostic_k * BUNDLE_CAPACITY,
                            sizeof(double) * chunk_size);
                    }
                    for (long int diagnostic_i = 0;
                         diagnostic_i < chunk_size;
                         ++diagnostic_i) {
                        // Diagnostic interpolation failures do not replace the
                        // photon's physical termination status.
                        d_status[0][diagnostic_i] = ACTIVE;
                    }
{TERMINAL_INTEGRATION_PARAM_SETUP}
{TERMINAL_MASK_SETUP}

                    numerical_interpolation(
                        commondata,
                        &numerical_params,
                        &spatial_context,
                        &numerical_window,
                        d_f_bundle[0],
                        d_status[0],
                        {INTERPOLATION_TERMINAL_ARGS}
                        d_metric_bundle[0],
                        NULL,
                        chunk_size,
                        0);

                    if (commondata->perform_normalization_check) {
                        {NORMALIZATION_KERNEL_NAME}(
                            d_f_bundle[0],
                            d_metric_bundle[0],
                            d_norm_bundle,
                            chunk_size,
                            0);
                    }
{TERMINAL_ANGULAR_MOMENTUM_KERNEL_CALL}
                    for (long int diagnostic_i = 0;
                         diagnostic_i < chunk_size;
                         ++diagnostic_i) {
                        const long int master_idx =
                            chunk_buffer[0][diagnostic_i];
                        if (d_status[0][diagnostic_i] ==
                                FAILURE_SPATIAL_INTERPOLATION ||
                            d_status[0][diagnostic_i] ==
                                FAILURE_TEMPORAL_INTERPOLATION) {
{TERMINAL_SKIP_ONE}
                            continue;
                        }
                        if (commondata->perform_normalization_check) {
                            const double current_norm_err =
                                {NORMALIZATION_ERROR_EXPR};
                            const double sidecar_norm_err =
                                {NORMALIZATION_SIDECAR_ERROR_EXPR};
                            if (!isfinite(current_norm_err) ||
                                !isfinite(sidecar_norm_err)) {
                                normalization_skipped_count++;
                            } else {
                                normalization_abs_by_ray[master_idx] =
                                    sidecar_norm_err;
                                if (terminal_crossing_results[master_idx].found) {
                                    terminal_crossing_results[master_idx].null_constraint_abs =
                                        sidecar_norm_err;
                                    terminal_crossing_results[master_idx].null_constraint_recorded = true;
                                }
                                if (current_norm_err > max_err_norm) {
                                    max_err_norm = current_norm_err;
                                    worst_ray_norm = master_idx;
                                }
                                if (all_photons_host.status[master_idx] !=
                                        STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED &&
                                    all_photons_host.status[master_idx] !=
                                        FAILURE_RKF45_REJECTION_LIMIT &&
                                    all_photons_host.status[master_idx] !=
                                        FAILURE_RKF45_MIN_STEP_LIMIT &&
                                    current_norm_err >
                                        max_err_norm_excluding_failures) {
                                    max_err_norm_excluding_failures =
                                        current_norm_err;
                                    worst_ray_norm_excluding_failures = master_idx;
                                }
                            }
                        } // END IF: evaluate terminal normalization
{TERMINAL_ANGULAR_MOMENTUM_EVALUATION}
                    } // END LOOP: evaluate terminal diagnostics in chunk
                } // END WHILE: terminal diagnostic slot has photons
            } // END LOOP: terminal diagnostic slots

            slot_manager_free(&diagnostic_tsm);

            if (commondata->perform_normalization_check) {
                printf("\n=================================================\n");
                printf(" NORMALIZATION DIAGNOSTIC REPORT\n");
                printf("=================================================\n");
                printf(
                    "  Max Absolute Error, all checked photons: %e (Ray %ld)\n",
                    max_err_norm,
                    worst_ray_norm);
                printf(
                    "  Max Absolute Error, excluding STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED; FAILURE_RKF45_REJECTION_LIMIT; FAILURE_RKF45_MIN_STEP_LIMIT: %e (Ray %ld)\n",
                    max_err_norm_excluding_failures,
                    worst_ray_norm_excluding_failures);
                printf(
                    "  Skipped terminal normalization diagnostics: %ld\n",
                    normalization_skipped_count);
            }
{TERMINAL_ANGULAR_MOMENTUM_REPORT}
        } // END IF: evaluate enabled terminal photon diagnostics
"""
    terminal_diagnostics_block = (
        terminal_diagnostics_block.replace(
            "{TERMINAL_GEOMETRY_SLOT}", terminal_geometry_slot_idx
        )
        .replace("{TERMINAL_DIAGNOSTIC_CONDITION}", terminal_diagnostic_condition)
        .replace(
            "{TERMINAL_DIAGNOSTIC_DECLARATIONS}",
            terminal_diagnostic_declarations,
        )
        .replace(
            "{TERMINAL_MASK_SAMPLING_TIME}",
            terminal_mask_sampling_time.replace("norm_ray", "diagnostic_ray"),
        )
        .replace("{TERMINAL_SKIP_ONE}", terminal_skip_one)
        .replace("{TERMINAL_SKIP_COUNT}", terminal_skip_count)
        .replace("{TERMINAL_INTEGRATION_PARAM_SETUP}", terminal_integration_param_setup)
        .replace(
            "{INTERPOLATION_TERMINAL_ARGS}",
            interpolation_terminal_integration_param_args,
        )
        .replace("{NORMALIZATION_KERNEL_NAME}", normalization_kernel_name)
        .replace(
            "{NORMALIZATION_ERROR_EXPR}",
            normalization_error_expr.replace("norm_i", "diagnostic_i"),
        )
        .replace(
            "{NORMALIZATION_SIDECAR_ERROR_EXPR}",
            normalization_sidecar_error_expr.replace("norm_i", "diagnostic_i"),
        )
        .replace(
            "{TERMINAL_ANGULAR_MOMENTUM_KERNEL_CALL}",
            terminal_angular_momentum_kernel_call,
        )
        .replace(
            "{TERMINAL_ANGULAR_MOMENTUM_EVALUATION}",
            terminal_angular_momentum_evaluation,
        )
        .replace(
            "{TERMINAL_ANGULAR_MOMENTUM_REPORT}",
            terminal_angular_momentum_report,
        )
    )

    angular_momentum_output = (
        r"""
        if (commondata->perform_Lz_check &&
            Lz_bin_path != NULL && Lz_bin_path[0] != '\0') {
            if (write_axial_angular_momentum_z(
                    Lz_bin_path,
                    initial_angular_momentum,
                    final_angular_momentum,
                    num_rays) != 0) {
                exit(1);
            }
        } // END IF: write initial and final L_z
"""
        if axisymmetric_about_z
        else ""
    )
    synthetic_slice_usage_output = r"""
        if (commondata->perform_synthetic_slice_check &&
            synthetic_slice_usage_bin_path != NULL &&
            synthetic_slice_usage_bin_path[0] != '\0') {
            FILE *synthetic_slice_usage_file =
                fopen(synthetic_slice_usage_bin_path, "wb");
            if (synthetic_slice_usage_file == NULL) {
                fprintf(stderr,
                        "ERROR: Could not open synthetic slice-use file '%s' for writing.\n",
                        synthetic_slice_usage_bin_path);
                exit(1);
            }
            for (long int photon_index = 0; photon_index < num_rays; ++photon_index) {
                const uint64_t serialized_photon_index = (uint64_t)photon_index;
                const uint8_t used_lower_endpoint =
                    synthetic_slice_usage_by_photon[photon_index].used_lower_endpoint ? 1U : 0U;
                const uint8_t used_upper_endpoint =
                    synthetic_slice_usage_by_photon[photon_index].used_upper_endpoint ? 1U : 0U;
                const double lower_request_time = used_lower_endpoint
                    ? synthetic_slice_usage_by_photon[photon_index].lower_request_time
                    : NAN;
                const double upper_request_time = used_upper_endpoint
                    ? synthetic_slice_usage_by_photon[photon_index].upper_request_time
                    : NAN;
                if (fwrite(&serialized_photon_index, sizeof(serialized_photon_index), 1,
                           synthetic_slice_usage_file) != 1 ||
                    fwrite(&used_lower_endpoint, sizeof(used_lower_endpoint), 1,
                           synthetic_slice_usage_file) != 1 ||
                    fwrite(&used_upper_endpoint, sizeof(used_upper_endpoint), 1,
                           synthetic_slice_usage_file) != 1 ||
                    fwrite(&lower_request_time, sizeof(lower_request_time), 1,
                           synthetic_slice_usage_file) != 1 ||
                    fwrite(&upper_request_time, sizeof(upper_request_time), 1,
                           synthetic_slice_usage_file) != 1) {
                    fprintf(stderr,
                            "ERROR: Could not write synthetic slice-use record for photon %ld to '%s'.\n",
                            photon_index,
                            synthetic_slice_usage_bin_path);
                    fclose(synthetic_slice_usage_file);
                    exit(1);
                }
            }
            if (fclose(synthetic_slice_usage_file) != 0) {
                fprintf(stderr,
                        "ERROR: Could not close synthetic slice-use file '%s'.\n",
                        synthetic_slice_usage_bin_path);
                exit(1);
            }
        } // END IF: write temporal stencil use
"""
    angular_momentum_cleanup = (
        r"""
        if (commondata->perform_Lz_check) {
            BHAH_FREE(d_angular_momentum_bundle);
            BHAH_FREE(initial_angular_momentum);
            BHAH_FREE(final_angular_momentum);
        } // END IF: free axial-angular-momentum arrays
"""
        if axisymmetric_about_z
        else ""
    )

    def memcpy_cpu(dest: str, src: str, size: str) -> str:
        return f"memcpy({dest}, {src}, {size});"

    def no_sync() -> str:
        return "// CPU memory operations above are complete at this point."

    stream_arg_current = ", current"
    stream_arg_next = ", next"

    def event_manager_call(slot: str) -> str:
        """
        Return the accepted-state event call for one CPU processing slot.

        :param slot: Name of the current or next CPU processing slot.
        :return: C call to the event-detection function for that slot.
        """
        return (
            f"event_detection_manager_kernel(commondata, d_f_bundle[{slot}], "
            f"d_f_start_bundle[{slot}], d_log_energy_bundle[{slot}], "
            f"integration_param_bridge[{slot}], d_integration_param_bundle[{slot}], "
            f"d_status[{slot}], "
            f"d_on_pos_non_terminal_plane_prev[{slot}], "
            f"d_on_pos_terminal_plane_prev[{slot}], "
            "plane_crossing_candidates, "
            f"d_chunk_buffer[{slot}], active_chunks[{slot}], {slot}, "
            f"&horizons, &numerical_window, photon_horizons, {str(fast_light).lower()});"
        )

    event_manager_current = event_manager_call("current")
    event_manager_next = event_manager_call("next")

    def _diagnostic_mask_setup(selection_expression: str) -> str:
        """
        Rebuild fixed spatial masks for a diagnostic photon chunk.

        :param selection_expression: C expression selecting the saved spatial and temporal nodes.
        :return: C loop filling the interpolation centers, masks, and diagnostic statuses.
        """
        return r"""
                    for (long int mask_i = 0; mask_i < chunk_size; ++mask_i) {
                        const long int master_idx = chunk_buffer[0][mask_i];
                        const batch_spatial_mask_selection_struct selection = @SELECTION@;
                        d_spatial_stencil_center_i0[0][mask_i] = selection.spatial_center.i0;
                        d_spatial_stencil_center_i2[0][mask_i] = selection.spatial_center.i2;
                        d_temporal_stencil_center[0][mask_i] = selection.temporal_center;
                        d_has_excluded_nodes[0][mask_i] = false;
                        if (time_window_manager_numerical_build_spatial_mask(
                                commondata, &numerical_params, &numerical_window,
                                &selection.spatial_center, selection.phi_plane,
                                selection.temporal_center, selection.sampling_time,
                                d_spatial_mask[0] + (size_t)mask_i * spatial_mask_node_count,
                                &d_has_excluded_nodes[0][mask_i]) !=
                            TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
                            d_status[0][mask_i] = FAILURE_SPATIAL_INTERPOLATION;
                        } // END IF: diagnostic mask unavailable
                    } // END LOOP: rebuild diagnostic spatial masks
""".replace("@SELECTION@", selection_expression)

    non_terminal_normalization_block = non_terminal_normalization_block.replace(
        "{NON_TERMINAL_MASK_SETUP}",
        _diagnostic_mask_setup(
            "batch_numerical_mask_selection(commondata, &numerical_params, "
            "&numerical_window, &plane_crossing_candidates[master_idx].non_terminal_candidate)"
        ),
    )
    terminal_diagnostics_block = terminal_diagnostics_block.replace(
        "{TERMINAL_MASK_SETUP}",
        _diagnostic_mask_setup("accepted_mask_selection[master_idx]"),
    )
    mask_diagnostic_condition = "commondata->perform_normalization_check"
    if axisymmetric_about_z:
        mask_diagnostic_condition += " || commondata->perform_Lz_check"
    mask_selection_helpers = r"""
typedef struct {
    NumericalSpatialStencilCenter spatial_center;
    uint64_t temporal_center;
    int phi_plane;
    double sampling_time;
} batch_spatial_mask_selection_struct;

/**
 * Select mask nodes at the validated observer or a refined crossing anchor.
 *
 * @param[in] commondata Photon evolution and interpolation parameters.
 * @param[in] numerical_params Numerical spatial grid parameters.
 * @param[in] numerical_window Retained numerical slice times.
 * @param[in] candidate Saved crossing anchor, or NULL for the observer.
 * @return Spatial center, temporal center, azimuthal plane, and sampling time.
 *
 * @pre Observer interpolation or crossing refinement has validated the anchor.
 */
static batch_spatial_mask_selection_struct batch_numerical_mask_selection(
    const commondata_struct *commondata, const params_struct *numerical_params,
    const NumericalTimeWindowManager *numerical_window,
    const photon_plane_crossing_candidate_t *candidate) {
    const REAL cartesian[3] = {
        candidate == NULL ? commondata->observer_x : candidate->state_before[1],
        candidate == NULL ? commondata->observer_y : candidate->state_before[2],
        candidate == NULL ? commondata->observer_z : candidate->state_before[3]};
    REAL native[3];
    int automatic_center[3], selected_center[3];
    (void)time_window_manager_numerical_resolve_spatial_target_and_stencil(
        numerical_params, cartesian, NULL, native, automatic_center, selected_center);
    const double sampling_time = candidate == NULL ? commondata->t_start : @ANCHOR_TIME@;
    const double midpoint = candidate == NULL ? commondata->t_start : @CROSSING_MIDPOINT@;
    const batch_spatial_mask_selection_struct selection = {
        {selected_center[0], selected_center[2]},
        time_window_manager_numerical_nearest_selected_slice_for_time(
            numerical_window, midpoint),
        native[1] < 0.0 || native[1] >= 3.14159265358979323846 - 1.0e-12 ? 0 : 1,
        sampling_time};
    return selection;
} // END FUNCTION: select diagnostic stencil nodes
""".replace(
        "@ANCHOR_TIME@",
        (
            "commondata->t_start"
            if fast_light
            else (
                "candidate->parameter_before"
                if normalized_eom
                else "candidate->state_before[0]"
            )
        ),
    ).replace(
        "@CROSSING_MIDPOINT@",
        (
            "commondata->t_start"
            if fast_light
            else (
                "0.5 * (candidate->parameter_before + candidate->parameter_after)"
                if normalized_eom
                else "candidate->state_before[0] + 0.5 * "
                "(candidate->parameter_after - candidate->parameter_before) * candidate->state_before[4]"
            )
        ),
    )

    puncture_helpers = puncture_diagnostics()
    body = rf"""
{gradient_validation}
    // Read requested horizon positions before evolving any photons in this tile.
    apparent_horizon_file_t horizons = {{0}};
    puncture_result_t *photon_horizons = NULL;
    puncture_diagnostics_t puncture_data = puncture_load_diagnostics(
        commondata->apparent_horizon_bin_path[0] ? "" : commondata->puncture_diagnostics_dir);

    // Initialize the caller's photon output records before event detection.
    for (long int i = 0; i < num_rays; ++i) {{
        results_buffer[i] = (blueprint_data_t){{0}};
        results_buffer[i].termination_type = FAILURE_GENERIC;
        results_buffer[i].non_terminal_plane_Lz = NAN;
        results_buffer[i].non_terminal_plane_Lz_time = NAN;
        results_buffer[i].non_terminal_plane_signed_distance = NAN;
    }} // END LOOP: initialize deterministic numerical-batch result records

    //==========================================
    // 1. HOST ALLOCATION
    //==========================================

    // The master host-side Structure of Arrays (SoA) tracking all photons $f^\mu$.
    PhotonStateSoA all_photons_host = {{0}};

    // Saved crossing steps and their later plane results, indexed by photon.
    photon_plane_crossing_candidates_t *plane_crossing_candidates = NULL;
    non_terminal_plane_crossing_result_t *non_terminal_crossing_results = NULL;
    terminal_plane_crossing_result_t *terminal_crossing_results = NULL;
    {malloc_pinned}(plane_crossing_candidates,
        sizeof(*plane_crossing_candidates) * num_rays);
    {malloc_pinned}(non_terminal_crossing_results,
        sizeof(*non_terminal_crossing_results) * num_rays);
    {malloc_pinned}(terminal_crossing_results,
        sizeof(*terminal_crossing_results) * num_rays);
    if (num_rays > 0 && (plane_crossing_candidates == NULL ||
                         non_terminal_crossing_results == NULL ||
                         terminal_crossing_results == NULL)) {{
        fprintf(stderr, "ERROR: could not allocate photon plane crossing records.\n");
        exit(EXIT_FAILURE);
    }}
    for (long int ray = 0; ray < num_rays; ++ray) {{
        plane_crossing_candidates[ray] = (photon_plane_crossing_candidates_t){{0}};
        non_terminal_crossing_results[ray] = (non_terminal_plane_crossing_result_t){{0}};
        terminal_crossing_results[ray] = (terminal_plane_crossing_result_t){{0}};
        plane_crossing_candidates[ray].non_terminal_candidate.photon_index = ray;
        plane_crossing_candidates[ray].non_terminal_candidate.plane_type =
            NON_TERMINAL_PLANE_EVENT;
        photon_plane_crossing_candidate_t *terminal_candidates[3] = {{
            &plane_crossing_candidates[ray].terminal_candidate_1,
            &plane_crossing_candidates[ray].terminal_candidate_2,
            &plane_crossing_candidates[ray].terminal_candidate_3
        }};
        for (int candidate = 0; candidate < 3; ++candidate) {{
            terminal_candidates[candidate]->photon_index = ray;
            terminal_candidates[candidate]->plane_type = TERMINAL_PLANE_EVENT;
        }}
        non_terminal_crossing_results[ray].null_constraint_abs = NAN;
        non_terminal_crossing_results[ray].Lz = NAN;
        terminal_crossing_results[ray].null_constraint_abs = NAN;
        terminal_crossing_results[ray].Lz = NAN;
    }} // END LOOP: initialize crossing records

    // {pin_comment} the state vector $f^\mu$.
    {malloc_pinned}(all_photons_host.f, sizeof(double) * 9 * num_rays);
    // {pin_comment} the integration parameter.
    {malloc_pinned}(all_photons_host.integration_param, sizeof(double) * num_rays);
    // {pin_comment} individual integration step sizes $h$.
    {malloc_pinned}(all_photons_host.h, sizeof(double) * num_rays);
    // {pin_comment} the trajectory termination status.
    {malloc_pinned}(all_photons_host.status, sizeof(termination_type_t) * num_rays);
    // {pin_comment} the number of step-size rejections.
    {malloc_pinned}(all_photons_host.rejection_retries, sizeof(int) * num_rays);
    {malloc_pinned}(all_photons_host.consecutive_min_steps, sizeof(int) * num_rays);
    if (all_photons_host.consecutive_min_steps == NULL) {{
        fprintf(stderr, "ERROR: could not allocate photon minimum-step counters.\n");
        exit(EXIT_FAILURE);
    }} // END IF: minimum-step allocation failed
    // {pin_comment} the previous nonterminal plane boundary state.
    {malloc_pinned}(all_photons_host.on_positive_side_of_non_terminal_plane_prev, sizeof(bool) * num_rays);
    // {pin_comment} the previous terminal-plane boundary state.
    {malloc_pinned}(all_photons_host.on_positive_side_of_terminal_plane_prev, sizeof(bool) * num_rays);
    // CPU-only numerical integration: direct commondata access and no extra execution-buffer setup.

    //==========================================
    // DOUBLE-BUFFERED BRIDGE ARRAYS
    //==========================================
    // Extraction buffer used by the TimeSlotManager to map sparse indices to contiguous execution blocks.
    long int *chunk_buffer[2];
    // Bridge array staging the state vector $f^\mu$ for memory transfers.
    double *f_bridge[2];
    // Bridge array staging the integration parameter for memory transfers.
    double *integration_param_bridge[2];
    // Bridge array staging the current integration step size $h$ for memory transfers.
    double *h_bridge[2];
    // Bridge array staging the current trajectory termination status for memory transfers.
    termination_type_t *status_bridge[2];
    // Incoming status for each trial, saved before retry classification.
    termination_type_t *trial_status[2];
    // Bridge array staging the number of step-size rejections for memory transfers.
    int *retries_bridge[2];
    int *min_steps_bridge[2];
    // Bridge array staging the previous nonterminal plane boundary side flag for memory transfers.
    bool *on_pos_non_terminal_plane_prev_bridge[2];
    // Bridge array staging the previous terminal-plane boundary side flag for memory transfers.
    bool *on_pos_terminal_plane_prev_bridge[2];
    // Bridge array staging the nonterminal plane event lock for memory transfers.
    // Bridge array staging the terminal-plane event lock for memory transfers.

    //==========================================
    // DOUBLE-BUFFERED CPU SCRATCHPADS
    //==========================================
    // Scratchpad tracking the current state vector $f^\mu$ bounding the RKF45 step.
    double *d_f_bundle[2];
    // Scratchpad locking the anchor state vector $f_{{start}}$ to calculate the final stage update.
    double *d_f_start_bundle[2];
    // Scratchpad tracking the intermediate cumulative RKF45 stage updates.
    double *d_f_temp_bundle[2];
    // Scratchpad persisting the symmetric metric tensor $g_{{\mu\nu}}$.
    double *d_metric_bundle[2];
    // Scratchpad carrying the common upper-only log-energy measure.
    double *d_log_energy_bundle[2];
    // Store the 40 geometry components used by the photon RHS.
    double *d_rhs_geometry_bundle[2];
    // Derivative tensor storing $\dot{{f}}^\mu$ across all 6 intermediate RKF45 stages.
    double *d_k_bundle[2];
    // Array regulating active integration step sizing $h$.
    double *d_h[2];
    // Array regulating integration-parameter progress.
    double *d_integration_param_bundle[2];
    // Array holding the current trajectory status limits.
    termination_type_t *d_status[2];
    // Array tracking sequential error rejections per photon.
    int *d_retries[2];
    int *d_min_steps[2];
    // Array flagging the previous nonterminal plane boundary side.
    bool *d_on_pos_non_terminal_plane_prev[2];
    // Array flagging the previous terminal-plane boundary side.
    bool *d_on_pos_terminal_plane_prev[2];
    // Array guarding nonterminal-plane intersection coordinates from multi-trigger overwrites.
    // Array guarding the terminal-plane intersection coordinates from multi-trigger overwrites.
    // Array carrying the absolute master indices $m_{{idx}}$ mapping the execution chunk.
    long int *d_chunk_buffer[2];
{spatial_center_declarations}
{gradient_declarations}

    // Loop iterator for instantiating the double-buffered operational arrays.
    for (int s = 0; s < 2; ++s) {{
        // {bridge_alloc_comment}
        {malloc_pinned}(chunk_buffer[s], sizeof(long int) * BUNDLE_CAPACITY); // Pin chunk buffers.
        {malloc_pinned}(f_bridge[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Pin $f^\mu$ bridges.
        {malloc_pinned}(integration_param_bridge[s], sizeof(double) * BUNDLE_CAPACITY); // Pin integration-parameter bridges.
        {malloc_pinned}(h_bridge[s], sizeof(double) * BUNDLE_CAPACITY); // Pin $h$ bridges.
        {malloc_pinned}(status_bridge[s], sizeof(termination_type_t) * BUNDLE_CAPACITY); // Pin status bridges.
        {malloc_pinned}(trial_status[s], sizeof(termination_type_t) * BUNDLE_CAPACITY);
        if (trial_status[s] == NULL) {{
            fprintf(stderr, "ERROR: could not allocate photon trial statuses.\n");
            exit(EXIT_FAILURE);
        }} // END IF: trial status allocation failed
        {malloc_pinned}(retries_bridge[s], sizeof(int) * BUNDLE_CAPACITY); // Pin retries bridges.
        {malloc_pinned}(min_steps_bridge[s], sizeof(int) * BUNDLE_CAPACITY);
        {malloc_pinned}(on_pos_non_terminal_plane_prev_bridge[s], sizeof(bool) * BUNDLE_CAPACITY); // Pin nonterminal-plane flag bridges.
        {malloc_pinned}(on_pos_terminal_plane_prev_bridge[s], sizeof(bool) * BUNDLE_CAPACITY); // Pin source flag bridges.

        // {scratch_alloc_comment}
        {malloc_device}(d_f_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate $f^\mu$ scratchpad.
        {malloc_device}(d_f_start_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate $f_{{start}}$ scratchpad.
        {malloc_device}(d_f_temp_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate temporary stage scratchpad.
        {malloc_device}(d_metric_bundle[s], sizeof(double) * 10 * BUNDLE_CAPACITY); // Allocate $g_{{\mu\nu}}$ scratchpad.
        {malloc_device}(d_log_energy_bundle[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate common log-energy scratchpad.
        {malloc_device}(d_rhs_geometry_bundle[s], sizeof(double) * 40 * BUNDLE_CAPACITY); // Allocate the geometry scratchpad.
        {malloc_device}(d_k_bundle[s], sizeof(double) * 6 * 9 * BUNDLE_CAPACITY); // Allocate derivative scratchpad.
        {malloc_device}(d_h[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate $h$ scratchpad.
        {malloc_device}(d_integration_param_bundle[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate integration-parameter scratchpad.
        {malloc_device}(d_status[s], sizeof(termination_type_t) * BUNDLE_CAPACITY); // Allocate status scratchpad.
        {malloc_device}(d_retries[s], sizeof(int) * BUNDLE_CAPACITY); // Allocate retries scratchpad.
        {malloc_device}(d_min_steps[s], sizeof(int) * BUNDLE_CAPACITY);
        if (min_steps_bridge[s] == NULL || d_min_steps[s] == NULL) {{
            fprintf(stderr, "ERROR: could not allocate minimum-step chunk counters.\n");
            exit(EXIT_FAILURE);
        }} // END IF: chunk counter allocation failed
        {malloc_device}(d_on_pos_non_terminal_plane_prev[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate nonterminal-plane flag scratchpad.
        {malloc_device}(d_on_pos_terminal_plane_prev[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate source flag scratchpad.
        {malloc_device}(d_chunk_buffer[s], sizeof(long int) * BUNDLE_CAPACITY); // Allocate chunk mapping scratchpad.
{spatial_center_allocations}
{gradient_allocations}
    }} // END LOOP: for s over 2

    // Scratchpad array holding the initialization constraint outputs. It is
    // always allocated because initialization must validate every ray before
    // RKF45 starts, independent of the optional terminal diagnostic.
    normalization_constraint_t *d_initial_norm_bundle = NULL;
    {malloc_device}(d_initial_norm_bundle, sizeof(normalization_constraint_t) * BUNDLE_CAPACITY);

    // Scratchpad array holding the terminal normalization diagnostic outputs.
    normalization_constraint_t *d_norm_bundle = NULL;
    // Host array storing one absolute direct-scale-equivalent null-constraint
    // error per photon in master-ray order.
    double *normalization_abs_by_ray = NULL;
    // Host array storing the nonterminal plane null-constraint error by photon.
    double *normalization_abs_non_terminal_by_ray = NULL;
    bool *non_terminal_plane_norm_recorded_by_ray = NULL;
    // Host array storing first synthetic temporal-stencil use times per photon.
    synthetic_slice_usage_t *synthetic_slice_usage_by_photon = NULL;

    if (commondata->perform_normalization_check) {{
        {malloc_device}(d_norm_bundle, sizeof(normalization_constraint_t) * BUNDLE_CAPACITY); // Allocate terminal normalization scratchpad.
        {malloc_pinned}(normalization_abs_by_ray, sizeof(double) * num_rays); // Allocate per-photon direct-scale-equivalent normalization sidecar buffer.
        {malloc_pinned}(normalization_abs_non_terminal_by_ray, sizeof(double) * num_rays); // Allocate sparse nonterminal normalization values.
        {malloc_pinned}(non_terminal_plane_norm_recorded_by_ray, sizeof(bool) * num_rays);
        for (long int norm_init_i = 0; norm_init_i < num_rays; ++norm_init_i) {{
            normalization_abs_by_ray[norm_init_i] = NAN; // Marks photons whose terminal normalization was not evaluated.
            normalization_abs_non_terminal_by_ray[norm_init_i] = NAN;
            non_terminal_plane_norm_recorded_by_ray[norm_init_i] = false;
        }} // END LOOP: for norm_init_i over num_rays
    }} // END IF: commondata->perform_normalization_check to allocate normalization scratchpad
    if (commondata->perform_synthetic_slice_check) {{
        {malloc_pinned}(synthetic_slice_usage_by_photon, sizeof(synthetic_slice_usage_t) * num_rays);
        for (long int usage_init_i = 0; usage_init_i < num_rays; ++usage_init_i) {{
            synthetic_slice_usage_by_photon[usage_init_i].used_lower_endpoint = false;
            synthetic_slice_usage_by_photon[usage_init_i].lower_request_time = NAN;
            synthetic_slice_usage_by_photon[usage_init_i].used_upper_endpoint = false;
            synthetic_slice_usage_by_photon[usage_init_i].upper_request_time = NAN;
        }} // END LOOP: initialize synthetic temporal-stencil usage records
    }} // END IF: allocate synthetic temporal-stencil usage records
{angular_momentum_allocation}

    // Event-detection kernels write final physical plane intersections directly to results_buffer.

    // Host-bound struct managing temporal binning of photon trajectories $x^\mu$.
    TimeSlotManager tsm;
    // The slot-manager upper bound is exclusive, so it must sit slightly above
    // t_start; exact equality would reject the initial photon state.
    const double slot_manager_t_max = commondata->t_start + 1.0e-5;
    // The main slot range is controlled by commondata, not by combined-file metadata.
    slot_manager_init(
        &tsm,
        commondata->slot_manager_t_min,
        slot_manager_t_max,
        commondata->slot_manager_delta_t,
        num_rays);

    // Numerical spacetime window manager and metadata-fed grid params.
    NumericalTimeWindowManager numerical_window;
    time_window_manager_numerical_set_inert(&numerical_window);

    // Seed the runtime params with generator-time defaults so coordinate-map
    // parameters like AMPL, SINHW, and any slope terms remain valid after
    // metadata overlay.
    commondata_struct commondata_for_params_defaults = *commondata;
    griddata_struct dummy_griddata[MAXNUMGRIDS];
    params_struct_set_to_default(&commondata_for_params_defaults, dummy_griddata);
    params_struct numerical_params = dummy_griddata[0].params;

    if (commondata->numerical_spacetime_bin_path[0] == '\0') {{
        fprintf(stderr,
                "ERROR: commondata->numerical_spacetime_bin_path is empty. "
                "The numerical photon batch integrator requires a validated combined .bin file path.\n");
        slot_manager_free(&tsm);
        exit(1);
    }} // END IF: numerical_spacetime_bin_path was empty

    if (time_window_manager_numerical_init(
            &numerical_window,
            commondata->numerical_spacetime_bin_path,
            commondata,
            commondata->numerical_spacetime_temporal_interp_half_width,
            &numerical_params) != TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
        fprintf(stderr,
                "ERROR: failed to initialize numerical time-window manager from '%s'.\n",
                commondata->numerical_spacetime_bin_path);
        slot_manager_free(&tsm);
        exit(1);
    }} // END IF: numerical time-window manager initialization failed
    const size_t spatial_mask_width =
        2 * (size_t)commondata->numerical_spacetime_spatial_interp_half_width + 1;
    const size_t spatial_mask_node_count = spatial_mask_width * spatial_mask_width;
    const size_t spatial_mask_ray_capacity = NRPYMIN((long int)BUNDLE_CAPACITY, num_rays);
    for (int s = 0; s < 2; ++s) {{
        {malloc_device}(d_spatial_mask[s], spatial_mask_node_count * spatial_mask_ray_capacity);
        if (d_spatial_mask[s] == NULL) {{
            fprintf(stderr, "ERROR: could not allocate photon spatial masks.\n");
            exit(EXIT_FAILURE);
        }} // END IF: spatial mask allocation failed
    }} // END LOOP: allocate photon spatial masks
    {temporal_cap_disable}
    printf("Numerical spacetime first stored slice time: %.15e\n",
           (double)commondata->t_numerical_initial);

    if (numerical_params.Nxx{phi_dim} != 2) {{
        fprintf(stderr,
                "ERROR: numerical spatial interpolation expects exactly two stored phi planes in native dimension {phi_dim}; got Nxx{phi_dim}=%d.\n",
                numerical_params.Nxx{phi_dim});
        time_window_manager_numerical_free(&numerical_window);
        slot_manager_free(&tsm);
        exit(1);
    }} // END IF: stored phi-plane count was incompatible
    // Cache metric outlier flags for spatial stencils wider than one node.
    numerical_window.prepare_spatial_outlier_flags = spatial_mask_width > 1;

    azimuthal_symmetry_spatial_lagrange_context_struct spatial_context;
    spatial_context.stored_phi_samples[0] =
        numerical_params.xxmin{phi_dim} + 0.5 * numerical_params.dxx{phi_dim};
    spatial_context.stored_phi_samples[1] =
        numerical_params.xxmin{phi_dim} + 1.5 * numerical_params.dxx{phi_dim};

    //==========================================
    // 2. INITIALIZATION PHASE
    //==========================================
    // Map the initial slot before any numerical interpolation. The mapped
    // window is shared by the temporary observer state and all later rays.
    const int initial_slot_idx = slot_get_index(&tsm, commondata->t_start);
    if (initial_slot_idx < 0) {{
        fprintf(stderr,
                "ERROR: initial photon time t_start=%e is outside the configured TimeSlotManager range.\n",
                (double)commondata->t_start);
        time_window_manager_numerical_free(&numerical_window);
        slot_manager_free(&tsm);
        exit(1);
    }} // END IF: initial photon time outside bounds

    if (time_window_manager_numerical_mmap_for_slot(
            &numerical_window, &tsm, initial_slot_idx) !=
        TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
        fprintf(stderr,
                "ERROR: failed to map numerical time window for initial slot %d at t_start=%e.\n",
                initial_slot_idx,
                (double)commondata->t_start);
        time_window_manager_numerical_free(&numerical_window);
        slot_manager_free(&tsm);
        exit(1);
    }} // END IF: initial slot time window mapping

    if (apparent_horizons_open(&horizons, &numerical_window, commondata->apparent_horizon_bin_path)) {{
        fprintf(stderr, "ERROR: cannot read synchronized apparent horizons from %s.\n",
                commondata->apparent_horizon_bin_path);
        exit(EXIT_FAILURE);
    }} // END IF: invalid apparent horizon file
    if (horizons.bytes != NULL) {{
        photon_horizons = calloc((size_t)num_rays, sizeof(*photon_horizons));
        if (photon_horizons == NULL) {{
            fprintf(stderr, "ERROR: cannot allocate photon horizon results.\n");
            exit(EXIT_FAILURE);
        }} // END IF: allocation failed
    }} // END IF: binary horizon checking enabled

    // Interpolate the metric once at the common observer event. This temporary
    // state uses the observer position and t_start, but does not represent a
    // renderable ray. The observer initializer consumes only this ten-component
    // covariant metric and constructs the tetrad once on the host.
    double observer_metric[10];
    double observer_tetrad[4][4];
    for (int observer_component = 0; observer_component < 9; ++observer_component) {{
        d_f_bundle[0][observer_component * BUNDLE_CAPACITY] = 0.0;
    }} // END LOOP: clear temporary observer state
    d_f_bundle[0][0 * BUNDLE_CAPACITY] = {observer_state_coordinate_time};
    d_f_bundle[0][1 * BUNDLE_CAPACITY] = commondata->observer_x;
    d_f_bundle[0][2 * BUNDLE_CAPACITY] = commondata->observer_y;
    d_f_bundle[0][3 * BUNDLE_CAPACITY] = commondata->observer_z;
{observer_state_integration_setup}

    // Observer interpolation uses a temporary status because failure here
    // prevents every ray in the batch from being initialized.
    termination_type_t observer_interpolation_status = ACTIVE;
    numerical_interpolation(
        commondata,
        &numerical_params,
        &spatial_context,
        &numerical_window,
        d_f_bundle[0],
        &observer_interpolation_status,
        {interpolation_initial_integration_param_args}
        d_metric_bundle[0],
        NULL,
        1,
        0);

    for (int observer_metric_component = 0; observer_metric_component < 10;
         ++observer_metric_component) {{
        observer_metric[observer_metric_component] =
            d_metric_bundle[0][observer_metric_component * BUNDLE_CAPACITY];
    }} // END LOOP: for observer_metric_component over metric

    if (observer_interpolation_status == FAILURE_SPATIAL_INTERPOLATION ||
        observer_interpolation_status == FAILURE_TEMPORAL_INTERPOLATION) {{
        fprintf(stderr,
                "ERROR: observer interpolation failed with status %d.\n",
                (int)observer_interpolation_status);
        time_window_manager_numerical_free(&numerical_window);
        slot_manager_free(&tsm);
        exit(1);
    }} // END IF: observer interpolation failed

    if ({mask_diagnostic_condition}) {{
        {malloc_pinned}(accepted_mask_selection, sizeof(*accepted_mask_selection) * num_rays);
        if (accepted_mask_selection == NULL) {{
            fprintf(stderr, "ERROR: could not allocate accepted photon stencil selections.\n");
            exit(EXIT_FAILURE);
        }} // END IF: accepted selection allocation failed
        const batch_spatial_mask_selection_struct observer_selection =
            batch_numerical_mask_selection(
                commondata, &numerical_params, &numerical_window, NULL);
        for (long int ray = 0; ray < num_rays; ++ray) {{
            accepted_mask_selection[ray] = observer_selection;
        }} // END LOOP: seed observer stencil selections
    }} // END IF: delayed metric diagnostics enabled

    // Validate the observer tetrad and fill the initial photon states with
    // direct p^mu values before the first RKF45 attempt.
    {set_initial_conditions_call}

    long int num_batches = (num_rays + BUNDLE_CAPACITY - 1) / BUNDLE_CAPACITY; // Total integer calculation defining total iterative blocks required to process all photon indices.
    double max_initial_constraint_error = 0.0;
    long int worst_initial_constraint_ray = -1;

    for (long int init_batch = 0; init_batch < num_batches; ++init_batch) {{ // Loop iterator $init_batch$ for evaluating the initialization constraint across sequential blocks.
        long int start_idx = init_batch * BUNDLE_CAPACITY; // Absolute starting index mapped to the master SoA for the current initialization batch.
        long int chunk_size = NRPYMIN((long int)BUNDLE_CAPACITY, num_rays - start_idx); // Dynamically sized operational boundary ensuring the active chunk does not exceed total trajectories.

        for (int init_i = 0; init_i < chunk_size; ++init_i) {{ // Loop index $init_i$ iterating over the specific initialization batch elements to pack the bridge.
            long int master_idx = start_idx + init_i; // Computes the absolute master index $m_{{idx}}$ tracking the photon within the global array.
            for (int init_k = 0; init_k < 9; ++init_k) {{ // Loop index $init_k$ iterating over the 9 tensor components of the state vector $f^\mu$.
                f_bridge[0][init_k * BUNDLE_CAPACITY + init_i] = all_photons_host.f[init_k * num_rays + master_idx]; // Assigns the active tensor state component to the primary bridge.
            }} // END LOOP: for init_k over 9
        }} // END LOOP: for init_i over chunk_size

        for (int c_k = 0; c_k < 9; ++c_k) {{ // Loop index $c_k$ orchestrating the memory transfer of the 9 state vector $f^\mu$ components.
            {memcpy_cpu("d_f_bundle[0] + c_k * BUNDLE_CAPACITY", "f_bridge[0] + c_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size")}
        }} // END LOOP: for c_k over 9
{initial_integration_param_setup}

        // Reuse the one observer metric for every ray in this chunk. No
        // per-ray observer interpolation or p^0 quadratic solve is performed.
        for (int metric_component = 0; metric_component < 10; ++metric_component) {{
            for (long int init_i = 0; init_i < chunk_size; ++init_i) {{
                d_metric_bundle[0][metric_component * BUNDLE_CAPACITY + init_i] =
                    observer_metric[metric_component];
            }} // END LOOP: for init_i over current chunk
        }} // END LOOP: replicate all ten metric components

        //==========================================
        // DIAGNOSTIC PROBE: REUSED METRIC INTEGRITY CHECK
        //==========================================
        // Confirms that the one observer metric was replicated without
        // corruption before momentum conversion and constraint evaluation.

        double *metric_diag_bridge; // Pointer storing temporary metric data to validate the interpolation sequence.
        // Memory allocation: temporary bridge for the replicated metric check.
        {malloc_pinned}(metric_diag_bridge, sizeof(double) * 10 * BUNDLE_CAPACITY);

        for (int m_k = 0; m_k < 10; ++m_k) {{
            // Loop index $m_k$ orchestrating memory transfer of the 10 metric tensor $g_{{\mu\nu}}$ components.
            {memcpy_cpu("metric_diag_bridge + m_k * BUNDLE_CAPACITY", "d_metric_bundle[0] + m_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size")}
        }} // END LOOP: for m_k over 10
        {no_sync()}

        long int metric_nan_count = 0; // Accumulator tracking the total number of metric tensor evaluations containing non-finite values.
        for (int m_diag_i = 0; m_diag_i < chunk_size; ++m_diag_i) {{ // Loop iterator $m_diag_i$ scanning each trajectory within the current initialization chunk.
            bool m_has_nan = false; // Boolean flag indicating if the specific metric tensor $g_{{\mu\nu}}$ contains a NaN or Inf value.
            for (int m_diag_k = 0; m_diag_k < 10; ++m_diag_k) {{ // Loop index $m_diag_k$ iterating over the 10 independent components of the symmetric metric tensor $g_{{\mu\nu}}$.
                if (isnan(metric_diag_bridge[m_diag_k * BUNDLE_CAPACITY + m_diag_i]) ||
                    isinf(metric_diag_bridge[m_diag_k * BUNDLE_CAPACITY + m_diag_i])) {{
                    m_has_nan = true; // Flags the trajectory metric state as invalid due to a non-finite value.
                    break; // Terminates the tensor component loop early to avoid unnecessary work upon detecting a failure.
                }} // END IF: check for NaN or Inf
            }} // END LOOP: for m_diag_k over 10
            if (m_has_nan) metric_nan_count++; // Increments the total accumulation of corrupted metric tensor evaluations.
        }} // END LOOP: for m_diag_i over chunk_size

        if (metric_nan_count > 0) {{
            fprintf(stderr,
                    "ERROR: Init Batch %ld: %ld rays have invalid numerical metric "
                    "G_mu_nu before momentum conversion.\n",
                    init_batch,
                    metric_nan_count);
            {free_pinned}(metric_diag_bridge);
            time_window_manager_numerical_free(&numerical_window);
            slot_manager_free(&tsm);
            exit(1);
        }} // END IF: metric_nan_count > 0 to fail
        // Memory Free: Purges the diagnostic bridge utilized for metric integrity checks.
        {free_pinned}(metric_diag_bridge);

        // The tetrad initializer already supplied the complete direct p^mu.
        {normalized_momentum_conversion}
{initial_angular_momentum_evaluation}

        // Evaluate the applicable initial constraint after normalized
        // conversion, if requested. For direct evolution this is g_mu_nu p^mu
        // p^nu and must be zero; for normalized evolution it is the spatial
        // normalized constraint and must be one.
        {normalization_kernel_name}(
            d_f_bundle[0],
            d_metric_bundle[0],
            d_initial_norm_bundle,
            chunk_size,
            0);

        const double initial_constraint_target = {initial_constraint_target};
        const double initial_constraint_tolerance =
            1.0e-9 * fmax(1.0, fabs(initial_constraint_target));
        for (long int constraint_i = 0; constraint_i < chunk_size; ++constraint_i) {{
            const double constraint_value = d_initial_norm_bundle[constraint_i].C;
            const double constraint_error =
                fabs(constraint_value - initial_constraint_target);
            const long int master_idx = start_idx + constraint_i;
            if (!isfinite(constraint_value) || !isfinite(constraint_error)) {{
                fprintf(
                    stderr,
                    "ERROR: initial photon constraint is non-finite for ray %ld.\n",
                    master_idx);
                time_window_manager_numerical_free(&numerical_window);
                slot_manager_free(&tsm);
                exit(1);
            }} // END IF: initial constraint was non-finite
            if (constraint_error > max_initial_constraint_error) {{
                max_initial_constraint_error = constraint_error;
                worst_initial_constraint_ray = master_idx;
            }} // END IF: update maximum initial constraint residual
        }} // END LOOP: for constraint_i over initial chunk
        if (max_initial_constraint_error > initial_constraint_tolerance) {{
            fprintf(
                stderr,
                "ERROR: initial photon constraint residual %e exceeds tolerance %e "
                "at ray %ld.\n",
                max_initial_constraint_error,
                initial_constraint_tolerance,
                worst_initial_constraint_ray);
            time_window_manager_numerical_free(&numerical_window);
            slot_manager_free(&tsm);
            exit(1);
        }} // END IF: initial constraint residual exceeded tolerance

        for (int c_k = 0; c_k < 9; ++c_k) {{ // Loop index $c_k$ orchestrating memory transfer of the 9 constrained state vector $f^\mu$ components.
            {memcpy_cpu("f_bridge[0] + c_k * BUNDLE_CAPACITY", "d_f_bundle[0] + c_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size")}
        }} // END LOOP: for c_k over 9
        {no_sync()}

        long int nonfinite_count = 0; // Accumulator tracking the total number of physical states $f^\mu$ containing non-finite values post-constraint solving.
        for (int gather_i = 0; gather_i < chunk_size; ++gather_i) {{ // Loop iterator $gather_i$ scanning each trajectory within the retrieved initialization chunk.
            long int master_idx = start_idx + gather_i; // Computes the absolute master index $m_{{idx}}$ mapping the localized chunk to the global master SoA.
            bool has_nonfinite = false; // Boolean flag indicating if the specific state vector $f^\mu$ contains a non-finite value.
            for (int gather_k = 0; gather_k < 9; ++gather_k) {{ // Loop index $gather_k$ iterating over the 9 tensor components of the state vector $f^\mu$.
                double val = f_bridge[0][gather_k * BUNDLE_CAPACITY + gather_i]; // Evaluates the updated numerical value of the specific tensor component.
                all_photons_host.f[gather_k * num_rays + master_idx] = val; // Maps the valid constrained tensor scalar back to the global Host SoA.
                if (!isfinite(val)) has_nonfinite = true; // Flags the physical state vector as invalid due to a non-finite evaluation.
            }} // END LOOP: for gather_k over 9
            if (has_nonfinite) nonfinite_count++; // Increments the total count of unresolved physical state vectors $f^\mu$.
        }} // END LOOP: for gather_i over chunk_size

        if (nonfinite_count > 0) {{
            fprintf(stderr,
                    "ERROR: Init Batch %ld: %ld rays contain nonfinite state values "
                    "after tetrad initialization. Aborting numerical batch integration.\n",
                    init_batch,
                    nonfinite_count);
            time_window_manager_numerical_free(&numerical_window);
            slot_manager_free(&tsm);
            exit(1);
        }} // END IF: nonfinite_count > 0 to abort
    }} // END LOOP: for init_batch over num_batches

    {free_device}(d_initial_norm_bundle); // Purges the initialization constraint scratchpad.

    puncture_result_t observer_horizon = {{0, HORIZON_OUTSIDE, NAN, NAN}};
    if (photon_horizons != NULL) {{
        const double point[3] = {{commondata->observer_x, commondata->observer_y, commondata->observer_z}};
        observer_horizon = apparent_horizons_classify(&horizons, &numerical_window, commondata->t_start, point);
    }} // END IF: evaluate the shared observer position
    long int sync_i; // Loop iterator spanning the photon states.
    for(sync_i = 0; sync_i < num_rays; ++sync_i) {{
        all_photons_host.status[sync_i] = ACTIVE; // Assigns the initial trajectory activity enum for the global physics engine.
        all_photons_host.integration_param[sync_i] = {initial_integration_param}; // Sets the initial integration parameter.
        all_photons_host.rejection_retries[sync_i] = 0; // Clears the error rejection scalar to initialize the step size convergence tracking.
        all_photons_host.consecutive_min_steps[sync_i] = 0;
        if (photon_horizons != NULL) {{
            photon_horizons[sync_i] = observer_horizon;
            if (observer_horizon.horizon_status == HORIZON_INSIDE) {{
                all_photons_host.status[sync_i] = STOP_CONDITION_INSIDE_APPARENT_HORIZON;
                continue;
            }} // END IF: photon starts inside a horizon
        }} // END IF: initialize horizon result
        all_photons_host.h[sync_i] = rkf45_cap_numerical_step(
            commondata, all_photons_host.h[sync_i],
            all_photons_host.f[4 * num_rays + sync_i],
            &all_photons_host.status[sync_i]);
        if (all_photons_host.status[sync_i] != ACTIVE)
            continue;

        int s_idx = slot_get_index(&tsm, {initial_coordinate_time}); // Maps coordinate time to a TimeSlotManager bin.
        if (s_idx != -1) {{
            slot_add_photon(&tsm, s_idx, sync_i); // Registers the active photon index to its corresponding temporal bin mapped by the orchestrator.
        }} // END IF: s_idx != -1 to add
    }} // END LOOP: for sync_i over num_rays

    // Hardware clock state marking the beginning of the active integration chunk.
    struct timespec batch_start_time;
    // Hardware clock state marking the conclusion of the active integration chunk.
    struct timespec batch_end_time;
    // Counter tracking the number of finalized full RKF45 attempts since the previous dashboard print.
    long int attempted_rkf45_steps_since_print = 0;
    // User-tunable exponential moving average window for the RKF45 throughput display.
    const double rkf45_speed_avg_window_sec = 1.0;
    // Smoothed RKF45 attempt throughput displayed in the live dashboard.
    double rkf45_attempts_per_sec_avg = 0.0;
    // Tracks whether the throughput average has received its first valid sample.
    bool rkf45_speed_avg_initialized = false;
    // Captures the initial hardware clock state prior to temporal loop execution.
    clock_gettime(CLOCK_MONOTONIC, &batch_start_time);

    // Allocates vertical terminal space for the dynamic multi-line progress dashboard.
    printf("\n\n\n\n\n\n\n\n");

    //==========================================
    // 3. TEMPORAL LOOP (The Engine)
    //==========================================
    // Integer tracking the global number of active photon trajectories to allow early loop termination.
    long int total_active_photons = num_rays;

    // Outer loop iterator for the physical time bins.
    for (int slot_idx = tsm.num_slots - 1; slot_idx >= 0; --slot_idx) {{
        // Evaluates the early exit condition to terminate the temporal engine if all geometric trajectories have concluded.
        if (total_active_photons <= 0) {{
            break; // Terminates the temporal engine early to avoid unnecessary work.
        }} // END IF: total_active_photons <= 0 to break

        if (tsm.slot_counts[slot_idx] > 0) {{
            if (time_window_manager_numerical_mmap_for_slot(
                    &numerical_window, &tsm, {geometry_slot_idx}) !=
                TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
                fprintf(stderr,
                        "ERROR: failed to map numerical time window for slot %d [%e, %e).\n",
                        {geometry_slot_idx},
                        (double)slot_lower_time(&tsm, {geometry_slot_idx}),
                        (double)slot_upper_time(&tsm, {geometry_slot_idx}));
                time_window_manager_numerical_free(&numerical_window);
                slot_manager_free(&tsm);
                exit(1);
            }} // END IF: slot time window mapping failed
        }} // END IF: slot_idx still contained active photons

        int current = 0; // Integer index tracking the primary active CPU buffer for execution.
        int next = 1; // Buffer index for the next photon chunk.
        long int active_chunks[2] = {{ 0, 0}}; // 1D array storing the total number of trajectories queued for each operational buffer.

        //==========================================
        // PHASE A: PRIME THE PUMP (Buffer 0)
        //==========================================
        // Populate the first bridge and launch synchronous integration on the primary buffer.

        active_chunks[current] = NRPYMIN((long int)BUNDLE_CAPACITY, tsm.slot_counts[slot_idx]); // Evaluates the active chunk size bounding the CPU staging to avoid CPU scratch overflow.

        if (active_chunks[current] > 0) {{
            slot_remove_chunk(&tsm, slot_idx, chunk_buffer[current], active_chunks[current]); // Extracts the execution chunk mapping from the Host-side temporal bin.

            // 1. Pack the 9-component tensors using cache-friendly forward sweeps
            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ iterating over the 9 tensor components of the state vectors.
                for (int bridge_i = 0; bridge_i < active_chunks[current]; ++bridge_i) {{  // Copy each photon's state into the active chunk.
                    long int m_idx = chunk_buffer[current][bridge_i]; // Index in the full photon state array.
                    f_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f[c_k * num_rays + m_idx]; // Packs the coordinate state vector $f^\mu$ into the transfer bridge.
                }} // END LOOP: for bridge_i over active_chunks[current]
            }} // END LOOP: for c_k over 9

            // 2. Pack the 1D arrays in a separate sequential loop
            for (int bridge_i = 0; bridge_i < active_chunks[current]; ++bridge_i) {{
                long int m_idx = chunk_buffer[current][bridge_i];
                h_bridge[current][bridge_i] = all_photons_host.h[m_idx]; // Packs the current integration step size $h$ into the transfer bridge.
                status_bridge[current][bridge_i] = all_photons_host.status[m_idx]; // Packs the trajectory status enum into the transfer bridge.
                retries_bridge[current][bridge_i] = all_photons_host.rejection_retries[m_idx]; // Packs the error rejection scalar into the transfer bridge.
                min_steps_bridge[current][bridge_i] = all_photons_host.consecutive_min_steps[m_idx];
                integration_param_bridge[current][bridge_i] = all_photons_host.integration_param[m_idx]; // Packs the integration parameter into the transfer bridge.
                on_pos_non_terminal_plane_prev_bridge[current][bridge_i] = all_photons_host.on_positive_side_of_non_terminal_plane_prev[m_idx]; // Packs the nonterminal plane boundary flag into the transfer bridge.
                on_pos_terminal_plane_prev_bridge[current][bridge_i] = all_photons_host.on_positive_side_of_terminal_plane_prev[m_idx]; // Packs the terminal-plane boundary flag into the transfer bridge.
            }} // END LOOP: for bridge_i over active_chunks[current]

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating CPU buffer copy of the 9 state vector components.
                // CPU buffer copy: Synchronously pushes bounded state vectors $f^\mu$ to CPU scratch strictly on buffer [current] to minimize latency.
                {memcpy_cpu("d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "f_bridge[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]")}
            }} // END LOOP: for c_k over 9
            // CPU buffer copy: Synchronously pushes step sizes $h$ to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_h[current]", "h_bridge[current]", "sizeof(double) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes status enums to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_status[current]", "status_bridge[current]", "sizeof(termination_type_t) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes rejection scalars to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_retries[current]", "retries_bridge[current]", "sizeof(int) * active_chunks[current]")}
            {memcpy_cpu("d_min_steps[current]", "min_steps_bridge[current]", "sizeof(int) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes integration parameters to CPU scratch on buffer [current].
            {memcpy_cpu("d_integration_param_bundle[current]", "integration_param_bridge[current]", "sizeof(double) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes nonterminal-plane boundary flags to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_on_pos_non_terminal_plane_prev[current]", "on_pos_non_terminal_plane_prev_bridge[current]", "sizeof(bool) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes terminal-plane boundary flags to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_on_pos_terminal_plane_prev[current]", "on_pos_terminal_plane_prev_bridge[current]", "sizeof(bool) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes nonterminal-plane intersection locks to CPU scratch strictly on buffer [current] to minimize latency.
            // CPU buffer copy: Synchronously pushes terminal-plane intersection locks to CPU scratch strictly on buffer [current] to minimize latency.
            // CPU buffer copy: Synchronously pushes chunk indices $m_{{ idx}}$ to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_chunk_buffer[current]", "chunk_buffer[current]", "sizeof(long int) * active_chunks[current]")}

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating CPU buffer baseline setup of the 9 state vector components.
                // CPU buffer copy: Duplicates the initial physical state vector $f^\mu$ to anchor the final RKF45 evaluation.
                {memcpy_cpu("d_f_start_bundle[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]")}
                // Initialize the RKF45 stage state from the current photon state.
                {memcpy_cpu("d_f_temp_bundle[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]")}
            }} // END LOOP: for c_k over 9
{trial_spatial_center_setup_current}

            for (int stage = 1; stage <= 6; ++stage) {{  // Loop iterator $stage$ executing the 6 discrete stages of the RKF45 Runge-Kutta numerical solver.
                // Interpolate the metric and RHS geometry for the active chunk.
                numerical_interpolation(
                    commondata,
                    &numerical_params,
                    &spatial_context,
                    &numerical_window,
                    d_f_temp_bundle[current],
                    d_status[current],
                    {interpolation_integration_param_args}
                    d_metric_bundle[current],
                    d_rhs_geometry_bundle[current],
                    active_chunks[current],
                    current);
                // RHS step: compute the geodesic equation derivatives $\dot{{ f}}^\mu$
                // on the active buffer.
{gradient_update_current}
                calculate_ode_rhs_kernel(d_f_temp_bundle[current], d_metric_bundle[current], d_rhs_geometry_bundle[current], {rhs_integration_param_args} d_k_bundle[current], stage, active_chunks[current]{stream_arg_current});
                // Stage 6 still computes its RHS; only intermediate state update is skipped.
                if (stage < 6) {{
                    // Stage update: accumulate the intermediate RKF45 state updates on the
                    // active buffer.
                    rkf45_stage_update(d_f_start_bundle[current], d_k_bundle[current], d_h[current], stage, active_chunks[current], d_f_temp_bundle[current]{stream_arg_current});
                }} // END IF: skip stage-6 intermediate state update
            }} // END LOOP: for stage over 6

{finalize_current}
{accepted_mask_save_current}
{accepted_metric_refresh_current}
{log_energy_evaluation_current}
            attempted_rkf45_steps_since_print += active_chunks[current];
            // Record accepted steps that change sides of either plane.
            {event_manager_current}

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating CPU buffer copy of the 9 state vector components.
                // CPU buffer copy: Retrieves updated coordinate states $f^\mu$ back to CPU RAM synchronously on the active buffer.
                {memcpy_cpu("f_bridge[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]")}
            }} // END LOOP: for c_k over 9
            // CPU buffer copy: Retrieves active step sizes $h$ back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("h_bridge[current]", "d_h[current]", "sizeof(double) * active_chunks[current]")}
            // CPU buffer copy: Retrieves updated status enums back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("status_bridge[current]", "d_status[current]", "sizeof(termination_type_t) * active_chunks[current]")}
            // CPU buffer copy: Retrieves active rejection counts back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("retries_bridge[current]", "d_retries[current]", "sizeof(int) * active_chunks[current]")}
            {memcpy_cpu("min_steps_bridge[current]", "d_min_steps[current]", "sizeof(int) * active_chunks[current]")}
            // CPU buffer copy: Retrieves integration-parameter progress from the active buffer.
            {memcpy_cpu("integration_param_bridge[current]", "d_integration_param_bundle[current]", "sizeof(double) * active_chunks[current]")}
            // CPU buffer copy: Retrieves updated nonterminal-plane boundary flags back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("on_pos_non_terminal_plane_prev_bridge[current]", "d_on_pos_non_terminal_plane_prev[current]", "sizeof(bool) * active_chunks[current]")}
            // CPU buffer copy: Retrieves updated terminal-plane boundary flags back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("on_pos_terminal_plane_prev_bridge[current]", "d_on_pos_terminal_plane_prev[current]", "sizeof(bool) * active_chunks[current]")}
        }} // END IF: active chunks available

        //==========================================
        // PHASE B: THE OVERLAP LOOP
        //==========================================
        // Pack the next photon chunk while processing the active chunk.

        while (active_chunks[current] > 0 || tsm.slot_counts[slot_idx] > 0) {{

            active_chunks[next] = NRPYMIN((long int)BUNDLE_CAPACITY, tsm.slot_counts[slot_idx]); // Evaluates the active chunk size bounding the upcoming CPU staging execution block.
            if (active_chunks[next] > 0) {{
                slot_remove_chunk(&tsm, slot_idx, chunk_buffer[next], active_chunks[next]); // Extracts the next execution chunk mapping from the Host-side temporal bin.

                // 1. Pack the 9-component tensors using cache-friendly forward sweeps
                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ iterating over the 9 tensor components of the state vectors.
                    for (int bridge_i = 0; bridge_i < active_chunks[next]; ++bridge_i) {{  // Copy each photon's state into the next chunk.
                        long int m_idx = chunk_buffer[next][bridge_i]; // Index in the full photon state array.
                        f_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f[c_k * num_rays + m_idx]; // Packs the coordinate state vector $f^\mu$ into the transfer bridge.
                    }} // END LOOP: for bridge_i over active_chunks[next]
                }} // END LOOP: for c_k over 9

                // 2. Pack the 1D arrays in a separate sequential loop
                for (int bridge_i = 0; bridge_i < active_chunks[next]; ++bridge_i) {{
                    long int m_idx = chunk_buffer[next][bridge_i];
                    h_bridge[next][bridge_i] = all_photons_host.h[m_idx]; // Packs the current integration step size $h$ into the transfer bridge.
                    status_bridge[next][bridge_i] = all_photons_host.status[m_idx]; // Packs the trajectory status enum into the transfer bridge.
                    retries_bridge[next][bridge_i] = all_photons_host.rejection_retries[m_idx]; // Packs the error rejection scalar into the transfer bridge.
                    min_steps_bridge[next][bridge_i] = all_photons_host.consecutive_min_steps[m_idx];
                    integration_param_bridge[next][bridge_i] = all_photons_host.integration_param[m_idx]; // Packs the integration parameter into the transfer bridge.
                    on_pos_non_terminal_plane_prev_bridge[next][bridge_i] = all_photons_host.on_positive_side_of_non_terminal_plane_prev[m_idx]; // Packs the nonterminal plane boundary flag into the transfer bridge.
                    on_pos_terminal_plane_prev_bridge[next][bridge_i] = all_photons_host.on_positive_side_of_terminal_plane_prev[m_idx]; // Packs the terminal-plane boundary flag into the transfer bridge.
                }} // END LOOP: for bridge_i over active_chunks[next]

                for (int c_k = 0; c_k < 9; ++c_k) {{  // Copy the nine state components for the next chunk.
                    // CPU buffer copy: Synchronously pushes bounded state vectors $f^\mu$ to CPU scratch strictly on buffer [next] to overlap execution.
                    {memcpy_cpu("d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "f_bridge[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]")}
                }} // END LOOP: for c_k over 9
                // CPU buffer copy: Synchronously pushes step sizes $h$ to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_h[next]", "h_bridge[next]", "sizeof(double) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes status enums to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_status[next]", "status_bridge[next]", "sizeof(termination_type_t) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes rejection scalars to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_retries[next]", "retries_bridge[next]", "sizeof(int) * active_chunks[next]")}
                {memcpy_cpu("d_min_steps[next]", "min_steps_bridge[next]", "sizeof(int) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes integration parameters to CPU scratch on buffer [next].
                {memcpy_cpu("d_integration_param_bundle[next]", "integration_param_bridge[next]", "sizeof(double) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes nonterminal-plane boundary flags to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_on_pos_non_terminal_plane_prev[next]", "on_pos_non_terminal_plane_prev_bridge[next]", "sizeof(bool) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes terminal-plane boundary flags to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_on_pos_terminal_plane_prev[next]", "on_pos_terminal_plane_prev_bridge[next]", "sizeof(bool) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes nonterminal-plane intersection locks to CPU scratch strictly on buffer [next] to overlap execution.
                // CPU buffer copy: Synchronously pushes terminal-plane intersection locks to CPU scratch strictly on buffer [next] to overlap execution.
                // CPU buffer copy: Synchronously pushes chunk indices $m_{{ idx}}$ to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_chunk_buffer[next]", "chunk_buffer[next]", "sizeof(long int) * active_chunks[next]")}

                for (int c_k = 0; c_k < 9; ++c_k) {{  // Initialize the RKF45 start and stage states.
                    // CPU buffer copy: Duplicates the initial physical state vector $f^\mu$ to anchor the upcoming RKF45 evaluation.
                    {memcpy_cpu("d_f_start_bundle[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]")}
                    // Initialize the RKF45 stage state from the current photon state.
                    {memcpy_cpu("d_f_temp_bundle[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]")}
                }} // END LOOP: for c_k over 9
{trial_spatial_center_setup_next}

                for (int stage = 1; stage <= 6; ++stage) {{  // Loop iterator $stage$ executing the 6 discrete stages of the upcoming RKF45 Runge-Kutta numerical solver.
                    // Interpolate the metric and RHS geometry for the next chunk.
                    numerical_interpolation(
                        commondata,
                        &numerical_params,
                        &spatial_context,
                        &numerical_window,
                        d_f_temp_bundle[next],
                        d_status[next],
                        {interpolation_integration_param_args_next}
                        d_metric_bundle[next],
                        d_rhs_geometry_bundle[next],
                        active_chunks[next],
                        next);
                    // RHS step: compute the geodesic equation derivatives
                    // $\dot{{ f}}^\mu$ on the alternate buffer.
{gradient_update_next}
                    calculate_ode_rhs_kernel(d_f_temp_bundle[next], d_metric_bundle[next], d_rhs_geometry_bundle[next], {rhs_integration_param_args_next} d_k_bundle[next], stage, active_chunks[next]{stream_arg_next});
                    // Stage 6 still computes its RHS; only intermediate state update is skipped.
                    if (stage < 6) {{
                        // Stage update: accumulate the intermediate RKF45 state updates on
                        // the alternate buffer.
                        rkf45_stage_update(d_f_start_bundle[next], d_k_bundle[next], d_h[next], stage, active_chunks[next], d_f_temp_bundle[next]{stream_arg_next});
                    }} // END IF: skip stage-6 intermediate state update
                }} // END LOOP: for stage over 6

{finalize_next}
{accepted_mask_save_next}
{accepted_metric_refresh_next}
{log_energy_evaluation_next}
                attempted_rkf45_steps_since_print += active_chunks[next];
                // Record accepted steps that change sides of either plane.
                {event_manager_next}
                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating CPU buffer copy of the 9 upcoming state vector components.
                    // CPU buffer copy: Retrieves updated coordinate states $f^\mu$ back to CPU RAM synchronously on the alternate buffer.
                    {memcpy_cpu("f_bridge[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]")}
                }} // END LOOP: for c_k over 9
                // CPU buffer copy: Retrieves upcoming active step sizes $h$ back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("h_bridge[next]", "d_h[next]", "sizeof(double) * active_chunks[next]")}
                // CPU buffer copy: Retrieves upcoming updated status enums back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("status_bridge[next]", "d_status[next]", "sizeof(termination_type_t) * active_chunks[next]")}
                // CPU buffer copy: Retrieves upcoming active rejection counts back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("retries_bridge[next]", "d_retries[next]", "sizeof(int) * active_chunks[next]")}
                {memcpy_cpu("min_steps_bridge[next]", "d_min_steps[next]", "sizeof(int) * active_chunks[next]")}
                // CPU buffer copy: Retrieves integration-parameter progress from the alternate buffer.
                {memcpy_cpu("integration_param_bridge[next]", "d_integration_param_bundle[next]", "sizeof(double) * active_chunks[next]")}
                // CPU buffer copy: Retrieves upcoming updated nonterminal-plane boundary flags back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("on_pos_non_terminal_plane_prev_bridge[next]", "d_on_pos_non_terminal_plane_prev[next]", "sizeof(bool) * active_chunks[next]")}
                // CPU buffer copy: Retrieves upcoming updated terminal-plane boundary flags back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("on_pos_terminal_plane_prev_bridge[next]", "d_on_pos_terminal_plane_prev[next]", "sizeof(bool) * active_chunks[next]")}
            }} // END IF: next chunks available

            if (active_chunks[current] > 0) {{
                // Finish the chunk before copying its states back to the full array.
                {no_sync()}

                // 1. Unpack 9-component tensors sequentially
                for (int fin_k = 0; fin_k < 9; ++fin_k) {{  // Loop index $fin_k$ retrieving the 9 tensor components back into the global Host matrix.
                    for (int fin_i = 0; fin_i < active_chunks[current]; ++fin_i) {{  // Loop iterator $fin_i$ unpacking the finalized physical data back to the global Host matrix.
                        long int m_idx = chunk_buffer[current][fin_i]; // Absolute master index $m_{{ idx}}$ retrieving the specific photon index from the execution chunk.
                        all_photons_host.f[fin_k * num_rays + m_idx] = f_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i]; // Unpacks the synchronized state vector $f^\mu$ into the global Host matrix.
                    }} // END LOOP: for fin_i over active_chunks[current]
                }} // END LOOP: for fin_k over 9

                // 2. Unpack 1D arrays sequentially
                for (int fin_i = 0; fin_i < active_chunks[current]; ++fin_i) {{
                    long int m_idx = chunk_buffer[current][fin_i];
                    all_photons_host.h[m_idx] = h_bridge[current][fin_i]; // Unpacks the synchronized step size $h$ into the global Host matrix.
                    all_photons_host.status[m_idx] = status_bridge[current][fin_i]; // Unpacks the synchronized trajectory status into the global Host matrix.
                    all_photons_host.rejection_retries[m_idx] = retries_bridge[current][fin_i]; // Unpacks the synchronized rejection count into the global Host matrix.
                    all_photons_host.consecutive_min_steps[m_idx] = min_steps_bridge[current][fin_i];
                    all_photons_host.integration_param[m_idx] = integration_param_bridge[current][fin_i]; // Unpacks the synchronized integration parameter into the global Host matrix.
                    all_photons_host.on_positive_side_of_non_terminal_plane_prev[m_idx] = on_pos_non_terminal_plane_prev_bridge[current][fin_i]; // Unpacks the synchronized nonterminal-plane boundary flag into the global Host matrix.
                    all_photons_host.on_positive_side_of_terminal_plane_prev[m_idx] = on_pos_terminal_plane_prev_bridge[current][fin_i]; // Unpacks the synchronized terminal-plane boundary flag into the global Host matrix.

                }} // END LOOP: for fin_i over active_chunks[current]

                // 3. TimeSlotManager State Update (Cache-hot, strictly sequential)
                for (int fin_i = 0; fin_i < active_chunks[current]; ++fin_i) {{
                    long int m_idx = chunk_buffer[current][fin_i];
                    // Preserve existing physical stops and time-domain errors.
                    const int min_step_limit = commondata->rkf45_max_consecutive_min_steps;
                    if (min_step_limit > 0 && min_steps_bridge[current][fin_i] >= min_step_limit &&
                        (status_bridge[current][fin_i] == REJECTED ||
                         (status_bridge[current][fin_i] == ACTIVE &&
                          slot_get_index(&tsm, {updated_coordinate_time}) != -1))) {{
                        status_bridge[current][fin_i] = FAILURE_RKF45_MIN_STEP_LIMIT;
                        all_photons_host.status[m_idx] = FAILURE_RKF45_MIN_STEP_LIMIT;
                    }} // END IF: consecutive minimum-step limit reached
                    if (status_bridge[current][fin_i] == ACTIVE) {{  // Evaluates the continuation logic if the trajectory remains within safe physical bounds.
                        int next_s_idx = slot_get_index(&tsm, {updated_coordinate_time}); // Maps updated coordinate time to the next operational bin.
                        if (next_s_idx != -1) {{  // Confirms the physical state has not exceeded the maximum simulation time bounds.
                            slot_add_photon(&tsm, next_s_idx, m_idx);  // Re-queues the updated physical state vector back into the host orchestrator.
                        }} else {{
                            all_photons_host.status[m_idx] = STOP_CONDITION_T_MAX_EXCEEDED; // Flags the physical state as stopped after excessive propagation time.
                            total_active_photons--; // Decrements the global counter as the physical trajectory has reached a terminal state.
                        }} // END ELSE: state flagged failed
                    }} // END IF: trajectory remains active
                    else if (status_bridge[current][fin_i] == REJECTED) {{   // RK error or interpolation failure remains retryable.
                        all_photons_host.status[m_idx] = ACTIVE; // Requeue status is active; retry count is preserved.
                        slot_add_photon(&tsm, slot_idx, m_idx); // Re-adds to the current bin to attempt integration with an adapted step-size scalar $h$.
                    }} else {{
                        total_active_photons--; // Decrements the global counter as the physical trajectory has reached a terminal state.
                    }} // END ELSE: trajectory reached terminal state
                }} // END LOOP: for fin_i over active_chunks[current]
            }} // END IF: active chunks completed

            //==========================================
            // PROGRESS DASHBOARD
            //==========================================
            // Evaluate computational throughput and temporal progress to update the terminal dashboard.

            // Captures the terminal hardware clock state for the current chunk.
            clock_gettime(CLOCK_MONOTONIC, &batch_end_time);

            // Evaluates the absolute wall-clock duration of the integration chunk in seconds.
            double elapsed_sec = (batch_end_time.tv_sec - batch_start_time.tv_sec) + (batch_end_time.tv_nsec - batch_start_time.tv_nsec) / 1e9;

            // Evaluates the finalized full-RKF45-attempt throughput across the current dashboard interval.
            double rkf45_attempts_per_sec =
                (elapsed_sec > 0.0) ? ((double)attempted_rkf45_steps_since_print / elapsed_sec) : 0.0;

            // Updates the smoothed RKF45 attempt throughput using an exponential moving average.
            if (elapsed_sec > 0.0 && attempted_rkf45_steps_since_print > 0) {{
                const double alpha = elapsed_sec / (rkf45_speed_avg_window_sec + elapsed_sec);

                if (!rkf45_speed_avg_initialized) {{
                    rkf45_attempts_per_sec_avg = rkf45_attempts_per_sec;
                    rkf45_speed_avg_initialized = true;
                }} else {{
                    rkf45_attempts_per_sec_avg =
                        (1.0 - alpha) * rkf45_attempts_per_sec_avg + alpha * rkf45_attempts_per_sec;
                }} // END ELSE: update the existing RKF45 throughput
            }} // END IF: elapsed interval had attempts

            // Evaluates the global completion ratio bounded between $0.0$ and $1.0$.
            double percent_done = 100.0 * (1.0 - ((double)total_active_photons / (double)num_rays));

            // Evaluates the physical coordinate time $t$ for the active temporal bin.
            double current_t = commondata->slot_manager_t_min + slot_idx * commondata->slot_manager_delta_t;

            // Defines the total character width of the dynamic loading bar visualization.
            int bar_width = 20;
            // Computes the integer index demarcating the active boundary within the loading bar.
            int pos = (int)(bar_width * percent_done / 100.0);
            // Character array storing the formatted loading bar string.
            char bar[21];

            // Loop iterator $bar_i$ constructing the ASCII loading bar visualizer.
            for (int bar_i = 0; bar_i < bar_width; ++bar_i) {{
                if (bar_i < pos) bar[bar_i] = '='; // Appends the completed progression character.
                else if (bar_i == pos) bar[bar_i] = '>'; // Appends the active vanguard character.
                else bar[bar_i] = ' '; // Appends the uncompleted progression character.
            }} // END LOOP: for bar_i over bar_width
            bar[bar_width] = '\0'; // Terminates the loading bar character array to prevent buffer overruns.

            // Count retryable rejections in this chunk, not consecutive retry
            // depths. Terminal failures and accepted-step caps are excluded.
            long int batch_rejections = 0;
            long int rkf45_rejections = 0;
            long int spatial_rejections = 0;
            long int temporal_rejections = 0;

            // Saved incoming statuses distinguish interpolation failures from
            // RKF45 error-control rejections after finalization.
            for (int sum_i = 0; sum_i < active_chunks[current]; ++sum_i) {{
                if (status_bridge[current][sum_i] == REJECTED) {{
                    batch_rejections++;
                    if (trial_status[current][sum_i] == FAILURE_SPATIAL_INTERPOLATION) {{
                        spatial_rejections++;
                    }} else if (trial_status[current][sum_i] == FAILURE_TEMPORAL_INTERPOLATION) {{
                        temporal_rejections++;
                    }} else {{
                        rkf45_rejections++;
                    }} // END ELSE: RKF45 error control rejected
                }} // END IF: rejected finalized status
            }} // END LOOP: count chunk rejection reasons

            // Every percentage uses all attempts in this chunk. Independent
            // rounding can make the displayed percentages differ by 0.1%.
            const double reject_percent = (active_chunks[current] > 0) ? (100.0 * (double)batch_rejections / (double)active_chunks[current]) : 0.0;
            const double rkf45_reject_percent = (active_chunks[current] > 0) ? (100.0 * (double)rkf45_rejections / (double)active_chunks[current]) : 0.0;
            const double spatial_reject_percent = (active_chunks[current] > 0) ? (100.0 * (double)spatial_rejections / (double)active_chunks[current]) : 0.0;
            const double temporal_reject_percent = (active_chunks[current] > 0) ? (100.0 * (double)temporal_rejections / (double)active_chunks[current]) : 0.0;

            // Restores the terminal cursor position vertically to overwrite the previous dashboard iteration.
            printf("\033[8A");
            printf("--------------------------------------------------\n"); // Prints the upper border of the diagnostic dashboard.
            printf(" Progress:   [%s] %5.1f%% \033[K\n", bar, percent_done); // Prints the global completion loading bar and percentage.
            printf(" Active:     %ld / %ld \033[K\n", total_active_photons, num_rays); // Prints the remaining active photon trajectories $x^\mu$.
            printf(" Slot Time:  Slot %d (t = %.1f) \033[K\n", slot_idx, current_t); // Prints the current physical temporal bin coordinate $t$.
            printf(" Speed:      %.2e RKF45 attempts/s, %.0fs avg \033[K\n", rkf45_attempts_per_sec_avg, rkf45_speed_avg_window_sec); // Print the mean RKF45 attempt rate.
            printf(" Rejects:    %ld (%.1f%%) \033[K\n", batch_rejections, reject_percent); // Prints the retryable rejection frequency.
            printf("             RKF45 %.1f%% | Spatial %.1f%% | Temporal %.1f%% \033[K\n", rkf45_reject_percent, spatial_reject_percent, temporal_reject_percent);
            printf("--------------------------------------------------\n"); // Prints the lower border of the diagnostic dashboard.
            fflush(stdout); // Flushes the standard output buffer to ensure instantaneous terminal rendering.

            // Resets the hardware clock state for the upcoming integration chunk.
            clock_gettime(CLOCK_MONOTONIC, &batch_start_time);
            attempted_rkf45_steps_since_print = 0;
            // --------------------------

            active_chunks[current] = 0; // Clears the execution queue tracker to indicate the active chunk has been fully processed.
            int temp = current; // Temporary integer scalar storing the primary buffer index for logical pointer swapping.
            current = next; // Shifts the primary execution tracker to the alternate buffer index.
            next = temp; // Reuse the cleared buffer for the next chunk.
        }} // END WHILE: alternate buffers process temporal slot

        //==========================================
        // PHASE C: THE TIME BARRIER
        //==========================================
        // CPU execution is synchronous here; no additional extra barrier is needed.

     }} // END LOOP: for slot_idx down to 0

    //==========================================
    // 4. RESULT WRITING, CLEANUP, AND FINALIZATION
    //==========================================
    // Process terminal photon trajectories, evaluate optional diagnostics, and extract final geometric intersections.

        printf("Refining recorded plane crossings with Brent's method...\n");
        fflush(stdout);
        find_plane_crossings_brent(
            commondata, &numerical_params, &spatial_context, &numerical_window,
            slot_manager_t_max, num_rays, plane_crossing_candidates,
            non_terminal_crossing_results, terminal_crossing_results);
        printf("Plane crossing refinement complete.\n");
        fflush(stdout);

        // Record selected plane positions before final blueprint calculation.
        for (long int ray = 0; ray < num_rays; ++ray) {{
            non_terminal_plane_crossing_result_t *non_terminal =
                &non_terminal_crossing_results[ray];
            terminal_plane_crossing_result_t *terminal =
                &terminal_crossing_results[ray];
            if (non_terminal->found) {{
                results_buffer[ray].y_nt = non_terminal->y_local;
                results_buffer[ray].z_nt = non_terminal->z_local;
                results_buffer[ray].non_terminal_plane_lambda =
                    non_terminal->physical_lambda;
                results_buffer[ray].non_terminal_plane_t =
                    non_terminal->coordinate_time;
                results_buffer[ray].non_terminal_plane_signed_distance =
                    non_terminal->plane_distance;
            }} // END IF: nonterminal crossing found
            if (terminal->found) {{
                terminal->trajectory_stop_reason = all_photons_host.status[ray];
                all_photons_host.status[ray] = STOP_CONDITION_TERMINAL_PLANE;
                results_buffer[ray].y_t = terminal->y_local;
                results_buffer[ray].z_t = terminal->z_local;
                results_buffer[ray].L_f = terminal->physical_lambda;
                results_buffer[ray].t_f = terminal->coordinate_time;
            }} else if (all_photons_host.status[ray] == STOP_CONDITION_INSIDE_APPARENT_HORIZON) {{
                results_buffer[ray].y_t = photon_horizons[ray].puncture_id;
            }} // END IF: record terminal geometry
        }} // END LOOP: copy refined plane crossings

{terminal_diagnostics_block}

{non_terminal_normalization_block}

        write_plane_crossings(non_terminal_crossings_path, num_rays,
            non_terminal_crossing_results, terminal_crossing_results, false);
        write_plane_crossings(terminal_crossings_path, num_rays,
            non_terminal_crossing_results, terminal_crossing_results, true);

        // Final output step: process escaped photons intersecting the celestial
        // sphere $r > r_{{escape}}$ and optionally write the direct-scale-equivalent
        // normalization-error sidecar.
        if (commondata->perform_normalization_check &&
            normalization_abs_by_ray != NULL &&
            norm_abs_bin_path != NULL &&
            norm_abs_bin_path[0] != '\0') {{
            {calc_blueprint_with_norm_abs}
        }} else {{
            {calc_blueprint}
        }} // END ELSE: blueprint post-processing without normalization sidecar
{angular_momentum_output}
{synthetic_slice_usage_output}

        write_puncture_diagnostics(&puncture_data, photon_horizons, commondata, &all_photons_host,
                                  results_buffer, num_rays, puncture_bin_path,
                                  {'true' if fast_light else 'false'});
        free(photon_horizons);
        apparent_horizons_close(&horizons);
        puncture_free_diagnostics(&puncture_data);

        if (commondata->perform_normalization_check) {{
            {free_device}(d_norm_bundle); // Purges the terminal normalization diagnostic scratchpad.
            {free_pinned}(normalization_abs_by_ray); // Purges the per-photon normalization sidecar buffer.
            {free_pinned}(normalization_abs_non_terminal_by_ray); // Purges the sparse nonterminal normalization values.
            {free_pinned}(non_terminal_plane_norm_recorded_by_ray);
        }} // END IF: commondata->perform_normalization_check to purge normalization scratchpad
        if (commondata->perform_synthetic_slice_check) {{
            {free_pinned}(synthetic_slice_usage_by_photon);
        }} // END IF: free synthetic temporal-stencil usage records
{angular_momentum_cleanup}

        // Loop iterator $s$ purging the double-buffered arrays across both CPU buffers.
        for (int s = 0; s < 2; ++s) {{
            // Host Memory Free: Purges bridge components supporting scatter logic mapped to CPU memory transfers.
            {free_pinned}(chunk_buffer[s]); // Purges the execution chunk mapping bridge.
            {free_pinned}(f_bridge[s]); // Purges the state vector $f^\mu$ bridge.
            {free_pinned}(integration_param_bridge[s]); // Purges the integration-parameter bridge.
            {free_pinned}(h_bridge[s]); // Purges the integration step size $h$ bridge.
            {free_pinned}(status_bridge[s]); // Purges the trajectory status bridge.
            {free_pinned}(trial_status[s]); // Release the saved trial statuses.
            {free_pinned}(retries_bridge[s]); // Purges the error rejection scalar bridge.
            {free_pinned}(min_steps_bridge[s]);
            {free_pinned}(on_pos_non_terminal_plane_prev_bridge[s]); // Purges the nonterminal plane boundary flag bridge.
            {free_pinned}(on_pos_terminal_plane_prev_bridge[s]); // Purges the terminal-plane boundary flag bridge.

            // Free the RKF45 work arrays for this buffer.
            {free_device}(d_f_bundle[s]); // Purges the state vector $f^\mu$ scratchpad.
            {free_device}(d_f_start_bundle[s]); // Purges the anchor state vector $f_{{start}}$ scratchpad.
            {free_device}(d_f_temp_bundle[s]); // Purges the temporary stage $f^\mu_{{temp}}$ scratchpad.
            {free_device}(d_metric_bundle[s]); // Purges the symmetric metric tensor $g_{{\mu\nu}}$ scratchpad.
            {free_device}(d_log_energy_bundle[s]); // Purges the common log-energy scratchpad.
            {free_device}(d_rhs_geometry_bundle[s]); // Purges the geometry scratchpad.
{gradient_cleanup}
            {free_device}(d_k_bundle[s]); // Purges the derivative tensor $\dot{{f}}^\mu$ scratchpad.
            {free_device}(d_h[s]); // Purges the active integration step sizing $h$ scratchpad.
            {free_device}(d_integration_param_bundle[s]); // Purges the integration-parameter scratchpad.
            {free_device}(d_status[s]); // Purges the current trajectory status limit scratchpad.
            {free_device}(d_retries[s]); // Purges the sequential error rejection scratchpad.
            {free_device}(d_min_steps[s]);
            {free_device}(d_on_pos_non_terminal_plane_prev[s]); // Purges the previous nonterminal plane boundary side scratchpad.
            {free_device}(d_on_pos_terminal_plane_prev[s]); // Purges the previous terminal-plane boundary side scratchpad.
            {free_device}(d_chunk_buffer[s]); // Purges the absolute master indices $m_{{idx}}$ mapping scratchpad.
{spatial_center_frees}
        }} // END LOOP: for s over 2


        // Host Memory Free: Purges the primary Host state and integration-parameter arrays.
        {free_pinned}(all_photons_host.f); // Purges the primary Host array state $f^\mu$.
        {free_pinned}(all_photons_host.integration_param); // Purges the primary Host integration-parameter array.
        {free_pinned}(all_photons_host.h); // Purges the primary Host array integration step size $h$.
        {free_pinned}(all_photons_host.status); // Purges the primary Host array trajectory status enum.
        {free_pinned}(all_photons_host.rejection_retries); // Purges the primary Host array error rejection scalar.
        {free_pinned}(all_photons_host.consecutive_min_steps);
        {free_pinned}(all_photons_host.on_positive_side_of_non_terminal_plane_prev); // Purges the primary Host array nonterminal plane boundary flag.
        {free_pinned}(all_photons_host.on_positive_side_of_terminal_plane_prev); // Purges the primary Host array terminal-plane boundary flag.
        {free_pinned}(plane_crossing_candidates);
        {free_pinned}(non_terminal_crossing_results);
        {free_pinned}(terminal_crossing_results);
        {free_pinned}(accepted_mask_selection);

        // Release the active numerical spacetime window before the slot lattice is destroyed.
        time_window_manager_numerical_free(&numerical_window);

        // Memory Free: Purges the temporal sorting struct mapping the Host-side execution grid.
        slot_manager_free(&tsm); // Purges the central time slot orchestrator.
    """

    cfc.register_CFunction(
        includes=includes,
        prefunc=plane_crossing_writer + puncture_helpers + mask_selection_helpers,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
