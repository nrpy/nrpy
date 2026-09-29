# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/event_detection_manager_kernel.py
"""
Generate delayed event detection for photon plane crossings.

After each accepted step, the generated kernel detects plane crossings and
retains them until enough accepted states exist for centered interpolation.
It reconstructs pending crossings before a physical stop and shifts accepted
state history for the next RKF45 step.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par


def event_detection_manager_kernel(
    maximum_degree: int,
    normalized_eom: bool = False,
    capture_event_state: bool = False,
    single_photon_step_limit: bool = False,
    numerical_time_limit: bool = False,
) -> None:
    """
    Register delayed plane detection with centered state interpolation.

    :param normalized_eom: Whether component zero stores affine parameter.
    :param capture_event_state: Whether to copy interpolated states to arrays.
    :param maximum_degree: Largest generated interpolation polynomial degree.
    :param single_photon_step_limit: Accept a final-step flag so the single photon
        can interpolate a pending crossing before its accepted-step limit.
    :param numerical_time_limit: Stop numerical batch photons at the coordinate-time
        slot boundary before processing pending plane crossings.
    :raises ValueError: If maximum_degree is below three.
    """
    if maximum_degree < 3:
        raise ValueError("Plane interpolation degree must be at least three")
    par.register_CodeParameters(
        "bool",
        __name__,
        ["non_terminal_plane_enabled", "terminal_plane_enabled"],
        [False, False],
        commondata=True,
        add_to_parfile=True,
        descriptions=[
            "Enable nonterminal-plane crossing detection.",
            "Enable terminal-plane crossing detection.",
        ],
    )
    parallelization = par.parval_from_str("parallelization")
    cd_access = parallel_utils.get_commondata_access(parallelization)
    commondata_arg = "" if parallelization == "cuda" else ", commondata"
    escape_statement = "return;" if parallelization == "cuda" else "continue;"
    required_post_steps = (maximum_degree + 1) // 2

    args = {
        "d_f_bundle": "const double *restrict",
        "d_log_energy_bundle": "const double *restrict",
        "d_f_history_bundle": "double *restrict",
        "d_integration_param": "const double *restrict",
        "d_integration_param_history": "double *restrict",
        "d_results_buffer": "blueprint_data_t *restrict",
        "d_status_bundle": "termination_type_t *restrict",
        "d_on_pos_non_terminal_plane_prev": "bool *restrict",
        "d_on_pos_terminal_plane_prev": "bool *restrict",
        "d_non_terminal_plane_event_found": "bool *restrict",
        "d_terminal_plane_event_found": "bool *restrict",
        "d_non_terminal_plane_crossing_pending": "bool *restrict",
        "d_terminal_plane_crossing_pending": "bool *restrict",
        "d_non_terminal_plane_steps_past": "int *restrict",
        "d_terminal_plane_steps_past": "int *restrict",
        "d_non_terminal_plane_event_degree": "int *restrict",
        "d_terminal_plane_event_degree": "int *restrict",
    }
    event_state_args = ""
    if capture_event_state:
        event_state_args = (
            "double *restrict d_non_terminal_plane_event_state_bundle, "
            "double *restrict d_terminal_plane_event_state_bundle, "
        )
        args["d_non_terminal_plane_event_state_bundle"] = "double *restrict"
        args["d_terminal_plane_event_state_bundle"] = "double *restrict"
    if single_photon_step_limit:
        args["force_pending_interpolation"] = "const bool"
    args["d_chunk_buffer"] = "const long int *restrict"
    args["chunk_size"] = "const int"
    force_arg = (
        "const bool force_pending_interpolation, " if single_photon_step_limit else ""
    )
    if parallelization != "cuda":
        args["commondata"] = "const commondata_struct *restrict"

    non_terminal_capture = ""
    terminal_capture = ""
    if capture_event_state:
        non_terminal_capture = """
                for (int component = 0; component < 9; ++component) {
                    WriteCUDA(&d_non_terminal_plane_event_state_bundle[
                        IDX_F(component, i)], f_int[component]);
                }
"""
        terminal_capture = """
                for (int component = 0; component < 9; ++component) {
                    WriteCUDA(&d_terminal_plane_event_state_bundle[
                        IDX_F(component, i)], f_int[component]);
                }
"""
    intersection_setup = (
        """
            double f_intersection[9];
            for (int component = 0; component < 9; ++component) {
                f_intersection[component] = f_int[component];
            }
            f_intersection[0] = event_integration_param;
            const double physical_lambda = f_int[0];
"""
        if normalized_eom
        else "const double physical_lambda = event_integration_param;"
    )
    intersection_state = "f_intersection" if normalized_eom else "f_int"
    time_limit_check = ""
    if numerical_time_limit:
        coordinate_time = (
            "integration_params[0]" if normalized_eom else "state_history[0]"
        )
        time_limit_check = f"""
    const double coordinate_time = {coordinate_time};
    if (d_status_bundle[i] == ACTIVE &&
        (!isfinite(coordinate_time) ||
         coordinate_time < {cd_access}slot_manager_t_min ||
         coordinate_time >= {cd_access}t_start + 1.0e-5)) {{
        d_status_bundle[i] = STOP_CONDITION_T_MAX_EXCEEDED;
    }}
"""
    if parallelization == "cuda":
        loop_start = """
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= chunk_size) return;
"""
        loop_end = ""
    else:
        loop_start = """
    #pragma omp parallel for
    for (long int i = 0; i < chunk_size; ++i) {
"""
        loop_end = "    } // END LOOP: accepted photon states"

    core = r"""
    #define IDX_F(c, ray) ((c) * BUNDLE_CAPACITY + (ray))
    #define IDX_HISTORY(node, c, ray) (((node) * 9 + (c)) * BUNDLE_CAPACITY + (ray))
    #define IDX_PARAM(node, ray) ((node) * BUNDLE_CAPACITY + (ray))

    // RK failures and rejected steps have no new accepted state.
    if (d_status_bundle[i] != ACTIVE) @ESCAPE@
    const long int master_idx = d_chunk_buffer[i];
    double state_history[@HISTORY_SIZE@ * 9];
    double integration_params[@HISTORY_SIZE@];
    for (int component = 0; component < 9; ++component) {
        state_history[component] = ReadCUDA(&d_f_bundle[IDX_F(component, i)]);
    }
    integration_params[0] = ReadCUDA(&d_integration_param[i]);
    for (int node = 1; node < @HISTORY_SIZE@; ++node) {
        for (int component = 0; component < 9; ++component) {
            state_history[9 * node + component] = ReadCUDA(
                &d_f_history_bundle[IDX_HISTORY(node - 1, component, i)]);
        }
        integration_params[node] = ReadCUDA(
            &d_integration_param_history[IDX_PARAM(node - 1, i)]);
    }
    const double x = state_history[1];
    const double y = state_history[2];
    const double z = state_history[3];

    // Detect both crossings before physical stop checks. A crossing and a
    // physical stop on this accepted step must still reconstruct the event.
    double w_normal[3] = {@CD@non_terminal_plane_normal_x,
                          @CD@non_terminal_plane_normal_y,
                          @CD@non_terminal_plane_normal_z};
    double w_dist = 0.0;
    if (@CD@non_terminal_plane_enabled &&
        !d_non_terminal_plane_event_found[i]) {
        const double w_sq = w_normal[0] * w_normal[0] +
                            w_normal[1] * w_normal[1] +
                            w_normal[2] * w_normal[2];
        if (!isfinite(w_sq) || w_sq <= 1.0e-28) {
            d_status_bundle[i] = FAILURE_GENERIC;
            @ESCAPE@
        }
        const double inverse = 1.0 / SqrtCUDA(w_sq);
        for (int axis = 0; axis < 3; ++axis) w_normal[axis] *= inverse;
        w_dist = @CD@non_terminal_plane_center_x * w_normal[0] +
                 @CD@non_terminal_plane_center_y * w_normal[1] +
                 @CD@non_terminal_plane_center_z * w_normal[2];
        const bool side = x * w_normal[0] + y * w_normal[1] +
                          z * w_normal[2] - w_dist > 0.0;
        if (d_non_terminal_plane_crossing_pending[i]) {
            ++d_non_terminal_plane_steps_past[i];
        } else if (side != d_on_pos_non_terminal_plane_prev[i]) {
            d_non_terminal_plane_crossing_pending[i] = true;
            d_non_terminal_plane_steps_past[i] = 1;
        }
        d_on_pos_non_terminal_plane_prev[i] = side;
    }

    double s_normal[3] = {@CD@terminal_plane_normal_x,
                          @CD@terminal_plane_normal_y,
                          @CD@terminal_plane_normal_z};
    double s_dist = 0.0;
    if (@CD@terminal_plane_enabled && !d_terminal_plane_event_found[i]) {
        const double s_sq = s_normal[0] * s_normal[0] +
                            s_normal[1] * s_normal[1] +
                            s_normal[2] * s_normal[2];
        if (!isfinite(s_sq) || s_sq <= 1.0e-28) {
            d_status_bundle[i] = FAILURE_GENERIC;
            @ESCAPE@
        }
        const double inverse = 1.0 / SqrtCUDA(s_sq);
        for (int axis = 0; axis < 3; ++axis) s_normal[axis] *= inverse;
        s_dist = @CD@terminal_plane_center_x * s_normal[0] +
                 @CD@terminal_plane_center_y * s_normal[1] +
                 @CD@terminal_plane_center_z * s_normal[2];
        const bool side = x * s_normal[0] + y * s_normal[1] +
                          z * s_normal[2] - s_dist > 0.0;
        if (d_terminal_plane_crossing_pending[i]) {
            ++d_terminal_plane_steps_past[i];
        } else if (side != d_on_pos_terminal_plane_prev[i]) {
            d_terminal_plane_crossing_pending[i] = true;
            d_terminal_plane_steps_past[i] = 1;
        }
        d_on_pos_terminal_plane_prev[i] = side;
    }

    const double log_energy_measure = ReadCUDA(&d_log_energy_bundle[i]);
    if (log_energy_measure > @CD@evolution_measure_max) {
        d_status_bundle[i] = STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED;
    }
    const double radius_squared = x * x + y * y + z * z;
    if (radius_squared > @CD@r_escape * @CD@r_escape &&
        d_status_bundle[i] == ACTIVE) {
        d_status_bundle[i] = STOP_CONDITION_COORD_RADIUS_EXCEEDED;
    }
@TIME_LIMIT_CHECK@
    if (d_status_bundle[i] != ACTIVE &&
        !d_non_terminal_plane_crossing_pending[i] &&
        !d_terminal_plane_crossing_pending[i]) @ESCAPE@

    // Process terminal plane first. If it stops the photon, reconstruct any
    // earlier pending nonterminal crossing using available accepted states.
    if (d_terminal_plane_crossing_pending[i] &&
        (d_terminal_plane_steps_past[i] >= @POST_STEPS@ ||
         d_status_bundle[i] != ACTIVE || @FORCE_PENDING@)) {
        double f_int[9];
        double event_integration_param;
        int actual_degree;
        const int result = find_event_time_and_state_centered(
            state_history, integration_params,
            d_terminal_plane_steps_past[i], s_normal, s_dist,
            &event_integration_param, f_int, &actual_degree);
        if (result != 1) {
            d_status_bundle[i] = result == 0
                ? FAILURE_PLANE_INTERPOLATION_HISTORY : FAILURE_GENERIC;
            @ESCAPE@
        }
        @INTERSECTION_SETUP@
        if (handle_terminal_plane_intersection(
                @INTERSECTION_STATE@, physical_lambda,
                &d_results_buffer[master_idx]@COMMONDATA_ARG@)) {
            d_status_bundle[i] = STOP_CONDITION_TERMINAL_PLANE;
            d_terminal_plane_event_found[i] = true;
            d_terminal_plane_event_degree[i] = actual_degree;
            @TERMINAL_CAPTURE@
        }
        d_terminal_plane_crossing_pending[i] = false;
        d_terminal_plane_steps_past[i] = 0;
    }

    if (d_non_terminal_plane_crossing_pending[i] &&
        (d_non_terminal_plane_steps_past[i] >= @POST_STEPS@ ||
         d_status_bundle[i] != ACTIVE || @FORCE_PENDING@)) {
        double f_int[9];
        double event_integration_param;
        int actual_degree;
        const int result = find_event_time_and_state_centered(
            state_history, integration_params,
            d_non_terminal_plane_steps_past[i], w_normal, w_dist,
            &event_integration_param, f_int, &actual_degree);
        if (result != 1) {
            d_status_bundle[i] = result == 0
                ? FAILURE_PLANE_INTERPOLATION_HISTORY : FAILURE_GENERIC;
            @ESCAPE@
        }
        @INTERSECTION_SETUP@
        if (handle_non_terminal_plane_intersection(
                @INTERSECTION_STATE@, physical_lambda,
                &d_results_buffer[master_idx]@COMMONDATA_ARG@)) {
            d_non_terminal_plane_event_found[i] = true;
            d_non_terminal_plane_event_degree[i] = actual_degree;
            @NON_TERMINAL_CAPTURE@
        }
        d_non_terminal_plane_crossing_pending[i] = false;
        d_non_terminal_plane_steps_past[i] = 0;
    }

    // Shift history only for photons that will take another RK step.
    if (d_status_bundle[i] == ACTIVE) {
        for (int node = @MAX_DEGREE@ - 1; node > 0; --node) {
            for (int component = 0; component < 9; ++component) {
                WriteCUDA(&d_f_history_bundle[IDX_HISTORY(node, component, i)],
                    state_history[9 * node + component]);
            }
            WriteCUDA(&d_integration_param_history[IDX_PARAM(node, i)],
                integration_params[node]);
        }
        for (int component = 0; component < 9; ++component) {
            WriteCUDA(&d_f_history_bundle[IDX_HISTORY(0, component, i)],
                state_history[component]);
        }
        WriteCUDA(&d_integration_param_history[IDX_PARAM(0, i)],
            integration_params[0]);
    }
    #undef IDX_PARAM
    #undef IDX_HISTORY
    #undef IDX_F
"""
    substitutions = {
        "@CD@": cd_access,
        "@ESCAPE@": escape_statement,
        "@HISTORY_SIZE@": str(maximum_degree + 1),
        "@MAX_DEGREE@": str(maximum_degree),
        "@POST_STEPS@": str(required_post_steps),
        "@FORCE_PENDING@": (
            "force_pending_interpolation" if single_photon_step_limit else "false"
        ),
        "@TIME_LIMIT_CHECK@": time_limit_check,
        "@INTERSECTION_SETUP@": intersection_setup,
        "@INTERSECTION_STATE@": intersection_state,
        "@COMMONDATA_ARG@": commondata_arg,
        "@NON_TERMINAL_CAPTURE@": non_terminal_capture,
        "@TERMINAL_CAPTURE@": terminal_capture,
    }
    for marker, replacement in substitutions.items():
        core = core.replace(marker, replacement)
    prefunc_kernel, body = parallel_utils.generate_kernel_and_launch_code(
        kernel_name="event_detection_manager_kernel",
        kernel_body=f"{loop_start}\n{core}\n{loop_end}",
        arg_dict_cuda=args.copy(),
        arg_dict_host=args.copy(),
        parallelization=parallelization,
        launch_dict={
            "threads_per_block": ["256", "1", "1"],
            "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
            "stream": "stream_idx",
        },
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
        thread_tiling_macro_suffix="RKF45",
    )
    prefunc = "\n\n".join(
        [
            cfc.CFunction_dict["find_event_time_and_state_centered"].full_function,
            cfc.CFunction_dict["handle_non_terminal_plane_intersection"].full_function,
            cfc.CFunction_dict["handle_terminal_plane_intersection"].full_function,
            prefunc_kernel,
        ]
    )
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if parallelization == "cuda":
        includes.extend(["cuda_intrinsics.h", "BHaH_device_defines.h"])
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc="Detect pending plane crossings and interpolate accepted photon states.",
        cfunc_type="void",
        name="event_detection_manager_kernel",
        params=(
            "const commondata_struct *restrict commondata, "
            "const double *restrict d_f_bundle, "
            "const double *restrict d_log_energy_bundle, "
            "double *restrict d_f_history_bundle, "
            "const double *restrict d_integration_param, "
            "double *restrict d_integration_param_history, "
            "blueprint_data_t *restrict d_results_buffer, "
            "termination_type_t *restrict d_status_bundle, "
            "bool *restrict d_on_pos_non_terminal_plane_prev, "
            "bool *restrict d_on_pos_terminal_plane_prev, "
            "bool *restrict d_non_terminal_plane_event_found, "
            "bool *restrict d_terminal_plane_event_found, "
            "bool *restrict d_non_terminal_plane_crossing_pending, "
            "bool *restrict d_terminal_plane_crossing_pending, "
            "int *restrict d_non_terminal_plane_steps_past, "
            "int *restrict d_terminal_plane_steps_past, "
            "int *restrict d_non_terminal_plane_event_degree, "
            "int *restrict d_terminal_plane_event_degree, "
            f"{event_state_args}{force_arg}"
            "const long int *restrict d_chunk_buffer, "
            "const int chunk_size, "
            "const int stream_idx"
        ),
        include_CodeParameters_h=False,
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
