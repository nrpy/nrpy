"""
Register saved crossing steps and fixed-anchor RKF45 Brent refinement.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par


def register_plane_crossing_candidates() -> None:
    """Register one nonterminal and three terminal candidates per photon."""
    if "photon_plane_crossing_candidates" in par.glb_extras_dict.get(
        "BHaH_defines", {}
    ):
        return

    Bdefines_h.register_BHaH_defines(
        "photon_plane_crossing_candidates",
        r"""
    // One accepted RKF45 step that changed sides of a plane.
    typedef struct {
        bool recorded;
        long int photon_index;
        event_type_t plane_type;
        double state_before[9];
        double parameter_before;
        double parameter_after;
    } photon_plane_crossing_candidate_t;

    // Terminal candidates retain their order along the evolved photon path.
    typedef struct {
        photon_plane_crossing_candidate_t non_terminal_candidate;
        photon_plane_crossing_candidate_t terminal_candidate_1;
        photon_plane_crossing_candidate_t terminal_candidate_2;
        photon_plane_crossing_candidate_t terminal_candidate_3;
        int terminal_candidate_count;
        bool additional_terminal_crossing_detected;
    } photon_plane_crossing_candidates_t;
""",
    )


def register_plane_crossing_root() -> None:
    """Register one refined plane state and its numerical result."""
    if "photon_plane_crossing_root" in par.glb_extras_dict.get("BHaH_defines", {}):
        return
    Bdefines_h.register_BHaH_defines(
        "photon_plane_crossing_root",
        r"""
    typedef enum {
        PLANE_CROSSING_ROOT_FOUND = 0,
        PLANE_CROSSING_ROOT_INVALID_INPUT = 1,
        PLANE_CROSSING_ROOT_INVALID_BRACKET = 2,
        PLANE_CROSSING_ROOT_RKF45_FAILURE = 3,
        PLANE_CROSSING_ROOT_ITERATION_LIMIT = 4
    } photon_plane_crossing_root_status_t;

    typedef struct {
        bool found;
        photon_plane_crossing_root_status_t status;
        double state[9];
        double parameter;
        double plane_distance;
        double parameter_interval_width;
        int iterations;
    } photon_plane_crossing_root_t;
""",
    )


def plane_crossing_brent_refine_candidate(
    normalized_eom: bool = False, spacetime_name: str = ""
) -> None:
    """
    Register fixed-anchor fifth-order RKF45 and Brent plane refinement.

    :param normalized_eom: Whether coordinate time is the integration parameter.
    :param spacetime_name: Analytical metric name; empty selects numerical data.
    """
    register_plane_crossing_candidates()
    register_plane_crossing_root()
    par.register_CodeParameter(
        "REAL",
        __name__,
        "plane_crossing_residual_tolerance",
        1.0e-12,
        commondata=True,
        add_to_parfile=True,
        description="Maximum absolute signed plane distance at a refined crossing.",
    )
    par.register_CodeParameter(
        "REAL",
        __name__,
        "plane_crossing_parameter_tolerance",
        1.0e-12,
        commondata=True,
        add_to_parfile=True,
        description="Maximum integration-parameter interval width at a crossing.",
    )
    parameter_array = (
        "double integration_parameter[1] = {candidate->parameter_before};"
        if normalized_eom
        else ""
    )
    interpolation_args = "integration_parameter, h, stage, " if normalized_eom else ""
    rhs_args = "integration_parameter, h, " if normalized_eom else ""
    analytical = bool(spacetime_name)
    synthetic_args = ""
    if not analytical:
        synthetic_args = (
            "NULL, NULL, "
            if "d_photon_indices"
            in cfc.CFunction_dict["numerical_interpolation"].params
            else ""
        )
    # Existing kernels index Structure-of-Arrays fields with BUNDLE_CAPACITY.
    # Generate one-ray copies with stride one so Brent trials use small arrays.
    scalar_sources = {
        name: cfc.CFunction_dict[name].full_function
        for name in ("calculate_ode_rhs_kernel", "rkf45_stage_update")
    }
    if analytical:
        scalar_sources[f"g4DD_metric_{spacetime_name}"] = cfc.CFunction_dict[
            f"g4DD_metric_{spacetime_name}"
        ].full_function
        scalar_sources[f"connections_{spacetime_name}"] = cfc.CFunction_dict[
            f"connections_{spacetime_name}"
        ].full_function
        if par.parval_from_str("parallelization") == "cuda":
            # Keep analytical CUDA imports local to this refinement path.
            # pylint: disable=import-outside-toplevel
            from nrpy.equations.general_relativity.geodesics.geodesics import (
                Geodesic_Equations,
            )
            from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import (
                calculate_ode_rhs_kernel,
                rkf45_stage_update,
            )

            prior_mode = par.parval_from_str("parallelization")
            prior_rhs = cfc.CFunction_dict["calculate_ode_rhs_kernel"]
            prior_stage = cfc.CFunction_dict["rkf45_stage_update"]
            try:
                par.set_parval_from_str("parallelization", "openmp")
                cfc.CFunction_dict.pop("calculate_ode_rhs_kernel")
                cfc.CFunction_dict.pop("rkf45_stage_update")
                geodesic_data = Geodesic_Equations[f"{spacetime_name}_photon"]
                rhs = (
                    geodesic_data.geodesic_eom_rhs_photon_normalized_christoffel()
                    if normalized_eom
                    else geodesic_data.geodesic_eom_rhs_photon_christoffel()
                )
                calculate_ode_rhs_kernel.calculate_ode_rhs_kernel(
                    rhs,
                    geodesic_data.xx,
                    rhs_uses_metric_derivatives=False,
                    normalized_eom=normalized_eom,
                )
                rkf45_stage_update.rkf45_stage_update()
                scalar_sources["calculate_ode_rhs_kernel"] = cfc.CFunction_dict[
                    "calculate_ode_rhs_kernel"
                ].full_function
                scalar_sources["rkf45_stage_update"] = cfc.CFunction_dict[
                    "rkf45_stage_update"
                ].full_function
            finally:
                cfc.CFunction_dict["calculate_ode_rhs_kernel"] = prior_rhs
                cfc.CFunction_dict["rkf45_stage_update"] = prior_stage
                par.set_parval_from_str("parallelization", prior_mode)
    else:
        scalar_sources["numerical_interpolation"] = cfc.CFunction_dict[
            "numerical_interpolation"
        ].full_function

    scalar_functions = []
    for function_name, source in scalar_sources.items():
        generated = source
        generated = generated.replace("BUNDLE_CAPACITY", "1")
        if analytical and par.parval_from_str("parallelization") == "cuda":
            if function_name == "rkf45_stage_update":
                for intrinsic in (
                    "ReadCUDA",
                    "WriteCUDA",
                    "MulCUDA",
                    "FusedMulAddCUDA",
                ):
                    generated = generated.replace(
                        f"{intrinsic}(", f"plane_crossing_host_{intrinsic}("
                    )
            if function_name in (
                f"g4DD_metric_{spacetime_name}",
                f"connections_{spacetime_name}",
            ):
                generated = generated.replace("BHAH_HD_INLINE", "static inline")
        generated = "".join(
            (
                line[len(line.rstrip("\r\n")) :]
                if line.strip() == "#pragma omp parallel for"
                else line
            )
            for line in generated.splitlines(keepends=True)
        )
        if function_name == "numerical_interpolation":
            identifier = "numerical_interpolation"
            renamed = []
            start = 0
            position = generated.find(identifier, start)
            while position != -1:
                before = generated[position - 1 : position]
                after = generated[
                    position + len(identifier) : position + len(identifier) + 1
                ]
                renamed.append(generated[start:position])
                renamed.append(
                    identifier
                    if before.isalnum()
                    or before == "_"
                    or after.isalnum()
                    or after == "_"
                    else "plane_crossing_numerical_interpolation"
                )
                start = position + len(identifier)
                position = generated.find(identifier, start)
            renamed.append(generated[start:])
            generated = "".join(renamed)
        elif function_name in ("calculate_ode_rhs_kernel", "rkf45_stage_update"):
            generated = generated.replace(
                function_name, f"plane_crossing_{function_name}"
            )
        scalar_functions.append(generated)
    trial_geometry_params = (
        ""
        if analytical
        else "const params_struct *numerical_params,\n"
        "    const azimuthal_symmetry_spatial_lagrange_context_struct *spatial_context,\n"
        "    const NumericalTimeWindowManager *numerical_window,\n    "
    )
    trial_stencil_params = (
        "" if analytical else "const int stencil_i0, const int stencil_i2,\n    "
    )
    trial_stencil_setup = (
        ""
        if analytical
        else "int selected_i0 = stencil_i0;\n    int selected_i2 = stencil_i2;"
    )
    trial_geometry_call = (
        f"g4DD_metric_{spacetime_name}(commondata, stage_state, metric);\n"
        f"        connections_{spacetime_name}(commondata, stage_state, geometry);"
        if analytical
        else "plane_crossing_numerical_interpolation(\n"
        "            commondata, numerical_params, spatial_context, numerical_window,\n"
        "            stage_state, &status, &selected_i0, &selected_i2,\n"
        f"            {interpolation_args}{synthetic_args}metric, geometry,\n"
        "            chunk_size, stream_idx);"
    )
    trial_status_declaration = (
        "" if analytical else "termination_type_t status = ACTIVE;"
    )
    trial_status_check = "" if analytical else "if (status != ACTIVE) return false;"
    host_intrinsics = ""
    if analytical and par.parval_from_str("parallelization") == "cuda":
        host_intrinsics = r"""
static double plane_crossing_host_ReadCUDA(const double *value) { return *value; }
static void plane_crossing_host_WriteCUDA(double *target, const double value) {
    *target = value;
}
static double plane_crossing_host_MulCUDA(const double a, const double b) {
    return a * b;
}
static double plane_crossing_host_FusedMulAddCUDA(
    const double a, const double b, const double c) { return fma(a, b, c); }
"""
    prefunc = host_intrinsics + "\n".join(scalar_functions) + rf"""
static double plane_crossing_signed_distance(
    const double *state, const double *normal, const double *center) {{
    return normal[0] * (state[1] - center[0]) +
           normal[1] * (state[2] - center[1]) +
           normal[2] * (state[3] - center[2]);
}}

static bool plane_crossing_rk5_trial(
    const commondata_struct *commondata,
    {trial_geometry_params}
    const photon_plane_crossing_candidate_t *candidate,
    {trial_stencil_params}
    const double trial_parameter, double *trial_state) {{
    const long int chunk_size = 1;
    const int stream_idx = 0;
    const double h_value = trial_parameter - candidate->parameter_before;
    if (!isfinite(h_value)) return false;
    double h[1] = {{h_value}};
    {parameter_array}
    double stage_state[9];
    double metric[10];
    double geometry[40];
    double k[54] = {{0.0}};
    {trial_status_declaration}
    {trial_stencil_setup}
    for (int component = 0; component < 9; ++component) {{
        stage_state[component] = candidate->state_before[component];
    }}
    for (int stage = 1; stage <= 6; ++stage) {{
        {trial_geometry_call}
        {trial_status_check}
        for (int component = 0; component < 10; ++component) {{
            if (!isfinite(metric[component])) return false;
        }}
        for (int component = 0; component < 40; ++component) {{
            if (!isfinite(geometry[component])) return false;
        }}
        plane_crossing_calculate_ode_rhs_kernel(
            stage_state, metric, geometry, {rhs_args}k,
            stage, chunk_size, stream_idx);
        for (int component = 0; component < 9; ++component) {{
            if (!isfinite(k[(stage - 1) * 9 + component])) return false;
        }}
        if (stage < 6) {{
            plane_crossing_rkf45_stage_update(
                candidate->state_before, k, h, stage,
                (int)chunk_size, stage_state, stream_idx);
        }}
    }} // END LOOP: evaluate RKF45 stages
    for (int component = 0; component < 9; ++component) {{
        const double update =
            (16.0 / 135.0) * k[component] +
            (6656.0 / 12825.0) * k[18 + component] +
            (28561.0 / 56430.0) * k[27 + component] -
            (9.0 / 50.0) * k[36 + component] +
            (2.0 / 55.0) * k[45 + component];
        trial_state[component] = candidate->state_before[component] +
                                 h_value * update;
        if (!isfinite(trial_state[component])) return false;
    }} // END LOOP: update trial photon state
    return true;
}} // END FUNCTION: plane_crossing_rk5_trial
"""
    body = r"""
    enum { MAX_BRENT_ITERATIONS = 75 };
    *root = (photon_plane_crossing_root_t){0};
    root->status = PLANE_CROSSING_ROOT_INVALID_INPUT;
    root->parameter = NAN;
    root->plane_distance = NAN;
    root->parameter_interval_width = NAN;
    for (int component = 0; component < 9; ++component) {
        root->state[component] = NAN;
    }
    if (!candidate->recorded ||
        !isfinite(candidate->parameter_before) ||
        !isfinite(candidate->parameter_after) ||
        candidate->parameter_before == candidate->parameter_after ||
        !isfinite(commondata->plane_crossing_residual_tolerance) ||
        commondata->plane_crossing_residual_tolerance <= 0.0 ||
        !isfinite(commondata->plane_crossing_parameter_tolerance) ||
        commondata->plane_crossing_parameter_tolerance <= 0.0) return;
    for (int component = 0; component < 9; ++component) {
        if (!isfinite(candidate->state_before[component])) return;
    }

    double normal[3], center[3];
    if (candidate->plane_type == NON_TERMINAL_PLANE_EVENT) {
        normal[0] = commondata->non_terminal_plane_normal_x;
        normal[1] = commondata->non_terminal_plane_normal_y;
        normal[2] = commondata->non_terminal_plane_normal_z;
        center[0] = commondata->non_terminal_plane_center_x;
        center[1] = commondata->non_terminal_plane_center_y;
        center[2] = commondata->non_terminal_plane_center_z;
    } else if (candidate->plane_type == TERMINAL_PLANE_EVENT) {
        normal[0] = commondata->terminal_plane_normal_x;
        normal[1] = commondata->terminal_plane_normal_y;
        normal[2] = commondata->terminal_plane_normal_z;
        center[0] = commondata->terminal_plane_center_x;
        center[1] = commondata->terminal_plane_center_y;
        center[2] = commondata->terminal_plane_center_z;
    } else return;
    const double normal_length = hypot(hypot(normal[0], normal[1]), normal[2]);
    if (!isfinite(normal_length) || normal_length <= 1.0e-14) return;
    for (int axis = 0; axis < 3; ++axis) {
        if (!isfinite(center[axis])) return;
        normal[axis] /= normal_length;
    }

    const REAL anchor_cartesian[3] = {
        (REAL)candidate->state_before[1],
        (REAL)candidate->state_before[2],
        (REAL)candidate->state_before[3]};
    REAL anchor_native[3];
    int automatic_center[3], selected_center[3];
    if (time_window_manager_numerical_resolve_spatial_target_and_stencil(
            numerical_params, anchor_cartesian, NULL, anchor_native,
            automatic_center, selected_center) !=
        TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        root->status = PLANE_CROSSING_ROOT_RKF45_FAILURE;
        return;
    }

    double endpoint_state[9];
    if (!plane_crossing_rk5_trial(
            commondata, numerical_params, spatial_context, numerical_window,
            candidate, selected_center[0], selected_center[2],
            candidate->parameter_after, endpoint_state)) {
        root->status = PLANE_CROSSING_ROOT_RKF45_FAILURE;
        return;
    }
    const double anchor_distance = plane_crossing_signed_distance(
        candidate->state_before, normal, center);
    const double endpoint_distance = plane_crossing_signed_distance(
        endpoint_state, normal, center);
    if (!isfinite(anchor_distance) || !isfinite(endpoint_distance)) {
        root->status = PLANE_CROSSING_ROOT_RKF45_FAILURE;
        return;
    }
    double low = fmin(candidate->parameter_before, candidate->parameter_after);
    double high = fmax(candidate->parameter_before, candidate->parameter_after);
    double f_low = candidate->parameter_before < candidate->parameter_after
                       ? anchor_distance : endpoint_distance;
    double f_high = candidate->parameter_before < candidate->parameter_after
                        ? endpoint_distance : anchor_distance;
    double previous_parameter = NAN, previous_distance = NAN;
    bool have_previous = false;

    // Endpoint roots have zero parameter interval width.
    if (anchor_distance == 0.0 || endpoint_distance == 0.0) {
        const bool at_anchor = anchor_distance == 0.0;
        root->found = true;
        root->status = PLANE_CROSSING_ROOT_FOUND;
        root->parameter = at_anchor ? candidate->parameter_before
                                    : candidate->parameter_after;
        root->plane_distance = 0.0;
        root->parameter_interval_width = 0.0;
        const double *state = at_anchor ? candidate->state_before : endpoint_state;
        for (int component = 0; component < 9; ++component) {
            root->state[component] = state[component];
        }
        return;
    }

    root->status = PLANE_CROSSING_ROOT_ITERATION_LIMIT;
    for (int iteration = 1;
         iteration <= MAX_BRENT_ITERATIONS;
         ++iteration) {
        // A crossing interval is required for each Brent update.
        if (signbit(f_low) == signbit(f_high)) {
            root->status = PLANE_CROSSING_ROOT_INVALID_BRACKET;
            return;
        }
        const double width = high - low;
        double trial_parameter = low + 0.5 * width;
        if (have_previous && f_low != f_high &&
            f_low != previous_distance && f_high != previous_distance) {
            trial_parameter =
                low * f_high * previous_distance /
                    ((f_low - f_high) * (f_low - previous_distance)) +
                high * f_low * previous_distance /
                    ((f_high - f_low) * (f_high - previous_distance)) +
                previous_parameter * f_low * f_high /
                    ((previous_distance - f_low) * (previous_distance - f_high));
        } else if (f_low != f_high) {
            trial_parameter = high - f_high * width / (f_high - f_low);
        }
        if (!isfinite(trial_parameter) ||
            trial_parameter <= low + 0.1 * width ||
            trial_parameter >= high - 0.1 * width) {
            trial_parameter = low + 0.5 * width;
        }
        if (trial_parameter <= low || trial_parameter >= high) return;

        double trial_state[9];
        if (!plane_crossing_rk5_trial(
                commondata, numerical_params, spatial_context, numerical_window,
                candidate, selected_center[0], selected_center[2],
                trial_parameter, trial_state)) {
            root->status = PLANE_CROSSING_ROOT_RKF45_FAILURE;
            return;
        }
        const double trial_distance = plane_crossing_signed_distance(
            trial_state, normal, center);
        if (!isfinite(trial_distance)) {
            root->status = PLANE_CROSSING_ROOT_RKF45_FAILURE;
            return;
        }
        root->iterations = iteration;
        root->parameter = trial_parameter;
        root->plane_distance = trial_distance;
        if (trial_distance == 0.0) {
            low = high = trial_parameter;
            f_low = f_high = 0.0;
        } else {
            const bool replace_low = signbit(trial_distance) == signbit(f_low);
            previous_parameter = replace_low ? low : high;
            previous_distance = replace_low ? f_low : f_high;
            have_previous = true;
            if (replace_low) {
                low = trial_parameter;
                f_low = trial_distance;
            } else {
                high = trial_parameter;
                f_high = trial_distance;
            }
        }
        root->parameter_interval_width = high - low;
        if (fabs(trial_distance) <= commondata->plane_crossing_residual_tolerance &&
            root->parameter_interval_width <=
                commondata->plane_crossing_parameter_tolerance) {
            for (int component = 0; component < 9; ++component) {
                root->state[component] = trial_state[component];
            }
            root->found = true;
            root->status = PLANE_CROSSING_ROOT_FOUND;
            return;
        }
    }
"""
    if analytical:
        spatial_setup = """    const REAL anchor_cartesian[3] = {
        (REAL)candidate->state_before[1],
        (REAL)candidate->state_before[2],
        (REAL)candidate->state_before[3]};
    REAL anchor_native[3];
    int automatic_center[3], selected_center[3];
    if (time_window_manager_numerical_resolve_spatial_target_and_stencil(
            numerical_params, anchor_cartesian, NULL, anchor_native,
            automatic_center, selected_center) !=
        TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        root->status = PLANE_CROSSING_ROOT_RKF45_FAILURE;
        return;
    }

"""
        body = body.replace(spatial_setup, "")
        body = body.replace(
            "commondata, numerical_params, spatial_context, numerical_window,\n"
            "            candidate, selected_center[0], selected_center[2],",
            "commondata, candidate,",
        ).replace(
            "commondata, numerical_params, spatial_context, numerical_window,\n"
            "                candidate, selected_center[0], selected_center[2],",
            "commondata, candidate,",
        )
    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        prefunc=prefunc,
        desc="Refine one photon plane crossing by fixed-anchor RKF45 and Brent trials.",
        cfunc_type="void",
        name=(
            "plane_crossing_brent_refine_candidate_analytical"
            if analytical
            else "plane_crossing_brent_refine_candidate"
        ),
        params=(
            "const commondata_struct *restrict commondata, "
            + (
                ""
                if analytical
                else "const params_struct *restrict numerical_params, "
                "const azimuthal_symmetry_spatial_lagrange_context_struct *restrict spatial_context, "
                "const NumericalTimeWindowManager *restrict numerical_window, "
            )
            + "const photon_plane_crossing_candidate_t *restrict candidate, "
            "photon_plane_crossing_root_t *restrict root"
        ),
        include_CodeParameters_h=False,
        body=body,
    )


def find_plane_crossings_brent(normalized_eom: bool = False) -> None:
    """
    Register slot-ordered refinement and selection of numerical plane crossings.

    :param normalized_eom: Whether coordinate time is the integration parameter.
    :raises ValueError: If the single-candidate Brent function is not registered.
    """
    if "plane_crossing_brent_refine_candidate" not in cfc.CFunction_dict:
        raise ValueError("Register plane_crossing_brent_refine_candidate first.")
    anchor_time = (
        "candidate->parameter_before"
        if normalized_eom
        else "candidate->state_before[0]"
    )
    non_terminal_precedes_terminal = (
        "non_terminal->integration_param >= selected->integration_param"
        if normalized_eom
        else "non_terminal->integration_param <= selected->integration_param"
    )
    body = r"""
    if (num_rays <= 0) return;
    if (num_rays > LONG_MAX / 4 || (size_t)num_rays > SIZE_MAX /
            (3 * sizeof(terminal_plane_crossing_result_t))) {
        fprintf(stderr, "ERROR: too many photon plane candidates.\n");
        exit(EXIT_FAILURE);
    }
    const long int candidate_capacity = 4 * num_rays;
    terminal_plane_crossing_result_t *terminal_candidates =
        (terminal_plane_crossing_result_t *)calloc(
            (size_t)(3 * num_rays), sizeof(*terminal_candidates));
    const long int chunk_capacity = NRPYMIN((long int)BUNDLE_CAPACITY,
                                            candidate_capacity);
    long int *chunk_indices =
        (long int *)malloc((size_t)chunk_capacity * sizeof(*chunk_indices));
    if (terminal_candidates == NULL || chunk_indices == NULL) {
        fprintf(stderr, "ERROR: could not allocate plane crossing results.\n");
        free(terminal_candidates);
        free(chunk_indices);
        exit(EXIT_FAILURE);
    }

    TimeSlotManager crossing_slots;
    slot_manager_init(&crossing_slots, commondata->slot_manager_t_min,
                      slot_manager_t_max, commondata->slot_manager_delta_t,
                      candidate_capacity);
    long int failures[5] = {0};
    long int outside_slots = 0;
    long int plane_geometry_failures = 0;
    long int terminal_overflow_count = 0;
    for (long int photon = 0; photon < num_rays; ++photon) {
        if (candidates[photon].additional_terminal_crossing_detected)
            ++terminal_overflow_count;
        non_terminal_results[photon] = (non_terminal_plane_crossing_result_t){0};
        terminal_results[photon] = (terminal_plane_crossing_result_t){0};
        non_terminal_results[photon].null_constraint_abs = NAN;
        non_terminal_results[photon].Lz = NAN;
        terminal_results[photon].null_constraint_abs = NAN;
        terminal_results[photon].Lz = NAN;
        for (int candidate_number = 0; candidate_number < 4;
             ++candidate_number) {
            const photon_plane_crossing_candidate_t *candidate = NULL;
            if (candidate_number == 0)
                candidate = &candidates[photon].non_terminal_candidate;
            else if (candidate_number == 1)
                candidate = &candidates[photon].terminal_candidate_1;
            else if (candidate_number == 2)
                candidate = &candidates[photon].terminal_candidate_2;
            else
                candidate = &candidates[photon].terminal_candidate_3;
            if (!candidate->recorded) continue;
            const double candidate_time = @ANCHOR_TIME@;
            const int slot = slot_get_index(&crossing_slots, candidate_time);
            if (slot < 0) {
                ++outside_slots;
                continue;
            }
            slot_add_photon(&crossing_slots, slot,
                            4 * photon + candidate_number);
        }
    }

    for (int slot = crossing_slots.num_slots - 1; slot >= 0; --slot) {
        if (crossing_slots.slot_counts[slot] == 0) continue;
        const bool window_available =
            time_window_manager_numerical_mmap_for_slot(
                numerical_window, &crossing_slots, slot) ==
            TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
        while (crossing_slots.slot_counts[slot] > 0) {
            const long int chunk_size = NRPYMIN(
                chunk_capacity, crossing_slots.slot_counts[slot]);
            slot_remove_chunk(&crossing_slots, slot, chunk_indices, chunk_size);
            #pragma omp parallel for
            for (long int lane = 0; lane < chunk_size; ++lane) {
                const long int candidate_index = chunk_indices[lane];
                const long int photon = candidate_index / 4;
                const int candidate_number = (int)(candidate_index % 4);
                const photon_plane_crossing_candidate_t *candidate = NULL;
                if (candidate_number == 0)
                    candidate = &candidates[photon].non_terminal_candidate;
                else if (candidate_number == 1)
                    candidate = &candidates[photon].terminal_candidate_1;
                else if (candidate_number == 2)
                    candidate = &candidates[photon].terminal_candidate_2;
                else
                    candidate = &candidates[photon].terminal_candidate_3;
                if (!window_available) {
                    #pragma omp atomic update
                    ++failures[PLANE_CROSSING_ROOT_RKF45_FAILURE];
                    continue;
                }
                photon_plane_crossing_root_t root;
                plane_crossing_brent_refine_candidate(
                    commondata, numerical_params, spatial_context,
                    numerical_window, candidate, &root);
                if (!root.found) {
                    #pragma omp atomic update
                    ++failures[root.status];
                    continue;
                }
                bool plane_valid;
                if (candidate_number == 0) {
                    plane_valid = handle_non_terminal_plane_intersection(
                        &root, &non_terminal_results[photon], commondata);
                } else {
                    terminal_plane_crossing_result_t *result =
                        &terminal_candidates[3 * photon + candidate_number - 1];
                    plane_valid = handle_terminal_plane_intersection(
                        &root, result, commondata);
                }
                if (!plane_valid) {
                    #pragma omp atomic update
                    ++plane_geometry_failures;
                }
            }
        }
    }

    for (long int photon = 0; photon < num_rays; ++photon) {
        terminal_plane_crossing_result_t *selected = &terminal_results[photon];
        for (int candidate_number = 1; candidate_number <= 3;
             ++candidate_number) {
            const terminal_plane_crossing_result_t *trial =
                &terminal_candidates[3 * photon + candidate_number - 1];
            if (trial->found && trial->inside_radial_bounds) {
                *selected = *trial;
                selected->selected_candidate_number = candidate_number;
                break;
            }
        }
        non_terminal_plane_crossing_result_t *non_terminal =
            &non_terminal_results[photon];
        if (selected->found && non_terminal->found &&
            !(@NON_TERMINAL_PRECEDES_TERMINAL@)) {
            *non_terminal = (non_terminal_plane_crossing_result_t){0};
            non_terminal->null_constraint_abs = NAN;
            non_terminal->Lz = NAN;
        }
    }
    if (terminal_overflow_count || outside_slots ||
        plane_geometry_failures || failures[1] ||
        failures[2] || failures[3] || failures[4]) {
        fprintf(stderr,
                "Plane crossing refinement: terminal crossing cap reached=%ld, "
                "outside time slots=%ld, "
                "invalid input=%ld, invalid bracket=%ld, "
                "RKF45 or numerical geometry failure=%ld, "
                "iteration limit=%ld, plane geometry failure=%ld.\n",
                terminal_overflow_count, outside_slots,
                failures[1], failures[2], failures[3],
                failures[4], plane_geometry_failures);
    }
    slot_manager_free(&crossing_slots);
    free(chunk_indices);
    free(terminal_candidates);
""".replace("@ANCHOR_TIME@", anchor_time).replace(
        "@NON_TERMINAL_PRECEDES_TERMINAL@", non_terminal_precedes_terminal
    )
    cfc.register_CFunction(
        includes=[
            "BHaH_defines.h",
            "BHaH_function_prototypes.h",
            "<limits.h>",
            "<stdint.h>",
            "<stdlib.h>",
            "<stdio.h>",
            "<math.h>",
        ],
        desc="Refine saved crossing steps in numerical time slots and select plane results.",
        cfunc_type="void",
        name="find_plane_crossings_brent",
        params=(
            "const commondata_struct *restrict commondata, "
            "const params_struct *restrict numerical_params, "
            "const azimuthal_symmetry_spatial_lagrange_context_struct *restrict spatial_context, "
            "NumericalTimeWindowManager *restrict numerical_window, "
            "const double slot_manager_t_max, "
            "const long int num_rays, "
            "const photon_plane_crossing_candidates_t *restrict candidates, "
            "non_terminal_plane_crossing_result_t *restrict non_terminal_results, "
            "terminal_plane_crossing_result_t *restrict terminal_results"
        ),
        include_CodeParameters_h=False,
        body=body,
    )


def find_plane_crossings_brent_analytical(
    spacetime_name: str, normalized_eom: bool = False
) -> None:
    """
    Register analytical Brent refinement and select plane crossings.

    :param spacetime_name: Analytical spacetime metric name.
    :param normalized_eom: Whether coordinate time is the integration parameter.
    """
    plane_crossing_brent_refine_candidate(normalized_eom, spacetime_name)
    projection_suffix = (
        "_host" if par.parval_from_str("parallelization") == "cuda" else ""
    )
    non_terminal_precedes_terminal = (
        "non_terminal->integration_param >= selected->integration_param"
        if normalized_eom
        else "non_terminal->integration_param <= selected->integration_param"
    )
    body = r"""
    if (num_rays <= 0) return;
    if ((size_t)num_rays > SIZE_MAX / (3 * sizeof(terminal_plane_crossing_result_t))) {
        fprintf(stderr, "ERROR: too many photon plane candidates.\n");
        exit(EXIT_FAILURE);
    }
    terminal_plane_crossing_result_t *terminal_candidates =
        (terminal_plane_crossing_result_t *)calloc(
            (size_t)(3 * num_rays), sizeof(*terminal_candidates));
    if (terminal_candidates == NULL) {
        fprintf(stderr, "ERROR: could not allocate plane crossing results.\n");
        exit(EXIT_FAILURE);
    }
    long int failures[5] = {0};
    long int plane_geometry_failures = 0;
    long int terminal_overflow_count = 0;
    #pragma omp parallel for reduction(+:plane_geometry_failures,terminal_overflow_count)
    for (long int photon = 0; photon < num_rays; ++photon) {
        if (candidates[photon].additional_terminal_crossing_detected)
            ++terminal_overflow_count;
        non_terminal_results[photon] = (non_terminal_plane_crossing_result_t){0};
        terminal_results[photon] = (terminal_plane_crossing_result_t){0};
        non_terminal_results[photon].null_constraint_abs = NAN;
        non_terminal_results[photon].Lz = NAN;
        terminal_results[photon].null_constraint_abs = NAN;
        terminal_results[photon].Lz = NAN;
        for (int candidate_number = 0; candidate_number < 4; ++candidate_number) {
            const photon_plane_crossing_candidate_t *candidate = NULL;
            if (candidate_number == 0)
                candidate = &candidates[photon].non_terminal_candidate;
            else if (candidate_number == 1)
                candidate = &candidates[photon].terminal_candidate_1;
            else if (candidate_number == 2)
                candidate = &candidates[photon].terminal_candidate_2;
            else
                candidate = &candidates[photon].terminal_candidate_3;
            if (!candidate->recorded) continue;
            photon_plane_crossing_root_t root;
            plane_crossing_brent_refine_candidate_analytical(
                commondata, candidate, &root);
            if (!root.found) {
                #pragma omp atomic update
                ++failures[root.status];
                continue;
            }
            bool plane_valid;
            if (candidate_number == 0) {
                plane_valid = handle_non_terminal_plane_intersection@HOST@(
                    &root, &non_terminal_results[photon], commondata);
            } else {
                terminal_plane_crossing_result_t *result =
                    &terminal_candidates[3 * photon + candidate_number - 1];
                plane_valid = handle_terminal_plane_intersection@HOST@(
                    &root, result, commondata);
            }
            if (!plane_valid) ++plane_geometry_failures;
        }
        terminal_plane_crossing_result_t *selected = &terminal_results[photon];
        for (int candidate_number = 1; candidate_number <= 3; ++candidate_number) {
            const terminal_plane_crossing_result_t *trial =
                &terminal_candidates[3 * photon + candidate_number - 1];
            if (trial->found && trial->inside_radial_bounds) {
                *selected = *trial;
                selected->selected_candidate_number = candidate_number;
                break;
            }
        }
        non_terminal_plane_crossing_result_t *non_terminal =
            &non_terminal_results[photon];
        if (selected->found && non_terminal->found &&
            !(@NON_TERMINAL_PRECEDES_TERMINAL@)) {
            *non_terminal = (non_terminal_plane_crossing_result_t){0};
            non_terminal->null_constraint_abs = NAN;
            non_terminal->Lz = NAN;
        }
    }
    if (terminal_overflow_count || plane_geometry_failures || failures[1] ||
        failures[2] || failures[3] || failures[4]) {
        fprintf(stderr,
                "Plane crossing refinement: terminal crossing cap reached=%ld, "
                "invalid input=%ld, invalid bracket=%ld, "
                "RKF45 or analytical geometry failure=%ld, "
                "iteration limit=%ld, plane geometry failure=%ld.\n",
                terminal_overflow_count, failures[1], failures[2], failures[3],
                failures[4], plane_geometry_failures);
    }
    free(terminal_candidates);
""".replace("@HOST@", projection_suffix).replace(
        "@NON_TERMINAL_PRECEDES_TERMINAL@", non_terminal_precedes_terminal
    )
    cfc.register_CFunction(
        includes=[
            "BHaH_defines.h",
            "BHaH_function_prototypes.h",
            "<stdint.h>",
            "<stdlib.h>",
            "<stdio.h>",
            "<math.h>",
        ],
        desc="Refine analytical photon plane crossings and select plane results.",
        cfunc_type="void",
        name="find_plane_crossings_brent_analytical",
        params=(
            "const commondata_struct *restrict commondata, "
            "const long int num_rays, "
            "const photon_plane_crossing_candidates_t *restrict candidates, "
            "non_terminal_plane_crossing_result_t *restrict non_terminal_results, "
            "terminal_plane_crossing_result_t *restrict terminal_results"
        ),
        include_CodeParameters_h=False,
        body=body,
    )
