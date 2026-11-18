# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/event_detection_manager_kernel.py
"""
Record accepted photon steps that cross either configured plane.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.puncture_diagnostics import (
    register_apparent_horizons,
)


def event_detection_manager_kernel(
    normalized_eom: bool = False, analytical: bool = False
) -> None:
    """
    Register accepted-step plane crossings and physical photon stops.

    :param normalized_eom: Whether coordinate time is the integration parameter.
    :param analytical: Whether analytic geometry has no numerical time window.
    """
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
    coordinate_time = (
        "d_integration_param_after[i]" if normalized_eom else "d_f_bundle[IDX_F(0, i)]"
    )
    horizon_params = ""
    horizon_check = ""
    horizon_desc = ""
    if not analytical:
        register_apparent_horizons()
        horizon_params = (
            ", const apparent_horizon_file_t *horizons, "
            "const NumericalTimeWindowManager *numerical_window, "
            "puncture_result_t *photon_horizons, const bool fast_light"
        )
        horizon_desc = """
        @param[in] horizons Apparent-horizon mapping, with an empty mapping disabling the check.
        @param[in] numerical_window Retained numerical spacetime slice times.
        @param[in,out] photon_horizons Latest horizon classification indexed by photon.
        @param fast_light Whether to classify at the fixed geometry sampling time.
"""
        horizon_check = r"""
        if (horizons->bytes != NULL) {
            const double point[3] = {x, y, z};
            const double lookup_time = fast_light ? commondata->t_start : coordinate_time;
            photon_horizons[photon_index] =
                apparent_horizons_classify(horizons, numerical_window, lookup_time, point);
            if (photon_horizons[photon_index].horizon_status == HORIZON_INSIDE)
                d_status_bundle[i] = STOP_CONDITION_INSIDE_APPARENT_HORIZON;
        } // END IF: apparent horizon checking enabled
"""
    body = r"""
    #define IDX_F(component, ray) ((component) * BUNDLE_CAPACITY + (ray))
    (void)stream_idx;
    #pragma omp parallel for
    for (long int i = 0; i < chunk_size; ++i) {
        // Only an accepted RKF45 step has a new state and crossing interval.
        if (d_status_bundle[i] != ACTIVE) continue;
        const long int photon_index = d_chunk_buffer[i];
        photon_plane_crossing_candidates_t *crossings =
            &plane_crossing_candidates[photon_index];
        const double x = d_f_bundle[IDX_F(1, i)];
        const double y = d_f_bundle[IDX_F(2, i)];
        const double z = d_f_bundle[IDX_F(3, i)];

        const double coordinate_time = @COORDINATE_TIME@;
        @HORIZON_CHECK@

        // Retain plane crossing candidates on the accepted stopping step.
        if (commondata->non_terminal_plane_enabled &&
            !crossings->non_terminal_candidate.recorded) {
            const double nx = commondata->non_terminal_plane_normal_x;
            const double ny = commondata->non_terminal_plane_normal_y;
            const double nz = commondata->non_terminal_plane_normal_z;
            const double normal_sq = nx * nx + ny * ny + nz * nz;
            if (!isfinite(normal_sq) || normal_sq <= 1.0e-28) {
                if (d_status_bundle[i] == ACTIVE) d_status_bundle[i] = FAILURE_GENERIC;
                continue;
            } // END IF: nonterminal plane normal invalid
            const double signed_distance =
                (x - commondata->non_terminal_plane_center_x) * nx +
                (y - commondata->non_terminal_plane_center_y) * ny +
                (z - commondata->non_terminal_plane_center_z) * nz;
            if (!isfinite(signed_distance)) {
                if (d_status_bundle[i] == ACTIVE) d_status_bundle[i] = FAILURE_GENERIC;
                continue;
            } // END IF: nonterminal plane distance nonfinite
            const bool side = signed_distance > 0.0;
            if (side != d_on_pos_non_terminal_plane_prev[i]) {
                photon_plane_crossing_candidate_t *candidate =
                    &crossings->non_terminal_candidate;
                for (int component = 0; component < 9; ++component) {
                    candidate->state_before[component] =
                        d_f_start_bundle[IDX_F(component, i)];
                } // END LOOP: component saves nonterminal anchor
                candidate->parameter_before = d_integration_param_before[i];
                candidate->parameter_after = d_integration_param_after[i];
                candidate->recorded = true;
            } // END IF: record nonterminal crossing step
            d_on_pos_non_terminal_plane_prev[i] = side;
        } // END IF: detect enabled nonterminal crossing

        if (commondata->terminal_plane_enabled &&
            !crossings->additional_terminal_crossing_detected) {
            const double nx = commondata->terminal_plane_normal_x;
            const double ny = commondata->terminal_plane_normal_y;
            const double nz = commondata->terminal_plane_normal_z;
            const double normal_sq = nx * nx + ny * ny + nz * nz;
            if (!isfinite(normal_sq) || normal_sq <= 1.0e-28) {
                if (d_status_bundle[i] == ACTIVE) d_status_bundle[i] = FAILURE_GENERIC;
                continue;
            } // END IF: terminal plane normal invalid
            const double signed_distance =
                (x - commondata->terminal_plane_center_x) * nx +
                (y - commondata->terminal_plane_center_y) * ny +
                (z - commondata->terminal_plane_center_z) * nz;
            if (!isfinite(signed_distance)) {
                if (d_status_bundle[i] == ACTIVE) d_status_bundle[i] = FAILURE_GENERIC;
                continue;
            } // END IF: terminal plane distance nonfinite
            const bool side = signed_distance > 0.0;
            if (side != d_on_pos_terminal_plane_prev[i]) {
                photon_plane_crossing_candidate_t *candidate = NULL;
                if (crossings->terminal_candidate_count == 0) {
                    candidate = &crossings->terminal_candidate_1;
                } // END IF: select first terminal candidate
                else if (crossings->terminal_candidate_count == 1) {
                    candidate = &crossings->terminal_candidate_2;
                } // END ELSE IF: select second terminal candidate
                else if (crossings->terminal_candidate_count == 2) {
                    candidate = &crossings->terminal_candidate_3;
                } // END ELSE IF: select third terminal candidate
                else {
                    crossings->additional_terminal_crossing_detected = true;
                } // END ELSE: terminal crossing capacity exceeded
                if (candidate != NULL) {
                    for (int component = 0; component < 9; ++component) {
                        candidate->state_before[component] =
                            d_f_start_bundle[IDX_F(component, i)];
                    } // END LOOP: component saves terminal anchor
                    candidate->parameter_before = d_integration_param_before[i];
                    candidate->parameter_after = d_integration_param_after[i];
                    candidate->recorded = true;
                    ++crossings->terminal_candidate_count;
                } // END IF: save available terminal candidate
            } // END IF: record terminal crossing step
            d_on_pos_terminal_plane_prev[i] = side;
        } // END IF: detect enabled terminal crossing

        const double log_energy_measure = d_log_energy_bundle[i];
        if (d_status_bundle[i] == ACTIVE && log_energy_measure > commondata->evolution_measure_max) {
            d_status_bundle[i] = STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED;
        } // END IF: photon energy limit exceeded
        const double radius_squared = x * x + y * y + z * z;
        if (radius_squared > commondata->r_escape * commondata->r_escape &&
            d_status_bundle[i] == ACTIVE) {
            d_status_bundle[i] = STOP_CONDITION_COORD_RADIUS_EXCEEDED;
        } // END IF: photon escape radius exceeded
        if (d_status_bundle[i] == ACTIVE &&
            (!isfinite(coordinate_time) ||
             coordinate_time < commondata->slot_manager_t_min ||
             coordinate_time >= commondata->t_start + 1.0e-5)) {
            d_status_bundle[i] = STOP_CONDITION_T_MAX_EXCEEDED;
        } // END IF: photon coordinate time invalid
    } // END LOOP: i processes accepted photon steps
    #undef IDX_F
""".replace("@COORDINATE_TIME@", coordinate_time).replace(
        "@HORIZON_CHECK@", horizon_check
    )
    if analytical:
        body = body.replace(
            "(!isfinite(coordinate_time) ||\n"
            "             coordinate_time < commondata->slot_manager_t_min ||\n"
            "             coordinate_time >= commondata->t_start + 1.0e-5)",
            "!isfinite(coordinate_time)",
        ).replace("STOP_CONDITION_T_MAX_EXCEEDED;", "FAILURE_GENERIC;")

    includes = ["BHaH_defines.h"]
    prefunc = ""
    if analytical and par.parval_from_str("parallelization") == "cuda":
        kernel_body = body.replace(
            "(void)stream_idx;\n    #pragma omp parallel for\n"
            "    for (long int i = 0; i < chunk_size; ++i) {",
            "const long int i = blockIdx.x * blockDim.x + threadIdx.x;\n"
            "    if (i >= chunk_size) return;",
        ).replace("continue;", "return;")
        kernel_body = kernel_body.replace(
            "    } // END LOOP: i processes accepted photon steps\n    #undef IDX_F",
            "    #undef IDX_F",
        ).replace("commondata->", "d_commondata.")
        cuda_args = {
            "d_f_bundle": "const double *restrict",
            "d_f_start_bundle": "const double *restrict",
            "d_log_energy_bundle": "const double *restrict",
            "d_integration_param_before": "const double *restrict",
            "d_integration_param_after": "const double *restrict",
            "d_status_bundle": "termination_type_t *restrict",
            "d_on_pos_non_terminal_plane_prev": "bool *restrict",
            "d_on_pos_terminal_plane_prev": "bool *restrict",
            "plane_crossing_candidates": "photon_plane_crossing_candidates_t *restrict",
            "d_chunk_buffer": "const long int *restrict",
            "chunk_size": "const int",
        }
        prefunc, body = parallel_utils.generate_kernel_and_launch_code(
            kernel_name="event_detection_manager_kernel",
            kernel_body=kernel_body,
            arg_dict_cuda=cuda_args,
            arg_dict_host=cuda_args,
            parallelization="cuda",
            launch_dict={
                "threads_per_block": ["256", "1", "1"],
                "blocks_per_grid": ["(chunk_size + 255) / 256", "1", "1"],
                "stream": "stream_idx",
            },
            cfunc_decorators="__global__",
        )
        includes.extend(["BHaH_device_defines.h", "cuda_intrinsics.h"])

    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc="""Record accepted crossing steps and apply photon stop conditions.

        State arrays use nine components with stride BUNDLE_CAPACITY.
        Integration parameters are affine parameters for direct evolution and
        coordinate times for normalized evolution. Rejected steps are skipped.
        CUDA uses the common parameters already stored in d_commondata.

        @param[in] commondata Plane geometry, stop thresholds, and numerical time limits; unused by the CUDA host launcher.
        @param[in] d_f_bundle Accepted photon states in component-major chunk storage.
        @param[in] d_f_start_bundle Photon states at the start of each accepted RKF45 step.
        @param[in] d_log_energy_bundle Logarithmic photon energy measures at the accepted states.
        @param[in] d_integration_param_before Integration parameters before each accepted step.
        @param[in] d_integration_param_after Integration parameters after each accepted step.
        @param[in,out] d_status_bundle Chunk termination statuses; only ACTIVE entries are processed.
        @param[in,out] d_on_pos_non_terminal_plane_prev Saved nonterminal plane sides, updated when the plane is checked.
        @param[in,out] d_on_pos_terminal_plane_prev Saved terminal plane sides, updated when the plane is checked.
        @param[in,out] plane_crossing_candidates Global photon crossing records, updated with accepted crossing steps.
        @param[in] d_chunk_buffer Mapping from chunk entries to global photon indices.
        @param chunk_size Number of photon entries in the chunk.
        @param stream_idx CUDA launch stream index; unused by OpenMP.
        """ + horizon_desc,
        cfunc_type="void",
        name="event_detection_manager_kernel",
        params=(
            "const commondata_struct *restrict commondata, "
            "const double *restrict d_f_bundle, "
            "const double *restrict d_f_start_bundle, "
            "const double *restrict d_log_energy_bundle, "
            "const double *restrict d_integration_param_before, "
            "const double *restrict d_integration_param_after, "
            "termination_type_t *restrict d_status_bundle, "
            "bool *restrict d_on_pos_non_terminal_plane_prev, "
            "bool *restrict d_on_pos_terminal_plane_prev, "
            "photon_plane_crossing_candidates_t *restrict plane_crossing_candidates, "
            "const long int *restrict d_chunk_buffer, "
            "const int chunk_size, "
            "const int stream_idx" + horizon_params
        ),
        include_CodeParameters_h=False,
        body=body,
    )
