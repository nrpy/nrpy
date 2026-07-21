# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/event_detection_manager_kernel.py
"""
Provides the C orchestrator for geometric event detection.

This module provides the high-level logic for detecting crossings of the observer
window and the source emission plane. It generates a C kernel that reads the current
and historical integration state bundles from global device memory into local arrays
to evaluate energy limits and coordinate-radius bounds before verifying
physical plane intersections. The geometric boundaries remain mathematically immutable
across all rendered tiles, ensuring consistent hit detection depth. The kernel calls
downstream interpolation routines to resolve precise boundary crossing coordinates
and outputs the filtered physical intersections to persistent blueprint structures
while shifting the valid trajectory history arrays to stage the next solver step.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par


def event_detection_manager_kernel(normalized_eom: bool = False) -> None:
    """
    Define the configuration and parameters for the event-detection kernel.

    :param normalized_eom: Whether the state stores affine parameter in ``f[0]``
        and coordinate time in the integration-parameter tracker.
    """
    par.register_CodeParameters(
        "REAL",
        __name__,
        ["energy_max"],
        [1e5],
        commondata=True,
        add_to_parfile=True,
    )

    parallelization = par.parval_from_str("parallelization")
    cd_access = parallel_utils.get_commondata_access(parallelization)

    find_event_c_code = cfc.CFunction_dict["find_event_time_and_state"].full_function
    window_c_code = cfc.CFunction_dict["handle_window_plane_intersection"].full_function
    source_c_code = cfc.CFunction_dict["handle_source_plane_intersection"].full_function

    arg_dict_cuda = {
        "d_f_bundle": "const double *restrict",
        "d_f_prev_bundle": "double *restrict",
        "d_f_pre_prev_bundle": "double *restrict",
        "d_integration_param": "const double *restrict",
        "d_integration_param_prev": "double *restrict",
        "d_integration_param_pre_prev": "double *restrict",
        "d_results_buffer": "blueprint_data_t *restrict",
        "d_status_bundle": "termination_type_t *restrict",
        "d_on_pos_window_prev": "bool *restrict",
        "d_on_pos_source_prev": "bool *restrict",
        "d_window_event_found": "bool *restrict",
        "d_source_event_found": "bool *restrict",
        "d_chunk_buffer": "const long int *restrict",
        "chunk_size": "const int",
    }

    arg_dict_host = {
        "d_f_bundle": "const double *restrict",
        "d_f_prev_bundle": "double *restrict",
        "d_f_pre_prev_bundle": "double *restrict",
        "d_integration_param": "const double *restrict",
        "d_integration_param_prev": "double *restrict",
        "d_integration_param_pre_prev": "double *restrict",
        "d_results_buffer": "blueprint_data_t *restrict",
        "d_status_bundle": "termination_type_t *restrict",
        "d_on_pos_window_prev": "bool *restrict",
        "d_on_pos_source_prev": "bool *restrict",
        "d_window_event_found": "bool *restrict",
        "d_source_event_found": "bool *restrict",
        "d_chunk_buffer": "const long int *restrict",
        "chunk_size": "const int",
    }

    # Pass commondata explicitly when not using CUDA's global memory
    if parallelization != "cuda":
        arg_dict_cuda["commondata"] = "const commondata_struct *restrict"
        arg_dict_host["commondata"] = "const commondata_struct *restrict"

    escape_statement = "return;" if parallelization == "cuda" else "continue;"

    # Variables to handle architecture differences dynamically
    commondata_arg = "" if parallelization == "cuda" else ", commondata"
    pragma_unroll = "#pragma unroll" if parallelization == "cuda" else ""
    intersection_state_setup = (
        """
            double f_intersection[9];
            for (int component = 0; component < 9; ++component) {
                f_intersection[component] = f_int[component];
            } // END LOOP: for component over event-state components
            f_intersection[0] = event_integration_param;
            const double physical_lambda = f_int[0];
        """
        if normalized_eom
        else "const double physical_lambda = event_integration_param;"
    )
    intersection_state_name = "f_intersection" if normalized_eom else "f_int"

    if parallelization == "cuda":
        loop_preamble = """
    //==========================================
    // CUDA THREAD IDENTIFICATION
    //==========================================
    // Thread ID maps to a unique photon index $i$ within the execution chunk.
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= chunk_size) return;
    """
        loop_postamble = ""
    else:
        loop_preamble = """
    //==========================================
    // OPENMP LOOP ARCHITECTURE
    //==========================================
    // Distribute photon rays across available CPU threads for parallel evaluation.
    #pragma omp parallel for
    for(long int i = 0; i < chunk_size; i++) {
    """
        loop_postamble = "    } // END LOOP: for i over chunk_size rays"

    core_math = rf"""
    // Resolves the absolute global memory index $m_{{idx}}$ of the trajectory to bypass local array overwriting.
    const long int master_idx = d_chunk_buffer[i];

    //==========================================
    // MACROS
    //==========================================
    // IDX_F maps a component to the flattened state bundle using SoA layout.
    // Memory Striding strictly uses BUNDLE_CAPACITY, preventing bounds failure on remainders.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))

    //==========================================
    // MODE-SPECIFIC ENERGY LIMIT CHECK
    //==========================================
    // Direct evolution stores $p^0$ in f[4]; normalized evolution stores its log-energy measure.
    const double energy_measure = ReadCUDA(&d_f_bundle[IDX_F(4, i)]);

    if (AbsCUDA(energy_measure) > {cd_access}energy_max) {{
        d_status_bundle[i] = FAILURE_ENERGY_LIMIT_EXCEEDED; // Stops a ray whose energy measure exceeded its limit.
        {escape_statement}
    }} // END IF: mode-specific energy limit exceeded

    // Terminated photons cleanly bypass the geometric evaluation logic.
    if (d_status_bundle[i] != ACTIVE) {escape_statement}

    //==========================================
    // LOCAL REGISTER HYDRATION
    //==========================================
    // 1-Pass read from global memory bundles into thread-local arrays to respect hardware register limits.
    double f_local[9], f_p_local[9], f_p_p_local[9]; // Thread-local arrays holding the state vector $f^\mu$ and its history.
    {pragma_unroll}
    for (int c = 0; c < 9; c++) {{ // Loops over the $9$ tensor components.
        f_local[c] = ReadCUDA(&d_f_bundle[IDX_F(c, i)]); // Hydrates the current step state $f^\mu_n$.
        f_p_local[c] = ReadCUDA(&d_f_prev_bundle[IDX_F(c, i)]); // Hydrates the previous step state $f^\mu_{{n-1}}$.
        f_p_p_local[c] = ReadCUDA(&d_f_pre_prev_bundle[IDX_F(c, i)]); // Hydrates the pre-previous step state $f^\mu_{{n-2}}$.
    }} // END LOOP: for c over 9 tensor components

    // Hydrate the event-interpolation parameter history.
    const double integration_param_local = ReadCUDA(&d_integration_param[i]);
    const double integration_param_prev_local =
        ReadCUDA(&d_integration_param_prev[i]);
    const double integration_param_pre_prev_local =
        ReadCUDA(&d_integration_param_pre_prev[i]);

    const double x = f_local[1]; // Extracts the local Cartesian coordinate $x$.
    const double y = f_local[2]; // Extracts the local Cartesian coordinate $y$.
    const double z = f_local[3]; // Extracts the local Cartesian coordinate $z$.

    //==========================================
    // CELESTIAL ESCAPE CHECK
    //==========================================
    // Evaluates if the photon has exceeded the coordinate escape radius $r_{{escape}}$.
    const double r_sq = x*x + y*y + z*z; // Computes the squared radial distance $r^2$ from the origin.
    if (r_sq > ({cd_access}r_escape * {cd_access}r_escape)) {{
        d_status_bundle[i] = TERMINATION_TYPE_COORD_RADIUS_EXCEEDED; // Marks the coordinate-radius limit termination.
        {escape_statement}
    }} // END IF: coordinate-radius escape limit exceeded

    //==========================================
    // EVENT DETECTION & TERMINATION CHECKS
    //==========================================
    // Evaluates physical plane intersections strictly using localized variables.

    // Window plane logic is guarded to lock the intersection coordinates permanently.
    if (!d_window_event_found[i]) {{
        //==========================================
        // GLOBAL WINDOW PLANE RECONSTRUCTION
        //==========================================

        double w_normal[3]; // 3D geometric unit normal vector $n_i$ of the global observer window.
        w_normal[0] = {cd_access}original_window_center_x - {cd_access}camera_pos_x; // Computes the $x$ component of the global window normal.
        w_normal[1] = {cd_access}original_window_center_y - {cd_access}camera_pos_y; // Computes the $y$ component of the global window normal.
        w_normal[2] = {cd_access}original_window_center_z - {cd_access}camera_pos_z; // Computes the $z$ component of the global window normal.

        const double mag_inv = 1.0 / SqrtCUDA(w_normal[0]*w_normal[0] + w_normal[1]*w_normal[1] + w_normal[2]*w_normal[2]); // Computes the inverse magnitude $1/|n|$ for normalization.
        w_normal[0] *= mag_inv; // Normalizes the $x$ component.
        w_normal[1] *= mag_inv; // Normalizes the $y$ component.
        w_normal[2] *= mag_inv; // Normalizes the $z$ component.

        // Calculates orthogonal distance $d_w$ from the origin to the global window plane.
        const double w_dist = {cd_access}original_window_center_x*w_normal[0] + {cd_access}original_window_center_y*w_normal[1] + {cd_access}original_window_center_z*w_normal[2];

        // Evaluates the global plane equation $E_w$ for the photon's current spatial position.
        const double w_val = x*w_normal[0] + y*w_normal[1] + z*w_normal[2] - w_dist;

        const bool on_pos_win_curr = (w_val > 1e-10); // Checks if the photon is on the positive side of the global window.
        const bool on_pos_win_prev = d_on_pos_window_prev[i]; // Retrieves the previous integration step's window side evaluation.

        if (on_pos_win_curr != on_pos_win_prev) {{ // Triggers intersection event if the physical plane was crossed.
            double f_int[9];  // Reconstructed $9$-component state vector $f^\mu$ at the intersection.
            double event_integration_param;
            find_event_time_and_state(f_local, f_p_local, f_p_p_local, integration_param_local, integration_param_prev_local, integration_param_pre_prev_local, w_normal, w_dist, &event_integration_param, f_int);
{intersection_state_setup}

            // Writes the physical intersection to the persistent master index array slot via global mapping.
            // The downstream function handle_window_plane_intersection safely handles mapping the global
            // spatial coordinates to the local tile offsets.
            // Window plane function call to pass commondata conditionally.
            if (handle_window_plane_intersection({intersection_state_name}, physical_lambda, &d_results_buffer[master_idx]{commondata_arg})) {{
                d_window_event_found[i] = true;
            }} // END IF: handle_window_plane_intersection succeeded
        }} // END IF: physical window plane was crossed
        d_on_pos_window_prev[i] = on_pos_win_curr; // Updates the window evaluation history for the next step.
    }} // END IF: window event not found

    // Source plane logic is guarded to lock the intersection coordinates permanently.
    if (!d_source_event_found[i]) {{
        const double s_normal[3] = {{{cd_access}source_plane_normal_x, {cd_access}source_plane_normal_y, {cd_access}source_plane_normal_z}}; // 3D geometric unit normal vector $n_i$ of the source plane.
        const double s_dist = {cd_access}source_plane_center_x*s_normal[0] + {cd_access}source_plane_center_y*s_normal[1] + {cd_access}source_plane_center_z*s_normal[2]; // Calculates orthogonal distance $d_s$ to the source plane.
        const double s_val = x*s_normal[0] + y*s_normal[1] + z*s_normal[2] - s_dist; // Evaluates the plane equation $E_s$ for the current position.

        const bool on_pos_src_curr = (s_val > 1e-10); // Checks if the photon is on the positive side of the source.
        const bool on_pos_src_prev = d_on_pos_source_prev[i]; // Retrieves the previous integration step's source side evaluation.

        if (on_pos_src_curr != on_pos_src_prev) {{ // Triggers intersection event if the plane was crossed.
            double f_int[9];  // Reconstructed $9$-component state vector $f^\mu$ at the intersection.
            double event_integration_param;
            find_event_time_and_state(f_local, f_p_local, f_p_p_local, integration_param_local, integration_param_prev_local, integration_param_pre_prev_local, s_normal, s_dist, &event_integration_param, f_int);
{intersection_state_setup}
            // Writes the physical intersection to the persistent master index array slot via global mapping.
            if (handle_source_plane_intersection({intersection_state_name}, physical_lambda, &d_results_buffer[master_idx]{commondata_arg})) {{
                d_status_bundle[i] = TERMINATION_TYPE_SOURCE_PLANE; // Marks the ray as terminated upon striking the source plane.
                d_source_event_found[i] = true; // Locks the source intersection to prevent future overwrites.
            }} // END IF: handle_source_plane_intersection succeeded
        }} // END IF: physical source plane was crossed
        d_on_pos_source_prev[i] = on_pos_src_curr; // Updates the source evaluation history for the next step.
    }} // END IF: source event not found

    //==========================================
    // HISTORY SHIFT
    //==========================================
    // Memory bundles are updated only for active trajectories to stage the next RKF45 step.
    if (d_status_bundle[i] == ACTIVE) {{
        WriteCUDA(
            &d_integration_param_pre_prev[i], integration_param_prev_local);
        WriteCUDA(&d_integration_param_prev[i], integration_param_local);
        {pragma_unroll}
        for (int c = 0; c < 9; c++) {{ // Loops over the $9$ tensor components to shift the history.
            WriteCUDA(&d_f_pre_prev_bundle[IDX_F(c, i)], f_p_local[c]); // Shifts the previous state $f^\mu_{{n-1}}$ to the pre-previous slot $f^\mu_{{n-2}}$.
            WriteCUDA(&d_f_prev_bundle[IDX_F(c, i)], f_local[c]); // Shifts the current state $f^\mu_{{n}}$ to the previous slot $f^\mu_{{n-1}}$.
        }} // END LOOP: for c over 9 tensor components to shift the history
    }} // END IF: trajectory is active

    #undef IDX_F
    """

    kernel_body = f"{loop_preamble}\n{core_math}\n{loop_postamble}"

    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
        "stream": "stream_idx",
    }

    prefunc_kernel, body = parallel_utils.generate_kernel_and_launch_code(
        kernel_name="event_detection_manager_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_host,
        parallelization=parallelization,
        launch_dict=launch_dict,
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
        thread_tiling_macro_suffix="RKF45",
    )

    prefunc = "\n\n".join(
        [find_event_c_code, window_c_code, source_c_code, prefunc_kernel]
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if parallelization == "cuda":
        includes.append("cuda_intrinsics.h")
        includes.append("BHaH_device_defines.h")

    desc = r""" Optimized detection of plane crossings using consolidated blueprints.

    @param d_f_bundle SoA pointer to the state array for step $f^\mu_{n}$.
    @param d_f_prev_bundle SoA pointer to the state array for step $f^\mu_{n-1}$.
    @param d_f_pre_prev_bundle SoA pointer to the state array for step $f^\mu_{n-2}$.
    @param d_integration_param Pointer to the current integration parameter.
    @param d_integration_param_prev Pointer to the preceding integration parameter.
    @param d_integration_param_pre_prev Pointer to the integration parameter two steps earlier.
    @param d_results_buffer Pointer to the flat array of blueprint data structures $b_i$.
    @param d_status_bundle Pointer to the array of termination statuses.
    @param d_on_pos_window_prev Array tracking the window plane side.
    @param d_on_pos_source_prev Array tracking the source plane side.
    @param d_window_event_found Array tracking if a window intersection has been locked.
    @param d_source_event_found Array tracking if a source intersection has been locked.
    @param d_chunk_buffer Array containing the absolute master mapping indices $m_{idx}$.
    @param chunk_size The number of active rays in the current bundle batch.
    @param stream_idx The hardware stream identifier."""
    cfunc_type = "void"
    name = "event_detection_manager_kernel"
    params = (
        "const commondata_struct *restrict commondata, "
        "const double *restrict d_f_bundle, "
        "double *restrict d_f_prev_bundle, "
        "double *restrict d_f_pre_prev_bundle, "
        "const double *restrict d_integration_param, "
        "double *restrict d_integration_param_prev, "
        "double *restrict d_integration_param_pre_prev, "
        "blueprint_data_t *restrict d_results_buffer, "
        "termination_type_t *restrict d_status_bundle, "
        "bool *restrict d_on_pos_window_prev, "
        "bool *restrict d_on_pos_source_prev, "
        "bool *restrict d_window_event_found, "
        "bool *restrict d_source_event_found, "
        "const long int *restrict d_chunk_buffer, "
        "const int chunk_size,"
        "const int stream_idx"
    )
    include_CodeParameters_h = False

    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
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
