"""
Generates the C orchestrator for geometric event detection.

This module provides the high-level logic for detecting crossings of the
observer window and the source emission plane. It operates strictly within 
the Split-Pipeline architecture, reading from VRAM scratchpad bundles and
incorporating absolute global mapping locks to prevent blueprint overwrites.
Author: Dalton J. Moone.
"""
import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code


def event_detection_manager_kernel() -> None:
    """
    Register the event_detection_manager global CUDA kernel.

    :raises SystemError: If C function registration fails during the pipeline compilation.
    """
    par.register_CodeParameters(
        "REAL",
        __name__,
        ["p_t_max"],
        [1e5],
        commondata=True,
        add_to_parfile=True,
    )

    find_event_c_code = cfc.CFunction_dict["find_event_time_and_state"].full_function
    window_c_code = cfc.CFunction_dict["handle_window_plane_intersection"].full_function
    source_c_code = cfc.CFunction_dict["handle_source_plane_intersection"].full_function

    kernel_body = r"""
    // --- THREAD IDENTIFICATION & BOUNDS ---
    // Thread ID maps to a unique photon index within the execution chunk.
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= chunk_size) return;

    // Resolves the absolute global memory index of the trajectory to bypass local array overwriting.
    const long int master_idx = d_chunk_buffer[i];

    // --- MACROS ---
    // IDX_F maps a component to the flattened state bundle using SoA layout.
    // Memory Striding strictly uses BUNDLE_CAPACITY, preventing bounds failure on remainders.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))

    // --- TEMPORAL EXPLOSION CHECK ---
    // Reads $p_t$ directly from VRAM to terminate doomed rays before register hydration.
    const double p_t = ReadCUDA(&d_f_bundle[IDX_F(4, i)]);

    if (AbsCUDA(p_t) > d_commondata.p_t_max) {
        d_status_bundle[i] = FAILURE_PT_TOO_BIG;
        return;
    }

    // Terminated photons cleanly bypass the geometric evaluation logic.
    if (d_status_bundle[i] != ACTIVE) return;

    // --- LOCAL REGISTER HYDRATION ---
    // 1-Pass read from global VRAM bundles into thread-local arrays to respect sm_86 limits.
    double f_local[9], f_p_local[9], f_p_p_local[9];
    #pragma unroll
    for (int c = 0; c < 9; c++) {
        f_local[c] = ReadCUDA(&d_f_bundle[IDX_F(c, i)]);
        f_p_local[c] = ReadCUDA(&d_f_prev_bundle[IDX_F(c, i)]);
        f_p_p_local[c] = ReadCUDA(&d_f_pre_prev_bundle[IDX_F(c, i)]);
    }

    // Hydrates the affine parameter $\lambda$ dependencies directly from the discrete VRAM trackers.
    const double lam_local = ReadCUDA(&d_affine[i]);
    const double lam_p_local = ReadCUDA(&d_affine_prev[i]);
    const double lam_p_p_local = ReadCUDA(&d_affine_pre_prev[i]);

    const double x = f_local[1]; 
    const double y = f_local[2]; 
    const double z = f_local[3];

    // --- CELESTIAL ESCAPE CHECK ---
    // Evaluates if the photon has exceeded the coordinate escape radius $r_{escape}$.
    const double r_sq = x*x + y*y + z*z;
    if (r_sq > (d_commondata.r_escape * d_commondata.r_escape)) {
        d_status_bundle[i] = TERMINATION_TYPE_CELESTIAL_SPHERE;
        return;
    }

    // --- EVENT DETECTION & TERMINATION CHECKS ---
    // Evaluates physical plane intersections strictly using localized variables.
    
    // Window plane logic is guarded to lock the intersection coordinates permanently.
    if (!d_window_event_found[i]) {
        double w_normal[3]; // 3D geometric unit normal vector $n_i$ of the observer window.
        w_normal[0] = d_commondata.window_center_x - d_commondata.camera_pos_x;
        w_normal[1] = d_commondata.window_center_y - d_commondata.camera_pos_y;
        w_normal[2] = d_commondata.window_center_z - d_commondata.camera_pos_z;
        
        const double mag_inv = 1.0 / SqrtCUDA(w_normal[0]*w_normal[0] + w_normal[1]*w_normal[1] + w_normal[2]*w_normal[2]);
        w_normal[0] *= mag_inv; w_normal[1] *= mag_inv; w_normal[2] *= mag_inv;
        
        const double w_dist = d_commondata.window_center_x*w_normal[0] + d_commondata.window_center_y*w_normal[1] + d_commondata.window_center_z*w_normal[2];
        const double w_val = x*w_normal[0] + y*w_normal[1] + z*w_normal[2] - w_dist;
        
        const bool on_pos_win_curr = (w_val > 1e-10);
        const bool on_pos_win_prev = d_on_pos_window_prev[i];

        if (on_pos_win_curr != on_pos_win_prev) {
            double f_int[9];  // Reconstructed 9-component state vector $f^\mu$ at the intersection.
            double lam_event; // Interpolated affine parameter $\lambda$ of the exact boundary crossing.
            find_event_time_and_state(f_local, f_p_local, f_p_p_local, lam_local, lam_p_local, lam_p_p_local, w_normal, w_dist, &lam_event, f_int);
            // Writes the physical intersection to the persistent master index array slot via global mapping.
            if (handle_window_plane_intersection(f_int, lam_event, &d_results_buffer[master_idx])) {
                d_window_event_found[i] = true;
            }
        }
        d_on_pos_window_prev[i] = on_pos_win_curr;
    }

    // Source plane logic is guarded to lock the intersection coordinates permanently.
    if (!d_source_event_found[i]) {
        const double s_normal[3] = {d_commondata.source_plane_normal_x, d_commondata.source_plane_normal_y, d_commondata.source_plane_normal_z};
        const double s_dist = d_commondata.source_plane_center_x*s_normal[0] + d_commondata.source_plane_center_y*s_normal[1] + d_commondata.source_plane_center_z*s_normal[2];
        const double s_val = x*s_normal[0] + y*s_normal[1] + z*s_normal[2] - s_dist;
        
        const bool on_pos_src_curr = (s_val > 1e-10);
        const bool on_pos_src_prev = d_on_pos_source_prev[i];

        if (on_pos_src_curr != on_pos_src_prev) {
            double f_int[9];  // Reconstructed 9-component state vector $f^\mu$ at the intersection.
            double lam_event; // Interpolated affine parameter $\lambda$ of the exact boundary crossing.
            find_event_time_and_state(f_local, f_p_local, f_p_p_local, lam_local, lam_p_local, lam_p_p_local, s_normal, s_dist, &lam_event, f_int);
            // Writes the physical intersection to the persistent master index array slot via global mapping.
            if (handle_source_plane_intersection(f_int, lam_event, &d_results_buffer[master_idx])) {
                d_status_bundle[i] = TERMINATION_TYPE_SOURCE_PLANE;
                d_source_event_found[i] = true;
            }
        }
        d_on_pos_source_prev[i] = on_pos_src_curr;
    }

    // --- HISTORY SHIFT ---
    // VRAM bundles are updated only for active trajectories to stage the next RKF45 step.
    if (d_status_bundle[i] == ACTIVE) {
        // Shifts the affine parameter tracking variables forward one discrete step.
        WriteCUDA(&d_affine_pre_prev[i], lam_p_local);
        WriteCUDA(&d_affine_prev[i], lam_local);
        #pragma unroll
        for (int c = 0; c < 9; c++) {
            WriteCUDA(&d_f_pre_prev_bundle[IDX_F(c, i)], f_p_local[c]);
            WriteCUDA(&d_f_prev_bundle[IDX_F(c, i)], f_local[c]);
        }
    }

    #undef IDX_F
    """

    arg_dict = {
        "d_f_bundle": "const double *restrict",
        "d_f_prev_bundle": "double *restrict",
        "d_f_pre_prev_bundle": "double *restrict",
        "d_affine": "const double *restrict",
        "d_affine_prev": "double *restrict",
        "d_affine_pre_prev": "double *restrict",
        "d_results_buffer": "blueprint_data_t *restrict",
        "d_status_bundle": "termination_type_t *restrict",
        "d_on_pos_window_prev": "bool *restrict",
        "d_on_pos_source_prev": "bool *restrict",
        "d_window_event_found": "bool *restrict",
        "d_source_event_found": "bool *restrict",
        "d_chunk_buffer": "const long int *restrict",
        "chunk_size": "const int"
    }

    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
    }

    prefunc_kernel, body = generate_kernel_and_launch_code(
        kernel_name="event_detection_manager_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict,
        arg_dict_host=arg_dict,
        parallelization="cuda",
        launch_dict=launch_dict,
        cfunc_decorators="__global__",
        thread_tiling_macro_suffix="RKF45"
    )

    # Use a safe join operation to prevent potential concatenation issues
    prefunc = "\n\n".join([find_event_c_code, window_c_code, source_c_code, prefunc_kernel])

    includes = [
        "BHaH_defines.h", 
        "BHaH_device_defines.h", 
        "BHaH_function_prototypes.h", 
        "<math.h>", 
        "<stdbool.h>",
        "cuda_intrinsics.h"
    ]
    
    desc = r"""@brief GPU-optimized detection of plane crossings using consolidated blueprints.

    @param d_f_bundle VRAM SoA pointer to the state array for step $f^\mu_{n}$.
    @param d_f_prev_bundle VRAM SoA pointer to the state array for step $f^\mu_{n-1}$.
    @param d_f_pre_prev_bundle VRAM SoA pointer to the state array for step $f^\mu_{n-2}$.
    @param d_affine VRAM pointer to the current affine parameter $\lambda_n$.
    @param d_affine_prev VRAM pointer to the history affine parameter $\lambda_{n-1}$.
    @param d_affine_pre_prev VRAM pointer to the history affine parameter $\lambda_{n-2}$.
    @param d_results_buffer Pointer to the flat array of blueprint data structures $b_i$.
    @param d_status_bundle Pointer to the array of termination statuses.
    @param d_on_pos_window_prev VRAM array tracking the window plane side.
    @param d_on_pos_source_prev VRAM array tracking the source plane side.
    @param d_window_event_found VRAM array tracking if a window intersection has been locked.
    @param d_source_event_found VRAM array tracking if a source intersection has been locked.
    @param d_chunk_buffer VRAM array containing the absolute master mapping indices $m_{idx}$.
    @param chunk_size The number of active rays in the current bundle batch."""
    
    cfunc_type = "void"
    name = "event_detection_manager_kernel"
    
    params = (
        "const double *restrict d_f_bundle, "
        "double *restrict d_f_prev_bundle, "
        "double *restrict d_f_pre_prev_bundle, "
        "const double *restrict d_affine, "
        "double *restrict d_affine_prev, "
        "double *restrict d_affine_pre_prev, "
        "blueprint_data_t *restrict d_results_buffer, "
        "termination_type_t *restrict d_status_bundle, "
        "bool *restrict d_on_pos_window_prev, "
        "bool *restrict d_on_pos_source_prev, "
        "bool *restrict d_window_event_found, "
        "bool *restrict d_source_event_found, "
        "const long int *restrict d_chunk_buffer, "
        "const int chunk_size"
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
        body=body
    )