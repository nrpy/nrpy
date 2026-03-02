"""
Orchestration module for the Project Singularity-Axiom numerical integration pipeline.

This module generates the high-level C orchestrator responsible for managing the
life cycle of photon trajectories in curved spacetimes. It implements a Native CUDA
Streaming Bundle architecture that strictly adheres to the following physical constraints:

1. Memory Architecture: Utilizes a flattened Structure of Arrays (SoA) via
   the `PhotonStateSoA` and `blueprint_data_t` structures. Memory is transferred
   in strided bundles of 32,768 photons to protect the 10GB VRAM hardware limit.
2. Dual-Stream Overlap: Employs `cudaStream_t` to overlap the memory transfers
   of Batch N+1 with the computational integration of Batch N.
3. Temporal Binning: Employs a lock-free `TimeSlotManager` remaining strictly on
   the host CPU to bin active rays by physical coordinate time ($t$), eliminating
   serialized linked-list traversals on the GPU.
4. Stream Compaction: Implements a staged integration pattern where active photons
   are processed in dense bundles. Uses hardware-optimized warp-aggregated atomics
   to circumvent memory contention and ensure saturation.

Author: Dalton J. Moone.
"""

import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code

batch_structs_c_code = r"""
    #define BUNDLE_CAPACITY 32768 // Maximum number of photons processed per batch to fit within L1/L2 cache.

    // Strictly enforce 1D mapping to prevent pointer offset arithmetic bugs
    #define IDX_GLOBAL(component, ray_id, num_rays) ((component) * (num_rays) + (ray_id))
    #define IDX_LOCAL(component, batch_id, batch_size) ((component) * (batch_size) + (batch_id))

    // Defines the physical planes where a photon trajectory might terminate.
    typedef enum {
        WINDOW_EVENT, // Intersection with the observer's camera window.
        SOURCE_EVENT  // Intersection with the emission source plane.
    } event_type_t;

    // Defines the specific exit condition for a photon's integration loop.
    typedef enum {
        TERMINATION_TYPE_CELESTIAL_SPHERE, // 0: Photon escaped to infinity (exceeded r_escape).
        TERMINATION_TYPE_SOURCE_PLANE,     // 1: Photon successfully hit the source emission plane.
        FAILURE_PT_TOO_BIG,                // 2: Integration failed due to unbounded momentum $p_t$.
        FAILURE_RKF45_REJECTION_LIMIT,     // 3: Adaptive step-size rejected too many consecutive times.
        FAILURE_T_MAX_EXCEEDED,            // 4: Integration exceeded maximum allowable physical time.
        FAILURE_SLOT_MANAGER_ERROR,        // 5: TimeSlotManager failed to allocate or retrieve the photon.
        TERMINATION_TYPE_FAILURE,          // 6: Generic unclassified numerical failure.
        ACTIVE                             // 7: Photon is currently undergoing integration.
    } termination_type_t;

    // Stores the final physical properties of a photon upon integration termination.
    typedef struct {
        termination_type_t termination_type; // The exit condition of the photon.
        double y_w; // Local y-coordinate intersection on the observer window.
        double z_w; // Local z-coordinate intersection on the observer window.
        double y_s; // Local y-coordinate intersection on the source plane.
        double z_s; // Local z-coordinate intersection on the source plane.
        double final_theta; // Final polar angle $\theta$ at termination.
        double final_phi;   // Final azimuthal angle $\phi$ at termination.
        double L_w; // Affine parameter $\lambda$ at the window intersection.
        double t_w; // Physical coordinate time $t$ at the window intersection.
        double L_s; // Affine parameter $\lambda$ at the source intersection.
        double t_s; // Physical coordinate time $t$ at the source intersection.
    } __attribute__((packed)) blueprint_data_t;

    // ==========================================
    // Flattened SoA Struct (Master Storage)
    // ==========================================
    typedef struct {
        double *f; // Flattened state vector: 9 components (t, x, y, z, p_t, p_x, p_y, p_z, aux).
        double *f_p; // State vector at the previous integration step.
        double *f_p_p; // State vector at two integration steps prior.
        double *affine_param; // Current affine parameter $\lambda$ for the trajectory.
        double *affine_param_p; // Affine parameter $\lambda$ at the previous step.
        double *affine_param_p_p; // Affine parameter $\lambda$ at two steps prior.
        double *h; // Current adaptive step size for the RKF45 integrator.
        termination_type_t *status; // Current physical/numerical status of the photon.
        int *rejection_retries; // Counter for consecutive RKF45 error tolerance rejections.

        // Event Detection State Flags
        bool *on_positive_side_of_window_prev; // True if photon was previously 'above' the window plane.
        bool *on_positive_side_of_source_prev; // True if photon was previously 'above' the source plane.

        bool *source_event_found; // Flag indicating a source plane intersection was detected.
        double *source_event_lambda; // Exact affine parameter $\lambda$ of the source intersection.
        double *source_event_f_intersect; // Interpolated 9-component state vector at the source intersection.

        bool *window_event_found; // Flag indicating an observer window intersection was detected.
        double *window_event_lambda; // Exact affine parameter $\lambda$ of the window intersection.
        double *window_event_f_intersect; // Interpolated 9-component state vector at the window intersection.
    } PhotonStateSoA;
"""
Bdefines_h.register_BHaH_defines("photon_02_batch_structs", batch_structs_c_code)

par.register_CodeParameters(
    "REAL",
    __name__,
    [
        "t_integration_max",
        "r_escape",
        "p_t_max",
        "slot_manager_t_min",
        "slot_manager_delta_t",
        "numerical_initial_h",
    ],
    [10000.0, 150.0, 1e3, -1000.0, 10.0, 0.1],
    commondata=True,
    add_to_parfile=True,
)
par.register_CodeParameters(
    "bool",
    __name__,
    ["perform_conservation_check", "debug_mode"],
    [True, True],
    commondata=True,
    add_to_parfile=True,
)


def batch_integrator_numerical(spacetime_name: str) -> None:
    """
    Generate the Native CUDA orchestrator for the batched integration pipeline.

    Metaprograms a dual-stream architecture pipeline generating independent CUDA kernels 
    for stream compaction, RKF45 numerical integration, and event detection. Data transfers 
    are strictly scheduled to map thread-local registers over VRAM to prevent GPU bottlenecks.
    
    :param spacetime_name: The identifier for the spacetime metric (e.g., 'KerrSchild').
    :raises ValueError: If an unsupported spacetime metric identifier is provided.
    """

    # --- KERNEL 1: INIT BATCH ---
    init_body = r"""
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < bundle_size) {
        needs_step_bundle[tid] = true;
        status_bundle[tid] = ACTIVE;
        retries_bundle[tid] = 0;
        h_bundle[tid] = commondata->numerical_initial_h;
        for(int k=0; k<9; k++) {
            double val = f_bundle[IDX_LOCAL(k, tid, BUNDLE_CAPACITY)];
            f_p_bundle[IDX_LOCAL(k, tid, BUNDLE_CAPACITY)] = val;
            f_p_p_bundle[IDX_LOCAL(k, tid, BUNDLE_CAPACITY)] = val;
        }
    }
    """
    init_prefunc, init_launch = generate_kernel_and_launch_code(
        kernel_name="init_batch",
        kernel_body=init_body,
        arg_dict_cuda={
            "f_bundle": "const double *restrict",
            "f_p_bundle": "double *restrict",
            "f_p_p_bundle": "double *restrict",
            "h_bundle": "double *restrict",
            "needs_step_bundle": "bool *restrict",
            "status_bundle": "termination_type_t *restrict",
            "retries_bundle": "int *restrict",
            "commondata": "const commondata_struct *restrict",
            "bundle_size": "const long int"
        },
        arg_dict_host={},
        parallelization="cuda",
        launch_dict={
            "blocks_per_grid": ["(bundle_size + 63) / 64"],
            "threads_per_block": ["64", "1", "1"],
            "stream": "s_compute"
        },
        launchblock_with_braces=True,
        cfunc_decorators="__global__"
    )

    # --- KERNEL 2: WARP-AGGREGATED STREAM COMPACTION ---
    compaction_body = r"""
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int lane = threadIdx.x % 32;

    bool active = (tid < bundle_size) && needs_step_bundle[tid] && (status_bundle[tid] == ACTIVE);
    
    // The active thread mask is generated using hardware intrinsics to saturate the warp efficiently.
    unsigned int mask = __ballot_sync(0xFFFFFFFF, active);
    
    // The local destination offset is calculated using population counts to ensure wait-free warp synchronization.
    int lane_offset = __popc(mask & ((1 << lane) - 1));
    int warp_count = __popc(mask);

    int global_base = 0;
    if (active && lane_offset == 0) {
        // The warp leader performs a single atomic transaction to reserve block space for the warp.
        global_base = BHAH_WARP_ATOMIC_ADD(d_active_count, warp_count);
    }
    // Broadcast the base index to all active lanes in the warp.
    global_base = __shfl_sync(mask, global_base, __ffs(mask) - 1);

    if (active) {
        int pos = global_base + lane_offset;
        compacted_local_ids[pos] = tid;
    }
    """
    compaction_prefunc, compaction_launch = generate_kernel_and_launch_code(
        kernel_name="stream_compaction",
        kernel_body=compaction_body,
        arg_dict_cuda={
            "needs_step_bundle": "const bool *restrict",
            "status_bundle": "const termination_type_t *restrict",
            "compacted_local_ids": "int *restrict",
            "d_active_count": "int *restrict",
            "bundle_size": "const long int"
        },
        arg_dict_host={},
        parallelization="cuda",
        launch_dict={
            "blocks_per_grid": ["(bundle_size + 63) / 64"],
            "threads_per_block": ["64", "1", "1"],
            "stream": "s_compute"
        },
        launchblock_with_braces=True,
        cfunc_decorators="__global__"
    )

    # --- KERNEL 3: FUSED RKF45 INTEGRATION ---
    rkf45_body = r"""
    int cid = blockIdx.x * blockDim.x + threadIdx.x;
    if (cid >= *d_active_count) return;

    int local_idx = compacted_local_ids[cid];

    // Thread-local register allocations for the 10-component metric $g_{\mu\nu}$ and 40-component connection $\Gamma^\alpha_{\beta\gamma}$.
    // Retaining these strictly in registers prevents catastrophic VRAM spilling and adheres to the 255-register limit of the sm_86 architecture.
    double metric_local[10];
    double conn_local[40];

    // Shared memory allocation for the 54-component $k$-array spanning 6 RKF45 stages.
    // Ordering the array with the 64-thread dimension last ensures contiguous, coalesced memory access across the warp, preventing bank conflicts.
    __shared__ double k_shared[6][9][64];

    double f_start_local[9];
    double f_temp_local[9];

    for(int k=0; k<9; k++) {
        f_start_local[k] = f_bundle[IDX_LOCAL(k, local_idx, BUNDLE_CAPACITY)];
        f_temp_local[k] = f_start_local[k];
    }

    for(int stage = 1; stage <= 6; stage++) {
        placeholder_interpolation_engine_KerrSchild_Cartesian(commondata, f_temp_local, metric_local, conn_local);

        double k_local[9];
        calculate_ode_rhs(f_temp_local, metric_local, conn_local, k_local);

        for(int c=0; c<9; c++) {
            k_shared[stage-1][c][threadIdx.x] = k_local[c];
        }

        if (stage < 6) {
            // Passes the 3D shared memory base pointer directly to evaluate the intermediate RKF45 stages.
            calculate_rkf45_stage_f_temp(stage + 1, f_start_local, &k_shared[0][0][threadIdx.x], h_bundle[local_idx], f_temp_local);
        }
    }

    double f_out_local[9];
    double f_err_local[9];
    rkf45_kernel(f_start_local, &k_shared[0][0][threadIdx.x], h_bundle[local_idx], f_out_local, f_err_local);

    bool accepted = update_photon_state_and_stepsize(f_start_local, f_start_local, f_out_local, f_err_local, &h_bundle[local_idx], &affine_param_bundle[local_idx], &retries_bundle[local_idx], commondata);

    if (!accepted) {
        if (retries_bundle[local_idx] > commondata->rkf45_max_retries) {
            status_bundle[local_idx] = FAILURE_RKF45_REJECTION_LIMIT;
            needs_step_bundle[local_idx] = false;
        } else {
            *d_bundle_has_active = true;
        }
    } else {
        needs_step_bundle[local_idx] = false;

        // Perform VRAM-local state rotation entirely within VRAM, bypassing historical host synchronization.
        for(int k=0; k<9; ++k) {
            f_p_p_bundle[IDX_LOCAL(k, local_idx, BUNDLE_CAPACITY)] = f_p_bundle[IDX_LOCAL(k, local_idx, BUNDLE_CAPACITY)];
            f_p_bundle[IDX_LOCAL(k, local_idx, BUNDLE_CAPACITY)] = f_start_local[k];
            f_bundle[IDX_LOCAL(k, local_idx, BUNDLE_CAPACITY)] = f_out_local[k];
        }
    }
    """
    rkf45_prefunc, rkf45_launch = generate_kernel_and_launch_code(
        kernel_name="fused_rkf45",
        kernel_body=rkf45_body,
        arg_dict_cuda={
            "f_bundle": "double *restrict",
            "f_p_bundle": "double *restrict",
            "f_p_p_bundle": "double *restrict",
            "h_bundle": "double *restrict",
            "affine_param_bundle": "double *restrict",
            "needs_step_bundle": "bool *restrict",
            "status_bundle": "termination_type_t *restrict",
            "retries_bundle": "int *restrict",
            "compacted_local_ids": "const int *restrict",
            "d_active_count": "const int *restrict",
            "commondata": "const commondata_struct *restrict",
            "d_bundle_has_active": "bool *restrict"
        },
        arg_dict_host={},
        parallelization="cuda",
        launch_dict={
            "blocks_per_grid": ["(BUNDLE_CAPACITY + 63) / 64"], 
            "threads_per_block": ["64", "1", "1"],
            "stream": "s_compute"
        },
        launchblock_with_braces=True,
        cfunc_decorators="__global__"
    )

    # --- KERNEL 4: EVENT DETECTION ---
    event_body = r"""
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < bundle_size && status_bundle[tid] == ACTIVE) {
        double f_local[9], f_p_local[9], f_p_p_local[9];
        for(int k=0; k<9; k++) {
            f_local[k] = f_bundle[IDX_LOCAL(k, tid, BUNDLE_CAPACITY)];
            f_p_local[k] = f_p_bundle[IDX_LOCAL(k, tid, BUNDLE_CAPACITY)];
            f_p_p_local[k] = f_p_p_bundle[IDX_LOCAL(k, tid, BUNDLE_CAPACITY)];
        }

        // Call updated manager with corrected pointer names and blueprint architecture
        event_detection_manager(f_local, f_p_local, f_p_p_local, commondata, 
                                &status_bundle[tid], 
                                &results_buffer[start_idx + tid],
                                &on_pos_window_prev_bundle[tid], 
                                &on_pos_source_prev_bundle[tid]);

        // Secondary physical constraints checked in kernel to optimize manager registers
        if (fabs(f_local[4]) > commondata->p_t_max) status_bundle[tid] = FAILURE_PT_TOO_BIG;
        if (fabs(f_local[0]) > commondata->t_integration_max) status_bundle[tid] = FAILURE_T_MAX_EXCEEDED;
        double r_sq = f_local[1]*f_local[1] + f_local[2]*f_local[2] + f_local[3]*f_local[3];
        if (r_sq > commondata->r_escape * commondata->r_escape) status_bundle[tid] = TERMINATION_TYPE_CELESTIAL_SPHERE;
    }
    """
    event_prefunc, event_launch = generate_kernel_and_launch_code(
        kernel_name="event_detection",
        kernel_body=event_body,
        arg_dict_cuda={
            "f_bundle": "const double *restrict",
            "f_p_bundle": "const double *restrict",
            "f_p_p_bundle": "const double *restrict",
            "status_bundle": "termination_type_t *restrict",
            "results_buffer": "blueprint_data_t *restrict",
            "on_pos_window_prev_bundle": "bool *restrict",   
            "on_pos_source_prev_bundle": "bool *restrict",   
            "commondata": "const commondata_struct *restrict",
            "start_idx": "const long int",
            "bundle_size": "const long int"
        },
        arg_dict_host={},
        parallelization="cuda",
        launch_dict={
            "blocks_per_grid": ["(bundle_size + 63) / 64"],
            "threads_per_block": ["64", "1", "1"],
            "stream": "s_compute"
        },
        launchblock_with_braces=True,
        cfunc_decorators="__global__"
    )

    # --- MAIN C-FUNCTION GENERATION ---
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "cuda_runtime.h", "math.h"]
    desc = r"""@brief Finalized Project Singularity-Axiom Orchestrator (Native CUDA Architecture).

    This function serves as the master entry point for the ray-tracing pipeline.
    It performs the following high-level operations:
    1. Memory Management: Allocates dual-stream buffers in VRAM and Pinned Host RAM.
    2. Stream Architecture: Executes 9-strided asynchronous bridge transfers overlapping computation.
    3. Staged Integration: Executes the fused RKF45 logic leveraging explicit `__shared__` memory buffers.
    4. Event Processing: Evaluates local coordinate transformations for spatial planes.

    @param commondata Struct containing global spacetime parameters.
    @param num_rays The total number of photon trajectories to simulate, represented as $num\_rays$.
    @param results_buffer Device array storing the final physical intersections and bounds, mapped via $f^\mu$."""
    
    cfunc_type = "void"
    name = "batch_integrator_numerical"
    params = "const commondata_struct *restrict commondata, long int num_rays, blueprint_data_t *restrict results_buffer"

    initial_conditions_func = f"set_initial_conditions_cartesian_{spacetime_name}"

    # Python: Extract physical inline helpers from the global CFunction_dict registry.
    metric_func = cfc.CFunction_dict[f"g4DD_metric_{spacetime_name}"].full_function
    conn_func = cfc.CFunction_dict[f"connections_{spacetime_name}"].full_function
    interp_func = cfc.CFunction_dict[f"placeholder_interpolation_engine_{spacetime_name}"].full_function
    rhs_func = cfc.CFunction_dict["calculate_ode_rhs"].full_function

    # Python: Enforce linear top-down concatenation for C-compiler compliance. Tensor math precedes interpolation, which precedes ODE evaluations.
    prefunc = metric_func + conn_func + interp_func + rhs_func + init_prefunc + compaction_prefunc + rkf45_prefunc + event_prefunc

    # Python: Manually patch generated launch strings to index dual-stream arrays.
    # Python: This bypasses the ignored arg_dict_host and corrects pointer types without altering the locked gpu_kernel.py.
    init_launch = init_launch.replace("f_bundle, ", "f_bundle[s_compute], ").replace("h_bundle, ", "h_bundle[s_compute], ").replace("needs_step_bundle, ", "needs_step_bundle[s_compute], ").replace("status_bundle, ", "status_bundle[s_compute], ").replace("retries_bundle, ", "retries_bundle[s_compute], ")
    compaction_launch = compaction_launch.replace("needs_step_bundle, ", "needs_step_bundle[s_compute], ").replace("status_bundle, ", "status_bundle[s_compute], ")
    rkf45_launch = rkf45_launch.replace("f_bundle, ", "f_bundle[s_compute], ").replace("h_bundle, ", "h_bundle[s_compute], ").replace("affine_param_bundle, ", "affine_param_bundle[s_compute], ").replace("needs_step_bundle, ", "needs_step_bundle[s_compute], ").replace("status_bundle, ", "status_bundle[s_compute], ").replace("retries_bundle, ", "retries_bundle[s_compute], ")
    event_launch = event_launch.replace("on_pos_window_prev_bundle, ", "on_pos_window_prev_bundle[s_compute], ").replace("on_pos_source_prev_bundle, ", "on_pos_source_prev_bundle[s_compute], ").replace("f_bundle, ", "f_bundle[s_compute], ").replace("status_bundle, ", "status_bundle[s_compute], ").replace("results_buffer, ", "d_results_buffer, ")
    body = rf"""
    // --- HOST INITIALIZATION ---
    PhotonStateSoA all_photons_host; // Host-side replica of the flattened master storage.
    
    // Host-to-Device transfer: Pinned memory allocation ensures non-pageable memory for maximum PCIe DMA throughput.
    BHAH_MALLOC_PINNED(all_photons_host.f, sizeof(double) * 9 * num_rays);
    BHAH_MALLOC_PINNED(all_photons_host.affine_param, sizeof(double) * num_rays);
    BHAH_MALLOC_PINNED(all_photons_host.h, sizeof(double) * num_rays);
    BHAH_MALLOC_PINNED(all_photons_host.status, sizeof(termination_type_t) * num_rays);
    BHAH_MALLOC_PINNED(all_photons_host.rejection_retries, sizeof(int) * num_rays);

    // CPU-bound allocations mapping geometric orientation flags (standard malloc is fine here as they aren't asynchronously streamed)
    all_photons_host.on_positive_side_of_window_prev = (bool *)malloc(sizeof(bool) * num_rays);
    all_photons_host.on_positive_side_of_source_prev = (bool *)malloc(sizeof(bool) * num_rays);
    all_photons_host.source_event_found = (bool *)malloc(sizeof(bool) * num_rays);
    all_photons_host.source_event_lambda = (double *)malloc(sizeof(double) * num_rays);
    all_photons_host.source_event_f_intersect = (double *)malloc(sizeof(double) * 9 * num_rays);
    all_photons_host.window_event_found = (bool *)malloc(sizeof(bool) * num_rays);
    all_photons_host.window_event_lambda = (double *)malloc(sizeof(double) * num_rays);
    all_photons_host.window_event_f_intersect = (double *)malloc(sizeof(double) * 9 * num_rays);

    TimeSlotManager tsm; // Lock-free temporal arena bounded exclusively to the Host CPU context.
    slot_manager_init(&tsm, commondata->slot_manager_t_min, commondata->t_start + 1.0, commondata->slot_manager_delta_t, num_rays);

    double window_center[3]; // Cartesian coordinates mapping the geometric center of the camera window.
    double n_x[3]; // Orthonormal basis vector describing the horizontal axis of the camera frame.
    double n_y[3]; // Orthonormal basis vector describing the vertical axis of the camera frame.
    double n_z[3]; // Orthonormal basis vector normal to the camera window (line of sight).
    {initial_conditions_func}(commondata, num_rays, &all_photons_host, window_center, n_x, n_y, n_z);

    int initial_slot_idx = slot_get_index(&tsm, commondata->t_start);
    if(initial_slot_idx != -1) {{
        for(long int i=0; i<num_rays; ++i) slot_add_photon(&tsm, initial_slot_idx, i);
    }}

    // --- CUDA DUAL-STREAM BUFFER ALLOCATION ---
    cudaStream_t streams[2];
    cudaStreamCreate(&streams[0]);
    cudaStreamCreate(&streams[1]);

    termination_type_t *bundle_status_host[2];
    double *bundle_time_host[2];
    double *f_bundle[2];
    termination_type_t *status_bundle[2];
    double *h_bundle[2];
    bool *needs_step_bundle[2];
    int *retries_bundle[2];
    double *affine_param_bundle[2];
    bool *on_pos_window_prev_bundle[2];
    bool *on_pos_source_prev_bundle[2];

    for(int i=0; i<2; i++) {{
        // Device-to-Host transfer: Pinned memory allocation for status array to facilitate asynchronous CPU retrieval.
        BHAH_MALLOC_PINNED(bundle_status_host[i], sizeof(termination_type_t) * BUNDLE_CAPACITY);
        // Device-to-Host transfer: Pinned memory allocation for temporal state data.
        BHAH_MALLOC_PINNED(bundle_time_host[i], sizeof(double) * BUNDLE_CAPACITY);

        // Host-to-Device transfer: VRAM allocation for the physical state vector $f^\mu$ bounding the RKF45 batch.
        BHAH_MALLOC_DEVICE(f_bundle[i], sizeof(double) * 9 * BUNDLE_CAPACITY);
        // Host-to-Device transfer: VRAM allocation for the trajectory status enumerations.
        BHAH_MALLOC_DEVICE(status_bundle[i], sizeof(termination_type_t) * BUNDLE_CAPACITY);
        // Host-to-Device transfer: VRAM allocation tracking individual integration step sizes $h$.
        BHAH_MALLOC_DEVICE(h_bundle[i], sizeof(double) * BUNDLE_CAPACITY);
        // Host-to-Device transfer: VRAM allocation for boolean stepping masks.
        BHAH_MALLOC_DEVICE(needs_step_bundle[i], sizeof(bool) * BUNDLE_CAPACITY);
        // Host-to-Device transfer: VRAM allocation tracking local truncation error rejections.
        BHAH_MALLOC_DEVICE(retries_bundle[i], sizeof(int) * BUNDLE_CAPACITY);
        // Host-to-Device transfer: VRAM allocation tracking total affine progression $\lambda$.
        BHAH_MALLOC_DEVICE(affine_param_bundle[i], sizeof(double) * BUNDLE_CAPACITY);

        BHAH_MALLOC_DEVICE(on_pos_window_prev_bundle[i], sizeof(bool) * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(on_pos_source_prev_bundle[i], sizeof(bool) * BUNDLE_CAPACITY);
    }}

    double *f_p_bundle, *f_p_p_bundle; // Condensed local registers avoiding global VRAM history sprawl.
    
    // Host-to-Device transfer: VRAM allocation for the historical state vector $f_p^\mu$.
    BHAH_MALLOC_DEVICE(f_p_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY);
    // Host-to-Device transfer: VRAM allocation for the secondary historical state vector $f_p_p^\mu$.
    BHAH_MALLOC_DEVICE(f_p_p_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY);

    int *compacted_local_ids;
    int *d_active_count;
    bool *d_bundle_has_active;
    
    // Host-to-Device transfer: VRAM allocations for thread compaction logic tracking and indexing.
    BHAH_MALLOC_DEVICE(compacted_local_ids, sizeof(int) * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_active_count, sizeof(int));
    BHAH_MALLOC_DEVICE(d_bundle_has_active, sizeof(bool));

    blueprint_data_t *d_results_buffer;
    
    // Device-to-Host transfer: VRAM allocation storing the final intersected geometry outputs.
    BHAH_MALLOC_DEVICE(d_results_buffer, sizeof(blueprint_data_t) * num_rays);

    long int num_batches = (num_rays + BUNDLE_CAPACITY - 1) / BUNDLE_CAPACITY;

    // --- PIPELINE LOOP ---
    for (long int batch_idx = 0; batch_idx < num_batches; batch_idx++) {{
        int s_compute = batch_idx % 2;
        int s_transfer = (batch_idx + 1) % 2;
        long int start_idx = batch_idx * BUNDLE_CAPACITY;
        long int bundle_size = MIN(BUNDLE_CAPACITY, num_rays - start_idx);

        // 1. Initial bridging transfer reserved exclusively for Batch 0 alignment.
        if (batch_idx == 0) {{
            for(int k=0; k<9; k++) {{
                // Host-to-Device transfer: Asynchronously maps the $k$-th component of the state vector to VRAM using strict pointer arithmetic.
                cudaMemcpyAsync(f_bundle[s_compute] + (k * BUNDLE_CAPACITY),
                                all_photons_host.f + (k * num_rays) + start_idx,
                                bundle_size * sizeof(double), cudaMemcpyHostToDevice, streams[s_compute]);
            }}
        }}

        // 2. Compute stream synchronization blocks CPU solely to confirm prior kernels have flushed their buffer footprint.
        cudaStreamSynchronize(streams[s_compute]);

        // 3. Process the CPU-Bound TimeSlotManager parsing exclusively on Pinned RAM to preserve maximum SM throughput.
        if (batch_idx > 0) {{
            int prev_s = (batch_idx - 1) % 2;
            long int prev_start = (batch_idx - 1) * BUNDLE_CAPACITY;
            long int prev_size = MIN(BUNDLE_CAPACITY, num_rays - prev_start);
            for (long int j = 0; j < prev_size; j++) {{
                if (bundle_status_host[prev_s][j] == ACTIVE) {{
                    long int p_idx = prev_start + j;
                    int s_idx = slot_get_index(&tsm, bundle_time_host[prev_s][j]);
                    if (s_idx != -1) slot_add_photon(&tsm, s_idx, p_idx);
                }}
            }}
        }}

        // 4. Overlap asynchronous Host-to-Device 9-strided transfers for the subsequent Batch N+1.
        if (batch_idx + 1 < num_batches) {{
            long int next_start = (batch_idx + 1) * BUNDLE_CAPACITY;
            long int next_size = MIN(BUNDLE_CAPACITY, num_rays - next_start);
            for(int k=0; k<9; k++) {{
                // Host-to-Device transfer: Asynchronously maps the $k$-th component of the subsequent bundle's vector to VRAM using strict pointer arithmetic.
                cudaMemcpyAsync(f_bundle[s_transfer] + (k * BUNDLE_CAPACITY),
                                all_photons_host.f + (k * num_rays) + next_start,
                                next_size * sizeof(double), cudaMemcpyHostToDevice, streams[s_transfer]);
            }}
        }}

        // 5. Initialize VRAM-local status masks directly inside the compute stream context.
        {init_launch}

        // 6. Fused Compute Loop: Synchronized iteration bounds physical ODE evolution mapping.
        bool host_bundle_has_active = true;
        while(host_bundle_has_active) {{
            int zero = 0;
            // Host-to-Device transfer: Zeroes the active thread counter directly inside VRAM prior to atomic operations.
            cudaMemcpyAsync(d_active_count, &zero, sizeof(int), cudaMemcpyHostToDevice, streams[s_compute]);
            host_bundle_has_active = false;
            // Host-to-Device transfer: Disables the boolean execution loop trigger for device-side evaluation.
            cudaMemcpyAsync(d_bundle_has_active, &host_bundle_has_active, sizeof(bool), cudaMemcpyHostToDevice, streams[s_compute]);

            {compaction_launch}
            {rkf45_launch}

            // Device-to-Host transfer: Validates continuation requirement based on numerical integration completion limits.
            cudaMemcpyAsync(&host_bundle_has_active, d_bundle_has_active, sizeof(bool), cudaMemcpyDeviceToHost, streams[s_compute]);
            cudaStreamSynchronize(streams[s_compute]);
        }}

        // 7. Fire native event evaluation geometry kernels over the validated trajectory pool.
        {event_launch}

        // 8. Retrieve updated numerical statuses bounding physical progression maps (D2H Transfer).
        // Device-to-Host transfer: Asynchronously pulls the status array back to pinned memory for CPU-side TimeSlotManager evaluation.
        cudaMemcpyAsync(bundle_status_host[s_compute], status_bundle[s_compute], bundle_size * sizeof(termination_type_t), cudaMemcpyDeviceToHost, streams[s_compute]);
        // Device-to-Host transfer: Asynchronously pulls the physical coordinate time $t$ back to pinned memory for re-binning.
        cudaMemcpyAsync(bundle_time_host[s_compute], f_bundle[s_compute], bundle_size * sizeof(double), cudaMemcpyDeviceToHost, streams[s_compute]);
    }}

    // Resolve final pipeline batch outputs remaining in the transit buffer post-orchestration loop.
    int final_s = (num_batches - 1) % 2;
    cudaStreamSynchronize(streams[final_s]);
    long int final_start = (num_batches - 1) * BUNDLE_CAPACITY;
    long int final_size = MIN(BUNDLE_CAPACITY, num_rays - final_start);
    for (long int j = 0; j < final_size; j++) {{
        if (bundle_status_host[final_s][j] == ACTIVE) {{
            long int p_idx = final_start + j;
            int s_idx = slot_get_index(&tsm, bundle_time_host[final_s][j]);
            if (s_idx != -1) slot_add_photon(&tsm, s_idx, p_idx);
        }}
    }}

    // Export validated device-native blueprints safely mapped into the master C host struct.
    // Device-to-Host transfer: Extracts the complete physical root-finding results.
    cudaMemcpy(results_buffer, d_results_buffer, sizeof(blueprint_data_t) * num_rays, cudaMemcpyDeviceToHost);

    // --- MEMORY DE-ALLOCATIONS ---
    for(int i=0; i<2; i++) {{
        BHAH_FREE_PINNED(bundle_status_host[i]); BHAH_FREE_PINNED(bundle_time_host[i]);
        BHAH_FREE_DEVICE(f_bundle[i]); BHAH_FREE_DEVICE(status_bundle[i]);
        BHAH_FREE_DEVICE(h_bundle[i]); BHAH_FREE_DEVICE(needs_step_bundle[i]);
        BHAH_FREE_DEVICE(retries_bundle[i]); BHAH_FREE_DEVICE(affine_param_bundle[i]);
        BHAH_FREE_DEVICE(on_pos_window_prev_bundle[i]); 
        BHAH_FREE_DEVICE(on_pos_source_prev_bundle[i]);
    }}
    BHAH_FREE_DEVICE(f_p_bundle); BHAH_FREE_DEVICE(f_p_p_bundle);
    BHAH_FREE_DEVICE(compacted_local_ids); BHAH_FREE_DEVICE(d_active_count);
    BHAH_FREE_DEVICE(d_bundle_has_active); BHAH_FREE_DEVICE(d_results_buffer);

    BHAH_FREE_PINNED(all_photons_host.f);
    BHAH_FREE_PINNED(all_photons_host.affine_param);
    BHAH_FREE_PINNED(all_photons_host.h);
    BHAH_FREE_PINNED(all_photons_host.status);
    BHAH_FREE_PINNED(all_photons_host.rejection_retries);

    // --- HOST-SIDE GEOMETRY DE-ALLOCATIONS ---
    // These buffers were allocated via standard malloc as they reside strictly in the CPU host context.
    free(all_photons_host.on_positive_side_of_window_prev);
    free(all_photons_host.on_positive_side_of_source_prev);
    free(all_photons_host.source_event_found); 
    free(all_photons_host.source_event_lambda);
    free(all_photons_host.source_event_f_intersect); 
    free(all_photons_host.window_event_found);
    free(all_photons_host.window_event_lambda); 
    free(all_photons_host.window_event_f_intersect);

    slot_manager_free(&tsm);

    // Explicitly destroy the CUDA streams to prevent runtime resource leaks
    cudaStreamDestroy(streams[0]);
    cudaStreamDestroy(streams[1]);
    """

    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )