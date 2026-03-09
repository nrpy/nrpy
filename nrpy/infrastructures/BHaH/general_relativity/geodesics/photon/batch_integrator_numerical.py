"""
Orchestration module for the Project Singularity-Axiom numerical integration pipeline.

This module generates the high-level C orchestrator responsible for managing the
life cycle of photon trajectories in numerical spacetimes. It strictly implements a 
Split-Pipeline architecture, decoupling the Runge-Kutta-Fehlberg 4(5) integration 
kernels. Intermediate tensors and state vectors are persisted in Global VRAM via 
flattened Structure of Arrays (SoA) scratchpad bundles. This design avoids register 
spilling and cache thrashing on the hardware by trading maximum VRAM bandwidth for 
register stability.
Author: Dalton J. Moone.
"""

import nrpy.c_function as cfc
import nrpy.params as par

def batch_integrator_numerical(spacetime_name: str) -> None:
    """
    Generate the Native CUDA orchestrator for the batched numerical integration pipeline.

    :param spacetime_name: The identifier for the spacetime metric (e.g., 'KerrSchild').
    :raises ValueError: If an unsupported spacetime metric identifier is provided.
    """
    
    # Register core physics and numerical simulation parameters to the global struct.
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

    prefunc = " "

    includes = [
        "BHaH_defines.h", 
        "BHaH_global_device_defines.h",
        "BHaH_function_prototypes.h", 
        "cuda_runtime.h", 
        "cuda_intrinsics.h"
    ]

    desc = r"""@brief Orchestrates the batched Split-Pipeline relativistic ray tracing loop.

    This function serves as the central Host-bound CPU loop for evaluating photon
    geodesics. It utilizes a TimeSlotManager to bin active rays by their physical 
    coordinate time $t$. The Split-Pipeline architecture maps mathematical tensors 
    like $g_{\mu\nu}$ and $\Gamma^\alpha_{\beta\gamma}$ to VRAM scratchpads to
    respect the 255-register hardware limit per thread on sm_86 architecture.

    @param commondata Struct containing global spacetime and numerical tolerances.
    @param num_rays The total number of photon trajectories to simulate.
    @param results_buffer Device array storing the final physical intersections."""

    cfunc_type = "void"

    name = "batch_integrator_numerical"

    params = "const commondata_struct *restrict commondata, long int num_rays, blueprint_data_t *restrict results_buffer"

    include_CodeParameters_h = True

    body = fr"""
    // --- 1. HOST & DEVICE ALLOCATION ---
    
    // The master host-side Structure of Arrays (SoA) tracking all photons $f^\mu$.
    PhotonStateSoA all_photons_host;

    // Host-to-Device transfer allocation: Pinned memory utilized to maximize PCIe DMA throughput for the state vector $f^\mu$.
    BHAH_MALLOC_PINNED(all_photons_host.f, sizeof(double) * 9 * num_rays);
    // Host-to-Device transfer allocation: Pinned memory utilized to maximize PCIe DMA throughput for the first derivative.
    BHAH_MALLOC_PINNED(all_photons_host.f_p, sizeof(double) * 9 * num_rays);
    // Host-to-Device transfer allocation: Pinned memory utilized to maximize PCIe DMA throughput for the second derivative.
    BHAH_MALLOC_PINNED(all_photons_host.f_p_p, sizeof(double) * 9 * num_rays);
    // Host-to-Device transfer allocation: Pinned memory tracking the physical affine parameter $\lambda$.
    BHAH_MALLOC_PINNED(all_photons_host.affine_param, sizeof(double) * num_rays);
    // Host-to-Device transfer allocation: Pinned memory tracking individual integration step sizes $h$.
    BHAH_MALLOC_PINNED(all_photons_host.h, sizeof(double) * num_rays);
    // Host-to-Device transfer allocation: Pinned memory tracking the current trajectory termination status.
    BHAH_MALLOC_PINNED(all_photons_host.status, sizeof(termination_type_t) * num_rays);
    // Host-to-Device transfer allocation: Pinned memory tracking the number of step-size rejections.
    BHAH_MALLOC_PINNED(all_photons_host.rejection_retries, sizeof(int) * num_rays);
    // Host-to-Device transfer allocation: Pinned memory tracking if the photon was previously on the positive side of the window.
    BHAH_MALLOC_PINNED(all_photons_host.on_positive_side_of_window_prev, sizeof(bool) * num_rays);
    // Host-to-Device transfer allocation: Pinned memory tracking if the photon was previously on the positive side of the source.
    BHAH_MALLOC_PINNED(all_photons_host.on_positive_side_of_source_prev, sizeof(bool) * num_rays);
    // Host-to-Device transfer allocation: Pinned memory tracking the history step $\lambda_{{n-1}}$.
    BHAH_MALLOC_PINNED(all_photons_host.affine_param_p, sizeof(double) * num_rays);
    // Host-to-Device transfer allocation: Pinned memory tracking the history step $\lambda_{{n-2}}$.
    BHAH_MALLOC_PINNED(all_photons_host.affine_param_p_p, sizeof(double) * num_rays);
    // Host-to-Device transfer allocation: Pinned memory locking the observer window intersection.
    BHAH_MALLOC_PINNED(all_photons_host.window_event_found, sizeof(bool) * num_rays);
    // Host-to-Device transfer allocation: Pinned memory locking the source emission plane intersection.
    BHAH_MALLOC_PINNED(all_photons_host.source_event_found, sizeof(bool) * num_rays);

    // Arrays of CUDA streams for asynchronous hardware orchestration.
    cudaStream_t streams[2];
    cudaStreamCreate(&streams[0]);
    cudaStreamCreate(&streams[1]);

    // Host-to-Device transfer: Maps the global spacetime constants to the device cache to ensure zero-latency read access.
    cudaMemcpyToSymbol(d_commondata, commondata, sizeof(commondata_struct));

    // --- DOUBLE-BUFFERED BRIDGE ARRAYS ---
    // Extraction buffer used by the TimeSlotManager to map sparse indices to contiguous execution blocks.
    long int *chunk_buffer[2];
    // Bridge array for the state vector $f^\mu$.
    double *f_bridge[2];
    // Bridge array for the first derivative of the state vector.
    double *f_p_bridge[2];
    // Bridge array for the second derivative of the state vector.
    double *f_p_p_bridge[2];
    // Bridge array storing the affine parameter $\lambda$.
    double *affine_bridge[2];
    // Bridge array storing the current integration step size $h$.
    double *h_bridge[2];
    // Bridge array storing the current trajectory termination status.
    termination_type_t *status_bridge[2];
    // Bridge array counting the number of step-size rejections.
    int *retries_bridge[2];
    // Bridge array tracking if the photon was previously on the positive side of the window.
    bool *on_pos_window_prev_bridge[2];
    // Bridge array tracking if the photon was previously on the positive side of the source plane.
    bool *on_pos_source_prev_bridge[2];
    // Bridge array storing the historical affine parameter $\lambda_{{n-1}}$ for chunked transfers.
    double *affine_p_bridge[2];
    // Bridge array storing the historical affine parameter $\lambda_{{n-2}}$ for chunked transfers.
    double *affine_p_p_bridge[2];
    // Bridge array tracking if the window event was previously found.
    bool *window_event_found_bridge[2];
    // Bridge array tracking if the source event was previously found.
    bool *source_event_found_bridge[2];

    // --- DOUBLE-BUFFERED VRAM SCRATCHPADS ---
    // VRAM scratchpad tracking the current state vector $f^\mu$ bounding the RKF45 step.
    double *d_f_bundle[2];
    // VRAM scratchpad locking the anchor state vector $f_{{start}}$ to calculate the final stage update.
    double *d_f_start_bundle[2];
    // VRAM scratchpad tracking the intermediate cumulative RKF45 stage updates.
    double *d_f_temp_bundle[2];
    // VRAM scratchpad tracking the history state $f^\mu_{{n-1}}$ for geometric intersection detection.
    double *d_f_prev_bundle[2];
    // VRAM scratchpad tracking the history state $f^\mu_{{n-2}}$ for geometric intersection detection.
    double *d_f_pre_prev_bundle[2];
    // VRAM scratchpad persisting the symmetric metric tensor $g_{{\mu\nu}}$.
    double *d_metric_bundle[2];
    // VRAM scratchpad persisting the Christoffel symbols $\Gamma^\alpha_{{\beta\gamma}}$.
    double *d_connection_bundle[2];
    // Derivative tensor storing $\dot{{f}}$ across all 6 intermediate RKF45 stages.
    double *d_k_bundle[2];
    // VRAM array regulating active integration step sizing $h$.
    double *d_h[2];
    // VRAM array regulating total affine parameter progress $\lambda$.
    double *d_affine[2];
    // VRAM array holding the current trajectory status limits.
    termination_type_t *d_status[2];
    // VRAM array tracking sequential error rejections per photon.
    int *d_retries[2];
    // VRAM array flagging the previous observer window boundary side.
    bool *d_on_pos_window_prev[2];
    // VRAM array flagging the previous source emission boundary side.
    bool *d_on_pos_source_prev[2];
    // VRAM array tracking historical affine parameter $\lambda_{{n-1}}$.
    double *d_affine_prev[2];
    // VRAM array tracking historical affine parameter $\lambda_{{n-2}}$.
    double *d_affine_pre_prev[2];
    // VRAM array guarding the window intersection coordinates from multi-trigger overwrites.
    bool *d_window_event_found[2];
    // VRAM array guarding the source intersection coordinates from multi-trigger overwrites.
    bool *d_source_event_found[2];
    // VRAM array carrying the absolute master indices $m_{{idx}}$ mapping the execution chunk.
    long int *d_chunk_buffer[2];

    // Loop iterator for instantiating the double-buffered operational arrays.
    for (int s = 0; s < 2; ++s) {{
        // Allocate Bridge arrays in Host Pinned Memory to allow fast, overlapped chunked transfers to the GPU.
        BHAH_MALLOC_PINNED(chunk_buffer[s], sizeof(long int) * BUNDLE_CAPACITY);
        BHAH_MALLOC_PINNED(f_bridge[s], sizeof(double) * 9 * BUNDLE_CAPACITY);
        BHAH_MALLOC_PINNED(f_p_bridge[s], sizeof(double) * 9 * BUNDLE_CAPACITY);
        BHAH_MALLOC_PINNED(f_p_p_bridge[s], sizeof(double) * 9 * BUNDLE_CAPACITY);
        BHAH_MALLOC_PINNED(affine_bridge[s], sizeof(double) * BUNDLE_CAPACITY);
        BHAH_MALLOC_PINNED(h_bridge[s], sizeof(double) * BUNDLE_CAPACITY);
        BHAH_MALLOC_PINNED(status_bridge[s], sizeof(termination_type_t) * BUNDLE_CAPACITY);
        BHAH_MALLOC_PINNED(retries_bridge[s], sizeof(int) * BUNDLE_CAPACITY);
        BHAH_MALLOC_PINNED(on_pos_window_prev_bridge[s], sizeof(bool) * BUNDLE_CAPACITY);
        BHAH_MALLOC_PINNED(on_pos_source_prev_bridge[s], sizeof(bool) * BUNDLE_CAPACITY);
        BHAH_MALLOC_PINNED(affine_p_bridge[s], sizeof(double) * BUNDLE_CAPACITY);
        BHAH_MALLOC_PINNED(affine_p_p_bridge[s], sizeof(double) * BUNDLE_CAPACITY);
        BHAH_MALLOC_PINNED(window_event_found_bridge[s], sizeof(bool) * BUNDLE_CAPACITY);
        BHAH_MALLOC_PINNED(source_event_found_bridge[s], sizeof(bool) * BUNDLE_CAPACITY);

        // Allocate 1D VRAM Scratchpad arrays ensuring strict adherence to the hardware bounds.
        BHAH_MALLOC_DEVICE(d_f_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_f_start_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_f_temp_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_f_prev_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_f_pre_prev_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_metric_bundle[s], sizeof(double) * 10 * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_connection_bundle[s], sizeof(double) * 40 * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_k_bundle[s], sizeof(double) * 6 * 9 * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_h[s], sizeof(double) * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_affine[s], sizeof(double) * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_status[s], sizeof(termination_type_t) * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_retries[s], sizeof(int) * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_on_pos_window_prev[s], sizeof(bool) * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_on_pos_source_prev[s], sizeof(bool) * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_affine_prev[s], sizeof(double) * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_affine_pre_prev[s], sizeof(double) * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_window_event_found[s], sizeof(bool) * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_source_event_found[s], sizeof(bool) * BUNDLE_CAPACITY);
        BHAH_MALLOC_DEVICE(d_chunk_buffer[s], sizeof(long int) * BUNDLE_CAPACITY);
    }}

    // Final intersection blueprint results buffer mapped directly to VRAM.
    blueprint_data_t *d_results_buffer;
    BHAH_MALLOC_DEVICE(d_results_buffer, sizeof(blueprint_data_t) * num_rays);

    // Lock-free temporal TimeSlotManager confined entirely to the Host CPU context.
    TimeSlotManager tsm;
    slot_manager_init(&tsm, commondata->slot_manager_t_min, commondata->t_start + 1.0, commondata->slot_manager_delta_t, num_rays);

    // --- DIAGNOSTIC MEMORY ALLOCATION ---
    // Host pointer tracking the initial conserved quantities prior to integration.
    conserved_quantities_t *initial_cq_host = NULL;
    // Host pointer tracking the terminal conserved quantities post integration.
    conserved_quantities_t *final_cq_host = NULL;

    if (commondata->perform_conservation_check) {{
        // Host-to-Device transfer allocation: Pinned memory utilized to maximize PCIe DMA throughput for the initial diagnostic data.
        BHAH_MALLOC_PINNED(initial_cq_host, sizeof(conserved_quantities_t) * num_rays);
        // Host-to-Device transfer allocation: Pinned memory utilized to maximize PCIe DMA throughput for the final diagnostic data.
        BHAH_MALLOC_PINNED(final_cq_host, sizeof(conserved_quantities_t) * num_rays);
    }}


    // --- 2. INITIALIZATION PHASE ---
    
    // 3D array storing the spatial Cartesian coordinates of the observer window center.
    double window_center_out[3]; 
    // 3D orthonormal basis vector pointing along the x-axis of the local window geometry.
    double n_x_out[3]; 
    // 3D orthonormal basis vector pointing along the y-axis of the local window geometry.
    double n_y_out[3]; 
    // 3D orthonormal basis vector pointing along the z-axis of the local window geometry.
    double n_z_out[3];

    // Kernel Launch: Evaluate initial conditions on the Host to populate the master $f^\mu$ state vector.
    // Hardware Justification: Operates synchronously as the Host array must be fully populated before VRAM offloading.
    set_initial_conditions_kernel_{spacetime_name}(commondata, num_rays, &all_photons_host, window_center_out, n_x_out, n_y_out, n_z_out, 0);

    // --- DIAGNOSTIC PROBE: DETAILED INITIAL POSITION & PLACEHOLDER ALIGNMENT ---
    // Algorithmic Step: Scans the master Host SoA immediately following the initialization kernel call.
    // Hardware Justification: This architectural step verifies that the Host-side logic has correctly populated the starting coordinates $x^\mu$ and zeroed the $p_t$ and $\lambda$ components.
    long int init_mismatch_count = 0; 
    long int mismatch_t = 0, mismatch_x = 0, mismatch_y = 0, mismatch_z = 0, mismatch_pt = 0, mismatch_lam = 0;

    for (long int p = 0; p < num_rays; p++) {{
        const double t_check = all_photons_host.f[0 * num_rays + p];
        const double x_check = all_photons_host.f[1 * num_rays + p];
        const double y_check = all_photons_host.f[2 * num_rays + p];
        const double z_check = all_photons_host.f[3 * num_rays + p];
        const double pt_check = all_photons_host.f[4 * num_rays + p];
        const double lam_check = all_photons_host.f[8 * num_rays + p];

        bool fail_t = fabs(t_check - commondata->t_start) > 1e-10;
        bool fail_x = fabs(x_check - commondata->camera_pos_x) > 1e-10;
        bool fail_y = fabs(y_check - commondata->camera_pos_y) > 1e-10;
        bool fail_z = fabs(z_check - commondata->camera_pos_z) > 1e-10;
        bool fail_pt = fabs(pt_check) > 1e-15;
        bool fail_lam = fabs(lam_check) > 1e-15;

        if (fail_t) mismatch_t++;
        if (fail_x) mismatch_x++;
        if (fail_y) mismatch_y++;
        if (fail_z) mismatch_z++;
        if (fail_pt) mismatch_pt++;
        if (fail_lam) mismatch_lam++;

        if (fail_t || fail_x || fail_y || fail_z || fail_pt || fail_lam) {{
            init_mismatch_count++;
        }}
    }}

    if (init_mismatch_count > 0) {{
        const double mismatch_percent = ((double)init_mismatch_count / (double)num_rays) * 100.0;
        printf("[DIAGNOSTIC] Initialization Alignment Check: %ld out of %ld rays (%.2f%%) fail coordinate/placeholder validation.\n", init_mismatch_count, num_rays, mismatch_percent);
    }}

    // Total integer calculation defining total iterative blocks required to process all photon indices.
    long int num_batches = (num_rays + BUNDLE_CAPACITY - 1) / BUNDLE_CAPACITY;

    // Loop iterator for evaluating the initialization constraint across sequential blocks.
    for (long int init_batch = 0; init_batch < num_batches; ++init_batch) {{
        long int start_idx = init_batch * BUNDLE_CAPACITY;
        long int chunk_size = MIN((long int)BUNDLE_CAPACITY, num_rays - start_idx);

        // Loop index iterating over the specific initialization batch elements to pack the bridge.
        for (int init_i = 0; init_i < chunk_size; ++init_i) {{
            long int master_idx = start_idx + init_i;
            for (int init_k = 0; init_k < 9; ++init_k) {{
                f_bridge[0][init_k * BUNDLE_CAPACITY + init_i] = all_photons_host.f[init_k * num_rays + master_idx];
            }}
        }}

        // H2D: Asynchronously pushes strictly bounded initial states to VRAM on stream 0.
        for (int c_k = 0; c_k < 9; ++c_k) {{
            cudaMemcpyAsync(d_f_bundle[0] + c_k * BUNDLE_CAPACITY, f_bridge[0] + c_k * BUNDLE_CAPACITY, sizeof(double) * chunk_size, cudaMemcpyHostToDevice, streams[0]);
        }}

        // Kernel Launch: Calculate $g_{{\mu\nu}}$ required for the Hamiltonian constraint strictly on stream 0.
        interpolation_kernel_{spacetime_name}(d_f_bundle[0], d_metric_bundle[0], NULL, chunk_size, 0);

        // --- DIAGNOSTIC PROBE: VRAM METRIC INTEGRITY CHECK ---
        double *metric_diag_bridge; 
        BHAH_MALLOC_PINNED(metric_diag_bridge, sizeof(double) * 10 * BUNDLE_CAPACITY); 
        
        // Device-to-Host transfer: Retrieves the symmetric metric tensor $g_{{\mu\nu}}$ for inspection. Only retrieves metric components for the active chunk.
        for (int m_k = 0; m_k < 10; ++m_k) {{
            cudaMemcpyAsync(metric_diag_bridge + m_k * BUNDLE_CAPACITY, d_metric_bundle[0] + m_k * BUNDLE_CAPACITY, sizeof(double) * chunk_size, cudaMemcpyDeviceToHost, streams[0]);
        }}
        // Hard synchronization barrier to inspect the payload prior to constraint solving.
        cudaStreamSynchronize(streams[0]);

        long int metric_nan_count = 0;
        for (int m_diag_i = 0; m_diag_i < chunk_size; ++m_diag_i) {{
            bool m_has_nan = false;
            for (int m_diag_k = 0; m_diag_k < 10; ++m_diag_k) {{
                if (isnan(metric_diag_bridge[m_diag_k * BUNDLE_CAPACITY + m_diag_i]) || 
                    isinf(metric_diag_bridge[m_diag_k * BUNDLE_CAPACITY + m_diag_i])) {{
                    m_has_nan = true;
                    break;
                }}
            }}
            if (m_has_nan) metric_nan_count++;
        }}

        if (metric_nan_count > 0) {{
            printf("[DIAGNOSTIC] Init Batch %ld: %ld rays have invalid Metric G_mu_nu before p_t solve.\n", init_batch, metric_nan_count);
        }}
        BHAH_FREE_PINNED(metric_diag_bridge); 
        
        // Kernel Launch: Solves the constraint $p_\mu p^\mu = 0$ to find the temporal momentum $p_t$ on stream 0.
        p0_reverse_kernel(d_f_bundle[0], d_metric_bundle[0], chunk_size, 0);

        // Device-to-Host transfer: Retrieves the mathematically constrained state vector back to CPU RAM.
        // NEW D2H: Retrieves strictly bounded constrained state vectors back to CPU RAM.
        for (int c_k = 0; c_k < 9; ++c_k) {{
            cudaMemcpyAsync(f_bridge[0] + c_k * BUNDLE_CAPACITY, d_f_bundle[0] + c_k * BUNDLE_CAPACITY, sizeof(double) * chunk_size, cudaMemcpyDeviceToHost, streams[0]);
        }}
        cudaStreamSynchronize(streams[0]);

        // Unpack the validated constraint states back to the master SoA array.
        long int nan_count = 0; 
        for (int gather_i = 0; gather_i < chunk_size; ++gather_i) {{
            long int master_idx = start_idx + gather_i;
            bool has_nan = false;
            for (int gather_k = 0; gather_k < 9; ++gather_k) {{
                double val = f_bridge[0][gather_k * BUNDLE_CAPACITY + gather_i];
                all_photons_host.f[gather_k * num_rays + master_idx] = val;
                if (isnan(val)) has_nan = true;
            }}
            if (has_nan) nan_count++;
        }}

        if (nan_count > 0) {{
            printf("[DIAGNOSTIC] Init Batch %ld: %ld rays contain NaN in state f^mu after p_t solve.\n", init_batch, nan_count);
        }}
    }}

    // --- BASELINE CONSERVED QUANTITIES ---
    // Algorithmic Step: Evaluate initial conserved quantities immediately after generating valid physical null states.
    // Hardware Justification: This baselines the data via chunked VRAM kernels before the Split-Pipeline begins mutating the state vectors.
    if (commondata->perform_conservation_check) {{
        calculate_conserved_quantities_universal_{spacetime_name}_photon(&all_photons_host, num_rays, initial_cq_host);
    }}


    // Loop iterator traversing the entire global ray count to synchronize starting properties.
    long int sync_i;
    for(sync_i = 0; sync_i < num_rays; ++sync_i) {{
        // Loop iterator updating tensor states.
        int sync_k;
        for (sync_k = 0; sync_k < 9; ++sync_k) {{
            all_photons_host.f_p[sync_k * num_rays + sync_i] = all_photons_host.f[sync_k * num_rays + sync_i];
            all_photons_host.f_p_p[sync_k * num_rays + sync_i] = all_photons_host.f[sync_k * num_rays + sync_i];
        }}
        all_photons_host.status[sync_i] = ACTIVE;
        all_photons_host.affine_param[sync_i] = 0.0;
        all_photons_host.rejection_retries[sync_i] = 0;

        // Initializes the affine parameter histories and sets the intersection locks to false.
        all_photons_host.affine_param_p[sync_i] = 0.0;
        all_photons_host.affine_param_p_p[sync_i] = 0.0;
        all_photons_host.window_event_found[sync_i] = false;
        all_photons_host.source_event_found[sync_i] = false;

        // Integer representing the assigned temporal slot for the current photon.
        int s_idx = slot_get_index(&tsm, all_photons_host.f[sync_i]);
        if (s_idx != -1) {{
            slot_add_photon(&tsm, s_idx, sync_i);
        }}
    }}


    // --- 3. TEMPORAL LOOP (The Engine) ---

    // Integer tracking the global number of active photon trajectories to allow early loop termination.
    long int total_active_photons = num_rays; 

    // Outer loop iterator for the physical time bins.
    for (int slot_idx = tsm.num_slots - 1; slot_idx >= 0; --slot_idx) {{ 
        // Evaluates the early exit condition to terminate the temporal engine if all geometric trajectories have concluded.
        if (total_active_photons <= 0) {{
            break;
        }}

        // Variables tracking the active stream state machine to overlap Host-to-Device latency.
        int current = 0;
        int next = 1;
        long int active_chunks[2] = {{0, 0}};

        // --- PHASE A: PRIME THE PUMP (Stream 0) ---
        // Algorithmic Step: Populate the first bridge and launch asynchronous integration on the primary stream.
        // Hardware Justification: Initializes the pipeline to allow subsequent Host-side packing to overlap with GPU compute.
        active_chunks[current] = MIN((long int)BUNDLE_CAPACITY, tsm.slot_counts[slot_idx]); 
        
        if (active_chunks[current] > 0) {{
            slot_remove_chunk(&tsm, slot_idx, chunk_buffer[current], active_chunks[current]); 

            // 1. Pack Host Data into Bridge Arrays [current]
            for (int bridge_i = 0; bridge_i < active_chunks[current]; ++bridge_i) {{
                long int m_idx = chunk_buffer[current][bridge_i];
                for (int c_k = 0; c_k < 9; ++c_k) {{
                    f_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f[c_k * num_rays + m_idx]; 
                    f_p_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p[c_k * num_rays + m_idx];
                    f_p_p_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p_p[c_k * num_rays + m_idx];
                }}
                h_bridge[current][bridge_i] = all_photons_host.h[m_idx]; 
                status_bridge[current][bridge_i] = all_photons_host.status[m_idx];
                retries_bridge[current][bridge_i] = all_photons_host.rejection_retries[m_idx]; 
                affine_bridge[current][bridge_i] = all_photons_host.affine_param[m_idx]; 
                on_pos_window_prev_bridge[current][bridge_i] = all_photons_host.on_positive_side_of_window_prev[m_idx]; 
                on_pos_source_prev_bridge[current][bridge_i] = all_photons_host.on_positive_side_of_source_prev[m_idx]; 
                affine_p_bridge[current][bridge_i] = all_photons_host.affine_param_p[m_idx]; 
                affine_p_p_bridge[current][bridge_i] = all_photons_host.affine_param_p_p[m_idx]; 
                window_event_found_bridge[current][bridge_i] = all_photons_host.window_event_found[m_idx]; 
                source_event_found_bridge[current][bridge_i] = all_photons_host.source_event_found[m_idx]; 
            }}

            // 2. Launch ASYNC Host-to-Device Transfers on streams[current]
            for (int c_k = 0; c_k < 9; ++c_k) {{
                cudaMemcpyAsync(d_f_bundle[current] + c_k * BUNDLE_CAPACITY, f_bridge[current] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[current], cudaMemcpyHostToDevice, streams[current]);
                cudaMemcpyAsync(d_f_prev_bundle[current] + c_k * BUNDLE_CAPACITY, f_p_bridge[current] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[current], cudaMemcpyHostToDevice, streams[current]); 
                cudaMemcpyAsync(d_f_pre_prev_bundle[current] + c_k * BUNDLE_CAPACITY, f_p_p_bridge[current] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[current], cudaMemcpyHostToDevice, streams[current]);
            }}
            cudaMemcpyAsync(d_h[current], h_bridge[current], sizeof(double) * active_chunks[current], cudaMemcpyHostToDevice, streams[current]);
            cudaMemcpyAsync(d_status[current], status_bridge[current], sizeof(termination_type_t) * active_chunks[current], cudaMemcpyHostToDevice, streams[current]); 
            cudaMemcpyAsync(d_retries[current], retries_bridge[current], sizeof(int) * active_chunks[current], cudaMemcpyHostToDevice, streams[current]); 
            cudaMemcpyAsync(d_affine[current], affine_bridge[current], sizeof(double) * active_chunks[current], cudaMemcpyHostToDevice, streams[current]);
            cudaMemcpyAsync(d_on_pos_window_prev[current], on_pos_window_prev_bridge[current], sizeof(bool) * active_chunks[current], cudaMemcpyHostToDevice, streams[current]); 
            cudaMemcpyAsync(d_on_pos_source_prev[current], on_pos_source_prev_bridge[current], sizeof(bool) * active_chunks[current], cudaMemcpyHostToDevice, streams[current]); 
            cudaMemcpyAsync(d_affine_prev[current], affine_p_bridge[current], sizeof(double) * active_chunks[current], cudaMemcpyHostToDevice, streams[current]);
            cudaMemcpyAsync(d_affine_pre_prev[current], affine_p_p_bridge[current], sizeof(double) * active_chunks[current], cudaMemcpyHostToDevice, streams[current]);
            cudaMemcpyAsync(d_window_event_found[current], window_event_found_bridge[current], sizeof(bool) * active_chunks[current], cudaMemcpyHostToDevice, streams[current]);
            cudaMemcpyAsync(d_source_event_found[current], source_event_found_bridge[current], sizeof(bool) * active_chunks[current], cudaMemcpyHostToDevice, streams[current]);
            cudaMemcpyAsync(d_chunk_buffer[current], chunk_buffer[current], sizeof(long int) * active_chunks[current], cudaMemcpyHostToDevice, streams[current]);

            // 3. Device-to-Device Baseline Synchronization 
            for (int c_k = 0; c_k < 9; ++c_k) {{
                cudaMemcpyAsync(d_f_start_bundle[current] + c_k * BUNDLE_CAPACITY, d_f_bundle[current] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[current], cudaMemcpyDeviceToDevice, streams[current]);
                cudaMemcpyAsync(d_f_temp_bundle[current] + c_k * BUNDLE_CAPACITY, d_f_bundle[current] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[current], cudaMemcpyDeviceToDevice, streams[current]);
            }}

            // 4. Launch the Compute Pipeline on streams[current]
            for (int stage = 1; stage <= 6; ++stage) {{
                interpolation_kernel_{spacetime_name}(d_f_temp_bundle[current], d_metric_bundle[current], d_connection_bundle[current], active_chunks[current], current); 
                calculate_ode_rhs_kernel(d_f_temp_bundle[current], d_metric_bundle[current], d_connection_bundle[current], d_k_bundle[current], stage, active_chunks[current], current); 
                rkf45_stage_update(d_f_start_bundle[current], d_k_bundle[current], d_h[current], stage, active_chunks[current], d_f_temp_bundle[current], current); 
            }}
            
            // Finalize step-size $h$ and detect geometric events
            rkf45_finalize_and_control(d_f_bundle[current], d_f_start_bundle[current], d_k_bundle[current], d_h[current], d_status[current], d_affine[current], d_retries[current], active_chunks[current], current); 
            event_detection_manager_kernel(d_f_bundle[current], d_f_prev_bundle[current], d_f_pre_prev_bundle[current], d_affine[current], d_affine_prev[current], d_affine_pre_prev[current], d_results_buffer, d_status[current],
                d_on_pos_window_prev[current], d_on_pos_source_prev[current], d_window_event_found[current], d_source_event_found[current], d_chunk_buffer[current], active_chunks[current], current);

            // 5. Launch ASYNC Device-to-Host Transfers on streams[current]
            for (int c_k = 0; c_k < 9; ++c_k) {{
                cudaMemcpyAsync(f_bridge[current] + c_k * BUNDLE_CAPACITY, d_f_bundle[current] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[current], cudaMemcpyDeviceToHost, streams[current]);
                cudaMemcpyAsync(f_p_bridge[current] + c_k * BUNDLE_CAPACITY, d_f_prev_bundle[current] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[current], cudaMemcpyDeviceToHost, streams[current]);
                cudaMemcpyAsync(f_p_p_bridge[current] + c_k * BUNDLE_CAPACITY, d_f_pre_prev_bundle[current] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[current], cudaMemcpyDeviceToHost, streams[current]);
            }}
            cudaMemcpyAsync(h_bridge[current], d_h[current], sizeof(double) * active_chunks[current], cudaMemcpyDeviceToHost, streams[current]);
            cudaMemcpyAsync(status_bridge[current], d_status[current], sizeof(termination_type_t) * active_chunks[current], cudaMemcpyDeviceToHost, streams[current]);
            cudaMemcpyAsync(retries_bridge[current], d_retries[current], sizeof(int) * active_chunks[current], cudaMemcpyDeviceToHost, streams[current]);
            cudaMemcpyAsync(affine_bridge[current], d_affine[current], sizeof(double) * active_chunks[current], cudaMemcpyDeviceToHost, streams[current]);
            cudaMemcpyAsync(on_pos_window_prev_bridge[current], d_on_pos_window_prev[current], sizeof(bool) * active_chunks[current], cudaMemcpyDeviceToHost, streams[current]);
            cudaMemcpyAsync(on_pos_source_prev_bridge[current], d_on_pos_source_prev[current], sizeof(bool) * active_chunks[current], cudaMemcpyDeviceToHost, streams[current]);
            cudaMemcpyAsync(affine_p_bridge[current], d_affine_prev[current], sizeof(double) * active_chunks[current], cudaMemcpyDeviceToHost, streams[current]);
            cudaMemcpyAsync(affine_p_p_bridge[current], d_affine_pre_prev[current], sizeof(double) * active_chunks[current], cudaMemcpyDeviceToHost, streams[current]);
            cudaMemcpyAsync(window_event_found_bridge[current], d_window_event_found[current], sizeof(bool) * active_chunks[current], cudaMemcpyDeviceToHost, streams[current]);
            cudaMemcpyAsync(source_event_found_bridge[current], d_source_event_found[current], sizeof(bool) * active_chunks[current], cudaMemcpyDeviceToHost, streams[current]);
        }}

        // --- PHASE B: THE OVERLAP LOOP ---
        // Algorithmic Step: Continuously alternate between streams, packing the next payload while syncing the current.
        // Hardware Justification: Completely hides PCIe DMA transfer latency behind the active RKF45 integration compute time.
        while (active_chunks[current] > 0 || tsm.slot_counts[slot_idx] > 0) {{ 
            
            // QUEUE NEXT: Prepare the upcoming payload on the CPU.
            active_chunks[next] = MIN((long int)BUNDLE_CAPACITY, tsm.slot_counts[slot_idx]); 
            if (active_chunks[next] > 0) {{
                slot_remove_chunk(&tsm, slot_idx, chunk_buffer[next], active_chunks[next]); 
                
                // Pack Host Data into Bridge Arrays using 'next' index
                for (int bridge_i = 0; bridge_i < active_chunks[next]; ++bridge_i) {{
                    long int m_idx = chunk_buffer[next][bridge_i];
                    for (int c_k = 0; c_k < 9; ++c_k) {{
                        f_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f[c_k * num_rays + m_idx];
                        f_p_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p[c_k * num_rays + m_idx];
                        f_p_p_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p_p[c_k * num_rays + m_idx];
                    }}
                    h_bridge[next][bridge_i] = all_photons_host.h[m_idx];
                    status_bridge[next][bridge_i] = all_photons_host.status[m_idx];
                    retries_bridge[next][bridge_i] = all_photons_host.rejection_retries[m_idx];
                    affine_bridge[next][bridge_i] = all_photons_host.affine_param[m_idx];
                    on_pos_window_prev_bridge[next][bridge_i] = all_photons_host.on_positive_side_of_window_prev[m_idx];
                    on_pos_source_prev_bridge[next][bridge_i] = all_photons_host.on_positive_side_of_source_prev[m_idx];
                    affine_p_bridge[next][bridge_i] = all_photons_host.affine_param_p[m_idx];
                    affine_p_p_bridge[next][bridge_i] = all_photons_host.affine_param_p_p[m_idx];
                    window_event_found_bridge[next][bridge_i] = all_photons_host.window_event_found[m_idx];
                    source_event_found_bridge[next][bridge_i] = all_photons_host.source_event_found[m_idx];
                }}

                // Launch ASYNC H2D -> Kernels -> D2H exactly as Phase A, strictly using streams[next]
                for (int c_k = 0; c_k < 9; ++c_k) {{
                    cudaMemcpyAsync(d_f_bundle[next] + c_k * BUNDLE_CAPACITY, f_bridge[next] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[next], cudaMemcpyHostToDevice, streams[next]);
                    cudaMemcpyAsync(d_f_prev_bundle[next] + c_k * BUNDLE_CAPACITY, f_p_bridge[next] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[next], cudaMemcpyHostToDevice, streams[next]);
                    cudaMemcpyAsync(d_f_pre_prev_bundle[next] + c_k * BUNDLE_CAPACITY, f_p_p_bridge[next] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[next], cudaMemcpyHostToDevice, streams[next]);
                }}
                cudaMemcpyAsync(d_h[next], h_bridge[next], sizeof(double) * active_chunks[next], cudaMemcpyHostToDevice, streams[next]);
                cudaMemcpyAsync(d_status[next], status_bridge[next], sizeof(termination_type_t) * active_chunks[next], cudaMemcpyHostToDevice, streams[next]);
                cudaMemcpyAsync(d_retries[next], retries_bridge[next], sizeof(int) * active_chunks[next], cudaMemcpyHostToDevice, streams[next]);
                cudaMemcpyAsync(d_affine[next], affine_bridge[next], sizeof(double) * active_chunks[next], cudaMemcpyHostToDevice, streams[next]);
                cudaMemcpyAsync(d_on_pos_window_prev[next], on_pos_window_prev_bridge[next], sizeof(bool) * active_chunks[next], cudaMemcpyHostToDevice, streams[next]);
                cudaMemcpyAsync(d_on_pos_source_prev[next], on_pos_source_prev_bridge[next], sizeof(bool) * active_chunks[next], cudaMemcpyHostToDevice, streams[next]);
                cudaMemcpyAsync(d_affine_prev[next], affine_p_bridge[next], sizeof(double) * active_chunks[next], cudaMemcpyHostToDevice, streams[next]);
                cudaMemcpyAsync(d_affine_pre_prev[next], affine_p_p_bridge[next], sizeof(double) * active_chunks[next], cudaMemcpyHostToDevice, streams[next]);
                cudaMemcpyAsync(d_window_event_found[next], window_event_found_bridge[next], sizeof(bool) * active_chunks[next], cudaMemcpyHostToDevice, streams[next]);
                cudaMemcpyAsync(d_source_event_found[next], source_event_found_bridge[next], sizeof(bool) * active_chunks[next], cudaMemcpyHostToDevice, streams[next]);
                cudaMemcpyAsync(d_chunk_buffer[next], chunk_buffer[next], sizeof(long int) * active_chunks[next], cudaMemcpyHostToDevice, streams[next]);

                for (int c_k = 0; c_k < 9; ++c_k) {{
                    cudaMemcpyAsync(d_f_start_bundle[next] + c_k * BUNDLE_CAPACITY, d_f_bundle[next] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[next], cudaMemcpyDeviceToDevice, streams[next]);
                    cudaMemcpyAsync(d_f_temp_bundle[next] + c_k * BUNDLE_CAPACITY, d_f_bundle[next] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[next], cudaMemcpyDeviceToDevice, streams[next]);
                }}

                for (int stage = 1; stage <= 6; ++stage) {{ 
                    interpolation_kernel_{spacetime_name}(d_f_temp_bundle[next], d_metric_bundle[next], d_connection_bundle[next], active_chunks[next], next);
                    calculate_ode_rhs_kernel(d_f_temp_bundle[next], d_metric_bundle[next], d_connection_bundle[next], d_k_bundle[next], stage, active_chunks[next], next);
                    rkf45_stage_update(d_f_start_bundle[next], d_k_bundle[next], d_h[next], stage, active_chunks[next], d_f_temp_bundle[next], next);
                }}

                rkf45_finalize_and_control(d_f_bundle[next], d_f_start_bundle[next], d_k_bundle[next], d_h[next], d_status[next], d_affine[next], d_retries[next], active_chunks[next], next);
                event_detection_manager_kernel(d_f_bundle[next], d_f_prev_bundle[next], d_f_pre_prev_bundle[next], d_affine[next], d_affine_prev[next], d_affine_pre_prev[next], d_results_buffer, d_status[next], d_on_pos_window_prev[next], d_on_pos_source_prev[next], d_window_event_found[next], d_source_event_found[next], d_chunk_buffer[next], active_chunks[next], next);

                for (int c_k = 0; c_k < 9; ++c_k) {{
                    cudaMemcpyAsync(f_bridge[next] + c_k * BUNDLE_CAPACITY, d_f_bundle[next] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[next], cudaMemcpyDeviceToHost, streams[next]);
                    cudaMemcpyAsync(f_p_bridge[next] + c_k * BUNDLE_CAPACITY, d_f_prev_bundle[next] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[next], cudaMemcpyDeviceToHost, streams[next]);
                    cudaMemcpyAsync(f_p_p_bridge[next] + c_k * BUNDLE_CAPACITY, d_f_pre_prev_bundle[next] + c_k * BUNDLE_CAPACITY, sizeof(double) * active_chunks[next], cudaMemcpyDeviceToHost, streams[next]);
                }}
                cudaMemcpyAsync(h_bridge[next], d_h[next], sizeof(double) * active_chunks[next], cudaMemcpyDeviceToHost, streams[next]);
                cudaMemcpyAsync(status_bridge[next], d_status[next], sizeof(termination_type_t) * active_chunks[next], cudaMemcpyDeviceToHost, streams[next]);
                cudaMemcpyAsync(retries_bridge[next], d_retries[next], sizeof(int) * active_chunks[next], cudaMemcpyDeviceToHost, streams[next]);
                cudaMemcpyAsync(affine_bridge[next], d_affine[next], sizeof(double) * active_chunks[next], cudaMemcpyDeviceToHost, streams[next]);
                cudaMemcpyAsync(on_pos_window_prev_bridge[next], d_on_pos_window_prev[next], sizeof(bool) * active_chunks[next], cudaMemcpyDeviceToHost, streams[next]);
                cudaMemcpyAsync(on_pos_source_prev_bridge[next], d_on_pos_source_prev[next], sizeof(bool) * active_chunks[next], cudaMemcpyDeviceToHost, streams[next]);
                cudaMemcpyAsync(affine_p_bridge[next], d_affine_prev[next], sizeof(double) * active_chunks[next], cudaMemcpyDeviceToHost, streams[next]);
                cudaMemcpyAsync(affine_p_p_bridge[next], d_affine_pre_prev[next], sizeof(double) * active_chunks[next], cudaMemcpyDeviceToHost, streams[next]);
                cudaMemcpyAsync(window_event_found_bridge[next], d_window_event_found[next], sizeof(bool) * active_chunks[next], cudaMemcpyDeviceToHost, streams[next]);
                cudaMemcpyAsync(source_event_found_bridge[next], d_source_event_found[next], sizeof(bool) * active_chunks[next], cudaMemcpyDeviceToHost, streams[next]);
            }}

            // SYNC CURRENT: CPU waits for the active compute stream to finish, then processes results.
            if (active_chunks[current] > 0) {{
                // Hard synchronization on the specific stream to guarantee memory stability before unpacking.
                cudaStreamSynchronize(streams[current]);

                // Unpack results from Bridge Arrays[current] back into all_photons_host
                for (int fin_i = 0; fin_i < active_chunks[current]; ++fin_i) {{
                    long int m_idx = chunk_buffer[current][fin_i];
                    for (int fin_k = 0; fin_k < 9; ++fin_k) {{
                        all_photons_host.f[fin_k * num_rays + m_idx] = f_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i];
                        all_photons_host.f_p[fin_k * num_rays + m_idx] = f_p_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i];
                        all_photons_host.f_p_p[fin_k * num_rays + m_idx] = f_p_p_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i];
                    }}
                    all_photons_host.h[m_idx] = h_bridge[current][fin_i];
                    all_photons_host.status[m_idx] = status_bridge[current][fin_i];
                    all_photons_host.rejection_retries[m_idx] = retries_bridge[current][fin_i];
                    all_photons_host.affine_param[m_idx] = affine_bridge[current][fin_i];
                    all_photons_host.on_positive_side_of_window_prev[m_idx] = on_pos_window_prev_bridge[current][fin_i];
                    all_photons_host.on_positive_side_of_source_prev[m_idx] = on_pos_source_prev_bridge[current][fin_i];
                    all_photons_host.affine_param_p[m_idx] = affine_p_bridge[current][fin_i];
                    all_photons_host.affine_param_p_p[m_idx] = affine_p_p_bridge[current][fin_i];
                    all_photons_host.window_event_found[m_idx] = window_event_found_bridge[current][fin_i];
                    all_photons_host.source_event_found[m_idx] = source_event_found_bridge[current][fin_i];

                    // Route trajectory status back to the TimeSlotManager.
                    if (status_bridge[current][fin_i] == ACTIVE) {{
                        int next_s_idx = slot_get_index(&tsm, all_photons_host.f[m_idx]); 
                        if (next_s_idx != -1) {{
                            slot_add_photon(&tsm, next_s_idx, m_idx); 
                        }} else {{
                            all_photons_host.status[m_idx] = FAILURE_T_MAX_EXCEEDED; 
                            total_active_photons--;
                        }}
                    }} else if (status_bridge[current][fin_i] == REJECTED) {{ 
                        // Re-add to current bin to attempt integration with an adapted step-size scalar $h$.
                        slot_add_photon(&tsm, slot_idx, m_idx); 
                    }} else {{
                        // Decrements the global counter as the photon has reached a terminal state.
                        total_active_photons--;
                    }}
                }}
                // Mark the current chunk as fully processed.
                active_chunks[current] = 0; 
            }}

            // SWAP POINTERS: Toggle the conveyor belts for the next cycle.
            int temp = current;
            current = next;
            next = temp;
        }}

        // --- PHASE C: THE TIME BARRIER ---
        // Algorithmic Step: Enforce a rigid hardware sync before advancing the physical time clock.
        // Hardware Justification: Prevents race conditions and ensures all $f^\mu$ states strictly adhere to the current temporal bin limits.
        cudaDeviceSynchronize(); 
    }}

    // --- 4. CLEANUP & FINALIZATION ---
    
    // Device-to-Host transfer: Extracts validated device-native blueprints $b_i$ containing geometric plane intersections.
    cudaMemcpy(results_buffer, d_results_buffer, sizeof(blueprint_data_t) * num_rays, cudaMemcpyDeviceToHost); 

    // Kernel Launch: Processes escaped photons intersecting the celestial sphere $r > r_{{escape}}$ 
    // appending 0 as the final argument since we are relying strictly on the primary hardware stream.
    calculate_and_fill_blueprint_data_universal(&all_photons_host, num_rays, results_buffer, 0);

    // Loop iterator purging the double-buffered arrays across both hardware streams.
    for (int s = 0; s < 2; ++s) {{
        // Host Memory Free: Purges bridge components supporting scatter logic mapped to PCIe DMA transfers.
        BHAH_FREE_PINNED(chunk_buffer[s]);
        BHAH_FREE_PINNED(f_bridge[s]);
        BHAH_FREE_PINNED(f_p_bridge[s]);
        BHAH_FREE_PINNED(f_p_p_bridge[s]);
        BHAH_FREE_PINNED(affine_bridge[s]);
        BHAH_FREE_PINNED(h_bridge[s]);
        BHAH_FREE_PINNED(status_bridge[s]);
        BHAH_FREE_PINNED(retries_bridge[s]);
        BHAH_FREE_PINNED(on_pos_window_prev_bridge[s]);
        BHAH_FREE_PINNED(on_pos_source_prev_bridge[s]);
        BHAH_FREE_PINNED(affine_p_bridge[s]);
        BHAH_FREE_PINNED(affine_p_p_bridge[s]);
        BHAH_FREE_PINNED(window_event_found_bridge[s]);
        BHAH_FREE_PINNED(source_event_found_bridge[s]);

        // Device Memory Free: Purges remaining VRAM operational pipeline scratchpads.
        BHAH_FREE_DEVICE(d_f_bundle[s]);
        BHAH_FREE_DEVICE(d_f_start_bundle[s]);
        BHAH_FREE_DEVICE(d_f_temp_bundle[s]);
        BHAH_FREE_DEVICE(d_f_prev_bundle[s]);
        BHAH_FREE_DEVICE(d_f_pre_prev_bundle[s]);
        BHAH_FREE_DEVICE(d_metric_bundle[s]);
        BHAH_FREE_DEVICE(d_connection_bundle[s]);
        BHAH_FREE_DEVICE(d_k_bundle[s]);
        BHAH_FREE_DEVICE(d_h[s]);
        BHAH_FREE_DEVICE(d_affine[s]);
        BHAH_FREE_DEVICE(d_status[s]);
        BHAH_FREE_DEVICE(d_retries[s]);
        BHAH_FREE_DEVICE(d_on_pos_window_prev[s]);
        BHAH_FREE_DEVICE(d_on_pos_source_prev[s]);
        BHAH_FREE_DEVICE(d_affine_prev[s]);
        BHAH_FREE_DEVICE(d_affine_pre_prev[s]);
        BHAH_FREE_DEVICE(d_window_event_found[s]);
        BHAH_FREE_DEVICE(d_source_event_found[s]);
        BHAH_FREE_DEVICE(d_chunk_buffer[s]);

        // Hardware Resource Free: Destroys the asynchronous orchestration streams.
        cudaStreamDestroy(streams[s]);
    }}


    // --- CPU CONSERVATION DRIFT EVALUATION ---
    if (commondata->perform_conservation_check) {{
        // Kernel Launch: Stream the final state vectors to VRAM to calculate terminal conserved quantities.
        calculate_conserved_quantities_universal_{spacetime_name}_photon(&all_photons_host, num_rays, final_cq_host);

        printf("\n=================================================\n");
        printf(" CONSERVED QUANTITIES DIAGNOSTIC REPORT\n");
        printf("=================================================\n");
        
        // Scalar variables tracking the maximum recorded relative drift for each physical quantity.
        double max_err_E = 0.0, max_err_Lz = 0.0, max_err_Q = 0.0;
        // Absolute master indices $m_{{idx}}$ identifying the trajectory responsible for the maximum numerical drift.
        long int worst_ray_E = -1, worst_ray_Lz = -1, worst_ray_Q = -1;

        // Loop iterator spanning the global dataset to calculate relative errors.
        for (long int i = 0; i < num_rays; i++) {{
            // Hardware Justification: Evaluate relative numerical drift natively on the CPU to prevent VRAM bottlenecks and leverage complex print formatting.
            double err_E = fabs((final_cq_host[i].E - initial_cq_host[i].E) / (initial_cq_host[i].E + 1e-15));
            double err_Lz = fabs((final_cq_host[i].Lz - initial_cq_host[i].Lz) / (initial_cq_host[i].Lz + 1e-15));
            double err_Q = fabs((final_cq_host[i].Q - initial_cq_host[i].Q) / (initial_cq_host[i].Q + 1e-15));

            if (err_E > max_err_E) {{ max_err_E = err_E; worst_ray_E = i; }}
            if (err_Lz > max_err_Lz) {{ max_err_Lz = err_Lz; worst_ray_Lz = i; }}
            if (err_Q > max_err_Q) {{ max_err_Q = err_Q; worst_ray_Q = i; }}
        }}

        printf("  Max Relative Error (Energy E): %e (Ray %ld)\n", max_err_E, worst_ray_E);
        printf("  Max Relative Error (Momentum Lz): %e (Ray %ld)\n", max_err_Lz, worst_ray_Lz);
        printf("  Max Relative Error (Carter Q): %e (Ray %ld)\n", max_err_Q, worst_ray_Q);
        printf("=================================================\n\n");

        // Host Memory Free: Purges pinned diagnostic buffers.
        BHAH_FREE_PINNED(initial_cq_host);
        BHAH_FREE_PINNED(final_cq_host);
    }}

    // Host Memory Free: Purges the primary Host array states f^mu and affine parameters \lambda.
    BHAH_FREE_PINNED(all_photons_host.f);
    BHAH_FREE_PINNED(all_photons_host.f_p);
    BHAH_FREE_PINNED(all_photons_host.f_p_p);
    BHAH_FREE_PINNED(all_photons_host.affine_param);
    BHAH_FREE_PINNED(all_photons_host.h);
    BHAH_FREE_PINNED(all_photons_host.status);
    BHAH_FREE_PINNED(all_photons_host.rejection_retries);
    BHAH_FREE_PINNED(all_photons_host.on_positive_side_of_window_prev);
    BHAH_FREE_PINNED(all_photons_host.on_positive_side_of_source_prev);
    BHAH_FREE_PINNED(all_photons_host.affine_param_p);
    BHAH_FREE_PINNED(all_photons_host.affine_param_p_p);
    BHAH_FREE_PINNED(all_photons_host.window_event_found);
    BHAH_FREE_PINNED(all_photons_host.source_event_found);
    
    // (The legacy un-indexed bridge and bundle frees have been removed from here)

    // Device Memory Free: Purges the final single-pointer intersection blueprint buffer.
    BHAH_FREE_DEVICE(d_results_buffer);

    // Memory Free: Purges the temporal sorting struct mapping the Host-side execution grid.
    slot_manager_free(&tsm);
    """


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