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

    prefunc = ""

    includes = [
        "BHaH_defines.h", 
        "BHaH_device_defines.h", 
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

    // Extraction buffer used by the TimeSlotManager to map sparse indices to contiguous execution blocks.
    long int *chunk_buffer;
    BHAH_MALLOC_PINNED(chunk_buffer, sizeof(long int) * BUNDLE_CAPACITY);

    // Bridge array for the state vector $f^\mu$.
    double *f_bridge;
    // Bridge array for the first derivative of the state vector.
    double *f_p_bridge;
    // Bridge array for the second derivative of the state vector.
    double *f_p_p_bridge;
    // Bridge array storing the affine parameter $\lambda$.
    double *affine_bridge;
    // Bridge array storing the current integration step size $h$.
    double *h_bridge;
    // Bridge array storing the current trajectory termination status.
    termination_type_t *status_bridge;
    // Bridge array counting the number of step-size rejections.
    int *retries_bridge;
    // Bridge array tracking if the photon was previously on the positive side of the window.
    bool *on_pos_window_prev_bridge;
    // Bridge array tracking if the photon was previously on the positive side of the source plane.
    bool *on_pos_source_prev_bridge;

    // Allocate Bridge arrays in Host Pinned Memory to allow fast chunked transfers to the GPU.
    BHAH_MALLOC_PINNED(f_bridge, sizeof(double) * 9 * BUNDLE_CAPACITY);
    BHAH_MALLOC_PINNED(f_p_bridge, sizeof(double) * 9 * BUNDLE_CAPACITY);
    BHAH_MALLOC_PINNED(f_p_p_bridge, sizeof(double) * 9 * BUNDLE_CAPACITY);
    BHAH_MALLOC_PINNED(affine_bridge, sizeof(double) * BUNDLE_CAPACITY);
    BHAH_MALLOC_PINNED(h_bridge, sizeof(double) * BUNDLE_CAPACITY);
    BHAH_MALLOC_PINNED(status_bridge, sizeof(termination_type_t) * BUNDLE_CAPACITY);
    BHAH_MALLOC_PINNED(retries_bridge, sizeof(int) * BUNDLE_CAPACITY);
    BHAH_MALLOC_PINNED(on_pos_window_prev_bridge, sizeof(bool) * BUNDLE_CAPACITY);
    BHAH_MALLOC_PINNED(on_pos_source_prev_bridge, sizeof(bool) * BUNDLE_CAPACITY);

    // Host-to-Device transfer: Maps the global spacetime constants to the device cache to ensure zero-latency read access.
    cudaMemcpyToSymbol(d_commondata, commondata, sizeof(commondata_struct));

    // VRAM scratchpad tracking the current state vector $f^\mu$ bounding the RKF45 step.
    double *d_f_bundle;
    // VRAM scratchpad locking the anchor state vector $f_{{start}}$ to calculate the final stage update.
    double *d_f_start_bundle;
    // VRAM scratchpad tracking the intermediate cumulative RKF45 stage updates.
    double *d_f_temp_bundle;
    // VRAM scratchpad tracking the history state $f^\mu_{{n-1}}$ for geometric intersection detection.
    double *d_f_prev_bundle;
    // VRAM scratchpad tracking the history state $f^\mu_{{n-2}}$ for geometric intersection detection.
    double *d_f_pre_prev_bundle;

    // Allocate 1D VRAM Scratchpad arrays ensuring strict adherence to the hardware bounds.
    BHAH_MALLOC_DEVICE(d_f_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_f_start_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_f_temp_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_f_prev_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_f_pre_prev_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY);

    // VRAM scratchpad persisting the symmetric metric tensor $g_{{\mu\nu}}$.
    double *d_metric_bundle;
    // VRAM scratchpad persisting the Christoffel symbols $\Gamma^\alpha_{{\beta\gamma}}$.
    double *d_connection_bundle;

    BHAH_MALLOC_DEVICE(d_metric_bundle, sizeof(double) * 10 * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_connection_bundle, sizeof(double) * 40 * BUNDLE_CAPACITY);

    // Derivative tensor storing $\dot{{f}}$ across all 6 intermediate RKF45 stages.
    double *d_k_bundle;
    
    // Allocate VRAM for the 6-stage $k$-bundle to strictly prevent register bleeding.
    BHAH_MALLOC_DEVICE(d_k_bundle, sizeof(double) * 6 * 9 * BUNDLE_CAPACITY);

    // VRAM array regulating active integration step sizing $h$.
    double *d_h;
    // VRAM array regulating total affine parameter progress $\lambda$.
    double *d_affine;
    // VRAM array holding the current trajectory status limits.
    termination_type_t *d_status;
    // VRAM array tracking sequential error rejections per photon.
    int *d_retries;
    // VRAM array flagging the previous observer window boundary side.
    bool *d_on_pos_window_prev;
    // VRAM array flagging the previous source emission boundary side.
    bool *d_on_pos_source_prev;

    BHAH_MALLOC_DEVICE(d_h, sizeof(double) * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_affine, sizeof(double) * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_status, sizeof(termination_type_t) * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_retries, sizeof(int) * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_on_pos_window_prev, sizeof(bool) * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_on_pos_source_prev, sizeof(bool) * BUNDLE_CAPACITY);

    // Final intersection blueprint results buffer mapped directly to VRAM.
    blueprint_data_t *d_results_buffer;
    BHAH_MALLOC_DEVICE(d_results_buffer, sizeof(blueprint_data_t) * num_rays);

    // Lock-free temporal TimeSlotManager confined entirely to the Host CPU context.
    TimeSlotManager tsm;
    slot_manager_init(&tsm, commondata->slot_manager_t_min, commondata->t_start + 1.0, commondata->slot_manager_delta_t, num_rays);


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
    set_initial_conditions_kernel_{spacetime_name}(commondata, num_rays, &all_photons_host, window_center_out, n_x_out, n_y_out, n_z_out);

    // Total integer calculation defining total iterative blocks required to process all photon indices.
    long int num_batches = (num_rays + BUNDLE_CAPACITY - 1) / BUNDLE_CAPACITY;

    // Loop iterator for evaluating the initialization constraint across sequential blocks.
    long int init_batch;
    for (init_batch = 0; init_batch < num_batches; ++init_batch) {{
        // The master index offset for the current batch initialization.
        long int start_idx = init_batch * BUNDLE_CAPACITY;
        // The constrained integer size of the active chunk bounded by the hardware capacity.
        long int chunk_size = MIN((long int)BUNDLE_CAPACITY, num_rays - start_idx);

        // Loop index iterating over the specific initialization batch elements.
        int init_i;
        for (init_i = 0; init_i < chunk_size; ++init_i) {{
            // The calculated global long integer index corresponding to the specific ray in the master SoA.
            long int master_idx = start_idx + init_i;
            // Loop index over the 9 tensor components of the state vector.
            int init_k;
            for (init_k = 0; init_k < 9; ++init_k) {{
                f_bridge[init_k * BUNDLE_CAPACITY + init_i] = all_photons_host.f[init_k * num_rays + master_idx];
            }}
        }}

        // Host-to-Device transfer: Pushes packed initial states to VRAM to exploit parallel geometric evaluations.
        cudaMemcpy(d_f_bundle, f_bridge, sizeof(double) * 9 * BUNDLE_CAPACITY, cudaMemcpyHostToDevice);

        // Kernel Launch: Calculate $g_{{\mu\nu}}$ required for the Hamiltonian constraint. Thread ID maps to photon index.
        interpolation_kernel_{spacetime_name}(d_f_bundle, d_metric_bundle, NULL, chunk_size);
        
        // Kernel Launch: Solves the constraint $p_\mu p^\mu = 0$ to find the temporal momentum $p_t$.
        p0_reverse_kernel(d_f_bundle, d_metric_bundle, chunk_size);

        // Device-to-Host transfer: Retrieves the mathematically constrained state vector back to CPU RAM.
        cudaMemcpy(f_bridge, d_f_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY, cudaMemcpyDeviceToHost);

        // Loop index for gathering post-constraint logic and re-scattering to the master array.
        int gather_i;
        for (gather_i = 0; gather_i < chunk_size; ++gather_i) {{
            // The calculated global index for the scatter update back to the main SoA.
            long int master_idx = start_idx + gather_i;
            // Loop index over the 9 tensor components to ensure complete state transfer.
            int gather_k;
            for (gather_k = 0; gather_k < 9; ++gather_k) {{
                all_photons_host.f[gather_k * num_rays + master_idx] = f_bridge[gather_k * BUNDLE_CAPACITY + gather_i];
            }}
        }}
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
        all_photons_host.on_positive_side_of_window_prev[sync_i] = true;
        all_photons_host.on_positive_side_of_source_prev[sync_i] = true;
        all_photons_host.status[sync_i] = ACTIVE;
        all_photons_host.h[sync_i] = commondata->numerical_initial_h;
        all_photons_host.affine_param[sync_i] = 0.0;
        all_photons_host.rejection_retries[sync_i] = 0;

        // Integer representing the assigned temporal slot for the current photon.
        int s_idx = slot_get_index(&tsm, all_photons_host.f[sync_i]);
        if (s_idx != -1) {{
            slot_add_photon(&tsm, s_idx, sync_i);
        }}
    }}


    // --- 3. TEMPORAL LOOP (The Engine) ---
    
    // Outer loop iterator for the physical time bins.
    int slot_idx;
    for (slot_idx = 0; slot_idx < tsm.num_slots; ++slot_idx) {{
        while (tsm.slot_counts[slot_idx] > 0) {{
            
            // Evaluated bounded integer size corresponding to active trajectories in this specific time bin.
            long int chunk_size = MIN((long int)BUNDLE_CAPACITY, tsm.slot_counts[slot_idx]);
            slot_remove_chunk(&tsm, slot_idx, chunk_buffer, chunk_size);

            // Inner loop iterator for populating contiguous memory bridges.
            int bridge_i;
            for (bridge_i = 0; bridge_i < chunk_size; ++bridge_i) {{
                // Master index retrieved from the extraction buffer array.
                long int m_idx = chunk_buffer[bridge_i];
                // Component loop iterator mapping values directly into the bridge transfer blocks.
                int c_k;
                for (c_k = 0; c_k < 9; ++c_k) {{
                    f_bridge[c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f[c_k * num_rays + m_idx];
                    f_p_bridge[c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p[c_k * num_rays + m_idx];
                    f_p_p_bridge[c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p_p[c_k * num_rays + m_idx];
                }}
                h_bridge[bridge_i] = all_photons_host.h[m_idx];
                status_bridge[bridge_i] = all_photons_host.status[m_idx];
                retries_bridge[bridge_i] = all_photons_host.rejection_retries[m_idx];
                affine_bridge[bridge_i] = all_photons_host.affine_param[m_idx];
                on_pos_window_prev_bridge[bridge_i] = all_photons_host.on_positive_side_of_window_prev[m_idx];
                on_pos_source_prev_bridge[bridge_i] = all_photons_host.on_positive_side_of_source_prev[m_idx];
            }}

            // Host-to-Device transfer: Pushes packed $f^\mu$ execution chunk to VRAM scratchpad.
            cudaMemcpy(d_f_bundle, f_bridge, sizeof(double) * 9 * BUNDLE_CAPACITY, cudaMemcpyHostToDevice);
            // Host-to-Device transfer: Pushes packed $f_p$ execution chunk to VRAM scratchpad.
            cudaMemcpy(d_f_prev_bundle, f_p_bridge, sizeof(double) * 9 * BUNDLE_CAPACITY, cudaMemcpyHostToDevice);
            // Host-to-Device transfer: Pushes packed $f_p_p$ execution chunk to VRAM scratchpad.
            cudaMemcpy(d_f_pre_prev_bundle, f_p_p_bridge, sizeof(double) * 9 * BUNDLE_CAPACITY, cudaMemcpyHostToDevice);
            // Host-to-Device transfer: Pushes packed integration step sizes $h$ to VRAM scratchpad.
            cudaMemcpy(d_h, h_bridge, sizeof(double) * BUNDLE_CAPACITY, cudaMemcpyHostToDevice);
            // Host-to-Device transfer: Pushes packed status flags to VRAM scratchpad.
            cudaMemcpy(d_status, status_bridge, sizeof(termination_type_t) * BUNDLE_CAPACITY, cudaMemcpyHostToDevice);
            // Host-to-Device transfer: Pushes packed rejection retries to VRAM scratchpad.
            cudaMemcpy(d_retries, retries_bridge, sizeof(int) * BUNDLE_CAPACITY, cudaMemcpyHostToDevice);
            // Host-to-Device transfer: Pushes packed affine progress $\lambda$ to VRAM scratchpad.
            cudaMemcpy(d_affine, affine_bridge, sizeof(double) * BUNDLE_CAPACITY, cudaMemcpyHostToDevice);
            // Host-to-Device transfer: Pushes packed window history side flags.
            cudaMemcpy(d_on_pos_window_prev, on_pos_window_prev_bridge, sizeof(bool) * BUNDLE_CAPACITY, cudaMemcpyHostToDevice);
            // Host-to-Device transfer: Pushes packed source history side flags.
            cudaMemcpy(d_on_pos_source_prev, on_pos_source_prev_bridge, sizeof(bool) * BUNDLE_CAPACITY, cudaMemcpyHostToDevice);

            // Device-to-Device transfer: Baselines the starting state $f_{{start}}$ to protect base vectors during RKF45 intermediate updates.
            cudaMemcpy(d_f_start_bundle, d_f_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY, cudaMemcpyDeviceToDevice);
            // Device-to-Device transfer: Mirrors initial state to the $f_{{temp}}$ integration accumulator.
            cudaMemcpy(d_f_temp_bundle, d_f_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY, cudaMemcpyDeviceToDevice);

            // Iterator for tracking the active RKF45 intermediate step mapping.
            int stage;
            for (stage = 1; stage <= 6; ++stage) {{
                // Kernel Launch: Evaluate $g_{{\mu\nu}}$ and $\Gamma^\alpha_{{\beta\gamma}}$ based on current $f_{{temp}}$.
                interpolation_kernel_{spacetime_name}(d_f_temp_bundle, d_metric_bundle, d_connection_bundle, chunk_size);
                
                // Kernel Launch: Compute the differential tensor $\dot{{f}}$ mapping strictly to Global VRAM $k$-bundles.
                calculate_ode_rhs_kernel(d_f_temp_bundle, d_metric_bundle, d_connection_bundle, d_k_bundle, stage, chunk_size);
                
                // Kernel Launch: Accumulate updates generating the intermediate vector $f_{{temp}}$ utilizing Butcher tableau coefficients.
                rkf45_stage_update(d_f_start_bundle, d_k_bundle, d_h, stage, chunk_size, d_f_temp_bundle);
            }}

            // Kernel Launch: Adjust step-size $h$ based on explicit local truncation limits bounding the trajectory errors.
            rkf45_finalize_and_control(d_f_bundle, d_f_start_bundle, d_k_bundle, d_h, d_status, d_affine, d_retries, chunk_size);
            
            // Kernel Launch: Check geometric intersection events across local orthonormal plane coordinates mapping spatial geometry.
            event_detection_manager_kernel(d_f_bundle, d_f_prev_bundle, d_f_pre_prev_bundle, d_results_buffer, d_status, d_on_pos_window_prev, d_on_pos_source_prev, chunk_size);

            // Device-to-Host transfer: Retrieves updated states $f^\mu$ back to the CPU Bridge arrays.
            cudaMemcpy(f_bridge, d_f_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY, cudaMemcpyDeviceToHost);
            // Device-to-Host transfer: Retrieves updated $f_p$ histories.
            cudaMemcpy(f_p_bridge, d_f_prev_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY, cudaMemcpyDeviceToHost);
            // Device-to-Host transfer: Retrieves updated $f_p_p$ histories.
            cudaMemcpy(f_p_p_bridge, d_f_pre_prev_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY, cudaMemcpyDeviceToHost);
            // Device-to-Host transfer: Retrieves modified step sizes $h$.
            cudaMemcpy(h_bridge, d_h, sizeof(double) * BUNDLE_CAPACITY, cudaMemcpyDeviceToHost);
            // Device-to-Host transfer: Retrieves updated termination flags.
            cudaMemcpy(status_bridge, d_status, sizeof(termination_type_t) * BUNDLE_CAPACITY, cudaMemcpyDeviceToHost);
            // Device-to-Host transfer: Retrieves current rejection retries.
            cudaMemcpy(retries_bridge, d_retries, sizeof(int) * BUNDLE_CAPACITY, cudaMemcpyDeviceToHost);
            // Device-to-Host transfer: Retrieves adapted affine parameter progression $\lambda$.
            cudaMemcpy(affine_bridge, d_affine, sizeof(double) * BUNDLE_CAPACITY, cudaMemcpyDeviceToHost);
            // Device-to-Host transfer: Retrieves persistent window boundary flags.
            cudaMemcpy(on_pos_window_prev_bridge, d_on_pos_window_prev, sizeof(bool) * BUNDLE_CAPACITY, cudaMemcpyDeviceToHost);
            // Device-to-Host transfer: Retrieves persistent source boundary flags.
            cudaMemcpy(on_pos_source_prev_bridge, d_on_pos_source_prev, sizeof(bool) * BUNDLE_CAPACITY, cudaMemcpyDeviceToHost);

            // Loop iterator traversing the finalized trajectory bundle to process statuses.
            int fin_i;
            for (fin_i = 0; fin_i < chunk_size; ++fin_i) {{
                // Master index mapped from the bounded extraction array.
                long int m_idx = chunk_buffer[fin_i];
                // Component loop iterator for unpacking bridge data.
                int fin_k;
                for (fin_k = 0; fin_k < 9; ++fin_k) {{
                    all_photons_host.f[fin_k * num_rays + m_idx] = f_bridge[fin_k * BUNDLE_CAPACITY + fin_i];
                    all_photons_host.f_p[fin_k * num_rays + m_idx] = f_p_bridge[fin_k * BUNDLE_CAPACITY + fin_i];
                    all_photons_host.f_p_p[fin_k * num_rays + m_idx] = f_p_p_bridge[fin_k * BUNDLE_CAPACITY + fin_i];
                }}
                all_photons_host.h[m_idx] = h_bridge[fin_i];
                all_photons_host.status[m_idx] = status_bridge[fin_i];
                all_photons_host.rejection_retries[m_idx] = retries_bridge[fin_i];
                all_photons_host.affine_param[m_idx] = affine_bridge[fin_i];
                all_photons_host.on_positive_side_of_window_prev[m_idx] = on_pos_window_prev_bridge[fin_i];
                all_photons_host.on_positive_side_of_source_prev[m_idx] = on_pos_source_prev_bridge[fin_i];

                if (status_bridge[fin_i] == ACTIVE) {{
                    // Extract $f^0$ (coordinate time) to evaluate the next discrete temporal bin.
                    int next_s_idx = slot_get_index(&tsm, all_photons_host.f[m_idx]); 
                    if (next_s_idx != -1) {{
                        slot_add_photon(&tsm, next_s_idx, m_idx);
                    }}
                }} else if (status_bridge[fin_i] == FAILURE_RKF45_REJECTION_LIMIT) {{
                    // Re-add to current bin to attempt integration with an adapted step-size scalar.
                    slot_add_photon(&tsm, slot_idx, m_idx);
                }}
            }}
        }}
    }}

    // --- 4. CLEANUP ---
    
    // Device-to-Host transfer: Extracts validated device-native blueprints mapped into the host pointer.
    cudaMemcpy(results_buffer, d_results_buffer, sizeof(blueprint_data_t) * num_rays, cudaMemcpyDeviceToHost);

    // Host Memory Free: Purges the primary Host array states.
    BHAH_FREE_PINNED(all_photons_host.f);
    BHAH_FREE_PINNED(all_photons_host.f_p);
    BHAH_FREE_PINNED(all_photons_host.f_p_p);
    BHAH_FREE_PINNED(all_photons_host.affine_param);
    BHAH_FREE_PINNED(all_photons_host.h);
    BHAH_FREE_PINNED(all_photons_host.status);
    BHAH_FREE_PINNED(all_photons_host.rejection_retries);
    BHAH_FREE_PINNED(all_photons_host.on_positive_side_of_window_prev);
    BHAH_FREE_PINNED(all_photons_host.on_positive_side_of_source_prev);
    
    // Host Memory Free: Purges bridge components supporting scatter logic.
    BHAH_FREE_PINNED(chunk_buffer);
    BHAH_FREE_PINNED(f_bridge);
    BHAH_FREE_PINNED(f_p_bridge);
    BHAH_FREE_PINNED(f_p_p_bridge);
    BHAH_FREE_PINNED(affine_bridge);
    BHAH_FREE_PINNED(h_bridge);
    BHAH_FREE_PINNED(status_bridge);
    BHAH_FREE_PINNED(retries_bridge);
    BHAH_FREE_PINNED(on_pos_window_prev_bridge);
    BHAH_FREE_PINNED(on_pos_source_prev_bridge);

    // Device Memory Free: Purges VRAM operational pipeline scratchpads.
    BHAH_FREE_DEVICE(d_f_bundle);
    BHAH_FREE_DEVICE(d_f_start_bundle);
    BHAH_FREE_DEVICE(d_f_temp_bundle);
    BHAH_FREE_DEVICE(d_f_prev_bundle);
    BHAH_FREE_DEVICE(d_f_pre_prev_bundle);
    BHAH_FREE_DEVICE(d_metric_bundle);
    BHAH_FREE_DEVICE(d_connection_bundle);
    BHAH_FREE_DEVICE(d_k_bundle);
    BHAH_FREE_DEVICE(d_h);
    BHAH_FREE_DEVICE(d_affine);
    BHAH_FREE_DEVICE(d_status);
    BHAH_FREE_DEVICE(d_retries);
    BHAH_FREE_DEVICE(d_on_pos_window_prev);
    BHAH_FREE_DEVICE(d_on_pos_source_prev);
    BHAH_FREE_DEVICE(d_results_buffer);

    // Explicitly free the temporal sorting struct mapping the execution grid.
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