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

    prefunc = r"""
    #define TIME_BLOCK(name, block) \
        { \
            cudaEvent_t start, stop; \
            cudaEventCreate(&start); cudaEventCreate(&stop); \
            cudaEventRecord(start); \
            block; \
            cudaEventRecord(stop); \
            cudaEventSynchronize(stop); \
            float ms = 0; \
            cudaEventElapsedTime(&ms, start, stop); \
            printf("  [TIMER] %-30s: %8.3f ms\n", name, ms); \
            cudaEventDestroy(start); cudaEventDestroy(stop); \
        }
    """

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

    // Host-to-Device transfer allocation: Pinned memory tracking the history step $\lambda_{{n-1}}$.
    BHAH_MALLOC_PINNED(all_photons_host.affine_param_p, sizeof(double) * num_rays);
    // Host-to-Device transfer allocation: Pinned memory tracking the history step $\lambda_{{n-2}}$.
    BHAH_MALLOC_PINNED(all_photons_host.affine_param_p_p, sizeof(double) * num_rays);
    // Host-to-Device transfer allocation: Pinned memory locking the observer window intersection.
    BHAH_MALLOC_PINNED(all_photons_host.window_event_found, sizeof(bool) * num_rays);
    // Host-to-Device transfer allocation: Pinned memory locking the source emission plane intersection.
    BHAH_MALLOC_PINNED(all_photons_host.source_event_found, sizeof(bool) * num_rays);

    // Bridge array storing the historical affine parameter $\lambda_{{n-1}}$ for chunked transfers.
    double *affine_p_bridge;
    // Bridge array storing the historical affine parameter $\lambda_{{n-2}}$ for chunked transfers.
    double *affine_p_p_bridge;
    // Bridge array tracking if the window event was previously found.
    bool *window_event_found_bridge;
    // Bridge array tracking if the source event was previously found.
    bool *source_event_found_bridge;

    // Allocate Bridge arrays in Host Pinned Memory to allow fast chunked transfers to the GPU.
    BHAH_MALLOC_PINNED(affine_p_bridge, sizeof(double) * BUNDLE_CAPACITY);
    BHAH_MALLOC_PINNED(affine_p_p_bridge, sizeof(double) * BUNDLE_CAPACITY);
    BHAH_MALLOC_PINNED(window_event_found_bridge, sizeof(bool) * BUNDLE_CAPACITY);
    BHAH_MALLOC_PINNED(source_event_found_bridge, sizeof(bool) * BUNDLE_CAPACITY);

    // VRAM array tracking historical affine parameter $\lambda_{{n-1}}$.
    double *d_affine_prev;
    // VRAM array tracking historical affine parameter $\lambda_{{n-2}}$.
    double *d_affine_pre_prev;
    // VRAM array guarding the window intersection coordinates from multi-trigger overwrites.
    bool *d_window_event_found;
    // VRAM array guarding the source intersection coordinates from multi-trigger overwrites.
    bool *d_source_event_found;
    // VRAM array carrying the absolute master indices $m_{{idx}}$ mapping the execution chunk.
    long int *d_chunk_buffer;

    BHAH_MALLOC_DEVICE(d_affine_prev, sizeof(double) * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_affine_pre_prev, sizeof(double) * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_window_event_found, sizeof(bool) * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_source_event_found, sizeof(bool) * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_chunk_buffer, sizeof(long int) * BUNDLE_CAPACITY);

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
    set_initial_conditions_kernel_{spacetime_name}(commondata, num_rays, &all_photons_host, window_center_out, n_x_out, n_y_out, n_z_out);



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    // --- DIAGNOSTIC PROBE: DETAILED INITIAL POSITION & PLACEHOLDER ALIGNMENT ---
    // Algorithmic Step: Scans the master Host SoA immediately following the initialization kernel call.
    // Hardware Justification: This architectural step verifies that the Host-side $set\_initial\_conditions$
    // routine has correctly populated the starting coordinates $x^\mu$ and zeroed the $p_t$ and $\lambda$ components.
    // Identifying coordinate misalignment here isolates initialization logic failures from downstream VRAM integration errors.
    long int init_mismatch_count = 0; // Tracks the total number of trajectories with unaligned initial states.
    long int mismatch_t = 0;          // Tracks failures in coordinate time $t$.
    long int mismatch_x = 0;          // Tracks failures in spatial coordinate $x$.
    long int mismatch_y = 0;          // Tracks failures in spatial coordinate $y$.
    long int mismatch_z = 0;          // Tracks failures in spatial coordinate $z$.
    long int mismatch_pt = 0;         // Tracks failures in temporal momentum $p_t$.
    long int mismatch_lam = 0;        // Tracks failures in affine parameter $\lambda$.

    // Loop iterator traversing the global ray count to validate Host-side coordinate hydration.
    for (long int p = 0; p < num_rays; p++) {{
        // Retrieving the initialized coordinate time $t$ and spatial coordinates $x^i$ from the master SoA.
        const double t_check = all_photons_host.f[0 * num_rays + p];
        const double x_check = all_photons_host.f[1 * num_rays + p];
        const double y_check = all_photons_host.f[2 * num_rays + p];
        const double z_check = all_photons_host.f[3 * num_rays + p];

        // Retrieving the placeholders for temporal momentum $p_t$ (index 4) and the affine parameter $\lambda$ (index 8).
        const double pt_check = all_photons_host.f[4 * num_rays + p];
        const double lam_check = all_photons_host.f[8 * num_rays + p];

        // Evaluates the physical state against the prescribed camera coordinates and numerical boundaries.
        // A tolerance of $10^{-10}$ is utilized to accommodate floating-point variance in basis transformations.
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

    // Reporting the fractional failure rate if any trajectories deviate from the intended camera parameters.
    if (init_mismatch_count > 0) {{
        const double mismatch_percent = ((double)init_mismatch_count / (double)num_rays) * 100.0;
        printf("[DIAGNOSTIC] Initialization Alignment Check: %ld out of %ld rays (%.2f%%) fail coordinate/placeholder validation.\n", 
            init_mismatch_count, num_rays, mismatch_percent);
        printf("             Component Breakdown of Failures:\n");
        printf("               $t$       (Target: %.2f) : %ld rays (%.2f%%)\n", commondata->t_start, mismatch_t, ((double)mismatch_t / num_rays) * 100.0);
        printf("               $x$       (Target: %.2f) : %ld rays (%.2f%%)\n", commondata->camera_pos_x, mismatch_x, ((double)mismatch_x / num_rays) * 100.0);
        printf("               $y$       (Target: %.2f) : %ld rays (%.2f%%)\n", commondata->camera_pos_y, mismatch_y, ((double)mismatch_y / num_rays) * 100.0);
        printf("               $z$       (Target: %.2f) : %ld rays (%.2f%%)\n", commondata->camera_pos_z, mismatch_z, ((double)mismatch_z / num_rays) * 100.0);
        printf("               $p_t$     (Target: 0.0)  : %ld rays (%.2f%%)\n", mismatch_pt, ((double)mismatch_pt / num_rays) * 100.0);
        printf("               $\\lambda$ (Target: 0.0)  : %ld rays (%.2f%%)\n", mismatch_lam, ((double)mismatch_lam / num_rays) * 100.0);
    }}


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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


        // --- DIAGNOSTIC PROBE: VRAM METRIC INTEGRITY CHECK ---
        // Algorithmic Step: Syncs the metric scratchpad to the Host to verify geometric stability.
        // Hardware Justification: This architectural step identifies coordinate singularities in VRAM 
        // before the Hamiltonian constraint $p_\mu p^\mu = 0$ is evaluated by the $p0\_reverse$ kernel.
        double *metric_diag_bridge; // Temporary bridge for metric tensor components.
        BHAH_MALLOC_PINNED(metric_diag_bridge, sizeof(double) * 10 * BUNDLE_CAPACITY); 

        // Device-to-Host transfer: Retrieves the symmetric metric tensor $g_{{\mu\nu}}$ for inspection.
        cudaMemcpy(metric_diag_bridge, d_metric_bundle, sizeof(double) * 10 * BUNDLE_CAPACITY, cudaMemcpyDeviceToHost); 

        long int metric_nan_count = 0; // Counter for non-finite metric components.
        int m_diag_i;
        for (m_diag_i = 0; m_diag_i < chunk_size; ++m_diag_i) {{
        bool m_has_nan = false;
        int m_diag_k;
        for (m_diag_k = 0; m_diag_k < 10; ++m_diag_k) {{
            // Evaluates $g_{{\mu\nu}}$ for $NaN$ or $Inf$ signifying a coordinate breakdown.
            if (isnan(metric_diag_bridge[m_diag_k * BUNDLE_CAPACITY + m_diag_i]) || 
                isinf(metric_diag_bridge[m_diag_k * BUNDLE_CAPACITY + m_diag_i])) {{
            m_has_nan = true;
            break;
            }}
        }}
        if (m_has_nan) {{
            metric_nan_count++;
        }}
        }}

        if (metric_nan_count > 0) {{
        double m_nan_percentage = ((double)metric_nan_count / (double)chunk_size) * 100.0;
        printf("[DIAGNOSTIC] Init Batch %ld: %ld rays (%.2f%%) have invalid Metric G_mu_nu before p_t solve.\n", 
                init_batch, metric_nan_count, m_nan_percentage);
        }}
        BHAH_FREE_PINNED(metric_diag_bridge); 


        
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

        // --- DIAGNOSTIC PROBE: INITIALIZATION NaN DETECTOR ---
      // Scans the newly constraint-solved batch for non-finite values.
      // Hardware Justification: This architectural step executes immediately after the Device-to-Host transfer
      // to verify the integrity of the initial state $f^\mu$ before committing it to the integration engine.
      long int nan_count = 0; // Tracks the total number of corrupted photon trajectories in the current batch.

      // Loop iterator scanning the bounded execution chunk.
      int diag_i;
      for (diag_i = 0; diag_i < chunk_size; ++diag_i) {{
        bool has_nan = false; // Persistent flag indicating a non-finite tensor component was identified.
        
        // Component loop iterator mapping values across the 9 tensor dimensions.
        int diag_k;
        for (diag_k = 0; diag_k < 9; ++diag_k) {{
          // Evaluates the bridged state vector $f^\mu$ for standard IEEE-754 NaN representation.
          if (isnan(f_bridge[diag_k * BUNDLE_CAPACITY + diag_i])) {{
            has_nan = true;
            break; // Terminates evaluation early to minimize redundant CPU cache line fetches.
          }}
        }}
        
        if (has_nan) {{
          nan_count++; // Increments the global error tracking counter.
        }}
      }}

      // Evaluates if any trajectories failed the Hamiltonian constraint solving phase.
      if (nan_count > 0) {{
        // Double precision scalar representing the fractional failure rate of the Hamiltonian constraint.
        double nan_percentage = ((double)nan_count / (double)chunk_size) * 100.0;
        printf("[DIAGNOSTIC] Init Batch %ld: %ld out of %ld rays (%.2f%%) contain NaN in state f^mu after p_t solve.\n", init_batch, nan_count, (long int)chunk_size, nan_percentage);
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
    int slot_idx;
    for (slot_idx = tsm.num_slots - 1; slot_idx >= 0; --slot_idx) {{
        // Evaluates the early exit condition to terminate the temporal engine if all geometric trajectories have concluded.
        if (total_active_photons <= 0) {{
        break;
        }}

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

                // Pack affine parameter progression and event guards into physically contiguous buffers for DMA efficiency.
                affine_p_bridge[bridge_i] = all_photons_host.affine_param_p[m_idx];
                affine_p_p_bridge[bridge_i] = all_photons_host.affine_param_p_p[m_idx];
                window_event_found_bridge[bridge_i] = all_photons_host.window_event_found[m_idx];
                source_event_found_bridge[bridge_i] = all_photons_host.source_event_found[m_idx];
            }}

            // --- HOST-TO-DEVICE PAYLOAD TRANSFER ---
            // Transfers the active execution chunk of the state vector $f^\mu$ and its history to the VRAM scratchpad.

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////
            TIME_BLOCK("H2D Bridge Transfers", {{
                // Transfers the active execution chunk of the state vector $f^\mu$ and its history to the VRAM scratchpad.
                // Hardware Justification: Component-wise Host-to-Device transfers bounded by $chunk\_size$ minimize PCIe bus saturation, while $BUNDLE\_CAPACITY$ enforces the rigid SoA stride.
                for (int c_k = 0; c_k < 9; ++c_k) {{
                    // Host-to-Device transfer: Pushes the active trajectory state $f^\mu$ to VRAM to minimize PCIe overhead.
                    cudaMemcpy(d_f_bundle + c_k * BUNDLE_CAPACITY, f_bridge + c_k * BUNDLE_CAPACITY, sizeof(double) * chunk_size, cudaMemcpyHostToDevice);
                    // Host-to-Device transfer: Pushes the history state $f^\mu_{{n-1}}$ to VRAM for geometric intersection detection.
                    cudaMemcpy(d_f_prev_bundle + c_k * BUNDLE_CAPACITY, f_p_bridge + c_k * BUNDLE_CAPACITY, sizeof(double) * chunk_size, cudaMemcpyHostToDevice);
                    // Host-to-Device transfer: Pushes the history state $f^\mu_{{n-2}}$ to VRAM.
                    cudaMemcpy(d_f_pre_prev_bundle + c_k * BUNDLE_CAPACITY, f_p_p_bridge + c_k * BUNDLE_CAPACITY, sizeof(double) * chunk_size, cudaMemcpyHostToDevice);
                }}

                // Host-to-Device transfer: Pushes the scalar integration step sizes $h$ to VRAM to utilize high-bandwidth device memory.
                cudaMemcpy(d_h, h_bridge, sizeof(double) * chunk_size, cudaMemcpyHostToDevice);
                // Host-to-Device transfer: Pushes the numerical status flags to VRAM to track active vs. terminated states per thread.
                cudaMemcpy(d_status, status_bridge, sizeof(termination_type_t) * chunk_size, cudaMemcpyHostToDevice);
                // Host-to-Device transfer: Pushes the consecutive rejection counters to VRAM to evaluate RKF45 adaptive limits.
                cudaMemcpy(d_retries, retries_bridge, sizeof(int) * chunk_size, cudaMemcpyHostToDevice);
                // Host-to-Device transfer: Pushes the affine parameter $\lambda$ to VRAM to maintain integration progress.
                cudaMemcpy(d_affine, affine_bridge, sizeof(double) * chunk_size, cudaMemcpyHostToDevice);
                // Host-to-Device transfer: Pushes the persistent window boundary flags to VRAM to detect geometric sign changes.
                cudaMemcpy(d_on_pos_window_prev, on_pos_window_prev_bridge, sizeof(bool) * chunk_size, cudaMemcpyHostToDevice);
                // Host-to-Device transfer: Pushes the persistent source boundary flags to VRAM to detect emission plane intersections.
                cudaMemcpy(d_on_pos_source_prev, on_pos_source_prev_bridge, sizeof(bool) * chunk_size, cudaMemcpyHostToDevice);

                // Host-to-Device transfer: Pushes the discrete affine history buffers tracking $\lambda_{{n-1}}$ and $\lambda_{{n-2}}$.
                cudaMemcpy(d_affine_prev, affine_p_bridge, sizeof(double) * chunk_size, cudaMemcpyHostToDevice);
                cudaMemcpy(d_affine_pre_prev, affine_p_p_bridge, sizeof(double) * chunk_size, cudaMemcpyHostToDevice);
                
                // Host-to-Device transfer: Pushes the event guards to lock geometric intersections upon detection.
                cudaMemcpy(d_window_event_found, window_event_found_bridge, sizeof(bool) * chunk_size, cudaMemcpyHostToDevice);
                cudaMemcpy(d_source_event_found, source_event_found_bridge, sizeof(bool) * chunk_size, cudaMemcpyHostToDevice);

                // Host-to-Device transfer: Pushes the master index keys $m_{{idx}}$ to map thread IDs to absolute global indices.
                cudaMemcpy(d_chunk_buffer, chunk_buffer, sizeof(long int) * chunk_size, cudaMemcpyHostToDevice);
            }});
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // --- DEVICE-TO-DEVICE BASELINE SYNCHRONIZATION ---

            // --- DEVICE-TO-DEVICE BASELINE SYNCHRONIZATION ---
            // Duplicates the freshly loaded initial state $f^\mu$ into the persistent integrator scratchpads.
            // Hardware Justification: Device-to-Device copies ensure zero-state corruption is prevented entirely within the high-bandwidth VRAM domain.
            for (int c_k = 0; c_k < 9; ++c_k) {{
                // Device-to-Device transfer: Baselines the anchor state $f_{{start}}$ within the VRAM boundary.
                cudaMemcpy(d_f_start_bundle + c_k * BUNDLE_CAPACITY, d_f_bundle + c_k * BUNDLE_CAPACITY, sizeof(double) * chunk_size, cudaMemcpyDeviceToDevice);
                // Device-to-Device transfer: Initializes the intermediate accumulator $f_{{temp}}$ within the VRAM boundary.
                cudaMemcpy(d_f_temp_bundle + c_k * BUNDLE_CAPACITY, d_f_bundle + c_k * BUNDLE_CAPACITY, sizeof(double) * chunk_size, cudaMemcpyDeviceToDevice);
            }}

            // Iterator for tracking the active RKF45 intermediate step mapping.
            int stage;

           ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            TIME_BLOCK("RKF45 6-Stage Compute", {{
                for (stage = 1; stage <= 6; ++stage) {{
                    // Kernel Launch: Evaluate $g_{{\mu\nu}}$ and $\Gamma^\alpha_{{\beta\gamma}}$ based on current $f_{{temp}}$.
                    interpolation_kernel_{spacetime_name}(d_f_temp_bundle, d_metric_bundle, d_connection_bundle, chunk_size);
                    
                    // Kernel Launch: Compute the differential tensor $\dot{{f}}$ mapping strictly to Global VRAM $k$-bundles.
                    calculate_ode_rhs_kernel(d_f_temp_bundle, d_metric_bundle, d_connection_bundle, d_k_bundle, stage, chunk_size);
                    
                    // Kernel Launch: Accumulate updates generating the intermediate vector $f_{{temp}}$ utilizing Butcher tableau coefficients.
                    rkf45_stage_update(d_f_start_bundle, d_k_bundle, d_h, stage, chunk_size, d_f_temp_bundle);
                }}
            }});
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


            // Kernel Launch: Adjust step-size $h$ based on explicit local truncation limits bounding the trajectory errors.
            rkf45_finalize_and_control(d_f_bundle, d_f_start_bundle, d_k_bundle, d_h, d_status, d_affine, d_retries, chunk_size);
            

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            TIME_BLOCK("Event Detection Kernel", {{
                // Kernel Launch: Check geometric intersection events across local orthonormal plane coordinates mapping spatial geometry.
                event_detection_manager_kernel(d_f_bundle, d_f_prev_bundle, d_f_pre_prev_bundle, d_affine, d_affine_prev, d_affine_pre_prev, d_results_buffer, d_status, d_on_pos_window_prev, d_on_pos_source_prev, d_window_event_found, d_source_event_found, d_chunk_buffer, chunk_size);
            }});
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            // --- DEVICE-TO-HOST STATE RETRIEVAL ---
            // Extracts the updated trajectory states and termination flags from VRAM back to the Host Pinned bridges.
            // Hardware Justification: Device-to-Host transfers bounded by $chunk\_size$ prevent buffer overruns and minimize PCIe latency.
            for (int c_k = 0; c_k < 9; ++c_k) {{
                // Device-to-Host transfer: Retrieves the integrated state $f^\mu$ to Host memory.
                cudaMemcpy(f_bridge + c_k * BUNDLE_CAPACITY, d_f_bundle + c_k * BUNDLE_CAPACITY, sizeof(double) * chunk_size, cudaMemcpyDeviceToHost);
                // Device-to-Host transfer: Retrieves the history state $f^\mu_{{n-1}}$ to Host memory.
                cudaMemcpy(f_p_bridge + c_k * BUNDLE_CAPACITY, d_f_prev_bundle + c_k * BUNDLE_CAPACITY, sizeof(double) * chunk_size, cudaMemcpyDeviceToHost);
                // Device-to-Host transfer: Retrieves the history state $f^\mu_{{n-2}}$ to Host memory.
                cudaMemcpy(f_p_p_bridge + c_k * BUNDLE_CAPACITY, d_f_pre_prev_bundle + c_k * BUNDLE_CAPACITY, sizeof(double) * chunk_size, cudaMemcpyDeviceToHost);
            }}

            // Device-to-Host transfer: Retrieves the adapted step sizes $h$ to Host Pinned memory to seed the next integration bin.
            cudaMemcpy(h_bridge, d_h, sizeof(double) * chunk_size, cudaMemcpyDeviceToHost);
            // Device-to-Host transfer: Retrieves the numerical termination statuses to Host Pinned memory for TimeSlotManager routing.
            cudaMemcpy(status_bridge, d_status, sizeof(termination_type_t) * chunk_size, cudaMemcpyDeviceToHost);
            // Device-to-Host transfer: Retrieves the integration rejection counters to Host Pinned memory.
            cudaMemcpy(retries_bridge, d_retries, sizeof(int) * chunk_size, cudaMemcpyDeviceToHost);
            // Device-to-Host transfer: Retrieves the updated affine parameter $\lambda$ to Host Pinned memory.
            cudaMemcpy(affine_bridge, d_affine, sizeof(double) * chunk_size, cudaMemcpyDeviceToHost);
            // Device-to-Host transfer: Retrieves the persistent window boundary flags to Host Pinned memory for the next trajectory cycle.
            cudaMemcpy(on_pos_window_prev_bridge, d_on_pos_window_prev, sizeof(bool) * chunk_size, cudaMemcpyDeviceToHost);
            // Device-to-Host transfer: Retrieves the persistent source boundary flags to Host Pinned memory.
            cudaMemcpy(on_pos_source_prev_bridge, d_on_pos_source_prev, sizeof(bool) * chunk_size, cudaMemcpyDeviceToHost);

            // Device-to-Host transfer: Retrieves the tracked affine history vectors back to Host Pinned memory for routing.
            cudaMemcpy(affine_p_bridge, d_affine_prev, sizeof(double) * chunk_size, cudaMemcpyDeviceToHost);
            cudaMemcpy(affine_p_p_bridge, d_affine_pre_prev, sizeof(double) * chunk_size, cudaMemcpyDeviceToHost);

            // Device-to-Host transfer: Retrieves the updated event locks back to Host Pinned memory.
            cudaMemcpy(window_event_found_bridge, d_window_event_found, sizeof(bool) * chunk_size, cudaMemcpyDeviceToHost);
            cudaMemcpy(source_event_found_bridge, d_source_event_found, sizeof(bool) * chunk_size, cudaMemcpyDeviceToHost);

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

                // Scatter memory assignment extracting history $\lambda$ and locked statuses to the master array.
                all_photons_host.affine_param_p[m_idx] = affine_p_bridge[fin_i];
                all_photons_host.affine_param_p_p[m_idx] = affine_p_p_bridge[fin_i];
                all_photons_host.window_event_found[m_idx] = window_event_found_bridge[fin_i];
                all_photons_host.source_event_found[m_idx] = source_event_found_bridge[fin_i];

                if (status_bridge[fin_i] == ACTIVE) {{
                    // Extract $f^0$ (coordinate time) to evaluate the next discrete temporal bin.
                    int next_s_idx = slot_get_index(&tsm, all_photons_host.f[m_idx]);
                    if (next_s_idx != -1) {{
                        slot_add_photon(&tsm, next_s_idx, m_idx);
                    }} else {{
                        // Terminates photons that fall outside the TimeSlotManager bounds to prevent persistent ACTIVE state.
                        all_photons_host.status[m_idx] = FAILURE_T_MAX_EXCEEDED;
                        total_active_photons--;
                    }}
                    }} else if (status_bridge[fin_i] == REJECTED) {{
                    // Re-add to current bin to attempt integration with an adapted step-size scalar $h$.
                    slot_add_photon(&tsm, slot_idx, m_idx);
                    }} else {{
                    // Decrements the global counter as the photon has reached a terminal physical state or numerical failure limit.
                    total_active_photons--;
                }}
            }}
        }}
    }}

    // --- 4. CLEANUP & FINALIZATION ---
    
    // Device-to-Host transfer: Extracts validated device-native blueprints $b_i$ containing geometric plane intersections to Host memory.
    cudaMemcpy(results_buffer, d_results_buffer, sizeof(blueprint_data_t) * num_rays, cudaMemcpyDeviceToHost);

    // VRAM Memory Free: Purges heavy integration scratchpads (e.g., the 6-stage $k$-bundle and $\Gamma^\alpha_{{\beta\gamma}}$ array) to maximize memory headroom for the finalization pass.
    BHAH_FREE_DEVICE(d_k_bundle);
    BHAH_FREE_DEVICE(d_connection_bundle);
    BHAH_FREE_DEVICE(d_metric_bundle);
    BHAH_FREE_DEVICE(d_f_temp_bundle);
    BHAH_FREE_DEVICE(d_f_start_bundle);
    
    // Kernel Launch: Processes escaped photons intersecting the celestial sphere $r > r_{{escape}}$ and synchronizes final coordinate mappings in 32k streaming bundles.
    calculate_and_fill_blueprint_data_universal(&all_photons_host, num_rays, results_buffer);


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

    // Host Memory Free: Purges the primary Host array states $f^\mu$ and affine parameters $\lambda$.
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
    
    // Host Memory Free: Purges bridge components supporting scatter logic mapped to PCIe DMA transfers.
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
    BHAH_FREE_PINNED(affine_p_bridge);
    BHAH_FREE_PINNED(affine_p_p_bridge);
    BHAH_FREE_PINNED(window_event_found_bridge);
    BHAH_FREE_PINNED(source_event_found_bridge);

    // Device Memory Free: Purges remaining VRAM operational pipeline scratchpads.
    BHAH_FREE_DEVICE(d_f_bundle);
    BHAH_FREE_DEVICE(d_f_prev_bundle);
    BHAH_FREE_DEVICE(d_f_pre_prev_bundle);
    BHAH_FREE_DEVICE(d_h);
    BHAH_FREE_DEVICE(d_affine);
    BHAH_FREE_DEVICE(d_status);
    BHAH_FREE_DEVICE(d_retries);
    BHAH_FREE_DEVICE(d_on_pos_window_prev);
    BHAH_FREE_DEVICE(d_on_pos_source_prev);
    BHAH_FREE_DEVICE(d_affine_prev);
    BHAH_FREE_DEVICE(d_affine_pre_prev);
    BHAH_FREE_DEVICE(d_window_event_found);
    BHAH_FREE_DEVICE(d_source_event_found);
    BHAH_FREE_DEVICE(d_chunk_buffer);
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