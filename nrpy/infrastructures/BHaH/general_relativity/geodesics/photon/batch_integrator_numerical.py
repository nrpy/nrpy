"""
Orchestration module for the Project Singularity-Axiom numerical integration pipeline.

This module structures the high-level C orchestrator responsible for managing the
life cycle of photon trajectories $x^\mu$ in numerical spacetimes. It implements a
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
    Construct the Native CUDA orchestrator for the batched numerical integration pipeline.

    :param spacetime_name: The identifier for the spacetime metric (e.g., 'KerrSchild').
    :raises ValueError: If an unsupported spacetime metric identifier is provided.
    """
    
    # Core physics and numerical simulation parameters for the global spacetime struct.
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
        ["perform_conservation_check"],
        [True],
        commondata=True,
        add_to_parfile=True,
    )

    parallelization = par.parval_from_str("parallelization")

    includes = [
        "BHaH_defines.h", 
        "BHaH_function_prototypes.h"
    ]
    
    if parallelization == "cuda":
        includes.extend(["cuda_runtime.h", "cuda_intrinsics.h","BHaH_global_device_defines.h"])

    desc = r"""@brief Central Host-bound CPU orchestrator for the batched Split-Pipeline relativistic ray tracing loop.

    This function acts as the primary loop for evaluating photon geodesics $x^\mu$.
    It utilizes a TimeSlotManager to bin active rays by their physical coordinate time $t$. 
    The Split-Pipeline architecture maps mathematical tensors like $g_{\mu\nu}$ and $\Gamma^\alpha_{\beta\gamma}$ to memory scratchpads.
    Mapping tensors to memory respects the hardware limit per thread on modern architectures.

    @param commondata Struct containing global spacetime and numerical tolerances.
    @param num_rays Total number of photon trajectories to simulate.
    @param results_buffer Device array storing the final physical intersections."""

    cfunc_type = "void"

    name = "batch_integrator_numerical"

    params = "const commondata_struct *restrict commondata, long int num_rays, blueprint_data_t *restrict results_buffer"

    include_CodeParameters_h = True

    # --- DYNAMIC MACRO GENERATION ---
    malloc_pinned = "BHAH_MALLOC_PINNED" if parallelization == "cuda" else "BHAH_MALLOC"
    malloc_device = "BHAH_MALLOC_DEVICE" if parallelization == "cuda" else "BHAH_MALLOC"

    if parallelization == "cuda":
        stream_setup_str = """
        // Array of CUDA streams for asynchronous hardware orchestration.
        cudaStream_t streams[2];
        // Initializes the primary CUDA stream mapped to the first double-buffer context.
        cudaStreamCreate(&streams[0]);
        // Initializes the secondary CUDA stream mapped to the second double-buffer context.
        cudaStreamCreate(&streams[1]);

        // Host-to-Device transfer: Maps the global spacetime constants to the device cache to ensure zero-latency read access.
        cudaMemcpyToSymbol(d_commondata, commondata, sizeof(commondata_struct));"""
        pin_comment = "Host-to-Device transfer allocation: Pinned memory utilized to maximize PCIe DMA throughput for"
        dev_comment = "Device memory allocation: Dedicated VRAM array ensuring strict adherence to the hardware bounds for"
        bridge_alloc_comment = "Allocate Bridge arrays in Host Pinned Memory to allow fast, overlapped chunked transfers to the GPU."
        scratch_alloc_comment = "Allocate 1D VRAM Scratchpad arrays ensuring strict adherence to the hardware bounds."
    else:
        stream_setup_str = """
        // Streams not required for CPU execution.
    
        // Direct struct access utilized; explicit device symbol copying bypassed."""
        pin_comment = "Host memory allocation: Primary CPU RAM mapped for"
        dev_comment = "Host memory allocation: Primary CPU RAM mapped for"
        bridge_alloc_comment = "Allocate memory arrays in Host RAM for the structural bridge payloads."
        scratch_alloc_comment = "Allocate 1D Host Scratchpad arrays for temporal data staging."


    
    # This block configures the hardware-specific teardown and memory synchronization sequences.
    if parallelization == "cuda":
        results_memcpy = "cudaMemcpy(results_buffer, d_results_buffer, sizeof(blueprint_data_t) * num_rays, cudaMemcpyDeviceToHost);"
        calc_blueprint = "calculate_and_fill_blueprint_data_universal(&all_photons_host, num_rays, results_buffer, 0);"
        set_intitial_con = f" set_initial_conditions_kernel_{spacetime_name}(commondata, num_rays, &all_photons_host, window_center_out, n_x_out, n_y_out, n_z_out,0);"
        stream_destroy = "cudaStreamDestroy(streams[s]); // Hardware Resource Free: Purges the hardware stream execution context."
        free_device = "BHAH_FREE_DEVICE"
        free_pinned = "BHAH_FREE_PINNED"
    else:
        results_memcpy = "memcpy(results_buffer, d_results_buffer, sizeof(blueprint_data_t) * num_rays);"
        calc_blueprint = "calculate_and_fill_blueprint_data_universal(&all_photons_host, num_rays, results_buffer, 0);"
        set_intitial_con = f" set_initial_conditions_kernel_{spacetime_name}(commondata, num_rays, &all_photons_host, window_center_out, n_x_out, n_y_out, n_z_out);"
        stream_destroy = "// Hardware Justification: Stream destruction natively omitted for synchronous CPU execution."
        free_device = "BHAH_FREE"
        free_pinned = "BHAH_FREE"



    # Establish memory transfer protocols based on the target parallelization architecture.
    if parallelization == "cuda":
        # Python: Generates asynchronous PCIe transfer commands for CUDA architectures.
        def memcpy_async(dest, src, size, direction, stream):
            return f"cudaMemcpyAsync({dest}, {src}, {size}, {direction}, {stream});"
        
        # Generates hardware synchronization barriers for CUDA streams.
        def stream_sync(stream):
            return f"cudaStreamSynchronize({stream});"
        
        stream_arg = ", 0"
        free_pinned = "BHAH_FREE_PINNED"
        
        h2d_f_comment = r"// Hardware Justification: Asynchronously pushes strictly bounded initial states to VRAM on stream 0 to minimize PCIe latency."
        d2h_metric_comment = r"// Hardware Justification: Retrieves the symmetric metric tensor $g_{\mu\nu}$ to Host memory for inspection restricted to the active chunk."
        sync_metric_comment = r"// Hardware Justification: Hard synchronization barrier on stream 0 to inspect the payload prior to constraint solving."
        d2h_f_comment = r"// Hardware Justification: Retrieves strictly bounded constrained state vectors back to CPU RAM on stream 0."
        sync_f_comment = r"// Hardware Justification: Hard synchronization barrier on stream 0 to guarantee PCIe transfer completion before Host unpacking."
    else:
        # Generates synchronous Host memory copies for native OpenMP architectures.
        def memcpy_async(dest, src, size, direction, stream):
            return f"memcpy({dest}, {src}, {size});"
        
        # Omits stream synchronization as CPU execution is inherently synchronous.
        def stream_sync(stream):
            return r"// Hardware Justification: Stream synchronization is natively omitted for synchronous CPU execution."
            
        stream_arg = ", 0"
        free_pinned = "BHAH_FREE"
        
        h2d_f_comment = r"// Hardware Justification: Duplicates strictly bounded initial states to localized pipeline arrays using standard memory copies."
        d2h_metric_comment = r"// Hardware Justification: Duplicates the symmetric metric tensor $g_{\mu\nu}$ to diagnostic bridges via standard memory copies."
        sync_metric_comment = r"// Hardware Justification: Barrier omitted for synchronous Host execution prior to constraint solving."
        d2h_f_comment = r"// Hardware Justification: Retrieves localized constrained state vectors back to the primary Host bridge using standard memory copies."
        sync_f_comment = r"// Hardware Justification: Barrier omitted for synchronous Host execution before unpacking."

    # Define stream args for CUDA
    stream_arg_current = ", current" if parallelization == "cuda" else ", current"
    stream_arg_next = ", next" if parallelization == "cuda" else ", next"

    body = fr"""
    // --- 1. HOST & DEVICE ALLOCATION ---
    
    // The master host-side Structure of Arrays (SoA) tracking all photons $f^\mu$.
    PhotonStateSoA all_photons_host;

    // {pin_comment} the state vector $f^\mu$.
    {malloc_pinned}(all_photons_host.f, sizeof(double) * 9 * num_rays);
    // {pin_comment} the first derivative $\dot{{f}}^\mu$.
    {malloc_pinned}(all_photons_host.f_p, sizeof(double) * 9 * num_rays);
    // {pin_comment} the second derivative $\ddot{{f}}^\mu$.
    {malloc_pinned}(all_photons_host.f_p_p, sizeof(double) * 9 * num_rays);
    // {pin_comment} the physical affine parameter $\lambda$.
    {malloc_pinned}(all_photons_host.affine_param, sizeof(double) * num_rays);
    // {pin_comment} individual integration step sizes $h$.
    {malloc_pinned}(all_photons_host.h, sizeof(double) * num_rays);
    // {pin_comment} the trajectory termination status.
    {malloc_pinned}(all_photons_host.status, sizeof(termination_type_t) * num_rays);
    // {pin_comment} the number of step-size rejections.
    {malloc_pinned}(all_photons_host.rejection_retries, sizeof(int) * num_rays);
    // {pin_comment} the previous observer window boundary state.
    {malloc_pinned}(all_photons_host.on_positive_side_of_window_prev, sizeof(bool) * num_rays);
    // {pin_comment} the previous source emission boundary state.
    {malloc_pinned}(all_photons_host.on_positive_side_of_source_prev, sizeof(bool) * num_rays);
    // {pin_comment} the history step $\lambda_{{n-1}}$.
    {malloc_pinned}(all_photons_host.affine_param_p, sizeof(double) * num_rays);
    // {pin_comment} the history step $\lambda_{{n-2}}$.
    {malloc_pinned}(all_photons_host.affine_param_p_p, sizeof(double) * num_rays);
    // {pin_comment} the observer window intersection lock.
    {malloc_pinned}(all_photons_host.window_event_found, sizeof(bool) * num_rays);
    // {pin_comment} the source emission plane intersection lock.
    {malloc_pinned}(all_photons_host.source_event_found, sizeof(bool) * num_rays);

    {stream_setup_str}

    // --- DOUBLE-BUFFERED BRIDGE ARRAYS ---
    // Extraction buffer used by the TimeSlotManager to map sparse indices to contiguous execution blocks.
    long int *chunk_buffer[2];
    // Bridge array staging the state vector $f^\mu$ for memory transfers.
    double *f_bridge[2];
    // Bridge array staging the first derivative $\dot{{f}}^\mu$ for memory transfers.
    double *f_p_bridge[2];
    // Bridge array staging the second derivative $\ddot{{f}}^\mu$ for memory transfers.
    double *f_p_p_bridge[2];
    // Bridge array staging the affine parameter $\lambda$ for memory transfers.
    double *affine_bridge[2];
    // Bridge array staging the current integration step size $h$ for memory transfers.
    double *h_bridge[2];
    // Bridge array staging the current trajectory termination status for memory transfers.
    termination_type_t *status_bridge[2];
    // Bridge array staging the number of step-size rejections for memory transfers.
    int *retries_bridge[2];
    // Bridge array staging the previous observer window boundary side flag for memory transfers.
    bool *on_pos_window_prev_bridge[2];
    // Bridge array staging the previous source emission boundary side flag for memory transfers.
    bool *on_pos_source_prev_bridge[2];
    // Bridge array staging the historical affine parameter $\lambda_{{n-1}}$ for chunked memory transfers.
    double *affine_p_bridge[2];
    // Bridge array staging the historical affine parameter $\lambda_{{n-2}}$ for chunked memory transfers.
    double *affine_p_p_bridge[2];
    // Bridge array staging the observer window event lock for memory transfers.
    bool *window_event_found_bridge[2];
    // Bridge array staging the source emission event lock for memory transfers.
    bool *source_event_found_bridge[2];

    // --- DOUBLE-BUFFERED VRAM SCRATCHPADS ---
    // Scratchpad tracking the current state vector $f^\mu$ bounding the RKF45 step.
    double *d_f_bundle[2];
    // Scratchpad locking the anchor state vector $f_{{start}}$ to calculate the final stage update.
    double *d_f_start_bundle[2];
    // Scratchpad tracking the intermediate cumulative RKF45 stage updates.
    double *d_f_temp_bundle[2];
    // Scratchpad tracking the history state $f^\mu_{{n-1}}$ for geometric intersection detection.
    double *d_f_prev_bundle[2];
    // Scratchpad tracking the history state $f^\mu_{{n-2}}$ for geometric intersection detection.
    double *d_f_pre_prev_bundle[2];
    // Scratchpad persisting the symmetric metric tensor $g_{{\mu\nu}}$.
    double *d_metric_bundle[2];
    // Scratchpad persisting the Christoffel symbols $\Gamma^\alpha_{{\beta\gamma}}$.
    double *d_connection_bundle[2];
    // Derivative tensor storing $\dot{{f}}^\mu$ across all 6 intermediate RKF45 stages.
    double *d_k_bundle[2];
    // Array regulating active integration step sizing $h$.
    double *d_h[2];
    // Array regulating total affine parameter progress $\lambda$.
    double *d_affine[2];
    // Array holding the current trajectory status limits.
    termination_type_t *d_status[2];
    // Array tracking sequential error rejections per photon.
    int *d_retries[2];
    // Array flagging the previous observer window boundary side.
    bool *d_on_pos_window_prev[2];
    // Array flagging the previous source emission boundary side.
    bool *d_on_pos_source_prev[2];
    // Array tracking historical affine parameter $\lambda_{{n-1}}$.
    double *d_affine_prev[2];
    // Array tracking historical affine parameter $\lambda_{{n-2}}$.
    double *d_affine_pre_prev[2];
    // Array guarding the window intersection coordinates from multi-trigger overwrites.
    bool *d_window_event_found[2];
    // Array guarding the source intersection coordinates from multi-trigger overwrites.
    bool *d_source_event_found[2];
    // Array carrying the absolute master indices $m_{{idx}}$ mapping the execution chunk.
    long int *d_chunk_buffer[2];

    // Loop iterator for instantiating the double-buffered operational arrays.
    for (int s = 0; s < 2; ++s) {{
        // {bridge_alloc_comment}
        {malloc_pinned}(chunk_buffer[s], sizeof(long int) * BUNDLE_CAPACITY); // Pin chunk buffers.
        {malloc_pinned}(f_bridge[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Pin $f^\mu$ bridges.
        {malloc_pinned}(f_p_bridge[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Pin $\dot{{f}}^\mu$ bridges.
        {malloc_pinned}(f_p_p_bridge[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Pin $\ddot{{f}}^\mu$ bridges.
        {malloc_pinned}(affine_bridge[s], sizeof(double) * BUNDLE_CAPACITY); // Pin $\lambda$ bridges.
        {malloc_pinned}(h_bridge[s], sizeof(double) * BUNDLE_CAPACITY); // Pin $h$ bridges.
        {malloc_pinned}(status_bridge[s], sizeof(termination_type_t) * BUNDLE_CAPACITY); // Pin status bridges.
        {malloc_pinned}(retries_bridge[s], sizeof(int) * BUNDLE_CAPACITY); // Pin retries bridges.
        {malloc_pinned}(on_pos_window_prev_bridge[s], sizeof(bool) * BUNDLE_CAPACITY); // Pin window flag bridges.
        {malloc_pinned}(on_pos_source_prev_bridge[s], sizeof(bool) * BUNDLE_CAPACITY); // Pin source flag bridges.
        {malloc_pinned}(affine_p_bridge[s], sizeof(double) * BUNDLE_CAPACITY); // Pin $\lambda_{{n-1}}$ bridges.
        {malloc_pinned}(affine_p_p_bridge[s], sizeof(double) * BUNDLE_CAPACITY); // Pin $\lambda_{{n-2}}$ bridges.
        {malloc_pinned}(window_event_found_bridge[s], sizeof(bool) * BUNDLE_CAPACITY); // Pin window lock bridges.
        {malloc_pinned}(source_event_found_bridge[s], sizeof(bool) * BUNDLE_CAPACITY); // Pin source lock bridges.

        // {scratch_alloc_comment}
        {malloc_device}(d_f_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate $f^\mu$ scratchpad.
        {malloc_device}(d_f_start_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate $f_{{start}}$ scratchpad.
        {malloc_device}(d_f_temp_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate temporary stage scratchpad.
        {malloc_device}(d_f_prev_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate $f^\mu_{{n-1}}$ scratchpad.
        {malloc_device}(d_f_pre_prev_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate $f^\mu_{{n-2}}$ scratchpad.
        {malloc_device}(d_metric_bundle[s], sizeof(double) * 10 * BUNDLE_CAPACITY); // Allocate $g_{{\mu\nu}}$ scratchpad.
        {malloc_device}(d_connection_bundle[s], sizeof(double) * 40 * BUNDLE_CAPACITY); // Allocate $\Gamma^\alpha_{{\beta\gamma}}$ scratchpad.
        {malloc_device}(d_k_bundle[s], sizeof(double) * 6 * 9 * BUNDLE_CAPACITY); // Allocate derivative scratchpad.
        {malloc_device}(d_h[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate $h$ scratchpad.
        {malloc_device}(d_affine[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate $\lambda$ scratchpad.
        {malloc_device}(d_status[s], sizeof(termination_type_t) * BUNDLE_CAPACITY); // Allocate status scratchpad.
        {malloc_device}(d_retries[s], sizeof(int) * BUNDLE_CAPACITY); // Allocate retries scratchpad.
        {malloc_device}(d_on_pos_window_prev[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate window flag scratchpad.
        {malloc_device}(d_on_pos_source_prev[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate source flag scratchpad.
        {malloc_device}(d_affine_prev[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate $\lambda_{{n-1}}$ scratchpad.
        {malloc_device}(d_affine_pre_prev[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate $\lambda_{{n-2}}$ scratchpad.
        {malloc_device}(d_window_event_found[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate window lock scratchpad.
        {malloc_device}(d_source_event_found[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate source lock scratchpad.
        {malloc_device}(d_chunk_buffer[s], sizeof(long int) * BUNDLE_CAPACITY); // Allocate chunk mapping scratchpad.
    }}

    // Device-native struct pointer storing the final physical plane intersections.
    blueprint_data_t *d_results_buffer;
    // {dev_comment} the blueprint results buffer to avoid mid-computation memory transfers.
    {malloc_device}(d_results_buffer, sizeof(blueprint_data_t) * num_rays);

    // Host-bound struct managing temporal binning of photon trajectories $x^\mu$.
    TimeSlotManager tsm;
    // Initializes the temporal manager strictly on the CPU to coordinate the Split-Pipeline batches.
    slot_manager_init(&tsm, commondata->slot_manager_t_min, commondata->t_start + 1.0, commondata->slot_manager_delta_t, num_rays);

    // --- DIAGNOSTIC MEMORY ALLOCATION ---
    // Host pointer tracking the initial conserved quantities prior to integration.
    conserved_quantities_t *initial_cq_host = NULL;
    // Host pointer tracking the terminal conserved quantities post integration.
    conserved_quantities_t *final_cq_host = NULL;

    if (commondata->perform_conservation_check) {{
        // {pin_comment} the initial diagnostic data.
        {malloc_pinned}(initial_cq_host, sizeof(conserved_quantities_t) * num_rays);
        // {pin_comment} the final diagnostic data.
        {malloc_pinned}(final_cq_host, sizeof(conserved_quantities_t) * num_rays);
    }}


    /* Algorithmic Step: Evaluate initial coordinate states and map global spacetime metrics to memory bounds. Hardware Justification: This pre-computation stage occurs prior to the temporal loop to maximize coalesced memory access during iterative integration. */
    // --- 2. INITIALIZATION PHASE ---
    
    double window_center_out[3]; // 3D array storing the spatial Cartesian coordinates $x^i$ of the observer window center.
    double n_x_out[3]; // 3D orthonormal basis vector pointing along the $x$-axis of the local window geometry.
    double n_y_out[3]; // 3D orthonormal basis vector pointing along the $y$-axis of the local window geometry.
    double n_z_out[3]; // 3D orthonormal basis vector pointing along the $z$-axis of the local window geometry.

    // Hardware Justification: Operates synchronously as the primary state array must be fully populated before pipeline dispatch.
    {set_intitial_con}

    /* Algorithmic Step: Scans the master Host SoA immediately following the initialization kernel call. Hardware Justification: This architectural step verifies that the Host-side logic has correctly populated the starting coordinates $x^\mu$ and zeroed the temporal momentum $p_t$ and affine parameter $\lambda$ components. */
    // --- DIAGNOSTIC PROBE: DETAILED INITIAL POSITION & PLACEHOLDER ALIGNMENT ---
    
    long int init_mismatch_count = 0; // Accumulator tracking the total number of physical state initializations that failed structural validation.
    long int mismatch_t = 0; // Counter tracking validation failures for the temporal coordinate $t$.
    long int mismatch_x = 0; // Counter tracking validation failures for the spatial coordinate $x$.
    long int mismatch_y = 0; // Counter tracking validation failures for the spatial coordinate $y$.
    long int mismatch_z = 0; // Counter tracking validation failures for the spatial coordinate $z$.
    long int mismatch_pt = 0; // Counter tracking validation failures for the temporal momentum $p_t$.
    long int mismatch_lam = 0; // Counter tracking validation failures for the affine parameter $\lambda$.

    for (long int p = 0; p < num_rays; p++) {{ // Loop iterator index $p$ mapping to a unique photon trajectory $x^\mu$ during diagnostic validation.
        const double t_check = all_photons_host.f[0 * num_rays + p]; // Evaluates the current temporal coordinate $t$ from the Host SoA.
        const double x_check = all_photons_host.f[1 * num_rays + p]; // Evaluates the current spatial coordinate $x$ from the Host SoA.
        const double y_check = all_photons_host.f[2 * num_rays + p]; // Evaluates the current spatial coordinate $y$ from the Host SoA.
        const double z_check = all_photons_host.f[3 * num_rays + p]; // Evaluates the current spatial coordinate $z$ from the Host SoA.
        const double pt_check = all_photons_host.f[4 * num_rays + p]; // Evaluates the initial temporal momentum $p_t$ from the Host SoA.
        const double lam_check = all_photons_host.f[8 * num_rays + p]; // Evaluates the initial affine parameter $\lambda$ from the Host SoA.

        bool fail_t = fabs(t_check - commondata->t_start) > 1e-10; // Boolean flag indicating temporal coordinate $t$ validation failure.
        bool fail_x = fabs(x_check - commondata->camera_pos_x) > 1e-10; // Boolean flag indicating spatial coordinate $x$ validation failure.
        bool fail_y = fabs(y_check - commondata->camera_pos_y) > 1e-10; // Boolean flag indicating spatial coordinate $y$ validation failure.
        bool fail_z = fabs(z_check - commondata->camera_pos_z) > 1e-10; // Boolean flag indicating spatial coordinate $z$ validation failure.
        bool fail_pt = fabs(pt_check) > 1e-15; // Boolean flag indicating temporal momentum $p_t$ validation failure.
        bool fail_lam = fabs(lam_check) > 1e-15; // Boolean flag indicating affine parameter $\lambda$ validation failure.

        if (fail_t) mismatch_t++; // Increments the validation failure counter for temporal coordinate $t$.
        if (fail_x) mismatch_x++; // Increments the validation failure counter for spatial coordinate $x$.
        if (fail_y) mismatch_y++; // Increments the validation failure counter for spatial coordinate $y$.
        if (fail_z) mismatch_z++; // Increments the validation failure counter for spatial coordinate $z$.
        if (fail_pt) mismatch_pt++; // Increments the validation failure counter for temporal momentum $p_t$.
        if (fail_lam) mismatch_lam++; // Increments the validation failure counter for affine parameter $\lambda$.

        if (fail_t || fail_x || fail_y || fail_z || fail_pt || fail_lam) {{
            init_mismatch_count++; // Increments the total accumulation of trajectory structural validation failures.
        }}
    }}

    if (init_mismatch_count > 0) {{
        const double mismatch_percent = ((double)init_mismatch_count / (double)num_rays) * 100.0; // Calculates the failure rate percentage for the initial state alignment.
        // Hardware Justification: This is a soft warning to surface initialization inconsistencies without halting execution.
        printf("[DIAGNOSTIC] Initialization Alignment Check: %ld out of %ld rays (%.2f%%) fail coordinate/placeholder validation.\n", init_mismatch_count, num_rays, mismatch_percent);
    }}

    long int num_batches = (num_rays + BUNDLE_CAPACITY - 1) / BUNDLE_CAPACITY; // Total integer calculation defining total iterative blocks required to process all photon indices.

    for (long int init_batch = 0; init_batch < num_batches; ++init_batch) {{ // Loop iterator $init_batch$ for evaluating the initialization constraint across sequential blocks.
        long int start_idx = init_batch * BUNDLE_CAPACITY; // Absolute starting index mapped to the master SoA for the current initialization batch.
        long int chunk_size = MIN((long int)BUNDLE_CAPACITY, num_rays - start_idx); // Dynamically sized operational boundary ensuring the active chunk does not exceed total trajectories.

        for (int init_i = 0; init_i < chunk_size; ++init_i) {{ // Loop index $init_i$ iterating over the specific initialization batch elements to pack the bridge.
            long int master_idx = start_idx + init_i; // Computes the absolute master index $m_{{idx}}$ tracking the photon within the global array.
            for (int init_k = 0; init_k < 9; ++init_k) {{ // Loop index $init_k$ iterating over the 9 tensor components of the state vector $f^\mu$.
                f_bridge[0][init_k * BUNDLE_CAPACITY + init_i] = all_photons_host.f[init_k * num_rays + master_idx]; // Assigns the active tensor state component to the primary bridge.
            }}
        }}

        for (int c_k = 0; c_k < 9; ++c_k) {{ // Loop index $c_k$ orchestrating the memory transfer of the 9 state vector $f^\mu$ components.
            {h2d_f_comment}
            {memcpy_async("d_f_bundle[0] + c_k * BUNDLE_CAPACITY", "f_bridge[0] + c_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size", "cudaMemcpyHostToDevice", "streams[0]")}
        }}

        // Hardware Justification: Calculates symmetric metric tensor $g_{{\mu\nu}}$ strictly on the primary operational pipeline for the Hamiltonian constraint.
        interpolation_kernel_{spacetime_name}(commondata,d_f_bundle[0], d_metric_bundle[0], NULL, chunk_size{stream_arg});

        /* Algorithmic Step: Extracts metric payload to confirm numerical stability prior to momentum solving. Hardware Justification: Ensures constraint solver convergence by intercepting unphysical spacetime regions immediately after calculation. */
        // --- DIAGNOSTIC PROBE: METRIC INTEGRITY CHECK ---
        
        double *metric_diag_bridge; // Pointer storing temporary metric data to validate the interpolation sequence.
        {malloc_pinned}(metric_diag_bridge, sizeof(double) * 10 * BUNDLE_CAPACITY); // Memory allocation: Temporary bridge mapped to maximize throughput for the metric $g_{{\mu\nu}}$ diagnostic sequence.
        
        for (int m_k = 0; m_k < 10; ++m_k) {{ 
            // Loop index $m_k$ orchestrating memory transfer of the 10 metric tensor $g_{{\mu\nu}}$ components.
            {d2h_metric_comment}
            {memcpy_async("metric_diag_bridge + m_k * BUNDLE_CAPACITY", "d_metric_bundle[0] + m_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size", "cudaMemcpyDeviceToHost", "streams[0]")}
        }}
        {sync_metric_comment}
        {stream_sync('streams[0]')}

        long int metric_nan_count = 0; // Accumulator tracking the total number of metric tensor evaluations containing non-finite values.
        for (int m_diag_i = 0; m_diag_i < chunk_size; ++m_diag_i) {{ // Loop iterator $m_diag_i$ scanning each trajectory within the current initialization chunk.
            bool m_has_nan = false; // Boolean flag indicating if the specific metric tensor $g_{{\mu\nu}}$ contains a NaN or Inf value.
            for (int m_diag_k = 0; m_diag_k < 10; ++m_diag_k) {{ // Loop index $m_diag_k$ iterating over the 10 independent components of the symmetric metric tensor $g_{{\mu\nu}}$.
                if (isnan(metric_diag_bridge[m_diag_k * BUNDLE_CAPACITY + m_diag_i]) || 
                    isinf(metric_diag_bridge[m_diag_k * BUNDLE_CAPACITY + m_diag_i])) {{
                    m_has_nan = true; // Flags the trajectory metric state as invalid due to a non-finite value.
                    break; // Terminates the tensor component loop early to conserve execution cycles upon detecting a failure.
                }}
            }}
            if (m_has_nan) metric_nan_count++; // Increments the total accumulation of corrupted metric tensor evaluations.
        }}

        if (metric_nan_count > 0) {{
            // Hardware Justification: This is a soft warning to highlight numerical metric singularities without aborting the physics engine.
            printf("[DIAGNOSTIC] Init Batch %ld: %ld rays have invalid Metric G_mu_nu before p_t solve.\n", init_batch, metric_nan_count);
        }}
        {free_pinned}(metric_diag_bridge); // Memory Free: Purges the diagnostic bridge utilized for metric integrity checks.
        
        // Hardware Justification: Solves the constraint $p_\mu p^\mu = 0$ to find the temporal momentum $p_t$ natively on the active pipeline.
        p0_reverse_kernel(d_f_bundle[0], d_metric_bundle[0], chunk_size{stream_arg});

        for (int c_k = 0; c_k < 9; ++c_k) {{ // Loop index $c_k$ orchestrating memory transfer of the 9 constrained state vector $f^\mu$ components.
            {d2h_f_comment}
            {memcpy_async("f_bridge[0] + c_k * BUNDLE_CAPACITY", "d_f_bundle[0] + c_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size", "cudaMemcpyDeviceToHost", "streams[0]")}
        }}
        {sync_f_comment}
        {stream_sync('streams[0]')}

        long int nan_count = 0; // Accumulator tracking the total number of physical states $f^\mu$ containing NaN values post-constraint solving.
        for (int gather_i = 0; gather_i < chunk_size; ++gather_i) {{ // Loop iterator $gather_i$ scanning each trajectory within the retrieved initialization chunk.
            long int master_idx = start_idx + gather_i; // Computes the absolute master index $m_{{idx}}$ mapping the localized chunk to the global master SoA.
            bool has_nan = false; // Boolean flag indicating if the specific state vector $f^\mu$ contains a NaN value.
            for (int gather_k = 0; gather_k < 9; ++gather_k) {{ // Loop index $gather_k$ iterating over the 9 tensor components of the state vector $f^\mu$.
                double val = f_bridge[0][gather_k * BUNDLE_CAPACITY + gather_i]; // Evaluates the updated numerical value of the specific tensor component.
                all_photons_host.f[gather_k * num_rays + master_idx] = val; // Maps the valid constrained tensor scalar back to the global Host SoA.
                if (isnan(val)) has_nan = true; // Flags the physical state vector as invalid due to a non-finite evaluation.
            }}
            if (has_nan) nan_count++; // Increments the total count of unresolved physical state vectors $f^\mu$.
        }}

        if (nan_count > 0) {{
            // Hardware Justification: This is a soft warning alerting to unresolved constraints $p_\mu p^\mu = 0$ for isolated trajectories.
            printf("[DIAGNOSTIC] Init Batch %ld: %ld rays contain NaN in state f^mu after p_t solve.\n", init_batch, nan_count);
        }}
    }}

    /* Algorithmic Step: Evaluate initial conserved quantities immediately after generating valid physical null states. Hardware Justification: This baselines the data via chunked pipeline kernels before the Split-Pipeline begins mutating the state vectors. */
    // --- BASELINE CONSERVED QUANTITIES ---
    
    if (commondata->perform_conservation_check) {{
        // Hardware Justification: Executes global conserved quantities natively via chunked device parameters before pipeline processing.
        calculate_conserved_quantities_universal_{spacetime_name}_photon(commondata, &all_photons_host, num_rays, initial_cq_host);
    }}

    long int sync_i; // Loop iterator index $sync_i$ spanning the entire global ray count to synchronize starting properties across history states.
    for(sync_i = 0; sync_i < num_rays; ++sync_i) {{
        int sync_k; // Loop index $sync_k$ iterating over the 9 tensor components to populate the historical derivatives $\dot{{f}}^\mu$ and $\ddot{{f}}^\mu$.
        for (sync_k = 0; sync_k < 9; ++sync_k) {{
            all_photons_host.f_p[sync_k * num_rays + sync_i] = all_photons_host.f[sync_k * num_rays + sync_i]; // Propagates the initial coordinate state vector $f^\mu$ to the first history derivative matrix.
            all_photons_host.f_p_p[sync_k * num_rays + sync_i] = all_photons_host.f[sync_k * num_rays + sync_i]; // Propagates the initial coordinate state vector $f^\mu$ to the second history derivative matrix.
        }}
        all_photons_host.status[sync_i] = ACTIVE; // Assigns the initial trajectory activity enum for the global physics engine.
        all_photons_host.affine_param[sync_i] = 0.0; // Sets the initial baseline progression scalar for the affine parameter $\lambda$.
        all_photons_host.rejection_retries[sync_i] = 0; // Clears the error rejection scalar to initialize the step size convergence tracking.

        all_photons_host.affine_param_p[sync_i] = 0.0; // Initializes the first historical affine parameter $\lambda_{{n-1}}$.
        all_photons_host.affine_param_p_p[sync_i] = 0.0; // Initializes the second historical affine parameter $\lambda_{{n-2}}$.
        all_photons_host.window_event_found[sync_i] = false; // Sets the observer window intersection logical lock to false.
        all_photons_host.source_event_found[sync_i] = false; // Sets the source emission intersection logical lock to false.

        int s_idx = slot_get_index(&tsm, all_photons_host.f[sync_i]); // Integer index $s_{{idx}}$ mapping the current photon's temporal coordinate $t$ to a discrete execution bin in the TimeSlotManager.
        if (s_idx != -1) {{
            slot_add_photon(&tsm, s_idx, sync_i); // Registers the active photon index to its corresponding temporal bin mapped by the orchestrator.
        }}
    }}

    // --- 3. TEMPORAL LOOP (The Engine) ---
    // Integer tracking the global number of active photon trajectories to allow early loop termination.
    long int total_active_photons = num_rays; 

    // Outer loop iterator for the physical time bins.
    for (int slot_idx = tsm.num_slots - 1; slot_idx >= 0; --slot_idx) {{  
        // Evaluates the early exit condition to terminate the temporal engine if all geometric trajectories have concluded.
        if (total_active_photons <= 0) {{ 
            break; // Terminates the temporal engine early to conserve hardware cycles.
        }}

        int current = 0; // Integer index tracking the primary active hardware stream for execution overlapping.
        int next = 1; // Integer index tracking the secondary hardware stream preparing the upcoming payload.
        long int active_chunks[2] = {{ 0, 0}}; // 1D array storing the total number of trajectories queued for each operational stream.

        /* Algorithmic Step: Populate the first bridge and launch asynchronous integration on the primary stream. Hardware Justification: Initializes the pipeline to allow subsequent Host-side packing to overlap with GPU compute. */
        // --- PHASE A: PRIME THE PUMP (Stream 0) ---
        
        active_chunks[current] = MIN((long int)BUNDLE_CAPACITY, tsm.slot_counts[slot_idx]); // Evaluates the active chunk size bounding the PCIe transfer to avoid VRAM overflow.
        
        if (active_chunks[current] > 0) {{ 
            slot_remove_chunk(&tsm, slot_idx, chunk_buffer[current], active_chunks[current]); // Extracts the execution chunk mapping from the Host-side temporal bin.

            for (int bridge_i = 0; bridge_i < active_chunks[current]; ++bridge_i) {{  // Loop iterator $bridge_i$ packing the physical state payloads into the Host-side bridge arrays.
                long int m_idx = chunk_buffer[current][bridge_i]; // Absolute master index $m_{{ idx}}$ mapping the active payload to the global trajectory matrix.
                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ iterating over the 9 tensor components of the state vectors.
                    f_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f[c_k * num_rays + m_idx]; // Packs the coordinate state vector $f^\mu$ into the transfer bridge.
                    f_p_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p[c_k * num_rays + m_idx]; // Packs the first derivative $\dot{{ f}}^\mu$ into the transfer bridge.
                    f_p_p_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p_p[c_k * num_rays + m_idx]; // Packs the second derivative $\ddot{{ f}}^\mu$ into the transfer bridge.
                }}
                h_bridge[current][bridge_i] = all_photons_host.h[m_idx]; // Packs the current integration step size $h$ into the transfer bridge.
                status_bridge[current][bridge_i] = all_photons_host.status[m_idx]; // Packs the trajectory status enum into the transfer bridge.
                retries_bridge[current][bridge_i] = all_photons_host.rejection_retries[m_idx]; // Packs the error rejection scalar into the transfer bridge.
                affine_bridge[current][bridge_i] = all_photons_host.affine_param[m_idx]; // Packs the affine parameter $\lambda$ into the transfer bridge.
                on_pos_window_prev_bridge[current][bridge_i] = all_photons_host.on_positive_side_of_window_prev[m_idx]; // Packs the observer window boundary flag into the transfer bridge.
                on_pos_source_prev_bridge[current][bridge_i] = all_photons_host.on_positive_side_of_source_prev[m_idx]; // Packs the source emission boundary flag into the transfer bridge.
                affine_p_bridge[current][bridge_i] = all_photons_host.affine_param_p[m_idx]; // Packs the historical affine parameter $\lambda_{{ n-1}}$ into the transfer bridge.
                affine_p_p_bridge[current][bridge_i] = all_photons_host.affine_param_p_p[m_idx]; // Packs the historical affine parameter $\lambda_{{ n-2}}$ into the transfer bridge.
                window_event_found_bridge[current][bridge_i] = all_photons_host.window_event_found[m_idx]; // Packs the observer window intersection lock into the transfer bridge.
                source_event_found_bridge[current][bridge_i] = all_photons_host.source_event_found[m_idx]; // Packs the source emission intersection lock into the transfer bridge.
            }}

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating Host-to-Device transfer of the 9 state vector components.
                // Host-to-Device transfer: Asynchronously pushes bounded state vectors $f^\mu$ to VRAM strictly on stream [current] to minimize latency.
                {memcpy_async("d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "f_bridge[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
                // Host-to-Device transfer: Asynchronously pushes first derivatives $\dot{{ f}}^\mu$ to VRAM strictly on stream [current] to minimize latency.
                {memcpy_async("d_f_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "f_p_bridge[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
                // Host-to-Device transfer: Asynchronously pushes second derivatives $\ddot{{ f}}^\mu$ to VRAM strictly on stream [current] to minimize latency.
                {memcpy_async("d_f_pre_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "f_p_p_bridge[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            }}
            // Host-to-Device transfer: Asynchronously pushes step sizes $h$ to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_h[current]", "h_bridge[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes status enums to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_status[current]", "status_bridge[current]", "sizeof(termination_type_t) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes rejection scalars to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_retries[current]", "retries_bridge[current]", "sizeof(int) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes affine parameters $\lambda$ to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_affine[current]", "affine_bridge[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes window boundary flags to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_on_pos_window_prev[current]", "on_pos_window_prev_bridge[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes source boundary flags to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_on_pos_source_prev[current]", "on_pos_source_prev_bridge[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes historical affine parameters $\lambda_{{ n-1}}$ to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_affine_prev[current]", "affine_p_bridge[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes historical affine parameters $\lambda_{{ n-2}}$ to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_affine_pre_prev[current]", "affine_p_p_bridge[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes window intersection locks to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_window_event_found[current]", "window_event_found_bridge[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes source intersection locks to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_source_event_found[current]", "source_event_found_bridge[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes chunk indices $m_{{ idx}}$ to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_chunk_buffer[current]", "chunk_buffer[current]", "sizeof(long int) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating Device-to-Device baseline setup of the 9 state vector components.
                // Device-to-Device transfer: Duplicates the initial physical state vector $f^\mu$ to anchor the final RKF45 evaluation.
                {memcpy_async("d_f_start_bundle[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToDevice", "streams[current]")}
                // Device-to-Device transfer: Primes the temporary state vector bundle $f^\mu_{{ temp}}$ for iterative stage accumulation.
                {memcpy_async("d_f_temp_bundle[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToDevice", "streams[current]")}
            }}

            for (int stage = 1; stage <= 6; ++stage) {{  // Loop iterator $stage$ executing the 6 discrete stages of the RKF45 Runge-Kutta numerical solver.
                // Kernel Launch: Evaluates the metric tensor $g_{{ \mu\nu}}$ and connection $\Gamma^\alpha_{{ \beta\gamma}}$ asynchronously on the active stream.
                interpolation_kernel_{ spacetime_name}(commondata, d_f_temp_bundle[current], d_metric_bundle[current], d_connection_bundle[current], active_chunks[current]{stream_arg_current}); 
                // Kernel Launch: Computes the geodesic equation right-hand-side derivatives $\dot{{ f}}^\mu$ asynchronously on the active stream.
                calculate_ode_rhs_kernel(d_f_temp_bundle[current], d_metric_bundle[current], d_connection_bundle[current], d_k_bundle[current], stage, active_chunks[current]{stream_arg_current});
                // Kernel Launch: Accumulates the intermediate RKF45 stage numerical updates asynchronously on the active stream.
                rkf45_stage_update(d_f_start_bundle[current], d_k_bundle[current], d_h[current], stage, active_chunks[current], d_f_temp_bundle[current]{stream_arg_current});
            }}
            
            // Kernel Launch: Applies Cash-Karp error control to finalize the step-size $h$ and update the integration baseline.
            rkf45_finalize_and_control(commondata, d_f_bundle[current], d_f_start_bundle[current], d_k_bundle[current], d_h[current], d_status[current], d_affine[current], d_retries[current], active_chunks[current]{stream_arg_current});
            // Kernel Launch: Detects geometric events and records intersection coordinate states asynchronously on the active stream.
            event_detection_manager_kernel(commondata, d_f_bundle[current], d_f_prev_bundle[current], d_f_pre_prev_bundle[current], d_affine[current], d_affine_prev[current], d_affine_pre_prev[current], d_results_buffer, d_status[current], d_on_pos_window_prev[current], d_on_pos_source_prev[current], d_window_event_found[current], d_source_event_found[current], d_chunk_buffer[current], active_chunks[current]{stream_arg_current});

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating Device-to-Host transfer of the 9 state vector components.
                // Device-to-Host transfer: Retrieves updated coordinate states $f^\mu$ back to CPU RAM asynchronously on the active stream.
                {memcpy_async("f_bridge[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
                // Device-to-Host transfer: Retrieves updated first derivatives $\dot{{ f}}^\mu$ back to CPU RAM asynchronously on the active stream.
                {memcpy_async("f_p_bridge[current] + c_k * BUNDLE_CAPACITY", "d_f_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
                // Device-to-Host transfer: Retrieves updated second derivatives $\ddot{{ f}}^\mu$ back to CPU RAM asynchronously on the active stream.
                {memcpy_async("f_p_p_bridge[current] + c_k * BUNDLE_CAPACITY", "d_f_pre_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            }}
            // Device-to-Host transfer: Retrieves active step sizes $h$ back to CPU RAM asynchronously on the active stream.
            {memcpy_async("h_bridge[current]", "d_h[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves updated status enums back to CPU RAM asynchronously on the active stream.
            {memcpy_async("status_bridge[current]", "d_status[current]", "sizeof(termination_type_t) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves active rejection counts back to CPU RAM asynchronously on the active stream.
            {memcpy_async("retries_bridge[current]", "d_retries[current]", "sizeof(int) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves total affine progression $\lambda$ back to CPU RAM asynchronously on the active stream.
            {memcpy_async("affine_bridge[current]", "d_affine[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves updated window boundary flags back to CPU RAM asynchronously on the active stream.
            {memcpy_async("on_pos_window_prev_bridge[current]", "d_on_pos_window_prev[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves updated source boundary flags back to CPU RAM asynchronously on the active stream.
            {memcpy_async("on_pos_source_prev_bridge[current]", "d_on_pos_source_prev[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves historical affine parameter $\lambda_{{ n-1}}$ back to CPU RAM asynchronously on the active stream.
            {memcpy_async("affine_p_bridge[current]", "d_affine_prev[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves historical affine parameter $\lambda_{{ n-2}}$ back to CPU RAM asynchronously on the active stream.
            {memcpy_async("affine_p_p_bridge[current]", "d_affine_pre_prev[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves active window locks back to CPU RAM asynchronously on the active stream.
            {memcpy_async("window_event_found_bridge[current]", "d_window_event_found[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves active source locks back to CPU RAM asynchronously on the active stream.
            {memcpy_async("source_event_found_bridge[current]", "d_source_event_found[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
        }}

        /* Algorithmic Step: Continuously alternate between streams, packing the next payload while syncing the current. Hardware Justification: Completely hides PCIe DMA transfer latency behind the active RKF45 integration compute time. */
        // --- PHASE B: THE OVERLAP LOOP ---
        
        while (active_chunks[current] > 0 || tsm.slot_counts[slot_idx] > 0) {{  
            
            active_chunks[next] = MIN((long int)BUNDLE_CAPACITY, tsm.slot_counts[slot_idx]); // Evaluates the active chunk size bounding the upcoming PCIe transfer execution block.
            if (active_chunks[next] > 0) {{ 
                slot_remove_chunk(&tsm, slot_idx, chunk_buffer[next], active_chunks[next]); // Extracts the next execution chunk mapping from the Host-side temporal bin.
                
                for (int bridge_i = 0; bridge_i < active_chunks[next]; ++bridge_i) {{  // Loop iterator $bridge_i$ packing the physical state payloads into the next Host-side bridge array.
                    long int m_idx = chunk_buffer[next][bridge_i]; // Absolute master index $m_{{ idx}}$ mapping the active payload to the global trajectory matrix.
                    for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ iterating over the 9 tensor components of the state vectors.
                        f_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f[c_k * num_rays + m_idx]; // Packs the coordinate state vector $f^\mu$ into the transfer bridge.
                        f_p_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p[c_k * num_rays + m_idx]; // Packs the first derivative $\dot{{ f}}^\mu$ into the transfer bridge.
                        f_p_p_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p_p[c_k * num_rays + m_idx]; // Packs the second derivative $\ddot{{ f}}^\mu$ into the transfer bridge.
                    }}
                    h_bridge[next][bridge_i] = all_photons_host.h[m_idx]; // Packs the current integration step size $h$ into the transfer bridge.
                    status_bridge[next][bridge_i] = all_photons_host.status[m_idx]; // Packs the trajectory status enum into the transfer bridge.
                    retries_bridge[next][bridge_i] = all_photons_host.rejection_retries[m_idx]; // Packs the error rejection scalar into the transfer bridge.
                    affine_bridge[next][bridge_i] = all_photons_host.affine_param[m_idx]; // Packs the affine parameter $\lambda$ into the transfer bridge.
                    on_pos_window_prev_bridge[next][bridge_i] = all_photons_host.on_positive_side_of_window_prev[m_idx]; // Packs the observer window boundary flag into the transfer bridge.
                    on_pos_source_prev_bridge[next][bridge_i] = all_photons_host.on_positive_side_of_source_prev[m_idx]; // Packs the source emission boundary flag into the transfer bridge.
                    affine_p_bridge[next][bridge_i] = all_photons_host.affine_param_p[m_idx]; // Packs the historical affine parameter $\lambda_{{ n-1}}$ into the transfer bridge.
                    affine_p_p_bridge[next][bridge_i] = all_photons_host.affine_param_p_p[m_idx]; // Packs the historical affine parameter $\lambda_{{ n-2}}$ into the transfer bridge.
                    window_event_found_bridge[next][bridge_i] = all_photons_host.window_event_found[m_idx]; // Packs the observer window intersection lock into the transfer bridge.
                    source_event_found_bridge[next][bridge_i] = all_photons_host.source_event_found[m_idx]; // Packs the source emission intersection lock into the transfer bridge.
                }}

                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating Host-to-Device transfer of the 9 state vector components for the upcoming payload.
                    // Host-to-Device transfer: Asynchronously pushes bounded state vectors $f^\mu$ to VRAM strictly on stream [next] to overlap execution.
                    {memcpy_async("d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "f_bridge[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                    // Host-to-Device transfer: Asynchronously pushes first derivatives $\dot{{ f}}^\mu$ to VRAM strictly on stream [next] to overlap execution.
                    {memcpy_async("d_f_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "f_p_bridge[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                    // Host-to-Device transfer: Asynchronously pushes second derivatives $\ddot{{ f}}^\mu$ to VRAM strictly on stream [next] to overlap execution.
                    {memcpy_async("d_f_pre_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "f_p_p_bridge[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                }}
                // Host-to-Device transfer: Asynchronously pushes step sizes $h$ to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_h[next]", "h_bridge[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes status enums to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_status[next]", "status_bridge[next]", "sizeof(termination_type_t) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes rejection scalars to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_retries[next]", "retries_bridge[next]", "sizeof(int) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes affine parameters $\lambda$ to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_affine[next]", "affine_bridge[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes window boundary flags to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_on_pos_window_prev[next]", "on_pos_window_prev_bridge[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes source boundary flags to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_on_pos_source_prev[next]", "on_pos_source_prev_bridge[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes historical affine parameters $\lambda_{{ n-1}}$ to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_affine_prev[next]", "affine_p_bridge[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes historical affine parameters $\lambda_{{ n-2}}$ to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_affine_pre_prev[next]", "affine_p_p_bridge[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes window intersection locks to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_window_event_found[next]", "window_event_found_bridge[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes source intersection locks to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_source_event_found[next]", "source_event_found_bridge[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes chunk indices $m_{{ idx}}$ to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_chunk_buffer[next]", "chunk_buffer[next]", "sizeof(long int) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}

                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating Device-to-Device baseline setup of the 9 state vector components for the upcoming payload.
                    // Device-to-Device transfer: Duplicates the initial physical state vector $f^\mu$ to anchor the upcoming RKF45 evaluation.
                    {memcpy_async("d_f_start_bundle[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToDevice", "streams[next]")}
                    // Device-to-Device transfer: Primes the temporary state vector bundle $f^\mu_{{ temp}}$ for the upcoming iterative stage accumulation.
                    {memcpy_async("d_f_temp_bundle[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToDevice", "streams[next]")}
                }}

                for (int stage = 1; stage <= 6; ++stage) {{  // Loop iterator $stage$ executing the 6 discrete stages of the upcoming RKF45 Runge-Kutta numerical solver.
                    // Kernel Launch: Evaluates the metric tensor $g_{{ \mu\nu}}$ and connection $\Gamma^\alpha_{{ \beta\gamma}}$ asynchronously on the alternate stream.
                    interpolation_kernel_{spacetime_name}(commondata, d_f_temp_bundle[next], d_metric_bundle[next], d_connection_bundle[next],active_chunks[next]{stream_arg_next});
                    // Kernel Launch: Computes the geodesic equation right-hand-side derivatives $\dot{{ f}}^\mu$ asynchronously on the alternate stream.
                    calculate_ode_rhs_kernel(d_f_temp_bundle[next], d_metric_bundle[next], d_connection_bundle[next], d_k_bundle[next], stage, active_chunks[next]{stream_arg_next});
                    // Kernel Launch: Accumulates the intermediate RKF45 stage numerical updates asynchronously on the alternate stream.
                    rkf45_stage_update(d_f_start_bundle[next], d_k_bundle[next], d_h[next], stage, active_chunks[next], d_f_temp_bundle[next]{stream_arg_next});
                }}

                // Kernel Launch: Applies Cash-Karp error control to finalize the step-size $h$ and update the upcoming integration baseline.
                rkf45_finalize_and_control(commondata, d_f_bundle[next], d_f_start_bundle[next], d_k_bundle[next], d_h[next], d_status[next], d_affine[next], d_retries[next], active_chunks[next]{stream_arg_next});
                // Kernel Launch: Detects geometric events and records intersection coordinate states asynchronously on the alternate stream.
                event_detection_manager_kernel(commondata, d_f_bundle[next], d_f_prev_bundle[next], d_f_pre_prev_bundle[next], d_affine[next], d_affine_prev[next], d_affine_pre_prev[next], d_results_buffer, d_status[next], d_on_pos_window_prev[next], d_on_pos_source_prev[next], d_window_event_found[next], d_source_event_found[next], d_chunk_buffer[next], active_chunks[next]{stream_arg_next});
                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating Device-to-Host transfer of the 9 upcoming state vector components.
                    // Device-to-Host transfer: Retrieves updated coordinate states $f^\mu$ back to CPU RAM asynchronously on the alternate stream.
                    {memcpy_async("f_bridge[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                    // Device-to-Host transfer: Retrieves updated first derivatives $\dot{{ f}}^\mu$ back to CPU RAM asynchronously on the alternate stream.
                    {memcpy_async("f_p_bridge[next] + c_k * BUNDLE_CAPACITY", "d_f_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                    // Device-to-Host transfer: Retrieves updated second derivatives $\ddot{{ f}}^\mu$ back to CPU RAM asynchronously on the alternate stream.
                    {memcpy_async("f_p_p_bridge[next] + c_k * BUNDLE_CAPACITY", "d_f_pre_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                }}
                // Device-to-Host transfer: Retrieves upcoming active step sizes $h$ back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("h_bridge[next]", "d_h[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming updated status enums back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("status_bridge[next]", "d_status[next]", "sizeof(termination_type_t) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming active rejection counts back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("retries_bridge[next]", "d_retries[next]", "sizeof(int) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming total affine progression $\lambda$ back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("affine_bridge[next]", "d_affine[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming updated window boundary flags back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("on_pos_window_prev_bridge[next]", "d_on_pos_window_prev[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming updated source boundary flags back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("on_pos_source_prev_bridge[next]", "d_on_pos_source_prev[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming historical affine parameter $\lambda_{{ n-1}}$ back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("affine_p_bridge[next]", "d_affine_prev[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming historical affine parameter $\lambda_{{ n-2}}$ back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("affine_p_p_bridge[next]", "d_affine_pre_prev[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming active window locks back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("window_event_found_bridge[next]", "d_window_event_found[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming active source locks back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("source_event_found_bridge[next]", "d_source_event_found[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
            }}

            if (active_chunks[current] > 0) {{ 
                // Device synchronization barrier strictly enforcing completion of the current stream before payload unpacking.
                {stream_sync("streams[current]")}

                for (int fin_i = 0; fin_i < active_chunks[current]; ++fin_i) {{  // Loop iterator $fin_i$ unpacking the finalized physical data back to the global Host matrix.
                    long int m_idx = chunk_buffer[current][fin_i]; // Absolute master index $m_{{ idx}}$ retrieving the specific photon index from the execution chunk.
                    for (int fin_k = 0; fin_k < 9; ++fin_k) {{  // Loop index $fin_k$ retrieving the 9 tensor components back into the global Host matrix.
                        all_photons_host.f[fin_k * num_rays + m_idx] = f_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i]; // Unpacks the synchronized state vector $f^\mu$ into the global Host matrix.
                        all_photons_host.f_p[fin_k * num_rays + m_idx] = f_p_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i]; // Unpacks the synchronized first derivative $\dot{{ f}}^\mu$ into the global Host matrix.
                        all_photons_host.f_p_p[fin_k * num_rays + m_idx] = f_p_p_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i]; // Unpacks the synchronized second derivative $\ddot{{ f}}^\mu$ into the global Host matrix.
                    }}
                    all_photons_host.h[m_idx] = h_bridge[current][fin_i]; // Unpacks the synchronized step size $h$ into the global Host matrix.
                    all_photons_host.status[m_idx] = status_bridge[current][fin_i]; // Unpacks the synchronized trajectory status into the global Host matrix.
                    all_photons_host.rejection_retries[m_idx] = retries_bridge[current][fin_i]; // Unpacks the synchronized rejection count into the global Host matrix.
                    all_photons_host.affine_param[m_idx] = affine_bridge[current][fin_i]; // Unpacks the synchronized affine parameter $\lambda$ into the global Host matrix.
                    all_photons_host.on_positive_side_of_window_prev[m_idx] = on_pos_window_prev_bridge[current][fin_i]; // Unpacks the synchronized window boundary flag into the global Host matrix.
                    all_photons_host.on_positive_side_of_source_prev[m_idx] = on_pos_source_prev_bridge[current][fin_i]; // Unpacks the synchronized source boundary flag into the global Host matrix.
                    all_photons_host.affine_param_p[m_idx] = affine_p_bridge[current][fin_i]; // Unpacks the synchronized historical affine parameter $\lambda_{{ n-1}}$ into the global Host matrix.
                    all_photons_host.affine_param_p_p[m_idx] = affine_p_p_bridge[current][fin_i]; // Unpacks the synchronized historical affine parameter $\lambda_{{ n-2}}$ into the global Host matrix.
                    all_photons_host.window_event_found[m_idx] = window_event_found_bridge[current][fin_i]; // Unpacks the synchronized window lock into the global Host matrix.
                    all_photons_host.source_event_found[m_idx] = source_event_found_bridge[current][fin_i]; // Unpacks the synchronized source lock into the global Host matrix.

                    if (status_bridge[current][fin_i] == ACTIVE) {{  // Evaluates the continuation logic if the trajectory remains within safe physical bounds.
                        int next_s_idx = slot_get_index(&tsm, all_photons_host.f[m_idx]); // Evaluates the updated temporal coordinate $t$ to determine the next operational bin.
                        if (next_s_idx != -1) {{  // Confirms the physical state has not exceeded the maximum simulation time bounds.
                            slot_add_photon(&tsm, next_s_idx, m_idx);  // Re-queues the updated physical state vector back into the host orchestrator.
                        }} else {{ 
                            all_photons_host.status[m_idx] = FAILURE_T_MAX_EXCEEDED; // Flags the physical state as permanently failed due to excessive propagation time.
                            total_active_photons--; // Decrements the global counter as the physical trajectory has reached a terminal state.
                        }}
                    }} else if (status_bridge[current][fin_i] == REJECTED) {{   // Evaluates the retry logic if the numerical step exceeded the requested tolerances.
                        slot_add_photon(&tsm, slot_idx, m_idx); // Re-adds to the current bin to attempt integration with an adapted step-size scalar $h$.
                    }} else {{ 
                        total_active_photons--; // Decrements the global counter as the physical trajectory has reached a terminal state.
                    }}
                }}
                active_chunks[current] = 0; // Clears the execution queue tracker to indicate the active chunk has been fully processed.
            }}

            int temp = current; // Temporary integer scalar storing the primary stream index for logical pointer swapping.
            current = next; // Shifts the primary execution tracker to the alternate stream index.
            next = temp; // Assigns the cleared stream index back to the upcoming payload queue.
        }}

       /* Algorithmic Step: Enforce a rigid hardware sync before advancing the physical time clock. 
        Hardware Justification: Prevents race conditions and ensures all coordinate states $f^\mu$ 
        strictly adhere to the current temporal bin limits. */
        // --- PHASE C: THE TIME BARRIER ---
        BHAH_DEVICE_SYNC(); // Synchronize global state convergence before advancing the central time engine. 

     }}

    /* Algorithmic Step: Process terminal photon trajectories and extract final geometric intersections. Hardware Justification: Memory transfers execute synchronously on the host to ensure all asynchronous execution pipelines have fully resolved. */
    // --- 4. CLEANUP & FINALIZATION ---
        
        // Device-to-Host transfer: Extracts validated device-native blueprints $b_i$ containing geometric plane intersections.
        {results_memcpy} 

        // Kernel Launch: Processes escaped photons intersecting the celestial sphere $r > r_{{escape}}$.
        {calc_blueprint}

        // Loop iterator $s$ purging the double-buffered arrays across both hardware streams.
        for (int s = 0; s < 2; ++s) {{
            // Host Memory Free: Purges bridge components supporting scatter logic mapped to PCIe DMA transfers.
            {free_pinned}(chunk_buffer[s]); // Purges the execution chunk mapping bridge.
            {free_pinned}(f_bridge[s]); // Purges the state vector $f^\mu$ bridge.
            {free_pinned}(f_p_bridge[s]); // Purges the first derivative $\dot{{f}}^\mu$ bridge.
            {free_pinned}(f_p_p_bridge[s]); // Purges the second derivative $\ddot{{f}}^\mu$ bridge.
            {free_pinned}(affine_bridge[s]); // Purges the affine parameter $\lambda$ bridge.
            {free_pinned}(h_bridge[s]); // Purges the integration step size $h$ bridge.
            {free_pinned}(status_bridge[s]); // Purges the trajectory status bridge.
            {free_pinned}(retries_bridge[s]); // Purges the error rejection scalar bridge.
            {free_pinned}(on_pos_window_prev_bridge[s]); // Purges the observer window boundary flag bridge.
            {free_pinned}(on_pos_source_prev_bridge[s]); // Purges the source emission boundary flag bridge.
            {free_pinned}(affine_p_bridge[s]); // Purges the historical affine parameter $\lambda_{{n-1}}$ bridge.
            {free_pinned}(affine_p_p_bridge[s]); // Purges the historical affine parameter $\lambda_{{n-2}}$ bridge.
            {free_pinned}(window_event_found_bridge[s]); // Purges the observer window intersection lock bridge.
            {free_pinned}(source_event_found_bridge[s]); // Purges the source emission intersection lock bridge.

            // Device Memory Free: Purges remaining VRAM operational pipeline scratchpads.
            {free_device}(d_f_bundle[s]); // Purges the state vector $f^\mu$ scratchpad.
            {free_device}(d_f_start_bundle[s]); // Purges the anchor state vector $f_{{start}}$ scratchpad.
            {free_device}(d_f_temp_bundle[s]); // Purges the temporary stage $f^\mu_{{temp}}$ scratchpad.
            {free_device}(d_f_prev_bundle[s]); // Purges the history state $f^\mu_{{n-1}}$ scratchpad.
            {free_device}(d_f_pre_prev_bundle[s]); // Purges the history state $f^\mu_{{n-2}}$ scratchpad.
            {free_device}(d_metric_bundle[s]); // Purges the symmetric metric tensor $g_{{\mu\nu}}$ scratchpad.
            {free_device}(d_connection_bundle[s]); // Purges the Christoffel symbols $\Gamma^\alpha_{{\beta\gamma}}$ scratchpad.
            {free_device}(d_k_bundle[s]); // Purges the derivative tensor $\dot{{f}}^\mu$ scratchpad.
            {free_device}(d_h[s]); // Purges the active integration step sizing $h$ scratchpad.
            {free_device}(d_affine[s]); // Purges the total affine parameter progress $\lambda$ scratchpad.
            {free_device}(d_status[s]); // Purges the current trajectory status limit scratchpad.
            {free_device}(d_retries[s]); // Purges the sequential error rejection scratchpad.
            {free_device}(d_on_pos_window_prev[s]); // Purges the previous observer window boundary side scratchpad.
            {free_device}(d_on_pos_source_prev[s]); // Purges the previous source emission boundary side scratchpad.
            {free_device}(d_affine_prev[s]); // Purges the historical affine parameter $\lambda_{{n-1}}$ scratchpad.
            {free_device}(d_affine_pre_prev[s]); // Purges the historical affine parameter $\lambda_{{n-2}}$ scratchpad.
            {free_device}(d_window_event_found[s]); // Purges the window intersection coordinate guard scratchpad.
            {free_device}(d_source_event_found[s]); // Purges the source intersection coordinate guard scratchpad.
            {free_device}(d_chunk_buffer[s]); // Purges the absolute master indices $m_{{idx}}$ mapping scratchpad.

            {stream_destroy}
        }}

        /* Algorithmic Step: Evaluate relative numerical drift on the CPU. Hardware Justification: Evaluating relative numerical drift natively on the CPU prevents VRAM bottlenecks and leverages complex print formatting. */
        // --- CPU CONSERVATION DRIFT EVALUATION ---
        if (commondata->perform_conservation_check) {{
            // Kernel Launch: Calculates terminal conserved quantities natively on the executing architecture.
            calculate_conserved_quantities_universal_{spacetime_name}_photon(commondata, &all_photons_host, num_rays, final_cq_host);

            printf("\n=================================================\n");
            printf(" CONSERVED QUANTITIES DIAGNOSTIC REPORT\n");
            printf("=================================================\n");
            
            // Scalar variables tracking the maximum recorded relative drift for energy $E$, angular momentum $L_z$, and Carter constant $Q$.
            double max_err_E = 0.0, max_err_Lz = 0.0, max_err_Q = 0.0;
            // Absolute master indices $m_{{idx}}$ identifying the trajectory responsible for the maximum numerical drift in each respective quantity.
            long int worst_ray_E = -1, worst_ray_Lz = -1, worst_ray_Q = -1;

            // Loop iterator $i$ spanning the global dataset to calculate relative errors natively on the CPU.
            for (long int i = 0; i < num_rays; i++) {{
                double err_E = fabs((final_cq_host[i].E - initial_cq_host[i].E) / (initial_cq_host[i].E + 1e-15)); // Evaluates the relative numerical drift for energy $E$.
                double err_Lz = fabs((final_cq_host[i].Lz - initial_cq_host[i].Lz) / (initial_cq_host[i].Lz + 1e-15)); // Evaluates the relative numerical drift for angular momentum $L_z$.
                double err_Q = fabs((final_cq_host[i].Q - initial_cq_host[i].Q) / (initial_cq_host[i].Q + 1e-15)); // Evaluates the relative numerical drift for Carter constant $Q$.

                if (err_E > max_err_E) {{ max_err_E = err_E; worst_ray_E = i; }} // Updates the maximum tracked relative error and associated index $i$ for energy $E$.
                if (err_Lz > max_err_Lz) {{ max_err_Lz = err_Lz; worst_ray_Lz = i; }} // Updates the maximum tracked relative error and associated index $i$ for angular momentum $L_z$.
                if (err_Q > max_err_Q) {{ max_err_Q = err_Q; worst_ray_Q = i; }} // Updates the maximum tracked relative error and associated index $i$ for Carter constant $Q$.
            }}

            printf("  Max Relative Error (Energy E): %e (Ray %ld)\n", max_err_E, worst_ray_E); // Output block printing the maximum relative error for energy $E$.
            printf("  Max Relative Error (Momentum Lz): %e (Ray %ld)\n", max_err_Lz, worst_ray_Lz); // Output block printing the maximum relative error for angular momentum $L_z$.
            printf("  Max Relative Error (Carter Q): %e (Ray %ld)\n", max_err_Q, worst_ray_Q); // Output block printing the maximum relative error for Carter constant $Q$.
            printf("=================================================\n\n");

            // Host Memory Free: Purges pinned diagnostic buffers.
            {free_pinned}(initial_cq_host); // Purges pinned initial diagnostic data buffer.
            {free_pinned}(final_cq_host); // Purges pinned final diagnostic data buffer.
        }}

        // Host Memory Free: Purges the primary Host array states $f^\mu$ and affine parameters $\lambda$.
        {free_pinned}(all_photons_host.f); // Purges the primary Host array state $f^\mu$.
        {free_pinned}(all_photons_host.f_p); // Purges the primary Host array first derivative $\dot{{f}}^\mu$.
        {free_pinned}(all_photons_host.f_p_p); // Purges the primary Host array second derivative $\ddot{{f}}^\mu$.
        {free_pinned}(all_photons_host.affine_param); // Purges the primary Host array affine parameter $\lambda$.
        {free_pinned}(all_photons_host.h); // Purges the primary Host array integration step size $h$.
        {free_pinned}(all_photons_host.status); // Purges the primary Host array trajectory status enum.
        {free_pinned}(all_photons_host.rejection_retries); // Purges the primary Host array error rejection scalar.
        {free_pinned}(all_photons_host.on_positive_side_of_window_prev); // Purges the primary Host array observer window boundary flag.
        {free_pinned}(all_photons_host.on_positive_side_of_source_prev); // Purges the primary Host array source emission boundary flag.
        {free_pinned}(all_photons_host.affine_param_p); // Purges the primary Host array historical affine parameter $\lambda_{{n-1}}$.
        {free_pinned}(all_photons_host.affine_param_p_p); // Purges the primary Host array historical affine parameter $\lambda_{{n-2}}$.
        {free_pinned}(all_photons_host.window_event_found); // Purges the primary Host array observer window intersection lock.
        {free_pinned}(all_photons_host.source_event_found); // Purges the primary Host array source emission intersection lock.
        
        // Device Memory Free: Purges the final single-pointer intersection blueprint buffer $b_i$.
        {free_device}(d_results_buffer); // Purges the intersection blueprint buffer $b_i$.

        // Memory Free: Purges the temporal sorting struct mapping the Host-side execution grid.
        slot_manager_free(&tsm); // Purges the central time slot orchestrator.
    """



    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body
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
