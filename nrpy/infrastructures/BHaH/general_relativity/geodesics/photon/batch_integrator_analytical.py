# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_analytical.py
r"""
Defines the orchestration module for the analytical integration pipeline.

This module structures the C orchestrator responsible for managing the life cycle of
photon trajectories $x^\mu$ in analytical spacetimes. It implements a split-pipeline
architecture, decoupling the Runge-Kutta-Fehlberg 4(5) integration kernels.
Intermediate tensors and state vectors are persisted in flattened Structure of Arrays
bundles. The orchestrator manages asynchronous data transfers and parallel execution
streams.

Double-buffered streams hide data transfer latency behind the integration compute
time. The primary execution stream dynamically computes the metric and checks the
normalization constraint without persistent memory allocation. Pinned memory increases
transfer bandwidth for diagnostic arrays. Evaluating relative numerical drift natively
prevents device memory bottlenecks. Pre-computing states prior to the temporal loop
enables coalesced memory access during iterative integration. Diagnostic probes verify
that the logic populated coordinates and complete metric-null initial momenta.
Intercepting unphysical spacetime regions after calculation
ensures constraint solver convergence. Evaluating conserved quantities establishes a
data baseline before the pipeline mutates the state vectors.

Each invocation is one tile batch: it receives the caller-computed ``num_rays``,
interpolates the analytic metric once at the observer event, constructs one
observer tetrad, and passes that tetrad to all rays in the invocation. Image
sample coordinates remain serialization metadata and never enter the evolving
``PhotonStateSoA``.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.time_slot_manager_helpers import (
    time_slot_manager_helpers,
)


def batch_integrator_analytical(spacetime_name: str) -> None:
    r"""
    Construct the Native CUDA orchestrator for the batched analytical integration pipeline.

    :param spacetime_name: The identifier for the spacetime metric (e.g., 'KerrSchild').
    """
    if "time_slot_manager" not in par.glb_extras_dict.get("BHaH_defines", {}):
        time_slot_manager_helpers()

    # Core physics and numerical simulation parameters for the global spacetime struct.
    par.register_CodeParameters(
        "REAL",
        __name__,
        ["r_escape"],
        [150.0],
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

    # Initialize the complete serialized record once before any event writes.
    init_arg_dict = {
        "num_rays": "const long int",
        "d_results_buffer": "blueprint_data_t *restrict",
    }
    if parallelization == "cuda":
        init_loop = """
        const long int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= num_rays) return;
        """
        init_loop_end = ""
        init_launch_dict = {
            "threads_per_block": ["256", "1", "1"],
            "blocks_per_grid": ["(num_rays + 255) / 256", "1", "1"],
        }
    else:
        init_loop = """
        #pragma omp parallel for
        for (long int i = 0; i < num_rays; ++i) {
        """
        init_loop_end = "        }\n"
        init_launch_dict = None
    init_kernel_body = f"""
    {init_loop}
        d_results_buffer[i].termination_type = FAILURE_GENERIC;
        d_results_buffer[i].y_nt = 0.0;
        d_results_buffer[i].z_nt = 0.0;
        d_results_buffer[i].y_t = 0.0;
        d_results_buffer[i].z_t = 0.0;
        d_results_buffer[i].final_theta = 0.0;
        d_results_buffer[i].final_phi = 0.0;
        d_results_buffer[i].non_terminal_plane_lambda = 0.0;
        d_results_buffer[i].non_terminal_plane_t = 0.0;
        d_results_buffer[i].L_f = 0.0;
        d_results_buffer[i].t_f = 0.0;
    {init_loop_end}
    """
    init_prefunc, init_launch = parallel_utils.generate_kernel_and_launch_code(
        kernel_name="initialize_blueprint_data",
        kernel_body=init_kernel_body,
        arg_dict_cuda=init_arg_dict,
        arg_dict_host=init_arg_dict,
        parallelization=parallelization,
        launch_dict=init_launch_dict,
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    if parallelization == "cuda":
        includes.extend(
            ["cuda_runtime.h", "cuda_intrinsics.h", "BHaH_global_device_defines.h"]
        )

    desc = r""" Central Host-bound CPU orchestrator for the batched Split-Pipeline relativistic ray tracing loop.

    This function acts as the primary loop for evaluating photon geodesics $x^\mu$.
    It utilizes a TimeSlotManager to bin active rays by their physical coordinate time $t$.
    Its integration parameter is the physical affine parameter $\lambda$.
    The Split-Pipeline architecture maps mathematical tensors like $g_{\mu\nu}$ and $\Gamma^\alpha_{\beta\gamma}$ to memory scratchpads.
    Mapping tensors to memory respects the hardware limit per thread on modern architectures.

    @param[in] commondata Struct containing global spacetime and integration tolerances.
    @param num_rays Number of rays in this tile batch. The caller derives this
                     from scan-density sampling; this integrator does not
                     calculate image placement or carry image metadata in the
                     evolving PhotonStateSoA.
    @param[out] results_buffer Device array storing the final physical intersections."""

    cfunc_type = "void"

    name = "batch_integrator_analytical"

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
        bridge_alloc_comment = (
            "Allocate memory arrays in Host RAM for the structural bridge payloads."
        )
        scratch_alloc_comment = (
            "Allocate 1D Host Scratchpad arrays for temporal data staging."
        )

    # This block configures the hardware-specific teardown and memory synchronization sequences.
    if parallelization == "cuda":
        results_memcpy = "cudaMemcpy(results_buffer, d_results_buffer, sizeof(blueprint_data_t) * num_rays, cudaMemcpyDeviceToHost);"
        calc_blueprint = "calculate_and_fill_blueprint_data_universal(&all_photons_host, num_rays, results_buffer, NULL, NULL, 0);"
        set_initial_conditions_call = " set_initial_conditions_kernel(commondata, num_rays, &all_photons_host, observer_metric, observer_tetrad);"
        stream_destroy = "cudaStreamDestroy(streams[s]); // Purges the hardware stream execution context."
        free_device = "BHAH_FREE_DEVICE"
        free_pinned = "BHAH_FREE_PINNED"
    else:
        results_memcpy = "memcpy(results_buffer, d_results_buffer, sizeof(blueprint_data_t) * num_rays);"
        calc_blueprint = "calculate_and_fill_blueprint_data_universal(&all_photons_host, num_rays, results_buffer, NULL, NULL, 0);"
        set_initial_conditions_call = " set_initial_conditions_kernel(commondata, num_rays, &all_photons_host, observer_metric, observer_tetrad);"
        stream_destroy = (
            "// Stream destruction natively omitted for synchronous CPU execution."
        )
        free_device = "BHAH_FREE"
        free_pinned = "BHAH_FREE"

    # Establish memory transfer protocols based on the target parallelization architecture.
    if parallelization == "cuda":
        # Generates asynchronous PCIe transfer commands for CUDA architectures.
        def memcpy_async(
            dest: str, src: str, size: str, direction: str, stream: str
        ) -> str:
            return f"cudaMemcpyAsync({dest}, {src}, {size}, {direction}, {stream});"

        # Generates hardware synchronization barriers for CUDA streams.
        def stream_sync(stream: str) -> str:
            return f"cudaStreamSynchronize({stream});"

        stream_arg = ", 0"
        free_pinned = "BHAH_FREE_PINNED"
    else:
        # Generates synchronous Host memory copies for native OpenMP architectures.
        def memcpy_async(
            dest: str, src: str, size: str, direction: str, stream: str
        ) -> str:
            # pylint: disable=unused-argument
            return f"memcpy({dest}, {src}, {size});"

        def stream_sync(stream: str) -> str:
            # pylint: disable=unused-argument
            return r"// Stream synchronization is natively omitted for synchronous CPU execution."

        stream_arg = ", 0"
        free_pinned = "BHAH_FREE"

    # Define stream args for CUDA
    stream_arg_current = ", current" if parallelization == "cuda" else ", current"
    stream_arg_next = ", next" if parallelization == "cuda" else ", next"
    initial_state_value = "commondata->t_start"

    if parallelization == "cuda":
        observer_state_transfer = "\n".join(
            [
                "        cudaMemcpy("
                "d_f_bundle[0] + observer_component * BUNDLE_CAPACITY, "
                "all_photons_host.f + observer_component * num_rays, "
                "sizeof(double), cudaMemcpyHostToDevice);"
            ]
        )
        observer_metric_transfer = "\n".join(
            [
                "        cudaMemcpy("
                "observer_metric + observer_component, "
                "d_metric_bundle[0] + observer_component * BUNDLE_CAPACITY, "
                "sizeof(double), cudaMemcpyDeviceToHost);"
            ]
        )
        observer_metric_wait = "        cudaStreamSynchronize(0);"
    else:
        observer_state_transfer = "\n".join(
            [
                "        memcpy("
                "d_f_bundle[0] + observer_component * BUNDLE_CAPACITY, "
                "all_photons_host.f + observer_component * num_rays, "
                "sizeof(double));"
            ]
        )
        observer_metric_transfer = "\n".join(
            [
                "        memcpy("
                "observer_metric + observer_component, "
                "d_metric_bundle[0] + observer_component * BUNDLE_CAPACITY, "
                "sizeof(double));"
            ]
        )
        observer_metric_wait = ""

    body = rf"""
    //==========================================
    // 1. HOST & DEVICE ALLOCATION
    //==========================================

    // The master host-side Structure of Arrays (SoA) tracking all photons $f^\mu$.
    PhotonStateSoA all_photons_host;

    // {pin_comment} the state vector $f^\mu$.
    {malloc_pinned}(all_photons_host.f, sizeof(double) * 9 * num_rays);
    // {pin_comment} the first derivative $\dot{{f}}^\mu$.
    {malloc_pinned}(all_photons_host.f_p, sizeof(double) * 9 * num_rays);
    // {pin_comment} the second derivative $\ddot{{f}}^\mu$.
    {malloc_pinned}(all_photons_host.f_p_p, sizeof(double) * 9 * num_rays);
    // {pin_comment} the physical integration parameter $\lambda$.
    {malloc_pinned}(all_photons_host.integration_param, sizeof(double) * num_rays);
    // {pin_comment} individual integration step sizes $h$.
    {malloc_pinned}(all_photons_host.h, sizeof(double) * num_rays);
    // {pin_comment} the trajectory termination status.
    {malloc_pinned}(all_photons_host.status, sizeof(termination_type_t) * num_rays);
    // {pin_comment} the number of step-size rejections.
    {malloc_pinned}(all_photons_host.rejection_retries, sizeof(int) * num_rays);
    // {pin_comment} the previous nonterminal plane boundary state.
    {malloc_pinned}(all_photons_host.on_positive_side_of_non_terminal_plane_prev, sizeof(bool) * num_rays);
    // {pin_comment} the previous terminal-plane boundary state.
    {malloc_pinned}(all_photons_host.on_positive_side_of_terminal_plane_prev, sizeof(bool) * num_rays);
    // {pin_comment} the history step $\lambda_{{n-1}}$.
    {malloc_pinned}(all_photons_host.integration_param_p, sizeof(double) * num_rays);
    // {pin_comment} the history step $\lambda_{{n-2}}$.
    {malloc_pinned}(all_photons_host.integration_param_p_p, sizeof(double) * num_rays);
    // {pin_comment} the nonterminal plane intersection lock.
    {malloc_pinned}(all_photons_host.non_terminal_plane_event_found, sizeof(bool) * num_rays);
    // {pin_comment} the terminal-plane intersection lock.
    {malloc_pinned}(all_photons_host.terminal_plane_event_found, sizeof(bool) * num_rays);

    {stream_setup_str}

    //==========================================
    // DOUBLE-BUFFERED BRIDGE ARRAYS
    //==========================================
    // Extraction buffer used by the TimeSlotManager to map sparse indices to contiguous execution blocks.
    long int *chunk_buffer[2];
    // Bridge array staging the state vector $f^\mu$ for memory transfers.
    double *f_bridge[2];
    // Bridge array staging the first derivative $\dot{{f}}^\mu$ for memory transfers.
    double *f_p_bridge[2];
    // Bridge array staging the second derivative $\ddot{{f}}^\mu$ for memory transfers.
    double *f_p_p_bridge[2];
    // Bridge array staging the integration parameter $\lambda$ for memory transfers.
    double *integration_param_bridge[2];
    // Bridge array staging the current integration step size $h$ for memory transfers.
    double *h_bridge[2];
    // Bridge array staging the current trajectory termination status for memory transfers.
    termination_type_t *status_bridge[2];
    // Bridge array staging the number of step-size rejections for memory transfers.
    int *retries_bridge[2];
    // Bridge array staging the previous nonterminal plane boundary side flag for memory transfers.
    bool *on_pos_non_terminal_plane_prev_bridge[2];
    // Bridge array staging the previous terminal-plane boundary side flag for memory transfers.
    bool *on_pos_terminal_plane_prev_bridge[2];
    // Bridge array staging the historical integration parameter $\lambda_{{n-1}}$ for chunked memory transfers.
    double *integration_param_p_bridge[2];
    // Bridge array staging the historical integration parameter $\lambda_{{n-2}}$ for chunked memory transfers.
    double *integration_param_p_p_bridge[2];
    // Bridge array staging the nonterminal plane event lock for memory transfers.
    bool *non_terminal_plane_event_found_bridge[2];
    // Bridge array staging the terminal-plane event lock for memory transfers.
    bool *terminal_plane_event_found_bridge[2];

    //==========================================
    // DOUBLE-BUFFERED VRAM SCRATCHPADS
    //==========================================
    // Scratchpad array holding the physical normalization diagnostic outputs.
    normalization_constraint_t *d_norm_bundle[2];
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
    // Array regulating total integration parameter progress $\lambda$.
    double *d_integration_param[2];
    // Array holding the current trajectory status limits.
    termination_type_t *d_status[2];
    // Array tracking sequential error rejections per photon.
    int *d_retries[2];
    // Array flagging the previous nonterminal plane boundary side.
    bool *d_on_pos_non_terminal_plane_prev[2];
    // Array flagging the previous terminal-plane boundary side.
    bool *d_on_pos_terminal_plane_prev[2];
    // Array tracking historical integration parameter $\lambda_{{n-1}}$.
    double *d_integration_param_prev[2];
    // Array tracking historical integration parameter $\lambda_{{n-2}}$.
    double *d_integration_param_pre_prev[2];
    // Array guarding nonterminal-plane intersection coordinates from multi-trigger overwrites.
    bool *d_non_terminal_plane_event_found[2];
    // Array guarding the terminal-plane intersection coordinates from multi-trigger overwrites.
    bool *d_terminal_plane_event_found[2];
    // Array carrying the absolute master indices $m_{{idx}}$ mapping the execution chunk.
    long int *d_chunk_buffer[2];

    // Loop iterator for instantiating the double-buffered operational arrays.
    for (int s = 0; s < 2; ++s) {{
        // {bridge_alloc_comment}
        {malloc_pinned}(chunk_buffer[s], sizeof(long int) * BUNDLE_CAPACITY); // Pin chunk buffers.
        {malloc_pinned}(f_bridge[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Pin $f^\mu$ bridges.
        {malloc_pinned}(f_p_bridge[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Pin $\dot{{f}}^\mu$ bridges.
        {malloc_pinned}(f_p_p_bridge[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Pin $\ddot{{f}}^\mu$ bridges.
        {malloc_pinned}(integration_param_bridge[s], sizeof(double) * BUNDLE_CAPACITY); // Pin $\lambda$ bridges.
        {malloc_pinned}(h_bridge[s], sizeof(double) * BUNDLE_CAPACITY); // Pin $h$ bridges.
        {malloc_pinned}(status_bridge[s], sizeof(termination_type_t) * BUNDLE_CAPACITY); // Pin status bridges.
        {malloc_pinned}(retries_bridge[s], sizeof(int) * BUNDLE_CAPACITY); // Pin retries bridges.
        {malloc_pinned}(on_pos_non_terminal_plane_prev_bridge[s], sizeof(bool) * BUNDLE_CAPACITY); // Pin nonterminal-plane flag bridges.
        {malloc_pinned}(on_pos_terminal_plane_prev_bridge[s], sizeof(bool) * BUNDLE_CAPACITY); // Pin terminal-plane side-flag bridges.
        {malloc_pinned}(integration_param_p_bridge[s], sizeof(double) * BUNDLE_CAPACITY); // Pin $\lambda_{{n-1}}$ bridges.
        {malloc_pinned}(integration_param_p_p_bridge[s], sizeof(double) * BUNDLE_CAPACITY); // Pin $\lambda_{{n-2}}$ bridges.
        {malloc_pinned}(non_terminal_plane_event_found_bridge[s], sizeof(bool) * BUNDLE_CAPACITY); // Pin nonterminal-plane lock bridges.
        {malloc_pinned}(terminal_plane_event_found_bridge[s], sizeof(bool) * BUNDLE_CAPACITY); // Pin terminal-plane lock bridges.

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
        {malloc_device}(d_integration_param[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate $\lambda$ scratchpad.
        {malloc_device}(d_status[s], sizeof(termination_type_t) * BUNDLE_CAPACITY); // Allocate status scratchpad.
        {malloc_device}(d_retries[s], sizeof(int) * BUNDLE_CAPACITY); // Allocate retries scratchpad.
        {malloc_device}(d_on_pos_non_terminal_plane_prev[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate nonterminal-plane flag scratchpad.
        {malloc_device}(d_on_pos_terminal_plane_prev[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate terminal-plane side-flag scratchpad.
        {malloc_device}(d_integration_param_prev[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate $\lambda_{{n-1}}$ scratchpad.
        {malloc_device}(d_integration_param_pre_prev[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate $\lambda_{{n-2}}$ scratchpad.
        {malloc_device}(d_non_terminal_plane_event_found[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate nonterminal-plane lock scratchpad.
        {malloc_device}(d_terminal_plane_event_found[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate terminal-plane lock scratchpad.
        {malloc_device}(d_chunk_buffer[s], sizeof(long int) * BUNDLE_CAPACITY); // Allocate chunk mapping scratchpad.
        {malloc_device}(d_norm_bundle[s], sizeof(normalization_constraint_t) * BUNDLE_CAPACITY); // Allocate diagnostic outputs scratchpad.
    }} // END LOOP: for s over 2

    // Device-native struct pointer storing the final physical plane intersections.
    blueprint_data_t *d_results_buffer;
    // {dev_comment} the blueprint results buffer to avoid mid-computation memory transfers.
    {malloc_device}(d_results_buffer, sizeof(blueprint_data_t) * num_rays);
    {init_launch}

    // Host-bound struct managing temporal binning of photon trajectories $x^\mu$.
    TimeSlotManager tsm;
    // Initializes the temporal manager strictly on the CPU to coordinate the Split-Pipeline batches.
    slot_manager_init(&tsm, commondata->slot_manager_t_min, commondata->t_start + 1.0, commondata->slot_manager_delta_t, num_rays);

    //==========================================
    // DIAGNOSTIC MEMORY ALLOCATION
    //==========================================
    // Host pointer tracking the initial conserved quantities prior to integration.
    conserved_quantities_t *initial_cq_host = NULL;
    // Host pointer tracking the terminal conserved quantities post integration.
    conserved_quantities_t *final_cq_host = NULL;

    if (commondata->perform_conservation_check) {{
        // {pin_comment} the initial diagnostic data.
        {malloc_pinned}(initial_cq_host, sizeof(conserved_quantities_t) * num_rays);
        // {pin_comment} the final diagnostic data.
        {malloc_pinned}(final_cq_host, sizeof(conserved_quantities_t) * num_rays);
    }} // END IF: commondata->perform_conservation_check to allocate diagnostic arrays

    //==========================================
    // 2. INITIALIZATION PHASE
    //==========================================
    // Evaluate initial coordinate states and map global spacetime metrics to memory bounds.

    double observer_metric[10]; // One covariant metric evaluated at the observer event.
    // One contravariant tetrad constructed for this tile-batch call and passed
    // unchanged to every ray in the current batch.
    double observer_tetrad[4][4];

    // Build a temporary observer state with zero momentum.  The analytic metric
    // interpolation reads only the event coordinates and is executed once.
    for (int observer_component = 0; observer_component < 9; ++observer_component) {{
        all_photons_host.f[observer_component * num_rays] = 0.0;
    }} // END LOOP: clear temporary observer state
    all_photons_host.f[0 * num_rays] = commondata->t_start;
    all_photons_host.f[1 * num_rays] = commondata->observer_x;
    all_photons_host.f[2 * num_rays] = commondata->observer_y;
    all_photons_host.f[3 * num_rays] = commondata->observer_z;

    for (int observer_component = 0; observer_component < 9; ++observer_component) {{
{observer_state_transfer}
    }} // END LOOP: transfer temporary observer state
    interpolation_kernel_{spacetime_name}(
        commondata,
        d_f_bundle[0],
        d_metric_bundle[0],
        NULL,
        1{stream_arg});
{observer_metric_wait}
    for (int observer_component = 0; observer_component < 10; ++observer_component) {{
{observer_metric_transfer}
    }} // END LOOP: retrieve one analytic observer metric
    {observer_metric_wait}
    for (int observer_component = 0; observer_component < 10; ++observer_component) {{
        if (!isfinite(observer_metric[observer_component])) {{
            fprintf(stderr, "ERROR: analytic observer metric component %d was not finite.\n", observer_component);
            exit(EXIT_FAILURE);
        }} // END IF: analytic observer metric invalid
    }} // END LOOP: validate analytic observer metric

    // Operates synchronously as the primary state array must be fully populated before pipeline dispatch.
    {set_initial_conditions_call}

    //==========================================
    // INITIAL-STATE FINITE AND NULL DIAGNOSTICS
    //==========================================
    double observer_metric_scale = 1.0;
    for (int metric_component = 0; metric_component < 10; ++metric_component) {{
        observer_metric_scale =
            fmax(observer_metric_scale, fabs(observer_metric[metric_component]));
    }} // END LOOP: scale one analytic observer metric
    long int invalid_initial_rays = 0;
    double max_initial_null_residual = 0.0;

    for (long int ray = 0; ray < num_rays; ++ray) {{
        const double t_value = all_photons_host.f[0 * num_rays + ray];
        const double x_value = all_photons_host.f[1 * num_rays + ray];
        const double y_value = all_photons_host.f[2 * num_rays + ray];
        const double z_value = all_photons_host.f[3 * num_rays + ray];
        const double p0_value = all_photons_host.f[4 * num_rays + ray];
        const double p1_value = all_photons_host.f[5 * num_rays + ray];
        const double p2_value = all_photons_host.f[6 * num_rays + ray];
        const double p3_value = all_photons_host.f[7 * num_rays + ray];
        const double affine_value = all_photons_host.f[8 * num_rays + ray];
        const bool finite_state =
            isfinite(t_value) && isfinite(x_value) && isfinite(y_value) &&
            isfinite(z_value) && isfinite(p0_value) && isfinite(p1_value) &&
            isfinite(p2_value) && isfinite(p3_value) && isfinite(affine_value);

        const double null_residual =
            observer_metric[0] * p0_value * p0_value +
            2.0 * observer_metric[1] * p0_value * p1_value +
            2.0 * observer_metric[2] * p0_value * p2_value +
            2.0 * observer_metric[3] * p0_value * p3_value +
            observer_metric[4] * p1_value * p1_value +
            2.0 * observer_metric[5] * p1_value * p2_value +
            2.0 * observer_metric[6] * p1_value * p3_value +
            observer_metric[7] * p2_value * p2_value +
            2.0 * observer_metric[8] * p2_value * p3_value +
            observer_metric[9] * p3_value * p3_value;
        const double momentum_scale =
            fmax(1.0, fabs(p0_value) + fabs(p1_value) +
                           fabs(p2_value) + fabs(p3_value));
        const double null_tolerance =
            1.0e-9 * observer_metric_scale * momentum_scale * momentum_scale;
        const double null_error = fabs(null_residual);
        max_initial_null_residual =
            fmax(max_initial_null_residual, null_error);

        const bool coordinate_state_matches =
            fabs(t_value - {initial_state_value}) <= 1.0e-10 &&
            fabs(x_value - commondata->observer_x) <= 1.0e-10 &&
            fabs(y_value - commondata->observer_y) <= 1.0e-10 &&
            fabs(z_value - commondata->observer_z) <= 1.0e-10 &&
            fabs(affine_value) <= 1.0e-15;
        if (!finite_state || !coordinate_state_matches ||
            !isfinite(null_residual) || null_error > null_tolerance) {{
            ++invalid_initial_rays;
        }}
    }} // END LOOP: validate analytic initial rays

    if (invalid_initial_rays > 0) {{
        fprintf(
            stderr,
            "ERROR: %ld analytic initial rays failed finite, coordinate, or null checks.\n",
            invalid_initial_rays);
        exit(EXIT_FAILURE);
    }}
    printf(
        "[DIAGNOSTIC] Analytic initialization maximum null residual: %.17e\n",
        max_initial_null_residual);


    //==========================================
    // BASELINE CONSERVED QUANTITIES
    //==========================================
    // Evaluate initial conserved quantities immediately after generating valid physical null states.

    if (commondata->perform_conservation_check) {{
        // Executes global conserved quantities natively via chunked device parameters before pipeline processing.
        calculate_conserved_quantities_universal_{spacetime_name}_photon(commondata, &all_photons_host, num_rays, initial_cq_host);
    }} // END IF: perform_conservation_check to evaluate baseline conserved

    long int sync_i; // Loop iterator index $sync_i$ spanning the entire global ray count to synchronize starting properties across history states.
    for(sync_i = 0; sync_i < num_rays; ++sync_i) {{
        int sync_k; // Loop index $sync_k$ iterating over the 9 tensor components to populate the historical derivatives $\dot{{f}}^\mu$ and $\ddot{{f}}^\mu$.
        for (sync_k = 0; sync_k < 9; ++sync_k) {{
            all_photons_host.f_p[sync_k * num_rays + sync_i] = all_photons_host.f[sync_k * num_rays + sync_i]; // Propagates the initial coordinate state vector $f^\mu$ to the first history derivative matrix.
            all_photons_host.f_p_p[sync_k * num_rays + sync_i] = all_photons_host.f[sync_k * num_rays + sync_i]; // Propagates the initial coordinate state vector $f^\mu$ to the second history derivative matrix.
        }} // END LOOP: for sync_k over 9
        all_photons_host.status[sync_i] = ACTIVE; // Assigns the initial trajectory activity enum for the global physics engine.
        all_photons_host.integration_param[sync_i] = 0.0; // Sets the initial baseline progression scalar for the integration parameter $\lambda$.
        all_photons_host.rejection_retries[sync_i] = 0; // Clears the error rejection scalar to initialize the step size convergence tracking.

        all_photons_host.integration_param_p[sync_i] = 0.0; // Initializes the first historical integration parameter $\lambda_{{n-1}}$.
        all_photons_host.integration_param_p_p[sync_i] = 0.0; // Initializes the second historical integration parameter $\lambda_{{n-2}}$.
        all_photons_host.non_terminal_plane_event_found[sync_i] = false; // Sets the nonterminal plane intersection logical lock to false.
        all_photons_host.terminal_plane_event_found[sync_i] = false; // Sets the terminal-plane intersection logical lock to false.

        int s_idx = slot_get_index(&tsm, all_photons_host.f[sync_i]); // Integer index $s_{{idx}}$ mapping the current photon's temporal coordinate $t$ to a discrete execution bin in the TimeSlotManager.
        if (s_idx != -1) {{
            slot_add_photon(&tsm, s_idx, sync_i); // Registers the active photon index to its corresponding temporal bin mapped by the orchestrator.
        }} // END IF: s_idx != -1 to add
    }} // END LOOP: for sync_i over num_rays

    // Hardware clock state marking the beginning of the active integration chunk.
    struct timespec batch_start_time;
    // Hardware clock state marking the conclusion of the active integration chunk.
    struct timespec batch_end_time;
    // Captures the initial hardware clock state prior to temporal loop execution.
    clock_gettime(CLOCK_MONOTONIC, &batch_start_time);

    // Allocates vertical terminal space for the dynamic multi-line progress dashboard.
    printf("\n\n\n\n\n\n\n");

    //==========================================
    // 3. TEMPORAL LOOP (The Engine)
    //==========================================
    // Integer tracking the global number of active photon trajectories to allow early loop termination.
    long int total_active_photons = num_rays;

    // Outer loop iterator for the physical time bins.
    for (int slot_idx = tsm.num_slots - 1; slot_idx >= 0; --slot_idx) {{
        // Evaluates the early exit condition to terminate the temporal engine if all geometric trajectories have concluded.
        if (total_active_photons <= 0) {{
            break; // Terminates the temporal engine early to conserve hardware cycles.
        }} // END IF: total_active_photons <= 0 to break

        int current = 0; // Integer index tracking the primary active hardware stream for execution overlapping.
        int next = 1; // Integer index tracking the secondary hardware stream preparing the upcoming payload.
        long int active_chunks[2] = {{ 0, 0}}; // 1D array storing the total number of trajectories queued for each operational stream.

        //==========================================
        // PHASE A: PRIME THE PUMP (Stream 0)
        //==========================================
        // Populate the first bridge and launch asynchronous integration on the primary stream.

        active_chunks[current] = NRPYMIN((long int)BUNDLE_CAPACITY, tsm.slot_counts[slot_idx]); // Evaluates the active chunk size bounding the PCIe transfer to avoid VRAM overflow.

        if (active_chunks[current] > 0) {{
            slot_remove_chunk(&tsm, slot_idx, chunk_buffer[current], active_chunks[current]); // Extracts the execution chunk mapping from the Host-side temporal bin.

            // 1. Pack the 9-component tensors using cache-friendly forward sweeps
            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ iterating over the 9 tensor components of the state vectors.
                for (int bridge_i = 0; bridge_i < active_chunks[current]; ++bridge_i) {{  // Loop iterator $bridge_i$ packing the physical state payloads into the Host-side bridge arrays.
                    long int m_idx = chunk_buffer[current][bridge_i]; // Absolute master index $m_{{ idx}}$ mapping the active payload to the global trajectory matrix.
                    f_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f[c_k * num_rays + m_idx]; // Packs the coordinate state vector $f^\mu$ into the transfer bridge.
                    f_p_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p[c_k * num_rays + m_idx]; // Packs the first derivative $\dot{{ f}}^\mu$ into the transfer bridge.
                    f_p_p_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p_p[c_k * num_rays + m_idx]; // Packs the second derivative $\ddot{{ f}}^\mu$ into the transfer bridge.
                }} // END LOOP: for bridge_i over active_chunks[current]
            }} // END LOOP: for c_k over 9

            // 2. Pack the 1D arrays in a separate sequential loop
            for (int bridge_i = 0; bridge_i < active_chunks[current]; ++bridge_i) {{
                long int m_idx = chunk_buffer[current][bridge_i];
                h_bridge[current][bridge_i] = all_photons_host.h[m_idx]; // Packs the current integration step size $h$ into the transfer bridge.
                status_bridge[current][bridge_i] = all_photons_host.status[m_idx]; // Packs the trajectory status enum into the transfer bridge.
                retries_bridge[current][bridge_i] = all_photons_host.rejection_retries[m_idx]; // Packs the error rejection scalar into the transfer bridge.
                integration_param_bridge[current][bridge_i] = all_photons_host.integration_param[m_idx]; // Packs the integration parameter $\lambda$ into the transfer bridge.
                on_pos_non_terminal_plane_prev_bridge[current][bridge_i] = all_photons_host.on_positive_side_of_non_terminal_plane_prev[m_idx]; // Packs the nonterminal plane boundary flag into the transfer bridge.
                on_pos_terminal_plane_prev_bridge[current][bridge_i] = all_photons_host.on_positive_side_of_terminal_plane_prev[m_idx]; // Packs the terminal-plane boundary flag into the transfer bridge.
                integration_param_p_bridge[current][bridge_i] = all_photons_host.integration_param_p[m_idx]; // Packs the historical integration parameter $\lambda_{{ n-1}}$ into the transfer bridge.
                integration_param_p_p_bridge[current][bridge_i] = all_photons_host.integration_param_p_p[m_idx]; // Packs the historical integration parameter $\lambda_{{ n-2}}$ into the transfer bridge.
                non_terminal_plane_event_found_bridge[current][bridge_i] = all_photons_host.non_terminal_plane_event_found[m_idx]; // Packs the nonterminal plane intersection lock into the transfer bridge.
                terminal_plane_event_found_bridge[current][bridge_i] = all_photons_host.terminal_plane_event_found[m_idx]; // Packs the terminal-plane intersection lock into the transfer bridge.
            }} // END LOOP: for bridge_i over active_chunks[current]

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating Host-to-Device transfer of the 9 state vector components.
                // Host-to-Device transfer: Asynchronously pushes bounded state vectors $f^\mu$ to VRAM strictly on stream [current] to minimize latency.
                {memcpy_async("d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "f_bridge[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
                // Host-to-Device transfer: Asynchronously pushes first derivatives $\dot{{ f}}^\mu$ to VRAM strictly on stream [current] to minimize latency.
                {memcpy_async("d_f_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "f_p_bridge[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
                // Host-to-Device transfer: Asynchronously pushes second derivatives $\ddot{{ f}}^\mu$ to VRAM strictly on stream [current] to minimize latency.
                {memcpy_async("d_f_pre_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "f_p_p_bridge[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            }} // END LOOP: for c_k over 9
            // Host-to-Device transfer: Asynchronously pushes step sizes $h$ to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_h[current]", "h_bridge[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes status enums to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_status[current]", "status_bridge[current]", "sizeof(termination_type_t) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes rejection scalars to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_retries[current]", "retries_bridge[current]", "sizeof(int) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes integration parameters $\lambda$ to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_integration_param[current]", "integration_param_bridge[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes nonterminal-plane boundary flags to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_on_pos_non_terminal_plane_prev[current]", "on_pos_non_terminal_plane_prev_bridge[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes terminal-plane boundary flags to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_on_pos_terminal_plane_prev[current]", "on_pos_terminal_plane_prev_bridge[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes historical integration parameters $\lambda_{{ n-1}}$ to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_integration_param_prev[current]", "integration_param_p_bridge[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes historical integration parameters $\lambda_{{ n-2}}$ to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_integration_param_pre_prev[current]", "integration_param_p_p_bridge[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes nonterminal-plane intersection locks to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_non_terminal_plane_event_found[current]", "non_terminal_plane_event_found_bridge[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes terminal-plane intersection locks to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_terminal_plane_event_found[current]", "terminal_plane_event_found_bridge[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Host-to-Device transfer: Asynchronously pushes chunk indices $m_{{ idx}}$ to VRAM strictly on stream [current] to minimize latency.
            {memcpy_async("d_chunk_buffer[current]", "chunk_buffer[current]", "sizeof(long int) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating Device-to-Device baseline setup of the 9 state vector components.
                // Device-to-Device transfer: Duplicates the initial physical state vector $f^\mu$ to anchor the final RKF45 evaluation.
                {memcpy_async("d_f_start_bundle[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToDevice", "streams[current]")}
                // Device-to-Device transfer: Primes the temporary state vector bundle $f^\mu_{{ temp}}$ for iterative stage accumulation.
                {memcpy_async("d_f_temp_bundle[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToDevice", "streams[current]")}
            }} // END LOOP: for c_k over 9

            for (int stage = 1; stage <= 6; ++stage) {{  // Loop iterator $stage$ executing the 6 discrete stages of the RKF45 Runge-Kutta numerical solver.
                // Kernel Launch: Evaluates the metric tensor $g_{{ \mu\nu}}$ and connection $\Gamma^\alpha_{{ \beta\gamma}}$ asynchronously on the active stream.
                interpolation_kernel_{ spacetime_name}(commondata, d_f_temp_bundle[current], d_metric_bundle[current], d_connection_bundle[current], active_chunks[current]{stream_arg_current});
                // Kernel Launch: Computes the geodesic equation right-hand-side derivatives $\dot{{ f}}^\mu$ asynchronously on the active stream.
                calculate_ode_rhs_kernel(d_f_temp_bundle[current], d_metric_bundle[current], d_connection_bundle[current], d_k_bundle[current], stage, active_chunks[current]{stream_arg_current});
                // Stage 6 still computes its RHS; only intermediate state update is skipped.
                if (stage < 6) {{
                    // Kernel Launch: Accumulates the intermediate RKF45 stage numerical updates asynchronously on the active stream.
                    rkf45_stage_update(d_f_start_bundle[current], d_k_bundle[current], d_h[current], stage, active_chunks[current], d_f_temp_bundle[current]{stream_arg_current});
                }} // END IF: skip stage-6 intermediate state update
            }} // END LOOP: for stage over 6

            // Kernel Launch: Applies RKF45 error control to finalize the step-size $h$ and update the integration baseline.
            rkf45_finalize_and_control(commondata, d_f_bundle[current], d_f_start_bundle[current], d_k_bundle[current], d_h[current], d_status[current], d_integration_param[current], d_retries[current], active_chunks[current]{stream_arg_current});
            // Kernel Launch: Detects geometric events and records intersection coordinate states asynchronously on the active stream.
            event_detection_manager_kernel(commondata, d_f_bundle[current], d_f_prev_bundle[current], d_f_pre_prev_bundle[current], d_integration_param[current], d_integration_param_prev[current], d_integration_param_pre_prev[current], d_results_buffer, d_status[current], d_on_pos_non_terminal_plane_prev[current], d_on_pos_terminal_plane_prev[current], d_non_terminal_plane_event_found[current], d_terminal_plane_event_found[current], d_chunk_buffer[current], active_chunks[current]{stream_arg_current});

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating Device-to-Host transfer of the 9 state vector components.
                // Device-to-Host transfer: Retrieves updated coordinate states $f^\mu$ back to CPU RAM asynchronously on the active stream.
                {memcpy_async("f_bridge[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
                // Device-to-Host transfer: Retrieves updated first derivatives $\dot{{ f}}^\mu$ back to CPU RAM asynchronously on the active stream.
                {memcpy_async("f_p_bridge[current] + c_k * BUNDLE_CAPACITY", "d_f_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
                // Device-to-Host transfer: Retrieves updated second derivatives $\ddot{{ f}}^\mu$ back to CPU RAM asynchronously on the active stream.
                {memcpy_async("f_p_p_bridge[current] + c_k * BUNDLE_CAPACITY", "d_f_pre_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            }} // END LOOP: for c_k over 9
            // Device-to-Host transfer: Retrieves active step sizes $h$ back to CPU RAM asynchronously on the active stream.
            {memcpy_async("h_bridge[current]", "d_h[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves updated status enums back to CPU RAM asynchronously on the active stream.
            {memcpy_async("status_bridge[current]", "d_status[current]", "sizeof(termination_type_t) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves active rejection counts back to CPU RAM asynchronously on the active stream.
            {memcpy_async("retries_bridge[current]", "d_retries[current]", "sizeof(int) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves total integration-parameter progression $\lambda$ back to CPU RAM asynchronously on the active stream.
            {memcpy_async("integration_param_bridge[current]", "d_integration_param[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves updated nonterminal-plane boundary flags back to CPU RAM asynchronously on the active stream.
            {memcpy_async("on_pos_non_terminal_plane_prev_bridge[current]", "d_on_pos_non_terminal_plane_prev[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves updated terminal-plane boundary flags back to CPU RAM asynchronously on the active stream.
            {memcpy_async("on_pos_terminal_plane_prev_bridge[current]", "d_on_pos_terminal_plane_prev[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves historical integration parameter $\lambda_{{ n-1}}$ back to CPU RAM asynchronously on the active stream.
            {memcpy_async("integration_param_p_bridge[current]", "d_integration_param_prev[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves historical integration parameter $\lambda_{{ n-2}}$ back to CPU RAM asynchronously on the active stream.
            {memcpy_async("integration_param_p_p_bridge[current]", "d_integration_param_pre_prev[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves active nonterminal-plane locks back to CPU RAM asynchronously on the active stream.
            {memcpy_async("non_terminal_plane_event_found_bridge[current]", "d_non_terminal_plane_event_found[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Device-to-Host transfer: Retrieves active terminal-plane locks back to CPU RAM asynchronously on the active stream.
            {memcpy_async("terminal_plane_event_found_bridge[current]", "d_terminal_plane_event_found[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
        }} // END IF: active chunks available

        //==========================================
        // PHASE B: THE OVERLAP LOOP
        //==========================================
        // Continuously alternate between streams, packing the next payload while syncing the current.

        while (active_chunks[current] > 0 || tsm.slot_counts[slot_idx] > 0) {{

            active_chunks[next] = NRPYMIN((long int)BUNDLE_CAPACITY, tsm.slot_counts[slot_idx]); // Evaluates the active chunk size bounding the upcoming PCIe transfer execution block.
            if (active_chunks[next] > 0) {{
                slot_remove_chunk(&tsm, slot_idx, chunk_buffer[next], active_chunks[next]); // Extracts the next execution chunk mapping from the Host-side temporal bin.

                // 1. Pack the 9-component tensors using cache-friendly forward sweeps
                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ iterating over the 9 tensor components of the state vectors.
                    for (int bridge_i = 0; bridge_i < active_chunks[next]; ++bridge_i) {{  // Loop iterator $bridge_i$ packing the physical state payloads into the next Host-side bridge array.
                        long int m_idx = chunk_buffer[next][bridge_i]; // Absolute master index $m_{{ idx}}$ mapping the active payload to the global trajectory matrix.
                        f_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f[c_k * num_rays + m_idx]; // Packs the coordinate state vector $f^\mu$ into the transfer bridge.
                        f_p_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p[c_k * num_rays + m_idx]; // Packs the first derivative $\dot{{ f}}^\mu$ into the transfer bridge.
                        f_p_p_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p_p[c_k * num_rays + m_idx]; // Packs the second derivative $\ddot{{ f}}^\mu$ into the transfer bridge.
                    }} // END LOOP: for bridge_i over active_chunks[next]
                }} // END LOOP: for c_k over 9

                // 2. Pack the 1D arrays in a separate sequential loop
                for (int bridge_i = 0; bridge_i < active_chunks[next]; ++bridge_i) {{
                    long int m_idx = chunk_buffer[next][bridge_i];
                    h_bridge[next][bridge_i] = all_photons_host.h[m_idx]; // Packs the current integration step size $h$ into the transfer bridge.
                    status_bridge[next][bridge_i] = all_photons_host.status[m_idx]; // Packs the trajectory status enum into the transfer bridge.
                    retries_bridge[next][bridge_i] = all_photons_host.rejection_retries[m_idx]; // Packs the error rejection scalar into the transfer bridge.
                    integration_param_bridge[next][bridge_i] = all_photons_host.integration_param[m_idx]; // Packs the integration parameter $\lambda$ into the transfer bridge.
                    on_pos_non_terminal_plane_prev_bridge[next][bridge_i] = all_photons_host.on_positive_side_of_non_terminal_plane_prev[m_idx]; // Packs the nonterminal plane boundary flag into the transfer bridge.
                    on_pos_terminal_plane_prev_bridge[next][bridge_i] = all_photons_host.on_positive_side_of_terminal_plane_prev[m_idx]; // Packs the terminal-plane boundary flag into the transfer bridge.
                    integration_param_p_bridge[next][bridge_i] = all_photons_host.integration_param_p[m_idx]; // Packs the historical integration parameter $\lambda_{{ n-1}}$ into the transfer bridge.
                    integration_param_p_p_bridge[next][bridge_i] = all_photons_host.integration_param_p_p[m_idx]; // Packs the historical integration parameter $\lambda_{{ n-2}}$ into the transfer bridge.
                    non_terminal_plane_event_found_bridge[next][bridge_i] = all_photons_host.non_terminal_plane_event_found[m_idx]; // Packs the nonterminal plane intersection lock into the transfer bridge.
                    terminal_plane_event_found_bridge[next][bridge_i] = all_photons_host.terminal_plane_event_found[m_idx]; // Packs the terminal-plane intersection lock into the transfer bridge.
                }} // END LOOP: for bridge_i over active_chunks[next]

                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating Host-to-Device transfer of the 9 state vector components for the upcoming payload.
                    // Host-to-Device transfer: Asynchronously pushes bounded state vectors $f^\mu$ to VRAM strictly on stream [next] to overlap execution.
                    {memcpy_async("d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "f_bridge[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                    // Host-to-Device transfer: Asynchronously pushes first derivatives $\dot{{ f}}^\mu$ to VRAM strictly on stream [next] to overlap execution.
                    {memcpy_async("d_f_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "f_p_bridge[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                    // Host-to-Device transfer: Asynchronously pushes second derivatives $\ddot{{ f}}^\mu$ to VRAM strictly on stream [next] to overlap execution.
                    {memcpy_async("d_f_pre_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "f_p_p_bridge[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                }} // END LOOP: for c_k over 9
                // Host-to-Device transfer: Asynchronously pushes step sizes $h$ to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_h[next]", "h_bridge[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes status enums to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_status[next]", "status_bridge[next]", "sizeof(termination_type_t) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes rejection scalars to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_retries[next]", "retries_bridge[next]", "sizeof(int) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes integration parameters $\lambda$ to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_integration_param[next]", "integration_param_bridge[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes nonterminal-plane boundary flags to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_on_pos_non_terminal_plane_prev[next]", "on_pos_non_terminal_plane_prev_bridge[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes terminal-plane boundary flags to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_on_pos_terminal_plane_prev[next]", "on_pos_terminal_plane_prev_bridge[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes historical integration parameters $\lambda_{{ n-1}}$ to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_integration_param_prev[next]", "integration_param_p_bridge[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes historical integration parameters $\lambda_{{ n-2}}$ to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_integration_param_pre_prev[next]", "integration_param_p_p_bridge[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
            // Host-to-Device transfer: Asynchronously pushes nonterminal-plane intersection locks to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_non_terminal_plane_event_found[next]", "non_terminal_plane_event_found_bridge[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes terminal-plane intersection locks to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_terminal_plane_event_found[next]", "terminal_plane_event_found_bridge[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Host-to-Device transfer: Asynchronously pushes chunk indices $m_{{ idx}}$ to VRAM strictly on stream [next] to overlap execution.
                {memcpy_async("d_chunk_buffer[next]", "chunk_buffer[next]", "sizeof(long int) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}

                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating Device-to-Device baseline setup of the 9 state vector components for the upcoming payload.
                    // Device-to-Device transfer: Duplicates the initial physical state vector $f^\mu$ to anchor the upcoming RKF45 evaluation.
                    {memcpy_async("d_f_start_bundle[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToDevice", "streams[next]")}
                    // Device-to-Device transfer: Primes the temporary state vector bundle $f^\mu_{{ temp}}$ for the upcoming iterative stage accumulation.
                    {memcpy_async("d_f_temp_bundle[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToDevice", "streams[next]")}
                }} // END LOOP: for c_k over 9

                for (int stage = 1; stage <= 6; ++stage) {{  // Loop iterator $stage$ executing the 6 discrete stages of the upcoming RKF45 Runge-Kutta numerical solver.
                    // Kernel Launch: Evaluates the metric tensor $g_{{ \mu\nu}}$ and connection $\Gamma^\alpha_{{ \beta\gamma}}$ asynchronously on the alternate stream.
                    interpolation_kernel_{spacetime_name}(commondata, d_f_temp_bundle[next], d_metric_bundle[next], d_connection_bundle[next],active_chunks[next]{stream_arg_next});
                    // Kernel Launch: Computes the geodesic equation right-hand-side derivatives $\dot{{ f}}^\mu$ asynchronously on the alternate stream.
                    calculate_ode_rhs_kernel(d_f_temp_bundle[next], d_metric_bundle[next], d_connection_bundle[next], d_k_bundle[next], stage, active_chunks[next]{stream_arg_next});
                    // Stage 6 still computes its RHS; only intermediate state update is skipped.
                    if (stage < 6) {{
                        // Kernel Launch: Accumulates the intermediate RKF45 stage numerical updates asynchronously on the alternate stream.
                        rkf45_stage_update(d_f_start_bundle[next], d_k_bundle[next], d_h[next], stage, active_chunks[next], d_f_temp_bundle[next]{stream_arg_next});
                    }} // END IF: skip stage-6 intermediate state update
                }} // END LOOP: for stage over 6

                // Kernel Launch: Applies RKF45 error control to finalize the step-size $h$ and update the upcoming integration baseline.
                rkf45_finalize_and_control(commondata, d_f_bundle[next], d_f_start_bundle[next], d_k_bundle[next], d_h[next], d_status[next], d_integration_param[next], d_retries[next], active_chunks[next]{stream_arg_next});
                // Kernel Launch: Detects geometric events and records intersection coordinate states asynchronously on the alternate stream.
                event_detection_manager_kernel(commondata, d_f_bundle[next], d_f_prev_bundle[next], d_f_pre_prev_bundle[next], d_integration_param[next], d_integration_param_prev[next], d_integration_param_pre_prev[next], d_results_buffer, d_status[next], d_on_pos_non_terminal_plane_prev[next], d_on_pos_terminal_plane_prev[next], d_non_terminal_plane_event_found[next], d_terminal_plane_event_found[next], d_chunk_buffer[next], active_chunks[next]{stream_arg_next});
                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating Device-to-Host transfer of the 9 upcoming state vector components.
                    // Device-to-Host transfer: Retrieves updated coordinate states $f^\mu$ back to CPU RAM asynchronously on the alternate stream.
                    {memcpy_async("f_bridge[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                    // Device-to-Host transfer: Retrieves updated first derivatives $\dot{{ f}}^\mu$ back to CPU RAM asynchronously on the alternate stream.
                    {memcpy_async("f_p_bridge[next] + c_k * BUNDLE_CAPACITY", "d_f_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                    // Device-to-Host transfer: Retrieves updated second derivatives $\ddot{{ f}}^\mu$ back to CPU RAM asynchronously on the alternate stream.
                    {memcpy_async("f_p_p_bridge[next] + c_k * BUNDLE_CAPACITY", "d_f_pre_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                }} // END LOOP: for c_k over 9
                // Device-to-Host transfer: Retrieves upcoming active step sizes $h$ back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("h_bridge[next]", "d_h[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming updated status enums back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("status_bridge[next]", "d_status[next]", "sizeof(termination_type_t) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming active rejection counts back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("retries_bridge[next]", "d_retries[next]", "sizeof(int) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming total integration-parameter progression $\lambda$ back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("integration_param_bridge[next]", "d_integration_param[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming updated nonterminal-plane boundary flags back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("on_pos_non_terminal_plane_prev_bridge[next]", "d_on_pos_non_terminal_plane_prev[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming updated terminal-plane boundary flags back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("on_pos_terminal_plane_prev_bridge[next]", "d_on_pos_terminal_plane_prev[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming historical integration parameter $\lambda_{{ n-1}}$ back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("integration_param_p_bridge[next]", "d_integration_param_prev[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming historical integration parameter $\lambda_{{ n-2}}$ back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("integration_param_p_p_bridge[next]", "d_integration_param_pre_prev[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming active nonterminal-plane locks back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("non_terminal_plane_event_found_bridge[next]", "d_non_terminal_plane_event_found[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Device-to-Host transfer: Retrieves upcoming active terminal-plane locks back to CPU RAM asynchronously on the alternate stream.
                {memcpy_async("terminal_plane_event_found_bridge[next]", "d_terminal_plane_event_found[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
            }} // END IF: next chunks available

            if (active_chunks[current] > 0) {{
                // Device synchronization barrier strictly enforcing completion of the current stream before payload unpacking.
                {stream_sync("streams[current]")}

                // 1. Unpack 9-component tensors sequentially
                for (int fin_k = 0; fin_k < 9; ++fin_k) {{  // Loop index $fin_k$ retrieving the 9 tensor components back into the global Host matrix.
                    for (int fin_i = 0; fin_i < active_chunks[current]; ++fin_i) {{  // Loop iterator $fin_i$ unpacking the finalized physical data back to the global Host matrix.
                        long int m_idx = chunk_buffer[current][fin_i]; // Absolute master index $m_{{ idx}}$ retrieving the specific photon index from the execution chunk.
                        all_photons_host.f[fin_k * num_rays + m_idx] = f_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i]; // Unpacks the synchronized state vector $f^\mu$ into the global Host matrix.
                        all_photons_host.f_p[fin_k * num_rays + m_idx] = f_p_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i]; // Unpacks the synchronized first derivative $\dot{{ f}}^\mu$ into the global Host matrix.
                        all_photons_host.f_p_p[fin_k * num_rays + m_idx] = f_p_p_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i]; // Unpacks the synchronized second derivative $\ddot{{ f}}^\mu$ into the global Host matrix.
                    }} // END LOOP: for fin_i over active_chunks[current]
                }} // END LOOP: for fin_k over 9

                // 2. Unpack 1D arrays sequentially
                for (int fin_i = 0; fin_i < active_chunks[current]; ++fin_i) {{
                    long int m_idx = chunk_buffer[current][fin_i];
                    all_photons_host.h[m_idx] = h_bridge[current][fin_i]; // Unpacks the synchronized step size $h$ into the global Host matrix.
                    all_photons_host.status[m_idx] = status_bridge[current][fin_i]; // Unpacks the synchronized trajectory status into the global Host matrix.
                    all_photons_host.rejection_retries[m_idx] = retries_bridge[current][fin_i]; // Unpacks the synchronized rejection count into the global Host matrix.
                    all_photons_host.integration_param[m_idx] = integration_param_bridge[current][fin_i]; // Unpacks the synchronized integration parameter $\lambda$ into the global Host matrix.
                    all_photons_host.on_positive_side_of_non_terminal_plane_prev[m_idx] = on_pos_non_terminal_plane_prev_bridge[current][fin_i]; // Unpacks the synchronized nonterminal-plane boundary flag into the global Host matrix.
                    all_photons_host.on_positive_side_of_terminal_plane_prev[m_idx] = on_pos_terminal_plane_prev_bridge[current][fin_i]; // Unpacks the synchronized terminal-plane boundary flag into the global Host matrix.
                    all_photons_host.integration_param_p[m_idx] = integration_param_p_bridge[current][fin_i]; // Unpacks the synchronized historical integration parameter $\lambda_{{ n-1}}$ into the global Host matrix.
                    all_photons_host.integration_param_p_p[m_idx] = integration_param_p_p_bridge[current][fin_i]; // Unpacks the synchronized historical integration parameter $\lambda_{{ n-2}}$ into the global Host matrix.
                    all_photons_host.non_terminal_plane_event_found[m_idx] = non_terminal_plane_event_found_bridge[current][fin_i]; // Unpacks the synchronized nonterminal-plane lock into the global Host matrix.
                    all_photons_host.terminal_plane_event_found[m_idx] = terminal_plane_event_found_bridge[current][fin_i]; // Unpacks the synchronized terminal-plane lock into the global Host matrix.
                }} // END LOOP: for fin_i over active_chunks[current]

                // 3. TimeSlotManager State Update (Cache-hot, strictly sequential)
                for (int fin_i = 0; fin_i < active_chunks[current]; ++fin_i) {{
                    long int m_idx = chunk_buffer[current][fin_i];
                    if (status_bridge[current][fin_i] == ACTIVE) {{  // Evaluates the continuation logic if the trajectory remains within safe physical bounds.
                        int next_s_idx = slot_get_index(&tsm, all_photons_host.f[m_idx]); // Evaluates the updated temporal coordinate $t$ to determine the next operational bin.
                        if (next_s_idx != -1) {{  // Confirms the physical state has not exceeded the maximum simulation time bounds.
                            slot_add_photon(&tsm, next_s_idx, m_idx);  // Re-queues the updated physical state vector back into the host orchestrator.
                        }} else {{
                            all_photons_host.status[m_idx] = STOP_CONDITION_T_MAX_EXCEEDED; // Flags the physical state as stopped after excessive propagation time.
                            total_active_photons--; // Decrements the global counter as the physical trajectory has reached a terminal state.
                        }} // END ELSE: state flagged failed
                    }} // END IF: trajectory remains active
                    else if (status_bridge[current][fin_i] == REJECTED) {{   // Evaluates the retry logic if the numerical step exceeded the requested tolerances.
                        slot_add_photon(&tsm, slot_idx, m_idx); // Re-adds to the current bin to attempt integration with an adapted step-size scalar $h$.
                    }} else {{
                        total_active_photons--; // Decrements the global counter as the physical trajectory has reached a terminal state.
                    }} // END ELSE: trajectory reached terminal state
                }} // END LOOP: for fin_i over active_chunks[current]
            }} // END IF: active chunks completed

            //==========================================
            // PROGRESS DASHBOARD
            //==========================================
            // Evaluate computational throughput and temporal progress to update the terminal dashboard.

            // Captures the terminal hardware clock state for the current chunk.
            clock_gettime(CLOCK_MONOTONIC, &batch_end_time);

            // Evaluates the absolute wall-clock duration of the integration chunk in seconds.
            double elapsed_sec = (batch_end_time.tv_sec - batch_start_time.tv_sec) + (batch_end_time.tv_nsec - batch_start_time.tv_nsec) / 1e9;

            // Evaluates the raw integration steps per second to monitor pipeline throughput.
            double steps_per_sec = (elapsed_sec > 0.0) ? ((double)active_chunks[current] / elapsed_sec) : 0.0;

            // Evaluates the global completion ratio bounded between $0.0$ and $1.0$.
            double percent_done = 100.0 * (1.0 - ((double)total_active_photons / (double)num_rays));

            // Evaluates the physical coordinate time $t$ for the active temporal bin.
            double current_t = commondata->slot_manager_t_min + slot_idx * commondata->slot_manager_delta_t;

            // Defines the total character width of the dynamic loading bar visualization.
            int bar_width = 20;
            // Computes the integer index demarcating the active boundary within the loading bar.
            int pos = (int)(bar_width * percent_done / 100.0);
            // Character array storing the formatted loading bar string.
            char bar[21];

            // Loop iterator $bar_i$ constructing the ASCII loading bar visualizer.
            for (int bar_i = 0; bar_i < bar_width; ++bar_i) {{
                if (bar_i < pos) bar[bar_i] = '='; // Appends the completed progression character.
                else if (bar_i == pos) bar[bar_i] = '>'; // Appends the active vanguard character.
                else bar[bar_i] = ' '; // Appends the uncompleted progression character.
            }} // END LOOP: for bar_i over bar_width
            bar[bar_width] = '\0'; // Terminates the loading bar character array to prevent buffer overruns.

            // Accumulator tracking the total number of adaptive step size $h$ rejections in the active chunk.
            long int batch_rejections = 0;

            // Loop iterator $sum_i$ scanning the finalized physical state bridge for error tolerance failures.
            for (int sum_i = 0; sum_i < active_chunks[current]; ++sum_i) {{
                batch_rejections += retries_bridge[current][sum_i]; // Accumulates the localized step rejection tally.
            }} // END LOOP: for sum_i over active_chunks[current]

            // Evaluates the relative frequency of adaptive step size $h$ rejections.
            double reject_percent = (active_chunks[current] > 0) ? (100.0 * (double)batch_rejections / (double)active_chunks[current]) : 0.0;

            // Restores the terminal cursor position vertically to overwrite the previous dashboard iteration.
            printf("\033[7A");
            printf("--------------------------------------------------\n"); // Prints the upper border of the diagnostic dashboard.
            printf(" Progress:   [%s] %5.1f%% \033[K\n", bar, percent_done); // Prints the global completion loading bar and percentage.
            printf(" Active:     %ld / %ld \033[K\n", total_active_photons, num_rays); // Prints the remaining active photon trajectories $x^\mu$.
            printf(" Slot Time:  Slot %d (t = %.1f) \033[K\n", slot_idx, current_t); // Prints the current physical temporal bin coordinate $t$.
            printf(" Speed:      %.2e integration steps/s \033[K\n", steps_per_sec); // Prints the pipeline execution throughput.
            printf(" Rejects:    %ld (%.1f%%) \033[K\n", batch_rejections, reject_percent); // Prints the adaptive step size $h$ rejection frequency.
            printf("--------------------------------------------------\n"); // Prints the lower border of the diagnostic dashboard.
            fflush(stdout); // Flushes the standard output buffer to ensure instantaneous terminal rendering.

            // Resets the hardware clock state for the upcoming integration chunk.
            clock_gettime(CLOCK_MONOTONIC, &batch_start_time);
            // --------------------------

            active_chunks[current] = 0; // Clears the execution queue tracker to indicate the active chunk has been fully processed.
            int temp = current; // Temporary integer scalar storing the primary stream index for logical pointer swapping.
            current = next; // Shifts the primary execution tracker to the alternate stream index.
            next = temp; // Assigns the cleared stream index back to the upcoming payload queue.
        }} // END WHILE: alternate streams process temporal slot

        //==========================================
        // PHASE C: THE TIME BARRIER
        //==========================================
        // Enforce a rigid hardware sync before advancing the physical time clock.
        BHAH_DEVICE_SYNC(); // Synchronize global state convergence before advancing the central time engine.

     }} // END LOOP: for slot_idx down to 0

    //==========================================
    // 4. Conserved Values & CLEANUP & FINALIZATION
    //==========================================
    // Process terminal photon trajectories and extract final geometric intersections.

        // Device-to-Host transfer: Extracts validated device-native blueprints $b_i$ containing geometric plane intersections.
        {results_memcpy}

        // Kernel Launch: Processes escaped photons intersecting the celestial sphere $r > r_{{escape}}$.
        {calc_blueprint}

        //==========================================
        // TERMINAL NORMALIZATION DIAGNOSTIC
        //==========================================
        // Evaluate terminal normalization constraint.
        if (commondata->perform_conservation_check) {{
            normalization_constraint_t *norm_diag_bridge; // Host pointer for retrieving diagnostic normalization records.
            {malloc_pinned}(norm_diag_bridge, sizeof(normalization_constraint_t) * BUNDLE_CAPACITY); // Pinned memory maximizes PCIe transfer bandwidth for diagnostic arrays.

            double max_err_norm = 0.0; // Scalar tracking the maximum absolute drift from the expected normalization constraint invariant $C$.
            long int worst_ray_norm = -1; // Absolute index identifying the trajectory $x^\mu$ with the highest constraint violation.

            long int norm_num_batches = (num_rays + BUNDLE_CAPACITY - 1) / BUNDLE_CAPACITY; // Integer calculation defining the total sequential blocks required to process all photon trajectories.

            for (long int norm_batch = 0; norm_batch < norm_num_batches; ++norm_batch) {{ // Loop iterator $norm_batch$ for evaluating the terminal normalization constraint across sequential chunks.
                long int start_idx = norm_batch * BUNDLE_CAPACITY; // Absolute starting index mapped to the master SoA for the current normalization batch.
                long int chunk_size = NRPYMIN((long int)BUNDLE_CAPACITY, num_rays - start_idx); // Dynamically sized operational boundary ensuring the active chunk does not exceed the total trajectory count.

                for (int norm_i = 0; norm_i < chunk_size; ++norm_i) {{ // Loop index $norm_i$ iterating over the specific normalization batch elements to pack the bridge.
                    long int master_idx = start_idx + norm_i; // Computes the absolute master index $m_{{idx}}$ tracking the specific photon within the global matrix.
                    for (int norm_k = 0; norm_k < 9; ++norm_k) {{ // Loop index $norm_k$ iterating over the 9 tensor components of the state vector $f^\mu$.
                        f_bridge[0][norm_k * BUNDLE_CAPACITY + norm_i] = all_photons_host.f[norm_k * num_rays + master_idx]; // Assigns the terminal tensor state vector component to the primary bridge array.
                    }} // END LOOP: for norm_k over 9
                }} // END LOOP: for norm_i over chunk_size

                for (int c_k = 0; c_k < 9; ++c_k) {{ // Loop index $c_k$ orchestrating the asynchronous memory transfer of the 9 state vector $f^\mu$ components.
                    {memcpy_async("d_f_bundle[0] + c_k * BUNDLE_CAPACITY", "f_bridge[0] + c_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size", "cudaMemcpyHostToDevice", "streams[0]")}
                }} // END LOOP: for c_k over 9

                // Kernel Launch: Calculates the symmetric metric tensor $g_{{\mu\nu}}$ strictly on the primary operational pipeline to supply the normalization constraint evaluator.
                interpolation_kernel_{spacetime_name}(commondata, d_f_bundle[0], d_metric_bundle[0], NULL, chunk_size{stream_arg});

                // Kernel Launch: Computes the normalization constraint $C = g_{{\mu\nu}} p^\mu p^\nu$ natively on the active hardware pipeline synchronously with the memory stream.
                normalization_constraint_photon(d_f_bundle[0], d_metric_bundle[0], d_norm_bundle[0], chunk_size{stream_arg});

                // Device-to-Host transfer: Retrieves active diagnostic normalization records back to the Host RAM asynchronously on the primary stream.
                {memcpy_async("norm_diag_bridge", "d_norm_bundle[0]", "sizeof(normalization_constraint_t) * chunk_size", "cudaMemcpyDeviceToHost", "streams[0]")}

                // Device synchronization barrier strictly enforcing completion of the primary stream before diagnostic payload unpacking.
                {stream_sync('streams[0]')}

                for (int norm_i = 0; norm_i < chunk_size; ++norm_i) {{ // Loop iterator $norm_i$ scanning each trajectory within the current diagnostic memory chunk.
                    double current_norm_err = fabs(norm_diag_bridge[norm_i].C); // Evaluates the absolute numerical drift for the physical scalar invariant $C$.
                    if (current_norm_err > max_err_norm) {{
                        max_err_norm = current_norm_err; // Updates the maximum tracked absolute error for the geometric normalization constraint.
                        worst_ray_norm = start_idx + norm_i; // Updates the absolute master index $m_{{idx}}$ associated with the maximum geometric constraint violation.
                    }} // END IF: current_norm_err > max_err_norm
                }} // END LOOP: for norm_i over chunk_size
            }} // END LOOP: for norm_batch over norm_num_batches

            printf("\n=================================================\n");
            printf(" NORMALIZATION DIAGNOSTIC REPORT\n");
            printf("=================================================\n");
            printf("  Max Absolute Error (Normalization Constraint |g_mu_nu p^mu p^nu|): %e (Ray %ld)\n", max_err_norm, worst_ray_norm);
            {free_pinned}(norm_diag_bridge); // Memory Free: Purges the diagnostic bridge utilized for the geometric normalization constraint checks.
        }} // END IF: commondata->perform_conservation_check to evaluate terminal normalization

        // Loop iterator $s$ purging the double-buffered arrays across both hardware streams.
        for (int s = 0; s < 2; ++s) {{
            // Host Memory Free: Purges bridge components supporting scatter logic mapped to PCIe DMA transfers.
            {free_pinned}(chunk_buffer[s]); // Purges the execution chunk mapping bridge.
            {free_pinned}(f_bridge[s]); // Purges the state vector $f^\mu$ bridge.
            {free_pinned}(f_p_bridge[s]); // Purges the first derivative $\dot{{f}}^\mu$ bridge.
            {free_pinned}(f_p_p_bridge[s]); // Purges the second derivative $\ddot{{f}}^\mu$ bridge.
            {free_pinned}(integration_param_bridge[s]); // Purges the integration parameter $\lambda$ bridge.
            {free_pinned}(h_bridge[s]); // Purges the integration step size $h$ bridge.
            {free_pinned}(status_bridge[s]); // Purges the trajectory status bridge.
            {free_pinned}(retries_bridge[s]); // Purges the error rejection scalar bridge.
            {free_pinned}(on_pos_non_terminal_plane_prev_bridge[s]); // Purges the nonterminal plane boundary flag bridge.
            {free_pinned}(on_pos_terminal_plane_prev_bridge[s]); // Purges the terminal-plane boundary flag bridge.
            {free_pinned}(integration_param_p_bridge[s]); // Purges the historical integration parameter $\lambda_{{n-1}}$ bridge.
            {free_pinned}(integration_param_p_p_bridge[s]); // Purges the historical integration parameter $\lambda_{{n-2}}$ bridge.
            {free_pinned}(non_terminal_plane_event_found_bridge[s]); // Purges the nonterminal plane intersection lock bridge.
            {free_pinned}(terminal_plane_event_found_bridge[s]); // Purges the terminal-plane intersection lock bridge.

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
            {free_device}(d_integration_param[s]); // Purges the total integration parameter progress $\lambda$ scratchpad.
            {free_device}(d_status[s]); // Purges the current trajectory status limit scratchpad.
            {free_device}(d_retries[s]); // Purges the sequential error rejection scratchpad.
            {free_device}(d_on_pos_non_terminal_plane_prev[s]); // Purges the previous nonterminal plane boundary side scratchpad.
            {free_device}(d_on_pos_terminal_plane_prev[s]); // Purges the previous terminal-plane boundary side scratchpad.
            {free_device}(d_integration_param_prev[s]); // Purges the historical integration parameter $\lambda_{{n-1}}$ scratchpad.
            {free_device}(d_integration_param_pre_prev[s]); // Purges the historical integration parameter $\lambda_{{n-2}}$ scratchpad.
            {free_device}(d_non_terminal_plane_event_found[s]); // Purges the nonterminal-plane intersection coordinate guard scratchpad.
            {free_device}(d_terminal_plane_event_found[s]); // Purges the terminal-plane intersection coordinate guard scratchpad.
            {free_device}(d_chunk_buffer[s]); // Purges the absolute master indices $m_{{idx}}$ mapping scratchpad.
            {free_device}(d_norm_bundle[s]); // Purges the diagnostic outputs scratchpad.

            {stream_destroy}
        }} // END LOOP: for s over 2

        //==========================================
        // CPU CONSERVATION DRIFT EVALUATION
        //==========================================
        // Evaluate relative numerical drift on the CPU.
        if (commondata->perform_conservation_check) {{
            // Kernel Launch: Calculates terminal conserved quantities natively on the executing architecture.
            calculate_conserved_quantities_universal_{spacetime_name}_photon(commondata, &all_photons_host, num_rays, final_cq_host);

            printf("\n=================================================\n");
            printf(" CONSERVED QUANTITIES DIAGNOSTIC REPORT\n");
            printf("=================================================\n");

            double max_err_E = 0.0; // Scalar variable tracking the maximum recorded relative drift for energy $E$.
            double max_err_Lz = 0.0; // Scalar variable tracking the maximum recorded relative drift for angular momentum $L_z$.
            double max_err_Q = 0.0; // Scalar variable tracking the maximum recorded relative drift for Carter constant $Q$.

            long int worst_ray_E = -1; // Absolute master index $m_{{idx}}$ identifying the trajectory responsible for the maximum relative numerical drift in energy $E$.
            long int worst_ray_Lz = -1; // Absolute master index $m_{{idx}}$ identifying the trajectory responsible for the maximum relative numerical drift in angular momentum $L_z$.
            long int worst_ray_Q = -1; // Absolute master index $m_{{idx}}$ identifying the trajectory responsible for the maximum relative numerical drift in Carter constant $Q$.

            double max_abs_err_E = 0.0; // Scalar variable tracking the maximum recorded absolute drift for energy $E$.
            double max_abs_err_Lz = 0.0; // Scalar variable tracking the maximum recorded absolute drift for angular momentum $L_z$.
            double max_abs_err_Q = 0.0; // Scalar variable tracking the maximum recorded absolute drift for Carter constant $Q$.

            long int worst_ray_abs_E = -1; // Absolute master index $m_{{idx}}$ identifying the trajectory responsible for the maximum absolute numerical drift in energy $E$.
            long int worst_ray_abs_Lz = -1; // Absolute master index $m_{{idx}}$ identifying the trajectory responsible for the maximum absolute numerical drift in angular momentum $L_z$.
            long int worst_ray_abs_Q = -1; // Absolute master index $m_{{idx}}$ identifying the trajectory responsible for the maximum absolute numerical drift in Carter constant $Q$.

            // Loop iterator $i$ spanning the global dataset to calculate errors natively on the CPU.
            for (long int i = 0; i < num_rays; i++) {{

                double err_E = fabs((final_cq_host[i].E - initial_cq_host[i].E) / (initial_cq_host[i].E + 1e-15)); // Evaluates the relative numerical drift for energy $E$.
                double err_Lz = fabs((final_cq_host[i].Lz - initial_cq_host[i].Lz) / (initial_cq_host[i].Lz + 1e-15)); // Evaluates the relative numerical drift for angular momentum $L_z$.
                double err_Q = fabs((final_cq_host[i].Q - initial_cq_host[i].Q) / (initial_cq_host[i].Q + 1e-15)); // Evaluates the relative numerical drift for Carter constant $Q$.

                double abs_err_E = fabs(final_cq_host[i].E - initial_cq_host[i].E); // Evaluates the absolute numerical drift for energy $E$.
                double abs_err_Lz = fabs(final_cq_host[i].Lz - initial_cq_host[i].Lz); // Evaluates the absolute numerical drift for angular momentum $L_z$.
                double abs_err_Q = fabs(final_cq_host[i].Q - initial_cq_host[i].Q); // Evaluates the absolute numerical drift for Carter constant $Q$.

                if (err_E > max_err_E) {{
                    max_err_E = err_E; // Updates the maximum tracked relative error for energy $E$.
                    worst_ray_E = i; // Updates the absolute master index $m_{{idx}}$ for the maximum relative energy drift.
                }} // END IF: err_E > max_err_E

                if (err_Lz > max_err_Lz) {{
                    max_err_Lz = err_Lz; // Updates the maximum tracked relative error for angular momentum $L_z$.
                    worst_ray_Lz = i; // Updates the absolute master index $m_{{idx}}$ for the maximum relative angular momentum drift.
                }} // END IF: err_Lz > max_err_Lz

                if (err_Q > max_err_Q) {{
                    max_err_Q = err_Q; // Updates the maximum tracked relative error for Carter constant $Q$.
                    worst_ray_Q = i; // Updates the absolute master index $m_{{idx}}$ for the maximum relative Carter constant drift.
                }} // END IF: err_Q > max_err_Q

                if (abs_err_E > max_abs_err_E) {{
                    max_abs_err_E = abs_err_E; // Updates the maximum tracked absolute error for energy $E$.
                    worst_ray_abs_E = i; // Updates the absolute master index $m_{{idx}}$ for the maximum absolute energy drift.
                }} // END IF: abs_err_E > max_abs_err_E

                if (abs_err_Lz > max_abs_err_Lz) {{
                    max_abs_err_Lz = abs_err_Lz; // Updates the maximum tracked absolute error for angular momentum $L_z$.
                    worst_ray_abs_Lz = i; // Updates the absolute master index $m_{{idx}}$ for the maximum absolute angular momentum drift.
                }} // END IF: abs_err_Lz > max_abs_err_Lz

                if (abs_err_Q > max_abs_err_Q) {{
                    max_abs_err_Q = abs_err_Q; // Updates the maximum tracked absolute error for Carter constant $Q$.
                    worst_ray_abs_Q = i; // Updates the absolute master index $m_{{idx}}$ for the maximum absolute Carter constant drift.
                }} // END IF: abs_err_Q > max_abs_err_Q
            }} // END LOOP: for i over num_rays

            printf("  Max Relative Error (Energy E): %e (Ray %ld)\n", max_err_E, worst_ray_E); // Output block printing the maximum relative error for energy $E$.
            printf("  Max Absolute Error (Energy E): %e (Ray %ld)\n\n", max_abs_err_E, worst_ray_abs_E); // Output block printing the maximum absolute error for energy $E$.

            printf("  Max Relative Error (Momentum Lz): %e (Ray %ld)\n", max_err_Lz, worst_ray_Lz); // Output block printing the maximum relative error for angular momentum $L_z$.
            printf("  Max Absolute Error (Momentum Lz): %e (Ray %ld)\n\n", max_abs_err_Lz, worst_ray_abs_Lz); // Output block printing the maximum absolute error for angular momentum $L_z$.

            printf("  Max Relative Error (Carter Q): %e (Ray %ld)\n", max_err_Q, worst_ray_Q); // Output block printing the maximum relative error for Carter constant $Q$.
            printf("  Max Absolute Error (Carter Q): %e (Ray %ld)\n", max_abs_err_Q, worst_ray_abs_Q); // Output block printing the maximum absolute error for Carter constant $Q$.

            printf("=================================================\n\n"); // Output block printing the terminal footer for the diagnostic sequence.

            // Host Memory Free: Purges pinned diagnostic buffers.
            {free_pinned}(initial_cq_host); // Purges pinned initial diagnostic data buffer.
            {free_pinned}(final_cq_host); // Purges pinned final diagnostic data buffer.
        }} // END IF: commondata->perform_conservation_check to evaluate numerical drift

        // Host Memory Free: Purges the primary Host array states $f^\mu$ and integration parameters $\lambda$.
        {free_pinned}(all_photons_host.f); // Purges the primary Host array state $f^\mu$.
        {free_pinned}(all_photons_host.f_p); // Purges the primary Host array first derivative $\dot{{f}}^\mu$.
        {free_pinned}(all_photons_host.f_p_p); // Purges the primary Host array second derivative $\ddot{{f}}^\mu$.
        {free_pinned}(all_photons_host.integration_param); // Purges the primary Host array integration parameter $\lambda$.
        {free_pinned}(all_photons_host.h); // Purges the primary Host array integration step size $h$.
        {free_pinned}(all_photons_host.status); // Purges the primary Host array trajectory status enum.
        {free_pinned}(all_photons_host.rejection_retries); // Purges the primary Host array error rejection scalar.
        {free_pinned}(all_photons_host.on_positive_side_of_non_terminal_plane_prev); // Purges the primary Host array nonterminal plane boundary flag.
        {free_pinned}(all_photons_host.on_positive_side_of_terminal_plane_prev); // Purges the primary Host array terminal-plane boundary flag.
        {free_pinned}(all_photons_host.integration_param_p); // Purges the primary Host array historical integration parameter $\lambda_{{n-1}}$.
        {free_pinned}(all_photons_host.integration_param_p_p); // Purges the primary Host array historical integration parameter $\lambda_{{n-2}}$.
        {free_pinned}(all_photons_host.non_terminal_plane_event_found); // Purges the primary Host array nonterminal plane intersection lock.
        {free_pinned}(all_photons_host.terminal_plane_event_found); // Purges the primary Host array terminal-plane intersection lock.

        // Device Memory Free: Purges the final single-pointer intersection blueprint buffer $b_i$.
        {free_device}(d_results_buffer); // Purges the intersection blueprint buffer $b_i$.

        // Memory Free: Purges the temporal sorting struct mapping the Host-side execution grid.
        slot_manager_free(&tsm); // Purges the central time slot orchestrator.
    """

    cfc.register_CFunction(
        includes=includes,
        prefunc=init_prefunc,
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
