# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_numerical.py
r"""
Generates the CUDA/OpenMP driver for numerical photon integration.

This module generates the C driver that initializes, evolves, and records photon
trajectories $x^\mu$ in numerical spacetimes. Its Split-Pipeline separates the
Runge-Kutta-Fehlberg 4(5) integration stages. Metric and Christoffel-symbol arrays
and photon state vectors use flattened Structure-of-Arrays layouts. Two work-array sets let CUDA overlap
data copies and integration; OpenMP uses the same indexing with host-memory copies.
Each ray chunk computes the metric and connection and advances the photon states.
Diagnostic probes verify
that initialization populated coordinates and zeroed the temporal momentum and
distance-traveled components. The initial momentum solve enforces the photon
normalization constraint. The driver records initial conserved quantities before
RKF45 updates the state vectors and evaluates terminal normalization afterward.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.time_slot_manager_helpers import (
    time_slot_manager_helpers,
)


def batch_integrator_numerical(spacetime_name: str) -> None:
    r"""
    Construct the CUDA/OpenMP driver for batched numerical photon integration.

    :param spacetime_name: The identifier for the spacetime metric (e.g., 'KerrSchild').
    """
    if "time_slot_manager" not in par.glb_extras_dict.get("BHaH_defines", {}):
        time_slot_manager_helpers()

    # Core physics and numerical simulation parameters for the global spacetime struct.
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "t_integration_max",
            "r_escape",
            "p_t_max",
            "numerical_initial_h",
        ],
        [10000.0, 150.0, 1e3, 0.1],
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

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    if parallelization == "cuda":
        includes.extend(
            ["cuda_runtime.h", "cuda_intrinsics.h", "BHaH_global_device_defines.h"]
        )

    desc = r""" Host-side driver for the batched Split-Pipeline relativistic ray tracing loop.

    This function acts as the primary loop for evaluating photon geodesics $x^\mu$.
    It utilizes a TimeSlotManager to bin active rays by their physical coordinate time $t$.
    The Split-Pipeline implementation stores $g_{\mu\nu}$ and $\Gamma^\alpha_{\beta\gamma}$ in scratch arrays.
    The scratch arrays hold tensor values for one bounded ray chunk.

    @param commondata Struct containing global spacetime and numerical tolerances.
    @param num_rays Total number of photon trajectories to simulate.
    @param results_buffer Caller-provided output array for final physical intersections."""

    cfunc_type = "void"

    name = "batch_integrator_numerical"

    params = "const commondata_struct *restrict commondata, long int num_rays, blueprint_data_t *restrict results_buffer"

    include_CodeParameters_h = True

    # --- DYNAMIC MACRO GENERATION ---
    malloc_pinned = "BHAH_MALLOC_PINNED" if parallelization == "cuda" else "BHAH_MALLOC"
    malloc_device = "BHAH_MALLOC_DEVICE" if parallelization == "cuda" else "BHAH_MALLOC"

    if parallelization == "cuda":
        stream_setup_str = """
        // CUDA streams for asynchronous kernel execution and memory transfers.
        cudaStream_t streams[2];
        // Initializes the primary CUDA stream mapped to the first double-buffer context.
        cudaStreamCreate(&streams[0]);
        // Initializes the secondary CUDA stream mapped to the second double-buffer context.
        cudaStreamCreate(&streams[1]);

        // Copy global spacetime and numerical parameters to CUDA constant memory.
        cudaMemcpyToSymbol(d_commondata, commondata, sizeof(commondata_struct));"""
        pin_comment = "Allocate pinned host memory for"
        dev_comment = "Allocate CUDA device memory for"
        bridge_alloc_comment = "Allocate pinned host arrays for CUDA ray-chunk copies."
        scratch_alloc_comment = "Allocate one-dimensional CUDA work arrays."
    else:
        stream_setup_str = """
        // OpenMP execution uses synchronous host-memory copies.

        // Functions read the common-data structure directly."""
        pin_comment = "Allocate host memory for"
        dev_comment = "Allocate host memory for"
        bridge_alloc_comment = "Allocate host-memory arrays for ray chunks."
        scratch_alloc_comment = "Allocate one-dimensional host work arrays for ray-chunk state, metric, RKF45, and status values."

    # Select cleanup and result-copy operations for the CUDA or OpenMP target.
    if parallelization == "cuda":
        results_memcpy = "cudaMemcpy(results_buffer, d_results_buffer, sizeof(blueprint_data_t) * num_rays, cudaMemcpyDeviceToHost);"
        calc_blueprint = "calculate_and_fill_blueprint_data_universal(&all_photons_host, num_rays, results_buffer, 0);"
        set_intitial_con = f" set_initial_conditions_kernel_{spacetime_name}(commondata, num_rays, &all_photons_host, window_center_out, n_x_out, n_y_out, n_z_out,0);"
        stream_destroy = "cudaStreamDestroy(streams[s]); // Destroy the CUDA stream."
        free_device = "BHAH_FREE_DEVICE"
        free_pinned = "BHAH_FREE_PINNED"
    else:
        results_memcpy = "memcpy(results_buffer, d_results_buffer, sizeof(blueprint_data_t) * num_rays);"
        calc_blueprint = "calculate_and_fill_blueprint_data_universal(&all_photons_host, num_rays, results_buffer, 0);"
        set_intitial_con = f" set_initial_conditions_kernel_{spacetime_name}(commondata, num_rays, &all_photons_host, window_center_out, n_x_out, n_y_out, n_z_out);"
        stream_destroy = "// OpenMP execution has no CUDA stream to destroy."
        free_device = "BHAH_FREE"
        free_pinned = "BHAH_FREE"

    # Select memory-copy operations for the CUDA or OpenMP target.
    if parallelization == "cuda":
        # Generate asynchronous PCIe transfers for CUDA builds.
        def memcpy_async(
            dest: str, src: str, size: str, direction: str, stream: str
        ) -> str:
            return f"cudaMemcpyAsync({dest}, {src}, {size}, {direction}, {stream});"

        # Generate synchronization calls for CUDA streams.
        def stream_sync(stream: str) -> str:
            return f"cudaStreamSynchronize({stream});"

        stream_arg = ", 0"
        free_pinned = "BHAH_FREE_PINNED"
    else:
        # Generate synchronous host-memory copies for OpenMP builds.
        def memcpy_async(
            dest: str, src: str, size: str, direction: str, stream: str
        ) -> str:
            # pylint: disable=unused-argument
            return f"memcpy({dest}, {src}, {size});"

        def stream_sync(stream: str) -> str:
            # pylint: disable=unused-argument
            return r"// OpenMP memory copies complete before the next statement."

        stream_arg = ", 0"
        free_pinned = "BHAH_FREE"

    # Select work-array indices for generated function calls.
    stream_arg_current = ", current" if parallelization == "cuda" else ", current"
    stream_arg_next = ", next" if parallelization == "cuda" else ", next"

    body = rf"""
    //==========================================
    // 1. MEMORY ALLOCATION
    //==========================================

    // The master host-side Structure of Arrays (SoA) tracking all photons $f^\mu$.
    PhotonStateSoA all_photons_host;

    // {pin_comment} the state vector $f^\mu$.
    {malloc_pinned}(all_photons_host.f, sizeof(double) * 9 * num_rays);
    // {pin_comment} the previous state $f^\mu_{{n-1}}$.
    {malloc_pinned}(all_photons_host.f_p, sizeof(double) * 9 * num_rays);
    // {pin_comment} the state from two steps earlier $f^\mu_{{n-2}}$.
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

    //==========================================
    // DOUBLE-BUFFERED TRANSFER ARRAYS
    //==========================================
    // Ray-index arrays filled from coordinate-time slots.
    long int *chunk_buffer[2];
    // Host transfer array holding the state vector $f^\mu$.
    double *f_bridge[2];
    // Host transfer array holding the previous state $f^\mu_{{n-1}}$.
    double *f_p_bridge[2];
    // Host transfer array holding the state from two steps earlier $f^\mu_{{n-2}}$.
    double *f_p_p_bridge[2];
    // Host transfer array holding the affine parameter $\lambda$.
    double *affine_bridge[2];
    // Host transfer array holding the current integration step size $h$.
    double *h_bridge[2];
    // Host transfer array holding the current trajectory termination status.
    termination_type_t *status_bridge[2];
    // Host transfer array holding the number of step-size rejections.
    int *retries_bridge[2];
    // Host transfer array holding the previous observer window boundary side flag.
    bool *on_pos_window_prev_bridge[2];
    // Host transfer array holding the previous source emission boundary side flag.
    bool *on_pos_source_prev_bridge[2];
    // Host transfer array holding the historical affine parameter $\lambda_{{n-1}}$.
    double *affine_p_bridge[2];
    // Host transfer array holding the historical affine parameter $\lambda_{{n-2}}$.
    double *affine_p_p_bridge[2];
    // Host transfer array holding the observer window event lock.
    bool *window_event_found_bridge[2];
    // Host transfer array holding the source emission event lock.
    bool *source_event_found_bridge[2];

    //==========================================
    // DOUBLE-BUFFERED WORK ARRAYS
    //==========================================
    // Work array holding the physical normalization diagnostic outputs.
    normalization_constraint_t *d_norm_bundle[2];
    // Work array holding the current state vector $f^\mu$ bounding the RKF45 step.
    double *d_f_bundle[2];
    // Work array holding the anchor state vector $f_{{start}}$ for the final stage update.
    double *d_f_start_bundle[2];
    // Work array holding the cumulative RKF45 stage updates.
    double *d_f_temp_bundle[2];
    // Work array holding the history state $f^\mu_{{n-1}}$ for geometric intersection detection.
    double *d_f_prev_bundle[2];
    // Work array holding the history state $f^\mu_{{n-2}}$ for geometric intersection detection.
    double *d_f_pre_prev_bundle[2];
    // Work array holding the symmetric metric tensor $g_{{\mu\nu}}$.
    double *d_metric_bundle[2];
    // Work array holding the Christoffel symbols $\Gamma^\alpha_{{\beta\gamma}}$.
    double *d_connection_bundle[2];
    // Array storing the nine state derivatives across all six RKF45 stages.
    double *d_k_bundle[2];
    // Array holding the current integration step size $h$.
    double *d_h[2];
    // Array holding the current affine parameter $\lambda$.
    double *d_affine[2];
    // Array holding each trajectory's current termination status.
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
    // Array of ray indices $m_{{idx}}$ in each photon chunk.
    long int *d_chunk_buffer[2];

    // Allocate both sets of double-buffered arrays.
    for (int s = 0; s < 2; ++s) {{
        // {bridge_alloc_comment}
        {malloc_pinned}(chunk_buffer[s], sizeof(long int) * BUNDLE_CAPACITY); // Allocate the ray-index array.
        {malloc_pinned}(f_bridge[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate the $f^\mu$ host-transfer array.
        {malloc_pinned}(f_p_bridge[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate the $f^\mu_{{n-1}}$ host-transfer array.
        {malloc_pinned}(f_p_p_bridge[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate the $f^\mu_{{n-2}}$ host-transfer array.
        {malloc_pinned}(affine_bridge[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate the $\lambda$ host-transfer array.
        {malloc_pinned}(h_bridge[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate the $h$ host-transfer array.
        {malloc_pinned}(status_bridge[s], sizeof(termination_type_t) * BUNDLE_CAPACITY); // Allocate the status host-transfer array.
        {malloc_pinned}(retries_bridge[s], sizeof(int) * BUNDLE_CAPACITY); // Allocate the retry-count host-transfer array.
        {malloc_pinned}(on_pos_window_prev_bridge[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate the window-side host-transfer array.
        {malloc_pinned}(on_pos_source_prev_bridge[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate the source-side host-transfer array.
        {malloc_pinned}(affine_p_bridge[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate the $\lambda_{{n-1}}$ host-transfer array.
        {malloc_pinned}(affine_p_p_bridge[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate the $\lambda_{{n-2}}$ host-transfer array.
        {malloc_pinned}(window_event_found_bridge[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate the window-intersection host-transfer array.
        {malloc_pinned}(source_event_found_bridge[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate the source-intersection host-transfer array.

        // {scratch_alloc_comment}
        {malloc_device}(d_f_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate the $f^\mu$ work array.
        {malloc_device}(d_f_start_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate the $f_{{start}}$ work array.
        {malloc_device}(d_f_temp_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate the temporary stage array.
        {malloc_device}(d_f_prev_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate the $f^\mu_{{n-1}}$ work array.
        {malloc_device}(d_f_pre_prev_bundle[s], sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate the $f^\mu_{{n-2}}$ work array.
        {malloc_device}(d_metric_bundle[s], sizeof(double) * 10 * BUNDLE_CAPACITY); // Allocate the $g_{{\mu\nu}}$ work array.
        {malloc_device}(d_connection_bundle[s], sizeof(double) * 40 * BUNDLE_CAPACITY); // Allocate the $\Gamma^\alpha_{{\beta\gamma}}$ work array.
        {malloc_device}(d_k_bundle[s], sizeof(double) * 6 * 9 * BUNDLE_CAPACITY); // Allocate the derivative work array.
        {malloc_device}(d_h[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate the $h$ work array.
        {malloc_device}(d_affine[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate the $\lambda$ work array.
        {malloc_device}(d_status[s], sizeof(termination_type_t) * BUNDLE_CAPACITY); // Allocate the status work array.
        {malloc_device}(d_retries[s], sizeof(int) * BUNDLE_CAPACITY); // Allocate the retries work array.
        {malloc_device}(d_on_pos_window_prev[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate the window-flag work array.
        {malloc_device}(d_on_pos_source_prev[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate the source-flag work array.
        {malloc_device}(d_affine_prev[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate the $\lambda_{{n-1}}$ work array.
        {malloc_device}(d_affine_pre_prev[s], sizeof(double) * BUNDLE_CAPACITY); // Allocate the $\lambda_{{n-2}}$ work array.
        {malloc_device}(d_window_event_found[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate the window-lock work array.
        {malloc_device}(d_source_event_found[s], sizeof(bool) * BUNDLE_CAPACITY); // Allocate the source-lock work array.
        {malloc_device}(d_chunk_buffer[s], sizeof(long int) * BUNDLE_CAPACITY); // Allocate the ray-index work array.
        {malloc_device}(d_norm_bundle[s], sizeof(normalization_constraint_t) * BUNDLE_CAPACITY); // Allocate the diagnostic-output work array.
    }} // END LOOP: for s over 2 to allocate double-buffered arrays

    // Work array storing the final physical plane intersections.
    blueprint_data_t *d_results_buffer;
    // {dev_comment} the blueprint results buffer to avoid mid-computation memory transfers.
    {malloc_device}(d_results_buffer, sizeof(blueprint_data_t) * num_rays);

    // Host structure for binning photon trajectories $x^\mu$ by coordinate time.
    TimeSlotManager tsm;
    // Initialize coordinate-time bins for Split-Pipeline ray chunks.
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
    // Evaluate initial photon states in bounded ray chunks.

    double window_center_out[3]; // 3D array storing the spatial Cartesian coordinates $x^i$ of the observer window center.
    double n_x_out[3]; // 3D orthonormal basis vector pointing along the $x$-axis of the local window geometry.
    double n_y_out[3]; // 3D orthonormal basis vector pointing along the $y$-axis of the local window geometry.
    double n_z_out[3]; // 3D orthonormal basis vector pointing along the $z$-axis of the local window geometry.

    // Operates synchronously because the primary state array must be complete before integration starts.
    {set_intitial_con}

    //==========================================
    // INITIAL-STATE CHECK
    //==========================================
    // Scans the master Host SoA immediately following the initialization kernel call.

    long int init_mismatch_count = 0; // Number of rays with an incorrect initial state.
    long int mismatch_t = 0; // Counter tracking validation failures for the temporal coordinate $t$.
    long int mismatch_x = 0; // Counter tracking validation failures for the spatial coordinate $x$.
    long int mismatch_y = 0; // Counter tracking validation failures for the spatial coordinate $y$.
    long int mismatch_z = 0; // Counter tracking validation failures for the spatial coordinate $z$.
    long int mismatch_pt = 0; // Counter tracking validation failures for the temporal momentum $p_t$.
    long int mismatch_lam = 0; // Counter tracking validation failures for the distance traveled.

    for (long int p = 0; p < num_rays; p++) {{ // Loop iterator index $p$ mapping to a unique photon trajectory $x^\mu$ during diagnostic validation.
        const double t_check = all_photons_host.f[0 * num_rays + p]; // Evaluates the current temporal coordinate $t$ from the Host SoA.
        const double x_check = all_photons_host.f[1 * num_rays + p]; // Evaluates the current spatial coordinate $x$ from the Host SoA.
        const double y_check = all_photons_host.f[2 * num_rays + p]; // Evaluates the current spatial coordinate $y$ from the Host SoA.
        const double z_check = all_photons_host.f[3 * num_rays + p]; // Evaluates the current spatial coordinate $z$ from the Host SoA.
        const double pt_check = all_photons_host.f[4 * num_rays + p]; // Evaluates the initial temporal momentum $p_t$ from the Host SoA.
        const double lam_check = all_photons_host.f[8 * num_rays + p]; // Evaluates the initial distance traveled from the Host SoA.

        bool fail_t = fabs(t_check - commondata->t_start) > 1e-10; // Boolean flag indicating temporal coordinate $t$ validation failure.
        bool fail_x = fabs(x_check - commondata->camera_pos_x) > 1e-10; // Boolean flag indicating spatial coordinate $x$ validation failure.
        bool fail_y = fabs(y_check - commondata->camera_pos_y) > 1e-10; // Boolean flag indicating spatial coordinate $y$ validation failure.
        bool fail_z = fabs(z_check - commondata->camera_pos_z) > 1e-10; // Boolean flag indicating spatial coordinate $z$ validation failure.
        bool fail_pt = fabs(pt_check) > 1e-15; // Boolean flag indicating temporal momentum $p_t$ validation failure.
        bool fail_lam = fabs(lam_check) > 1e-15; // Boolean flag indicating distance traveled validation failure.

        if (fail_t) mismatch_t++; // Increments the validation failure counter for temporal coordinate $t$.
        if (fail_x) mismatch_x++; // Increments the validation failure counter for spatial coordinate $x$.
        if (fail_y) mismatch_y++; // Increments the validation failure counter for spatial coordinate $y$.
        if (fail_z) mismatch_z++; // Increments the validation failure counter for spatial coordinate $z$.
        if (fail_pt) mismatch_pt++; // Increments the validation failure counter for temporal momentum $p_t$.
        if (fail_lam) mismatch_lam++; // Increments the validation failure counter for the distance traveled.

        if (fail_t || fail_x || fail_y || fail_z || fail_pt || fail_lam) {{
            init_mismatch_count++; // Count this incorrect initial state.
        }} // END IF: validate initialization coordinates
    }} // END LOOP: for p over num_rays to validate initialization

    if (init_mismatch_count > 0) {{
        const double mismatch_percent = ((double)init_mismatch_count / (double)num_rays) * 100.0; // Percentage of rays with an incorrect initial state.
        // This warning reports initialization inconsistencies without halting execution.
        printf("[DIAGNOSTIC] Initial-state check: %ld out of %ld rays (%.2f%%) have incorrect coordinates, p_t, or path length.\n", init_mismatch_count, num_rays, mismatch_percent);
    }} // END IF: init_mismatch_count > 0 to print diagnostic

    long int num_batches = (num_rays + BUNDLE_CAPACITY - 1) / BUNDLE_CAPACITY; // Number of ray chunks needed for the initial constraint calculation.

    for (long int init_batch = 0; init_batch < num_batches; ++init_batch) {{ // Loop over ray chunks for the initial constraint calculation.
        long int start_idx = init_batch * BUNDLE_CAPACITY; // First ray index in the current chunk.
        long int chunk_size = NRPYMIN((long int)BUNDLE_CAPACITY, num_rays - start_idx); // Number of active trajectories in this chunk.

        for (int init_i = 0; init_i < chunk_size; ++init_i) {{ // Loop over photons in the current chunk.
            long int master_idx = start_idx + init_i; // Ray index selecting this photon in the global Structure-of-Arrays state.
            for (int init_k = 0; init_k < 9; ++init_k) {{ // Loop index $init_k$ iterating over the nine components of the state vector $f^\mu$.
                f_bridge[0][init_k * BUNDLE_CAPACITY + init_i] = all_photons_host.f[init_k * num_rays + master_idx]; // Assign the active state-vector component to the primary host-transfer array.
            }} // END LOOP: for init_k over nine state-vector components
        }} // END LOOP: for init_i over chunk_size to pack the host-transfer array

        for (int c_k = 0; c_k < 9; ++c_k) {{ // Copy the 9 state-vector components into the first work array.
            {memcpy_async("d_f_bundle[0] + c_k * BUNDLE_CAPACITY", "f_bridge[0] + c_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size", "cudaMemcpyHostToDevice", "streams[0]")}
        }} // END LOOP: for c_k over 9 to fill first work array

        // Calculate $g_{{\mu\nu}}$ for the Hamiltonian constraint.
        interpolation_kernel_{spacetime_name}(commondata,d_f_bundle[0], d_metric_bundle[0], NULL, chunk_size{stream_arg});

        //==========================================
        // DIAGNOSTIC PROBE: METRIC INTEGRITY CHECK
        //==========================================
        // Extract metric values for the current ray chunk to check numerical stability before solving for momentum.

        double *metric_diag_bridge; // Pointer storing temporary metric data to validate the interpolation sequence.
        // Temporary host array for metric $g_{{\mu\nu}}$ values.
        {malloc_pinned}(metric_diag_bridge, sizeof(double) * 10 * BUNDLE_CAPACITY);

        for (int m_k = 0; m_k < 10; ++m_k) {{
            // Loop index $m_k$ copying the 10 metric tensor $g_{{\mu\nu}}$ components.
            {memcpy_async("metric_diag_bridge + m_k * BUNDLE_CAPACITY", "d_metric_bundle[0] + m_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size", "cudaMemcpyDeviceToHost", "streams[0]")}
        }} // END LOOP: for m_k over 10 to copy metric tensor
        {stream_sync('streams[0]')}

        long int metric_nan_count = 0; // Accumulator tracking the total number of metric tensor evaluations containing non-finite values.
        for (int m_diag_i = 0; m_diag_i < chunk_size; ++m_diag_i) {{ // Loop iterator $m_diag_i$ scanning each trajectory within the current initialization chunk.
            bool m_has_nan = false; // Boolean flag indicating if the specific metric tensor $g_{{\mu\nu}}$ contains a NaN or Inf value.
            for (int m_diag_k = 0; m_diag_k < 10; ++m_diag_k) {{ // Loop index $m_diag_k$ iterating over the 10 independent components of the symmetric metric tensor $g_{{\mu\nu}}$.
                if (isnan(metric_diag_bridge[m_diag_k * BUNDLE_CAPACITY + m_diag_i]) ||
                    isinf(metric_diag_bridge[m_diag_k * BUNDLE_CAPACITY + m_diag_i])) {{
                    m_has_nan = true; // Flags the trajectory metric state as invalid due to a non-finite value.
                    break; // Terminates the tensor component loop early to conserve execution cycles upon detecting a failure.
                }} // END IF: check for NaN or Inf in metric
            }} // END LOOP: for m_diag_k over 10 to check metric tensor components
            if (m_has_nan) metric_nan_count++; // Increments the total accumulation of corrupted metric tensor evaluations.
        }} // END LOOP: for m_diag_i over chunk_size to scan for metric integrity

        if (metric_nan_count > 0) {{
            // Warn about invalid metric values without aborting photon integration.
            printf("[DIAGNOSTIC] Init Batch %ld: %ld rays have invalid Metric G_mu_nu before p_t solve.\n", init_batch, metric_nan_count);
        }} // END IF: metric_nan_count > 0 to print diagnostic
        // Free the diagnostic host-transfer array used for metric checks.
        {free_pinned}(metric_diag_bridge);

        // Solve $p_\mu p^\mu = 0$ to find temporal momentum $p_t$.
        p0_reverse_kernel(d_f_bundle[0], d_metric_bundle[0], chunk_size{stream_arg});

        for (int c_k = 0; c_k < 9; ++c_k) {{ // Copy the 9 constrained state-vector components $f^\mu$ from the first work array.
            {memcpy_async("f_bridge[0] + c_k * BUNDLE_CAPACITY", "d_f_bundle[0] + c_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size", "cudaMemcpyDeviceToHost", "streams[0]")}
        }} // END LOOP: for c_k over 9 to retrieve constrained state vector
        {stream_sync('streams[0]')}

        long int nan_count = 0; // Number of photon state vectors containing NaN after solving the null constraint.
        for (int gather_i = 0; gather_i < chunk_size; ++gather_i) {{ // Loop over photons in the retrieved chunk.
            long int master_idx = start_idx + gather_i; // Ray index selecting this photon in the global Structure-of-Arrays state.
            bool has_nan = false; // Whether this photon state vector contains a NaN value.
            for (int gather_k = 0; gather_k < 9; ++gather_k) {{ // Loop index $gather_k$ iterating over the nine components of the state vector $f^\mu$.
                double val = f_bridge[0][gather_k * BUNDLE_CAPACITY + gather_i]; // Evaluate the updated value of this state-vector component.
                all_photons_host.f[gather_k * num_rays + master_idx] = val; // Store the constrained state-vector component in the host Structure-of-Arrays state.
                if (isnan(val)) has_nan = true; // Flags the physical state vector as invalid due to a non-finite evaluation.
            }} // END LOOP: for gather_k over nine state-vector components
            if (has_nan) nan_count++; // Increments the total count of unresolved physical state vectors $f^\mu$.
        }} // END LOOP: for gather_i over chunk_size to retrieve updated constrained state vectors

        if (nan_count > 0) {{
            // This is a soft warning alerting to unresolved constraints $p_\mu p^\mu = 0$ for isolated trajectories.
            printf("[DIAGNOSTIC] Init Batch %ld: %ld rays contain NaN in state f^mu after p_t solve.\n", init_batch, nan_count);
        }} // END IF: nan_count > 0 to print diagnostic
    }} // END LOOP: for init_batch over num_batches to evaluate initialization constraints

    //==========================================
    // BASELINE CONSERVED QUANTITIES
    //==========================================
    // Evaluate initial conserved quantities immediately after generating valid physical null states.

    if (commondata->perform_conservation_check) {{
        // Compute conserved quantities from initialized photon states before RKF45 integration.
        calculate_conserved_quantities_universal_{spacetime_name}_photon(commondata, &all_photons_host, num_rays, initial_cq_host);
    }} // END IF: perform_conservation_check to evaluate baseline conserved quantities

    long int sync_i; // Loop iterator index $sync_i$ spanning the entire global ray count to synchronize starting properties across history states.
    for(sync_i = 0; sync_i < num_rays; ++sync_i) {{
        int sync_k; // Loop index $sync_k$ over the 9 components copied into the two history states.
        for (sync_k = 0; sync_k < 9; ++sync_k) {{
            all_photons_host.f_p[sync_k * num_rays + sync_i] = all_photons_host.f[sync_k * num_rays + sync_i]; // Set $f^\mu_{{n-1}}$ to the initial state.
            all_photons_host.f_p_p[sync_k * num_rays + sync_i] = all_photons_host.f[sync_k * num_rays + sync_i]; // Set $f^\mu_{{n-2}}$ to the initial state.
        }} // END LOOP: for sync_k over 9 to initialize both history states
        all_photons_host.status[sync_i] = ACTIVE; // Mark the photon trajectory active.
        all_photons_host.affine_param[sync_i] = 0.0; // Sets the initial baseline progression scalar for the affine parameter $\lambda$.
        all_photons_host.rejection_retries[sync_i] = 0; // Clear the rejected-step count.

        all_photons_host.affine_param_p[sync_i] = 0.0; // Initializes the first historical affine parameter $\lambda_{{n-1}}$.
        all_photons_host.affine_param_p_p[sync_i] = 0.0; // Initializes the second historical affine parameter $\lambda_{{n-2}}$.
        all_photons_host.window_event_found[sync_i] = false; // Sets the observer window intersection logical lock to false.
        all_photons_host.source_event_found[sync_i] = false; // Sets the source emission intersection logical lock to false.

        int s_idx = slot_get_index(&tsm, all_photons_host.f[sync_i]); // Integer index $s_{{idx}}$ mapping the current photon's temporal coordinate $t$ to a discrete execution bin in the TimeSlotManager.
        if (s_idx != -1) {{
            slot_add_photon(&tsm, s_idx, sync_i); // Registers the active photon index to its corresponding temporal bin mapped by the driver.
        }} // END IF: s_idx != -1 to add photon to slot
    }} // END LOOP: for sync_i over num_rays to synchronize starting properties

    // Monotonic-clock value at the start of an integration chunk.
    struct timespec batch_start_time;
    // Monotonic-clock value at the end of an integration chunk.
    struct timespec batch_end_time;
    // Read the monotonic clock before temporal-loop execution.
    clock_gettime(CLOCK_MONOTONIC, &batch_start_time);

    // Reserve terminal lines for the progress display.
    printf("\n\n\n\n\n\n\n");

    //==========================================
    // 3. COORDINATE-TIME LOOP
    //==========================================
    // Integer tracking the global number of active photon trajectories to allow early loop termination.
    long int total_active_photons = num_rays;

    // Outer loop iterator for the physical time bins.
    for (int slot_idx = tsm.num_slots - 1; slot_idx >= 0; --slot_idx) {{
        // Stop advancing coordinate-time slots after all photon trajectories terminate.
        if (total_active_photons <= 0) {{
            break; // No active photon remains in later coordinate-time slots.
        }} // END IF: total_active_photons <= 0 to break

        int current = 0; // Index of the work arrays holding the active photon chunk.
        int next = 1; // Index of the work arrays prepared for the next photon chunk.
        long int active_chunks[2] = {{ 0, 0}}; // Number of trajectories in each work-array set.

        //==========================================
        // PHASE A: PROCESS FIRST RAY CHUNK
        //==========================================
        // Populate the first transfer array and process it with the current work arrays.

        active_chunks[current] = NRPYMIN((long int)BUNDLE_CAPACITY, tsm.slot_counts[slot_idx]); // Bound the active chunk by the allocated work-array capacity.

        if (active_chunks[current] > 0) {{
            slot_remove_chunk(&tsm, slot_idx, chunk_buffer[current], active_chunks[current]); // Remove ray indices from the current coordinate-time slot.

            // 1. Pack the nine-component state vectors using forward sweeps
            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ iterating over the nine components of the state vectors.
                for (int bridge_i = 0; bridge_i < active_chunks[current]; ++bridge_i) {{  // Loop iterator $bridge_i$ packing each photon state into the Host-side transfer arrays.
                    long int m_idx = chunk_buffer[current][bridge_i]; // Ray index $m_{{idx}}$ selecting this photon in the host Structure-of-Arrays state.
                    f_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f[c_k * num_rays + m_idx]; // Pack the photon state vector $f^\mu$ into the host-transfer array.
                    f_p_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p[c_k * num_rays + m_idx]; // Pack $f^\mu_{{n-1}}$ into the host-transfer array.
                    f_p_p_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p_p[c_k * num_rays + m_idx]; // Pack $f^\mu_{{n-2}}$ into the host-transfer array.
                }} // END LOOP: for bridge_i over active_chunks[current] to pack photon states
            }} // END LOOP: for c_k over nine state-vector components

            // 2. Pack the 1D arrays in a separate sequential loop
            for (int bridge_i = 0; bridge_i < active_chunks[current]; ++bridge_i) {{
                long int m_idx = chunk_buffer[current][bridge_i];
                h_bridge[current][bridge_i] = all_photons_host.h[m_idx]; // Packs the current integration step size $h$ into the host-transfer array.
                status_bridge[current][bridge_i] = all_photons_host.status[m_idx]; // Packs the trajectory status enum into the host-transfer array.
                retries_bridge[current][bridge_i] = all_photons_host.rejection_retries[m_idx]; // Packs the error rejection scalar into the host-transfer array.
                affine_bridge[current][bridge_i] = all_photons_host.affine_param[m_idx]; // Packs the affine parameter $\lambda$ into the host-transfer array.
                on_pos_window_prev_bridge[current][bridge_i] = all_photons_host.on_positive_side_of_window_prev[m_idx]; // Packs the observer window boundary flag into the host-transfer array.
                on_pos_source_prev_bridge[current][bridge_i] = all_photons_host.on_positive_side_of_source_prev[m_idx]; // Packs the source emission boundary flag into the host-transfer array.
                affine_p_bridge[current][bridge_i] = all_photons_host.affine_param_p[m_idx]; // Packs the historical affine parameter $\lambda_{{ n-1}}$ into the host-transfer array.
                affine_p_p_bridge[current][bridge_i] = all_photons_host.affine_param_p_p[m_idx]; // Packs the historical affine parameter $\lambda_{{ n-2}}$ into the host-transfer array.
                window_event_found_bridge[current][bridge_i] = all_photons_host.window_event_found[m_idx]; // Packs the observer window intersection lock into the host-transfer array.
                source_event_found_bridge[current][bridge_i] = all_photons_host.source_event_found[m_idx]; // Packs the source emission intersection lock into the host-transfer array.
            }} // END LOOP: for bridge_i over active_chunks[current] to pack 1D arrays

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Copy the 9 state-vector components into the current work arrays.
                // Copy $f^\mu$ to the current work array.
                {memcpy_async("d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "f_bridge[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
                // Copy $f^\mu_{{n-1}}$ to the current work array.
                {memcpy_async("d_f_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "f_p_bridge[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
                // Copy $f^\mu_{{n-2}}$ to the current work array.
                {memcpy_async("d_f_pre_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "f_p_p_bridge[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            }} // END LOOP: for c_k over 9 to fill current work arrays
            // Copy step sizes $h$ to the current work array.
            {memcpy_async("d_h[current]", "h_bridge[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Copy trajectory statuses to the current work array.
            {memcpy_async("d_status[current]", "status_bridge[current]", "sizeof(termination_type_t) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Copy rejection counts to the current work array.
            {memcpy_async("d_retries[current]", "retries_bridge[current]", "sizeof(int) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Copy affine parameters $\lambda$ to the current work array.
            {memcpy_async("d_affine[current]", "affine_bridge[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Copy previous window-side flags to the current work array.
            {memcpy_async("d_on_pos_window_prev[current]", "on_pos_window_prev_bridge[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Copy previous source-side flags to the current work array.
            {memcpy_async("d_on_pos_source_prev[current]", "on_pos_source_prev_bridge[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Copy $\lambda_{{n-1}}$ to the current work array.
            {memcpy_async("d_affine_prev[current]", "affine_p_bridge[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Copy $\lambda_{{n-2}}$ to the current work array.
            {memcpy_async("d_affine_pre_prev[current]", "affine_p_p_bridge[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Copy window-intersection flags to the current work array.
            {memcpy_async("d_window_event_found[current]", "window_event_found_bridge[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Copy source-intersection flags to the current work array.
            {memcpy_async("d_source_event_found[current]", "source_event_found_bridge[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}
            // Copy global ray indices $m_{{idx}}$ to the current work array.
            {memcpy_async("d_chunk_buffer[current]", "chunk_buffer[current]", "sizeof(long int) * active_chunks[current]", "cudaMemcpyHostToDevice", "streams[current]")}

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Set the two RKF45 state arrays from the current photon state.
                // Copy $f^\mu$ into the RKF45 starting state.
                {memcpy_async("d_f_start_bundle[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToDevice", "streams[current]")}
                // Copy $f^\mu$ into the RKF45 stage state.
                {memcpy_async("d_f_temp_bundle[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToDevice", "streams[current]")}
            }} // END LOOP: for c_k over 9 to set current RKF45 states

            for (int stage = 1; stage <= 6; ++stage) {{  // Loop iterator $stage$ executing the 6 discrete stages of the RKF45 Runge-Kutta numerical solver.
                // Evaluate metric tensor $g_{{ \mu\nu}}$ and connection $\Gamma^\alpha_{{ \beta\gamma}}$ for this stage.
                interpolation_kernel_{ spacetime_name}(commondata, d_f_temp_bundle[current], d_metric_bundle[current], d_connection_bundle[current], active_chunks[current]{stream_arg_current});
                // Compute geodesic-equation right-hand-side derivatives $\dot{{ f}}^\mu$.
                calculate_ode_rhs_kernel(d_f_temp_bundle[current], d_metric_bundle[current], d_connection_bundle[current], d_k_bundle[current], stage, active_chunks[current]{stream_arg_current});
                // Accumulate this RKF45 stage.
                rkf45_stage_update(d_f_start_bundle[current], d_k_bundle[current], d_h[current], stage, active_chunks[current], d_f_temp_bundle[current]{stream_arg_current});
            }} // END LOOP: for stage over 6 to execute RKF45 stages

            // Apply Cash-Karp error control and update step size $h$.
            rkf45_finalize_and_control(commondata, d_f_bundle[current], d_f_start_bundle[current], d_k_bundle[current], d_h[current], d_status[current], d_affine[current], d_retries[current], active_chunks[current]{stream_arg_current});
            // Detect geometric events and record intersection coordinates.
            event_detection_manager_kernel(commondata, d_f_bundle[current], d_f_prev_bundle[current], d_f_pre_prev_bundle[current], d_affine[current], d_affine_prev[current], d_affine_pre_prev[current], d_results_buffer, d_status[current], d_on_pos_window_prev[current], d_on_pos_source_prev[current], d_window_event_found[current], d_source_event_found[current], d_chunk_buffer[current], active_chunks[current]{stream_arg_current});

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Copy the 9 updated state-vector components from the current work arrays.
                // Copy updated $f^\mu$ from the current work array.
                {memcpy_async("f_bridge[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
                // Copy updated $f^\mu_{{n-1}}$ from the current work array.
                {memcpy_async("f_p_bridge[current] + c_k * BUNDLE_CAPACITY", "d_f_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
                // Copy updated $f^\mu_{{n-2}}$ from the current work array.
                {memcpy_async("f_p_p_bridge[current] + c_k * BUNDLE_CAPACITY", "d_f_pre_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            }} // END LOOP: for c_k over 9 to retrieve current states
            // Copy updated step sizes $h$ from the current work array.
            {memcpy_async("h_bridge[current]", "d_h[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Copy updated trajectory statuses from the current work array.
            {memcpy_async("status_bridge[current]", "d_status[current]", "sizeof(termination_type_t) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Copy updated rejection counts from the current work array.
            {memcpy_async("retries_bridge[current]", "d_retries[current]", "sizeof(int) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Copy updated affine parameters $\lambda$ from the current work array.
            {memcpy_async("affine_bridge[current]", "d_affine[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Copy updated window-side flags from the current work array.
            {memcpy_async("on_pos_window_prev_bridge[current]", "d_on_pos_window_prev[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Copy updated source-side flags from the current work array.
            {memcpy_async("on_pos_source_prev_bridge[current]", "d_on_pos_source_prev[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Copy $\lambda_{{n-1}}$ from the current work array.
            {memcpy_async("affine_p_bridge[current]", "d_affine_prev[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Copy $\lambda_{{n-2}}$ from the current work array.
            {memcpy_async("affine_p_p_bridge[current]", "d_affine_pre_prev[current]", "sizeof(double) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Copy window-intersection flags from the current work array.
            {memcpy_async("window_event_found_bridge[current]", "d_window_event_found[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
            // Copy source-intersection flags from the current work array.
            {memcpy_async("source_event_found_bridge[current]", "d_source_event_found[current]", "sizeof(bool) * active_chunks[current]", "cudaMemcpyDeviceToHost", "streams[current]")}
        }} // END IF: active_chunks[current] > 0 to process first ray chunk

        //==========================================
        // PHASE B: THE OVERLAP LOOP
        //==========================================
        // Alternate work-array sets while processing consecutive ray chunks.

        while (active_chunks[current] > 0 || tsm.slot_counts[slot_idx] > 0) {{

            active_chunks[next] = NRPYMIN((long int)BUNDLE_CAPACITY, tsm.slot_counts[slot_idx]); // Bound the next ray chunk by work-array capacity.
            if (active_chunks[next] > 0) {{
                slot_remove_chunk(&tsm, slot_idx, chunk_buffer[next], active_chunks[next]); // Remove the next ray indices from the current coordinate-time slot.

                // 1. Pack the nine-component state vectors using forward sweeps
                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ iterating over the nine components of the state vectors.
                    for (int bridge_i = 0; bridge_i < active_chunks[next]; ++bridge_i) {{  // Loop iterator $bridge_i$ packing each photon state into the next host-transfer array.
                        long int m_idx = chunk_buffer[next][bridge_i]; // Ray index $m_{{idx}}$ selecting this photon in the host Structure-of-Arrays state.
                        f_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f[c_k * num_rays + m_idx]; // Pack the photon state vector $f^\mu$ into the host-transfer array.
                        f_p_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p[c_k * num_rays + m_idx]; // Pack $f^\mu_{{n-1}}$ into the host-transfer array.
                        f_p_p_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p_p[c_k * num_rays + m_idx]; // Pack $f^\mu_{{n-2}}$ into the host-transfer array.
                    }} // END LOOP: for bridge_i over active_chunks[next] to pack photon states
                }} // END LOOP: for c_k over nine state-vector components

                // 2. Pack the 1D arrays in a separate sequential loop
                for (int bridge_i = 0; bridge_i < active_chunks[next]; ++bridge_i) {{
                    long int m_idx = chunk_buffer[next][bridge_i];
                    h_bridge[next][bridge_i] = all_photons_host.h[m_idx]; // Packs the current integration step size $h$ into the host-transfer array.
                    status_bridge[next][bridge_i] = all_photons_host.status[m_idx]; // Packs the trajectory status enum into the host-transfer array.
                    retries_bridge[next][bridge_i] = all_photons_host.rejection_retries[m_idx]; // Packs the error rejection scalar into the host-transfer array.
                    affine_bridge[next][bridge_i] = all_photons_host.affine_param[m_idx]; // Packs the affine parameter $\lambda$ into the host-transfer array.
                    on_pos_window_prev_bridge[next][bridge_i] = all_photons_host.on_positive_side_of_window_prev[m_idx]; // Packs the observer window boundary flag into the host-transfer array.
                    on_pos_source_prev_bridge[next][bridge_i] = all_photons_host.on_positive_side_of_source_prev[m_idx]; // Packs the source emission boundary flag into the host-transfer array.
                    affine_p_bridge[next][bridge_i] = all_photons_host.affine_param_p[m_idx]; // Packs the historical affine parameter $\lambda_{{ n-1}}$ into the host-transfer array.
                    affine_p_p_bridge[next][bridge_i] = all_photons_host.affine_param_p_p[m_idx]; // Packs the historical affine parameter $\lambda_{{ n-2}}$ into the host-transfer array.
                    window_event_found_bridge[next][bridge_i] = all_photons_host.window_event_found[m_idx]; // Packs the observer window intersection lock into the host-transfer array.
                    source_event_found_bridge[next][bridge_i] = all_photons_host.source_event_found[m_idx]; // Packs the source emission intersection lock into the host-transfer array.
                }} // END LOOP: for bridge_i over active_chunks[next] to pack 1D arrays

                for (int c_k = 0; c_k < 9; ++c_k) {{  // Copy the next chunk's 9 state-vector components into the next work arrays.
                    // Copy $f^\mu$ to the next work array.
                    {memcpy_async("d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "f_bridge[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                    // Copy $f^\mu_{{n-1}}$ to the next work array.
                    {memcpy_async("d_f_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "f_p_bridge[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                    // Copy $f^\mu_{{n-2}}$ to the next work array.
                    {memcpy_async("d_f_pre_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "f_p_p_bridge[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                }} // END LOOP: for c_k over 9 to fill next work arrays
                // Copy step sizes $h$ to the next work array.
                {memcpy_async("d_h[next]", "h_bridge[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Copy trajectory statuses to the next work array.
                {memcpy_async("d_status[next]", "status_bridge[next]", "sizeof(termination_type_t) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Copy rejection counts to the next work array.
                {memcpy_async("d_retries[next]", "retries_bridge[next]", "sizeof(int) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Copy affine parameters $\lambda$ to the next work array.
                {memcpy_async("d_affine[next]", "affine_bridge[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Copy previous window-side flags to the next work array.
                {memcpy_async("d_on_pos_window_prev[next]", "on_pos_window_prev_bridge[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Copy previous source-side flags to the next work array.
                {memcpy_async("d_on_pos_source_prev[next]", "on_pos_source_prev_bridge[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Copy $\lambda_{{n-1}}$ to the next work array.
                {memcpy_async("d_affine_prev[next]", "affine_p_bridge[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Copy $\lambda_{{n-2}}$ to the next work array.
                {memcpy_async("d_affine_pre_prev[next]", "affine_p_p_bridge[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Copy window-intersection flags to the next work array.
                {memcpy_async("d_window_event_found[next]", "window_event_found_bridge[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Copy source-intersection flags to the next work array.
                {memcpy_async("d_source_event_found[next]", "source_event_found_bridge[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}
                // Copy global ray indices $m_{{idx}}$ to the next work array.
                {memcpy_async("d_chunk_buffer[next]", "chunk_buffer[next]", "sizeof(long int) * active_chunks[next]", "cudaMemcpyHostToDevice", "streams[next]")}

                for (int c_k = 0; c_k < 9; ++c_k) {{  // Set the two RKF45 state arrays from the next photon state.
                    // Copy $f^\mu$ into the RKF45 starting state.
                    {memcpy_async("d_f_start_bundle[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToDevice", "streams[next]")}
                    // Copy $f^\mu$ into the RKF45 stage state.
                    {memcpy_async("d_f_temp_bundle[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToDevice", "streams[next]")}
                }} // END LOOP: for c_k over 9 to set next RKF45 states

                for (int stage = 1; stage <= 6; ++stage) {{  // Loop iterator $stage$ executing the 6 discrete stages of the upcoming RKF45 Runge-Kutta numerical solver.
                    // Evaluate metric tensor $g_{{ \mu\nu}}$ and connection $\Gamma^\alpha_{{ \beta\gamma}}$ for this stage.
                    interpolation_kernel_{spacetime_name}(commondata, d_f_temp_bundle[next], d_metric_bundle[next], d_connection_bundle[next],active_chunks[next]{stream_arg_next});
                    // Compute geodesic-equation right-hand-side derivatives $\dot{{ f}}^\mu$.
                    calculate_ode_rhs_kernel(d_f_temp_bundle[next], d_metric_bundle[next], d_connection_bundle[next], d_k_bundle[next], stage, active_chunks[next]{stream_arg_next});
                    // Accumulate this RKF45 stage.
                    rkf45_stage_update(d_f_start_bundle[next], d_k_bundle[next], d_h[next], stage, active_chunks[next], d_f_temp_bundle[next]{stream_arg_next});
                }} // END LOOP: for stage over 6 to execute RKF45 stages

                // Apply Cash-Karp error control and update step size $h$.
                rkf45_finalize_and_control(commondata, d_f_bundle[next], d_f_start_bundle[next], d_k_bundle[next], d_h[next], d_status[next], d_affine[next], d_retries[next], active_chunks[next]{stream_arg_next});
                // Detect geometric events and record intersection coordinates.
                event_detection_manager_kernel(commondata, d_f_bundle[next], d_f_prev_bundle[next], d_f_pre_prev_bundle[next], d_affine[next], d_affine_prev[next], d_affine_pre_prev[next], d_results_buffer, d_status[next], d_on_pos_window_prev[next], d_on_pos_source_prev[next], d_window_event_found[next], d_source_event_found[next], d_chunk_buffer[next], active_chunks[next]{stream_arg_next});
                for (int c_k = 0; c_k < 9; ++c_k) {{  // Copy the next chunk's 9 updated state-vector components from the next work arrays.
                    // Copy updated $f^\mu$ from the next work array.
                    {memcpy_async("f_bridge[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                    // Copy updated $f^\mu_{{n-1}}$ from the next work array.
                    {memcpy_async("f_p_bridge[next] + c_k * BUNDLE_CAPACITY", "d_f_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                    // Copy updated $f^\mu_{{n-2}}$ from the next work array.
                    {memcpy_async("f_p_p_bridge[next] + c_k * BUNDLE_CAPACITY", "d_f_pre_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                }} // END LOOP: for c_k over 9 to retrieve next states
                // Copy updated step sizes $h$ from the next work array.
                {memcpy_async("h_bridge[next]", "d_h[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Copy updated trajectory statuses from the next work array.
                {memcpy_async("status_bridge[next]", "d_status[next]", "sizeof(termination_type_t) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Copy updated rejection counts from the next work array.
                {memcpy_async("retries_bridge[next]", "d_retries[next]", "sizeof(int) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Copy updated affine parameters $\lambda$ from the next work array.
                {memcpy_async("affine_bridge[next]", "d_affine[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Copy updated window-side flags from the next work array.
                {memcpy_async("on_pos_window_prev_bridge[next]", "d_on_pos_window_prev[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Copy updated source-side flags from the next work array.
                {memcpy_async("on_pos_source_prev_bridge[next]", "d_on_pos_source_prev[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Copy $\lambda_{{n-1}}$ from the next work array.
                {memcpy_async("affine_p_bridge[next]", "d_affine_prev[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Copy $\lambda_{{n-2}}$ from the next work array.
                {memcpy_async("affine_p_p_bridge[next]", "d_affine_pre_prev[next]", "sizeof(double) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Copy window-intersection flags from the next work array.
                {memcpy_async("window_event_found_bridge[next]", "d_window_event_found[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
                // Copy source-intersection flags from the next work array.
                {memcpy_async("source_event_found_bridge[next]", "d_source_event_found[next]", "sizeof(bool) * active_chunks[next]", "cudaMemcpyDeviceToHost", "streams[next]")}
            }} // END IF: active_chunks[next] > 0 to process upcoming ray chunk

            if (active_chunks[current] > 0) {{
                // Wait for the current work set before unpacking this ray chunk.
                {stream_sync("streams[current]")}

                // 1. Unpack nine-component state vectors sequentially
                for (int fin_k = 0; fin_k < 9; ++fin_k) {{  // Copy the nine state-vector components into the host Structure-of-Arrays state.
                    for (int fin_i = 0; fin_i < active_chunks[current]; ++fin_i) {{  // Copy each completed photon state into its host-array entry.
                        long int m_idx = chunk_buffer[current][fin_i]; // Ray index $m_{{idx}}$ selecting the photon in the host Structure-of-Arrays state.
                        all_photons_host.f[fin_k * num_rays + m_idx] = f_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i]; // Copy the synchronized state-vector component $f^\mu$ into the host array.
                        all_photons_host.f_p[fin_k * num_rays + m_idx] = f_p_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i]; // Unpack synchronized $f^\mu_{{n-1}}$ into the global host array.
                        all_photons_host.f_p_p[fin_k * num_rays + m_idx] = f_p_p_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i]; // Unpack synchronized $f^\mu_{{n-2}}$ into the global host array.
                    }} // END LOOP: for fin_i over active_chunks[current] to unpack finalized data
                }} // END LOOP: for fin_k over nine state-vector components

                // 2. Unpack 1D arrays sequentially
                for (int fin_i = 0; fin_i < active_chunks[current]; ++fin_i) {{
                    long int m_idx = chunk_buffer[current][fin_i];
                    all_photons_host.h[m_idx] = h_bridge[current][fin_i]; // Copy the synchronized step size $h$ into the host array.
                    all_photons_host.status[m_idx] = status_bridge[current][fin_i]; // Copy the synchronized trajectory status into the host array.
                    all_photons_host.rejection_retries[m_idx] = retries_bridge[current][fin_i]; // Copy the synchronized rejection count into the host array.
                    all_photons_host.affine_param[m_idx] = affine_bridge[current][fin_i]; // Copy the synchronized affine parameter $\lambda$ into the host array.
                    all_photons_host.on_positive_side_of_window_prev[m_idx] = on_pos_window_prev_bridge[current][fin_i]; // Copy the synchronized window-side flag into the host array.
                    all_photons_host.on_positive_side_of_source_prev[m_idx] = on_pos_source_prev_bridge[current][fin_i]; // Copy the synchronized source-side flag into the host array.
                    all_photons_host.affine_param_p[m_idx] = affine_p_bridge[current][fin_i]; // Copy the synchronized affine parameter $\lambda_{{n-1}}$ into the host array.
                    all_photons_host.affine_param_p_p[m_idx] = affine_p_p_bridge[current][fin_i]; // Copy the synchronized affine parameter $\lambda_{{n-2}}$ into the host array.
                    all_photons_host.window_event_found[m_idx] = window_event_found_bridge[current][fin_i]; // Copy the window-intersection flag into the host array.
                    all_photons_host.source_event_found[m_idx] = source_event_found_bridge[current][fin_i]; // Copy the source-intersection flag into the host array.
                }} // END LOOP: for fin_i over active_chunks[current] to unpack 1D arrays

                // 3. TimeSlotManager State Update (Cache-hot, strictly sequential)
                for (int fin_i = 0; fin_i < active_chunks[current]; ++fin_i) {{
                    long int m_idx = chunk_buffer[current][fin_i];
                    if (status_bridge[current][fin_i] == ACTIVE) {{  // Evaluates the continuation logic if the trajectory remains within safe physical bounds.
                        int next_s_idx = slot_get_index(&tsm, all_photons_host.f[m_idx]); // Use the updated coordinate time $t$ to select the next time slot.
                        if (next_s_idx != -1) {{  // Confirms the physical state has not exceeded the maximum simulation time bounds.
                            slot_add_photon(&tsm, next_s_idx, m_idx);  // Add the ray index to its next coordinate-time slot.
                        }} else {{
                            all_photons_host.status[m_idx] = FAILURE_T_MAX_EXCEEDED; // Flags the physical state as permanently failed due to excessive propagation time.
                            total_active_photons--; // Decrements the global counter as the physical trajectory has reached a terminal state.
                        }} // END ELSE: flag state as failed and decrement total active photons
                    }} // END IF: trajectory remains active
                    else if (status_bridge[current][fin_i] == REJECTED) {{   // Evaluates the retry logic if the numerical step exceeded the requested tolerances.
                        slot_add_photon(&tsm, slot_idx, m_idx); // Re-adds to the current bin to attempt integration with an adapted step-size scalar $h$.
                    }} else {{
                        total_active_photons--; // Decrements the global counter as the physical trajectory has reached a terminal state.
                    }} // END ELSE: trajectory reached terminal state
                }} // END LOOP: for fin_i over active_chunks[current] to update TimeSlotManager state
            }} // END IF: active_chunks[current] > 0 to complete and unpack current chunk

            //==========================================
            // PROGRESS DISPLAY
            //==========================================
            // Calculate integration rate and completed-ray fraction.

            // Read the monotonic clock after the current chunk.
            clock_gettime(CLOCK_MONOTONIC, &batch_end_time);

            // Evaluates the absolute wall-clock duration of the integration chunk in seconds.
            double elapsed_sec = (batch_end_time.tv_sec - batch_start_time.tv_sec) + (batch_end_time.tv_nsec - batch_start_time.tv_nsec) / 1e9;

            // Calculate raw integration steps per second.
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
                else if (bar_i == pos) bar[bar_i] = '>'; // Mark the current progress position.
                else bar[bar_i] = ' '; // Appends the uncompleted progression character.
            }} // END LOOP: for bar_i over bar_width to construct loading bar
            bar[bar_width] = '\0'; // Terminates the loading bar character array to prevent buffer overruns.

            // Accumulator tracking the total number of adaptive step size $h$ rejections in the active chunk.
            long int batch_rejections = 0;

            // Loop iterator $sum_i$ scanning the finalized physical state host-transfer array for error tolerance failures.
            for (int sum_i = 0; sum_i < active_chunks[current]; ++sum_i) {{
                batch_rejections += retries_bridge[current][sum_i]; // Accumulates the localized step rejection tally.
            }} // END LOOP: for sum_i over active_chunks[current] to calculate rejection count

            // Evaluates the relative frequency of adaptive step size $h$ rejections.
            double reject_percent = (active_chunks[current] > 0) ? (100.0 * (double)batch_rejections / (double)active_chunks[current]) : 0.0;

            // Move the terminal cursor up to overwrite the previous progress display.
            printf("\033[7A");
            printf("--------------------------------------------------\n"); // Print upper border of the progress display.
            printf(" Progress:   [%s] %5.1f%% \033[K\n", bar, percent_done); // Prints the global completion loading bar and percentage.
            printf(" Active:     %ld / %ld \033[K\n", total_active_photons, num_rays); // Prints the remaining active photon trajectories $x^\mu$.
            printf(" Slot Time:  Slot %d (t = %.1f) \033[K\n", slot_idx, current_t); // Prints the current physical temporal bin coordinate $t$.
            printf(" Speed:      %.2e integration steps/s \033[K\n", steps_per_sec); // Print integration throughput.
            printf(" Rejects:    %ld (%.1f%%) \033[K\n", batch_rejections, reject_percent); // Prints the adaptive step size $h$ rejection frequency.
            printf("--------------------------------------------------\n"); // Print lower border of the progress display.
            fflush(stdout); // Flushes the standard output buffer to ensure instantaneous terminal rendering.

            // Read the monotonic clock before the next integration chunk.
            clock_gettime(CLOCK_MONOTONIC, &batch_start_time);
            // --------------------------

            active_chunks[current] = 0; // Mark the current work-array set empty.
            int temp = current; // Save the current work-array index during the swap.
            current = next; // Use the prepared work-array set for the next ray chunk.
            next = temp; // Reuse the completed work-array set for a later ray chunk.
        }} // END WHILE: alternate work arrays to process temporal bin

        //==========================================
        // PHASE C: ADVANCE THE TIME SLOT
        //==========================================
        // Complete all parallel calls and copies before advancing the coordinate-time slot.
        BHAH_DEVICE_SYNC(); // Wait for all photon states in this slot.

     }} // END LOOP: for slot_idx down to 0 to process all temporal bins

    //==========================================
    // 4. Conserved Values & CLEANUP & FINALIZATION
    //==========================================
    // Process terminal photon trajectories and extract final geometric intersections.

        // Copy final geometric intersections $b_i$ into the caller's output array.
        {results_memcpy}

        // Process escaped photons intersecting the celestial sphere $r > r_{{escape}}$.
        {calc_blueprint}

        //==========================================
        // TERMINAL NORMALIZATION DIAGNOSTIC
        //==========================================
        // Evaluate terminal normalization constraint.
        if (commondata->perform_conservation_check) {{
            normalization_constraint_t *norm_diag_bridge; // Host array for diagnostic normalization values.
            {malloc_pinned}(norm_diag_bridge, sizeof(normalization_constraint_t) * BUNDLE_CAPACITY); // Allocate one ray chunk of normalization values.

            double max_err_norm = 0.0; // Scalar tracking the maximum absolute drift from the expected normalization constraint invariant $C$.
            long int worst_ray_norm = -1; // Absolute index identifying the trajectory $x^\mu$ with the highest constraint violation.

            long int norm_num_batches = (num_rays + BUNDLE_CAPACITY - 1) / BUNDLE_CAPACITY; // Integer calculation defining the total sequential blocks required to process all photon trajectories.

            for (long int norm_batch = 0; norm_batch < norm_num_batches; ++norm_batch) {{ // Loop iterator $norm_batch$ for evaluating the terminal normalization constraint across sequential chunks.
                long int start_idx = norm_batch * BUNDLE_CAPACITY; // First ray index in this normalization chunk of the host Structure-of-Arrays state.
                long int chunk_size = NRPYMIN((long int)BUNDLE_CAPACITY, num_rays - start_idx); // Number of active trajectories in this chunk.

                for (int norm_i = 0; norm_i < chunk_size; ++norm_i) {{ // Loop index $norm_i$ iterating over the specific normalization batch elements to pack the host-transfer array.
                    long int master_idx = start_idx + norm_i; // Ray index $m_{{idx}}$ selecting the photon in the host Structure-of-Arrays state.
                    for (int norm_k = 0; norm_k < 9; ++norm_k) {{ // Loop index $norm_k$ iterating over the nine components of the state vector $f^\mu$.
                        f_bridge[0][norm_k * BUNDLE_CAPACITY + norm_i] = all_photons_host.f[norm_k * num_rays + master_idx]; // Assign the terminal state-vector component to the primary host-transfer array.
                    }} // END LOOP: for norm_k over nine terminal state-vector components
                }} // END LOOP: for norm_i over chunk_size to pack host-transfer array

                for (int c_k = 0; c_k < 9; ++c_k) {{ // Copy the nine state-vector components $f^\mu$ into the first work array.
                    {memcpy_async("d_f_bundle[0] + c_k * BUNDLE_CAPACITY", "f_bridge[0] + c_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size", "cudaMemcpyHostToDevice", "streams[0]")}
                }} // END LOOP: for c_k over 9 state-vector components

                // Calculate $g_{{\mu\nu}}$ for the normalization constraint.
                interpolation_kernel_{spacetime_name}(commondata, d_f_bundle[0], d_metric_bundle[0], NULL, chunk_size{stream_arg});

                // Compute $C = g_{{\mu\nu}} p^\mu p^\nu$.
                normalization_constraint_photon(d_f_bundle[0], d_metric_bundle[0], d_norm_bundle[0], chunk_size{stream_arg});

                // Copy diagnostic normalization values from the first work array.
                {memcpy_async("norm_diag_bridge", "d_norm_bundle[0]", "sizeof(normalization_constraint_t) * chunk_size", "cudaMemcpyDeviceToHost", "streams[0]")}

                // Wait for normalization values before checking this ray chunk.
                {stream_sync('streams[0]')}

                for (int norm_i = 0; norm_i < chunk_size; ++norm_i) {{ // Loop iterator $norm_i$ scanning each trajectory within the current diagnostic memory chunk.
                    double current_norm_err = fabs(norm_diag_bridge[norm_i].C); // Evaluates the absolute numerical drift for the physical scalar invariant $C$.
                    if (current_norm_err > max_err_norm) {{
                        max_err_norm = current_norm_err; // Updates the maximum tracked absolute error for the geometric normalization constraint.
                        worst_ray_norm = start_idx + norm_i; // Updates the absolute master index $m_{{idx}}$ associated with the maximum geometric constraint violation.
                    }} // END IF: current_norm_err > max_err_norm
                }} // END LOOP: for norm_i over chunk_size to scan trajectory and evaluate constraint violation
            }} // END LOOP: for norm_batch over norm_num_batches to evaluate terminal normalization constraint

            printf("\n=================================================\n");
            printf(" NORMALIZATION DIAGNOSTIC REPORT\n");
            printf("=================================================\n");
            printf("  Max Absolute Error (Normalization Constraint |g_mu_nu p^mu p^nu|): %e (Ray %ld)\n", max_err_norm, worst_ray_norm);
            {free_pinned}(norm_diag_bridge); // Free the host-transfer array used for normalization checks.
        }} // END IF: commondata->perform_conservation_check to evaluate terminal normalization constraint

        // Free both sets of double-buffered arrays.
        for (int s = 0; s < 2; ++s) {{
            // Free host-transfer arrays.
            {free_pinned}(chunk_buffer[s]); // Free the ray-index host-transfer array.
            {free_pinned}(f_bridge[s]); // Free the state vector $f^\mu$ host-transfer array.
            {free_pinned}(f_p_bridge[s]); // Free the $f^\mu_{{n-1}}$ host-transfer array.
            {free_pinned}(f_p_p_bridge[s]); // Free the $f^\mu_{{n-2}}$ host-transfer array.
            {free_pinned}(affine_bridge[s]); // Free the affine-parameter host-transfer array.
            {free_pinned}(h_bridge[s]); // Free the step-size host-transfer array.
            {free_pinned}(status_bridge[s]); // Free the trajectory-status host-transfer array.
            {free_pinned}(retries_bridge[s]); // Free the step-rejection host-transfer array.
            {free_pinned}(on_pos_window_prev_bridge[s]); // Free the observer-window-side host-transfer array.
            {free_pinned}(on_pos_source_prev_bridge[s]); // Free the source-plane-side host-transfer array.
            {free_pinned}(affine_p_bridge[s]); // Free the $\lambda_{{n-1}}$ host-transfer array.
            {free_pinned}(affine_p_p_bridge[s]); // Free the $\lambda_{{n-2}}$ host-transfer array.
            {free_pinned}(window_event_found_bridge[s]); // Free the observer-window-intersection host-transfer array.
            {free_pinned}(source_event_found_bridge[s]); // Free the source-plane-intersection host-transfer array.

            // Free work arrays.
            {free_device}(d_f_bundle[s]); // Free the state vector $f^\mu$ array.
            {free_device}(d_f_start_bundle[s]); // Free the anchor state vector $f_{{start}}$ array.
            {free_device}(d_f_temp_bundle[s]); // Free the temporary stage $f^\mu_{{temp}}$ array.
            {free_device}(d_f_prev_bundle[s]); // Free the history state $f^\mu_{{n-1}}$ array.
            {free_device}(d_f_pre_prev_bundle[s]); // Free the history state $f^\mu_{{n-2}}$ array.
            {free_device}(d_metric_bundle[s]); // Free the symmetric metric tensor $g_{{\mu\nu}}$ array.
            {free_device}(d_connection_bundle[s]); // Free the Christoffel-symbol array.
            {free_device}(d_k_bundle[s]); // Free the RKF45 stage-derivative array.
            {free_device}(d_h[s]); // Free the integration-step array $h$.
            {free_device}(d_affine[s]); // Free the affine-parameter array $\lambda$.
            {free_device}(d_status[s]); // Free the trajectory-status array.
            {free_device}(d_retries[s]); // Free the step-rejection array.
            {free_device}(d_on_pos_window_prev[s]); // Free the previous window-side flags.
            {free_device}(d_on_pos_source_prev[s]); // Free the previous source-side flags.
            {free_device}(d_affine_prev[s]); // Free the historical affine-parameter array $\lambda_{{n-1}}$.
            {free_device}(d_affine_pre_prev[s]); // Free the historical affine-parameter array $\lambda_{{n-2}}$.
            {free_device}(d_window_event_found[s]); // Free the window-intersection flags.
            {free_device}(d_source_event_found[s]); // Free the source-intersection flags.
            {free_device}(d_chunk_buffer[s]); // Free the global ray-index array.
            {free_device}(d_norm_bundle[s]); // Free the diagnostic-output array.

            {stream_destroy}
        }} // END LOOP: for s over 2 to free double-buffered arrays

        //==========================================
        // CPU CONSERVATION DRIFT EVALUATION
        //==========================================
        // Evaluate relative numerical drift on the CPU.
        if (commondata->perform_conservation_check) {{
            // Calculate terminal conserved quantities on the selected CPU or GPU.
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

            // Loop over all rays to calculate errors on the host.
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
            }} // END LOOP: for i over num_rays to calculate errors on host

            printf("  Max Relative Error (Energy E): %e (Ray %ld)\n", max_err_E, worst_ray_E); // Output block printing the maximum relative error for energy $E$.
            printf("  Max Absolute Error (Energy E): %e (Ray %ld)\n\n", max_abs_err_E, worst_ray_abs_E); // Output block printing the maximum absolute error for energy $E$.

            printf("  Max Relative Error (Momentum Lz): %e (Ray %ld)\n", max_err_Lz, worst_ray_Lz); // Output block printing the maximum relative error for angular momentum $L_z$.
            printf("  Max Absolute Error (Momentum Lz): %e (Ray %ld)\n\n", max_abs_err_Lz, worst_ray_abs_Lz); // Output block printing the maximum absolute error for angular momentum $L_z$.

            printf("  Max Relative Error (Carter Q): %e (Ray %ld)\n", max_err_Q, worst_ray_Q); // Output block printing the maximum relative error for Carter constant $Q$.
            printf("  Max Absolute Error (Carter Q): %e (Ray %ld)\n", max_abs_err_Q, worst_ray_abs_Q); // Output block printing the maximum absolute error for Carter constant $Q$.

            printf("=================================================\n\n"); // Output block printing the terminal footer for the diagnostic sequence.

            // Free the host diagnostic arrays.
            {free_pinned}(initial_cq_host); // Free initial conserved quantities.
            {free_pinned}(final_cq_host); // Free final conserved quantities.
        }} // END IF: commondata->perform_conservation_check to evaluate numerical drift on the CPU

        // Free the host arrays holding ray states and affine parameters.
        {free_pinned}(all_photons_host.f); // Free the state vector $f^\mu$.
        {free_pinned}(all_photons_host.f_p); // Free the previous state $f^\mu_{{n-1}}$.
        {free_pinned}(all_photons_host.f_p_p); // Free the state from two steps earlier $f^\mu_{{n-2}}$.
        {free_pinned}(all_photons_host.affine_param); // Free the affine parameter $\lambda$.
        {free_pinned}(all_photons_host.h); // Free the integration step size $h$.
        {free_pinned}(all_photons_host.status); // Free the trajectory-status array.
        {free_pinned}(all_photons_host.rejection_retries); // Free the step-rejection array.
        {free_pinned}(all_photons_host.on_positive_side_of_window_prev); // Free the observer-window-side flags.
        {free_pinned}(all_photons_host.on_positive_side_of_source_prev); // Free the source-plane-side flags.
        {free_pinned}(all_photons_host.affine_param_p); // Free the historical affine parameter $\lambda_{{n-1}}$.
        {free_pinned}(all_photons_host.affine_param_p_p); // Free the historical affine parameter $\lambda_{{n-2}}$.
        {free_pinned}(all_photons_host.window_event_found); // Free the observer-window-intersection flags.
        {free_pinned}(all_photons_host.source_event_found); // Free the source-plane-intersection flags.

        // Free the intersection-result array $b_i$.
        {free_device}(d_results_buffer);

        // Free the time-slot arrays.
        slot_manager_free(&tsm);
    """

    cfc.register_CFunction(
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
