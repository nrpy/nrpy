# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_numerical.py
r"""
Emit the CPU/OpenMP numerical-spacetime photon batch integrator.

This module registers the C orchestrator for batched photon geodesics in a
numerical spacetime. The generated integrator consumes a validated combined
numerical raytracing ``.bin`` path stored on ``commondata``, initializes a
``NumericalTimeWindowManager``, maps numerical time windows per time slot, and
calls ``numerical_interpolation()`` for metric and Christoffel data.

The numerical time-window logic assumes the generator also registers
``rkf45_finalize_and_control_kernel(enable_numerical_time_window_step_cap=True)``
so accepted RKF45 steps remain inside the mapped numerical window.

The script keeps the broad photon orchestration, RKF45 stepping, and
``TimeSlotManager`` structure used by the analytic photon batch integrator, but
all metric/connection evaluations go through the numerical ``.bin`` data path.
It does not compute analytic conserved quantities.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.time_slot_manager_helpers import (
    time_slot_manager_helpers,
)


def batch_integrator_numerical(
    spacetime_name: str,
    dataset_coord_system: str,
) -> None:
    r"""
    Construct the CPU numerical-spacetime photon batch integrator.

    :param spacetime_name: Spacetime identifier used by photon-specific initial-condition helpers.
    :param dataset_coord_system: Coordinate system used by the numerical dataset.
    :raises ValueError: If the configured parallelization mode is not OpenMP.

    Doctests:
    >>> import os
    >>> import nrpy.c_function as cfc
    >>> os.environ["XDG_CACHE_HOME"] = "/tmp"
    >>> cfc.CFunction_dict.clear()
    >>> batch_integrator_numerical("Schwarzschild", "Spherical")
    >>> generated = cfc.CFunction_dict["batch_integrator_numerical"].full_function
    >>> "@param[in] commondata" in generated
    True
    >>> "time_window_manager_numerical_mmap_for_slot" in generated
    True
    """
    if "time_slot_manager" not in par.glb_extras_dict.get("BHaH_defines", {}):
        time_slot_manager_helpers()

    parallelization = par.parval_from_str("parallelization")
    if parallelization != "openmp":
        raise ValueError(
            "batch_integrator_numerical currently supports only "
            "parallelization='openmp'."
        )
    if dataset_coord_system in ("Spherical", "SinhSpherical"):
        phi_dim = 2
    elif dataset_coord_system in ("Cylindrical", "SinhCylindrical"):
        phi_dim = 1
    else:
        raise ValueError(
            "batch_integrator_numerical currently supports only "
            "dataset_coord_system in ('Spherical', 'SinhSpherical', "
            f"'Cylindrical', 'SinhCylindrical'); found '{dataset_coord_system}'."
        )

    # Supplied-name mapping: the checklist aliases numerical_time_window_manager_set_inert,
    # numerical_time_window_manager_init, numerical_time_window_manager_mmap_for_slot,
    # and numerical_time_window_manager_free correspond to the supplied helper names
    # time_window_manager_numerical_set_inert/init/mmap_for_slot/free.
    # Core physics and numerical simulation parameters for the global spacetime struct.
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "r_escape",
            "p_t_max",
            "numerical_initial_h",
        ],
        [150.0, 1e3, 0.1],
        commondata=True,
        add_to_parfile=True,
    )
    par.register_CodeParameter(
        "char[4096]",
        __name__,
        "numerical_spacetime_bin_path",
        "",
        commondata=True,
        add_to_parfile=True,
        description=(
            "Path to the validated combined numerical raytracing .bin file used by "
            "the numerical photon batch integrator."
        ),
    )
    par.register_CodeParameters(
        "bool",
        __name__,
        ["perform_normalization_check"],
        [False],
        commondata=True,
        add_to_parfile=True,
    )

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<stdio.h>",
        "<stdlib.h>",
        "<string.h>",
        "<time.h>",
    ]

    desc = r"""CPU numerical-spacetime photon batch integrator.

    This function bins active rays by coordinate time using TimeSlotManager,
    maps combined numerical-spacetime .bin time windows through
    NumericalTimeWindowManager, interpolates metric and Christoffels through
    numerical_interpolation(), advances photons with the existing RKF45 kernels,
    and writes final blueprint results.

    @param[in] commondata Struct containing global numerical-spacetime, camera,
                          and integration parameters.
    @param num_rays Total number of photon trajectories to simulate.
    @param[out] results_buffer Host array storing the final physical
                               intersections."""

    cfunc_type = "void"

    name = "batch_integrator_numerical"

    params = "const commondata_struct *restrict commondata, long int num_rays, blueprint_data_t *restrict results_buffer"

    include_CodeParameters_h = True

    malloc_pinned = "BHAH_MALLOC"
    malloc_device = "BHAH_MALLOC"
    free_pinned = "BHAH_FREE"
    free_device = "BHAH_FREE"
    pin_comment = "Host memory allocation: CPU RAM mapped for"
    bridge_alloc_comment = (
        "Allocate memory arrays in Host RAM for the structural bridge payloads."
    )
    scratch_alloc_comment = (
        "Allocate 1D Host scratchpad arrays for temporal data staging."
    )
    results_memcpy = (
        "// Event outputs are written directly to results_buffer on the CPU."
    )
    calc_blueprint = "calculate_and_fill_blueprint_data_universal(&all_photons_host, num_rays, results_buffer, 0);"
    set_initial_conditions_call = (
        f" set_initial_conditions_kernel_{spacetime_name}(commondata, num_rays, "
        "&all_photons_host, window_center_out, n_x_out, n_y_out, n_z_out);"
    )

    def memcpy_cpu(dest: str, src: str, size: str) -> str:
        return f"memcpy({dest}, {src}, {size});"

    def no_sync() -> str:
        return "// CPU memory operations above are complete at this point."

    stream_arg = ", 0"
    stream_arg_current = ", current"
    stream_arg_next = ", next"

    body = rf"""
    //==========================================
    // 1. HOST ALLOCATION
    //==========================================

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

    // CPU-only numerical integration: direct commondata access and no extra execution-buffer setup.

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

    //==========================================
    // DOUBLE-BUFFERED CPU SCRATCHPADS
    //==========================================
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
    }} // END LOOP: for s over 2 to instantiate the double-buffered operational arrays

    // Scratchpad array holding the terminal normalization diagnostic outputs.
    normalization_constraint_t *d_norm_bundle = NULL;

    if (commondata->perform_normalization_check) {{
        {malloc_device}(d_norm_bundle, sizeof(normalization_constraint_t) * BUNDLE_CAPACITY); // Allocate terminal normalization scratchpad.
    }} // END IF: commondata->perform_normalization_check to allocate normalization scratchpad

    // Event-detection kernels write final physical plane intersections directly to results_buffer.

    // Host-bound struct managing temporal binning of photon trajectories $x^\mu$.
    TimeSlotManager tsm;
    // The slot-manager upper bound is exclusive, so it must sit slightly above
    // t_start; exact equality would reject the initial photon state.
    const double slot_manager_t_max = commondata->t_start + 1.0e-5;
    // The main slot range is controlled by commondata, not by combined-file metadata.
    slot_manager_init(
        &tsm,
        commondata->slot_manager_t_min,
        slot_manager_t_max,
        commondata->slot_manager_delta_t,
        num_rays);

    // Numerical spacetime window manager and metadata-fed grid params.
    NumericalTimeWindowManager numerical_window;
    time_window_manager_numerical_set_inert(&numerical_window);

    // Seed the runtime params with generator-time defaults so coordinate-map
    // parameters like AMPL and SINHW remain valid after metadata overlay.
    commondata_struct commondata_for_params_defaults = *commondata;
    griddata_struct dummy_griddata[MAXNUMGRIDS];
    params_struct_set_to_default(&commondata_for_params_defaults, dummy_griddata);
    params_struct numerical_params = dummy_griddata[0].params;

    if (commondata->numerical_spacetime_bin_path[0] == '\0') {{
        fprintf(stderr,
                "ERROR: commondata->numerical_spacetime_bin_path is empty. "
                "The numerical photon batch integrator requires a validated combined .bin file path.\n");
        slot_manager_free(&tsm);
        exit(1);
    }} // END IF: numerical_spacetime_bin_path was empty

    if (time_window_manager_numerical_init(
            &numerical_window,
            commondata->numerical_spacetime_bin_path,
            commondata,
            commondata->numerical_spacetime_temporal_interp_order,
            &numerical_params) != TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
        fprintf(stderr,
                "ERROR: failed to initialize numerical time-window manager from '%s'.\n",
                commondata->numerical_spacetime_bin_path);
        slot_manager_free(&tsm);
        exit(1);
    }} // END IF: numerical time-window manager initialization failed

    if (numerical_params.Nxx{phi_dim} != 2) {{
        fprintf(stderr,
                "ERROR: numerical spatial interpolation expects exactly two stored phi planes in native dimension {phi_dim}; got Nxx{phi_dim}=%d.\n",
                numerical_params.Nxx{phi_dim});
        time_window_manager_numerical_free(&numerical_window);
        slot_manager_free(&tsm);
        exit(1);
    }} // END IF: stored phi-plane count was incompatible with azimuthal symmetry

    azimuthal_symmetry_spatial_lagrange_context_struct spatial_context;
    spatial_context.stored_phi_samples[0] =
        numerical_params.xxmin{phi_dim} + 0.5 * numerical_params.dxx{phi_dim};
    spatial_context.stored_phi_samples[1] =
        numerical_params.xxmin{phi_dim} + 1.5 * numerical_params.dxx{phi_dim};

    //==========================================
    // 2. INITIALIZATION PHASE
    //==========================================
    // Evaluate initial coordinate states and map global spacetime metrics to memory bounds.

    double window_center_out[3]; // 3D array storing the spatial Cartesian coordinates $x^i$ of the observer window center.
    double n_x_out[3]; // 3D orthonormal basis vector pointing along the $x$-axis of the local window geometry.
    double n_y_out[3]; // 3D orthonormal basis vector pointing along the $y$-axis of the local window geometry.
    double n_z_out[3]; // 3D orthonormal basis vector pointing along the $z$-axis of the local window geometry.

    // Operates synchronously as the primary state array must be fully populated before pipeline dispatch.
    {set_initial_conditions_call}

    //==========================================
    // DIAGNOSTIC PROBE: INITIAL-STATE ALIGNMENT
    //==========================================
    // Scan the master Host SoA immediately after the initial-condition kernel call.

    long int init_mismatch_count = 0; // Accumulator tracking the total number of physical state initializations that failed structural validation.

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

        if (fail_t || fail_x || fail_y || fail_z || fail_pt || fail_lam) {{
            init_mismatch_count++; // Increments the total accumulation of trajectory structural validation failures.
        }} // END IF: validate initialization coordinates
    }} // END LOOP: for p over num_rays to validate initialization

    if (init_mismatch_count > 0) {{
        const double mismatch_percent = ((double)init_mismatch_count / (double)num_rays) * 100.0; // Calculates the failure rate percentage for the initial state alignment.
        // This is a soft warning to surface initialization inconsistencies without halting execution.
        printf("[DIAGNOSTIC] Initialization Alignment Check: %ld out of %ld rays (%.2f%%) fail initial-state validation.\n", init_mismatch_count, num_rays, mismatch_percent);
    }} // END IF: init_mismatch_count > 0 to print diagnostic

    const int initial_slot_idx = slot_get_index(&tsm, commondata->t_start);
    if (initial_slot_idx < 0) {{
        fprintf(stderr,
                "ERROR: initial photon time t_start=%e is outside the configured TimeSlotManager range.\n",
                (double)commondata->t_start);
        time_window_manager_numerical_free(&numerical_window);
        slot_manager_free(&tsm);
        exit(1);
    }} // END IF: initial photon time was outside the configured slot range

    if (time_window_manager_numerical_mmap_for_slot(
            &numerical_window, &tsm, initial_slot_idx) !=
        TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
        fprintf(stderr,
                "ERROR: failed to map numerical time window for initial slot %d at t_start=%e.\n",
                initial_slot_idx,
                (double)commondata->t_start);
        time_window_manager_numerical_free(&numerical_window);
        slot_manager_free(&tsm);
        exit(1);
    }} // END IF: initial slot time window mapping failed

    long int num_batches = (num_rays + BUNDLE_CAPACITY - 1) / BUNDLE_CAPACITY; // Total integer calculation defining total iterative blocks required to process all photon indices.

    for (long int init_batch = 0; init_batch < num_batches; ++init_batch) {{ // Loop iterator $init_batch$ for evaluating the initialization constraint across sequential blocks.
        long int start_idx = init_batch * BUNDLE_CAPACITY; // Absolute starting index mapped to the master SoA for the current initialization batch.
        long int chunk_size = NRPYMIN((long int)BUNDLE_CAPACITY, num_rays - start_idx); // Dynamically sized operational boundary ensuring the active chunk does not exceed total trajectories.

        for (int init_i = 0; init_i < chunk_size; ++init_i) {{ // Loop index $init_i$ iterating over the specific initialization batch elements to pack the bridge.
            long int master_idx = start_idx + init_i; // Computes the absolute master index $m_{{idx}}$ tracking the photon within the global array.
            for (int init_k = 0; init_k < 9; ++init_k) {{ // Loop index $init_k$ iterating over the 9 tensor components of the state vector $f^\mu$.
                f_bridge[0][init_k * BUNDLE_CAPACITY + init_i] = all_photons_host.f[init_k * num_rays + master_idx]; // Assigns the active tensor state component to the primary bridge.
            }} // END LOOP: for init_k over 9 to iterate over tensor components
        }} // END LOOP: for init_i over chunk_size to pack the bridge

        for (int c_k = 0; c_k < 9; ++c_k) {{ // Loop index $c_k$ orchestrating the memory transfer of the 9 state vector $f^\mu$ components.
            {memcpy_cpu("d_f_bundle[0] + c_k * BUNDLE_CAPACITY", "f_bridge[0] + c_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size")}
        }} // END LOOP: for c_k over 9 to orchestrate memory transfer of state vector

        // Calculates symmetric metric tensor $g_{{\mu\nu}}$ strictly on the primary CPU buffer for the Hamiltonian constraint.
        numerical_interpolation(
            commondata,
            &numerical_params,
            &spatial_context,
            &numerical_window,
            d_f_bundle[0],
            d_metric_bundle[0],
            NULL,
            chunk_size,
            0);

        //==========================================
        // DIAGNOSTIC PROBE: METRIC INTEGRITY CHECK
        //==========================================
        // Extracts metric payload to confirm numerical stability prior to momentum solving.

        double *metric_diag_bridge; // Pointer storing temporary metric data to validate the interpolation sequence.
        // Memory allocation: Temporary bridge mapped to maximize throughput for the metric $g_{{\mu\nu}}$ diagnostic sequence.
        {malloc_pinned}(metric_diag_bridge, sizeof(double) * 10 * BUNDLE_CAPACITY);

        for (int m_k = 0; m_k < 10; ++m_k) {{
            // Loop index $m_k$ orchestrating memory transfer of the 10 metric tensor $g_{{\mu\nu}}$ components.
            {memcpy_cpu("metric_diag_bridge + m_k * BUNDLE_CAPACITY", "d_metric_bundle[0] + m_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size")}
        }} // END LOOP: for m_k over 10 to orchestrate memory transfer of metric tensor
        {no_sync()}

        long int metric_nan_count = 0; // Accumulator tracking the total number of metric tensor evaluations containing non-finite values.
        for (int m_diag_i = 0; m_diag_i < chunk_size; ++m_diag_i) {{ // Loop iterator $m_diag_i$ scanning each trajectory within the current initialization chunk.
            bool m_has_nan = false; // Boolean flag indicating if the specific metric tensor $g_{{\mu\nu}}$ contains a NaN or Inf value.
            for (int m_diag_k = 0; m_diag_k < 10; ++m_diag_k) {{ // Loop index $m_diag_k$ iterating over the 10 independent components of the symmetric metric tensor $g_{{\mu\nu}}$.
                if (isnan(metric_diag_bridge[m_diag_k * BUNDLE_CAPACITY + m_diag_i]) ||
                    isinf(metric_diag_bridge[m_diag_k * BUNDLE_CAPACITY + m_diag_i])) {{
                    m_has_nan = true; // Flags the trajectory metric state as invalid due to a non-finite value.
                    break; // Terminates the tensor component loop early to avoid unnecessary work upon detecting a failure.
                }} // END IF: check for NaN or Inf in metric
            }} // END LOOP: for m_diag_k over 10 to check metric tensor components
            if (m_has_nan) metric_nan_count++; // Increments the total accumulation of corrupted metric tensor evaluations.
        }} // END LOOP: for m_diag_i over chunk_size to scan for metric integrity

        if (metric_nan_count > 0) {{
            fprintf(stderr,
                    "ERROR: Init Batch %ld: %ld rays have invalid numerical metric "
                    "G_mu_nu before p_t solve.\n",
                    init_batch,
                    metric_nan_count);
            {free_pinned}(metric_diag_bridge);
            time_window_manager_numerical_free(&numerical_window);
            slot_manager_free(&tsm);
            exit(1);
        }} // END IF: metric_nan_count > 0 to fail loudly on invalid numerical metric
        // Memory Free: Purges the diagnostic bridge utilized for metric integrity checks.
        {free_pinned}(metric_diag_bridge);

        // Solves the constraint $p_\mu p^\mu = 0$ to find the temporal momentum $p_t$ natively on the active pipeline.
        p0_reverse_kernel(d_f_bundle[0], d_metric_bundle[0], chunk_size{stream_arg});

        for (int c_k = 0; c_k < 9; ++c_k) {{ // Loop index $c_k$ orchestrating memory transfer of the 9 constrained state vector $f^\mu$ components.
            {memcpy_cpu("f_bridge[0] + c_k * BUNDLE_CAPACITY", "d_f_bundle[0] + c_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size")}
        }} // END LOOP: for c_k over 9 to orchestrate memory transfer of constrained state vector
        {no_sync()}

        long int nonfinite_count = 0; // Accumulator tracking the total number of physical states $f^\mu$ containing non-finite values post-constraint solving.
        for (int gather_i = 0; gather_i < chunk_size; ++gather_i) {{ // Loop iterator $gather_i$ scanning each trajectory within the retrieved initialization chunk.
            long int master_idx = start_idx + gather_i; // Computes the absolute master index $m_{{idx}}$ mapping the localized chunk to the global master SoA.
            bool has_nonfinite = false; // Boolean flag indicating if the specific state vector $f^\mu$ contains a non-finite value.
            for (int gather_k = 0; gather_k < 9; ++gather_k) {{ // Loop index $gather_k$ iterating over the 9 tensor components of the state vector $f^\mu$.
                double val = f_bridge[0][gather_k * BUNDLE_CAPACITY + gather_i]; // Evaluates the updated numerical value of the specific tensor component.
                all_photons_host.f[gather_k * num_rays + master_idx] = val; // Maps the valid constrained tensor scalar back to the global Host SoA.
                if (!isfinite(val)) has_nonfinite = true; // Flags the physical state vector as invalid due to a non-finite evaluation.
            }} // END LOOP: for gather_k over 9 to iterate over tensor components
            if (has_nonfinite) nonfinite_count++; // Increments the total count of unresolved physical state vectors $f^\mu$.
        }} // END LOOP: for gather_i over chunk_size to retrieve updated constrained state vectors

        if (nonfinite_count > 0) {{
            fprintf(stderr,
                    "ERROR: Init Batch %ld: %ld rays contain nonfinite state values "
                    "after p_t solve. Aborting numerical batch integration.\n",
                    init_batch,
                    nonfinite_count);
            time_window_manager_numerical_free(&numerical_window);
            slot_manager_free(&tsm);
            exit(1);
        }} // END IF: nonfinite_count > 0 to abort on invalid post-p_t states
    }} // END LOOP: for init_batch over num_batches to evaluate initialization constraints


    long int sync_i; // Loop iterator index $sync_i$ spanning the entire global ray count to synchronize starting properties across history states.
    for(sync_i = 0; sync_i < num_rays; ++sync_i) {{
        int sync_k; // Loop index $sync_k$ iterating over the 9 tensor components to populate the historical derivatives $\dot{{f}}^\mu$ and $\ddot{{f}}^\mu$.
        for (sync_k = 0; sync_k < 9; ++sync_k) {{
            all_photons_host.f_p[sync_k * num_rays + sync_i] = all_photons_host.f[sync_k * num_rays + sync_i]; // Propagates the initial coordinate state vector $f^\mu$ to the first history derivative matrix.
            all_photons_host.f_p_p[sync_k * num_rays + sync_i] = all_photons_host.f[sync_k * num_rays + sync_i]; // Propagates the initial coordinate state vector $f^\mu$ to the second history derivative matrix.
        }} // END LOOP: for sync_k over 9 to propagate historical derivatives
        all_photons_host.status[sync_i] = ACTIVE; // Assigns the initial trajectory activity enum for the global physics engine.
        all_photons_host.affine_param[sync_i] = 0.0; // Sets the initial baseline progression scalar for the affine parameter $\lambda$.
        all_photons_host.rejection_retries[sync_i] = 0; // Clears the error rejection scalar to initialize the step size convergence tracking.

        all_photons_host.affine_param_p[sync_i] = 0.0; // Initializes the first historical affine parameter $\lambda_{{n-1}}$.
        all_photons_host.affine_param_p_p[sync_i] = 0.0; // Initializes the second historical affine parameter $\lambda_{{n-2}}$.
        all_photons_host.window_event_found[sync_i] = false; // Sets the observer window intersection logical lock to false.
        all_photons_host.source_event_found[sync_i] = false; // Sets the source emission intersection logical lock to false.

        int s_idx = slot_get_index(&tsm, all_photons_host.f[0 * num_rays + sync_i]); // Integer index $s_{{idx}}$ mapping the current photon's temporal coordinate $t$ to a discrete execution bin in the TimeSlotManager.
        if (s_idx != -1) {{
            slot_add_photon(&tsm, s_idx, sync_i); // Registers the active photon index to its corresponding temporal bin mapped by the orchestrator.
        }} // END IF: s_idx != -1 to add photon to slot
    }} // END LOOP: for sync_i over num_rays to synchronize starting properties

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
            break; // Terminates the temporal engine early to avoid unnecessary work.
        }} // END IF: total_active_photons <= 0 to break

        if (tsm.slot_counts[slot_idx] > 0) {{
            if (time_window_manager_numerical_mmap_for_slot(
                    &numerical_window, &tsm, slot_idx) !=
                TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
                fprintf(stderr,
                        "ERROR: failed to map numerical time window for slot %d [%e, %e).\n",
                        slot_idx,
                        (double)slot_lower_time(&tsm, slot_idx),
                        (double)slot_upper_time(&tsm, slot_idx));
                time_window_manager_numerical_free(&numerical_window);
                slot_manager_free(&tsm);
                exit(1);
            }} // END IF: slot time window mapping failed
        }} // END IF: slot_idx still contained active photons

        int current = 0; // Integer index tracking the primary active CPU buffer for execution.
        int next = 1; // Integer index tracking the secondary CPU buffer preparing the upcoming payload.
        long int active_chunks[2] = {{ 0, 0}}; // 1D array storing the total number of trajectories queued for each operational buffer.

        //==========================================
        // PHASE A: PRIME THE PUMP (Buffer 0)
        //==========================================
        // Populate the first bridge and launch synchronous integration on the primary buffer.

        active_chunks[current] = NRPYMIN((long int)BUNDLE_CAPACITY, tsm.slot_counts[slot_idx]); // Evaluates the active chunk size bounding the CPU staging to avoid CPU scratch overflow.

        if (active_chunks[current] > 0) {{
            slot_remove_chunk(&tsm, slot_idx, chunk_buffer[current], active_chunks[current]); // Extracts the execution chunk mapping from the Host-side temporal bin.

            // 1. Pack the 9-component tensors using cache-friendly forward sweeps
            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ iterating over the 9 tensor components of the state vectors.
                for (int bridge_i = 0; bridge_i < active_chunks[current]; ++bridge_i) {{  // Loop iterator $bridge_i$ packing the physical state payloads into the Host-side bridge arrays.
                    long int m_idx = chunk_buffer[current][bridge_i]; // Absolute master index $m_{{ idx}}$ mapping the active payload to the global trajectory matrix.
                    f_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f[c_k * num_rays + m_idx]; // Packs the coordinate state vector $f^\mu$ into the transfer bridge.
                    f_p_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p[c_k * num_rays + m_idx]; // Packs the first derivative $\dot{{ f}}^\mu$ into the transfer bridge.
                    f_p_p_bridge[current][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p_p[c_k * num_rays + m_idx]; // Packs the second derivative $\ddot{{ f}}^\mu$ into the transfer bridge.
                }} // END LOOP: for bridge_i over active_chunks[current] to pack payloads
            }} // END LOOP: for c_k over 9 to pack tensor components

            // 2. Pack the 1D arrays in a separate sequential loop
            for (int bridge_i = 0; bridge_i < active_chunks[current]; ++bridge_i) {{
                long int m_idx = chunk_buffer[current][bridge_i];
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
            }} // END LOOP: for bridge_i over active_chunks[current] to pack 1D arrays

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating CPU buffer copy of the 9 state vector components.
                // CPU buffer copy: Synchronously pushes bounded state vectors $f^\mu$ to CPU scratch strictly on buffer [current] to minimize latency.
                {memcpy_cpu("d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "f_bridge[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]")}
                // CPU buffer copy: Synchronously pushes first derivatives $\dot{{ f}}^\mu$ to CPU scratch strictly on buffer [current] to minimize latency.
                {memcpy_cpu("d_f_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "f_p_bridge[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]")}
                // CPU buffer copy: Synchronously pushes second derivatives $\ddot{{ f}}^\mu$ to CPU scratch strictly on buffer [current] to minimize latency.
                {memcpy_cpu("d_f_pre_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "f_p_p_bridge[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]")}
            }} // END LOOP: for c_k over 9 to orchestrate CPU buffer copy
            // CPU buffer copy: Synchronously pushes step sizes $h$ to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_h[current]", "h_bridge[current]", "sizeof(double) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes status enums to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_status[current]", "status_bridge[current]", "sizeof(termination_type_t) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes rejection scalars to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_retries[current]", "retries_bridge[current]", "sizeof(int) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes affine parameters $\lambda$ to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_affine[current]", "affine_bridge[current]", "sizeof(double) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes window boundary flags to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_on_pos_window_prev[current]", "on_pos_window_prev_bridge[current]", "sizeof(bool) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes source boundary flags to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_on_pos_source_prev[current]", "on_pos_source_prev_bridge[current]", "sizeof(bool) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes historical affine parameters $\lambda_{{ n-1}}$ to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_affine_prev[current]", "affine_p_bridge[current]", "sizeof(double) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes historical affine parameters $\lambda_{{ n-2}}$ to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_affine_pre_prev[current]", "affine_p_p_bridge[current]", "sizeof(double) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes window intersection locks to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_window_event_found[current]", "window_event_found_bridge[current]", "sizeof(bool) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes source intersection locks to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_source_event_found[current]", "source_event_found_bridge[current]", "sizeof(bool) * active_chunks[current]")}
            // CPU buffer copy: Synchronously pushes chunk indices $m_{{ idx}}$ to CPU scratch strictly on buffer [current] to minimize latency.
            {memcpy_cpu("d_chunk_buffer[current]", "chunk_buffer[current]", "sizeof(long int) * active_chunks[current]")}

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating CPU buffer baseline setup of the 9 state vector components.
                // CPU buffer copy: Duplicates the initial physical state vector $f^\mu$ to anchor the final RKF45 evaluation.
                {memcpy_cpu("d_f_start_bundle[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]")}
                // CPU buffer copy: Primes the temporary state vector bundle $f^\mu_{{ temp}}$ for iterative stage accumulation.
                {memcpy_cpu("d_f_temp_bundle[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]")}
            }} // END LOOP: for c_k over 9 to setup CPU buffer baseline

            for (int stage = 1; stage <= 6; ++stage) {{  // Loop iterator $stage$ executing the 6 discrete stages of the RKF45 Runge-Kutta numerical solver.
                // Interpolation step: evaluate the metric tensor $g_{{ \mu\nu}}$
                // and connection $\Gamma^\alpha_{{ \beta\gamma}}$ on the active buffer.
                numerical_interpolation(
                    commondata,
                    &numerical_params,
                    &spatial_context,
                    &numerical_window,
                    d_f_temp_bundle[current],
                    d_metric_bundle[current],
                    d_connection_bundle[current],
                    active_chunks[current],
                    current);
                long int bad_interp_current = 0;
                for (int interp_i = 0; interp_i < active_chunks[current]; ++interp_i) {{
                    bool bad_interp = false;
                    for (int interp_c = 0; interp_c < 10; ++interp_c) {{
                        const double val = d_metric_bundle[current][interp_c * BUNDLE_CAPACITY + interp_i];
                        if (!isfinite(val)) {{
                            bad_interp = true;
                            break;
                        }} // END IF: one metric component was not finite
                    }} // END LOOP: for interp_c over metric components
                    if (!bad_interp) {{
                        for (int interp_c = 0; interp_c < 40; ++interp_c) {{
                            const double val = d_connection_bundle[current][interp_c * BUNDLE_CAPACITY + interp_i];
                            if (!isfinite(val)) {{
                                bad_interp = true;
                                break;
                            }} // END IF: one connection component was not finite
                        }} // END LOOP: for interp_c over connection components
                    }} // END IF: metric components were finite before checking the connection
                    if (bad_interp)
                        bad_interp_current++;
                }} // END LOOP: for interp_i over active chunks on the active buffer
                if (bad_interp_current > 0) {{
                    fprintf(stderr,
                            "ERROR: Slot %d stage %d buffer %d: %ld rays had "
                            "nonfinite numerical interpolation output. "
                            "Aborting numerical batch integration.\n",
                            slot_idx,
                            stage,
                            current,
                            bad_interp_current);
                    time_window_manager_numerical_free(&numerical_window);
                    slot_manager_free(&tsm);
                    exit(1);
                }} // END IF: bad_interp_current > 0 to abort before RHS on invalid interpolation
                // RHS step: compute the geodesic equation derivatives $\dot{{ f}}^\mu$
                // on the active buffer.
                calculate_ode_rhs_kernel(d_f_temp_bundle[current], d_metric_bundle[current], d_connection_bundle[current], d_k_bundle[current], stage, active_chunks[current]{stream_arg_current});
                // Stage update: accumulate the intermediate RKF45 state updates on the
                // active buffer.
                rkf45_stage_update(d_f_start_bundle[current], d_k_bundle[current], d_h[current], stage, active_chunks[current], d_f_temp_bundle[current]{stream_arg_current});
            }} // END LOOP: for stage over 6 to execute RKF45 stages

            // Finalize step: apply Cash-Karp error control to update the integration
            // baseline and step size $h$.
            rkf45_finalize_and_control(commondata, d_f_bundle[current], d_f_start_bundle[current], d_k_bundle[current], d_h[current], d_status[current], d_affine[current], d_retries[current], active_chunks[current]{stream_arg_current});
            // Event step: detect geometric events and record intersection coordinate
            // states on the active buffer.
            event_detection_manager_kernel(commondata, d_f_bundle[current], d_f_prev_bundle[current], d_f_pre_prev_bundle[current], d_affine[current], d_affine_prev[current], d_affine_pre_prev[current], results_buffer, d_status[current], d_on_pos_window_prev[current], d_on_pos_source_prev[current], d_window_event_found[current], d_source_event_found[current], d_chunk_buffer[current], active_chunks[current]{stream_arg_current});

            for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating CPU buffer copy of the 9 state vector components.
                // CPU buffer copy: Retrieves updated coordinate states $f^\mu$ back to CPU RAM synchronously on the active buffer.
                {memcpy_cpu("f_bridge[current] + c_k * BUNDLE_CAPACITY", "d_f_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]")}
                // CPU buffer copy: Retrieves updated first derivatives $\dot{{ f}}^\mu$ back to CPU RAM synchronously on the active buffer.
                {memcpy_cpu("f_p_bridge[current] + c_k * BUNDLE_CAPACITY", "d_f_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]")}
                // CPU buffer copy: Retrieves updated second derivatives $\ddot{{ f}}^\mu$ back to CPU RAM synchronously on the active buffer.
                {memcpy_cpu("f_p_p_bridge[current] + c_k * BUNDLE_CAPACITY", "d_f_pre_prev_bundle[current] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[current]")}
            }} // END LOOP: for c_k over 9 to orchestrate CPU buffer copy
            // CPU buffer copy: Retrieves active step sizes $h$ back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("h_bridge[current]", "d_h[current]", "sizeof(double) * active_chunks[current]")}
            // CPU buffer copy: Retrieves updated status enums back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("status_bridge[current]", "d_status[current]", "sizeof(termination_type_t) * active_chunks[current]")}
            // CPU buffer copy: Retrieves active rejection counts back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("retries_bridge[current]", "d_retries[current]", "sizeof(int) * active_chunks[current]")}
            // CPU buffer copy: Retrieves total affine progression $\lambda$ back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("affine_bridge[current]", "d_affine[current]", "sizeof(double) * active_chunks[current]")}
            // CPU buffer copy: Retrieves updated window boundary flags back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("on_pos_window_prev_bridge[current]", "d_on_pos_window_prev[current]", "sizeof(bool) * active_chunks[current]")}
            // CPU buffer copy: Retrieves updated source boundary flags back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("on_pos_source_prev_bridge[current]", "d_on_pos_source_prev[current]", "sizeof(bool) * active_chunks[current]")}
            // CPU buffer copy: Retrieves historical affine parameter $\lambda_{{ n-1}}$ back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("affine_p_bridge[current]", "d_affine_prev[current]", "sizeof(double) * active_chunks[current]")}
            // CPU buffer copy: Retrieves historical affine parameter $\lambda_{{ n-2}}$ back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("affine_p_p_bridge[current]", "d_affine_pre_prev[current]", "sizeof(double) * active_chunks[current]")}
            // CPU buffer copy: Retrieves active window locks back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("window_event_found_bridge[current]", "d_window_event_found[current]", "sizeof(bool) * active_chunks[current]")}
            // CPU buffer copy: Retrieves active source locks back to CPU RAM synchronously on the active buffer.
            {memcpy_cpu("source_event_found_bridge[current]", "d_source_event_found[current]", "sizeof(bool) * active_chunks[current]")}
        }} // END IF: active_chunks[current] > 0 to prime the pump on Buffer 0

        //==========================================
        // PHASE B: THE OVERLAP LOOP
        //==========================================
        // Continuously alternate between buffers, packing the next payload while syncing the current.

        while (active_chunks[current] > 0 || tsm.slot_counts[slot_idx] > 0) {{

            active_chunks[next] = NRPYMIN((long int)BUNDLE_CAPACITY, tsm.slot_counts[slot_idx]); // Evaluates the active chunk size bounding the upcoming CPU staging execution block.
            if (active_chunks[next] > 0) {{
                slot_remove_chunk(&tsm, slot_idx, chunk_buffer[next], active_chunks[next]); // Extracts the next execution chunk mapping from the Host-side temporal bin.

                // 1. Pack the 9-component tensors using cache-friendly forward sweeps
                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ iterating over the 9 tensor components of the state vectors.
                    for (int bridge_i = 0; bridge_i < active_chunks[next]; ++bridge_i) {{  // Loop iterator $bridge_i$ packing the physical state payloads into the next Host-side bridge array.
                        long int m_idx = chunk_buffer[next][bridge_i]; // Absolute master index $m_{{ idx}}$ mapping the active payload to the global trajectory matrix.
                        f_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f[c_k * num_rays + m_idx]; // Packs the coordinate state vector $f^\mu$ into the transfer bridge.
                        f_p_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p[c_k * num_rays + m_idx]; // Packs the first derivative $\dot{{ f}}^\mu$ into the transfer bridge.
                        f_p_p_bridge[next][c_k * BUNDLE_CAPACITY + bridge_i] = all_photons_host.f_p_p[c_k * num_rays + m_idx]; // Packs the second derivative $\ddot{{ f}}^\mu$ into the transfer bridge.
                    }} // END LOOP: for bridge_i over active_chunks[next] to pack payloads
                }} // END LOOP: for c_k over 9 to pack tensor components

                // 2. Pack the 1D arrays in a separate sequential loop
                for (int bridge_i = 0; bridge_i < active_chunks[next]; ++bridge_i) {{
                    long int m_idx = chunk_buffer[next][bridge_i];
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
                }} // END LOOP: for bridge_i over active_chunks[next] to pack 1D arrays

                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating CPU buffer copy of the 9 state vector components for the upcoming payload.
                    // CPU buffer copy: Synchronously pushes bounded state vectors $f^\mu$ to CPU scratch strictly on buffer [next] to overlap execution.
                    {memcpy_cpu("d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "f_bridge[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]")}
                    // CPU buffer copy: Synchronously pushes first derivatives $\dot{{ f}}^\mu$ to CPU scratch strictly on buffer [next] to overlap execution.
                    {memcpy_cpu("d_f_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "f_p_bridge[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]")}
                    // CPU buffer copy: Synchronously pushes second derivatives $\ddot{{ f}}^\mu$ to CPU scratch strictly on buffer [next] to overlap execution.
                    {memcpy_cpu("d_f_pre_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "f_p_p_bridge[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]")}
                }} // END LOOP: for c_k over 9 to orchestrate CPU buffer copy
                // CPU buffer copy: Synchronously pushes step sizes $h$ to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_h[next]", "h_bridge[next]", "sizeof(double) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes status enums to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_status[next]", "status_bridge[next]", "sizeof(termination_type_t) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes rejection scalars to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_retries[next]", "retries_bridge[next]", "sizeof(int) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes affine parameters $\lambda$ to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_affine[next]", "affine_bridge[next]", "sizeof(double) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes window boundary flags to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_on_pos_window_prev[next]", "on_pos_window_prev_bridge[next]", "sizeof(bool) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes source boundary flags to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_on_pos_source_prev[next]", "on_pos_source_prev_bridge[next]", "sizeof(bool) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes historical affine parameters $\lambda_{{ n-1}}$ to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_affine_prev[next]", "affine_p_bridge[next]", "sizeof(double) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes historical affine parameters $\lambda_{{ n-2}}$ to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_affine_pre_prev[next]", "affine_p_p_bridge[next]", "sizeof(double) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes window intersection locks to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_window_event_found[next]", "window_event_found_bridge[next]", "sizeof(bool) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes source intersection locks to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_source_event_found[next]", "source_event_found_bridge[next]", "sizeof(bool) * active_chunks[next]")}
                // CPU buffer copy: Synchronously pushes chunk indices $m_{{ idx}}$ to CPU scratch strictly on buffer [next] to overlap execution.
                {memcpy_cpu("d_chunk_buffer[next]", "chunk_buffer[next]", "sizeof(long int) * active_chunks[next]")}

                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating CPU buffer baseline setup of the 9 state vector components for the upcoming payload.
                    // CPU buffer copy: Duplicates the initial physical state vector $f^\mu$ to anchor the upcoming RKF45 evaluation.
                    {memcpy_cpu("d_f_start_bundle[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]")}
                    // CPU buffer copy: Primes the temporary state vector bundle $f^\mu_{{ temp}}$ for the upcoming iterative stage accumulation.
                    {memcpy_cpu("d_f_temp_bundle[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]")}
                }} // END LOOP: for c_k over 9 to setup CPU buffer baseline

                for (int stage = 1; stage <= 6; ++stage) {{  // Loop iterator $stage$ executing the 6 discrete stages of the upcoming RKF45 Runge-Kutta numerical solver.
                    // Interpolation step: evaluate the metric tensor
                    // $g_{{ \mu\nu}}$ and connection $\Gamma^\alpha_{{ \beta\gamma}}$
                    // on the alternate buffer.
                    numerical_interpolation(
                        commondata,
                        &numerical_params,
                        &spatial_context,
                        &numerical_window,
                        d_f_temp_bundle[next],
                        d_metric_bundle[next],
                        d_connection_bundle[next],
                        active_chunks[next],
                        next);
                    long int bad_interp_next = 0;
                    for (int interp_i = 0; interp_i < active_chunks[next]; ++interp_i) {{
                        bool bad_interp = false;
                        for (int interp_c = 0; interp_c < 10; ++interp_c) {{
                            const double val = d_metric_bundle[next][interp_c * BUNDLE_CAPACITY + interp_i];
                            if (!isfinite(val)) {{
                                bad_interp = true;
                                break;
                            }} // END IF: one metric component was not finite
                        }} // END LOOP: for interp_c over metric components
                        if (!bad_interp) {{
                            for (int interp_c = 0; interp_c < 40; ++interp_c) {{
                                const double val = d_connection_bundle[next][interp_c * BUNDLE_CAPACITY + interp_i];
                                if (!isfinite(val)) {{
                                    bad_interp = true;
                                    break;
                                }} // END IF: one connection component was not finite
                            }} // END LOOP: for interp_c over connection components
                        }} // END IF: metric components were finite before checking the connection
                        if (bad_interp)
                            bad_interp_next++;
                    }} // END LOOP: for interp_i over active chunks on the alternate buffer
                    if (bad_interp_next > 0) {{
                        fprintf(stderr,
                                "ERROR: Slot %d stage %d buffer %d: %ld rays had "
                                "nonfinite numerical interpolation output. "
                                "Aborting numerical batch integration.\n",
                                slot_idx,
                                stage,
                                next,
                                bad_interp_next);
                        time_window_manager_numerical_free(&numerical_window);
                        slot_manager_free(&tsm);
                        exit(1);
                    }} // END IF: bad_interp_next > 0 to abort before RHS on invalid interpolation
                    // RHS step: compute the geodesic equation derivatives
                    // $\dot{{ f}}^\mu$ on the alternate buffer.
                    calculate_ode_rhs_kernel(d_f_temp_bundle[next], d_metric_bundle[next], d_connection_bundle[next], d_k_bundle[next], stage, active_chunks[next]{stream_arg_next});
                    // Stage update: accumulate the intermediate RKF45 state updates on
                    // the alternate buffer.
                    rkf45_stage_update(d_f_start_bundle[next], d_k_bundle[next], d_h[next], stage, active_chunks[next], d_f_temp_bundle[next]{stream_arg_next});
                }} // END LOOP: for stage over 6 to execute RKF45 stages

                // Finalize step: apply Cash-Karp error control to update the upcoming
                // integration baseline and step size $h$.
                rkf45_finalize_and_control(commondata, d_f_bundle[next], d_f_start_bundle[next], d_k_bundle[next], d_h[next], d_status[next], d_affine[next], d_retries[next], active_chunks[next]{stream_arg_next});
                // Event step: detect geometric events and record intersection
                // coordinate states on the alternate buffer.
                event_detection_manager_kernel(commondata, d_f_bundle[next], d_f_prev_bundle[next], d_f_pre_prev_bundle[next], d_affine[next], d_affine_prev[next], d_affine_pre_prev[next], results_buffer, d_status[next], d_on_pos_window_prev[next], d_on_pos_source_prev[next], d_window_event_found[next], d_source_event_found[next], d_chunk_buffer[next], active_chunks[next]{stream_arg_next});
                for (int c_k = 0; c_k < 9; ++c_k) {{  // Loop index $c_k$ orchestrating CPU buffer copy of the 9 upcoming state vector components.
                    // CPU buffer copy: Retrieves updated coordinate states $f^\mu$ back to CPU RAM synchronously on the alternate buffer.
                    {memcpy_cpu("f_bridge[next] + c_k * BUNDLE_CAPACITY", "d_f_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]")}
                    // CPU buffer copy: Retrieves updated first derivatives $\dot{{ f}}^\mu$ back to CPU RAM synchronously on the alternate buffer.
                    {memcpy_cpu("f_p_bridge[next] + c_k * BUNDLE_CAPACITY", "d_f_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]")}
                    // CPU buffer copy: Retrieves updated second derivatives $\ddot{{ f}}^\mu$ back to CPU RAM synchronously on the alternate buffer.
                    {memcpy_cpu("f_p_p_bridge[next] + c_k * BUNDLE_CAPACITY", "d_f_pre_prev_bundle[next] + c_k * BUNDLE_CAPACITY", "sizeof(double) * active_chunks[next]")}
                }} // END LOOP: for c_k over 9 to orchestrate CPU buffer copy
                // CPU buffer copy: Retrieves upcoming active step sizes $h$ back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("h_bridge[next]", "d_h[next]", "sizeof(double) * active_chunks[next]")}
                // CPU buffer copy: Retrieves upcoming updated status enums back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("status_bridge[next]", "d_status[next]", "sizeof(termination_type_t) * active_chunks[next]")}
                // CPU buffer copy: Retrieves upcoming active rejection counts back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("retries_bridge[next]", "d_retries[next]", "sizeof(int) * active_chunks[next]")}
                // CPU buffer copy: Retrieves upcoming total affine progression $\lambda$ back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("affine_bridge[next]", "d_affine[next]", "sizeof(double) * active_chunks[next]")}
                // CPU buffer copy: Retrieves upcoming updated window boundary flags back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("on_pos_window_prev_bridge[next]", "d_on_pos_window_prev[next]", "sizeof(bool) * active_chunks[next]")}
                // CPU buffer copy: Retrieves upcoming updated source boundary flags back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("on_pos_source_prev_bridge[next]", "d_on_pos_source_prev[next]", "sizeof(bool) * active_chunks[next]")}
                // CPU buffer copy: Retrieves upcoming historical affine parameter $\lambda_{{ n-1}}$ back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("affine_p_bridge[next]", "d_affine_prev[next]", "sizeof(double) * active_chunks[next]")}
                // CPU buffer copy: Retrieves upcoming historical affine parameter $\lambda_{{ n-2}}$ back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("affine_p_p_bridge[next]", "d_affine_pre_prev[next]", "sizeof(double) * active_chunks[next]")}
                // CPU buffer copy: Retrieves upcoming active window locks back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("window_event_found_bridge[next]", "d_window_event_found[next]", "sizeof(bool) * active_chunks[next]")}
                // CPU buffer copy: Retrieves upcoming active source locks back to CPU RAM synchronously on the alternate buffer.
                {memcpy_cpu("source_event_found_bridge[next]", "d_source_event_found[next]", "sizeof(bool) * active_chunks[next]")}
            }} // END IF: active_chunks[next] > 0 to process upcoming payload

            if (active_chunks[current] > 0) {{
                // CPU memory operations are complete before payload unpacking.
                {no_sync()}

                // 1. Unpack 9-component tensors sequentially
                for (int fin_k = 0; fin_k < 9; ++fin_k) {{  // Loop index $fin_k$ retrieving the 9 tensor components back into the global Host matrix.
                    for (int fin_i = 0; fin_i < active_chunks[current]; ++fin_i) {{  // Loop iterator $fin_i$ unpacking the finalized physical data back to the global Host matrix.
                        long int m_idx = chunk_buffer[current][fin_i]; // Absolute master index $m_{{ idx}}$ retrieving the specific photon index from the execution chunk.
                        all_photons_host.f[fin_k * num_rays + m_idx] = f_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i]; // Unpacks the synchronized state vector $f^\mu$ into the global Host matrix.
                        all_photons_host.f_p[fin_k * num_rays + m_idx] = f_p_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i]; // Unpacks the synchronized first derivative $\dot{{ f}}^\mu$ into the global Host matrix.
                        all_photons_host.f_p_p[fin_k * num_rays + m_idx] = f_p_p_bridge[current][fin_k * BUNDLE_CAPACITY + fin_i]; // Unpacks the synchronized second derivative $\ddot{{ f}}^\mu$ into the global Host matrix.
                    }} // END LOOP: for fin_i over active_chunks[current] to unpack finalized data
                }} // END LOOP: for fin_k over 9 to retrieve tensor components

                // 2. Unpack 1D arrays sequentially
                for (int fin_i = 0; fin_i < active_chunks[current]; ++fin_i) {{
                    long int m_idx = chunk_buffer[current][fin_i];
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
                }} // END LOOP: for fin_i over active_chunks[current] to unpack 1D arrays

                // 3. TimeSlotManager State Update (Cache-hot, strictly sequential)
                for (int fin_i = 0; fin_i < active_chunks[current]; ++fin_i) {{
                    long int m_idx = chunk_buffer[current][fin_i];
                    if (status_bridge[current][fin_i] == ACTIVE) {{  // Evaluates the continuation logic if the trajectory remains within safe physical bounds.
                        int next_s_idx = slot_get_index(&tsm, all_photons_host.f[0 * num_rays + m_idx]); // Evaluates the updated temporal coordinate $t$ to determine the next operational bin.
                        if (next_s_idx != -1) {{  // Confirms the physical state has not exceeded the maximum simulation time bounds.
                            slot_add_photon(&tsm, next_s_idx, m_idx);  // Re-queues the updated physical state vector back into the host orchestrator.
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
            }} // END IF: active_chunks[current] > 0 to complete buffer and unpack

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
            }} // END LOOP: for bar_i over bar_width to construct loading bar
            bar[bar_width] = '\0'; // Terminates the loading bar character array to prevent buffer overruns.

            // Accumulator tracking the total number of adaptive step size $h$ rejections in the active chunk.
            long int batch_rejections = 0;

            // Loop iterator $sum_i$ scanning the finalized physical state bridge for error tolerance failures.
            for (int sum_i = 0; sum_i < active_chunks[current]; ++sum_i) {{
                batch_rejections += retries_bridge[current][sum_i]; // Accumulates the localized step rejection tally.
            }} // END LOOP: for sum_i over active_chunks[current] to calculate rejection count

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
            int temp = current; // Temporary integer scalar storing the primary buffer index for logical pointer swapping.
            current = next; // Shifts the primary execution tracker to the alternate buffer index.
            next = temp; // Assigns the cleared buffer index back to the upcoming payload queue.
        }} // END WHILE: alternating buffers to process temporal bin

        //==========================================
        // PHASE C: THE TIME BARRIER
        //==========================================
        // CPU execution is synchronous here; no additional extra barrier is needed.

     }} // END LOOP: for slot_idx down to 0 to process all temporal bins

    //==========================================
    // 4. RESULT WRITING, CLEANUP, AND FINALIZATION
    //==========================================
    // Process terminal photon trajectories and extract final geometric intersections.

        // CPU buffer copy: Extracts validated CPU-side blueprints $b_i$ containing geometric plane intersections.
        {results_memcpy}

        // Final output step: process escaped photons intersecting the celestial
        // sphere $r > r_{{escape}}$.
        {calc_blueprint}

        //==========================================
        // TERMINAL NORMALIZATION DIAGNOSTIC
        //==========================================
        if (commondata->perform_normalization_check) {{
            TimeSlotManager norm_tsm;
            int normalization_failure_mode = 0;
            long int normalization_failure_ray = -1;
            int normalization_failure_slot = -1;
            slot_manager_init(
                &norm_tsm,
                commondata->slot_manager_t_min,
                slot_manager_t_max,
                commondata->slot_manager_delta_t,
                num_rays);

            double max_err_norm = 0.0;
            long int worst_ray_norm = -1;

            for (long int norm_ray = 0; norm_ray < num_rays; ++norm_ray) {{
                const int norm_slot_idx = slot_get_index(
                    &norm_tsm, all_photons_host.f[0 * num_rays + norm_ray]);
                if (norm_slot_idx < 0) {{
                    continue;
                }} // END IF: norm_slot_idx < 0 to skip terminal states outside the numerical time domain
                slot_add_photon(&norm_tsm, norm_slot_idx, norm_ray);
            }} // END LOOP: for norm_ray over num_rays to populate the terminal slot manager

            for (int norm_slot_idx = norm_tsm.num_slots - 1;
                 norm_slot_idx >= 0 && normalization_failure_mode == 0;
                 --norm_slot_idx) {{
                if (norm_tsm.slot_counts[norm_slot_idx] <= 0) {{
                    continue;
                }} // END IF: norm_tsm.slot_counts[norm_slot_idx] <= 0 to skip an empty terminal slot

                if (time_window_manager_numerical_mmap_for_slot(
                        &numerical_window, &norm_tsm, norm_slot_idx) !=
                    TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
                    normalization_failure_mode = 2;
                    normalization_failure_slot = norm_slot_idx;
                    break;
                }} // END IF: terminal normalization mmap fails

                while (norm_tsm.slot_counts[norm_slot_idx] > 0 &&
                       normalization_failure_mode == 0) {{
                    const long int chunk_size = NRPYMIN(
                        (long int)BUNDLE_CAPACITY, norm_tsm.slot_counts[norm_slot_idx]);
                    slot_remove_chunk(
                        &norm_tsm, norm_slot_idx, chunk_buffer[0], chunk_size);

                    for (int norm_k = 0; norm_k < 9; ++norm_k) {{
                        for (long int norm_i = 0; norm_i < chunk_size; ++norm_i) {{
                            const long int master_idx = chunk_buffer[0][norm_i];
                            f_bridge[0][norm_k * BUNDLE_CAPACITY + norm_i] =
                                all_photons_host.f[norm_k * num_rays + master_idx];
                        }} // END LOOP: for norm_i over chunk_size to pack the terminal state bundle
                    }} // END LOOP: for norm_k over 9 to pack the terminal state components

                    for (int norm_k = 0; norm_k < 9; ++norm_k) {{
                        {memcpy_cpu("d_f_bundle[0] + norm_k * BUNDLE_CAPACITY", "f_bridge[0] + norm_k * BUNDLE_CAPACITY", "sizeof(double) * chunk_size")}
                    }} // END LOOP: for norm_k over 9 to copy the terminal state bundle to the CPU scratchpad

                    numerical_interpolation(
                        commondata,
                        &numerical_params,
                        &spatial_context,
                        &numerical_window,
                        d_f_bundle[0],
                        d_metric_bundle[0],
                        NULL,
                        chunk_size,
                        0);

                    normalization_constraint_photon(
                        d_f_bundle[0],
                        d_metric_bundle[0],
                        d_norm_bundle,
                        chunk_size,
                        0);

                    for (long int norm_i = 0; norm_i < chunk_size; ++norm_i) {{
                        const double current_norm_err = fabs(d_norm_bundle[norm_i].C);
                        if (!isfinite(current_norm_err)) {{
                            normalization_failure_mode = 3;
                            normalization_failure_ray = chunk_buffer[0][norm_i];
                            break;
                        }} // END IF: current_norm_err is non-finite
                        if (current_norm_err > max_err_norm) {{
                            max_err_norm = current_norm_err;
                            worst_ray_norm = chunk_buffer[0][norm_i];
                        }} // END IF: current_norm_err > max_err_norm
                    }} // END LOOP: for norm_i over chunk_size to scan constraint values
                }} // END WHILE: norm_tsm.slot_counts[norm_slot_idx] > 0 to process the terminal slot
            }} // END LOOP: for norm_slot_idx down to 0 to process all terminal slots

            slot_manager_free(&norm_tsm);

            if (normalization_failure_mode != 0) {{
                if (normalization_failure_mode == 1) {{
                    fprintf(stderr,
                            "ERROR: terminal normalization diagnostic could not bin photon %ld.\n",
                            normalization_failure_ray);
                }} else if (normalization_failure_mode == 2) {{
                    fprintf(stderr,
                            "ERROR: failed to map numerical time window for terminal normalization slot %d.\n",
                            normalization_failure_slot);
                }} else if (normalization_failure_mode == 3) {{
                    fprintf(stderr,
                            "ERROR: terminal normalization diagnostic produced a non-finite constraint for photon %ld.\n",
                            normalization_failure_ray);
                }} // END ELSE IF: normalization_failure_mode == 3 to report a non-finite terminal constraint
                time_window_manager_numerical_free(&numerical_window);
                slot_manager_free(&tsm);
                exit(1);
            }} // END IF: normalization_failure_mode != 0 to abort terminal normalization diagnostics

            printf("\n=================================================\n");
            printf(" NORMALIZATION DIAGNOSTIC REPORT\n");
            printf("=================================================\n");
            printf(
                "  Max Absolute Error (Normalization Constraint |g_mu_nu p^mu p^nu|): %e (Ray %ld)\n",
                max_err_norm,
                worst_ray_norm);
        }} // END IF: commondata->perform_normalization_check to evaluate terminal normalization constraint

        if (commondata->perform_normalization_check) {{
            {free_device}(d_norm_bundle); // Purges the terminal normalization diagnostic scratchpad.
        }} // END IF: commondata->perform_normalization_check to purge normalization scratchpad

        // Loop iterator $s$ purging the double-buffered arrays across both CPU buffers.
        for (int s = 0; s < 2; ++s) {{
            // Host Memory Free: Purges bridge components supporting scatter logic mapped to CPU memory transfers.
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

            // Host Memory Free: Purges remaining CPU scratch operational pipeline scratchpads.
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
        }} // END LOOP: for s over 2 to purge double-buffered arrays


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

        // Release the active numerical spacetime window before the slot lattice is destroyed.
        time_window_manager_numerical_free(&numerical_window);

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
