# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/rkf45_finalize_and_control_kernel.py
r"""
Defines the native kernel and host-side orchestrator for the RKF45 Finalize step.

This module provides the computational kernel responsible for calculating the final
5th-order solution, estimating the local truncation error, and executing the adaptive
step-size controller logic. The kernel bounds local memory usage by calculating the
5th-order candidate and error component-by-component for all 9 tensor components. The
L1 momentum floor is calculated using the initial state to prevent redundant reads
during the main loop, while scalar loads map directly to local memory. It calculates
physical state updates using exact double-precision coefficients, explicitly checking
for numerical singularities to reject corrupted trajectory steps. Additionally,
truncation errors are calculated directly via coefficient deltas to avoid catastrophic
floating-point cancellation against the anchor state. Upon step acceptance, the kernel
commits the updated state, advances the integration-parameter tracker, resets
the retry counter, and sets the ray status to active; otherwise, it increments
the retries and sets a rejected or failure status. Direct geodesic evolution
tracks affine parameter, while normalized evolution tracks coordinate time and
uses its log-energy tolerance in the embedded error norm.

When requested by `enable_numerical_time_window_step_cap`, the generated kernel
also caps accepted next-step sizes using `rkf45_max_delta_t` so backward
numerical-spacetime ray tracing remains inside the mapped time window.
The numerical photon example registers the companion time-window manager before
registering this kernel, so that the shared `rkf45_max_delta_t` parameter exists
when the accepted-step cap is emitted.

Together, the time-window manager and this kernel enforce one shared contract:
the manager maps enough lower-time numerical data for one photon slot assuming
the next accepted RKF45 step will not move farther backward in coordinate time
than the promised lookahead, and this kernel makes that assumption true by
limiting the next accepted integration-parameter step when necessary.

For the full mapped-window reasoning, including why the lower side of the
window carries extra backward-time lookahead while the upper side only needs
the centered interpolation halo, see
`time_window_manager_numerical_required_grid_range()` and
`time_window_manager_numerical_stencil_for_time()` in the companion
`time_window_manager_numerical` helper.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

from typing import List, Union

import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par


def rkf45_finalize_and_control_kernel(
    enable_numerical_time_window_step_cap: bool = False,
    normalized_eom: bool = False,
) -> None:
    r"""
    Global kernel for RKF45 finalization and error control.

    The kernel computes the 4th and 5th order solutions, calculates the error
    norm, and updates the photon's status (ACTIVE/REJECTED) and step size $h$
    in global memory.

    :param enable_numerical_time_window_step_cap: Whether to cap accepted RKF45
        step sizes so numerical-spacetime interpolation remains inside the
        currently mapped time window. The caller must register the numerical
        time-window manager separately before enabling this option.
    :param normalized_eom: Whether the nine-component state uses normalized
        photon variables with coordinate time as its integration parameter.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> import os
    >>> import tempfile
    >>> cfc.CFunction_dict.clear()
    >>> with tempfile.TemporaryDirectory(dir=os.getcwd()) as temp_dir:
    ...     old_cache_home = os.environ.get("XDG_CACHE_HOME")
    ...     _ = os.environ.__setitem__("XDG_CACHE_HOME", temp_dir)
    ...     rkf45_finalize_and_control_kernel(enable_numerical_time_window_step_cap=True)
    ...     generated = cfc.CFunction_dict["rkf45_finalize_and_control"].full_function
    ...     if old_cache_home is None:
    ...         _ = os.environ.pop("XDG_CACHE_HOME", None)
    ...     else:
    ...         _ = os.environ.__setitem__("XDG_CACHE_HOME", old_cache_home)
    >>> "rkf45_checked_floor_to_long" in generated
    True
    >>> "rkf45_checked_floor_to_long(slot_position, &slot_idx)" in generated
    True
    >>> cfc.CFunction_dict.clear()
    >>> with tempfile.TemporaryDirectory(dir=os.getcwd()) as temp_dir:
    ...     old_cache_home = os.environ.get("XDG_CACHE_HOME")
    ...     _ = os.environ.__setitem__("XDG_CACHE_HOME", temp_dir)
    ...     rkf45_finalize_and_control_kernel(
    ...         enable_numerical_time_window_step_cap=True, normalized_eom=True
    ...     )
    ...     generated = cfc.CFunction_dict["rkf45_finalize_and_control"].full_function
    ...     if old_cache_home is None:
    ...         _ = os.environ.pop("XDG_CACHE_HOME", None)
    ...     else:
    ...         _ = os.environ.__setitem__("XDG_CACHE_HOME", old_cache_home)
    >>> "rkf45_log_energy_tolerance" in generated
    True
    >>> "rkf45_constraint_tolerance" not in generated
    True
    """
    real_param_names: List[str] = [
        "rkf45_error_tolerance",
        "rkf45_absolute_error_tolerance",
        "rkf45_h_min",
        "rkf45_h_max",
    ]
    real_param_defaults: List[Union[str, int, float]] = [1e-8, 1e-8, 1e-10, 10.0]
    if normalized_eom:
        real_param_names.append("rkf45_log_energy_tolerance")
        real_param_defaults.append(1e-8)
    # The accepted-step cap is emitted only for numerical-spacetime builds.
    # The caller owns registration of the shared lookahead and slot-lattice
    # CodeParameters before this kernel consumes them.
    real_param_names.append("numerical_initial_h")
    real_param_defaults.append(1.0)
    par.register_CodeParameters(
        "REAL",
        __name__,
        real_param_names,
        real_param_defaults,
        commondata=True,
        add_to_parfile=True,
    )
    par.register_CodeParameter(
        "int", __name__, "rkf45_max_retries", 10, commondata=True, add_to_parfile=True
    )

    parallelization = par.parval_from_str("parallelization")
    cd_access = parallel_utils.get_commondata_access(parallelization)
    pragma_unroll = "#pragma unroll" if parallelization == "cuda" else ""

    arg_dict_cuda = {
        "d_f_persistent": "double *restrict",
        "d_f_start": "const double *restrict",
        "d_k_bundle": "const double *restrict",
        "d_h": "double *restrict",
        "d_status": "termination_type_t *restrict",
        "d_integration_param": "double *restrict",
        "d_retries": "int *restrict",
        "chunk_size": "const long int",
    }

    arg_dict_host = {
        "d_f_persistent": "double *restrict",
        "d_f_start": "const double *restrict",
        "d_k_bundle": "const double *restrict",
        "d_h": "double *restrict",
        "d_status": "termination_type_t *restrict",
        "d_integration_param": "double *restrict",
        "d_retries": "int *restrict",
        "chunk_size": "const long int",
    }
    # Pass commondata explicitly when not using CUDA's global memory
    if parallelization != "cuda":
        arg_dict_cuda["commondata"] = "const commondata_struct *restrict"
        arg_dict_host["commondata"] = "const commondata_struct *restrict"

    if parallelization == "cuda":
        loop_preamble = """
    //==========================================
    // CUDA THREAD IDENTIFICATION
    //==========================================
    // The identifier $i$ represents the global thread index mapped to a specific photon ray.
    const long int i = blockIdx.x * blockDim.x + threadIdx.x;

    // Guard prevents out-of-bounds VRAM access for threads exceeding the active chunk.
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

    prefunc = r"""
static inline int rkf45_checked_floor_to_long(
    const double value,
    long int *out) {
  if (out == NULL || !isfinite(value) ||
      value < (double)LONG_MIN || value > (double)LONG_MAX) {
    return 1;
  } // END IF: floor() source value was not representable as long int
  *out = (long int)floor(value);
  return 0;
} // END FUNCTION: rkf45_checked_floor_to_long
"""

    accepted_time_window_step_cap = ""
    proposed_time_window_step_cap = ""
    if enable_numerical_time_window_step_cap and normalized_eom:
        proposed_time_window_step_cap = rf"""
    if ({cd_access}rkf45_max_delta_t > 0.0 &&
        h_new_abs > {cd_access}rkf45_max_delta_t) {{
        h_new_abs = {cd_access}rkf45_max_delta_t;
        h_new = h_sign * h_new_abs;
    }} // END IF: next coordinate-time step exceeded rkf45_max_delta_t
"""
        accepted_time_window_step_cap = rf"""
        //==========================================
        // NUMERICAL TIME-WINDOW ACCEPTED-STEP CAP
        //==========================================
        {{
            const double t_accepted = ReadCUDA(&d_integration_param[i]);
            const double slot_position =
                (t_accepted - {cd_access}slot_manager_t_min) /
                {cd_access}slot_manager_delta_t;
            long int slot_idx = 0L;

            if (rkf45_checked_floor_to_long(slot_position, &slot_idx) != 0) {{
                h_new_abs = {cd_access}rkf45_h_min;
                h_new = h_sign * h_new_abs;
                WriteCUDA(&d_status[i], TERMINATION_TYPE_FAILURE);
            }} else {{
                const double slot_lower =
                    {cd_access}slot_manager_t_min +
                    (double)slot_idx * {cd_access}slot_manager_delta_t;
                const double delta_to_lower_slot_edge =
                    fmax(0.0, t_accepted - slot_lower);
                const double allowed_delta_t =
                    delta_to_lower_slot_edge + {cd_access}rkf45_max_delta_t;

                if (allowed_delta_t > 0.0 && h_new_abs > allowed_delta_t) {{
                    h_new_abs = allowed_delta_t;
                    h_new = h_sign * h_new_abs;
                }} // END IF: next coordinate-time step exceeded mapped-window lookahead
            }} // END ELSE: slot position was representable as long int
        }} // END BLOCK: normalized accepted-step time-window cap
"""
    elif enable_numerical_time_window_step_cap:
        accepted_time_window_step_cap = rf"""
        //==========================================
        // NUMERICAL TIME-WINDOW ACCEPTED-STEP CAP
        //==========================================
        // Reverse ray tracing moves toward lower coordinate time. The numerical
        // time-window manager mapped enough lower-time data under the
        // assumption that the next accepted RKF45 step will stay within the
        // current slot plus rkf45_max_delta_t of extra lookahead. Convert that
        // allowed coordinate-time motion into an affine-parameter cap using
        // the accepted fifth-order estimate of p^0 = dt/dlambda.
        // For the full cross-file contract, see
        // time_window_manager_numerical_required_grid_range() and
        // time_window_manager_numerical_stencil_for_time() in the companion
        // time_window_manager_numerical helper.
        {{
            const double t_accepted = f_5th_cache[0];
            // f_5th_cache[4] is used instead of the old p^0 because it is the
            // accepted fifth-order value after the just-completed step, hence
            // the best local estimate for converting the next affine step into
            // coordinate-time motion.
            const double p0_accepted_abs = fabs(f_5th_cache[4]);
            const double slot_position =
                (t_accepted - {cd_access}slot_manager_t_min) /
                {cd_access}slot_manager_delta_t;
            long int slot_idx = 0L;

            if (rkf45_checked_floor_to_long(slot_position, &slot_idx) != 0) {{
                h_new = {cd_access}rkf45_h_min;
                WriteCUDA(&d_status[i], TERMINATION_TYPE_FAILURE);
            }} else {{
                const double slot_lower =
                    {cd_access}slot_manager_t_min +
                    (double)slot_idx * {cd_access}slot_manager_delta_t;
                const double delta_to_lower_slot_edge =
                    fmax(0.0, t_accepted - slot_lower);
                const double allowed_delta_t =
                    delta_to_lower_slot_edge + {cd_access}rkf45_max_delta_t;

                if (p0_accepted_abs > 0.0 && allowed_delta_t > 0.0) {{
                    // Local estimate:
                    //   |Delta t| approx |p^0| Delta lambda
                    // so require:
                    //   Delta lambda <= 0.9 * allowed_delta_t / |p^0|
                    // Use a fixed 0.9 geometric safety margin here, independent
                    // of the RKF45 adaptive error-control safety factor, so the
                    // step stays away from the mmap boundary even if p^0 changes
                    // during the next trial step.
                    const double time_window_h_cap =
                        0.9 * allowed_delta_t / p0_accepted_abs;
                    if (time_window_h_cap < {cd_access}rkf45_h_min) {{
                        // This configuration should be unreachable in practice.
                        // Fail the ray rather than stepping outside the mapped
                        // numerical time window that the companion manager
                        // promised to keep available.
                        h_new = time_window_h_cap;
                        WriteCUDA(&d_status[i], TERMINATION_TYPE_FAILURE);
                    }} else {{
                        h_new = fmin(h_new, time_window_h_cap);
                        h_new = fmax(h_new, {cd_access}rkf45_h_min);
                    }} // END ELSE: time-window cap remained compatible with h_min
                }} // END IF: accepted temporal derivative supports a finite cap
            }} // END ELSE: slot position was representable as long int
        }} // END BLOCK: numerical time-window accepted-step cap
"""

    if normalized_eom:
        tolerance_load = rf"""
    const double log_energy_tol =
        {cd_access}rkf45_log_energy_tolerance;
"""
        error_normalization = r"""
        (void)p_L1;
        if (comp > 0 && comp < 8) {
            const double state_mag = fmax(AbsCUDA(f_n), AbsCUDA(f_5th_val));
            double scale = 0.0;

            if (comp <= 3) {
                scale = AddCUDA(atol, MulCUDA(rtol, state_mag));
            } else if (comp == 4) {
                scale = log_energy_tol;
            } else {
                scale = AddCUDA(atol, MulCUDA(rtol, fmax(1.0, state_mag)));
            } // END ELSE: normalized-state tolerance selection

            const double current_err = DivCUDA(err_abs, scale);
            if (isnan(current_err) || isinf(current_err)) {
                err_norm = 1e30;
            } else {
                err_norm = fmax(err_norm, current_err);
            } // END ELSE: normalized-error accumulation
        } // END IF: exclude integration parameter and Eulerian length from error norm
"""
        adaptive_step_control = rf"""
    const double h_sign = (h_local < 0.0) ? -1.0 : 1.0;
    double h_new_abs = MulCUDA(safety, MulCUDA(AbsCUDA(h_local), factor));
    h_new_abs = fmax(h_new_abs, {cd_access}rkf45_h_min);
    h_new_abs = fmin(h_new_abs, {cd_access}rkf45_h_max);
    double h_new = h_sign * h_new_abs;
"""
        integration_parameter_update = r"""
        // Normalized evolution advances coordinate time by the accepted step.
        const double old_integration_param = ReadCUDA(&d_integration_param[i]);
        WriteCUDA(
            &d_integration_param[i], AddCUDA(old_integration_param, h_local));
"""
    else:
        tolerance_load = ""
        error_normalization = r"""
        if (comp < 8) {
            double scale = 0.0;

            if (comp == 0) {
                scale = atol;
            } else if (comp <= 3) {
                scale = AddCUDA(atol, MulCUDA(rtol, AbsCUDA(f_n)));
            } else if (comp == 4) {
                scale = AddCUDA(atol, MulCUDA(rtol, AbsCUDA(f_n)));
            } else {
                scale = AddCUDA(atol, MulCUDA(rtol, p_L1));
            } // END ELSE: direct-geodesic tolerance selection

            const double current_err = DivCUDA(err_abs, scale);
            if (isnan(current_err) || isinf(current_err)) {
                err_norm = 1e30;
            } else {
                err_norm = fmax(err_norm, current_err);
            } // END ELSE: direct-geodesic error accumulation
        } // END IF: exclude Eulerian length from error norm
"""
        adaptive_step_control = rf"""
    double h_new = MulCUDA(safety, MulCUDA(h_local, factor));
    h_new = fmax(h_new, {cd_access}rkf45_h_min);
    h_new = fmin(h_new, {cd_access}rkf45_h_max);
"""
        integration_parameter_update = r"""
        // Direct geodesic evolution advances affine parameter by the accepted step.
        const double old_integration_param = ReadCUDA(&d_integration_param[i]);
        WriteCUDA(
            &d_integration_param[i], AddCUDA(old_integration_param, h_local));
"""

    core_math = rf"""
    //==========================================
    // MACRO DEFINITIONS FOR BUNDLE ACCESS
    //==========================================
    // IDX_F maps a component to the flattened state bundle using SoA layout.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
    // IDX_K maps a component to the flattened derivative bundle.
    #define IDX_K(s, c, ray_id) ((s) * 9 * BUNDLE_CAPACITY + (c) * BUNDLE_CAPACITY + (ray_id))

    //==========================================
    // PARAMETER LOAD
    //==========================================
    // Load tolerances and state variables from global memory and constant memory structs.
    const double rtol = {cd_access}rkf45_error_tolerance; // Relative tolerance bounds.
    const double atol = {cd_access}rkf45_absolute_error_tolerance; // Absolute tolerance bounds.
{tolerance_load}
    double h_local = ReadCUDA(&d_h[i]); // Local step size $h$.
    const int retries = ReadCUDA(&d_retries[i]); // Current count of failed adaptation attempts.

    //==========================================
    // COMPUTE & CACHE LOOP
    //==========================================
    // We compute the 5th order candidate and error component-by-component to bound register usage.
    double f_5th_cache[9]; // Thread-local register cache for the evaluated 5th-order state $f^{{\mu}}$.
    double err_norm = 0.0; // Accumulator for the maximum normalized truncation error norm $L_\infty$.

    //==========================================
    // L1 MOMENTUM FLOOR
    //==========================================
    // Evaluates the $L_1$ momentum floor using the initial state to avoid redundant memory reads inside the main loop.
    double p_L1 = 0.0; // Intermediate accumulator for the $L_1$ momentum floor.
    {{
        const double px = ReadCUDA(&d_f_start[IDX_F(5, i)]); // Spatial momentum $p_x$.
        const double py = ReadCUDA(&d_f_start[IDX_F(6, i)]); // Spatial momentum $p_y$.
        const double pz = ReadCUDA(&d_f_start[IDX_F(7, i)]); // Spatial momentum $p_z$.
        p_L1 = AddCUDA(AbsCUDA(px), AddCUDA(AbsCUDA(py), AbsCUDA(pz))); // The scalar $L_1$ momentum norm.
    }} // END BLOCK: L1 momentum floor calculation

    for (int comp = 0; comp < 9; ++comp) {{
        //==========================================
        // BASE STATE AND DERIVATIVE COMPONENT LOAD
        //==========================================
        // Scalar loads mapped directly to registers for the current tensor component.
        const double f_n = ReadCUDA(&d_f_start[IDX_F(comp, i)]);     // Base state tensor component $f_n$.
        const double k0  = ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]); // Stage 0 derivative vector $k_0$.
        // Stage 1 (k1) is mathematically zeroed out in RKF45, so it is skipped.
        const double k2  = ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]); // Stage 2 derivative vector $k_2$.
        const double k3  = ReadCUDA(&d_k_bundle[IDX_K(3, comp, i)]); // Stage 3 derivative vector $k_3$.
        const double k4  = ReadCUDA(&d_k_bundle[IDX_K(4, comp, i)]); // Stage 4 derivative vector $k_4$.
        const double k5  = ReadCUDA(&d_k_bundle[IDX_K(5, comp, i)]); // Stage 5 derivative vector $k_5$.

        //==========================================
        // 5TH ORDER CANDIDATE EVALUATION & STATE CORRUPTION SAFEGUARD
        //==========================================
        // Evaluates the physical state update $f^\mu_{{n+1}}$ utilizing exact double-precision Runge-Kutta coefficients to catch non-linear tensor interactions near event horizons before memory persistence.
        double update = MulCUDA(16.0 / 135.0, k0); // Intermediate accumulator for the 5th order step update.
        update = FusedMulAddCUDA(6656.0 / 12825.0, k2, update); // Accumulates the stage 2 derivative vector $k_2$.
        update = FusedMulAddCUDA(28561.0 / 56430.0, k3, update); // Accumulates the stage 3 derivative vector $k_3$.
        update = FusedMulAddCUDA(-9.0 / 50.0, k4, update); // Accumulates the stage 4 derivative vector $k_4$.
        update = FusedMulAddCUDA(2.0 / 55.0, k5, update); // Accumulates the stage 5 derivative vector $k_5$.

        const double f_5th_val = FusedMulAddCUDA(h_local, update, f_n); // The 5th order candidate state $f^\mu_{{n+1}}$.
        f_5th_cache[comp] = f_5th_val; // Caches the evaluated component to thread-local registers.

        // Evaluates the physical state for numerical singularities to guarantee the rejection of corrupted trajectory steps.
        if (isnan(f_5th_val) || isinf(f_5th_val)) {{
            err_norm = 1e30; // Forces an artificially massive error norm $L_\infty$ to guarantee step rejection.
        }} // END IF: numerical singularities check

        //==========================================
        // TRUNCATION ERROR EVALUATION
        //==========================================
        // Evaluates the truncation error directly via coefficient deltas ($C_5 - C_4$) to prevent catastrophic floating-point cancellation against the anchor state $f_n$.
        double err_val = MulCUDA(1.0 / 360.0, k0); // Intermediate accumulator for the truncation error.
        err_val = FusedMulAddCUDA(-128.0 / 4275.0, k2, err_val); // Accumulates the stage 2 error coefficient delta.
        err_val = FusedMulAddCUDA(-2197.0 / 75240.0, k3, err_val); // Accumulates the stage 3 error coefficient delta.
        err_val = FusedMulAddCUDA(1.0 / 50.0, k4, err_val); // Accumulates the stage 4 error coefficient delta.
        err_val = FusedMulAddCUDA(2.0 / 55.0, k5, err_val); // Accumulates the stage 5 error coefficient delta.

        const double err_abs = AbsCUDA(MulCUDA(h_local, err_val)); // The absolute magnitude of the local truncation error.

        //==========================================
        // ERROR NORMALIZATION & L-INFINITY NORM
        //==========================================
{error_normalization}
    }} // END LOOP: for comp over 9 tensor components

    //==========================================
    // UNIFIED ADAPTIVE CONTROL LOGIC
    //==========================================
    // Evaluates the mathematically optimal adaptive step size $h$ for subsequent integration.
    const double safety = 0.9; // Fixed RKF45 damping factor for next-step scaling.
    double factor = (err_norm > 1e-15) ? pow(DivCUDA(1.0, err_norm), 0.2) : 2.0; // Growth or shrink factor for the adaptive step $h$.

{adaptive_step_control}
{proposed_time_window_step_cap}

    if (err_norm <= 1.0) {{
        //==========================================
        // ACCEPTED STEP MEMORY COMMIT
        //==========================================
        // Commits the local register cache to persistent global memory.
        {pragma_unroll}
        for (int comp = 0; comp < 9; ++comp) {{
            WriteCUDA(&d_f_persistent[IDX_F(comp, i)], f_5th_cache[comp]); // Commits the cached state component to persistent memory.
        }} // END LOOP: for comp over 9 tensor components to commit state

        // Advance the affine-parameter or coordinate-time tracker.
{integration_parameter_update}

        // Reset Retries & Set Status.
        WriteCUDA(&d_retries[i], 0); // Resets the retry counter to zero for the active ray.
        WriteCUDA(&d_status[i], ACTIVE); // Updates the ray status flag to active.

{accepted_time_window_step_cap}
        // Commit the newly adapted step size $h$ to memory.
        WriteCUDA(&d_h[i], h_new); // Commits the newly adapted step size $h$ to memory.
    }} // END IF: err_norm <= 1.0
    else {{
        //==========================================
        // REJECTED STEP HANDLING
        //==========================================
        // Increment Retry Counter.
        const int new_retries = retries + 1; // Updated retry count for the current step.
        WriteCUDA(&d_retries[i], new_retries); // Commits the incremented retry count to memory.

        // Update Status for fatal unrecoverable errors vs recoverable rejections.
        if (new_retries > {cd_access}rkf45_max_retries) {{
            WriteCUDA(&d_status[i], FAILURE_RKF45_REJECTION_LIMIT); // Updates the ray status flag to a fatal rejection failure.
        }} else {{
            WriteCUDA(&d_status[i], REJECTED); // Updates the ray status flag to a recoverable rejection.
        }} // END ELSE: recoverable rejection status

        // Commit the scaled retry step size $h$ to memory without overwriting persistent state vectors.
        WriteCUDA(&d_h[i], h_new); // Commits the scaled retry step size $h$ to memory.
    }} // END ELSE: rejected step handling

    //==========================================
    // MACRO CLEANUP
    //==========================================
    // Undefine macros to ensure hermetic compilation and prevent redefinition errors.
    #undef IDX_F
    #undef IDX_K
    """

    kernel_body = f"{loop_preamble}\n{core_math}\n{loop_postamble}"

    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
        "stream": "stream_idx",
    }

    generated_prefunc, launch_code = parallel_utils.generate_kernel_and_launch_code(
        kernel_name="rkf45_finalize_and_control_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_host,
        parallelization=parallelization,
        launch_dict=launch_dict,
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
    )
    prefunc += generated_prefunc

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "limits.h"]
    if parallelization == "cuda":
        includes.append("cuda_intrinsics.h")

    desc = r""" Finalizes the RKF45 step, checks errors, and updates state/stepsize.

    @param d_f_persistent Pointer to the persistent nine-component state (updated on acceptance).
    @param d_f_start Pointer to the base nine-component state (read-only).
    @param d_k_bundle Pointer to all 6 derivative vectors $k_n$.
    @param d_h Pointer to the step size $h$.
    @param d_status Pointer to the ray status flag.
    @param d_integration_param Pointer to the affine-parameter or coordinate-time tracker.
    @param d_retries Pointer to the retry counter.
    @param chunk_size The number of rays in the current bundle.

    @note When numerical-spacetime support requests the optional accepted-step
    cap, this routine becomes the runtime enforcement layer for the
    slot-based numerical time window mapped by the companion manager.
    """

    cfunc_type = "void"

    name = "rkf45_finalize_and_control"

    params = (
        "const commondata_struct *restrict commondata, "
        "double *restrict d_f_persistent, "
        "const double *restrict d_f_start, "
        "const double *restrict d_k_bundle, "
        "double *restrict d_h, "
        "termination_type_t *restrict d_status, "
        "double *restrict d_integration_param, "
        "int *restrict d_retries, "
        "const long int chunk_size, "
        "const int stream_idx"
    )

    include_CodeParameters_h = False

    body = launch_code

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
