r"""
Provides the native kernel and host-side orchestrator for the RKF45 Finalize & Control step (Kernel 6).

Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
This module provides the global memory kernel responsible for computing the final 5th-order
solution, estimating the local truncation error, and executing the adaptive step-size
controller logic on VRAM bundles. To prevent register spilling on the NVIDIA RTX 3080 (which has 
a hard limit of 255 registers per thread), we strictly enforce a Split-Pipeline architecture, 
reading and writing intermediate stages via global VRAM scratchpads rather than fused registers.
Author: Dalton J. Moone.
"""

import nrpy.params as par
import nrpy.c_function as cfc
from nrpy.helpers.parallelization.utilities import (
    generate_kernel_and_launch_code,
    get_commondata_access,
)

def rkf45_finalize_and_control_kernel() -> None:
    r"""
    This function provides the global kernel for RKF45 finalization and error control.

    The kernel computes the 4th and 5th order solutions, calculates the error norm,
    and updates the photon's status (ACTIVE/REJECTED) and step size $h$ in global memory.

    :raises TypeError: If incorrect parameters are passed to the code generation functions.
    """
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "rkf45_error_tolerance",
            "rkf45_absolute_error_tolerance",
            "rkf45_h_min",
            "rkf45_h_max",
            "rkf45_safety_factor",
            "numerical_initial_h",
        ],
        [1e-8, 1e-8, 1e-10, 10.0, 0.9, 1.0],
        commondata=True,
        add_to_parfile=True,
    )
    par.register_CodeParameter(
        "int", __name__, "rkf45_max_retries", 10, commondata=True, add_to_parfile=True
    )

    parallelization = par.parval_from_str("parallelization")
    cd_access = get_commondata_access(parallelization)

    arg_dict_cuda = {
        "d_f_persistent": "double *restrict",
        "d_f_start": "const double *restrict",
        "d_k_bundle": "const double *restrict",
        "d_h": "double *restrict",
        "d_status": "termination_type_t *restrict",
        "d_affine": "double *restrict",
        "d_retries": "int *restrict",
        "chunk_size": "const long int"
    }

    arg_dict_host = {
        "d_f_persistent": "double *restrict",
        "d_f_start": "const double *restrict",
        "d_k_bundle": "const double *restrict",
        "d_h": "double *restrict",
        "d_status": "termination_type_t *restrict",
        "d_affine": "double *restrict",
        "d_retries": "int *restrict",
        "chunk_size": "const long int"
    }

    if parallelization == "cuda":
        loop_preamble = """
    // --- CUDA THREAD IDENTIFICATION ---
    // The identifier $i$ represents the global thread index mapped to a specific photon ray.
    const long int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Guard prevents out-of-bounds VRAM access for threads exceeding the active chunk.
    if (i >= chunk_size) return;
    """
        loop_postamble = ""
    else:
        loop_preamble = """
    // --- OPENMP LOOP ARCHITECTURE ---
    // Distribute photon rays across available CPU threads for parallel evaluation.
    #pragma omp parallel for
    for(long int i = 0; i < chunk_size; i++) {
    """
        loop_postamble = "    } // End OpenMP loop"

    core_math = rf"""
    // --- MACRO DEFINITIONS FOR BUNDLE ACCESS ---
    // IDX_F maps a component to the flattened state bundle using SoA layout.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
    // IDX_K maps a component to the flattened derivative bundle.
    #define IDX_K(s, c, ray_id) ((s) * 9 * BUNDLE_CAPACITY + (c) * BUNDLE_CAPACITY + (ray_id))

    // --- PARAMETER LOAD ---
    // Load tolerances and state variables from global memory and constant memory structs.
    const double rtol = {cd_access}rkf45_error_tolerance; // Relative tolerance bounds.
    const double atol = {cd_access}rkf45_absolute_error_tolerance; // Absolute tolerance bounds.
    double h_local = ReadCUDA(&d_h[i]); // Local step size $h$.
    const int retries = ReadCUDA(&d_retries[i]); // Current count of failed adaptation attempts.

    // --- COMPUTE & CACHE LOOP ---
    // We compute the 5th order candidate and error component-by-component to bound register usage.
    double f_5th_cache[9]; // Thread-local register cache for the evaluated 5th-order state $f^{{\mu}}$.
    double err_norm = 0.0; // Accumulator for the maximum normalized truncation error norm $L_\infty$.

    // --- L1 MOMENTUM FLOOR ---
    // Evaluates the $L_1$ momentum floor using the initial state to avoid redundant memory reads inside the main loop.
    double p_L1 = 0.0; // Intermediate accumulator for the $L_1$ momentum floor.
    {{
        const double px = ReadCUDA(&d_f_start[IDX_F(5, i)]); // Spatial momentum $p_x$.
        const double py = ReadCUDA(&d_f_start[IDX_F(6, i)]); // Spatial momentum $p_y$.
        const double pz = ReadCUDA(&d_f_start[IDX_F(7, i)]); // Spatial momentum $p_z$.
        p_L1 = AddCUDA(AbsCUDA(px), AddCUDA(AbsCUDA(py), AbsCUDA(pz))); // The scalar $L_1$ momentum norm.
    }}

    for (int comp = 0; comp < 9; ++comp) {{
        // --- BASE STATE AND DERIVATIVE COMPONENT LOAD ---
        // Scalar loads mapped directly to registers for the current tensor component.
        const double f_n = ReadCUDA(&d_f_start[IDX_F(comp, i)]);     // Base state tensor component $f_n$.
        const double k0  = ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]); // Stage 0 derivative vector $k_0$.
        const double k1  = ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]); // Stage 1 derivative vector $k_1$.
        const double k2  = ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]); // Stage 2 derivative vector $k_2$.
        const double k3  = ReadCUDA(&d_k_bundle[IDX_K(3, comp, i)]); // Stage 3 derivative vector $k_3$.
        const double k4  = ReadCUDA(&d_k_bundle[IDX_K(4, comp, i)]); // Stage 4 derivative vector $k_4$.
        const double k5  = ReadCUDA(&d_k_bundle[IDX_K(5, comp, i)]); // Stage 5 derivative vector $k_5$.

        // --- 5TH ORDER CANDIDATE EVALUATION & STATE CORRUPTION SAFEGUARD ---
        /* Evaluates the physical state update $f^\mu_{{n+1}}$ utilizing exact double-precision Runge-Kutta coefficients to catch non-linear tensor interactions near event horizons before memory persistence. */
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
        }}

        // --- TRUNCATION ERROR EVALUATION --- 
        /* Evaluates the truncation error directly via coefficient deltas ($C_5 - C_4$) to prevent catastrophic floating-point cancellation against the anchor state $f_n$. */
        double err_val = MulCUDA(1.0 / 360.0, k0); // Intermediate accumulator for the truncation error.
        err_val = FusedMulAddCUDA(-128.0 / 4275.0, k2, err_val); // Accumulates the stage 2 error coefficient delta.
        err_val = FusedMulAddCUDA(-2197.0 / 75240.0, k3, err_val); // Accumulates the stage 3 error coefficient delta.
        err_val = FusedMulAddCUDA(1.0 / 50.0, k4, err_val); // Accumulates the stage 4 error coefficient delta.
        err_val = FusedMulAddCUDA(2.0 / 55.0, k5, err_val); // Accumulates the stage 5 error coefficient delta.

        const double err_abs = AbsCUDA(MulCUDA(h_local, err_val)); // The absolute magnitude of the local truncation error.

        // --- ERROR NORMALIZATION & L-INFINITY NORM ---
        // Evaluates the normalized error for each tensor component to dictate the RKF45 adaptive step size $h$.
        if (comp < 8) {{ // Excludes the affine parameter $\lambda$ (index 8) from the error check.
            double scale = 0.0; // The bounded tolerance scale limit for the current tensor component.

            if (comp == 0) {{ 
                // Coordinate Time $t$: Applies pure absolute tolerance to prevent secular drift.
                scale = atol; // Applies pure absolute tolerance bounds.
            }} else if (comp <= 3) {{ 
                // Spatial Position $x^i$: Mixed tolerance scaling bounds spatial position variation.
                scale = AddCUDA(atol, MulCUDA(rtol, AbsCUDA(f_n))); // Mixed tolerance scaling bounds positional variation.
            }} else if (comp == 4) {{ 
                // Temporal Momentum $p_t$: Mixed tolerance scaling enforces energy conservation bounds.
                scale = AddCUDA(atol, MulCUDA(rtol, AbsCUDA(f_n))); // Mixed tolerance scaling bounds temporal momentum variation.
            }} else {{ 
                // Spatial Momentum $p_i$: Mixed tolerance bounded by the $L_1$ momentum floor.
                scale = AddCUDA(atol, MulCUDA(rtol, p_L1)); // Mixed tolerance scaling bounds spatial momentum variation.
            }}

            double current_err = DivCUDA(err_abs, scale); // The normalized error for the specific tensor component.
            
            // Evaluates the calculated error norm for numerical singularities to guard against IEEE-754 fmax behavior.
            if (isnan(current_err) || isinf(current_err)) {{
                err_norm = 1e30; // Forces an artificially massive error norm $L_\infty$ to guarantee step rejection.
            }} else {{
                err_norm = fmax(err_norm, current_err); // Accumulates the maximum normalized error equivalent to the $L_\infty$ norm.
            }}
        }}
    }}

    // --- UNIFIED ADAPTIVE CONTROL LOGIC ---
    // Evaluates the mathematically optimal adaptive step size $h$ for subsequent integration.
    double safety = {cd_access}rkf45_safety_factor; // Safety factor for step-size scaling.
    double factor = (err_norm > 1e-15) ? pow(DivCUDA(1.0, err_norm), 0.2) : 2.0; // Growth or shrink factor for the adaptive step $h$.
    
    double h_new = MulCUDA(safety, MulCUDA(h_local, factor)); // The candidate new step size $h$.
    h_new = fmax(h_new, {cd_access}rkf45_h_min); // Enforces the minimum step size $h_{{min}}$ bound.
    h_new = fmin(h_new, {cd_access}rkf45_h_max); // Enforces the maximum step size $h_{{max}}$ bound.

    if (err_norm <= 1.0) {{
        // --- ACCEPTED STEP MEMORY COMMIT ---
        // Commits the local register cache to persistent global memory.
        #pragma unroll
        for (int comp = 0; comp < 9; ++comp) {{
            WriteCUDA(&d_f_persistent[IDX_F(comp, i)], f_5th_cache[comp]); // Commits the cached state component to persistent memory.
        }}

        // Update Affine Parameter $\lambda$.
        const double old_affine = ReadCUDA(&d_affine[i]); // Previous affine parameter $\lambda$.
        WriteCUDA(&d_affine[i], AddCUDA(old_affine, h_local)); // Commits the updated affine parameter $\lambda$ to memory.

        // Reset Retries & Set Status.
        WriteCUDA(&d_retries[i], 0); // Resets the retry counter to zero for the active ray.
        WriteCUDA(&d_status[i], ACTIVE); // Updates the ray status flag to active.
        
        // Commit the newly adapted step size $h$ to memory.
        WriteCUDA(&d_h[i], h_new); // Commits the newly adapted step size $h$ to memory.
    }} else {{
        // --- REJECTED STEP HANDLING ---
        // Increment Retry Counter.
        const int new_retries = retries + 1; // Updated retry count for the current step.
        WriteCUDA(&d_retries[i], new_retries); // Commits the incremented retry count to memory.

        // Update Status for fatal unrecoverable errors vs recoverable rejections.
        if (new_retries > {cd_access}rkf45_max_retries) {{
            WriteCUDA(&d_status[i], FAILURE_RKF45_REJECTION_LIMIT); // Updates the ray status flag to a fatal rejection failure.
        }} else {{
            WriteCUDA(&d_status[i], REJECTED); // Updates the ray status flag to a recoverable rejection.
        }}
        
        // Commit the scaled retry step size $h$ to memory without overwriting persistent state vectors.
        WriteCUDA(&d_h[i], h_new); // Commits the scaled retry step size $h$ to memory.
    }}

    // --- MACRO CLEANUP ---
    // Undefine macros to ensure hermetic compilation and prevent redefinition errors.
    #undef IDX_F
    #undef IDX_K
    """

    kernel_body = f"{loop_preamble}\n{core_math}\n{loop_postamble}"

    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
        "stream": "stream_idx"
    }

    prefunc, launch_code = generate_kernel_and_launch_code(
        kernel_name="rkf45_finalize_and_control_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_host,
        parallelization=parallelization,
        launch_dict=launch_dict,
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "cuda_intrinsics.h", "<math.h>"]
    
    desc = r"""@brief Finalizes the RKF45 step, checks errors, and updates state/stepsize.
    
    @param d_f_persistent Pointer to the persistent state $f^{\mu}$ (updated on acceptance).
    @param d_f_start Pointer to the base state $f^{\mu}$ (read-only).
    @param d_k_bundle Pointer to all 6 derivative vectors $k_n$.
    @param d_h Pointer to the step size $h$.
    @param d_status Pointer to the ray status flag.
    @param d_affine Pointer to the affine parameter $\lambda$.
    @param d_retries Pointer to the retry counter.
    @param chunk_size The number of rays in the current bundle.
    """

    cfunc_type = "void"
    
    name = "rkf45_finalize_and_control"
    
    params = (
        "double *restrict d_f_persistent, "
        "const double *restrict d_f_start, "
        "const double *restrict d_k_bundle, "
        "double *restrict d_h, "
        "termination_type_t *restrict d_status, "
        "double *restrict d_affine, "
        "int *restrict d_retries, "
        "const long int chunk_size, "
        "const int stream_idx"
    )

    params = (
        "double *restrict d_f_persistent, "
        "const double *restrict d_f_start, "
        "const double *restrict d_k_bundle, "
        "double *restrict d_h, "
        "termination_type_t *restrict d_status, "
        "double *restrict d_affine, "
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