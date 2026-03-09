"""
Generates the native CUDA kernel and host-side orchestrator for the RKF45 Finalize & Control step (Kernel 6).

Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
This module provides the global memory kernel responsible for computing the final 5th-order
solution, estimating the local truncation error, and executing the adaptive step-size
controller logic on VRAM bundles. It implements a "Split-Pipeline" architecture to minimize
register pressure by operating on persistent VRAM scratchpads rather than fused registers.

Author: Dalton J. Moone.
"""

import nrpy.params as par
import nrpy.c_function as cfc
from nrpy.helpers.parallelization.utilities import (
    generate_kernel_and_launch_code,
    get_commondata_access,
)

def rkf45_finalize_and_control_kernel() -> None:
    """
    Register the global CUDA kernel for RKF45 finalization and error control.

    Generates the C function `rkf45_finalize_and_control_kernel` and its host launcher.
    The kernel computes the 4th and 5th order solutions, calculates the error norm,
    and updates the photon's status (ACTIVE/REJECTED) and step size $h$ in global VRAM.

    :raises TypeError: If incorrect parameters are passed to the code generation functions.
    """

    # Python: Registers RKF45 control parameters to CodeParameters for global access.
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

    # Python: Retrieve the accessor for constant memory (e.g., "d_commondata.").
    cd_access = get_commondata_access("cuda")

    # Python: Define the argument dictionary for the CUDA kernel generation utility.
    arg_dict_cuda = {
        "d_f_persistent": "double *restrict",        # Read/Write: Persistent state bundle
        "d_f_start": "const double *restrict",       # Input: Base state at t_n
        "d_k_bundle": "const double *restrict",      # Input: All 6 k derivative vectors
        "d_h": "double *restrict",                   # Read/Write: Adaptive step size
        "d_status": "termination_type_t *restrict",  # Write: Ray status flags
        "d_affine": "double *restrict",              # Read/Write: Affine parameter lambda
        "d_retries": "int *restrict",                # Read/Write: Retry counter
        "chunk_size": "const long int"        # Input: Active bundle size
    }

    # Python: Define the GPU kernel body using raw strings.
    # Note: We strictly avoid declaring double k[6][9] to prevent Register Bleeding.
    kernel_body = rf"""
    // --- THREAD IDENTIFICATION ---
    // Map the global thread index to the unique ray identifier within the bundle.
    const long int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Architectural Guard: Prevent out-of-bounds VRAM access if grid size > chunk_size.
    if (i >= chunk_size) return;

    // --- MACRO DEFINITIONS FOR BUNDLE ACCESS ---
    // IDX_F: Maps (component, ray) to the flattened State-of-Arrays (SoA) layout.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
    // IDX_K: Maps (stage, component, ray) to the flattened derivative bundle.
    #define IDX_K(s, c, ray_id) ((s) * 9 * BUNDLE_CAPACITY + (c) * BUNDLE_CAPACITY + (ray_id))

    // --- PARAMETER LOAD ---
    // Load tolerances and state variables from Global VRAM and Constant Memory.
    const double rtol = {cd_access}rkf45_error_tolerance;
    const double atol = {cd_access}rkf45_absolute_error_tolerance;
    double h_local = ReadCUDA(&d_h[i]);
    const int retries = ReadCUDA(&d_retries[i]);

    // --- COMPUTE & CACHE LOOP ---
    // We compute the 5th order candidate and error component-by-component.
    // We cache the result in 'f_5th_cache' (18 registers) instead of storing all K-vectors (108 registers).
    double f_5th_cache[9];
    double err_norm = 0.0;

    // --- L1 MOMENTUM FLOOR ---
    // Evaluates the L1 momentum floor using the initial state to avoid re-reading inside the loop.
    double p_L1 = 0.0;
    {{
        const double px = ReadCUDA(&d_f_start[IDX_F(5, i)]);
        const double py = ReadCUDA(&d_f_start[IDX_F(6, i)]);
        const double pz = ReadCUDA(&d_f_start[IDX_F(7, i)]);
        p_L1 = AddCUDA(AbsCUDA(px), AddCUDA(AbsCUDA(py), AbsCUDA(pz)));
    }}

    for (int comp = 0; comp < 9; ++comp) {{
        // 1. Load Base State & Derivative Components
        // Scalar loads mapped directly to registers for the current tensor component.
        const double f_n = ReadCUDA(&d_f_start[IDX_F(comp, i)]);
        const double k0  = ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]);
        const double k1  = ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]);
        const double k2  = ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]);
        const double k3  = ReadCUDA(&d_k_bundle[IDX_K(3, comp, i)]);
        const double k4  = ReadCUDA(&d_k_bundle[IDX_K(4, comp, i)]);
        const double k5  = ReadCUDA(&d_k_bundle[IDX_K(5, comp, i)]);

        // --- 5TH ORDER CANDIDATE EVALUATION & STATE CORRUPTION SAFEGUARD ---
        // Evaluates the physical state update $f^\mu_{{n+1}}$ utilizing exact double-precision Runge-Kutta coefficients.
        // This architectural step occurs here to catch non-linear tensor interactions near event horizons before VRAM persistence.
        double update = MulCUDA(16.0 / 135.0, k0); // Intermediate accumulator for the 5th order step update.
        update = FusedMulAddCUDA(6656.0 / 12825.0, k2, update);
        update = FusedMulAddCUDA(28561.0 / 56430.0, k3, update);
        update = FusedMulAddCUDA(-9.0 / 50.0, k4, update);
        update = FusedMulAddCUDA(2.0 / 55.0, k5, update);

        const double f_5th_val = FusedMulAddCUDA(h_local, update, f_n); // The 5th order candidate state $f^\mu_{{n+1}}$.
        f_5th_cache[comp] = f_5th_val; // Caches the evaluated component to thread-local registers.

        // Evaluates the physical state for numerical singularities to guarantee the rejection of corrupted trajectory steps.
        if (isnan(f_5th_val) || isinf(f_5th_val)) {{
        err_norm = 1e30; // Forces an artificially massive error norm $L_\infty$ to guarantee step rejection.
        }}
        // 3. Compute Truncation Error 
        // Evaluates the truncation error directly via coefficient deltas ($C_5 - C_4$) to prevent 
        // catastrophic floating-point cancellation against the anchor state $f_n$.
        // Delta Coeffs: (16/135 - 25/216), 0, (6656/12825 - 1408/2565), (28561/56430 - 2197/4104), (-9/50 - -1/5), 2/55
        double err_val = MulCUDA(1.0 / 360.0, k0);
        err_val = FusedMulAddCUDA(-128.0 / 4275.0, k2, err_val);
        err_val = FusedMulAddCUDA(-2197.0 / 75240.0, k3, err_val);
        err_val = FusedMulAddCUDA(1.0 / 50.0, k4, err_val);
        err_val = FusedMulAddCUDA(2.0 / 55.0, k5, err_val);

        const double err_abs = AbsCUDA(MulCUDA(h_local, err_val));

        // --- ERROR NORMALIZATION & L-INFINITY NORM ---
        // Evaluates the normalized error for each tensor component to dictate the RKF45 adaptive step size $h$.
        if (comp < 8) {{ // Excludes the affine parameter (index 8) from the error check.
        double scale = 0.0;

        if (comp == 0) {{ 
            // Coordinate Time $t$: Applies pure absolute tolerance to prevent secular drift.
            scale = atol; 
        }} else if (comp <= 3) {{ 
            // Spatial Position $x^i$: Mixed tolerance scaling bounds spatial position variation.
            scale = AddCUDA(atol, MulCUDA(rtol, AbsCUDA(f_n)));
        }} else if (comp == 4) {{ 
            // Temporal Momentum $p_t$: Mixed tolerance scaling enforces energy conservation bounds.
            scale = AddCUDA(atol, MulCUDA(rtol, AbsCUDA(f_n)));
        }} else {{ 
            // Spatial Momentum $p_i$: Mixed tolerance bounded by the $L_1$ momentum floor.
            scale = AddCUDA(atol, MulCUDA(rtol, p_L1));
        }}

        // Accumulates the maximum normalized error equivalent to the $L_\infty$ norm.
        double current_err = DivCUDA(err_abs, scale);
        
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
    double safety = d_commondata.rkf45_safety_factor;
    double factor = (err_norm > 1e-15) ? pow(DivCUDA(1.0, err_norm), 0.2) : 2.0;
    
    double h_new = MulCUDA(safety, MulCUDA(h_local, factor));
    h_new = fmax(h_new, d_commondata.rkf45_h_min);
    h_new = fmin(h_new, d_commondata.rkf45_h_max);

    if (err_norm <= 1.0) {{
        // === ACCEPTED STEP ===
        // Flush Cached State to Persistent VRAM.
        // This coalesced write commits the local register cache to global memory.
        #pragma unroll
        for (int comp = 0; comp < 9; ++comp) {{
        WriteCUDA(&d_f_persistent[IDX_F(comp, i)], f_5th_cache[comp]);
        }}

        // Update Affine Parameter $\lambda$.
        const double old_affine = ReadCUDA(&d_affine[i]);
        WriteCUDA(&d_affine[i], AddCUDA(old_affine, h_local));

        // Reset Retries & Set Status.
        WriteCUDA(&d_retries[i], 0);
        WriteCUDA(&d_status[i], ACTIVE);
        
        // Commit the newly adapted step size $h$ to VRAM.
        WriteCUDA(&d_h[i], h_new); 
    }} else {{
        // === REJECTED STEP ===
        // Increment Retry Counter.
        const int new_retries = retries + 1;
        WriteCUDA(&d_retries[i], new_retries);

        // Update Status.
        if (new_retries > d_commondata.rkf45_max_retries) {{
        WriteCUDA(&d_status[i], FAILURE_RKF45_REJECTION_LIMIT);
        }} else {{
        WriteCUDA(&d_status[i], REJECTED);
        }}
        
        // Commit the scaled retry step size $h$ to VRAM.
        // Note: d_f_persistent is NOT updated, preserving the valid state.
        WriteCUDA(&d_h[i], h_new); 
    }}
    """

    # Python: Generate the CUDA kernel and the C host wrapper.
    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
        "stream": "stream_idx"
    }

    prefunc, body = generate_kernel_and_launch_code(
        kernel_name="rkf45_finalize_and_control_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_cuda,
        parallelization="cuda",
        launch_dict=launch_dict,
        cfunc_decorators="__global__",
    )

    # Python: Define arguments for C-Function registration strictly before the call.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "cuda_intrinsics.h", "<math.h>"]
    
    desc = r"""@brief Finalizes the RKF45 step, checks errors, and updates state/stepsize.
    
    @param d_f_persistent Device pointer to the persistent state (updated on acceptance).
    @param d_f_start Device pointer to the base state (read-only).
    @param d_k_bundle Device pointer to all 6 derivative vectors.
    @param d_h Device pointer to the step size.
    @param d_status Device pointer to the ray status (ACTIVE, REJECTED, FAILURE).
    @param d_affine Device pointer to the affine parameter.
    @param d_retries Device pointer to the retry counter.
    @param chunk_size The number of rays in the current bundle."""

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

    # Python: Register the complete C function (Host Launcher + Device Kernel).
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body
    )