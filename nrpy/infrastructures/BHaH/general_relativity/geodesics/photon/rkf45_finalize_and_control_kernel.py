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
        "chunk_size": "const long int"               # Input: Active bundle size
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
    #define IDX_F(c, ray_id) ((c) * chunk_size + (ray_id))
    // IDX_K: Maps (stage, component, ray) to the flattened derivative bundle.
    #define IDX_K(s, c, ray_id) ((s) * 9 * chunk_size + (c) * chunk_size + (ray_id))

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
    
    // NOTE: For L1 momentum floor, we need the sum of absolute momenta from the START state.
    // We load them specifically here to avoid re-reading inside the loop.
    double p_L1 = 0.0;
    {{
        const double px = ReadCUDA(&d_f_start[IDX_F(5, i)]);
        const double py = ReadCUDA(&d_f_start[IDX_F(6, i)]);
        const double pz = ReadCUDA(&d_f_start[IDX_F(7, i)]);
        p_L1 = AddCUDA(AbsCUDA(px), AddCUDA(AbsCUDA(py), AbsCUDA(pz)));
    }}

    for (int comp = 0; comp < 9; ++comp) {{
        // 1. Load Base State
        const double f_n = ReadCUDA(&d_f_start[IDX_F(comp, i)]);

        // 2. Load Derivative Components (Scalar Loads -> Registers)
        // We only hold these 6 doubles in registers for the duration of this loop iteration.
        const double k0 = ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]);
        const double k1 = ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]);
        const double k2 = ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]);
        const double k3 = ReadCUDA(&d_k_bundle[IDX_K(3, comp, i)]);
        const double k4 = ReadCUDA(&d_k_bundle[IDX_K(4, comp, i)]);
        const double k5 = ReadCUDA(&d_k_bundle[IDX_K(5, comp, i)]);

        // 3. Compute 4th Order Baseline (Error Estimator)
        // Coeffs: 25/216, 0, 1408/2565, 2197/4104, -1/5, 0
        double f_4th = MulCUDA(0.1157407407407407, k0);
        f_4th = FusedMulAddCUDA(0.5489278752436647, k2, f_4th);
        f_4th = FusedMulAddCUDA(0.5353313840155946, k3, f_4th);
        f_4th = FusedMulSubCUDA(0.2, k4, f_4th);
        f_4th = FusedMulAddCUDA(h_local, f_4th, f_n);

        // 4. Compute 5th Order Candidate (Physical Update)
        // Coeffs: 16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55
        double update = MulCUDA(0.1185185185185185, k0);
        update = FusedMulAddCUDA(0.5189863547758285, k2, update);
        update = FusedMulAddCUDA(0.5061137692716641, k3, update);
        update = FusedMulSubCUDA(0.18, k4, update);
        update = FusedMulAddCUDA(0.0363636363636364, k5, update);
        
        const double f_5th_val = FusedMulAddCUDA(h_local, update, f_n);
        
        // Cache the result for potential write-back to Global VRAM.
        f_5th_cache[comp] = f_5th_val;

        // 5. Error Normalization
        if (comp < 8) {{ // Exclude affine param (index 8) from error check
            const double err_abs = AbsCUDA(SubCUDA(f_5th_val, f_4th));
            double scale = 0.0;

            if (comp == 0) {{ // Time t
                scale = atol;
            }} else if (comp <= 3) {{ // Position x, y, z
                scale = AddCUDA(atol, MulCUDA(rtol, AbsCUDA(f_n)));
            }} else if (comp == 4) {{ // Energy p_t
                scale = AddCUDA(atol, MulCUDA(rtol, AbsCUDA(f_n)));
            }} else {{ // Spatial Momentum p_i (Uses L1 Floor)
                scale = AddCUDA(atol, MulCUDA(rtol, p_L1));
            }}
            
            // Accumulate Maximum Normalized Error
            double current_err = DivCUDA(err_abs, scale);
            
            // EXPLICIT NaN REJECTION: Guard against IEEE 754 fmax behavior
            if (isnan(current_err)) {{
                err_norm = 1e30; // Force an artificially massive error to guarantee rejection
            }} else {{
                err_norm = fmax(err_norm, current_err);
            }}
        }}
    }}

    // --- CONTROL LOGIC ---
    if (err_norm <= 1.0) {{
        // === ACCEPTED STEP ===
        
        // 1. Flush Cached State to Persistent VRAM
        // This coalesced write commits the local register cache to global memory.
        #pragma unroll
        for (int comp = 0; comp < 9; ++comp) {{
            WriteCUDA(&d_f_persistent[IDX_F(comp, i)], f_5th_cache[comp]);
        }}

        // 2. Update Affine Parameter lambda
        const double old_affine = ReadCUDA(&d_affine[i]);
        WriteCUDA(&d_affine[i], AddCUDA(old_affine, h_local));

        // 3. Reset Retries & Set Status
        WriteCUDA(&d_retries[i], 0);
        WriteCUDA(&d_status[i], ACTIVE);

        // 4. Grow Step Size (Safety Factor)
        // Scale h based on error norm, bounded by h_max.
        double safety = {cd_access}rkf45_safety_factor;
        double factor = (err_norm > 1e-15) ? pow(DivCUDA(1.0, err_norm), 0.2) : 2.0;
        double h_new = MulCUDA(safety, MulCUDA(h_local, factor));
        h_new = fmin(h_new, {cd_access}rkf45_h_max);
        
        WriteCUDA(&d_h[i], h_new);

    }} else {{
        // === REJECTED STEP ===

        // 1. Increment Retry Counter
        const int new_retries = retries + 1;
        WriteCUDA(&d_retries[i], new_retries);

        // 2. Shrink Step Size
        // Reduce h by half, bounded by h_min.
        double h_new = MulCUDA(0.5, h_local);
        h_new = fmax(h_new, {cd_access}rkf45_h_min);
        WriteCUDA(&d_h[i], h_new);

        // 3. Update Status
        if (new_retries > {cd_access}rkf45_max_retries) {{
            WriteCUDA(&d_status[i], FAILURE_RKF45_REJECTION_LIMIT);
        }} else {{
            WriteCUDA(&d_status[i], REJECTED); 
        }}
        
        // Note: d_f_persistent is NOT updated, preserving the valid state for retry.
    }}
    """

    # Python: Generate the CUDA kernel and the C host wrapper.
    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
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
        "const long int chunk_size"
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