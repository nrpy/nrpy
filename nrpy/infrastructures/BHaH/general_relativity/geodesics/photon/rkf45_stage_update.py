"""
Generates the native CUDA kernel and host-side orchestrator for the RKF45 Stage Update (Kernel 5).

Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
This module provides the global memory kernel responsible for evaluating the intermediate
stages of the RKF45 algorithm. It operates on "Streaming Bundles" in VRAM to minimize
bandwidth overhead on the RTX 3080.

Author: Dalton J. Moone.
"""

import nrpy.c_function as cfc
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code

def rkf45_stage_update() -> None:
    """
    Register the global CUDA kernel for RKF45 intermediate stage updates.

    Generates the C function `rkf45_stage_update_kernel` and its host launcher.
    The kernel reads the base state $f_{start}$ and the computed derivative vectors $k$
    from VRAM bundles, applies the Butcher Tableau coefficients, and writes the
    resulting temporary state $f_{temp}$ back to VRAM for the next interpolation step.

    :raises TypeError: If incorrect parameters are passed to the code generation functions.
    """
    
    # Python: Define the argument dictionary for the CUDA kernel.
    # Note: d_k_bundle is a flattened 1D array of size [6 * 9 * bundle_capacity].
    arg_dict = {
        "d_f_start": "const double *restrict",
        "d_k_bundle": "const double *restrict",
        "d_h": "const double *restrict",
        "stage": "const int",
        "chunk_size": "const int",
        "d_f_temp": "double *restrict"
    }

    # Python: Define the GPU kernel body using raw strings.
    # We strictly use cuda_intrinsics.h macros for all arithmetic.
    # WE EXPLICITLY DEFINE MACROS HERE TO ENSURE ATOMICITY.
    kernel_body = r"""
    // --- THREAD IDENTIFICATION ---
    // The identifier $i$ represents the local thread index within the current bundle batch.
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    // Architectural Justification: Truncation prevents out-of-bounds VRAM access.
    if (i >= chunk_size) return;

    // --- MACRO DEFINITIONS FOR BUNDLE ACCESS ---
    // IDX_F maps a (component, ray) pair to the flattened state bundle.
    // Layout: [Component][RayID] (Structure of Arrays for coalescing).
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))

    // IDX_K maps a (stage, component, ray) triplet to the flattened derivative bundle.
    // Layout: [Stage][Component][RayID].
    #define IDX_K(s, c, ray_id) ((s) * 9 * BUNDLE_CAPACITY + (c) * BUNDLE_CAPACITY + (ray_id))

    // --- STATE LOADING ---
    // Load the local step size $h$ from VRAM.
    const double h = ReadCUDA(&d_h[i]);

    // --- BUTCHER TABLEAU EVALUATION ---
    // Iterate over all 9 physical components of the state vector $f^\mu$.
    for (int comp = 0; comp < 9; ++comp) {
        // Load the base state component $f_{start}$ from VRAM.
        const double f_n = ReadCUDA(&d_f_start[IDX_F(comp, i)]);
        
        // Initialize the accumulator for the update step.
        double update_val = 0.0;

        // Apply coefficients based on the current RKF45 stage.
        // Math operations use intrinsics to ensure precise rounding modes.
        switch (stage) {
        case 1:
        // We just finished k1. Prepare f_temp for k2: f_n + (1/4) * h * k1
        update_val = MulCUDA(0.25, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]));
        break;
        case 2:
        // We just finished k2. Prepare f_temp for k3: f_n + h * (3/32*k1 + 9/32*k2)
        update_val = FusedMulAddCUDA(0.09375, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]), 
                                    MulCUDA(0.28125, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)])));
        break;
        case 3:
        // Prepare for k4.
        update_val = FusedMulAddCUDA(1932.0 / 2197.0, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]),
                                    FusedMulAddCUDA(-7200.0 / 2197.0, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]),
                                                    MulCUDA(7296.0 / 2197.0, ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]))));
        break;
        case 4:
        // Prepare for k5.
        update_val = FusedMulAddCUDA(439.0 / 216.0, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]),
                                    FusedMulAddCUDA(-8.0, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]),
                                                    FusedMulAddCUDA(3680.0 / 513.0, ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]),
                                                                    MulCUDA(-845.0 / 4104.0, ReadCUDA(&d_k_bundle[IDX_K(3, comp, i)])))));
        break;
        case 5:
        // Prepare for k6.
        update_val = FusedMulAddCUDA(-8.0 / 27.0, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]),
                                    FusedMulAddCUDA(2.0, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]),
                                                    FusedMulAddCUDA(-3544.0 / 2565.0, ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]),
                                                                    FusedMulAddCUDA(1859.0 / 4104.0, ReadCUDA(&d_k_bundle[IDX_K(3, comp, i)]),
                                                                                    MulCUDA(-0.275, ReadCUDA(&d_k_bundle[IDX_K(4, comp, i)]))))));
        break;
        case 6:
        // k6 is the final derivative; we don't need to update f_temp anymore 
        // because rkf45_finalize_and_control will compute the final f_n+1.
        return; 
        }

        // --- GLOBAL VRAM WRITE ---
        // Compute $f_{temp} = f_n + h * update\_val$.
        const double f_result = FusedMulAddCUDA(h, update_val, f_n);
        
        // Write the intermediate state to the destination bundle.
        WriteCUDA(&d_f_temp[IDX_F(comp, i)], f_result);
    }
    """

    # Python: Generate the kernel and the C host wrapper.
    # We use a standard 1D block size. The orchestrator determines the grid size.
    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
    }

    prefunc, body = generate_kernel_and_launch_code(
        kernel_name="rkf45_stage_update_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict,
        arg_dict_host=arg_dict,
        parallelization="cuda",
        launch_dict=launch_dict,
        cfunc_decorators="__global__",
    )

    # Python: Define arguments for C-Function registration strictly before the call.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "cuda_intrinsics.h"]
    
    desc = r"""@brief Orchestrates the CUDA kernel for RKF45 intermediate stage updates.
    
    @param d_f_start Pointer to the base state bundle ($f_{start}$) in VRAM.
    @param d_k_bundle Pointer to the flattened derivative array in VRAM.
    @param d_h Pointer to the step size array in VRAM.
    @param stage The current RKF45 stage index (1-6).
    @param chunk_size The number of active rays in the current bundle.
    @param d_f_temp Pointer to the destination bundle for the intermediate state ($f_{temp}$).
    """

    cfunc_type = "void"
    name = "rkf45_stage_update"
    
    params = (
        "const double *restrict d_f_start, "
        "const double *restrict d_k_bundle, "
        "const double *restrict d_h, "
        "const int stage, "
        "const int chunk_size, "
        "double *restrict d_f_temp"
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
        body=body,
    )