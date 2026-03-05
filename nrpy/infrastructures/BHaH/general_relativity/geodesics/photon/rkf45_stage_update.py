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
                // Stage 1 is simply the base state itself (no $k$ contribution).
                update_val = 0.0; 
                break;
            case 2:
                // $k_1$ contribution: 1/4
                update_val = MulCUDA(0.25, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]));
                break;
            case 3:
                // $k_1$: 3/32, $k_2$: 9/32
                update_val = FusedMulAddCUDA(0.09375, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]),
                             MulCUDA(0.28125, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)])));
                break;
            case 4:
                // $k_1$: 1932/2197, $k_2$: -7200/2197, $k_3$: 7296/2197
                update_val = FusedMulAddCUDA(1932.0/2197.0, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]),
                             FusedMulAddCUDA(-7200.0/2197.0, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]),
                             MulCUDA(7296.0/2197.0, ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]))));
                break;
            case 5:
                // $k_1$: 439/216, $k_2$: -8, $k_3$: 3680/513, $k_4$: -845/4104
                update_val = FusedMulAddCUDA(439.0/216.0, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]),
                             FusedMulAddCUDA(-8.0, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]),
                             FusedMulAddCUDA(3680.0/513.0, ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]),
                             MulCUDA(-845.0/4104.0, ReadCUDA(&d_k_bundle[IDX_K(3, comp, i)])))));
                break;
            case 6:
                // $k_1$: -8/27, $k_2$: 2, $k_3$: -3544/2565, $k_4$: 1859/4104, $k_5$: -11/40
                update_val = FusedMulAddCUDA(-8.0/27.0, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]),
                             FusedMulAddCUDA(2.0, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]),
                             FusedMulAddCUDA(-3544.0/2565.0, ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]),
                             FusedMulAddCUDA(1859.0/4104.0, ReadCUDA(&d_k_bundle[IDX_K(3, comp, i)]),
                             MulCUDA(-0.275, ReadCUDA(&d_k_bundle[IDX_K(4, comp, i)]))))));
                break;
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