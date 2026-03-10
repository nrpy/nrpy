"""
Generates the native CUDA kernel and host-side orchestrator for the RKF45 Stage Update (Kernel 5).

Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
This module provides the global memory kernel responsible for evaluating the intermediate
stages of the Runge-Kutta-Fehlberg 4(5) algorithm for relativistic ray tracing on Numerical Spacetimes.
It strictly enforces a Split-Pipeline Architecture, operating on "Streaming Bundles" in VRAM to minimize
bandwidth overhead and avoid fused kernels that cause register spilling on the RTX 3080 target hardware.
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

    # Define the argument dictionary for the CUDA device pointers.
    # Note: d_k_bundle is a flattened 1D array of size [6 * 9 * bundle_capacity].
    arg_dict_cuda = {
        "d_f_start": "const double *restrict",
        "d_k_bundle": "const double *restrict",
        "d_h": "const double *restrict",
        "stage": "const int",
        "chunk_size": "const int",
        "d_f_temp": "double *restrict"
    }

    # Define the argument dictionary for the host pointers and structs.
    arg_dict_host = {
        "d_f_start": "const double *restrict",
        "d_k_bundle": "const double *restrict",
        "d_h": "const double *restrict",
        "stage": "const int",
        "chunk_size": "const int",
        "d_f_temp": "double *restrict"
    }

    # Define the GPU kernel body using raw strings.
    kernel_body = r"""
    // --- THREAD IDENTIFICATION ---
    // Establishes the thread execution index for the 1D computational grid.
    // Architectural Justification: Thread ID maps to a unique photon index to ensure coalesced memory access.
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    // Boundary check for the active bundle size.
    // Architectural Justification: Truncation prevents out-of-bounds VRAM access.
    if (i >= chunk_size) return;

    // --- MACRO DEFINITIONS FOR BUNDLE ACCESS ---
    // Mapping function for the state bundle layout $f^mu$.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))

    // Mapping function for the derivative bundle layout $k^mu$.
    #define IDX_K(s, c, ray_id) ((s) * 9 * BUNDLE_CAPACITY + (c) * BUNDLE_CAPACITY + (ray_id))

    // --- STATE LOADING ---
    // Loads the step size $h$ into local registers prior to tensor evaluation.
    // Architectural Justification: Minimizes repeated global memory accesses during the 9-component loop.
    const double h = ReadCUDA(&d_h[i]);

    // --- BUTCHER TABLEAU EVALUATION ---
    // Evaluates the intermediate Runge-Kutta stages using pre-computed derivative bundles.
    // Architectural Justification: Fused multiply-add intrinsics are utilized to ensure exact IEEE 754 rounding behavior across stages.
    for (int comp = 0; comp < 9; ++comp) {
        
        // Load the base state component $f_{start}$ from VRAM.
        const double f_n = ReadCUDA(&d_f_start[IDX_F(comp, i)]);
        
        // Accumulator for the intermediate update step $f_{temp}$.
        double update_val = 0.0;

        // Apply coefficients based on the current RKF45 stage.
        switch (stage) {
        case 1:
        // Compute intermediate state for $k_2$.
        update_val = MulCUDA(0.25, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]));
        break;
        case 2:
        // Compute intermediate state for $k_3$.
        update_val = FusedMulAddCUDA(0.09375, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]), 
                                    MulCUDA(0.28125, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)])));
        break;
        case 3:
        // Compute intermediate state for $k_4$.
        update_val = FusedMulAddCUDA(1932.0 / 2197.0, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]),
                                    FusedMulAddCUDA(-7200.0 / 2197.0, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]),
                                                    MulCUDA(7296.0 / 2197.0, ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]))));
        break;
        case 4:
        // Compute intermediate state for $k_5$.
        update_val = FusedMulAddCUDA(439.0 / 216.0, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]),
                                    FusedMulAddCUDA(-8.0, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]),
                                                    FusedMulAddCUDA(3680.0 / 513.0, ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]),
                                                                    MulCUDA(-845.0 / 4104.0, ReadCUDA(&d_k_bundle[IDX_K(3, comp, i)])))));
        break;
        case 5:
        // Compute intermediate state for $k_6$.
        update_val = FusedMulAddCUDA(-8.0 / 27.0, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]),
                                    FusedMulAddCUDA(2.0, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]),
                                                    FusedMulAddCUDA(-3544.0 / 2565.0, ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]),
                                                                    FusedMulAddCUDA(1859.0 / 4104.0, ReadCUDA(&d_k_bundle[IDX_K(3, comp, i)]),
                                                                                    MulCUDA(-0.275, ReadCUDA(&d_k_bundle[IDX_K(4, comp, i)]))))));
        break;
        case 6:
        // Finalize derivative computation without updating $f_{temp}$.
        return; 
        }

        // --- GLOBAL VRAM WRITE ---
        // Computes the final intermediate state $f_{temp} = f_n + h * update_val$.
        // Architectural Justification: Output must be written to global VRAM to communicate with the subsequent evaluation kernel due to the split-pipeline constraint.
        const double f_result = FusedMulAddCUDA(h, update_val, f_n);
        
        // Write the intermediate state $f_{temp}$ to the destination bundle in VRAM.
        WriteCUDA(&d_f_temp[IDX_F(comp, i)], f_result);
    }
    """

    # Generate the kernel and the C host wrapper code blocks.
    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
        "stream": "stream_idx"
    }

    prefunc, launch_code = generate_kernel_and_launch_code(
        kernel_name="rkf45_stage_update_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_host,
        parallelization="cuda",
        launch_dict=launch_dict,
        cfunc_decorators="__global__",
    )

    # Define arguments for C-Function registration sequentially.
    prefunc = prefunc
    
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "cuda_intrinsics.h"]
    
    desc = r"""@brief Orchestrates the CUDA kernel for RKF45 intermediate stage updates.
    
    @param d_f_start Pointer to the base state bundle ($f_{start}$) in VRAM.
    @param d_k_bundle Pointer to the flattened derivative array in VRAM.
    @param d_h Pointer to the step size array ($h$) in VRAM.
    @param stage The current RKF45 stage index (1-6).
    @param chunk_size The number of active rays in the current bundle.
    @param d_f_temp Pointer to the destination bundle for the intermediate state ($f_{temp}$).
    @param stream_idx The active CUDA execution stream identifier.
    """

    cfunc_type = "void"
    
    name = "rkf45_stage_update"
    
    params = (
        "const double *restrict d_f_start, "
        "const double *restrict d_k_bundle, "
        "const double *restrict d_h, "
        "const int stage, "
        "const int chunk_size, "
        "double *restrict d_f_temp, "
        "const int stream_idx"
    )
    
    include_CodeParameters_h = False

    body = f"""
    // --- HOST-SIDE ORCHESTRATION ---
    // Wraps the generated launch code to initiate the GPU kernel.
    {launch_code}
    """

    # Register the complete C function (Host Launcher + Device Kernel).
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