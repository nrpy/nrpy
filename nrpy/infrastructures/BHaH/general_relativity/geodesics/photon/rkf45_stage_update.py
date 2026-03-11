r"""
Provides the native kernel and host-side orchestrator for the RKF45 Stage Update (Kernel 5).

Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
This module provides the global memory kernel responsible for evaluating the intermediate
stages of the Runge-Kutta-Fehlberg 4(5) algorithm for relativistic ray tracing on Numerical Spacetimes.
It strictly enforces a Split-Pipeline Architecture, operating on "Streaming Bundles" in VRAM to minimize
bandwidth overhead and avoid fused kernels that cause register spilling on the RTX 3080 target hardware.

Author: Dalton J. Moone.
"""

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code


def rkf45_stage_update() -> None:
    r"""
    Orchestrates the global memory kernel for RKF45 intermediate stage updates.

    The kernel reads the base state $f_{start}$ and the computed derivative vectors $k^{\mu}$
    from memory bundles, applies the Butcher Tableau coefficients, and writes the
    resulting temporary state $f_{temp}$ back to memory for the next interpolation step.
    """
    parallelization = par.parval_from_str("parallelization")

    arg_dict_cuda = {
        "d_f_start": "const double *restrict",
        "d_k_bundle": "const double *restrict",
        "d_h": "const double *restrict",
        "stage": "const int",
        "chunk_size": "const int",
        "d_f_temp": "double *restrict",
    }

    arg_dict_host = {
        "d_f_start": "const double *restrict",
        "d_k_bundle": "const double *restrict",
        "d_h": "const double *restrict",
        "stage": "const int",
        "chunk_size": "const int",
        "d_f_temp": "double *restrict",
    }

    if parallelization == "cuda":
        loop_preamble = """
    // --- CUDA THREAD IDENTIFICATION ---
    // Thread index $i$ maps to a unique photon ray to ensure coalesced VRAM access.
    const long int i = blockIdx.x * blockDim.x + threadIdx.x; // Global thread index $i$.

    // Guard prevents out-of-bounds VRAM access for threads exceeding the active chunk.
    if (i >= chunk_size) return; // Exits if thread index $i$ exceeds $chunk_size$.
    """
        loop_postamble = ""
    else:
        loop_preamble = """
    // --- OPENMP LOOP ARCHITECTURE ---
    // Distribute photon rays across available CPU threads for parallel evaluation.
    #pragma omp parallel for
    for(long int i = 0; i < chunk_size; i++) { // Thread index $i$ maps to a unique photon ray.
    """
        loop_postamble = "    } // End OpenMP loop"

    core_math = r"""
    // --- MACRO DEFINITIONS FOR BUNDLE ACCESS ---
    // Mapping function for the state bundle layout $f^{\mu}$.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))

    // Mapping function for the derivative bundle layout $k^{\mu}$.
    #define IDX_K(s, c, ray_id) ((s) * 9 * BUNDLE_CAPACITY + (c) * BUNDLE_CAPACITY + (ray_id))

    // --- STATE LOADING ---
    // Loading step size $h$ into local registers minimizes repeated global memory accesses during the 9-component loop.
    const double h = ReadCUDA(&d_h[i]); // Local register for step size $h$.

    // --- BUTCHER TABLEAU EVALUATION ---
    // Fused multiply-add intrinsics evaluate the intermediate Runge-Kutta stages to ensure exact IEEE 754 rounding behavior.

    // Bypass the computation entirely for stage 6 to ensure OpenMP compliance.
    if (stage != 6) {
        int comp; // Loop index for iterating over the tensor components.
        for (comp = 0; comp < 9; ++comp) {

            // Load the base state component $f_{start}$ from memory.
            const double f_n = ReadCUDA(&d_f_start[IDX_F(comp, i)]); // Component of the base state $f_{start}$.

            // Accumulator for the intermediate update step $f_{temp}$.
            double update_val = 0.0; // Accumulates the stage update $k^{\mu}$ contributions.

            // Apply coefficients based on the current RKF45 stage.
            switch (stage) {
            case 1:
            // Compute intermediate state for $k_2$.
            update_val = MulCUDA(0.25, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)])); // Applies the $k_1$ coefficient.
            break;
            case 2:
            // Compute intermediate state for $k_3$.
            update_val = FusedMulAddCUDA(0.09375, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]),
                                        MulCUDA(0.28125, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]))); // Applies the $k_1$ and $k_2$ coefficients.
            break;
            case 3:
            // Compute intermediate state for $k_4$.
            update_val = FusedMulAddCUDA(1932.0 / 2197.0, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]),
                                        FusedMulAddCUDA(-7200.0 / 2197.0, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]),
                                                        MulCUDA(7296.0 / 2197.0, ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)])))); // Applies the $k_1$, $k_2$, and $k_3$ coefficients.
            break;
            case 4:
            // Compute intermediate state for $k_5$.
            update_val = FusedMulAddCUDA(439.0 / 216.0, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]),
                                        FusedMulAddCUDA(-8.0, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]),
                                                        FusedMulAddCUDA(3680.0 / 513.0, ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]),
                                                                        MulCUDA(-845.0 / 4104.0, ReadCUDA(&d_k_bundle[IDX_K(3, comp, i)]))))); // Applies coefficients up to $k_4$.
            break;
            case 5:
            // Compute intermediate state for $k_6$.
            update_val = FusedMulAddCUDA(-8.0 / 27.0, ReadCUDA(&d_k_bundle[IDX_K(0, comp, i)]),
                                        FusedMulAddCUDA(2.0, ReadCUDA(&d_k_bundle[IDX_K(1, comp, i)]),
                                                        FusedMulAddCUDA(-3544.0 / 2565.0, ReadCUDA(&d_k_bundle[IDX_K(2, comp, i)]),
                                                                        FusedMulAddCUDA(1859.0 / 4104.0, ReadCUDA(&d_k_bundle[IDX_K(3, comp, i)]),
                                                                                        MulCUDA(-0.275, ReadCUDA(&d_k_bundle[IDX_K(4, comp, i)])))))); // Applies coefficients up to $k_5$.
            break;
            }

            // --- GLOBAL MEMORY WRITE ---
            // Writing the computed update $f_{temp} = f_n + h \times update_val$ to global memory strictly enforces the split-pipeline communication constraint.
            const double f_result = FusedMulAddCUDA(h, update_val, f_n); // Computes the step update and stores it in $f_{result}$.

            // Write the intermediate state $f_{temp}$ to the destination bundle in memory.
            WriteCUDA(&d_f_temp[IDX_F(comp, i)], f_result); // Writes $f_{temp}$ to global memory.
        }
    }

    // --- MACRO CLEANUP ---
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

    prefunc, launch_code = generate_kernel_and_launch_code(
        kernel_name="rkf45_stage_update_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_host,
        parallelization=parallelization,
        launch_dict=launch_dict,
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if parallelization == "cuda":
        includes.append("cuda_intrinsics.h")

    desc = r"""@brief Orchestrates the memory kernel for RKF45 intermediate stage updates.

    @param d_f_start Pointer to the base state bundle ($f_{start}$) in memory.
    @param d_k_bundle Pointer to the flattened derivative array $k^{\mu}$ in memory.
    @param d_h Pointer to the step size array $h$ in memory.
    @param stage The current RKF45 stage index ($1-6$).
    @param chunk_size The number of active rays in the current bundle.
    @param d_f_temp Pointer to the destination bundle for the intermediate state ($f_{temp}$).
    @param stream_idx The active execution stream identifier.
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
    // Wraps the generated launch code to initiate the execution kernel.
    {launch_code}
    """

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
