r"""
Evaluates the temporal momentum component and provides the host-side orchestrator.

This module evaluates the quadratic Hamiltonian constraint to enforce physical null
trajectories for a batch of photons. It operates to find the negative root of the
constraint equation for $p^t$.

Author: Dalton J. Moone.
"""

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code


def p0_reverse_kernel(p0_expr: sp.Expr) -> None:
    r"""
    Orchestrate the global kernel for the initial temporal momentum calculation.

    The kernel loads the photon's spatial momenta and metric components
    from global memory, executes the Hamiltonian constraint solver,
    and writes the resulting temporal momentum $p^0$ back to the state vector bundle.

    :param p0_expr: The SymPy expression representing the negative root of the Hamiltonian constraint.
    :raises ValueError: If the symbolic expression fails to parse during code generation.
    """
    if not p0_expr:
        raise ValueError("p0_expr must contain a valid symbolic expression.")

    parallelization = par.parval_from_str("parallelization")

    arg_dict = {
        "d_f_bundle": "double *restrict",
        "d_metric_bundle": "const double *restrict",
        "chunk_size": "const int",
    }

    # Dynamically format read/write macros and SIMD generation based on target architecture.
    if parallelization == "cuda":
        read_fmt = "ReadCUDA(&{0})"
        write_fmt = "WriteCUDA(&{0}, {1})"
        enable_simd = True
        loop_preamble = """
    // --- CUDA THREAD IDENTIFICATION ---
    // The identifier $i$ maps directly to a unique photon ray index in the global VRAM batch.
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    // Guard prevents out-of-bounds VRAM access for threads exceeding the active chunk.
    if (i >= chunk_size) return;
    """
        loop_postamble = ""
    else:
        read_fmt = "{0}"
        write_fmt = "{0} = {1}"
        enable_simd = False
        loop_preamble = """
    // --- OPENMP LOOP ARCHITECTURE ---
    // Distribute photon rays across available CPU threads for parallel evaluation.
    #pragma omp parallel for
    for(int i = 0; i < chunk_size; i++) {
    """
        loop_postamble = "    } // End OpenMP loop"

    # Generate the raw C math string from the SymPy expression.
    body_math = ccg.c_codegen(
        [p0_expr],
        ["p0_val"],
        enable_cse=True,
        enable_simd=enable_simd,
        verbose=False,
        include_braces=False,
    )

    # Translate SIMD macro signatures to native CUDA hardware intrinsics if targeting GPU.
    if parallelization == "cuda":
        body_math = body_math.replace("SIMD", "CUDA")

    # Map the 2D symmetric metric components to explicit local register loads.
    metric_loads = []
    k = 0
    for m in range(4):
        for n in range(m, 4):
            comp_name = f"metric_g4DD{m}{n}"
            load_str = read_fmt.format(f"d_metric_bundle[IDX_METRIC({k}, i)]")
            metric_loads.append(
                f"// Covariant metric component $g_{{{m}{n}}}$.\n    const double {comp_name} = {load_str};"
            )
            k += 1
    metric_load_str = "\n    ".join(metric_loads)

    # Dynamically inject access patterns based on target architecture.
    pU1_load = read_fmt.format("d_f_bundle[IDX_F(5, i)]")
    pU2_load = read_fmt.format("d_f_bundle[IDX_F(6, i)]")
    pU3_load = read_fmt.format("d_f_bundle[IDX_F(7, i)]")
    p0_write = write_fmt.format("d_f_bundle[IDX_F(4, i)]", "p0_val")

    core_math = f"""
    // --- MACRO DEFINITIONS FOR BUNDLE ACCESS ---
    // IDX_F maps a component to the flattened state bundle using SoA layout.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
    // IDX_METRIC maps a component to the flattened symmetric metric bundle.
    #define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))

    // --- SPATIAL MOMENTUM UNPACKING ---
    // Unpacking variables from device memory to local registers minimizes global memory transactions.
    // Contravariant spatial momentum component $p^x$.
    const double pU1 = {pU1_load};
    // Contravariant spatial momentum component $p^y$.
    const double pU2 = {pU2_load};
    // Contravariant spatial momentum component $p^z$.
    const double pU3 = {pU3_load};

    // --- METRIC TENSOR UNPACKING ---
    // Load the 10 independent metric components $g_{{\\mu\\nu}}$ from memory into explicitly named registers.
    {metric_load_str}

    // --- HAMILTONIAN CONSTRAINT ROOT FINDING ---
    // Evaluate the algebraic solution for the temporal momentum $p^0$.
    // Local register storing the final evaluated negative root.
    double p0_val = 0.0;
    {body_math}

    // --- GLOBAL MEMORY WRITE ---
    // Write the resulting temporal momentum back to component index 4 of the state bundle in memory.
    {p0_write};

    // --- MACRO CLEANUP ---
    // Undefine macros to ensure hermetic compilation and prevent redefinition errors.
    #undef IDX_F
    #undef IDX_METRIC
    """

    kernel_body = f"{loop_preamble}\n{core_math}\n{loop_postamble}"

    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
        "stream": "stream_idx",
    }

    prefunc, launch_code = generate_kernel_and_launch_code(
        kernel_name="p0_reverse_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict,
        arg_dict_host=arg_dict,
        parallelization=parallelization,
        launch_dict=launch_dict,
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
        thread_tiling_macro_suffix="RKF45",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if parallelization == "cuda":
        includes.append("cuda_intrinsics.h")

    desc = r"""@brief Orchestrates the kernel for the initial temporal momentum calculation.

    @param d_f_bundle Pointer to the state vector bundle $f^{\mu}$ in memory.
    @param d_metric_bundle Pointer to the pre-calculated metric bundle $g_{\mu\nu}$ in memory.
    @param chunk_size The number of active rays in the current bundle batch.
    @param stream_idx The hardware stream index for asynchronous execution.
    """

    cfunc_type = "void"

    name = "p0_reverse_kernel"

    params = (
        "double *restrict d_f_bundle, "
        "const double *restrict d_metric_bundle, "
        "const int chunk_size, "
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
