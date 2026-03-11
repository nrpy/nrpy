r"""
Module defines the global kernel and host-side orchestrator for computing the photon geodesic ODE right-hand sides.

Module evaluates the spatial and temporal derivatives $\dot{f}$ required during
the RKF45 integration step. It reads pre-calculated metric and connection tensors
from global memory bundles to minimize register pressure on the target architecture.

Author: Dalton J. Moone.
"""

from typing import List

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code


def calculate_ode_rhs_kernel(
    geodesic_rhs_expressions: List[sp.Expr], coordinate_symbols: List[sp.Symbol]
) -> None:
    r"""
    Register the global kernel to compute the ODE right-hand side.

    The generated kernel maps memory tensor data into thread-local registers matching
    the symbols expected by the generated geodesic equations. It computes the $9$
    derivative components and writes them to the specific stage offset in the
    RKF45 derivative bundle.

    :param geodesic_rhs_expressions: The mathematical right-hand side evaluations representing the geodesic equations.
    :param coordinate_symbols: The spatial and temporal coordinate variables in order.
    :raises ValueError: If the provided geodesic expression list is empty.
    """
    if not geodesic_rhs_expressions:
        raise ValueError(
            "geodesic_rhs_expressions must contain at least one mathematical expression."
        )

    parallelization = par.parval_from_str("parallelization")

    # Identify all unique mathematical symbols used in the generated RHS expressions.
    used_symbol_names = {
        str(sym) for expr in geodesic_rhs_expressions for sym in expr.free_symbols
    }

    # Define the argument dictionary for the hardware kernel generation.
    arg_dict_cuda = {
        "d_f_temp_bundle": "const double *restrict",
        "d_metric_bundle": "const double *restrict",
        "d_connection_bundle": "const double *restrict",
        "d_k_bundle": "double *restrict",
        "stage": "const int",
        "chunk_size": "const long int",
    }

    arg_dict_host = {
        "d_f_temp_bundle": "const double *restrict",
        "d_metric_bundle": "const double *restrict",
        "d_connection_bundle": "const double *restrict",
        "d_k_bundle": "double *restrict",
        "stage": "const int",
        "chunk_size": "const long int",
    }

    # Build the memory data unpacking block.
    # This maps global memory directly to the local scalar registers expected by ccg.c_codegen.
    preamble_lines = [
        "// --- STATE VECTOR & COORDINATE UNPACKING ---",
        "// Load spacetime coordinates $x^{\\mu}$ from the global state bundle.",
    ]

    for j, sym in enumerate(coordinate_symbols):
        if str(sym) in used_symbol_names:
            preamble_lines.append(
                f"const double {str(sym)} = d_f_temp_bundle[IDX_F({j}, i)]; // Maps the coordinate ${str(sym)}$ from global memory to a thread-local register."
            )

    preamble_lines.extend(
        [
            "\n    // --- MOMENTUM UNPACKING ---",
            "// Load contravariant four-momenta $p^{\\mu}$ from the global state bundle.",
        ]
    )
    for j in range(4):
        if f"pU{j}" in used_symbol_names:
            preamble_lines.append(
                f"const double pU{j} = d_f_temp_bundle[IDX_F({j+4}, i)]; // Maps the momentum component $p^{{{j}}}$ from global memory to a thread-local register."
            )

    preamble_lines.extend(
        [
            "\n    // --- METRIC TENSOR UNPACKING ---",
            "// Load the symmetric covariant metric $g_{\\mu\\nu}$ from the pre-calculated memory bundle.",
        ]
    )
    curr_idx = 0
    for m in range(4):
        for n in range(m, 4):
            comp_name = f"metric_g4DD{m}{n}"
            if comp_name in used_symbol_names:
                preamble_lines.append(
                    f"const double {comp_name} = d_metric_bundle[IDX_METRIC({curr_idx}, i)]; // Maps the metric component $g_{{{m}{n}}}$ from global memory to a thread-local register."
                )
            curr_idx += 1

    preamble_lines.extend(
        [
            "\n    // --- CHRISTOFFEL CONNECTION UNPACKING ---",
            "// Load Christoffel symbols $\\Gamma^{\\alpha}_{\\mu\\nu}$ from the pre-calculated memory bundle.",
        ]
    )
    curr_idx = 0
    for a in range(4):
        for m in range(4):
            for n in range(m, 4):
                comp_name = f"conn_Gamma4UDD{a}{m}{n}"
                if comp_name in used_symbol_names:
                    preamble_lines.append(
                        f"const double {comp_name} = d_connection_bundle[IDX_CONN({curr_idx}, i)]; // Maps the connection component $\\Gamma^{{{a}}}_{{{m}{n}}}$ from global memory to a thread-local register."
                    )
                curr_idx += 1

    preamble_unpacking_str = "\n    ".join(preamble_lines)

    # Generate the raw C math string from the SymPy expressions.
    # Output targets are local scalar registers k_out_0 through k_out_8.
    k_array_outputs = [f"k_out_{j}" for j in range(9)]

    enable_simd = parallelization == "cuda"
    raw_c_code = ccg.c_codegen(
        geodesic_rhs_expressions,
        k_array_outputs,
        enable_cse=True,
        enable_simd=enable_simd,
        include_braces=False,
        verbose=False,
    )

    if parallelization == "cuda":
        loop_preamble = """
    // --- CUDA THREAD IDENTIFICATION ---
    // The identifier $i$ represents the global thread index mapped to a specific photon ray.
    const long int i = blockIdx.x * blockDim.x + threadIdx.x; // Thread ID maps to unique photon index.

    // Hardware Justification: Guard prevents out-of-bounds memory access for threads exceeding the active chunk.
    if (i >= chunk_size) return;
    """
        loop_postamble = ""
        # Translate SIMD macro signatures to native CUDA hardware intrinsics.
        body_math = raw_c_code.replace("SIMD", "CUDA")
    else:
        loop_preamble = """
    // --- OPENMP LOOP ARCHITECTURE ---
    // Distribute photon rays across available CPU threads for parallel evaluation.
    #pragma omp parallel for
    for(long int i = 0; i < chunk_size; i++) {
    """
        loop_postamble = "    } // End OpenMP loop"
        body_math = raw_c_code

    core_math = rf"""
    // --- MACRO DEFINITIONS FOR BUNDLE ACCESS ---
    // IDX_F maps a component to the flattened state bundle using SoA layout.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id)) // Computes the 1D index for the state bundle.
    // IDX_METRIC maps a component to the flattened symmetric metric bundle.
    #define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id)) // Computes the 1D index for the metric bundle.
    // IDX_CONN maps a component to the flattened Christoffel connection bundle.
    #define IDX_CONN(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id)) // Computes the 1D index for the connection bundle.
    // IDX_K maps a stage and component triplet to the flattened derivative bundle.
    // Hardware Justification: The stage index dictates the massive offset $stage \\times 9 \\times Capacity$.
    #define IDX_K(s, c, ray_id) (((s) - 1) * 9 * BUNDLE_CAPACITY + (c) * BUNDLE_CAPACITY + (ray_id)) // Computes the 1D index for the RKF45 derivative bundle.

    {preamble_unpacking_str}

    // --- GEODESIC RHS EVALUATION ---
    // Local register declarations to capture the evaluated derivatives $\dot{{f}}$.
    double k_out_0, k_out_1, k_out_2, k_out_3, k_out_4, k_out_5, k_out_6, k_out_7, k_out_8; // Thread-local registers allocate memory for the $9$ derivative outputs.

    // Evaluate the derivatives $dx^{{\mu}}/d\lambda$ and $dp^{{\mu}}/d\lambda$ using hardware FMA instructions.
    {body_math}

    // --- GLOBAL MEMORY WRITE ---
    // Write the computed derivatives to the correct RKF45 stage offset within the massive derivative bundle.
    d_k_bundle[IDX_K(stage, 0, i)] = k_out_0; // Write derivative component $0$ to memory.
    d_k_bundle[IDX_K(stage, 1, i)] = k_out_1; // Write derivative component $1$ to memory.
    d_k_bundle[IDX_K(stage, 2, i)] = k_out_2; // Write derivative component $2$ to memory.
    d_k_bundle[IDX_K(stage, 3, i)] = k_out_3; // Write derivative component $3$ to memory.
    d_k_bundle[IDX_K(stage, 4, i)] = k_out_4; // Write derivative component $4$ to memory.
    d_k_bundle[IDX_K(stage, 5, i)] = k_out_5; // Write derivative component $5$ to memory.
    d_k_bundle[IDX_K(stage, 6, i)] = k_out_6; // Write derivative component $6$ to memory.
    d_k_bundle[IDX_K(stage, 7, i)] = k_out_7; // Write derivative component $7$ to memory.
    d_k_bundle[IDX_K(stage, 8, i)] = k_out_8; // Write derivative component $8$ to memory.

    // --- MACRO CLEANUP ---
    #undef IDX_F
    #undef IDX_METRIC
    #undef IDX_CONN
    #undef IDX_K
    """

    kernel_body = f"{loop_preamble}\n{core_math}\n{loop_postamble}"

    # Generate the kernel and the C host wrapper.
    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
        "stream": "stream_idx",
    }

    prefunc_kernel, launch_code = generate_kernel_and_launch_code(
        kernel_name="calculate_ode_rhs_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_host,
        parallelization=parallelization,
        launch_dict=launch_dict,
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
        thread_tiling_macro_suffix="RKF45",
    )

    # --- CANONICAL MASTER ORDER SEQUENCE ---
    prefunc = prefunc_kernel

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if parallelization == "cuda":
        includes.append("cuda_intrinsics.h")

    desc = r"""@brief Orchestrates the memory kernel for computing the photon geodesic ODE right-hand sides.

    @param d_f_temp_bundle Pointer to the intermediate state bundle $f^{\mu}$ in memory.
    @param d_metric_bundle Pointer to the pre-calculated metric bundle $g_{\mu\nu}$ in memory.
    @param d_connection_bundle Pointer to the pre-calculated connection bundle $\Gamma^{\alpha}_{\beta\gamma}$ in memory.
    @param d_k_bundle Pointer to the massive derivative bundle array in memory.
    @param stage The current RKF45 stage index used to offset the write location.
    @param chunk_size The number of active rays in the current bundle batch.
    """

    cfunc_type = "void"

    name = "calculate_ode_rhs_kernel"

    params = (
        "const double *restrict d_f_temp_bundle, "
        "const double *restrict d_metric_bundle, "
        "const double *restrict d_connection_bundle, "
        "double *restrict d_k_bundle, "
        "const int stage, "
        "const long int chunk_size,"
        "const int stream_idx"
    )

    include_CodeParameters_h = False

    body = launch_code

    # Register the complete C function using the canonical Master Order sequence.
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
