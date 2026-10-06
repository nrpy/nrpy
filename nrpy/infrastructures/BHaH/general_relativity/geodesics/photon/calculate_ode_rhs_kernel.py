# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/calculate_ode_rhs_kernel.py
r"""
Provides the kernel generation logic for computing photon geodesic derivatives.

This module constructs and registers the CUDA or OpenMP function used during each
RKF45 stage of photon-geodesic integration. It translates the SymPy expressions for
the coordinate and momentum derivatives into C code. For every photon, the generated
function reads the state, metric, and Christoffel-symbol arrays, evaluates the nine
geodesic right-hand sides, and writes them to the selected RKF45 stage.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

from typing import List

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par


def calculate_ode_rhs_kernel(
    geodesic_rhs_expressions: List[sp.Expr], coordinate_symbols: List[sp.Symbol]
) -> None:
    r"""
    Register the CUDA or OpenMP function for computing the geodesic ODE right-hand sides.

    The generated function assigns array components to the local scalar names used by
    the geodesic equations. It computes the nine derivative components and writes them
    to the selected stage in the RKF45 derivative array.

    :param geodesic_rhs_expressions: The mathematical right-hand side evaluations
        representing the geodesic equations.
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

    # Define arguments for the CUDA and OpenMP functions.
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

    # Assign array components to the local scalar names expected by ccg.c_codegen.
    preamble_lines = [
        "//==========================================",
        "// STATE VECTOR & COORDINATE UNPACKING",
        "//==========================================",
        "// Load spacetime coordinates $x^{\\mu}$ from the state array.",
    ]

    for j, sym in enumerate(coordinate_symbols):
        if str(sym) in used_symbol_names:
            preamble_lines.append(
                f"const double {str(sym)} = d_f_temp_bundle[IDX_F({j}, i)]; // Coordinate ${str(sym)}$ for this photon."
            )

    preamble_lines.extend(
        [
            "\n    //==========================================\n    // MOMENTUM UNPACKING\n    //==========================================",
            "// Load contravariant four-momenta $p^{\\mu}$ from the state array.",
        ]
    )
    for j in range(4):
        if f"pU{j}" in used_symbol_names:
            preamble_lines.append(
                f"const double pU{j} = d_f_temp_bundle[IDX_F({j+4}, i)]; // Momentum component $p^{{{j}}}$ for this photon."
            )

    preamble_lines.extend(
        [
            "\n    //==========================================\n    // METRIC TENSOR UNPACKING\n    //==========================================",
            "// Load the symmetric covariant metric $g_{\\mu\\nu}$ from the pre-calculated memory array.",
        ]
    )
    curr_idx = 0
    for m in range(4):
        for n in range(m, 4):
            comp_name = f"metric_g4DD{m}{n}"
            if comp_name in used_symbol_names:
                preamble_lines.append(
                    f"const double {comp_name} = d_metric_bundle[IDX_METRIC({curr_idx}, i)]; // Metric component $g_{{{m}{n}}}$ for this photon."
                )
            curr_idx += 1

    preamble_lines.extend(
        [
            "\n    //==========================================\n    // CHRISTOFFEL CONNECTION UNPACKING\n    //==========================================",
            "// Load Christoffel symbols $\\Gamma^{\\alpha}_{\\mu\\nu}$ from the pre-calculated memory array.",
        ]
    )
    curr_idx = 0
    for a in range(4):
        for m in range(4):
            for n in range(m, 4):
                comp_name = f"conn_Gamma4UDD{a}{m}{n}"
                if comp_name in used_symbol_names:
                    preamble_lines.append(
                        f"const double {comp_name} = d_connection_bundle[IDX_CONN({curr_idx}, i)]; // Christoffel symbol $\\Gamma^{{{a}}}_{{{m}{n}}}$ for this photon."
                    )
                curr_idx += 1

    preamble_unpacking_str = "\n    ".join(preamble_lines)

    # Generate the raw C math string from the SymPy expressions.
    # Output targets are local scalars k_out_0 through k_out_8.
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
    //==========================================
    // CUDA THREAD IDENTIFICATION
    //==========================================
    // The identifier $i$ represents the global thread index mapped to a specific photon ray.
    const long int i = blockIdx.x * blockDim.x + threadIdx.x; // Thread ID maps to unique photon index.

    // Guard prevents out-of-bounds memory access for threads exceeding the active chunk.
    if (i >= chunk_size) return;
    """
        loop_postamble = ""
        # Translate SIMD macro signatures to native CUDA hardware intrinsics.
        body_math = raw_c_code.replace("SIMD", "CUDA")
    else:
        loop_preamble = """
    //==========================================
    // OPENMP PARALLEL LOOP
    //==========================================
    // Distribute photon rays across available CPU threads for parallel evaluation.
    #pragma omp parallel for
    for(long int i = 0; i < chunk_size; i++) {
    """
        loop_postamble = "    } // END LOOP: for i over chunk_size rays"
        body_math = raw_c_code

    core_math = rf"""
    //==========================================
    // MACRO DEFINITIONS FOR ARRAY ACCESS
    //==========================================
    // IDX_F maps a component to the flattened state array using SoA layout.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id)) // Computes the 1D index for the state array.
    // IDX_METRIC maps a component to the flattened symmetric metric array.
    #define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id)) // Computes the 1D index for the metric array.
    // IDX_CONN maps a component to the flattened Christoffel connection array.
    #define IDX_CONN(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id)) // Computes the 1D index for the connection array.
    // IDX_K maps a stage and component triplet to the flattened derivative array.
    #define IDX_K(s, c, ray_id) (((s) - 1) * 9 * BUNDLE_CAPACITY + (c) * BUNDLE_CAPACITY + (ray_id)) // Computes the 1D index for the RKF45 derivative array.

    {preamble_unpacking_str}

    //==========================================
    // GEODESIC RHS EVALUATION
    //==========================================
    // Local scalars for the evaluated derivatives $\dot{{f}}$.
    double k_out_0, k_out_1, k_out_2, k_out_3, k_out_4, k_out_5, k_out_6, k_out_7, k_out_8; // Nine geodesic right-hand sides for this photon.

    // Evaluate the derivatives $dx^{{\mu}}/d\lambda$ and $dp^{{\mu}}/d\lambda$.
    {body_math}

    //==========================================
    // RKF45 DERIVATIVE ARRAY WRITE
    //==========================================
    // Write the computed derivatives to the correct RKF45 stage offset within the RKF45 derivative array.
    d_k_bundle[IDX_K(stage, 0, i)] = k_out_0; // Write derivative component $0$ to memory.
    d_k_bundle[IDX_K(stage, 1, i)] = k_out_1; // Write derivative component $1$ to memory.
    d_k_bundle[IDX_K(stage, 2, i)] = k_out_2; // Write derivative component $2$ to memory.
    d_k_bundle[IDX_K(stage, 3, i)] = k_out_3; // Write derivative component $3$ to memory.
    d_k_bundle[IDX_K(stage, 4, i)] = k_out_4; // Write derivative component $4$ to memory.
    d_k_bundle[IDX_K(stage, 5, i)] = k_out_5; // Write derivative component $5$ to memory.
    d_k_bundle[IDX_K(stage, 6, i)] = k_out_6; // Write derivative component $6$ to memory.
    d_k_bundle[IDX_K(stage, 7, i)] = k_out_7; // Write derivative component $7$ to memory.
    d_k_bundle[IDX_K(stage, 8, i)] = k_out_8; // Write derivative component $8$ to memory.

    //==========================================
    // MACRO CLEANUP
    //==========================================
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

    prefunc_kernel, launch_code = parallel_utils.generate_kernel_and_launch_code(
        kernel_name="calculate_ode_rhs_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_host,
        parallelization=parallelization,
        launch_dict=launch_dict,
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
        thread_tiling_macro_suffix="RKF45",
    )

    # Step 1: Canonical sequence
    prefunc = prefunc_kernel

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if parallelization == "cuda":
        includes.append("cuda_intrinsics.h")

    desc = r""" Runs the CUDA or OpenMP function for the photon geodesic ODE right-hand sides.

    @param d_f_temp_bundle Pointer to the intermediate state array $f^{\mu}$ in memory.
    @param d_metric_bundle Pointer to the pre-calculated metric array $g_{\mu\nu}$ in memory.
    @param d_connection_bundle Pointer to the pre-calculated connection array $\Gamma^{\alpha}_{\beta\gamma}$ in memory.
    @param d_k_bundle Pointer to the RKF45 derivative array in memory.
    @param stage The current RKF45 stage index used to offset the write location.
    @param chunk_size The number of active rays in the current ray chunk.
    @param stream_idx Work-array index; CUDA uses the corresponding stream.
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

    # Register the complete C function
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
