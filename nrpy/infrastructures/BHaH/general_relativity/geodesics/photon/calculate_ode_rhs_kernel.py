# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/calculate_ode_rhs_kernel.py
r"""
Provides the kernel generation logic for computing photon geodesic derivatives.

This module defines the Python function that constructs and registers the host-side
orchestrator and global kernel required during the RKF45 integration step for photon
geodesics. It translates SymPy expressions for spatial and temporal derivatives into
C code and incorporates parallelism macros. The generated kernel maps global memory
bundle data into thread-local registers to evaluate the right-hand sides of the
geodesic equations. A boundary guard prevents out-of-bounds memory access for
threads exceeding the active chunk. The stage index dictates the offset for memory
alignment.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

from typing import List, Optional

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par

SUPPORTED_INTERPOLATION_METHODS = ("g4DD", "g4DD_d0", "GammaUDD")


def calculate_ode_rhs_kernel(
    geodesic_rhs_expressions: List[sp.Expr],
    coordinate_symbols: List[sp.Symbol],
    use_metric_derivative_rhs: Optional[bool] = None,
    normalized_eom: bool = False,
    interpolation_method: Optional[str] = None,
) -> None:
    r"""
    Provide the global kernel registration for computing the ODE right-hand side.

    The generated kernel maps memory tensor data into thread-local registers matching
    the symbols expected by the generated geodesic equations. Direct geodesic
    evolution receives Christoffel symbols, while numerical-spacetime evolution
    receives first metric derivatives. Normalized evolution uses coordinate time from
    the integration-parameter bundle at each Cash-Karp stage. The kernel computes the
    nine derivative components and writes them to the stage offset in the RKF45
    derivative bundle.

    :param geodesic_rhs_expressions: The mathematical right-hand side evaluations
        representing the geodesic equations.
    :param coordinate_symbols: The spatial and temporal coordinate variables in order.
    :param use_metric_derivative_rhs: Legacy selection for whether the RHS
        consumes first metric derivatives instead of Christoffel symbols. If
        `interpolation_method` is supplied, it determines the selection.
    :param normalized_eom: Whether the state stores normalized photon variables
        ``u`` and ``PiD_i`` instead of direct four-momentum ``pU_i``.
    :param interpolation_method: Numerical data method selected at Python code
        generation time. ``GammaUDD`` selects the Christoffel RHS; the other
        methods select the metric-derivative RHS.
    :raises ValueError: If the provided geodesic expression list is empty or if
        `interpolation_method` is unsupported.
    """
    if not geodesic_rhs_expressions:
        raise ValueError(
            "geodesic_rhs_expressions must contain at least one mathematical expression."
        )
    if interpolation_method is not None:
        if interpolation_method not in SUPPORTED_INTERPOLATION_METHODS:
            raise ValueError(
                "interpolation_method must be one of "
                f"{SUPPORTED_INTERPOLATION_METHODS}; found '{interpolation_method}'."
            )
        selected_use_metric_derivative_rhs = interpolation_method != "GammaUDD"
    else:
        selected_use_metric_derivative_rhs = bool(use_metric_derivative_rhs)

    parallelization = par.parval_from_str("parallelization")

    # Select the generated RHS input contract once for this registration.
    geometry_bundle_arg_name = "d_rhs_geometry_bundle"
    geometry_index_macro = "IDX_RHS_GEOMETRY"
    if selected_use_metric_derivative_rhs:
        geometry_bundle_description = "first metric-derivative bundle"
        geometry_symbol_prefix = "metric_g4DD_dD"
        geometry_components = [
            (mu, nu, sigma)
            for mu in range(4)
            for nu in range(mu, 4)
            for sigma in range(4)
        ]
        geometry_index_comment = (
            "IDX_RHS_GEOMETRY maps a symmetric metric component and its "
            "derivative direction to the flattened metric-derivative bundle."
        )
    else:
        geometry_bundle_description = "Christoffel-symbol bundle"
        geometry_symbol_prefix = "conn_Gamma4UDD"
        geometry_components = [
            (alpha, mu, nu)
            for alpha in range(4)
            for mu in range(4)
            for nu in range(mu, 4)
        ]
        geometry_index_comment = (
            "IDX_RHS_GEOMETRY maps a Christoffel component to the flattened "
            "Christoffel-symbol bundle."
        )

    if normalized_eom:
        integration_parameter_args = {
            "d_integration_param_bundle": "const double *restrict",
            "d_h": "const double *restrict",
        }
        coordinate_time_preamble = r"""
    // Coordinate time is the independent integration parameter in normalized evolution.
    double rkf45_stage_time_fraction = 0.0;
    switch (stage) {
      case 1:
        rkf45_stage_time_fraction = 0.0;
        break;
      case 2:
        rkf45_stage_time_fraction = 1.0 / 4.0;
        break;
      case 3:
        rkf45_stage_time_fraction = 3.0 / 8.0;
        break;
      case 4:
        rkf45_stage_time_fraction = 12.0 / 13.0;
        break;
      case 5:
        rkf45_stage_time_fraction = 1.0;
        break;
      case 6:
        rkf45_stage_time_fraction = 1.0 / 2.0;
        break;
      default:
        break;
    } // END SWITCH: select Cash-Karp stage time fraction
    const double t = d_integration_param_bundle[i] +
                     rkf45_stage_time_fraction * d_h[i];
    """
    else:
        integration_parameter_args = {}
        coordinate_time_preamble = ""

    # Identify all unique mathematical symbols used in the generated RHS expressions.
    used_symbol_names = {
        str(sym) for expr in geodesic_rhs_expressions for sym in expr.free_symbols
    }

    # Define the argument dictionary for the hardware kernel generation.
    arg_dict_cuda = {
        "d_f_temp_bundle": "const double *restrict",
        "d_metric_bundle": "const double *restrict",
        geometry_bundle_arg_name: "const double *restrict",
        **integration_parameter_args,
        "d_k_bundle": "double *restrict",
        "stage": "const int",
        "chunk_size": "const long int",
    }

    arg_dict_host = {
        "d_f_temp_bundle": "const double *restrict",
        "d_metric_bundle": "const double *restrict",
        geometry_bundle_arg_name: "const double *restrict",
        **integration_parameter_args,
        "d_k_bundle": "double *restrict",
        "stage": "const int",
        "chunk_size": "const long int",
    }

    # Build the memory data unpacking block.
    # This maps global memory directly to the local scalar registers expected by ccg.c_codegen.
    preamble_lines = [
        "//==========================================",
        "// STATE VECTOR & COORDINATE UNPACKING",
        "//==========================================",
        "// Load spatial coordinates from the global state bundle.",
    ]

    for j, sym in enumerate(coordinate_symbols):
        if j == 0:
            if not normalized_eom and str(sym) in used_symbol_names:
                preamble_lines.append(
                    f"const double {str(sym)} = d_f_temp_bundle[IDX_F({j}, i)];"
                )
        elif str(sym) in used_symbol_names:
            preamble_lines.append(
                f"const double {str(sym)} = d_f_temp_bundle[IDX_F({j}, i)];"
            )

    if coordinate_time_preamble and str(coordinate_symbols[0]) in used_symbol_names:
        preamble_lines.append(coordinate_time_preamble)

    if normalized_eom:
        preamble_lines.extend(
            [
                "\n    //==========================================\n    // NORMALIZED MOMENTUM UNPACKING\n    //==========================================",
                "// Load the normalized photon energy and covariant spatial momentum.",
                "const double u = d_f_temp_bundle[IDX_F(4, i)];",
            ]
        )
        for j in range(3):
            preamble_lines.append(
                f"const double PiD{j} = d_f_temp_bundle[IDX_F({j + 5}, i)];"
            )
    else:
        preamble_lines.extend(
            [
                "\n    //==========================================\n    // MOMENTUM UNPACKING\n    //==========================================",
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
            "\n    //==========================================\n    // METRIC TENSOR UNPACKING\n    //==========================================",
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
            "\n    //==========================================\n    // GEOMETRY UNPACKING\n    //==========================================",
            f"// Load symbols from the {geometry_bundle_description}.",
        ]
    )
    for geometry_index, component_indices in enumerate(geometry_components):
        comp_name = f"{geometry_symbol_prefix}{''.join(map(str, component_indices))}"
        if comp_name in used_symbol_names:
            preamble_lines.append(
                f"const double {comp_name} = {geometry_bundle_arg_name}[{geometry_index_macro}({geometry_index}, i)];"
            )

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
    // OPENMP LOOP ARCHITECTURE
    //==========================================
    // Distribute photon rays across available CPU threads for parallel evaluation.
    #pragma omp parallel for
    for(long int i = 0; i < chunk_size; i++) {
    """
        loop_postamble = "    } // END LOOP: for i over chunk_size rays"
        body_math = raw_c_code

    core_math = rf"""
    //==========================================
    // MACRO DEFINITIONS FOR BUNDLE ACCESS
    //==========================================
    // IDX_F maps a component to the flattened state bundle using SoA layout.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id)) // Computes the 1D index for the state bundle.
    // IDX_METRIC maps a component to the flattened symmetric metric bundle.
    #define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id)) // Computes the 1D index for the metric bundle.
    // {geometry_index_comment}
    #define {geometry_index_macro}(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
    // IDX_K maps a stage and component triplet to the flattened derivative bundle.
    #define IDX_K(s, c, ray_id) (((s) - 1) * 9 * BUNDLE_CAPACITY + (c) * BUNDLE_CAPACITY + (ray_id)) // Computes the 1D index for the RKF45 derivative bundle.

    {preamble_unpacking_str}

    //==========================================
    // GEODESIC RHS EVALUATION
    //==========================================
    // Local register declarations to capture the evaluated derivatives $\dot{{f}}$.
    double k_out_0, k_out_1, k_out_2, k_out_3, k_out_4, k_out_5, k_out_6, k_out_7, k_out_8; // Thread-local registers allocate memory for the $9$ derivative outputs.

    // Evaluate the derivatives $dx^{{\mu}}/d\lambda$ and $dp^{{\mu}}/d\lambda$ using hardware FMA instructions.
    {body_math}

    //==========================================
    // GLOBAL MEMORY WRITE
    //==========================================
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

    //==========================================
    // MACRO CLEANUP
    //==========================================
    #undef IDX_F
    #undef IDX_METRIC
    #undef {geometry_index_macro}
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

    desc_lines = [
        r"""Orchestrates the memory kernel for computing photon geodesic ODE right-hand sides.

@param d_f_temp_bundle Pointer to the intermediate state bundle $f^{\mu}$ in memory.
@param d_metric_bundle Pointer to the pre-calculated metric bundle $g_{\mu\nu}$ in memory.""",
        f"@param {geometry_bundle_arg_name} Pointer to the pre-calculated "
        f"{geometry_bundle_description} in memory.",
    ]
    if normalized_eom:
        desc_lines.extend(
            [
                "@param d_integration_param_bundle Pointer to the base coordinate-time bundle.",
                "@param d_h Pointer to the per-ray RKF45 step-size bundle.",
            ]
        )
    desc_lines.extend(
        [
            "@param d_k_bundle Pointer to the flattened derivative bundle.",
            "@param stage Current RKF45 stage index used to select the write offset.",
            "@param chunk_size Number of active rays in the bundle batch.",
            "@param stream_idx Active execution stream identifier.",
        ]
    )
    desc = "\n".join(desc_lines)

    cfunc_type = "void"

    name = "calculate_ode_rhs_kernel"

    params_list = [
        "const double *restrict d_f_temp_bundle",
        "const double *restrict d_metric_bundle",
        f"const double *restrict {geometry_bundle_arg_name}",
    ]
    if normalized_eom:
        params_list.extend(
            [
                "const double *restrict d_integration_param_bundle",
                "const double *restrict d_h",
            ]
        )
    params_list.extend(
        [
            "double *restrict d_k_bundle",
            "const int stage",
            "const long int chunk_size",
            "const int stream_idx",
        ]
    )
    params = ", ".join(params_list)

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
