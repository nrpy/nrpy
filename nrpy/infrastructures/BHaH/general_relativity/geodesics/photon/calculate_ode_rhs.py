"""
Registers C function for computing photon geodesic ODE right-hand sides.

Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
Optimized with a Preamble pattern to minimize register pressure and ensure 
safe hardware-agnostic compilation for global GPU architectures.
"""

from typing import List

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc


def calculate_ode_rhs(
    geodesic_rhs_expressions: List[sp.Expr], coordinate_symbols: List[sp.Symbol]
) -> None:
    """
    Registers the portable GPU-ready derivatives for the geodesic ODE system.

    Implements the preamble pattern to map local registers efficiently before 
    computing the right-hand-side evaluations of the integration step.

    Args:
        geodesic_rhs_expressions (List[sp.Expr]): The mathematical right-hand side evaluations.
        coordinate_symbols (List[sp.Symbol]): The spatial/temporal coordinate variables in order.

    >>> import sympy as sp
    >>> calculate_ode_rhs([], [])
    """
    # 1. Define C-Function metadata in order of appearance
    includes = ["BHaH_defines.h"]
    desc = """@brief Portable GPU-ready derivatives for the geodesic ODE system.
    Evaluates the Christoffel connections and spacetime metric mappings to formulate 
    the right-hand side of the geodesic integration step."""
    cfunc_type = "void"
    name = "calculate_ode_rhs"
    params = """
                  const double *restrict f_temp,
                  const double *restrict metric_g4DD,
                  const double *restrict conn_GammaUDD,
                  double *restrict k_array,
                  const int bundle_capacity,
                  const int stage,
                  const int j"""

    used_symbol_names = {
        str(sym) for expr in geodesic_rhs_expressions for sym in expr.free_symbols
    }

    preamble_lines = [
        "// Ensure hardware-agnostic local array indexing",
        "#ifndef IDX_LOCAL",
        "#define IDX_LOCAL(component, batch_id, batch_size) ((component) * (batch_size) + (batch_id))",
        "#endif\n",
        "// Step 1: Unpack coordinate position data from the 9-component intermediate state vector",
    ]

    # Map the first 4 components logically as spacetime coordinates (t, x, y, z)
    for i, sym in enumerate(coordinate_symbols):
        if str(sym) in used_symbol_names:
            preamble_lines.append(
                f"const double {str(sym)} = f_temp[IDX_LOCAL({i}, j, bundle_capacity)]; // The spacetime coordinate corresponding to {str(sym)}."
            )

    preamble_lines.append(
        "\n  // Step 2: Unpack physical momenta from the intermediate state vector"
    )
    for i in range(4):
        if f"pU{i}" in used_symbol_names:
            # Indices 4, 5, 6, 7 map to p_t, p_x, p_y, p_z respectively
            preamble_lines.append(
                f"const double pU{i} = f_temp[IDX_LOCAL({i+4}, j, bundle_capacity)]; // The momentum component representing contravariant p^{i}."
            )

    preamble_lines.append(
        "\n  // Step 3: Unpack local spacetime metric tensor components"
    )
    curr_idx = 0
    for m in range(4):
        for n in range(m, 4):
            comp_name = f"metric_g4DD{m}{n}"
            if comp_name in used_symbol_names:
                preamble_lines.append(
                    f"const double {comp_name} = metric_g4DD[IDX_LOCAL({curr_idx}, j, bundle_capacity)]; // The symmetric covariant metric tensor component g_{m}{n}."
                )
            curr_idx += 1

    preamble_lines.append(
        "\n  // Step 4: Unpack Christoffel connection coefficients"
    )
    curr_idx = 0
    for a in range(4):
        for m in range(4):
            for n in range(m, 4):
                comp_name = f"conn_Gamma4UDD{a}{m}{n}"
                if comp_name in used_symbol_names:
                    preamble_lines.append(
                        f"const double {comp_name} = conn_GammaUDD[IDX_LOCAL({curr_idx}, j, bundle_capacity)]; // The symmetric Christoffel symbol of the second kind \\Gamma^{a}_{m}{n}."
                    )
                curr_idx += 1

    body = "\n  ".join(preamble_lines) + "\n\n"

    k_array_outputs = [
        f"k_array[IDX_LOCAL((stage - 1) * 9 + {i}, j, bundle_capacity)]"
        for i in range(9)
    ]

    body += ccg.c_codegen(
        geodesic_rhs_expressions,
        k_array_outputs,
        enable_cse=True,
        include_braces=False,
        verbose=False,
    )

    prefunc = """
    #ifdef USE_GPU
    #pragma omp declare target
    #endif
    """
    postfunc = """
    #ifdef USE_GPU
    #pragma omp end declare target
    #endif
    """

    # Register the C function adhering to the exact order expected by the framework.
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        postfunc=postfunc,
    )