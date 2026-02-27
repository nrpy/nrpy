"""
Register C function for computing massive particle geodesic ODE right-hand sides.

This module registers the 'calculate_ode_rhs_massive' C function, which evaluates
the right-hand sides (RHS) of the coupled system of 8 first-order ODEs governing
massive particle geodesics.

It utilizes the flattened SoA architecture via local batch indexing (j) and
supports GPU offloading via OpenMP.

The system is derived from:
    d(x^mu)/d(tau) = u^mu
    d(u^mu)/d(tau) = -Gamma^mu_{alpha beta} u^alpha u^beta

Author: Dalton J. Moone
"""

import logging
import os
import sys
from typing import List

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc


def calculate_ode_rhs_massive(
    geodesic_rhs_expressions: List[sp.Expr], coordinate_symbols: List[sp.Symbol]
) -> None:
    """
    Generate the C engine for the massive particle geodesic ODE right-hand sides.

    :param geodesic_rhs_expressions: List of 8 symbolic expressions for the RHS.
    :param coordinate_symbols: List of 4 symbolic coordinates (e.g., [t, x, y, z]).
    """
    # Step 1: Define metadata and GPU pragmas
    includes = ["BHaH_defines.h"]
    desc = """@brief Portable GPU-ready derivatives for the massive geodesic ODE system.

    Calculates dx^mu/dtau and du^mu/dtau using the provided state vector and
    pre-computed Christoffel symbols. Uses SoA indexing for SIMD/GPU efficiency.

    Input:
        f_temp: Local state vector [t, x, y, z, u^t, u^x, u^y, u^z]
        metric_g4DD: Flattened metric array.
        conn_GammaUDD: Flattened Christoffel symbol array.
        bundle_capacity: Total capacity of the SIMD/GPU bundle.
        stage: RK stage index.
        j: Local index within the bundle.
    Output:
        k_array: The computed derivatives stored in SoA format."""

    name = "calculate_ode_rhs_massive"
    cfunc_type = "void"

    params = """
                  const double *restrict f_temp,
                  const double *restrict metric_g4DD,
                  const double *restrict conn_GammaUDD,
                  double *restrict k_array,
                  const int bundle_capacity,
                  const int stage,
                  const int j"""

    # Step 2: Analyze symbols used in the expressions to avoid unnecessary unpacking
    used_symbol_names = {
        str(sym) for expr in geodesic_rhs_expressions for sym in expr.free_symbols
    }

    # Step 3: Generate the Dynamic Preamble
    preamble_lines = ["// Step 1: Unpack coordinates from state vector [t, x, y, z]"]
    for i, sym in enumerate(coordinate_symbols):
        if str(sym) in used_symbol_names:
            preamble_lines.append(
                f"const double {str(sym)} = f_temp[IDX_LOCAL({i}, j, bundle_capacity)];"
            )

    preamble_lines.append(
        "\n  // Step 2: Unpack 4-velocity components [u^0, u^1, u^2, u^3]"
    )
    for i in range(4):
        if f"uU{i}" in used_symbol_names:
            preamble_lines.append(
                f"const double uU{i} = f_temp[IDX_LOCAL({i+4}, j, bundle_capacity)];"
            )

    preamble_lines.append("\n  // Step 3: Unpack metric components")
    curr_idx = 0
    for m in range(4):
        for n in range(m, 4):
            comp_name = f"metric_g4DD{m}{n}"
            if comp_name in used_symbol_names:
                preamble_lines.append(
                    f"const double {comp_name} = metric_g4DD[IDX_LOCAL({curr_idx}, j, bundle_capacity)];"
                )
            curr_idx += 1

    preamble_lines.append("\n  // Step 4: Unpack Christoffel symbols (Connections)")
    curr_idx = 0
    for a in range(4):
        for m in range(4):
            for n in range(m, 4):
                comp_name = f"conn_Gamma4UDD{a}{m}{n}"
                if comp_name in used_symbol_names:
                    preamble_lines.append(
                        f"const double {comp_name} = conn_GammaUDD[IDX_LOCAL({curr_idx}, j, bundle_capacity)];"
                    )
                curr_idx += 1

    # Step 4: Map SymPy RHS results to the flattened k_array
    # massive system has 8 components (4 position derivatives, 4 velocity derivatives)
    k_array_outputs = [
        f"k_array[IDX_LOCAL((stage - 1) * 8 + {i}, j, bundle_capacity)]"
        for i in range(8)
    ]

    # Step 5: Generate Calculation Body
    print(f" -> Generating C worker function: {name} ...")
    kernel = ccg.c_codegen(
        geodesic_rhs_expressions, k_array_outputs, enable_cse=True, include_braces=False
    )

    body = "\n  ".join(preamble_lines) + "\n\n" + kernel

    # Step 6: Register C Function with GPU target declarations
    cfc.register_CFunction(
        prefunc="#ifdef USE_GPU\n#pragma omp declare target\n#endif",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=False,
        postfunc="#ifdef USE_GPU\n#pragma omp end declare target\n#endif",
    )
    print(f"    ... {name}() registration complete.")


if __name__ == "__main__":
    # Ensure local modules can be imported
    sys.path.append(os.getcwd())

    # Configure logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("TestCalculateODERHS")

    try:
        from nrpy.equations.general_relativity.geodesics.geodesics import (
            Geodesic_Equations,
        )
    except ImportError:
        logger.error("Error: Could not import Geodesic_Equations for testing.")
        sys.exit(1)

    SPACETIME = "KerrSchild_Cartesian"
    GEO_KEY = f"{SPACETIME}_massive"

    logger.info("Test: Generating ODE RHS C-code for %s...", SPACETIME)

    try:
        # 1. Acquire Symbolic Data
        logger.info(" -> Acquiring symbolic geodesic RHS expressions...")
        geodesic_data = Geodesic_Equations[GEO_KEY]

        # 2. Run the Generator
        calculate_ode_rhs_massive(geodesic_data.geodesic_rhs, geodesic_data.xx)

        # 3. Validation and Output
        func_name = "calculate_ode_rhs_massive"
        if func_name in cfc.CFunction_dict:
            filename = f"{func_name}.c"
            with open(filename, "w", encoding="utf-8") as f:
                f.write(cfc.CFunction_dict[func_name].full_function)
            logger.info(" -> PASS: Wrote %s", filename)
        else:
            raise RuntimeError(
                f"Function {func_name} not registered in CFunction_dict."
            )

    except RuntimeError as e:
        logger.error(" -> FAIL: Test failed with error: %s", e)
        import traceback

        traceback.print_exc()
        sys.exit(1)
