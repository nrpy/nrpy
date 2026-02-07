"""
Register C function for computing massive particle geodesic ODE right-hand sides.

This module registers the 'calculate_ode_rhs_massive' C function, which evaluates
the right-hand sides (RHS) of the coupled system of 8 first-order ODEs governing
massive particle geodesics.

It generates a C-code preamble to unpack the state vector f[8] into local
coordinate (x^mu) and 4-velocity (u^mu) variables before computing the derivatives.

The system is derived from:
    d(u^mu)/d(tau) = -Gamma^mu_{alpha beta} u^alpha u^beta
    d(x^mu)/d(tau) = u^mu

Author: Dalton J. Moone
"""

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
    :param coordinate_symbols: List of 4 symbolic coordinates (e.g., [t, x, y, z]) mapping to f[0]..f[3].
    """
    # 1. Define metadata
    includes = ["BHaH_defines.h"]
    desc = (
        "@brief Computes the 8 derivatives required for the geodesic ODE system.\n"
        "\n"
        "Calculates dx^mu/dtau and du^mu/dtau using the provided state vector and "
        "pre-computed Christoffel symbols.\n"
        "\n"
        "Input:\n"
        "    f[8]: State vector [t, x, y, z, u^t, u^x, u^y, u^z].\n"
        "    conn: Struct containing analytic Christoffel symbols.\n"
        "Output:\n"
        "    rhs_out[8]: The computed derivatives [dt, dx, dy, dz, du^t, du^x, du^y, du^z]."
    )
    name = "calculate_ode_rhs_massive"

    c_params_str = """
                  const double f[8],
                  const connection_struct *restrict conn,
                  double rhs_out[8]"""

    # 2. Analyze symbols used in the expressions
    used_symbol_names = set()
    for expr in geodesic_rhs_expressions:
        for sym in expr.free_symbols:
            used_symbol_names.add(str(sym))

    preamble_lines = []

    # 3. Unpack y[0]..y[3] (Coordinates) using the provided coordinate_symbols list
    # We iterate through the list to map the specific spacetime coordinates
    # to the state vector indices 0-3.
    preamble_lines.append("// Unpack position coordinates from f[0]..f[3]")
    for i, sym in enumerate(coordinate_symbols):
        sym_name = str(sym)
        # Only unpack if the coordinate is actually used in the RHS expressions
        if sym_name in used_symbol_names:
            preamble_lines.append(f"const double {sym_name} = f[{i}];")

    if preamble_lines:
        preamble_lines.append("")

    # 4. Unpack y[4]..y[7] (4-Velocity)
    # We check for standard velocity names used in GeodesicEquations (massive)
    preamble_lines.append("// Unpack 4-velocity components from f[4]..f[7]")
    for i in range(4):
        coord_name = f"uU{i}"
        if coord_name in used_symbol_names:
            preamble_lines.append(f"const double {coord_name} = f[{i+4}];")

    preamble = "\n  ".join(preamble_lines)

    # 6. Generate Calculation Body
    # Use include_braces=False to cleanly prepend our preamble
    kernel = ccg.c_codegen(
        geodesic_rhs_expressions,
        [f"rhs_out[{i}]" for i in range(8)],
        enable_cse=True,
        verbose=False,
        include_braces=False,
    )

    body = preamble + "\n\n" + kernel

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=c_params_str,
        include_CodeParameters_h=False,
        body=body,
    )


if __name__ == "__main__":
    import logging
    import os
    import sys

    # Ensure local modules can be imported
    sys.path.append(os.getcwd())

    try:
        from nrpy.equations.general_relativity.geodesics.geodesics import (
            Geodesic_Equations,
        )
    except ImportError:
        print("Error: Could not import Geodesic_Equations for testing.")
        sys.exit(1)

    # Configure logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("TestCalculateODERHS")

    SPACETIME = "KerrSchild_Cartesian"
    GEO_KEY = f"{SPACETIME}_massive"

    logger.info("Test: Generating ODE RHS C-code...")

    try:
        # 1. Acquire Symbolic Data
        geodesic_data = Geodesic_Equations[GEO_KEY]

        # 2. Run the Generator (Passing both RHS expressions and Coordinate list)
        calculate_ode_rhs_massive(geodesic_data.geodesic_rhs, geodesic_data.xx)

        # 3. Output
        func_name = "calculate_ode_rhs_massive"
        if func_name in cfc.CFunction_dict:
            filename = f"{func_name}.c"
            with open(filename, "w", encoding="utf-8") as f:
                f.write(cfc.CFunction_dict[func_name].full_function)
            logger.info(" -> Success! Wrote %s", filename)
        else:
            raise RuntimeError(f"Function {func_name} not registered.")

    except Exception as e:  # pylint: disable=broad-exception-caught
        logger.error("Test failed: %s", e)
        sys.exit(1)
