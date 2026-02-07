"""
Register C function for computing the initial time component of 4-velocity.

This module registers the 'u0_massive_{spacetime_name}' C function. It enforces
the 4-velocity normalization constraint for massive particles (u.u = -1) by solving
the quadratic Hamiltonian constraint for the time component u^0.

It generates a preamble to unpack the state vector f[8] into local coordinates
and spatial velocity components (u^i) required for the calculation.

Author: Dalton J. Moone
"""

import logging

# Step 0.a: Import standard Python modules
import sys
from typing import List

# Step 0.b: Import third-party modules
import sympy as sp

# Step 0.c: Import NRPy core modules
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc


def u0_massive(
    u0_expr: sp.Expr, spacetime_name: str, coord_symbols: List[sp.Symbol]
) -> None:
    """
    Generate and register the C function to compute u^0 for a massive particle.

    :param u0_expr: The SymPy expression for u^0.
    :param spacetime_name: Name of the spacetime (used for function naming).
    :param coord_symbols: A list of SymPy symbols representing the coordinates
                          (e.g., [t, x, y, z]).
    """
    # Step 3: Define C function metadata
    includes = ["BHaH_defines.h"]
    desc = (
        "@brief Computes the initial time-component of the 4-velocity (u^0).\n"
        "\n"
        "Solves the quadratic Hamiltonian constraint equation:\n"
        "    g_munu u^mu u^nu = -1\n"
        "for the positive root of u^0, given the spatial velocity components.\n"
        "\n"
        "Input:\n"
        "    commondata: Simulation parameters (mass, spin, etc.).\n"
        "    metric: The metric tensor components at the current location.\n"
        "    f[8]: The state vector (specifically spatial velocities f[5]..f[7]).\n"
        "Output:\n"
        "    u0_out: The computed u^0 component."
    )
    name = f"u0_massive_{spacetime_name}"

    params = (
        "const commondata_struct *restrict commondata, "
        "const metric_struct *restrict metric, "
        "const double f[8], "
        "double *restrict u0_out"
    )

    # Step 4: Generate C body
    print(f" -> Generating C worker function: {name} (Spacetime: {spacetime_name})...")

    # 4a. Generate the Math Body (using CSE)
    body_math = ccg.c_codegen(
        [u0_expr], ["*u0_out"], enable_cse=True, verbose=False, include_braces=False
    )

    # 4b. Generate the Dynamic Preamble
    preamble_lines = ["// Unpack position coordinates from f[0]..f[3]"]

    # Unpack coordinates: y[0] -> t, y[1] -> x, etc.
    for i, symbol in enumerate(coord_symbols):
        preamble_lines.append(f"const double {str(symbol)} = f[{i}];")

    preamble_lines.append("")
    preamble_lines.append("// Unpack spatial velocity components from f[5]..f[7]")
    preamble_lines.append(
        "// Note: f[4] is u^0 (which we are computing), so we skip it."
    )

    # Standard 3-velocity components usually map to uU1, uU2, uU3
    # We map y[5]->uU1, y[6]->uU2, y[7]->uU3
    preamble_lines.append("const double uU1 = f[5];")
    preamble_lines.append("const double uU2 = f[6];")
    preamble_lines.append("const double uU3 = f[7];")
    preamble_lines.append("")

    preamble = "\n  ".join(preamble_lines)

    # Combine preamble and math
    full_body = preamble + body_math

    # Step 5: Register the C function
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=full_body,
    )
    print(f"    ... {name}() registration complete.")


if __name__ == "__main__":
    import os

    # Ensure local modules can be imported
    sys.path.append(os.getcwd())

    try:
        from nrpy.equations.general_relativity.geodesics.analytic_spacetimes import (
            Analytic_Spacetimes,
        )
        from nrpy.equations.general_relativity.geodesics.geodesics import (
            Geodesic_Equations,
        )
    except ImportError as e:
        print(f"Error: Could not import required modules: {e}")
        sys.exit(1)

    # Configure logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("TestU0MassiveAnalytic")

    SPACETIME = "KerrSchild_Cartesian"
    GEO_KEY = f"{SPACETIME}_massive"

    logger.info("Test: Generating u0_massive C-code for %s...", SPACETIME)

    try:
        # 1. Acquire Symbolic Data
        logger.info(" -> Acquiring symbolic Hamiltonian constraint (u^0)...")
        geodesic_data = Geodesic_Equations[GEO_KEY]

        # 2. Acquire Coordinate Symbols
        logger.info(" -> Acquiring coordinate symbols from AnalyticSpacetimes...")
        spacetime_data = Analytic_Spacetimes[SPACETIME]

        if geodesic_data.u0_massive is None:
            raise ValueError(f"u0_massive is None for key {GEO_KEY}")

        # 3. Run the Generator
        logger.info(" -> Calling u0_massive_analytic()...")
        u0_massive(geodesic_data.u0_massive, SPACETIME, spacetime_data.xx)

        # 4. Validation
        cfunc_name = f"u0_massive_{SPACETIME}"

        if cfunc_name not in cfc.CFunction_dict:
            raise RuntimeError(
                f"FAIL: '{cfunc_name}' was not registered in cfc.CFunction_dict."
            )
        logger.info(" -> PASS: '%s' function registered successfully.", cfunc_name)

        # 5. Output Files
        filename = f"{cfunc_name}.c"
        cfunc = cfc.CFunction_dict[cfunc_name]
        with open(filename, "w", encoding="utf-8") as f:
            f.write(cfunc.full_function)
        logger.info(" -> Written to %s", filename)

    except Exception as e:  # pylint: disable=broad-exception-caught
        logger.error(" -> FAIL: u0_massive_analytic test failed with error: %s", e)
        import traceback

        traceback.print_exc()
        sys.exit(1)
