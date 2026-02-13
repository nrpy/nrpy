"""
Register C function for computing the initial time component of 4-momentum.

This module registers the 'p0_reverse' C function. It enforces
the 4-momentum normalization constraint for photons (p.p = 0) by solving
the quadratic Hamiltonian constraint for the time component p^0.

It generates a preamble to unpack the state vector f[9] into spatial momentum components (p^i) required for the calculation.

Author: Dalton J. Moone
"""

import logging

# Step 0.a: Import standard Python modules
import sys

# Step 0.b: Import third-party modules
import sympy as sp

# Step 0.c: Import NRPy core modules
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc


def p0_reverse(p0_expr: sp.Expr) -> None:
    """
    Generate and register the C function to compute p^0 for a photon.

    :param p0_expr: The SymPy expression for p^0.
    """
    # Step 3: Define C function metadata
    includes = ["BHaH_defines.h"]
    desc = """@brief Computes the initial time-component of the 4-momentum (p^0).

        Solves the quadratic Hamiltonian constraint equation:
            g_munu p^mu p^nu = 0
        for the negative root of p^0, given the spatial momentum components.

        Input:
            metric: The metric tensor components at the current location.
            f[9]: The state vector (specifically spatial momentum f[5]..f[7]).
        Output:
            p0_out: The computed p^0 component."""
    name = "p0_reverse"

    params = (
        "const metric_struct *restrict metric, "
        "const double f[9], "
        "double *restrict p0_out"
    )

    # Step 4: Generate C body
    print(f" -> Generating C worker function: {name} ...")

    # 4a. Generate the Math Body (using CSE)
    body_math = ccg.c_codegen(
        [p0_expr], ["*p0_out"], enable_cse=True, verbose=False, include_braces=False
    )

    # 4b. Generate the Dynamic Preamble
    preamble_lines = ["// Unpack spatial momentum components from f[5]..f[7]"]
    preamble_lines.append(
        "// Note: f[4] is p^0 (which we are computing), so we skip it."
    )

    # Standard 3-momentum components usually map to pU1, pU2, pU3
    # We map f[5]->pU1, f[6]->pU2, f[7]->pU3
    preamble_lines.append("const double pU1 = f[5];")
    preamble_lines.append("const double pU2 = f[6];")
    preamble_lines.append("const double pU3 = f[7];")
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
        include_CodeParameters_h=False,
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
    logger = logging.getLogger("Testp0ReverseAnalytic")

    SPACETIME = "KerrSchild_Cartesian"
    GEO_KEY = f"{SPACETIME}_photon"

    logger.info("Test: Generating p0_reverse C-code for %s...", SPACETIME)

    try:
        # 1. Acquire Symbolic Data
        logger.info(" -> Acquiring symbolic Hamiltonian constraint (p^0)...")
        geodesic_data = Geodesic_Equations[GEO_KEY]

        # 2. Acquire Coordinate Symbols
        logger.info(" -> Acquiring coordinate symbols from AnalyticSpacetimes...")
        spacetime_data = Analytic_Spacetimes[SPACETIME]

        if geodesic_data.p0_photon is None:
            raise ValueError(f"p0_reverse is None for key {GEO_KEY}")

        # 3. Run the Generator
        logger.info(" -> Calling p0_reverse()...")
        p0_reverse(geodesic_data.p0_photon)

        # 4. Validation
        cfunc_name = "p0_reverse"

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
        logger.error(" -> FAIL: p0_reverse test failed with error: %s", e)
        import traceback

        traceback.print_exc()
        sys.exit(1)
