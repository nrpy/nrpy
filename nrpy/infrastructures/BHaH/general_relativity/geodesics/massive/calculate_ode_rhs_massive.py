"""
Generate the C function for the massive particle geodesic ODE right-hand sides.

Author: Dalton J. Moone
"""
from typing import List
import sympy as sp
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc


def calculate_ode_rhs_massive(geodesic_rhs_expressions: List[sp.Expr]) -> None:
    """
    Generate the C engine for the massive particle geodesic ODE right-hand sides.

    :param geodesic_rhs_expressions: List of 8 symbolic expressions for the RHS.
    """
    name = "calculate_ode_rhs_massive"
    includes = ["BHaH_defines.h"]
    desc = "Calculates the RHS of the 8 massive particle geodesic ODEs."
    params = """const double y[8],
                  const connection_struct *restrict conn,
                  double rhs_out[8]"""

    # Generate Calculation Body directly from expressions
    body = ccg.c_codegen(
        geodesic_rhs_expressions,
        [f"rhs_out[{i}]" for i in range(8)],
        enable_cse=True,
        verbose=False,
    )

    cfc.register_CFunction(
        name=name,
        includes=includes,
        desc=desc,
        params=params,
        body=body,
    )


if __name__ == "__main__":
    import logging
    import sys
    import os

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

        # 2. Run the Generator (Generic, no spacetime string passed)
        calculate_ode_rhs_massive(geodesic_data.geodesic_rhs)

        # 3. Output
        func_name = "calculate_ode_rhs_massive"
        if func_name in cfc.CFunction_dict:
            filename = f"{func_name}.c"
            with open(filename, "w", encoding="utf-8") as f:
                f.write(cfc.CFunction_dict[func_name].full_function)
            logger.info(f" -> Success! Wrote {filename}")
        else:
            raise RuntimeError(f"Function {func_name} not registered.")

    except Exception as e:
        logger.error(f"Test failed: {e}")
        sys.exit(1)
