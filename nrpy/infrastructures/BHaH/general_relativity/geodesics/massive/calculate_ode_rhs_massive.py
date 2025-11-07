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

    Args:
        geodesic_rhs_expressions: List of 8 symbolic expressions for the RHS.
    """
    name = "calculate_ode_rhs_massive"
    includes = ["BHaH_defines.h"]
    desc = "Calculates the RHS of the 8 massive particle geodesic ODEs."
    params = """const double y[8],
                  const connection_struct *restrict conn,
                  double rhs_out[8]"""
    body = ccg.c_codegen(
        geodesic_rhs_expressions,
        [f"rhs_out[{i}]" for i in range(8)],
        enable_cse=True,
    )

    cfc.register_CFunction(
        name=name,
        includes=includes,
        desc=desc,
        params=params,
        body=body,
    )