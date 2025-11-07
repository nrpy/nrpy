"""
Generate the C function for the Kerr metric in Cartesian Kerr-Schild coordinates.

Author: Dalton J. Moone
"""
from typing import List
import sympy as sp
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc


def g4dd_kerr_schild(g4DD_expressions: List[List[sp.Expr]]) -> None:
    """
    Generate and register the C worker function for the Kerr-Schild metric.

    This function takes the symbolic expressions for the 10 unique components
    of the Kerr-Schild metric tensor and generates an optimized C function
    that computes them at a given position.

    Args:
        g4DD_expressions: A list of lists containing the symbolic expressions
                          for the Kerr-Schild metric components.
    """
    name = "g4DD_kerr_schild"
    includes = ["BHaH_defines.h"]
    desc = "Computes the 10 unique components of the Kerr metric in Cartesian Kerr-Schild coords."
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const double y[4],
                  metric_struct *restrict metric"""

    list_of_g4DD_syms = [g4DD_expressions[i][j] for i in range(4) for j in range(i, 4)]
    list_of_g4DD_C_vars = [f"metric->g{i}{j}" for i in range(4) for j in range(i, 4)]

    body = ccg.c_codegen(
        list_of_g4DD_syms, list_of_g4DD_C_vars, enable_cse=True, cse_varprefix="KS"
    )

    cfc.register_CFunction(
        name=name,
        includes=includes,
        desc=desc,
        params=params,
        body=body,
        include_CodeParameters_h=True,
    )