"""
Generate the C function for the Kerr-Schild Christoffel symbols.

Author: Dalton J. Moone
"""
from typing import List
import sympy as sp
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.indexedexp as ixp


def con_kerr_schild(Gamma4UDD_expressions: List[List[List[sp.Expr]]]) -> None:
    """
    Generate and register the C worker for the Kerr-Schild Christoffel symbols.

    Args:
        Gamma4UDD_expressions: A rank-3 list of lists of lists containing the
                               symbolic expressions for the Christoffel symbols.
    """
    name = "con_kerr_schild"
    includes = ["BHaH_defines.h"]
    desc = "Computes the 40 unique Christoffel symbols for the Kerr metric in Kerr-Schild coords."
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const double y[4],
                  connection_struct *restrict conn"""

    list_of_Gamma_syms = [
        Gamma4UDD_expressions[i][j][k]
        for i in range(4)
        for j in range(4)
        for k in range(j, 4)
    ]
    conn_Gamma4UDD = ixp.declarerank3("conn->Gamma4UDD", dimension=4)
    list_of_Gamma_C_vars = [
        str(conn_Gamma4UDD[i][j][k])
        for i in range(4)
        for j in range(4)
        for k in range(j, 4)
    ]

    body = ccg.c_codegen(
        list_of_Gamma_syms,
        list_of_Gamma_C_vars,
        enable_cse=True,
        cse_varprefix="conn",
    )

    cfc.register_CFunction(
        name=name,
        includes=includes,
        desc=desc,
        params=params,
        body=body,
        include_CodeParameters_h=True,
    )