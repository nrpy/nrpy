"""
Generate the C helper for computing u^t and u^phi from a radius for a circular orbit.

Author: Dalton J. Moone
"""
import sympy as sp
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc


def calculate_ut_uphi_from_r(ut_expr: sp.Expr, uphi_expr: sp.Expr) -> None:
    """
    Generate the C helper for computing u^t and u^phi from a radius.

    Args:
        ut_expr: The symbolic expression for u^t.
        uphi_expr: The symbolic expression for u^phi.
    """
    name = "calculate_ut_uphi_from_r"
    includes = ["BHaH_defines.h"]
    desc = "Computes u^t and u^phi for a circular orbit at a given radius using a numerically stable method."
    params = """const double r_initial,
                  const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double *ut, double *uphi"""
    body = ccg.c_codegen(
        [ut_expr, uphi_expr], ["*ut", "*uphi"], enable_cse=True, cse_varprefix="ic"
    )

    cfc.register_CFunction(
        name=name,
        includes=includes,
        desc=desc,
        params=params,
        body=body,
        include_CodeParameters_h=True,
    )