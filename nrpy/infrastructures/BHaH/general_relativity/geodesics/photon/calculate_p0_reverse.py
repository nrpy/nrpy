"""
Generates the C function to compute the time component of the reverse 4-momentum.

Author: Dalton J. Moone
"""

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import sympy as sp


def calculate_p0_reverse(p0_expr: sp.Expr) -> None:
    """
    Generate and register the C function to compute p^0 from the null condition.

    This function takes the pre-computed symbolic expression for the time component
    of the reverse-time 4-momentum (p^0), which solves the null geodesic
    condition g_μν p^μ p^ν = 0, and generates the corresponding C function.

    Args:
        p0_expr: The symbolic sympy expression for p^0.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = ["BHaH_defines.h"]
    desc = r"""@brief Computes reverse-time p^0 from the null condition g_munu p^mu p^nu = 0.
    @param[in]  metric The metric components at the evaluation point.
    @param[in]  y      The 9-component state vector (only momentum is used).
    @return The calculated value of p^0.
    """
    name = "calculate_p0_reverse"
    cfunc_type = "double"
    params = "const metric_struct *restrict metric, const double y[9]"

    # The core logic: call c_codegen to get the C code string.
    # Note the use of keyword arguments as required by the standard.
    #
    # Step 1: Implement the null geodesic condition for p^0.
    # The null condition g_μν p^μ p^ν = 0 is a quadratic equation for p^0.
    # The physically correct (negative) root for reverse-time integration is:
    # p^0 = [-(g_0i p^i) + sqrt((g_0i p^i)^2 - g_00 (g_ij p^i p^j))] / g_00
    p0_C_code = ccg.c_codegen(
        sympyexpr=p0_expr,
        output_varname_str="double p0_val",
        enable_cse=True,
        include_braces=False,
        cse_varprefix="p0_intermed",
    )

    # Construct the final function body with a return statement.
    body = f"{{\n{p0_C_code}\nreturn p0_val;\n}}"

    # Register the C function.
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )