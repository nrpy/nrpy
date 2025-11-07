"""
Generate the C function for computing conserved quantities (E, L_i, Q).

Author: Dalton J. Moone
"""
from typing import List
import sympy as sp
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc


def check_conservation_massive(
    kerr_E_expr: sp.Expr,
    kerr_L_exprs: List[sp.Expr],
    kerr_Q_expr: sp.Expr,
    schw_E_expr: sp.Expr,
    schw_L_exprs: List[sp.Expr],
    schw_Q_expr: sp.Expr,
) -> None:
    """
    Generate the C engine for computing conserved quantities (E, L_i, Q).

    This engine contains two distinct C code blocks, one for Kerr (using the
    full Carter constant) and one for Schwarzschild (where Q simplifies to L^2),
    and dispatches to the correct one at runtime.

    Args:
        kerr_E_expr: Symbolic expression for energy in Kerr.
        kerr_L_exprs: List of symbolic expressions for [Lx, Ly, Lz] in Kerr.
        kerr_Q_expr: Symbolic expression for the Carter Constant in Kerr.
        schw_E_expr: Symbolic expression for energy in Schwarzschild.
        schw_L_exprs: List of symbolic expressions for [Lx, Ly, Lz] in Schwarzschild.
        schw_Q_expr: Symbolic expression for L^2 in Schwarzschild.
    """
    name = "check_conservation_massive"
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "<stdlib.h>"]
    desc = "Computes conserved quantities (E, L_i, Q) for a given massive particle state vector."
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric_params_in,
                  const double y[8],
                  double *E, double *Lx, double *Ly, double *Lz, double *Q"""

    kerr_expressions = [kerr_E_expr] + kerr_L_exprs + [kerr_Q_expr]
    schw_expressions = [schw_E_expr] + schw_L_exprs + [schw_Q_expr]
    output_vars = ["*E", "*Lx", "*Ly", "*Lz", "*Q"]

    body_kerr = ccg.c_codegen(
        kerr_expressions, output_vars, enable_cse=True, include_braces=False
    )
    body_schw = ccg.c_codegen(
        schw_expressions, output_vars, enable_cse=True, include_braces=False
    )

    body = f"""
    // Allocate a metric struct on the HEAP to get a pointer.
    // This is necessary because the symbolic expressions for conserved
    // quantities were generated assuming pointer access (metric->g...).
    metric_struct* metric = (metric_struct*)malloc(sizeof(metric_struct));
    if (metric == NULL) {{
        fprintf(stderr, "Error: Failed to allocate memory for metric_struct in check_conservation.n");
        exit(1);
    }}
    
    // Call the dispatcher to fill the struct with metric components.
    // We pass the pointer 'metric' directly.
    g4DD_metric(commondata, params, metric_params_in, y, metric);

    // Dispatch to the correct block of code based on metric type.
    if (metric_params_in->type == Kerr) {{
        {body_kerr}
    }} else {{ // Both Schwarzschild and Schwarzschild_Standard
        {body_schw}
    }}
    
    // Free the heap-allocated memory to prevent memory leaks.
    free(metric);
    """

    cfc.register_CFunction(
        name=name,
        includes=includes,
        desc=desc,
        params=params,
        body=body,
        include_CodeParameters_h=True,
    )