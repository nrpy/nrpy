"""
Generates the C function to check conserved quantities for photon geodesics.

Author: Dalton J. Moone
"""

from typing import List
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import sympy as sp


def check_conservation(
    kerr_expressions: List[sp.Expr], schw_expressions: List[sp.Expr]
) -> None:
    """
    Generate and register the C function to check conserved quantities.

    This function creates a C dispatcher that computes Energy (E), the three
    components of angular momentum (L_i), and the Carter Constant (Q) for a
    given state vector. It generates two separate C code blocks, one for the
    Kerr metric and one for the Schwarzschild metric (where Q simplifies to L^2),
    and uses a switch statement to execute the correct one at runtime.

    Args:
        kerr_expressions: List of 5 sympy expressions [E, Lx, Ly, Lz, Q_kerr].
        schw_expressions: List of 5 sympy expressions [E, Lx, Ly, Lz, Q_schw].
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "stdlib.h"]
    desc = r"""@brief Computes conserved quantities (E, L_i, Q/L^2) for a given state vector.
    @param[in]  commondata        Pointer to commondata struct with runtime parameters.
    @param[in]  params            Pointer to params struct (unused, for signature compatibility).
    @param[in]  metric_params_in  Pointer to the metric_params struct specifying the metric type.
    @param[in]  y                 The 9-component state vector.
    @param[out] E                 Pointer to store the calculated energy.
    @param[out] Lx                Pointer to store the calculated x-angular momentum.
    @param[out] Ly                Pointer to store the calculated y-angular momentum.
    @param[out] Lz                Pointer to store the calculated z-angular momentum.
    @param[out] Q                 Pointer to store the calculated Carter Constant.
    """
    name = "check_conservation"
    params = """const commondata_struct *restrict commondata,
        const params_struct *restrict params,
        const metric_params *restrict metric_params_in,
        const double y[9],
        double *E, double *Lx, double *Ly, double *Lz, double *Q"""

    # Step 1: Implement the formulas for conserved quantities.
    # Energy: E = -p_t = -g_{tμ} p^μ
    # Angular Momentum (z-comp): L_z = p_φ = g_{φμ} p^μ  (or x*p_y - y*p_x in Cartesian)
    # Carter Constant (null): Q = p_θ^2 + cos^2(θ) * (a^2(E^2 - p_t^2) - L_z^2/sin^2(θ))
    # For Schwarzschild (a=0), Q simplifies to L_x^2 + L_y^2 + L_z^2.
    output_vars = ["*E", "*Lx", "*Ly", "*Lz", "*Q"]
    body_C_code_kerr = ccg.c_codegen(
        sympyexpr=kerr_expressions,
        output_varname_str=output_vars,
        enable_cse=True,
        include_braces=False,
        cse_varprefix="kerr_cons_intermed",
    )
    body_C_code_schw = ccg.c_codegen(
        sympyexpr=schw_expressions,
        output_varname_str=output_vars,
        enable_cse=True,
        include_braces=False,
        cse_varprefix="schw_cons_intermed",
    )

    # Construct the final function body, which includes logic to get the metric
    # and a switch statement to dispatch to the correct kernel.
    body = f"""
    // Unpack parameters from commondata struct that are needed symbolically.
    const REAL a_spin = commondata->a_spin;

    // Allocate a metric_struct on the heap to prevent stack overflow in parallel loops.
    metric_struct* metric = (metric_struct*)malloc(sizeof(metric_struct));
    if (metric == NULL) {{
        fprintf(stderr, "Error: Failed to allocate memory for metric_struct in check_conservation.\\n");
        exit(1);
    }}

    // Call the dispatcher to fill the allocated struct with metric components.
    g4DD_metric(commondata, params, metric_params_in, y, metric);

    // Dispatch to the correct conservation law implementation.
    if (metric_params_in->type == Kerr) {{
        {body_C_code_kerr}
    }} else {{ // Both Schwarzschild types use the same Cartesian conservation laws.
        {body_C_code_schw}
    }}

    free(metric);
    """

    # Register the C function.
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type="void",
        name=name,
        params=params,
        body=body,
    )