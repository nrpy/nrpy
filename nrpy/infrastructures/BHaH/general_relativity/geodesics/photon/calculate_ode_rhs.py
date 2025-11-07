"""
Generates the C function for the right-hand side of the geodesic ODEs for photons.

Author:Dalton J. Moone
"""

from typing import List
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import sympy as sp


def calculate_ode_rhs(all_rhs_expressions: List[sp.Expr]) -> None:
    """
    Generate and register the C function for the photon geodesic ODE RHS.

    This function takes the pre-computed symbolic expressions for the right-hand
    sides (RHS) of the 9 geodesic ODEs (4 for position, 4 for momentum, 1 for
    path length) and generates a single, optimized C function to compute them.

    Args:
        all_rhs_expressions: A list of 9 sympy expressions for the RHSs.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = ["BHaH_defines.h"]
    desc = r"""@brief Calculates the right-hand sides (RHS) of the 9 photon geodesic ODEs.

    This function implements the generic geodesic equation using pre-computed
    Christoffel symbols. It is a pure "engine" function that does not depend
    on any specific metric's parameters, only on the geometric values passed
    to it via the metric and connection structs.

    @param[in]  y         The 9-component state vector [t, x, y, z, p^t, p^x, p^y, p^z, L].
    @param[in]  metric    A pointer to the metric_struct holding g_μν.
    @param[in]  conn      A pointer to the connection_struct holding Γ^α_μν.
    @param[out] rhs_out   A pointer to the 9-component output array for the RHS results.
    """
    name = "calculate_ode_rhs"
    params = "const double y[9], const metric_struct *restrict metric, const connection_struct *restrict conn, double rhs_out[9]"

    # The core logic: call c_codegen to generate the body of the function.
    # The output variable names are set to match the C array `rhs_out`.
    #
    # Step 1: Implement the system of 9 first-order ODEs for reverse-time ray-tracing.
    # The system is:
    # 1. Position ODE: dx^α/dκ = p^α
    # 2. Momentum ODE: dp^α/dκ = -Γ^α_μν p^μ p^ν
    # 3. Path Length ODE: dL/dκ = sqrt(γ_ij p^i p^j)
    rhs_output_vars = [f"rhs_out[{i}]" for i in range(9)]
    body = ccg.c_codegen(
        sympyexpr=all_rhs_expressions,
        output_varname_str=rhs_output_vars,
        enable_cse=True,
        cse_varprefix="ode_rhs_intermed",
    )

    # Register the C function.
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )