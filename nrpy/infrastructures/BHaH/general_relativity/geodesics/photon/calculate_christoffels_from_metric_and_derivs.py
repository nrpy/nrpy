"""
Generates the C worker for computing Christoffel symbols algebraically.

Author: Dalton J. Moone
"""

from typing import List
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.indexedexp as ixp
import sympy as sp


def calculate_christoffels_from_metric_and_derivs(Gamma4UDD_num_recipe: List[List[List[sp.Expr]]]) -> None:
    """
    Generate and register the C worker for algebraic Christoffel calculation.

    This function takes the pure symbolic recipe for the Christoffel symbols,
    which is expressed in terms of abstract metric and metric derivative
    placeholders. It generates a C function that takes these values as inputs
    (from structs) and computes the Christoffel symbols. This is a core
    component of the numerical metric pipeline.

    Args:
        Gamma4UDD_num_recipe: The symbolic recipe for the Christoffel symbols.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = ["BHaH_defines.h"]
    desc = r"""@brief Computes Christoffel symbols from pre-interpolated metric and derivative values.

    This is a pure algebraic engine. It takes structs containing the values of
    the metric and its first derivatives at a point and computes the
    Christoffel symbols using the standard formula.

    @param[in]  g4DD_in   Pointer to a struct containing the 10 unique metric components.
    @param[in]  g4DDdD_in Pointer to a struct containing the 40 unique metric derivative components.
    @param[out] conn_out  Pointer to the connection_struct to be filled with the results.
    """
    name = "calculate_christoffels_from_metric_and_derivs"
    params = """const metric_struct *restrict g4DD_in,
                const g4DD_deriv_struct *restrict g4DDdD_in,
                connection_struct *restrict conn_out"""

    # Step 1: Prepare lists of symbolic expressions and C variable names for code generation.
    list_of_Gamma_C_vars: List[str] = []
    list_of_Gamma_syms: List[sp.Expr] = []
    conn_Gamma4UDD = ixp.declarerank3("conn_out->Gamma4UDD", dimension=4)
    for alpha in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                list_of_Gamma_C_vars.append(str(conn_Gamma4UDD[alpha][mu][nu]))
                list_of_Gamma_syms.append(Gamma4UDD_num_recipe[alpha][mu][nu])

    # Step 2: Generate a C code preamble to unpack input structs into local variables.
    # This is crucial because the symbolic recipe was built using these exact local variable names.
    preamble = "// Unpack input structs into local variables that match symbolic recipe.\\n"
    for i in range(4):
        for j in range(i, 4):
            preamble += f"    const double g4DD{i}{j} = g4DD_in->g{i}{j};\\n"

    for i in range(4):
        for j in range(i, 4):
            for k in range(4):
                preamble += f"    const double g4DDdD{i}{j}d{k} = g4DDdD_in->g{i}{j}d{k};\\n"

    # Step 3: Generate the core computational kernel from the symbolic recipe.
    # The formula for the Christoffel symbols of the second kind is:
    # Γ^α_μν = (1/2) g^αδ (g_νδ,μ + g_μδ,ν - g_μν,δ)
    kernel_C_code = ccg.c_codegen(
        sympyexpr=list_of_Gamma_syms,
        output_varname_str=list_of_Gamma_C_vars,
        enable_cse=True,
        cse_varprefix="num_conn_intermed",        
        include_braces=False,
        verbose=False   
    )
   

    # Step 4: Assemble the full function body.
    body = preamble + kernel_C_code

    # Register the C function.
    cfc.register_CFunction(
        includes=includes, desc=desc, name=name, params=params, body=body
    )
    
