"""
Register C function for computing analytic Christoffel symbols.

This module registers the 'connections_{spacetime_name}' C function, which calculates
the 40 unique components of the Christoffel symbols (Gamma^alpha_mu_nu)
for a specific spacetime using a flattened SoA architecture.

It generates a preamble to unpack coordinate variables directly from the flattened 
state vector using macro indexing. 

Author: Dalton J. Moone
"""

# Step 0.a: Import standard Python modules
from typing import List, Set

# Step 0.b: Import third-party modules
import sympy as sp

# Step 0.c: Import NRPy core modules
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par
from nrpy.equations.general_relativity.geodesics.analytic_spacetimes import (
    Analytic_Spacetimes,
)


def connections(
    Gamma4UDD_exprs: List[List[List[sp.Expr]]], spacetime_name: str, PARTICLE: str
) -> None:
    """
    Generate and register the C function to compute Christoffel symbols.

    The C function computes the unique components of the Christoffel symbols
    (Gamma^alpha_{mu nu}), exploiting symmetry in the lower indices (nu >= mu),
    and writes them to a flattened 1D batch array.

    :param Gamma4UDD_exprs: A 4x4x4 list of SymPy expressions for connections.
    :param spacetime_name: Name of the spacetime (used for documentation).
    :param PARTICLE: The type of particle ("massive" or "photon").
                     Determines state vector documentation logic.
    :raises ValueError: If PARTICLE is not "massive" or "photon".
    """
    # Step 1: Specific setup based on particle type
    if PARTICLE == "massive":
        array_size = 8
    elif PARTICLE == "photon":
        array_size = 9
    else:
        raise ValueError(f"Unsupported PARTICLE: {PARTICLE}")

    # Step 2: Prepare symbolic expressions and 1D flattened C variable names
    list_of_Gamma_syms = []
    list_of_Gamma_C_vars = []

    k = 0
    for alpha in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                list_of_Gamma_syms.append(Gamma4UDD_exprs[alpha][mu][nu])
                # Map 3D tensor components to 1D flat array using macro
                list_of_Gamma_C_vars.append(f"conn_GammaUDD[IDX_LOCAL({k}, batch_id, batch_size)]")
                k += 1

    # Step 3: Generate the Dynamic Preamble
    # 3.a: Analyze symbols used in the expressions first
    used_symbol_names: Set[str] = set()
    for expr in list_of_Gamma_syms:
        for sym in expr.free_symbols:
            used_symbol_names.add(str(sym))

    # 3.b: Retrieve coordinate symbols from the Analytic Spacetime registry
    xx_symbols = Analytic_Spacetimes[spacetime_name].xx
    preamble_lines = [
        f"// Unpack position coordinates from flattened state vector (State vector size: {array_size})"
    ]

    for i, symbol in enumerate(xx_symbols):
        sym_name = str(symbol)
        # Check if symbol is used before unpacking to avoid unused variable warnings
        if sym_name in used_symbol_names:
            preamble_lines.append(f"const double {sym_name} = f[IDX_LOCAL({i}, batch_id, batch_size)];")
    preamble = "\n  ".join(preamble_lines)

    # Step 4: Define C function metadata
    includes = ["BHaH_defines.h"]
    desc = (
        f"@brief Computes the 40 unique Christoffel symbols for the {spacetime_name} metric "
        f"for a {PARTICLE} particle.\n"
    )
    name = f"connections_{spacetime_name}"
    
    # Updated SoA compatible signature
    params = (
        "const commondata_struct *restrict commondata, "
        "const double *restrict f, "
        "double *restrict conn_GammaUDD, "
        "const int batch_size, "
        "const int batch_id"
    )

    # Step 5: Generate C Body
    print(
        f" -> Generating C worker function: {name} (Spacetime: {spacetime_name}, Particle: {PARTICLE})..."
    )

    # Generate kernel without braces so we can wrap preamble + kernel together
    kernel = ccg.c_codegen(
        list_of_Gamma_syms,
        list_of_Gamma_C_vars,
        enable_cse=True,
        verbose=False,
        include_braces=False,
    )

    # Wrap the entire body so the GPU can call it from the ODE RHS kernel
    body = f"""
        #ifdef USE_GPU
        #pragma omp declare target
        #endif
        {preamble}

        {kernel}
        #ifdef USE_GPU
        #pragma omp end declare target
        #endif
        """

    # Step 6: Register the C function
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=True,
    )
    print(f"    ... {name}() registration complete.")


if __name__ == "__main__":
    import logging
    import os
    import sys

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
    logger = logging.getLogger("TestConnectionsAnalytic")

    SPACETIME = "KerrSchild_Cartesian"
    PARTICLE_test = "massive"  # Can be "massive" or "photon"
    GEO_KEY = f"{SPACETIME}_{PARTICLE_test}"

    logger.info(
        "Test: Generating Connections C-code for %s (%s)...", SPACETIME, PARTICLE_test
    )

    try:
        # 1. Acquire Symbolic Data
        logger.info(" -> Acquiring symbolic Christoffel symbols...")
        geodesic_data = Geodesic_Equations[GEO_KEY]

        # 2. Run the Generator
        logger.info(" -> Calling connections()...")
        connections(geodesic_data.Gamma4UDD, SPACETIME, PARTICLE_test)

        # 3. Validation
        cfunc_name = f"connections_{SPACETIME}"

        # Check C Function Registration
        if cfunc_name not in cfc.CFunction_dict:
            raise RuntimeError(
                f"FAIL: '{cfunc_name}' was not registered in cfc.CFunction_dict."
            )

        logger.info(" -> PASS: '%s' function registered successfully.", cfunc_name)

        # 4. Output Files to Current Directory
        Bdefines_h.output_BHaH_defines_h(project_dir=".")
        for func_name, c_function in cfc.CFunction_dict.items():
            filename = f"{func_name}.c"
            with open(filename, "w", encoding="utf-8") as f:
                f.write(c_function.full_function)
            logger.info("    ... Wrote %s", filename)

        logger.info(" -> Success! All files generated.")

    except (ImportError, RuntimeError, ValueError) as e:
        logger.error(" -> FAIL: connections_analytic test failed with error: %s", e)
        import traceback

        traceback.print_exc()
        sys.exit(1)
