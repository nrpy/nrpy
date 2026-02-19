"""
Register C function for computing analytic Christoffel symbols.

This module registers the 'connections_{spacetime_name}' C function, which calculates
the 40 unique components of the Christoffel symbols (Gamma^alpha_mu_nu)
for a specific spacetime. It also registers the required 'connection_struct'.

It generates a preamble to unpack f[0]..f[3] from the state vector into coordinate
variables (e.g., t, x, y, z). Spacetime parameters stored in `commondata` are passed
as a separate argument and are accessed via the infrastructure (e.g., macros in
`BHaH_defines.h`), not by additional unpacking in this preamble.

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

    Registers the 'connection_struct' and the 'connections_{spacetime_name}' C function.
    The C function computes the unique components of the Christoffel symbols
    (Gamma^alpha_{mu nu}), exploiting symmetry in the lower indices (nu >= mu).

    :param Gamma4UDD_exprs: A 4x4x4 list of SymPy expressions for connections.
    :param spacetime_name: Name of the spacetime (used for documentation).
    :param PARTICLE: The type of particle ("massive" or "photon").
                     Determines array size for the state vector f.
    :raises ValueError: If PARTICLE is not "massive" or "photon".
    """
    # Step 1: Specific setup based on particle type
    if PARTICLE == "massive":
        array_size = 8
    elif PARTICLE == "photon":
        array_size = 9
    else:
        raise ValueError(f"Unsupported PARTICLE: {PARTICLE}")

    # Step 2: Register connection_struct in BHaH_defines.h
    # We use list comprehension to generate components: Gamma4UDD000, Gamma4UDD001...
    connection_components = [
        f"Gamma4UDD{i}{j}{k}" for i in range(4) for j in range(4) for k in range(j, 4)
    ]
    connections_struct_str = (
        "typedef struct { double "
        + ", ".join(connection_components)
        + "; } connection_struct;"
    )

    # Register under the "after_general" section of BHaH_defines.h
    Bdefines_h.register_BHaH_defines("after_general", connections_struct_str)

    # Step 3: Prepare symbolic expressions and C variable names
    list_of_Gamma_syms = []
    list_of_Gamma_C_vars = []

    for alpha in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                list_of_Gamma_syms.append(Gamma4UDD_exprs[alpha][mu][nu])
                list_of_Gamma_C_vars.append(f"conn->Gamma4UDD{alpha}{mu}{nu}")

    # Step 4: Generate the Dynamic Preamble
    # 4.a: Analyze symbols used in the expressions first
    used_symbol_names: Set[str] = set()
    for expr in list_of_Gamma_syms:
        for sym in expr.free_symbols:
            used_symbol_names.add(str(sym))

    # 4.b: Retrieve coordinate symbols from the Analytic Spacetime registry
    xx_symbols = Analytic_Spacetimes[spacetime_name].xx
    preamble_lines = [
        f"// Unpack position coordinates from f[0]..f[3] (State vector size: {array_size})"
    ]

    for i, symbol in enumerate(xx_symbols):
        sym_name = str(symbol)
        # Check if symbol is used before unpacking to avoid unused variable warnings
        if sym_name in used_symbol_names:
            preamble_lines.append(f"const double {sym_name} = f[{i}];")
    preamble = "\n  ".join(preamble_lines)

    # Step 5: Define C function metadata
    includes = ["BHaH_defines.h"]
    desc = (
        f"@brief Computes the 40 unique Christoffel symbols for the {spacetime_name} metric "
        f"for a {PARTICLE} particle.\n"
    )
    name = f"connections_{spacetime_name}"
    params = (
        "const commondata_struct *restrict commondata, "
        f"const double f[{array_size}], "
        "connection_struct *restrict conn"
    )

    # Step 6: Generate C Body
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

    # Combine Preamble + Kernel
    body = preamble + "\n\n" + kernel

    # Step 7: Register the C function
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

        # Check Struct Registration
        bhah_defines_dict = par.glb_extras_dict.get("BHaH_defines", {})
        if "after_general" not in bhah_defines_dict:
            raise RuntimeError("FAIL: BHaH_defines_h 'after_general' section is empty.")

        defines_str = bhah_defines_dict["after_general"]
        if "typedef struct { double Gamma4UDD000" not in defines_str:
            raise RuntimeError(
                "FAIL: connection_struct definition not found in BHaH_defines."
            )
        logger.info(
            " -> PASS: connection_struct registered successfully in BHaH_defines."
        )

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
