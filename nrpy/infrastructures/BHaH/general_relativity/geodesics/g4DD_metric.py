"""
Register C function for computing metric components.

This module registers the 'g4DD_metric_{spacetime_name}' C function, which calculates
the 10 unique components of the metric tensor for a specific spacetime.
It also registers the required 'metric_struct' in BHaH_defines.h using
the specific infrastructure hook "after_general".

It generates a preamble to unpack f[8] -> coordinates and commondata -> parameters.

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


def g4DD_metric(g4DD_exprs: List[List[sp.Expr]], spacetime_name: str) -> None:
    """
    Generate and register the C function to compute metric components.

    Registers the 'metric_struct' and the 'g4DD_metric' C function.
    The C function computes the upper-triangular components of the metric
    tensor (g_munu) based on the provided symbolic expressions.

    :param g4DD_exprs: A 4x4 list of SymPy expressions representing the metric.
    :param spacetime_name: Name of the spacetime (used for documentation).
    """
    # Step 1: Register metric_struct in BHaH_defines.h
    # The struct stores the 10 unique components of the symmetric metric tensor.
    # We use list comprehension to generate components: g4DD00, g4DD01, g4DD02, g4DD03, g4DD11...
    metric_components = [f"g4DD{i}{j}" for i in range(4) for j in range(i, 4)]
    metric_struct_str = (
        "typedef struct { double " + ", ".join(metric_components) + "; } metric_struct;"
    )

    # Register under the "after_general" section of BHaH_defines.h
    Bdefines_h.register_BHaH_defines("after_general", metric_struct_str)

    # Step 2: Define C function metadata
    includes = ["BHaH_defines.h"]
    desc = f"@brief Computes the 10 unique components of the {spacetime_name} metric.\n"
    name = f"g4DD_metric_{spacetime_name}"
    params = (
        "const commondata_struct *restrict commondata, "
        "const double f[8], "
        "metric_struct *restrict metric"
    )

    # Step 3: Prepare symbolic expressions and C variable names
    # We only compute the upper triangle (mu <= nu) due to symmetry.
    list_of_g4DD_syms = []
    list_of_g4DD_C_vars = []

    for mu in range(4):
        for nu in range(mu, 4):
            list_of_g4DD_syms.append(g4DD_exprs[mu][nu])
            list_of_g4DD_C_vars.append(f"metric->g4DD{mu}{nu}")

    # Step 4: Generate the Dynamic Preamble
    # This unpacks y[] -> (t,x,y,z) and commondata -> (a_spin, M_scale, etc)

    # 4.a: Analyze symbols used in the expressions
    used_symbol_names: Set[str] = set()
    for expr in list_of_g4DD_syms:
        for sym in expr.free_symbols:
            used_symbol_names.add(str(sym))

    # 4.b: Retrieve coordinate symbols from the Analytic Spacetime registry
    # Note: Removed "if not in" check because Analytic_Spacetimes uses lazy loading.
    # Accessing it directly triggers the load/validation.
    xx_symbols = Analytic_Spacetimes[spacetime_name].xx
    preamble_lines = ["// Unpack position coordinates from f[0]..f[3]"]

    for i, symbol in enumerate(xx_symbols):
        sym_name = str(symbol)
        # Check if symbol is used before unpacking to avoid unused variable warnings
        if sym_name in used_symbol_names:
            preamble_lines.append(f"const double {sym_name} = f[{i}];")

    # Use 2 spaces for indentation to match existing style
    preamble = "\n  ".join(preamble_lines)

    # Step 5: Generate C body
    print(f" -> Generating C worker function: {name} (Spacetime: {spacetime_name})...")

    # Generate kernel without braces to combine with preamble
    kernel = ccg.c_codegen(
        list_of_g4DD_syms,
        list_of_g4DD_C_vars,
        enable_cse=True,
        verbose=False,
        include_braces=False,
    )

    # Combine Preamble + Kernel
    body = preamble + "\n\n" + kernel

    # Step 6: Register the C function
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    print(f"    ... {name}() registration complete.")


if __name__ == "__main__":
    import logging
    import os
    import sys

    # Ensure local modules can be imported
    sys.path.append(os.getcwd())

    # Configure logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("TestG4DDAnalytic")

    SPACETIME = "KerrSchild_Cartesian"
    logger.info("Test: Generating Metric C-code for %s...", SPACETIME)

    try:
        # 1. Acquire Symbolic Data
        logger.info(" -> Acquiring symbolic metric expressions...")
        metric_data = Analytic_Spacetimes[SPACETIME]

        # 2. Run the Generator
        logger.info(" -> Calling g4DD_metric()...")
        g4DD_metric(metric_data.g4DD, SPACETIME)

        # 3. Validation
        # Check C Function Registration
        cfunc_name = f"g4DD_metric_{SPACETIME}"
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
        if "typedef struct { double g4DD00" not in defines_str:
            raise RuntimeError(
                "FAIL: metric_struct definition not found in BHaH_defines."
            )
        logger.info(" -> PASS: metric_struct registered successfully in BHaH_defines.")

        # 4. Output Files to Current Directory
        logger.info(" -> Writing C files to current directory...")

        # A. Generate BHaH_defines.h
        # This file contains the typedef for 'metric_struct' and other core definitions.
        Bdefines_h.output_BHaH_defines_h(project_dir=".")
        logger.info("    ... Wrote BHaH_defines.h (contains metric_struct definition)")

        # B. Generate the C function files manually
        # Iterates through registered functions and writes them to .c files.
        for func_name, c_function in cfc.CFunction_dict.items():
            filename = f"{func_name}.c"
            with open(filename, "w", encoding="utf-8") as f:
                f.write(c_function.full_function)
            logger.info("    ... Wrote %s", filename)

        logger.info(" -> Success! All files generated.")

    except Exception as e:  # pylint: disable=broad-exception-caught
        logger.error(" -> FAIL: g4dd_analytic test failed with error: %s", e)
        sys.exit(1)
