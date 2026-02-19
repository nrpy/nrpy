"""
Register C function for computing metric components.

This module registers the 'g4DD_metric_{spacetime_name}' C function,
which calculates the 10 unique components of the metric tensor for a specific spacetime.
It also registers the required 'metric_struct' in BHaH_defines.h using
the specific infrastructure hook "after_general".

It generates a preamble to unpack f[0]..f[3] from the state vector into coordinate
variables. Spacetime parameters stored in `commondata` are provided through the
function arguments and accessed via the BHaH infrastructure, rather than being
explicitly unpacked into local variables in this preamble.

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
from nrpy.equations.general_relativity.geodesics.analytic_spacetimes import (
    Analytic_Spacetimes,
)


def g4DD_metric(
    g4DD_exprs: List[List[sp.Expr]], spacetime_name: str, PARTICLE: str
) -> None:
    """
    Generate and register the C function to compute metric components.

    Registers the 'metric_struct' and the 'g4DD_metric' C function.
    The C function computes the upper-triangular components of the metric
    tensor (g_munu) based on the provided symbolic expressions.

    :param g4DD_exprs: A 4x4 list of SymPy expressions representing the metric.
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

    # Step 2: Register metric_struct in BHaH_defines.h
    # The struct stores the 10 unique components of the symmetric metric tensor.
    metric_components = [f"g4DD{i}{j}" for i in range(4) for j in range(i, 4)]
    metric_struct_str = (
        "typedef struct { double " + ", ".join(metric_components) + "; } metric_struct;"
    )

    # Register under the "after_general" section of BHaH_defines.h
    Bdefines_h.register_BHaH_defines("after_general", metric_struct_str)

    # Step 3: Define C function metadata
    includes = ["BHaH_defines.h"]
    desc = (
        f"@brief Computes the 10 unique components of the {spacetime_name} metric "
        f"for a {PARTICLE} particle.\n"
    )
    name = f"g4DD_metric_{spacetime_name}"
    params = (
        "const commondata_struct *restrict commondata, "
        f"const double f[{array_size}], "
        "metric_struct *restrict metric"
    )

    # Step 4: Prepare symbolic expressions and C variable names
    list_of_g4DD_syms = []
    list_of_g4DD_C_vars = []

    for mu in range(4):
        for nu in range(mu, 4):
            list_of_g4DD_syms.append(g4DD_exprs[mu][nu])
            list_of_g4DD_C_vars.append(f"metric->g4DD{mu}{nu}")

    # Step 5: Generate the Dynamic Preamble
    # 5.a: Analyze symbols used in the expressions
    used_symbol_names: Set[str] = set()
    for expr in list_of_g4DD_syms:
        for sym in expr.free_symbols:
            used_symbol_names.add(str(sym))

    # 5.b: Retrieve coordinate symbols from the Analytic Spacetime registry
    xx_symbols = Analytic_Spacetimes[spacetime_name].xx
    preamble_lines = [
        f"// Unpack position coordinates from f[0]..f[3] (State vector size: {array_size})"
    ]

    for i, symbol in enumerate(xx_symbols):
        sym_name = str(symbol)
        if sym_name in used_symbol_names:
            preamble_lines.append(f"const double {sym_name} = f[{i}];")

    preamble = "\n  ".join(preamble_lines)

    # Step 6: Generate C body
    print(
        f" -> Generating C worker function: {name} (Spacetime: {spacetime_name}, Particle: {PARTICLE})..."
    )

    kernel = ccg.c_codegen(
        list_of_g4DD_syms,
        list_of_g4DD_C_vars,
        enable_cse=True,
        verbose=False,
        include_braces=False,
    )

    body = preamble + "\n\n" + kernel

    # Step 7: Register the C function
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

    sys.path.append(os.getcwd())

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("TestG4DDMetric")

    SPACETIME = "KerrSchild_Cartesian"
    PARTICLE_test = "photon"  # Can be "massive" or "photon"

    logger.info(
        "Test: Generating Metric C-code for %s (%s)...", SPACETIME, PARTICLE_test
    )

    try:
        # 1. Acquire Symbolic Data
        logger.info(" -> Acquiring symbolic metric expressions...")
        metric_data = Analytic_Spacetimes[SPACETIME]

        # 2. Run the Generator
        logger.info(" -> Calling g4DD_metric()...")
        g4DD_metric(metric_data.g4DD, SPACETIME, PARTICLE_test)

        # 3. Validation
        cfunc_name = f"g4DD_metric_{SPACETIME}"
        if cfunc_name not in cfc.CFunction_dict:
            raise RuntimeError(f"FAIL: '{cfunc_name}' was not registered.")

        logger.info(" -> PASS: '%s' function registered successfully.", cfunc_name)

        # 4. Output Files
        Bdefines_h.output_BHaH_defines_h(project_dir=".")
        for func_name, c_function in cfc.CFunction_dict.items():
            filename = f"{func_name}.c"
            with open(filename, "w", encoding="utf-8") as f:
                f.write(c_function.full_function)
            logger.info("    ... Wrote %s", filename)

        logger.info(" -> Success! All files generated.")

    except (RuntimeError, ValueError) as e:
        logger.error(" -> FAIL: g4DD_metric.py test failed with error: %s", e)
        sys.exit(1)
