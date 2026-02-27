"""
Register C function for computing metric components.

This module registers the 'g4DD_metric_{spacetime_name}' C function,
which calculates the 10 unique components of the metric tensor for a specific spacetime.
It utilizes a flattened SoA architecture for SIMT-compatible batch processing on GPUs.

It generates a preamble to unpack coordinate variables directly from the flattened
state vector using macro indexing. Spacetime parameters stored in `commondata` are
provided through the function arguments.

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

    The C function computes the upper-triangular components of the metric
    tensor (g_munu) and writes them to a flattened 1D batch array.

    :param g4DD_exprs: A 4x4 list of SymPy expressions representing the metric.
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

    # Step 2: Define C function metadata
    includes = ["BHaH_defines.h"]
    desc = (
        f"@brief Computes the 10 unique components of the {spacetime_name} metric "
        f"for a {PARTICLE} particle.\n"
    )
    name = f"g4DD_metric_{spacetime_name}"

    # Updated SoA compatible signature
    params = (
        "const commondata_struct *restrict commondata, "
        "const double *restrict f, "
        "double *restrict metric_g4DD, "
        "const int batch_size, "
        "const int batch_id"
    )

    # Step 3: Prepare symbolic expressions and 1D flattened C variable names
    list_of_g4DD_syms = []
    list_of_g4DD_C_vars = []

    k = 0
    for mu in range(4):
        for nu in range(mu, 4):
            list_of_g4DD_syms.append(g4DD_exprs[mu][nu])
            # Map 2D symmetric tensor components to 1D flat array using macro
            list_of_g4DD_C_vars.append(f"metric_g4DD[IDX_LOCAL({k}, batch_id, batch_size)]")
            k += 1

    # Step 4: Generate the Dynamic Preamble
    # 4.a: Analyze symbols used in the expressions
    used_symbol_names: Set[str] = set()
    for expr in list_of_g4DD_syms:
        for sym in expr.free_symbols:
            used_symbol_names.add(str(sym))

    # 4.b: Retrieve coordinate symbols from the Analytic Spacetime registry
    xx_symbols = Analytic_Spacetimes[spacetime_name].xx
    preamble_lines = [
        f"// Unpack position coordinates from flattened state vector (State vector size: {array_size})"
    ]

    for i, symbol in enumerate(xx_symbols):
        sym_name = str(symbol)
        if sym_name in used_symbol_names:
            # Unpack coordinates directly from flattened input array
            preamble_lines.append(f"const double {sym_name} = f[IDX_LOCAL({i}, batch_id, batch_size)];")

    preamble = "\n  ".join(preamble_lines)

    # Step 5: Generate C body
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

    # Combine Preamble + Kernel with GPU portability wrappers
    body = f"""
        {preamble}

        {kernel} """

    prefunc = """
    #ifdef USE_GPU
    #pragma omp declare target
    #endif
    """

    postfunc = """
    #ifdef USE_GPU
    #pragma omp end declare target
    #endif
    """

    # Step 6: Register the C function
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=True,
        postfunc=postfunc
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
