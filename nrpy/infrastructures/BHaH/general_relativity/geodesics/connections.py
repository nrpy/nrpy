"""
Register C function for computing analytic Christoffel symbols.

This module registers the C function, which calculates the 40 unique components
of the Christoffel symbols for a specific spacetime. It generates a preamble
to unpack coordinate variables directly from a thread-local state vector to
minimize VRAM accesses and prevent register spilling on hardware architectures.
Author: Dalton J. Moone.
"""

# Python: Import standard modules
from typing import List, Set

# Python: Import third-party modules
import sympy as sp

# Python: Import NRPy core modules
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
from nrpy.equations.general_relativity.geodesics.analytic_spacetimes import (
    Analytic_Spacetimes,
)

def connections(
    Gamma4UDD_exprs: List[List[List[sp.Expr]]], spacetime_name: str, PARTICLE: str
) -> None:
    """
    Generate and register the C function to compute Christoffel symbols.

    The C function computes the unique components of the Christoffel symbols
    (Gamma^alpha_{mu nu}$) and writes them to a thread-local 1D array.

    :param Gamma4UDD_exprs: A 4x4x4 list of SymPy expressions for connections.
    :param spacetime_name: Name of the spacetime.
    :param PARTICLE: The type of particle ("massive" or "photon").
    :raises ValueError: If PARTICLE is not "massive" or "photon".
    """
    # Python: Specific setup based on particle type
    if PARTICLE == "massive":
        array_size = 8
    elif PARTICLE == "photon":
        array_size = 9
    else:
        raise ValueError(f"Unsupported PARTICLE: {PARTICLE}")

    list_of_Gamma_syms = []
    list_of_Gamma_C_vars = []

    k = 0
    for alpha in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                list_of_Gamma_syms.append(Gamma4UDD_exprs[alpha][mu][nu])
                # Python: Map 3D tensor components to the 1D thread-local C array
                list_of_Gamma_C_vars.append(f"Gamma_local[{k}]")
                k += 1

    # Python: Extract coordinates dynamically to populate the preamble
    used_symbol_names: Set[str] = set()
    for expr in list_of_Gamma_syms:
        for sym in expr.free_symbols:
            used_symbol_names.add(str(sym))

    xx_symbols = Analytic_Spacetimes[spacetime_name].xx
    preamble_lines = [
        "// Unpack position coordinates $x^i$ from the thread-local state vector.",
        f"// Evaluated at compile time for state vector size: {array_size}"
    ]

    for i, symbol in enumerate(xx_symbols):
        sym_name = str(symbol)
        if sym_name in used_symbol_names:
            # Python: Unpack coordinates strictly into fast local registers
            preamble_lines.append(f"const double {sym_name} = f_local[{i}];")

    preamble = "\n  ".join(preamble_lines)

    kernel = ccg.c_codegen(
        list_of_Gamma_syms,
        list_of_Gamma_C_vars,
        enable_cse=True,
        verbose=False,
        include_braces=False,
    )

    # Python: Define C-Function metadata in strict chronological order
    includes = ["BHaH_defines.h"]
    desc = rf"""@brief Computes the 40 unique Christoffel symbols $Gamma^alpha_{{mu nu}}$ for the {spacetime_name} metric.
    @param commondata Struct containing global spacetime parameters.
    @param f_local Thread-local array containing the 1D flattened state vector.
    @param Gamma_local Thread-local array where the connection components are stored."""
    cfunc_type = "static BHAH_HD_INLINE void"
    name = f"connections_{spacetime_name}"
    params = (
        "const commondata_struct *restrict commondata, "
        "const double *restrict f_local, "
        "double *restrict Gamma_local"
    )
    include_CodeParameters_h = True
    body = rf"""
    // --- CONNECTION EVALUATION & THREAD-LOCAL UNPACKING ---
    // Algorithmic Step: Extract spatial coordinates $x^i$ and compute $Gamma^alpha_{{mu nu}}$.
    // Hardware Justification: Thread-local execution prevents latency-heavy reads from VRAM 
    // and guarantees mathematical intermediates stay in hardware registers.
    {preamble}

    {kernel}
    """

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
    )

if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")