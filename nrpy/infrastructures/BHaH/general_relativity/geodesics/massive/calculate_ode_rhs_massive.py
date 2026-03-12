r"""
Define C function for computing massive particle geodesic ODE right-hand sides.

This module defines the 'calculate_ode_rhs_massive' C function, which evaluates
the right-hand sides (RHS) of the coupled system of 8 first-order ODEs governing
massive particle geodesics.

It utilizes a thread-local architecture for maximum CPU register efficiency,
discarding the previous SoA block indexing methodology.

The system is derived from:
    $dx^\mu/d\tau = u^\mu$
    $du^\mu/d\tau = -\Gamma^\mu_{\alpha\beta} u^\alpha u^\beta$
Author: Dalton J. Moone.
"""

import sys
from typing import List

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc


def calculate_ode_rhs_massive(
    geodesic_rhs_expressions: List[sp.Expr], coordinate_symbols: List[sp.Symbol]
) -> None:
    """
    Define the C engine for the massive particle geodesic ODE right-hand sides.

    This function processes the SymPy expressions into highly optimized, thread-local
    C code for execution on CPU architectures.

    :param geodesic_rhs_expressions: List of 8 symbolic expressions for the RHS.
    :param coordinate_symbols: List of 4 symbolic coordinates (e.g., [t, x, y, z]).
    """
    # Extract unique symbols used in the mathematical expressions
    used_symbol_names = {
        str(sym) for expr in geodesic_rhs_expressions for sym in expr.free_symbols
    }

    # Construct the thread-local unpacking preamble
    preamble_lines = [
        "  // --- THREAD-LOCAL STATE UNPACKING ---",
        "  // Algorithmic Step: Extract spatial coordinates $x^i$ and 4-velocity $u^\\mu$ from the local state vector.",
        "  // Hardware Justification: Caching these values into local scalars forces the compiler to allocate fast hardware registers.",
    ]

    for i, sym in enumerate(coordinate_symbols):
        if str(sym) in used_symbol_names:
            preamble_lines.append(
                f"  const double {str(sym)} = f_local[{i}]; // Unpack coordinate ${str(sym)}$."
            )

    for i in range(4):
        if f"uU{i}" in used_symbol_names:
            preamble_lines.append(
                f"  const double uU{i} = f_local[{i+4}]; // Unpack 4-velocity component $u^{i}$."
            )

    preamble_lines.append("\n  // --- METRIC AND CONNECTION UNPACKING ---")
    preamble_lines.append(
        "  // Algorithmic Step: Extract pre-computed metric $g_{\\mu\\nu}$ and Christoffel symbols $\\Gamma^\\alpha_{\\mu\\nu}$."
    )
    preamble_lines.append(
        "  // Hardware Justification: Reading from the thread-local arrays guarantees L1 cache hits."
    )

    curr_idx = 0
    for m in range(4):
        for n in range(m, 4):
            comp_name = f"metric_g4DD{m}{n}"
            if comp_name in used_symbol_names:
                preamble_lines.append(
                    f"  const double {comp_name} = metric_local[{curr_idx}]; // Unpack metric component $g_{{{m}{n}}}$."
                )
            curr_idx += 1

    curr_idx = 0
    for a in range(4):
        for m in range(4):
            for n in range(m, 4):
                comp_name = f"conn_Gamma4UDD{a}{m}{n}"
                if comp_name in used_symbol_names:
                    preamble_lines.append(
                        f"  const double {comp_name} = Gamma_local[{curr_idx}]; // Unpack Christoffel symbol $\\Gamma^{a}_{{{m}{n}}}$."
                    )
                curr_idx += 1

    # Map SymPy RHS results directly to the thread-local output array
    k_array_outputs = [f"k_local[{i}]" for i in range(8)]

    # Generate the highly optimized CSE core math block
    kernel = ccg.c_codegen(
        geodesic_rhs_expressions, k_array_outputs, enable_cse=True, include_braces=False
    )

    # Define C-function arguments and metadata in Master Order
    includes = ["BHaH_defines.h"]

    desc = r"""@brief Portable CPU derivatives for the massive geodesic ODE system.

    Calculates $dx^\mu/d\tau$ and $du^\mu/d\tau$ using the provided state vector and
    pre-computed Christoffel symbols $\Gamma^\mu_{\alpha\beta}$. Uses thread-local arrays to guarantee
    intermediate values remain in hardware registers.

    @param f_local Local state vector $[t, x, y, z, u^t, u^x, u^y, u^z]$.
    @param metric_local Thread-local flattened metric array $g_{\mu\nu}$.
    @param Gamma_local Thread-local flattened Christoffel symbol array $\Gamma^\mu_{\alpha\beta}$.
    @param k_local The computed derivatives $dx^\mu/d\tau$ and $du^\mu/d\tau$ stored in thread-local format."""

    cfunc_type = "BHAH_HD_INLINE void"

    name = "calculate_ode_rhs_massive"

    params = (
        "const double *restrict f_local, "
        "const double *restrict metric_local, "
        "const double *restrict Gamma_local, "
        "double *restrict k_local"
    )

    include_CodeParameters_h = False

    body = "\n".join(preamble_lines) + "\n\n"
    body += "  // --- GEODESIC EQUATION EVALUATION ---\n"
    body += "  // Algorithmic Step: Evaluate the SymPy-generated Common Subexpression Elimination (CSE) block for the ODE RHS.\n"
    body += "  // Hardware Justification: CSE minimizes total FLOPs and maximizes register reuse during derivative calculation.\n"
    body += f"  {kernel}"

    # Register the function, omitting empty kwargs to prevent bloat
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

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
