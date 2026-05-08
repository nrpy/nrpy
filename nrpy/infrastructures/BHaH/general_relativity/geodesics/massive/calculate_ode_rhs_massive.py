# nrpy/infrastructures/BHaH/general_relativity/geodesics/massive/calculate_ode_rhs_massive.py
r"""
Provides the infrastructure for computing massive particle geodesic right-hand sides.

This module defines the Python function that constructs and registers the C kernel
evaluating the right-hand sides (RHS) of the eight first-order ordinary differential
equations governing massive particle geodesics. The mathematical system derives from
the geodesic equations $dx^\mu/d\tau = u^\mu$ and
$du^\mu/d\tau = -\Gamma^\mu_{\alpha\beta} u^\alpha u^\beta$. The code generation
process translates SymPy expressions into C code. Caching values into local scalars
and reading from local arrays organizes memory access. Common Subexpression Elimination
(CSE) reduces floating-point operations during the derivative calculation.

Author: Dalton J. Moone
"""

from typing import List

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc


def calculate_ode_rhs_massive(
    geodesic_rhs_expressions: List[sp.Expr], coordinate_symbols: List[sp.Symbol]
) -> None:
    r"""
    Define the C engine for the massive particle geodesic ODE right-hand sides.

    This function processes the SymPy expressions into thread-local C code. Caching
    values into local scalars and reading from thread-local arrays organizes data
    transfer and memory allocation. Common Subexpression Elimination (CSE) reduces
    floating-point operations and increases register reuse during the derivative
    calculation.

    :param geodesic_rhs_expressions: List of 8 symbolic expressions for the RHS.
    :param coordinate_symbols: List of 4 symbolic coordinates (e.g., [t, x, y, z]).
    """
    # Extract unique symbols used in the mathematical expressions
    used_symbol_names = {
        str(sym) for expr in geodesic_rhs_expressions for sym in expr.free_symbols
    }

    # Construct the thread-local unpacking preamble
    preamble_lines = [
        "  //==========================================",
        "  // THREAD-LOCAL STATE UNPACKING",
        "  //==========================================",
        "  // Extract spatial coordinates $x^i$ and 4-velocity $u^\\mu$ from the local state vector.",
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

    preamble_lines.append(
        "\n  //==========================================\n  // METRIC AND CONNECTION UNPACKING\n  //=========================================="
    )
    preamble_lines.append(
        "  // Extract pre-computed metric $g_{\\mu\\nu}$ and Christoffel symbols $\\Gamma^\\alpha_{\\mu\\nu}$."
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

    # Define C-function arguments and metadata
    includes = ["BHaH_defines.h"]

    desc = r""" Portable derivatives for the massive geodesic ODE system.

    Calculates $dx^\mu/d\tau$ and $du^\mu/d\tau$ using the provided state vector and
    pre-computed Christoffel symbols $\Gamma^\mu_{\alpha\beta}$. Thread-local arrays
    organize intermediate values during the calculation.

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
    body += "  //==========================================\n"
    body += "  // GEODESIC EQUATION EVALUATION\n"
    body += "  //==========================================\n"
    body += "  // Evaluate the SymPy-generated Common Subexpression Elimination (CSE) block for the ODE RHS.\n"
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
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
