"""
Module for generating the C code that computes the temporal momentum component.

This constraint solver forces the four-momentum vector to remain null
by solving the quadratic Hamiltonian constraint for the negative root.
Author: Dalton J. Moone.
"""

import sympy as sp
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc

def p0_reverse(p0_expr: sp.Expr) -> None:
    """
    Generate and register the C function to compute the initial temporal momentum component.

    :param p0_expr: The SymPy expression representing the negative root of the Hamiltonian constraint.
    :raises ValueError: If the symbolic expression fails to parse during code generation.
    """
    # Python: Map the 2D symmetric metric components to a 1D thread-local array index.
    metric_map = {}
    k = 0
    for i in range(4):
        for j in range(i, 4):
            old_symbol = f"metric_g4DD{i}{j}"
            new_access = f"metric_local[{k}]"
            metric_map[old_symbol] = new_access
            k += 1

    # Python: Generate the raw C math string from the SymPy expression.
    body_math = ccg.c_codegen(
        [p0_expr], ["*p0_out"], enable_cse=True, verbose=False, include_braces=False
    )

    # Python: Replace symbolic metric variables with the thread-local array accesses.
    for old, new in metric_map.items():
        body_math = body_math.replace(old, new)

    # Python: Define the ordered variables for CFunction registration.
    includes = ["BHaH_defines.h"]

    desc = r"""@brief Computes the initial temporal component of the 4-momentum $p^0$.

    Detailed algorithm: Solves the quadratic Hamiltonian constraint $p_\mu p^\mu = 0$
    specifically for the negative root to enforce physical null photon trajectories.
    This operation executes strictly within thread-local registers to satisfy the sm_86
    architecture limit of 255 registers per thread, bypassing VRAM bottlenecks.

    @param metric_local Thread-local array containing the 10 independent metric components $g_{\mu\nu}$.
    @param f_local Thread-local array containing the 9-component state vector $f^\mu$.
    @param p0_out Pointer to the thread-local scalar storing the resulting temporal momentum $p^0$."""

    cfunc_type = "BHAH_HD_INLINE void"
    
    name = "p0_reverse"

    params = (
        "const double *restrict metric_local, "
        "const double *restrict f_local, "
        "double *restrict p0_out"
    )

    include_CodeParameters_h = False

    body = rf"""
    // --- SPATIAL MOMENTUM EXTRACTION ---
    /* * Algorithmic Step: Extract contravariant spatial momentum components to solve the Hamiltonian constraint.
     * Architectural Justification: Thread-local reads bypass global memory latency during the root-finding phase,
     * ensuring immediate access to the $f^\mu$ state vector components mapped directly to physical registers.
     */
    const double pU1 = f_local[5]; // Contravariant spatial momentum component $p^x$.
    const double pU2 = f_local[6]; // Contravariant spatial momentum component $p^y$.
    const double pU3 = f_local[7]; // Contravariant spatial momentum component $p^z$.

    // --- HAMILTONIAN CONSTRAINT ROOT FINDING ---
    /* * Algorithmic Step: Evaluate the generated algebraic solution for $p^0$.
     * Architectural Justification: The symbolic generation maps algebraic operations directly to hardware 
     * Fused Multiply-Add (FMA) instructions, maximizing floating-point execution throughput.
     */
    {body_math}
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