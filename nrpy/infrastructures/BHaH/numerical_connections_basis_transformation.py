"""
Generate a C function that transforms numerical grid-basis Christoffels to Cartesian basis.

This module registers a standalone BHaH geodesics infrastructure C function named
``numerical_connections_basis_transformation``. The generated C function assumes
that an upstream interpolation or reconstruction stage has already evaluated the
covariant 4-metric, its first derivatives with respect to the grid-basis spacetime
coordinates, the spatial Jacobian from grid coordinates to Cartesian coordinates,
and the spatial Jacobian derivatives at the photon position. From those local
inputs it computes the grid-basis Levi-Civita connection, applies the spatial
basis-transformation correction appropriate for Christoffel symbols, and packs the
40 upper-triangular output components in the same order expected by the geodesics
RHS infrastructure.

Author: OpenAI Codex
        support **at** openai **dot** com
"""

# Import standard modules
from typing import List

# Import third-party modules
import sympy as sp

# Import NRPy core modules
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.indexedexp as ixp


def _metric_component_pairs() -> List[str]:
    """
    Return the canonical upper-triangular 4-metric component suffixes.

    :return: The packed component suffix ordering.
    """
    return ["00", "01", "02", "03", "11", "12", "13", "22", "23", "33"]


def numerical_connections_basis_transformation() -> None:
    r"""
    Register the numerical basis-transformation Christoffel helper.

    The generated C function consumes local geometric data at a photon position:
    the covariant metric \(g_{\mu\nu}\), the derivatives \(\partial_A g_{\mu\nu}\),
    the spatial Jacobian \(J^i{}_a\), and the Jacobian derivatives
    \(\partial_b J^i{}_c\). It computes the grid-basis Levi-Civita connection,
    transforms the spatial basis to Cartesian coordinates, and writes the 40
    upper-triangular connection components \(\Gamma^\alpha_{\mu\nu}\) to the
    output buffer in the same order already used by the geodesics RHS kernel.
    """
    metric_suffixes = _metric_component_pairs()

    # Step 1: Declare symbolic local-input arrays using the packed storage contract.
    g4DD_packed = [sp.symbols(f"g4DD{suffix}", real=True) for suffix in metric_suffixes]
    dg4DD_packed = [
        sp.symbols(f"dg4DD{dirn}{suffix}", real=True)
        for dirn in range(4)
        for suffix in metric_suffixes
    ]
    J_packed = [sp.symbols(f"J{i}{a}", real=True) for i in range(3) for a in range(3)]
    dJ_packed = [
        sp.symbols(f"dJ{b}{i}{c}", real=True)
        for b in range(3)
        for i in range(3)
        for c in range(3)
    ]

    # Step 2: Rebuild symmetric tensor storage from packed inputs.
    g4DD = ixp.zerorank2(dimension=4)
    packed_idx = 0
    for mu in range(4):
        for nu in range(mu, 4):
            g4DD[mu][nu] = g4DD_packed[packed_idx]
            g4DD[nu][mu] = g4DD_packed[packed_idx]
            packed_idx += 1

    g4DD_dD = ixp.zerorank3(dimension=4)
    for dirn in range(4):
        packed_idx = dirn * len(metric_suffixes)
        for mu in range(4):
            for nu in range(mu, 4):
                g4DD_dD[mu][nu][dirn] = dg4DD_packed[packed_idx]
                g4DD_dD[nu][mu][dirn] = dg4DD_packed[packed_idx]
                packed_idx += 1

    J = [[sp.sympify(0) for _ in range(3)] for __ in range(3)]
    packed_idx = 0
    for i in range(3):
        for a in range(3):
            J[i][a] = J_packed[packed_idx]
            packed_idx += 1

    dJ = [[[sp.sympify(0) for _ in range(3)] for __ in range(3)] for ___ in range(3)]
    packed_idx = 0
    for b in range(3):
        for i in range(3):
            for c in range(3):
                dJ[b][i][c] = dJ_packed[packed_idx]
                packed_idx += 1

    # Step 3: Compute the grid-basis 4D Christoffel symbols.
    g4UU, _ = ixp.symm_matrix_inverter4x4(g4DD)
    Gamma_grid = ixp.zerorank3(dimension=4)
    for alpha in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                term = sp.sympify(0)
                for beta in range(4):
                    term += (
                        sp.Rational(1, 2)
                        * g4UU[alpha][beta]
                        * (
                            g4DD_dD[nu][beta][mu]
                            + g4DD_dD[mu][beta][nu]
                            - g4DD_dD[mu][nu][beta]
                        )
                    )
                Gamma_grid[alpha][mu][nu] = term
                Gamma_grid[alpha][nu][mu] = term

    # Step 4: Compute the inverse spatial Jacobian K^a_i.
    K_matrix = sp.Matrix(J).inv()
    K = [[sp.sympify(0) for _ in range(3)] for __ in range(3)]
    for a in range(3):
        for i in range(3):
            K[a][i] = K_matrix[a, i]

    # Step 5: Transform the connection from (t, u^a) to (t, x^i).
    Gamma_cart = ixp.zerorank3(dimension=4)

    # Step 5.a: Upper time index, no non-tensor correction term, but spatial lower
    #           indices still transform with the inverse Jacobian factors.
    Gamma_cart[0][0][0] = Gamma_grid[0][0][0]
    for j in range(3):
        term = sp.sympify(0)
        for c in range(3):
            term += K[c][j] * Gamma_grid[0][0][c + 1]
        Gamma_cart[0][0][j + 1] = term
        Gamma_cart[0][j + 1][0] = term

    for j in range(3):
        for k in range(j, 3):
            term = sp.sympify(0)
            for b in range(3):
                for c in range(3):
                    term += K[b][j] * K[c][k] * Gamma_grid[0][b + 1][c + 1]
            Gamma_cart[0][j + 1][k + 1] = term
            Gamma_cart[0][k + 1][j + 1] = term

    # Step 5.b: Upper spatial index, mixed lower indices with time unchanged.
    for i in range(3):
        term_00 = sp.sympify(0)
        for a in range(3):
            term_00 += J[i][a] * Gamma_grid[a + 1][0][0]
        Gamma_cart[i + 1][0][0] = term_00

        for j in range(3):
            term_0j = sp.sympify(0)
            for a in range(3):
                for c in range(3):
                    term_0j += J[i][a] * K[c][j] * Gamma_grid[a + 1][0][c + 1]
            Gamma_cart[i + 1][0][j + 1] = term_0j
            Gamma_cart[i + 1][j + 1][0] = term_0j

    # Step 5.c: Upper spatial index with two spatial lower indices, including the
    #           non-tensor correction term.
    for i in range(3):
        for j in range(3):
            for k in range(j, 3):
                tensor_term = sp.sympify(0)
                correction_term = sp.sympify(0)
                for a in range(3):
                    for b in range(3):
                        for c in range(3):
                            tensor_term += (
                                J[i][a]
                                * K[b][j]
                                * K[c][k]
                                * Gamma_grid[a + 1][b + 1][c + 1]
                            )
                for b in range(3):
                    for c in range(3):
                        correction_term += K[b][j] * K[c][k] * dJ[b][i][c]
                Gamma_cart[i + 1][j + 1][k + 1] = tensor_term - correction_term
                Gamma_cart[i + 1][k + 1][j + 1] = tensor_term - correction_term

    # Step 6: Pack only the 40 upper-triangular lower-index outputs.
    output_exprs: List[sp.Expr] = []
    output_vars: List[str] = []
    packed_idx = 0
    for alpha in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                output_exprs.append(Gamma_cart[alpha][mu][nu])
                output_vars.append(f"Gamma_local[{packed_idx}]")
                packed_idx += 1

    kernel = ccg.c_codegen(
        output_exprs,
        output_vars,
        enable_cse=True,
        verbose=False,
        include_braces=False,
    )

    # Step 7: Build the packed-array unpacking preamble for the generated C body.
    preamble_lines = [
        "(void)commondata;",
        "",
        "//==========================================",
        "// LOCAL INPUT UNPACKING",
        "//==========================================",
        "// Unpack the metric, derivatives, Jacobian, and Jacobian derivatives from",
        "// their packed local-array storage into scalar registers.",
    ]
    for packed_idx, suffix in enumerate(metric_suffixes):
        preamble_lines.append(f"const double g4DD{suffix} = g4DD_local[{packed_idx}];")

    for dirn in range(4):
        for packed_idx, suffix in enumerate(metric_suffixes):
            preamble_lines.append(
                f"const double dg4DD{dirn}{suffix} = "
                f"dg4DD_local[{dirn * len(metric_suffixes) + packed_idx}];"
            )

    packed_idx = 0
    for i in range(3):
        for a in range(3):
            preamble_lines.append(f"const double J{i}{a} = J_local[{packed_idx}];")
            packed_idx += 1

    packed_idx = 0
    for b in range(3):
        for i in range(3):
            for c in range(3):
                preamble_lines.append(
                    f"const double dJ{b}{i}{c} = dJ_local[{packed_idx}];"
                )
                packed_idx += 1

    preamble = "\n    ".join(preamble_lines)

    # Step 8: Register the generated C function.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r""" Compute the 40 upper-triangular Christoffel symbols from local numerical geometry data.

    @param commondata Struct containing global infrastructure parameters.
    @param g4DD_local Packed local 4-metric values $g_{\mu\nu}(u_p)$ in the order
      {00,01,02,03,11,12,13,22,23,33}.
    @param dg4DD_local Packed local metric derivatives $\partial_A g_{\mu\nu}(u_p)$.
      The outer index is the derivative direction $A \in \{0,1,2,3\}$ and the inner
      metric packing order matches g4DD_local.
    @param J_local Packed spatial Jacobian values $J^i{}_a(u_p)$ with index order
      $(i,a)$ and $i$ outermost.
    @param dJ_local Packed spatial Jacobian derivatives $\partial_b J^i{}_c(u_p)$
      with index order $(b,i,c)$ and $b$ outermost.
    @param Gamma_local Output array storing the 40 packed upper-triangular
      Christoffel components $\Gamma^\alpha_{\mu\nu}$ in the existing geodesics order.
    """
    cfunc_type = "void"
    name = "numerical_connections_basis_transformation"
    params = (
        "const commondata_struct *restrict commondata, "
        "const double *restrict g4DD_local, "
        "const double *restrict dg4DD_local, "
        "const double *restrict J_local, "
        "const double *restrict dJ_local, "
        "double *restrict Gamma_local"
    )
    body = rf"""
    //==========================================
    // LOCAL DATA INGEST
    //==========================================
    {preamble}

    //==========================================
    // CHRISTOFFEL EVALUATION AND BASIS TRANSFORMATION
    //==========================================
    // Compute the grid-basis 4D Levi-Civita connection, transform the spatial
    // basis to Cartesian coordinates, and pack the 40 upper-triangular outputs.
    {kernel}
    """

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
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
