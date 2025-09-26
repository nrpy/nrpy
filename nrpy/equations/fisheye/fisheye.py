r"""
Core Jacobian (basis) transformation functions.
For Fisheye <-> Cartesian basis transformations (no spherical, no rfm).

We define Jacobians for the purely radial Fisheye map at a point with
primed Cartesian coordinates x'^a:

Given r = sqrt(δ_ij x^i x^j) and r' = r'(r), the map is
    x'^i = x^i * (r'/r).

Define the Jacobian and its inverse by
    J^i{}_a        = ∂x^i / ∂x'^a,
    (J^{-1})^a{}_i = ∂x'^a / ∂x^i.

For the fisheye map:
    J^i{}_a        = δ^i{}_a (r/r') + (x'^i x'^a / r'^2) * (dr/dr' - r/r'),
    (J^{-1})^a{}_i = δ^a{}_i (r'/r) + (x'^a x'^i / r'^2) * (dr'/dr - r'/r).

These formulas follow the “Fisheye Coordinates” appendix of:
J. A. Faber, T. W. Baumgarte, Z. B. Etienne, S. L. Shapiro, and K. Taniguchi,
"Relativistic hydrodynamics in the presence of puncture black holes," Phys. Rev. D (2017),
Appendix: Fisheye Coordinates.

Author: Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from typing import Sequence, Tuple
import sympy as sp
import nrpy.indexedexp as ixp


# ================================================================
# Fisheye profile r'(r) and its radial derivative dr'/dr (helpers)
# ================================================================


def fisheye_rprime_of_r(
    r: sp.Expr,
    a: Sequence[sp.Expr],
    R: Sequence[sp.Expr],
    s: Sequence[sp.Expr],
) -> sp.Expr:
    r"""
    Compute r'(r) for the radial fisheye profile.

    :param r: Physical radius \(r\).
    :param a: Coefficients \([a_0, a_1, \ldots, a_n]\), where \(a_n\) is the outermost scale.
    :param R: Transition centers \([R_1, \ldots, R_n]\).
    :param s: Transition widths \([s_1, \ldots, s_n]\).
    :return: The fisheye radius \(r'(r)\) as a SymPy expression.

    Doctests:
    >>> r = sp.Integer(3)
    >>> # Pure global scale: n=0, a=[2] -> r' = 2*r
    >>> fisheye_rprime_of_r(r, a=[sp.Integer(2)], R=[], s=[])
    6
    """
    n = len(R)
    assert len(s) == n and len(a) == n + 1
    rp = a[-1] * r
    for i in range(n):
        num = (a[i] - a[i + 1]) * s[i]
        denom = 2 * sp.tanh(R[i] / s[i])
        term = sp.log(sp.cosh((r + R[i]) / s[i]) / sp.cosh((r - R[i]) / s[i]))
        rp += num / denom * term
    return sp.simplify(rp)


def fisheye_drprime_dr(
    r: sp.Expr,
    a: Sequence[sp.Expr],
    R: Sequence[sp.Expr],
    s: Sequence[sp.Expr],
) -> sp.Expr:
    r"""
    Compute \(dr'/dr\) for the radial fisheye profile.

    :param r: Physical radius \(r\) (symbolic allowed).
    :param a: Coefficients \([a_0, a_1, \ldots, a_n]\), where \(a_n\) is the outermost scale.
    :param R: Transition centers \([R_1, \ldots, R_n]\).
    :param s: Transition widths \([s_1, \ldots, s_n]\).
    :return: The radial derivative \(dr'/dr\) as a SymPy expression.

    Doctests:
    >>> r = sp.symbols("r", positive=True)
    >>> # Pure global scale: n=0, a=[2] -> dr'/dr = 2
    >>> fisheye_drprime_dr(r, a=[sp.Integer(2)], R=[], s=[])
    2
    """
    n = len(R)
    assert len(s) == n and len(a) == n + 1
    drp = a[-1]
    for i in range(n):
        coeff = (a[i] - a[i + 1]) / sp.tanh(R[i] / s[i])
        bracket = (sp.tanh((r + R[i]) / s[i]) - sp.tanh((r - R[i]) / s[i])) / 2
        drp += coeff * bracket
    return sp.simplify(drp)


# ================================================================
# Core: build J and J^{-1} from {x'^i, r', r, dr'/dr}
# ================================================================


def build_fisheye_jacobians_from_xprime(
    xprimeU: Sequence[sp.Expr],
    rprime: sp.Expr,
    r: sp.Expr,
    drprime_dr: sp.Expr,
) -> Tuple[Sequence[Sequence[sp.Expr]], Sequence[Sequence[sp.Expr]]]:
    r"""
    Construct spatial Jacobians \(J\) and \(J^{-1}\) for the fisheye map.

    Construct \(J^i{}_a = \partial x^i / \partial x'^a\) and
    \((J^{-1})^a{}_i = \partial x'^a / \partial x^i\) using the exact radial formulas.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\) at the point.
    :param rprime: Fisheye radius \(r' = \sqrt{x'^i x'^i}\).
    :param r: Physical radius \(r = \sqrt{x^i x^i}\) (consistent with \(r'\)).
    :param drprime_dr: Radial derivative \(dr'/dr\) evaluated at \(r\).
    :return: Tuple \((J, Jinv)\) where \(J\) is \(3\times 3\) with entries \(J[i][a]\),
             and \(Jinv\) is \(3\times 3\) with entries \(Jinv[a][i]\).

    Doctests:
    >>> # Global scaling test: x'^i = 2 x^i => r' = 2r, dr'/dr = 2.
    >>> x0p, x1p, x2p = 2, 0, 0
    >>> xprimeU = [sp.Integer(x0p), sp.Integer(x1p), sp.Integer(x2p)]
    >>> rprime = sp.Integer(2)     # since point corresponds to x=(1,0,0)
    >>> r = sp.Integer(1)
    >>> drpdr = sp.Integer(2)
    >>> J, Jinv = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drpdr)
    >>> [sp.simplify(J[i][i]) for i in range(3)]
    [1/2, 1/2, 1/2]
    >>> [sp.simplify(Jinv[i][i]) for i in range(3)]
    [2, 2, 2]
    >>> # Check JJ^{-1} = I at this sample:
    >>> M = [[sp.simplify(sum(J[i][k]*Jinv[k][j] for k in range(3))) for j in range(3)] for i in range(3)]
    >>> M
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    """
    dim = 3
    J = ixp.zerorank2(dimension=dim)  # rows i, cols a
    Jinv = ixp.zerorank2(dimension=dim)  # rows a, cols i

    dr_drprime = sp.simplify(1 / drprime_dr)

    s1 = sp.simplify(r / rprime)
    s2 = sp.simplify(dr_drprime - r / rprime)

    s1i = sp.simplify(rprime / r)
    s2i = sp.simplify(drprime_dr - rprime / r)

    rprime2 = sp.simplify(rprime**2)
    fact = sp.simplify(1 / rprime2)

    # J^i{}_a
    for i in range(dim):
        for aidx in range(dim):
            delta = sp.Integer(1) if i == aidx else sp.Integer(0)
            J[i][aidx] = delta * s1 + fact * xprimeU[i] * xprimeU[aidx] * s2

    # (J^{-1})^a{}_i
    for aidx in range(dim):
        for i in range(dim):
            delta = sp.Integer(1) if aidx == i else sp.Integer(0)
            Jinv[aidx][i] = delta * s1i + fact * xprimeU[aidx] * xprimeU[i] * s2i

    return J, Jinv


# ================================================================
# Transforms: fisheye <-> Cartesian
# ================================================================


def basis_transform_vectorU_from_fisheye_to_Cartesian(
    xprimeU: Sequence[sp.Expr],
    rprime: sp.Expr,
    r: sp.Expr,
    drprime_dr: sp.Expr,
    src_vectorU_prime: Sequence[sp.Expr],
) -> Sequence[sp.Expr]:
    r"""
    Transform a contravariant vector from fisheye basis to Cartesian.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\).
    :param rprime: Fisheye radius \(r'\).
    :param r: Physical radius \(r\).
    :param drprime_dr: Radial derivative \(dr'/dr\) at \(r\).
    :param src_vectorU_prime: Source contravariant vector \(v'^a\) in fisheye basis.
    :return: Transformed contravariant vector \(v^i\) in Cartesian basis.

    Doctests:
    >>> # With global r'=2r (J = 1/2 I), a contravariant vector halves.
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> rprime, r, drpdr = sp.Integer(2), sp.Integer(1), sp.Integer(2)
    >>> vprime = [sp.Integer(3), sp.Integer(4), sp.Integer(5)]
    >>> v = basis_transform_vectorU_from_fisheye_to_Cartesian(xprimeU, rprime, r, drpdr, vprime)
    >>> [sp.simplify(comp) for comp in v]
    [3/2, 2, 5/2]
    """
    J, _ = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drprime_dr)
    dim = 3
    dst = ixp.zerorank1(dimension=dim)
    for i in range(dim):
        for aidx in range(dim):
            dst[i] += J[i][aidx] * src_vectorU_prime[aidx]
    return dst


def basis_transform_vectorD_from_fisheye_to_Cartesian(
    xprimeU: Sequence[sp.Expr],
    rprime: sp.Expr,
    r: sp.Expr,
    drprime_dr: sp.Expr,
    src_vectorD_prime: Sequence[sp.Expr],
) -> Sequence[sp.Expr]:
    r"""
    Transform a covariant vector from fisheye basis to Cartesian.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\).
    :param rprime: Fisheye radius \(r'\).
    :param r: Physical radius \(r\).
    :param drprime_dr: Radial derivative \(dr'/dr\) at \(r\).
    :param src_vectorD_prime: Source covariant vector \(\omega'_a\) in fisheye basis.
    :return: Transformed covariant vector \(\omega_i\) in Cartesian basis.

    Doctests:
    >>> # With global r'=2r (Jinv = 2 I), a covector doubles.
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> rprime, r, drpdr = sp.Integer(2), sp.Integer(1), sp.Integer(2)
    >>> wprime = [sp.Integer(1), sp.Integer(-2), sp.Integer(3)]
    >>> w = basis_transform_vectorD_from_fisheye_to_Cartesian(xprimeU, rprime, r, drpdr, wprime)
    >>> [sp.simplify(comp) for comp in w]
    [2, -4, 6]
    """
    _, Jinv = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drprime_dr)
    dim = 3
    dst = ixp.zerorank1(dimension=dim)
    for i in range(dim):
        for aidx in range(dim):
            dst[i] += Jinv[aidx][i] * src_vectorD_prime[aidx]
    return dst


def basis_transform_tensorDD_from_fisheye_to_Cartesian(
    xprimeU: Sequence[sp.Expr],
    rprime: sp.Expr,
    r: sp.Expr,
    drprime_dr: sp.Expr,
    src_tensorDD_prime: Sequence[Sequence[sp.Expr]],
) -> Sequence[Sequence[sp.Expr]]:
    r"""
    Transform a rank-2 covariant tensor from fisheye basis to Cartesian.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\).
    :param rprime: Fisheye radius \(r'\).
    :param r: Physical radius \(r\).
    :param drprime_dr: Radial derivative \(dr'/dr\) at \(r\).
    :param src_tensorDD_prime: Source covariant tensor \(T'_{ab}\) in fisheye basis.
    :return: Transformed covariant tensor \(T_{ij}\) in Cartesian basis.

    Doctests:
    >>> # With global r'=2r, Jinv = 2I, so T_{ij} = 4 T'_{ij} for covariant rank-2.
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> rprime, r, drpdr = sp.Integer(2), sp.Integer(1), sp.Integer(2)
    >>> Tprime = [[sp.Integer(1), sp.Integer(2), sp.Integer(0)],
    ...           [sp.Integer(2), sp.Integer(5), sp.Integer(-1)],
    ...           [sp.Integer(0), sp.Integer(-1), sp.Integer(4)]]
    >>> T = basis_transform_tensorDD_from_fisheye_to_Cartesian(xprimeU, rprime, r, drpdr, Tprime)
    >>> [[sp.simplify(T[i][j]) for j in range(3)] for i in range(3)]
    [[4, 8, 0], [8, 20, -4], [0, -4, 16]]
    """
    _, Jinv = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drprime_dr)
    dim = 3
    dst = ixp.zerorank2(dimension=dim)
    for i in range(dim):
        for j in range(dim):
            for aidx in range(dim):
                for bidx in range(dim):
                    dst[i][j] += (
                        Jinv[aidx][i] * Jinv[bidx][j] * src_tensorDD_prime[aidx][bidx]
                    )
    return dst


def basis_transform_vectorU_from_Cartesian_to_fisheye(
    xprimeU: Sequence[sp.Expr],
    rprime: sp.Expr,
    r: sp.Expr,
    drprime_dr: sp.Expr,
    src_vectorU_cart: Sequence[sp.Expr],
) -> Sequence[sp.Expr]:
    r"""
    Transform a contravariant vector from Cartesian to fisheye basis.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\).
    :param rprime: Fisheye radius \(r'\).
    :param r: Physical radius \(r\).
    :param drprime_dr: Radial derivative \(dr'/dr\) at \(r\).
    :param src_vectorU_cart: Source contravariant vector \(v^i\) in Cartesian basis.
    :return: Transformed contravariant vector \(v'^a\) in fisheye basis.

    Doctests:
    >>> # Inverse of previous: vectorU scales up by factor 2.
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> rprime, r, drpdr = sp.Integer(2), sp.Integer(1), sp.Integer(2)
    >>> vcart = [sp.Rational(3,2), sp.Integer(2), sp.Rational(5,2)]
    >>> vprime = basis_transform_vectorU_from_Cartesian_to_fisheye(xprimeU, rprime, r, drpdr, vcart)
    >>> [sp.simplify(comp) for comp in vprime]
    [3, 4, 5]
    """
    _, Jinv = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drprime_dr)
    dim = 3
    dst = ixp.zerorank1(dimension=dim)
    for aidx in range(dim):
        for i in range(dim):
            dst[aidx] += Jinv[aidx][i] * src_vectorU_cart[i]
    return dst


def basis_transform_vectorD_from_Cartesian_to_fisheye(
    xprimeU: Sequence[sp.Expr],
    rprime: sp.Expr,
    r: sp.Expr,
    drprime_dr: sp.Expr,
    src_vectorD_cart: Sequence[sp.Expr],
) -> Sequence[sp.Expr]:
    r"""
    Transform a covariant vector from Cartesian to fisheye basis.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\).
    :param rprime: Fisheye radius \(r'\).
    :param r: Physical radius \(r\).
    :param drprime_dr: Radial derivative \(dr'/dr\) at \(r\).
    :param src_vectorD_cart: Source covariant vector \(\omega_i\) in Cartesian basis.
    :return: Transformed covariant vector \(\omega'_a\) in fisheye basis.

    Doctests:
    >>> # Inverse of previous covector: halves.
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> rprime, r, drpdr = sp.Integer(2), sp.Integer(1), sp.Integer(2)
    >>> wcart = [sp.Integer(2), sp.Integer(-4), sp.Integer(6)]
    >>> wprime = basis_transform_vectorD_from_Cartesian_to_fisheye(xprimeU, rprime, r, drpdr, wcart)
    >>> [sp.simplify(comp) for comp in wprime]
    [1, -2, 3]
    """
    J, _ = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drprime_dr)
    dim = 3
    dst = ixp.zerorank1(dimension=dim)
    for aidx in range(dim):
        for i in range(dim):
            dst[aidx] += J[i][aidx] * src_vectorD_cart[i]
    return dst


def basis_transform_tensorDD_from_Cartesian_to_fisheye(
    xprimeU: Sequence[sp.Expr],
    rprime: sp.Expr,
    r: sp.Expr,
    drprime_dr: sp.Expr,
    src_tensorDD_cart: Sequence[Sequence[sp.Expr]],
) -> Sequence[Sequence[sp.Expr]]:
    r"""
    Transform a rank-2 covariant tensor from Cartesian to fisheye basis.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\).
    :param rprime: Fisheye radius \(r'\).
    :param r: Physical radius \(r\).
    :param drprime_dr: Radial derivative \(dr'/dr\) at \(r\).
    :param src_tensorDD_cart: Source covariant tensor \(T_{ij}\) in Cartesian basis.
    :return: Transformed covariant tensor \(T'_{ab}\) in fisheye basis.

    Doctests:
    >>> # Inverse of previous tensor test: divides by 4.
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> rprime, r, drpdr = sp.Integer(2), sp.Integer(1), sp.Integer(2)
    >>> Tcart = [[sp.Integer(4), sp.Integer(8), sp.Integer(0)],
    ...          [sp.Integer(8), sp.Integer(20), sp.Integer(-4)],
    ...          [sp.Integer(0), sp.Integer(-4), sp.Integer(16)]]
    >>> Tprime = basis_transform_tensorDD_from_Cartesian_to_fisheye(xprimeU, rprime, r, drpdr, Tcart)
    >>> [[sp.simplify(Tprime[i][j]) for j in range(3)] for i in range(3)]
    [[1, 2, 0], [2, 5, -1], [0, -1, 4]]
    """
    J, _ = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drprime_dr)
    dim = 3
    dst = ixp.zerorank2(dimension=dim)
    for aidx in range(dim):
        for bidx in range(dim):
            for i in range(dim):
                for j in range(dim):
                    dst[aidx][bidx] += J[i][aidx] * J[j][bidx] * src_tensorDD_cart[i][j]
    return dst


# ================================================================
# 4-tensor (UU), time-independent spatial map
# ================================================================


def basis_transform_4tensorUU_from_time_indep_fisheye_to_Cartesian(
    xprimeU: Sequence[sp.Expr],
    rprime: sp.Expr,
    r: sp.Expr,
    drprime_dr: sp.Expr,
    T4UU_prime: Sequence[Sequence[sp.Expr]],
) -> Sequence[Sequence[sp.Expr]]:
    r"""
    Transform a contravariant 4-tensor (time-independent spatial map) from fisheye to Cartesian.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\).
    :param rprime: Fisheye radius \(r'\).
    :param r: Physical radius \(r\).
    :param drprime_dr: Radial derivative \(dr'/dr\) at \(r\).
    :param T4UU_prime: Source 4-tensor \(T'^{\mu\nu}\) (contravariant) in fisheye basis.
    :return: Transformed 4-tensor \(T^{\mu\nu}\) (contravariant) in Cartesian basis.

    Doctests:
    >>> # For global r'=2r, each spatial contravariant index contributes a factor 1/2.
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> rprime, r, drpdr = sp.Integer(2), sp.Integer(1), sp.Integer(2)
    >>> # Spatial block diag [2,4,6] -> after map becomes [1,2,3] per single index, i.e., 1/2 each.
    >>> T4p = [[sp.Integer(7),0,0,0],[0,sp.Integer(2),0,0],[0,0,sp.Integer(4),0],[0,0,0,sp.Integer(6)]]
    >>> T4 = basis_transform_4tensorUU_from_time_indep_fisheye_to_Cartesian(xprimeU, rprime, r, drpdr, T4p)
    >>> [sp.simplify(T4[i][i]) for i in range(4)]
    [7, 1/2, 1, 3/2]
    """
    dim = 3
    J, _ = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drprime_dr)

    Jac4 = ixp.zerorank2(dimension=4)
    Jac4[0][0] = sp.Integer(1)
    for i in range(dim):
        for j in range(dim):
            Jac4[i + 1][j + 1] = J[i][j]

    dst = ixp.zerorank2(dimension=4)
    for mu in range(4):
        for nu in range(4):
            for a in range(4):
                for b in range(4):
                    dst[mu][nu] += Jac4[mu][a] * Jac4[nu][b] * T4UU_prime[a][b]
    return dst


def basis_transform_4tensorUU_from_Cartesian_to_time_indep_fisheye(
    xprimeU: Sequence[sp.Expr],
    rprime: sp.Expr,
    r: sp.Expr,
    drprime_dr: sp.Expr,
    T4UU_cart: Sequence[Sequence[sp.Expr]],
) -> Sequence[Sequence[sp.Expr]]:
    r"""
    Transform a contravariant 4-tensor (time-independent spatial map) from Cartesian to fisheye.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\).
    :param rprime: Fisheye radius \(r'\).
    :param r: Physical radius \(r\).
    :param drprime_dr: Radial derivative \(dr'/dr\) at \(r\).
    :param T4UU_cart: Source 4-tensor \(T^{\mu\nu}\) (contravariant) in Cartesian basis.
    :return: Transformed 4-tensor \(T'^{\mu\nu}\) (contravariant) in fisheye basis.

    Doctests:
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> rprime, r, drpdr = sp.Integer(2), sp.Integer(1), sp.Integer(2)
    >>> # Spatial block diag [1,2,3] -> after inverse map becomes [4,8,12] (factor 2 per index).
    >>> T4c = [[sp.Integer(7),0,0,0],[0,sp.Integer(1),0,0],[0,0,sp.Integer(2),0],[0,0,0,sp.Integer(3)]]
    >>> T4p = basis_transform_4tensorUU_from_Cartesian_to_time_indep_fisheye(xprimeU, rprime, r, drpdr, T4c)
    >>> [sp.simplify(T4p[i][i]) for i in range(4)]
    [7, 4, 8, 12]
    """
    dim = 3
    _, Jinv = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drprime_dr)

    Jac4inv = ixp.zerorank2(dimension=4)
    Jac4inv[0][0] = sp.Integer(1)
    for a in range(dim):
        for i in range(dim):
            Jac4inv[a + 1][i + 1] = Jinv[a][i]

    dst = ixp.zerorank2(dimension=4)
    for mu in range(4):
        for nu in range(4):
            for a in range(4):
                for b in range(4):
                    dst[mu][nu] += Jac4inv[mu][a] * Jac4inv[nu][b] * T4UU_cart[a][b]
    return dst


# ================================================================
# Doctest runner (style as in your original snippet)
# ================================================================
if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
