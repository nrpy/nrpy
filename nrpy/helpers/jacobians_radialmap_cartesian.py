r"""
Core Jacobian (basis) transformation functions.
For Cartesian <-> radially-mapped Cartesian bases (purely radial map).

We consider a radial map on R^3:
Given r = sqrt(δ_ij x^i x^j) and a profile r'(r), define x'^i = x^i * (r'/r).

Define the Jacobian and its inverse by
    J^i{}_a        = ∂x^i / ∂x'^a,
    (J^{-1})^a{}_i = ∂x'^a / ∂x^i,

with closed-form expressions
    J^i{}_a        = δ^i{}_a (r/r') + (x'^i x'^a / r'^2) * (dr/dr' - r/r'),
    (J^{-1})^a{}_i = δ^a{}_i (r'/r) + (x'^a x'^i / r'^2) * (dr'/dr - r'/r),

where dr/dr' = 1 / (dr'/dr).

This helper:
- Accepts any radial profile via callables rprime_of_r(r) and drprime_dr(r).
- Does not register CodeParameters.
- Avoids inverting r'(r) unless you pass r or a callable r_of_rprime(r').

Reference (fisheye profile example):
J. A. Faber, T. W. Baumgarte, Z. B. Etienne, S. L. Shapiro, and K. Taniguchi,
"Relativistic hydrodynamics in the presence of puncture black holes,"
Phys. Rev. D 96, 104035 (2017), Appendix: Fisheye Coordinates.
"""

from __future__ import annotations
from typing import Callable, Optional, Sequence, Tuple

import sympy as sp
import nrpy.indexedexp as ixp


def _resolve_r_rprime_and_drprime(
    xprimeU: Sequence[sp.Expr],
    rprime_of_r: Callable[[sp.Expr], sp.Expr],
    drprime_dr: Callable[[sp.Expr], sp.Expr],
    r: Optional[sp.Expr] = None,
    r_of_rprime: Optional[Callable[[sp.Expr], sp.Expr]] = None,
) -> Tuple[sp.Expr, sp.Expr, sp.Expr]:
    """
    Resolve (r', r, dr'/dr) at a point given x'^i and the radial profile.

    :param xprimeU: Primed Cartesian coordinates [x'^0, x'^1, x'^2].
    :param rprime_of_r: Callable returning r'(r).
    :param drprime_dr: Callable returning dr'/dr evaluated at r.
    :param r: Optional known r at the point (bypasses inversion).
    :param r_of_rprime: Optional callable r(r') for exact inversion when available.
    :return: (rprime, r, drprime_dr_at_r).

    Doctest (scaling r' = 2 r):
    >>> rprime_of_r = lambda rv: 2*rv
    >>> drprime     = lambda rv: sp.Integer(2)
    >>> rp, rr, drp = _resolve_r_rprime_and_drprime([sp.Integer(2), 0, 0], rprime_of_r, drprime)
    >>> (sp.simplify(rp), sp.simplify(rr), sp.simplify(drp))
    (2, 1, 2)

    Doctest (nonlinear r' = r + r^3 at r=1, so r' = 2, dr'/dr = 4):
    >>> rprime_of_r = lambda rv: rv + rv**3
    >>> drprime     = lambda rv: 1 + 3*rv**2
    >>> rp, rr, drp = _resolve_r_rprime_and_drprime([sp.Integer(2), 0, 0], rprime_of_r, drprime, r=sp.Integer(1))
    >>> (sp.simplify(rp), sp.simplify(rr), sp.simplify(drp))
    (2, 1, 4)
    """
    rprime = sp.sqrt(sum(comp * comp for comp in xprimeU))

    if r is None:
        if r_of_rprime is not None:
            r = sp.simplify(r_of_rprime(rprime))
        else:
            # Fallback for pure global scaling: r' = a0 * r => r = r'/a0, with a0 = r'(1)
            a0 = sp.simplify(rprime_of_r(sp.Integer(1)))
            if a0 != 0:
                r = sp.simplify(rprime / a0)
            else:
                r = sp.symbols("r", positive=True)

    drp = sp.simplify(drprime_dr(r))
    return rprime, r, drp


def build_radialmap_jacobians_from_xprime(
    xprimeU: Sequence[sp.Expr],
    rprime_of_r: Callable[[sp.Expr], sp.Expr],
    drprime_dr: Callable[[sp.Expr], sp.Expr],
    r: Optional[sp.Expr] = None,
    r_of_rprime: Optional[Callable[[sp.Expr], sp.Expr]] = None,
) -> Tuple[Sequence[Sequence[sp.Expr]], Sequence[Sequence[sp.Expr]]]:
    """
    Construct J^i{}_a = ∂x^i/∂x'^a and (J^{-1})^a{}_i = ∂x'^a/∂x^i for a radial map.

    :param xprimeU: Primed Cartesian coordinates [x'^0, x'^1, x'^2] at the point.
    :param rprime_of_r: Callable returning r'(r).
    :param drprime_dr: Callable returning dr'/dr(r).
    :param r: Optional r at the point (bypasses inversion).
    :param r_of_rprime: Optional callable r(r') for exact inversion when available.
    :return: (J, Jinv) as two 3x3 arrays with entries J[i][a], Jinv[a][i].

    Doctest (r' = 2 r => J = (1/2)I, Jinv = 2I):
    >>> rprime_of_r = lambda rv: 2*rv
    >>> drprime     = lambda rv: sp.Integer(2)
    >>> J, Jinv = build_radialmap_jacobians_from_xprime([sp.Integer(2), 0, 0], rprime_of_r, drprime)
    >>> [sp.simplify(J[i][i]) for i in range(3)]
    [1/2, 1/2, 1/2]
    >>> [sp.simplify(Jinv[i][i]) for i in range(3)]
    [2, 2, 2]

    Doctest (nonlinear r' = r + r^3 at r=1; expect J=diag(1/4,1/2,1/2), Jinv=diag(4,2,2)):
    >>> rprime_of_r = lambda rv: rv + rv**3
    >>> drprime     = lambda rv: 1 + 3*rv**2
    >>> J, Jinv = build_radialmap_jacobians_from_xprime([sp.Integer(2), 0, 0], rprime_of_r, drprime, r=sp.Integer(1))
    >>> [sp.simplify(J[i][i]) for i in range(3)]
    [1/4, 1/2, 1/2]
    >>> [sp.simplify(Jinv[i][i]) for i in range(3)]
    [4, 2, 2]
    """
    dim = 3
    rprime, r_phys, drp_dr = _resolve_r_rprime_and_drprime(
        xprimeU, rprime_of_r, drprime_dr, r=r, r_of_rprime=r_of_rprime
    )

    J = ixp.zerorank2(dimension=dim)     # rows i, cols a
    Jinv = ixp.zerorank2(dimension=dim)  # rows a, cols i

    dr_drprime = sp.simplify(sp.Integer(1) / drp_dr)
    s1 = sp.simplify(r_phys / rprime)
    s2 = sp.simplify(dr_drprime - r_phys / rprime)
    s1i = sp.simplify(rprime / r_phys)
    s2i = sp.simplify(drp_dr - rprime / r_phys)
    inv_rp2 = sp.simplify(sp.Integer(1) / (rprime**2))

    for i in range(dim):
        for a in range(dim):
            delta = sp.Integer(1) if i == a else sp.Integer(0)
            J[i][a] = delta * s1 + inv_rp2 * xprimeU[i] * xprimeU[a] * s2

    for a in range(dim):
        for i in range(dim):
            delta = sp.Integer(1) if a == i else sp.Integer(0)
            Jinv[a][i] = delta * s1i + inv_rp2 * xprimeU[a] * xprimeU[i] * s2i

    return J, Jinv


def basis_transform_vectorU_from_radialmap_to_Cartesian(
    xprimeU: Sequence[sp.Expr],
    src_vectorU_prime: Sequence[sp.Expr],
    rprime_of_r: Callable[[sp.Expr], sp.Expr],
    drprime_dr: Callable[[sp.Expr], sp.Expr],
    r: Optional[sp.Expr] = None,
    r_of_rprime: Optional[Callable[[sp.Expr], sp.Expr]] = None,
) -> Sequence[sp.Expr]:
    """
    Transform a contravariant vector from the radial-mapped basis to the Cartesian basis.

    :param xprimeU: Primed Cartesian coordinates [x'^0, x'^1, x'^2] at the point.
    :param src_vectorU_prime: Source contravariant vector v'^a (radial-mapped basis).
    :param rprime_of_r: Callable returning r'(r).
    :param drprime_dr: Callable returning dr'/dr(r).
    :param r: Optional r at the point.
    :param r_of_rprime: Optional callable r(r') for exact inversion when available.
    :return: Transformed contravariant vector v^i in Cartesian basis.

    Doctest (r' = 2 r):
    >>> rp = lambda rv: 2*rv
    >>> drp = lambda rv: sp.Integer(2)
    >>> v = basis_transform_vectorU_from_radialmap_to_Cartesian([sp.Integer(2),0,0],[3,4,5], rp, drp)
    >>> [sp.simplify(comp) for comp in v]
    [3/2, 2, 5/2]

    Doctest (nonlinear r' = r + r^3 at r=1; J=diag(1/4,1/2,1/2)):
    >>> rp = lambda rv: rv + rv**3
    >>> drp = lambda rv: 1 + 3*rv**2
    >>> v = basis_transform_vectorU_from_radialmap_to_Cartesian([sp.Integer(2),0,0],[3,4,5], rp, drp, r=sp.Integer(1))
    >>> [sp.simplify(comp) for comp in v]
    [3/4, 2, 5/2]
    """
    J, _ = build_radialmap_jacobians_from_xprime(xprimeU, rprime_of_r, drprime_dr, r, r_of_rprime)
    dim = 3
    dst = ixp.zerorank1(dimension=dim)
    for i in range(dim):
        for a in range(dim):
            dst[i] += J[i][a] * src_vectorU_prime[a]
    return dst


def basis_transform_vectorD_from_radialmap_to_Cartesian(
    xprimeU: Sequence[sp.Expr],
    src_vectorD_prime: Sequence[sp.Expr],
    rprime_of_r: Callable[[sp.Expr], sp.Expr],
    drprime_dr: Callable[[sp.Expr], sp.Expr],
    r: Optional[sp.Expr] = None,
    r_of_rprime: Optional[Callable[[sp.Expr], sp.Expr]] = None,
) -> Sequence[sp.Expr]:
    """
    Transform a covariant vector from the radial-mapped basis to the Cartesian basis.

    :param xprimeU: Primed Cartesian coordinates [x'^0, x'^1, x'^2] at the point.
    :param src_vectorD_prime: Source covariant vector ω'_a (radial-mapped basis).
    :param rprime_of_r: Callable returning r'(r).
    :param drprime_dr: Callable returning dr'/dr(r).
    :param r: Optional r at the point.
    :param r_of_rprime: Optional callable r(r') for exact inversion when available.
    :return: Transformed covariant vector ω_i in Cartesian basis.

    Doctest (r' = 2 r):
    >>> rp = lambda rv: 2*rv
    >>> drp = lambda rv: sp.Integer(2)
    >>> w = basis_transform_vectorD_from_radialmap_to_Cartesian([sp.Integer(2),0,0],[1,-2,3], rp, drp)
    >>> [sp.simplify(comp) for comp in w]
    [2, -4, 6]

    Doctest (nonlinear r' = r + r^3 at r=1; Jinv=diag(4,2,2)):
    >>> rp = lambda rv: rv + rv**3
    >>> drp = lambda rv: 1 + 3*rv**2
    >>> w = basis_transform_vectorD_from_radialmap_to_Cartesian([sp.Integer(2),0,0],[1,-2,3], rp, drp, r=sp.Integer(1))
    >>> [sp.simplify(comp) for comp in w]
    [4, -4, 6]
    """
    _, Jinv = build_radialmap_jacobians_from_xprime(xprimeU, rprime_of_r, drprime_dr, r, r_of_rprime)
    dim = 3
    dst = ixp.zerorank1(dimension=dim)
    for i in range(dim):
        for a in range(dim):
            dst[i] += Jinv[a][i] * src_vectorD_prime[a]
    return dst


def basis_transform_tensorDD_from_radialmap_to_Cartesian(
    xprimeU: Sequence[sp.Expr],
    src_tensorDD_prime: Sequence[Sequence[sp.Expr]],
    rprime_of_r: Callable[[sp.Expr], sp.Expr],
    drprime_dr: Callable[[sp.Expr], sp.Expr],
    r: Optional[sp.Expr] = None,
    r_of_rprime: Optional[Callable[[sp.Expr], sp.Expr]] = None,
) -> Sequence[Sequence[sp.Expr]]:
    """
    Transform a rank-2 covariant tensor from the radial-mapped basis to the Cartesian basis.

    :param xprimeU: Primed Cartesian coordinates [x'^0, x'^1, x'^2] at the point.
    :param src_tensorDD_prime: Source covariant tensor T'_{ab} (radial-mapped basis).
    :param rprime_of_r: Callable returning r'(r).
    :param drprime_dr: Callable returning dr'/dr(r).
    :param r: Optional r at the point.
    :param r_of_rprime: Optional callable r(r') for exact inversion when available.
    :return: Transformed covariant tensor T_{ij} in Cartesian basis.

    Doctest (r' = 2 r):
    >>> rp = lambda rv: 2*rv
    >>> drp = lambda rv: sp.Integer(2)
    >>> Tprime = [[sp.Integer(1), 2, 0],[2, 5, -1],[0, -1, 4]]
    >>> T = basis_transform_tensorDD_from_radialmap_to_Cartesian([sp.Integer(2),0,0], Tprime, rp, drp)
    >>> [[sp.simplify(T[i][j]) for j in range(3)] for i in range(3)]
    [[4, 8, 0], [8, 20, -4], [0, -4, 16]]

    Doctest (nonlinear r' = r + r^3 at r=1; Jinv=diag(4,2,2)):
    >>> rp = lambda rv: rv + rv**3
    >>> drp = lambda rv: 1 + 3*rv**2
    >>> Tprime = [[sp.Integer(1), 2, 0],[2, 5, -1],[0, -1, 4]]
    >>> T = basis_transform_tensorDD_from_radialmap_to_Cartesian([sp.Integer(2),0,0], Tprime, rp, drp, r=sp.Integer(1))
    >>> [[sp.simplify(T[i][j]) for j in range(3)] for i in range(3)]
    [[16, 16, 0], [16, 20, -4], [0, -4, 16]]
    """
    _, Jinv = build_radialmap_jacobians_from_xprime(xprimeU, rprime_of_r, drprime_dr, r, r_of_rprime)
    dim = 3
    dst = ixp.zerorank2(dimension=dim)
    for i in range(dim):
        for j in range(dim):
            for a in range(dim):
                for b in range(dim):
                    dst[i][j] += Jinv[a][i] * Jinv[b][j] * src_tensorDD_prime[a][b]
    return dst


def basis_transform_vectorU_from_Cartesian_to_radialmap(
    xprimeU: Sequence[sp.Expr],
    src_vectorU_cart: Sequence[sp.Expr],
    rprime_of_r: Callable[[sp.Expr], sp.Expr],
    drprime_dr: Callable[[sp.Expr], sp.Expr],
    r: Optional[sp.Expr] = None,
    r_of_rprime: Optional[Callable[[sp.Expr], sp.Expr]] = None,
) -> Sequence[sp.Expr]:
    """
    Transform a contravariant vector from the Cartesian basis to the radial-mapped basis.

    :param xprimeU: Primed Cartesian coordinates [x'^0, x'^1, x'^2] at the point.
    :param src_vectorU_cart: Source contravariant vector v^i (Cartesian basis).
    :param rprime_of_r: Callable returning r'(r).
    :param drprime_dr: Callable returning dr'/dr(r).
    :param r: Optional r at the point.
    :param r_of_rprime: Optional callable r(r') for exact inversion when available.
    :return: Transformed contravariant vector v'^a (radial-mapped basis).

    Doctest (r' = 2 r):
    >>> rp = lambda rv: 2*rv
    >>> drp = lambda rv: sp.Integer(2)
    >>> vprime = basis_transform_vectorU_from_Cartesian_to_radialmap([sp.Integer(2),0,0],[sp.Rational(3,2),2,sp.Rational(5,2)], rp, drp)
    >>> [sp.simplify(comp) for comp in vprime]
    [3, 4, 5]

    Doctest (nonlinear r' = r + r^3 at r=1; Jinv=diag(4,2,2)):
    >>> rp = lambda rv: rv + rv**3
    >>> drp = lambda rv: 1 + 3*rv**2
    >>> vprime = basis_transform_vectorU_from_Cartesian_to_radialmap([sp.Integer(2),0,0],[sp.Rational(3,4), 2, sp.Rational(5,2)], rp, drp, r=sp.Integer(1))
    >>> [sp.simplify(comp) for comp in vprime]
    [3, 4, 5]
    """
    _, Jinv = build_radialmap_jacobians_from_xprime(xprimeU, rprime_of_r, drprime_dr, r, r_of_rprime)
    dim = 3
    dst = ixp.zerorank1(dimension=dim)
    for a in range(dim):
        for i in range(dim):
            dst[a] += Jinv[a][i] * src_vectorU_cart[i]
    return dst


def basis_transform_vectorD_from_Cartesian_to_radialmap(
    xprimeU: Sequence[sp.Expr],
    src_vectorD_cart: Sequence[sp.Expr],
    rprime_of_r: Callable[[sp.Expr], sp.Expr],
    drprime_dr: Callable[[sp.Expr], sp.Expr],
    r: Optional[sp.Expr] = None,
    r_of_rprime: Optional[Callable[[sp.Expr], sp.Expr]] = None,
) -> Sequence[sp.Expr]:
    """
    Transform a covariant vector from the Cartesian basis to the radial-mapped basis.

    :param xprimeU: Primed Cartesian coordinates [x'^0, x'^1, x'^2] at the point.
    :param src_vectorD_cart: Source covariant vector ω_i (Cartesian basis).
    :param rprime_of_r: Callable returning r'(r).
    :param drprime_dr: Callable returning dr'/dr(r).
    :param r: Optional r at the point.
    :param r_of_rprime: Optional callable r(r') for exact inversion when available.
    :return: Transformed covariant vector ω'_a (radial-mapped basis).

    Doctest (r' = 2 r):
    >>> rp = lambda rv: 2*rv
    >>> drp = lambda rv: sp.Integer(2)
    >>> wprime = basis_transform_vectorD_from_Cartesian_to_radialmap([sp.Integer(2),0,0],[2,-4,6], rp, drp)
    >>> [sp.simplify(comp) for comp in wprime]
    [1, -2, 3]

    Doctest (nonlinear r' = r + r^3 at r=1; J=diag(1/4,1/2,1/2)):
    >>> rp = lambda rv: rv + rv**3
    >>> drp = lambda rv: 1 + 3*rv**2
    >>> wprime = basis_transform_vectorD_from_Cartesian_to_radialmap([sp.Integer(2),0,0],[4,-4,6], rp, drp, r=sp.Integer(1))
    >>> [sp.simplify(comp) for comp in wprime]
    [1, -2, 3]
    """
    J, _ = build_radialmap_jacobians_from_xprime(xprimeU, rprime_of_r, drprime_dr, r, r_of_rprime)
    dim = 3
    dst = ixp.zerorank1(dimension=dim)
    for a in range(dim):
        for i in range(dim):
            dst[a] += J[i][a] * src_vectorD_cart[i]
    return dst


def basis_transform_tensorDD_from_Cartesian_to_radialmap(
    xprimeU: Sequence[sp.Expr],
    src_tensorDD_cart: Sequence[Sequence[sp.Expr]],
    rprime_of_r: Callable[[sp.Expr], sp.Expr],
    drprime_dr: Callable[[sp.Expr], sp.Expr],
    r: Optional[sp.Expr] = None,
    r_of_rprime: Optional[Callable[[sp.Expr], sp.Expr]] = None,
) -> Sequence[Sequence[sp.Expr]]:
    """
    Transform a rank-2 covariant tensor from the Cartesian basis to the radial-mapped basis.

    :param xprimeU: Primed Cartesian coordinates [x'^0, x'^1, x'^2] at the point.
    :param src_tensorDD_cart: Source covariant tensor T_{ij} (Cartesian basis).
    :param rprime_of_r: Callable returning r'(r).
    :param drprime_dr: Callable returning dr'/dr(r).
    :param r: Optional r at the point.
    :param r_of_rprime: Optional callable r(r') for exact inversion when available.
    :return: Transformed covariant tensor T'_{ab} (radial-mapped basis).

    Doctest (r' = 2 r):
    >>> rp = lambda rv: 2*rv
    >>> drp = lambda rv: sp.Integer(2)
    >>> Tcart = [[sp.Integer(4), 8, 0],[8, 20, -4],[0, -4, 16]]
    >>> Tprime = basis_transform_tensorDD_from_Cartesian_to_radialmap([sp.Integer(2),0,0], Tcart, rp, drp)
    >>> [[sp.simplify(Tprime[i][j]) for j in range(3)] for i in range(3)]
    [[1, 2, 0], [2, 5, -1], [0, -1, 4]]

    Doctest (nonlinear r' = r + r^3 at r=1; J=diag(1/4,1/2,1/2)):
    >>> rp = lambda rv: rv + rv**3
    >>> drp = lambda rv: 1 + 3*rv**2
    >>> Tcart = [[sp.Integer(16), 16, 0],[16, 20, -4],[0, -4, 16]]
    >>> Tprime = basis_transform_tensorDD_from_Cartesian_to_radialmap([sp.Integer(2),0,0], Tcart, rp, drp, r=sp.Integer(1))
    >>> [[sp.simplify(Tprime[i][j]) for j in range(3)] for i in range(3)]
    [[1, 2, 0], [2, 5, -1], [0, -1, 4]]
    """
    J, _ = build_radialmap_jacobians_from_xprime(xprimeU, rprime_of_r, drprime_dr, r, r_of_rprime)
    dim = 3
    dst = ixp.zerorank2(dimension=dim)
    for a in range(dim):
        for b in range(dim):
            for i in range(dim):
                for j in range(dim):
                    dst[a][b] += J[i][a] * J[j][b] * src_tensorDD_cart[i][j]
    return dst


def basis_transform_4tensorUU_from_time_indep_radialmap_to_Cartesian(
    xprimeU: Sequence[sp.Expr],
    T4UU_prime: Sequence[Sequence[sp.Expr]],
    rprime_of_r: Callable[[sp.Expr], sp.Expr],
    drprime_dr: Callable[[sp.Expr], sp.Expr],
    r: Optional[sp.Expr] = None,
    r_of_rprime: Optional[Callable[[sp.Expr], sp.Expr]] = None,
) -> Sequence[Sequence[sp.Expr]]:
    """
    Transform a contravariant 4-tensor (time-independent spatial map) from the radial-mapped
    basis to the Cartesian basis.

    :param xprimeU: Primed Cartesian coordinates [x'^0, x'^1, x'^2] at the point.
    :param T4UU_prime: Source 4-tensor T'^{μν} (contravariant) in the radial-mapped basis.
    :param rprime_of_r: Callable returning r'(r).
    :param drprime_dr: Callable returning dr'/dr(r).
    :param r: Optional r at the point.
    :param r_of_rprime: Optional callable r(r') for exact inversion when available.
    :return: Transformed 4-tensor T^{μν} (contravariant) in Cartesian basis.

    Doctest (r' = 2 r; spatial J = (1/2)I):
    >>> rp = lambda rv: 2*rv
    >>> drp = lambda rv: sp.Integer(2)
    >>> T4p = [[sp.Integer(7),0,0,0],[0,sp.Integer(2),0,0],[0,0,sp.Integer(4),0],[0,0,0,sp.Integer(6)]]
    >>> T4  = basis_transform_4tensorUU_from_time_indep_radialmap_to_Cartesian([sp.Integer(2),0,0], T4p, rp, drp)
    >>> [sp.simplify(T4[i][i]) for i in range(4)]
    [7, 1/2, 1, 3/2]

    Doctest (nonlinear r' = r + r^3 at r=1; spatial J = diag(1/4,1/2,1/2)):
    >>> rp = lambda rv: rv + rv**3
    >>> drp = lambda rv: 1 + 3*rv**2
    >>> T4p = [[sp.Integer(7),0,0,0],[0,sp.Integer(2),0,0],[0,0,sp.Integer(4),0],[0,0,0,sp.Integer(6)]]
    >>> T4  = basis_transform_4tensorUU_from_time_indep_radialmap_to_Cartesian([sp.Integer(2),0,0], T4p, rp, drp, r=sp.Integer(1))
    >>> [sp.nsimplify(T4[i][i]) for i in range(4)]
    [7, 1/8, 1, 3/2]
    """
    J, _ = build_radialmap_jacobians_from_xprime(xprimeU, rprime_of_r, drprime_dr, r, r_of_rprime)

    dim = 3
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


def basis_transform_4tensorUU_from_Cartesian_to_time_indep_radialmap(
    xprimeU: Sequence[sp.Expr],
    T4UU_cart: Sequence[Sequence[sp.Expr]],
    rprime_of_r: Callable[[sp.Expr], sp.Expr],
    drprime_dr: Callable[[sp.Expr], sp.Expr],
    r: Optional[sp.Expr] = None,
    r_of_rprime: Optional[Callable[[sp.Expr], sp.Expr]] = None,
) -> Sequence[Sequence[sp.Expr]]:
    """
    Transform a contravariant 4-tensor (time-independent spatial map) from the Cartesian basis
    to the radial-mapped basis.

    :param xprimeU: Primed Cartesian coordinates [x'^0, x'^1, x'^2] at the point.
    :param T4UU_cart: Source 4-tensor T^{μν} (contravariant) in Cartesian basis.
    :param rprime_of_r: Callable returning r'(r).
    :param drprime_dr: Callable returning dr'/dr(r).
    :param r: Optional r at the point.
    :param r_of_rprime: Optional callable r(r') for exact inversion when available.
    :return: Transformed 4-tensor T'^{μν} (contravariant) in the radial-mapped basis.

    Doctest (r' = 2 r; spatial Jinv = 2I):
    >>> rp = lambda rv: 2*rv
    >>> drp = lambda rv: sp.Integer(2)
    >>> T4c = [[sp.Integer(7),0,0,0],[0,sp.Integer(1),0,0],[0,0,sp.Integer(2),0],[0,0,0,sp.Integer(3)]]
    >>> T4p = basis_transform_4tensorUU_from_Cartesian_to_time_indep_radialmap([sp.Integer(2),0,0], T4c, rp, drp)
    >>> [sp.simplify(T4p[i][i]) for i in range(4)]
    [7, 4, 8, 12]

    Doctest (nonlinear r' = r + r^3 at r=1; spatial Jinv = diag(4,2,2)):
    >>> rp = lambda rv: rv + rv**3
    >>> drp = lambda rv: 1 + 3*rv**2
    >>> T4c = [[sp.Integer(7),0,0,0],[0,sp.Rational(1,8),0,0],[0,0,sp.Integer(1),0],[0,0,0,sp.Rational(3,2)]]
    >>> T4p = basis_transform_4tensorUU_from_Cartesian_to_time_indep_radialmap([sp.Integer(2),0,0], T4c, rp, drp, r=sp.Integer(1))
    >>> [sp.nsimplify(T4p[i][i]) for i in range(4)]
    [7, 2, 4, 6]
    """
    _, Jinv = build_radialmap_jacobians_from_xprime(xprimeU, rprime_of_r, drprime_dr, r, r_of_rprime)

    dim = 3
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


##################################################
if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
