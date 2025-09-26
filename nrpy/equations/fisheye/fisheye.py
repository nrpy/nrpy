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
import nrpy.params as par  # NRPy CodeParameters

# ---------------------------------------------------------------------
# Register Fisheye CodeParameters (max length 5 for each array)
# n = number of transitions (0..4). We use the first n+1 entries of a,
# and the first n entries of R and s.
# ---------------------------------------------------------------------
_FISH_MAX = 5

par.register_CodeParameter(
    cparam_type="INT",
    module=__name__,
    name="fisheye_n",
    defaultvalue=0,
    description="Number of fisheye transition zones n (0..4).",
    add_to_parfile=True,
    commondata=True,
    add_to_set_CodeParameters_h=False,
    assumption="IntNonNegative",
)

par.register_CodeParameter(
    cparam_type=f"REAL[{_FISH_MAX}]",
    module=__name__,
    name="fisheye_a",
    defaultvalue=[1, 1, 1, 1, 1],
    description=(
        "Fisheye coefficients [a_0, ..., a_k]; use entries 0..n inclusive "
        "for algorithm that needs a_0..a_n (with n <= 4)."
    ),
    add_to_parfile=True,
    commondata=True,
    add_to_set_CodeParameters_h=False,
    assumption="Real",
)

par.register_CodeParameter(
    cparam_type=f"REAL[{_FISH_MAX}]",
    module=__name__,
    name="fisheye_R",
    defaultvalue=[0, 0, 0, 0, 0],
    description="Fisheye transition centers [R_1, ..., R_n] (use first n entries).",
    add_to_parfile=True,
    commondata=True,
    add_to_set_CodeParameters_h=False,
    assumption="Real",
)

par.register_CodeParameter(
    cparam_type=f"REAL[{_FISH_MAX}]",
    module=__name__,
    name="fisheye_s",
    defaultvalue=[1, 1, 1, 1, 1],
    description="Fisheye transition widths [s_1, ..., s_n] (use first n entries).",
    add_to_parfile=True,
    commondata=True,
    add_to_set_CodeParameters_h=False,
    assumption="RealPositive",
)


# ================================================================
# Fisheye profile r'(r) and its radial derivative dr'/dr (helpers)
# ================================================================


def fisheye_rprime_of_r(r: sp.Expr) -> sp.Expr:
    r"""
    Compute r'(r) for the radial fisheye profile using NRPy CodeParameters.

    :param r: Physical radius \(r\).
    :return: The fisheye radius \(r'(r)\) as a SymPy expression.

    Doctests:
    >>> # Configure a pure global scale: n=0, a_0=2 -> r' = 2*r
    >>> par.set_parval_from_str("fisheye_n", 0)
    >>> par.set_parval_from_str("fisheye_a", [2,1,1,1,1])
    >>> par.set_parval_from_str("fisheye_R", [0,0,0,0,0])
    >>> par.set_parval_from_str("fisheye_s", [1,1,1,1,1])
    >>> fisheye_rprime_of_r(sp.Integer(3))
    6
    """
    n = int(par.parval_from_str("fisheye_n"))
    if not (0 <= n <= _FISH_MAX - 1):
        raise ValueError("fisheye_n must be in [0, 4].")
    a_arr = list(par.parval_from_str("fisheye_a"))
    R_arr = list(par.parval_from_str("fisheye_R"))
    s_arr = list(par.parval_from_str("fisheye_s"))

    rp = sp.sympify(a_arr[n]) * r  # a_n * r
    for i in range(1, n + 1):
        a_im1 = sp.sympify(a_arr[i - 1])
        a_i = sp.sympify(a_arr[i])
        R_i = sp.sympify(R_arr[i - 1])
        s_i = sp.sympify(s_arr[i - 1])
        denom = 2 * sp.tanh(R_i / s_i)
        term = sp.log(sp.cosh((r + R_i) / s_i) / sp.cosh((r - R_i) / s_i))
        rp += (a_im1 - a_i) * s_i / denom * term
    return sp.simplify(rp)


def fisheye_drprime_dr(r: sp.Expr) -> sp.Expr:
    r"""
    Compute \(dr'/dr\) for the radial fisheye profile using NRPy CodeParameters.

    :param r: Physical radius \(r\) (symbolic allowed).
    :return: The radial derivative \(dr'/dr\) as a SymPy expression.

    Doctests:
    >>> # Pure global scale: n=0, a_0=2 -> dr'/dr = 2
    >>> par.set_parval_from_str("fisheye_n", 0)
    >>> par.set_parval_from_str("fisheye_a", [2,1,1,1,1])
    >>> fisheye_drprime_dr(sp.symbols("r", positive=True))
    2
    """
    n = int(par.parval_from_str("fisheye_n"))
    if not (0 <= n <= _FISH_MAX - 1):
        raise ValueError("fisheye_n must be in [0, 4].")
    a_arr = list(par.parval_from_str("fisheye_a"))
    R_arr = list(par.parval_from_str("fisheye_R"))
    s_arr = list(par.parval_from_str("fisheye_s"))

    drp = sp.sympify(a_arr[n])
    for i in range(1, n + 1):
        a_im1 = sp.sympify(a_arr[i - 1])
        a_i = sp.sympify(a_arr[i])
        R_i = sp.sympify(R_arr[i - 1])
        s_i = sp.sympify(s_arr[i - 1])
        coeff = (a_im1 - a_i) / sp.tanh(R_i / s_i)
        bracket = (sp.tanh((r + R_i) / s_i) - sp.tanh((r - R_i) / s_i)) / 2
        drp += coeff * bracket
    return sp.simplify(drp)


# ================================================================
# Core helpers
# ================================================================

def _compute_r_rprime_and_drprime(xprimeU: Sequence[sp.Expr]) -> Tuple[sp.Expr, sp.Expr, sp.Expr]:
    r"""
    Given fisheye coordinates \(x'^i\), compute \(r'=\sqrt{x'^i x'^i}\),
    the corresponding physical radius \(r\), and \(dr'/dr\).

    Notes
    -----
    - For general fisheye with transitions (n>0), inverting \(r'(r)\) may
      require a numeric solve; here we keep doctests to n=0 where
      \(r' = a_0 r\) and inversion is analytic.
    - If symbolic inversion is required, we return a positive symbol `r`
      and keep expressions in terms of `r`.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\).
    :return: Tuple \((rprime, r, drprime\_dr)\).
    """
    rprime = sp.sqrt(sum(comp * comp for comp in xprimeU))
    n = int(par.parval_from_str("fisheye_n"))
    a_arr = list(par.parval_from_str("fisheye_a"))

    # If n=0, r' = a0 * r  => r = r'/a0
    if n == 0:
        a0 = sp.sympify(a_arr[0])
        r = sp.simplify(rprime / a0)
        drp = sp.sympify(a_arr[0])
        return rprime, r, drp

    # General case: attempt symbolic fallback
    r_sym = sp.symbols("r", positive=True)
    drp = fisheye_drprime_dr(r_sym)
    return rprime, r_sym, drp


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
    >>> # Global scaling: n=0, a0=2 -> r' = 2 r, J = (1/2) I, Jinv = 2 I.
    >>> par.set_parval_from_str("fisheye_n", 0)
    >>> par.set_parval_from_str("fisheye_a", [2,1,1,1,1])
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]  # r' = 2
    >>> rprime, r, drpdr = _compute_r_rprime_and_drprime(xprimeU)
    >>> J, Jinv = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drpdr)
    >>> [sp.simplify(J[i][i]) for i in range(3)]
    [1/2, 1/2, 1/2]
    >>> [sp.simplify(Jinv[i][i]) for i in range(3)]
    [2, 2, 2]
    """
    dim = 3
    J = ixp.zerorank2(dimension=dim)   # rows i, cols a
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
    src_vectorU_prime: Sequence[sp.Expr],
) -> Sequence[sp.Expr]:
    r"""
    Transform a contravariant vector from fisheye basis to Cartesian.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\) of the evaluation point.
    :param src_vectorU_prime: Source contravariant vector \(v'^a\) in fisheye basis.
    :return: Transformed contravariant vector \(v^i\) in Cartesian basis.

    Doctests:
    >>> # Set fisheye: n=0, a0=2 -> r' = 2 r (pure scaling)
    >>> par.set_parval_from_str("fisheye_n", 0)
    >>> par.set_parval_from_str("fisheye_a", [2,1,1,1,1])
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> vprime = [sp.Integer(3), sp.Integer(4), sp.Integer(5)]
    >>> v = basis_transform_vectorU_from_fisheye_to_Cartesian(xprimeU, vprime)
    >>> [sp.simplify(comp) for comp in v]
    [3/2, 2, 5/2]
    """
    rprime, r, drpdr = _compute_r_rprime_and_drprime(xprimeU)
    J, _ = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drpdr)
    dim = 3
    dst = ixp.zerorank1(dimension=dim)
    for i in range(dim):
        for aidx in range(dim):
            dst[i] += J[i][aidx] * src_vectorU_prime[aidx]
    return dst


def basis_transform_vectorD_from_fisheye_to_Cartesian(
    xprimeU: Sequence[sp.Expr],
    src_vectorD_prime: Sequence[sp.Expr],
) -> Sequence[sp.Expr]:
    r"""
    Transform a covariant vector from fisheye basis to Cartesian.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\) of the evaluation point.
    :param src_vectorD_prime: Source covariant vector \(\omega'_a\) in fisheye basis.
    :return: Transformed covariant vector \(\omega_i\) in Cartesian basis.

    Doctests:
    >>> par.set_parval_from_str("fisheye_n", 0)
    >>> par.set_parval_from_str("fisheye_a", [2,1,1,1,1])
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> wprime = [sp.Integer(1), sp.Integer(-2), sp.Integer(3)]
    >>> w = basis_transform_vectorD_from_fisheye_to_Cartesian(xprimeU, wprime)
    >>> [sp.simplify(comp) for comp in w]
    [2, -4, 6]
    """
    rprime, r, drpdr = _compute_r_rprime_and_drprime(xprimeU)
    _, Jinv = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drpdr)
    dim = 3
    dst = ixp.zerorank1(dimension=dim)
    for i in range(dim):
        for aidx in range(dim):
            dst[i] += Jinv[aidx][i] * src_vectorD_prime[aidx]
    return dst


def basis_transform_tensorDD_from_fisheye_to_Cartesian(
    xprimeU: Sequence[sp.Expr],
    src_tensorDD_prime: Sequence[Sequence[sp.Expr]],
) -> Sequence[Sequence[sp.Expr]]:
    r"""
    Transform a rank-2 covariant tensor from fisheye basis to Cartesian.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\) of the evaluation point.
    :param src_tensorDD_prime: Source covariant tensor \(T'_{ab}\) in fisheye basis.
    :return: Transformed covariant tensor \(T_{ij}\) in Cartesian basis.

    Doctests:
    >>> par.set_parval_from_str("fisheye_n", 0)
    >>> par.set_parval_from_str("fisheye_a", [2,1,1,1,1])
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> Tprime = [[sp.Integer(1), sp.Integer(2), sp.Integer(0)],
    ...           [sp.Integer(2), sp.Integer(5), sp.Integer(-1)],
    ...           [sp.Integer(0), sp.Integer(-1), sp.Integer(4)]]
    >>> T = basis_transform_tensorDD_from_fisheye_to_Cartesian(xprimeU, Tprime)
    >>> [[sp.simplify(T[i][j]) for j in range(3)] for i in range(3)]
    [[4, 8, 0], [8, 20, -4], [0, -4, 16]]
    """
    rprime, r, drpdr = _compute_r_rprime_and_drprime(xprimeU)
    _, Jinv = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drpdr)
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
    src_vectorU_cart: Sequence[sp.Expr],
) -> Sequence[sp.Expr]:
    r"""
    Transform a contravariant vector from Cartesian to fisheye basis.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\) of the evaluation point.
    :param src_vectorU_cart: Source contravariant vector \(v^i\) in Cartesian basis.
    :return: Transformed contravariant vector \(v'^a\) in fisheye basis.

    Doctests:
    >>> par.set_parval_from_str("fisheye_n", 0)
    >>> par.set_parval_from_str("fisheye_a", [2,1,1,1,1])
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> vcart = [sp.Rational(3,2), sp.Integer(2), sp.Rational(5,2)]
    >>> vprime = basis_transform_vectorU_from_Cartesian_to_fisheye(xprimeU, vcart)
    >>> [sp.simplify(comp) for comp in vprime]
    [3, 4, 5]
    """
    rprime, r, drpdr = _compute_r_rprime_and_drprime(xprimeU)
    _, Jinv = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drpdr)
    dim = 3
    dst = ixp.zerorank1(dimension=dim)
    for aidx in range(dim):
        for i in range(dim):
            dst[aidx] += Jinv[aidx][i] * src_vectorU_cart[i]
    return dst


def basis_transform_vectorD_from_Cartesian_to_fisheye(
    xprimeU: Sequence[sp.Expr],
    src_vectorD_cart: Sequence[sp.Expr],
) -> Sequence[sp.Expr]:
    r"""
    Transform a covariant vector from Cartesian to fisheye basis.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\) of the evaluation point.
    :param src_vectorD_cart: Source covariant vector \(\omega_i\) in Cartesian basis.
    :return: Transformed covariant vector \(\omega'_a\) in fisheye basis.

    Doctests:
    >>> par.set_parval_from_str("fisheye_n", 0)
    >>> par.set_parval_from_str("fisheye_a", [2,1,1,1,1])
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> wcart = [sp.Integer(2), sp.Integer(-4), sp.Integer(6)]
    >>> wprime = basis_transform_vectorD_from_Cartesian_to_fisheye(xprimeU, wcart)
    >>> [sp.simplify(comp) for comp in wprime]
    [1, -2, 3]
    """
    rprime, r, drpdr = _compute_r_rprime_and_drprime(xprimeU)
    J, _ = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drpdr)
    dim = 3
    dst = ixp.zerorank1(dimension=dim)
    for aidx in range(dim):
        for i in range(dim):
            dst[aidx] += J[i][aidx] * src_vectorD_cart[i]
    return dst


def basis_transform_tensorDD_from_Cartesian_to_fisheye(
    xprimeU: Sequence[sp.Expr],
    src_tensorDD_cart: Sequence[Sequence[sp.Expr]],
) -> Sequence[Sequence[sp.Expr]]:
    r"""
    Transform a rank-2 covariant tensor from Cartesian to fisheye basis.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\) of the evaluation point.
    :param src_tensorDD_cart: Source covariant tensor \(T_{ij}\) in Cartesian basis.
    :return: Transformed covariant tensor \(T'_{ab}\) in fisheye basis.

    Doctests:
    >>> par.set_parval_from_str("fisheye_n", 0)
    >>> par.set_parval_from_str("fisheye_a", [2,1,1,1,1])
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> Tcart = [[sp.Integer(4), sp.Integer(8), sp.Integer(0)],
    ...          [sp.Integer(8), sp.Integer(20), sp.Integer(-4)],
    ...          [sp.Integer(0), sp.Integer(-4), sp.Integer(16)]]
    >>> Tprime = basis_transform_tensorDD_from_Cartesian_to_fisheye(xprimeU, Tcart)
    >>> [[sp.simplify(Tprime[i][j]) for j in range(3)] for i in range(3)]
    [[1, 2, 0], [2, 5, -1], [0, -1, 4]]
    """
    rprime, r, drpdr = _compute_r_rprime_and_drprime(xprimeU)
    J, _ = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drpdr)
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
    T4UU_prime: Sequence[Sequence[sp.Expr]],
) -> Sequence[Sequence[sp.Expr]]:
    r"""
    Transform a contravariant 4-tensor (time-independent spatial map) from fisheye to Cartesian.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\) of the evaluation point.
    :param T4UU_prime: Source 4-tensor \(T'^{\mu\nu}\) (contravariant) in fisheye basis.
    :return: Transformed 4-tensor \(T^{\mu\nu}\) (contravariant) in Cartesian basis.

    Doctests:
    >>> par.set_parval_from_str("fisheye_n", 0)
    >>> par.set_parval_from_str("fisheye_a", [2,1,1,1,1])
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> T4p = [[sp.Integer(7),0,0,0],[0,sp.Integer(2),0,0],[0,0,sp.Integer(4),0],[0,0,0,sp.Integer(6)]]
    >>> T4 = basis_transform_4tensorUU_from_time_indep_fisheye_to_Cartesian(xprimeU, T4p)
    >>> [sp.simplify(T4[i][i]) for i in range(4)]
    [7, 1/2, 1, 3/2]
    """
    rprime, r, drpdr = _compute_r_rprime_and_drprime(xprimeU)
    J, _ = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drpdr)

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


def basis_transform_4tensorUU_from_Cartesian_to_time_indep_fisheye(
    xprimeU: Sequence[sp.Expr],
    T4UU_cart: Sequence[Sequence[sp.Expr]],
) -> Sequence[Sequence[sp.Expr]]:
    r"""
    Transform a contravariant 4-tensor (time-independent spatial map) from Cartesian to fisheye.

    :param xprimeU: Primed Cartesian coordinates \([x'^0, x'^1, x'^2]\) of the evaluation point.
    :param T4UU_cart: Source 4-tensor \(T^{\mu\nu}\) (contravariant) in Cartesian basis.
    :return: Transformed 4-tensor \(T'^{\mu\nu}\) (contravariant) in fisheye basis.

    Doctests:
    >>> par.set_parval_from_str("fisheye_n", 0)
    >>> par.set_parval_from_str("fisheye_a", [2,1,1,1,1])
    >>> xprimeU = [sp.Integer(2), sp.Integer(0), sp.Integer(0)]
    >>> T4c = [[sp.Integer(7),0,0,0],[0,sp.Integer(1),0,0],[0,0,sp.Integer(2),0],[0,0,0,sp.Integer(3)]]
    >>> T4p = basis_transform_4tensorUU_from_Cartesian_to_time_indep_fisheye(xprimeU, T4c)
    >>> [sp.simplify(T4p[i][i]) for i in range(4)]
    [7, 4, 8, 12]
    """
    rprime, r, drpdr = _compute_r_rprime_and_drprime(xprimeU)
    _, Jinv = build_fisheye_jacobians_from_xprime(xprimeU, rprime, r, drpdr)

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
