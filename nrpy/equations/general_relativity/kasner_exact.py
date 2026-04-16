"""
Exact Kasner expressions used by BHaH initial-data and diagnostics code.

Author: Nishita Jadoo
        njadoo **at** uidaho **dot* com

Metric (synchronous coordinates):
  ds^2 = -dt^2 + t^(2 p1) dx^2 + t^(2 p2) dy^2 + t^(2 p3) dz^2

Kasner constraints:
  p1 + p2 + p3 = 1
  p1^2 + p2^2 + p3^2 = 1

Reference:
  Edward Kasner, "Geometrical Theorems on Einstein's Cosmological Equations",
  American Journal of Mathematics 43 (1921), 217-221.

ADM fields in Cartesian coordinates:
  gamma_ij = diag(t^(2 p1), t^(2 p2), t^(2 p3))
  alpha = 1
  beta^i = 0, B^i = 0
  K_ij = -(1/2) d_t gamma_ij  =>  K_ii = -p_i t^(2 p_i - 1)

BSSN construction pipeline:
  1) Build ADM fields above.
  2) Convert ADM -> BSSN in Cartesian basis via ADM_to_BSSN.
  3) Transform tensors to the raw reference-metric basis.
  4) Compute h_ij = (gammabar_ij - gammahat_ij) / ReDD_ij.
  5) Compute Lambdabar^i from Christoffels and rescale to lambdaU.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN
from nrpy.equations.basis_transforms.jacobians import BasisTransforms
import nrpy.reference_metric as refmetric


def kasner_adm_quantities(
    t_phys: sp.Expr, p1: sp.Expr, p2: sp.Expr, p3: sp.Expr
) -> Tuple[List[List[sp.Expr]], List[List[sp.Expr]], sp.Expr, List[sp.Expr], List[sp.Expr]]:
    """
    Construct exact ADM fields for the Kasner metric.

    :param t_phys: Physical time coordinate ``t`` in the Kasner metric.
    :param p1: Kasner exponent along ``x``.
    :param p2: Kasner exponent along ``y``.
    :param p3: Kasner exponent along ``z``.
    :return: ``(gammaDD, KDD, alpha, betaU, BU)`` in Cartesian basis.

    :Example:
    >>> t, p1, p2, p3 = sp.symbols("t p1 p2 p3", real=True)
    >>> gammaDD, KDD, alpha, betaU, BU = kasner_adm_quantities(t, p1, p2, p3)
    >>> alpha
    1
    >>> gammaDD[0][0]
    t**(2*p1)
    >>> KDD[2][2]
    -p3*t**(2*p3 - 1)
    >>> all(gammaDD[i][j] == 0 for i in range(3) for j in range(3) if i != j)
    True
    >>> all(KDD[i][j] == 0 for i in range(3) for j in range(3) if i != j)
    True
    """
    gammaDD = ixp.zerorank2()
    gammaDD[0][0] = t_phys ** (2 * p1)
    gammaDD[1][1] = t_phys ** (2 * p2)
    gammaDD[2][2] = t_phys ** (2 * p3)

    KDD = ixp.zerorank2()
    KDD[0][0] = -p1 * t_phys ** (2 * p1 - 1)
    KDD[1][1] = -p2 * t_phys ** (2 * p2 - 1)
    KDD[2][2] = -p3 * t_phys ** (2 * p3 - 1)

    alpha = sp.sympify(1)
    betaU = ixp.zerorank1()
    BU = ixp.zerorank1()
    return gammaDD, KDD, alpha, betaU, BU


def kasner_adm_initial_data(
    t0_default: float = 1.0,
    p1_default: float = -1.0 / 3.0,
    p2_default: float = 2.0 / 3.0,
    p3_default: float = 2.0 / 3.0,
) -> Tuple[List[List[sp.Expr]], List[List[sp.Expr]], sp.Expr, List[sp.Expr], List[sp.Expr]]:
    """
    Construct exact Kasner ADM initial data in Cartesian coordinates.

    Registers ``KASNER_t0``, ``KASNER_p1``, ``KASNER_p2``, and ``KASNER_p3`` as commondata
    CodeParameters, and evaluates the exact ADM fields at ``t = KASNER_t0``.

    :return: ``(gammaDD, KDD, alpha, betaU, BU)``.

    :Example:
    >>> gammaDD, KDD, alpha, betaU, BU = kasner_adm_initial_data()
    >>> alpha
    1
    >>> all(gammaDD[i][j] == 0 for i in range(3) for j in range(3) if i != j)
    True
    >>> all(KDD[i][j] == 0 for i in range(3) for j in range(3) if i != j)
    True
    """
    t0, p1, p2, p3 = par.register_CodeParameters(
        "REAL",
        __name__,
        ["KASNER_t0", "KASNER_p1", "KASNER_p2", "KASNER_p3"],
        [t0_default, p1_default, p2_default, p3_default],
        commondata=True,
    )
    return kasner_adm_quantities(t0, p1, p2, p3)


def kasner_exact_bssn_exprs(CoordSystem: str) -> Dict[str, sp.Expr | List[List[sp.Expr]]]:
    """
    Build exact Kasner BSSN expressions in the raw reference-metric basis.

    :param CoordSystem: Coordinate system used for the BSSN variables.
    :return: Dict with exact ``cf``, ``trK``, ``hDD``, ``AbarDD``, and ``lambdaU`` expressions.

    :Example:
    >>> exprs = kasner_exact_bssn_exprs("Cartesian")
    Setting up reference_metric[Cartesian]...
    >>> sorted(exprs.keys())
    ['AbarDD', 'cf', 'hDD', 'lambdaU', 'trK']
    >>> len(exprs["hDD"]) == 3 and len(exprs["hDD"][0]) == 3
    True
    >>> len(exprs["lambdaU"]) == 3
    True
    """
    t_phys = sp.Symbol("KASNER_t_phys", real=True)
    p1 = sp.Symbol("KASNER_p1", real=True)
    p2 = sp.Symbol("KASNER_p2", real=True)
    p3 = sp.Symbol("KASNER_p3", real=True)
    gammaDD, KDD, _alpha, betaU, BU = kasner_adm_quantities(t_phys, p1, p2, p3)
    adm2bssn = ADM_to_BSSN(
        gammaDD=gammaDD,
        KDD=KDD,
        betaU=betaU,
        BU=BU,
        CoordSystem="Cartesian",
    )

    basis_transforms = BasisTransforms(CoordSystem)
    gammabarDD = basis_transforms.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
        adm2bssn.gammabarDD
    )
    AbarDD = basis_transforms.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
        adm2bssn.AbarDD
    )

    rfm = refmetric.reference_metric[CoordSystem]
    hDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            hDD[i][j] = (gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[i][j]

    gammabarUU, _detgammabar = ixp.symm_matrix_inverter3x3(gammabarDD)

    gammabarDD_dD = ixp.zerorank3()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                gammabarDD_dD[i][j][k] = sp.diff(gammabarDD[i][j], rfm.xx[k])

    GammabarUDD = ixp.zerorank3()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for m in range(3):
                    GammabarUDD[i][j][k] += sp.Rational(1, 2) * gammabarUU[i][m] * (
                        gammabarDD_dD[m][j][k]
                        + gammabarDD_dD[m][k][j]
                        - gammabarDD_dD[j][k][m]
                    )

    lambdaU = ixp.zerorank1()
    for i in range(3):
        LambdabarU = sp.sympify(0)
        for j in range(3):
            for k in range(3):
                LambdabarU += gammabarUU[j][k] * (
                    GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k]
                )
        lambdaU[i] = LambdabarU / rfm.ReU[i]

    return {
        "cf": adm2bssn.cf,
        "trK": adm2bssn.trK,
        "hDD": hDD,
        "AbarDD": AbarDD,
        "lambdaU": lambdaU,
    }


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
