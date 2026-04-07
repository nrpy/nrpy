"""
Symbolic helper for exact Kasner diagnostics.

This module builds codegen-ready SymPy expressions for exact Kasner BSSN
quantities in the raw reference-metric basis at physical time
`KASNER_t0 + time`.
"""

from __future__ import annotations

from typing import Dict, List

import sympy as sp

import nrpy.indexedexp as ixp
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN
import nrpy.equations.basis_transforms.jacobians as jac


def kasner_exact_bssn_exprs(CoordSystem: str) -> Dict[str, sp.Expr | List[List[sp.Expr]]]:
    """
    Build exact Kasner BSSN expressions in the raw reference-metric basis.

    :param CoordSystem: Coordinate system used for the BSSN variables.
    :return: Dict with exact ``cf``, ``trK``, ``hDD``, ``AbarDD``, and ``lambdaU`` expressions.
    """
    t_phys = sp.Symbol("KASNER_t_phys", real=True)
    p1 = sp.Symbol("KASNER_p1", real=True)
    p2 = sp.Symbol("KASNER_p2", real=True)
    p3 = sp.Symbol("KASNER_p3", real=True)

    gammaDD = ixp.zerorank2()
    gammaDD[0][0] = t_phys ** (2 * p1)
    gammaDD[1][1] = t_phys ** (2 * p2)
    gammaDD[2][2] = t_phys ** (2 * p3)

    KDD = ixp.zerorank2()
    KDD[0][0] = -p1 * t_phys ** (2 * p1 - 1)
    KDD[1][1] = -p2 * t_phys ** (2 * p2 - 1)
    KDD[2][2] = -p3 * t_phys ** (2 * p3 - 1)

    betaU = ixp.zerorank1()
    BU = ixp.zerorank1()
    adm2bssn = ADM_to_BSSN(
        gammaDD=gammaDD,
        KDD=KDD,
        betaU=betaU,
        BU=BU,
        CoordSystem="Cartesian",
    )

    gammabarDD = jac.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
        CoordSystem, adm2bssn.gammabarDD
    )
    AbarDD = jac.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
        CoordSystem, adm2bssn.AbarDD
    )
    import nrpy.reference_metric as refmetric

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
