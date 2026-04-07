"""Helpers for exact Kasner initial data."""

from __future__ import annotations

from typing import List, Tuple

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.params as par


def kasner_adm_initial_data(
    t0_default: float = 1.0,
    p1_default: float = -1.0 / 3.0,
    p2_default: float = 2.0 / 3.0,
    p3_default: float = 2.0 / 3.0,
) -> Tuple[List[List[sp.Expr]], List[List[sp.Expr]], sp.Expr, List[sp.Expr], List[sp.Expr]]:
    """
    Construct exact Kasner ADM initial data in Cartesian coordinates.

    :return: ``(gammaDD, KDD, alpha, betaU, BU)``.
    """
    t0, p1, p2, p3 = par.register_CodeParameters(
        "REAL",
        __name__,
        ["KASNER_t0", "KASNER_p1", "KASNER_p2", "KASNER_p3"],
        [t0_default, p1_default, p2_default, p3_default],
        commondata=True,
    )

    gammaDD = ixp.zerorank2()
    gammaDD[0][0] = t0 ** (2 * p1)
    gammaDD[1][1] = t0 ** (2 * p2)
    gammaDD[2][2] = t0 ** (2 * p3)

    KDD = ixp.zerorank2()
    KDD[0][0] = -p1 * t0 ** (2 * p1 - 1)
    KDD[1][1] = -p2 * t0 ** (2 * p2 - 1)
    KDD[2][2] = -p3 * t0 ** (2 * p3 - 1)

    alpha = sp.sympify(1)
    betaU = ixp.zerorank1()
    BU = ixp.zerorank1()

    return gammaDD, KDD, alpha, betaU, BU
