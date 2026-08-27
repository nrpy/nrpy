"""
Project BSSN variables onto the algebraic determinant and trace constraints.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List, Tuple

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.reference_metric as refmetric
from nrpy.equations.general_relativity.BSSN_quantities import (
    BSSN_quantities,
    BSSNQuantities,
)
from nrpy.reference_metric import ReferenceMetric


def _project_general_rfm(
    Bq: BSSNQuantities, rfm: ReferenceMetric
) -> Tuple[List[List[sp.Expr]], List[List[sp.Expr]]]:
    """
    Project BSSN variables for a fully general reference metric.

    :param Bq: BSSN quantities to project.
    :param rfm: Fully general reference metric used by the BSSN quantities.
    :return: Projected rescaled conformal metric and extrinsic-curvature tensors.
    """
    hprimeDD_enforced = ixp.zerorank2()
    aprimeDD_enforced = ixp.zerorank2()
    _, detgammabar = ixp.symm_matrix_inverter3x3(Bq.gammabarDD)
    nrpyAbs = sp.Function("nrpyAbs")
    q = (nrpyAbs(rfm.detgammahat) / detgammabar) ** sp.Rational(1, 3)
    gammabarprimeDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            gammabarprimeDD[i][j] = q * Bq.gammabarDD[i][j]
    gammabarprimeUU, _ = ixp.symm_matrix_inverter3x3(gammabarprimeDD)
    trAbarprime = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            trAbarprime += gammabarprimeUU[i][j] * Bq.AbarDD[i][j]
    for i in range(3):
        for j in range(3):
            hprimeDD_enforced[i][j] = (
                gammabarprimeDD[i][j] - rfm.ghatDD[i][j]
            ) / rfm.ReDD[i][j]
            aprimeDD_enforced[i][j] = (
                Bq.AbarDD[i][j]
                - sp.Rational(1, 3) * gammabarprimeDD[i][j] * trAbarprime
            ) / rfm.ReDD[i][j]
    return hprimeDD_enforced, aprimeDD_enforced


def _project_orthogonal_rfm(
    Bq: BSSNQuantities,
) -> Tuple[List[List[sp.Expr]], List[List[sp.Expr]]]:
    """
    Project BSSN variables for an orthogonal reference metric.

    :param Bq: BSSN quantities to project.
    :return: Projected rescaled conformal metric and extrinsic-curvature tensors.
    """
    hprimeDD_enforced = ixp.zerorank2()
    aprimeDD_enforced = ixp.zerorank2()
    # Every standard RFM is orthogonal. With diagonal scale-factor matrix D,
    # gammabar = D H D, gammahat = D I D, and Abar = D a D, where H = I+h.
    # Thus det(gammabar)/det(gammahat) = det(H), so q = det(H)^(-1/3)
    # and h' = q H - I. Also gammabar'^(-1) = q^(-1) D^(-1) H^(-1)
    # D^(-1), hence tr(gammabar'^(-1) Abar) = q^(-1) (H^(-1):a).
    # The q factors cancel in gammabar' tr/3, leaving a' = a-H(H^(-1):a)/3.
    HDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            HDD[i][j] = Bq.hDD[i][j] + sp.KroneckerDelta(i, j)
    HUU, detH = ixp.symm_matrix_inverter3x3(HDD)
    q = detH ** sp.Rational(-1, 3)
    trA = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            trA += HUU[i][j] * Bq.aDD[i][j]
    for i in range(3):
        for j in range(3):
            hprimeDD_enforced[i][j] = q * HDD[i][j] - sp.KroneckerDelta(i, j)
            aprimeDD_enforced[i][j] = Bq.aDD[i][j] - sp.Rational(1, 3) * HDD[i][j] * trA
    return hprimeDD_enforced, aprimeDD_enforced


def BSSN_algebraic_constraints(
    CoordSystem: str, enable_rfm_precompute: bool = False
) -> Tuple[List[List[sp.Expr]], List[List[sp.Expr]]]:
    """
    Project BSSN variables onto their determinant and trace constraints.

    :param CoordSystem: Coordinate system used by BSSN quantities.
    :param enable_rfm_precompute: Whether to enable reference-metric precomputation.
    :return: Projected rescaled conformal metric and extrinsic-curvature tensors.
    """
    cache_key = CoordSystem + ("_rfm_precompute" if enable_rfm_precompute else "")
    Bq = BSSN_quantities[cache_key]
    if CoordSystem.startswith("GeneralRFM"):
        rfm = refmetric.reference_metric[cache_key]
        return _project_general_rfm(Bq, rfm)
    return _project_orthogonal_rfm(Bq)


if __name__ == "__main__":
    import os

    import nrpy.validate_expressions.validate_expressions as ve

    for validation_coord, validation_precompute in (
        ("Cartesian", False),
        ("SinhSpherical", True),
        ("GeneralRFM", True),
    ):
        hprimeDD, aprimeDD = BSSN_algebraic_constraints(
            validation_coord, validation_precompute
        )
        hprimeDD = [
            [expr.subs(sp.Function("nrpyAbs"), sp.Abs) for expr in row]
            for row in hprimeDD
        ]
        aprimeDD = [
            [expr.subs(sp.Function("nrpyAbs"), sp.Abs) for expr in row]
            for row in aprimeDD
        ]
        results_dict = ve.process_dictionary_of_expressions(
            {"hprimeDD": hprimeDD, "aprimeDD": aprimeDD},
            fixed_mpfs_for_free_symbols=True,
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{validation_coord}"
            + ("_rfm_precompute" if validation_precompute else ""),
            results_dict,
        )
