"""
Core Jacobian (basis) transformation functions.
For reference metric basis to/from the Cartesian and Spherical bases.

We define Jacobians relative to the reference metric
basis at a point x^j_rfm=(xx0,xx1,xx2)_rfm on the source grid:

Jac_dUCart_dDrfmUD[i][j] = dx^i_Cart / dx^j_rfm

via exact differentiation (courtesy SymPy), and the inverse Jacobian

Jac_dUrfm_dDCartUD[i][j] = dx^i_rfm / dx^j_Cart

using NRPy's generic_matrix_inverter3x3() function.

Author: Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

from typing import Sequence

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.reference_metric as refmetric


def basis_transform_vectorU_from_rfmbasis_to_Cartesian(
    CoordSystem: str, src_vectorU: Sequence[sp.Expr]
) -> Sequence[sp.Expr]:
    """
    Transform a vector from the reference metric basis to the Cartesian basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param src_vectorU: The source vector in the reference metric basis.
    :return: Transformed vector in Cartesian basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    Cart_dst_vectorU = ixp.zerorank1()
    for i in range(3):
        for l in range(3):
            Cart_dst_vectorU[i] += rfm.Jac_dUCart_dDrfmUD[i][l] * src_vectorU[l]
    return Cart_dst_vectorU


def basis_transform_vectorU_from_rfmbasis_to_Spherical(
    CoordSystem: str, src_vectorU: Sequence[sp.Expr]
) -> Sequence[sp.Expr]:
    """
    Transform a vector from the reference metric basis to the Spherical basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param src_vectorU: The source vector in the reference metric basis.
    :return: Transformed vector in Spherical basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    Sph_dst_vectorU = ixp.zerorank1()
    for i in range(3):
        for l in range(3):
            Sph_dst_vectorU[i] += rfm.Jac_dUSph_dDrfmUD[i][l] * src_vectorU[l]
    return Sph_dst_vectorU


def basis_transform_vectorD_from_rfmbasis_to_Cartesian(
    CoordSystem: str, src_vectorD: Sequence[sp.Expr]
) -> Sequence[sp.Expr]:
    """
    Transform a covariant vector from the reference metric basis to the Cartesian basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param src_vectorD: The source covariant vector in the reference metric basis.
    :return: Transformed covariant vector in Cartesian basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    Cart_dst_vectorD = ixp.zerorank1()
    for i in range(3):
        for l in range(3):
            Cart_dst_vectorD[i] += rfm.Jac_dUrfm_dDCartUD[l][i] * src_vectorD[l]
    return Cart_dst_vectorD


def basis_transform_vectorD_from_rfmbasis_to_Spherical(
    CoordSystem: str, src_vectorD: Sequence[sp.Expr]
) -> Sequence[sp.Expr]:
    """
    Transform a covariant vector from the reference metric basis to the Spherical basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param src_vectorD: The source covariant vector in the reference metric basis.
    :return: Transformed covariant vector in Spherical basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    Sph_dst_vectorD = ixp.zerorank1()
    for i in range(3):
        for l in range(3):
            Sph_dst_vectorD[i] += rfm.Jac_dUrfm_dDSphUD[l][i] * src_vectorD[l]
    return Sph_dst_vectorD


def basis_transform_tensorDD_from_rfmbasis_to_Cartesian(
    CoordSystem: str, src_tensorDD: Sequence[Sequence[sp.Expr]]
) -> Sequence[Sequence[sp.Expr]]:
    """
    Transform a rank-2 tensor from the reference metric basis to the Cartesian basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param src_tensorDD: The source tensor in the reference metric basis.
    :return: Transformed tensor in Cartesian basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    Cart_dst_tensorDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            for l in range(3):
                for m in range(3):
                    Cart_dst_tensorDD[i][j] += (
                        rfm.Jac_dUrfm_dDCartUD[l][i]
                        * rfm.Jac_dUrfm_dDCartUD[m][j]
                        * src_tensorDD[l][m]
                    )
    return Cart_dst_tensorDD


def basis_transform_tensorDD_from_rfmbasis_to_Spherical(
    CoordSystem: str, src_tensorDD: Sequence[Sequence[sp.Expr]]
) -> Sequence[Sequence[sp.Expr]]:
    """
    Transform a rank-2 tensor from the reference metric basis to the Spherical basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param src_tensorDD: The source tensor in the reference metric basis.
    :return: Transformed tensor in Spherical basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    Sph_dst_tensorDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            for l in range(3):
                for m in range(3):
                    Sph_dst_tensorDD[i][j] += (
                        rfm.Jac_dUrfm_dDSphUD[l][i]
                        * rfm.Jac_dUrfm_dDSphUD[m][j]
                        * src_tensorDD[l][m]
                    )
    return Sph_dst_tensorDD


def basis_transform_vectorU_from_Cartesian_to_rfmbasis(
    CoordSystem: str, Cart_src_vectorU: Sequence[sp.Expr]
) -> Sequence[sp.Expr]:
    """
    Transform a vector from the Cartesian basis to the reference metric basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param Cart_src_vectorU: The source vector in Cartesian basis.
    :return: Transformed vector in reference metric basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    rfm_dst_vectorU = ixp.zerorank1()
    for i in range(3):
        for l in range(3):
            rfm_dst_vectorU[i] += rfm.Jac_dUrfm_dDCartUD[i][l] * Cart_src_vectorU[l]
    return rfm_dst_vectorU


def basis_transform_vectorU_from_Spherical_to_rfmbasis(
    CoordSystem: str, Sph_src_vectorU: Sequence[sp.Expr]
) -> Sequence[sp.Expr]:
    """
    Transform a vector from the Spherical basis to the reference metric basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param Sph_src_vectorU: The source vector in Spherical basis.
    :return: Transformed vector in reference metric basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    rfm_dst_vectorU = ixp.zerorank1()
    for i in range(3):
        for l in range(3):
            rfm_dst_vectorU[i] += rfm.Jac_dUrfm_dDSphUD[i][l] * Sph_src_vectorU[l]
    return rfm_dst_vectorU


def basis_transform_vectorD_from_Cartesian_to_rfmbasis(
    CoordSystem: str, Cart_src_vectorD: Sequence[sp.Expr]
) -> Sequence[sp.Expr]:
    """
    Transform a covariant vector from the Cartesian basis to the reference metric basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param Cart_src_vectorD: The source covariant vector in Cartesian basis.
    :return: Transformed covariant vector in reference metric basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    rfm_dst_vectorD = ixp.zerorank1()
    for i in range(3):
        for l in range(3):
            rfm_dst_vectorD[i] += rfm.Jac_dUCart_dDrfmUD[l][i] * Cart_src_vectorD[l]
    return rfm_dst_vectorD


def basis_transform_vectorD_from_Spherical_to_rfmbasis(
    CoordSystem: str, Sph_src_vectorD: Sequence[sp.Expr]
) -> Sequence[sp.Expr]:
    """
    Transform a covariant vector from the Spherical basis to the reference metric basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param Sph_src_vectorD: The source covariant vector in Spherical basis.
    :return: Transformed covariant vector in reference metric basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    rfm_dst_vectorD = ixp.zerorank1()
    for i in range(3):
        for l in range(3):
            rfm_dst_vectorD[i] += rfm.Jac_dUSph_dDrfmUD[l][i] * Sph_src_vectorD[l]
    return rfm_dst_vectorD


def basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
    CoordSystem: str, Cart_src_tensorDD: Sequence[Sequence[sp.Expr]]
) -> Sequence[Sequence[sp.Expr]]:
    """
    Transform a rank-2 tensor from the Cartesian basis to the reference metric basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param Cart_src_tensorDD: The source tensor in the Cartesian basis.
    :return: Transformed tensor in reference metric basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    rfm_dst_tensorDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            for l in range(3):
                for m in range(3):
                    rfm_dst_tensorDD[i][j] += (
                        rfm.Jac_dUCart_dDrfmUD[l][i]
                        * rfm.Jac_dUCart_dDrfmUD[m][j]
                        * Cart_src_tensorDD[l][m]
                    )
    return rfm_dst_tensorDD


def basis_transform_tensorDD_from_Spherical_to_rfmbasis(
    CoordSystem: str, Sph_src_tensorDD: Sequence[Sequence[sp.Expr]]
) -> Sequence[Sequence[sp.Expr]]:
    """
    Transform a rank-2 tensor from the Spherical basis to the reference metric basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param Sph_src_tensorDD: The source tensor in the Spherical basis.
    :return: Transformed tensor in reference metric basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    rfm_dst_tensorDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            for l in range(3):
                for m in range(3):
                    rfm_dst_tensorDD[i][j] += (
                        rfm.Jac_dUSph_dDrfmUD[l][i]
                        * rfm.Jac_dUSph_dDrfmUD[m][j]
                        * Sph_src_tensorDD[l][m]
                    )
    return rfm_dst_tensorDD


# Useful for stress-energy tensor; assumes reference metric is time-independent.
def basis_transform_4tensorUU_from_time_indep_rfmbasis_to_Cartesian(
    CoordSystem: str, T4UU: Sequence[Sequence[sp.Expr]]
) -> Sequence[Sequence[sp.Expr]]:
    """
    Transform a 4-tensor from time-independent reference metric basis to Cartesian basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param T4UU: The source 4-tensor in the reference metric basis.
    :return: Transformed 4-tensor in Cartesian basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]

    Jac4_dUCart_dDrfmUD = ixp.zerorank2(dimension=4)
    Jac4_dUCart_dDrfmUD[0][0] = sp.sympify(1)
    for i in range(3):
        for j in range(3):
            Jac4_dUCart_dDrfmUD[i + 1][j + 1] = rfm.Jac_dUCart_dDrfmUD[i][j]

    Cart_dst_T4UU = ixp.zerorank2(dimension=4)
    for mu in range(4):
        for nu in range(4):
            for delta in range(4):
                for sigma in range(4):
                    Cart_dst_T4UU[mu][nu] += (
                        Jac4_dUCart_dDrfmUD[mu][delta]
                        * Jac4_dUCart_dDrfmUD[nu][sigma]
                        * T4UU[delta][sigma]
                    )
    return Cart_dst_T4UU


# Useful for stress-energy tensor; assumes reference metric is time-independent.
def basis_transform_4tensorUU_from_time_indep_rfmbasis_to_Spherical(
    CoordSystem: str, T4UU: Sequence[Sequence[sp.Expr]]
) -> Sequence[Sequence[sp.Expr]]:
    """
    Transform a 4-tensor from time-independent reference metric basis to Spherical basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param T4UU: The source 4-tensor in the reference metric basis.
    :return: Transformed 4-tensor in Spherical basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]

    Jac4_dUSph_dDrfmUD = ixp.zerorank2(dimension=4)
    Jac4_dUSph_dDrfmUD[0][0] = sp.sympify(1)
    for i in range(3):
        for j in range(3):
            Jac4_dUSph_dDrfmUD[i + 1][j + 1] = rfm.Jac_dUSph_dDrfmUD[i][j]

    Sph_dst_T4UU = ixp.zerorank2(dimension=4)
    for mu in range(4):
        for nu in range(4):
            for delta in range(4):
                for sigma in range(4):
                    Sph_dst_T4UU[mu][nu] += (
                        Jac4_dUSph_dDrfmUD[mu][delta]
                        * Jac4_dUSph_dDrfmUD[nu][sigma]
                        * T4UU[delta][sigma]
                    )
    return Sph_dst_T4UU


def basis_transform_4tensorUU_from_Cartesian_to_time_indep_rfmbasis(
    CoordSystem: str, T4UU: Sequence[Sequence[sp.Expr]]
) -> Sequence[Sequence[sp.Expr]]:
    """
    Transform a 4-tensor from Cartesian basis to time-independent reference metric basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param T4UU: The source 4-tensor in Cartesian basis.
    :return: Transformed 4-tensor in reference metric basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]

    Jac4_dUrfm_dDCartUD = ixp.zerorank2(dimension=4)
    Jac4_dUrfm_dDCartUD[0][0] = sp.sympify(1)
    for i in range(3):
        for j in range(3):
            Jac4_dUrfm_dDCartUD[i + 1][j + 1] = rfm.Jac_dUrfm_dDCartUD[i][j]

    rfm_dst_T4UU = ixp.zerorank2(dimension=4)
    for mu in range(4):
        for nu in range(4):
            for delta in range(4):
                for sigma in range(4):
                    rfm_dst_T4UU[mu][nu] += (
                        Jac4_dUrfm_dDCartUD[mu][delta]
                        * Jac4_dUrfm_dDCartUD[nu][sigma]
                        * T4UU[delta][sigma]
                    )
    return rfm_dst_T4UU


def basis_transform_4tensorUU_from_Spherical_to_time_indep_rfmbasis(
    CoordSystem: str, T4UU: Sequence[Sequence[sp.Expr]]
) -> Sequence[Sequence[sp.Expr]]:
    """
    Transform a 4-tensor from Spherical basis to time-independent reference metric basis.

    :param CoordSystem: Coordinate system to use for the reference metric.
    :param T4UU: The source 4-tensor in Spherical basis.
    :return: Transformed 4-tensor in reference metric basis.
    """
    rfm = refmetric.reference_metric[CoordSystem]

    Jac4_dUrfm_dDSphUD = ixp.zerorank2(dimension=4)
    Jac4_dUrfm_dDSphUD[0][0] = sp.sympify(1)
    for i in range(3):
        for j in range(3):
            Jac4_dUrfm_dDSphUD[i + 1][j + 1] = rfm.Jac_dUrfm_dDSphUD[i][j]

    rfm_dst_T4UU = ixp.zerorank2(dimension=4)
    for mu in range(4):
        for nu in range(4):
            for delta in range(4):
                for sigma in range(4):
                    rfm_dst_T4UU[mu][nu] += (
                        Jac4_dUrfm_dDSphUD[mu][delta]
                        * Jac4_dUrfm_dDSphUD[nu][sigma]
                        * T4UU[delta][sigma]
                    )
    return rfm_dst_T4UU


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
