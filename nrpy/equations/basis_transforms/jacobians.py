# nrpy/equations/basis_transforms/jacobians.py
"""
Core basis-transformation functions for reference-metric, Cartesian, and Spherical bases.

These routines use Jacobian and inverse-Jacobian matrices that are computed and
stored on each ``ReferenceMetric`` instance in ``nrpy/reference_metric.py``:

``Jac_dUCart_dDrfmUD``, ``Jac_dUrfm_dDCartUD``, ``Jac_dUSph_dDrfmUD``,
and ``Jac_dUrfm_dDSphUD``.

This module applies those precomputed matrices to vectors/tensors; it does not
construct Jacobians itself.

Author: Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

from typing import Dict, List

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.reference_metric as refmetric


class BasisTransforms:
    """
    Construct and apply basis transformations for a specific coordinate system.

    :param CoordSystem: Coordinate system used to fetch the corresponding ReferenceMetric.
    """

    def __init__(self, CoordSystem: str) -> None:
        self.CoordSystem = CoordSystem
        self.rfm = refmetric.reference_metric[CoordSystem]

    def basis_transform_vectorU_from_rfmbasis_to_Cartesian(
        self, src_vectorU: List[sp.Expr]
    ) -> List[sp.Expr]:
        """
        Transform a vector from the reference metric basis to the Cartesian basis.

        :param src_vectorU: Source vector in the reference metric basis.
        :return: Transformed vector in Cartesian basis.
        """
        Cart_dst_vectorU = ixp.zerorank1()
        for i in range(3):
            for l in range(3):
                Cart_dst_vectorU[i] += (
                    self.rfm.Jac_dUCart_dDrfmUD[i][l] * src_vectorU[l]
                )
        return Cart_dst_vectorU

    def basis_transform_vectorU_from_rfmbasis_to_Spherical(
        self, src_vectorU: List[sp.Expr]
    ) -> List[sp.Expr]:
        """
        Transform a vector from the reference metric basis to the Spherical basis.

        :param src_vectorU: Source vector in the reference metric basis.
        :return: Transformed vector in Spherical basis.
        """
        Sph_dst_vectorU = ixp.zerorank1()
        for i in range(3):
            for l in range(3):
                Sph_dst_vectorU[i] += self.rfm.Jac_dUSph_dDrfmUD[i][l] * src_vectorU[l]
        return Sph_dst_vectorU

    def basis_transform_vectorD_from_rfmbasis_to_Cartesian(
        self, src_vectorD: List[sp.Expr]
    ) -> List[sp.Expr]:
        """
        Transform a covariant vector from the reference metric basis to the Cartesian basis.

        :param src_vectorD: Source covariant vector in the reference metric basis.
        :return: Transformed covariant vector in Cartesian basis.
        """
        Cart_dst_vectorD = ixp.zerorank1()
        for i in range(3):
            for l in range(3):
                Cart_dst_vectorD[i] += (
                    self.rfm.Jac_dUrfm_dDCartUD[l][i] * src_vectorD[l]
                )
        return Cart_dst_vectorD

    def basis_transform_vectorD_from_rfmbasis_to_Spherical(
        self, src_vectorD: List[sp.Expr]
    ) -> List[sp.Expr]:
        """
        Transform a covariant vector from the reference metric basis to the Spherical basis.

        :param src_vectorD: Source covariant vector in the reference metric basis.
        :return: Transformed covariant vector in Spherical basis.
        """
        Sph_dst_vectorD = ixp.zerorank1()
        for i in range(3):
            for l in range(3):
                Sph_dst_vectorD[i] += self.rfm.Jac_dUrfm_dDSphUD[l][i] * src_vectorD[l]
        return Sph_dst_vectorD

    def basis_transform_tensorDD_from_rfmbasis_to_Cartesian(
        self, src_tensorDD: List[List[sp.Expr]]
    ) -> List[List[sp.Expr]]:
        """
        Transform a rank-2 tensor from the reference metric basis to the Cartesian basis.

        :param src_tensorDD: Source tensor in the reference metric basis.
        :return: Transformed tensor in Cartesian basis.
        """
        Cart_dst_tensorDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                for l in range(3):
                    for m in range(3):
                        Cart_dst_tensorDD[i][j] += (
                            self.rfm.Jac_dUrfm_dDCartUD[l][i]
                            * self.rfm.Jac_dUrfm_dDCartUD[m][j]
                            * src_tensorDD[l][m]
                        )
        return Cart_dst_tensorDD

    def basis_transform_tensorDD_from_rfmbasis_to_Spherical(
        self, src_tensorDD: List[List[sp.Expr]]
    ) -> List[List[sp.Expr]]:
        """
        Transform a rank-2 tensor from the reference metric basis to the Spherical basis.

        :param src_tensorDD: Source tensor in the reference metric basis.
        :return: Transformed tensor in Spherical basis.
        """
        Sph_dst_tensorDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                for l in range(3):
                    for m in range(3):
                        Sph_dst_tensorDD[i][j] += (
                            self.rfm.Jac_dUrfm_dDSphUD[l][i]
                            * self.rfm.Jac_dUrfm_dDSphUD[m][j]
                            * src_tensorDD[l][m]
                        )
        return Sph_dst_tensorDD

    def basis_transform_vectorU_from_Cartesian_to_rfmbasis(
        self, Cart_src_vectorU: List[sp.Expr]
    ) -> List[sp.Expr]:
        """
        Transform a vector from the Cartesian basis to the reference metric basis.

        :param Cart_src_vectorU: Source vector in Cartesian basis.
        :return: Transformed vector in reference metric basis.
        """
        rfm_dst_vectorU = ixp.zerorank1()
        for i in range(3):
            for l in range(3):
                rfm_dst_vectorU[i] += (
                    self.rfm.Jac_dUrfm_dDCartUD[i][l] * Cart_src_vectorU[l]
                )
        return rfm_dst_vectorU

    def basis_transform_vectorU_from_Spherical_to_rfmbasis(
        self, Sph_src_vectorU: List[sp.Expr]
    ) -> List[sp.Expr]:
        """
        Transform a vector from the Spherical basis to the reference metric basis.

        :param Sph_src_vectorU: Source vector in Spherical basis.
        :return: Transformed vector in reference metric basis.
        """
        rfm_dst_vectorU = ixp.zerorank1()
        for i in range(3):
            for l in range(3):
                rfm_dst_vectorU[i] += (
                    self.rfm.Jac_dUrfm_dDSphUD[i][l] * Sph_src_vectorU[l]
                )
        return rfm_dst_vectorU

    def basis_transform_vectorD_from_Cartesian_to_rfmbasis(
        self, Cart_src_vectorD: List[sp.Expr]
    ) -> List[sp.Expr]:
        """
        Transform a covariant vector from the Cartesian basis to the reference metric basis.

        :param Cart_src_vectorD: Source covariant vector in Cartesian basis.
        :return: Transformed covariant vector in reference metric basis.
        """
        rfm_dst_vectorD = ixp.zerorank1()
        for i in range(3):
            for l in range(3):
                rfm_dst_vectorD[i] += (
                    self.rfm.Jac_dUCart_dDrfmUD[l][i] * Cart_src_vectorD[l]
                )
        return rfm_dst_vectorD

    def basis_transform_vectorD_from_Spherical_to_rfmbasis(
        self, Sph_src_vectorD: List[sp.Expr]
    ) -> List[sp.Expr]:
        """
        Transform a covariant vector from the Spherical basis to the reference metric basis.

        :param Sph_src_vectorD: Source covariant vector in Spherical basis.
        :return: Transformed covariant vector in reference metric basis.
        """
        rfm_dst_vectorD = ixp.zerorank1()
        for i in range(3):
            for l in range(3):
                rfm_dst_vectorD[i] += (
                    self.rfm.Jac_dUSph_dDrfmUD[l][i] * Sph_src_vectorD[l]
                )
        return rfm_dst_vectorD

    def basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
        self, Cart_src_tensorDD: List[List[sp.Expr]]
    ) -> List[List[sp.Expr]]:
        """
        Transform a rank-2 tensor from the Cartesian basis to the reference metric basis.

        :param Cart_src_tensorDD: Source tensor in Cartesian basis.
        :return: Transformed tensor in reference metric basis.
        """
        rfm_dst_tensorDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                for l in range(3):
                    for m in range(3):
                        rfm_dst_tensorDD[i][j] += (
                            self.rfm.Jac_dUCart_dDrfmUD[l][i]
                            * self.rfm.Jac_dUCart_dDrfmUD[m][j]
                            * Cart_src_tensorDD[l][m]
                        )
        return rfm_dst_tensorDD

    def basis_transform_tensorDD_from_Spherical_to_rfmbasis(
        self, Sph_src_tensorDD: List[List[sp.Expr]]
    ) -> List[List[sp.Expr]]:
        """
        Transform a rank-2 tensor from the Spherical basis to the reference metric basis.

        :param Sph_src_tensorDD: Source tensor in Spherical basis.
        :return: Transformed tensor in reference metric basis.
        """
        rfm_dst_tensorDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                for l in range(3):
                    for m in range(3):
                        rfm_dst_tensorDD[i][j] += (
                            self.rfm.Jac_dUSph_dDrfmUD[l][i]
                            * self.rfm.Jac_dUSph_dDrfmUD[m][j]
                            * Sph_src_tensorDD[l][m]
                        )
        return rfm_dst_tensorDD

    # Useful for stress-energy tensor; assumes reference metric is time-independent.
    def basis_transform_4tensorUU_from_time_indep_rfmbasis_to_Cartesian(
        self, T4UU: List[List[sp.Expr]]
    ) -> List[List[sp.Expr]]:
        """
        Transform a 4-tensor from time-independent reference metric basis to Cartesian basis.

        :param T4UU: Source 4-tensor in reference metric basis.
        :return: Transformed 4-tensor in Cartesian basis.
        """
        Jac4_dUCart_dDrfmUD = ixp.zerorank2(dimension=4)
        Jac4_dUCart_dDrfmUD[0][0] = sp.sympify(1)
        for i in range(3):
            for j in range(3):
                Jac4_dUCart_dDrfmUD[i + 1][j + 1] = self.rfm.Jac_dUCart_dDrfmUD[i][j]

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
        self, T4UU: List[List[sp.Expr]]
    ) -> List[List[sp.Expr]]:
        """
        Transform a 4-tensor from time-independent reference metric basis to Spherical basis.

        :param T4UU: Source 4-tensor in reference metric basis.
        :return: Transformed 4-tensor in Spherical basis.
        """
        Jac4_dUSph_dDrfmUD = ixp.zerorank2(dimension=4)
        Jac4_dUSph_dDrfmUD[0][0] = sp.sympify(1)
        for i in range(3):
            for j in range(3):
                Jac4_dUSph_dDrfmUD[i + 1][j + 1] = self.rfm.Jac_dUSph_dDrfmUD[i][j]

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
        self, T4UU: List[List[sp.Expr]]
    ) -> List[List[sp.Expr]]:
        """
        Transform a 4-tensor from Cartesian basis to time-independent reference metric basis.

        :param T4UU: Source 4-tensor in Cartesian basis.
        :return: Transformed 4-tensor in reference metric basis.
        """
        Jac4_dUrfm_dDCartUD = ixp.zerorank2(dimension=4)
        Jac4_dUrfm_dDCartUD[0][0] = sp.sympify(1)
        for i in range(3):
            for j in range(3):
                Jac4_dUrfm_dDCartUD[i + 1][j + 1] = self.rfm.Jac_dUrfm_dDCartUD[i][j]

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
        self, T4UU: List[List[sp.Expr]]
    ) -> List[List[sp.Expr]]:
        """
        Transform a 4-tensor from Spherical basis to time-independent reference metric basis.

        :param T4UU: Source 4-tensor in Spherical basis.
        :return: Transformed 4-tensor in reference metric basis.
        """
        Jac4_dUrfm_dDSphUD = ixp.zerorank2(dimension=4)
        Jac4_dUrfm_dDSphUD[0][0] = sp.sympify(1)
        for i in range(3):
            for j in range(3):
                Jac4_dUrfm_dDSphUD[i + 1][j + 1] = self.rfm.Jac_dUrfm_dDSphUD[i][j]

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


class BasisTransforms_dict(Dict[str, BasisTransforms]):
    """Custom dictionary for storing BasisTransforms objects."""

    def __getitem__(self, CoordSystem_in: str) -> BasisTransforms:
        if CoordSystem_in not in self:
            print(f"Setting up basis_transforms[{CoordSystem_in}]...")
            self.__setitem__(CoordSystem_in, BasisTransforms(CoordSystem_in))
        return dict.__getitem__(self, CoordSystem_in)

    def __setitem__(self, CoordSystem: str, value: BasisTransforms) -> None:
        dict.__setitem__(self, CoordSystem, value)

    def __delitem__(self, CoordSystem: str) -> None:
        dict.__delitem__(self, CoordSystem)


basis_transforms = BasisTransforms_dict()


if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.validate_expressions.validate_expressions as ve

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    src_vectorU_rfm = ixp.declarerank1("src_vectorU_rfm")
    src_vectorD_rfm = ixp.declarerank1("src_vectorD_rfm")
    src_tensorDD_rfm = ixp.declarerank2("src_tensorDD_rfm")

    src_vectorU_Cart = ixp.declarerank1("src_vectorU_Cart")
    src_vectorD_Cart = ixp.declarerank1("src_vectorD_Cart")
    src_tensorDD_Cart = ixp.declarerank2("src_tensorDD_Cart")

    src_vectorU_Sph = ixp.declarerank1("src_vectorU_Sph")
    src_vectorD_Sph = ixp.declarerank1("src_vectorD_Sph")
    src_tensorDD_Sph = ixp.declarerank2("src_tensorDD_Sph")

    src_T4UU_rfm = ixp.declarerank2("src_T4UU_rfm", dimension=4)
    src_T4UU_Cart = ixp.declarerank2("src_T4UU_Cart", dimension=4)
    src_T4UU_Sph = ixp.declarerank2("src_T4UU_Sph", dimension=4)

    for Coord in refmetric.supported_CoordSystems:
        bt = basis_transforms[Coord]
        validation_dict = {
            "vectorU_rfm_to_Cart": bt.basis_transform_vectorU_from_rfmbasis_to_Cartesian(
                src_vectorU_rfm
            ),
            "vectorU_rfm_to_Spherical": bt.basis_transform_vectorU_from_rfmbasis_to_Spherical(
                src_vectorU_rfm
            ),
            "vectorD_rfm_to_Cart": bt.basis_transform_vectorD_from_rfmbasis_to_Cartesian(
                src_vectorD_rfm
            ),
            "vectorD_rfm_to_Spherical": bt.basis_transform_vectorD_from_rfmbasis_to_Spherical(
                src_vectorD_rfm
            ),
            "tensorDD_rfm_to_Cart": bt.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(
                src_tensorDD_rfm
            ),
            "tensorDD_rfm_to_Spherical": bt.basis_transform_tensorDD_from_rfmbasis_to_Spherical(
                src_tensorDD_rfm
            ),
            "vectorU_Cart_to_rfm": bt.basis_transform_vectorU_from_Cartesian_to_rfmbasis(
                src_vectorU_Cart
            ),
            "vectorU_Spherical_to_rfm": bt.basis_transform_vectorU_from_Spherical_to_rfmbasis(
                src_vectorU_Sph
            ),
            "vectorD_Cart_to_rfm": bt.basis_transform_vectorD_from_Cartesian_to_rfmbasis(
                src_vectorD_Cart
            ),
            "vectorD_Spherical_to_rfm": bt.basis_transform_vectorD_from_Spherical_to_rfmbasis(
                src_vectorD_Sph
            ),
            "tensorDD_Cart_to_rfm": bt.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
                src_tensorDD_Cart
            ),
            "tensorDD_Spherical_to_rfm": bt.basis_transform_tensorDD_from_Spherical_to_rfmbasis(
                src_tensorDD_Sph
            ),
            "T4UU_rfm_to_Cart": bt.basis_transform_4tensorUU_from_time_indep_rfmbasis_to_Cartesian(
                src_T4UU_rfm
            ),
            "T4UU_rfm_to_Spherical": bt.basis_transform_4tensorUU_from_time_indep_rfmbasis_to_Spherical(
                src_T4UU_rfm
            ),
            "T4UU_Cart_to_rfm": bt.basis_transform_4tensorUU_from_Cartesian_to_time_indep_rfmbasis(
                src_T4UU_Cart
            ),
            "T4UU_Spherical_to_rfm": bt.basis_transform_4tensorUU_from_Spherical_to_time_indep_rfmbasis(
                src_T4UU_Sph
            ),
        }

        results_dict = ve.process_dictionary_of_expressions(
            validation_dict, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{Coord}",
            results_dict,
        )
