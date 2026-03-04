"""
Construct symbolic SO(3) matrix-rotation equations.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import sympy as sp

import nrpy.indexedexp as ixp


class SO3Expressions:
    """
    Construct and store symbolic SO(3) matrix-rotation expressions.

    This class provides a class-backed representation of the rotation equations
    used by BHaH rotation infrastructure.
    """

    def __init__(self) -> None:
        """Set up all core SO(3) symbolic expressions and store them as attributes."""
        # Step 1: Declare core symbols for cumulative hats and rotation matrix expressions.
        self.xhatU: List[sp.Expr] = ixp.declarerank1("xhatU")
        self.yhatU: List[sp.Expr] = ixp.declarerank1("yhatU")
        self.zhatU: List[sp.Expr] = ixp.declarerank1("zhatU")
        self.R_from_hats: List[List[sp.Expr]] = (
            SO3Expressions.build_rotation_matrix_from_hats(
                self.xhatU, self.yhatU, self.zhatU
            )
        )

        # Step 2: Declare generic symbols for vector/tensor rotations.
        self.R: List[List[sp.Expr]] = ixp.declarerank2("R", symmetry="nosym")
        self.R_dst: List[List[sp.Expr]] = ixp.declarerank2("R_dst", symmetry="nosym")
        self.R_src: List[List[sp.Expr]] = ixp.declarerank2("R_src", symmetry="nosym")
        self.vU: List[sp.Expr] = ixp.declarerank1("vU")
        self.tDD: List[List[sp.Expr]] = ixp.declarerank2("tDD", symmetry="nosym")

        self.vU_rotated: List[sp.Expr] = SO3Expressions.apply_R_to_vector(
            self.R, self.vU
        )
        self.vU_unrotated: List[sp.Expr] = SO3Expressions.apply_RT_to_vector(
            self.R, self.vU
        )
        self.tDD_rotated: List[List[sp.Expr]] = SO3Expressions.apply_R_to_tensorDD(
            self.R, self.tDD
        )
        self.tDD_unrotated: List[List[sp.Expr]] = SO3Expressions.apply_RT_to_tensorDD(
            self.R, self.tDD
        )
        self.deltaR_dst_from_src: List[List[sp.Expr]] = (
            SO3Expressions.relative_rotation_dst_from_src(self.R_dst, self.R_src)
        )

        # Step 3: Declare symbols for scalar/vector algebra helpers.
        self.aU: List[sp.Expr] = ixp.declarerank1("aU")
        self.bU: List[sp.Expr] = ixp.declarerank1("bU")
        self.dot_a_b: sp.Expr = SO3Expressions.dot_product3(self.aU, self.bU)
        self.norm_a: sp.Expr = SO3Expressions.norm3(self.aU)
        self.cross_a_b: List[sp.Expr] = SO3Expressions.cross_product3(self.aU, self.bU)

        self.det_from_hats: sp.Expr = SO3Expressions.det_from_columns(
            self.xhatU, self.yhatU, self.zhatU
        )
        self.hat_invariants: Dict[str, sp.Expr] = (
            SO3Expressions.hat_validation_invariants(self.xhatU, self.yhatU, self.zhatU)
        )

        # Step 4: Declare symbols for axis-angle/Rodrigues helpers.
        self.nU: List[sp.Expr] = ixp.declarerank1("nU")
        self.dphi: sp.Expr = sp.Symbol("dphi")
        self.nU_unit: List[sp.Expr]
        self.n_norm: sp.Expr
        self.nU_unit, self.n_norm = SO3Expressions.normalize_axis(self.nU)
        self.R_rodrigues_unit_axis: List[List[sp.Expr]] = (
            SO3Expressions.rodrigues_matrix_from_unit_axis(self.nU_unit, self.dphi)
        )
        self.R_rodrigues_axis_angle: List[List[sp.Expr]] = (
            SO3Expressions.rodrigues_matrix_from_axis_angle(self.nU, self.dphi)
        )

        self.trace_R: sp.Expr = SO3Expressions.matrix_trace3(self.R)
        self.denom: sp.Expr = sp.Symbol("denom")
        self.axis_general_branch: List[sp.Expr] = (
            SO3Expressions.axis_angle_general_branch_axis(self.R, self.denom)
        )
        self.pi_branch_diagonal_terms: List[sp.Expr] = (
            SO3Expressions.axis_angle_pi_branch_diagonal_terms(self.R)
        )

    @staticmethod
    def build_rotation_matrix_from_hats(
        xhatU: List[sp.Expr], yhatU: List[sp.Expr], zhatU: List[sp.Expr]
    ) -> List[List[sp.Expr]]:
        """
        Build R from cumulative hats, with columns R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat.

        :param xhatU: xhat components in fixed basis.
        :param yhatU: yhat components in fixed basis.
        :param zhatU: zhat components in fixed basis.
        :return: Rotation matrix mapping rotating-frame components to fixed-frame components.
        """
        R = ixp.zerorank2()
        for i in range(3):
            R[i][0] = xhatU[i]
            R[i][1] = yhatU[i]
            R[i][2] = zhatU[i]
        return R

    @staticmethod
    def apply_R_to_vector(R: List[List[sp.Expr]], vU: List[sp.Expr]) -> List[sp.Expr]:
        """
        Compute v_dst = R v_src.

        :param R: Rotation matrix.
        :param vU: Source vector.
        :return: Destination vector.
        """
        v_dstU = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                v_dstU[i] += R[i][j] * vU[j]
        return v_dstU

    @staticmethod
    def apply_RT_to_vector(R: List[List[sp.Expr]], vU: List[sp.Expr]) -> List[sp.Expr]:
        """
        Compute v_dst = R^T v_src.

        :param R: Rotation matrix.
        :param vU: Source vector.
        :return: Destination vector.
        """
        v_dstU = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                v_dstU[i] += R[j][i] * vU[j]
        return v_dstU

    @staticmethod
    def apply_R_to_tensorDD(
        R: List[List[sp.Expr]], tDD: List[List[sp.Expr]]
    ) -> List[List[sp.Expr]]:
        """
        Compute T_dst = R T_src R^T.

        :param R: Rotation matrix.
        :param tDD: Source covariant rank-2 tensor.
        :return: Destination covariant rank-2 tensor.
        """
        t_dstDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        t_dstDD[i][j] += R[i][k] * tDD[k][l] * R[j][l]
        return t_dstDD

    @staticmethod
    def apply_RT_to_tensorDD(
        R: List[List[sp.Expr]], tDD: List[List[sp.Expr]]
    ) -> List[List[sp.Expr]]:
        """
        Compute T_dst = R^T T_src R.

        :param R: Rotation matrix.
        :param tDD: Source covariant rank-2 tensor.
        :return: Destination covariant rank-2 tensor.
        """
        t_dstDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        t_dstDD[i][j] += R[k][i] * tDD[k][l] * R[l][j]
        return t_dstDD

    @staticmethod
    def relative_rotation_dst_from_src(
        R_dst: List[List[sp.Expr]], R_src: List[List[sp.Expr]]
    ) -> List[List[sp.Expr]]:
        """
        Compute DeltaR_dst_from_src = R_dst^T R_src.

        :param R_dst: Absolute destination rotation.
        :param R_src: Absolute source rotation.
        :return: Relative rotation from source rotating basis to destination rotating basis.
        """
        deltaR = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    deltaR[i][j] += R_dst[k][i] * R_src[k][j]
        return deltaR

    @staticmethod
    def dot_product3(aU: List[sp.Expr], bU: List[sp.Expr]) -> sp.Expr:
        """
        Compute Euclidean dot product in 3D.

        :param aU: First vector.
        :param bU: Second vector.
        :return: Dot product a dot b.
        """
        out = sp.sympify(0)
        for i in range(3):
            out += aU[i] * bU[i]
        return out

    @staticmethod
    def norm3(aU: List[sp.Expr]) -> sp.Expr:
        """
        Compute Euclidean norm in 3D.

        :param aU: Vector.
        :return: Norm |a|.
        """
        return sp.sqrt(SO3Expressions.dot_product3(aU, aU))

    @staticmethod
    def cross_product3(aU: List[sp.Expr], bU: List[sp.Expr]) -> List[sp.Expr]:
        """
        Compute Euclidean cross product in 3D.

        :param aU: First vector.
        :param bU: Second vector.
        :return: Cross product a cross b.
        """
        outU = ixp.zerorank1()
        outU[0] = aU[1] * bU[2] - aU[2] * bU[1]
        outU[1] = aU[2] * bU[0] - aU[0] * bU[2]
        outU[2] = aU[0] * bU[1] - aU[1] * bU[0]
        return outU

    @staticmethod
    def det_from_columns(
        xhatU: List[sp.Expr], yhatU: List[sp.Expr], zhatU: List[sp.Expr]
    ) -> sp.Expr:
        """
        Compute determinant of a matrix whose columns are xhat, yhat, zhat.

        :param xhatU: First column.
        :param yhatU: Second column.
        :param zhatU: Third column.
        :return: Determinant.
        """
        cross_yz = SO3Expressions.cross_product3(yhatU, zhatU)
        return SO3Expressions.dot_product3(xhatU, cross_yz)

    @staticmethod
    def hat_validation_invariants(
        xhatU: List[sp.Expr], yhatU: List[sp.Expr], zhatU: List[sp.Expr]
    ) -> Dict[str, sp.Expr]:
        """
        Construct scalar invariants used to validate orthonormal right-handed hats.

        :param xhatU: Candidate xhat.
        :param yhatU: Candidate yhat.
        :param zhatU: Candidate zhat.
        :return: Dictionary of invariant expressions.
        """
        return {
            "x_dot_y": SO3Expressions.dot_product3(xhatU, yhatU),
            "x_dot_z": SO3Expressions.dot_product3(xhatU, zhatU),
            "y_dot_z": SO3Expressions.dot_product3(yhatU, zhatU),
            "x_norm": SO3Expressions.norm3(xhatU),
            "y_norm": SO3Expressions.norm3(yhatU),
            "z_norm": SO3Expressions.norm3(zhatU),
            "detR": SO3Expressions.det_from_columns(xhatU, yhatU, zhatU),
        }

    @staticmethod
    def normalize_axis(nU: List[sp.Expr]) -> Tuple[List[sp.Expr], sp.Expr]:
        """
        Normalize a 3-vector axis.

        :param nU: Axis components.
        :return: Tuple (unit_axis, norm).
        """
        nnorm = SO3Expressions.norm3(nU)
        nU_unit = ixp.zerorank1()
        for i in range(3):
            nU_unit[i] = nU[i] / nnorm
        return nU_unit, nnorm

    @staticmethod
    def rodrigues_matrix_from_unit_axis(
        nU_unit: List[sp.Expr], dphi: sp.Expr
    ) -> List[List[sp.Expr]]:
        """
        Build Rodrigues rotation matrix from a unit axis and angle.

        :param nU_unit: Unit rotation axis.
        :param dphi: Rotation angle.
        :return: Rodrigues rotation matrix.
        """
        nx = nU_unit[0]
        ny = nU_unit[1]
        nz = nU_unit[2]
        c = sp.cos(dphi)
        s = sp.sin(dphi)
        one_minus_c = 1 - c

        R = ixp.zerorank2()
        R[0][0] = c + one_minus_c * nx * nx
        R[0][1] = one_minus_c * nx * ny - s * nz
        R[0][2] = one_minus_c * nx * nz + s * ny
        R[1][0] = one_minus_c * ny * nx + s * nz
        R[1][1] = c + one_minus_c * ny * ny
        R[1][2] = one_minus_c * ny * nz - s * nx
        R[2][0] = one_minus_c * nz * nx - s * ny
        R[2][1] = one_minus_c * nz * ny + s * nx
        R[2][2] = c + one_minus_c * nz * nz
        return R

    @staticmethod
    def rodrigues_matrix_from_axis_angle(
        nU: List[sp.Expr], dphi: sp.Expr
    ) -> List[List[sp.Expr]]:
        """
        Build Rodrigues rotation matrix from axis-angle input.

        :param nU: Rotation axis.
        :param dphi: Rotation angle.
        :return: Rodrigues rotation matrix.
        """
        nU_unit, _ = SO3Expressions.normalize_axis(nU)
        return SO3Expressions.rodrigues_matrix_from_unit_axis(nU_unit, dphi)

    @staticmethod
    def matrix_trace3(R: List[List[sp.Expr]]) -> sp.Expr:
        """
        Compute trace of a 3x3 matrix.

        :param R: Matrix.
        :return: trace(R).
        """
        return R[0][0] + R[1][1] + R[2][2]

    @staticmethod
    def axis_angle_general_branch_axis(
        R: List[List[sp.Expr]], denom: sp.Expr
    ) -> List[sp.Expr]:
        """
        Recover axis components from skew(R)/(2 sin(phi)) in the general branch.

        :param R: Rotation matrix.
        :param denom: Denominator 2 sin(phi).
        :return: Axis vector components.
        """
        nU = ixp.zerorank1()
        nU[0] = (R[2][1] - R[1][2]) / denom
        nU[1] = (R[0][2] - R[2][0]) / denom
        nU[2] = (R[1][0] - R[0][1]) / denom
        return nU

    @staticmethod
    def axis_angle_pi_branch_diagonal_terms(R: List[List[sp.Expr]]) -> List[sp.Expr]:
        """
        Build diagonal-derived axis-magnitude terms used in the pi branch.

        :param R: Rotation matrix.
        :return: Terms [n0, n1, n2] before clamping/square-root safeguards.
        """
        return [
            sp.Rational(1, 2) * (R[0][0] + 1),
            sp.Rational(1, 2) * (R[1][1] + 1),
            sp.Rational(1, 2) * (R[2][2] + 1),
        ]


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

    SO3_eqs = SO3Expressions()
    results_dict = ve.process_dictionary_of_expressions(
        SO3_eqs.__dict__, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        # File basename. If this is set to "trusted_module_test1", then
        #   trusted results_dict will be stored in tests/trusted_module_test1.py
        f"{os.path.splitext(os.path.basename(__file__))[0]}",
        results_dict,
    )
