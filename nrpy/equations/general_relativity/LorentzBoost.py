"""
Infrastructure for performing Lorentz boosts in four dimensions.

Author: Thiago Assumpção
        assumpcaothiago **at** gmail **dot* com

License: BSD 2-Clause
"""

# Import needed modules
from typing import List, cast

import sympy as sp  # For symbolic computations

import nrpy.indexedexp as ixp  # NRPy: Symbolic indexed expression (e.g., tensors, vectors, etc.) support


class LorentzBoost:
    """Sets up Lorentz boost in Cartesian coordinates."""

    def __init__(self, vBoost: List[sp.Expr]) -> None:
        """
        Set up Lorentz boost in Cartesian coordinates.

        :param vBoost: 3-velocity vector in Cartesian coordinates for boost direction.

        .. note::
            Class variables that are set in this function: LorentzMatrix and InverseLorentzMatrix
                Lorentz transformation matrix in Cartesian coordinates. See, e.g.,
                https://en.wikipedia.org/wiki/Lorentz_transformation#Proper_transformations:
        """
        # Compute Lorentz matrix
        self.LorentzMatrix = self._compute_LorentzMatrix(vBoost)

        # Compute inverse Lorentz matrix using the 3-velocity vector with a negative sign
        inverse_vBoost = [-vB for vB in vBoost]
        self.InverseLorentzMatrix = self._compute_LorentzMatrix(inverse_vBoost)

    def _compute_LorentzMatrix(self, vBoost: List[sp.Expr]) -> List[List[sp.Expr]]:
        """
        Create Lorentz boost matrix in Cartesian coordinates.

        :param vBoost: 3-velocity vector in Cartesian coordinates for boost direction.

        :return: A 4x4 Lorentz transformation matrix in Cartesian coordinates.
        """
        # Step 1: Compute Lorentz factor from 3-velocity
        vBoost2 = sp.sympify(0)
        for i in range(3):
            vBoost2 += vBoost[i] ** 2
        gamma_lorentz = 1 / sp.sqrt(1 - vBoost2)

        # Step 2: Declare indexed expression for Lorentz matrix
        LorentzMatrix = ixp.zerorank2(dimension=4)

        # Step 2.a: Compute LorentzMatrix[0][0]
        LorentzMatrix[0][0] = gamma_lorentz

        # Step 2.b: Compute LorentzMatrix[0][i+1]
        for i in range(3):
            LorentzMatrix[0][i + 1] = -gamma_lorentz * vBoost[i]

        # Step 2.c: Compute LorentzMatrix[i + 1][j + 1]
        for i in range(3):
            for j in range(i, 3):
                LorentzMatrix[i + 1][j + 1] = (
                    (gamma_lorentz - 1) * vBoost[i] * vBoost[j] / vBoost2
                )
                if i == j:  # Add 1 to diagonal term
                    LorentzMatrix[i + 1][j + 1] += sp.sympify(1)

        # Step 2.d: Symmetrize matrix
        _LorentzMatrix = ixp.symmetrize_rank2(
            indexedexp=LorentzMatrix, symmetry="sym01", dimension=4
        )
        # Since ixp.symmetrize_rank2 returns a List[List[sp.Expr]] type,
        # we need to recast LorentzMatrix into a List[List[sp.Expr]] again
        LorentzMatrix = cast(List[List[sp.Expr]], _LorentzMatrix)

        return LorentzMatrix

    def boost_vectorU(self, vectorU: List[sp.Expr]) -> List[sp.Expr]:
        """
        Boost a contravariant 4-vector (upper index).

        :param vectorU: 4-vector in Cartesian coordinates to be boosted.

        :return: Boosted 4-vector in Cartesian coordinates.
        """
        # Declare Lorentz matrix locally
        LorentzMatrix = self.LorentzMatrix

        # Declare indexed expression for vector to be boosted
        boosted_vectorU = ixp.zerorank1(dimension=4)

        # Perform Lorentz transformation on vectorU
        for mu in range(4):
            for nu in range(4):
                boosted_vectorU[mu] += LorentzMatrix[mu][nu] * vectorU[nu]
        return boosted_vectorU

    def inverse_boost_vectorU(self, vectorU: List[sp.Expr]) -> List[sp.Expr]:
        """
        Perform inverse boost on a contravariant 4-vector (upper index).

        :param vectorU: 4-vector in Cartesian coordinates to be inversely boosted.

        :return: Inversely boosted 4-vector in Cartesian coordinates.
        """
        # Declare inverse Lorentz matrix locally
        InverseLorentzMatrix = self.InverseLorentzMatrix

        # Declare indexed expression for vector to be inversely boosted
        inverse_boosted_vectorU = ixp.zerorank1(dimension=4)

        # Perform Lorentz transformation on vectorU
        for mu in range(4):
            for nu in range(4):
                inverse_boosted_vectorU[mu] += (
                    InverseLorentzMatrix[mu][nu] * vectorU[nu]
                )
        return inverse_boosted_vectorU

    def boost_vectorD(self, vectorD: List[sp.Expr]) -> List[sp.Expr]:
        """
        Boost a covariant 4-vector (lower index).

        :param vectorD: Vector in Cartesian coordinates to be boosted.

        :return: Boosted vector in Cartesian coordinates.
        """
        # Declare inverse Lorentz matrix locally
        InverseLorentzMatrix = self.InverseLorentzMatrix

        # Declare indexed expression for tensor to be boosted
        boosted_vectorD = ixp.zerorank1(dimension=4)

        # Perform Lorentz transformation
        for mu in range(4):
            for nu in range(4):
                boosted_vectorD[mu] += InverseLorentzMatrix[nu][mu] * vectorD[nu]
        return boosted_vectorD

    def boost_tensorDD(self, tensorDD: List[List[sp.Expr]]) -> List[List[sp.Expr]]:
        """
        Boost a tensor with two lower indices.

        :param tensorDD: tensor in Cartesian coordinates to be boosted.

        :return: Boosted tensor in Cartesian coordinates.
        """
        # Declare inverse Lorentz matrix locally
        InverseLorentzMatrix = self.InverseLorentzMatrix

        # Declare indexed expression for tensor to be boosted
        boosted_tensorDD = ixp.zerorank2(dimension=4)

        # Perform Lorentz transformation
        for mu1 in range(4):
            for mu2 in range(4):
                for nu1 in range(4):
                    for nu2 in range(4):
                        boosted_tensorDD[mu1][mu2] += (
                            InverseLorentzMatrix[nu1][mu1]
                            * InverseLorentzMatrix[nu2][mu2]
                            * tensorDD[nu1][nu2]
                        )
        return boosted_tensorDD

    def boost_tensorDDD(
        self, tensorDDD: List[List[List[sp.Expr]]]
    ) -> List[List[List[sp.Expr]]]:
        """
        Boost a tensor with three lower indices.

        :param tensorDDD: tensor in Cartesian coordinates to be boosted.

        :return: Boosted tensor in Cartesian coordinates.
        """
        # Declare inverse Lorentz matrix locally
        InverseLorentzMatrix = self.InverseLorentzMatrix

        # Declare indexed expression for tensor to be boosted
        boosted_tensorDDD = ixp.zerorank3(dimension=4)

        # Perform Lorentz transformation
        for mu1 in range(4):
            for mu2 in range(4):
                for mu3 in range(4):
                    for nu1 in range(4):
                        for nu2 in range(4):
                            for nu3 in range(4):
                                boosted_tensorDDD[mu1][mu2][mu3] += (
                                    InverseLorentzMatrix[nu1][mu1]
                                    * InverseLorentzMatrix[nu2][mu2]
                                    * InverseLorentzMatrix[nu3][mu3]
                                    * tensorDDD[nu1][nu2][nu3]
                                )
        return boosted_tensorDDD

    def boost_tensorDDDD(
        self, tensorDDDD: List[List[List[List[sp.Expr]]]]
    ) -> List[List[List[List[sp.Expr]]]]:
        """
        Boost a tensor with three four indices.

        :param tensorDDDD: tensor in Cartesian coordinates to be boosted.

        :return: Boosted tensor in Cartesian coordinates.
        """
        # Declare inverse Lorentz matrix locally
        InverseLorentzMatrix = self.InverseLorentzMatrix

        # Declare indexed expression for tensor to be boosted
        boosted_tensorDDDD = ixp.zerorank4(dimension=4)

        # Perform Lorentz transformation
        for mu1 in range(4):
            for mu2 in range(4):
                for mu3 in range(4):
                    for mu4 in range(4):
                        for nu1 in range(4):
                            for nu2 in range(4):
                                for nu3 in range(4):
                                    for nu4 in range(4):
                                        boosted_tensorDDDD[mu1][mu2][mu3][mu4] += (
                                            InverseLorentzMatrix[nu1][mu1]
                                            * InverseLorentzMatrix[nu2][mu2]
                                            * InverseLorentzMatrix[nu3][mu3]
                                            * InverseLorentzMatrix[nu4][mu4]
                                            * tensorDDDD[nu1][nu2][nu3][nu4]
                                        )
        return boosted_tensorDDDD


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

    # Declare symbolic variables for the 3-velocity vector
    vB = cast(List[sp.Expr], ixp.declarerank1("vB", dimension=3))

    # Initialize LorentzBoost with the symbolic 3-velocity vector
    lb = LorentzBoost(vB)

    # Prepare input_dict by initially setting it with lb.__dict__
    input_dict = lb.__dict__.copy()

    # Declare symbolic tensors for testing the boost functions
    vU = cast(List[sp.Expr], ixp.declarerank1("vU", dimension=4))
    vD = cast(List[sp.Expr], ixp.declarerank1("vD", dimension=4))
    tDD = cast(
        List[List[sp.Expr]], ixp.declarerank2("tDD", dimension=4, symmetry="sym01")
    )
    tDDD = cast(
        List[List[List[sp.Expr]]],
        ixp.declarerank3("tDDD", dimension=4, symmetry="sym012"),
    )
    tDDDD = cast(
        List[List[List[List[sp.Expr]]]],
        ixp.declarerank4("tDDDD", dimension=4, symmetry="sym0123"),
    )

    # Extend input_dict with the symbolic expressions for the boosted quantities
    input_dict["boosted_vU"] = lb.boost_vectorU(vU)
    input_dict["boosted_vD"] = lb.boost_vectorU(vD)
    input_dict["inverse_boosted_vU"] = lb.inverse_boost_vectorU(vU)
    input_dict["boosted_tDD"] = lb.boost_tensorDD(tDD)
    input_dict["boosted_tDDD"] = lb.boost_tensorDDD(tDDD)
    input_dict["boosted_tDDDD"] = lb.boost_tensorDDDD(tDDDD)

    # Process the dictionary of expressions
    results_dict = ve.process_dictionary_of_expressions(
        input_dict, fixed_mpfs_for_free_symbols=True
    )

    # Validate expressions
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        # File basename. If this is set to "trusted_module_test1", then
        #   trusted results_dict will be stored in tests/trusted_module_test1.py
        f"{os.path.splitext(os.path.basename(__file__))[0]}",
        results_dict,
    )
