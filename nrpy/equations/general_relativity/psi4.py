"""
Constructs psi_4, with respect to an arbitrary, unspecified tetrad.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

# Step 1.a: import needed modules
from typing import List, cast
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
from nrpy.equations.general_relativity.psi4_tetrads import Psi4Tetrads


class Psi4:
    """
    Class responsible for constructing psi_4 with respect to an arbitrary, unspecified tetrad.

    :param CoordSystem: The coordinate system to be used. Default is 'Cartesian'.
    :param enable_rfm_precompute: Flag to enable/disable reference metric precomputation. Default is False.
    :param tetrad: Specify the tetrad explicitly? Default is True. False means to leave it symbolic.
    """

    def __init__(
        self,
        CoordSystem: str = "Cartesian",
        enable_rfm_precompute: bool = False,
        tetrad: str = "leave_symbolic",
    ):
        # Step 1.b: Import all ADM quantities as written in terms of BSSN quantities
        BtoA = BSSN_to_ADM(
            CoordSystem=CoordSystem,
            enable_rfm_precompute=enable_rfm_precompute,
        )

        # Step 1.c: Set up tetrad vectors
        mre4U: List[sp.Expr]
        mim4U: List[sp.Expr]
        n4U: List[sp.Expr]
        if tetrad == "leave_symbolic":
            # For code validation against NRPy+ psi_4 tutorial module (Tutorial-Psi4.ipynb);
            #   ensures a more complete code validation.
            mre4U = cast(List[sp.Expr], ixp.declarerank1("mre4U", dimension=4))
            mim4U = cast(List[sp.Expr], ixp.declarerank1("mim4U", dimension=4))
            n4U = cast(List[sp.Expr], ixp.declarerank1("n4U", dimension=4))
        else:
            BP4t = Psi4Tetrads(
                CoordSystem=CoordSystem,
                enable_rfm_precompute=enable_rfm_precompute,
                tetrad=tetrad,
            )
            mre4U = BP4t.mre4U
            mim4U = BP4t.mim4U
            n4U = BP4t.n4U

        # Step 2: Construct the (rank-4) Riemann curvature tensor associated with the ADM 3-metric:
        RDDDD = ixp.zerorank4()
        gammaDDdDD = BtoA.gammaDDdDD

        for i in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        RDDDD[i][k][l][m] = sp.Rational(1, 2) * (
                            +gammaDDdDD[i][m][k][l]
                            + gammaDDdDD[k][l][i][m]
                            - gammaDDdDD[i][l][k][m]
                            - gammaDDdDD[k][m][i][l]
                        )

        # ... then we add the term on the right:
        gammaDD = BtoA.gammaDD
        GammaUDD = BtoA.GammaUDD

        for i in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            for p in range(3):
                                RDDDD[i][k][l][m] += gammaDD[n][p] * (
                                    GammaUDD[n][k][l] * GammaUDD[p][i][m]
                                    - GammaUDD[n][k][m] * GammaUDD[p][i][l]
                                )

        # Step 3: Construct the (rank-4) tensor in term 1 of psi_4 (referring to Eq 5.1 in
        #   Baker, Campanelli, Lousto (2001); https://arxiv.org/pdf/gr-qc/0104063.pdf
        rank4term1DDDD = ixp.zerorank4()
        KDD = BtoA.KDD

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        rank4term1DDDD[i][j][k][l] = (
                            RDDDD[i][j][k][l]
                            + KDD[i][k] * KDD[l][j]
                            - KDD[i][l] * KDD[k][j]
                        )

        # Step 4: Construct the (rank-3) tensor in term 2 of psi_4 (referring to Eq 5.1 in
        #   Baker, Campanelli, Lousto (2001); https://arxiv.org/pdf/gr-qc/0104063.pdf
        rank3term2DDD = ixp.zerorank3()
        KDDdD = BtoA.KDDdD

        for j in range(3):
            for k in range(3):
                for l in range(3):
                    rank3term2DDD[j][k][l] = sp.Rational(1, 2) * (
                        KDDdD[j][k][l] - KDDdD[j][l][k]
                    )

        # ... then we construct the second term in this sum:
        #  \Gamma^{p}_{j[k} K_{l]p} = \frac{1}{2} (\Gamma^{p}_{jk} K_{lp}-\Gamma^{p}_{jl} K_{kp}):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for p in range(3):
                        rank3term2DDD[j][k][l] += sp.Rational(1, 2) * (
                            GammaUDD[p][j][k] * KDD[l][p]
                            - GammaUDD[p][j][l] * KDD[k][p]
                        )

        # Finally, we multiply the term by $-8$:
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    rank3term2DDD[j][k][l] *= sp.sympify(-8)

        # Step 5: Construct the (rank-2) tensor in term 3 of psi_4 (referring to Eq 5.1 in
        #   Baker, Campanelli, Lousto (2001); https://arxiv.org/pdf/gr-qc/0104063.pdf

        # Step 5.1: Construct 3-Ricci tensor R_{ij} = gamma^{im} R_{ijml}
        RDD = ixp.zerorank2()
        gammaUU = BtoA.gammaUU
        for j in range(3):
            for l in range(3):
                for i in range(3):
                    for m in range(3):
                        RDD[j][l] += gammaUU[i][m] * RDDDD[i][j][m][l]

        # Step 5.2: Construct K^p_l = gamma^{pi} K_{il}
        KUD = ixp.zerorank2()
        for p in range(3):
            for l in range(3):
                for i in range(3):
                    KUD[p][l] += gammaUU[p][i] * KDD[i][l]

        # Step 5.3: Construct trK = gamma^{ij} K_{ij}
        trK = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                trK += gammaUU[i][j] * KDD[i][j]

        # Next we put these terms together to construct the entire term in parentheses:
        # +4 \left(R_{jl} - K_{jp} K^p_l + K K_{jl} \right),
        rank2term3DD = ixp.zerorank2()
        for j in range(3):
            for l in range(3):
                rank2term3DD[j][l] = RDD[j][l] + trK * KDD[j][l]
                for p in range(3):
                    rank2term3DD[j][l] += -KDD[j][p] * KUD[p][l]
        # Finally we multiply by +4:
        for j in range(3):
            for l in range(3):
                rank2term3DD[j][l] *= sp.sympify(4)

        # Step 6: Construct real & imaginary parts of psi_4
        #         by contracting constituent rank 2, 3, and 4
        #         tensors with input tetrads mre4U, mim4U, & n4U.

        def tetrad_product__Real_psi4(
            n: List[sp.Expr],
            Mre: List[sp.Expr],
            Mim: List[sp.Expr],
            mu: int,
            nu: int,
            eta: int,
            delta: int,
        ) -> sp.Expr:
            return cast(
                sp.Expr,
                +n[mu] * Mre[nu] * n[eta] * Mre[delta]
                - n[mu] * Mim[nu] * n[eta] * Mim[delta],
            )

        def tetrad_product__Imag_psi4(
            n: List[sp.Expr],
            Mre: List[sp.Expr],
            Mim: List[sp.Expr],
            mu: int,
            nu: int,
            eta: int,
            delta: int,
        ) -> sp.Expr:
            return cast(
                sp.Expr,
                -n[mu] * Mre[nu] * n[eta] * Mim[delta]
                - n[mu] * Mim[nu] * n[eta] * Mre[delta],
            )

        # We split psi_4 into three pieces, to expedite & possibly parallelize C code generation.
        self.psi4_re_pt = [sp.sympify(0), sp.sympify(0), sp.sympify(0)]
        self.psi4_im_pt = [sp.sympify(0), sp.sympify(0), sp.sympify(0)]
        # First term:
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        self.psi4_re_pt[0] += rank4term1DDDD[i][j][k][
                            l
                        ] * tetrad_product__Real_psi4(
                            n4U, mre4U, mim4U, i + 1, j + 1, k + 1, l + 1
                        )
                        self.psi4_im_pt[0] += rank4term1DDDD[i][j][k][
                            l
                        ] * tetrad_product__Imag_psi4(
                            n4U, mre4U, mim4U, i + 1, j + 1, k + 1, l + 1
                        )

        # Second term:
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    self.psi4_re_pt[1] += (
                        rank3term2DDD[j][k][l]
                        * sp.Rational(1, 2)
                        * (
                            +tetrad_product__Real_psi4(
                                n4U, mre4U, mim4U, 0, j + 1, k + 1, l + 1
                            )
                            - tetrad_product__Real_psi4(
                                n4U, mre4U, mim4U, j + 1, 0, k + 1, l + 1
                            )
                        )
                    )
                    self.psi4_im_pt[1] += (
                        rank3term2DDD[j][k][l]
                        * sp.Rational(1, 2)
                        * (
                            +tetrad_product__Imag_psi4(
                                n4U, mre4U, mim4U, 0, j + 1, k + 1, l + 1
                            )
                            - tetrad_product__Imag_psi4(
                                n4U, mre4U, mim4U, j + 1, 0, k + 1, l + 1
                            )
                        )
                    )
        # Third term:
        for j in range(3):
            for l in range(3):
                self.psi4_re_pt[2] += rank2term3DD[j][l] * (
                    sp.Rational(1, 4)
                    * (
                        +tetrad_product__Real_psi4(
                            n4U, mre4U, mim4U, 0, j + 1, 0, l + 1
                        )
                        - tetrad_product__Real_psi4(
                            n4U, mre4U, mim4U, j + 1, 0, 0, l + 1
                        )
                        - tetrad_product__Real_psi4(
                            n4U, mre4U, mim4U, 0, j + 1, l + 1, 0
                        )
                        + tetrad_product__Real_psi4(
                            n4U, mre4U, mim4U, j + 1, 0, l + 1, 0
                        )
                    )
                )
                self.psi4_im_pt[2] += rank2term3DD[j][l] * (
                    sp.Rational(1, 4)
                    * (
                        +tetrad_product__Imag_psi4(
                            n4U, mre4U, mim4U, 0, j + 1, 0, l + 1
                        )
                        - tetrad_product__Imag_psi4(
                            n4U, mre4U, mim4U, j + 1, 0, 0, l + 1
                        )
                        - tetrad_product__Imag_psi4(
                            n4U, mre4U, mim4U, 0, j + 1, l + 1, 0
                        )
                        + tetrad_product__Imag_psi4(
                            n4U, mre4U, mim4U, j + 1, 0, l + 1, 0
                        )
                    )
                )


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

    for in_tetrad in ["quasiKinnersley", "leave_symbolic"]:
        for Coord in [
            "Spherical",
            "SinhSpherical",
            "SinhSpherical_rfm_precompute",
        ]:
            enable_rfm_pre = False
            Coord_in = Coord
            if "rfm_precompute" in Coord:
                Coord_in = Coord.replace("_rfm_precompute", "")
                enable_rfm_pre = True
            psi4 = Psi4(
                CoordSystem=Coord_in,
                enable_rfm_precompute=enable_rfm_pre,
                tetrad=in_tetrad,
            )
            results_dict = ve.process_dictionary_of_expressions(
                psi4.__dict__, fixed_mpfs_for_free_symbols=True
            )
            ve.compare_or_generate_trusted_results(
                os.path.abspath(__file__),
                os.getcwd(),
                # File basename. If this is set to "trusted_module_test1", then
                #   trusted results_dict will be stored in tests/trusted_module_test1.py
                f"{os.path.splitext(os.path.basename(__file__))[0]}_{in_tetrad}_{Coord}",
                results_dict,
            )
