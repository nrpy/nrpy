"""
This module constructs tetrads needed to
compute psi_4 (as well as other Weyl scalars
and invariants if desired).

Authors: Zachariah B. Etienne
         (zachetie **at** gmail **dot* com),
         and Patrick Nelson
"""

# Step 1.a: import all needed modules
from typing import List, cast, Union
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import nrpy.params as par  # NRPy+: Parameter interface
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
from nrpy.helpers.cached_functions import cached_simplify

# Step 1.b: Declare free parameter
# use_metric_to_construct_unit_normal=False: consistent with WeylScal4 ETK thorn.
par.register_param(bool, __name__, "use_metric_to_construct_unit_normal", False)


class Psi4Tetrads:
    """
    Class responsible for constructing tetrads needed for the computation
    of Weyl scalars like psi_4, and other GR invariants.

    :param CoordSystem: The coordinate system to be used. Default is 'Cartesian'.
    :param enable_rfm_precompute: Flag to enable/disable reference metric precomputation. Default is False.
    :param tetrad: quasiKinnersley = choice made in Baker, Campanelli, and Lousto. PRD 65, 044001 (2002)

    :raises ValueError: If an unsupported tetrad choice is made.
    """

    def __init__(
        self,
        CoordSystem: str = "Cartesian",
        enable_rfm_precompute: bool = False,
        tetrad: str = "quasiKinnersley",
    ):
        # Step 1.c: Check if tetrad choice is implemented:
        self.tetrad_choice = tetrad
        if self.tetrad_choice != "quasiKinnersley":
            raise ValueError(
                f"ERROR: tetrad = {self.tetrad_choice} currently unsupported!"
            )

        # Step 1.d: Given the chosen coordinate system, set up
        #           corresponding reference metric and needed
        #           reference metric quantities
        # The following function call sets up the reference metric
        #    and related quantities, including rescaling matrices ReDD,
        #    ReU, and hatted quantities.
        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]

        # Step 1.f: Import all ADM quantities as written in terms of BSSN quantities
        BtoA = BSSN_to_ADM(
            CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
        )

        # Step 2.a: Declare the Cartesian x,y,z in terms of
        #           xx0,xx1,xx2.
        x = rfm.xx_to_Cart[0]
        y = rfm.xx_to_Cart[1]
        z = rfm.xx_to_Cart[2]

        # Step 2.b: Declare detgamma and gammaUU from
        #           BSSN.ADM_in_terms_of_BSSN;
        #           simplify detgamma & gammaUU expressions,
        #           which expedites Psi4 codegen.
        detgamma = cached_simplify(BtoA.detgamma)
        gammaUU = ixp.zerorank2()
        for i in range(3):
            for j in range(i, 3):
                # The simplify() here is SLOW, so we try to use it sparingly.
                gammaUU[i][j] = gammaUU[j][i] = cached_simplify(BtoA.gammaUU[i][j])

        # Step 2.c: Define v1U and v2U
        v1UCart = [-y, x, sp.sympify(0)]
        v2UCart = [x, y, z]

        # Step 2.d: Construct the Jacobian d x_Cart^i / d xx^j
        # Step 2.e: Invert above Jacobian to get needed d xx^j / d x_Cart^i
        Jac_dUrfm_dDCartUD = rfm.Jac_dUrfm_dDCartUD

        # Step 2.e.i: Simplify expressions for d xx^j / d x_Cart^i:
        for i in range(3):
            for j in range(3):
                Jac_dUrfm_dDCartUD[i][j] = cached_simplify(Jac_dUrfm_dDCartUD[i][j])

        # Step 2.f: Transform v1U and v2U from the Cartesian to the xx^i basis
        v1U = ixp.zerorank1()
        v2U = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                v1U[i] += Jac_dUrfm_dDCartUD[i][j] * v1UCart[j]
                v2U[i] += Jac_dUrfm_dDCartUD[i][j] * v2UCart[j]

        # Step 2.g: Define v3U
        v3U = ixp.zerorank1()
        LeviCivitaSymbolDDD = ixp.LeviCivitaSymbol_dim3_rank3()
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    for d in range(3):
                        v3U[a] += (
                            sp.sqrt(detgamma)
                            * gammaUU[a][d]
                            * LeviCivitaSymbolDDD[d][b][c]
                            * v1U[b]
                            * v2U[c]
                        )

        # Step 2.g.i: Simplify expressions for v1U,v2U,v3U. This greatly expedites the C code generation (~10x faster)
        #             Drat. Simplification with certain versions of SymPy & coord systems results in a hang. Let's just
        #             evaluate the expressions so the most trivial optimizations can be performed.
        for a in range(3):
            v1U[a] = v1U[a].doit()  # cached_simplify(v1U[a])
            v2U[a] = v2U[a].doit()  # cached_simplify(v2U[a])
            v3U[a] = v3U[a].doit()  # cached_simplify(v3U[a])

        # Step 2.h: Define omega_{ij}
        omegaDD = ixp.zerorank2()
        gammaDD = BtoA.gammaDD

        def v_vectorDU(
            v1U: List[sp.Expr], v2U: List[sp.Expr], v3U: List[sp.Expr], i: int, a: int
        ) -> sp.Expr:
            if i == 0:
                return v1U[a]
            if i == 1:
                return v2U[a]
            if i == 2:
                return v3U[a]
            raise ValueError("ERROR: unknown vector!")

        def update_omega(
            omegaDD: List[List[sp.Expr]],
            i: int,
            j: int,
            v1U: List[sp.Expr],
            v2U: List[sp.Expr],
            v3U: List[sp.Expr],
            gammaDD: List[List[sp.Expr]],
        ) -> None:
            omegaDD[i][j] = sp.sympify(0)
            for a in range(3):
                for b in range(3):
                    omegaDD[i][j] += (
                        v_vectorDU(v1U, v2U, v3U, i, a)
                        * v_vectorDU(v1U, v2U, v3U, j, b)
                        * gammaDD[a][b]
                    )

        # Step 2.i: Define e^a_i. Note that:
        #           omegaDD[0][0] = \omega_{11} above;
        #           omegaDD[1][1] = \omega_{22} above, etc.
        # First e_1^a: Orthogonalize & normalize:
        e1U = ixp.zerorank1()
        update_omega(omegaDD, 0, 0, v1U, v2U, v3U, gammaDD)
        for a in range(3):
            e1U[a] = v1U[a] / sp.sqrt(omegaDD[0][0])

        # Next e_2^a: First orthogonalize:
        e2U = ixp.zerorank1()
        update_omega(omegaDD, 0, 1, e1U, v2U, v3U, gammaDD)
        for a in range(3):
            e2U[a] = v2U[a] - omegaDD[0][1] * e1U[a]
        # Then normalize:
        update_omega(omegaDD, 1, 1, e1U, e2U, v3U, gammaDD)
        for a in range(3):
            e2U[a] /= sp.sqrt(omegaDD[1][1])

        # Next e_3^a: First orthogonalize:
        e3U = ixp.zerorank1()
        update_omega(omegaDD, 0, 2, e1U, e2U, v3U, gammaDD)
        update_omega(omegaDD, 1, 2, e1U, e2U, v3U, gammaDD)
        for a in range(3):
            e3U[a] = v3U[a] - omegaDD[0][2] * e1U[a] - omegaDD[1][2] * e2U[a]
        # Then normalize:
        update_omega(omegaDD, 2, 2, e1U, e2U, e3U, gammaDD)
        for a in range(3):
            e3U[a] /= sp.sqrt(omegaDD[2][2])

        # Step 2.j: Construct l^mu, n^mu, and m^mu, based on r^mu, theta^mu, phi^mu, and u^mu:
        r4U = cast(List[sp.Expr], ixp.zerorank1(dimension=4))
        # u4U = cast(List[Union[int, sp.Expr]], ixp.zerorank1(dimension=4))
        u4U = cast(List[Union[int, sp.Expr]], ixp.zerorank1(dimension=4))
        theta4U = cast(List[sp.Expr], ixp.zerorank1(dimension=4))
        phi4U = cast(List[sp.Expr], ixp.zerorank1(dimension=4))

        for a in range(3):
            r4U[a + 1] = e2U[a]
            theta4U[a + 1] = e3U[a]
            phi4U[a + 1] = e1U[a]

        # FIXME? assumes alpha=1, beta^i = 0
        if par.parval_from_str("use_metric_to_construct_unit_normal"):
            # Eq. 2.116 in Baumgarte & Shapiro:
            #  n^mu = {1/alpha, -beta^i/alpha}. Note that n_mu = {alpha,0}, so n^mu n_mu = -1.
            Bq = BSSN_quantities[
                CoordSystem + "_rfm_precompute"
                if enable_rfm_precompute
                else CoordSystem
            ]
            u4U[0] = 1 / Bq.alpha
            for i in range(3):
                u4U[i + 1] = -Bq.betaU[i] / Bq.alpha
        else:
            u4U[0] = 1

        self.l4U = ixp.zerorank1(dimension=4)
        self.n4U = ixp.zerorank1(dimension=4)
        self.mre4U = ixp.zerorank1(dimension=4)
        self.mim4U = ixp.zerorank1(dimension=4)

        # M_SQRT1_2 = 1 / sqrt(2) (defined in math.h on Linux)
        M_SQRT1_2 = par.register_CodeParameter(
            "#define",
            __name__,
            "M_SQRT1_2",
            "",
            commondata=False,
            add_to_parfile=False,
            add_to_set_CodeParameters_h=False,
            add_to_glb_code_params_dict=False,
        )
        isqrt2 = M_SQRT1_2  # 1/sp.sqrt(2) <- SymPy drops precision to 15 sig. digits in unit tests
        for mu in range(4):
            self.l4U[mu] = isqrt2 * (u4U[mu] + r4U[mu])
            self.n4U[mu] = isqrt2 * (u4U[mu] - r4U[mu])
            self.mre4U[mu] = isqrt2 * theta4U[mu]
            self.mim4U[mu] = isqrt2 * phi4U[mu]


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
        psi4 = Psi4Tetrads(CoordSystem=Coord_in, enable_rfm_precompute=enable_rfm_pre)
        results_dict = ve.process_dictionary_of_expressions(
            psi4.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{Coord}",
            results_dict,
        )
