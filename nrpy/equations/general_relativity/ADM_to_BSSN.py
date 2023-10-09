"""
Construct expressions for BSSN quantities in terms of ADM quantities.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

# Step 1: Import needed core NRPy+ modules
from typing import Sequence, Tuple, List
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import nrpy.params as par  # NRPy+: Parameter interface
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support


# NRPy+: This module depends on the parameter EvolvedConformalFactor_cf,
#        which is defined in BSSN.BSSN_quantities
import nrpy.equations.general_relativity.BSSN_quantities  # pylint: disable=unused-import


class ADM_to_BSSN:
    """Sets up and stores expressions for BSSN variables in terms of ADM quantities."""

    def __init__(
        self,
        gammaDD: Sequence[Sequence[sp.Expr]],
        KDD: Sequence[Sequence[sp.Expr]],
        betaU: Sequence[sp.Expr],
        BU: Sequence[sp.Expr],
        CoordSystem: str = "Cartesian",
        enable_rfm_precompute: bool = False,
        compute_cf_only: bool = False,
    ):
        """
        Initialize the ADM_to_BSSN class.

        :param gammaDD: Input 3x3 metric tensor.
        :param KDD: Input 3x3 extrinsic curvature tensor.
        :param betaU: Input 3-vector shift.
        :param BU: Input 3-vector B.
        :param CoordSystem: String describing the coordinate system of the inputs.
        :param enable_rfm_precompute: Boolean flag to enable reference metric precomputation
        :param compute_cf_only: Boolean flag to compute only the conformal factor.
        """
        # Step 2: All ADM quantities were input into this function in the Spherical or Cartesian
        #         basis, as functions of r,th,ph or x,y,z, respectively. In Steps 1 and 2 above,
        #         we converted them to the xx0,xx1,xx2 basis, and as functions of xx0,xx1,xx2.
        #         Here we convert ADM quantities to their BSSN Curvilinear counterparts:

        # Step 2.a: Convert ADM $\gamma_{ij}$ to BSSN $\bar{gamma}_{ij}$:
        #           We have (Eqs. 2 and 3 of [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):
        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]
        if gammaDD is None:
            gammaDD = ixp.declarerank2("gammaDD", symmetry="sym01")
        if (
            KDD is None
        ):  # Use "is None" instead of "==None", as the former is more correct.
            KDD = ixp.declarerank2("KDD", symmetry="sym01")
        if (
            betaU is None
        ):  # Use "is None" instead of "==None", as the former is more correct.
            betaU = ixp.declarerank1("betaU")
        if (
            BU is None
        ):  # Use "is None" instead of "==None", as the former is more correct.
            BU = ixp.declarerank1("BU")

        # Step 2.b: Set the conformal factor variable cf, which is set
        #           by the "BSSN_quantities::EvolvedConformalFactor_cf" parameter. For example if
        #           "EvolvedConformalFactor_cf" is set to "phi", we can use Eq. 3 of
        #           [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf),
        #           which in arbitrary coordinates is written:

        # \bar{gamma}_{ij} = (\frac{\bar{gamma}}{gamma})^{1/3}*gamma_{ij}.
        gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)  # _gammaUU unused.
        self.gammabarDD = ixp.zerorank2()
        self.hDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self.gammabarDD[i][j] = (rfm.detgammahat / gammaDET) ** (
                    sp.Rational(1, 3)
                ) * gammaDD[i][j]
                self.hDD[i][j] = (self.gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[
                    i
                ][j]
        gammabarUU, gammabarDET = ixp.symm_matrix_inverter3x3(self.gammabarDD)

        self.cf = sp.sympify(0)

        EvolvedConformalFactor_cf = par.parval_from_str("EvolvedConformalFactor_cf")
        if EvolvedConformalFactor_cf == "phi":
            # phi = \frac{1}{12} log(\frac{gamma}{\bar{gamma}}).
            self.cf = sp.Rational(1, 12) * sp.log(gammaDET / gammabarDET)
        elif EvolvedConformalFactor_cf == "chi":
            # chi = exp(-4*phi) = exp(-4*\frac{1}{12}*(\frac{gamma}{\bar{gamma}}))
            #      = exp(-\frac{1}{3}*log(\frac{gamma}{\bar{gamma}})) = (\frac{gamma}{\bar{gamma}})^{-1/3}.
            #
            self.cf = (gammaDET / gammabarDET) ** (-sp.Rational(1, 3))
        elif EvolvedConformalFactor_cf == "W":
            # W = exp(-2*phi) = exp(-2*\frac{1}{12}*log(\frac{gamma}{\bar{gamma}}))
            #   = exp(-\frac{1}{6}*log(\frac{gamma}{\bar{gamma}})) = (\frac{gamma}{bar{gamma}})^{-1/6}.
            self.cf = (gammaDET / gammabarDET) ** (-sp.Rational(1, 6))
        else:
            raise ValueError(
                f"Error EvolvedConformalFactor_cf type = {EvolvedConformalFactor_cf} unknown."
            )
        if compute_cf_only:
            return

        # Step 2.c: Convert the extrinsic curvature K_{ij} to the trace-free extrinsic
        #           curvature \bar{A}_{ij}, plus the trace of the extrinsic curvature K,
        #           where (Eq. 3 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)):

        # \bar{gamma}_{ij} = (\frac{\bar{gamma}}{gamma})^{1/3}*gamma_{ij}.
        # K = gamma^{ij} K_{ij}, and
        # \bar{A}_{ij} &= (\frac{\bar{gamma}}{gamma})^{1/3}*(K_{ij} - \frac{1}{3}*gamma_{ij}*K)
        self.trK = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                self.trK += gammaUU[i][j] * KDD[i][j]

        self.AbarDD = ixp.zerorank2()
        self.aDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self.AbarDD[i][j] = (rfm.detgammahat / gammaDET) ** (
                    sp.Rational(1, 3)
                ) * (KDD[i][j] - sp.Rational(1, 3) * gammaDD[i][j] * self.trK)
                self.aDD[i][j] = self.AbarDD[i][j] / rfm.ReDD[i][j]

        # Step 2.d: \bar{Lambda}^i
        # \bar{Lambda}^i = \bar{gamma}^{jk}(\bar{Gamma}^i_{jk} - \hat{Gamma}^i_{jk}).
        # First compute Christoffel symbols \bar{Gamma}^i_{jk}, with respect to barred metric:
        GammabarUDD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        GammabarUDD[i][j][k] += (
                            sp.Rational(1, 2)
                            * gammabarUU[i][l]
                            * (
                                sp.diff(self.gammabarDD[l][j], rfm.xx[k])
                                + sp.diff(self.gammabarDD[l][k], rfm.xx[j])
                                - sp.diff(self.gammabarDD[j][k], rfm.xx[l])
                            )
                        )
        # Next evaluate \bar{Lambda}^i, based on GammabarUDD above and GammahatUDD
        #       (from the reference metric):
        LambdabarU = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    LambdabarU[i] += gammabarUU[j][k] * (
                        GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k]
                    )
        for i in range(3):
            # We evaluate LambdabarU[i] here to ensure proper cancellations. If these cancellations
            #   are not applied, certain expressions (e.g., lambdaU[0] in StaticTrumpet) will
            #   cause SymPy's (v1.5+) CSE algorithm to hang
            LambdabarU[i] = LambdabarU[i].doit()
        lambdaU = ixp.zerorank1()
        for i in range(3):
            lambdaU[i] = LambdabarU[i] / rfm.ReU[i]

        # Step 2.e: Rescale beta^i and B^i according to the prescription described in
        #         the [BSSN in curvilinear coordinates tutorial module](Tutorial-BSSNCurvilinear.ipynb)
        #         (also [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):
        #
        # \mathcal{V}^i &= beta^i/(ReU[i])
        # \mathcal{B}^i &= B^i/(ReU[i])
        self.vetU = ixp.zerorank1()
        self.betU = ixp.zerorank1()
        for i in range(3):
            self.vetU[i] = betaU[i] / rfm.ReU[i]
            self.betU[i] = BU[i] / rfm.ReU[i]


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

    # Cannot import InitialData_Spherical here, as that would result in a circular dependency.
    def StaticTrumpet_standalone() -> (
        Tuple[
            List[List[sp.Expr]],
            List[List[sp.Expr]],
            sp.Expr,
            List[sp.Expr],
            List[sp.Expr],
        ]
    ):
        """
        Set up Static Trumpet initial data.

        Cannot import InitialData_Spherical here, as that would result in a circular dependency.
        """
        M, r, th = sp.symbols("M r th", real=True)
        psi0 = sp.sqrt(1 + M / r)
        IDgammaDD = ixp.zerorank2()
        IDgammaDD[0][0] = psi0**4
        IDgammaDD[1][1] = psi0**4 * r**2
        IDgammaDD[2][2] = psi0**4 * r**2 * sp.sin(th) ** 2
        IDKDD = ixp.zerorank2()
        IDKDD[0][0] = -M / r**2
        IDKDD[1][1] = M
        IDKDD[2][2] = M * sp.sin(th) ** 2
        IDalpha = r / (r + M)
        IDbetaU = ixp.zerorank1()
        IDbetaU[0] = M * r / (r + M) ** 2
        IDBU = ixp.zerorank1()
        return IDgammaDD, IDKDD, IDalpha, IDbetaU, IDBU

    ingammaDD, inKDD, inalpha, inbetaU, inBU = StaticTrumpet_standalone()
    adm2bssn = ADM_to_BSSN(ingammaDD, inKDD, inbetaU, inBU, "Spherical")
    results_dict = ve.process_dictionary_of_expressions(
        adm2bssn.__dict__, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        # File basename. If this is set to "trusted_module_test1", then
        #   trusted results_dict will be stored in tests/trusted_module_test1.py
        f"{os.path.splitext(os.path.basename(__file__))[0]}_StaticTrumpet",
        results_dict,
    )
