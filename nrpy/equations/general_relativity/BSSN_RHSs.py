"""
This module constructs the right-hand sides (RHS)
  expressions of the BSSN time evolution equations.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

# Step 1.a: import all needed modules from NRPy+:
from typing import Dict, List
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends

import nrpy.params as par  # NRPy+: Parameter interface
import nrpy.grid as gri
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity import T4munu

have_already_called_BSSN_RHSs_function = False

# Step 1.b: Set the coordinate system for the numerical grid:
#  DO NOT SET IN STANDALONE PYTHON MODULE
# par.set_parval_from_str("reference_metric::CoordSystem","Spherical")


class BSSNRHSs:
    """
    Constructs the right-hand side (RHS) expressions for the BSSN (Baumgarte-Shapiro-Shibata-Nakamura) time
    evolution equations. These equations are used for numerical relativity simulations and
    are derived from Einstein's field equations.
    """

    def __init__(
        self, CoordSystem: str = "Cartesian", enable_rfm_precompute: bool = False
    ):
        enable_T4munu = par.parval_from_str("enable_T4munu")

        # Step 1.c: Given the chosen coordinate system, set up
        #           corresponding reference metric and needed
        #           reference metric quantities
        # The following function call sets up the reference metric
        #    and related quantities, including rescaling matrices ReDD,
        #    ReU, and hatted quantities.
        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]

        # Step 1.e: Import all basic (unrescaled) BSSN scalars & tensors
        Bq = BSSN_quantities[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]
        gammabarDD = Bq.gammabarDD
        AbarDD = Bq.AbarDD
        LambdabarU = Bq.LambdabarU
        trK = Bq.trK
        alpha = Bq.alpha
        betaU = Bq.betaU

        # Step 1.f: Import all needed rescaled BSSN tensors:
        cf = Bq.cf
        lambdaU = Bq.lambdaU

        # Step 2.a.i: Import derivative expressions for betaU defined in the BSSN_quantities module:
        betaU_dD = Bq.betaU_dD
        betaU_dDD = Bq.betaU_dDD
        # Step 2.a.ii: Import derivative expression for gammabarDD
        gammabarDD_dupD = Bq.gammabarDD_dupD

        # fmt: off

        # Step 2.a.iii: First term of \partial_t \bar{\gamma}_{i j} right-hand side:
        # \beta^k \bar{\gamma}_{ij,k} + \beta^k_{,i} \bar{\gamma}_{kj} + \beta^k_{,j} \bar{\gamma}_{ik}
        gammabar_rhsDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    gammabar_rhsDD[i][j] += betaU[k] * gammabarDD_dupD[i][j][k] + betaU_dD[k][i] * gammabarDD[k][j] \
                                            + betaU_dD[k][j] * gammabarDD[i][k]

        # Step 2.b.i: First import \bar{A}_{ij}'s contraction trAbar = \bar{A}^k_k
        #           from BSSN_quantities
        trAbar = Bq.trAbar

        # Step 2.b.ii: Import detgammabar from BSSN_quantities:
        detgammabar = Bq.detgammabar
        detgammabar_dD = Bq.detgammabar_dD

        # Step 2.b.ii: Compute the contraction \bar{D}_k \beta^k = \beta^k_{,k} + \frac{\beta^k \bar{\gamma}_{,k}}{2 \bar{\gamma}}
        Dbarbetacontraction = sp.sympify(0)
        for k in range(3):
            Dbarbetacontraction += betaU_dD[k][k] + betaU[k] * detgammabar_dD[k] / (2 * detgammabar)

        # Step 2.b.iii: Second term of \partial_t \bar{\gamma}_{i j} right-hand side:
        # \frac{2}{3} \bar{\gamma}_{i j} \left (\alpha \bar{A}_{k}^{k} - \bar{D}_{k} \beta^{k}\right )
        for i in range(3):
            for j in range(3):
                gammabar_rhsDD[i][j] += sp.Rational(2, 3) * gammabarDD[i][j] * (alpha * trAbar - Dbarbetacontraction)

        # Step 2.c: Third term of \partial_t \bar{\gamma}_{i j} right-hand side:
        # -2 \alpha \bar{A}_{ij}
        for i in range(3):
            for j in range(3):
                gammabar_rhsDD[i][j] += -2 * alpha * AbarDD[i][j]

        # Step 3.a: First term of \partial_t \bar{A}_{i j}:
        # \beta^k \partial_k \bar{A}_{ij} + \partial_i \beta^k \bar{A}_{kj} + \partial_j \beta^k \bar{A}_{ik}

        # First define AbarDD_dupD:
        AbarDD_dupD = Bq.AbarDD_dupD # From Bq.AbarUU_AbarUD_trAbar_AbarDD_dD()

        Abar_rhsDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    Abar_rhsDD[i][j] += betaU[k] * AbarDD_dupD[i][j][k] + betaU_dD[k][i] * AbarDD[k][j] \
                                        + betaU_dD[k][j] * AbarDD[i][k]

        # Step 3.b: Second term of \partial_t \bar{A}_{i j}:
        # - (2/3) \bar{A}_{i j} \bar{D}_{k} \beta^{k} - 2 \alpha \bar{A}_{i k} {\bar{A}^{k}}_{j} + \alpha \bar{A}_{i j} K
        gammabarUU = Bq.gammabarUU  # From Bq.gammabar__inverse_and_derivs()
        AbarUD = Bq.AbarUD  # From Bq.AbarUU_AbarUD_trAbar()
        for i in range(3):
            for j in range(3):
                Abar_rhsDD[i][j] += -sp.Rational(2, 3) * AbarDD[i][j] * Dbarbetacontraction + alpha * AbarDD[i][j] * trK
                for k in range(3):
                    Abar_rhsDD[i][j] += -2 * alpha * AbarDD[i][k] * AbarUD[k][j]

        # Step 3.c.i: Define partial derivatives of \phi in terms of evolved quantity "cf":
        phi_dD = Bq.phi_dD
        phi_dupD = Bq.phi_dupD
        exp_m4phi = Bq.exp_m4phi
        phi_dBarD = Bq.phi_dBarD  # phi_dBarD = Dbar_i phi = phi_dD (since phi is a scalar)
        phi_dBarDD = Bq.phi_dBarDD  # phi_dBarDD = Dbar_i Dbar_j phi (covariant derivative)

        # Step 3.c.ii: Define RbarDD
        RbarDD = Bq.RbarDD

        # Step 3.c.iii: Define first and second derivatives of \alpha, as well as
        #         \bar{D}_i \bar{D}_j \alpha, which is defined just like phi
        alpha_dD = ixp.declarerank1("alpha_dD")
        alpha_dDD = ixp.declarerank2("alpha_dDD", symmetry="sym01")
        alpha_dBarD = alpha_dD
        alpha_dBarDD = ixp.zerorank2()
        GammabarUDD = Bq.GammabarUDD  # Defined in Bq.gammabar__inverse_and_derivs()
        for i in range(3):
            for j in range(3):
                alpha_dBarDD[i][j] = alpha_dDD[i][j]
                for k in range(3):
                    alpha_dBarDD[i][j] += - GammabarUDD[k][i][j] * alpha_dD[k]

        # Step 3.c.iv: Define the terms in curly braces:
        curlybrackettermsDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                curlybrackettermsDD[i][j] = -2 * alpha * phi_dBarDD[i][j] + 4 * alpha * phi_dBarD[i] * phi_dBarD[j] \
                                            + 2 * alpha_dBarD[i] * phi_dBarD[j] \
                                            + 2 * alpha_dBarD[j] * phi_dBarD[i] \
                                            - alpha_dBarDD[i][j] + alpha * RbarDD[i][j]

        # Step 3.c.v: Compute the trace:
        curlybracketterms_trace = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                curlybracketterms_trace += gammabarUU[i][j] * curlybrackettermsDD[i][j]

        # Step 3.c.vi: Third and final term of Abar_rhsDD[i][j]:
        for i in range(3):
            for j in range(3):
                Abar_rhsDD[i][j] += exp_m4phi * (curlybrackettermsDD[i][j] -
                                                 sp.Rational(1, 3) * gammabarDD[i][j] * curlybracketterms_trace)

        # Step 4: Right-hand side of conformal factor variable "cf". Supported
        #          options include: cf=phi, cf=W=e^(-2*phi) (default), and cf=chi=e^(-4*phi)
        # \partial_t phi = \left[\beta^k \partial_k \phi \right] <- TERM 1
        #                  + \frac{1}{6} \left (\bar{D}_{k} \beta^{k} - \alpha K \right ) <- TERM 2
        self.cf_rhs = sp.Rational(1, 6) * (Dbarbetacontraction - alpha * trK)  # Term 2
        for k in range(3):
            self.cf_rhs += betaU[k] * phi_dupD[k]  # Term 1

        # Next multiply to convert phi_rhs to cf_rhs.
        EvolvedConformalFactor_cf = par.parval_from_str("EvolvedConformalFactor_cf")
        if EvolvedConformalFactor_cf == "phi":
            pass  # do nothing; cf_rhs = phi_rhs
        elif EvolvedConformalFactor_cf == "W":
            self.cf_rhs *= -2 * cf  # cf_rhs = -2*cf*phi_rhs
        elif EvolvedConformalFactor_cf == "chi":
            self.cf_rhs *= -4 * cf  # cf_rhs = -4*cf*phi_rhs
        else:
            raise ValueError(f"EvolvedConformalFactor_cf == {EvolvedConformalFactor_cf} unsupported!")

        # Step 5: right-hand side of trK (trace of extrinsic curvature):
        # \partial_t K = \beta^k \partial_k K <- TERM 1
        #           + \frac{1}{3} \alpha K^{2} <- TERM 2
        #           + \alpha \bar{A}_{i j} \bar{A}^{i j} <- TERM 3
        #           - - e^{-4 \phi} (\bar{D}_{i} \bar{D}^{i} \alpha + 2 \bar{D}^{i} \alpha \bar{D}_{i} \phi ) <- TERM 4
        # TERM 2:
        self.trK_rhs = sp.Rational(1, 3) * alpha * trK * trK
        trK_dupD = ixp.declarerank1("trK_dupD")
        for i in range(3):
            # TERM 1:
            self.trK_rhs += betaU[i] * trK_dupD[i]
        for i in range(3):
            for j in range(3):
                # TERM 4:
                self.trK_rhs += -exp_m4phi * gammabarUU[i][j] * (alpha_dBarDD[i][j] + 2 * alpha_dBarD[j] * phi_dBarD[i])
        AbarUU = Bq.AbarUU
        for i in range(3):
            for j in range(3):
                # TERM 3:
                self.trK_rhs += alpha * AbarDD[i][j] * AbarUU[i][j]

        # Step 6: right-hand side of \partial_t \bar{\Lambda}^i:
        # \partial_t \bar{\Lambda}^i = \beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k <- TERM 1
        #                            + \bar{\gamma}^{j k} \hat{D}_{j} \hat{D}_{k} \beta^{i} <- TERM 2
        #                            + \frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j} <- TERM 3
        #                            + \frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j} <- TERM 4
        #                            - 2 \bar{A}^{i j} (\partial_{j} \alpha - 6 \partial_{j} \phi) <- TERM 5
        #                            + 2 \alpha \bar{A}^{j k} \Delta_{j k}^{i} <- TERM 6
        #                            - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K <- TERM 7

        # Step 6.a: Term 1 of \partial_t \bar{\Lambda}^i: \beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k
        # First we declare \bar{\Lambda}^i and \bar{\Lambda}^i_{,j} in terms of \lambda^i and \lambda^i_{,j}
        self.LambdabarU_dupD = ixp.zerorank2()  # Used on the RHS of the Gamma-driving shift conditions
        lambdaU_dupD = ixp.declarerank2("lambdaU_dupD", symmetry="nosym")
        for i in range(3):
            for j in range(3):
                self.LambdabarU_dupD[i][j] = lambdaU_dupD[i][j] * rfm.ReU[i] + lambdaU[i] * rfm.ReUdD[i][j]

        self.Lambdabar_rhsU = ixp.zerorank1()  # Used on the RHS of the Gamma-driving shift conditions
        for i in range(3):
            for k in range(3):
                self.Lambdabar_rhsU[i] += betaU[k] * self.LambdabarU_dupD[i][k] - betaU_dD[i][k] * LambdabarU[k]  # Term 1

        # Step 6.b: Term 2 of \partial_t \bar{\Lambda}^i = \bar{\gamma}^{jk} (Term 2a + Term 2b + Term 2c)
        # Term 2a: \bar{\gamma}^{jk} \beta^i_{,kj}
        Term2aUDD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    Term2aUDD[i][j][k] += betaU_dDD[i][k][j]
        # Term 2b: \hat{\Gamma}^i_{mk,j} \beta^m + \hat{\Gamma}^i_{mk} \beta^m_{,j}
        #          + \hat{\Gamma}^i_{dj}\beta^d_{,k} - \hat{\Gamma}^d_{kj} \beta^i_{,d}
        Term2bUDD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for m in range(3):
                        Term2bUDD[i][j][k] += rfm.GammahatUDDdD[i][m][k][j] * betaU[m] \
                                              + rfm.GammahatUDD[i][m][k] * betaU_dD[m][j] \
                                              + rfm.GammahatUDD[i][m][j] * betaU_dD[m][k] \
                                              - rfm.GammahatUDD[m][k][j] * betaU_dD[i][m]
        # Term 2c: \hat{\Gamma}^i_{dj}\hat{\Gamma}^d_{mk} \beta^m - \hat{\Gamma}^d_{kj} \hat{\Gamma}^i_{md} \beta^m
        Term2cUDD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for m in range(3):
                        for d in range(3):
                            Term2cUDD[i][j][k] += (rfm.GammahatUDD[i][d][j] * rfm.GammahatUDD[d][m][k]
                                                   - rfm.GammahatUDD[d][k][j] * rfm.GammahatUDD[i][m][d]) * betaU[m]

        Lambdabar_rhsUpieceU = ixp.zerorank1()

        # Put it all together to get Term 2:
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.Lambdabar_rhsU[i] += gammabarUU[j][k] * (Term2aUDD[i][j][k] + Term2bUDD[i][j][k] + Term2cUDD[i][j][k])
                    Lambdabar_rhsUpieceU[i] += gammabarUU[j][k] * (
                                Term2aUDD[i][j][k] + Term2bUDD[i][j][k] + Term2cUDD[i][j][k])

        # Step 6.c: Term 3 of \partial_t \bar{\Lambda}^i:
        #    \frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j}
        DGammaU = Bq.DGammaU  # From Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
        for i in range(3):
            self.Lambdabar_rhsU[i] += sp.Rational(2, 3) * DGammaU[i] * Dbarbetacontraction  # Term 3

        # Step 6.d: Term 4 of \partial_t \bar{\Lambda}^i:
        #           \frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j}
        detgammabar_dDD = Bq.detgammabar_dDD  # From Bq.detgammabar_and_derivs()
        Dbarbetacontraction_dBarD = ixp.zerorank1()
        for k in range(3):
            for m in range(3):
                Dbarbetacontraction_dBarD[m] += betaU_dDD[k][k][m] + \
                                                (betaU_dD[k][m] * detgammabar_dD[k] +
                                                 betaU[k] * detgammabar_dDD[k][m]) / (2 * detgammabar) \
                                                - betaU[k] * detgammabar_dD[k] * detgammabar_dD[m] / (
                                                            2 * detgammabar * detgammabar)
        for i in range(3):
            for m in range(3):
                self.Lambdabar_rhsU[i] += sp.Rational(1, 3) * gammabarUU[i][m] * Dbarbetacontraction_dBarD[m]

        # Step 6.e: Term 5 of \partial_t \bar{\Lambda}^i:
        #           - 2 \bar{A}^{i j} (\partial_{j} \alpha - 6 \alpha \partial_{j} \phi)
        for i in range(3):
            for j in range(3):
                self.Lambdabar_rhsU[i] += -2 * AbarUU[i][j] * (alpha_dD[j] - 6 * alpha * phi_dD[j])

        # Step 6.f: Term 6 of \partial_t \bar{\Lambda}^i:
        #           2 \alpha \bar{A}^{j k} \Delta^{i}_{j k}
        DGammaUDD = Bq.DGammaUDD  # From RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.Lambdabar_rhsU[i] += 2 * alpha * AbarUU[j][k] * DGammaUDD[i][j][k]

        # Step 6.g: Term 7 of \partial_t \bar{\Lambda}^i:
        #           -\frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K
        trK_dD = ixp.declarerank1("trK_dD")
        for i in range(3):
            for j in range(3):
                self.Lambdabar_rhsU[i] += -sp.Rational(4, 3) * alpha * gammabarUU[i][j] * trK_dD[j]

        # Step 7: Add T4munu if enable_T4munu
        if enable_T4munu:
            if "T4UU00" not in gri.glb_gridfcs_dict:
                _ = gri.register_gridfunctions_for_single_rank2(
                    "T4UU",
                    symmetry="sym01",
                    dimension=4,
                    group="AUXEVOL",
                )
            (
                sourceterm_trK_rhs,
                sourceterm_Abar_rhsDD,
                sourceterm_Lambdabar_rhsU
            ) = T4munu.BSSN_RHSs_T4UU_source_terms(
                CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
            )
            self.trK_rhs += sourceterm_trK_rhs
            for i in range(3):
                self.Lambdabar_rhsU[i] += sourceterm_Lambdabar_rhsU[i]
                for j in range(3):
                    Abar_rhsDD[i][j] += sourceterm_Abar_rhsDD[i][j]


        # Step 8: Rescale the RHS quantities so that the evolved
        #         variables are smooth across coord singularities
        self.h_rhsDD = ixp.zerorank2()
        self.a_rhsDD = ixp.zerorank2()
        self.lambda_rhsU = ixp.zerorank1()
        for i in range(3):
            self.lambda_rhsU[i] = self.Lambdabar_rhsU[i] / rfm.ReU[i]
            for j in range(3):
                self.h_rhsDD[i][j] = gammabar_rhsDD[i][j] / rfm.ReDD[i][j]
                self.a_rhsDD[i][j] = Abar_rhsDD[i][j] / rfm.ReDD[i][j]
        # print(str(Abar_rhsDD[2][2]).replace("**","^").replace("_","").replace("xx","x").replace("sin(x2)","Sin[x2]").replace("sin(2*x2)","Sin[2*x2]").replace("cos(x2)","Cos[x2]").replace("detgbaroverdetghat","detg"))
        # print(str(Dbarbetacontraction).replace("**","^").replace("_","").replace("xx","x").replace("sin(x2)","Sin[x2]").replace("detgbaroverdetghat","detg"))
        # print(betaU_dD)
        # print(str(trK_rhs).replace("xx2","xx3").replace("xx1","xx2").replace("xx0","xx1").replace("**","^").replace("_","").replace("sin(xx2)","Sinx2").replace("xx","x").replace("sin(2*x2)","Sin2x2").replace("cos(x2)","Cosx2").replace("detgbaroverdetghat","detg"))
        # print(str(bet_rhsU[0]).replace("xx2","xx3").replace("xx1","xx2").replace("xx0","xx1").replace("**","^").replace("_","").replace("sin(xx2)","Sinx2").replace("xx","x").replace("sin(2*x2)","Sin2x2").replace("cos(x2)","Cosx2").replace("detgbaroverdetghat","detg"))

        # fmt: on

        self.BSSN_RHSs_varnames: List[str] = []
        self.BSSN_RHSs_exprs: List[sp.Expr] = []
        self.BSSN_RHSs_varnames += ["cf_rhs", "trK_rhs"]
        self.BSSN_RHSs_exprs += [self.cf_rhs, self.trK_rhs]
        for i in range(3):
            self.BSSN_RHSs_varnames += [f"lambda_rhsU{i}"]
            self.BSSN_RHSs_exprs += [self.lambda_rhsU[i]]
            for j in range(i, 3):
                self.BSSN_RHSs_varnames += [f"a_rhsDD{i}{j}", f"h_rhsDD{i}{j}"]
                self.BSSN_RHSs_exprs += [self.a_rhsDD[i][j], self.h_rhsDD[i][j]]
        # Sort the lists alphabetically by varname:
        sorted_list = sorted(zip(self.BSSN_RHSs_varnames, self.BSSN_RHSs_exprs))
        self.BSSN_RHSs_varnames, self.BSSN_RHSs_exprs = [
            list(t) for t in zip(*sorted_list)
        ]


class BSSNRHSs_dict(Dict[str, BSSNRHSs]):
    """Custom dictionary for storing BSSNRHSs objects."""

    def __getitem__(self, CoordSystem_in: str) -> BSSNRHSs:
        if CoordSystem_in not in self:
            enable_T4munu = par.parval_from_str("enable_T4munu")

            # In case [CoordSystem]_rfm_precompute is passed:
            CoordSystem = CoordSystem_in.replace("_rfm_precompute", "")
            print(
                f"Setting up BSSN_RHSs for CoordSystem = {CoordSystem}, enable_T4munu={enable_T4munu}."
            )
            self.__setitem__(
                CoordSystem, BSSNRHSs(CoordSystem, enable_rfm_precompute=False)
            )
            self.__setitem__(
                CoordSystem + "_rfm_precompute",
                BSSNRHSs(CoordSystem, enable_rfm_precompute=True),
            )
        return dict.__getitem__(self, CoordSystem_in)

    def __setitem__(self, CoordSystem: str, value: BSSNRHSs) -> None:
        dict.__setitem__(self, CoordSystem, value)

    def __delitem__(self, CoordSystem: str) -> None:
        dict.__delitem__(self, CoordSystem)


BSSN_RHSs = BSSNRHSs_dict()


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

    # enable T4munu to ensure a maximally comprehensive test.
    par.set_parval_from_str("enable_T4munu", True)

    for Coord in [
        "Spherical",
        "SinhSpherical_rfm_precompute",
        "Cartesian",
        "SinhCartesian",
        "SinhCylindrical",
        "SinhSymTP",
    ]:
        brhs = BSSN_RHSs[Coord]
        results_dict = ve.process_dictionary_of_expressions(
            brhs.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{Coord}",
            results_dict,
        )
