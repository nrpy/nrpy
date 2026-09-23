"""
Construct the physical spacetime Christoffel symbols from BSSN quantities.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.reference_metric as refmetric
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
from nrpy.equations.general_relativity.g4munu_conversions import (
    ADM_to_g4DD,
    ADM_to_g4UU,
)


class BSSN_to_g4Christoffel:
    r"""
    Construct the physical spacetime Christoffel symbols from BSSN quantities.

    This class reuses the existing BSSN-to-ADM reconstruction already present in
    NRPy, then assembles the missing time-derivative sector needed to build the
    physical four-metric derivatives and Christoffel symbols. The K_{ij}
    reconstruction in this class always projects out any residual conformal
    trace from \bar{A}_{ij} so that the physical trace remains trK.

    :param CoordSystem: Coordinate system to use.
    :param enable_rfm_precompute: Whether to enable reference-metric precomputation.
    """

    def __init__(
        self,
        CoordSystem: str = "Cartesian",
        enable_rfm_precompute: bool = False,
    ) -> None:
        r"""
        Initialize the BSSN_to_g4Christoffel class.

        :param CoordSystem: Coordinate system to use.
        :param enable_rfm_precompute: Whether to enable reference-metric precomputation.
        :raises ValueError: If the coordinate system is unsupported.
        :raises ValueError: If enable_rfm_precompute is not a bool.
        """
        # Step 1.a: Validate constructor inputs before triggering cached dictionary lookups.
        if CoordSystem not in refmetric.supported_CoordSystems:
            raise ValueError(
                f"Unsupported CoordSystem '{CoordSystem}'. "
                f"Supported coordinate systems: {refmetric.supported_CoordSystems}"
            )
        if not isinstance(enable_rfm_precompute, bool):
            raise ValueError(
                "enable_rfm_precompute must be a bool, "
                f"got {type(enable_rfm_precompute).__name__}"
            )

        # Step 1.b: Given the chosen coordinate system, set up the reference metric.
        coord_system_key = CoordSystem + (
            "_rfm_precompute" if enable_rfm_precompute else ""
        )
        rfm = refmetric.reference_metric[coord_system_key]

        # Step 1.c: Import all needed BSSN and ADM quantities already defined elsewhere in NRPy.
        Bq = BSSN_quantities[coord_system_key]
        BtoA = BSSN_to_ADM(
            CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
        )

        # Step 1.d: Define gauge quantities and their derivatives.
        self.alpha = Bq.alpha
        self.alpha_dD = ixp.declarerank1("alpha_dD")
        self.alpha_d0 = sp.symbols("alpha_d0", real=True)
        self.vetU_d0 = ixp.declarerank1("vetU_d0")

        # Step 1.e: Import the ADM and BSSN tensors needed for the 4-Christoffel construction.
        self.betaU = Bq.betaU
        self.betaU_dD = Bq.betaU_dD
        self.gammaDD = BtoA.gammaDD
        self.gammaDDdD = BtoA.gammaDDdD
        self.gammaUU = BtoA.gammaUU
        self.GammaUDD = BtoA.GammaUDD

        # Step 1.f: Reconstruct K_{ij}, always enforcing the trace-free projection
        #           so that the physical trace remains trK even if numerical error
        #           introduces a small conformal trace in \bar{A}_{ij}.
        self.KDD = ixp.zerorank2()
        self.trAbar = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                self.trAbar += Bq.gammabarUU[i][j] * Bq.AbarDD[i][j]

        self.AbarDD_for_KDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self.AbarDD_for_KDD[i][j] = (
                    Bq.AbarDD[i][j]
                    - sp.Rational(1, 3) * Bq.gammabarDD[i][j] * self.trAbar
                )
                self.KDD[i][j] = (
                    self.AbarDD_for_KDD[i][j] / Bq.exp_m4phi
                    + sp.Rational(1, 3) * self.gammaDD[i][j] * Bq.trK
                )

        # Step 2.a: Construct the covariant shift \beta_i = \gamma_{ij} \beta^j.
        self.betaD = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                self.betaD[i] += self.gammaDD[i][j] * self.betaU[j]

        # Step 2.b: Construct the spatial derivatives
        #           \beta_{i,k} = \gamma_{ij,k} \beta^j + \gamma_{ij} \beta^j_{,k}.
        self.betaDdD = ixp.zerorank2()
        for i in range(3):
            for k in range(3):
                for j in range(3):
                    self.betaDdD[i][k] += (
                        self.gammaDDdD[i][j][k] * self.betaU[j]
                        + self.gammaDD[i][j] * self.betaU_dD[j][k]
                    )

        # Step 3.a: Construct D_i \beta_j = \beta_{j,i} - \Gamma^k_{ij} \beta_k.
        self.DbetaDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                # D_i beta_j starts from partial_i beta_j = beta_j,i.
                self.DbetaDD[i][j] = self.betaDdD[j][i]
                for k in range(3):
                    self.DbetaDD[i][j] += -self.GammaUDD[k][i][j] * self.betaD[k]

        # Step 3.b: Construct
        #           \gamma_{ij,t} = D_i \beta_j + D_j \beta_i - 2 \alpha K_{ij}.
        self.gammaDDd0 = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self.gammaDDd0[i][j] = (
                    self.DbetaDD[i][j]
                    + self.DbetaDD[j][i]
                    - 2 * self.alpha * self.KDD[i][j]
                )

        # Step 4.a: Construct
        #           \beta^i_{,t} = \partial_t \mathrm{vet}^i \operatorname{ReU}[i].
        self.betaU_d0 = ixp.zerorank1()
        for i in range(3):
            self.betaU_d0[i] = self.vetU_d0[i] * rfm.ReU[i]

        # Step 4.b: Construct
        #           \beta_{i,t} = \gamma_{ij,t} \beta^j + \gamma_{ij} \beta^j_{,t}.
        self.betaDd0 = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                self.betaDd0[i] += (
                    self.gammaDDd0[i][j] * self.betaU[j]
                    + self.gammaDD[i][j] * self.betaU_d0[j]
                )

        # Step 5.a: Construct the spacetime metric g_{\mu\nu}
        #           from \gamma_{ij}, \beta^i, and \alpha.
        self.g4DD = ADM_to_g4DD(self.gammaDD, self.betaU, self.alpha)

        # Step 5.b: Construct the inverse spacetime metric g^{\mu\nu}.
        self.g4UU = ADM_to_g4UU(
            self.gammaDD, self.betaU, self.alpha, gammaUU=self.gammaUU
        )

        # Step 6: Assemble all first derivatives g_{\mu\nu,\rho} of the spacetime metric.
        self.g4DD_dD = ixp.zerorank3(dimension=4)

        # Step 6.a: Construct
        #           g_{00,k} = -2 \alpha \alpha_{,k}
        #                    + \gamma_{ij,k} \beta^i \beta^j
        #                    + 2 \beta_i \beta^i_{,k}.
        for k in range(3):
            self.g4DD_dD[0][0][k + 1] = -2 * self.alpha * self.alpha_dD[k]
            for i in range(3):
                self.g4DD_dD[0][0][k + 1] += 2 * self.betaD[i] * self.betaU_dD[i][k]

                for j in range(3):
                    # Step 6.c: Construct g_{ij,k} = \gamma_{ij,k}.
                    self.g4DD_dD[i + 1][j + 1][k + 1] = self.gammaDDdD[i][j][k]
                    self.g4DD_dD[0][0][k + 1] += (
                        self.gammaDDdD[i][j][k] * self.betaU[i] * self.betaU[j]
                    )

                # Step 6.b: Construct g_{0i,k} = \beta_{i,k}.
                self.g4DD_dD[0][i + 1][k + 1] = self.betaDdD[i][k]
                self.g4DD_dD[i + 1][0][k + 1] = self.betaDdD[i][k]

        # Step 6.d: Construct
        #           g_{00,t} = -2 \alpha \alpha_{,t}
        #                    + \gamma_{ij,t} \beta^i \beta^j
        #                    + 2 \beta_i \beta^i_{,t}.
        self.g4DD_dD[0][0][0] = -2 * self.alpha * self.alpha_d0
        for i in range(3):
            self.g4DD_dD[0][0][0] += 2 * self.betaD[i] * self.betaU_d0[i]

            for j in range(3):
                # Step 6.f: Construct g_{ij,t} = \gamma_{ij,t}.
                self.g4DD_dD[i + 1][j + 1][0] = self.gammaDDd0[i][j]
                self.g4DD_dD[0][0][0] += (
                    self.gammaDDd0[i][j] * self.betaU[i] * self.betaU[j]
                )

            # Step 6.e: Construct g_{0i,t} = \beta_{i,t}.
            self.g4DD_dD[0][i + 1][0] = self.betaDd0[i]
            self.g4DD_dD[i + 1][0][0] = self.betaDd0[i]

        # Step 6.g: Construct the static-spacetime metric derivatives by
        #           retaining all spatial derivatives and removing only the
        #           coordinate-time derivatives.
        self.g4DD_dD_static = ixp.zerorank3(dimension=4)
        for mu in range(4):
            for nu in range(4):
                for derivative_direction in range(4):
                    self.g4DD_dD_static[mu][nu][derivative_direction] = self.g4DD_dD[
                        mu
                    ][nu][derivative_direction]
                self.g4DD_dD_static[mu][nu][0] = sp.sympify(0)

        # Step 7.a: Construct ADM-reduced helper contractions for the
        #           four-Christoffel symbols.
        self.KUD = ixp.zerorank2()
        self.DalphaU = ixp.zerorank1()
        self.DbetaUD = ixp.zerorank2()
        self.KbetaD = ixp.zerorank1()
        self.Kbetabeta = sp.sympify(0)
        self.betaDalpha = sp.sympify(0)
        for i in range(3):
            self.betaDalpha += self.betaU[i] * self.alpha_dD[i]
            for j in range(3):
                self.KbetaD[i] += self.KDD[i][j] * self.betaU[j]
                self.Kbetabeta += self.KDD[i][j] * self.betaU[i] * self.betaU[j]
                self.DalphaU[i] += self.gammaUU[i][j] * self.alpha_dD[j]
                self.DbetaUD[i][j] = self.betaU_dD[i][j]
                for k in range(3):
                    self.KUD[i][j] += self.gammaUU[i][k] * self.KDD[k][j]
                    self.DbetaUD[i][j] += self.GammaUDD[i][j][k] * self.betaU[k]

        # Step 7.b: Construct the spacetime Christoffels using the ADM-reduced
        #           formulas, avoiding the much larger generic four-metric
        #           contraction while preserving the same physical connection.
        self.Gamma4UDD = ixp.zerorank3(dimension=4)
        self.Gamma4UDD[0][0][0] = (
            self.alpha_d0 + self.betaDalpha - self.Kbetabeta
        ) / self.alpha
        for i in range(3):
            self.Gamma4UDD[0][0][i + 1] = (
                self.alpha_dD[i] - self.KbetaD[i]
            ) / self.alpha
            self.Gamma4UDD[0][i + 1][0] = self.Gamma4UDD[0][0][i + 1]

            betaD_Db = sp.sympify(0)
            KUDbeta = sp.sympify(0)
            for j in range(3):
                betaD_Db += self.betaU[j] * self.DbetaUD[i][j]
                KUDbeta += self.KUD[i][j] * self.betaU[j]
            self.Gamma4UDD[i + 1][0][0] = (
                self.betaU_d0[i]
                + betaD_Db
                + self.alpha * self.DalphaU[i]
                - 2 * self.alpha * KUDbeta
                + self.betaU[i]
                / self.alpha
                * (self.Kbetabeta - self.alpha_d0 - self.betaDalpha)
            )

            for j in range(3):
                self.Gamma4UDD[0][i + 1][j + 1] = -self.KDD[i][j] / self.alpha

                self.Gamma4UDD[i + 1][0][j + 1] = (
                    self.DbetaUD[i][j]
                    - self.alpha * self.KUD[i][j]
                    + self.betaU[i] / self.alpha * (self.KbetaD[j] - self.alpha_dD[j])
                )
                self.Gamma4UDD[i + 1][j + 1][0] = self.Gamma4UDD[i + 1][0][j + 1]

                for k in range(3):
                    self.Gamma4UDD[i + 1][j + 1][k + 1] = (
                        self.GammaUDD[i][j][k]
                        + self.betaU[i] * self.KDD[j][k] / self.alpha
                    )

        # Step 7.c: Construct static-spacetime Christoffel symbols directly
        #           from the static metric derivatives. The dynamic symbols
        #           above intentionally retain their compact ADM expressions.
        self.Gamma4UDD_static = ixp.zerorank3(dimension=4)
        for alpha in range(4):
            for mu in range(4):
                for nu in range(mu, 4):
                    term = sp.sympify(0)
                    for beta in range(4):
                        term += (
                            sp.Rational(1, 2)
                            * self.g4UU[alpha][beta]
                            * (
                                self.g4DD_dD_static[beta][nu][mu]
                                + self.g4DD_dD_static[beta][mu][nu]
                                - self.g4DD_dD_static[mu][nu][beta]
                            )
                        )
                    self.Gamma4UDD_static[alpha][mu][nu] = term
                    self.Gamma4UDD_static[alpha][nu][mu] = term


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    import os

    import nrpy.validate_expressions.validate_expressions as ve

    for Coord in [
        "Spherical",
        "SinhSpherical",
        "SinhSpherical_rfm_precompute",
        "Cartesian",
        "SinhCartesian",
        "SinhCylindrical",
        "SinhSymTP",
    ]:
        enable_rfm_pre = "rfm_precompute" in Coord
        bssn_to_g4christoffel = BSSN_to_g4Christoffel(
            Coord.replace("_rfm_precompute", "", enable_rfm_pre),
            enable_rfm_precompute=enable_rfm_pre,
        )
        results_dict = {
            "g4DD": bssn_to_g4christoffel.g4DD,
            "g4UU": bssn_to_g4christoffel.g4UU,
            "g4DD_dD": bssn_to_g4christoffel.g4DD_dD,
            "g4DD_dD_static": bssn_to_g4christoffel.g4DD_dD_static,
            "gammaDDd0": bssn_to_g4christoffel.gammaDDd0,
            "Gamma4UDD": bssn_to_g4christoffel.Gamma4UDD,
            "Gamma4UDD_static": bssn_to_g4christoffel.Gamma4UDD_static,
        }
        results_dict = ve.process_dictionary_of_expressions(
            results_dict, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{Coord}",
            results_dict,
        )
