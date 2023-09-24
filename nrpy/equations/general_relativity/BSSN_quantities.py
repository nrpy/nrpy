"""
This module provides functions that declare and define useful BSSN quantities

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from typing import Dict, List
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends

# Step 1: Import all needed modules from NRPy+:
import nrpy.params as par  # NRPy+: Parameter interface
import nrpy.grid as gri  # NRPy+: Functions having to do with numerical grids
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support

#  Declare/initialize parameters for this module
par.register_param(str, __name__, "EvolvedConformalFactor_cf", "W")
par.register_param(bool, __name__, "detgbarOverdetghat_equals_one", True)
par.register_param(bool, __name__, "LeaveRicciSymbolic", False)
par.register_param(bool, __name__, "enable_T4munu", False)


class BSSNQuantities:
    """
    This class handles the BSSN quantities involved in General Relativity simulations.
    It takes care of registering the necessary grid functions for these quantities.

    :param CoordSystem: (string) The coordinate system being used, defaults to "Cartesian"
    :param enable_rfm_precompute: (bool) Whether or not to enable precomputation for reference metric, defaults to False
    """

    def __init__(
        self, CoordSystem: str = "Cartesian", enable_rfm_precompute: bool = False
    ) -> None:
        """
        Set up all BSSN quantities & store to the class object.

        :param CoordSystem: (string) The coordinate system being used, defaults to "Cartesian"
        :param enable_rfm_precompute: (bool) Whether or not to enable precomputation for reference metric, defaults to False
        """
        # Step 2: Register all needed BSSN gridfunctions if needed.

        #   Check to see if this function has already been called.
        #   If so, do not register the gridfunctions again!
        LeaveRicciSymbolic = par.parval_from_str("LeaveRicciSymbolic")

        if any("hDD00" in gf.name for gf in gri.glb_gridfcs_dict.values()):
            self.hDD = ixp.declarerank2("hDD", symmetry="sym01")
            self.aDD = ixp.declarerank2("aDD", symmetry="sym01")
            self.lambdaU = ixp.declarerank1("lambdaU")
            self.vetU = ixp.declarerank1("vetU")
            self.betU = ixp.declarerank1("betU")
            self.trK, self.cf, self.alpha = sp.symbols(
                "trK cf alpha"
            )  # , real=True)  # no tensor gridfunction is real=True. Do this for consistency.
        else:
            # Step 2.a: Register indexed quantities, using ixp.register_... functions
            self.hDD = gri.register_gridfunctions_for_single_rank2(
                "hDD", symmetry="sym01", f_infinity=0.0, wavespeed=1.0
            )
            self.aDD = gri.register_gridfunctions_for_single_rank2(
                "aDD", symmetry="sym01", f_infinity=0.0, wavespeed=1.0
            )
            self.lambdaU = gri.register_gridfunctions_for_single_rank1(
                "lambdaU", f_infinity=0.0, wavespeed=1.0
            )
            self.vetU = gri.register_gridfunctions_for_single_rank1(
                "vetU", f_infinity=0.0, wavespeed=1.0
            )
            self.betU = gri.register_gridfunctions_for_single_rank1(
                "betU", f_infinity=0.0, wavespeed=1.0
            )

            # Step 2.b: Register scalar quantities, using gri.register_gridfunctions()
            self.trK, self.cf, self.alpha = gri.register_gridfunctions(
                ["trK", "cf", "alpha"],
                f_infinity=[0.0, 1.0, 1.0],
                wavespeed=[1.0, 1.0, sp.sqrt(2.0)],
            )

        if LeaveRicciSymbolic and not any(
            "RbarDD00" in gf.name for gf in gri.glb_gridfcs_dict.values()
        ):
            self.RbarDD = gri.register_gridfunctions_for_single_rank2(
                "RbarDD",
                symmetry="sym01",
                group="AUXEVOL",
                gf_array_name="auxevol_gfs",
            )
        else:
            self.RbarDD = ixp.declarerank2("RbarDD", symmetry="sym01")

        # fmt: off
        # Step 3: Define all basic conformal BSSN tensors
        #        gammabarDD,AbarDD,LambdabarU,betaU,BU
        #        in terms of BSSN gridfunctions.
        # Step 3.a: Defines gammabarDD, AbarDD, LambdabarU, betaU, BU
        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]

        # Step 3.a.i: gammabarDD and AbarDD:
        self.gammabarDD = ixp.zerorank2()
        self.AbarDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                # gammabar_{ij}  = h_{ij}*ReDD[i][j] + gammahat_{ij}
                self.gammabarDD[i][j] = self.hDD[i][j] * rfm.ReDD[i][j] + rfm.ghatDD[i][j]
                # Abar_{ij}      = a_{ij}*ReDD[i][j]
                self.AbarDD[i][j] = self.aDD[i][j] * rfm.ReDD[i][j]

        # Step 3.a.ii: LambdabarU, betaU, and BU:
        self.LambdabarU = ixp.zerorank1()
        self.betaU = ixp.zerorank1()
        self.BU = ixp.zerorank1()
        for i in range(3):
            self.LambdabarU[i] = self.lambdaU[i] * rfm.ReU[i]
            self.betaU[i] = self.vetU[i] * rfm.ReU[i]
            self.BU[i] = self.betU[i] * rfm.ReU[i]


        # Step 4: gammabarUU and spatial derivatives of gammabarDD,
        #         including GammabarUDD
        # Step 4.a: Defines gammabarUU, gammabarDD_dD, gammabarDD_dupD, gammabarDD_dDD, GammabarUDD
        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]

        # Step 4.a.i: gammabarUU:
        self.gammabarUU, dummydet = ixp.symm_matrix_inverter3x3(self.gammabarDD)

        # Step 4.b.i: gammabarDDdD[i][j][k]
        #               = \hat{\gamma}_{ij,k} + h_{ij,k} \text{ReDD[i][j]} + h_{ij} \text{ReDDdD[i][j][k]}.
        self.gammabarDD_dD = ixp.zerorank3()
        self.gammabarDD_dupD = ixp.zerorank3()
        hDD_dD = ixp.declarerank3("hDD_dD", symmetry="sym01")
        hDD_dupD = ixp.declarerank3("hDD_dupD", symmetry="sym01")
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.gammabarDD_dD[i][j][k] = rfm.ghatDDdD[i][j][k] + \
                                             hDD_dD[i][j][k] * rfm.ReDD[i][j] + self.hDD[i][j] * rfm.ReDDdD[i][j][k]

                    # Compute associated upwinded derivative, needed for the \bar{\gamma}_{ij} RHS
                    self.gammabarDD_dupD[i][j][k] = rfm.ghatDDdD[i][j][k] + \
                                               hDD_dupD[i][j][k] * rfm.ReDD[i][j] + self.hDD[i][j] * rfm.ReDDdD[i][j][k]

        # Step 4.b.ii: Compute gammabarDD_dDD in terms of the rescaled BSSN quantity hDD
        #      and its derivatives, as well as the reference metric and rescaling
        #      matrix, and its derivatives (expression given below):
        hDD_dDD = ixp.declarerank4("hDD_dDD", symmetry="sym01_sym23")
        self.gammabarDD_dDD = ixp.zerorank4()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        # gammabar_{ij,kl} = gammahat_{ij,kl}
                        #                  + h_{ij,kl} ReDD[i][j]
                        #                  + h_{ij,k} ReDDdD[i][j][l] + h_{ij,l} ReDDdD[i][j][k]
                        #                  + h_{ij} ReDDdDD[i][j][k][l]
                        self.gammabarDD_dDD[i][j][k][l] = rfm.ghatDDdDD[i][j][k][l]
                        self.gammabarDD_dDD[i][j][k][l] += hDD_dDD[i][j][k][l] * rfm.ReDD[i][j]
                        self.gammabarDD_dDD[i][j][k][l] += hDD_dD[i][j][k] * rfm.ReDDdD[i][j][l] + \
                                                      hDD_dD[i][j][l] * rfm.ReDDdD[i][j][k]
                        self.gammabarDD_dDD[i][j][k][l] += self.hDD[i][j] * rfm.ReDDdDD[i][j][k][l]

        # Step 4.b.iii: Define barred Christoffel symbol \bar{\Gamma}^{i}_{kl} = GammabarUDD[i][k][l] (see expression below)
        self.GammabarUDD = ixp.zerorank3()
        for i in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        # Gammabar^i_{kl} = 1/2 * gammabar^{im} ( gammabar_{mk,l} + gammabar_{ml,k} - gammabar_{kl,m}):
                        self.GammabarUDD[i][k][l] += sp.Rational(1, 2) * self.gammabarUU[i][m] * \
                                                (self.gammabarDD_dD[m][k][l] + self.gammabarDD_dD[m][l][k] - self.gammabarDD_dD[k][l][m])

        # Step 5: det(gammabarDD) and its derivatives
        # Step 5.a: Defines detgammabar, detgammabar_dD, detgammabar_dDD
        # Ignore return values of declare_BSSN_gridfunctions_if_not_declared_already() here, as they are unused

        detgbarOverdetghat = sp.sympify(1)
        detgbarOverdetghat_dD = ixp.zerorank1()
        detgbarOverdetghat_dDD = ixp.zerorank2()

        if not par.parval_from_str("detgbarOverdetghat_equals_one"):
            raise ValueError(
                "Error: detgbarOverdetghat_equals_one=False is not fully implemented yet."
            )

        ## Approach for implementing detgbarOverdetghat_equals_one=False:
        #     detgbarOverdetghat = gri.register_gridfunctions("AUX", ["detgbarOverdetghat"])
        #     detgbarOverdetghatInitial = gri.register_gridfunctions("AUX", ["detgbarOverdetghatInitial"])
        #     detgbarOverdetghat_dD = ixp.declarerank1("detgbarOverdetghat_dD")
        #     detgbarOverdetghat_dDD = ixp.declarerank2("detgbarOverdetghat_dDD", "sym01")

        # Step 5.b: Define detgammabar, detgammabar_dD, and detgammabar_dDD (needed for \partial_t \bar{\Lambda}^i below)
        self.detgammabar = detgbarOverdetghat * rfm.detgammahat

        self.detgammabar_dD = ixp.zerorank1()
        for i in range(3):
            self.detgammabar_dD[i] = (
                detgbarOverdetghat_dD[i] * rfm.detgammahat
                + detgbarOverdetghat * rfm.detgammahatdD[i]
            )

        self.detgammabar_dDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self.detgammabar_dDD[i][j] = (
                    detgbarOverdetghat_dDD[i][j] * rfm.detgammahat
                    + detgbarOverdetghat_dD[i] * rfm.detgammahatdD[j]
                    + detgbarOverdetghat_dD[j] * rfm.detgammahatdD[i]
                    + detgbarOverdetghat * rfm.detgammahatdDD[i][j]
                )

        # Step 6: Quantities related to conformal traceless
        #         extrinsic curvature AbarDD:
        #         AbarUU, AbarUD, and trAbar
        # Step 6.a: Defines AbarUU, AbarUD, trAbar, AbarDD_dD, AbarDD_dupD

        # Step 6.a.i: Compute Abar^{ij} in terms of Abar_{ij} and gammabar^{ij}
        self.AbarUU = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        # Abar^{ij} = gammabar^{ik} gammabar^{jl} Abar_{kl}
                        self.AbarUU[i][j] += (
                            self.gammabarUU[i][k]
                            * self.gammabarUU[j][l]
                            * self.AbarDD[k][l]
                        )

        # Step 6.a.ii: Compute Abar^i_j in terms of Abar_{ij} and gammabar^{ij}
        self.AbarUD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    # Abar^i_j = gammabar^{ik} Abar_{kj}
                    self.AbarUD[i][j] += self.gammabarUU[i][k] * self.AbarDD[k][j]

        # Step 6.a.iii: Compute Abar^k_k = trace of Abar:
        self.trAbar = sp.sympify(0)
        for k in range(3):
            for j in range(3):
                # Abar^k_k = gammabar^{kj} Abar_{jk}
                self.trAbar += self.gammabarUU[k][j] * self.AbarDD[j][k]

        # Step 6.a.iv: Compute Abar_{ij,k}
        self.AbarDD_dD = ixp.zerorank3()
        self.AbarDD_dupD = ixp.zerorank3()
        aDD_dD = ixp.declarerank3("aDD_dD", symmetry="sym01")
        aDD_dupD = ixp.declarerank3("aDD_dupD", symmetry="sym01")
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.AbarDD_dupD[i][j][k] = (
                        rfm.ReDDdD[i][j][k] * self.aDD[i][j]
                        + rfm.ReDD[i][j] * aDD_dupD[i][j][k]
                    )
                    self.AbarDD_dD[i][j][k] = (
                        rfm.ReDDdD[i][j][k] * self.aDD[i][j]
                        + rfm.ReDD[i][j] * aDD_dD[i][j][k]
                    )

        # Step 7: The conformal ("barred") Ricci tensor RbarDD
        #         and associated quantities
        # Step 7.a: Defines RbarDD, DGammaUDD, gammabarDD_dHatD, DGammaU

        # Step 7.a.i: Define \varepsilon_{ij} = epsDD[i][j]
        epsDD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                epsDD[i][j] = self.hDD[i][j] * rfm.ReDD[i][j]

        # Step 7.a.ii: Define epsDD_dD[i][j][k]
        hDD_dD = ixp.declarerank3("hDD_dD", symmetry="sym01")
        epsDD_dD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    epsDD_dD[i][j][k] = (
                        hDD_dD[i][j][k] * rfm.ReDD[i][j]
                        + self.hDD[i][j] * rfm.ReDDdD[i][j][k]
                    )

        # Step 7.a.iii: Define epsDD_dDD[i][j][k][l]
        hDD_dDD = ixp.declarerank4("hDD_dDD", symmetry="sym01_sym23")
        epsDD_dDD = ixp.zerorank4()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        epsDD_dDD[i][j][k][l] = (
                            hDD_dDD[i][j][k][l] * rfm.ReDD[i][j]
                            + hDD_dD[i][j][k] * rfm.ReDDdD[i][j][l]
                            + hDD_dD[i][j][l] * rfm.ReDDdD[i][j][k]
                            + self.hDD[i][j] * rfm.ReDDdDD[i][j][k][l]
                        )

        # Step 7.a.iv: DhatgammabarDDdD[i][j][l] = \bar{\gamma}_{ij;\hat{l}}
        # \bar{\gamma}_{ij;\hat{l}} = \varepsilon_{i j,l}
        #                           - \hat{\Gamma}^m_{i l} \varepsilon_{m j}
        #                           - \hat{\Gamma}^m_{j l} \varepsilon_{i m}
        self.gammabarDD_dHatD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                for l in range(3):
                    self.gammabarDD_dHatD[i][j][l] = epsDD_dD[i][j][l]
                    for m in range(3):
                        self.gammabarDD_dHatD[i][j][l] += (
                            -rfm.GammahatUDD[m][i][l] * epsDD[m][j]
                            - rfm.GammahatUDD[m][j][l] * epsDD[i][m]
                        )

        # Step 7.a.v: \bar{\gamma}_{ij;\hat{l},k} = DhatgammabarDD_dHatD_dD[i][j][l][k]:
        #        \bar{\gamma}_{ij;\hat{l},k} = \varepsilon_{ij,lk}
        #                                      - \hat{\Gamma}^m_{i l,k} \varepsilon_{m j}
        #                                      - \hat{\Gamma}^m_{i l} \varepsilon_{m j,k}
        #                                      - \hat{\Gamma}^m_{j l,k} \varepsilon_{i m}
        #                                      - \hat{\Gamma}^m_{j l} \varepsilon_{i m,k}
        gammabarDD_dHatD_dD = ixp.zerorank4()
        for i in range(3):
            for j in range(3):
                for l in range(3):
                    for k in range(3):
                        gammabarDD_dHatD_dD[i][j][l][k] = epsDD_dDD[i][j][l][k]
                        for m in range(3):
                            gammabarDD_dHatD_dD[i][j][l][k] += (
                                -rfm.GammahatUDDdD[m][i][l][k] * epsDD[m][j]
                                - rfm.GammahatUDD[m][i][l] * epsDD_dD[m][j][k]
                                - rfm.GammahatUDDdD[m][j][l][k] * epsDD[i][m]
                                - rfm.GammahatUDD[m][j][l] * epsDD_dD[i][m][k]
                            )

        # Step 7.a.vi: \bar{\gamma}_{ij;\hat{l}\hat{k}} = DhatgammabarDD_dHatDD[i][j][l][k]
        #          \bar{\gamma}_{ij;\hat{l}\hat{k}} = \partial_k \hat{D}_{l} \varepsilon_{i j}
        #                                           - \hat{\Gamma}^m_{lk} \left(\hat{D}_{m} \varepsilon_{i j}\right)
        #                                           - \hat{\Gamma}^m_{ik} \left(\hat{D}_{l} \varepsilon_{m j}\right)
        #                                           - \hat{\Gamma}^m_{jk} \left(\hat{D}_{l} \varepsilon_{i m}\right)
        gammabarDD_dHatDD = ixp.zerorank4()
        for i in range(3):
            for j in range(3):
                for l in range(3):
                    for k in range(3):
                        gammabarDD_dHatDD[i][j][l][k] = gammabarDD_dHatD_dD[i][j][l][k]
                        for m in range(3):
                            gammabarDD_dHatDD[i][j][l][k] += (
                                -rfm.GammahatUDD[m][l][k]
                                * self.gammabarDD_dHatD[i][j][m]
                                - rfm.GammahatUDD[m][i][k]
                                * self.gammabarDD_dHatD[m][j][l]
                                - rfm.GammahatUDD[m][j][k]
                                * self.gammabarDD_dHatD[i][m][l]
                            )

        # Step 7.b: Second term of RhatDD: compute \hat{D}_{j} \bar{\Lambda}^{k} = LambarU_dHatD[k][j]
        lambdaU_dD = ixp.declarerank2("lambdaU_dD", symmetry="nosym")
        LambarU_dHatD = ixp.zerorank2()
        for j in range(3):
            for k in range(3):
                LambarU_dHatD[k][j] = (
                    lambdaU_dD[k][j] * rfm.ReU[k] + self.lambdaU[k] * rfm.ReUdD[k][j]
                )
                for m in range(3):
                    LambarU_dHatD[k][j] += (
                        rfm.GammahatUDD[k][m][j] * self.lambdaU[m] * rfm.ReU[m]
                    )

        # Step 7.c: Conformal Ricci tensor, part 3: The \Delta^{k} \Delta_{(i j) k}
        #           + \bar{\gamma}^{k l}*(2 \Delta_{k(i}^{m} \Delta_{j) m l}
        #           + \Delta_{i k}^{m} \Delta_{m j l}) terms

        # Step 7.c.i: Define \Delta^i_{jk} = \bar{\Gamma}^i_{jk} - \hat{\Gamma}^i_{jk} = DGammaUDD[i][j][k]
        self.DGammaUDD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.DGammaUDD[i][j][k] = (
                        self.GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k]
                    )

        # Step 7.c.ii: Define \Delta^i = \bar{\gamma}^{jk} \Delta^i_{jk}
        self.DGammaU = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.DGammaU[i] += self.gammabarUU[j][k] * self.DGammaUDD[i][j][k]

        # If we wish to leave Ricci symbolic, there's no point continuing to define RbarDD symbolically:
        if LeaveRicciSymbolic:
            pass
        else:
            # Step 7.c.iii: Define \Delta_{ijk} = \bar{\gamma}_{im} \Delta^m_{jk}
            DGammaDDD = ixp.zerorank3()
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        for m in range(3):
                            DGammaDDD[i][j][k] += (
                                self.gammabarDD[i][m] * self.DGammaUDD[m][j][k]
                            )

            # Step 7.d: Summing the terms and defining \bar{R}_{ij}
            # Step 7.d.i: Add the first term to RbarDD:
            #         Rbar_{ij} += - \frac{1}{2} \bar{\gamma}^{k l} \hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j}
            self.RbarDD = ixp.zerorank2()
            RbarDDpiece = ixp.zerorank2()
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        for l in range(3):
                            self.RbarDD[i][j] += (
                                -sp.Rational(1, 2)
                                * self.gammabarUU[k][l]
                                * gammabarDD_dHatDD[i][j][l][k]
                            )
                            RbarDDpiece[i][j] += (
                                -sp.Rational(1, 2)
                                * self.gammabarUU[k][l]
                                * gammabarDD_dHatDD[i][j][l][k]
                            )

            # Step 7.d.ii: Add the second term to RbarDD:
            #         Rbar_{ij} += (1/2) * (gammabar_{ki} Lambar^k_{;\hat{j}} + gammabar_{kj} Lambar^k_{;\hat{i}})
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        self.RbarDD[i][j] += sp.Rational(1, 2) * (
                            self.gammabarDD[k][i] * LambarU_dHatD[k][j]
                            + self.gammabarDD[k][j] * LambarU_dHatD[k][i]
                        )

            # Step 7.d.iii: Add the remaining term to RbarDD:
            #      Rbar_{ij} += \Delta^{k} \Delta_{(i j) k} = 1/2 \Delta^{k} (\Delta_{i j k} + \Delta_{j i k})
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        self.RbarDD[i][j] += (
                            sp.Rational(1, 2)
                            * self.DGammaU[k]
                            * (DGammaDDD[i][j][k] + DGammaDDD[j][i][k])
                        )

            # Step 7.d.iv: Add the final term to RbarDD:
            #      Rbar_{ij} += \bar{\gamma}^{k l} (\Delta^{m}_{k i} \Delta_{j m l}
            #                   + \Delta^{m}_{k j} \Delta_{i m l}
            #                   + \Delta^{m}_{i k} \Delta_{m j l})
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        for l in range(3):
                            for m in range(3):
                                self.RbarDD[i][j] += self.gammabarUU[k][l] * (
                                    self.DGammaUDD[m][k][i] * DGammaDDD[j][m][l]
                                    + self.DGammaUDD[m][k][j] * DGammaDDD[i][m][l]
                                    + self.DGammaUDD[m][i][k] * DGammaDDD[m][j][l]
                                )

        # Step 8: The unrescaled shift vector betaU spatial derivatives:
        #         betaUdD & betaUdDD, written in terms of the
        #         rescaled shift vector vetU
        # Step 8.i: Defines betaU_dD, betaU_dupD, betaU_dDD

        # Step 8.ii: Compute the unrescaled shift vector beta^i = ReU[i]*vet^i
        vetU_dD = ixp.declarerank2("vetU_dD", symmetry="nosym")
        vetU_dupD = ixp.declarerank2(
            "vetU_dupD", symmetry="nosym"
        )  # Needed for upwinded \beta^i_{,j}
        vetU_dDD = ixp.declarerank3(
            "vetU_dDD", symmetry="sym12"
        )  # Needed for \beta^i_{,j}
        self.betaU_dD = ixp.zerorank2()
        self.betaU_dupD = ixp.zerorank2()  # Needed for, e.g., \beta^i RHS
        self.betaU_dDD = ixp.zerorank3()  # Needed for, e.g., \bar{\Lambda}^i RHS
        for i in range(3):
            for j in range(3):
                self.betaU_dD[i][j] = (
                    vetU_dD[i][j] * rfm.ReU[i] + self.vetU[i] * rfm.ReUdD[i][j]
                )
                self.betaU_dupD[i][j] = (
                    vetU_dupD[i][j] * rfm.ReU[i] + self.vetU[i] * rfm.ReUdD[i][j]
                )  # Needed for \beta^i RHS
                for k in range(3):
                    # Needed for, e.g., \bar{\Lambda}^i RHS:
                    self.betaU_dDD[i][j][k] = (
                        vetU_dDD[i][j][k] * rfm.ReU[i]
                        + vetU_dD[i][j] * rfm.ReUdD[i][k]
                        + vetU_dD[i][k] * rfm.ReUdD[i][j]
                        + self.vetU[i] * rfm.ReUdDD[i][j][k]
                    )

        # Step 9: Standard BSSN conformal factor phi,
        #         and its partial and covariant derivatives,
        #         all in terms of BSSN gridfunctions like cf
        # Step 9.a: Defines phi_dD, phi_dupD, phi_dDD, exp_m4phi, phi_dBarD, phi_dBarDD

        # Step 9.a.i: Define partial derivatives of \phi in terms of evolved quantity "cf":
        cf_dD = ixp.declarerank1("cf_dD")
        cf_dupD = ixp.declarerank1("cf_dupD")  # Needed for \partial_t \phi next.
        cf_dDD = ixp.declarerank2("cf_dDD", symmetry="sym01")
        self.phi_dD = ixp.zerorank1()
        self.phi_dupD = ixp.zerorank1()
        self.phi_dDD = ixp.zerorank2()

        # Step 9.a.ii: Assuming cf=phi, define exp_m4phi, phi_dD,
        #              phi_dupD (upwind finite-difference version of phi_dD), and phi_DD
        if par.parval_from_str("EvolvedConformalFactor_cf") == "phi":
            for i in range(3):
                self.phi_dD[i] = cf_dD[i]
                self.phi_dupD[i] = cf_dupD[i]
                for j in range(3):
                    self.phi_dDD[i][j] = cf_dDD[i][j]
            self.exp_m4phi = sp.exp(-4 * self.cf)

        # Step 9.a.iii: Assuming cf=W=e^{-2 phi}, define exp_m4phi, phi_dD,
        #               phi_dupD (upwind finite-difference version of phi_dD), and phi_DD
        if par.parval_from_str("EvolvedConformalFactor_cf") == "W":
            # \partial_i W = \partial_i (e^{-2 phi}) = -2 e^{-2 phi} \partial_i phi
            # -> \partial_i phi = -\partial_i cf / (2 cf)
            for i in range(3):
                self.phi_dD[i] = -cf_dD[i] / (2 * self.cf)
                self.phi_dupD[i] = -cf_dupD[i] / (2 * self.cf)
                for j in range(3):
                    # \partial_j \partial_i phi = - \partial_j [\partial_i cf / (2 cf)]
                    #                           = - cf_{,ij} / (2 cf) + \partial_i cf \partial_j cf / (2 cf^2)
                    self.phi_dDD[i][j] = (
                        -cf_dDD[i][j] + cf_dD[i] * cf_dD[j] / self.cf
                    ) / (2 * self.cf)
            self.exp_m4phi = self.cf * self.cf

        # Step 9.a.iv: Assuming cf=chi=e^{-4 phi}, define exp_m4phi, phi_dD,
        #              phi_dupD (upwind finite-difference version of phi_dD), and phi_DD
        if par.parval_from_str("EvolvedConformalFactor_cf") == "chi":
            # \partial_i chi = \partial_i (e^{-4 phi}) = -4 e^{-4 phi} \partial_i phi
            # -> \partial_i phi = -\partial_i cf / (4 cf)
            for i in range(3):
                self.phi_dD[i] = -cf_dD[i] / (4 * self.cf)
                self.phi_dupD[i] = -cf_dupD[i] / (4 * self.cf)
                for j in range(3):
                    # \partial_j \partial_i phi = - \partial_j [\partial_i cf / (4 cf)]
                    #                           = - cf_{,ij} / (4 cf) + \partial_i cf \partial_j cf / (4 cf^2)
                    self.phi_dDD[i][j] = (
                        -cf_dDD[i][j] + cf_dD[i] * cf_dD[j] / self.cf
                    ) / (4 * self.cf)
            self.exp_m4phi = self.cf

        # Step 9.a.v: Error out if unsupported EvolvedConformalFactor_cf choice is made:
        cf_choice = par.parval_from_str("EvolvedConformalFactor_cf")
        if cf_choice not in ("phi", "W", "chi"):
            raise ValueError(
                "EvolvedConformalFactor_cf == "
                + par.parval_from_str("EvolvedConformalFactor_cf")
                + " unsupported!"
            )

        # Step 9.b: Define phi_dBarD = phi_dD (since phi is a scalar) and phi_dBarDD (covariant derivative)
        #          \bar{D}_i \bar{D}_j \phi = \phi_{;\bar{i}\bar{j}} = \bar{D}_i \phi_{,j}
        #                                   = \phi_{,ij} - \bar{\Gamma}^k_{ij} \phi_{,k}
        self.phi_dBarD = self.phi_dD
        self.phi_dBarDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self.phi_dBarDD[i][j] = self.phi_dDD[i][j]
                for k in range(3):
                    self.phi_dBarDD[i][j] += -self.GammabarUDD[k][i][j] * self.phi_dD[k]

        self.Ricci_varnames: List[str] = []
        self.Ricci_exprs: List[sp.Expr] = []
        for i in range(3):
            for j in range(i, 3):
                self.Ricci_varnames += [f"RbarDD{i}{j}"]
                self.Ricci_exprs += [self.RbarDD[i][j]]
        # Sort the lists alphabetically by varname:
        sorted_list = sorted(zip(self.Ricci_varnames, self.Ricci_exprs))
        self.Ricci_varnames, self.Ricci_exprs = [list(t) for t in zip(*sorted_list)]

    # fmt: on


class BSSNQuantities_dict(Dict[str, BSSNQuantities]):
    """Custom dictionary for storing BSSNQuantities objects."""

    def __getitem__(self, CoordSystem_in: str) -> BSSNQuantities:
        if CoordSystem_in not in self:
            # In case [CoordSystem]_rfm_precompute is passed:
            CoordSystem = CoordSystem_in.replace("_rfm_precompute", "")
            enable_T4munu = par.parval_from_str("enable_T4munu")
            print(
                f"Setting up BSSN_quantities for CoordSystem = {CoordSystem}, enable_T4munu={enable_T4munu}."
            )
            self.__setitem__(
                CoordSystem, BSSNQuantities(CoordSystem, enable_rfm_precompute=False)
            )
            self.__setitem__(
                CoordSystem + "_rfm_precompute",
                BSSNQuantities(CoordSystem, enable_rfm_precompute=True),
            )
        return dict.__getitem__(self, CoordSystem_in)

    def __setitem__(self, CoordSystem: str, value: BSSNQuantities) -> None:
        dict.__setitem__(self, CoordSystem, value)

    def __delitem__(self, CoordSystem: str) -> None:
        dict.__delitem__(self, CoordSystem)


BSSN_quantities = BSSNQuantities_dict()

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

    # enable T4munu to ensure a maximally comprehensive test (has no impact here).
    par.set_parval_from_str("enable_T4munu", True)

    for Coord in [
        "Spherical",
        "SinhSpherical",
        "SinhSpherical_rfm_precompute",
        "Cartesian",
        "SinhCartesian",
        "SinhCylindrical",
        "SinhSymTP",
    ]:
        bq = BSSN_quantities[Coord]
        results_dict = ve.process_dictionary_of_expressions(
            bq.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{Coord}",
            results_dict,
        )
