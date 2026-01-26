"""
Construct all expressions required for the general relativistic hydrodynamics equations.

These equations are used for numerical relativity simulations where the dynamics of matter are coupled to those of the spacetime.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from typing import List

# Step 1.a: import all needed modules from NRPy:
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy depends

import nrpy.indexedexp as ixp  # NRPy: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.reference_metric as refmetric  # NRPy: Reference metric support
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
from nrpy.equations.general_relativity.g4munu_conversions import (
    ADM_to_g4DD,
    ADM_to_g4UU,
)


class GRHD_Equations:
    """Construct and store expressions for the GRHD equations."""

    def __init__(
        self,
        CoordSystem: str = "Cartesian",
        enable_rfm_precompute: bool = False,
    ) -> None:
        """
        Initialize and set up all GRHD quantities, storing them within the class object.

        :param enable_rfm_precompute: Whether to enable reference-metric
            precomputation, defaults to False.
        :param CoordSystem: The coordinate system being used, defaults
            to "Cartesian".

        """
        # Step 1.c: Given the chosen coordinate system, set up
        #           corresponding reference metric and needed
        #           reference metric quantities
        # The following function call sets up the reference metric
        #    and related quantities, including rescaling matrices ReDD,
        #    ReU, and hatted quantities.
        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]
        self.ReU = rfm.ReU
        self.GammahatUDD = rfm.GammahatUDD

        # Step 1.e: Import all basic (unrescaled) BSSN scalars & tensors
        self.Bq = BSSN_quantities[
            CoordSystem + ("_rfm_precompute" if enable_rfm_precompute else "")
        ]

        # First define hydrodynamical quantities
        self.rescaledvU = ixp.declarerank1("rescaledvU", dimension=3)
        self.u4Ut, self.rho_b, self.P, self.h, self.Ye, self.S = sp.symbols(
            "u4Ut rhob P h Ye S", real=True
        )

        self.VU = ixp.zerorank1()
        self.u4U = ixp.zerorank1(dimension=4)
        self.u4U[0] = self.u4Ut
        for i in range(3):
            self.VU[i] = self.rescaledvU[i] * self.ReU[i]
            self.u4U[i + 1] = self.VU[i] * self.u4U[0]

        self.betaU = self.Bq.betaU
        self.betaU_dD = self.Bq.betaU_dD
        self.alpha = self.Bq.alpha
        self.alpha_dD = ixp.declarerank1("alpha_dD", dimension=3)

        # ADM in terms of BSSN
        self.AitoB = BSSN_to_ADM(
            CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
        )

        self.gammaDD = self.AitoB.gammaDD
        self.e6phi = (self.Bq.exp_m4phi ** (0.5)) ** (-3.0)
        self.KDD = self.AitoB.KDD

        for i in range(3):
            for j in range(3):
                self.gammaDD[i][j] = sp.together(self.gammaDD[i][j])

        # Declare class variables that will be defined later
        self.T4UU: List[List[sp.Expr]]
        self.rescaledT4UU: List[List[sp.Expr]]

        self.T4UD: List[List[sp.Expr]]
        self.rescaledT4UD: List[List[sp.Expr]]

        self.rho_star: sp.Expr
        self.Ye_star: sp.Expr
        self.S_star: sp.Expr
        self.tau_tilde: sp.Expr
        self.S_tildeD: List[sp.Expr]
        self.rescaledS_tildeD: List[sp.Expr]

        self.VU_from_u4U: List[sp.Expr]

        self.rho_star_fluxU: List[sp.Expr]
        self.rescaled_rho_star_fluxU: List[sp.Expr]

        self.Ye_star_fluxU: List[sp.Expr]
        self.rescaled_Ye_star_fluxU: List[sp.Expr]

        self.S_star_fluxU: List[sp.Expr]
        self.rescaled_S_star_fluxU: List[sp.Expr]

        self.tau_tilde_fluxU: List[sp.Expr]
        self.rescaled_tau_tilde_fluxU: List[sp.Expr]

        self.S_tilde_fluxUD: List[List[sp.Expr]]
        self.rescaled_S_tilde_fluxUD: List[List[sp.Expr]]

        self.tau_source_term: sp.Expr
        self.rho_star_connection_term: sp.Expr
        self.Ye_star_connection_term: sp.Expr
        self.S_star_connection_term: sp.Expr
        self.tau_connection_term: sp.Expr

        self.S_tilde_connection_termsD: List[sp.Expr]

        self.S_tilde_source_termD: List[sp.Expr]

    def compute_T4UU(self) -> None:
        """Define T^{mu nu} (a 4-dimensional tensor)."""
        gammaDD = self.gammaDD
        betaU = self.betaU
        alpha = self.alpha
        rho_b = self.rho_b
        P = self.P
        h = self.h
        u4U = self.u4U

        # define g^{mu nu} in terms of the ADM quantities:
        g4UU_tmp = ADM_to_g4UU(gammaDD, betaU, alpha)

        g4UU = ixp.zerorank2(dimension=4)

        for mu in range(4):
            for nu in range(4):
                g4UU[mu][nu] = sp.together(g4UU_tmp[mu][nu])

        # compute T^{mu nu}
        self.T4UU = ixp.zerorank2(dimension=4)
        for mu in range(4):
            for nu in range(4):
                self.T4UU[mu][nu] = sp.together(
                    (rho_b * h) * u4U[mu] * u4U[nu] + P * g4UU[mu][nu]
                )

        self.rescaledT4UU = ixp.zerorank2(dimension=4)
        self.rescaledT4UU[0][0] = self.T4UU[0][0]
        for mu in range(3):
            self.rescaledT4UU[0][mu + 1] = sp.together(
                self.T4UU[0][mu + 1] / (self.ReU[mu])
            )
            self.rescaledT4UU[mu + 1][0] = self.rescaledT4UU[0][mu + 1]
        for mu in range(3):
            for nu in range(3):
                self.rescaledT4UU[mu + 1][nu + 1] = sp.together(
                    self.T4UU[mu + 1][nu + 1] / (self.ReU[mu] * self.ReU[nu])
                )

    def compute_T4UD(self) -> None:
        """Define T^{mu}_{nu} (a 4-dimensional tensor)."""
        gammaDD = self.gammaDD
        betaU = self.betaU
        alpha = self.alpha

        # compute T^mu_nu = T^{mu delta} g_{delta nu}, needed for S_tilde flux.
        # we'll need g_{alpha nu} in terms of ADM quantities:
        g4DD_tmp = ADM_to_g4DD(gammaDD, betaU, alpha)

        g4DD = ixp.zerorank2(dimension=4)

        for mu in range(4):
            for nu in range(4):
                g4DD[mu][nu] = sp.together(g4DD_tmp[mu][nu])

        self.T4UD = ixp.zerorank2(dimension=4)
        for mu in range(4):
            for nu in range(4):
                for delta in range(4):
                    self.T4UD[mu][nu] += sp.together(
                        self.T4UU[mu][delta] * g4DD[delta][nu]
                    )

        self.rescaledT4UD = ixp.zerorank2(dimension=4)
        self.rescaledT4UD[0][0] = self.T4UD[0][0]
        for mu in range(3):
            self.rescaledT4UD[0][mu + 1] = sp.together(
                self.T4UD[0][mu + 1] * (self.ReU[mu])
            )
            self.rescaledT4UD[mu + 1][0] = sp.together(
                self.T4UD[mu + 1][0] / (self.ReU[mu])
            )
        for mu in range(3):
            for nu in range(3):
                self.rescaledT4UD[mu + 1][nu + 1] = sp.together(
                    self.T4UD[mu + 1][nu + 1] * (self.ReU[nu] / self.ReU[mu])
                )

    def compute_rho_star(self) -> None:
        """Compute densitized conserved density."""
        alpha = self.alpha
        e6phi = self.e6phi
        rho_b = self.rho_b
        u4U = self.u4U

        # Compute rho_star:
        self.rho_star = alpha * e6phi * rho_b * u4U[0]

    def compute_Ye_star(self) -> None:
        """Compute densitized conserved electron fraction."""
        alpha = self.alpha
        e6phi = self.e6phi
        rho_b = self.rho_b
        u4U = self.u4U
        Ye = self.Ye

        # Compute Ye_star:
        self.Ye_star = Ye * alpha * e6phi * rho_b * u4U[0]

    def compute_S_star(self) -> None:
        """Compute densitized conserved entropy."""
        alpha = self.alpha
        e6phi = self.e6phi
        S = self.S
        u4U = self.u4U

        # Compute S_star:
        self.S_star = alpha * e6phi * S * u4U[0]

    def compute_tau_tilde(self) -> None:
        """Compute densitized conserved energy."""
        alpha = self.alpha
        e6phi = self.e6phi
        T4UU = self.T4UU
        rho_star = self.rho_star

        self.tau_tilde = alpha**2 * e6phi * T4UU[0][0] - rho_star

    def compute_S_tildeD(self) -> None:
        """Compute densitized conserved momentum."""
        alpha = self.alpha
        e6phi = self.e6phi
        T4UD = self.T4UD

        self.S_tildeD = ixp.zerorank1(dimension=3)
        for i in range(3):
            self.S_tildeD[i] = alpha * e6phi * T4UD[0][i + 1]

        self.rescaledS_tildeD = ixp.zerorank1(dimension=3)
        for i in range(3):
            self.rescaledS_tildeD[i] = alpha * e6phi * self.rescaledT4UD[0][i + 1]

    # Define the fluxes for the GRHD equations
    # vU from u4U may be needed for computing rho_star_flux from u4U
    def compute_vU_from_u4U__no_speed_limit(self) -> None:
        """Compute fluid three velocity from the four-velocity: v^i = u^i/u^0."""
        u4U = self.u4U

        self.VU_from_u4U = ixp.zerorank1(dimension=3)
        for j in range(3):
            self.VU_from_u4U[j] = u4U[j + 1] / u4U[0]

    def compute_rho_star_fluxU(self) -> None:
        """Density flux term."""
        VU = self.VU
        rho_star = self.rho_star

        self.rho_star_fluxU = ixp.zerorank1(dimension=3)
        self.rescaled_rho_star_fluxU = ixp.zerorank1(dimension=3)
        for j in range(3):
            self.rho_star_fluxU[j] = rho_star * VU[j]
            self.rescaled_rho_star_fluxU[j] = sp.together(
                self.rho_star_fluxU[j] / self.ReU[j]
            )

    def compute_Ye_star_fluxU(self) -> None:
        """Electron fraction flux term."""
        VU = self.VU
        Ye_star = self.Ye_star

        self.Ye_star_fluxU = ixp.zerorank1(dimension=3)
        self.rescaled_Ye_star_fluxU = ixp.zerorank1(dimension=3)
        for j in range(3):
            self.Ye_star_fluxU[j] = Ye_star * VU[j]
            self.rescaled_Ye_star_fluxU[j] = sp.together(
                self.Ye_star_fluxU[j] / self.ReU[j]
            )

    def compute_S_star_fluxU(self) -> None:
        """Entropy flux term."""
        VU = self.VU
        S_star = self.S_star

        self.S_star_fluxU = ixp.zerorank1(dimension=3)
        self.rescaled_S_star_fluxU = ixp.zerorank1(dimension=3)
        for j in range(3):
            self.S_star_fluxU[j] = S_star * VU[j]
            self.rescaled_S_star_fluxU[j] = sp.together(
                self.S_star_fluxU[j] / self.ReU[j]
            )

    def compute_tau_tilde_fluxU(self) -> None:
        """Energy flux term."""
        alpha = self.alpha
        e6phi = self.e6phi
        VU = self.VU
        T4UU = self.T4UU
        rho_star = self.rho_star

        self.tau_tilde_fluxU = ixp.zerorank1(dimension=3)
        self.rescaled_tau_tilde_fluxU = ixp.zerorank1(dimension=3)
        for j in range(3):
            self.tau_tilde_fluxU[j] = (
                alpha**2 * e6phi * T4UU[0][j + 1] - rho_star * VU[j]
            )
            self.rescaled_tau_tilde_fluxU[j] = sp.together(
                alpha**2 * e6phi * self.rescaledT4UU[0][j + 1]
                - rho_star * VU[j] / self.ReU[j]
            )

    def compute_S_tilde_fluxUD(self) -> None:
        """Momentum flux term."""
        alpha = self.alpha
        e6phi = self.e6phi
        T4UD = self.T4UD

        self.S_tilde_fluxUD = ixp.zerorank2(dimension=3)
        self.rescaled_S_tilde_fluxUD = ixp.zerorank2(dimension=3)
        for j in range(3):
            for i in range(3):
                self.S_tilde_fluxUD[j][i] = alpha * e6phi * T4UD[j + 1][i + 1]
                self.rescaled_S_tilde_fluxUD[j][i] = (
                    alpha * e6phi * self.rescaledT4UD[j + 1][i + 1]
                )

    def compute_tau_source_term(self) -> None:
        """Source terms for energy equation."""
        KDD = self.KDD
        betaU = self.betaU
        alpha = self.alpha
        e6phi = self.e6phi
        alpha_dD = self.alpha_dD
        T4UU = self.T4UU

        self.tau_source_term = sp.sympify(0)
        # Term 1:
        for i in range(3):
            for j in range(3):
                self.tau_source_term += (
                    T4UU[0][0] * betaU[i] * betaU[j]
                    + 2 * T4UU[0][i + 1] * betaU[j]
                    + T4UU[i + 1][j + 1]
                ) * KDD[i][j]
        # Term 2:
        for i in range(3):
            self.tau_source_term += (
                -(T4UU[0][0] * betaU[i] + T4UU[0][i + 1]) * alpha_dD[i]
            )
        # Term 3:
        self.tau_source_term *= alpha * e6phi

    def compute_all_connection_terms(self) -> None:
        """Source terms from connection coefficients for density, electron fraction, entropy, and energy equations."""
        self.rho_star_connection_term = sp.sympify(0)
        self.Ye_star_connection_term = sp.sympify(0)
        self.S_star_connection_term = sp.sympify(0)
        self.tau_connection_term = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                self.rho_star_connection_term += (
                    self.GammahatUDD[i][i][j] * self.rho_star_fluxU[j]
                )
                self.Ye_star_connection_term += (
                    self.GammahatUDD[i][i][j] * self.Ye_star_fluxU[j]
                )
                self.S_star_connection_term += (
                    self.GammahatUDD[i][i][j] * self.S_star_fluxU[j]
                )
                self.tau_connection_term += (
                    self.GammahatUDD[i][i][j] * self.tau_tilde_fluxU[j]
                )

    def compute_S_tilde_connection_termsD(self) -> None:
        """Source terms from connection coefficients for momentum equation."""
        self.S_tilde_connection_termsD = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                for l in range(3):
                    self.S_tilde_connection_termsD[i] += (
                        self.GammahatUDD[j][j][l] * self.S_tilde_fluxUD[l][i]
                        - self.GammahatUDD[l][j][i] * self.S_tilde_fluxUD[j][l]
                    )

    def compute_S_tilde_source_termD(self) -> None:
        """Source terms for momentum equation."""
        alpha = self.alpha
        alpha_dD = self.alpha_dD
        betaU = self.betaU
        betaU_dD = self.betaU_dD
        e6phi = self.e6phi
        T4UU = self.T4UU
        T4UD = self.T4UD

        first_termD = ixp.zerorank1()
        second_termD = ixp.zerorank1()
        third_termD = ixp.zerorank1()
        self.S_tilde_source_termD = ixp.zerorank1()

        covhatdD_gammaDD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    covhatdD_gammaDD[i][j][k] += (self.Bq.exp_m4phi ** (-1)) * (
                        sp.sympify(4) * self.Bq.gammabarDD[j][k] * self.Bq.phi_dD[i]
                        + self.Bq.gammabarDD_dD[j][k][i]
                    )
                    for l in range(3):
                        covhatdD_gammaDD[i][j][k] -= (self.Bq.exp_m4phi ** (-1)) * (
                            self.GammahatUDD[l][i][j] * self.Bq.gammabarDD[l][k]
                            + self.GammahatUDD[l][i][k] * self.Bq.gammabarDD[j][l]
                        )

        for i in range(3):
            first_termD[i] -= T4UU[0][0] * alpha * alpha_dD[i]
            for j in range(3):
                second_termD[i] += T4UD[0][j + 1] * (betaU_dD[j][i])
                for k in range(3):
                    second_termD[i] += (
                        T4UD[0][j + 1] * self.GammahatUDD[j][i][k] * betaU[k]
                    )
                    third_termD[i] += (
                        sp.sympify(1.0 / 2.0)
                        * covhatdD_gammaDD[i][j][k]
                        * (
                            T4UU[0][0] * betaU[j] * betaU[k]
                            + sp.sympify(2) * T4UU[0][j + 1] * betaU[k]
                            + T4UU[j + 1][k + 1]
                        )
                    )

        for i in range(3):
            self.S_tilde_source_termD[i] += (
                alpha * e6phi * (first_termD[i] + second_termD[i] + third_termD[i])
            )

    def construct_all_equations(self) -> None:
        """Run all class functions, in the appropriate order, to define all quantities required for the GRHD evolution equations."""
        self.compute_T4UU()
        self.compute_T4UD()

        self.compute_rho_star()
        self.compute_rho_star_fluxU()

        self.compute_Ye_star()
        self.compute_Ye_star_fluxU()

        self.compute_S_star()
        self.compute_S_star_fluxU()

        self.compute_tau_tilde()
        self.compute_tau_tilde_fluxU()

        self.compute_S_tildeD()
        self.compute_S_tilde_fluxUD()

        self.compute_tau_source_term()

        self.compute_all_connection_terms()

        self.compute_S_tilde_source_termD()
        self.compute_S_tilde_connection_termsD()

        print("Finished setting up GRHD equations.")


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
        "Cartesian",
        "SinhCartesian",
        "SinhCylindrical",
        "SinhSymTP",
    ]:
        enable_rfm_pre = "rfm_precompute" in Coord
        grhd_eqs = GRHD_Equations(Coord.replace("_rfm_precompute", "", enable_rfm_pre))
        grhd_eqs.construct_all_equations()
        results_dict = ve.process_dictionary_of_expressions(
            grhd_eqs.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{Coord}",
            results_dict,
        )
