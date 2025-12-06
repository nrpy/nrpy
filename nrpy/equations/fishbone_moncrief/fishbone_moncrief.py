"""
Construct Fishbone-Moncrief initial data for GRMHD simulations.

Implements the Fishbone-Moncrief torus (Fishbone & Moncrief, 1976, ApJ 207)
around a spinning black hole, following:
- $l(r)$ from Eq. 3.8,
- specific enthalpy $h$ from Eq. 3.6 (with integration constant as described below Eq. 3.6),
- 4-velocity components using Eqs. 3.3, 2.12 (corrected typo), 2.13, and 2.16,
- BL metric and inverse using Eq. 2.1 with shorthands from Eq. 3.5,
- BL to Kerr-Schild transformation following the HARM implementation,
- Kerr-Schild metric using Cook (Living Reviews, 2000) Eq. 79 and Appendix of Etienne et al. (2017, arXiv:1704.00599),
- magnetic field seed with $B^i = nabla times A$ (Eq. 15 of Porth et al., 2015, arXiv:1501.07276),
- ADM quantities from Eq. 4.49 of Gourgoulhon (arXiv:0703035), consistent with B&S Eq. 2.122,
- comoving magnetic field using Eqs. 23, 24, 27, and 31 of Noble et al. (2006, astro-ph/0503420),
- Valencia velocity using Eq. 11 of Porth et al. (2015, arXiv:1501.07276).

Final data are mapped from spherical to the chosen coordinate system (default Cartesian)
via the NRPy reference-metric framework.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot** com
"""

from typing import List

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.params as par

thismodule = __name__


class FishboneMoncriefID:
    """Construct and store Fishbone-Moncrief initial data for GRMHD simulations."""

    def __init__(self) -> None:
        """
        Initialize and compute Fishbone-Moncrief initial data in spherical coordinates.

        Sets up ADM quantities (gammaDD, KDD, alpha, betaU), hydrodynamic quantities
        (rho_initial, Pressure_initial, LorentzFactor, Valencia3velocityU),
        and magnetic field quantities (BtildeU, smallb2).
        """
        # Register code parameters
        self.r_in = par.register_CodeParameter(
            "REAL", thismodule, "r_in", defaultvalue=6.0, commondata=True
        )
        self.r_at_max_density = par.register_CodeParameter(
            "REAL", thismodule, "r_at_max_density", defaultvalue=12.0, commondata=True
        )
        self.a = par.register_CodeParameter(
            "REAL", thismodule, "a", defaultvalue=0.9375, commondata=True
        )
        self.M = par.register_CodeParameter(
            "REAL", thismodule, "M", defaultvalue=1.0, commondata=True
        )
        self.kappa = par.register_CodeParameter(
            "REAL", thismodule, "kappa", defaultvalue=1.0e-3, commondata=True
        )
        self.Gamma = par.register_CodeParameter(
            "REAL", thismodule, "Gamma", defaultvalue=4.0 / 3.0, commondata=True
        )
        self.A_b = par.register_CodeParameter(
            "REAL", thismodule, "A_b", defaultvalue=1.0, commondata=True
        )

        # Define spherical coordinates consistent with other NRPy 2 modules
        self.r, self.th, self.ph = sp.symbols("r th ph", real=True)

        # Initialize output quantities
        self.gammaDD: List[List[sp.Expr]] = ixp.zerorank2()
        self.KDD: List[List[sp.Expr]] = ixp.zerorank2()
        self.alpha: sp.Expr = sp.sympify(0)
        self.betaU: List[sp.Expr] = ixp.zerorank1()
        self.BU: List[sp.Expr] = ixp.zerorank1()

        # Hydrodynamic quantities
        self.rho_initial: sp.Expr = sp.sympify(0)
        self.Pressure_initial: sp.Expr = sp.sympify(0)
        self.LorentzFactor: sp.Expr = sp.sympify(0)
        self.Valencia3velocityU: List[sp.Expr] = ixp.zerorank1()

        # Magnetic field quantities
        self.BtildeU: List[sp.Expr] = ixp.zerorank1()
        self.smallb2: sp.Expr = sp.sympify(0)

        # Compute all initial data
        self._compute_initial_data()

    def _calculate_l_at_r(self, r_val: sp.Expr) -> sp.Expr:
        """
        Compute the specific angular momentum $l(r)$ at radius $r$ using Eq. 3.8 of Fishbone and Moncrief (1976).

        :param r_val: Boyer-Lindquist radius.
        :return: $l(r)$.
        """
        M = self.M
        a = self.a
        l_val = sp.sqrt(M / r_val**3) * (
            r_val**4
            + r_val**2 * a**2
            - 2 * M * r_val * a**2
            - a * sp.sqrt(M * r_val) * (r_val**2 - a**2)
        )
        l_val /= r_val**2 - 3 * M * r_val + 2 * a * sp.sqrt(M * r_val)
        return l_val

    def _compute_initial_data(self) -> None:
        """Compute and store all initial-data expressions."""
        # Retrieve parameters
        M = self.M
        a = self.a
        r_in = self.r_in
        r_at_max_density = self.r_at_max_density
        kappa = self.kappa
        Gamma = self.Gamma
        A_b = self.A_b

        # Spherical coordinates
        r = self.r
        th = self.th
        ph = self.ph

        # Step 1: Angular momentum at r_at_max_density (co-rotating disk; positive root) [Eq. 3.8]
        l_val = self._calculate_l_at_r(r_at_max_density)

        # Step 2: Log enthalpy terms (Eq. 3.6)
        Delta = r**2 - 2 * M * r + a**2
        Sigma = r**2 + a**2 * sp.cos(th) ** 2
        Aexpr = (r**2 + a**2) ** 2 - Delta * a**2 * sp.sin(th) ** 2

        tmp3 = sp.sqrt(1 + 4 * l_val**2 * Sigma**2 * Delta / (Aexpr * sp.sin(th)) ** 2)
        ln_h = sp.Rational(1, 2) * sp.log((1 + tmp3) / (Sigma * Delta / Aexpr))
        ln_h -= sp.Rational(1, 2) * tmp3
        ln_h -= 2 * a * M * r * l_val / Aexpr

        # Radially independent integration constant part of ln h (text below Eq. 3.6)
        Delin = r_in**2 - 2 * M * r_in + a**2
        Sigin = r_in**2 + a**2 * sp.cos(sp.pi / 2) ** 2
        Ain = (r_in**2 + a**2) ** 2 - Delin * a**2 * sp.sin(sp.pi / 2) ** 2
        tmp3in = sp.sqrt(
            1 + 4 * l_val**2 * Sigin**2 * Delin / (Ain * sp.sin(sp.pi / 2)) ** 2
        )
        mln_h_in = -sp.Rational(1, 2) * sp.log((1 + tmp3in) / (Sigin * Delin / Ain))
        mln_h_in += sp.Rational(1, 2) * tmp3in
        mln_h_in += 2 * a * M * r_in * l_val / Ain

        hm1 = sp.exp(ln_h + mln_h_in) - 1

        # Step 3: Density and pressure (power via exp/log workaround)
        self.rho_initial = sp.exp(
            (1 / (Gamma - 1)) * sp.log(hm1 * (Gamma - 1) / (kappa * Gamma))
        )
        self.Pressure_initial = kappa * self.rho_initial**Gamma

        # Step 4: Covariant velocity components u_mu in BL coordinates
        # exp(-2 chi) uses Eq. 2.16 with exp(2 nu) and exp(2 psi) from Eq. 3.5
        exp2nu = Sigma * Delta / Aexpr
        exp2psi = Aexpr * sp.sin(th) ** 2 / Sigma
        expm2chi = exp2nu / exp2psi
        # u_(phi) from Eq. 3.3
        u_pphip = sp.sqrt((-1 + sp.sqrt(1 + 4 * l_val**2 * expm2chi)) / 2)
        # u_(t) from Eq. 2.13
        u_ptp = -sp.sqrt(1 + u_pphip**2)

        # u_phi from Eq. 2.12 (corrected typo) and u_t from Eq. 2.13
        uBL4D = ixp.zerorank1(dimension=4)
        uBL4D[3] = sp.sqrt(exp2psi) * u_pphip
        omega = 2 * a * M * r / Aexpr
        uBL4D[0] = u_ptp * sp.sqrt(exp2nu) - omega * uBL4D[3]

        # Step 5: Inverse BL metric g^{mu nu} from Eq. 2.1 with Eq. 3.5 shorthands
        gPhys4BLUU = ixp.zerorank2(dimension=4)
        gPhys4BLUU[0][0] = -Aexpr / (Delta * Sigma)
        gPhys4BLUU[0][3] = gPhys4BLUU[3][0] = -2 * a * M * r / (Delta * Sigma)
        gPhys4BLUU[3][3] = -4 * a**2 * M**2 * r**2 / (
            Delta * Aexpr * Sigma
        ) + Sigma**2 / (Aexpr * Sigma * sp.sin(th) ** 2)

        uBL4U = ixp.zerorank1(dimension=4)
        for i in range(4):
            for j in range(4):
                uBL4U[i] += gPhys4BLUU[i][j] * uBL4D[j]

        # Step 6: Transform u^mu to Kerr-Schild basis (HARM implementation)
        #         https://github.com/atchekho/harmpi/blob/master/init.c
        transformBLtoKS = ixp.zerorank2(dimension=4)
        for i in range(4):
            transformBLtoKS[i][i] = 1
        transformBLtoKS[0][1] = 2 * M * r / (r**2 - 2 * M * r + a**2)
        transformBLtoKS[3][1] = a / (r**2 - 2 * M * r + a**2)
        uKS4U = ixp.zerorank1(dimension=4)
        for i in range(4):
            for j in range(4):
                uKS4U[i] += transformBLtoKS[i][j] * uBL4U[j]

        # Step 7: Kerr-Schild metric components (Cook Eq. 79; Etienne et al. 2017 Appendix)
        rhoKS2 = r**2 + a**2 * sp.cos(th) ** 2
        alphaKS = 1 / sp.sqrt(1 + 2 * M * r / rhoKS2)
        betaKSU = ixp.zerorank1()
        betaKSU[0] = alphaKS**2 * 2 * M * r / rhoKS2
        gammaKSDD = ixp.zerorank2()
        gammaKSDD[0][0] = 1 + 2 * M * r / rhoKS2
        gammaKSDD[0][2] = gammaKSDD[2][0] = (
            -(1 + 2 * M * r / rhoKS2) * a * sp.sin(th) ** 2
        )
        gammaKSDD[1][1] = rhoKS2
        gammaKSDD[2][2] = (
            r**2 + a**2 + 2 * M * r / rhoKS2 * a**2 * sp.sin(th) ** 2
        ) * sp.sin(th) ** 2

        # Extrinsic curvature components (Cook Eq. 79 form; indices 0..2 are spatial)
        AA = a**2 * sp.cos(2 * th) + a**2 + 2 * r**2
        BB = AA + 4 * M * r
        DD = sp.sqrt(2 * M * r / (a**2 * sp.cos(th) ** 2 + r**2) + 1)
        KDD = ixp.zerorank2()
        KDD[0][0] = (
            DD
            * (AA + 2 * M * r)
            / (AA**2 * BB)
            * (4 * M * (a**2 * sp.cos(2 * th) + a**2 - 2 * r**2))
        )
        KDD[0][1] = KDD[1][0] = (
            DD / (AA * BB) * 8 * a**2 * M * r * sp.sin(th) * sp.cos(th)
        )
        KDD[0][2] = KDD[2][0] = (
            DD
            / AA**2
            * (-2 * a * M * sp.sin(th) ** 2 * (a**2 * sp.cos(2 * th) + a**2 - 2 * r**2))
        )
        KDD[1][1] = DD / BB * 4 * M * r**2
        KDD[1][2] = KDD[2][1] = (
            DD / (AA * BB) * (-8 * a**3 * M * r * sp.sin(th) ** 3 * sp.cos(th))
        )
        KDD[2][2] = (
            DD
            / (AA**2 * BB)
            * (
                2
                * M
                * r
                * sp.sin(th) ** 2
                * (
                    a**4 * (r - M) * sp.cos(4 * th)
                    + a**4 * (M + 3 * r)
                    + 4 * a**2 * r**2 * (2 * r - M)
                    + 4 * a**2 * r * sp.cos(2 * th) * (a**2 + r * (M + 2 * r))
                    + 8 * r**5
                )
            )
        )

        # Step 8: Build gPhys4UU from alpha, beta, and gamma (3-metric) [Eq. 4.49 of Gourgoulhon; B&S Eq. 2.122]
        gammaKSUU, _gammaKSDET = ixp.symm_matrix_inverter3x3(gammaKSDD)
        gPhys4UU = ixp.zerorank2(dimension=4)
        gPhys4UU[0][0] = -1 / alphaKS**2
        for i in range(1, 4):
            gPhys4UU[0][i] = gPhys4UU[i][0] = betaKSU[i - 1] / alphaKS**2
        for i in range(1, 4):
            for j in range(1, 4):
                gPhys4UU[i][j] = (
                    gammaKSUU[i - 1][j - 1]
                    - betaKSU[i - 1] * betaKSU[j - 1] / alphaKS**2
                )

        # Step 9: Vector potential and Btilde^i in spherical basis [Eq. 15 of Porth et al. 2015]
        A_3vecpotentialD = ixp.zerorank1()
        A_3vecpotentialD[2] = -A_b * self.rho_initial

        BtildeU = ixp.zerorank1()
        BtildeU[0] = sp.diff(A_3vecpotentialD[2], th) - sp.diff(A_3vecpotentialD[1], ph)
        BtildeU[1] = sp.diff(A_3vecpotentialD[0], ph) - sp.diff(A_3vecpotentialD[2], r)
        BtildeU[2] = sp.diff(A_3vecpotentialD[1], r) - sp.diff(A_3vecpotentialD[0], th)

        # Step 10: 3+1 quantities alpha, beta^i, gamma_{ij}, gamma^{ij}, gamma (Eq. 4.49 of Gourgoulhon)
        alpha = sp.sqrt(1 / (-gPhys4UU[0][0]))
        betaU = ixp.zerorank1()
        for i in range(3):
            betaU[i] = alpha**2 * gPhys4UU[0][i + 1]
        gammaUU = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                gammaUU[i][j] = gPhys4UU[i + 1][j + 1] + betaU[i] * betaU[j] / alpha**2

        gammaDD, igammaDET = ixp.symm_matrix_inverter3x3(gammaUU)
        gammaDET = 1 / igammaDET

        betaD = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                betaD[i] += gammaDD[i][j] * betaU[j]

        beta2 = sp.sympify(0)
        for i in range(3):
            beta2 += betaU[i] * betaD[i]

        gPhys4DD = ixp.zerorank2(dimension=4)
        gPhys4DD[0][0] = -(alpha**2) + beta2
        for i in range(3):
            gPhys4DD[0][i + 1] = gPhys4DD[i + 1][0] = betaD[i]
            for j in range(3):
                gPhys4DD[i + 1][j + 1] = gammaDD[i][j]

        # Step 11: Magnetic field in comoving frame and Lorentz factor [Eqs. 23, 24, 27, 31 of Noble et al. 2006]
        uKS4D = ixp.zerorank1(dimension=4)
        for i in range(4):
            for j in range(4):
                uKS4D[i] += gPhys4DD[i][j] * uKS4U[j]

        BU = ixp.zerorank1()
        for i in range(3):
            BU[i] = BtildeU[i] / sp.sqrt(gammaDET)

        BU0_u = sp.sympify(0)
        for i in range(3):
            BU0_u += uKS4D[i + 1] * BU[i] / alpha

        smallbU = ixp.zerorank1(dimension=4)
        smallbU[0] = BU0_u / sp.sqrt(4 * sp.pi)
        for i in range(3):
            smallbU[i + 1] = (BU[i] / alpha + BU0_u * uKS4U[i + 1]) / (
                sp.sqrt(4 * sp.pi) * uKS4U[0]
            )

        smallbD = ixp.zerorank1(dimension=4)
        for i in range(4):
            for j in range(4):
                smallbD[i] += gPhys4DD[i][j] * smallbU[j]

        smallb2 = sp.sympify(0)
        for i in range(4):
            smallb2 += smallbU[i] * smallbD[i]

        LorentzFactor = alpha * uKS4U[0]
        # Valencia 3-velocity from Eq. 11 of Porth et al. (2015)
        Valencia3velocityU = ixp.zerorank1()
        for i in range(3):
            Valencia3velocityU[i] = uKS4U[i + 1] / (alpha * uKS4U[0]) + betaU[i] / alpha

        # Step 12: Store results as class attributes
        self.alpha = alpha
        for i in range(3):
            self.betaU[i] = betaU[i]
            self.Valencia3velocityU[i] = Valencia3velocityU[i]
            self.BtildeU[i] = BtildeU[i]
            for j in range(3):
                self.gammaDD[i][j] = gammaDD[i][j]
                self.KDD[i][j] = KDD[i][j]
        self.LorentzFactor = LorentzFactor
        self.smallb2 = smallb2
        # BU is the time derivative of shift, set to zero for this initial data
        self.BU = ixp.zerorank1()


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

    ID = FishboneMoncriefID()
    results_dict = ve.process_dictionary_of_expressions(
        ID.__dict__, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        f"{os.path.splitext(os.path.basename(__file__))[0]}",
        results_dict,
    )
