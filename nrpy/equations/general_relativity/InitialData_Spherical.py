# nrpy/equations/general_relativity/InitialData_Spherical.py
"""
Set up initial data for solving Einstein's equations of general relativity, for data most naturally specified in Spherical coordinates.

Outputs ADM quantities.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List, Tuple, cast

import sympy as sp

# NRPy: This module depends on the parameter EvolvedConformalFactor_cf,
#        which is defined in BSSN.BSSN_quantities
import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN


class InitialData_Spherical:
    """Construct and store Spherical initial data for Einstein's equations of general relativity, as ADM quantities."""

    def __init__(self, IDtype: str, override_gauge_with_standard: bool = False) -> None:
        self.IDtype = IDtype

        self.gammaDD = ixp.zerorank2()
        self.KDD = ixp.zerorank2()
        self.alpha = cast(sp.Expr, sp.sympify(0))
        self.betaU = ixp.zerorank1()
        self.BU = ixp.zerorank1()

        self.r, self.th, self.ph = sp.symbols("r th ph", real=True)

        if IDtype == "UIUCBlackHole":
            ID_defines_gauge_quantities = False
            self.gammaDD, self.KDD = self.UIUCBlackHole()
        elif IDtype == "StaticTrumpet":
            ID_defines_gauge_quantities = True
            (
                self.gammaDD,
                self.KDD,
                self.alpha,
                self.betaU,
                self.BU,
            ) = self.StaticTrumpet()
        elif IDtype == "OffsetKerrSchild":
            ID_defines_gauge_quantities = True
            (
                self.gammaDD,
                self.KDD,
                self.alpha,
                self.betaU,
                self.BU,
            ) = self.OffsetKerrSchild()
        elif IDtype in ("TeukolskyWave", "TeukolskyWave21", "TeukolskyWave22"):
            # Teukolsky waves are defined in the alpha=1, beta^i=0 gauge; overriding it is an error.
            if override_gauge_with_standard:
                raise ValueError(f"{IDtype} requires alpha=1 and betaU=0")
            ID_defines_gauge_quantities = True
            mode = 20 if IDtype == "TeukolskyWave" else int(IDtype[-2:])
            (
                self.gammaDD,
                self.KDD,
                self.alpha,
                self.betaU,
                self.BU,
            ) = self.TeukolskyWave(mode=mode)
        else:
            raise ValueError(f"IDtype = {IDtype} is not supported.")

        if not ID_defines_gauge_quantities or override_gauge_with_standard:
            # Set shift betaU and time derivative of shift BU to zero -- the standard choice.
            self.betaU = ixp.zerorank1()
            self.BU = ixp.zerorank1()
            # Next compute alpha. Standard choice is alpha = 1/psi**2 (psi = BSSN conformal factor),
            #   where psi = exp(phi); chi = 1/psi**4; W = 1/psi**2
            adm2bssn = ADM_to_BSSN(
                self.gammaDD,
                self.KDD,
                self.betaU,
                self.BU,
                "Spherical",
                compute_cf_only=True,
            )
            EvolvedConformalFactor_cf = par.parval_from_str("EvolvedConformalFactor_cf")
            cf = adm2bssn.cf
            if EvolvedConformalFactor_cf == "phi":
                self.alpha = sp.exp(-2 * cf)
            elif EvolvedConformalFactor_cf == "chi":
                self.alpha = sp.sqrt(cf)
            elif EvolvedConformalFactor_cf == "W":
                self.alpha = cf
            else:
                raise ValueError(
                    f"Error EvolvedConformalFactor_cf type = {EvolvedConformalFactor_cf} unknown."
                )

    # fmt: off
    def UIUCBlackHole(self) -> Tuple[List[List[sp.Expr]], List[List[sp.Expr]]]:
        """
        Set ADM quantities for a spinning black hole in the UIUC initial data slicing, Liu, Etienne, & Shapiro (2009) https://arxiv.org/pdf/1001.4077.pdf.

        :return: A tuple containing the 3-metric gammaDD and the extrinsic curvature KDD.
        """
        M = par.register_CodeParameter("REAL", __name__, "M", 1.0, commondata=True)
        chi = par.register_CodeParameter(
            "REAL", __name__, "chi", 0.99, commondata=True, add_to_parfile=False
        )
        r = self.r
        th = self.th

        # Step 1: Set psi, the conformal factor:
        # Spin per unit mass
        a = M*chi

        # Defined under equation 1 in Liu, Etienne, & Shapiro (2009) https://arxiv.org/pdf/1001.4077.pdf
        # Boyer - Lindquist outer horizon
        rp = M + sp.sqrt(M**2 - a**2)
        # Boyer - Lindquist inner horizon
        rm = M - sp.sqrt(M**2 - a**2)

        # Boyer - Lindquist radius in terms of UIUC radius
        # Eq. 11
        # r_{BL} = r * ( 1 + r_+ / 4r )^2
        rBL = r*(1 + rp / (4*r))**2

        # Expressions found below Eq. 2
        # Sigma = r_{BL}^2 + a^2 cos^2 theta
        SIG = rBL**2 + a**2*sp.cos(th)**2

        # Delta = r_{BL}^2 - 2Mr_{BL} + a^2
        DEL = rBL**2 - 2*M*rBL + a**2

        # A = (r_{BL}^2 + a^2)^2 - Delta a^2 sin^2 theta
        AA = (rBL**2 + a**2)**2 - DEL*a**2*sp.sin(th)**2

        # *** The ADM 3-metric in spherical basis ***
        gammaDD = ixp.zerorank2()
        # Declare the nonzero components of the 3-metric
        # (Eq. 13 of Liu, Etienne, & Shapiro, https://arxiv.org/pdf/1001.4077.pdf):

        # ds^2 = Sigma (r + r_+/4)^2 / ( r^3 (r_{BL} - r_- ) * dr^2 +
        # Sigma d theta^2  +  (A sin^2 theta) / Sigma  *  d\phi^2

        gammaDD[0][0] = ((SIG*(r + rp/4)**2)/(r**3*(rBL - rm)))
        gammaDD[1][1] = SIG
        gammaDD[2][2] = AA/SIG*sp.sin(th)**2

        # *** The physical trace-free extrinsic curvature in spherical basis ***
        # Nonzero components of the extrinsic curvature K, given by
        # Eq. 14 of Liu, Etienne, & Shapiro, https://arxiv.org/pdf/1001.4077.pdf:
        KDD = ixp.zerorank2()


        # K_{r phi} = K_{phi r} = (Ma sin^2 theta) / (Sigma sqrt{A Sigma}) *
        #     [3r^4_{BL} + 2a^2 r^2_{BL} - a^4 - a^2 (r^2_{BL} - a^2) sin^2 theta] *
        #     (1 + r_+ / 4r) (1 / sqrt{r(r_{BL} - r_-)})

        KDD[0][2] = KDD[2][0] = (M*a*sp.sin(th)**2)/(SIG*sp.sqrt(AA*SIG))* \
                                (3*rBL**4 + 2*a**2*rBL**2 - a**4- a**2*(rBL**2 - a**2)*
                                 sp.sin(th)**2)*(1 + rp/(4*r))*1/sp.sqrt(r*(rBL - rm))

        # Components of the extrinsic curvature K, given by
        # Eq. 15 of Liu, Etienne, & Shapiro, https://arxiv.org/pdf/1001.4077.pdf:

        # K_{theta phi} = K_{phi theta} = -(2a^3 Mr_{BL} cos theta sin^3 theta) /
        #         (Sigma sqrt{A Sigma}) x (r - r_+ / 4) sqrt{(r_{BL} - r_-) / r }

        KDD[1][2] = KDD[2][1] = -((2*a**3*M*rBL*sp.cos(th)*sp.sin(th)**3)/
                                  (SIG*sp.sqrt(AA*SIG)))*(r - rp/4)*sp.sqrt((rBL - rm)/r)

        return gammaDD, KDD

    def StaticTrumpet(self) -> Tuple[List[List[sp.Expr]],List[List[sp.Expr]], sp.Expr, List[sp.Expr], List[sp.Expr]]:
        """
        Set ADM quantities for trumpet black hole initial data, Dennison and Baumgarte (2014) https://arxiv.org/pdf/1403.5484.pdf.

        :return: A tuple containing the 3-metric (gammaDD), the extrinsic curvature (KDD), the lapse function (alpha), the shift vector (betaU), and its time derivative (BU).
        """
        M = par.register_CodeParameter("REAL", __name__, "M", 1.0, commondata=True)
        r = self.r
        th = self.th

        # Step 1: Set psi, the StaticTrumpet conformal factor
        # Dennison and Baumgarte (2014) Eq. 13
        # https://arxiv.org/pdf/1403.5484.pdf

        # psi = sqrt{1 + M/r }
        psi0 = sp.sqrt(1 + M/r)

        # *** The physical spatial metric in spherical basis ***
        # Set the upper-triangle of the matrix...
        # Eq. 15
        # gamma_{ij} = psi^4 * eta_{ij}
        # eta_00 = 1, eta_11 = r^2, eta_22 = r^2 * sin^2 (theta)
        gammaDD = ixp.zerorank2()
        gammaDD[0][0] = psi0**4
        gammaDD[1][1] = psi0**4 * r**2
        gammaDD[2][2] = psi0**4 * r**2*sp.sin(th)**2
        # ... then apply symmetries to get the other components

        # *** The physical trace-free extrinsic curvature in spherical basis ***
        # Set the upper-triangle of the matrix...

        # Eq.19 and 20
        KDD = ixp.zerorank2()

        # K_{rr} = M / r^2
        KDD[0][0] = -M / r**2

        # K_{theta theta} = K_{phi phi} / sin^2 theta = M
        KDD[1][1] = M

        KDD[2][2] = M * sp.sin(th)**2
        # ... then apply symmetries to get the other components

        # Lapse function and shift vector
        # Eq. 15
        # alpha = r / (r+M)
        alpha = r / (r + M)

        betaU = ixp.zerorank1()
        # beta^r = Mr / (r + M)^2
        betaU[0] = M*r / (r + M)**2

        BU = ixp.zerorank1()

        return gammaDD, KDD, alpha, betaU, BU

    def OffsetKerrSchild(self) -> Tuple[List[List[sp.Expr]],List[List[sp.Expr]], sp.Expr, List[sp.Expr], List[sp.Expr]]:
        """
        Set ADM quantities for a spinning black hole in Kerr-Schild coordinates, with a radial offset r0 removing r<r0 region.
        See e.g., Etienne et al (2017) https://arxiv.org/pdf/1704.00599.pdf

        :return: A tuple containing the 3-metric (gammaDD), the extrinsic curvature (KDD), the lapse function (alpha), the shift vector (betaU), and its time derivative (BU).
        """
        M, a, r0 = par.register_CodeParameters("REAL", __name__, ["M", "a", "r0"], [1.0, 0.9, 1.0], commondata=True)
        r = self.r
        th = self.th

        # Step 1: Define rho^2, alpha, beta^(r_{KS}), beta^(theta), beta^(phi), gamma_{r_{KS}theta}, gamma_{theta\phi}

        # r_{KS} = r + r0
        rKS = r+r0

        # rho^2 = rKS^2 + a^2*cos^2(theta)
        rho2 = rKS*rKS + a*a*sp.cos(th)**2

        # alpha = 1/sqrt{1 + M*rKS/rho^2}
        alpha = 1/(sp.sqrt(1 + 2*M*rKS/rho2))

        # Initialize the shift vector, \beta^i, to zero.
        betaU = ixp.zerorank1()
        # beta^r = alpha^2*2Mr/rho^2
        betaU[0] = alpha*alpha*2*M*rKS/rho2

        # Time derivative of shift vector beta^i, B^i, is zero.
        BU = ixp.zerorank1()

        # Step 2: Define and construct nonzero components gamma_{r_{KS}r_{KS}}$, gamma_{r_{KS}phi},
        #         gamma_{thetatheta}, gamma_{phiphi}

        # Initialize \gamma_{ij} to zero.
        gammaDD = ixp.zerorank2()

        # gammaDD{rKS rKS} = 1 +2M*rKS/rho^2
        gammaDD[0][0] = 1 + 2*M*rKS/rho2

        # gammaDD{rKS phi} = -a*gammaDD{r r}*sin^2(theta)
        gammaDD[0][2] = gammaDD[2][0] = -a*gammaDD[0][0]*sp.sin(th)**2

        # gammaDD{theta theta} = rho^2
        gammaDD[1][1] = rho2

        # gammaDD{phi phi} = (rKS^2 + a^2 + 2Mr/rho^2*a^2*sin^2(theta))*sin^2(theta)
        gammaDD[2][2] = (rKS*rKS + a*a + 2*M*rKS*a*a*sp.sin(th)**2/rho2)*sp.sin(th)**2

        # Step 3: Define useful quantities A, B, C
        # A = (a^2*cos^2(2theta) + a^2 + 2r^2)
        A = a*a*sp.cos(2*th) + a*a + 2*rKS*rKS

        # B = A + 4M*rKS
        B = A + 4*M*rKS

        # D = \sqrt(2M*rKS/(a^2cos^2(theta) + rKS^2) + 1)
        D = sp.sqrt(2*M*rKS/(a*a*sp.cos(th)**2 + rKS*rKS) + 1)


        # Step 4: Define the extrinsic curvature in spherical polar coordinates

        # Establish the 3x3 zero-matrix
        KDD = ixp.zerorank2()

        # *** Fill in the nonzero components ***
        # *** This will create an upper-triangular matrix ***
        # K_{r r} = D(A+2Mr)/(A^2*B)[4M(a^2*cos(2theta) + a^2 - 2r^2)]
        KDD[0][0] = D*(A+2*M*rKS)/(A*A*B)*(4*M*(a*a*sp.cos(2*th)+a*a-2*rKS*rKS))

        # K_{r theta} = D/(AB)[8a^2*Mr*sin(theta)cos(theta)]
        KDD[0][1] = KDD[1][0] = D/(A*B)*(8*a*a*M*rKS*sp.sin(th)*sp.cos(th))

        # K_{r phi} = D/A^2[-2aMsin^2(theta)(a^2cos(2theta)+a^2-2r^2)]
        KDD[0][2] = KDD[2][0] =  D/(A*A)*(-2*a*M*sp.sin(th)**2*(a*a*sp.cos(2*th)+a*a-2*rKS*rKS))

        # K_{theta theta} = D/B[4Mr^2]
        KDD[1][1] = D/B*(4*M*rKS*rKS)

        # K_{theta phi} = D/(AB)*(-8*a^3*Mr*sin^3(theta)cos(theta))
        KDD[1][2] = KDD[2][1] = D/(A*B)*(-8*a**3*M*rKS*sp.sin(th)**3*sp.cos(th))

        # K_{phi phi} = D/(A^2*B)[2Mr*sin^2(theta)(a^4(M+3r)
        #   +4a^2r^2(2r-M)+4a^2r*cos(2theta)(a^2+r(M+2r))+8r^5)]
        KDD[2][2] = D/(A*A*B)*(2*M*rKS*sp.sin(th)**2*(a**4*(rKS-M)*sp.cos(4*th)
                                + a**4*(M+3*rKS)+4*a*a*rKS*rKS*(2*rKS-M)
                                + 4*a*a*rKS*sp.cos(2*th)*(a*a + rKS*(M + 2*rKS)) + 8*rKS**5))

        return gammaDD, KDD, alpha, betaU, BU

    # fmt: on

    def TeukolskyWave(
        self, branch: str = "regular", mode: int = 20
    ) -> Tuple[
        List[List[sp.Expr]], List[List[sp.Expr]], sp.Expr, List[sp.Expr], List[sp.Expr]
    ]:
        """
        Set ADM quantities for a linearized (l=2) Teukolsky wave in flat space, Teukolsky (1982) https://doi.org/10.1103/PhysRevD.26.745.
        Explicit formulas follow Rinne (2009) https://arxiv.org/pdf/0809.1761.pdf, with pulse F(q) = q exp(-lambda q^2).

        Conventions: gamma_{ij} = eta_{ij} + h_{ij}, where eta_{ij} = diag(1, r^2, r^2 sin^2 theta) is the flat
        spherical metric, with both indices down in (r, theta, phi) order. K_{ij} = -1/2 partial_t gamma_{ij} in the
        gauge alpha = 1, beta^i = 0, B^i = 0. Angular functions are Rinne's unnormalized real functions (cosine
        branch for m > 0). The outgoing wave depends on u = t - r, the ingoing wave on v = t + r.

        Note: near r = 0 the regular branch stays finite only through cancellation between large terms, so
        evaluating these expressions at very small r loses floating-point precision. A grid kernel that samples
        r -> 0 needs a series expansion there.

        :param branch: Wave branch: "outgoing", "ingoing", or "regular" (outgoing + ingoing; regular at the origin).
        :param mode: Teukolsky-wave mode: 20 for even-parity (2,0), 21 for odd-parity (2,1), or 22 for even-parity (2,2).
        :return: A tuple containing the 3-metric (gammaDD), the extrinsic curvature (KDD), the lapse function (alpha), the shift vector (betaU), and its time derivative (BU).
        :raises ValueError: If the requested mode or branch is unsupported.
        """
        if mode not in (20, 21, 22):
            raise ValueError("TeukolskyWave mode must be 20, 21 or 22")
        if branch not in ("outgoing", "ingoing", "regular"):
            raise ValueError("branch must be 'outgoing', 'ingoing', or 'regular'")

        Amp, lam, t0 = par.register_CodeParameters(
            "REAL", __name__, ["TW_Amp", "TW_lambda", "TW_t0"], [1.0e-6, 1.0, 0.0]
        )
        r = self.r
        th = self.th
        ph = self.ph

        # Step 1: Set the pulse derivatives Fn = d^n F/dq^n at the retarded time
        #         u = t - r (outgoing) and pn = Fn(v) at the advanced time v = t + r
        #         (ingoing). Since du/dt = dv/dt = 1, d^n/dt^n gives the n-th derivative.
        #         t stays symbolic until Step 4.
        t = sp.symbols("t", real=True)
        u = t - r
        v = t + r
        F = [sp.diff(u * sp.exp(-lam * u**2), t, n) for n in range(5)]
        p = [sp.diff(v * sp.exp(-lam * v**2), t, n) for n in range(5)]

        # Step 2: Set the radial coefficients: A, B, C (even parity; modes 20, 22) and
        #         K_rad, L_rad (odd parity; mode 21; not the extrinsic curvature K_{ij}).
        # Step 2.a: Outgoing wave.
        A_out = 24 * Amp * (F[2] / r**3 + 3 * F[1] / r**4 + 3 * F[0] / r**5)
        B_out = (
            -4
            * Amp
            * (F[3] / r**2 + 3 * F[2] / r**3 + 6 * F[1] / r**4 + 6 * F[0] / r**5)
        )
        C_out = (
            2
            * Amp
            * (
                F[4] / r
                + 2 * F[3] / r**2
                + 3 * F[2] / r**3
                + 3 * F[1] / r**4
                + 3 * F[0] / r**5
            )
        )
        K_rad_out = -4 * Amp * (F[2] / r**2 + 3 * F[1] / r**3 + 3 * F[0] / r**4)
        L_rad_out = (
            2 * Amp * (F[3] / r + 2 * F[2] / r**2 + 3 * F[1] / r**3 + 3 * F[0] / r**4)
        )
        # Step 2.b: Ingoing wave.
        A_in = 24 * Amp * (-p[2] / r**3 + 3 * p[1] / r**4 - 3 * p[0] / r**5)
        B_in = (
            4
            * Amp
            * (-p[3] / r**2 + 3 * p[2] / r**3 - 6 * p[1] / r**4 + 6 * p[0] / r**5)
        )
        C_in = (
            2
            * Amp
            * (
                -p[4] / r
                + 2 * p[3] / r**2
                - 3 * p[2] / r**3
                + 3 * p[1] / r**4
                - 3 * p[0] / r**5
            )
        )
        K_rad_in = 4 * Amp * (p[2] / r**2 - 3 * p[1] / r**3 + 3 * p[0] / r**4)
        L_rad_in = (
            2 * Amp * (p[3] / r - 2 * p[2] / r**2 + 3 * p[1] / r**3 - 3 * p[0] / r**4)
        )
        # Step 2.c: Select the branch; the regular wave is outgoing + ingoing.
        if branch == "outgoing":
            A, B, C, K_rad, L_rad = A_out, B_out, C_out, K_rad_out, L_rad_out
        elif branch == "ingoing":
            A, B, C, K_rad, L_rad = A_in, B_in, C_in, K_rad_in, L_rad_in
        else:
            A, B, C = A_out + A_in, B_out + B_in, C_out + C_in
            K_rad, L_rad = K_rad_out + K_rad_in, L_rad_out + L_rad_in

        # Step 3: Set the 3-metric gamma_{ij} = eta_{ij} + h_{ij},
        #         where eta_{ij} = diag(1, r^2, r^2 sin^2 theta).
        s = sp.sin(th)
        c = sp.cos(th)
        gammaDD = ixp.zerorank2()
        if mode == 20:
            # Step 3.a: Mode 20, even parity, axisymmetric: Y = 3 cos^2 theta - 1
            Y = 3 * c**2 - 1
            gammaDD[0][0] = 1 + A * Y
            gammaDD[0][1] = gammaDD[1][0] = -6 * r * B * s * c
            gammaDD[1][1] = r**2 * (1 - A * Y / 2 + 3 * C * s**2)
            gammaDD[2][2] = r**2 * s**2 * (1 - A * Y / 2 - 3 * C * s**2)
        elif mode == 22:
            # Step 3.b: Mode 22, even parity: P = cos(2 phi), Q = sin(2 phi)
            P = sp.cos(2 * ph)
            Q = sp.sin(2 * ph)
            gammaDD[0][0] = 1 + A * s**2 * P
            gammaDD[0][1] = gammaDD[1][0] = 2 * r * B * s * c * P
            gammaDD[0][2] = gammaDD[2][0] = -2 * r * B * s**2 * Q
            gammaDD[1][1] = r**2 * (1 - A * s**2 * P / 2 + C * (1 + c**2) * P)
            gammaDD[1][2] = gammaDD[2][1] = -2 * r**2 * C * s * c * Q
            gammaDD[2][2] = r**2 * s**2 * (1 - A * s**2 * P / 2 - C * (1 + c**2) * P)
        else:
            # Step 3.c: Mode 21, odd parity
            gammaDD[0][0] = sp.sympify(1)
            gammaDD[0][1] = gammaDD[1][0] = r * K_rad * c * sp.sin(ph)
            gammaDD[0][2] = gammaDD[2][0] = r * s * K_rad * (c**2 - s**2) * sp.cos(ph)
            gammaDD[1][1] = r**2 * (1 - L_rad * s * sp.sin(ph))
            gammaDD[1][2] = gammaDD[2][1] = -(r**2) * L_rad * s**2 * c * sp.cos(ph)
            gammaDD[2][2] = r**2 * s**2 * (1 + L_rad * s * sp.sin(ph))

        # Step 4: Set the extrinsic curvature K_{ij} = -1/2 partial_t gamma_{ij}
        #         (alpha = 1, beta^i = 0), then evaluate both on the slice t = t0.
        KDD = [[(-sp.diff(g, t) / 2).subs(t, t0) for g in row] for row in gammaDD]
        gammaDD = [[g.subs(t, t0) for g in row] for row in gammaDD]

        # Step 5: Set the gauge: unit lapse, zero shift, zero shift time derivative.
        alpha = sp.sympify(1)
        betaU = ixp.zerorank1()
        BU = ixp.zerorank1()

        return gammaDD, KDD, alpha, betaU, BU


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

    for ID_type in [
        "UIUCBlackHole",
        "StaticTrumpet",
        "OffsetKerrSchild",
        "TeukolskyWave",
        "TeukolskyWave21",
        "TeukolskyWave22",
    ]:
        ID = InitialData_Spherical(ID_type)
        results_dict = ve.process_dictionary_of_expressions(
            ID.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{ID_type}",
            results_dict,
        )

    # Step 1: Validate the Teukolsky-wave data against the companion equation packet
    #         (packet sections in parentheses).
    TW = InitialData_Spherical("TeukolskyWave")
    xx = [TW.r, TW.th, TW.ph]
    TW_Amp, TW_lambda, TW_t0 = sp.symbols("TW_Amp TW_lambda TW_t0", real=True)
    etaDD = ixp.zerorank2()
    etaDD[0][0] = sp.sympify(1)
    etaDD[1][1] = xx[0] ** 2
    etaDD[2][2] = xx[0] ** 2 * sp.sin(xx[1]) ** 2

    # Step 1.a: External benchmark (Sec. 13): outgoing mode 20 at t0 = r = 1,
    #           theta = pi/4, phi = pi/6, lambda = 1, Amp = 10^-6; compared exactly.
    point = {
        xx[0]: 1,
        xx[1]: sp.pi / 4,
        xx[2]: sp.pi / 6,
        TW_lambda: 1,
        TW_t0: 1,
        TW_Amp: sp.Rational(1, 10**6),
    }
    gammaDD_bench = ixp.zerorank2()
    gammaDD_bench[0][0] = sp.Rational(250009, 250000)
    gammaDD_bench[1][1] = sp.Rational(199991, 200000)
    gammaDD_bench[2][2] = sp.Rational(1000009, 2000000)
    KDD_bench = ixp.zerorank2()
    KDD_bench[0][0] = sp.Rational(18, 10**6)
    KDD_bench[0][1] = KDD_bench[1][0] = sp.Rational(72, 10**6)
    KDD_bench[1][1] = sp.Rational(-765, 10**7)
    KDD_bench[2][2] = sp.Rational(2925, 10**8)
    gammaDD_out, KDD_out, _, _, _ = TW.TeukolskyWave(branch="outgoing", mode=20)
    for i in range(3):
        for j in range(3):
            assert sp.simplify(gammaDD_out[i][j].subs(point) - gammaDD_bench[i][j]) == 0
            assert sp.simplify(KDD_out[i][j].subs(point) - KDD_bench[i][j]) == 0

    # Step 1.b: Structural checks for each mode, regular branch, evaluated with the same fixed
    #           numerical inputs as the trusted results (randomly sampled inputs can lose
    #           too much precision to these expressions' 1/r^5 terms). With
    #           h_{ij} = gamma_{ij} - eta_{ij} and the parity matrix P = diag(1, -1, 1):
    #           gauge (Sec. 4), symmetry, zero amplitude and amplitude linearity (Sec. 15),
    #           regular = outgoing + ingoing and K_{ij}(t0 = 0) = 0 (Sec. 8),
    #           trace-free eta^{ij} h_{ij} = 0 (Sec. 9), and the mode's angular symmetry:
    #           mode 20 independent of phi (Sec. 10); mode 22 unchanged under
    #           phi -> phi + pi and odd under phi -> phi + pi/2 (Sec. 11); tensor parity
    #           P h(r, pi - theta, phi + pi) P = +h (modes 20, 22) or -h (mode 21) (Sec. 12).
    # Christoffel symbols of the diagonal eta_{ij}, used in Step 1.c:
    # Gamma^i_{jk} = 1/2 eta^{ii} (eta_{ij,k} + eta_{ik,j} - eta_{jk,i}).
    GammaUDD = ixp.zerorank3()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                GammaUDD[i][j][k] = (
                    sp.diff(etaDD[i][j], xx[k])
                    + sp.diff(etaDD[i][k], xx[j])
                    - sp.diff(etaDD[j][k], xx[i])
                ) / (2 * etaDD[i][i])
    div_point = {
        xx[0]: sp.Rational(13, 10),
        xx[1]: sp.Rational(7, 10),
        xx[2]: sp.Rational(11, 10),
        TW_lambda: sp.Rational(9, 10),
        TW_t0: sp.Rational(3, 10),
        TW_Amp: 1,
    }
    parity = [1, -1, 1]
    flip = {xx[1]: sp.pi - xx[1], xx[2]: xx[2] + sp.pi}
    for mode_TW in (20, 21, 22):
        gammaDD_reg, KDD_reg, alpha_TW, betaU_TW, BU_TW = TW.TeukolskyWave(mode=mode_TW)
        assert alpha_TW == 1 and betaU_TW == [0, 0, 0] and BU_TW == [0, 0, 0]
        gammaDD_out, KDD_out, _, _, _ = TW.TeukolskyWave(
            branch="outgoing", mode=mode_TW
        )
        gammaDD_in, KDD_in, _, _, _ = TW.TeukolskyWave(branch="ingoing", mode=mode_TW)
        parity_sign = -1 if mode_TW == 21 else 1
        trace_h = sp.sympify(0)
        for i in range(3):
            trace_h += (gammaDD_reg[i][i] - etaDD[i][i]) / etaDD[i][i]
            for j in range(3):
                hDD_ij = gammaDD_reg[i][j] - etaDD[i][j]
                KDD_ij = KDD_reg[i][j]
                assert ve.check_zero(
                    hDD_ij - (gammaDD_reg[j][i] - etaDD[j][i]),
                    fixed_mpfs_for_free_symbols=True,
                )
                assert ve.check_zero(
                    KDD_ij - KDD_reg[j][i], fixed_mpfs_for_free_symbols=True
                )
                for expr in (hDD_ij, KDD_ij):
                    assert ve.check_zero(
                        expr.subs(TW_Amp, 0), fixed_mpfs_for_free_symbols=True
                    )
                    assert ve.check_zero(
                        expr.subs(TW_Amp, 2 * TW_Amp) - 2 * expr,
                        fixed_mpfs_for_free_symbols=True,
                    )
                    assert ve.check_zero(
                        parity[i] * parity[j] * expr.subs(flip, simultaneous=True)
                        - parity_sign * expr,
                        fixed_mpfs_for_free_symbols=True,
                    )
                    if mode_TW == 20:
                        assert ve.check_zero(
                            sp.diff(expr, xx[2]), fixed_mpfs_for_free_symbols=True
                        )
                    if mode_TW == 22:
                        assert ve.check_zero(
                            expr.subs(xx[2], xx[2] + sp.pi) - expr,
                            fixed_mpfs_for_free_symbols=True,
                        )
                        assert ve.check_zero(
                            expr.subs(xx[2], xx[2] + sp.pi / 2) + expr,
                            fixed_mpfs_for_free_symbols=True,
                        )
                assert ve.check_zero(
                    gammaDD_reg[i][j]
                    - gammaDD_out[i][j]
                    - gammaDD_in[i][j]
                    + etaDD[i][j],
                    fixed_mpfs_for_free_symbols=True,
                )
                assert ve.check_zero(
                    KDD_ij - KDD_out[i][j] - KDD_in[i][j],
                    fixed_mpfs_for_free_symbols=True,
                )
                assert ve.check_zero(
                    KDD_ij.subs(TW_t0, 0), fixed_mpfs_for_free_symbols=True
                )
        assert ve.check_zero(trace_h, fixed_mpfs_for_free_symbols=True)

        # Step 1.c: Check that h_{ij} is divergence-free, eta^{jj} D_j h_{ij} = 0 (Sec. 20).
        #           Check the outgoing and ingoing branches separately (their sum is the
        #           regular branch), at an exact rational point with 50 significant digits;
        #           ve.check_zero() cannot resolve these derivative-heavy expressions to zero.
        for gammaDD_branch in (gammaDD_out, gammaDD_in):
            for i in range(3):
                div_hD = sp.sympify(0)
                for j in range(3):
                    # D_j h_{ij} = h_{ij,j} - Gamma^l_{ji} h_{lj} - Gamma^l_{jj} h_{il}
                    hDD_covdD = sp.diff(gammaDD_branch[i][j] - etaDD[i][j], xx[j])
                    for l in range(3):
                        hDD_covdD -= GammaUDD[l][j][i] * (
                            gammaDD_branch[l][j] - etaDD[l][j]
                        ) + GammaUDD[l][j][j] * (gammaDD_branch[i][l] - etaDD[i][l])
                    div_hD += hDD_covdD / etaDD[j][j]
                assert abs(sp.N(div_hD.subs(div_point), 50)) < 1e-40

    # Step 2: Check that the standard-gauge override is rejected for every mode.
    for ID_type in ["TeukolskyWave", "TeukolskyWave21", "TeukolskyWave22"]:
        try:
            InitialData_Spherical(ID_type, override_gauge_with_standard=True)
        except ValueError as exc:
            assert f"{ID_type} requires alpha=1 and betaU=0" in str(exc)
        else:
            raise AssertionError(f"{ID_type} must reject the standard-gauge override.")

    print(
        "*** PASS *** TeukolskyWave packet benchmark and analytic identities validated."
    )
