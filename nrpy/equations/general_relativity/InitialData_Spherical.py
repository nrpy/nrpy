"""
Set up initial data for solving Einstein's equations of general relativity, for data most naturally specified in Spherical coordinates.

Outputs ADM quantities.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from typing import Tuple, List
import sympy as sp
import nrpy.indexedexp as ixp
import nrpy.params as par

from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN

# NRPy+: This module depends on the parameter EvolvedConformalFactor_cf,
#        which is defined in BSSN.BSSN_quantities
import nrpy.equations.general_relativity.BSSN_quantities  # pylint: disable=unused-import


class InitialData_Spherical:
    """Construct and store Spherical initial data for Einstein's equations of general relativity, as ADM quantities."""

    def __init__(self, IDtype: str, override_gauge_with_standard: bool = False) -> None:
        self.IDtype = IDtype

        self.gammaDD = ixp.zerorank2()
        self.KDD = ixp.zerorank2()
        self.alpha = sp.sympify(0)
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
        """Set ADM quantities for a spinning black hole in the UIUC initial data slicing, Liu, Etienne, & Shapiro (2009) https://arxiv.org/pdf/1001.4077.pdf."""
        M, chi = par.register_CodeParameters("REAL", __name__, ["M", "chi"], [1.0, 0.99], commondata=True)
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
        """Set ADM quantities for trumpet black hole initial data, Dennison and Baumgarte (2014) https://arxiv.org/pdf/1403.5484.pdf."""
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
