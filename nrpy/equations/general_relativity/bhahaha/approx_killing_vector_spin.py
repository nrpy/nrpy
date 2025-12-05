# equations/general_relativity/bhahaha/approx_killing_vector_spin.py
<<<<<<< HEAD

=======
>>>>>>> origin/main
"""
Single-pass integrands for AKV-based spin diagnostics on apparent horizons.

 (1) Approximate Killing Vector (AKV) spin axis on an AH
 (2) Spin magnitude components J_m once the axis is chosen

Public API:
 - sqrtq : scalar surface element sqrt(det q_AB)
 - Hmn_integrand[3][3] : AKV quadratic-form integrands (symmetric)
 - Nmn_integrand[3][3] : AKV normalization-form integrands (symmetric)
 - Jm_integrand[3]     : spin integrands per rigid-rotation basis direction

Note: The global 1/(8*pi) prefactor for J is applied by the caller during surface integration.

Sources:
 - BHaHAHA (https://arxiv.org/abs/2505.15912)
 - QuasiLocalMeasures docs (https://einsteintoolkit.org/thornguide/EinsteinAnalysis/QuasiLocalMeasures/documentation.html)
 - QuasiLocalMeasures source code (https://bitbucket.org/einsteintoolkit/einsteinanalysis)

 (Foundations of AKV on S²; H_{mn}, N_{mn}, and shear minimization)
 - Cook & Whiting (2007), *Phys. Rev. D 76, 041501*, arXiv:0706.0199
 (Eigenvalue formulation of approximate Killing fields)
 - Beetle (2008), *Phys. Rev. D 78, 084043*, arXiv:0808.1745
 (Practical AKV spin on apparent horizons; axis determination and normalization)
 - Owen (2009), *Phys. Rev. D 80, 084012*, arXiv:0907.0280
 (Hands-on recipe used widely in NR implementations of AKV spin)
 - Lovelace, Owen, Pfeiffer & Chu (2008), *Phys. Rev. D 78, 084017*, arXiv:0805.4192
 (Quasilocal angular momentum integrand used here: (K_ij − K γ_ij) φ^i s^j)
 - Brown & York (1993), *Phys. Rev. D 47, 1407*, arXiv:gr-qc/9209012
 (Reference-metric / pole-regular formalism used in NRPy+ code generation)
 - Ruchlin, Etienne & Baumgarte (2018), *Phys. Rev. D 97, 064036*, arXiv:1712.07658

Author: Wesley Inselman
"""

from typing import Dict, List, cast

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.reference_metric as refmetric


class ApproxKillingSpinClass:
    """
    Build symbolic expressions for AKV spin and axis on an apparent horizon.

    Coordinates: spherical (r,theta,phi), AH given by F=r-h(theta,phi)=0.
    Inputs follow NRPy 'var_dD' naming; no .subs or simplify used.
    The caller applies the overall 1/(8*pi) prefactor when integrating J.
    """

    def __init__(self, CoordSystem: str = "Spherical"):
        # Accept either "Spherical" or "Spherical_rfm_precompute" and parse locally.
        if CoordSystem not in ["Spherical", "Spherical_rfm_precompute"]:
            raise ValueError("Use 'Spherical' or 'Spherical_rfm_precompute'.")
        self.CoordSystem = CoordSystem
        self._rfm = refmetric.reference_metric[self.CoordSystem]

        # Public outputs
        self.sqrtq: sp.Expr
        self.Hmn_integrand: List[List[sp.Expr]]
        self.Nmn_integrand: List[List[sp.Expr]]
        self.Jm_integrand: List[sp.Expr]

        # Build pipeline
        self._build_3metric_and_derivs_()
        self._build_surface_geometry_()
        self._build_l1_basis_and_forms_()
        self._build_spin_integrands_()

    # ------------------------
    # 3-metric & K_ij in spherical-basis
    # ------------------------
    def _build_3metric_and_derivs_(self) -> None:
        # Conformal factor W and derivatives
        W = sp.Symbol("WW", real=True)
        WdD = ixp.declarerank1("partial_D_WW")

        # h_ij (barred-metric deviation) and derivatives
        hDD = ixp.declarerank2("hDD", symmetry="sym01")
        partial_D_hDD = cast(
            List[List[List[sp.Expr]]],
            ixp.declarerank3("partial_D_hDD", symmetry="sym12"),
        )

        # gammabar_ij
        self._gammabarDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self._gammabarDD[i][j] = (
                    hDD[i][j] * self._rfm.ReDD[i][j] + self._rfm.ghatDD[i][j]
                )

        # Physical 3-metric gamma_ij = gammabar_ij / W^2
        self._gammaDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self._gammaDD[i][j] = self._gammabarDD[i][j] / (W * W)

        # gamma_ij,k = -2/W^3 * gammabar_ij W_{,k} + gammabar_ij,k / W^2
        gammabarDDdD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    gammabarDDdD[i][j][k] = (
                        self._rfm.ghatDDdD[i][j][k]
                        + partial_D_hDD[k][i][j] * self._rfm.ReDD[i][j]
                        + hDD[i][j] * self._rfm.ReDDdD[i][j][k]
                    )
        self._gammaDDdD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self._gammaDDdD[i][j][k] = (
                        -2 * self._gammabarDD[i][j] * WdD[k]
                    ) / (W * W * W) + gammabarDDdD[i][j][k] / (W * W)

        # gamma^ij and det(gamma)
        self._gammaUU, _detgamma_unused = ixp.symm_matrix_inverter3x3(self._gammaDD)

        # Extrinsic curvature from tracefree piece aDD and trace trK
        aDD = ixp.declarerank2("aDD", symmetry="sym01")
        self._trK = sp.Symbol("trK", real=True)
        self._AbarDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self._AbarDD[i][j] = aDD[i][j] * self._rfm.ReDD[i][j]
        self._KDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self._KDD[i][j] = (
                    self._AbarDD[i][j] / (W * W)
                    + sp.Rational(1, 3) * self._gammaDD[i][j] * self._trK
                )

    # ------------------------
    # 2-surface geometry on r=h(theta,phi)
    # ------------------------
    def _build_surface_geometry_(self) -> None:
        # Horizon shape and derivatives
        self._h = sp.Symbol("hh", real=True)
        self._h_dD = ixp.declarerank1(
            "hh_dD"
        )  # [h_r, h_theta, h_phi]; angular components nonzero in practice
        self._h_dDD = ixp.declarerank2("hh_dDD", symmetry="sym01")

        # Level-set gradient F_{,i} in spherical basis: [1, -h_,theta, -h_,phi]
        F_dD = [sp.sympify(1), -self._h_dD[1], -self._h_dD[2]]

        # Pullback Jacobian p^i{}_A with y^A=(theta,phi) and x^i=(r,theta,phi) on r=h(theta,phi)
        self._pU_i_A = ixp.zerorank2()  # pU[i][A]
        self._pU_i_A[0][1] = self._h_dD[1]  # dr/dtheta
        self._pU_i_A[0][2] = self._h_dD[2]  # dr/dphi
        self._pU_i_A[1][1] = sp.sympify(1)  # dtheta/dtheta
        self._pU_i_A[2][2] = sp.sympify(1)  # dphi/dphi

        # Induced 2-metric q_AB = p^i_A p^j_B gamma_ij
        self._qDD = ixp.zerorank2()
        for A in range(1, 3):
            for B in range(1, 3):
                s = sp.sympify(0)
                for i in range(3):
                    for j in range(3):
                        s += (
                            self._pU_i_A[i][A]
                            * self._pU_i_A[j][B]
                            * self._gammaDD[i][j]
                        )
                self._qDD[A][B] = s

        # det(q_AB) and q^{AB}
        q2x2 = [[self._qDD[1][1], self._qDD[1][2]], [self._qDD[2][1], self._qDD[2][2]]]
        self._qUU2, self._detq2 = ixp.symm_matrix_inverter2x2(q2x2)
        self._qUU = ixp.zerorank2()
        self._qUU[1][1] = self._qUU2[0][0]
        self._qUU[1][2] = self._qUU2[0][1]
        self._qUU[2][1] = self._qUU2[1][0]
        self._qUU[2][2] = self._qUU2[1][1]
        self.sqrtq = sp.sqrt(self._detq2)

        # ∂_C q_AB via chain rule
        self._qDDdD = ixp.zerorank3()  # q_AB,C
        dxk_dyC = ixp.zerorank2()
        dxk_dyC[0][1] = self._h_dD[1]
        dxk_dyC[0][2] = self._h_dD[2]
        dxk_dyC[1][1] = sp.sympify(1)
        dxk_dyC[2][2] = sp.sympify(1)

        dp_iA_dC = ixp.zerorank3()
        dp_iA_dC[0][1][1] = self._h_dDD[1][1]
        dp_iA_dC[0][1][2] = self._h_dDD[1][2]
        dp_iA_dC[0][2][1] = self._h_dDD[2][1]
        dp_iA_dC[0][2][2] = self._h_dDD[2][2]

        for A in range(1, 3):
            for B in range(1, 3):
                for C in range(1, 3):
                    s = sp.sympify(0)
                    for i in range(3):
                        for j in range(3):
                            for k in range(3):
                                s += (
                                    self._pU_i_A[i][A]
                                    * self._pU_i_A[j][B]
                                    * self._gammaDDdD[i][j][k]
                                    * dxk_dyC[k][C]
                                )
                    for i in range(3):
                        for j in range(3):
                            s += (
                                dp_iA_dC[i][A][C]
                                * self._pU_i_A[j][B]
                                * self._gammaDD[i][j]
                            )
                            s += (
                                self._pU_i_A[i][A]
                                * dp_iA_dC[j][B][C]
                                * self._gammaDD[i][j]
                            )
                    self._qDDdD[A][B][C] = s

        # 2D Christoffels Γ^C_{AB}
        self._Gamma2D = ixp.zerorank3()
        for A in range(1, 3):
            for B in range(1, 3):
                for C in range(1, 3):
                    val = sp.sympify(0)
                    for D in range(1, 3):
                        val += (
                            sp.Rational(1, 2)
                            * self._qUU[C][D]
                            * (
                                self._qDDdD[A][D][B]
                                + self._qDDdD[B][D][A]
                                - self._qDDdD[A][B][D]
                            )
                        )
                    self._Gamma2D[A][B][C] = val

        # Unit normal s^i = gamma^{ij} F_{,j} / sqrt(gamma^{mn} F_{,m} F_{,n})
        unnorm_sU = ixp.zerorank1()
        for i in range(3):
            tmp = sp.sympify(0)
            for j in range(3):
                tmp += self._gammaUU[i][j] * F_dD[j]
            unnorm_sU[i] = tmp
        lam2 = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                lam2 += self._gammaUU[i][j] * F_dD[i] * F_dD[j]
        lam = sp.sqrt(lam2)
        self._sU = ixp.zerorank1()
        for i in range(3):
            self._sU[i] = unnorm_sU[i] / lam

    # ------------------------
    # l=1 rigid-rotation basis on (theta,phi) and AKV forms
    # ------------------------
    def _build_l1_basis_and_forms_(self) -> None:
        # Standard SO(3) generators in spherical angles:
        #   xi_(z) = ∂_phi
        #   xi_(x) =  sin(phi) ∂_theta + cot(theta) cos(phi) ∂_phi
        #   xi_(y) = -cos(phi) ∂_theta + cot(theta) sin(phi) ∂_phi
        theta = self._rfm.xx[1]
        phi = self._rfm.xx[2]
        sin = sp.sin
        cos = sp.cos
        # IMPORTANT (pole regularization):
        # cot(theta) = cos(theta)/sin(theta). In generated code, the rfm approach
        # used by BHaHAHA keeps regular parts finite at the poles; this symbolic
        # representation is compatible with that regularization.
        cot = cos(theta) / sin(theta)

        self._xiU = [
            ixp.zerorank1() for _ in range(3)
        ]  # contravariant components over A∈{1,2}
        # z
        self._xiU[2][1] = sp.sympify(0)
        self._xiU[2][2] = sp.sympify(1)
        # x
        self._xiU[0][1] = sin(phi)
        self._xiU[0][2] = cot * cos(phi)
        # y
        self._xiU[1][1] = -cos(phi)
        self._xiU[1][2] = cot * sin(phi)

        # Covariant xi_A = q_AB xi^B
        self._xiD = [ixp.zerorank1() for _ in range(3)]
        for m in range(3):
            for A in range(1, 3):
                tmp = sp.sympify(0)
                for B in range(1, 3):
                    tmp += self._qDD[A][B] * self._xiU[m][B]
                self._xiD[m][A] = tmp

        # Angular derivatives ∂_C xi^A (analytic)
        d_xiU = [ixp.zerorank2() for _ in range(3)]  # [m][A][C]
        d_cot_dtheta = -(
            sp.sympify(1) + cot * cot
        )  # d/dtheta(cot)= -csc^2 = -(1+cot^2)
        # m=2 (z): zeros already
        # m=0 (x)
        d_xiU[0][1][1] = sp.sympify(0)
        d_xiU[0][1][2] = cos(phi)
        d_xiU[0][2][1] = d_cot_dtheta * cos(phi)
        d_xiU[0][2][2] = -cot * sin(phi)
        # m=1 (y)
        d_xiU[1][1][1] = sp.sympify(0)
        d_xiU[1][1][2] = sin(phi)
        d_xiU[1][2][1] = d_cot_dtheta * sin(phi)
        d_xiU[1][2][2] = cot * cos(phi)

        # Divergence D_A xi^A
        self._div_xi = [sp.sympify(0) for _ in range(3)]
        for m in range(3):
            divm = sp.sympify(0)
            for A in range(1, 3):
                divm += d_xiU[m][A][A]
                for C in range(1, 3):
                    divm += self._Gamma2D[A][C][A] * self._xiU[m][C]
            self._div_xi[m] = divm

        # 2D covariant derivative of xi_B: D_A xi_B = ∂_A xi_B - Γ^C_{AB} xi_C
        # with ∂_A xi_B = (∂_A q_BD) xi^D + q_BD (∂_A xi^D)
        self._D_xiD = [ixp.zerorank2() for _ in range(3)]  # [m][A][B]
        for m in range(3):
            for A in range(1, 3):
                for B in range(1, 3):
                    term = sp.sympify(0)
                    for Dd in range(1, 3):
                        term += self._qDDdD[B][Dd][A] * self._xiU[m][Dd]
                    for Dd in range(1, 3):
                        term += self._qDD[B][Dd] * d_xiU[m][Dd][A]
                    for Cc in range(1, 3):
                        term += -self._Gamma2D[A][B][Cc] * self._xiD[m][Cc]
                    self._D_xiD[m][A][B] = term

        # Conformal Killing (shear) operator: S_AB[m] = D_A xi_B + D_B xi_A - (D_C xi^C) q_AB
        self._SDD = [ixp.zerorank2() for _ in range(3)]
        for m in range(3):
            for A in range(1, 3):
                for B in range(1, 3):
                    self._SDD[m][A][B] = (
                        self._D_xiD[m][A][B]
                        + self._D_xiD[m][B][A]
                        - self._div_xi[m] * self._qDD[A][B]
                    )

        # Quadratic-form integrand H_mn = S_AB[m] S_CD[n] q^{AC} q^{BD} * sqrt(q)
        H = ixp.zerorank2()
        for m in range(3):
            for n in range(3):
                s = sp.sympify(0)
                for A in range(1, 3):
                    for B in range(1, 3):
                        for C in range(1, 3):
                            for D in range(1, 3):
                                s += (
                                    self._SDD[m][A][B]
                                    * self._SDD[n][C][D]
                                    * self._qUU[A][C]
                                    * self._qUU[B][D]
                                )
                H[m][n] = s * self.sqrtq

        # Normalization-form integrand N_mn = q_AB xi^A_(m) xi^B_(n) * sqrt(q)
        N = ixp.zerorank2()
        for m in range(3):
            for n in range(3):
                t = sp.sympify(0)
                for A in range(1, 3):
                    for B in range(1, 3):
                        t += self._qDD[A][B] * self._xiU[m][A] * self._xiU[n][B]
                N[m][n] = t * self.sqrtq

        # Optional symmetrization for finite-precision robustness
        self.Hmn_integrand = ixp.zerorank2()
        self.Nmn_integrand = ixp.zerorank2()
        for m in range(3):
            for n in range(3):
                self.Hmn_integrand[m][n] = sp.Rational(1, 2) * (H[m][n] + H[n][m])
                self.Nmn_integrand[m][n] = sp.Rational(1, 2) * (N[m][n] + N[n][m])

    # ------------------------
    # Spin magnitude integrands per basis direction
    # ------------------------
    def _build_spin_integrands_(self) -> None:
        # J[phi] = ∮ (K_ij - K gamma_ij) phi^i s^j * sqrt(q) dtheta dphi
        # The global factor 1/(8*pi) is applied by the caller.
        KminusKg = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                KminusKg[i][j] = self._KDD[i][j] - self._trK * self._gammaDD[i][j]

        # Basis phi^i_(m) = p^i_A xi^A_(m)
        phiU_m = [ixp.zerorank1() for _ in range(3)]
        for m in range(3):
            for i in range(3):
                tmp = sp.sympify(0)
                for A in range(1, 3):
                    tmp += self._pU_i_A[i][A] * self._xiU[m][A]
                phiU_m[m][i] = tmp

        # Jm_integrand[m] = (K_ij - K gamma_ij) phi^i s^j * sqrt(q)
        self.Jm_integrand = ixp.zerorank1()
        for m in range(3):
            val = sp.sympify(0)
            for i in range(3):
                for j in range(3):
                    val += KminusKg[i][j] * phiU_m[m][i] * self._sU[j]
            self.Jm_integrand[m] = val * self.sqrtq


class ApproxKillingSpinClass_dict(Dict[str, ApproxKillingSpinClass]):
    """Lazy mapping from coordinate-system name to ApproxKillingSpinClass instances."""

    def __getitem__(self, key: str) -> ApproxKillingSpinClass:
        if key not in self:
            if key not in ["Spherical", "Spherical_rfm_precompute"]:
                raise KeyError("Use 'Spherical' or 'Spherical_rfm_precompute'.")
            print(f"Setting up ApproxKillingSpin[{key}]...")
            self.__setitem__(key, ApproxKillingSpinClass(CoordSystem=key))
        return dict.__getitem__(self, key)

    def __setitem__(self, key: str, value: ApproxKillingSpinClass) -> None:
        if key not in ["Spherical", "Spherical_rfm_precompute"]:
            raise KeyError("Use 'Spherical' or 'Spherical_rfm_precompute'.")
        dict.__setitem__(self, key, value)


ApproxKillingSpin = ApproxKillingSpinClass_dict()


if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.validate_expressions.validate_expressions as ve

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted}")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    for validation_key in ["Spherical", "Spherical_rfm_precompute"]:
        akv = ApproxKillingSpin[validation_key]

        # Export only the final, public AKV quantities needed by the implementation.
        export_only = {}

        # Scalar surface element
        export_only["sqrtq"] = akv.sqrtq

        # Quadratic-form integrands H_mn
        for m_idx in range(3):
            for n_idx in range(3):
                export_only[f"Hmn_integrand_{m_idx}{n_idx}"] = akv.Hmn_integrand[m_idx][
                    n_idx
                ]

        # Normalization-form integrands N_mn
        for m_idx in range(3):
            for n_idx in range(3):
                export_only[f"Nmn_integrand_{m_idx}{n_idx}"] = akv.Nmn_integrand[m_idx][
                    n_idx
                ]

        # Spin integrands J_m
        for m_idx in range(3):
            export_only[f"Jm_integrand_{m_idx}"] = akv.Jm_integrand[m_idx]

        results_dict = ve.process_dictionary_of_expressions(
            export_only, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{validation_key}",
            results_dict,
        )
