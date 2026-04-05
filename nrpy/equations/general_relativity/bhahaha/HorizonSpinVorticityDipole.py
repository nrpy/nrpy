# equations/general_relativity/bhahaha/HorizonSpinVorticityDipole.py
"""
NRPy module for Owen et al's vorticity-dipole-based black hole spin diagnostics.

Reference: Owen et al, arXiv:1708.07325v2

Goals:
  * Provide symbolic, single-pass integrands on the apparent horizon for:
      - Vorticity-dipole spin-axis vector:
            N^i = ∮ Omega * x^i dA
      - Quasilocal angular-momentum vector:
            J^i = (1/(8*pi)) ∮ omega_i * phi^i[X^a] dA
    where:
      - Omega is the curvature (2D curl) of the normal-bundle connection,
      - omega_i is constructed from K_ij and s^i following Owen et al.,
      - phi^i[X^a] are rotation generators built from background Cartesian
        scalars X^a, tangent to the horizon.

  * Enforce:
      - No SymPy .subs() or simplify() calls.
      - All indexed quantities via nrpy.indexedexp.
      - All FD-like derivatives named using "*_dD".
      - Use reference_metric.reference_metric.
      - Public interface minimal and focused on C code generation.
      - All per-point integrands depend only on:
            { h, h_dD, h_dDD, hDD, hDD_dD,
              W, W_dD, aDD, trK }
        plus reference-metric data and standard constants.
      - Spin axis integrals computable in a single pass.
      - Spin magnitude integrals computable in a single pass (may be separate).
      - No requirement for additional horizon-only gridfunctions computed
        in separate preprocessing passes; any derivatives beyond the given
        list appear only as symbolic "var_dD" placeholders that the
        C layer may fill directly from volume data as needed.

Key choices (physics and implementation):

  * Level set:
        F(r,theta,phi) = r - h(theta,phi) = 0
    in Spherical reference-metric coordinates xx[0]=r, xx[1]=theta, xx[2]=phi.

  * Physical 3-metric:
        gammabar_ij = hDD_ij * ReDD_ij + ghat_ij
        gamma_ij    = gammabar_ij / W^2       (W = e^{-2phi})

  * Horizon unit normal:
        s_i ∝ ∂_i F,
        s^i = gamma^{ij} ∂_j F / lambda,
        lambda^2 = gamma^{ij} ∂_i F ∂_j F.
    Orientation "outward" vs "inward" is configurable; the final
    omega_i * phi^i and Omega * x^i are invariant under s^i -> -s^i.

  * Extrinsic curvature (BSSN input):
        Abar_ij = aDD_ij * ReDD_ij
        K_ij    = (1/W^2) Abar_ij + (1/3) gamma_ij trK

  * Horizon vorticity one-form:
        omega_i = q_i^p K_pq s^q
    with q_i^j = delta_i^j - s_i s^j the tangential projector.
    The trace term drops out after projection; we encode an explicit
    K_sign to reconcile global sign conventions of K_ij.

  * Rotation generators (Owen et al.):
    For each Cartesian scalar X^a in {x,y,z}, define:
        phi^i[X^a] = epsilon^{ijk} s_j ∂_k X^a
    using the 3D Levi-Civita tensor built from det(gamma), ensuring
    tangentiality and covariance. The quasilocal spin component is:
        J_a = (1/(8*pi)) ∮ omega_i phi^i[X^a] dA.

  * Normal-bundle curvature scalar:
        Omega = epsilon^{AB} D_A omega_B
    is best computed from 2D derivatives on the surface. To keep this
    module single-pass compatible and independent of a specific FD
    stencil, we:
      - Define Omega through a 3D-projected covariant curl expression,
        expecting inputs nabla_p omega_q as "omegaD_covD_pq" if users
        wish to build it explicitly; and
      - Provide a direct Omega-x^i integrand helper that can be used once
        Omega is available as a gridfunction built from volume data.

    This keeps the interface consistent with the allowed primitive
    variables and avoids forcing extra horizon-surface passes here.

Public interface (minimal):

  * JCart_densU[3]:
        Per-point densities:
            JCart_densU[a] = (1/(8*pi)) * omega_i * phi^i[X^a]
        for a = 0:x,1:y,2:z. These are to be integrated with the
        existing area element and quadrature weights in BHaHAHA.

  * Omega_xCart_densU(Omega_sym)[3]:
        Method returning:
            (Omega_sym * X^a)
        per axis, to be multiplied by dA externally.
        Omega_sym is a SymPy symbol or expression representing Omega
        in terms of primitives. This keeps Omega construction flexible
        while supporting single-pass dipole integrals.

  * J_about_axis(nU)[0]:
        Given a constant axis n^a, returns the scalar density
            (1/(8*pi)) * omega_i * phi^i[zeta]
        with zeta = n·X, normalized n.

  * Omega_from_covDomega(covD_omegaD):
        Helper constructing scalar Omega from supplied covariant
        derivative nabla_p omega_q via 3D-to-2D projection:
            Omega = epsilon^{ijk} s_i q_j^p q_k^q (nabla_p omega_q).

Usage notes:

  * All expressions are built in Spherical reference-metric coordinates.
  * This module does not perform any FD differentiation itself; symbols
    like "hDD_dD" or "omegaD_covD" are declared and may be populated
    by generated C code using volume-level data in a single pass.

Authors: Ralston Graves
         Zachariah B. Etienne
         zachetie **at** gmail **dot* com
"""

from typing import Dict, List, Optional, Tuple

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.reference_metric as refmetric


class HorizonSpinVorticityDipoleClass:
    """
    Symbolic construction of vorticity-dipole and quasilocal-spin integrands.

    Assumptions: a horizon defined by F(r,theta,phi) = r - h(theta,phi) = 0
    in a Spherical reference-metric chart.

    Configuration parameters:
      CoordSystem        : must be "Spherical".
      enable_rfm_precompute : use precompute variant of reference metric.
      K_input            : "BSSN" or "external".
                            - "BSSN": K_ij built from aDD, trK, W.
                            - "external": K_ij provided directly as KDD.
      normal_orientation : "outward" or "inward" (flip s^i).
      K_sign             : +1 or -1; reconciles K_ij sign conventions.
    """

    def __init__(
        self,
        CoordSystem: str = "Spherical",
        enable_rfm_precompute: bool = False,
        K_input: str = "BSSN",
        normal_orientation: str = "outward",
        K_sign: int = +1,
    ):
        # Basic validation
        if CoordSystem != "Spherical":
            raise ValueError(
                f"Unsupported CoordSystem '{CoordSystem}'. Only 'Spherical' is supported."
            )
        if K_input not in ("BSSN", "external"):
            raise ValueError("K_input must be 'BSSN' or 'external'.")
        if normal_orientation not in ("outward", "inward"):
            raise ValueError("normal_orientation must be 'outward' or 'inward'.")
        if K_sign not in (+1, -1):
            raise ValueError("K_sign must be +1 or -1.")

        self._CoordSystem = CoordSystem
        self._K_input = K_input
        self._normal_orientation = normal_orientation
        self._K_sign = sp.Integer(1) if K_sign == +1 else sp.Integer(-1)

        # Reference metric object
        self._rfm = refmetric.reference_metric[
            (
                self._CoordSystem + "_rfm_precompute"
                if enable_rfm_precompute
                else self._CoordSystem
            )
        ]

        # Core build (no .subs(), no simplify())
        self._build_all()

    # -------------------------------------------------------------------------
    # Internal: construct all needed quantities
    # -------------------------------------------------------------------------

    def _build_all(self) -> None:
        """
        Build gamma_ij, s^i, omega_i, rotation generators, and J densities.

        Only minimal, necessary attributes are made public at the end.
        """
        # Short-hands
        # r = self._rfm.xx[0] # Unused; h serves as surface
        th = self._rfm.xx[1]
        ph = self._rfm.xx[2]

        # ------------------------------------------------------------------
        # 0. Horizon shape and level-set gradient
        # ------------------------------------------------------------------
        hh = sp.Symbol("hh", real=True)

        # h_dD: partial derivatives of h with respect to (r,theta,phi) coords.
        # h depends only on angles, but we do not enforce that algebraically.
        hh_dD = ixp.declarerank1("hh_dD")

        # F = r - h(theta,phi); F_{,i}
        F_dD = ixp.zerorank1()
        F_dD[0] = sp.sympify(1)
        F_dD[1] = -hh_dD[1]
        F_dD[2] = -hh_dD[2]

        # ------------------------------------------------------------------
        # 1. Physical 3-metric gamma_ij from conformal + reference metric
        # ------------------------------------------------------------------
        hDD = ixp.declarerank2("hDD", symmetry="sym01")
        WW = sp.Symbol("WW", real=True)

        gammabarDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                gammabarDD[i][j] = (
                    hDD[i][j] * self._rfm.ReDD[i][j] + self._rfm.ghatDD[i][j]
                )

        self._gammaDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self._gammaDD[i][j] = gammabarDD[i][j] / (WW * WW)

        # Inverse metric and determinant
        self._gammaUU, self._detgamma = ixp.symm_matrix_inverter3x3(self._gammaDD)
        self._sqrt_detgamma = sp.sqrt(self._detgamma)

        # ------------------------------------------------------------------
        # 2. Horizon unit normal s^i
        # ------------------------------------------------------------------
        # lambda^2 = gamma^{ij} F_i F_j
        lamb2 = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                lamb2 += self._gammaUU[i][j] * F_dD[i] * F_dD[j]
        lamb = sp.sqrt(lamb2)

        # s_i = F_i / lambda
        self._sD = ixp.zerorank1()
        for i in range(3):
            self._sD[i] = F_dD[i] / lamb

        if self._normal_orientation == "inward":
            for i in range(3):
                self._sD[i] = -self._sD[i]

        # s^i = gamma^{ij} s_j
        self._sU = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                self._sU[i] += self._gammaUU[i][j] * self._sD[j]

        # ------------------------------------------------------------------
        # 3. Extrinsic curvature K_ij
        # ------------------------------------------------------------------
        if self._K_input == "external":
            KDD = ixp.declarerank2("KDD", symmetry="sym01")
        else:
            aDD = ixp.declarerank2("aDD", symmetry="sym01")
            trK = sp.Symbol("trK", real=True)

            AbarDD = ixp.zerorank2()
            for i in range(3):
                for j in range(3):
                    AbarDD[i][j] = aDD[i][j] * self._rfm.ReDD[i][j]

            exp4phi = 1 / (WW * WW)
            KDD = ixp.zerorank2()
            for i in range(3):
                for j in range(3):
                    KDD[i][j] = (
                        exp4phi * AbarDD[i][j]
                        + sp.Rational(1, 3) * self._gammaDD[i][j] * trK
                    )

        # ------------------------------------------------------------------
        # 4. Tangential projector q_i^j
        # ------------------------------------------------------------------
        deltaDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                deltaDD[i][j] = sp.sympify(1) if i == j else sp.sympify(0)

        # q_i^j = delta_i^j - s_i s^j
        self._qDU = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self._qDU[i][j] = deltaDD[i][j] - self._sD[i] * self._sU[j]

        # ------------------------------------------------------------------
        # 5. Horizon vorticity one-form omega_i = K_sign * q_i^p K_pq s^q
        # ------------------------------------------------------------------
        K_dot_sU_D = ixp.zerorank1()
        for p in range(3):
            for q in range(3):
                K_dot_sU_D[p] += KDD[p][q] * self._sU[q]

        self._omegaD = ixp.zerorank1()
        for i in range(3):
            for p in range(3):
                self._omegaD[i] += self._K_sign * self._qDU[i][p] * K_dot_sU_D[p]

        # ------------------------------------------------------------------
        # 6. Levi-Civita tensor epsilon^{ijk} and rotation generators phi^i[X^a]
        # ------------------------------------------------------------------
        # epsilon^{ijk} with upper indices from sqrt(detgamma)
        self._epsUUU = ixp.LeviCivitaTensorUUU_dim3_rank3(self._sqrt_detgamma)

        # On-surface Cartesian scalars X^a = (x,y,z):
        # Use standard spherical relations with r->hh (no .subs; define directly).
        s_th = sp.sin(th)
        c_th = sp.cos(th)
        s_ph = sp.sin(ph)
        c_ph = sp.cos(ph)

        xCart = hh * s_th * c_ph
        yCart = hh * s_th * s_ph
        zCart = hh * c_th

        self._xCartU = [xCart, yCart, zCart]

        # dX^a/dx^k in (r,theta,phi), evaluated with r=hh:
        dX_dxx = [[sp.sympify(0) for _ in range(3)] for __ in range(3)]

        # a=0 -> x = r sin th cos ph -> with r=hh:
        dX_dxx[0][0] = s_th * c_ph  # dx/dr
        dX_dxx[0][1] = hh * c_th * c_ph  # dx/dtheta
        dX_dxx[0][2] = -hh * s_th * s_ph  # dx/dphi

        # a=1 -> y
        dX_dxx[1][0] = s_th * s_ph
        dX_dxx[1][1] = hh * c_th * s_ph
        dX_dxx[1][2] = hh * s_th * c_ph

        # a=2 -> z
        dX_dxx[2][0] = c_th
        dX_dxx[2][1] = -hh * s_th
        dX_dxx[2][2] = sp.sympify(0)

        # phi^i[X^a] = epsilon^{ijk} s_j ∂_k X^a
        self._phiU_by_cart = [[sp.sympify(0) for _ in range(3)] for __ in range(3)]
        for a in range(3):
            for i in range(3):
                acc = sp.sympify(0)
                for j in range(3):
                    for k in range(3):
                        acc += self._epsUUU[i][j][k] * self._sD[j] * dX_dxx[a][k]
                self._phiU_by_cart[a][i] = acc

        # ------------------------------------------------------------------
        # 7. Quasilocal J^i densities: JCart_densU[a] = (1/(8*pi)) omega_i phi^i[X^a]
        # ------------------------------------------------------------------
        one_over_8pi = sp.Integer(1) / (8 * sp.pi)

        omega_dot_phi = ixp.zerorank1()
        for a in range(3):
            for i in range(3):
                omega_dot_phi[a] += self._omegaD[i] * self._phiU_by_cart[a][i]

        self.JCart_densU = ixp.zerorank1()
        for a in range(3):
            self.JCart_densU[a] = one_over_8pi * omega_dot_phi[a]

    # -------------------------------------------------------------------------
    # Public helpers
    # -------------------------------------------------------------------------

    def J_about_axis(
        self, nU: Optional[List[sp.Expr]] = None
    ) -> Tuple[sp.Expr, List[sp.Expr]]:
        """
        Scalar spin density about a constant Cartesian axis n^a.

        Implements:
          J_density[n] = (1/(8*pi)) * omega_i * phi^i[zeta],
        where zeta = n·X, with X^a = (x,y,z) and phi^i[zeta] defined as:
          phi^i[zeta] = epsilon^{ijk} s_j ∂_k (n·X).

        :param nU: optional [n^x, n^y, n^z]; if None, symbolic components used.
        :return: (J_density_expr, nU_unit)
        """
        if nU is None:
            nU = [
                sp.Symbol("nU0", real=True),
                sp.Symbol("nU1", real=True),
                sp.Symbol("nU2", real=True),
            ]

        # Normalize axis
        n2 = sp.sympify(0)
        for a in range(3):
            n2 += nU[a] * nU[a]
        nmag = sp.sqrt(n2)
        nU_unit = [nU[a] / nmag for a in range(3)]

        # d/dx^k of zeta = n·X:
        # ∂_k zeta = sum_a n^a ∂_k X^a
        Dzeta_dD = ixp.zerorank1()
        # Rebuild local dX_dxx as above (no state mutation)
        th = self._rfm.xx[1]
        ph = self._rfm.xx[2]
        hh = sp.Symbol("hh", real=True)
        s_th = sp.sin(th)
        c_th = sp.cos(th)
        s_ph = sp.sin(ph)
        c_ph = sp.cos(ph)
        dX_dxx = [[sp.sympify(0) for _ in range(3)] for __ in range(3)]
        dX_dxx[0][0] = s_th * c_ph
        dX_dxx[0][1] = hh * c_th * c_ph
        dX_dxx[0][2] = -hh * s_th * s_ph
        dX_dxx[1][0] = s_th * s_ph
        dX_dxx[1][1] = hh * c_th * s_ph
        dX_dxx[1][2] = hh * s_th * c_ph
        dX_dxx[2][0] = c_th
        dX_dxx[2][1] = -hh * s_th
        dX_dxx[2][2] = sp.sympify(0)

        for k in range(3):
            for a in range(3):
                Dzeta_dD[k] += nU_unit[a] * dX_dxx[a][k]

        # phi^i[zeta] = eps^{ijk} s_j ∂_k zeta
        phiU = ixp.zerorank1()
        for i in range(3):
            acc = sp.sympify(0)
            for j in range(3):
                for k in range(3):
                    acc += self._epsUUU[i][j][k] * self._sD[j] * Dzeta_dD[k]
            phiU[i] = acc

        # J_density = (1/(8*pi)) * omega_i phi^i[zeta]
        one_over_8pi = sp.Integer(1) / (8 * sp.pi)
        Jdens = sp.sympify(0)
        for i in range(3):
            Jdens += self._omegaD[i] * phiU[i]

        return one_over_8pi * Jdens, nU_unit

    def Omega_from_covDomega(self, covD_omegaD: List[List[sp.Expr]]) -> sp.Expr:
        """
        Build boost-gauge-invariant scalar Omega = epsilon^{AB} D_A omega_B.
        This uses 3D projection:
          Omega = epsilon^{ijk} s_i q_j^p q_k^q (nabla_p omega_q),

        where covD_omegaD[p][q] = nabla_p omega_q is supplied
        (e.g., built from volume data and Christoffels in C).

        :param covD_omegaD: 3x3 nested list for nabla_p omega_q.
        :return: SymPy expression for Omega.
        """
        Omega = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for p in range(3):
                        for q in range(3):
                            Omega += (
                                self._epsUUU[i][j][k]
                                * self._sD[i]
                                * self._qDU[j][p]
                                * self._qDU[k][q]
                                * covD_omegaD[p][q]
                            )
        return Omega

    def Omega_xCart_densU(self, Omega_expr: sp.Expr) -> List[sp.Expr]:
        """
        Given Omega ito primitive variables, return per-axis integrands Omega * X^a.

        These are intended to be multiplied by the existing area element
        and quadrature weights in a single pass:

          N^a ≈ sum_grid Omega_xCart_densU[a] * dA * weights.

        :param Omega_expr: SymPy expression for Omega.
        :return: List [Omega * x, Omega * y, Omega * z].
        """
        densU = ixp.zerorank1()
        for a in range(3):
            densU[a] = Omega_expr * self._xCartU[a]
        return densU


class HorizonSpinVorticityDipoleClass_dict(Dict[str, HorizonSpinVorticityDipoleClass]):
    """
    Dictionary-style accessor for HorizonSpinVorticityDipoleClass.

      * "Spherical"
      * "Spherical_rfm_precompute"
    """

    def __getitem__(self, key: str) -> HorizonSpinVorticityDipoleClass:
        if key not in self:
            if key == "Spherical":
                obj = HorizonSpinVorticityDipoleClass(
                    CoordSystem="Spherical",
                    enable_rfm_precompute=False,
                )
            elif key == "Spherical_rfm_precompute":
                obj = HorizonSpinVorticityDipoleClass(
                    CoordSystem="Spherical",
                    enable_rfm_precompute=True,
                )
            else:
                raise KeyError(
                    "Supported keys: 'Spherical', 'Spherical_rfm_precompute'."
                )
            dict.__setitem__(self, key, obj)
        return dict.__getitem__(self, key)

    def __setitem__(self, key: str, value: HorizonSpinVorticityDipoleClass) -> None:
        if key not in ["Spherical", "Spherical_rfm_precompute"]:
            raise KeyError("Supported keys: 'Spherical', 'Spherical_rfm_precompute'.")
        dict.__setitem__(self, key, value)

    def __delitem__(self, key: str) -> None:
        dict.__delitem__(self, key)


# Public accessor
HorizonSpinVorticityDipole = HorizonSpinVorticityDipoleClass_dict()


if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.validate_expressions.validate_expressions as ve

    # Run doctests (none currently defined)
    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Validation: export only the quasilocal JCart_densU components
    for validation_key in ["Spherical", "Spherical_rfm_precompute"]:
        modobj = HorizonSpinVorticityDipole[validation_key]

        export_only = {}
        for axis_idx in range(3):
            export_only[f"JCart_densU_{axis_idx}"] = modobj.JCart_densU[axis_idx]
            # print(modobj.JCart_densU[axis_idx].free_symbols)

        results_dict = ve.process_dictionary_of_expressions(
            export_only, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{validation_key}",
            results_dict,
        )
