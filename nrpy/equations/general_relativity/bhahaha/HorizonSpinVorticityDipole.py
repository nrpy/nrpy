# equations/general_relativity/bhahaha/HorizonSpinVorticityDipole.py
"""
HorizonSpinVorticityDipole.py.
---------------------------------------

Symbolic NRPy construction of single-pass integrands on a closed 2-surface S
for the Cartesian spin components J^a using the Owen identity and
reference-metric infrastructure. Everything is built on the surface F=0 from
the outset (no post-hoc substitutions).

This version keeps runtime behavior identical but improves and clarifies
documentation, with explicit links to the notation used in Owen et al.,
"Black Hole Spin Axis in Numerical Relativity" (2017): arXiv:1708.07325v2

Mathematical summary
---------------------------------
Objects live on a spatial slice with physical 3-metric gamma_ij and
extrinsic curvature K_ij (ADM/BSSN). Let S be the closed 2-surface defined
as the level set F=0 with unit normal s^i (spatial). Define the tangential
projector:
    q_i^j = delta_i^j - s_i s^j .

Owen et al. (and Brown-York) introduce a tangential "momentum density"
one-form omega_A on S. Restricted to spatial indices and projected tangent
to S, the trace term drops out and we obtain (see their Sec. 2):
    omega_i = q_i^p K_pq s^q .

Because different codes adopt different sign conventions for K_ij, we make
that choice explicit and configurable:
    omega_i = (K_sign) * q_i^p K_pq s^q ,
where K_sign = +1 matches the Owen/Brown-York convention if your K_ij
matches theirs; choose K_sign = -1 to compensate for the opposite global
sign in your stored K_ij.

Let epsilon^{ijk} be the 3D Levi-Civita tensor density with upper indices
(constructed here from sqrt(det gamma)). Given a Cartesian scalar X^a in
{ x, y, z }, the associated coordinate-rotation generator tangent to S is
computed without Christoffels as
    phi^i[X^a] = epsilon^{ijk} s_j d_k X^a ,
where d_k X^a are partial derivatives of the Cartesian scalars with respect
to the reference-metric coordinates, all evaluated directly on S.

The quasilocal angular momentum densities for the three Cartesian axes are
then
    J_a_density = (1/(8*pi)) * omega_i * phi^i[X^a] ,  a in {x,y,z} .

Optional Omega curvature of the normal bundle
---------------------------------------------
Owen et al. also use
    Omega := epsilon^{AB} D_A omega_B ,
which equals (in 3D index notation with projections)
    Omega = epsilon^{ijk} s_i q_j^p q_k^q (nabla_p omega_q) .
On a closed S and for any scalar zeta with phi_A[zeta] = epsilon_A^B D_B zeta,
the following identity holds (integration by parts, no boundary term):
    integral_S (omega_A phi^A[zeta]) dA = integral_S (zeta * Omega) dA .
This module provides a helper to assemble Omega from a supplied covariant
derivative tensor nabla_p omega_q, so you can compare (1/(8*pi)) * omega.phi
densities to (1/(8*pi)) * zeta * Omega without constructing Christoffels here.

Boost gauge and orientation notes
---------------------------------
Under a boost of the null normals adapted to S, omega_A -> omega_A - D_A a.
The scalar Omega is boost invariant. Flipping the surface normal s^i -> -s^i
flips both omega and phi, leaving omega.phi unchanged, so J_a_density is
independent of the choice "outward" vs "inward" up to the overall sign choice
you control via normal_orientation.

Construction strategy
---------------------
This module does not build expressions in the 3D bulk and then substitute
r -> h(theta,phi). Instead, background quantities with r-dependence (flat
spherical reference metric ghat_ij, "rescaling" Re_ij, and d(Cart^a)/d(rfm^k))
are constructed directly on the surface r == h(theta,phi). Thus every
contraction is intrinsically on S and we never call .subs.

Public members
--------------
* JCart_densU[a] : (1/(8*pi)) * (omega_i * phi^i[X^a]) for a = 0:x,1:y,2:z.
* J_about_axis(nU): scalar spin density about a constant Cartesian axis n^a.
* Omega_scalar_from_covDomega(covD_omegaD): builds Omega from nabla_p omega_q.
* omega_moment_density_from_zeta(zeta, covD_omegaD): returns (1/(8*pi)) * zeta * Omega.

Configuration parameters
------------------------
* CoordSystem: reference-metric chart; must be "Spherical".
* enable_rfm_precompute: reuse the rfm precompute variant (coords only).
* K_input: "BSSN" (rebuild K_ij from aDD, trK, W) or "external" (supply K_ij).
* normal_orientation: "outward" or "inward".
* K_sign: +1 or -1 to reconcile extrinsic-curvature sign conventions.

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
    Build symbolic single-pass integrands for J on S: {F=0}.

    Parameters
    ----------
    CoordSystem : str
        Reference-metric chart; must be "Spherical".
    enable_rfm_precompute : bool
        If True, use the precompute variant of the reference metric (for coords only).
    K_input : {"BSSN","external"}
        If "BSSN", reconstruct K_ij from aDD, trK, and W;
        if "external", accept physical K_ij directly via KDD.
    normal_orientation : {"outward","inward"}
        Orientation choice for s^i; inward flips signs of s^i and s_i.
    K_sign : {+1,-1}
        Global sign to reconcile K_ij convention. We implement
            omega_i = (K_sign) * q_i^p K_pq s^q .
        Set K_sign=+1 to match Owen/Brown-York if your K_ij matches their sign.
        Set K_sign=-1 if your stored K_ij has the opposite sign.

    Public interface
    ----------------
    * JCart_densU[a] : (1/(8*pi)) * (omega_i phi^i[X^a]) for a in {0,1,2}.
    * J_about_axis(nU) -> (J_density_expr, nU_unit)
      Scalar spin density for a constant axis n^a (normalized).
    * Omega_scalar_from_covDomega(covD_omegaD) -> Omega
      Given covariant derivatives nabla_p omega_q, form
          Omega = epsilon^{ijk} s_i q_j^p q_k^q (nabla_p omega_q) .
      No Christoffels are built internally.
    * omega_moment_density_from_zeta(zeta, covD_omegaD) -> (1/(8*pi)) * zeta * Omega
    """

    # -------------------- ctor --------------------

    def __init__(
        self,
        CoordSystem: str = "Spherical",
        enable_rfm_precompute: bool = False,
        K_input: str = "BSSN",
        normal_orientation: str = "outward",
        K_sign: int = +1,
    ):
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
        # Store the chosen K-sign convention as a SymPy integer for exactness.
        self._K_sign = sp.Integer(1) if K_sign == +1 else sp.Integer(-1)

        # Coordinate interface (theta,phi symbols etc.). We do not use ghat/Re/Jac arrays from here.
        self._rfm = refmetric.reference_metric[
            (
                (self._CoordSystem + "_rfm_precompute")
                if enable_rfm_precompute
                else self._CoordSystem
            )
        ]

        # Core build: everything is assembled intrinsically on the surface (no .subs ever).
        self._build_symbolic_expressions()

    # -------------------- internal helpers (private) --------------------

    def _spherical_surface_background(
        self,
    ) -> Tuple[
        List[List[sp.Expr]], List[List[sp.Expr]], List[List[sp.Expr]], sp.Symbol
    ]:
        """
        Construct the on-surface flat-spherical background and Cartesian Jacobians.

        This function assumes ``r == h(theta, phi)`` from the outset and builds the
        on-surface reference metric, rescaling tensor, Cartesian Jacobians, and the
        horizon-shape symbol ``hh``.

        :return: A 4-tuple ``(ghatDD_surf, ReDD_surf, dUCart_surf, hh)`` where
                 each entry is evaluated on the surface ``r=h(theta,phi)``.

        Doctests:
        TBD
        """
        th = self._rfm.xx[1]  # theta
        ph = self._rfm.xx[2]  # phi
        hh = sp.Symbol("hh", real=True)  # scalar h(theta,phi) on the surface

        # Flat spherical metric on the surface: diag(1, r^2, r^2 sin^2 theta), with r->hh.
        ghatDD_surf = [[sp.sympify(0) for _ in range(3)] for __ in range(3)]
        ghatDD_surf[0][0] = sp.Integer(1)
        ghatDD_surf[1][1] = hh**2
        ghatDD_surf[2][2] = (hh**2) * sp.sin(th) ** 2

        # In NRPy conformal machinery, ReDD plays the role of reference scaling factors.
        # Using ReDD = ghatDD on-surface ensures consistency with the standard split.
        ReDD_surf = [[ghatDD_surf[i][j] for j in range(3)] for i in range(3)]

        # Cartesian scalars: x = r sin theta cos phi, y = r sin theta sin phi, z = r cos theta, with r->hh.
        s_th, c_th = sp.sin(th), sp.cos(th)
        s_ph, c_ph = sp.sin(ph), sp.cos(ph)

        # Build d(Cart^a)/d(rfm^k) evaluated at r=hh.
        dUCart_surf = [[sp.sympify(0) for _ in range(3)] for __ in range(3)]

        # a = 0 -> x
        dUCart_surf[0][0] = s_th * c_ph  # dx/dr
        dUCart_surf[0][1] = hh * c_th * c_ph  # dx/dtheta
        dUCart_surf[0][2] = -hh * s_th * s_ph  # dx/dphi

        # a = 1 -> y
        dUCart_surf[1][0] = s_th * s_ph
        dUCart_surf[1][1] = hh * c_th * s_ph
        dUCart_surf[1][2] = hh * s_th * c_ph

        # a = 2 -> z
        dUCart_surf[2][0] = c_th
        dUCart_surf[2][1] = -hh * s_th
        dUCart_surf[2][2] = sp.Integer(0)

        return ghatDD_surf, ReDD_surf, dUCart_surf, hh

    def _build_symbolic_expressions(self) -> None:
        """
        Assemble required tensors using NRPy idioms without constructing Christoffels.

        This function constructs all background objects (flat metric, rescalings,
        Cartesian Jacobians) directly on the surface ``r=h(theta,phi)`` and then
        builds the projected quantities needed for the spin densities.

        Doctests:
        TBD
        """
        # ---- On-surface background: flat spherical metric, rescalings, and dCart/drfm ----
        ghatDD_surf, ReDD_surf, dUCart_surf, _ = self._spherical_surface_background()

        # ---- Level set F = r - h(theta,phi) and its first derivatives (in r,theta,phi basis) ----
        # h depends only on (theta,phi); we keep a symbolic hh_dD for clarity.
        hh_dD = ixp.declarerank1("hh_dD")
        F_dD = ixp.zerorank1()
        F_dD[0] = sp.sympify(1)
        F_dD[1] = -hh_dD[1]
        F_dD[2] = -hh_dD[2]

        # ---- Physical 3-metric gamma_ij from conformal split, evaluated on the surface ----
        #   gammabar_ij = ghat_ij + h_ij * Re_ij   (Re_ij chosen as ghat_ij on-surface)
        #   gamma_ij    = gammabar_ij / W^2        (with W = e^{-2 phi})
        hDD = ixp.declarerank2("hDD", symmetry="sym01")
        WW = sp.Symbol("WW", real=True)

        gammabarDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                gammabarDD[i][j] = ghatDD_surf[i][j] + hDD[i][j] * ReDD_surf[i][j]

        self._gammaDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self._gammaDD[i][j] = gammabarDD[i][j] / (WW * WW)

        # Invert gamma_ij -> gamma^ij and compute det gamma, sqrt(det gamma)
        self._gammaUU, self._detgamma = ixp.symm_matrix_inverter3x3(self._gammaDD)
        self._sqrt_detgamma = sp.sqrt(self._detgamma)

        # ---- Unit normal s_i and s^i to the surface F=0 (within the spatial slice) ----
        # s_i proportional to d_i F, normalized by lambda = sqrt(gamma^{ij} d_iF d_jF).
        lamb2 = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                lamb2 += self._gammaUU[i][j] * F_dD[i] * F_dD[j]
        lamb = sp.sqrt(lamb2)

        self._sD = ixp.zerorank1()
        for i in range(3):
            self._sD[i] = F_dD[i] / lamb

        if self._normal_orientation == "inward":
            for i in range(3):
                self._sD[i] = -self._sD[i]

        self._sU = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                self._sU[i] += self._gammaUU[i][j] * self._sD[j]

        # ---- Extrinsic curvature K_ij (either external or built from BSSN), on the surface ----
        if self._K_input == "external":
            KDD = ixp.declarerank2("KDD", symmetry="sym01")
            trK = sp.Symbol(
                "trK", real=True
            )  # not needed in omega_i after tangential projection
        else:
            aDD = ixp.declarerank2("aDD", symmetry="sym01")
            trK = sp.Symbol("trK", real=True)

            AbarDD = ixp.zerorank2()
            for i in range(3):
                for j in range(3):
                    AbarDD[i][j] = aDD[i][j] * ReDD_surf[i][j]

            exp4phi = 1 / (WW * WW)  # since W = e^{-2 phi} implies e^{4 phi} = 1/W^2
            KDD = ixp.zerorank2()
            for i in range(3):
                for j in range(3):
                    KDD[i][j] = (
                        exp4phi * AbarDD[i][j]
                        + sp.Rational(1, 3) * self._gammaDD[i][j] * trK
                    )

        # ---- Tangential projector q_i^j = delta_i^j - s_i s^j ----
        deltaDU = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                deltaDU[i][j] = sp.sympify(1) if i == j else sp.sympify(0)

        qDU = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                qDU[i][j] = deltaDU[i][j] - self._sD[i] * self._sU[j]

        # ---- Horizon vorticity one-form omega_i (stored with a lower index) ----
        # Derivation sketch:
        # Starting from omega_i = h_i^p (K_pq - K gamma_pq) s^q with h_i^p = delta_i^p - s_i s^p,
        # the trace term cancels after tangential projection:
        #   h_i^p gamma_pq s^q = (delta_i^p - s_i s^p) s_p = s_i - s_i (s^p s_p) = 0 .
        # Hence omega_i = h_i^p K_pq s^q = q_i^p K_pq s^q, up to the global K_sign.
        K_dot_sU_D = ixp.zerorank1()
        for p in range(3):
            for q in range(3):
                K_dot_sU_D[p] += KDD[p][q] * self._sU[q]

        self._omegaD = ixp.zerorank1()
        for i in range(3):
            for p in range(3):
                # omega_i = (K_sign) * q_i^p K_pq s^q
                self._omegaD[i] += self._K_sign * qDU[i][p] * K_dot_sU_D[p]

        # ---- Curved-space Levi-Civita with upper indices and phi^i[X^a] generators ----
        # epsilon^{ijk} = hat_epsilon^{ijk}/sqrt(det gamma); phi^i[X^a] = epsilon^{ijk} s_j d_k X^a .
        epsUUU = ixp.LeviCivitaTensorUUU_dim3_rank3(self._sqrt_detgamma)

        self._phiU_by_cart = [[sp.sympify(0) for _ in range(3)] for __ in range(3)]
        for a in range(3):
            for i in range(3):
                acc = sp.sympify(0)
                for j in range(3):
                    for k in range(3):
                        acc += epsUUU[i][j][k] * self._sD[j] * dUCart_surf[a][k]
                self._phiU_by_cart[a][i] = acc

        # ---- Public J densities (scaled by 1/(8*pi)); all are on-surface by construction ----
        one_over_8pi = sp.Integer(1) / (8 * sp.pi)

        omega_dot_phi_no_pref_U = ixp.zerorank1()
        for a in range(3):
            for i in range(3):
                omega_dot_phi_no_pref_U[a] += self._omegaD[i] * self._phiU_by_cart[a][i]

        self.JCart_densU = ixp.zerorank1()
        for a in range(3):
            self.JCart_densU[a] = one_over_8pi * omega_dot_phi_no_pref_U[a]

        # ---- Local omega.phi helper (unscaled) ----
        # Useful for comparing against (1/(8*pi)) * zeta * Omega identity.
        self._omega_dot_phiU = ixp.zerorank1()
        for a in range(3):
            self._omega_dot_phiU[a] = omega_dot_phi_no_pref_U[a]

        # Cache objects needed for optional Omega computations
        self._epsUUU = epsUUU
        self._qDU = qDU

    # -------------------- public spin-direction helper --------------------

    def J_about_axis(
        self, nU: Optional[List[sp.Expr]] = None
    ) -> Tuple[sp.Expr, List[sp.Expr]]:
        """
        Compute the scalar spin density about a constant Cartesian axis.

        This function evaluates
        ``J[n]_dens = (1/(8*pi)) * omega_i * epsilon^{ijk} s_j D_k (n·X)``
        using on-surface Jacobians and the precomputed normal and vorticity.

        :param nU: Optional list ``[n^x, n^y, n^z]`` of axis components. If ``None``,
                   symbolic components are used and normalized internally.
        :return: A tuple ``(J_density_expr, nU_unit)`` where ``nU_unit`` is the normalized axis.

        Doctests:
        TBD
        """
        if nU is None:
            nU = [
                sp.Symbol("nU0", real=True),
                sp.Symbol("nU1", real=True),
                sp.Symbol("nU2", real=True),
            ]

        # Normalize the input axis to make the returned unit direction explicit.
        n2 = sp.sympify(0)
        for cart_axis in range(3):
            n2 += nU[cart_axis] * nU[cart_axis]
        nmag = sp.sqrt(n2 + sp.sympify(0))
        nU_unit = [nU[cart_axis] / nmag for cart_axis in range(3)]

        # Recreate on-surface Jacobians for a self-contained method (no shared state mutation).
        _, _, dUCart_surf, _ = self._spherical_surface_background()

        Dzeta_dD = ixp.zerorank1()
        for k in range(3):
            for cart_axis in range(3):
                Dzeta_dD[k] += nU_unit[cart_axis] * dUCart_surf[cart_axis][k]

        epsUUU = self._epsUUU  # already built from sqrt(det gamma) on-surface

        phiU = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    phiU[i] += epsUUU[i][j][k] * self._sD[j] * Dzeta_dD[k]

        one_over_8pi = sp.Integer(1) / (8 * sp.pi)
        Jdens = sp.sympify(0)
        for i in range(3):
            Jdens += self._omegaD[i] * phiU[i]
        return one_over_8pi * Jdens, nU_unit

    # -------------------- optional Omega helpers (no Christoffels built internally) --------------------

    def Omega_scalar_from_covDomega(self, covD_omegaD: List[List[sp.Expr]]) -> sp.Expr:
        """
        Build the scalar ``Omega = epsilon^{AB} D_A omega_B`` from a supplied covariant derivative.

        The 3D-to-2D projection used is
        ``Omega = epsilon^{ijk} s_i q_j^p q_k^q (nabla_p omega_q)``.

        :param covD_omegaD: 3×3 nested list ``(p,q)`` with entries ``nabla_p omega_q``.
        :return: SymPy expression for ``Omega``.

        Doctests:
        TBD
        """
        Omega = sp.sympify(0)
        # epsilon^{ijk} s_i q_j^p q_k^q (nabla_p omega_q)
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

    def omega_moment_density_from_zeta(
        self, zeta: sp.Expr, covD_omegaD: List[List[sp.Expr]]
    ) -> sp.Expr:
        """
        Build the moment-of-vorticity density ``(1/(8*pi)) * zeta * Omega``.

        ``Omega`` is computed via :py:meth:`Omega_scalar_from_covDomega`.

        :param zeta: Scalar potential (e.g., a Cartesian scalar ``X^a`` or ``n·X``).
        :param covD_omegaD: 3×3 nested list representing ``nabla_p omega_q``.
        :return: SymPy expression equal to ``(1/(8*pi)) * zeta * Omega``.

        Doctests:
        TBD
        """
        Omega = self.Omega_scalar_from_covDomega(covD_omegaD)
        return (sp.Integer(1) / (8 * sp.pi)) * zeta * Omega


class HorizonSpinVorticityDipoleClass_dict(Dict[str, HorizonSpinVorticityDipoleClass]):
    """
    Dictionary-like accessor.

    * "Spherical" -> reference metric without precompute
    * "Spherical_rfm_precompute" -> reference metric with precompute enabled

    Note:
        The on-surface build is independent of precompute; the key only alters
        which ``rfm`` entry provides coordinate symbols.
    """

    def __getitem__(self, key: str) -> HorizonSpinVorticityDipoleClass:
        if key not in self:
            if key == "Spherical":
                self.__setitem__(
                    key,
                    HorizonSpinVorticityDipoleClass(
                        CoordSystem="Spherical", enable_rfm_precompute=False
                    ),
                )
            elif key == "Spherical_rfm_precompute":
                self.__setitem__(
                    key,
                    HorizonSpinVorticityDipoleClass(
                        CoordSystem="Spherical", enable_rfm_precompute=True
                    ),
                )
            else:
                raise KeyError(
                    "Supported keys: 'Spherical', 'Spherical_rfm_precompute'."
                )
        return dict.__getitem__(self, key)

    def __setitem__(self, key: str, value: HorizonSpinVorticityDipoleClass) -> None:
        if key not in ["Spherical", "Spherical_rfm_precompute"]:
            raise KeyError("Supported keys: 'Spherical', 'Spherical_rfm_precompute'.")
        dict.__setitem__(self, key, value)


# Instantiate accessor
HorizonSpinVorticityDipole = HorizonSpinVorticityDipoleClass_dict()


# -------------------- validation harness (exports EXACTLY 3 entries) --------------------
if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.validate_expressions.validate_expressions as ve

    # Doctests
    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Export exactly the 3 J components
    for validation_key in ["Spherical", "Spherical_rfm_precompute"]:
        modobj = HorizonSpinVorticityDipole[validation_key]

        export_only: Dict[str, sp.Expr] = {}
        for axis_idx in range(3):
            export_only[f"JCart_densU_{axis_idx}"] = modobj.JCart_densU[axis_idx]

        # Process and compare/generate trusted results for exactly these 3 keys
        results_dict = ve.process_dictionary_of_expressions(
            export_only, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{validation_key}",
            results_dict,
        )
