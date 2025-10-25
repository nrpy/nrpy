# -*- coding: utf-8 -*-
r"""
Omega-based quasilocal spin diagnostics on an apparent horizon surface.

Overview
--------
Constructs symbolic expressions for spin_vector(), as described here:
https://spectre-code.org/group__SurfacesGroup.html#ga03595bdaf4f20da98d151470998a22bf

Specifically, this module builds the intrinsic two-geometry on a horizon surface
``r = h(θ, φ)`` and assembles the single-pass **integrands** needed by an
Ω-based quasilocal spin diagnostic. It also provides reduction helpers that
convert accumulated surface integrals ("RunSums") into a spin **vector** with a
well-specified near-zero policy.


What this module constructs
---------------------------
Per surface point (no integration performed here):

* ``area_integrand`` (= 1) and ``area_density`` (= ``sqrt(q)``).
* ``ricci_scalar`` (= R) via the **intrinsic** 2D Christoffel route.
* ``spin_function`` (= Ω) using
  ``Ω = ε^{AB}(∂_A X_B − Γ^C_{AB} X_C)`` with ``X_B = e_B^i K_{ij} s^j``.
* Moment ingredients in an arbitrary **MeasurementFrame**:
  ``measurement_frame_xU[3]``, ``xR_momentU[3]`` (= ``x^i R``),
  ``xOmega_momentU[3]`` (= ``x^i Ω``), and ``zOmegaU[3]`` (= ``z_α Ω``).
* Eigen-operator helpers used by external solvers for the spin potential ``z``:
  ``laplacian_of_z`` (= ``Δz``), ``laplacian_of_y`` (= ``Δy``), and
  ``div_R_grad_z`` (= ``∇·(R∇z)``).
* Validation integrands: ``gauss_bonnet_integrand`` (= R) and
  ``omega_constraint_integrand`` (= Ω).

Integration policy (explicit & uniform)
---------------------------------------
**Nothing is premultiplied by ``sqrt(q)`` in this module.** To integrate a
per-point scalar ``f`` you must multiply by ``area_density = sqrt(q)`` and any
quadrature weights & angular steps externally:

``∮ f dA  ≃  Σ f_ab * area_density_ab * weights * Δθ * Δφ``

All **moment** integrands are **unweighted**: ``measurement_frame_xU``,
``xR_momentU``, ``xOmega_momentU``, and ``zOmegaU`` are not premultiplied by
``sqrt(q)``.

Frames
------
* **MetricDataFrame (MDF)**: where geometry lives (``γ_ij``, ``K_ij``, ``s^i``,
  ``e_A^i``, ``q_AB``, R, Ω).
* **MeasurementFrame (MF)**: used **only** for the position vector ``x^i`` that
  appears inside moment integrals; set via ``set_measurement_frame_coords(...)``.
  MF may differ from MDF.

Conventions & scope
-------------------
* Coordinates: ambient reference-metric infrastructure in spherical coordinates.
* Intrinsic (2D Christoffel) route for R; max derivative order: second.
* Orientation: ``ε^{θφ} = σ / sqrt(q)`` with configurable ``σ ∈ {+1, −1}``.
* A uniform substitution ``r→h`` (and ``f0_of_xx0→h``) is applied at the end of
  expression construction so all outputs are evaluated **on the surface**.

Outputs and reductions at a glance
----------------------------------
Use ``get_public_integrands()`` to retrieve all per-point fields. Then, after
performing your own quadratures, feed the resulting integrals into the helpers:

* ``reduce_centroids_and_direction(sums)`` → centroids ``x0U``, Ricci shift
  ``xRcorrU``, direction integral ``IU``, its norm ``normI``, and unit vector
  ``nU``.
* ``Salpha_and_magnitude_from_zOmega(zOmega_integrals)`` → ``S_α`` and ``S``
  with the ``1/(8π)`` normalization applied.
* ``magnitude_from_zOmega(zOmega_integrals)`` → legacy convenience returning
  ``S`` only.
* ``compute_spin_vector(sums, S, eps, ...)`` → final spin vector ``SU`` with a
  documented near-zero fallback (zero vector or ``z_α`` axis).
* ``validation_residuals_from_sums(sums)`` → cheap checks for
  Gauss–Bonnet (``∮R dA − 8π``) and the Ω-constraint (``∮Ω dA``).

Centroid note
-------------
This implementation uses the **area-weighted** centroid
``x0^i = (∮ x^i dA)/A``. Some external codes use an unweighted point-average of
collocation coordinates; do not mix the two when comparing results.

Authors: Ralston Graves
         Zachariah B. Etienne
         zachetie **at** gmail **dot* com
------
This file provides documentation-heavy wrappers around symbolic expressions; it
does not perform numerical integration or eigen-solves.
"""

from typing import Any, Dict, Iterable, List, Optional, Tuple, Union, cast

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.reference_metric as refmetric


def _sympify_int(x: int) -> sp.Integer:
    r"""
    Return the input as a SymPy ``Integer``.

    This helper ensures that integers used inside symbolic expressions are
    normalized to SymPy's own integer type. Doing so avoids accidental Python
    ``int`` propagation inside ``Piecewise`` and other SymPy constructs.

    :param x: A Python integer to normalize.
    :return: ``sympy.Integer`` representing ``x``.
    """
    return sp.Integer(int(sp.Integer(x)))


class SpECTRESpinEstimateClass:
    r"""
    Build Ω-based quasilocal spin diagnostic pieces on ``r = h(θ, φ)``.

    Parameters
    ----------
    CoordSystem : str, optional
        Ambient coordinate system. Only ``"Spherical"`` is supported.
    enable_rfm_precompute : bool, optional
        If ``True``, use the precompute-enabled reference-metric infrastructure.
        This affects only how coordinate symbols are provided; mathematics is
        unchanged.
    orientation_sign : int, optional
        Orientation σ entering ``ε^{θφ} = σ / sqrt(q)``. Choose ``+1`` (default)
        or ``−1``. Changing σ flips the sign of Ω and all Ω-weighted moments.

    Public methods (selected)
    -------------------------
    bind_surface_derivs(surface_derivs, *, enforce_shared_stencil=False)
        Validate & record external derivative placeholders for the 2D operators
        and Ω. Shapes and optional shared-stencil checks are enforced.
    set_measurement_frame_coords(xMeasU)
        Set the three MeasurementFrame coordinates ``x^i`` used **only** inside
        the moment integrands.
    set_orientation_sign(sign)
        Flip / set the sign convention used to define Ω.
    get_public_integrands()
        Return all per-point integrands and scalars (see module docstring).
    reduce_centroids_and_direction(sums)
        Single-pass RunSums reduction to centroids, Ricci shift, and direction.
    near_zero_policy(A, normI, integral_abs_Omega, eps)
        Dimensionless near-zero criterion as a SymPy inequality.
    Salpha_and_magnitude_from_zOmega(zOmega_integrals)
        Return ``{"SalphaU": [S0, S1, S2], "S": S}`` with ``1/(8π)`` factor.
    magnitude_from_zOmega(zOmega_integrals)
        Legacy convenience returning ``S`` only.
    compute_spin_vector(sums, S, eps, fallback_choice="zero", zOmega_integrals=None)
        Turnkey spin-vector constructor with near-zero fallback.
    validation_residuals_from_sums(sums)
        Cheap checks: ``{"gauss_bonnet": ∮R dA − 8π, "omega_constraint": ∮Ω dA}``.

    Notes
    -----
    * All expressions are constructed from fundamental variables without
      symbolic differentiation beyond second derivatives.
    * A final substitution ``r→h`` and ``f0_of_xx0→h`` is applied uniformly.
    * ``area_integrand`` is ``1`` by construction; always multiply by
      ``area_density`` externally when integrating.
    """

    # -------------------------------------------------------------------------
    # Construction
    # -------------------------------------------------------------------------
    def __init__(
        self,
        CoordSystem: str = "Spherical",
        enable_rfm_precompute: bool = False,
        orientation_sign: int = +1,
    ):
        if CoordSystem != "Spherical":
            raise ValueError(
                f"Unsupported CoordSystem '{CoordSystem}'. Only 'Spherical' is supported."
            )
        if int(orientation_sign) not in (+1, -1):
            raise ValueError("orientation_sign must be +1 or -1.")
        self.CoordSystem = CoordSystem
        self._orientation_sign = _sympify_int(orientation_sign)

        self._rfm = refmetric.reference_metric[
            (
                (self.CoordSystem + "_rfm_precompute")
                if enable_rfm_precompute
                else self.CoordSystem
            )
        ]

        # --- Ambient metric inputs (BSSN-like split) ---
        # Shape h and its angular derivatives on the surface
        self._h = sp.Symbol("hh", real=True)
        self._h_dD = ixp.declarerank1("hh_dD")  # components used: 1 (θ) and 2 (φ)
        self._h_dDD = ixp.declarerank2(
            "hh_dDD", symmetry="sym01"
        )  # angular second derivs

        # Conformal metric deformation h_ij and its FD derivatives
        self._hDD = ixp.declarerank2("hDD", symmetry="sym01")
        self._partial_D_hDD = cast(
            List[List[List[sp.Expr]]],
            ixp.declarerank3("partial_D_hDD", symmetry="sym12"),
        )

        # Conformal factor W and its FD gradient
        self._W = sp.Symbol("WW", real=True)
        self._WdD = ixp.declarerank1("partial_D_WW")

        # Tracefree A_ij piece and trK
        self._aDD = ixp.declarerank2("aDD", symmetry="sym01")
        self._trK = sp.Symbol("trK", real=True)

        # --- Build gamma_ij, gamma^ij, detgamma ---
        gammabarDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                gammabarDD[i][j] = (
                    self._hDD[i][j] * self._rfm.ReDD[i][j] + self._rfm.ghatDD[i][j]
                )
        self._gammaDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self._gammaDD[i][j] = gammabarDD[i][j] / (self._W**2)
        self._gammaUU, self._detgamma = ixp.symm_matrix_inverter3x3(self._gammaDD)

        # --- Level set and unit normal s^i in ambient coords (r,θ,φ) ---
        # F_,i = [1, -h_θ, -h_φ]
        F_dD = [sp.sympify(1), -self._h_dD[1], -self._h_dD[2]]
        nU = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                nU[i] += self._gammaUU[i][j] * F_dD[j]
        lambda_sq = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                lambda_sq += self._gammaUU[i][j] * F_dD[i] * F_dD[j]
        self._lamb = sp.sqrt(lambda_sq)
        self._sU = ixp.zerorank1()
        for i in range(3):
            self._sU[i] = nU[i] / self._lamb

        # --- Tangent basis vectors e_A^i (A=θ,φ; i=r,θ,φ) ---
        # e_θ^i = (h_θ, 1, 0), e_φ^i = (h_φ, 0, 1)
        self._eADU = [[sp.sympify(0) for _ in range(3)] for __ in range(2)]
        self._eADU[0][0] = self._h_dD[1]
        self._eADU[0][1] = sp.sympify(1)
        self._eADU[0][2] = sp.sympify(0)
        self._eADU[1][0] = self._h_dD[2]
        self._eADU[1][1] = sp.sympify(0)
        self._eADU[1][2] = sp.sympify(1)

        # --- Induced 2-metric q_AB, inverse, det, and sqrt(q) ---
        self._q2DD = ixp.zerorank2(dimension=2)
        for A in range(2):
            for B in range(2):
                for i in range(3):
                    for j in range(3):
                        self._q2DD[A][B] += (
                            self._gammaDD[i][j] * self._eADU[A][i] * self._eADU[B][j]
                        )
        self._q2UU, self._detq2 = ixp.symm_matrix_inverter2x2(self._q2DD)
        self._sqrtq = sp.sqrt(self._detq2)

        # --- Orientation tensor ε^{AB} (unit convention uses +1/√q); σ handled later ---
        self._eps2UU_unit = ixp.zerorank2(dimension=2)
        self._eps2UU_unit[0][1] = +1 / self._sqrtq
        self._eps2UU_unit[1][0] = -1 / self._sqrtq

        # --- Extrinsic curvature K_ij and spin one-form components X_B ---
        AbarDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                AbarDD[i][j] = self._aDD[i][j] * self._rfm.ReDD[i][j]
        self._KDD = ixp.zerorank2()
        exp4phi = 1 / (self._W**2)
        for i in range(3):
            for j in range(3):
                self._KDD[i][j] = (
                    exp4phi * AbarDD[i][j]
                    + sp.Rational(1, 3) * self._gammaDD[i][j] * self._trK
                )

        # X_B = e_B^i K_ij s^j
        self._XB = [sp.sympify(0), sp.sympify(0)]
        for B in range(2):
            for i in range(3):
                for j in range(3):
                    self._XB[B] += self._eADU[B][i] * self._KDD[i][j] * self._sU[j]

        # --- Declare all surface-derivative inputs; bind in bind_surface_derivs() ---
        # Derivatives of q_AB on the surface
        self._q2DD_dD = ixp.declarerank3("q2DD_dD", dimension=2)  # q_AB,C
        self._q2DD_dDD = ixp.declarerank4(
            "q2DD_dDD", symmetry="sym01_sym23", dimension=2
        )  # q_AB,CD

        # Raw partial of X_B on the surface
        self._XB_dD = ixp.declarerank2("XB_dD", dimension=2)  # (X_B)_,A

        # Eigen-operator derivative inputs
        self._zeta_dD = ixp.declarerank1("zeta_dD", dimension=2)  # z_,A
        self._zeta_dDD = ixp.declarerank2(
            "zeta_dDD", symmetry="sym01", dimension=2
        )  # z_,AB
        self._y_aux_dD = ixp.declarerank1("y_aux_dD", dimension=2)  # y_,A
        self._y_aux_dDD = ixp.declarerank2(
            "y_aux_dDD", symmetry="sym01", dimension=2
        )  # y_,AB

        # Conservative divergence flux-density derivative: d_A[ sqrt(q) R q^{AB} z_,B ]
        self._Rgradz_flux_density_dD = ixp.declarerank1(
            "Rgradz_flux_density_dD", dimension=2
        )

        # Three scalar test modes z_alpha for zOmegaU
        self._zU = ixp.declarerank1("zU", dimension=3)

        # --- Measurement-frame x^i for moments (defaults to Cartesian coords of the point) ---
        X = self._rfm.xx_to_Cart[0]
        Y = self._rfm.xx_to_Cart[1]
        Z = self._rfm.xx_to_Cart[2]
        self._xMeasU = [X, Y, Z]

        # --- Build all 2D connections, curvature, Ω_base (σ=+1), and helpers ---
        self._build_intrinsic_ops_and_omega_base()
        self._xR_momentU = [self._xMeasU[i] * self._R for i in range(3)]

        # --- Apply the horizon substitution uniformly (r -> h and f0_of_xx0 -> h) ---
        self._apply_surface_substitutions_all()

        # --- Apply orientation sign σ to define Ω and rebuild dependent moments ---
        self._update_omega_with_orientation()

        # Record of bound (external) surface derivatives for diagnostics
        self._bound_surface_derivs: Optional[Dict[str, object]] = None

    # -------------------------------------------------------------------------
    # Public: bind inputs provided on the surface by FD kernels (validate & record)
    # -------------------------------------------------------------------------
    def bind_surface_derivs(
        self, surface_derivs: Dict[str, object], *, enforce_shared_stencil: bool = False
    ) -> None:
        r"""
        Validate and record horizon-surface derivative placeholders.

        This function **does not** alter already-built symbolic expressions; it
        only checks and records the provided derivative objects for auditing and
        tooling, enforcing shapes and (optionally) a shared stencil ID.

        Required keys in ``surface_derivs`` (fixed shapes implied by declarations):
          - ``"q2DD_dD"``                 (2×2×2)    ``q_{AB,C}``
          - ``"q2DD_dDD"``                (2×2×2×2)  ``q_{AB,CD}``
          - ``"XB_dD"``                   (2×2)      ``(X_B)_,A``
          - ``"zeta_dD"``, ``"zeta_dDD"`` (2), (2×2) derivatives for ``z``
          - ``"y_aux_dD"``, ``"y_aux_dDD"`` (2), (2×2) derivatives for ``y``
          - ``"Rgradz_flux_density_dD"``  (2)        ``∂_A[ √q R q^{AB} z_,B ]``
          - ``"zU"``                      (3)        three scalar test modes ``z_α``

        :param surface_derivs: Dictionary of externally computed surface
            derivatives.
        :param enforce_shared_stencil: If ``True``, validate that all derivative
            inputs share a common ``stencil_id`` (either global or per-array attribute).
        :raises KeyError: A required key is missing from ``surface_derivs``.
        :raises ValueError: A derivative object has an incorrect shape, or
            stencil IDs are inconsistent when ``enforce_shared_stencil=True``.
        """
        required = [
            "q2DD_dD",
            "q2DD_dDD",
            "XB_dD",
            "zeta_dD",
            "zeta_dDD",
            "y_aux_dD",
            "y_aux_dDD",
            "Rgradz_flux_density_dD",
            "zU",
        ]
        for key in required:
            if key not in surface_derivs:
                raise KeyError(f"Missing required surface derivative input: {key}")

        # Minimal shape validation helpers
        def _len1(x: Any) -> int:
            return len(x)

        def _len2(x: Any) -> Tuple[int, int]:
            return (len(x), len(x[0]))

        def _len3(x: Any) -> Tuple[int, int, int]:
            return (len(x), len(x[0]), len(x[0][0]))

        def _len4(x: Any) -> Tuple[int, int, int, int]:
            return (len(x), len(x[0]), len(x[0][0]), len(x[0][0][0]))

        # Shapes (expect (2), (2×2), (2×2×2), (2×2×2×2), and (3))
        if _len3(surface_derivs["q2DD_dD"]) != (2, 2, 2):
            raise ValueError("q2DD_dD must have shape (2,2,2)")
        if _len4(surface_derivs["q2DD_dDD"]) != (2, 2, 2, 2):
            raise ValueError("q2DD_dDD must have shape (2,2,2,2)")
        if _len2(surface_derivs["XB_dD"]) != (2, 2):
            raise ValueError("XB_dD must have shape (2,2)")
        if _len1(surface_derivs["zeta_dD"]) != 2 or _len2(
            surface_derivs["zeta_dDD"]
        ) != (2, 2):
            raise ValueError("zeta_dD must be (2); zeta_dDD must be (2,2)")
        if _len1(surface_derivs["y_aux_dD"]) != 2 or _len2(
            surface_derivs["y_aux_dDD"]
        ) != (2, 2):
            raise ValueError("y_aux_dD must be (2); y_aux_dDD must be (2,2)")
        if _len1(surface_derivs["Rgradz_flux_density_dD"]) != 2:
            raise ValueError("Rgradz_flux_density_dD must be (2)")
        if _len1(surface_derivs["zU"]) != 3:
            raise ValueError("zU must be (3)")

        # Enforce shared stencil ID if provided
        if enforce_shared_stencil:
            # Accept either a global "stencil_id" or per-array "stencil_id" fields.
            stencil_id = surface_derivs.get("stencil_id", None)
            to_check = [
                "q2DD_dD",
                "q2DD_dDD",
                "XB_dD",
                "zeta_dD",
                "zeta_dDD",
                "y_aux_dD",
                "y_aux_dDD",
                "Rgradz_flux_density_dD",
            ]
            # If per-array attributes exist, require equality:
            found_ids = []
            for key in to_check:
                arr = surface_derivs[key]
                sid = getattr(arr, "stencil_id", stencil_id)
                if sid is not None:
                    found_ids.append(sid)
            if len(found_ids) > 0 and any(s != found_ids[0] for s in found_ids):
                raise ValueError(
                    "Shared-stencil requirement violated: non-uniform stencil_id across inputs."
                )

        # Record only (for tooling & audits)
        self._bound_surface_derivs = dict(surface_derivs)

    # -------------------------------------------------------------------------
    # Public: set measurement-frame coordinates x^i used in moment integrands
    # -------------------------------------------------------------------------
    def set_measurement_frame_coords(self, xMeasU: List[sp.Expr]) -> None:
        r"""
        Set measurement-frame coordinates ``x^i`` for moment integrands.

        This affects ``measurement_frame_xU``, ``xR_momentU``, and
        ``xOmega_momentU`` (all unweighted by design).

        :param xMeasU: Three expressions representing the MeasurementFrame
            coordinates ``x^i`` of the surface point.
        :raises ValueError: ``xMeasU`` does not have length 3.
        """
        if len(xMeasU) != 3:
            raise ValueError("xMeasU must have length 3.")
        self._xMeasU = list(xMeasU)
        # Rebuild dependent (unweighted) moment outputs:
        self._xR_momentU = [self._xMeasU[i] * self._R for i in range(3)]
        # Ω-dependent moments updated through `_update_omega_with_orientation()` to keep σ in sync:
        self._update_omega_with_orientation()

    # -------------------------------------------------------------------------
    # Public: set orientation sign σ used in ε^{θφ} = σ / sqrt(q)
    # -------------------------------------------------------------------------
    def set_orientation_sign(self, sign: int) -> None:
        r"""
        Set the surface-orientation sign ``σ ∈ {+1, −1}``.

        This affects Ω and all Ω-weighted moments. Use this to switch between
        alternate conventions on the 2-surface.

        :param sign: The orientation sign, ``+1`` or ``−1``.
        :raises ValueError: ``sign`` is not ``+1`` or ``−1``.
        """
        if int(sign) not in (+1, -1):
            raise ValueError("orientation_sign must be +1 or -1.")
        self._orientation_sign = _sympify_int(sign)
        self._update_omega_with_orientation()

    # -------------------------------------------------------------------------
    # Public: get all per-point outputs needed by the diagnostic sweep
    # -------------------------------------------------------------------------
    def get_public_integrands(self) -> Dict[str, object]:
        r"""
        Return all per-point integrands and scalars.

        Integration policy: **No field is premultiplied by ``sqrt(q)``**. Multiply
        by ``area_density`` and quadrature weights externally when forming
        surface integrals.

        :return: Dictionary with the following keys (values are SymPy expressions
            or lists of expressions):
            * ``area_integrand``: unit scalar ``1`` (by policy).
            * ``area_density``: ``sqrt(q)``.
            * ``ricci_scalar``: the 2D Ricci scalar ``R``.
            * ``spin_function``: the spin function ``Ω``.
            * ``measurement_frame_xU``: list[3] of MeasurementFrame coordinates.
            * ``xR_momentU``: list[3] with components ``x^i * R``.
            * ``xOmega_momentU``: list[3] with components ``x^i * Ω``.
            * ``zOmegaU``: list[3] with components ``z_α * Ω``.
            * ``laplacian_of_z``: ``Δz`` (for eigen-operators).
            * ``laplacian_of_y``: ``Δy`` (auxiliary mode).
            * ``div_R_grad_z``: conservative ``∇·(R∇z)`` divided by ``sqrt(q)``.
            * ``gauss_bonnet_integrand``: same as ``R``.
            * ``omega_constraint_integrand``: same as ``Ω``.
        """
        out: Dict[str, object] = {}
        out["area_integrand"] = sp.sympify(1)  # unit scalar by policy
        out["area_density"] = self._sqrtq
        out["ricci_scalar"] = self._R
        out["spin_function"] = self._Omega
        out["measurement_frame_xU"] = self._xMeasU
        out["xR_momentU"] = self._xR_momentU
        out["xOmega_momentU"] = self._xOmega_momentU
        out["zOmegaU"] = self._zOmegaU
        out["laplacian_of_z"] = self._laplacian_of_z
        out["laplacian_of_y"] = self._laplacian_of_y
        out["div_R_grad_z"] = self._div_R_grad_z
        out["gauss_bonnet_integrand"] = self._R
        out["omega_constraint_integrand"] = self._Omega
        return out

    # -------------------------------------------------------------------------
    # Public: RunSums single-pass reduction to centroids & direction
    # -------------------------------------------------------------------------
    def reduce_centroids_and_direction(
        self, sums: Dict[str, object]
    ) -> Dict[str, object]:
        r"""
        Reduce single-pass RunSums to centroids and direction vector.

        Expect ``sums`` in either RunSums style
        ``{"A", "XU"[3], "R0", "XRU"[3], "O0", "XOU"[3], (optional) "Oabs"}``
        or in the legacy flat style with ``XU0..2``, ``XRU0..2``, ``XOU0..2``.

        :param sums: Integrated quantities (RunSums dictionary).
        :return: ``{"x0U", "xRcorrU", "IU", "normI", "nU"}`` where
            * ``x0U`` is the area-weighted centroid ``(∮ x^i dA)/A``.
            * ``xRcorrU`` are the Ricci-shifted centers
              ``(XR^i − x0^i R0)/(8π)``.
            * ``IU`` is the direction integral
              ``∮ Ω (x^i − x0^i − x_R^i) dA``.
            * ``normI`` is the Euclidean norm of ``IU``.
            * ``nU`` are components of the unit vector ``IU/‖IU‖`` with a
              safe ``Piecewise`` when ``‖IU‖=0``.
        """
        # Normalize input to RunSums form
        if "XU" in sums:
            A = sums["A"]
            XU = list(cast(Iterable[Any], sums["XU"]))
            R0 = sums["R0"]
            XRU = list(cast(Iterable[Any], sums["XRU"]))
            O0 = sums["O0"]
            XOU = list(cast(Iterable[Any], sums["XOU"]))
        else:
            A = sums["A"]
            XU = [sums["XU0"], sums["XU1"], sums["XU2"]]
            R0 = sums["R0"]
            XRU = [sums["XRU0"], sums["XRU1"], sums["XRU2"]]
            O0 = sums["O0"]
            XOU = [sums["XOU0"], sums["XOU1"], sums["XOU2"]]

        x0U = [XU[i] / A for i in range(3)]

        # x_R^i = (XR^i - x0^i * R0) / (8π)
        xRcorrU = [(XRU[i] - x0U[i] * R0) / (8 * sp.pi) for i in range(3)]

        # I^i = XO^i - (x0^i + x_R^i) * O0
        IU = [XOU[i] - (x0U[i] + xRcorrU[i]) * O0 for i in range(3)]

        # Euclidean norm in measurement frame
        normI = sp.sqrt(IU[0] * IU[0] + IU[1] * IU[1] + IU[2] * IU[2])

        # n^i = I^i / ‖I‖ (left symbolic; policy handles near-zero)
        nU = [
            sp.Piecewise((sp.sympify(0), sp.Eq(normI, 0)), (IU[i] / normI, True))
            for i in range(3)
        ]

        return {"x0U": x0U, "xRcorrU": xRcorrU, "IU": IU, "normI": normI, "nU": nU}

    # -------------------------------------------------------------------------
    # Public: robustness policy near zero spin-direction signal
    # -------------------------------------------------------------------------
    def near_zero_policy(
        self,
        A: sp.Expr,
        normI: sp.Expr,
        integral_abs_Omega: sp.Expr,
        eps: Union[float, sp.Expr],
    ) -> Dict[str, sp.Expr]:
        r"""
        Implement the dimensionless near-zero spin policy.

        The policy triggers when ``‖I‖ ≤ ε * R_char * ∮|Ω| dA``, with
        ``R_char = sqrt(A/(4π))``.

        :param A: Total surface area ``∮ dA``.
        :param normI: Euclidean norm of the direction integral ``I``.
        :param integral_abs_Omega: Integral of the absolute value ``∮|Ω| dA``.
        :param eps: Dimensionless tolerance ``ε``.
        :return: ``{"trigger": <inequality>}``, where the value is a SymPy
            boolean suitable for use in ``Piecewise``.
        """
        R_char = sp.sqrt(A / (4 * sp.pi))
        condition = sp.Le(normI, sp.sympify(eps) * R_char * integral_abs_Omega)
        return {"trigger": condition}

    # -------------------------------------------------------------------------
    # Public: compute S_α and S from ∮ z_α Ω dA (both returned)
    # -------------------------------------------------------------------------
    def Salpha_and_magnitude_from_zOmega(
        self, zOmega_integrals: List[sp.Expr]
    ) -> Dict[str, object]:
        r"""
        Compute ``S_α`` and ``S`` from the integrals ``∮ z_α Ω dA``.

        Given the three surface integrals ``ZΩ_α := ∮ z_α Ω dA``, this function
        returns ``S_α = (1/(8π)) ZΩ_α`` and the magnitude
        ``S = sqrt(Σ_α S_α²)``.

        :param zOmega_integrals: ``[∮ z_0 Ω dA, ∮ z_1 Ω dA, ∮ z_2 Ω dA]``.
        :return: ``{"SalphaU": [S_0, S_1, S_2], "S": S}``.
        :raises ValueError: ``zOmega_integrals`` does not have length 3.
        """
        if len(zOmega_integrals) != 3:
            raise ValueError("zOmega_integrals must contain three integrals.")
        SalphaU = [zOmega_integrals[a] / (8 * sp.pi) for a in range(3)]
        S = sp.sqrt(SalphaU[0] ** 2 + SalphaU[1] ** 2 + SalphaU[2] ** 2)
        return {"SalphaU": SalphaU, "S": S}

    def magnitude_from_zOmega(self, zOmega_integrals: List[sp.Expr]) -> sp.Expr:
        r"""
        Compute spin magnitude ``S`` from the integrals ``∮ z_α Ω dA`` (legacy).

        ``S = (1/(8π)) * sqrt( \sum_α (ZΩ_α)^2 )``.

        :param zOmega_integrals: ``[∮ z_0 Ω dA, ∮ z_1 Ω dA, ∮ z_2 Ω dA]``.
        :return: The spin magnitude ``S``.
        :raises ValueError: ``zOmega_integrals`` does not have length 3.
        """
        if len(zOmega_integrals) != 3:
            raise ValueError("zOmega_integrals must contain three integrals.")
        ssum = sp.sympify(0)
        for a in range(3):
            ssum += zOmega_integrals[a] * zOmega_integrals[a]
        return sp.sqrt(ssum) / (8 * sp.pi)

    # -------------------------------------------------------------------------
    # Public: turnkey spin-vector constructor with near-zero fallback
    # -------------------------------------------------------------------------
    def compute_spin_vector(
        self,
        sums: Dict[str, object],
        S: sp.Expr,
        eps: Union[float, sp.Expr],
        fallback_choice: str = "zero",
        zOmega_integrals: Optional[List[sp.Expr]] = None,
    ) -> Dict[str, object]:
        r"""
        Compose reductions and policies to return a spin vector ``S^i``.

        :param sums: RunSums dictionary (see ``reduce_centroids_and_direction``).
            May optionally contain ``"Oabs" = ∮|Ω| dA`` for a concrete trigger.
        :param S: Spin magnitude (dimensionful or dimensionless). The returned
            vector has Euclidean norm ``S`` by construction, unless the near-zero
            policy triggers and the fallback returns the zero vector.
        :param eps: Tolerance ``ε`` in the near-zero criterion.
        :param fallback_choice: One of ``{"zero", "zalpha"}``. If ``"zero"``,
            return the zero vector when the policy triggers. If ``"zalpha"``,
            fall back to the direction defined by ``S_α = (1/8π) ∮ z_α Ω dA``.
        :param zOmega_integrals: Required when ``fallback_choice='zalpha'``;
            ignored otherwise.
        :return: ``{"SU": [Sx, Sy, Sz], "near_zero_trigger": <bool>,
            "x0U", "xRcorrU", "IU", "nU", "normI"}``.
        :raises ValueError: ``fallback_choice`` is not one of ``{"zero","zalpha"}``.
        """
        red = self.reduce_centroids_and_direction(sums)

        # Optional absolute-Ω integral for the policy
        Oabs = sums.get("Oabs", None)
        if Oabs is None:
            # Keep symbolic placeholder to make the Piecewise explicit.
            Oabs = sp.Symbol("integral_abs_Omega", real=True, nonnegative=True)

        policy = self.near_zero_policy(
            A=sums["A"], normI=red["normI"], integral_abs_Omega=Oabs, eps=eps
        )
        trigger = policy["trigger"]

        # Nominal direction
        nU = cast(List[sp.Expr], red["nU"])
        nominal_SU = [S * nU[i] for i in range(3)]

        # Fallback
        if fallback_choice not in ("zero", "zalpha"):
            raise ValueError("fallback_choice must be one of {'zero', 'zalpha'}.")

        if (
            fallback_choice == "zalpha"
            and zOmega_integrals is not None
            and len(zOmega_integrals) == 3
        ):
            smag = self.Salpha_and_magnitude_from_zOmega(zOmega_integrals)
            SalphaU = cast(List[sp.Expr], smag["SalphaU"])  # already 1/(8π) scaled
            normSalpha = sp.sqrt(SalphaU[0] ** 2 + SalphaU[1] ** 2 + SalphaU[2] ** 2)
            # Unit from S_α; handle norm=0 safely:
            nz_dir = [
                sp.Piecewise(
                    (sp.sympify(0), sp.Eq(normSalpha, 0)),
                    (SalphaU[i] / normSalpha, True),
                )
                for i in range(3)
            ]
            fallback_SU = [S * nz_dir[i] for i in range(3)]
        else:
            fallback_SU = [sp.sympify(0), sp.sympify(0), sp.sympify(0)]

        # Piecewise assemble
        SU = [
            sp.Piecewise((fallback_SU[i], trigger), (nominal_SU[i], True))
            for i in range(3)
        ]

        out = {"SU": SU, "near_zero_trigger": trigger, "normI": red["normI"]}
        out.update({k: red[k] for k in ("x0U", "xRcorrU", "IU", "nU")})
        return out

    # -------------------------------------------------------------------------
    # Public: validation residuals using RunSums (cheap, consistent operators)
    # -------------------------------------------------------------------------
    def validation_residuals_from_sums(
        self, sums: Dict[str, object]
    ) -> Dict[str, sp.Expr]:
        r"""
        Return cheap integral-level validation residuals from RunSums.

        Uses the same surface-derivative operators that feed ``R`` and ``Ω``.

        :param sums: Integrated quantities (see ``reduce_centroids_and_direction``).
        :return: ``{"gauss_bonnet": ∮R dA − 8π, "omega_constraint": ∮Ω dA}``.
        """
        R0 = sums["R0"]
        O0 = sums["O0"]
        return {"gauss_bonnet": R0 - 8 * sp.pi, "omega_constraint": O0}

    # =========================================================================
    # Internal: build 2D connections, Ricci scalar, Ω_base (σ=+1), eigen-ops
    # =========================================================================
    def _build_intrinsic_ops_and_omega_base(self) -> None:
        r"""
        Construct intrinsic geometric operators and ``Ω_base``.

        Builds (i) derivatives of ``q^{AB}``, (ii) 2D Christoffel symbols and
        their derivatives, (iii) the 2D Ricci scalar ``R``, (iv) the base spin
        function ``Ω_base`` with unit orientation (``σ=+1``), and (v) Laplacians
        and conservative divergence terms used by eigen-operators.
        """
        # q^{AB} derivatives from identity: (qUU)_,C = -q^{AE} q^{BF} q_{EF,C}
        self._q2UUdD = ixp.zerorank3(dimension=2)
        for A in range(2):
            for B in range(2):
                for C in range(2):
                    temp = sp.sympify(0)
                    for E in range(2):
                        for F in range(2):
                            temp += (
                                -self._q2UU[A][E]
                                * self._q2UU[B][F]
                                * self._q2DD_dD[E][F][C]
                            )
                    self._q2UUdD[A][B][C] = temp

        # 2D Christoffel symbols: Γ^C_{AB} = 1/2 q^{CD} (q_{BD,A} + q_{AD,B} - q_{AB,D})
        self._GammaU2DD = ixp.zerorank3(dimension=2)
        for C in range(2):
            for A in range(2):
                for B in range(2):
                    val = sp.sympify(0)
                    for D in range(2):
                        val += (
                            sp.Rational(1, 2)
                            * self._q2UU[C][D]
                            * (
                                self._q2DD_dD[B][D][A]
                                + self._q2DD_dD[A][D][B]
                                - self._q2DD_dD[A][B][D]
                            )
                        )
                    self._GammaU2DD[C][A][B] = val

        # Derivatives of Γ using product rule with qUUdD and q2DD_dDD:
        self._GammaU2DD_dD = ixp.zerorank4(dimension=2)
        for E in range(2):
            for C in range(2):
                for A in range(2):
                    for B in range(2):
                        term = sp.sympify(0)
                        for D in range(2):
                            bracket = (
                                self._q2DD_dD[B][D][A]
                                + self._q2DD_dD[A][D][B]
                                - self._q2DD_dD[A][B][D]
                            )
                            d_bracket = (
                                self._q2DD_dDD[B][D][A][E]
                                + self._q2DD_dDD[A][D][B][E]
                                - self._q2DD_dDD[A][B][D][E]
                            )
                            term += sp.Rational(1, 2) * (
                                self._q2UUdD[C][D][E] * bracket
                                + self._q2UU[C][D] * d_bracket
                            )
                        self._GammaU2DD_dD[C][A][B][E] = term

        # 2D Ricci scalar:
        # R = q^{AB} * ( ∂_C Γ^C_{AB} − ∂_B Γ^C_{AC} + Γ^C_{AB} Γ^D_{CD} − Γ^C_{AD} Γ^D_{BC} )
        dC_GammaCAB = ixp.zerorank2(dimension=2)  # AB
        dB_GammaCAC = ixp.zerorank2(dimension=2)  # AB with B acting on Γ^C_{AC}
        for A in range(2):
            for B in range(2):
                sumC = sp.sympify(0)
                sumB = sp.sympify(0)
                for C in range(2):
                    sumC += self._GammaU2DD_dD[C][A][B][C]
                    sumB += self._GammaU2DD_dD[C][A][C][B]
                dC_GammaCAB[A][B] = sumC
                dB_GammaCAC[A][B] = sumB

        GammaCAB_GammaDCD = ixp.zerorank2(dimension=2)
        GammaCAD_GammaDBC = ixp.zerorank2(dimension=2)
        for A in range(2):
            for B in range(2):
                g1 = sp.sympify(0)
                g2 = sp.sympify(0)
                for C in range(2):
                    for D in range(2):
                        g1 += self._GammaU2DD[C][A][B] * self._GammaU2DD[D][C][D]
                        g2 += self._GammaU2DD[C][A][D] * self._GammaU2DD[D][B][C]
                GammaCAB_GammaDCD[A][B] = g1
                GammaCAD_GammaDBC[A][B] = g2

        R_raw = ixp.zerorank2(dimension=2)
        for A in range(2):
            for B in range(2):
                R_raw[A][B] = (
                    dC_GammaCAB[A][B]
                    - dB_GammaCAC[A][B]
                    + GammaCAB_GammaDCD[A][B]
                    - GammaCAD_GammaDBC[A][B]
                )
        self._R = sp.sympify(0)
        for A in range(2):
            for B in range(2):
                self._R += self._q2UU[A][B] * R_raw[A][B]

        # Ω_base (σ=+1): ε_unit^{AB} * ( X_B,_A − Γ^C_{AB} X_C )
        Omega_sum = sp.sympify(0)
        for A in range(2):
            for B in range(2):
                covXB = self._XB_dD[B][A]
                for C in range(2):
                    covXB += -self._GammaU2DD[C][A][B] * self._XB[C]
                Omega_sum += self._eps2UU_unit[A][B] * covXB
        self._Omega_base = Omega_sum  # σ factor applied later

        # zΩ and xΩ moments are rebuilt with σ in `_update_omega_with_orientation()`.

        # Laplacians: q^{AB} ( z_,AB − Γ^C_{AB} z_,C ) and same for y
        self._laplacian_of_z = sp.sympify(0)
        self._laplacian_of_y = sp.sympify(0)
        for A in range(2):
            for B in range(2):
                term_z = self._zeta_dDD[A][B]
                term_y = self._y_aux_dDD[A][B]
                for C in range(2):
                    term_z += -self._GammaU2DD[C][A][B] * self._zeta_dD[C]
                    term_y += -self._GammaU2DD[C][A][B] * self._y_aux_dD[C]
                self._laplacian_of_z += self._q2UU[A][B] * term_z
                self._laplacian_of_y += self._q2UU[A][B] * term_y

        # Conservative divergence: (1/sqrt(q)) * ∂_A[ sqrt(q) R q^{AB} z_,B ]
        flux_div = sp.sympify(0)
        for A in range(2):
            flux_div += self._Rgradz_flux_density_dD[A]
        self._div_R_grad_z = flux_div / self._sqrtq

    # -------------------------------------------------------------------------
    # Internal: uniformly apply r->h and f0_of_xx0->h to stored expressions
    # -------------------------------------------------------------------------
    def _apply_surface_substitutions_all(self) -> None:
        """Apply the on-surface substitutions uniformly to cached expressions."""
        subs_map = {
            self._rfm.xx[0]: self._h,
            sp.sympify("f0_of_xx0"): self._h,
        }

        def apply(expr: Union[sp.Expr, List[Any]]) -> Union[sp.Expr, List[Any]]:
            if isinstance(expr, list):
                return [apply(e) for e in expr]
            return sp.sympify(expr).subs(subs_map)

        # Scalars
        self._detgamma = apply(self._detgamma)
        self._detq2 = apply(self._detq2)
        self._sqrtq = apply(self._sqrtq)
        self._R = apply(self._R)
        self._Omega_base = apply(self._Omega_base)
        self._laplacian_of_z = apply(self._laplacian_of_z)
        self._laplacian_of_y = apply(self._laplacian_of_y)
        self._div_R_grad_z = apply(self._div_R_grad_z)

        # Vectors / lists
        self._xMeasU = apply(self._xMeasU)
        self._xR_momentU = apply(self._xR_momentU)

    # -------------------------------------------------------------------------
    # Internal: apply orientation sign σ to Ω and rebuild Ω-dependent moments
    # -------------------------------------------------------------------------
    def _update_omega_with_orientation(self) -> None:
        """Update ``Ω`` and all Ω-dependent moments to respect the current ``σ``."""
        # Ω = σ * Ω_base
        self._Omega = self._orientation_sign * self._Omega_base
        self._xOmega_momentU = [self._xMeasU[i] * self._Omega for i in range(3)]
        self._zOmegaU = [self._zU[a] * self._Omega for a in range(3)]

    # =========================================================================
    # Dict-style cached factory
    # =========================================================================


class SpECTRESpinEstimateClass_dict(Dict[str, SpECTRESpinEstimateClass]):
    r"""
    Dictionary-like accessor for :class:`SpECTRESpinEstimateClass` instances.

    Keys
    ----
    * ``"Spherical"`` → reference metric without precompute.
    * ``"Spherical_rfm_precompute"`` → reference metric with precompute enabled.

    Notes
    -----
    The key only alters which ``rfm`` entry provides coordinate symbols; the
    resulting symbolic expressions are otherwise identical.
    """

    def __getitem__(self, key: str) -> SpECTRESpinEstimateClass:
        if key not in self:
            if key == "Spherical":
                enable_rfm_precompute = False
            elif key == "Spherical_rfm_precompute":
                enable_rfm_precompute = True
            else:
                raise KeyError(
                    "Supported keys are 'Spherical' and 'Spherical_rfm_precompute'."
                )
            print(f"Setting up SpECTRESpinEstimateClass[{key}]...")
            self.__setitem__(
                key,
                SpECTRESpinEstimateClass(
                    CoordSystem="Spherical",
                    enable_rfm_precompute=enable_rfm_precompute,
                    orientation_sign=+1,  # default convention ε^{θφ}=+1/√q
                ),
            )
        return dict.__getitem__(self, key)

    def __setitem__(self, key: str, value: SpECTRESpinEstimateClass) -> None:
        if key not in ["Spherical", "Spherical_rfm_precompute"]:
            raise KeyError(
                "Supported keys are 'Spherical' and 'Spherical_rfm_precompute'."
            )
        dict.__setitem__(self, key, value)


# Public handle
SpECTRESpinEstimate = SpECTRESpinEstimateClass_dict()

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

    # Sanity validation of expressions for both factory options
    for validation_key in ["Spherical", "Spherical_rfm_precompute"]:
        omega_calc = SpECTRESpinEstimate[validation_key]
        results_dict = ve.process_dictionary_of_expressions(
            omega_calc.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{validation_key}",
            results_dict,
        )
