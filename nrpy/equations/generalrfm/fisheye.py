# nrpy/equations/generalrfm/fisheye.py
"""
N transition fisheye mapping and induced reference metric.

The map is radial: CartU[i] = lam(r) * xx[i] with lam(r) = rbar(r) / r.
This module constructs CartU, the Jacobian, the induced flat reference metric, and its
first and second derivatives in 3D raw Cartesian coordinates.

This file also includes detailed derivations to justify the tensor expressions used
in the implementation.

1. Raw radius

Given raw Cartesian coordinates xx[0], xx[1], xx[2], define

    r(xx) = sqrt( (xx[0])^2 + (xx[1])^2 + (xx[2])^2 )

and note the standard identity

    d r / d xx[j] = xx[j] / r.

2. N transition fisheye radius map

The N transition fisheye map is defined by:
- plateau stretch factors a_0, a_1, ..., a_N
- transition centers R_1, R_2, ..., R_N
- width parameters s_1, s_2, ..., s_N
- differences delta_a_i = a_{i-1} - a_i

Single transition kernel:

    G(r; R, s) =
        s / (2 * tanh(R / s)) *
        log( cosh((r + R)/s) / cosh((r - R)/s) ).

Unscaled radius map:

    rbar_unscaled(r) =
        a_N * r + sum_{i=1..N} (a_{i-1} - a_i) * G(r; R_i, s_i).

Scaled physical radius map:

    rbar(r) = c * rbar_unscaled(r).

3. Forward map

The physical Cartesian coordinates are obtained by radial rescaling:

    CartU[i] = (rbar(r) / r) * xx[i] = lam(r) * xx[i],

where

    lam(r) = rbar(r) / r.

4. Jacobian derivation

Differentiate CartU[i] = lam(r) * xx[i] with respect to xx[j]. Using
d r / d xx[j] = xx[j] / r and the product rule:

    d CartU[i] / d xx[j]
      = lam(r) * delta_ij + xx[i] * d lam / d xx[j]

and since lam depends on xx only through r:

    d lam / d xx[j] = (d lam / d r) * (d r / d xx[j]) = (d lam / d r) * xx[j] / r.

Therefore the Jacobian has the compact form used in the code:

    d CartU[i] / d xx[j]
      = lam * delta_ij + (d lam / d r) * (xx[i] * xx[j] / r).

The implementation stores the factor (d lam / d r) / r and multiplies by xx[i] * xx[j].

5. Induced reference metric and the A(r), B(r) decomposition

The physical space is flat with metric delta_mn in the physical Cartesian coordinates.
The induced metric in the raw coordinates is

    ghatDD[i][j] = delta_mn
                   * (d CartU[m] / d xx[i])
                   * (d CartU[n] / d xx[j]).

By spherical symmetry, any rank-2 tensor built from xx[i] can be decomposed as

    ghatDD[i][j] = A(r) * delta_ij + B(r) * xx[i] * xx[j].

For the fisheye map:

- Tangential (angular) scaling is set by the factor multiplying directions orthogonal
  to xx[i], which is lam(r). Thus the tangential eigenvalue is lam(r)^2, so

      A(r) = lam(r)^2 = (rbar(r) / r)^2.

- Radial scaling is determined by how radial line elements transform. Along the radial
  direction, drbar = (d rbar / d r) * dr, so the radial metric component is

      ghat_rr = (d rbar / d r)^2.

In the A,B decomposition, the radial eigenvalue is A(r) + B(r) * r^2, so

    A(r) + B(r) * r^2 = (d rbar / d r)^2,

which yields

    B(r) = ( (d rbar / d r)^2 - A(r) ) / r^2.

Defining u = rbar(r) and p = d u / d r, this becomes the algebraic form used in code:

    A(r) = u^2 / r^2
    B(r) = p^2 / r^2 - u^2 / r^4.

6. First derivatives of the metric

Write

    ghatDD[i][j] = A(r) * delta_ij + B(r) * xx[i] * xx[j].

Using
- d r / d xx[k] = xx[k] / r
- d(xx[i] * xx[j]) / d xx[k] = delta_ik * xx[j] + xx[i] * delta_jk

we obtain

    ghatDDdD[i][j][k]
      = (A'(r) / r) * xx[k] * delta_ij
        + (B'(r) / r) * xx[k] * xx[i] * xx[j]
        + B(r) * (delta_ik * xx[j] + xx[i] * delta_jk).

This is exactly the structure implemented for ghatDDdD.

7. Second derivatives of the metric

Differentiate ghatDDdD[i][j][k] with respect to xx[l]. A useful identity is

    d (xx[k] / r) / d xx[l] = delta_kl / r - xx[k] * xx[l] / r^3.

When differentiating (A'(r)/r) * xx[k] and (B'(r)/r) * xx[k], the terms combine into
compact coefficients

    coeff_A = A''/r^2 - A'/r^3
    coeff_B = B''/r^2 - B'/r^3

that multiply xx[k] * xx[l]. The remaining pieces come from differentiating the
explicit xx factors in xx[i] * xx[j] and the deltas. The final assembled expression
for ghatDDdDD is implemented directly in tensor-component form.

8. Closed-form radial derivatives

To avoid expensive Cartesian differentiation, this module builds all radial derivatives
in 1D in terms of a single symbol r_sym, and only at the end substitutes r_sym -> r(xx).

When requested, the kernel derivatives G', G'', G''' are computed in closed form, and
the map derivatives are accumulated in a single pass.

The helper _radial_metric_scalars uses u = rbar, p = du/dr, q = d2u/dr2, t = d3u/dr3
to write lam, dlam/dr, A, B and their first and second radial derivatives in explicit
algebraic form, which is then converted to full Cartesian expressions.
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple, Union

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.params as par

thismodule = __name__

__all__ = ["GeneralRFMFisheye", "make_fisheye_rfm"]

_CodeParameterDefault = Union[str, int, float]

_KroneckerDelta3D: Tuple[
    Tuple[int, int, int], Tuple[int, int, int], Tuple[int, int, int]
] = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
)


class GeneralRFMFisheye:
    """
    Build an N transition fisheye map and the induced flat reference metric in 3D.

    This class constructs and stores:
    * xx_to_CartU: the forward map CartU[i] = lam(r) * xx[i]
    * dCart_dxxUD: the Jacobian d CartU[i] / d xx[j]
    * ghatDD: the induced reference metric ghatDD[i][j]
    * ghatDDdD: the first derivatives ghatDDdD[i][j][k]
    * ghatDDdDD: the second derivatives ghatDDdDD[i][j][k][l]

    The derivations that justify the tensor expressions are documented in the module
    docstring and in helper docstrings.
    """

    def __init__(
        self,
        num_transitions: int,
        *,
        a_defaults: Optional[Sequence[float]] = None,
        R_defaults: Optional[Sequence[float]] = None,
        s_defaults: Optional[Sequence[float]] = None,
        c_default: float = 1.0,
        use_closed_form_derivs: bool = True,
        use_power_map_subs: bool = True,
        power_map_max: int = 6,
    ) -> None:
        """
        Initialize and compute the fisheye map, Jacobian, and induced reference metric tensors.

        Fisheye parameters are registered as NRPy CodeParameters so they can be overridden
        at code generation time.

        :param num_transitions: Number of fisheye transitions N. Must be at least 1.
        :param a_defaults: Default plateau stretch factors [a0, ..., aN].
                           If None, all entries are set to 1.0. Length must be
                           num_transitions + 1.
        :param R_defaults: Default transition centers [R1, ..., RN].
                           If None, they are set to [1.0, 2.0, ..., float(N)].
                           Length must be num_transitions.
        :param s_defaults: Default transition widths [s1, ..., sN].
                           If None, all entries are set to 0.5. Length must be
                           num_transitions.
        :param c_default: Default global scaling factor c. Defaults to 1.0.
        :param use_closed_form_derivs: If True, use closed form radial derivatives for
                                       the single transition kernel and the radius map.
                                       If False, compute derivatives using sp.diff in 1D.
        :param use_power_map_subs: If True, substitute powers of r_sym using r and r^2
                                   via an xreplace mapping (codegen shape optimization).
        :param power_map_max: Maximum absolute power used in the power map substitutions.
        :raises ValueError: If num_transitions is less than 1.
        :raises ValueError: If any default list has an inconsistent length.
        """
        n_transitions = int(num_transitions)
        if n_transitions < 1:
            raise ValueError(f"num_transitions must be >= 1; got {num_transitions}.")

        a_default_floats = (
            [1.0] * (n_transitions + 1)
            if a_defaults is None
            else [float(v) for v in a_defaults]
        )
        if len(a_default_floats) != n_transitions + 1:
            raise ValueError(
                f"a_defaults must have length {n_transitions + 1}; got {len(a_default_floats)}."
            )

        R_default_floats = (
            [float(i + 1) for i in range(n_transitions)]
            if R_defaults is None
            else [float(v) for v in R_defaults]
        )
        if len(R_default_floats) != n_transitions:
            raise ValueError(
                f"R_defaults must have length {n_transitions}; got {len(R_default_floats)}."
            )

        s_default_floats = (
            [0.5] * n_transitions
            if s_defaults is None
            else [float(v) for v in s_defaults]
        )
        if len(s_default_floats) != n_transitions:
            raise ValueError(
                f"s_defaults must have length {n_transitions}; got {len(s_default_floats)}."
            )

        self.num_transitions = n_transitions

        self.a_list = par.register_CodeParameters(
            "REAL",
            thismodule,
            [f"fisheye_a{i}" for i in range(n_transitions + 1)],
            _as_par_default_list(a_default_floats),
            commondata=True,
        )
        self.R_list = par.register_CodeParameters(
            "REAL",
            thismodule,
            [f"fisheye_R{i + 1}" for i in range(n_transitions)],
            _as_par_default_list(R_default_floats),
            commondata=True,
        )
        self.s_list = par.register_CodeParameters(
            "REAL",
            thismodule,
            [f"fisheye_s{i + 1}" for i in range(n_transitions)],
            _as_par_default_list(s_default_floats),
            commondata=True,
        )
        self.c = par.register_CodeParameter(
            "REAL",
            thismodule,
            "fisheye_c",
            defaultvalue=float(c_default),
            commondata=True,
        )

        self.xx = list(ixp.declarerank1("xx", dimension=3))

        # r_squared = xx[0]^2 + xx[1]^2 + xx[2]^2
        r_squared = sum(self.xx[i] ** 2 for i in range(3))

        # r = sqrt(r_squared)
        self.r = sp.sqrt(r_squared)

        # Build 1D radial expressions in r_sym, then substitute r_sym -> r(xx)
        r_sym = sp.Symbol("r_fisheye", real=True, positive=True)

        # (rb0_u, rb1_u, rb2_u, rb3_u) = (rbar_unscaled, d/dr, d2/dr2, d3/dr3)
        rb0_u, rb1_u, rb2_u, rb3_u = _unscaled_radius_map_and_derivs(
            r_sym,
            self.a_list,
            self.R_list,
            self.s_list,
            use_closed_form_derivs=use_closed_form_derivs,
        )

        # u = rbar = c * rbar_unscaled
        u = self.c * rb0_u
        # p = du/dr
        p = self.c * rb1_u
        # q = d2u/dr2
        q = self.c * rb2_u
        # t = d3u/dr3
        t = self.c * rb3_u

        # Build radial scalars:
        # lam(r), dlam/dr, A(r), B(r), A'(r), B'(r), A''(r), B''(r)
        lam_r, dlam_dr_r, A_r, B_r, A1_r, B1_r, A2_r, B2_r = _radial_metric_scalars(
            r_sym, u, p, q, t
        )

        # Replace r_sym with the full Cartesian radius r(xx)
        if use_power_map_subs:
            sub = _rsym_substitution_dict(
                r_sym, self.r, r_squared, max_power=int(power_map_max)
            )

            def replace(expr: sp.Expr) -> sp.Expr:
                return expr.xreplace(sub)

        else:

            def replace(expr: sp.Expr) -> sp.Expr:
                return expr.subs({r_sym: self.r})

        # lam = lam(r(xx)), dlam_dr = (dlam/dr)(r(xx)), etc.
        lam = replace(lam_r)
        dlam_dr = replace(dlam_dr_r)
        A = replace(A_r)
        B = replace(B_r)
        A1 = replace(A1_r)
        B1 = replace(B1_r)
        A2 = replace(A2_r)
        B2 = replace(B2_r)

        # rbar_unscaled(xx) = rb0_u with r_sym -> r(xx)
        self.rbar_unscaled = replace(rb0_u)
        # rbar(xx) = u with r_sym -> r(xx)
        self.rbar = replace(u)

        # inv_r = 1/r, inv_r2 = 1/r^2, inv_r3 = 1/r^3
        inv_r = 1 / self.r
        inv_r2 = 1 / r_squared
        inv_r3 = inv_r2 * inv_r

        # A1_over_r = A'(r)/r, B1_over_r = B'(r)/r
        A1_over_r = A1 * inv_r
        B1_over_r = B1 * inv_r

        # coeff_A = A''/r^2 - A'/r^3
        coeff_A = A2 * inv_r2 - A1 * inv_r3
        # coeff_B = B''/r^2 - B'/r^3
        coeff_B = B2 * inv_r2 - B1 * inv_r3

        # dlam_over_r = (dlam/dr)/r
        dlam_over_r = dlam_dr * inv_r

        xx = self.xx
        xx_prod = [[xx[i] * xx[j] for j in range(3)] for i in range(3)]
        KD = _KroneckerDelta3D

        # Forward map: CartU[i] = lam * xx[i]
        self.xx_to_CartU = [lam * xx[i] for i in range(3)]

        # Jacobian:
        # d CartU[mu] / d xx[i] = lam * delta_{mu i} + (dlam/dr)/r * xx[mu] * xx[i]
        self.dCart_dxxUD = [
            [lam * KD[mu][i] + dlam_over_r * xx_prod[mu][i] for i in range(3)]
            for mu in range(3)
        ]

        # Reference metric:
        # ghatDD[i][j] = A * delta_ij + B * xx[i] * xx[j]
        self.ghatDD = ixp.zerorank2(dimension=3)
        for i in range(3):
            for j in range(3):
                self.ghatDD[i][j] = A * KD[i][j] + B * xx_prod[i][j]

        # First derivatives:
        # ghatDDdD[i][j][k] =
        #   (A'(r)/r) * xx[k] * delta_ij
        # + (B'(r)/r) * xx[k] * xx[i] * xx[j]
        # + B * (delta_{ik} * xx[j] + xx[i] * delta_{jk})
        self.ghatDDdD = ixp.zerorank3(dimension=3)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.ghatDDdD[i][j][k] = (
                        A1_over_r * xx[k] * KD[i][j]
                        + B1_over_r * xx[k] * xx_prod[i][j]
                        + B * (KD[i][k] * xx[j] + xx[i] * KD[j][k])
                    )

        # Second derivatives:
        # Implemented in component form using coeff_A, coeff_B and product rule expansions.
        self.ghatDDdDD = ixp.zerorank4(dimension=3)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        term1 = (coeff_A * xx_prod[k][l] + A1_over_r * KD[k][l]) * KD[
                            i
                        ][j]
                        term2 = coeff_B * xx_prod[k][l] * xx_prod[i][j] + B1_over_r * (
                            KD[k][l] * xx_prod[i][j]
                            + xx[k] * (KD[i][l] * xx[j] + xx[i] * KD[j][l])
                        )
                        term3 = (
                            B1_over_r * xx[l] * KD[i][k] * xx[j]
                            + B * KD[i][k] * KD[j][l]
                        )
                        term4 = (
                            B1_over_r * xx[l] * xx[i] * KD[j][k]
                            + B * KD[i][l] * KD[j][k]
                        )
                        self.ghatDDdDD[i][j][k][l] = term1 + term2 + term3 + term4


def make_fisheye_rfm(
    num_transitions: int,
    *,
    a_defaults: Optional[Sequence[float]] = None,
    R_defaults: Optional[Sequence[float]] = None,
    s_defaults: Optional[Sequence[float]] = None,
    c_default: float = 1.0,
    use_closed_form_derivs: bool = False,
    use_power_map_subs: bool = False,
    power_map_max: int = 6,
) -> GeneralRFMFisheye:
    """
    Construct and return a GeneralRFMFisheye instance.

    :param num_transitions: Number of fisheye transitions. Must be at least 1.
    :param a_defaults: Default plateau stretch factors. Length must be num_transitions + 1.
    :param R_defaults: Default transition centers. Length must be num_transitions.
    :param s_defaults: Default transition widths. Length must be num_transitions.
    :param c_default: Default global scaling factor.
    :param use_closed_form_derivs: If True, use closed form radial derivatives.
    :param use_power_map_subs: If True, substitute powers of r_sym using r and r^2.
    :param power_map_max: Maximum power for power map substitutions.
    :return: A newly constructed GeneralRFMFisheye object.
    """
    return GeneralRFMFisheye(
        num_transitions,
        a_defaults=a_defaults,
        R_defaults=R_defaults,
        s_defaults=s_defaults,
        c_default=c_default,
        use_closed_form_derivs=use_closed_form_derivs,
        use_power_map_subs=use_power_map_subs,
        power_map_max=power_map_max,
    )


def _as_par_default_list(values: Sequence[float]) -> list[_CodeParameterDefault]:
    """
    Convert a float sequence into a list suitable for NRPy parameter registration.

    NRPy register_CodeParameters expects default values of type str, int, float, or a list
    of those types. This helper converts a sequence of floats into a list of
    _CodeParameterDefault entries.

    :param values: Sequence of float defaults.
    :return: List of defaults as _CodeParameterDefault.
    """
    defaults: list[_CodeParameterDefault] = []
    for val in values:
        defaults.append(float(val))
    return defaults


def _kernel_G(r: sp.Expr, R: sp.Expr, s: sp.Expr) -> sp.Expr:
    """
    Single transition kernel G(r; R, s) used in the multi transition map.

    The kernel is

        G(r; R, s) =
            s / (2 * tanh(R / s)) *
            log( cosh((r + R)/s) / cosh((r - R)/s) ).

    :param r: Raw radius r (symbol or expression).
    :param R: Transition center R.
    :param s: Transition width s.
    :return: Kernel value G(r; R, s).
    """
    # G(r; R, s) = s / (2*tanh(R/s)) * log( cosh((r+R)/s) / cosh((r-R)/s) )
    tanh_R_over_s = sp.tanh(R / s)
    return (s / (2 * tanh_R_over_s)) * sp.log(
        sp.cosh((r + R) / s) / sp.cosh((r - R) / s)
    )


def _kernel_G_and_derivs(
    r: sp.Expr, R: sp.Expr, s: sp.Expr
) -> Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """
    Return the single transition kernel G(r; R, s) and its first three radial derivatives.

    The kernel is

        G(r; R, s) =
            s / (2 * tanh(R / s)) *
            log( cosh((r + R)/s) / cosh((r - R)/s) ).

    Closed form radial derivatives (no sp.diff calls needed):

        G'(r)  = (tanh((r+R)/s) - tanh((r-R)/s)) / (2*tanh(R/s))

        G''(r) = (sech((r+R)/s)^2 - sech((r-R)/s)^2) / (2*s*tanh(R/s))

        G'''(r)= (sech((r-R)/s)^2*tanh((r-R)/s) - sech((r+R)/s)^2*tanh((r+R)/s))
                 / (s^2*tanh(R/s))

    In the implementation, sech(x)^2 is written as 1/cosh(x)^2 to avoid introducing
    a separate SymPy sech function.

    :param r: Raw radius r (symbol or expression).
    :param R: Transition center R.
    :param s: Transition width s.
    :return: Tuple (G, G1, G2, G3) where Gk = d^k G / d r^k.
    """
    # tanh_R_over_s = tanh(R/s)
    tanh_R_over_s = sp.tanh(R / s)
    # inv_2_tanh_R_over_s = 1 / (2*tanh(R/s))
    inv_2_tanh_R_over_s = sp.Integer(1) / (2 * tanh_R_over_s)

    # arg_plus = (r + R)/s, arg_minus = (r - R)/s
    arg_plus = (r + R) / s
    arg_minus = (r - R) / s

    # cosh_plus = cosh((r + R)/s), cosh_minus = cosh((r - R)/s)
    cosh_plus = sp.cosh(arg_plus)
    cosh_minus = sp.cosh(arg_minus)

    # G = s/(2*tanh(R/s)) * log(cosh_plus/cosh_minus)
    G = (s * inv_2_tanh_R_over_s) * sp.log(cosh_plus / cosh_minus)

    # G1 = (tanh(arg_plus) - tanh(arg_minus)) / (2*tanh(R/s))
    tanh_plus = sp.tanh(arg_plus)
    tanh_minus = sp.tanh(arg_minus)
    G1 = inv_2_tanh_R_over_s * (tanh_plus - tanh_minus)

    # sech2(x) = 1/cosh(x)^2, so:
    # G2 = (sech2(arg_plus) - sech2(arg_minus)) / (2*s*tanh(R/s))
    sech2_plus = sp.Integer(1) / (cosh_plus * cosh_plus)
    sech2_minus = sp.Integer(1) / (cosh_minus * cosh_minus)
    G2 = inv_2_tanh_R_over_s * (sech2_plus - sech2_minus) / s

    # G3 = (sech2(arg_minus)*tanh(arg_minus) - sech2(arg_plus)*tanh(arg_plus)) / (s^2*tanh(R/s))
    inv_tanh_R_over_s = sp.Integer(1) / tanh_R_over_s
    G3 = (
        inv_tanh_R_over_s
        * (sech2_minus * tanh_minus - sech2_plus * tanh_plus)
        / (s * s)
    )

    return G, G1, G2, G3


def _unscaled_radius_map_and_derivs(
    r: sp.Expr,
    a_list: Sequence[sp.Expr],
    R_list: Sequence[sp.Expr],
    s_list: Sequence[sp.Expr],
    *,
    use_closed_form_derivs: bool,
) -> Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """
    Construct the unscaled N transition radius map and its first three radial derivatives.

    The unscaled radius map is

        rbar_unscaled(r) =
            a_N * r + sum_{i=1..N} (a_{i-1} - a_i) * G(r; R_i, s_i)

    where G(r; R, s) is the single transition kernel.

    This helper returns:
    - rb0 = rbar_unscaled(r)
    - rb1 = d rb0 / d r
    - rb2 = d^2 rb0 / d r^2
    - rb3 = d^3 rb0 / d r^3

    There are two computation modes:

    1. If use_closed_form_derivs is True, each transition contributes G, G', G'', G'''
       from _kernel_G_and_derivs, and the map derivatives are accumulated in a single pass.

    2. If use_closed_form_derivs is False, the map rb0 is built from _kernel_G, then
       derivatives are computed using sp.diff in 1D.

    :param r: Radial symbol r used for 1D differentiation.
    :param a_list: Plateau stretch factors [a0, ..., aN] as SymPy expressions.
    :param R_list: Transition centers [R1, ..., RN] as SymPy expressions.
    :param s_list: Transition widths [s1, ..., sN] as SymPy expressions.
    :param use_closed_form_derivs: Select closed form or sp.diff derivatives.
    :return: Tuple (rb0, rb1, rb2, rb3).
    :raises ValueError: If the parameter list lengths are inconsistent.
    """
    n_transitions = len(R_list)
    if len(a_list) != n_transitions + 1 or len(s_list) != n_transitions:
        raise ValueError("Inconsistent fisheye parameter list lengths.")

    a_last = a_list[-1]

    if use_closed_form_derivs:
        # rb0 = a_N * r
        rb0 = a_last * r
        # rb1 = d(rb0)/dr = a_N
        rb1 = a_last
        # rb2 = d2(rb0)/dr2 = 0, rb3 = d3(rb0)/dr3 = 0
        rb2 = sp.Integer(0)
        rb3 = sp.Integer(0)

        for i in range(n_transitions):
            # delta_a = a_{i} - a_{i+1}
            delta_a = a_list[i] - a_list[i + 1]
            # Add delta_a * (G, G1, G2, G3) to (rb0, rb1, rb2, rb3)
            G0, G1, G2, G3 = _kernel_G_and_derivs(r, R_list[i], s_list[i])
            rb0 += delta_a * G0
            rb1 += delta_a * G1
            rb2 += delta_a * G2
            rb3 += delta_a * G3

        return rb0, rb1, rb2, rb3

    # rb0 = a_N * r + sum_i (a_i - a_{i+1}) * G(r; R_i, s_i)
    rb0 = a_last * r
    for i in range(n_transitions):
        delta_a = a_list[i] - a_list[i + 1]
        rb0 += delta_a * _kernel_G(r, R_list[i], s_list[i])

    # rb1 = d(rb0)/dr, rb2 = d(rb1)/dr, rb3 = d(rb2)/dr
    rb1 = sp.diff(rb0, r)
    rb2 = sp.diff(rb1, r)
    rb3 = sp.diff(rb2, r)
    return rb0, rb1, rb2, rb3


def _radial_metric_scalars(
    r: sp.Expr, u: sp.Expr, p: sp.Expr, q: sp.Expr, t: sp.Expr
) -> Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """
    Build radial scalar functions for the map and for the induced reference metric.

    Inputs are radial quantities:
    - u = rbar(r)
    - p = d u / d r
    - q = d^2 u / d r^2
    - t = d^3 u / d r^3

    The map is CartU[i] = lam(r) * xx[i] with

        lam(r) = u / r

    and

        d lam / d r = (p * r - u) / r^2.

    The induced metric is written in the radial-tangential decomposition

        ghatDD[i][j] = A(r) * delta_ij + B(r) * xx[i] * xx[j].

    For a spherically symmetric radius map u = rbar(r):
    - the tangential eigenvalue is A(r) = (u/r)^2 = u^2 / r^2
    - the radial eigenvalue is (d rbar / d r)^2 = p^2
      and since the radial eigenvalue in the A,B decomposition is A + B * r^2, we get

          A + B * r^2 = p^2

      so

          B(r) = p^2 / r^2 - u^2 / r^4.

    Radial derivatives are obtained by differentiating with respect to r and using
    u' = p, p' = q, q' = t. The closed forms used in the code are:

        A1 = dA/dr = 2*u*p/r^2 - 2*u^2/r^3

        B1 = dB/dr = 2*p*q/r^2 - 2*p^2/r^3 - 2*u*p/r^4 + 4*u^2/r^5

        A2 = d^2A/dr^2 = 2*(p^2 + u*q)/r^2 - 8*u*p/r^3 + 6*u^2/r^4

        B2 = d^2B/dr^2 = 2*(q^2 + p*t)/r^2 - 8*p*q/r^3 + (4*p^2 - 2*u*q)/r^4
                         + 16*u*p/r^5 - 20*u^2/r^6

    :param r: Radial coordinate symbol r.
    :param u: Physical radius u = rbar(r).
    :param p: First derivative p = d u / d r.
    :param q: Second derivative q = d^2 u / d r^2.
    :param t: Third derivative t = d^3 u / d r^3.
    :return: Tuple (lam, dlam_dr, A, B, A1, B1, A2, B2).
    """
    # lam = u / r
    lam = u / r
    # dlam_dr = (p*r - u) / r^2
    dlam_dr = (p * r - u) / (r**2)

    # A = u^2 / r^2
    A = u**2 / (r**2)
    # B = p^2 / r^2 - u^2 / r^4
    B = p**2 / (r**2) - u**2 / (r**4)

    # A1 = dA/dr = 2*u*p/r^2 - 2*u^2/r^3
    A1 = 2 * u * p / (r**2) - 2 * u**2 / (r**3)
    # B1 = dB/dr = 2*p*q/r^2 - 2*p^2/r^3 - 2*u*p/r^4 + 4*u^2/r^5
    B1 = 2 * p * q / (r**2) - 2 * p**2 / (r**3) - 2 * u * p / (r**4) + 4 * u**2 / (r**5)

    # A2 = d2A/dr2 = 2*(p^2 + u*q)/r^2 - 8*u*p/r^3 + 6*u^2/r^4
    A2 = 2 * (p**2 + u * q) / (r**2) - 8 * u * p / (r**3) + 6 * u**2 / (r**4)
    # B2 = d2B/dr2 = 2*(q^2 + p*t)/r^2 - 8*p*q/r^3 + (4*p^2 - 2*u*q)/r^4
    #               + 16*u*p/r^5 - 20*u^2/r^6
    B2 = (
        2 * (q**2 + p * t) / (r**2)
        - 8 * p * q / (r**3)
        + (4 * p**2 - 2 * u * q) / (r**4)
        + 16 * u * p / (r**5)
        - 20 * u**2 / (r**6)
    )

    return lam, dlam_dr, A, B, A1, B1, A2, B2


def _rsym_substitution_dict(
    r_sym: sp.Symbol, r: sp.Expr, r2: sp.Expr, *, max_power: int
) -> Dict[sp.Expr, sp.Expr]:
    """
    Build an xreplace substitution dictionary mapping powers of r_sym into expressions of r and r^2.

    This is a codegen-shape optimization: even powers can be expressed without an extra
    square root. The mapping includes at least:
    - r_sym -> r
    - r_sym^2 -> r2
    - r_sym^-1 -> 1/r
    - r_sym^-2 -> 1/r2

    and extends to higher powers up to max_power.

    :param r_sym: The 1D radial symbol used in intermediate expressions.
    :param r: The Cartesian radius expression r(xx).
    :param r2: The squared radius expression r^2(xx).
    :param max_power: Maximum absolute power to include in the mapping.
    :return: Dictionary suitable for expr.xreplace(...).
    """
    # inv_r = 1/r, inv_r2 = 1/r2
    inv_r = 1 / r
    inv_r2 = 1 / r2

    # Base substitutions:
    # r_sym -> r
    # r_sym^2 -> r2
    # r_sym^-1 -> 1/r
    # r_sym^-2 -> 1/r2
    sub: Dict[sp.Expr, sp.Expr] = {
        r_sym: r,
        r_sym**2: r2,
        r_sym**-1: inv_r,
        r_sym**-2: inv_r2,
    }

    # Extend to higher powers using r2 and r to reduce sqrt usage
    for n in range(3, max_power + 1):
        half = n // 2
        if n % 2 == 0:
            # r_sym^n -> (r2)^(n/2)
            sub[r_sym**n] = r2**half
            # r_sym^-n -> (1/r2)^(n/2)
            sub[r_sym**-n] = inv_r2**half
        else:
            # r_sym^n -> (r2)^((n-1)/2) * r
            sub[r_sym**n] = r2**half * r
            # r_sym^-n -> (1/r2)^((n-1)/2) * (1/r)
            sub[r_sym**-n] = inv_r2**half * inv_r

    return sub


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

    for N in (1, 2):
        print(f"Setting up GeneralRFMFisheye[N={N}]...")
        _fisheye = GeneralRFMFisheye(
            num_transitions=N,
            use_closed_form_derivs=True,
            use_power_map_subs=True,
            power_map_max=6,
        )
        results_dict = ve.process_dictionary_of_expressions(
            _fisheye.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_N{N}",
            results_dict,
        )
