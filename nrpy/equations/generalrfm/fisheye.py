# nrpy/equations/generalrfm/fisheye.py
"""
General N transition fisheye raw to physical (xx -> Cart) mapping and reference metric.

This module implements a radially symmetric, arbitrary N transition fisheye mapping
from raw Cartesian coordinates xx^i to physical Cartesian coordinates Cart^i,
together with the associated flat space reference metric in the raw coordinates and
its first and second derivatives.

Compared to the original implementation, this version is optimized for performance:

* It exploits the purely radial form of the map, Cart^i = λ(r) xx^i, where
  λ(r) = rbar(r) / r.
* All "expensive" derivatives are taken only with respect to a 1D radial symbol r,
  never with respect to the Cartesian coordinates xx^i.
* The reference metric and its derivatives are expressed in terms of radial scalar
  functions of r and simple Cartesian tensors (δ_ij and xx^i), avoiding repeated
  large sympy.diff calls with respect to xx^i.
* All radial derivatives of the transition kernel and radius map are built in closed
  form, avoiding SymPy differentiation (sp.diff) entirely.

This module is the *single source of truth* for the **equations** of the fisheye map:
* forward mapping xx -> Cart,
* the radial map r -> rbar(r),
* the inverse-map 1D root equation f(r) = rbar(r) - rCart (used by Newton-Raphson),
* and the reference metric induced by the map.

Infrastructure modules (e.g. BHaH codegen) should consume the SymPy expressions from
this module and handle C codegen + numerical algorithms externally.

Summary of the geometry:

* Raw to physical map: xx^i -> Cart^i(xx).
* Reference metric: ghat_ij = δ_mn (∂ Cart^m / ∂ xx^i) (∂ Cart^n / ∂ xx^j).
* First derivatives: ghat_ij,k computed from analytic radial formulas.
* Second derivatives: ghat_ij,kl computed from analytic radial formulas.

The raw radius is

    r = sqrt( sum_i (xx^i)^2 ).

The N transition fisheye map is defined by

* Plateau stretch (zoom out) factors a_0, ..., a_N.
* Transition centers R_1, ..., R_N.
* Width parameters s_1, ..., s_N.
* Differences Δ a_i = a_{i-1} - a_i.

The single transition kernel is

    G(r; R, s) =
        s / (2 * tanh(R / s)) *
        log( cosh( (r + R) / s ) / cosh( (r - R) / s ) ).

The unscaled radius map is

    rbar_unscaled(r) =
        a_N * r + sum_{i=1}^N (a_{i-1} - a_i) * G(r; R_i, s_i).

A global scale factor c produces the final radius map

    rbar(r) = c * rbar_unscaled(r).

The physical Cartesian coordinates are obtained by a purely radial rescaling,

    Cart^i = (rbar(r) / r) * xx^i,

which leaves the angular coordinates unchanged.

Inverse-map equation (Cart -> xx):

Let rCart = ||Cart||. Taking norms of Cart^i = (rbar(r)/r) xx^i gives

    rCart = rbar(r).

So the inverse map reduces to solving the 1D root equation

    f(r) = rbar(r) - rCart = 0,   with f'(r) = d rbar / dr,

and then recovering

    xx^i = (r / rCart) Cart^i    (for rCart > 0), with xx^i = 0 at rCart = 0.

From this mapping, the flat reference metric in raw coordinates can be written in the
standard "radial-tangential" decomposition

    ghat_ij(r, xx) = A(r) δ_ij + B(r) xx^i xx^j,

where A(r) encodes the tangential (angular) scaling and A(r) + B(r) r^2 encodes the
radial scaling. For a spherically symmetric map of the physical radius rbar(r) we have

    A(r)     = (rbar(r) / r)^2,
    A(r) + B(r) r^2 = (d rbar / d r)^2,

so that

    B(r) = ( (d rbar / d r)^2 - A(r) ) / r^2.

The first and second derivatives ghat_ij,k and ghat_ij,kl are then built from
A(r), B(r) and their radial derivatives using the chain rule and the simple
Cartesian identities

    ∂_k r = xx^k / r,
    ∂_k xx^i = δ^i_k,
    ∂_l ∂_k xx^i = 0,

without ever asking SymPy to differentiate huge expressions with respect to xx^i.

Fisheye parameters are stored as NRPy code parameters so that they can be
configured at code generation time.
"""

from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, List, Optional, Tuple, Union, cast

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.params as par

CodeParameterDefaultList = List[Union[str, int, float]]

# Kronecker delta in 3D as a plain Python constant for fast tensor assembly.
_KroneckerDelta3D: Tuple[Tuple[int, int, int], Tuple[int, int, int], Tuple[int, int, int]] = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
)


@dataclass(frozen=True)
class FisheyeRadialExprs:
    """
    Container for 1D radial SymPy expressions for the N-transition fisheye map.

    These expressions are intended for infrastructure code that needs:
      * rbar(r) and drbar/dr for forward mapping codegen, and/or
      * f(r) = rbar(r) - rCart and f'(r) for inverse-map Newton-Raphson.

    All expressions are built in terms of a single radial symbol r_sym (1D), not
    the full Cartesian radius sqrt(xx^i xx^i).

    :ivar num_transitions: Number of fisheye transitions N (>= 1).
    :ivar r_sym: 1D raw-radius symbol used to build all radial expressions.
    :ivar rCart_sym: 1D physical-radius symbol rCart = ||Cart|| used in f(r).
    :ivar rbar_unscaled_of_r: Unscaled radius map rbar_unscaled(r_sym).
    :ivar rbar_of_r: Scaled radius map rbar(r_sym) = c * rbar_unscaled(r_sym).
    :ivar drbar_dr_of_r: First derivative d rbar / d r.
    :ivar lam_of_r: Radial scaling λ(r) = rbar(r)/r.
    :ivar dlam_dr_of_r: dλ/dr = (r rbar'(r) - rbar(r))/r^2.
    :ivar lam_at_origin: Limiting λ(0) = c * a0 (useful for safe r->0 handling).
    :ivar f_of_r: Root function f(r) = rbar(r) - rCart.
    :ivar fprime_of_r: Root derivative f'(r) = d rbar / d r.
    """

    num_transitions: int
    r_sym: sp.Symbol
    rCart_sym: sp.Symbol

    rbar_unscaled_of_r: sp.Expr
    rbar_of_r: sp.Expr
    drbar_dr_of_r: sp.Expr

    lam_of_r: sp.Expr
    dlam_dr_of_r: sp.Expr
    lam_at_origin: sp.Expr

    f_of_r: sp.Expr
    fprime_of_r: sp.Expr


# -----------------------------------------------------------------------------
# Internal helpers: parameter registration and purely-radial kernel/map builders.
# -----------------------------------------------------------------------------
@dataclass(frozen=True)
class _FisheyeParams:
    """
    Container for registered fisheye CodeParameters.

    This exists so that radial-only expression builders (e.g. fisheye_radial_exprs)
    can register and access the canonical fisheye parameters without instantiating
    the full GeneralRFMFisheye object (which also builds Jacobians/metrics/derivs).
    """

    num_transitions: int
    a_list: List[sp.Expr]
    R_list: List[sp.Expr]
    s_list: List[sp.Expr]
    c: sp.Expr


def _register_fisheye_params(
    num_transitions: int,
    a_default: Optional[List[float]] = None,
    R_default: Optional[List[float]] = None,
    s_default: Optional[List[float]] = None,
    c_default: float = 1.0,
) -> _FisheyeParams:
    """
    Register the fisheye CodeParameters and return the associated SymPy symbols.

    The registration logic is shared between:
      * GeneralRFMFisheye (full metric + derivatives), and
      * fisheye_radial_exprs() (radial-only expressions for forward/inverse maps).

    :param num_transitions: Number of fisheye transitions N. Must be at least 1.
    :param a_default: Default plateau stretch factors [a0, ..., aN].
                      If None, all entries are set to 1.0. Length must be num_transitions + 1.
    :param R_default: Default raw transition centers [R1, ..., RN].
                      If None, they are set to [1.0, 2.0, ..., float(N)].
                      Length must be num_transitions.
    :param s_default: Default raw transition widths [s1, ..., sN].
                      If None, all entries are set to 0.5. Length must be num_transitions.
    :param c_default: Default global scaling factor c. Defaults to 1.0.
    :return: _FisheyeParams containing the registered CodeParameter symbols.
    :raises ValueError: If num_transitions is less than 1.
    :raises ValueError: If any of the default lists have inconsistent lengths.
    """
    if num_transitions < 1:
        raise ValueError(
            f"num_transitions must be >= 1; got num_transitions = {num_transitions}."
        )

    # ---------------------------------------------------------------------
    # Step 1: Set up default parameter values and register NRPy CodeParameters.
    # ---------------------------------------------------------------------
    if a_default is None:
        # Default: all plateaus have unit stretch before global scaling.
        a_default_list: List[float] = [1.0 for _ in range(num_transitions + 1)]
    else:
        a_default_list = [float(val) for val in a_default]
    if len(a_default_list) != num_transitions + 1:
        raise ValueError(
            "a_default must have length num_transitions + 1; "
            f"got len(a_default) = {len(a_default_list)} while "
            f"num_transitions + 1 = {num_transitions + 1}."
        )

    if R_default is None:
        # Simple monotonically increasing centers as a fallback.
        R_default_list: List[float] = [float(i + 1) for i in range(num_transitions)]
    else:
        R_default_list = [float(val) for val in R_default]
    if len(R_default_list) != num_transitions:
        raise ValueError(
            "R_default must have length num_transitions; "
            f"got len(R_default) = {len(R_default_list)} while "
            f"num_transitions = {num_transitions}."
        )

    if s_default is None:
        # Default: moderate transition widths.
        s_default_list: List[float] = [0.5 for _ in range(num_transitions)]
    else:
        s_default_list = [float(val) for val in s_default]
    if len(s_default_list) != num_transitions:
        raise ValueError(
            "s_default must have length num_transitions; "
            f"got len(s_default) = {len(s_default_list)} while "
            f"num_transitions = {num_transitions}."
        )

    # Plateau stretch factors a_0..a_N:
    a_names = [f"fisheye_a{i}" for i in range(num_transitions + 1)]
    a_list = list(
        par.register_CodeParameters(
            "REAL",
            __name__,
            a_names,
            cast(CodeParameterDefaultList, a_default_list),
            commondata=True,
        )
    )

    # Raw transition centers R_1..R_N:
    R_names = [f"fisheye_R{i + 1}" for i in range(num_transitions)]
    R_list = list(
        par.register_CodeParameters(
            "REAL",
            __name__,
            R_names,
            cast(CodeParameterDefaultList, R_default_list),
            commondata=True,
        )
    )

    # Raw transition widths s_1..s_N:
    s_names = [f"fisheye_s{i + 1}" for i in range(num_transitions)]
    s_list = list(
        par.register_CodeParameters(
            "REAL",
            __name__,
            s_names,
            cast(CodeParameterDefaultList, s_default_list),
            commondata=True,
        )
    )

    # Global scale factor c:
    c = par.register_CodeParameter(
        "REAL",
        __name__,
        "fisheye_c",
        c_default,
        commondata=True,
    )

    return _FisheyeParams(
        num_transitions=num_transitions,
        a_list=a_list,
        R_list=R_list,
        s_list=s_list,
        c=cast(sp.Expr, c),
    )


def _fisheye_G_kernel_and_derivs(
    r: sp.Expr, R: sp.Expr, s: sp.Expr, max_deriv: int
) -> Tuple[sp.Expr, ...]:
    """
    Return the single-transition kernel G(r; R, s) and (optionally) its radial derivatives.

    The kernel is

        G(r; R, s) =
            s / (2 * tanh(R / s)) *
            log( cosh( (r + R) / s ) / cosh( (r - R) / s ) ).

    Closed-form radial derivatives (no SymPy sp.diff calls needed):

        G'(r)  = [tanh((r+R)/s) - tanh((r-R)/s)] / [2*tanh(R/s)],

        G''(r) = [sech^2((r+R)/s) - sech^2((r-R)/s)] / [2*s*tanh(R/s)],

        G'''(r)= [sech^2((r-R)/s)*tanh((r-R)/s) - sech^2((r+R)/s)*tanh((r+R)/s)]
                 / [s^2*tanh(R/s)].

    :param r: Raw radius r (may be a symbol or an expression).
    :param R: Raw transition center R.
    :param s: Raw transition width s, assumed to be strictly positive.
    :param max_deriv: Maximum derivative order to return (0..3).
    :return: Tuple (G,) or (G, G1) or (G, G1, G2) or (G, G1, G2, G3).
    :raises ValueError: If max_deriv is outside 0..3.
    """
    if max_deriv < 0 or max_deriv > 3:
        raise ValueError(f"max_deriv must be in [0,3]; got max_deriv = {max_deriv}.")

    tanh_R_over_s = sp.tanh(R / s)
    inv_2_tanh_R_over_s = sp.Integer(1) / (2 * tanh_R_over_s)

    arg_plus = (r + R) / s
    arg_minus = (r - R) / s

    cosh_plus = sp.cosh(arg_plus)
    cosh_minus = sp.cosh(arg_minus)

    G = cast(
        sp.Expr,
        (s * inv_2_tanh_R_over_s) * sp.log(cosh_plus / cosh_minus),
    )
    if max_deriv == 0:
        return (G,)

    tanh_plus = sp.tanh(arg_plus)
    tanh_minus = sp.tanh(arg_minus)
    G1 = cast(sp.Expr, inv_2_tanh_R_over_s * (tanh_plus - tanh_minus))
    if max_deriv == 1:
        return (G, G1)

    # sech^2(x) = 1/cosh(x)^2. Using cosh(x) avoids introducing sech() as a separate
    # SymPy function (useful for some codegen backends).
    sech2_plus = cast(sp.Expr, sp.Integer(1) / (cosh_plus * cosh_plus))
    sech2_minus = cast(sp.Expr, sp.Integer(1) / (cosh_minus * cosh_minus))
    G2 = cast(sp.Expr, inv_2_tanh_R_over_s * (sech2_plus - sech2_minus) / s)
    if max_deriv == 2:
        return (G, G1, G2)

    inv_tanh_R_over_s = sp.Integer(1) / tanh_R_over_s
    G3 = cast(
        sp.Expr,
        inv_tanh_R_over_s
        * (sech2_minus * tanh_minus - sech2_plus * tanh_plus)
        / (s * s),
    )
    return (G, G1, G2, G3)


def _build_unscaled_radius_map_and_derivs(
    r: sp.Expr,
    num_transitions: int,
    a_list: List[sp.Expr],
    R_list: List[sp.Expr],
    s_list: List[sp.Expr],
    max_deriv: int,
) -> Tuple[sp.Expr, ...]:
    """
    Construct the unscaled N transition radius map rbar_unscaled(r) and (optionally) its first three radial derivatives.

    The map is

        rbar_unscaled(r) =
            a_N * r + sum_{i=1}^N (a_{i-1} - a_i) * G(r; R_i, s_i),

    where G(r; R, s) is the single transition kernel.

    :param r: Raw radius r (may be a symbol or an expression).
    :param num_transitions: Number of fisheye transitions N (>= 1).
    :param a_list: Plateau stretch factors [a0, ..., aN].
    :param R_list: Transition centers [R1, ..., RN].
    :param s_list: Transition widths [s1, ..., sN].
    :param max_deriv: Maximum derivative order to return (0..3).
    :return: Tuple (rb0,) or (rb0, rb1) or (rb0, rb1, rb2) or (rb0, rb1, rb2, rb3),
             where rbk = d^k rbar_unscaled / d r^k.
    :raises ValueError: If max_deriv is outside 0..3.
    """
    if max_deriv < 0 or max_deriv > 3:
        raise ValueError(f"max_deriv must be in [0,3]; got max_deriv = {max_deriv}.")

    a_N = a_list[-1]
    rbar0 = cast(sp.Expr, a_N * r)

    # Fast path: only build rbar itself.
    if max_deriv == 0:
        for i in range(num_transitions):
            delta_a_i = a_list[i] - a_list[i + 1]
            R_i = R_list[i]
            s_i = s_list[i]
            G = _fisheye_G_kernel_and_derivs(r, R_i, s_i, max_deriv=0)[0]
            rbar0 += cast(sp.Expr, delta_a_i * G)
        return (cast(sp.Expr, rbar0),)

    rbar1 = cast(sp.Expr, a_N)
    rbar2 = cast(sp.Expr, sp.Integer(0)) if max_deriv >= 2 else None
    rbar3 = cast(sp.Expr, sp.Integer(0)) if max_deriv >= 3 else None

    for i in range(num_transitions):
        delta_a_i = a_list[i] - a_list[i + 1]
        R_i = R_list[i]
        s_i = s_list[i]

        G_terms = _fisheye_G_kernel_and_derivs(r, R_i, s_i, max_deriv=max_deriv)

        rbar0 += cast(sp.Expr, delta_a_i * G_terms[0])
        rbar1 += cast(sp.Expr, delta_a_i * G_terms[1])

        if max_deriv >= 2:
            assert rbar2 is not None
            rbar2 += cast(sp.Expr, delta_a_i * G_terms[2])
        if max_deriv >= 3:
            assert rbar3 is not None
            rbar3 += cast(sp.Expr, delta_a_i * G_terms[3])

    out: List[sp.Expr] = [cast(sp.Expr, rbar0), cast(sp.Expr, rbar1)]
    if max_deriv >= 2:
        assert rbar2 is not None
        out.append(cast(sp.Expr, rbar2))
    if max_deriv >= 3:
        assert rbar3 is not None
        out.append(cast(sp.Expr, rbar3))
    return tuple(out)


def _make_rsym_to_r_xreplace_dict(
    r_sym: sp.Symbol, r: sp.Expr, r2: sp.Expr, max_power: int = 6
) -> Dict[sp.Expr, sp.Expr]:
    """
    Build an xreplace dictionary mapping powers of r_sym into expressions built from r = sqrt(r2) and r2 = r^2.

    This is purely a performance and codegen-shape optimization: even powers of r can
    be expressed without an extra sqrt, e.g. r^2 -> r2, r^4 -> r2^2, etc.

    :param r_sym: The 1D radial symbol used in intermediate expressions.
    :param r: The Cartesian radius expression sqrt(r2).
    :param r2: The squared Cartesian radius expression (sum_i xx_i^2).
    :param max_power: Maximum |power| to include in the mapping.
    :return: Dictionary suitable for SymPy expr.xreplace().
    """
    inv_r = cast(sp.Expr, 1 / r)
    inv_r2 = cast(sp.Expr, 1 / r2)

    sub: Dict[sp.Expr, sp.Expr] = {
        r_sym: r,
        r_sym**2: r2,
        r_sym**-1: inv_r,
        r_sym**-2: inv_r2,
    }

    for n in range(3, max_power + 1):
        if n % 2 == 0:
            sub[r_sym**n] = cast(sp.Expr, r2 ** (n // 2))
            sub[r_sym**-n] = cast(sp.Expr, inv_r2 ** (n // 2))
        else:
            sub[r_sym**n] = cast(sp.Expr, (r2 ** ((n - 1) // 2)) * r)
            sub[r_sym**-n] = cast(sp.Expr, (inv_r2 ** ((n - 1) // 2)) * inv_r)

    return sub


@lru_cache(maxsize=None)
def fisheye_radial_exprs(num_transitions: int) -> FisheyeRadialExprs:
    """
    Build and return the canonical 1D radial SymPy expressions for an N-transition fisheye map.

    This function is intentionally cached so that repeated codegen paths (e.g. multiple
    infrastructure functions needing rbar and/or f(r)) do not repeatedly rebuild the
    same symbolic expressions.

    NOTE: Calling this function registers the required fisheye CodeParameters:
      * fisheye_a0..fisheye_aN
      * fisheye_R1..fisheye_RN
      * fisheye_s1..fisheye_sN
      * fisheye_c

    :param num_transitions: Number of fisheye transitions N (>= 1).
    :return: A FisheyeRadialExprs container holding rbar(r), derivatives, and f(r).
    :raises ValueError: If num_transitions is less than 1.
    """
    if num_transitions < 1:
        raise ValueError(
            f"num_transitions must be >= 1 for fisheye; got {num_transitions}."
        )

    # Register the required CodeParameters without instantiating the full
    # GeneralRFMFisheye object (which would also build Jacobians/metrics/derivs).
    params = _register_fisheye_params(num_transitions=num_transitions)

    r_sym = sp.Symbol("r_fisheye", real=True, positive=True)
    rCart_sym = sp.Symbol("rCart", real=True, nonnegative=True)

    # Build rbar_unscaled(r_sym) and its first radial derivative in closed form.
    rbar_unscaled_of_r, drbar_unscaled_dr_of_r = _build_unscaled_radius_map_and_derivs(
        r=r_sym,
        num_transitions=num_transitions,
        a_list=params.a_list,
        R_list=params.R_list,
        s_list=params.s_list,
        max_deriv=1,
    )

    # Scaled physical radius and radial derivative.
    rbar_of_r = cast(sp.Expr, params.c * rbar_unscaled_of_r)
    drbar_dr_of_r = cast(sp.Expr, params.c * drbar_unscaled_dr_of_r)

    lam_of_r = cast(sp.Expr, rbar_of_r / r_sym)
    dlam_dr_of_r = cast(sp.Expr, (drbar_dr_of_r * r_sym - rbar_of_r) / (r_sym**2))

    # Safe r->0 limit: rbar(r) ~ c*a0*r => lam(0) = c*a0.
    lam_at_origin = cast(sp.Expr, params.c * params.a_list[0])

    # Inverse-map 1D root equation: rCart = rbar(r).
    f_of_r = cast(sp.Expr, rbar_of_r - rCart_sym)
    fprime_of_r = cast(sp.Expr, drbar_dr_of_r)

    return FisheyeRadialExprs(
        num_transitions=num_transitions,
        r_sym=r_sym,
        rCart_sym=rCart_sym,
        rbar_unscaled_of_r=rbar_unscaled_of_r,
        rbar_of_r=rbar_of_r,
        drbar_dr_of_r=drbar_dr_of_r,
        lam_of_r=lam_of_r,
        dlam_dr_of_r=dlam_dr_of_r,
        lam_at_origin=lam_at_origin,
        f_of_r=f_of_r,
        fprime_of_r=fprime_of_r,
    )


class GeneralRFMFisheye:
    """
    Construct an N transition fisheye raw to physical map and reference metric.

    The mapping is radially symmetric in the raw Cartesian coordinates xx^i and
    is defined by a multi-transition kernel. Given plateau stretch factors a_i,
    transition centers R_i, widths s_i, and a global scale c, the raw radius

        r = sqrt( sum_i (xx^i)^2 )

    is mapped to a physical radius rbar(r). The physical Cartesian coordinates are

        Cart^i = (rbar(r) / r) * xx^i,

    for r > 0, with the expression understood in the limiting sense at r = 0.

    From this transformation, the flat space reference metric ghat_ij in the raw
    coordinates admits the decomposition

        ghat_ij = A(r) δ_ij + B(r) xx^i xx^j,

    where A(r) = (rbar / r)^2 encodes the tangential scaling and
    A(r) + B(r) r^2 = (d rbar / d r)^2 encodes the radial scaling. The functions
    A(r), B(r) and their radial derivatives are computed symbolically in 1D using
    a dedicated radial symbol r_fisheye, then converted into full expressions in
    (xx^0, xx^1, xx^2) by the substitution r_fisheye -> r.

    This approach is significantly faster than directly differentiating the full
    metric with respect to xx^i using SymPy, especially for the second derivatives
    ghat_ij,kl, which would otherwise require many large second-order Cartesian
    derivatives of complicated expressions.

    :ivar num_transitions: Number of fisheye transitions N.
    :ivar a_list: Plateau stretch factors a_0, ..., a_N, stored as NRPy code
        parameters.
    :ivar R_list: Raw transition centers R_1, ..., R_N, stored as NRPy code
        parameters.
    :ivar s_list: Raw transition widths s_1, ..., s_N, stored as NRPy code
        parameters.
    :ivar c: Global scaling factor c, stored as an NRPy code parameter.
    :ivar xx: Raw Cartesian coordinate symbols [xx0, xx1, xx2].
    :ivar r: Raw radius r = sqrt( sum_i xx_i^2 ).
    :ivar rbar_unscaled: Unscaled physical radius rbar_unscaled(r(xx)).
    :ivar rbar: Scaled physical radius rbar(r(xx)) = c * rbar_unscaled(r(xx)).
    :ivar xx_to_CartU: Raw to physical Cartesian map Cart^i(xx^j).
    :ivar dCart_dxxUD: Jacobian ∂ Cart^i / ∂ xx^j, built from radial formulas.
    :ivar ghatDD: Reference metric ghat_ij in raw coordinates.
    :ivar ghatDDdD: First derivatives ∂_k ghat_ij, built from radial formulas.
    :ivar ghatDDdDD: Second derivatives ∂_l ∂_k ghat_ij, built from radial formulas.
    """

    def __init__(
        self,
        num_transitions: int,
        a_default: Optional[List[float]] = None,
        R_default: Optional[List[float]] = None,
        s_default: Optional[List[float]] = None,
        c_default: float = 1.0,
    ) -> None:
        """
        Initialize the N transition fisheye mapping and reference metric in 3D.

        Fisheye parameters are registered as NRPy CodeParameters so that they
        can be overridden at code generation time.

        :param num_transitions: Number of fisheye transitions N. Must be at least 1.
        :param a_default: Default plateau stretch factors [a0, ..., aN].
                          If None, all entries are set to 1.0. Length must be
                          num_transitions + 1.
        :param R_default: Default raw transition centers [R1, ..., RN].
                          If None, they are set to [1.0, 2.0, ..., float(N)].
                          Length must be num_transitions.
        :param s_default: Default raw transition widths [s1, ..., sN].
                          If None, all entries are set to 0.5. Length must be
                          num_transitions.
        :param c_default: Default global scaling factor c. Defaults to 1.0.

        :raises ValueError: If num_transitions is less than 1.
        :raises ValueError: If any of the default lists have inconsistent lengths.
        """
        if num_transitions < 1:
            raise ValueError(
                f"num_transitions must be >= 1; got num_transitions = {num_transitions}."
            )

        self.num_transitions = num_transitions

        # ---------------------------------------------------------------------
        # Step 1: Set up default parameter values and register NRPy CodeParameters.
        # ---------------------------------------------------------------------
        params = _register_fisheye_params(
            num_transitions=num_transitions,
            a_default=a_default,
            R_default=R_default,
            s_default=s_default,
            c_default=c_default,
        )
        self.a_list = params.a_list
        self.R_list = params.R_list
        self.s_list = params.s_list
        self.c = params.c

        # ---------------------------------------------------------------------
        # Step 2: Define raw coordinates, radial symbol, and radial map.
        # ---------------------------------------------------------------------
        # Raw Cartesian coordinates xx^i (3D).
        self.xx = list(ixp.declarerank1("xx", dimension=3))

        # Raw radius r = sqrt(xx^2 + yy^2 + zz^2).
        # This is the actual radius expression in terms of xx^i.
        r2 = cast(sp.Expr, sum(self.xx[i] ** 2 for i in range(3)))
        self.r = cast(sp.Expr, sp.sqrt(r2))

        # For performance, we build the radial map and its derivatives in terms of
        # a *symbolic* radial variable r_sym, and only at the end substitute
        # r_sym -> r(xx). This keeps all expensive symbolic work strictly 1D in r_sym.
        #
        # In addition, all radial derivatives of the kernel and map are built in
        # closed form, so no SymPy differentiation is performed.
        r_sym = sp.Symbol("r_fisheye", real=True, positive=True)

        # Unscaled radius map and its first three derivatives (closed form).
        (
            rbar_unscaled_of_r,
            drbar_unscaled_dr_of_r,
            d2rbar_unscaled_dr2_of_r,
            d3rbar_unscaled_dr3_of_r,
        ) = _build_unscaled_radius_map_and_derivs(
            r=r_sym,
            num_transitions=self.num_transitions,
            a_list=self.a_list,
            R_list=self.R_list,
            s_list=self.s_list,
            max_deriv=3,
        )

        # Scaled physical radius and its radial derivatives:
        rbar_of_r = cast(sp.Expr, self.c * rbar_unscaled_of_r)
        drbar_dr = cast(sp.Expr, self.c * drbar_unscaled_dr_of_r)
        d2rbar_dr2 = cast(sp.Expr, self.c * d2rbar_unscaled_dr2_of_r)
        d3rbar_dr3 = cast(sp.Expr, self.c * d3rbar_unscaled_dr3_of_r)

        # Substitution dictionary for converting radial expressions to Cartesian.
        # This maps even powers of r_sym into r2 = r^2 to reduce sqrt usage in downstream codegen.
        sub = _make_rsym_to_r_xreplace_dict(r_sym=r_sym, r=self.r, r2=r2, max_power=6)

        # Now build the actual coordinate-dependent radii by substituting r_sym -> r(xx).
        self.rbar_unscaled = cast(sp.Expr, rbar_unscaled_of_r.xreplace(sub))
        self.rbar = cast(sp.Expr, rbar_of_r.xreplace(sub))

        # ---------------------------------------------------------------------
        # Radial scaling factor λ(r) = rbar(r) / r and its radial derivative.
        #     λ(r)     = rbar(r) / r,
        #     λ'(r)    = (d rbar / d r * r - rbar) / r^2.
        # ---------------------------------------------------------------------
        lam_of_r = cast(sp.Expr, rbar_of_r / r_sym)
        dlam_dr_of_r = cast(sp.Expr, (drbar_dr * r_sym - rbar_of_r) / (r_sym**2))

        # ---------------------------------------------------------------------
        # Radial metric functions A(r) and B(r) and their derivatives.
        #
        # In spherical coordinates built from the raw radius r:
        #     g_rr = (d rbar / d r)^2,
        #     g_T  = (rbar / r)^2,
        #
        # where g_T is the tangential (angular) part. In the Cartesian decomposition
        #     g_ij = A(r) δ_ij + B(r) xx^i xx^j,
        # the eigenvalues are
        #     g_T     = A(r),
        #     g_rr    = A(r) + B(r) r^2.
        # Therefore:
        #     A(r) = (rbar / r)^2,
        #     B(r) = (g_rr - A(r)) / r^2
        #          = ( (d rbar / d r)^2 - A(r) ) / r^2.
        #
        # We compute A, B and their derivatives using closed-form formulas
        # in terms of u=rbar, p=u', q=u'', t=u''' to avoid sp.diff entirely.
        # ---------------------------------------------------------------------
        u, p, q, t = rbar_of_r, drbar_dr, d2rbar_dr2, d3rbar_dr3
        r = r_sym  # shorthand for formulas below

        # A(r) = u² / r²
        A_of_r = cast(sp.Expr, u**2 / r**2)

        # B(r) = p²/r² - u²/r⁴
        B_of_r = cast(sp.Expr, p**2 / r**2 - u**2 / r**4)

        # A'(r) = 2up/r² - 2u²/r³
        A1_of_r = cast(sp.Expr, 2 * u * p / r**2 - 2 * u**2 / r**3)

        # B'(r) = 2pq/r² - 2p²/r³ - 2up/r⁴ + 4u²/r⁵
        B1_of_r = cast(
            sp.Expr,
            2 * p * q / r**2 - 2 * p**2 / r**3 - 2 * u * p / r**4 + 4 * u**2 / r**5,
        )

        # A''(r) = 2(p² + uq)/r² - 8up/r³ + 6u²/r⁴
        A2_of_r = cast(
            sp.Expr,
            2 * (p**2 + u * q) / r**2 - 8 * u * p / r**3 + 6 * u**2 / r**4,
        )

        # B''(r) = 2(q² + pt)/r² - 8pq/r³ + (4p² - 2uq)/r⁴ + 16up/r⁵ - 20u²/r⁶
        B2_of_r = cast(
            sp.Expr,
            2 * (q**2 + p * t) / r**2
            - 8 * p * q / r**3
            + (4 * p**2 - 2 * u * q) / r**4
            + 16 * u * p / r**5
            - 20 * u**2 / r**6,
        )

        # Convert all radial functions to full expressions in the Cartesian coordinates
        # by substituting r_sym -> r(xx). Use xreplace for efficiency.
        lam, dlam_dr, A, B, A1, B1, A2, B2 = [
            cast(sp.Expr, expr.xreplace(sub))
            for expr in (
                lam_of_r,
                dlam_dr_of_r,
                A_of_r,
                B_of_r,
                A1_of_r,
                B1_of_r,
                A2_of_r,
                B2_of_r,
            )
        ]

        # ---------------------------------------------------------------------
        # Pre-compute common factors for the tensor loops.
        # ---------------------------------------------------------------------
        inv_r = cast(sp.Expr, 1 / self.r)
        inv_r2 = cast(sp.Expr, 1 / r2)
        inv_r3 = cast(sp.Expr, inv_r2 * inv_r)

        A1_over_r = cast(sp.Expr, A1 * inv_r)
        B1_over_r = cast(sp.Expr, B1 * inv_r)

        coeff_A = cast(sp.Expr, A2 * inv_r2 - A1 * inv_r3)
        coeff_B = cast(sp.Expr, B2 * inv_r2 - B1 * inv_r3)

        dlam_over_r = cast(sp.Expr, dlam_dr * inv_r)

        # Pre-compute xx products.
        xx = self.xx
        xx_prod = [[xx[i] * xx[j] for j in range(3)] for i in range(3)]

        # ---------------------------------------------------------------------
        # Step 3: Raw (xx) -> physical Cartesian (Cart) map.
        #         Cart^i = λ(r(xx)) * xx^i
        #
        # We still store Cart^i explicitly for completeness, but note that the
        # Jacobian is built below from the radial identities rather than by
        # differentiating these expressions with respect to xx^i.
        # ---------------------------------------------------------------------
        self.xx_to_CartU = [lam * xx[i] for i in range(3)]

        # ---------------------------------------------------------------------
        # Step 4: Jacobian d Cart^i / d xx^j.
        #
        # For a purely radial map Cart^i = λ(r) xx^i, with r = sqrt(xx^k xx^k),
        # we have
        #
        #     ∂_j r = xx^j / r,
        #     ∂_j λ = λ'(r) ∂_j r = λ'(r) xx^j / r,
        #
        # so
        #
        #     ∂_j Cart^i = ∂_j (λ xx^i)
        #                 = λ δ^i_j + xx^i ∂_j λ
        #                 = λ δ^i_j + (λ'(r) / r) xx^i xx^j.
        #
        # This avoids any SymPy differentiation with respect to xx^i and is
        # much faster than calling sp.diff on Cart^i(xx^j).
        # ---------------------------------------------------------------------
        KD = _KroneckerDelta3D
        self.dCart_dxxUD = [
            [lam * KD[mu][i] + dlam_over_r * xx_prod[mu][i] for i in range(3)]
            for mu in range(3)
        ]

        # ---------------------------------------------------------------------
        # Step 5: Reference metric ghat_ij = A(r) δ_ij + B(r) xx^i xx^j.
        # ---------------------------------------------------------------------
        self.ghatDD = ixp.zerorank2(dimension=3)
        for i in range(3):
            for j in range(3):
                self.ghatDD[i][j] = A * KD[i][j] + B * xx_prod[i][j]

        # ---------------------------------------------------------------------
        # Step 6: First derivatives ghat_ij,k via radial formulas.
        # ---------------------------------------------------------------------
        self.ghatDDdD = ixp.zerorank3(dimension=3)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.ghatDDdD[i][j][k] = (
                        A1_over_r * xx[k] * KD[i][j]
                        + B1_over_r * xx[k] * xx_prod[i][j]
                        + B * (KD[i][k] * xx[j] + xx[i] * KD[j][k])
                    )

        # ---------------------------------------------------------------------
        # Step 7: Second derivatives ghat_ij,kl via radial formulas.
        # ---------------------------------------------------------------------
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

    # -------------------------------------------------------------------------
    # Internal helper methods.
    # -------------------------------------------------------------------------

    def _G_kernel(self, r: sp.Expr, R: sp.Expr, s: sp.Expr) -> sp.Expr:
        """
        Single transition kernel G(r; R, s) used in the multi transition map.

        The kernel is

            G(r; R, s) =
                s / (2 * tanh(R / s)) *
                log( cosh( (r + R) / s ) / cosh( (r - R) / s ) ).

        :param r: Raw radius r (may be a symbol or an expression).
        :param R: Raw transition center R.
        :param s: Raw transition width s, assumed to be strictly positive.
        :return: The kernel value G(r; R, s).
        """
        return cast(
            sp.Expr, _fisheye_G_kernel_and_derivs(r=r, R=R, s=s, max_deriv=0)[0]
        )

    def _build_unscaled_radius_map(self, r: sp.Expr) -> sp.Expr:
        """
        Construct the unscaled N transition radius map rbar_unscaled(r).

        The map is

            rbar_unscaled(r) =
                a_N * r + sum_{i=1}^N (a_{i-1} - a_i) * G(r; R_i, s_i),

        where G(r; R, s) is the single transition kernel defined in _G_kernel.

        :param r: Raw radius r (may be a symbol or an expression).
        :return: The unscaled physical radius rbar_unscaled(r).
        """
        return cast(
            sp.Expr,
            _build_unscaled_radius_map_and_derivs(
                r=r,
                num_transitions=self.num_transitions,
                a_list=self.a_list,
                R_list=self.R_list,
                s_list=self.s_list,
                max_deriv=0,
            )[0],
        )


def build_generalrfm_fisheye(
    n_fisheye_transitions: int,
    a_default: Optional[List[float]] = None,
    R_default: Optional[List[float]] = None,
    s_default: Optional[List[float]] = None,
    c_default: float = 1.0,
) -> GeneralRFMFisheye:
    """
    Construct a GeneralRFMFisheye instance.

    This function instantiates GeneralRFMFisheye using the input
    n_fisheye_transitions value and passes through the default plateau
    stretches, transition centers, widths, and global scale, which in turn
    define the default values of the underlying NRPy CodeParameters.

    :param n_fisheye_transitions: Number of fisheye transitions N. Must be at least 1.
    :param a_default: Default plateau stretch factors [a0, ..., aN].
                      If None, all entries are set to 1.0. Length must be
                      n_fisheye_transitions + 1.
    :param R_default: Default raw transition centers [R1, ..., RN].
                      If None, they are set to [1.0, 2.0, ..., float(N)].
                      Length must be n_fisheye_transitions.
    :param s_default: Default raw transition widths [s1, ..., sN].
                      If None, all entries are set to 0.5. Length must be
                      n_fisheye_transitions.
    :param c_default: Default global scaling factor c. Defaults to 1.0.
    :return: A newly constructed GeneralRFMFisheye instance.

    :raises ValueError: If n_fisheye_transitions is less than 1.
    :raises ValueError: If any of the default lists have inconsistent lengths.
    """
    if n_fisheye_transitions < 1:
        raise ValueError(
            "n_fisheye_transitions must be >= 1; "
            f"got n_fisheye_transitions = {n_fisheye_transitions}."
        )

    return GeneralRFMFisheye(
        num_transitions=n_fisheye_transitions,
        a_default=a_default,
        R_default=R_default,
        s_default=s_default,
        c_default=c_default,
    )


if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.validate_expressions.validate_expressions as ve

    # Run doctests for this module (none are defined explicitly, but this keeps
    # behavior consistent with other NRPy modules).
    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Basic validation for a small set of configurations (3D only).
    # Note: thanks to the radial formulas used for ghat_ij, ghat_ij,k, and
    # ghat_ij,kl, this validation is significantly faster than an approach
    # based on repeatedly calling sp.diff with respect to xx^i.
    for N in (1, 2):
        print(f"Setting up GeneralRFMFisheye[N={N}]...")
        _fisheye = GeneralRFMFisheye(num_transitions=N)
        results_dict = ve.process_dictionary_of_expressions(
            _fisheye.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_N{N}",
            results_dict,
        )
