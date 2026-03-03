# nrpy/equations/generalrfm/fisheye.py
"""
General N transition fisheye raw to physical (xx to Cart) mapping and reference metric.

This module implements a radially symmetric, arbitrary N transition fisheye mapping
from raw Cartesian coordinates xx[i] to physical Cartesian coordinates Cart[i],
together with the induced flat reference metric in the raw coordinates and its
first and second derivatives.

Key equations used by this module:

Coordinate map and radii
Define the raw radius:
    r = sqrt(xx[0]**2 + xx[1]**2 + xx[2]**2)

The map is a purely radial rescaling:
    Cart[i](xx) = lam(r) * xx[i]
    lam(r)      = rbar(r) / r

N transition fisheye radius map

Parameters
- plateau stretch factors a0, a1, ..., aN
- transition centers R1, R2, ..., RN
- transition widths s1, s2, ..., sN
- global scaling factor c

Single transition kernel:

    G(r; R, s) = s / (2*tanh(R/s)) * log( cosh((r+R)/s) / cosh((r-R)/s) )

Unscaled radius map:

    rbar_unscaled(r) = aN*r + sum_{i=1..N} (a_{i-1} - a_i) * G(r; R_i, s_i)

Scaled physical radius:

    rbar(r) = c * rbar_unscaled(r)

Induced reference metric

Physical space is flat in Cart coordinates. The induced metric in raw coordinates is:

    ghat_ij = delta_mn * (dCart[m]/dxx[i]) * (dCart[n]/dxx[j])

Jacobian for a radial map

From Cart[i] = lam(r) * xx[i] and dr/dxx[j] = xx[j]/r:

    dCart[i]/dxx[j] = lam * delta^i_j + (dlam/dr)/r * xx[i] * xx[j]

Radial-tangential decomposition of ghat_ij

For any radial map, the induced metric can be written:

    ghat_ij = A(r) * delta_ij + B(r) * xx[i] * xx[j]

Tangential scaling:

    A(r) = lam(r)**2 = (rbar(r)/r)**2

Radial scaling is given by the radial eigenvalue A + B*r**2. For a radial map r -> rbar(r):

    A(r) + B(r) * r**2 = (drbar/dr)**2

Therefore:

    B(r) = ( (drbar/dr)**2 - A(r) ) / r**2

First derivatives of ghat_ij

Using dr/dxx[k] = xx[k]/r and dxx[i]/dxx[k] = delta^i_k:

    d_k ghat_ij =
        (A'(r)/r) * xx[k] * delta_ij
      + (B'(r)/r) * xx[k] * xx[i] * xx[j]
      + B(r) * (delta^i_k * xx[j] + xx[i] * delta^j_k)

Second derivatives of ghat_ij

Differentiate again and use:

    d_l (xx[k]/r) = delta^k_l / r - xx[k]*xx[l] / r**3

This produces coefficient combinations:

    (A''/r**2 - A'/r**3) * xx[k]*xx[l] + (A'/r) * delta_kl

(and similarly for B), which appear in the implementation through precomputed
scalar factors.
"""

from typing import Dict, List, Tuple, cast

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.params as par

thismodule = __name__

__all__ = ["GeneralRFMFisheye", "build_fisheye"]

_KRONECKER_DELTA_3D: Tuple[
    Tuple[int, int, int], Tuple[int, int, int], Tuple[int, int, int]
] = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
)


class GeneralRFMFisheye:
    """
    Construct and store an N transition fisheye map and induced reference metric.

    This class registers fisheye parameters as NRPy code parameters and builds:
    - the forward map Cart[i](xx)
    - the Jacobian dCart[i]/dxx[j]
    - the induced reference metric ghat_ij and its first and second derivatives

    The implementation uses the radial scalar functions A(r) and B(r) from the
    decomposition ghat_ij = A(r) delta_ij + B(r) xx[i] xx[j], so that derivatives
    can be computed using simple Cartesian identities instead of repeated SymPy
    differentiation with respect to xx[i].
    """

    def __init__(
        self,
        num_transitions: int,
    ) -> None:
        """
        Initialize the fisheye map and induced reference metric expressions.

        :param num_transitions: Number of fisheye transitions. Must be at least 1.
        :raises ValueError: If num_transitions is less than 1 or default list lengths are inconsistent.

        Derivation summary (high level)

        The map is Cart[i] = lam(r) * xx[i] with lam(r) = rbar(r)/r. The Jacobian is:

            dCart[i]/dxx[j] = lam * delta^i_j + (dlam/dr)/r * xx[i] * xx[j]

        The induced metric can be written:

            ghat_ij = A(r) * delta_ij + B(r) * xx[i] * xx[j]

        with:

            A(r) = lam(r)**2

        and from the radial eigenvalue:

            A(r) + B(r) * r**2 = (drbar/dr)**2

        so:

            B(r) = ( (drbar/dr)**2 - A(r) ) / r**2

        First and second derivatives follow by differentiating A(r) and B(r) with
        respect to r and applying:

            dr/dxx[k] = xx[k]/r
            d_l (xx[k]/r) = delta^k_l / r - xx[k]*xx[l] / r**3
        """
        if num_transitions < 1:
            raise ValueError(f"num_transitions must be >= 1; got {num_transitions}.")

        self.num_transitions = num_transitions

        # Step 1: Register fisheye CodeParameters
        # - a0..aN, R1..RN, s1..sN, c
        self.a_list = par.register_CodeParameters(
            "REAL",
            thismodule,
            [f"fisheye_a{i}" for i in range(num_transitions + 1)],
            [1.0] * (num_transitions + 1),
            commondata=False,
        )
        self.R_list = par.register_CodeParameters(
            "REAL",
            thismodule,
            [f"fisheye_R{i + 1}" for i in range(num_transitions)],
            [float(i + 1) for i in range(num_transitions)],
            commondata=False,
        )
        self.s_list = par.register_CodeParameters(
            "REAL",
            thismodule,
            [f"fisheye_s{i + 1}" for i in range(num_transitions)],
            [0.5] * num_transitions,
            commondata=False,
        )
        self.c = par.register_CodeParameter(
            "REAL", thismodule, "fisheye_c", defaultvalue=1.0, commondata=False
        )

        # Step 2: Define raw Cartesian coordinates and radius
        # xx[i] are raw Cartesian coordinates
        self.xx = list(ixp.declarerank1("xx", dimension=3))

        # r2 = xx[0]**2 + xx[1]**2 + xx[2]**2
        r2 = sum(self.xx[i] ** 2 for i in range(3))

        # r = sqrt(r2)
        self.r = sp.sqrt(r2)

        # Step 3: Build the radial map in 1D using a dedicated radial symbol
        # r_sym is used only for 1D radial expressions
        r_sym = sp.Symbol("r_fisheye", real=True, positive=True)

        # rbar_unscaled(r) = aN*r + sum_i (a[i] - a[i+1]) * G(r; R_i, s_i)
        # and derivatives through third order: d/dr, d2/dr2, d3/dr3
        (
            rbar_unscaled_sym,
            drbar_unscaled_sym,
            d2rbar_unscaled_sym,
            d3rbar_unscaled_sym,
        ) = _radius_map_unscaled_and_derivs_closed_form(
            r=r_sym, a_list=self.a_list, R_list=self.R_list, s_list=self.s_list
        )

        # rbar(r) = c * rbar_unscaled(r)
        u = self.c * rbar_unscaled_sym
        # p = drbar/dr, q = d2rbar/dr2, t = d3rbar/dr3
        p = self.c * drbar_unscaled_sym
        q = self.c * d2rbar_unscaled_sym
        t = self.c * d3rbar_unscaled_sym

        # lam(r) = rbar(r) / r
        lam_sym = u / r_sym
        # dlam/dr = (p*r - u) / r**2
        dlam_sym = (p * r_sym - u) / (r_sym**2)

        # Metric scalars from the radial map:
        # A(r)  = (rbar/r)**2
        # B(r)  = ( (drbar/dr)**2 - A(r) ) / r**2 = p**2/r**2 - u**2/r**4
        # plus A1, B1, A2, B2
        A_sym, B_sym, A1_sym, B1_sym, A2_sym, B2_sym = _A_B_and_derivs_from_rbar(
            r=r_sym, u=u, p=p, q=q, t=t
        )

        # Step 4: Convert 1D radial expressions into full Cartesian expressions
        # Optionally reduce explicit sqrt usage by rewriting powers of r_sym using r and r2.
        sub = _rpow_xreplace_dict(
            r_sym=r_sym,
            r=self.r,
            r2=cast(sp.Expr, r2),
        )
        self.rbar_unscaled = rbar_unscaled_sym.xreplace(sub)
        self.rbar = u.xreplace(sub)
        lam, dlam_dr, A, B, A1, B1, A2, B2 = [
            expr.xreplace(sub)
            for expr in (
                lam_sym,
                dlam_sym,
                A_sym,
                B_sym,
                A1_sym,
                B1_sym,
                A2_sym,
                B2_sym,
            )
        ]

        # Step 5: Precompute common radial factors used in tensor assembly
        # inv_r  = 1/r
        # inv_r2 = 1/r**2  (implemented as 1/r2)
        # inv_r3 = 1/r**3
        inv_r = sp.Integer(1) / self.r
        inv_r2 = sp.Integer(1) / r2
        inv_r3 = inv_r2 * inv_r

        # A1_over_r = A'(r)/r
        # B1_over_r = B'(r)/r
        A1_over_r = A1 * inv_r
        B1_over_r = B1 * inv_r

        # coeff_A = A''/r**2 - A'/r**3
        # coeff_B = B''/r**2 - B'/r**3
        coeff_A = A2 * inv_r2 - A1 * inv_r3
        coeff_B = B2 * inv_r2 - B1 * inv_r3

        # dlam_over_r = (dlam/dr) / r
        dlam_over_r = dlam_dr * inv_r

        # Step 6: Build map, Jacobian, metric, and derivatives

        # Precompute xx products: xx_prod[i][j] = xx[i] * xx[j]
        xx = self.xx
        xx_prod = [[xx[i] * xx[j] for j in range(3)] for i in range(3)]
        delta = _KRONECKER_DELTA_3D

        # Forward map:
        # Cart[i] = lam(r(xx)) * xx[i]
        self.xx_to_CartU = [lam * xx[i] for i in range(3)]

        # Jacobian:
        # dCart[mu]/dxx[j] = lam * delta[mu][j] + (dlam/dr)/r * xx[mu] * xx[j]
        self.dCart_dxxUD = [
            [lam * delta[mu][j] + dlam_over_r * xx_prod[mu][j] for j in range(3)]
            for mu in range(3)
        ]

        # Inverse Jacobian:
        # dxx[i]/dCart[j] = (1/lam) * delta[i][j]
        #                   - (dlam_over_r / (lam * (lam + dlam_over_r * r^2)))
        #                     * xx[i] * xx[j]
        #
        # For radial maps, lam + dlam_over_r * r^2 = drbar/dr.
        p_local = p.xreplace(sub)
        inv_lam = sp.Integer(1) / lam
        inv_corr = -dlam_over_r / (lam * p_local)
        self.dxx_dCartUD = [
            [inv_lam * delta[i][j] + inv_corr * xx_prod[i][j] for j in range(3)]
            for i in range(3)
        ]

        # Metric decomposition:
        # ghat_ij = A(r) * delta_ij + B(r) * xx[i] * xx[j]
        self.ghatDD = ixp.zerorank2(dimension=3)
        for i in range(3):
            for j in range(3):
                self.ghatDD[i][j] = A * delta[i][j] + B * xx_prod[i][j]

        # First derivatives:
        # d_k ghat_ij =
        #   (A'(r)/r) * xx[k] * delta_ij
        # + (B'(r)/r) * xx[k] * xx[i] * xx[j]
        # + B(r) * (delta[i][k] * xx[j] + xx[i] * delta[j][k])
        self.ghatDDdD = ixp.zerorank3(dimension=3)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.ghatDDdD[i][j][k] = (
                        A1_over_r * xx[k] * delta[i][j]
                        + B1_over_r * xx[k] * xx_prod[i][j]
                        + B * (delta[i][k] * xx[j] + xx[i] * delta[j][k])
                    )

        # Second derivatives:
        # This is the result of differentiating the first-derivative expression and using:
        # d_l (xx[k]/r) = delta[k][l]/r - xx[k]*xx[l]/r**3
        # The implementation uses coeff_A and coeff_B, where:
        # coeff_A = A''/r**2 - A'/r**3
        # coeff_B = B''/r**2 - B'/r**3
        self.ghatDDdDD = ixp.zerorank4(dimension=3)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        # termA corresponds to the part from differentiating the A term
                        termA = (coeff_A * xx_prod[k][l] + A1_over_r * delta[k][l]) * (
                            delta[i][j]
                        )

                        # termB corresponds to the part from differentiating the B prefactor
                        termB = coeff_B * xx_prod[k][l] * xx_prod[i][j] + B1_over_r * (
                            delta[k][l] * xx_prod[i][j]
                            + xx[k] * (delta[i][l] * xx[j] + xx[i] * delta[j][l])
                        )

                        # termC corresponds to differentiating the explicit xx factors in B * (delta*x + x*delta)
                        termC = B1_over_r * xx[l] * (
                            delta[i][k] * xx[j] + xx[i] * delta[j][k]
                        ) + B * (delta[i][k] * delta[j][l] + delta[i][l] * delta[j][k])

                        self.ghatDDdDD[i][j][k][l] = termA + termB + termC

    def radius_map_unscaled_and_derivs_closed_form(
        self, r: sp.Expr
    ) -> Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
        """
        Return the unscaled fisheye radial map and first 3 derivatives.

        :param r: Radial coordinate.
        :return: (rbar_unscaled, drbar_unscaled, d2rbar_unscaled, d3rbar_unscaled)
        """
        return _radius_map_unscaled_and_derivs_closed_form(
            r=r,
            a_list=self.a_list,
            R_list=self.R_list,
            s_list=self.s_list,
        )


def build_fisheye(num_transitions: int) -> GeneralRFMFisheye:
    """
    Construct a GeneralRFMFisheye instance.

    :param num_transitions: Number of fisheye transitions. Must be at least 1.
    :return: A newly constructed GeneralRFMFisheye instance.
    """
    return GeneralRFMFisheye(num_transitions=num_transitions)


def _G_kernel(r: sp.Expr, R: sp.Expr, s: sp.Expr) -> sp.Expr:
    """
    Compute the single transition kernel G(r; R, s).

    Kernel definition:

        G(r; R, s) = s / (2*tanh(R/s)) * log( cosh((r+R)/s) / cosh((r-R)/s) )

    :param r: Raw radius expression.
    :param R: Transition center parameter.
    :param s: Transition width parameter.
    :return: The kernel value G(r; R, s).
    """
    # G(r; R, s) = s/(2*tanh(R/s)) * log( cosh((r+R)/s) / cosh((r-R)/s) )
    return cast(
        sp.Expr,
        (s / (2 * sp.tanh(R / s)))
        * sp.log(sp.cosh((r + R) / s) / sp.cosh((r - R) / s)),
    )


def _radius_map_unscaled(
    r: sp.Expr, a_list: List[sp.Expr], R_list: List[sp.Expr], s_list: List[sp.Expr]
) -> sp.Expr:
    """
    Compute the unscaled N transition radius map rbar_unscaled(r).

    With delta_a_i = a_{i-1} - a_i, the unscaled map is:

        rbar_unscaled(r) = aN*r + sum_{i=1..N} delta_a_i * G(r; R_i, s_i)

    :param r: Raw radius expression.
    :param a_list: Plateau stretch factors a0..aN.
    :param R_list: Transition centers R1..RN.
    :param s_list: Transition widths s1..sN.
    :return: The unscaled radius map rbar_unscaled(r).
    """
    # rb = aN * r
    rb = a_list[-1] * r

    # rb += sum_i (a[i] - a[i+1]) * G(r; R_i, s_i)
    for i, (R_i, s_i) in enumerate(zip(R_list, s_list)):
        delta_a = a_list[i] - a_list[i + 1]
        rb += delta_a * _G_kernel(r=r, R=R_i, s=s_i)

    return cast(sp.Expr, rb)


def _G_and_derivs_closed_form(
    r: sp.Expr, R: sp.Expr, s: sp.Expr
) -> Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """
    Compute G(r; R, s) and its first three r derivatives in closed form.

    Kernel:

        G(r; R, s) = s / (2*tanh(R/s)) * log( cosh((r+R)/s) / cosh((r-R)/s) )

    Closed-form derivatives:

        G1 = (tanh((r+R)/s) - tanh((r-R)/s)) / (2*tanh(R/s))

        G2 = (sech((r+R)/s)**2 - sech((r-R)/s)**2) / (2*s*tanh(R/s))

        G3 = (sech((r-R)/s)**2 * tanh((r-R)/s)
              - sech((r+R)/s)**2 * tanh((r+R)/s)) / (s**2*tanh(R/s))

    Implementation note: sech(x)**2 is represented as 1/cosh(x)**2.

    :param r: Raw radius expression.
    :param R: Transition center parameter.
    :param s: Transition width parameter.
    :return: Tuple (G0, G1, G2, G3) where Gk is the kth derivative with respect to r.
    """
    th = sp.tanh(R / s)
    inv2th = sp.Integer(1) / (2 * th)

    # ap = (r + R)/s, am = (r - R)/s
    ap = (r + R) / s
    am = (r - R) / s

    cp = sp.cosh(ap)
    cm = sp.cosh(am)
    tp = sp.tanh(ap)
    tm = sp.tanh(am)

    # G0 = s/(2*tanh(R/s)) * log( cosh(ap) / cosh(am) )
    G0 = (s * inv2th) * sp.log(cp / cm)

    # G1 = (tanh(ap) - tanh(am)) / (2*tanh(R/s))
    G1 = inv2th * (tp - tm)

    # sech2p = 1/cosh(ap)**2, sech2m = 1/cosh(am)**2
    sech2p = sp.Integer(1) / (cp * cp)
    sech2m = sp.Integer(1) / (cm * cm)

    # G2 = (sech2p - sech2m) / (2*s*tanh(R/s))
    G2 = inv2th * (sech2p - sech2m) / s

    # G3 = (sech2m*tm - sech2p*tp) / (s**2 * tanh(R/s))
    G3 = (sech2m * tm - sech2p * tp) / (s * s * th)

    return G0, G1, G2, G3


def _radius_map_unscaled_and_derivs_closed_form(
    r: sp.Expr, a_list: List[sp.Expr], R_list: List[sp.Expr], s_list: List[sp.Expr]
) -> Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """
    Compute rbar_unscaled(r) and its first three r derivatives in closed form.

    The unscaled map is:

        rbar_unscaled(r) = aN*r + sum_{i=1..N} (a_{i-1} - a_i) * G(r; R_i, s_i)

    This function builds rbar_unscaled and its derivatives through third order
    by summing the corresponding kernel derivatives.

    :param r: Raw radius expression.
    :param a_list: Plateau stretch factors a0..aN.
    :param R_list: Transition centers R1..RN.
    :param s_list: Transition widths s1..sN.
    :return: Tuple (rb0, rb1, rb2, rb3) where rbk is the kth derivative with respect to r.
    """
    # rb0 = rbar_unscaled(r)
    # rb1 = d(rbar_unscaled)/dr
    # rb2 = d2(rbar_unscaled)/dr2
    # rb3 = d3(rbar_unscaled)/dr3
    aN = a_list[-1]
    rb0 = aN * r
    rb1 = aN
    rb2 = sp.Integer(0)
    rb3 = sp.Integer(0)

    for i, (R_i, s_i) in enumerate(zip(R_list, s_list)):
        # delta_a = a[i] - a[i+1]
        delta_a = a_list[i] - a_list[i + 1]

        # Add delta_a * G and its derivatives
        G0, G1, G2, G3 = _G_and_derivs_closed_form(r=r, R=R_i, s=s_i)
        rb0 += delta_a * G0
        rb1 += delta_a * G1
        rb2 += delta_a * G2
        rb3 += delta_a * G3

    return rb0, rb1, rb2, rb3


def _A_B_and_derivs_from_rbar(
    r: sp.Expr, u: sp.Expr, p: sp.Expr, q: sp.Expr, t: sp.Expr
) -> Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """
    Compute A, B and first and second r derivatives from rbar and its derivatives.

    Definitions:
    - u = rbar(r)
    - p = du/dr
    - q = d2u/dr2
    - t = d3u/dr3

    Metric decomposition:
    - A(r) = (u/r)**2 = u**2 / r**2
    - A(r) + B(r)*r**2 = p**2
      so B(r) = (p**2 - A(r)) / r**2 = p**2/r**2 - u**2/r**4

    Differentiating A and B with respect to r and simplifying yields the closed-form
    expressions implemented for A1, B1, A2, and B2.

    :param r: Raw radius symbol.
    :param u: rbar(r).
    :param p: First derivative drbar/dr.
    :param q: Second derivative d2rbar/dr2.
    :param t: Third derivative d3rbar/dr3.
    :return: Tuple (A, B, A1, B1, A2, B2).
    """
    # A(r) = u**2 / r**2
    A = u**2 / r**2

    # B(r) = p**2/r**2 - u**2/r**4
    B = p**2 / r**2 - u**2 / r**4

    # A1(r) = dA/dr = 2*u*p/r**2 - 2*u**2/r**3
    A1 = 2 * u * p / r**2 - 2 * u**2 / r**3

    # B1(r) = dB/dr = 2*p*q/r**2 - 2*p**2/r**3 - 2*u*p/r**4 + 4*u**2/r**5
    B1 = 2 * p * q / r**2 - 2 * p**2 / r**3 - 2 * u * p / r**4 + 4 * u**2 / r**5

    # A2(r) = d2A/dr2 = 2*(p**2 + u*q)/r**2 - 8*u*p/r**3 + 6*u**2/r**4
    A2 = 2 * (p**2 + u * q) / r**2 - 8 * u * p / r**3 + 6 * u**2 / r**4

    # B2(r) =
    #   2*(q**2 + p*t)/r**2 - 8*p*q/r**3 + (4*p**2 - 2*u*q)/r**4 + 16*u*p/r**5 - 20*u**2/r**6
    B2 = (
        2 * (q**2 + p * t) / r**2
        - 8 * p * q / r**3
        + (4 * p**2 - 2 * u * q) / r**4
        + 16 * u * p / r**5
        - 20 * u**2 / r**6
    )
    return A, B, A1, B1, A2, B2


def _rpow_xreplace_dict(
    r_sym: sp.Symbol, r: sp.Expr, r2: sp.Expr
) -> Dict[sp.Expr, sp.Expr]:
    """
    Build an xreplace dictionary mapping powers of r_sym into expressions using r and r2.

    Motivation

    After building 1D radial expressions in terms of r_sym, the final expressions are
    converted into full Cartesian form by replacing r_sym with r = sqrt(r2), where:

        r2 = xx[0]**2 + xx[1]**2 + xx[2]**2

    A direct substitution can introduce many explicit sqrt calls. Even powers of r can
    be represented without sqrt:

    - r**2 -> r2
    - r**4 -> r2**2
    - r**6 -> r2**3

    and similarly for negative powers using 1/r2.

    :param r_sym: 1D radial symbol appearing in intermediate expressions.
    :param r: Cartesian radius expression sqrt(r2).
    :param r2: Cartesian squared radius expression.
    :return: Dictionary suitable for SymPy expr.xreplace().
    """
    inv_r = sp.Integer(1) / r
    inv_r2 = sp.Integer(1) / r2

    # Base substitutions:
    # r_sym -> r
    # r_sym**2 -> r2
    # r_sym**-1 -> 1/r
    # r_sym**-2 -> 1/r2
    sub: Dict[sp.Expr, sp.Expr] = {
        r_sym: r,
        r_sym**2: r2,
        r_sym**-1: inv_r,
        r_sym**-2: inv_r2,
    }

    # Higher powers:
    # even n: r_sym**n -> r2**(n/2)
    # odd  n: r_sym**n -> r2**((n-1)/2) * r
    for n in range(3, 7):  # max power replacement = 6
        if n % 2 == 0:
            sub[r_sym**n] = r2 ** (n // 2)
            sub[r_sym**-n] = inv_r2 ** (n // 2)
        else:
            sub[r_sym**n] = (r2 ** ((n - 1) // 2)) * r
            sub[r_sym**-n] = (inv_r2 ** ((n - 1) // 2)) * inv_r

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

    for N in (1, 2):
        fisheye = GeneralRFMFisheye(num_transitions=N)
        results_dict = ve.process_dictionary_of_expressions(
            fisheye.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_N{N}",
            results_dict,
        )
