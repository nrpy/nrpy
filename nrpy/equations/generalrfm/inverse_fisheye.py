"""
Inverse fisheye specified in physical Cartesian coordinates.

Equal polar and equatorial centers give spherical transitions; distinct centers
select spheroidal transitions within the same map.

Unlike the forward fisheye provider, this construction specifies computational
coordinates explicitly in terms of physical coordinates. It is not generally
the inverse of the forward fisheye map with unchanged parameter values.

Key equations used by this module and the generated BHaH geometry functions:

Coordinates and parameters

    x[i] = Cart[i]             physical Cartesian coordinates
    q[i] = xx[i] = xhat[i]     computational Cartesian coordinates
    r    = sqrt(x[0]**2 + x[1]**2 + x[2]**2)

For layers j = 1, ..., N:

    a0 = 1 <= a1 <= ... <= aN
    R_j  = R_{j,phys}          positive physical equatorial transition center
    Rz_j = R_{z,j,phys}        positive physical polar transition center
    w_j  = w_{j,phys}          positive physical equatorial transition width
    delta_j = 1/a_{j-1} - 1/a_j
    c > 0                     internally computed normalization

The width w_j is the distance between the locations where the centered tanh
argument is -1 and +1. The notation s_j = w_j/2 below is only an algebraic
shorthand; no physical-to-internal transition-parameter conversion is needed.

Spherical inverse map

    drhat/dr = c * [1/aN + sum_{j=1..N} delta_j
                   * (tanh(2*(r+R_j)/w_j) - tanh(2*(r-R_j)/w_j))
                   / (2*tanh(2*R_j/w_j))]

Define the single-transition integral kernel:

    G(rho; R, w) = w/(4*tanh(2*R/w))
                  * log(cosh(2*(rho+R)/w) / cosh(2*(rho-R)/w))

Integrating with rhat(0) = 0 gives:

    rhat(r) = c * [r/aN + sum_{j=1..N} delta_j * G(r; R_j, w_j)]
    q[i]    = (rhat(r)/r) * x[i]

Away from transitions, drhat/dr is approximately c/a_j, hence a fixed
computational radial spacing gives physical spacing:

    Delta r approximately (a_j/c) * Delta rhat

Spheroidal extension

Each layer has its own ellipsoidal radius:

    rho_j**2 = Q_j = x[0]**2 + x[1]**2 + (R_j/Rz_j)**2 * x[2]**2

Replace r by rho_j separately in each G(r)/r term, leaving 1/aN unchanged:

    S(x)   = 1/aN + sum_{j=1..N} delta_j * G(rho_j; R_j, w_j)/rho_j
    lam(x) = c * S(x)
    q[i]   = lam(x) * x[i]

Here S is scale_unscaled, lam is scale, and q is Cart_to_xxU. The ratio
rhat/r = lam depends on direction as well as physical radius. The surface
rho_j = R_j has equatorial semiaxis R_j and polar semiaxis Rz_j. Its polar
width is w_z,j = w_j*Rz_j/R_j, not an independent parameter. Rz_j < R_j is
oblate; Rz_j = R_j for every layer recovers the spherical inverse map.

At the origin, the removable limits telescope:

    lim_{rho->0} G(rho)/rho = 1
    S(0) = 1/aN + sum_j delta_j = 1/a0 = 1
    lam(0) = c,    q(0) = 0

Physical-to-computational conversion is explicit. For the reverse direction,
let qnorm = sqrt(sum_i q[i]**2) and n[i] = q[i]/qnorm. Solve along that ray:

    r * lam(r*n) = qnorm,    x[i] = r*n[i]

For admissible parameters this scalar map is strictly increasing. Since
c/aN <= lam <= c, its root satisfies:

    qnorm/c <= r <= aN*qnorm/c

The zero-radius case is handled directly. InverseFisheyeComponent represents
the numerical position and Jacobian evaluations used by generated BHaH C.

Normalization and grid extents

For a prescribed physical axis extent L > 0, evaluate the map with c = 1:

    Lperp_unscaled = L * S(L, 0, 0)
    Lz_unscaled   = L * S(0, 0, L)
    c             = L/Lperp_unscaled = 1/S(L, 0, 0)
    Lhat_z        = c*Lz_unscaled = L*S(0, 0, L)/S(L, 0, 0)

normalization(L) returns (c, Lhat_z). The computational domain is:

    q[0], q[1] in [-L, L],    q[2] in [-Lhat_z, Lhat_z]

The boundary intersects each physical Cartesian axis at +/-L; its entire
physical shape is not prescribed. For a cell-centered grid with Nx = Ny:

    Delta q[0] = Delta q[1] = 2*L/Nx
    Delta q[2]             = 2*Lhat_z/Nz

Equal computational and central physical spacings require Nz = Nx*Lhat_z/L.
Round Nz to an admissible integer while retaining Lhat_z and the actual
spacing. For grids including both endpoints, replace Nx and Nz by Nx-1 and
Nz-1 in the spacing formulas. The spherical limit has Lhat_z = L.

Stable kernel value and radial derivatives

For a single layer, write s = w/2, t = tanh(R/s), and:

    tp = tanh((rho+R)/s),    tm = tanh((rho-R)/s)
    G1 = dG/drho       = (tp-tm)/(2*t)
    S(u) = sech(u)**2 = 4*exp(-2*abs(u))/(1+exp(-2*abs(u)))**2
    G2 = d2G/drho2     = (S((rho+R)/s)-S((rho-R)/s))/(2*s*t)
    G3 = d3G/drho3     = (S((rho-R)/s)*tm-S((rho+R)/s)*tp)/(s**2*t)

The imported _G_kernel evaluates G without overflowing cosh or subtracting
nearly equal log-cosh values. For rho >= 0 and R, s > 0, define:

    m = min(rho, R)/s,    M = max(rho, R)/s
    E = exp(-2*(M-m))
    D = 2*m + log1p(E*expm1(-4*m)/(1+E))
    G = s*D/(2*t)

D equals log(cosh((rho+R)/s)/cosh((rho-R)/s)). Derivatives use the smooth
closed forms above rather than differentiating the Min/Max representation.

kernel_coefficients returns f0, f1, f2, f3, where Q = rho**2 and:

    f0 = G/rho
    fk = d^k f0/dQ^k,    d/dQ = (1/(2*rho))*d/drho
    f1 = (rho*G1-G)/(2*rho**3)
    f2 = (rho**2*G2-3*rho*G1+3*G)/(4*rho**5)
    f3 = (rho**3*G3-6*rho**2*G2+15*rho*G1-15*G)/(8*rho**7)

For rho below min(2*s, max(7*s/20, R/4)), use the even Taylor polynomial
through rho**24 and its Q derivatives. With y(u)=exp(-2*u), the tanh
derivatives have the form:

    d^n tanh(u)/du^n = y*A_n(y)/(1+y)**(n+1)
    A_1(y) = 4
    A_{n+1}(y) = -2*((1-n*y)*A_n(y)+y*(1+y)*dA_n(y)/dy)

For even n, A_n(y) contains (1-y). Dividing it out cancels the tanh(R/s)
denominator analytically, so no small sech**2 factor is formed by subtraction.

Between the Taylor cutoff and R, compute G-rho and G1-1 before forming
f1, f2, and f3. With z=exp(-2*R/s), a=exp(-2*(R+rho)/s), and
b=exp(-2*(R-rho)/s):

    G-rho = (2*rho*z+s*(1+z)*(log1p(a)-log1p(b))/2)/(1-z)
    G1-1 = (2*z-(1+z)*(a/(1+a)+b/(1+b)))/(1-z)

These differences retain derivatives when G rounds to rho and G1 rounds to 1.

Cartesian scale derivatives consumed by BHaH

For each layer define b_j = (1, 1, (R_j/Rz_j)**2), h_j,i = b_j,i*x[i],
and f_{j,k} = fk evaluated at rho_j. Kronecker deltas are denoted delta_ik;
layer weights delta_j remain the reciprocal-stretch differences above.
There is no sum in the definition of h_j,i. Then:

    v_i   = d lam/dx[i]
          = c * sum_j delta_j * 2*f_{j,1}*h_j,i
    W_ik  = d2 lam/(dx[i] dx[k])
          = c * sum_j delta_j * (4*f_{j,2}*h_j,i*h_j,k
                                + 2*f_{j,1}*b_j,i*delta_ik)
    U_ikl = d3 lam/(dx[i] dx[k] dx[l])
          = c * sum_j delta_j * (8*f_{j,3}*h_j,i*h_j,k*h_j,l
              + 4*f_{j,2}*(delta_ik*b_j,i*h_j,l
                           + delta_il*b_j,i*h_j,k + delta_kl*b_j,k*h_j,i))

Jacobians and induced reference metric

The following tensors are assembled by the generated BHaH geometry functions,
not stored as symbolic metric arrays by this class. Their components refer
to the computational basis, evaluated at the recovered physical position.

    J_ai = dq[a]/dx[i] = lam*delta_ai + x[a]*v_i
    d    = lam + sum_b x[b]*v_b
    K_ai = dx[a]/dq[i] = (delta_ai-x[a]*v_i/d)/lam
    ghat_ij = sum_a K_ai*K_aj
    ghat^ij = sum_a J_ia*J_ja
    det(ghat) = det(K)**2 = 1/(lam**4*d**2)

In particular J(0) = c*I, K(0) = I/c, and ghat_ij(0) = delta_ij/c**2.
Metric derivatives must be with respect to q, not x. Define physical
map derivatives F2, F3 and computational inverse-map derivatives H, T:

    F2_abc  = d2 q[a]/(dx[b] dx[c])
            = delta_ab*v_c + delta_ac*v_b + x[a]*W_bc
    F3_abcd = d3 q[a]/(dx[b] dx[c] dx[d])
            = delta_ab*W_cd + delta_ac*W_bd + delta_ad*W_bc + x[a]*U_bcd
    H_aij   = d2 x[a]/(dq[i] dq[j])
            = -sum_{b,c,d} K_ab*F2_bcd*K_ci*K_dj
    T_aijk  = d3 x[a]/(dq[i] dq[j] dq[k])
            = -sum_{b,c,d,e} K_ab*F3_bcde*K_ci*K_dj*K_ek
              -sum_{b,c,d} K_ab*F2_bcd*(H_cij*K_dk + H_cik*K_dj + H_cjk*K_di)

Differentiating the metric gives the quantities stored in AUXEVOL:

    ghatDDdD[i][j][k] = sum_a (H_aik*K_aj + K_ai*H_ajk)
    ghatDDdDD[i][j][k][l] = sum_a (T_aikl*K_aj + H_aik*H_ajl
                                  + H_ail*H_ajk + K_ai*T_ajkl)

The physical coordinate-line spacing used by the BHaH spacing helper is:

    ds_i = sqrt(ghat_ii)*abs(Delta q[i]),    ds_min = min_i ds_i

This is a coordinate-line spacing criterion, not an assumption that the
spheroidal computational coordinates are orthogonal.

Optional conservative timestep bound
-----------------------------------
Let delta_j = 1/a_(j-1) - 1/a_j, k_j = R_j/Rz_j, and

    t_j = (k_j + 1/k_j)/2
    C = c*(1/a_n + sum_j delta_j*t_j)
    ds_bound = min_i abs(Delta q[i])/C
    dt_bound = CFL_FACTOR*ds_bound

This avoids sampling the whole grid. To see why it is conservative, set
B_j = diag(1, 1, k_j^2), f_j = G_j(rho_j)/rho_j, and
P_j = x*(B_j*x)^T/rho_j^2. The radial kernel satisfies
0 <= G'_j <= f_j <= 1. Since P_j is a rank-one projection,
||I-P_j||_2 <= t_j. The layer Jacobian can be written as
G'_j*I + (f_j-G'_j)*(I-P_j), and its norm is bounded by t_j.
Consequently ||J||_2 <= C, so ||K*e_i||_2 >= 1/C and ds_bound <= ds_i
at every point, with the origin covered by continuity. The bound may choose
a smaller timestep than the sampled spacing rule; it does not replace a
stability analysis of the evolution scheme.

Author: Nishita Jadoo
        njadoo **at** uidaho **dot* edu

Doctests:
>>> model = GeneralRFMInverseFisheye(1)
>>> x, y, z = model.Cart
>>> unit = {model.a_list[1]: sp.Integer(1)}
>>> [sp.simplify(v.xreplace(unit)) for v in model.Cart_to_xxU]
[Cart0*inverse_fisheye_c, Cart1*inverse_fisheye_c, Cart2*inverse_fisheye_c]
>>> c, Lz = model.normalization(sp.Integer(10))
>>> sp.simplify(c.xreplace(unit)), sp.simplify(Lz.xreplace(unit))
(1, 10)
>>> spherical = {model.Rz_list[0]: model.R_list[0]}
>>> sp.simplify((Lz - 10).xreplace(spherical))
0
>>> [v.xreplace({x: 0, y: 0, z: 0}) for v in model.Cart_to_xxU]
[0, 0, 0]
>>> q = sp.symbols("q0:3", real=True)
>>> coordinate = InverseFisheyeComponent(2, 0, 1, 0, *q)
>>> all(sp.diff(coordinate, q[j]) == InverseFisheyeComponent(2, 1, 1, j, *q) for j in range(3))
True
>>> GeneralRFMInverseFisheye(0)
Traceback (most recent call last):
...
ValueError: num_transitions must be a positive integer.
>>> import math
>>> import mpmath as mp
>>> rho = sp.Symbol("rho", real=True, nonnegative=True)
>>> f3 = model.kernel_coefficients(rho, 0)[3]
>>> evaluated = sp.lambdify((rho, model.R_list[0], model.w_list[0]), f3,
...                         [{"log1p": math.log1p, "expm1": math.expm1}, "math"])
>>> mp.mp.dps = 70
>>> def reference(q):
...     r = mp.sqrt(q)
...     return mp.mpf('.25') * mp.log(mp.cosh(4 * (r + 5)) / mp.cosh(4 * (r - 5))) / (2 * mp.tanh(20) * r)
>>> expected = mp.diff(reference, mp.mpf('.25'), 3)
>>> abs((mp.mpf(evaluated(.5, 5, .5)) - expected) / expected) < mp.mpf('1e-12')
True
"""

from typing import List, Tuple, cast

import sympy as sp
from sympy.codegen.cfunctions import log1p
from sympy.printing.c import C99CodePrinter

import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.equations.generalrfm.fisheye import _G_kernel


class InverseFisheyeComponent(sp.Function):
    """
    Represent a numerically evaluated map or Jacobian component in generated C.

    Arguments are N, kind (0: position, 1: dx/dq, 2: dq/dx), i, j, q0, q1, q2.
    Infrastructure registers the named helper before emitting a BHaH project.
    This avoids pretending the implicit forward map has an analytic inverse.
    """

    def fdiff(self, argindex: int = 1) -> sp.Expr:
        """
        Differentiate a physical coordinate using the numerical inverse Jacobian.

        :param argindex: One-based argument number selected by SymPy.
        :return: The corresponding Jacobian component or an unevaluated derivative.
        """
        n, kind, i, _, *xx = self.args
        if kind == 0 and argindex in (5, 6, 7):
            return cast(sp.Expr, InverseFisheyeComponent(n, 1, i, argindex - 5, *xx))
        return cast(sp.Expr, super().fdiff(argindex))

    def _ccode(self, printer: C99CodePrinter) -> str:
        """
        Print the call to the provider's generated numerical component helper.

        :param printer: Active SymPy C printer.
        :return: C expression evaluating the requested component.
        """
        n, kind, i, j, *xx = self.args
        args = ", ".join(printer.doprint(arg) for arg in [*xx, kind, i, j])
        return f"inverse_fisheye_component__N{n}(params, {args})"


class GeneralRFMInverseFisheye:
    """Build the explicit physical-to-computational map and its local geometry."""

    def __init__(self, num_transitions: int) -> None:
        """
        Register physical transition parameters and build the inverse map.

        :param num_transitions: Positive number of spheroidal transitions.
        :raises ValueError: If the transition count is not a positive integer.
        """
        if (
            not isinstance(num_transitions, int)
            or isinstance(num_transitions, bool)
            or num_transitions < 1
        ):
            raise ValueError("num_transitions must be a positive integer.")
        self.num_transitions = num_transitions
        self.a_list = [sp.Integer(1)] + par.register_CodeParameters(
            "REAL",
            __name__,
            [f"inverse_fisheye_a{i+1}" for i in range(num_transitions)],
            [1.0] * num_transitions,
            commondata=False,
            add_to_parfile=False,
        )
        self.R_list = par.register_CodeParameters(
            "REAL",
            __name__,
            [f"inverse_fisheye_R{i+1}" for i in range(num_transitions)],
            [float(i + 1) for i in range(num_transitions)],
            commondata=False,
            add_to_parfile=False,
        )
        self.Rz_list = par.register_CodeParameters(
            "REAL",
            __name__,
            [f"inverse_fisheye_Rz{i+1}" for i in range(num_transitions)],
            [float(i + 1) for i in range(num_transitions)],
            commondata=False,
            add_to_parfile=False,
        )
        self.w_list = par.register_CodeParameters(
            "REAL",
            __name__,
            [f"inverse_fisheye_w{i+1}" for i in range(num_transitions)],
            [0.5] * num_transitions,
            commondata=False,
            add_to_parfile=False,
        )
        self.c = par.register_CodeParameter(
            "REAL",
            __name__,
            "inverse_fisheye_c",
            defaultvalue=1.0,
            commondata=False,
            add_to_parfile=False,
        )
        self.Cart = list(ixp.declarerank1("Cart", dimension=3))
        self.scale = self.c * self.scale_unscaled(self.Cart)
        self.Cart_to_xxU = [self.scale * x for x in self.Cart]

    def scale_unscaled(self, Cart: List[sp.Expr]) -> sp.Expr:
        """
        Return the unscaled factor multiplying physical Cartesian coordinates.

        :param Cart: Physical Cartesian position relative to the patch center.
        :return: Scale with its analytic value at the origin.
        """
        scale: sp.Expr = 1 / self.a_list[-1]
        for j, (R, Rz, w) in enumerate(zip(self.R_list, self.Rz_list, self.w_list)):
            rho = sp.sqrt(Cart[0] ** 2 + Cart[1] ** 2 + (R / Rz) ** 2 * Cart[2] ** 2)
            kernel = (
                w
                * sp.log(sp.cosh(2 * (rho + R) / w) / sp.cosh(2 * (rho - R) / w))
                / (4 * sp.tanh(2 * R / w))
            )
            ratio = sp.Piecewise((sp.Integer(1), sp.Eq(rho, 0)), (kernel / rho, True))
            scale += (1 / self.a_list[j] - 1 / self.a_list[j + 1]) * ratio
        return scale

    def normalization(self, physical_extent: sp.Expr) -> Tuple[sp.Expr, sp.Expr]:
        """
        Compute c and the polar computational half-width from physical extent L.

        With S(x) = scale_unscaled(x), enforce q[0](L, 0, 0) = L:

            c = 1/S(L, 0, 0)
            Lhat_z = c*L*S(0, 0, L)

        Thus the computational half-widths (L, L, Lhat_z) map to physical
        axis endpoints at +/-L in all three directions.

        :param physical_extent: Positive physical axis extent L.
        :return: Normalization c and computational polar extent Lz.
        """
        zero = sp.Integer(0)
        equatorial = self.scale_unscaled([physical_extent, zero, zero])
        polar = self.scale_unscaled([zero, zero, physical_extent])
        return 1 / equatorial, physical_extent * polar / equatorial

    def kernel_coefficients(self, radius: sp.Symbol, layer: int) -> List[sp.Expr]:
        """
        Return G(rho)/rho and three derivatives with respect to rho squared.

        A Taylor branch handles the removable origin singularity and cancellation
        in derivatives. Below the transition center, the outer branch computes
        G-rho and G1-1 directly; beyond it, the branch uses a stable log-cosh
        difference. The Taylor polynomial includes powers through rho**24.

        :param radius: Nonnegative ellipsoidal radius symbol.
        :param layer: Zero-based transition index.
        :return: Four coefficients for assembling Cartesian scale derivatives.
        """
        R, w = self.R_list[layer], self.w_list[layer]
        s = w / 2
        y = sp.Symbol("y", positive=True)
        z = sp.exp(-2 * R / s)
        th = sp.tanh(R / s)
        # d^n tanh(u)/du^n = y*A_n(y)/(1+y)^(n+1), y=exp(-2*u).
        # For even n, A_n is divisible by (1-y), which cancels tanh(u).
        polynomial = sp.Poly(4, y)
        coefficients = [sp.Integer(1)]
        for order in range(1, 25):
            if order % 2 == 0:
                quotient = polynomial.exquo(sp.Poly(1 - y, y))
                normalized = quotient.set_domain(sp.QQ).mul_ground(
                    1 / sp.factorial(order + 1)
                )
                coefficients.append(
                    z
                    * sp.horner(normalized.as_expr(), y).xreplace({y: z})
                    / (1 + z) ** order
                )
            polynomial = -2 * (
                sp.Poly(1 - order * y, y) * polynomial
                + sp.Poly(y * (1 + y), y) * polynomial.diff()
            )
        # Stable value; unlike differentiating Min/Max this keeps derivatives smooth.
        G = _G_kernel(radius, R, s)
        tp, tm = sp.tanh((radius + R) / s), sp.tanh((radius - R) / s)
        G1 = (tp - tm) / (2 * th)

        def sech_squared(argument: sp.Expr) -> sp.Expr:
            e = sp.exp(-2 * sp.Abs(argument))
            value: sp.Expr = 4 * e / (1 + e) ** 2
            return value

        sp_ = sech_squared((radius + R) / s)
        sm = sech_squared((radius - R) / s)
        G2 = (sp_ - sm) / (2 * s * th)
        G3 = (sm * tm - sp_ * tp) / (s * s * th)
        a = sp.exp(-2 * (R + radius) / s)
        b = sp.exp(-2 * (R - radius) / s)
        G_minus_radius = (2 * radius * z + s * (1 + z) * (log1p(a) - log1p(b)) / 2) / (
            1 - z
        )
        G1_minus_one = (2 * z - (1 + z) * (a / (1 + a) + b / (1 + b))) / (1 - z)
        exact = [
            G / radius,
            (radius * G1 - G) / (2 * radius**3),
            (radius**2 * G2 - 3 * radius * G1 + 3 * G) / (4 * radius**5),
            (radius**3 * G3 - 6 * radius**2 * G2 + 15 * radius * G1 - 15 * G)
            / (8 * radius**7),
        ]
        inside = [
            1 + G_minus_radius / radius,
            (radius * G1_minus_one - G_minus_radius) / (2 * radius**3),
            (radius**2 * G2 - 3 * radius * G1_minus_one + 3 * G_minus_radius)
            / (4 * radius**5),
            (
                radius**3 * G3
                - 6 * radius**2 * G2
                + 15 * radius * G1_minus_one
                - 15 * G_minus_radius
            )
            / (8 * radius**7),
        ]
        cutoff = sp.Min(2 * s, sp.Max(7 * s / 20, R / 4))
        q = sp.Symbol("q", nonnegative=True)
        result: List[sp.Expr] = []
        for order in range(4):
            series: sp.Expr = sp.Integer(0)
            for index in reversed(range(order, len(coefficients))):
                coefficient = (
                    coefficients[index]
                    * sp.factorial(index)
                    / sp.factorial(index - order)
                    / s ** (2 * order)
                )
                series = series * q + coefficient
            series = series.xreplace({q: (radius / s) ** 2})
            result.append(
                sp.Piecewise(
                    (series, radius < cutoff),
                    (inside[order], radius < R),
                    (exact[order], True),
                )
            )
        return result


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: {results.attempted} test(s)")
