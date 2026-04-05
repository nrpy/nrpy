# BHaHAHA/area.py
"""
Calculate various geometric properties of a surface embedded in a 3D spherical space.

Include functions to compute the area, differential area elements, induced metrics, and differential arclength elements for a surface with a radius function h(theta, phi). Express geometric properties symbolically using SymPy for computations.

Functions:
- area: Calculates a symbolic expression for the area element of the surface.
- area2: Computes an alternate symbolic expression for the area element.
- compute_q2DD: Computes the 2D induced metric tensor q2DD for the surface parameterized by (theta, phi).
- area3: Calculates the differential area element using the induced metric.
- spin_NewtonRaphson: Computes the Newton-Raphson iteration for spin calculations.
- spin_HalleysMethod: Computes Halley's method iteration for spin calculations.
- circumferential_arclength: Calculates the differential arclength element along a specified direction (theta or phi).

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List, cast

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.reference_metric as refmetric
from nrpy.equations.general_relativity.bhahaha.ExpansionFunctionTheta import (
    ExpansionFunctionTheta,
)


def area() -> sp.Expr:
    """
    Calculate the area element for a surface with a radius function h(theta, phi).

    This is Eq 33 of Lin & Novak, https://arxiv.org/pdf/gr-qc/0702038v2

    DocTests:
    >>> expr = area()
    Setting up reference_metric[Spherical]...
    Setting up ExpansionFunctionThetaClass[Spherical]...
    >>> isinstance(expr, sp.Expr)
    True

    :return: The symbolic expression for the area element.
    """
    rfm = refmetric.reference_metric["Spherical"]
    Th = ExpansionFunctionTheta["Spherical"]
    physgammahatDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            physgammahatDD[i][j] = Th.gammaDD[i][j] / rfm.ReDD[i][j]
    h = sp.Symbol("hh", real=True)
    h_dD = ixp.declarerank1("hh_dD")
    term1a = (
        physgammahatDD[0][0] * h_dD[1] ** 2
        + 2 * physgammahatDD[0][1] * h * h_dD[1]
        + physgammahatDD[1][1] * h**2
    )
    term1b = (
        physgammahatDD[0][0] * h_dD[2] ** 2
        + 2 * physgammahatDD[0][2] * h * h_dD[2] * sp.sin(rfm.xx[1])
        + physgammahatDD[2][2] * h**2 * sp.sin(rfm.xx[1]) ** 2
    )
    term2 = (
        physgammahatDD[0][0] * h_dD[1] * h_dD[2]
        + physgammahatDD[0][1] * h * h_dD[2]
        + physgammahatDD[0][2] * h * h_dD[1] * sp.sin(rfm.xx[1])
        + physgammahatDD[1][2] * h**2 * sp.sin(rfm.xx[1])
    )
    area_element = sp.sqrt(term1a * term1b - term2**2).replace(sp.sympify("xx0"), h)
    return cast(sp.Expr, area_element)


def area2() -> sp.Expr:
    """
    Calculate an alternate symbolic expression for the area element of the surface.

    This is Eq 24 of Hui & Lin, https://arxiv.org/pdf/2404.16511

    DocTests:
    >>> expr = area2()
    >>> isinstance(expr, sp.Expr)
    True

    :return: The symbolic expression for the area element.
    """
    Th = ExpansionFunctionTheta["Spherical"]
    h = sp.Symbol("hh", real=True)
    h_dD = ixp.declarerank1("hh_dD")
    term1a = (
        Th.gammaDD[0][0] * h_dD[1] ** 2
        + 2 * Th.gammaDD[0][1] * h_dD[1]
        + Th.gammaDD[1][1]
    )
    term1b = (
        Th.gammaDD[0][0] * h_dD[2] ** 2
        + 2 * Th.gammaDD[0][2] * h_dD[2]
        + Th.gammaDD[2][2]
    )
    term2 = (
        Th.gammaDD[0][0] * h_dD[1] * h_dD[2]
        + Th.gammaDD[0][1] * h_dD[2]
        + Th.gammaDD[0][2] * h_dD[1]
        + Th.gammaDD[1][2]
    )
    area_element = sp.sqrt((term1a * term1b - term2**2).replace(sp.sympify("xx0"), h))
    return cast(sp.Expr, area_element)


def compute_q2DD() -> List[List[sp.Expr]]:
    """
    Compute the 2D induced metric tensor q2DD for the surface parameterized by (theta, phi).

    The induced metric q2DD is constructed by considering how the radius function h(theta, phi)
    induces a metric on the 2D surface within the original 3D metric gamma_{ij} (in the spherical basis).

    The surface is embedded in 3D space with coordinates (r, theta, phi), where r = h(theta, phi).
    The position vector on the surface is parameterized by:

        p^k = (h(theta, phi), theta, phi)

    To derive the induced metric, we need to compute how the changes in the parameters theta and phi
    are translated into changes on the surface.

    The partial derivatives of the position vector p^k with respect to the parameters (theta, phi) are:

        partial p^k / partialtheta and partial p^k / partialphi

    These partial derivatives are used to compute the induced metric tensor q_{ij} on the surface, which is given by:

        q_ij = (partial p^k / partial x^i) * (partial p^l / partial x^j) * gamma_kl

    where x^i represents the surface coordinates (theta, phi), and gamma_{kl} are the components of the 3D spherical metric.

    In this code, the induced metric tensor qDD is first computed in the full 3D space (including r, theta, phi).
    Since we are only interested in the surface parameterized by (theta, phi), we extract the 2x2 sub-matrix q2DD
    that corresponds to the angular components (theta, phi). This induced 2-metric q2DD is constructed as follows:

    1. Calculate the derivatives of the position vector pU with respect to the coordinates (r, theta, phi).
    2. Use these derivatives to construct the induced metric qDD by summing over the ambient metric components gamma_{kl}.
    3. Extract the 2x2 sub-matrix q2DD from qDD, which includes only the angular components (theta, phi).

    DocTests:
    >>> q2 = compute_q2DD()
    >>> isinstance(q2, list)
    True
    >>> len(q2) == 2 and all(len(row) == 2 for row in q2)
    True

    :return: The 2x2 induced metric tensor for the surface.
    """
    rfm = refmetric.reference_metric["Spherical"]
    Th = ExpansionFunctionTheta["Spherical"]
    h = sp.Symbol("hh", real=True)
    h_dD = ixp.declarerank1("hh_dD")
    # Position vector on the surface, in spherical coordinates
    # pU = p^k = (h(theta, phi), theta, phi), where h(theta, phi) is the radius
    pU = [h, rfm.xx[1], rfm.xx[2]]

    # Derivatives of the position vector with respect to the coordinates
    # pU_dD[i][j] represents the derivative of pU[i] with respect to coordinate j
    pU_dD = ixp.zerorank2()  # Initialize a 2D array (rank-2 tensor) with all zeros

    # Setting derivatives of the radial coordinate (h) with respect to theta and phi
    pU_dD[0][0] = sp.sympify(0)  # dh/dr = 0, since h is a function of (theta, phi)
    pU_dD[0][1] = h_dD[1]  # dh/dtheta = partial derivative of h with respect to theta
    pU_dD[0][2] = h_dD[2]  # dh/dphi = partial derivative of h with respect to phi

    # Calculating the derivatives of the angular components (theta, phi)
    # with respect to the spherical coordinate basis (r, theta, phi)
    for i in range(1, 3):  # Loop over theta and phi components (i = 1, 2)
        for j in range(3):  # Loop over r, theta, phi (j = 0, 1, 2)
            # Derivative of angular component i with respect to coordinate j
            pU_dD[i][j] = sp.diff(pU[i], rfm.xx[j])

    # Initializing the induced metric tensor q_DD on the surface (3x3 tensor)
    qDD = ixp.zerorank2()

    # Calculating the induced metric tensor by summing over the ambient metric components
    # qDD[i][j] represents the induced metric on the surface
    # q_{ij} = p^k_i p^l_j gamma_{kl}
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    # Summing contributions from the ambient metric and position vector derivatives
                    qDD[i][j] += pU_dD[k][i] * pU_dD[l][j] * Th.gammaDD[k][l]

    # Extracting the 2x2 induced metric tensor q2DD for the surface parameterized by (theta, phi)
    # The induced 2-metric q2DD represents the metric on the 2D surface, obtained by focusing
    # on the angular components (theta, phi) of the full 3D induced metric qDD.
    q2DD = ixp.zerorank2(dimension=2)  # Initialize a 2x2 tensor

    # Extracting the sub-matrix corresponding to the angular components
    # Since qDD is a 3x3 matrix that includes r, theta, and phi, we need to extract
    # the components that describe the angular directions only (i.e., ignoring the radial part).
    for i in range(2):  # Loop over theta and phi components (i = 0, 1)
        for j in range(2):  # Loop over theta and phi components (j = 0, 1)
            q2DD[i][j] = qDD[i + 1][j + 1]

    return q2DD


def area3() -> sp.Expr:
    """
    Calculate the differential area element for a surface with a radius function h(theta, phi).

    The function returns the area element expressed in terms of symbolic variables.

    The area element dA is given by:

    dA = sqrt(det(q2DD)) dtheta dphi in spherical coordinates.

    where q2DD is the induced 2x2 metric tensor on the surface parameterized by (theta, phi).
    This metric tensor is derived from the embedding of the surface into the 3D spherical space.

    DocTests:
    >>> expr = area3()
    >>> isinstance(expr, sp.Expr)
    True

    :return: The symbolic expression for the differential area element.
    """
    q2DD = compute_q2DD()
    # Computing the determinant of the 2x2 induced metric tensor q2DD, the _ variable is q2UU, ignored.
    _, q2det = ixp.symm_matrix_inverter2x2(q2DD)

    # Returning the square root of the determinant, which represents the area element
    # Replace the placeholder variable "xx0" with the actual radius function h
    area_element = sp.sqrt(q2det.replace(sp.sympify("xx0"), sp.Symbol("hh")))
    return cast(sp.Expr, area_element)


def spin_NewtonRaphson() -> sp.Expr:
    """
    Compute the Newton-Raphson iteration for spin calculations based on Eq 5.2 from Alcubierre et al arXiv:gr-qc/0411149.

    DocTests:
    >>> expr = spin_NewtonRaphson()
    >>> isinstance(expr, sp.Expr)
    True

    :return: The Newton-Raphson iteration formula as a symbolic expression.
    """
    # Define symbols for calculation
    x, Cr, E, K = sp.symbols("x C_r E K")

    # Define the function f(x) based on Eq 5.2 of Alcubierre et al arXiv:gr-qc/0411149,
    # where x = (a/m)^2. Given Cr, the spin x is found when this expression equals zero.
    f_of_x = ((1 + sp.sqrt(1 - x**2)) / sp.pi) * sp.elliptic_e(
        -(x**2) / (1 + sp.sqrt(1 - x**2)) ** 2
    ) - Cr

    # Differentiate the function f(x) with respect to x and simplify
    fprime_of_x = sp.diff(f_of_x, x).simplify()

    # Substitute elliptic functions with symbols E and K.
    # In the C code, E(x) & K(x) are computed and stored in symbols E & K, respectively.
    fprime_of_x = fprime_of_x.subs(
        {
            sp.elliptic_e(-(x**2) / (1 + sp.sqrt(1 - x**2)) ** 2): E,
            sp.elliptic_k(-(x**2) / (1 + sp.sqrt(1 - x**2)) ** 2): K,
        }
    )
    f_of_x = f_of_x.subs({sp.elliptic_e(-(x**2) / (1 + sp.sqrt(1 - x**2)) ** 2): E})

    # Return the Newton-Raphson iteration formula
    return cast(sp.Expr, x - f_of_x / fprime_of_x)


def spin_HalleysMethod() -> sp.Expr:
    """
    Compute Halley's method iteration for spin calculations based on Eq 5.2 from Alcubierre et al arXiv:gr-qc/0411149.

    DocTests:
    >>> expr = spin_HalleysMethod()
    >>> isinstance(expr, sp.Expr)
    True

    :return: Halley's iteration formula as a symbolic expression.
    """
    # Define symbols for calculation
    x, Cr, E, K = sp.symbols("x C_r E K")

    # Define the function f(x) based on Eq 5.2 of Alcubierre et al arXiv:gr-qc/0411149,
    # where x = (a/m)^2. Given Cr, the spin x is found when this expression equals zero.
    f_of_x = ((1 + sp.sqrt(1 - x**2)) / sp.pi) * sp.elliptic_e(
        -(x**2) / (1 + sp.sqrt(1 - x**2)) ** 2
    ) - Cr

    # Differentiate the function f(x) with respect to x and simplify
    fp_of_x = sp.diff(f_of_x, x)

    fpp_of_x = sp.diff(fp_of_x, x)

    # Substitute elliptic functions with symbols E and K.
    # In the C code, E(x) & K(x) are computed and stored in symbols E & K, respectively.
    fpp_of_x = fpp_of_x.subs(
        {
            sp.elliptic_e(-(x**2) / (1 + sp.sqrt(1 - x**2)) ** 2): E,
            sp.elliptic_k(-(x**2) / (1 + sp.sqrt(1 - x**2)) ** 2): K,
        }
    ).simplify()
    fp_of_x = fp_of_x.subs(
        {
            sp.elliptic_e(-(x**2) / (1 + sp.sqrt(1 - x**2)) ** 2): E,
            sp.elliptic_k(-(x**2) / (1 + sp.sqrt(1 - x**2)) ** 2): K,
        }
    ).simplify()
    f_of_x = f_of_x.subs(
        {sp.elliptic_e(-(x**2) / (1 + sp.sqrt(1 - x**2)) ** 2): E}
    ).simplify()

    # Return the Halley's iteration formula
    return cast(
        sp.Expr, x - 2 * f_of_x * fp_of_x / (2 * fp_of_x * fp_of_x - f_of_x * fpp_of_x)
    )


def circumferential_arclength(direction: str) -> sp.Expr:
    """
    Calculate the differential arclength element for a surface with a radius function h(theta, phi).

    The calculation is performed in spherical coordinates, along a specific circumferential direction.

    DocTests:
    >>> arclength_theta = circumferential_arclength("theta")
    >>> isinstance(arclength_theta, sp.Expr)
    True
    >>> arclength_phi = circumferential_arclength("phi")
    >>> isinstance(arclength_phi, sp.Expr)
    True
    >>> circumferential_arclength("invalid")
    Traceback (most recent call last):
    ...
    ValueError: Invalid direction specified. Use 'theta' or 'phi'.

    :param direction: The direction for the circumferential arclength. Either "theta" or "phi".
    :return: The differential arclength element in the specified direction.
    :raises ValueError: If an invalid direction is specified.
    """
    q2DD = compute_q2DD()

    if direction == "theta":
        # Differential arclength in the theta direction
        q_theta_theta = q2DD[0][0]
        arclength_expr = sp.sqrt(
            q_theta_theta.replace(sp.sympify("xx0"), sp.Symbol("hh"))
        )
    elif direction == "phi":
        # Differential arclength in the phi direction
        q_phi_phi = q2DD[1][1]
        arclength_expr = sp.sqrt(q_phi_phi.replace(sp.sympify("xx0"), sp.Symbol("hh")))
    else:
        raise ValueError("Invalid direction specified. Use 'theta' or 'phi'.")

    return cast(sp.Expr, arclength_expr)


if __name__ == "__main__":
    import doctest
    import sys

    import nrpy.validate_expressions.validate_expressions as ve

    # Run doctests
    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Validate expressions for the area functions
    ve.assert_equal(area(), area2())
    ve.assert_equal(area(), area3())
    ve.assert_equal(area2(), area3())

    # Validation expression in Mathematica:
    # f= (1 + Sqrt[1 - x ^ 2]) / Pi * EllipticE[-x^2 / (1 + Sqrt[1 - x^2])^2] - Cr; CForm[FullSimplify[x - f/D[f, x]]]
    # so
    x_main, Cr_main = sp.symbols("x C_r")
    trusted_spin_expr = (
        x_main
        * (
            -(Cr_main * sp.pi * sp.sqrt(1 - x_main**2))
            + (1 + sp.sqrt(1 - x_main**2)) * sp.Symbol("K")
        )
    ) / (
        (-1 + x_main**2 - sp.sqrt(1 - x_main**2)) * sp.Symbol("E")
        + (1 + sp.sqrt(1 - x_main**2)) * sp.Symbol("K")
    )
    ve.assert_equal(spin_NewtonRaphson(), trusted_spin_expr)
