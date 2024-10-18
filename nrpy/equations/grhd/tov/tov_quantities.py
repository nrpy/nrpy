"""
Get ADM quantities and energy-momentum tensor for a Tolman-Oppenheimer-Volkoff (TOV) star in terms of the quantities used when solving the TOV equations.

A good resource to learn more about the TOV equation is its Wikipedia entry:
https://en.wikipedia.org/wiki/Tolman-Oppenheimer-Volkoff_equation.

Author: Leonardo Rosa Werneck
        wernecklr **at** gmail **dot* com
"""

from typing import List, Tuple, Union

from sympy import Expr, Symbol, sin, sqrt

from nrpy.indexedexp import zerorank1, zerorank2


def tov_ADM_quantities(
    rbar: Union[Symbol, Expr] = Symbol("rbar", real=True),
    theta: Union[Symbol, Expr] = Symbol("theta", real=True),
    expnu: Union[Symbol, Expr] = Symbol("expnu", real=True),
    exp4phi: Union[Symbol, Expr] = Symbol("exp4phi", real=True),
) -> Tuple[
    Union[Symbol, Expr],
    List[Union[Symbol, Expr]],
    List[Union[Symbol, Expr]],
    List[List[Union[Symbol, Expr]]],
    List[List[Union[Symbol, Expr]]],
]:
    """
    Set ADM quantities in terms of the TOV solution quantities.

    :param rbar: Isotropic radius, rbar = e^{-2phi} r. Defaults to Symbol("rbar", real=True).
    :param theta: Polar angle. Defaults to Symbol("theta", real=True).
    :param expnu:  g_{tt} in the TOV solution. Defaults to Symbol("expnu", real=True).
    :param exp4phi: Conformal factor, psi = e^{4phi}. Defaults to Symbol("exp4phi", real=True).

    :return: List with ADM quantities: alpha, beta^i, B^i, gamma_{ij}, K_{ij}.
    """
    alpha = sqrt(expnu)
    betaU = zerorank1()
    BU = zerorank1()
    KDD = zerorank2()

    gammaDD = zerorank2()
    gammaDD[0][0] = exp4phi
    gammaDD[1][1] = exp4phi * rbar**2
    gammaDD[2][2] = exp4phi * rbar**2 * sin(theta) ** 2

    return alpha, betaU, BU, gammaDD, KDD


def tov_T4UU(
    rho: Union[Symbol, Expr] = Symbol("rho", real=True),
    P: Union[Symbol, Expr] = Symbol("P", real=True),
    rbar: Union[Symbol, Expr] = Symbol("rbar", real=True),
    theta: Union[Symbol, Expr] = Symbol("theta", real=True),
    expnu: Union[Symbol, Expr] = Symbol("expnu", real=True),
    exp4phi: Union[Symbol, Expr] = Symbol("exp4phi", real=True),
) -> List[List[Union[Symbol, Expr]]]:
    """
    Set the stress-energy tensor in terms of the TOV solution quantities.
    Here, T^{mu nu} = diag(rho/e^{nu}, PP, PP/r^2, PP/r^2/sin^2(theta)).

    :param rho: Fluid density. Defaults to Symbol("rho", real=True).
    :param P: Fluid pressure. Defaults to Symbol("P", real=True).
    :param rbar: Isotropic radius, rbar = e^{-2phi} r. Defaults to Symbol("rbar", real=True).
    :param theta: Polar angle. Defaults to Symbol("theta", real=True).
    :param expnu:  g_{tt} in the TOV solution. Defaults to Symbol("expnu", real=True).
    :param exp4phi: Conformal factor, psi = e^{4phi}. Defaults to Symbol("exp4phi", real=True).

    :return: Stress-energy tensor as rank-2 indexed expression.
    """
    T4UU = zerorank2(dimension=4)
    T4UU[0][0] = rho / expnu
    T4UU[1][1] = P / exp4phi
    T4UU[2][2] = P / (exp4phi * rbar**2)
    T4UU[3][3] = P / (exp4phi * rbar**2 * sin(theta) ** 2)

    return T4UU
