"""
Construct symbolic expressions for GRHD shock-test initial data.

These expressions are intended for high-resolution shock-capturing test
problems and provide symbolic density and pressure profiles that can be fed
into NRPy-generated infrastructure code.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from typing import Dict, Tuple

import sympy as sp

import nrpy.equations.grhd.Min_Max_and_Piecewise_Expressions as noif


def balsara0(x: sp.Expr, bound: sp.Expr = sp.Integer(0)) -> Tuple[sp.Expr, sp.Expr]:
    """
    Construct the Balsara 0 shock-tube profile.

    :param x: Coordinate aligned with the discontinuity.
    :param bound: Position of the shock interface.
    :return: Density and pressure expressions.
    """
    # Step 1: Set the left and right primitive states.
    rho_left = sp.Integer(1)
    rho_right = sp.Rational(1, 8)
    press_left = sp.Integer(1)
    press_right = sp.Rational(1, 10)

    # Step 2: Build the symbolic piecewise profile.
    rho = (
        noif.coord_less_bound(x, bound) * rho_left
        + noif.coord_geq_bound(x, bound) * rho_right
    )
    press = (
        noif.coord_less_bound(x, bound) * press_left
        + noif.coord_geq_bound(x, bound) * press_right
    )
    return rho, press


def hydro1(x: sp.Expr, bound: sp.Expr = sp.Rational(1, 2)) -> Tuple[sp.Expr, sp.Expr]:
    """
    Construct the Hydro 1 shock-tube profile.

    :param x: Coordinate aligned with the discontinuity.
    :param bound: Position of the shock interface.
    :return: Density and pressure expressions.
    """
    # Step 1: Set the left and right primitive states.
    rho_left = sp.Integer(10)
    rho_right = sp.Integer(1)
    press_left = sp.Rational(40, 3)
    press_right = sp.Integer(0)

    # Step 2: Build the symbolic piecewise profile.
    rho = (
        noif.coord_less_bound(x, bound) * rho_left
        + noif.coord_geq_bound(x, bound) * rho_right
    )
    press = (
        noif.coord_less_bound(x, bound) * press_left
        + noif.coord_geq_bound(x, bound) * press_right
    )
    return rho, press


def hydro2(x: sp.Expr, bound: sp.Expr = sp.Rational(1, 2)) -> Tuple[sp.Expr, sp.Expr]:
    """
    Construct the Hydro 2 shock-tube profile.

    :param x: Coordinate aligned with the discontinuity.
    :param bound: Position of the shock interface.
    :return: Density and pressure expressions.
    """
    # Step 1: Set the left and right primitive states.
    rho = sp.Integer(1)
    press_left = sp.Integer(1000)
    press_right = sp.Rational(1, 100)

    # Step 2: Build the symbolic piecewise profile.
    press = (
        noif.coord_less_bound(x, bound) * press_left
        + noif.coord_geq_bound(x, bound) * press_right
    )
    return rho, press


def TMM(
    theta: sp.Expr,
    bound: sp.Expr = sp.pi / 2,
    kappa: sp.Expr = sp.Integer(1),
    gamma: sp.Expr = sp.Rational(6, 5),
) -> Tuple[sp.Expr, sp.Expr]:
    """
    Construct the TMM angular shock profile.

    :param theta: Polar angle used to locate the discontinuity.
    :param bound: Position of the angular interface.
    :param kappa: Polytropic constant used to set the pressure.
    :param gamma: Polytropic exponent used to set the pressure.
    :return: Density and pressure expressions.
    """
    # Step 1: Set the left and right primitive states.
    rho_left = sp.Rational(1, 100000000)
    rho_right = sp.Rational(1, 10000000)

    # Step 2: Build the symbolic piecewise profile.
    rho = (
        noif.coord_less_bound(theta, bound) * rho_left
        + noif.coord_greater_bound(theta, bound) * rho_right
    )
    press = kappa * rho**gamma
    return rho, press


def cylindrical_explosion(
    r: sp.Expr,
    r_in: sp.Expr = sp.Rational(4, 5),
    r_out: sp.Expr = sp.Integer(1),
) -> Tuple[sp.Expr, sp.Expr]:
    """
    Construct the cylindrical explosion profile.

    :param r: Radial coordinate used to define the explosion profile.
    :param r_in: Inner transition radius.
    :param r_out: Outer transition radius.
    :return: Density and pressure expressions.
    """
    # Step 1: Precompute the radial interpolation factors.
    r_out_minus_r_in = r_out - r_in
    r_out_minus_r = r_out - r
    r_minus_r_in = r - r_in

    # Step 2: Interpolate the density profile through the transition region.
    rho_in = sp.Rational(1, 100)
    rho_out = sp.Rational(1, 10000)
    rho_mid = sp.exp(
        (r_out_minus_r * sp.log(rho_in) + r_minus_r_in * sp.log(rho_out))
        / r_out_minus_r_in
    )
    rho = (
        noif.coord_leq_bound(r, r_in) * rho_in
        + noif.coord_greater_bound(r, r_in) * noif.coord_less_bound(r, r_out) * rho_mid
        + noif.coord_geq_bound(r, r_out) * rho_out
    )

    # Step 3: Interpolate the pressure profile through the transition region.
    press_in = sp.Integer(1)
    press_out = sp.Rational(3, 100000)
    press_mid = sp.exp(
        (r_out_minus_r * sp.log(press_in) + r_minus_r_in * sp.log(press_out))
        / r_out_minus_r_in
    )
    press = (
        noif.coord_leq_bound(r, r_in) * press_in
        + noif.coord_greater_bound(r, r_in)
        * noif.coord_less_bound(r, r_out)
        * press_mid
        + noif.coord_geq_bound(r, r_out) * press_out
    )
    return rho, press


class GRHD_ShockTests:
    """
    Construct symbolic primitive data for supported GRHD shock tests.

    :param IDType: Shock-test initial data type.
    :param x: Coordinate aligned with one-dimensional shock interfaces.
    :param r: Radial coordinate for cylindrical-explosion data.
    :param bound: Location of the discontinuity.
    :param theta: Angular coordinate for the TMM profile.
    :param kappa: Polytropic constant for the TMM profile.
    :param gamma: Polytropic exponent for the TMM profile.
    :param r_in: Inner transition radius for the cylindrical explosion.
    :param r_out: Outer transition radius for the cylindrical explosion.
    :raises ValueError: If the requested shock-test type is unsupported.
    """

    def __init__(
        self,
        IDType: str,
        x: sp.Expr,
        r: sp.Expr,
        bound: sp.Expr = sp.Integer(1),
        theta: sp.Expr = sp.Integer(0),
        kappa: sp.Expr = sp.Integer(1),
        gamma: sp.Expr = sp.Rational(6, 5),
        r_in: sp.Expr = sp.Rational(4, 5),
        r_out: sp.Expr = sp.Integer(1),
    ) -> None:
        # Step 1: Select the requested shock-test profile.
        if IDType == "balsara0":
            rho, press = balsara0(x=x, bound=bound)
        elif IDType == "hydro1":
            rho, press = hydro1(x=x, bound=bound)
        elif IDType == "hydro2":
            rho, press = hydro2(x=x, bound=bound)
        elif IDType == "TMM":
            rho, press = TMM(
                theta=theta,
                bound=bound,
                kappa=kappa,
                gamma=gamma,
            )
        elif IDType == "cylindrical_explosion":
            rho, press = cylindrical_explosion(r=r, r_in=r_in, r_out=r_out)
        else:
            raise ValueError(
                f"Shock-test initial data type IDType={IDType} is not supported."
            )

        # Step 2: Store the symbolic primitive variables.
        self.rho = rho
        self.press = press


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

    # Step 1: Build representative symbolic expressions for each supported test.
    exprs_dict: Dict[str, sp.Expr] = {}
    x_sym, r_sym, theta_sym = sp.symbols("x r theta", real=True)
    shock_tests: Dict[str, GRHD_ShockTests] = {
        "balsara0": GRHD_ShockTests(
            IDType="balsara0",
            x=x_sym,
            r=r_sym,
            bound=sp.Rational(1, 3),
        ),
        "hydro1": GRHD_ShockTests(
            IDType="hydro1",
            x=x_sym,
            r=r_sym,
            bound=sp.Rational(2, 5),
        ),
        "hydro2": GRHD_ShockTests(
            IDType="hydro2",
            x=x_sym,
            r=r_sym,
            bound=sp.Rational(3, 5),
        ),
        "TMM": GRHD_ShockTests(
            IDType="TMM",
            x=x_sym,
            r=r_sym,
            bound=sp.pi / 4,
            theta=theta_sym,
            kappa=sp.Rational(7, 5),
            gamma=sp.Rational(4, 3),
        ),
        "cylindrical_explosion": GRHD_ShockTests(
            IDType="cylindrical_explosion",
            x=x_sym,
            r=r_sym,
            r_in=sp.Rational(3, 4),
            r_out=sp.Rational(9, 8),
        ),
    }

    # Step 2: Replace the symbolic no-if absolute-value helper for validation.
    for shock_test_name, shock_test in shock_tests.items():
        exprs_dict[f"{shock_test_name}_rho"] = shock_test.rho.subs(
            sp.Function("nrpyAbs"), sp.Abs
        )
        exprs_dict[f"{shock_test_name}_press"] = shock_test.press.subs(
            sp.Function("nrpyAbs"), sp.Abs
        )

    # Step 3: Compare against trusted expressions, or generate them if needed.
    results_dict = ve.process_dictionary_of_expressions(
        exprs_dict, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        f"{os.path.splitext(os.path.basename(__file__))[0]}",
        results_dict,
    )
