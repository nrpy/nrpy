"""
The Goldberg formula for computing spin-weighted spherical harmonics is implemented here.
  (https://aip.scitation.org/doi/10.1063/1.1705135)
Wikipedia also has an article on Spin-Weighted Spherical Hamronics:
  (https://en.wikipedia.org/w/index.php?title=Spin-weighted_spherical_harmonics&oldid=853425244)

Authors: Brandon Clark
         Zachariah B. Etienne
         zachetie **at** gmail **dot* com
"""
from typing import cast, Dict
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends

from nrpy.helpers.cached_functions import cached_simplify


def Y(
    s: int,
    l: int,
    m: int,
    th: sp.Symbol,
    ph: sp.Symbol,
    GenerateMathematicaCode: bool = False,
) -> sp.Expr:
    """
    Compute the spin-weighted spherical harmonics using the Goldberg formula.
    Refer to https://aip.scitation.org/doi/10.1063/1.1705135 and Wikipedia's article on
    Spin-weighted spherical harmonics:
    https://en.wikipedia.org/w/index.php?title=Spin-weighted_spherical_harmonics&oldid=853425244

    :param s: Spin weight of the spherical harmonic.
    :param l: Degree of the spherical harmonic.
    :param m: Order of the spherical harmonic.
    :param th: Theta symbol.
    :param ph: Phi symbol.
    :param GenerateMathematicaCode: Boolean to control whether to generate code compatible with Mathematica.

    :return: Expression for the spin-weighted spherical harmonic.
    """
    Sum: sp.Expr = sp.sympify(0)
    # Not sympifying -1 causes SymPy to revert to double-precision arithmetic
    negative_one = sp.sympify(-1)
    for r in range(l - s + 1):
        if GenerateMathematicaCode:
            # Mathematica needs expression to be in terms of cotangent, so that code validation below
            #    yields identity with existing Mathematica notebook on spin-weighted spherical harmonics.
            Sum += (
                sp.binomial(l - s, r)
                * sp.binomial(l + s, r + s - m)
                * negative_one ** (l - r - s)
                * sp.exp(sp.I * m * ph)
                * sp.cot(th / 2) ** (2 * r + s - m)
            )
        else:
            # SymPy C code generation cannot handle the cotangent function, so define cot(th/2) as 1/tan(th/2):
            Sum += (
                sp.binomial(l - s, r)
                * sp.binomial(l + s, r + s - m)
                * negative_one ** (l - r - s)
                * sp.exp(sp.I * m * ph)
                / sp.tan(th / 2) ** (2 * r + s - m)
            )
    return cast(
        sp.Expr,
        sp.sqrt(
            sp.factorial(l + m)
            * sp.factorial(l - m)
            * (2 * l + 1)
            / (4 * sp.pi * sp.factorial(l + s) * sp.factorial(l - s))
        )
        * cached_simplify(negative_one**m * sp.sin(th / 2) ** (2 * l) * Sum),
    )


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

    expr_dict: Dict[str, sp.Expr] = {}
    _th, _ph = sp.symbols("th ph", real=True)
    for _l in range(2, 9):
        for _m in range(-_l, +_l):
            expr_dict[f"Y_{{s=-2, l={_l}, m={_m}}}"] = Y(-2, _l, _m, _th, _ph)
    results_dict = ve.process_dictionary_of_expressions(
        expr_dict, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        # File basename. If this is set to "trusted_module_test1", then
        #   trusted results_dict will be stored in tests/trusted_module_test1.py
        f"{os.path.splitext(os.path.basename(__file__))[0]}",
        results_dict,
    )
