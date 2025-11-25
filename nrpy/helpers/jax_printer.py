"""
Custom JAX printer to handle custom power simplification.
This module performs the same power simplification as custom_c_codegen_functions.py for double

Author: Siddharth Mahesh; sm0193 **at** mix **dot* wvu **dot* edu
"""

from typing import Any, Dict, Union

import sympy as sp

try:
    from sympy.printing.numpy import JaxPrinter as Printer
    from sympy.printing.numpy import _known_constants_jax as known_constants
    from sympy.printing.numpy import _known_functions_jax as known_functions
except ImportError:
    # Fallback for older SymPy versions
    from sympy.printing.numpy import NumPyPrinter as Printer
    from sympy.printing.numpy import _known_constants_numpy, _known_functions_numpy

    known_functions = {k: "jnp." + v for k, v in _known_functions_numpy.items()}
    known_constants = {k: "jnp." + v for k, v in _known_constants_numpy.items()}


# adding a type ignore as mypy does not let me inherit from JaxPrinter.
# Disable specific pylint errors owing to sympy
# pylint: disable=too-many-ancestors, abstract-method
class NRPyJaxPrinter(Printer):  # type: ignore
    """Custom JAX printer to handle custom power simplification."""

    _module = "jnp"
    _kf = known_functions
    _kc = known_constants

    def __init__(self, settings: Union[None, Dict[str, Any]] = None) -> None:
        """
        Initialize the NRPyJaxPrinter.

        :param settings: Settings for the printer.(defaults to None)
        """
        super().__init__(settings=settings)

    def _print_Pow(self, expr: sp.Basic, rational: bool = False) -> Union[str, Any]:
        """
        Print a power expression.

        :param expr: Power expression to print.
        :param rational: Boolean indicating whether to use rational exponents.
        :return: String representation of the power expression.
        """
        base, exp = expr.as_base_exp()
        b = self._print(base)
        retval = None

        def muln(n: int) -> str:
            """
            Return a string of repeated multiplication.

            :param n: Number of times to repeat multiplication.
            :return: String of repeated multiplication.
            """
            return "(" + "*".join([f"({b})"] * n) + ")"

        # Fractional exponents mapped to sqrt/cbrt
        if exp == sp.Rational(1, 2) or (
            getattr(exp, "is_Float", False) and float(exp) == 0.5
        ):
            retval = f"{self._module}.sqrt({b})"
        if exp == -sp.Rational(1, 2) or (
            getattr(exp, "is_Float", False) and float(exp) == -0.5
        ):
            retval = f"(1.0/{self._module}.sqrt({b}))"
        if exp == sp.Rational(1, 3):
            retval = f"{self._module}.cbrt({b})"
        if exp == -sp.Rational(1, 3):
            retval = f"(1.0/{self._module}.cbrt({b}))"
        if exp == sp.Rational(1, 6):
            retval = f"{self._module}.sqrt({self._module}.cbrt({b}))"
        if exp == -sp.Rational(1, 6):
            retval = f"(1.0/({self._module}.sqrt({self._module}.cbrt({b}))))"
        # Small integer powers mapped to repeated multiplication
        if isinstance(exp, sp.Integer):
            n = int(exp)
            if n in (2, 3, 4, 5):
                retval = muln(n)
            if n in (-1, -2, -3, -4, -5):
                retval = f"(1.0/({muln(-n)}))"

        # Fallback to default JAX handling
        if retval is None:
            return super()._print_Pow(expr, rational=rational)
        return retval
