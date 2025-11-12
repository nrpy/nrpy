"""
Custom JAX printer to handle custom power simplification.
This module performs the same power simplification as custom_c_codegen_functions.py for double

Author: Siddharth Mahesh; sm0193 **at** mix **dot* wvu **dot* edu
"""

from typing import Any, Union

import sympy as sp
from sympy.printing.numpy import JaxPrinter


# adding a type ignore as mypy does not let me inherit from JaxPrinter
class NRPyJaxPrinter(JaxPrinter):  # type: ignore
    """Custom JAX printer to handle custom power simplification."""

    def __init__(self) -> None:
        """Initialize the NRPyJaxPrinter."""
        super().__init__()

    def _print_Pow(self, expr: sp.Basic, rational: bool = False) -> Union[str, Any]:
        """
        Print a power expression.

        :param expr: Power expression to print.
        :param rational: Boolean indicating whether to use rational exponents.
        :return: String representation of the power expression.
        """
        base, exp = expr.as_base_exp()
        b = self._print(base)

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
            return f"jax.numpy.sqrt({b})"
        if exp == -sp.Rational(1, 2) or (
            getattr(exp, "is_Float", False) and float(exp) == -0.5
        ):
            return f"(1.0/jax.numpy.sqrt({b}))"
        if exp == sp.Rational(1, 3):
            return f"jax.numpy.cbrt({b})"
        if exp == -sp.Rational(1, 3):
            return f"(1.0/jax.numpy.cbrt({b}))"
        if exp == sp.Rational(1, 6):
            return f"jax.numpy.sqrt(jax.numpy.cbrt({b}))"
        if exp == -sp.Rational(1, 6):
            return f"(1.0/(jax.numpy.sqrt(jax.numpy.cbrt({b}))))"

        # Small integer powers mapped to repeated multiplication
        if isinstance(exp, sp.Integer):
            n = int(exp)
            if n in (2, 3, 4, 5):
                return muln(n)
            if n in (-1, -2, -3, -4, -5):
                return f"(1.0/({muln(-n)}))"

        # Fallback to default JAX handling
        return super()._print_Pow(expr, rational=rational)
