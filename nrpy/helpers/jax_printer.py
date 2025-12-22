"""
Custom JAX printer to handle custom power simplification.
This module performs the same power simplification as custom_c_codegen_functions.py for double

Author: Siddharth Mahesh; sm0193 **at** mix **dot* wvu **dot* edu
"""

from typing import Any, Dict, Union, cast

import sympy as sp

try:
    from sympy.printing.numpy import JaxPrinter as Printer
    from sympy.printing.numpy import _jax_known_constants as known_constants
    from sympy.printing.numpy import _jax_known_functions as known_functions

    # SymPy's default JaxPrinter uses 'jax.numpy', but we want 'jnp'
    known_functions = {
        k: v.replace("jax.numpy", "jnp") for k, v in known_functions.items()
    }
    known_constants = {
        k: v.replace("jax.numpy", "jnp") for k, v in known_constants.items()
    }
except ImportError:
    # Fallback for older SymPy versions
    from sympy.printing.numpy import NumPyPrinter as Printer  # type: ignore
    from sympy.printing.numpy import _known_constants_numpy, _known_functions_numpy

    known_functions = {k: "jnp." + v for k, v in _known_functions_numpy.items()}
    known_constants = {k: "jnp." + v for k, v in _known_constants_numpy.items()}


# adding a type ignore as mypy does not let me inherit from JaxPrinter.
# Disable specific pylint errors owing to sympy
# pylint: disable=too-many-ancestors, abstract-method
class NRPyJaxPrinter(Printer):
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
        base, exp = cast(sp.Expr, expr).as_base_exp()
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

    def _print_ArrayElementwiseApplyFunc(self, expr: sp.Basic) -> str:
        """
        Print a SymPy ArrayElementwiseApplyFunc expression by inlining the lambda body.

        This lowers elementwise array-application nodes into JAX-broadcastable scalar
        expressions over arrays (e.g., `jnp.sin(A)`, `jnp.abs(A)`, etc.), avoiding
        unsupported SymPy printer nodes during code generation.

        :param expr: ArrayElementwiseApplyFunc expression to print.
        :return: String representation of the elementwise-applied expression.
        """
        # mypy note:
        # SymPy printer hooks accept `Basic`, but `Basic` doesn't declare `.function` / `.expr`.
        # These are runtime attributes on this node, so we cast to Any for safe inspection.
        e = cast(Any, expr)

        # ArrayElementwiseApplyFunc commonly provides `.function` and `.expr`. If not,
        # fall back to the conventional `(function, element)` structure in `expr.args`.
        if hasattr(e, "function") and hasattr(e, "expr"):
            func = e.function  # scalar function to apply (often a sympy.Lambda)
            arr = e.expr  # array operand/expression
        else:
            func, arr = expr.args  # expected: (function, element)

        # Print the array operand once. JAX elementwise ops generally broadcast over arrays.
        arr_str = self._print(arr)

        # Most commonly, SymPy wraps this as a unary Lambda(var, body).
        if isinstance(func, sp.Lambda):
            # If a multi-argument lambda appears, fall back to vectorization so codegen proceeds.
            if len(func.variables) != 1:
                return f"jnp.vectorize({self._print(func)})({arr_str})"

            var = func.variables[0]
            body = func.expr

            # Use a placeholder + string replacement to avoid SymPy rewriting back into the
            # original ArrayElementwiseApplyFunc (or similar array-expression nodes).
            placeholder = sp.Symbol("__AEAF_PH__")
            scalar_body = body.xreplace({var: placeholder})

            body_str = self._print(scalar_body)
            ph_str = self._print(placeholder)

            # Parenthesize the injected array to preserve operator precedence.
            return body_str.replace(ph_str, f"({arr_str})")

        # Non-Lambda case: print as a callable applied to the array.
        # Parentheses ensure correct precedence if `func` prints as an expression.
        return f"({self._print(func)})({arr_str})"
