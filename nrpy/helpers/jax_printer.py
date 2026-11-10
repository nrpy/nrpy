"""
Custom JAX printer used by py_codegen to print SymPy expressions as JAX code.

The printer performs the same power simplification as
custom_c_codegen_functions.py for double precision values, prints integers and
rationals that do not fit in a signed 32-bit integer as floats, and prints Max
and Min as nested jnp.maximum and jnp.minimum calls.

Author: Siddharth Mahesh
        sm0193 **at** mix **dot* wvu **dot* edu
"""

from typing import TYPE_CHECKING, Any, Dict, Tuple, Union, cast

import sympy as sp

if TYPE_CHECKING:
    # SymPy printer classes may be untyped (Any) on older SymPy installs, which
    # makes mypy complain about subclassing. Provide a small typed shim instead.
    # pylint: disable=unused-argument
    class Printer:  # pragma: no cover
        """Dummy class definition for typehinting."""

        _module: str
        _kf: Dict[Any, Any]
        _kc: Dict[Any, Any]

        def __init__(self, settings: Any = ...) -> None: ...
        def _print(self, expr: Any) -> Any: ...
        def _print_Pow(self, expr: sp.Basic, rational: bool = ...) -> Any: ...
        def _print_Integer(self, expr: sp.Integer) -> str: ...
        def _print_Rational(self, expr: sp.Rational) -> str: ...
        def doprint(
            self,
            expr: Any,
            assign_to: Any = ...,
            **kwargs: Any,
        ) -> str:
            """
             Typing-only stub for SymPy printers.

             :param expr: Expression to print.
             :param assign_to: Optional assignment target (e.g., a symbol/name).
            :param kwargs: Additional printer options.
             :return: Printed representation of `expr`.
            """
            return ""

else:
    try:
        from sympy.printing.numpy import JaxPrinter as Printer
    except ImportError:
        # Fallback for older SymPy versions
        from sympy.printing.numpy import NumPyPrinter as Printer

try:
    # SymPy-dev / newer SymPy
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
    from sympy.printing.numpy import _known_constants_numpy, _known_functions_numpy

    known_functions = {k: "jnp." + v for k, v in _known_functions_numpy.items()}
    known_constants = {k: "jnp." + v for k, v in _known_constants_numpy.items()}


# In operations with arrays, JAX rejects Python integers that do not fit in a
# signed 32-bit integer (signed 64-bit with jax_enable_x64), so integers of
# magnitude INT32_BOUND or larger are printed as floats.
INT32_BOUND = 2**31


# Disable specific pylint errors owing to sympy
# pylint: disable=too-many-ancestors, abstract-method
class NRPyJaxPrinter(Printer):
    """
    Print SymPy expressions as JAX code, with NRPy power simplification, large integers as floats, and Max/Min as jnp calls.

    Doctests:
    >>> x, y = sp.symbols("x y", real=True)
    >>> p = NRPyJaxPrinter()
    >>> p.doprint(3*x + 2**40*y)
    '3*x + 1099511627776.0*y'
    >>> p.doprint(sp.Rational(10**30 + 1, 2)*x)
    '(5e+29)*x'
    >>> p.doprint(sp.Rational(1, 2**32)*x)
    '(2.3283064365386963e-10)*x'
    >>> p.doprint(sp.Rational(1, 3)*x)
    '(1/3)*x'
    >>> p.doprint(sp.Max(0, 1 - 4*x))
    'jnp.maximum(0, 1 - 4*x)'
    >>> p.doprint(sp.Min(x, y, 1))
    'jnp.minimum(1, jnp.minimum(x, y))'
    """

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

    def _print_Integer(self, expr: sp.Integer) -> str:
        """
        Print an integer, as a float if it does not fit in a signed 32-bit integer.
        :param expr: Integer to print.
        :return: String representation of the integer.
        """
        if abs(expr.p) < INT32_BOUND:
            return super()._print_Integer(expr)
        return repr(float(expr.p))

    def _print_Rational(self, expr: sp.Rational) -> str:
        """
        Print a rational number, as a float if its numerator or denominator does not fit in a signed 32-bit integer.
        :param expr: Rational number to print.
        :return: String representation of the rational number.
        """
        if abs(expr.p) < INT32_BOUND and expr.q < INT32_BOUND:
            return super()._print_Rational(expr)
        return repr(float(expr))

    def _print_nested(self, func: str, args: Tuple[sp.Basic, ...]) -> str:
        """
        Print a function of several arguments as nested calls of a two-argument JAX function.
        :param func: Name of the two-argument function in the jnp module, e.g., "maximum".
        :param args: Arguments of the function.
        :return: String representation of the nested calls.
        """
        result = self._print(args[-1])
        for arg in reversed(args[:-1]):
            result = f"{self._module}.{func}({self._print(arg)}, {result})"
        return str(result)

    def _print_Max(self, expr: sp.Basic) -> str:
        """
        Print Max as nested jnp.maximum calls, which need no import of functools.
        :param expr: Max expression to print.
        :return: String representation of the Max expression.
        """
        return self._print_nested("maximum", expr.args)

    def _print_Min(self, expr: sp.Basic) -> str:
        """
        Print Min as nested jnp.minimum calls, which need no import of functools.
        :param expr: Min expression to print.
        :return: String representation of the Min expression.
        """
        return self._print_nested("minimum", expr.args)

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

            # Ensure strict string typing for mypy (Base _print returns Any)
            body_str = str(self._print(scalar_body))
            ph_str = str(self._print(placeholder))

            # Parenthesize the injected array to preserve operator precedence.
            return body_str.replace(ph_str, f"({arr_str})")

        # Non-Lambda case: print as a callable applied to the array.
        # Parentheses ensure correct precedence if `func` prints as an expression.
        return f"({self._print(func)})({arr_str})"


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
