"""
Assert SymPy Expression Equality.

Authors: Ken Sible, Zachariah Etienne
Emails:  ksible *at* outlook *dot* com
         zachetie *at** gmail *dot** com
"""

import copy
import hashlib
import importlib
import random
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union, cast

import black
import sympy as sp
from mpmath import fabs, mp, mpc, mpf  # type: ignore

# Contains the constants to be shared throughout all unittests.
# Typical value for precision is 30 significant digits.
precision = 30


def flatten_tensor(tensor: List[Any]) -> List[Any]:
    """
    Recursively flatten a tensor (list of lists of lists of...) to a simple list.

    :param tensor: The tensor to flatten
    :return: A flattened list
    """
    flat_list = []
    for element in tensor:
        if isinstance(element, list):
            flat_list.extend(flatten_tensor(element))
        else:
            flat_list.append(element)
    return flat_list


def _expression_structure(expression: object) -> Tuple[object, ...]:
    """
    Describe supported expression/list nesting without including values.

    :param expression: SymPy expression or recursively nested list to inspect.
    :return: A tuple that preserves every list boundary and expression leaf.
    :raises TypeError: If a leaf is not a SymPy expression.
    """
    if isinstance(expression, sp.Expr):
        return ("expression",)
    if isinstance(expression, list):
        return (
            "list",
            tuple(_expression_structure(element) for element in expression),
        )
    raise TypeError(
        "assert_equal dictionary values must be SymPy expressions or nested lists "
        f"of SymPy expressions, not {type(expression).__name__}"
    )


def _nonfinite_values_match(
    value_1: Union[mpc, mpf], value_2: Union[mpc, mpf]
) -> Optional[bool]:
    """
    Compare non-finite components explicitly.

    :param value_1: First numerical value.
    :param value_2: Second numerical value.
    :return: None when both values are finite, otherwise whether all components match.

    Doctests:
    >>> _nonfinite_values_match(mpf("1.0"), mpf("2.0")) is None
    True
    >>> _nonfinite_values_match(mpc(mp.nan, 1), mpc(mp.nan, 1))
    True
    >>> _nonfinite_values_match(mpc(mp.nan, 1), mpc(mp.nan, 2))
    False
    >>> _nonfinite_values_match(mpc(mp.nan, 0), mpc(0, mp.nan))
    False
    >>> _nonfinite_values_match(mpf("+inf"), mpf("+inf"))
    True
    >>> _nonfinite_values_match(mpf("+inf"), mpf("-inf"))
    False
    """
    # Normalize real and complex inputs so component placement remains part of
    # the sentinel contract.
    value_1_mpc = mpc(value_1)
    value_2_mpc = mpc(value_2)
    # None tells callers to retain their existing finite tolerance comparison.
    if mp.isfinite(value_1_mpc) and mp.isfinite(value_2_mpc):
        return None

    # Compare real and imaginary parts independently: NaN is equal only to NaN
    # in the same component, while infinity must also keep its sign.
    for component_1, component_2 in (
        (value_1_mpc.real, value_2_mpc.real),
        (value_1_mpc.imag, value_2_mpc.imag),
    ):
        if mp.isnan(component_1) or mp.isnan(component_2):
            if not (mp.isnan(component_1) and mp.isnan(component_2)):
                return False
        elif mp.isinf(component_1) or mp.isinf(component_2):
            if component_1 != component_2:
                return False
        elif component_1 != component_2:
            return False
    return True


# Dictionary values use Any to match heterogeneous expression dictionaries
# throughout NRPy and avoid mutable-list invariance for typed tensors. Runtime
# structure validation below still enforces recursively nested SymPy leaves.
def assert_equal(
    vardict_1: Union[Dict[str, Any], sp.Expr, str],
    vardict_2: Union[Dict[str, Any], sp.Expr, int],
    suppress_message: bool = False,
) -> None:
    """
    Assert the equality of SymPy expressions or dictionaries containing SymPy expressions.

    :param vardict_1: A scalar expression input or dictionary whose values are SymPy expressions or recursively nested lists of them.
    :param vardict_2: A scalar expression input or dictionary of supported expression values to compare with vardict_1.
    :param suppress_message: If False, prints a success message when the assertion passes.

    Dictionary inputs must have identical raw key sets and list nesting. Names
    generated for numerically processed leaves must not collide. Matching NaN
    components and same-signed infinities compare equal; one-sided or differently
    structured non-finite values do not. Keys containing ``funcform`` undergo
    structural and leaf validation but are intentionally omitted from numerical
    comparison and processed-name collision checks.

    :raises AssertionError: If dictionary structure or numerical values differ.
    :raises TypeError: If a dictionary leaf is not a SymPy expression.

    Darglint cannot infer the ``TypeError`` propagated by recursive leaf validation.
    # noqa: DAR402

    Doctests:
    >>> from sympy import sin, cos
    >>> from sympy.abc import x

    >>> assert_equal(sin(2*x), 2*sin(x)*cos(x))
    Assertion Passed!
    >>> assert_equal(cos(2*x), cos(x)**2 - sin(x)**2)
    Assertion Passed!
    >>> assert_equal(cos(2*x), 1 - 2*sin(x)**2)
    Assertion Passed!
    >>> assert_equal(cos(2*x), 1 + 2*sin(x)**2)
    Traceback (most recent call last):
    ...
    AssertionError
    >>> vardict_1 = {'A': sin(2*x), 'B': cos(2*x)}
    >>> vardict_2 = {'A': 2*sin(x)*cos(x), 'B': cos(x)**2 - sin(x)**2}
    >>> assert_equal(vardict_1, vardict_2)
    Assertion Passed!
    >>> vardict_1b = {'C': sin(2*x), 'D': cos(2*x)}
    >>> vardict_2b = {'C': 2*sin(x)*cos(x), 'D': cos(x)**2 - sin(x)**2 + 1}  # The +1 messes up the equality.
    >>> assert_equal(vardict_1b, vardict_2b)
    Traceback (most recent call last):
    ...
    AssertionError
    >>> assert_equal({'A': sp.Integer(1)}, {'A': sp.Integer(1), 'B': sp.Integer(2)}, suppress_message=True)
    Traceback (most recent call last):
    ...
    AssertionError
    >>> assert_equal({'A': sp.Integer(1)}, {'B': sp.Integer(1)}, suppress_message=True)
    Traceback (most recent call last):
    ...
    AssertionError
    >>> assert_equal({'A': []}, {'B': []}, suppress_message=True)
    Traceback (most recent call last):
    ...
    AssertionError
    >>> assert_equal({'A': [[x, x]]}, {'A': [x, x]}, suppress_message=True)
    Traceback (most recent call last):
    ...
    AssertionError
    >>> assert_equal({'A': [[[x]]]}, {'A': [[[x]]]}, suppress_message=True)
    >>> assert_equal({'A': [x], 'A_0': x}, {'A': [x], 'A_0': x}, suppress_message=True)
    Traceback (most recent call last):
    ...
    AssertionError
    >>> assert_equal({'A_funcform': sp.Integer(1)}, {'A_funcform': sp.Integer(2)}, suppress_message=True)
    >>> assert_equal({'A_funcform': [sp.Integer(1)], 'A_funcform_0': sp.Integer(2)}, {'A_funcform': [sp.Integer(9)], 'A_funcform_0': sp.Integer(8)}, suppress_message=True)
    >>> assert_equal({'A_funcform': 1}, {'A_funcform': 1}, suppress_message=True)
    Traceback (most recent call last):
    ...
    TypeError: assert_equal dictionary values must be SymPy expressions or nested lists of SymPy expressions, not int
    >>> import contextlib
    >>> import io
    >>> rejected_nonfinite_mismatch = False
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     try:
    ...         assert_equal(sp.nan, sp.Integer(1), suppress_message=True)
    ...     except AssertionError:
    ...         rejected_nonfinite_mismatch = True
    >>> rejected_nonfinite_mismatch
    True
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     assert_equal(sp.nan, sp.nan, suppress_message=True)
    >>> rejected_infinity_mismatch = False
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     try:
    ...         assert_equal(sp.oo, -sp.oo, suppress_message=True)
    ...     except AssertionError:
    ...         rejected_infinity_mismatch = True
    >>> rejected_infinity_mismatch
    True
    >>> assert_equal(sp.oo, sp.oo, suppress_message=True)
    >>> assert_equal(sp.oo + sp.I, sp.oo + sp.I, suppress_message=True)
    >>> assert_equal(sp.zoo, sp.zoo, suppress_message=True)
    >>> assert_equal('(a^2 - b^2) - (a + b)*(a - b)', 0)
    Assertion Passed!
    """
    # Scalar inputs share the dictionary comparison path after sympification.
    if not isinstance(vardict_1, dict):
        vardict_1 = {"": cast(sp.Expr, sp.sympify(vardict_1))}
    if not isinstance(vardict_2, dict):
        vardict_2 = {"": sp.sympify(vardict_2)}

    if vardict_1.keys() != vardict_2.keys():
        raise AssertionError

    # Compare raw nesting before flattening. Otherwise empty tensors disappear,
    # and differently shaped tensors can produce the same numbered leaf names.
    for key in vardict_1:
        if _expression_structure(vardict_1[key]) != _expression_structure(
            vardict_2[key]
        ):
            raise AssertionError

    # Among numerically processed entries, a scalar key such as A_0 can collide
    # with the generated name for A[0]. Reject such dictionaries instead of
    # allowing a later assignment to hide one compared expression.
    for dictionary in (vardict_1, vardict_2):
        processed_names = []
        for key, expression in dictionary.items():
            if "funcform" in key:
                continue
            if isinstance(expression, sp.Expr):
                processed_names.append(str(key))
            else:
                # Raw structure validation above guarantees every remaining
                # dictionary value is a recursively valid list.
                processed_names.extend(
                    f"{key}_{index}"
                    for index, _ in enumerate(flatten_tensor(expression))
                )
        if len(processed_names) != len(set(processed_names)):
            raise AssertionError

    vardict_1_results_dict = process_dictionary_of_expressions(
        vardict_1, fixed_mpfs_for_free_symbols=True
    )
    vardict_2_results_dict = process_dictionary_of_expressions(
        vardict_2, fixed_mpfs_for_free_symbols=True
    )
    # The raw precheck protects empty tensors and shape. This processed check is
    # a second guard against future changes in filtering or name generation.
    if vardict_1_results_dict.keys() != vardict_2_results_dict.keys():
        raise AssertionError

    for var, value_1 in vardict_1_results_dict.items():
        value_2 = vardict_2_results_dict[var]
        # Resolve non-finite sentinels before subtraction, whose NaN result
        # cannot participate safely in an ordered tolerance predicate.
        nonfinite_values_match = _nonfinite_values_match(value_1, value_2)
        if nonfinite_values_match is False:
            print(
                f"Error on {var}: non-finite values do not match: {value_1} != {value_2}"
            )
            raise AssertionError
        if nonfinite_values_match is True:
            continue

        tolerance = 10 ** (-4.0 / 5.0 * precision)
        abs_diff = fabs(value_1 - value_2)
        avg = (value_1 + value_2) / 2
        if abs_diff > fabs(avg) * tolerance:
            print(f"Error on {var}: {abs_diff} > {fabs(avg) * tolerance}")
            raise AssertionError

    if not suppress_message:
        print("Assertion Passed!")


def inject_mpfs_into_cse_expression(
    free_symbols_dict: Dict[sp.Symbol, Any],
    replaced: List[Tuple[mpf, mpf]],
    reduced: List[sp.Expr],
) -> Union[mpc, mpf]:
    """
    Calculate a numerical value for a given expression using the values in free_symbols_dict.

    :param free_symbols_dict: Dictionary of free symbols and their corresponding values.
    :param replaced: List of tuples containing simplified expressions from SymPy's cse algorithm.
    :param reduced: List of reduced expressions from SymPy's cse algorithm.
    :return: Numerical value of the expression as either mpf or mpc.

    Doctests:
    >>> inject_mpfs_into_cse_expression({}, [], [sp.oo])
    mpf('+inf')
    >>> inject_mpfs_into_cse_expression({}, [], [-sp.oo])
    mpf('-inf')
    """
    # Original logic, assumes reduced list only has one element
    reduced_expr = reduced[0]
    free_symbols_dict_orig = copy.deepcopy(free_symbols_dict)
    for lhs, rhs in replaced:
        # Using .evalf(n=mp.dps) at the end ensures that purely numerical expressions are simplified
        # while maintaining the desired precision. This significantly speeds up the validation algorithm.
        free_symbols_dict[lhs] = rhs.xreplace(free_symbols_dict).evalf(n=mp.dps)

    reduced_expr = reduced_expr.xreplace(free_symbols_dict)

    try:
        res = mpf(reduced_expr)
    except TypeError:
        if reduced_expr == sp.nan:
            res = mp.nan
            print(
                "inject_mpfs_into_cse_expression warning: after making replacements, found NaN.\n"
                "   The caller must validate it against an explicit matching sentinel.\n",
                free_symbols_dict_orig,
                replaced,
                reduced,
            )
            partial_env = {}  # type: ignore
            for lhs, rhs in replaced:
                rhs_numeric = rhs.xreplace(partial_env).xreplace(free_symbols_dict)
                print(lhs, "=", rhs_numeric, rhs_numeric.evalf())  # helpful debug print
                partial_env[lhs] = rhs_numeric.evalf(n=mp.dps)
        # SymPy infinities do not pass through mpf() on every supported stack.
        elif reduced_expr == sp.oo:
            res = mpf("+inf")
        elif reduced_expr == -sp.oo:
            res = mpf("-inf")
        else:
            # Convert real and imaginary components separately. mpc() cannot
            # consume SymPy expressions such as oo + I or zoo directly.
            real_expr, imag_expr = reduced_expr.as_real_imag()

            def sympy_component_to_mpf(component: sp.Expr) -> mpf:
                if component == sp.nan:
                    return mp.nan
                if component == sp.oo:
                    return mpf("+inf")
                if component == -sp.oo:
                    return mpf("-inf")
                return mpf(sp.N(component, mp.dps))

            res = mpc(
                sympy_component_to_mpf(real_expr),
                sympy_component_to_mpf(imag_expr),
            )
    return res


def convert_free_symbols_set_to_mpf_dict(
    free_symbols_set: Set[sp.Symbol],
    fixed_mpfs_for_free_symbols: bool = False,
    hex_offset: int = 0,
) -> Dict[sp.Symbol, mpf]:
    """
    Convert a set of free symbols into a dictionary with mpf values.

    :param free_symbols_set: Set of free symbols to be converted.
    :param fixed_mpfs_for_free_symbols: Flag to indicate if the mpf values should be fixed.
    :param hex_offset: Offset for random number; zero by default.

    :return: Dictionary with free symbols as keys and their corresponding mpf values as values.
    """
    free_symbols_dict: Dict[sp.Symbol, Any] = {}

    # Setting each variable in free_symbols_set to a random number in [0, 1) according to the hashed string
    # representation of each variable.
    for var in free_symbols_set:
        # Make sure PI is set to its correct value, pi
        if str(var) in ("PI", "M_PI"):
            free_symbols_dict[var] = mp.mpf(mp.pi)
        # Then make sure M_SQRT1_2 is set to its correct value, 1/sqrt(2)
        elif str(var) in ("SQRT1_2", "M_SQRT1_2"):
            free_symbols_dict[var] = mp.mpf(1.0 / mp.sqrt(2.0))
        # All other free variables are set to random numbers
        else:
            if fixed_mpfs_for_free_symbols:
                random.seed(
                    int(hashlib.md5(str(var).encode()).hexdigest(), 16 + hex_offset)
                )
            free_symbols_dict[var] = mp.mpf(random.random())
    return free_symbols_dict


def convert_one_expression_to_mpfmpc(
    expr: sp.Expr,
    fixed_mpfs_for_free_symbols: bool = False,
    hex_offset: int = 0,
    verbose: bool = True,
) -> Union[mpc, mpf]:
    """
    Convert a given SymPy expression to either mpf or mpc form.

    :param expr: SymPy expression to convert
    :param fixed_mpfs_for_free_symbols: Whether to fix mpf values for free symbols
    :param hex_offset: Offset for random number; zero by default.
    :param verbose: Whether to print verbose messages
    :return: The expression in either mpf or mpc form

    >>> convert_one_expression_to_mpfmpc(sp.sympify(0))
    mpf('0.0')
    """
    mp.dps = precision
    # Speeds up evaluation of Spherical rfm quantities by about 5%
    if expr == sp.sympify(0):
        return mpf(0.0)

    # Using SymPy's cse algorithm to optimize our value substitution
    replaced, reduced = sp.cse(expr, order="none")

    free_symbols_set = set(reduced[0].free_symbols)
    for subexpr in replaced:
        free_symbols_set.update(subexpr[1].free_symbols)

    mpf_symbols_dict = convert_free_symbols_set_to_mpf_dict(
        free_symbols_set,
        fixed_mpfs_for_free_symbols=fixed_mpfs_for_free_symbols,
        hex_offset=hex_offset,
    )

    # Calculate our result_value
    result_value = inject_mpfs_into_cse_expression(mpf_symbols_dict, replaced, reduced)

    if mp.fabs(result_value) != mp.mpf("0.0") and mp.fabs(result_value) < 10 ** (
        (-4.0 / 5.0) * mp.dps
    ):
        if verbose:
            print(
                f"Found |result| ({mp.fabs(result_value)}) close to zero. "
                f"Checking if indeed it should be zero. dps={mp.dps}"
            )
        mp.dps = 2 * precision
        mpf_symbols_dict = convert_free_symbols_set_to_mpf_dict(
            free_symbols_set,
            fixed_mpfs_for_free_symbols=fixed_mpfs_for_free_symbols,
        )
        new_result_value = inject_mpfs_into_cse_expression(
            mpf_symbols_dict, replaced, reduced
        )
        if mp.fabs(new_result_value) < 10 ** (-(4.0 / 5.0) * mp.dps):
            if verbose:
                print(
                    f"After re-evaluating with twice the digits of precision, |result| dropped to "
                    f"{mp.fabs(new_result_value)}. Setting value to zero. dps={mp.dps}"
                )
            result_value = mp.mpf("0.0")
        else:
            if verbose:
                print(
                    f"After re-evaluating with twice the digits of precision, |result| didn't change: "
                    f"{mp.fabs(new_result_value)}. NOT setting value to zero. dps={mp.dps}"
                )

        mp.dps = precision
    return result_value


def process_dictionary_of_expressions(
    dictionary: Dict[Any, Any],
    fixed_mpfs_for_free_symbols: bool = False,
    verbose: bool = True,
) -> Dict[str, Union[mpc, mpf]]:
    """
    Process a dictionary of symbolic expressions and convert them to a dictionary of numerical expressions.

    :param dictionary: The input dictionary with symbolic expressions.
    :param fixed_mpfs_for_free_symbols: Whether to fix mpf values for free symbols
    :param verbose: Flag to indicate if the function should print verbose output.

    :return: A dictionary with the numerical evaluation of the symbolic expressions.
    """
    results_dict = {}
    for lhs, rhs in sorted(dictionary.items()):
        if "funcform" not in lhs:
            if isinstance(rhs, sp.Expr):
                results_dict[str(lhs)] = convert_one_expression_to_mpfmpc(
                    rhs,
                    fixed_mpfs_for_free_symbols=fixed_mpfs_for_free_symbols,
                    verbose=verbose,
                )

            if isinstance(rhs, list):
                for element, item in enumerate(flatten_tensor(rhs)):
                    if isinstance(item, sp.Expr):
                        results_dict[f"{lhs}_{element}"] = (
                            convert_one_expression_to_mpfmpc(
                                item,
                                fixed_mpfs_for_free_symbols=fixed_mpfs_for_free_symbols,
                                verbose=verbose,
                            )
                        )
    return results_dict


def check_zero(
    expression: sp.Expr,
    fixed_mpfs_for_free_symbols: bool = False,
    hex_offset: int = 0,
    verbose: bool = False,
) -> bool:
    """
    Check if a given expression evaluates to zero.

    :param expression: The SymPy expression to check.
    :param fixed_mpfs_for_free_symbols: Whether to fix mpf values for free symbols
    :param hex_offset: Offset for random number; zero by default.
    :param verbose: Flag for additional output.
    :return: True if the expression evaluates to zero, else False.

    >>> from sympy import sin, cos
    >>> from sympy.abc import a,b,x

    >>> check_zero(sin(2*x) - 2*sin(x)*cos(x))
    True
    >>> check_zero(cos(2*x) - (cos(x)**2 - sin(x)**2))
    True
    >>> check_zero(cos(2*x) - (1 - 2*sin(x)**2))
    True
    >>> check_zero(cos(2*x) - (1 + 2*sin(x)**2))
    False
    >>> check_zero((a**2 - b**2) - (a + b)*(a - b))
    True
    >>> check_zero(sp.sympify(0))
    True
    """
    result_value = convert_one_expression_to_mpfmpc(
        expression,
        fixed_mpfs_for_free_symbols=fixed_mpfs_for_free_symbols,
        hex_offset=hex_offset,
        verbose=verbose,
    )

    return bool(result_value == mp.mpf("0.0"))  # Explicitly returning a boolean


def output_trusted(
    os_path_abspath: str,
    os_getcwd: str,
    trusted_file_basename: str,
    results_dict: Dict[Any, Union[mpc, mpf]],
) -> None:
    """
    Write a trusted dictionary to a file with appropriate import headers.

    :param os_path_abspath: The absolute path to the current module. Use os.path.abspath(__file__)
    :param os_getcwd: The current working directory. Use os.getcwd()
    :param trusted_file_basename: The base name for the trusted file.
    :param results_dict: The dictionary containing the results to output as trusted.

    :raises RuntimeError: If an error occurs during file writing.
    """
    output_str = f"trusted_dict = {results_dict}\n"

    # Determine necessary import headers based on the presence of mpf and mpc in the output string
    if "mpc(" in output_str and "mpf(" in output_str:
        header = "from mpmath import mpc, mpf  # type: ignore\n"
    elif "mpf(" in output_str:
        header = "from mpmath import mpf  # type: ignore\n"
    elif "mpc(" in output_str:
        header = "from mpmath import mpc  # type: ignore\n"
    else:
        header = ""

    # Calculate the relative path to the directory where the trusted file will be stored
    trusted_file_relpath = Path(os_path_abspath).relative_to(os_getcwd)
    outdir = trusted_file_relpath.parent / "tests"

    # Determine the output directory and create it if it doesn't exist
    outdir.mkdir(parents=True, exist_ok=True)

    outfile_path = outdir / f"{trusted_file_basename}.py"

    try:
        with outfile_path.open("w", encoding="utf-8") as f:
            formatted_str = black.format_str(header + output_str, mode=black.FileMode())
            f.write(formatted_str)
    except Exception as e:
        raise RuntimeError(f"Error writing to file: {e}") from e

    print(f"*** Outputted trusted dictionary to {outfile_path} . ***")


def compare_against_trusted(
    os_path_abspath: str,
    os_getcwd: str,
    trusted_file_basename: str,
    results_dict: Dict[Any, Union[mpc, mpf]],
) -> None:
    """
    Compare the results dictionary against a trusted dictionary.

    :param os_path_abspath: The absolute path to the current module. Use os.path.abspath(__file__)
    :param os_getcwd: The current working directory. Use os.getcwd()
    :param trusted_file_basename: The base name for the trusted file.
    :param results_dict: The dictionary containing the results to compare against trusted.

    Matching NaN components and same-signed infinities are accepted as explicit
    sentinels. A non-finite value never matches a finite value or a differently
    structured non-finite value.

    :raises ImportError: If the trusted dictionary module cannot be imported.
    :raises ValueError: If there's a mismatch in the number of tested expressions or their values.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> from nrpy.tests.reference_metric_GeneralRFM_fisheyeN2 import trusted_dict
    >>> module_path = str(Path.cwd() / "nrpy/reference_metric.py")
    >>> cwd = str(Path.cwd())
    >>> matching_results = dict(trusted_dict)
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     compare_against_trusted(module_path, cwd, "reference_metric_GeneralRFM_fisheyeN2", matching_results)
    >>> finite_for_trusted_nan = dict(trusted_dict)
    >>> finite_for_trusted_nan["Cart_to_xx_0"] = mpf("1.0")
    >>> rejected_finite_for_trusted_nan = False
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     try:
    ...         compare_against_trusted(module_path, cwd, "reference_metric_GeneralRFM_fisheyeN2", finite_for_trusted_nan)
    ...     except ValueError:
    ...         rejected_finite_for_trusted_nan = True
    >>> rejected_finite_for_trusted_nan
    True
    >>> nan_for_trusted_finite = dict(trusted_dict)
    >>> nan_for_trusted_finite["Cartx"] = mp.nan
    >>> rejected_nan_for_trusted_finite = False
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     try:
    ...         compare_against_trusted(module_path, cwd, "reference_metric_GeneralRFM_fisheyeN2", nan_for_trusted_finite)
    ...     except ValueError:
    ...         rejected_nan_for_trusted_finite = True
    >>> rejected_nan_for_trusted_finite
    True
    """
    # Calculate the relative path to the directory where the trusted file will be stored
    trusted_file_relpath = Path(os_path_abspath).relative_to(os_getcwd).parent
    trusted_file_relpath_dots = (
        str(trusted_file_relpath).replace("/", ".").replace("\\", ".")
    )

    modulename = f"{trusted_file_relpath_dots}.tests.{trusted_file_basename}".lstrip(
        "."
    )
    try:
        trusted = importlib.import_module(modulename)
    except ImportError as exc:
        raise ImportError(
            f"It seems the trusted output dictionary file has not been generated for: {modulename} (file open failed)"
        ) from exc

    trusted_results_dict = trusted.trusted_dict

    if len(trusted_results_dict) != len(results_dict):
        missing_keys: Set["str"] = set()
        if len(results_dict) > len(trusted_results_dict):
            missing_keys = set(results_dict.keys()) - set(trusted_results_dict.keys())
        if len(results_dict) < len(trusted_results_dict):
            missing_keys = set(trusted_results_dict.keys()) - set(results_dict.keys())
        raise ValueError(
            "\n*** FAILURE: NUMBER OF TESTED EXPRESSIONS DOES NOT MATCH ***\n"
            f"{len(results_dict)} = number of expressions tested in results dictionary.\n"
            f"{len(trusted_results_dict)} = number of expressions tested in trusted results dictionary.\n"
            f"Here are the missing expressions: {missing_keys}\n"
            f"If you trust the new version, then delete {trusted_file_relpath}/tests/{trusted_file_basename}.py and rerun to generate a new version.\n"
            "BE SURE TO INDICATE THE REASON FOR UPDATING THE TRUSTED FILE IN THE COMMIT MESSAGE.\n"
            "***************\n"
        )

    for key in trusted_results_dict:
        computed_value = results_dict[key]
        trusted_value = mpc(trusted_results_dict[key])
        # A matching sentinel is intentional evidence; every one-sided or
        # differently structured non-finite result is a hard mismatch.
        nonfinite_values_match = _nonfinite_values_match(computed_value, trusted_value)
        if nonfinite_values_match is True:
            continue
        if nonfinite_values_match is False:
            disagreement = "non-finite values do not match"
        else:
            # We compare relative error here. Pass if abs_diff <= |trusted value| * tolerance
            abs_diff = computed_value - trusted_value
            tolerance = 10 ** (-4.0 / 5.0 * mp.dps)
            if fabs(abs_diff) <= fabs(trusted_value) * tolerance:
                continue
            disagreement = (
                f"{fabs(abs_diff)} > {fabs(trusted_value)} * {tolerance} disagreement"
            )

        raise ValueError(
            "\n*** FAILURE ***\n"
            f"****{key}**** mismatch in {trusted_file_basename}! {disagreement}\n"
            f"{computed_value} <- computed result,\n{trusted_results_dict[key]} <- trusted value in file {trusted_file_relpath}/tests/{trusted_file_basename}.py\n"
            f"If you trust the new version, then delete {trusted_file_relpath}/tests/{trusted_file_basename}.py and rerun to generate a new version.\n"
            "BE SURE TO INDICATE THE REASON FOR UPDATING THE TRUSTED FILE IN THE COMMIT MESSAGE.\n"
            "***************\n"
        )

    print(f"*** PASS *** {trusted_file_basename} expressions successfully validated.")


def compare_or_generate_trusted_results(
    os_path_abspath: str,
    os_getcwd: str,
    trusted_file_basename: str,
    results_dict: Dict[Any, Union[mpc, mpf]],
) -> None:
    """
    Compare the results dictionary against a trusted dictionary or (if trusted dict does not exist) generate a new one.

    :param os_path_abspath: The absolute path to the current module. Use os.path.abspath(__file__)
    :param os_getcwd: The current working directory. Use os.getcwd()
    :param trusted_file_basename: The base name for the trusted file.
    :param results_dict: The dictionary containing the results to compare against trusted.
    """
    # Calculate the relative path to the directory where the trusted file will be stored
    trusted_file_relpath = Path(os_path_abspath).relative_to(os_getcwd)
    outdir = trusted_file_relpath.parent / "tests"

    # Determine the output directory and create it if it doesn't exist
    if not outdir.is_dir():
        outdir.mkdir(parents=True, exist_ok=True)

    outfile_path = outdir / f"{trusted_file_basename}.py"

    if outfile_path.exists():
        compare_against_trusted(
            os_path_abspath,
            os_getcwd,
            trusted_file_basename,
            results_dict,
        )
    else:
        output_trusted(
            os_path_abspath,
            os_getcwd,
            trusted_file_basename,
            results_dict,
        )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
