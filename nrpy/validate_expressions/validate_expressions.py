"""
Assert SymPy Expression Equality.

Authors: Ken Sible, Zachariah Etienne
Emails:  ksible *at* outlook *dot* com
         zachetie *at** gmail *dot** com
"""

from pathlib import Path
import random
from typing import Dict, Union, List, Tuple, Any, Set
import hashlib

import importlib
import black
import sympy as sp
from mpmath import mp, mpf, mpc, fabs  # type: ignore

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


def assert_equal(
    vardict_1: Union[
        Dict[str, Union[sp.Expr, List[sp.Expr], List[List[sp.Expr]]]], sp.Expr, str
    ],
    vardict_2: Union[
        Dict[str, Union[sp.Expr, List[sp.Expr], List[List[sp.Expr]]]], sp.Expr, int
    ],
    suppress_message: bool = False,
) -> None:
    """
    Assert the equality of SymPy expressions or dictionaries containing SymPy expressions.

    :param vardict_1: A SymPy expression or dictionary of SymPy expressions.
    :param vardict_2: A SymPy expression or dictionary of SymPy expressions to compare with vardict_1.
    :param suppress_message: If False, prints a success message when the assertion passes.

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
    >>> assert_equal('(a^2 - b^2) - (a + b)*(a - b)', 0)
    Assertion Passed!
    """
    if not isinstance(vardict_1, dict):
        vardict_1 = {"": sp.sympify(vardict_1)}
    if not isinstance(vardict_2, dict):
        vardict_2 = {"": sp.sympify(vardict_2)}

    vardict_1_results_dict = process_dictionary_of_expressions(
        vardict_1, fixed_mpfs_for_free_symbols=True
    )
    vardict_2_results_dict = process_dictionary_of_expressions(
        vardict_2, fixed_mpfs_for_free_symbols=True
    )
    for var_1, var_2 in zip(vardict_1_results_dict, vardict_2_results_dict):
        tolerance = 10 ** (-4.0 / 5.0 * precision)
        abs_diff = fabs(vardict_1_results_dict[var_1] - vardict_2_results_dict[var_2])
        avg = (vardict_1_results_dict[var_1] + vardict_2_results_dict[var_2]) / 2
        if abs_diff > fabs(avg) * tolerance:
            print(f"Error on {var_1}: {abs_diff} > {fabs(avg)*tolerance}")
            raise AssertionError

    if not suppress_message:
        print("Assertion Passed!")


def inject_mpfs_into_cse_expression(
    free_symbols_dict: Dict[sp.Symbol, Any],
    replaced: List[Tuple[mpf, mpf]],
    reduced: List[sp.Expr],
) -> Union[mpf, mpc]:
    """
    Calculate a numerical value for a given expression using the values in free_symbols_dict.

    :param free_symbols_dict: Dictionary of free symbols and their corresponding values.
    :param replaced: List of tuples containing simplified expressions from SymPy's cse algorithm.
    :param reduced: List of reduced expressions from SymPy's cse algorithm.
    :return: Numerical value of the expression as either mpf or mpc.
    """
    # Original logic, assumes reduced list only has one element
    reduced_expr = reduced[0]

    for lhs, rhs in replaced:
        free_symbols_dict[lhs] = rhs.xreplace(free_symbols_dict)

    reduced_expr = reduced_expr.xreplace(free_symbols_dict)

    # nrpyAbs = sp.Function("nrpyAbs")
    # reduced_expr = reduced_expr.subs(nrpyAbs, sp.Abs)

    try:
        res = mpf(reduced_expr)
    except TypeError:
        if reduced_expr == sp.nan:
            res = mp.nan
            print(
                "inject_mpfs_into_cse_expression warning: after making replacements, found NaN.\n"
                "   Seems to happen in SymTP Jacobians: rfm.Jac_dUrfm_dDSphUD[i][0]"
            )
        else:
            res = mpc(sp.N(reduced_expr, mp.dps))
    return res


def convert_free_symbols_set_to_mpf_dict(
    free_symbols_set: Set[sp.Symbol],
    fixed_mpfs_for_free_symbols: bool = False,
) -> Dict[sp.Symbol, mpf]:
    """
    Convert a set of free symbols into a dictionary with mpf values.

    :param free_symbols_set: Set of free symbols to be converted.
    :param fixed_mpfs_for_free_symbols: Flag to indicate if the mpf values should be fixed.

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
        elif str(var) == "M_SQRT1_2":
            free_symbols_dict[var] = mp.mpf(1.0 / mp.sqrt(2.0))
        # All other free variables are set to random numbers
        else:
            if fixed_mpfs_for_free_symbols:
                random.seed(int(hashlib.md5(str(var).encode()).hexdigest(), 16))
            free_symbols_dict[var] = mp.mpf(random.random())
    return free_symbols_dict


def convert_one_expression_to_mpfmpc(
    expr: sp.Expr,
    fixed_mpfs_for_free_symbols: bool = False,
    verbose: bool = True,
) -> Union[mpf, mpc]:
    """
    Convert a given SymPy expression to either mpf or mpc form.

    :param expr: SymPy expression to convert
    :param fixed_mpfs_for_free_symbols: Whether to fix mpf values for free symbols
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
        free_symbols_set, fixed_mpfs_for_free_symbols=fixed_mpfs_for_free_symbols
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
) -> Dict[str, Union[mpf, mpc]]:
    """
    Process a dictionary of symbolic expressions and convert them to a dictionary of numerical expressions.

    :param dictionary: The input dictionary with symbolic expressions.
    :param fixed_mpfs_for_free_symbols: Flag to indicate if the mpf values should be fixed.
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
                        results_dict[
                            f"{lhs}_{element}"
                        ] = convert_one_expression_to_mpfmpc(
                            item,
                            fixed_mpfs_for_free_symbols=fixed_mpfs_for_free_symbols,
                            verbose=verbose,
                        )
    return results_dict


def check_zero(
    expression: sp.Expr,
    fixed_mpfs_for_free_symbols: bool = False,
    verbose: bool = False,
) -> bool:
    """
    Check if a given expression evaluates to zero.

    :param expression: The SymPy expression to check.
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
        verbose=verbose,
    )

    return bool(result_value == mp.mpf("0.0"))  # Explicitly returning a boolean


def output_trusted(
    os_path_abspath: str,
    os_getcwd: str,
    trusted_file_basename: str,
    results_dict: Dict[Any, Union[mpf, mpc]],
) -> None:
    """
    Write a trusted dictionary to a file with appropriate import headers.

    :param os_path_abspath: The absolute path to the current module. Use os.path.abspath(__file__)
    :param os_getcwd: The current working directory. Use os.getcwd()
    :param trusted_file_basename: The base name for the trusted file.
    :param results_dict: The dictionary containing the results to output as trusted.
    """
    output_str = f"trusted_dict = {results_dict}\n"

    # Determine necessary import headers based on the presence of mpf and mpc in the output string
    if "mpc(" in output_str and "mpf(" in output_str:
        header = "from mpmath import mpf, mpc  # type: ignore\n"
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
    results_dict: Dict[Any, Union[mpf, mpc]],
) -> None:
    """
    Compare the results dictionary against a trusted dictionary.

    :param os_path_abspath: The absolute path to the current module. Use os.path.abspath(__file__)
    :param os_getcwd: The current working directory. Use os.getcwd()
    :param trusted_file_basename: The base name for the trusted file.
    :param results_dict: The dictionary containing the results to compare against trusted.
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
        # We compare relative error here. Pass if abs_diff <= |trusted value| * tolerance
        tolerance = 10 ** (-4.0 / 5.0 * mp.dps)
        abs_diff = results_dict[key] - mpc(trusted_results_dict[key])
        if fabs(abs_diff) > fabs(mpc(trusted_results_dict[key])) * tolerance:
            raise ValueError(
                "\n*** FAILURE ***\n"
                f"****{key}**** mismatch in {trusted_file_basename}! {fabs(abs_diff)} >= {fabs(mpc(trusted_results_dict[key]))} * {tolerance} disagreement\n"
                f"{results_dict[key]} <- computed result,\n{trusted_results_dict[key]} <- trusted value in file {trusted_file_relpath}/tests/{trusted_file_basename}.py\n"
                f"If you trust the new version, then delete {trusted_file_relpath}/tests/{trusted_file_basename}.py and rerun to generate a new version.\n"
                "BE SURE TO INDICATE THE REASON FOR UPDATING THE TRUSTED FILE IN THE COMMIT MESSAGE.\n"
                "***************\n"
            )

    print(f"*** PASS *** {trusted_file_basename} expressions successfully validated.")


def compare_or_generate_trusted_results(
    os_path_abspath: str,
    os_getcwd: str,
    trusted_file_basename: str,
    results_dict: Dict[Any, Union[mpf, mpc]],
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
