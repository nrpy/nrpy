# nrpy/py_codegen.py
"""
Core NRPy module used for generating Python code kernels.
This module is a stripped down version of c_codegen.py
for the case of JAX-optimized applications for SEBOB,
which does not use FD-based kernels.

Authors: Zachariah B. Etienne; zachetie **at** gmail **dot* com
         Ken Sible; ksible **at** outlook **dot* com
         Siddharth Mahesh; sm0193 **at** mix **dot* wvu **dot* edu
"""

import sys
from typing import Any, Dict, List, Optional, Sequence, Union

import sympy as sp
from typing_extensions import Literal

import nrpy.params as par
from nrpy.c_codegen import apply_substitution_dict
from nrpy.helpers.cse_preprocess_postprocess import (  # NRPy: CSE preprocessing and postprocessing
    cse_postprocess,
)
from nrpy.helpers.jax_printer import NRPyJaxPrinter
from nrpy.helpers.type_annotation_utilities import (
    generate_class_representation,
    validate_literal_arguments,
)

printer = NRPyJaxPrinter()


class PyCodeGen:
    """Store and process input parameters to c_codegen() below."""

    def __init__(
        self,
        prestring: str = "",
        poststring: str = "",
        verbose: bool = True,
        enable_cse: bool = True,
        cse_sorting: Literal["canonical", "none"] = "canonical",
        cse_varprefix: str = "",
        postproc_substitution_dict: Optional[Dict[str, str]] = None,
    ) -> None:
        """
        Initialize the PyCodeGen class with provided options for generating JAX-compatible Python code.

        :param prestring: String to be included before the code.
        :param poststring: String to be included after the code.
        :param verbose: Boolean to enable verbose output.
        :param enable_cse: Boolean to enable common subexpression elimination.
        :param cse_sorting: Sorting method for common subexpression elimination.
        :param cse_varprefix: Variable prefix for common subexpression elimination.
        :param postproc_substitution_dict: Dictionary for postprocessing substitutions.

        :raises ValueError: If 'Infrastructure' is not 'jax'.
        """
        validate_literal_arguments()
        Infrastructure = par.parval_from_str("Infrastructure")
        if Infrastructure != "JAX":
            raise ValueError("Infrastructure must be 'jax' for py_codegen")
        self.prestring = prestring
        self.poststring = poststring
        self.verbose = verbose
        self.enable_cse = enable_cse
        self.cse_sorting = cse_sorting
        self.cse_varprefix = cse_varprefix
        self.postproc_substitution_dict = (
            postproc_substitution_dict if postproc_substitution_dict is not None else {}
        )

    def __repr__(self) -> str:
        """
        Create a human-readable representation of the PyCodeGen object and what's in it.

        Generates a string representation of the PyCodeGen instance, detailing its configuration
        settings and options in a format that is easy to read and understand. This can be useful
        for debugging or logging purposes, providing quick insights into the instance's state.

        :return: A string that represents the current state and configuration of the PyCodeGen object.
        """
        return generate_class_representation()


def py_codegen(
    sympyexpr: Union[
        Sequence[Union[sp.Basic, sp.Expr, sp.Symbol]],
        Union[sp.Basic, sp.Expr, sp.Symbol],
    ],
    output_varname_str: Union[List[str], str],
    **kwargs: Any,
) -> str:
    """
    Output JAX-compatible Python code given SymPy expressions and variable names.

    :param sympyexpr: A SymPy expression or list of SymPy expressions to be converted.
    :param output_varname_str: A string or list of strings representing the variable name(s) in the output.
    :param kwargs: Additional keyword arguments for customization. They are used to initialize a PyCodeGen object.

    :raises TypeError: If either `sympyexpr` or `output_varname_str` is a tuple, as tuples can cause issues in this function.
    :raises ValueError: If the length of `sympyexpr` and `output_varname_str` do not match, indicating a mismatch between the number of expressions and output variable names.
    :raises ValueError: Under various conditions within the function that are specific to the configuration provided via `kwargs`, such as incompatible options or unsupported configurations.

    :return: A string containing the generated JAX-compatible Python code.

    >>> x, y, z = sp.symbols("x y z", real=True)
    >>> print(py_codegen(x**2 + sp.sqrt(y) - sp.sin(x*z), "blah", verbose=False))
    blah = ((x)*(x)) + jnp.sqrt(y) - jnp.sin(x*z)
    <BLANKLINE>
    >>> print(py_codegen(x**5 + x**3 + sp.sin(x**3), "blah"))
    #
    #  Original SymPy expression:
    #  "blah = x**5 + x**3 + sin(x**3)"
    tmp0 = ((x)*(x)*(x))
    blah = tmp0 + ((x)*(x)*(x)*(x)*(x)) + jnp.sin(tmp0)
    <BLANKLINE>
    """
    # Injected tuples wreak havoc in this function, so check for them & error out if spotted.
    if isinstance(sympyexpr, tuple):
        raise TypeError("sympyexpr cannot be a tuple")
    if isinstance(output_varname_str, tuple):
        raise TypeError("output_varname_str cannot be a tuple")
    PCGParams = PyCodeGen(**kwargs)

    # Step 1: Initialize
    #  commentblock: comment block containing the input SymPy string,
    #                set only if verbose==True
    #  outstring:    the output Python code string
    commentblock = outstring = ""

    # Step 2a: If sympyexpr and output_varname_str are not lists,
    #          convert them to lists of one element each, to
    #          simplify proceeding code.
    output_varname_str = (
        output_varname_str
        if isinstance(output_varname_str, list)
        else [output_varname_str]
    )
    sympyexpr_list = sympyexpr if isinstance(sympyexpr, list) else [sympyexpr]
    sympyexpr_list = sympyexpr_list[
        :
    ]  # Make a copy of sympyexpr_list to safeguard against the input expressions being changed.

    # Step 2b: Check that output_varname_str and sympyexpr_list lists
    #          are the same length
    if len(output_varname_str) != len(sympyexpr_list):
        raise ValueError(
            f"Length of SymPy expressions list ({len(sympyexpr_list)}) != Length of corresponding output variable name list ({len(output_varname_str)})"
        )
    output_varname_str = output_varname_str[
        :
    ]  # Make a copy of output_varname_str to safeguard against the input strings being changed.

    # Step 3: If PCGParams.verbose, then output the original SymPy
    #         expression(s) in code comments prior to actual Python code
    if PCGParams.verbose:
        plural = "s" if len(output_varname_str) > 1 else ""
        commentblock += f"#\n#  Original SymPy expression{plural}:\n"

        expressions = "\n".join(
            (
                f'#  "[{varname} = {expr}]"'
                if len(output_varname_str) > 1
                else f'#  "{varname} = {expr}"'
            )
            for varname, expr in zip(output_varname_str, sympyexpr_list)
        )

        commentblock += f"{expressions}\n"

    # Step 3a: If common subexpression elimination (CSE) disabled, then
    #         just output the SymPy string in the most boring way,
    if not PCGParams.enable_cse:
        # If CSE is disabled:
        for i, expr in enumerate(sympyexpr_list):
            if PCGParams.postproc_substitution_dict:
                expr = apply_substitution_dict(
                    expr, PCGParams.postproc_substitution_dict
                )
            processed_code = printer.doprint(
                expr,
                output_varname_str[i],
            )
            outstring += f"{processed_code}\n"
    # Step 3b: If CSE enabled, then perform CSE using SymPy and then
    #          resulting JAX code.
    else:
        # If CSE is enabled:

        varprefix = PCGParams.cse_varprefix
        varnames: List[str] = []
        sympyexprs: List[sp.Basic] = []

        # Iterate over sympy expressions
        for idx, expr in enumerate(sympyexpr_list):
            var_name = output_varname_str[idx]
            varnames.append(var_name)
            sympyexprs.append(expr)

        # Check sympy version and process the main group
        sp_version = (
            sp.__version__.lower().replace("rc", ".").replace("b", ".")
        )  # turn 1.2rc1 -> 1.2.1, 1.11b1 -> 1.11.1
        parts = sp_version.split(".")
        major = int(parts[0])
        minor = int(parts[1])
        if (major, minor) < (1, 3):
            print(
                f"Warning: SymPy version {sp.__version__} does not support CSE postprocessing."
            )
            cse_results = sp.cse(
                sympyexprs,
                sp.numbered_symbols(varprefix + "tmp"),
                order=PCGParams.cse_sorting,
            )
        else:
            cse_results = cse_postprocess(
                sp.cse(
                    sympyexprs,
                    sp.numbered_symbols(varprefix + "tmp"),
                    order=PCGParams.cse_sorting,
                )
            )

        # Processing common subexpressions and results from cse_postprocess
        # cse_results[0] contains common subexpression definitions, does not specify varnames.
        for common_subexpression in cse_results[0]:
            if PCGParams.postproc_substitution_dict:
                common_subexpression[1] = apply_substitution_dict(
                    common_subexpression[1], PCGParams.postproc_substitution_dict
                )
            outstring += (
                printer.doprint(
                    sp.expand(common_subexpression[1]),
                    common_subexpression[0],
                )
                + "\n"
            )

        # cse_results[1] specifies the varnames in terms of CSE variables.
        for i, result in enumerate(cse_results[1]):
            if PCGParams.postproc_substitution_dict:
                result = apply_substitution_dict(
                    result, PCGParams.postproc_substitution_dict
                )
            outstring += (
                printer.doprint(
                    sp.expand(result),
                    varnames[i],
                )
                + "\n"
            )

        # End of group processing
    # Step 4: Construct final output string
    final_Pycode_output_str = commentblock
    final_Pycode_output_str += (
        f"{PCGParams.prestring}" f"{outstring}{PCGParams.poststring}"
    )

    # Step 5: Return result string
    return final_Pycode_output_str


if __name__ == "__main__":
    par.set_parval_from_str("Infrastructure", "JAX")
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
