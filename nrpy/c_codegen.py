"""
This core NRPy+ module is used for
generating C code kernels.

Authors: Zachariah B. Etienne; zachetie **at** gmail **dot* com
         Ken Sible; ksible **at** outlook **dot* com
"""

import logging

import re  # Regular expressions can be toxic due to edge cases -- we use them sparingly
import sys
from typing import List, Union, Dict, Any, Optional, Sequence, Tuple
import sympy as sp
import nrpy.finite_difference as fin
import nrpy.params as par

from nrpy.helpers.simd import expr_convert_to_simd_intrins
from nrpy.helpers.generic import superfast_uniq, clang_format
from nrpy.helpers.cse_preprocess_postprocess import (
    cse_preprocess,
    cse_postprocess,
)  # NRPy+: CSE preprocessing and postprocessing


class CCodeGen:
    """
    Stores and processes input parameters to c_codegen() below
    """

    def __init__(
        self,
        prestring: str = "",
        poststring: str = "",
        include_braces: bool = True,
        c_type: str = "double",
        c_type_alias: str = "",
        verbose: bool = True,
        enable_cse: bool = True,
        cse_sorting: str = "canonical",
        cse_varprefix: str = "",
        enable_cse_preprocess: bool = False,
        enable_simd: bool = False,
        simd_find_more_subs: bool = False,
        simd_find_more_FMAsFMSs: bool = True,
        simd_debug: bool = False,
        enable_GoldenKernels: bool = False,
        SCALAR_TMP_varnames: Optional[List[str]] = None,
        SCALAR_TMP_sympyexprs: Optional[List[str]] = None,
        postproc_substitution_dict: Optional[Dict[str, str]] = None,
        enable_fd_codegen: bool = False,
        automatically_read_gf_data_from_memory: bool = False,
        enforce_c_parameters_must_be_defined: bool = False,
        enable_fd_functions: bool = False,
        mem_alloc_style: str = "210",
        upwind_control_vec: Union[List[sp.Symbol], sp.Symbol] = sp.Symbol("unset"),
        symbol_to_Rational_dict: Optional[Dict[sp.Basic, sp.Rational]] = None,
        clang_format_enable: bool = False,
        clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 200}",
    ) -> None:
        """
        Initializes the CCodeGen class with provided options for generating C code.

        :param prestring: String to be included before the code.
        :param poststring: String to be included after the code.
        :param include_braces: Boolean to decide whether to include braces.
        :param c_type: C data type, such as 'double'.
        :param c_type_alias: Alias for the C data type, such as 'REAL' or 'CCTK_REAL'.
        :param verbose: Boolean for verbosity.
        :param enable_cse: Boolean to enable common subexpression elimination.
        :param cse_sorting: Sorting method for common subexpression elimination.
        :param cse_varprefix: Variable prefix for common subexpression elimination.
        :param enable_cse_preprocess: Boolean for preprocessing common subexpression elimination.
        :param enable_simd: Boolean to enable SIMD.
        :param simd_find_more_subs: Boolean for SIMD optimization.
        :param simd_find_more_FMAsFMSs: Boolean for SIMD optimization.
        :param simd_debug: Boolean for SIMD debugging.
        :param enable_GoldenKernels: Boolean to enable Golden Kernels.
        :param SCALAR_TMP_varnames: List of temporary scalar variable names.
        :param SCALAR_TMP_sympyexprs: List of temporary scalar sympy expressions.
        :param postproc_substitution_dict: Dictionary for postprocessing substitutions.
        :param enable_fd_codegen: Boolean to enable finite difference code generation.
        :param automatically_read_gf_data_from_memory: Boolean for automatic memory read.
        :param enforce_c_parameters_must_be_defined: Boolean to enforce parameter definitions.
        :param enable_fd_functions: Boolean to enable finite difference functions.
        :param mem_alloc_style: Memory allocation style.
        :param upwind_control_vec: Upwind control vector as a symbol or list of symbols.
        :param clang_format_enable: Boolean to enable clang formatting.
        :param clang_format_options: Options for clang formatting.
        """
        self.prestring = prestring
        self.poststring = poststring
        self.include_braces = include_braces
        self.c_type = c_type
        self.c_type_alias = c_type_alias
        self.verbose = verbose
        self.enable_cse = enable_cse
        self.cse_sorting = cse_sorting
        self.cse_varprefix = cse_varprefix
        self.enable_cse_preprocess = enable_cse_preprocess
        self.enable_simd = enable_simd
        self.simd_find_more_subs = simd_find_more_subs
        self.simd_find_more_FMAsFMSs = simd_find_more_FMAsFMSs
        self.simd_debug = simd_debug
        self.enable_GoldenKernels = enable_GoldenKernels
        self.SCALAR_TMP_varnames = (
            SCALAR_TMP_varnames if SCALAR_TMP_varnames is not None else []
        )
        self.SCALAR_TMP_sympyexprs = (
            SCALAR_TMP_sympyexprs if SCALAR_TMP_sympyexprs is not None else []
        )
        self.postproc_substitution_dict = (
            postproc_substitution_dict if postproc_substitution_dict is not None else {}
        )
        self.enable_fd_codegen = enable_fd_codegen
        self.automatically_read_gf_data_from_memory = (
            automatically_read_gf_data_from_memory
        )
        self.enforce_c_parameters_must_be_defined = enforce_c_parameters_must_be_defined
        self.enable_fd_functions = enable_fd_functions
        self.mem_alloc_style = mem_alloc_style
        self.upwind_control_vec = upwind_control_vec
        self.symbol_to_Rational_dict = symbol_to_Rational_dict
        self.clang_format_enable = clang_format_enable
        self.clang_format_options = clang_format_options

        self.fd_order = par.parval_from_str("finite_difference::fd_order")

        # Now, process input!

        # Set c_type and c_type_alias
        if self.c_type not in ("float", "double", "long double"):
            raise ValueError(
                "c_type must be a standard C type for floating point numbers: float, double, or long double,\n"
                "as it is used to find the appropriate transcendental function; e.g., sin(), sinf(), or sinl().\n"
                "c_type_alias will appear in the generated code, and will be set according to infrastructure.\n"
                f"You chose c_type={self.c_type}"
            )
        Infrastructure = par.parval_from_str("Infrastructure")
        if self.enable_simd:
            if self.c_type not in "double":
                raise ValueError(
                    "SIMD output currently only supports double precision. Sorry!"
                )
            # If enable_simd==True, then check if c_type=="double". If not, error out.
            #         Otherwise set c_type="REAL_SIMD_ARRAY", which should be #define'd
            #         within the C code. For example for AVX-256, the C code should have
            #         #define REAL_SIMD_ARRAY __m256d
            if Infrastructure in ("NRPy", "BHaH", "BaseETK"):
                self.c_type = self.c_type_alias = "REAL_SIMD_ARRAY"
            else:
                raise ValueError("FIXME: Please specify the c_type for SIMD")
        else:
            if Infrastructure == "NRPy":
                self.c_type_alias = self.c_type
            elif Infrastructure == "BHaH":
                self.c_type_alias = "REAL"
            elif Infrastructure in ("BaseETK", "CarpetX"):
                self.c_type_alias = "CCTK_REAL"
            else:
                self.c_type_alias = ""

        if self.enable_GoldenKernels:
            self.enable_cse_preprocess = True
            self.simd_find_more_subs = True
            self.simd_find_more_FMAsFMSs = True

        if self.enable_cse_preprocess:
            sympy_version = sp.__version__.replace("rc", "...").replace("b", "...")
            sympy_major_version = int(sympy_version.split(".")[0])
            sympy_minor_version = int(sympy_version.split(".")[1])
            if sympy_major_version < 1 or (
                sympy_major_version == 1 and sympy_minor_version < 3
            ):
                logging.warning(
                    "SymPy version %s does not support CSE preprocessing. Disabling...",
                    sympy_version,
                )
                self.enable_cse_preprocess = False

        if (
            not isinstance(self.upwind_control_vec, sp.Symbol)
            and str(self.upwind_control_vec) == "unset"
        ):
            if not isinstance(self.upwind_control_vec, list) or not isinstance(
                self.upwind_control_vec[0], sp.Symbol
            ):
                raise ValueError(
                    "When specifying upwind_control_vec, it must be a list of sympy expressions"
                )

        if self.fd_order % 2 != 0 or self.fd_order <= 0:
            raise ValueError(
                "Due to the ambiguity of centering a lopsided stencil,"
                "we choose fd_order to be an even, positive integer by convention."
            )

        if self.enable_fd_codegen:
            self.automatically_read_gf_data_from_memory = True

        if self.enable_fd_functions and not self.enable_fd_codegen:
            raise ValueError("enable_fd_functions=True requires enable_fd_codegen=True")


def c_codegen(
    sympyexpr: Union[
        Sequence[Union[sp.Basic, sp.Expr, sp.Symbol]],
        Union[sp.Basic, sp.Expr, sp.Symbol],
    ],
    output_varname_str: Union[List[str], str],
    **kwargs: Any,
) -> str:
    """
    Outputs C code given SymPy expressions and variable names.

    :param sympyexpr: A SymPy expression or list of SymPy expressions to be converted.
    :param output_varname_str: A string or list of strings representing the variable name(s) in the output.
    :param kwargs: Additional keyword arguments for customization.
    :return: A string containing the generated C code.

    >>> x, y, z = sp.symbols("x y z", real=True)
    >>> print(c_codegen(x**2 + sp.sqrt(y) - sp.sin(x*z), "double blah", include_braces=False, verbose=False))
    double blah = ((x)*(x)) + sqrt(y) - sin(x*z);
    <BLANKLINE>
    >>> print(c_codegen(x**5 + x**3 + x - 1/x, "REAL_SIMD_ARRAY blah", include_braces=False, verbose=False, enable_simd=True))
    const double dbl_Integer_1 = 1.0;
    const REAL_SIMD_ARRAY _Integer_1 = ConstSIMD(dbl_Integer_1);
    <BLANKLINE>
    const double dbl_NegativeOne_ = -1.0;
    const REAL_SIMD_ARRAY _NegativeOne_ = ConstSIMD(dbl_NegativeOne_);
    <BLANKLINE>
    REAL_SIMD_ARRAY blah = FusedMulAddSIMD(MulSIMD(x, x), x, FusedMulAddSIMD(MulSIMD(MulSIMD(MulSIMD(x, x), x), x), x, SubSIMD(x, DivSIMD(_Integer_1, x))));
    <BLANKLINE>
    >>> print(c_codegen(x**5 + x**3 + sp.sin(x**3), "REAL_SIMD_ARRAY blah", enable_simd=True))
    /*
     *  Original SymPy expression:
     *  "REAL_SIMD_ARRAY blah = x**5 + x**3 + sin(x**3)"
     */
    {
    const REAL_SIMD_ARRAY tmp0 = MulSIMD(MulSIMD(x, x), x);
    REAL_SIMD_ARRAY blah = AddSIMD(tmp0, FusedMulAddSIMD(MulSIMD(MulSIMD(MulSIMD(x, x), x), x), x, SinSIMD(tmp0)));
    }
    <BLANKLINE>
    >>> print(c_codegen(x**5 + x**3 + sp.sin(x**3), "REAL_SIMD_ARRAY blah", enable_simd=True))
    /*
     *  Original SymPy expression:
     *  "REAL_SIMD_ARRAY blah = x**5 + x**3 + sin(x**3)"
     */
    {
    const REAL_SIMD_ARRAY tmp0 = MulSIMD(MulSIMD(x, x), x);
    REAL_SIMD_ARRAY blah = AddSIMD(tmp0, FusedMulAddSIMD(MulSIMD(MulSIMD(MulSIMD(x, x), x), x), x, SinSIMD(tmp0)));
    }
    <BLANKLINE>
    """
    CCGParams = CCodeGen(**kwargs)

    # Step 1: Initialize
    #  commentblock: comment block containing the input SymPy string,
    #                set only if verbose==True
    #  outstring:    the output C code string
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

    # Step 3: Process code generation for reading/writing gridfunctions and evaluating finite differences
    if CCGParams.automatically_read_gf_data_from_memory or CCGParams.enable_fd_codegen:
        free_symbols_list: List[sp.Basic] = []
        for expr in sympyexpr_list:
            free_symbols_list.extend(expr.free_symbols)
        list_of_deriv_vars = fin.extract_list_of_deriv_var_strings_from_sympyexpr_list(
            free_symbols_list, sp.Symbol("unset")
        )
        (
            list_of_base_gridfunction_names_in_derivs,
            list_of_deriv_operators,
        ) = fin.extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars(
            list_of_deriv_vars
        )

        # Store fdcoeffs and stencil (tuple) for each unique
        #   derivative operator to deriv_operator_dict
        deriv_operator_dict = {}
        for deriv_op in superfast_uniq(list_of_deriv_operators):
            deriv_operator_dict[deriv_op] = fin.compute_fdcoeffs_fdstencl(
                deriv_op, CCGParams.fd_order
            )

        fdcoeffs: List[List[sp.Rational]] = [
            [] for _ in range(len(list_of_deriv_operators))
        ]
        fdstencl: List[List[List[int]]] = [
            [[] for _ in range(4)] for __ in range(len(list_of_deriv_operators))
        ]
        for i, deriv_op in enumerate(list_of_deriv_operators):
            fdcoeffs[i], fdstencl[i] = deriv_operator_dict[deriv_op]

        read_from_memory_C_code = fin.read_gfs_from_memory(
            list_of_base_gridfunction_names_in_derivs,
            fdstencl,
            free_symbols_list,
            mem_alloc_style=CCGParams.mem_alloc_style,
            enable_simd=CCGParams.enable_simd,
        )

        # This calls outputC as needed to construct a C kernel that does gridfunction management with or without FDs
        return gridfunction_management_and_FD_codegen(
            sympyexpr_list,
            output_varname_str,
            list_of_deriv_vars,
            list_of_base_gridfunction_names_in_derivs,
            list_of_deriv_operators,
            deriv_operator_dict,
            read_from_memory_Ccode=read_from_memory_C_code,
            upwind_control_vec=CCGParams.upwind_control_vec,
            enable_fd_functions=CCGParams.enable_fd_functions,
            enable_simd=CCGParams.enable_simd,
        )

    # Step 4: If CCGParams.verbose, then output the original SymPy
    #         expression(s) in code comments prior to actual C code
    if CCGParams.verbose:
        plural = "s" if len(output_varname_str) > 1 else ""
        commentblock += f"/*\n *  Original SymPy expression{plural}:\n"

        expressions = "\n".join(
            f'*  "[{varname} = {expr}]"'
            if len(output_varname_str) > 1
            else f'*  "{varname} = {expr}"'
            for varname, expr in zip(output_varname_str, sympyexpr_list)
        )

        commentblock += f" {expressions}\n */\n"

    # Step 4a: If common subexpression elimination (CSE) disabled, then
    #         just output the SymPy string in the most boring way,
    #         nearly consistent with SymPy's ccode() function,
    #         though with support for float & long double types
    #         as well.
    simd_RATIONAL_decls = RATIONAL_decls = ""
    if not CCGParams.enable_cse:
        # If CSE is disabled:
        for i, expr in enumerate(sympyexpr_list):
            if CCGParams.postproc_substitution_dict:
                expr = apply_substitution_dict(
                    expr, CCGParams.postproc_substitution_dict
                )
            processed_code = sp.ccode(
                expr,
                output_varname_str[i],
                user_functions=custom_functions_for_SymPy_ccode,
            )
            outstring += f"{ccode_postproc(processed_code, CCGParams)}\n"
    # Step 4b: If CSE enabled, then perform CSE using SymPy and then
    #          resulting C code.
    else:
        # If CSE is enabled:
        simd_const_varnms = []
        simd_const_values = []

        varprefix = CCGParams.cse_varprefix
        if CCGParams.enable_cse_preprocess or CCGParams.enable_simd:
            if not CCGParams.symbol_to_Rational_dict:
                # If enable_cse_preprocess == True, then perform partial factorization
                # If enable_simd == True, then declare _NegativeOne_ in preprocessing
                factor_negative = (
                    CCGParams.enable_simd and CCGParams.simd_find_more_subs
                )

                sympyexpr_list, CCGParams.symbol_to_Rational_dict = cse_preprocess(
                    sympyexpr_list,
                    prefix=varprefix,
                    declare_neg1_as_symbol=CCGParams.enable_simd,
                    negative=factor_negative,
                    factor=CCGParams.enable_cse_preprocess,
                )
        if CCGParams.symbol_to_Rational_dict:
            for v in CCGParams.symbol_to_Rational_dict:
                p, q = float(CCGParams.symbol_to_Rational_dict[v].p), float(
                    CCGParams.symbol_to_Rational_dict[v].q
                )
                if not CCGParams.enable_simd:
                    RATIONAL_decls += f"const {CCGParams.c_type_alias} {str(v)} = "

                    # Since Integer is a subclass of Rational in SymPy, we need only check whether
                    # the denominator q = 1 to determine if a rational is an integer.
                    RATIONAL_decls += (
                        f"{str(p)}/{str(q)};\n" if q != 1 else f"{str(p)};\n"
                    )

        #####
        # Prior to the introduction of the SCALAR_TMP type, NRPy+
        # did not insert sympy.Eq objects into the list of functions
        # that were passed to the CSE. SCALAR_TMP places certain
        # ordering requirements on the system which can only be
        # untangled if sympy.Eq objects are used.
        ###

        # This set of symbols includes the SCALAR_TMPs,
        #    but as sp.Eq objects.
        sympyexprs_SCALAR_TMPs_are_sp_Eq: List[Union[sp.Basic, sp.Equality]] = []

        # This set of names does not include the SCALAR_TMPs
        varnames_excluding_SCALAR_TMPs = []

        # Iterate over sympy expressions
        for idx, expr in enumerate(sympyexpr_list):
            var_name = output_varname_str[idx]
            # Use sp.Eq for specified variable names
            if (
                CCGParams.SCALAR_TMP_varnames
                and var_name in CCGParams.SCALAR_TMP_varnames
            ):
                if output_varname_str[idx] != var_name:
                    raise ValueError(
                        f"Error: where set, CCGParams.SCALAR_TMP_varnames[{idx}] must equal output_varname_str[{idx}]"
                    )
                tmp_sympy_expr = CCGParams.SCALAR_TMP_sympyexprs[idx]
                sympyexprs_SCALAR_TMPs_are_sp_Eq.append(sp.Eq(tmp_sympy_expr, expr))
            else:
                varnames_excluding_SCALAR_TMPs.append(var_name)
                sympyexprs_SCALAR_TMPs_are_sp_Eq.append(expr)

        # Check sympy version and process the main group
        sympy_version = tuple(map(int, re.findall(r"\d+", sp.__version__)))
        if sympy_version < (1, 3):
            print(
                f"Warning: SymPy version {sp.__version__} does not support CSE postprocessing."
            )
            cse_results = sp.cse(
                sympyexprs_SCALAR_TMPs_are_sp_Eq,
                sp.numbered_symbols(CCGParams.cse_varprefix + "tmp"),
                order=CCGParams.cse_sorting,
            )
        else:
            cse_results = cse_postprocess(
                sp.cse(
                    sympyexprs_SCALAR_TMPs_are_sp_Eq,
                    sp.numbered_symbols(CCGParams.cse_varprefix + "tmp"),
                    order=CCGParams.cse_sorting,
                )
            )

        # Important: cse_postprocess() has just removed SCALAR_TMPs from the sympyexprs_SCALAR_TMPs_are_sp_Eq.

        # Processing common subexpressions and results from cse_postprocess
        # cse_results[0] contains common subexpression definitions, does not specify varnames.
        for common_subexpression in cse_results[0]:
            full_type_string = f"const {CCGParams.c_type_alias} "
            if CCGParams.enable_simd:
                outstring += (
                    f"{full_type_string}{common_subexpression[0]} = "
                    f"{expr_convert_to_simd_intrins(common_subexpression[1], CCGParams.symbol_to_Rational_dict, varprefix, CCGParams.simd_find_more_FMAsFMSs)};\n"
                )
            else:
                if CCGParams.postproc_substitution_dict:
                    common_subexpression[1] = apply_substitution_dict(
                        common_subexpression[1], CCGParams.postproc_substitution_dict
                    )
                outstring += (
                    full_type_string
                    + ccode_postproc(
                        sp.ccode(
                            common_subexpression[1],
                            common_subexpression[0],
                            user_functions=custom_functions_for_SymPy_ccode,
                        ),
                        CCGParams,
                    )
                    + "\n"
                )

        # cse_results[1] specifies the varnames in terms of CSE variables.
        for i, result in enumerate(cse_results[1]):
            if CCGParams.enable_simd:
                outstring += (
                    f"{varnames_excluding_SCALAR_TMPs[i]} = "
                    f"{expr_convert_to_simd_intrins(result, CCGParams.symbol_to_Rational_dict, varprefix, CCGParams.simd_find_more_FMAsFMSs)};\n"
                )
            else:
                if CCGParams.postproc_substitution_dict:
                    result = apply_substitution_dict(
                        result, CCGParams.postproc_substitution_dict
                    )
                outstring += (
                    ccode_postproc(
                        sp.ccode(
                            result,
                            varnames_excluding_SCALAR_TMPs[i],
                            user_functions=custom_functions_for_SymPy_ccode,
                        ),
                        CCGParams,
                    )
                    + "\n"
                )

        # End of group processing

        # Complication: SIMD functions require numerical constants to be stored in SIMD arrays
        # Resolution: This function extends lists "simd_const_varnms" and "simd_const_values",
        #             which store the name of each constant SIMD array (e.g., _Integer_1) and
        #             the value of each variable (e.g., 1.0).
        if CCGParams.enable_simd and CCGParams.symbol_to_Rational_dict:
            for v in CCGParams.symbol_to_Rational_dict:
                p, q = float(CCGParams.symbol_to_Rational_dict[v].p), float(
                    CCGParams.symbol_to_Rational_dict[v].q
                )
                simd_const_varnms.append(str(v))
                simd_const_values.append(str(p) if q == 1 else f"{p}/{q}")

        # Step 5b.i: If enable_simd is True, and there is at least
        #            one SIMD const variable, then declare the
        #            simd_const_varnms and simd_const_values arrays
        if CCGParams.enable_simd and simd_const_varnms:
            # Sort the lists and remove duplicates
            # MyPy cannot figure this out. Generally it has problems with zip(). Lame.
            simd_const_varnms, simd_const_values = zip(  # type: ignore
                *sorted(set(zip(simd_const_varnms, simd_const_values)))
            )
            if len(simd_const_varnms) != len(simd_const_values):
                raise ValueError(
                    "Error: SIMD constant declaration arrays simd_const_varnms[] and simd_const_values[] have inconsistent sizes!"
                )

            for i, varname in enumerate(simd_const_varnms):
                simd_RATIONAL_decls += (
                    f"const double dbl{varname} = {simd_const_values[i]};\n"
                )
                simd_RATIONAL_decls += (
                    f"const REAL_SIMD_ARRAY {varname} = ConstSIMD(dbl{varname});\n"
                )
                simd_RATIONAL_decls += "\n"

    # Step 6: Construct final output string
    final_Ccode_output_str = commentblock
    if CCGParams.include_braces:
        final_Ccode_output_str += "{\n"
    final_Ccode_output_str += (
        f"{CCGParams.prestring}{RATIONAL_decls}{simd_RATIONAL_decls}"
        f"{outstring}{CCGParams.poststring}"
    )
    if CCGParams.include_braces:
        final_Ccode_output_str += "}\n"

    # Step 7: Return result string
    return final_Ccode_output_str


# Sometimes SymPy has problems evaluating complicated expressions involving absolute
#    values, resulting in hangs. So instead of using sp.Abs(), if we instead use
#    nrpyAbs, we can sidestep the internal SymPy evaluation and force the C
#    codegen to output our desired fabs().
nrpyAbs = sp.Function("nrpyAbs")
custom_functions_for_SymPy_ccode = {
    "nrpyAbs": "fabs",
    "Pow": [
        (lambda b, e: e == 0.5, lambda b, e: f"sqrt({b})"),
        (lambda b, e: e == -0.5, lambda b, e: f"(1.0/sqrt({b}))"),
        (lambda b, e: e == sp.S.One / 3, lambda b, e: f"cbrt({b})"),
        (lambda b, e: e == -sp.S.One / 3, lambda b, e: f"(1.0/cbrt({b}))"),
        (lambda b, e: e == 2, lambda b, e: f"(({b})*({b}))"),
        (lambda b, e: e == 3, lambda b, e: f"(({b})*({b})*({b}))"),
        (lambda b, e: e == 4, lambda b, e: f"(({b})*({b})*({b})*({b}))"),
        (lambda b, e: e == 5, lambda b, e: f"(({b})*({b})*({b})*({b})*({b}))"),
        (lambda b, e: e == -1, lambda b, e: f"(1.0/({b}))"),
        (lambda b, e: e == -2, lambda b, e: f"(1.0/(({b})*({b})))"),
        (lambda b, e: e == -3, lambda b, e: f"(1.0/(({b})*({b})*({b})))"),
        (lambda b, e: e == -4, lambda b, e: f"(1.0/(({b})*({b})*({b})*({b})))"),
        (lambda b, e: e == -5, lambda b, e: f"(1.0/(({b})*({b})*({b})*({b})*({b})))"),
        (lambda b, e: e != -5, "pow"),
    ],
}


def ccode_postproc(string: str, CCGParams: CCodeGen) -> str:
    """
    Processes the generated C code string for functions related to specific data types
    and precision. Appends the appropriate suffix to standard C math library functions
    based on the data c_type (e.g., pow -> powf in single precision), and removes
    the "L" suffix on floating point numbers when not in long double precision.

    :param string: The original C code string.
    :param CCGParams: The CCodeGen object containing c_type information.
    :return: The processed C code string.
    """

    # Define the dictionary to map the c_type to corresponding cmath function suffix
    cmath_suffixes = {
        "float": "f",
        "double": "",
        "long double": "l",
    }

    # If the c_type is not one of the known keys, raise an error
    if CCGParams.c_type not in cmath_suffixes:
        raise ValueError(f"{__name__}::c_type = '{CCGParams.c_type}' not supported")

    # Get the corresponding cmath function suffix from the dictionary
    cmath_suffix = cmath_suffixes[CCGParams.c_type]

    # Append the cmath function suffix to standard C math library functions:
    c_funcs = [
        "pow",
        "sqrt",
        "cbrt",
        "sin",
        "cos",
        "tan",
        "sinh",
        "cosh",
        "tanh",
        "exp",
        "log",
        "fabs",
        "fmin",
        "fmax",
    ]

    # Add "(" to the end of each function name and join them with '|' to create a pattern that matches any of them
    pattern = "|".join([f"{func}\\(" for func in c_funcs])

    # Use a lambda function to add the suffix to the matched function name
    string = re.sub(
        pattern, lambda match: f"{match.group()[:-1]}{cmath_suffix}(", string
    )

    # If c_type is not 'long double', get rid of the "L" suffix on floating point numbers:
    if CCGParams.c_type != "long double":
        string = re.sub(r"([0-9.]+)L/([0-9.]+)L", "(\\1 / \\2)", string)

    return string


def apply_substitution_dict(
    expr: sp.Basic, postproc_substitution_dict: Dict[str, str]
) -> sp.Basic:
    """
    Given a SymPy expression and a substitution dictionary, this function
    applies the substitutions defined in the dictionary to the expression.

    :param expr: A SymPy expression
    :param postproc_substitution_dict: A dictionary with the original variable names
                                       as keys and the substrings to append as values.

    :return: Modified SymPy expression after substitutions.
    :rtype: sp.Basic

    :Example:

    >>> import sympy as sp
    >>> x = sp.Symbol('x')
    >>> apply_substitution_dict(x**2 + 1, {'x': '_new'})
    x_new**2 + 1
    """

    # Loop through all the free symbols in the given SymPy expression.
    for sym in list(expr.free_symbols):
        ss = str(sym)
        # If the current symbol's string representation is in the substitution dictionary...
        if ss in postproc_substitution_dict:
            # ... then substitute the symbol in the expression with the new symbol formed by appending the string from the dictionary.
            expr = expr.subs(sym, sp.symbols(ss + postproc_substitution_dict[ss]))

    # Return the modified expression.
    return expr


def gridfunction_management_and_FD_codegen(
    sympyexpr_list: List[sp.Basic],
    output_varname_str: Union[List[str], str],
    list_of_deriv_vars: List[Union[sp.Symbol, sp.Basic]],
    list_of_base_gridfunction_names_in_derivs: List[str],
    list_of_deriv_operators: List[str],
    deriv_operator_dict: Dict[str, Tuple[List[sp.Rational], List[List[int]]]],
    read_from_memory_Ccode: str,
    **kwargs: Any,
) -> str:
    """
    This function generates C code kernels for reading/writing gridfunctions
    and performing finite-differences with gridfunction data.

    :param List[sp.Basic] sympyexpr_list: List of sympy expressions.
    :param Union[List[str], str] output_varname_str: Output variable name(s) as string or list of strings.
    :param List[Union[sp.Symbol, sp.Basic]] list_of_deriv_vars: List of variables for derivative operations.
    :param List[str] list_of_base_gridfunction_names_in_derivs: List of base grid function names used in derivatives.
    :param List[str] list_of_deriv_operators: List of derivative operators.
    :param str read_from_memory_Ccode: C code string to read from memory.
    :param Any kwargs: Optional additional parameters.
    :return str: A string of code generated based on input parameters.

    >>> import nrpy.indexedexp as ixp
    >>> import nrpy.grid as gri
    >>> par.set_parval_from_str("Infrastructure", "BHaH")
    >>> gri.glb_gridfcs_dict.clear()
    >>> mem_alloc_style = "210"
    >>> enable_simd=False
    >>> c_type="double"
    >>> c_type_alias="REAL"
    >>> enable_fd_functions=False
    >>> fd_order = 2
    >>> upwind_control_vec = gri.register_gridfunctions_for_single_rank1("vetU", group="EVOL")
    >>> vU       = gri.register_gridfunctions_for_single_rank1("vU", group="EVOL")
    >>> hDD      = gri.register_gridfunctions_for_single_rank2("hDD", group="EVOL", symmetry="sym01")
    >>> hDD_dD   = ixp.declarerank3("hDD_dD", symmetry="sym01")
    >>> hDD_dupD = ixp.declarerank3("hDD_dupD", symmetry="sym01")
    >>> a0, a1, b, c = par.register_CodeParameters(c_type_alias=c_type_alias, module=__name__, names=["a0", "a1", "b", "c"], defaultvalues=1)
    >>> lhs_list = ["REAL a0",                       "REAL a1"]
    >>> exprlist = [b*hDD[1][0] + c*hDD_dD[0][1][1], c*hDD_dupD[0][2][2] + b*hDD_dupD[0][2][0] + a1*a0*vU[1]]
    >>> print(exprlist)
    [b*hDD01 + c*hDD_dD011, a0*a1*vU1 + b*hDD_dupD020 + c*hDD_dupD022]
    >>> free_symbols_list = []
    >>> for expr in exprlist:
    ...     free_symbols_list.extend(expr.free_symbols)
    >>> list_of_deriv_vars = fin.extract_list_of_deriv_var_strings_from_sympyexpr_list(free_symbols_list, sp.Symbol("unset"))
    >>> print(list_of_deriv_vars)
    [hDD_dD011, hDD_ddnD020, hDD_ddnD022, hDD_dupD020, hDD_dupD022]
    >>> list_of_base_gridfunction_names_in_derivs, list_of_deriv_operators = fin.extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars(list_of_deriv_vars)
    >>> print(list_of_base_gridfunction_names_in_derivs)
    ['hDD01', 'hDD02', 'hDD02', 'hDD02', 'hDD02']
    >>> fdcoeffs = [[] for _ in list_of_deriv_operators]
    >>> fdstencl = [[[] for _ in range(4)] for __ in list_of_deriv_operators]
    >>> for i, deriv_op in enumerate(list_of_deriv_operators): fdcoeffs[i], fdstencl[i] = fin.compute_fdcoeffs_fdstencl(deriv_op, fd_order)
    >>> deriv_operator_dict = {}
    >>> for deriv_op in superfast_uniq(list_of_deriv_operators):
    ...     deriv_operator_dict[deriv_op] = fin.compute_fdcoeffs_fdstencl(
    ...         deriv_op, fd_order
    ...     )
    >>> memread_Ccode = fin.read_gfs_from_memory(
    ...    list_of_base_gridfunction_names_in_derivs,
    ...    fdstencl,
    ...    free_symbols_list,
    ...    mem_alloc_style=mem_alloc_style,
    ...    enable_simd=enable_simd,
    ... )
    >>> print(gridfunction_management_and_FD_codegen(exprlist, ["a0", "a1"], list_of_deriv_vars,
    ...          list_of_base_gridfunction_names_in_derivs, list_of_deriv_operators, deriv_operator_dict,
    ...          memread_Ccode, upwind_control_vec=upwind_control_vec,
    ...          enable_fd_functions=enable_fd_functions, enable_simd=enable_simd))
    /*
     * NRPy+-Generated GF Access/FD Code, Step 1 of 3:
     * Read gridfunction(s) from main memory and compute FD stencils as needed.
     */
    const REAL hDD01_i1m1 = in_gfs[IDX4(HDD01GF, i0, i1-1, i2)];
    const REAL hDD01 = in_gfs[IDX4(HDD01GF, i0, i1, i2)];
    const REAL hDD01_i1p1 = in_gfs[IDX4(HDD01GF, i0, i1+1, i2)];
    const REAL hDD02_i2m2 = in_gfs[IDX4(HDD02GF, i0, i1, i2-2)];
    const REAL hDD02_i2m1 = in_gfs[IDX4(HDD02GF, i0, i1, i2-1)];
    const REAL hDD02_i0m2 = in_gfs[IDX4(HDD02GF, i0-2, i1, i2)];
    const REAL hDD02_i0m1 = in_gfs[IDX4(HDD02GF, i0-1, i1, i2)];
    const REAL hDD02 = in_gfs[IDX4(HDD02GF, i0, i1, i2)];
    const REAL hDD02_i0p1 = in_gfs[IDX4(HDD02GF, i0+1, i1, i2)];
    const REAL hDD02_i0p2 = in_gfs[IDX4(HDD02GF, i0+2, i1, i2)];
    const REAL hDD02_i2p1 = in_gfs[IDX4(HDD02GF, i0, i1, i2+1)];
    const REAL hDD02_i2p2 = in_gfs[IDX4(HDD02GF, i0, i1, i2+2)];
    const REAL vU1 = in_gfs[IDX4(VU1GF, i0, i1, i2)];
    const REAL FDPart1_Rational_1_2 = 1.0/2.0;
    const REAL FDPart1_Integer_2 = 2.0;
    const REAL FDPart1_Rational_3_2 = 3.0/2.0;
    const REAL FDPart1tmp0 = FDPart1_Rational_3_2*hDD02;
    const REAL UpwindAlgInputhDD_ddnD020 = invdxx0*(-FDPart1_Integer_2*hDD02_i0m1 + FDPart1_Rational_1_2*hDD02_i0m2 + FDPart1tmp0);
    const REAL UpwindAlgInputhDD_ddnD022 = invdxx2*(-FDPart1_Integer_2*hDD02_i2m1 + FDPart1_Rational_1_2*hDD02_i2m2 + FDPart1tmp0);
    const REAL UpwindAlgInputhDD_dupD020 = invdxx0*(FDPart1_Integer_2*hDD02_i0p1 - FDPart1_Rational_1_2*hDD02_i0p2 - FDPart1tmp0);
    const REAL UpwindAlgInputhDD_dupD022 = invdxx2*(FDPart1_Integer_2*hDD02_i2p1 - FDPart1_Rational_1_2*hDD02_i2p2 - FDPart1tmp0);
    const REAL UpwindControlVectorU0 = vetU0;
    const REAL UpwindControlVectorU2 = vetU2;
    const REAL hDD_dD011 = FDPart1_Rational_1_2*invdxx1*(-hDD01_i1m1 + hDD01_i1p1);
    <BLANKLINE>
    /*
     * NRPy+-Generated GF Access/FD Code, Step 2 of 3:
     * Implement upwinding algorithm.
     */
    const REAL Upwind0 = UPWIND_ALG(UpwindControlVectorU0);
    const REAL Upwind2 = UPWIND_ALG(UpwindControlVectorU2);
    const REAL hDD_dupD020 = Upwind0*(-UpwindAlgInputhDD_ddnD020 + UpwindAlgInputhDD_dupD020) + UpwindAlgInputhDD_ddnD020;
    const REAL hDD_dupD022 = Upwind2*(-UpwindAlgInputhDD_ddnD022 + UpwindAlgInputhDD_dupD022) + UpwindAlgInputhDD_ddnD022;
    <BLANKLINE>
    /*
     * NRPy+-Generated GF Access/FD Code, Step 3 of 3:
     * Evaluate SymPy expressions and write to main memory.
     */
    a0 = b*hDD01 + c*hDD_dD011;
    a1 = a0*a1*vU1 + b*hDD_dupD020 + c*hDD_dupD022;
    <BLANKLINE>
    """

    CCGParams = CCodeGen(**kwargs)

    # Step 5.a.i: Read gridfunctions from memory at needed pts.
    # *** No need to do anything here; already set in
    #     string "read_from_memory_Ccode". ***

    # Step 5.a.ii: Perform arithmetic needed for finite differences
    #              associated with input expressions provided in
    #              sympyexpr_list[].rhs.
    #           Note that FDexprs and FDlhsvarnames contain
    #          A) Finite difference expressions (constructed
    #             in steps above) and associated variable names,
    #             and
    #          B) Input expressions sympyexpr_list[], which
    #             in general depend on finite difference
    #             variables.
    FDexprs: List[sp.Basic] = []
    FDlhsvarnames: List[str] = []
    symbol_to_Rational_dict: Dict[sp.Basic, sp.Rational] = {}
    if CCGParams.enable_fd_functions:
        # Compute finite differences using function calls (instead of inlined calculations)?
        # Fdd FD functions to c_function's outC_function_dict (C function dictionary),
        #   AND return the full set of needed calls to these functions (to funccall_list)
        # FIXME: complete this.
        pass
        # funccall_list = fin.add_FD_func_to_outC_function_dict(
        #     list_of_deriv_vars,
        #     list_of_base_gridfunction_names_in_derivs,
        #     list_of_deriv_operators,
        #     fdcoeffs,
        #     fdstencl,
        # )
    else:

        def construct_deriv_prototypes() -> (
            Tuple[Dict[sp.Basic, sp.Rational], List[sp.Basic], List[str]]
        ):
            deriv_var_list = []
            deriv_op_list = []
            fdcoeffs_list = []
            fdstencl_list = []
            for deriv_op, deriv_op_tuple in sorted(deriv_operator_dict.items()):
                deriv_var_list += [sp.Symbol(f"FDPROTO_{deriv_op}")]
                deriv_op_list += [deriv_op]
                fdcoeffs_list += [deriv_op_tuple[0]]
                fdstencl_list += [deriv_op_tuple[1]]

            (
                proto_FDexprs,
                proto_FDlhsvarnames,
            ) = fin.proto_FD_operators_to_sympy_expressions(
                deriv_var_list,
                fdcoeffs_list,
                fdstencl_list,
                enable_simd=CCGParams.enable_simd,
            )

            # Factorize proto_FDexprs and convert Rationals to symbols, store symbol->Rational dictionary.
            proto_FDexprs, symbol_to_Rational_dict = cse_preprocess(
                proto_FDexprs,
                prefix="FDPart1",
                declare_neg1_as_symbol=CCGParams.enable_simd,
                negative=CCGParams.enable_simd,
                factor=True,
            )

            FDexprs = []
            FDlhsvarnames = []
            for i, deriv_var_symbol in enumerate(list_of_deriv_vars):
                # unpack
                operator = list_of_deriv_operators[i]
                proto_idx = deriv_op_list.index(operator)
                gf_name = list_of_base_gridfunction_names_in_derivs[i]

                FDlhsvarnames += [
                    proto_FDlhsvarnames[proto_idx].replace(
                        str(deriv_var_list[proto_idx]), str(deriv_var_symbol)
                    )
                ]

                replace_dict = {}
                for symb in proto_FDexprs[proto_idx].free_symbols:
                    if "FDPROTO" in str(symb):
                        replace_dict[symb] = sp.Symbol(
                            str(symb).replace("FDPROTO", gf_name)
                        )
                FDexprs += [proto_FDexprs[proto_idx].xreplace(replace_dict)]

            return symbol_to_Rational_dict, FDexprs, FDlhsvarnames

        (
            symbol_to_Rational_dict,
            FDexprs,
            FDlhsvarnames,
        ) = construct_deriv_prototypes()

    # Step 5.b.i: (Upwinded derivatives algorithm, part 1):
    # If an upwinding control vector is specified, determine
    #    which of the elements of the vector will be required.
    #    This ensures that those elements are read from memory.
    # For example, if a symmetry axis is specified,
    #     upwind derivatives with respect to only
    #     two of the three dimensions are used. Here
    #     we find all directions used for upwinding.
    for deriv_op in list_of_deriv_operators:
        if "dupD" in deriv_op and not isinstance(CCGParams.upwind_control_vec, list):
            print(
                "Warning: found an upwinded derivative, but upwind_control_vec was not set to a list."
            )
    upwind_directions = []
    if isinstance(CCGParams.upwind_control_vec, list):
        upwind_directions_unsorted_withdups = []
        for deriv_op in list_of_deriv_operators:
            if "dupD" in deriv_op:
                if deriv_op[-1].isdigit():
                    dirn = int(deriv_op[-1])
                    upwind_directions_unsorted_withdups.append(dirn)
                else:
                    raise ValueError(
                        f"Error: Derivative operator {deriv_op} does not contain a valid direction (must be an integer)."
                    )
        if len(upwind_directions_unsorted_withdups) > 0:
            upwind_directions = superfast_uniq(upwind_directions_unsorted_withdups)
            upwind_directions = sorted(upwind_directions, key=sp.default_sort_key)
        #   If upwind control vector is specified,
        #        add upwind control vectors to the
        #        derivative expression list, so its
        #        needed elements are read from memory.
        for dirn in upwind_directions:
            FDexprs += [CCGParams.upwind_control_vec[dirn]]
            FDlhsvarnames += [
                f"const {CCGParams.c_type_alias} UpwindControlVectorU{dirn}"
            ]

    if FDlhsvarnames:
        # Zip the lists, sort by FDlhsvarnames, then unzip
        sorted_pairs = sorted(zip(FDlhsvarnames, FDexprs))
        sorted_varnames, sorted_exprs = zip(*sorted_pairs)

        # Convert the sorted tuples back to lists
        FDlhsvarnames = list(sorted_varnames)
        FDexprs = list(sorted_exprs)

    # Step 5.x: Output useful code comment regarding
    #           which step we are on. *At most* this
    #           is a 3-step process:
    #        1. Read from memory & compute FD stencils,
    #        2. Perform upwinding, and
    #        3. Evaluate remaining expressions+write
    #           results to main memory.
    NRPy_FD_StepNumber = 1
    NRPy_FD__Number_of_Steps = 1
    if len(read_from_memory_Ccode) > 0:
        NRPy_FD__Number_of_Steps += 1
    if not isinstance(CCGParams.upwind_control_vec, str) and len(upwind_directions) > 0:
        NRPy_FD__Number_of_Steps += 1

    Coutput = ""

    # Copy kwargs
    kwargs_FDPart1 = kwargs.copy()
    kwargs_FDPart2 = kwargs.copy()
    if len(read_from_memory_Ccode) > 0:
        Coutput += f"""/*\n * NRPy+-Generated GF Access/FD Code, Step {NRPy_FD_StepNumber} of {NRPy_FD__Number_of_Steps}:
 * Read gridfunction(s) from main memory and compute FD stencils as needed.\n */
"""

        NRPy_FD_StepNumber = NRPy_FD_StepNumber + 1
        # We choose the CSE temporary variable prefix "FDpart1" for the finite difference coefficients:
        kwargs_FDPart1.update(
            {
                "cse_varprefix": "FDPart1",
                "verbose": False,
                "include_braces": False,
                "enable_cse_preprocess": False,
                "simd_find_more_subs": False,
                "automatically_read_gf_data_from_memory": False,
                "enable_fd_codegen": False,
                "symbol_to_Rational_dict": symbol_to_Rational_dict,
            }
        )

        if CCGParams.enable_fd_functions:
            # Compute finite differences using function calls (instead of inlined calculations)
            # FIXME:
            pass
            # Coutput += indent_Ccode(read_from_memory_Ccode, indent=indent)
            # for funccall in funccall_list:
            #     Coutput += indent_Ccode(funccall, indent=indent)
            # if not isinstance(CCGParams.upwind_control_vec, str):
            #     # # Compute finite differences using inlined calculations
            #     Coutput += indent_Ccode(
            #         c_codegen(FDexprs, FDlhsvarnames, **kwargs_FDPart1), indent=indent
            #     )

        else:
            Coutput += read_from_memory_Ccode + c_codegen(
                FDexprs,
                FDlhsvarnames,
                **kwargs_FDPart1,
            )

    # Step 5.b.ii: Implement control-vector upwinding algorithm.
    if not isinstance(CCGParams.upwind_control_vec, str):
        if len(upwind_directions) > 0:
            Coutput += f"""\n/*\n * NRPy+-Generated GF Access/FD Code, Step {NRPy_FD_StepNumber} of {NRPy_FD__Number_of_Steps}:
 * Implement upwinding algorithm.\n */
"""
            NRPy_FD_StepNumber += 1
            if CCGParams.enable_simd:
                for n in ["0", "1"]:
                    Coutput += f"""const double tmp_upwind_Integer_{n} = {n}.000000000000000000000000000000000;\n
const REAL_SIMD_ARRAY upwind_Integer_{n} = ConstSIMD(tmp_upwind_Integer_{n});
"""
            for dirn in upwind_directions:
                Coutput += f"const {CCGParams.c_type_alias} Upwind{dirn} = UPWIND_ALG(UpwindControlVectorU{dirn});\n"

        upwindU = [sp.sympify(0) for _ in range(3)]
        # Populate upwindU with symbolic expressions based on direction
        for direction in upwind_directions:
            upwindU[direction] = sp.sympify(f"Upwind{direction}")

        upwind_expr_list, var_list = [], []

        # Iterate over the list of derivative variables
        for i, deriv_var in enumerate(list_of_deriv_vars):
            operator = list_of_deriv_operators[i]

            # Check if the operator is a 5-length string and contains "dupD"
            if len(operator) == 5 and "dupD" in operator:
                var_dupD = sp.sympify(f"UpwindAlgInput{deriv_var}")
                var_ddnD = sp.sympify(
                    f"UpwindAlgInput{str(deriv_var).replace('_dupD', '_ddnD')}"
                )

                # Extract direction for upwind operation
                upwind_direction = int(operator[-1])

                # Calculate upwind expression
                upwind_expr = (
                    upwindU[upwind_direction] * (var_dupD - var_ddnD) + var_ddnD
                )

                # Update expression and variable lists
                upwind_expr_list.append(upwind_expr)
                var_list.append(f"const {CCGParams.c_type_alias} {str(deriv_var)}")

        # Copy kwargs
        kwargs_FDPart2 = kwargs_FDPart1.copy()
        kwargs_FDPart2["symbol_to_Rational_dict"] = {}
        # We choose the CSE temporary variable prefix "FDpart1" for the finite difference coefficients:
        kwargs_FDPart2.update(
            {
                "enable_cse_preprocess": CCGParams.enable_cse_preprocess,
                "simd_find_more_subs": CCGParams.simd_find_more_subs,
                "cse_varprefix": "FDPart2",
            }
        )
        Coutput += c_codegen(
            upwind_expr_list,
            var_list,
            **kwargs_FDPart2,
        )

    # Step 5.c.i: Add input RHS & LHS expressions from
    #             sympyexpr_list[]
    Coutput += f"""
/*\n * NRPy+-Generated GF Access/FD Code, Step {NRPy_FD_StepNumber} of {NRPy_FD__Number_of_Steps}:
 * Evaluate SymPy expressions and write to main memory.
 */
"""

    exprs = []
    lhsvarnames = []
    for i, expr in enumerate(sympyexpr_list):
        exprs += [expr]
        if CCGParams.enable_simd:
            lhsvarnames += [f"const REAL_SIMD_ARRAY __RHS_exp_{str(i)}"]
        else:
            lhsvarnames += [output_varname_str[i]]

    # Step 5.c.ii: Write output to gridfunctions specified in
    #              sympyexpr_list[].lhs.
    write_to_mem_string = ""
    if CCGParams.enable_simd:
        for i, output_varname in enumerate(output_varname_str):
            write_to_mem_string += f"WriteSIMD(&{output_varname}, __RHS_exp_{i});\n"

    # outputC requires as its second argument a list of strings.
    #   Sometimes when the lhs's are simple constants, but the inputs
    #   contain gridfunctions, it is necessary to convert the lhs's
    #   to strings:
    lhsvarnamestrings = []
    for lhs in lhsvarnames:
        lhsvarnamestrings.append(str(lhs))

    # Copy kwargs
    kwargs_FDPart3 = kwargs_FDPart2.copy()
    # We choose the CSE temporary variable prefix "FDpart2" for the finite difference coefficients:
    kwargs_FDPart3.update(
        {
            "cse_varprefix": "FDPart3",
        }
    )

    codegen = c_codegen(exprs, lhsvarnamestrings, **kwargs_FDPart3)
    Coutput += codegen
    if write_to_mem_string != "":
        Coutput += f"\n{write_to_mem_string}"

    if CCGParams.clang_format_enable:
        Coutput = clang_format(
            Coutput, clang_format_options=CCGParams.clang_format_options
        )

    return Coutput


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
