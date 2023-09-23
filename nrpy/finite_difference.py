"""
This module provides helper functions for c_codegen
 to generate finite-difference C-code kernels.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import sys  # Standard Python module for multiplatform OS-level functions
from typing import Union, List, Tuple, Any, Dict
from operator import itemgetter
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import nrpy.params as par  # NRPy+: parameter interface
import nrpy.grid as gri  # NRPy+: Functions having to do with numerical grids
from nrpy.helpers.generic import superfast_uniq

par.register_param(py_type=int, module=__name__, name="fd_order", value=4)


#######################################################
#  FINITE-DIFFERENCE COEFFICIENT ALGORITHM


#  Define the to-be-inverted matrix, A.
#  We define A row-by-row, according to the prescription
#  derived in notes/notes.pdf, via the following pattern
#  that applies for arbitrary order.
#
#  As an example, consider a 5-point finite difference
#  stencil (4th-order accurate), where we wish to compute
#  some derivative at the center point.
#
#  Then A is given by:
#
#  -2^0  -1^0  1  1^0   2^0
#  -2^1  -1^1  0  1^1   2^1
#  -2^2  -1^2  0  1^2   2^2
#  -2^3  -1^3  0  1^3   2^3
#  -2^4  -1^4  0  1^4   2^4
#
#  Then right-multiplying A^{-1}
#  by (1 0 0 0 0)^T will yield 0th deriv. stencil
#  by (0 1 0 0 0)^T will yield 1st deriv. stencil
#  by (0 0 1 0 0)^T will yield 2nd deriv. stencil
#  etc.
#
#  Next suppose we want an upwinded, 4th-order accurate
#  stencil. For this case, A is given by:
#
#  -1^0  1  1^0   2^0   3^0
#  -1^1  0  1^1   2^1   3^1
#  -1^2  0  1^2   2^2   3^2
#  -1^3  0  1^3   2^3   3^3
#  -1^4  0  1^4   2^4   3^4
#
#  ... and similarly for the downwinded derivative.
#
#  Finally, let's consider a 3rd-order accurate
#  stencil. This would correspond to an in-place
#  upwind stencil with stencil radius of 2 gridpoints,
#  where other, centered derivatives are 4th-order
#  accurate. For this case, A is given by:
#
#  -1^0  1  1^0   2^0
#  -1^1  0  1^1   2^1
#  -1^2  0  1^2   2^2
#  -1^3  0  1^3   2^3
#  -1^4  0  1^4   2^4
#
#  ... and similarly for the downwinded derivative.
#
#  The general pattern is as follows:
#
#  1) The top row is all 1's,
#  2) If the second row has N elements (N must be odd),
#  .... then the radius of the stencil is rs = (N-1)/2
#  .... and the j'th row e_j = j-rs-1. For example,
#  .... for 4th order, we have rs = 2
#  .... j  | element
#  .... 1  | -2
#  .... 2  | -1
#  .... 3  |  0
#  .... 4  |  1
#  .... 5  |  2
#  3) The L'th row, L>2 will be the same as the second
#  .... row, but with each element e_j -> e_j^(L-1)
#  A1 is used later to validate the inverted
#  matrix.
def setup_FD_matrix__return_inverse_lowlevel(
    stencil_width: int, UPDOWNWIND_stencil_shift: int
) -> Any:
    """
    Function to set up finite difference matrix and return its inverse.

    :param stencil_width: Width of the stencil.
    :param UPDOWNWIND_stencil_shift: Shift in the stencil for upwind or downwind.

    Returns
    -------
    sp.Matrix
        Inverse of the finite difference matrix.

    Example
    -------
    >>> setup_FD_matrix__return_inverse_lowlevel(3, 0)
    Matrix([
    [0, -1/2, 1/2],
    [1,    0,  -1],
    [0,  1/2, 1/2]])
    """
    M = sp.zeros(stencil_width, stencil_width)
    for i in range(stencil_width):
        for j in range(stencil_width):
            if i == 0:
                M[
                    (i, j)
                ] = 1  # Setting n^0 = 1 for all n, including n=0, because this matches the pattern
            else:
                dist_from_xeq0_col = (
                    j - sp.Rational((stencil_width - 1), 2) + UPDOWNWIND_stencil_shift
                )
                if dist_from_xeq0_col == 0:
                    M[(i, j)] = 0
                else:
                    M[(i, j)] = dist_from_xeq0_col**i
    # LU decomposition here speeds up the inversion by about 10-20%, esp for higher order FD computations.
    L, U, _ = M.LUdecomposition()
    L_inv = L.inv()
    U_inv = U.inv()
    # M's inverse is then U_inv * L_inv
    return U_inv * L_inv


class Matrix_dict(Dict[Tuple[int, int], sp.Matrix]):
    """Custom dictionary for storing FD coeff matrices, as the inversion is expensive"""

    def __getitem__(self, key: Tuple[int, int]) -> sp.Matrix:
        if key not in self:
            stencil_width = key[0]
            UPDOWNWIND_stencil_shift = key[1]
            self.__setitem__(
                key,
                setup_FD_matrix__return_inverse_lowlevel(
                    stencil_width=stencil_width,
                    UPDOWNWIND_stencil_shift=UPDOWNWIND_stencil_shift,
                ),
            )
        return dict.__getitem__(self, key)

    def __setitem__(self, key: Tuple[int, int], value: sp.Matrix) -> None:
        dict.__setitem__(self, key, value)


FD_Matrix_dict = Matrix_dict()


def setup_FD_matrix__return_inverse(
    stencil_width: int, UPDOWNWIND_stencil_shift: int
) -> Any:
    """
    Wrapper function to set up finite difference matrix and return its inverse.
    If the inputs have already been processed through
    setup_FD_matrix__return_inverse_lowlevel() (i.e., the matrix has already
    been inverted), this function will simply return the inverse matrix stored within
    FD_Matrix_dict((stencil_width, UPDOWNWIND_stencil_shift). Otherwise, it will
    populate the FD_Matrix_dict with this matrix.

    :param stencil_width: Width of the stencil.
    :param UPDOWNWIND_stencil_shift: Shift in the stencil for upwind or downwind.

    Returns
    -------
    sp.Matrix
        Inverse of the finite difference matrix.

    Example
    -------
    >>> setup_FD_matrix__return_inverse(3, 0)
    Matrix([
    [0, -1/2, 1/2],
    [1,    0,  -1],
    [0,  1/2, 1/2]])
    """
    return FD_Matrix_dict[(stencil_width, UPDOWNWIND_stencil_shift)]


def compute_fdcoeffs_fdstencl(
    derivstring: str, fd_order: int
) -> Tuple[List[sp.Rational], List[List[int]]]:
    """
    Construct finite difference coefficients and stencils for a given derivative type.
    """
    # Step 0: Set finite differencing order, stencil size, and up/downwinding
    if "dKOD" in derivstring:
        fd_order += 2  # par.parval_from_str("FD_KO_ORDER__CENTDERIVS_PLUS")

    stencil_width = fd_order + 1
    UPDOWNWIND_stencil_shift = 0
    # dup/dnD = single-point-offset upwind/downwinding.
    if "dupD" in derivstring:
        UPDOWNWIND_stencil_shift = 1
    elif "ddnD" in derivstring:
        UPDOWNWIND_stencil_shift = -1
    # dfullup/dnD = full upwind/downwinding.
    elif "dfullupD" in derivstring:
        UPDOWNWIND_stencil_shift = int(fd_order / 2)
    elif "dfulldnD" in derivstring:
        UPDOWNWIND_stencil_shift = -int(fd_order / 2)

    # Step 1: Set up FD matrix and return the inverse, as documented above.
    Minv = setup_FD_matrix__return_inverse(stencil_width, UPDOWNWIND_stencil_shift)

    # Step 2:
    #     Based on the input derivative string,
    #     pick out the relevant row of the matrix
    #     inverse, as outlined in the detailed code
    #     comments prior to this function definition.
    derivtype = "FirstDeriv"
    matrixrow = 1
    if "DDD" in derivstring:
        raise ValueError(
            "Error: Only derivatives up to second order currently supported."
        )
    if "DD" in derivstring:
        if derivstring[len(derivstring) - 1] == derivstring[len(derivstring) - 2]:
            # Assuming i==j, we call \partial_i \partial_j gf an "unmixed" second derivative,
            #     or more simply, just "SecondDeriv":
            derivtype = "SecondDeriv"
            matrixrow = 2
        else:
            # Assuming i!=j, we call \partial_i \partial_j gf a MIXED second derivative,
            #     which is computed using a composite of first derivative operations.
            derivtype = "MixedSecondDeriv"
    elif "dKOD" in derivstring:
        derivtype = "KreissOligerDeriv"
        matrixrow = stencil_width - 1
    else:
        # Up/downwinded and first derivs are all of "FirstDeriv" type
        pass

    # Step 3:
    #     Set finite difference coefficients
    #     and stencil points corresponding to
    #     each finite difference coefficient.
    fdcoeffs = []
    fdstencl = []
    if derivtype != "MixedSecondDeriv":
        for i in range(stencil_width):
            idx3 = [0, 0, 0]
            # First compute finite difference coefficient.
            fdcoeff = sp.factorial(matrixrow) * Minv[(i, matrixrow)]
            # Do not store fdcoeff or fdstencil if
            # finite difference coefficient is zero.
            if fdcoeff != 0:
                fdcoeffs.append(fdcoeff)
                if derivtype == "KreissOligerDeriv":
                    fdcoeffs[i] *= (-1) ** (
                        sp.Rational((stencil_width + 1), 2)
                    ) / 2**matrixrow

                # Next store finite difference stencil point
                # corresponding to coefficient.
                gridpt_posn = (
                    i - int((stencil_width - 1) / 2) + UPDOWNWIND_stencil_shift
                )
                if gridpt_posn != 0:
                    dirn = int(derivstring[len(derivstring) - 1])
                    idx3[dirn] = gridpt_posn
                fdstencl.append(idx3)
    else:
        # Mixed second derivative finite difference coeffs
        #     consist of products of first deriv coeffs,
        #     defined in first Minv matrix row.
        for i in range(stencil_width):
            for j in range(stencil_width):
                idx3 = [0, 0, 0]

                # First compute finite difference coefficient.
                fdcoeff = (sp.factorial(matrixrow) * Minv[(i, matrixrow)]) * (
                    sp.factorial(matrixrow) * Minv[(j, matrixrow)]
                )

                # Do not store fdcoeff or fdstencil if
                # finite difference coefficient is zero.
                if fdcoeff != 0:
                    fdcoeffs.append(fdcoeff)

                    # Next store finite difference stencil point
                    # corresponding to coefficient.
                    gridpt_posn1 = i - int((stencil_width - 1) / 2)
                    gridpt_posn2 = j - int((stencil_width - 1) / 2)
                    dirn1 = int(derivstring[len(derivstring) - 1])
                    dirn2 = int(derivstring[len(derivstring) - 2])
                    idx3[dirn1] = gridpt_posn1
                    idx3[dirn2] = gridpt_posn2
                    fdstencl.append(idx3)
    return fdcoeffs, fdstencl


#####################################
# STEP 0: DECLARE FD HELPER FUNCTIONS
def symbol_is_gridfunction_Cparameter_or_other(var: sp.Basic) -> str:
    """
    Determines a given variable's type.

    This function checks whether a given variable is in the list of global grid functions or C parameters,
    and returns a string indicating its type: either 'gridfunction', 'Cparameter', or 'other'.

    Parameters:
    var (sp.Symbol): The variable whose type is to be determined.

    Returns:
    str: A string indicating the type of the variable. Can be 'gridfunction', 'Cparameter', or 'other'.

    Raises:
    ValueError: If the variable is both a gridfunction and a Cparameter, or if it is neither and
    cannot be identified, a ValueError is raised.

    >>> vetU = gri.register_gridfunctions_for_single_rank1("vetU")
    >>> symbol_is_gridfunction_Cparameter_or_other(vetU[0])
    'gridfunction'
    >>> symbol_is_gridfunction_Cparameter_or_other(sp.Symbol('x'))
    'other'
    >>> a, b = par.register_CodeParameters(c_type_alias="double", module="gridtest", names=["a", "b"], defaultvalues=1)
    >>> symbol_is_gridfunction_Cparameter_or_other(a)
    'Cparameter'
    """

    var_is_gf = str(var) in gri.glb_gridfcs_dict
    var_is_Cparameter = any(
        str(var) == param.name for key, param in par.glb_code_params_dict.items()
    )

    if var_is_Cparameter and var_is_gf:
        raise ValueError(
            f"Error: variable {var} is registered both as a gridfunction and as a Cparameter."
        )

    if not var_is_Cparameter and not var_is_gf:
        return "other"
    if var_is_Cparameter:
        return "Cparameter"
    if var_is_gf:
        return "gridfunction"

    raise ValueError("grid.py: Could not find variable_type.")


#########################################
# STEP 1: EXTRACT DERIVATIVES TO COMPUTE
#         FROM LIST OF SYMPY EXPRESSIONS
def extract_list_of_deriv_var_strings_from_sympyexpr_list(
    list_of_free_symbols: List[sp.Basic],
    upwind_control_vec: Union[List[sp.Basic], sp.Basic],
) -> List[sp.Basic]:
    """
    Extract derivative expressions from SymPy expressions' free symbols.

    Args:
    list_of_free_symbols (List[sp.Basic]): List of free symbols from SymPy expressions.
    upwind_control_vec (Union[List[sp.Symbol], str]): Upwind control vector.

    Returns:
    List[sp.Symbol]: List of derivative variables, creating _ddnD in case upwinding is enabled with control vector.

    Raises:
    ValueError: If a variable in the SymPy expression isn't registered as a gridfunction or Cparameter.

    Doctest:
    >>> import nrpy.indexedexp as ixp
    >>> from typing import cast
    >>> gri.glb_gridfcs_dict.clear()
    >>> vU = gri.register_gridfunctions_for_single_rank1("vU")
    >>> vU_dD = cast(List[List[sp.Symbol]], ixp.declarerank2("vU_dD"))
    >>> expr = vU_dD[0][1] * vU[2] + vU_dD[2][0] * vU[1]
    >>> list_of_free_symbols = expr.free_symbols
    >>> upwind_control_vec = "unset"
    >>> extract_list_of_deriv_var_strings_from_sympyexpr_list(list_of_free_symbols, upwind_control_vec)
    [vU_dD01, vU_dD20]
    """
    # Create a list of free symbols that are classified neither as gridfunctions nor
    # as C parameters. These must be derivatives, hence we call the list "list_of_deriv_vars"
    list_of_deriv_vars_with_duplicates = []

    for var in list_of_free_symbols:
        vartype = symbol_is_gridfunction_Cparameter_or_other(var)
        if vartype == "other":
            if any(s in str(var) for s in ["_dD", "_dKOD", "_dupD", "_ddnD"]):
                list_of_deriv_vars_with_duplicates.append(var)
            else:
                pass
                # err_msg = (
                #     f'Error: Unregistered variable "{var}" in SymPy expression. '
                #     "All variables in SymPy expressions passed to FD_c_codegen() must be registered "
                #     "in NRPy+ as either a gridfunction or Cparameter, by calling "
                #     f'{var} = register_gridfunctions...() (in ixp/grid) if "{var}" is a gridfunction, or '
                #     f"{var} = Cparameters() (in par) otherwise (e.g., if it is a free parameter set at C runtime)."
                # )
                # raise ValueError(err_msg)

    list_of_deriv_vars = superfast_uniq(list_of_deriv_vars_with_duplicates)
    # For upwinding with respect to a control vector (e.g., BSSN shift vector), for each variable
    # with suffix _dupD, append to the list_of_deriv_vars the corresponding _ddnD.
    if not isinstance(upwind_control_vec, str):
        list_of_deriv_vars.extend(
            sp.sympify(str(var).replace("_dupD", "_ddnD"))
            for var in list_of_deriv_vars
            if "_dupD" in str(var)
        )

    # Sort the list_of_deriv_vars for consistency in the C code output and potential cache miss reduction.
    return sorted(list_of_deriv_vars, key=sp.default_sort_key)


########################################
# STEP 2: EXTRACT BASE GRIDFUNCTIONS AND
#         DERIVATIVE OPERATORS FROM LIST
#         OF DERIVATIVES
def fd_temp_variable_name(
    gf_basename: str,
    i0_offset: int,
    i1_offset: int,
    i2_offset: int,
) -> str:
    """Generate a unique variable name for use in finite-difference code. This function takes into account
    the gridfunction name, and the offsets in the three dimensions.

    Example: If a grid function is named hDD00, and we want to read from memory data at i0+1,i1,i2-1,
    we store the value of this grid function as hDD00_i0p1_i1_i2m1; this function generates this name.
    varname = hDD00, i0_offset = 0, i1_offset = 2, i2_offset = 1 results in
    output: "hDD00_i1p2_i2p1"

    :param gf_basename: The name of the grid function
    :param i0_offset: The offset for i0
    :param i1_offset: The offset for i1
    :param i2_offset: The offset for i2
    :return: Returns a unique name for a point of data for a grid function

    >>> fd_temp_variable_name("hDD00", -2, 0, -1)
    'hDD00_i0m2_i2m1'
    >>> fd_temp_variable_name("hDD00", 0, 4, 3)
    'hDD00_i1p4_i2p3'
    >>> fd_temp_variable_name("hDD00", 0, 0, 0)
    'hDD00'
    >>> fd_temp_variable_name("", -2, 0, -1)
    Traceback (most recent call last):
    ...
    ValueError: gf_basename must be provided
    """

    if not gf_basename:
        raise ValueError("gf_basename must be provided")

    offsets = [i0_offset, i1_offset, i2_offset]
    suffixes = ["i0", "i1", "i2"]

    suffix_list = []
    for i in range(3):
        if offsets[i] < 0:
            suffix_list.append(suffixes[i] + "m" + str(abs(offsets[i])))
        elif offsets[i] > 0:
            suffix_list.append(suffixes[i] + "p" + str(offsets[i]))

    if len(suffix_list) == 0:  # All offsets are 0
        return gf_basename
    return f"{gf_basename}_{'_'.join(suffix_list)}"


def extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars(
    list_of_deriv_vars: List[sp.Basic],
) -> Tuple[List[str], List[str]]:
    """
    Extract from list_of_deriv_vars a list of base gridfunctions
    and a list of derivative operators.

    Args:
    list_of_deriv_vars (List[sp.Basic]): List of derivative variables.

    Returns:
    Tuple[List[str], List[str]]: Tuple containing lists of base gridfunctions and derivative operators.

    Doctest:
    >>> import nrpy.indexedexp as ixp
    >>> c_dD = ixp.declarerank1("c_dD")
    >>> aDD_dD = ixp.declarerank3("aDD_dD")
    >>> aDD_dKOD = ixp.declarerank3("aDD_dD")
    >>> vetU_dKOD = ixp.declarerank2("vetU_dKOD")
    >>> hDD_dDD = ixp.declarerank4("hDD_dDD")
    >>> extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars(
    ...    [c_dD[0], aDD_dD[0][1][2], aDD_dKOD[0][1][2], vetU_dKOD[2][1], hDD_dDD[0][1][1][2]])
    (['c', 'aDD01', 'aDD01', 'vetU2', 'hDD01'], ['dD0', 'dD2', 'dD2', 'dKOD1', 'dDD12'])
    """
    list_of_base_gridfunction_names_in_derivs = []
    list_of_deriv_operators = []

    for var in list_of_deriv_vars:
        varstr = str(var)

        # Step 2a.1: Check that the number of integers appearing
        #            in the suffix of a variable name matches the
        #            number of U's + D's in the variable name:
        num_UDs = varstr.count("D") + varstr.count("U")
        # length of string with digits stripped from the right---v
        num_digits_at_end = len(varstr) - len(varstr.rstrip("0123456789"))

        if num_UDs != num_digits_at_end:
            error_msg = (
                f"Error: {varstr} has {num_UDs} U's and D's, but "
                f"{num_digits_at_end} integers at the end. These must be equal. "
                "Please rename your gridfunction."
            )
            raise TypeError(error_msg)

        # Step 2a.2: Based on the variable name, find the rank of
        #            the underlying gridfunction of which we're
        #            trying to take the derivative.
        # rank = "number of juxtaposed Us and Ds before the underscore in a derivative expression"
        rank = 0
        underscore_position = varstr.rfind("_")  # Find the last occurrence of "_"
        if underscore_position != -1:
            # count contiguous "U"s and "D"s before underscore
            i = underscore_position - 1
            while i >= 0 and varstr[i] in ["U", "D"]:
                rank += 1
                i -= 1

        # Step 2a.3: Based on the variable name, find the order
        #            of the derivative we're trying to take.
        deriv_order = 0  # deriv_order = "number of Ds after the underscore in a derivative expression"
        for i in range(underscore_position + 1, len(varstr)):
            if varstr[i] == "D":
                deriv_order += 1

        # Step 2a.4: Based on derivative order and rank,
        #            store the base gridfunction name and
        #            derivative operator in the respective lists.
        base_gridfunction = (
            varstr[:underscore_position]
            + varstr[len(varstr) - deriv_order - rank : len(varstr) - deriv_order]
        )
        deriv_operator = (
            varstr[underscore_position + 1 : len(varstr) - deriv_order - rank]
            + varstr[len(varstr) - deriv_order :]
        )

        list_of_base_gridfunction_names_in_derivs.append(base_gridfunction)
        list_of_deriv_operators.append(deriv_operator)

    return list_of_base_gridfunction_names_in_derivs, list_of_deriv_operators


def read_gfs_from_memory(
    list_of_base_gridfunction_names_in_derivs: List[str],
    fdstencl: List[List[List[int]]],
    free_symbols_list: List[sp.Basic],
    mem_alloc_style: str,
    enable_simd: bool,
) -> str:
    """
    Generates C code to read grid functions from memory for finite difference derivative
    calculations. Memory reads are sorted to minimize cache misses.

    :param list_of_base_gridfunction_names_in_derivs: A list of base grid function names
        which are involved in derivative computations.
    :param fdstencl: A list of lists representing the finite difference stencil to be used.
    :param free_symbols_list: A list of symbols present in all SymPy expressions that are
        passed to the parent function. If there's only one symbol, you can pass it directly.
    :param mem_alloc_style: A string that denotes the memory allocation style. Options are
        '012' or '210'. The default '210' means the innermost loop is the "x" direction.
    :param enable_simd: A boolean flag that indicates whether SIMD (Single Instruction,
        Multiple Data) should be enabled or not.
    :param c_type: As string that specifies the C data type (e.g., float, double, REAL_SIMD_ARRAY, etc)
    :return: A string containing C code that represents the reading of grid functions at the
        necessary points in memory.

    >>> import nrpy.indexedexp as ixp
    >>> gri.glb_gridfcs_dict.clear()
    >>> par.set_parval_from_str("Infrastructure", "BHaH")
    >>> vU       = gri.register_gridfunctions_for_single_rank1("vU", group="EVOL")
    >>> hDD      = gri.register_gridfunctions_for_single_rank2("hDD", group="EVOL", symmetry="sym01")
    >>> hDD_dD   = ixp.declarerank3("hDD_dD", symmetry="sym01")
    >>> hDD_dupD = ixp.declarerank3("hDD_dupD", symmetry="sym01")
    >>> a0, a1, b, c = par.register_CodeParameters(c_type_alias="REAL", module=__name__, names=["a0", "a1", "b", "c"], defaultvalues=1)
    >>> exprlist = [b*hDD[1][0] + c*hDD_dD[0][1][1], c*hDD_dupD[0][2][2] + b*hDD_dupD[0][2][0] + a1*a0*vU[1]]
    >>> print(exprlist)
    [b*hDD01 + c*hDD_dD011, a0*a1*vU1 + b*hDD_dupD020 + c*hDD_dupD022]
    >>> free_symbols_list = []
    >>> for expr in exprlist:
    ...     free_symbols_list.extend(expr.free_symbols)
    >>> list_of_deriv_vars = extract_list_of_deriv_var_strings_from_sympyexpr_list(free_symbols_list, sp.Symbol("unset"))
    >>> print(list_of_deriv_vars)
    [hDD_dD011, hDD_ddnD020, hDD_ddnD022, hDD_dupD020, hDD_dupD022]
    >>> list_of_base_gridfunction_names_in_derivs, list_of_deriv_operators = extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars(list_of_deriv_vars)
    >>> print(list_of_base_gridfunction_names_in_derivs)
    ['hDD01', 'hDD02', 'hDD02', 'hDD02', 'hDD02']
    >>> fdcoeffs = [[] for _ in list_of_deriv_operators]
    >>> fdstencl = [[[] for _ in range(4)] for __ in list_of_deriv_operators]
    >>> for i, deriv_op in enumerate(list_of_deriv_operators): fdcoeffs[i], fdstencl[i] = compute_fdcoeffs_fdstencl(deriv_op, 2)
    >>> print(fdstencl)
    [[[0, -1, 0], [0, 1, 0]], [[-2, 0, 0], [-1, 0, 0], [0, 0, 0]], [[0, 0, -2], [0, 0, -1], [0, 0, 0]], [[0, 0, 0], [1, 0, 0], [2, 0, 0]], [[0, 0, 0], [0, 0, 1], [0, 0, 2]]]
    >>> print(list_of_deriv_operators)
    ['dD1', 'ddnD0', 'ddnD2', 'dupD0', 'dupD2']
    >>> print(read_gfs_from_memory(list_of_base_gridfunction_names_in_derivs, fdstencl, free_symbols_list, "210", enable_simd=False))
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
    <BLANKLINE>
    >>> gri.glb_gridfcs_dict["vU1"].gf_array_name = "simd_in_gfs_vU1"
    >>> print(read_gfs_from_memory(list_of_base_gridfunction_names_in_derivs, fdstencl, free_symbols_list, "012", enable_simd=True))
    const REAL_SIMD_ARRAY hDD01_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1-1, i2)]);
    const REAL_SIMD_ARRAY hDD01 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2)]);
    const REAL_SIMD_ARRAY hDD01_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1+1, i2)]);
    const REAL_SIMD_ARRAY hDD02_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0-2, i1, i2)]);
    const REAL_SIMD_ARRAY hDD02_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0-1, i1, i2)]);
    const REAL_SIMD_ARRAY hDD02_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2-2)]);
    const REAL_SIMD_ARRAY hDD02_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2-1)]);
    const REAL_SIMD_ARRAY hDD02 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2)]);
    const REAL_SIMD_ARRAY hDD02_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2+1)]);
    const REAL_SIMD_ARRAY hDD02_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2+2)]);
    const REAL_SIMD_ARRAY hDD02_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0+1, i1, i2)]);
    const REAL_SIMD_ARRAY hDD02_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0+2, i1, i2)]);
    const REAL_SIMD_ARRAY vU1 = ReadSIMD(&simd_in_gfs_vU1[IDX4(VU1GF, i0, i1, i2)]);
    <BLANKLINE>
    """

    # Step 4a: Compile list of points to read from memory
    #          for each gridfunction i, based on list
    #          provided in fdstencil[i][].
    list_of_points_read_from_memory_with_duplicates: Dict[str, List[str]] = {
        k: [] for k in gri.glb_gridfcs_dict
    }
    for j, derivgfname in enumerate(list_of_base_gridfunction_names_in_derivs):
        if derivgfname in gri.glb_gridfcs_dict:
            list_of_points_read_from_memory_with_duplicates[derivgfname].extend(
                [f"{coord[0]},{coord[1]},{coord[2]}" for coord in fdstencl[j]]
            )

    # Step 4b: "Zeroth derivative" case:
    #     If gridfunction appears in expression not
    #     as derivative (i.e., by itself), it must
    #     be read from memory as well.
    for var in free_symbols_list:
        if symbol_is_gridfunction_Cparameter_or_other(var) == "gridfunction":
            if str(var) in gri.glb_gridfcs_dict:
                list_of_points_read_from_memory_with_duplicates[str(var)].append(
                    "0,0,0"
                )

    # Step 4c: Remove duplicates when reading from memory;
    #     do not needlessly read the same variable
    #     from memory twice.
    list_of_points_read_from_memory = [
        superfast_uniq(v)
        for v in list_of_points_read_from_memory_with_duplicates.values()
    ]

    # Step 4d: Minimize cache misses:
    #      Sort the list of points read from
    #      main memory by how they are stored
    #      in memory.

    # Step 4d.i: Define a function that maps a gridpoint
    #     index (i,j,k) to a unique memory "address",
    #     which will correspond to the correct ordering
    #     of actual memory addresses.
    #
    def unique_idx(idx3: List[int], mem_alloc_style: str) -> int:
        """
        Generates a unique index that ensures the indices are ordered in memory.
        This function uses a given offset and size, and supports two memory allocation styles: '210' and '012'.
        Note: The values of offset and size are used solely for ordering purposes and should not be modified.

        :param idx3: A list of three indices.
        :param mem_alloc_style: The memory allocation style, either '210' or '012'.
        :return: A unique index based on the input indices and memory allocation style.
        :raises ValueError: If the provided memory allocation style is unsupported.
        """
        offset = 50  # Offset for memory ordering. Do not modify.
        size = 100  # Assumed size in each direction. Do not modify.

        if mem_alloc_style == "210":
            return (
                idx3[0] + offset + size * (idx3[1] + offset + size * idx3[2] + offset)
            )
        if mem_alloc_style == "012":
            return (
                idx3[2] + offset + size * (idx3[1] + offset + size * idx3[0] + offset)
            )
        raise ValueError(f"Error: mem_alloc_style = {mem_alloc_style} is unsupported.")

    # Step 4d.ii: For each gridfunction and
    #      point read from memory, call unique_idx,
    #      then sort according to memory "address"
    def sort_memory_accesses(
        list_of_points_read_from_memory: List[List[str]], mem_alloc_style: str
    ) -> List[List[str]]:
        """
        Process and sort the list of memory access points for each grid function.

        This function sorts a given list of memory access points for each grid function
        by transforming them into unique indices and ordering them. The transformation
        is performed by the nested helper function memory_idx_from_str. The sorted list
        of memory access points enables optimal access pattern during data processing.

        Args:
            list_of_points_read_from_memory: A nested list where each sublist corresponds
                to a grid function and contains points read from memory for that grid function.
            mem_alloc_style: A string that denotes the memory allocation style used in
                transforming the memory address to a unique index.

        Returns:
            A nested list where each sublist corresponds to a grid function and contains
            the sorted points read from memory for that grid function.
        """

        def memory_idx_from_str(idx_str: str, mem_alloc_style: str) -> int:
            """Helper function to get unique index from a string of indices"""
            idx3 = [int(x) for x in idx_str.split(",")]
            return unique_idx(idx3, mem_alloc_style)

        sorted_list_of_points_read_from_memory: List[List[str]] = [
            [] for _ in gri.glb_gridfcs_dict
        ]

        for gfidx, points in enumerate(list_of_points_read_from_memory):
            if points:  # Continue only if reading at least one point from memory.
                # Convert the point indices to memory indices based on a memory allocation style.
                memory_indices = [
                    memory_idx_from_str(idx, mem_alloc_style) for idx in points
                ]

                # Zip together memory indices and points, then sort based on memory indices.
                # Finally, unzip them to get the sorted points.
                _, sorted_points = zip(
                    *sorted(zip(memory_indices, points), key=itemgetter(0))
                )
                sorted_list_of_points_read_from_memory[gfidx] = list(sorted_points)
        return sorted_list_of_points_read_from_memory

    sorted_list_of_points_read_from_memory = sort_memory_accesses(
        list_of_points_read_from_memory, mem_alloc_style
    )

    read_gf_from_memory_Ccode = ""

    # Initialize these to be unsorted.
    name_sorted_gfs = list(gri.glb_gridfcs_dict.keys())
    name_sorted_gfs_idxs = list(range(len(name_sorted_gfs)))

    # Pair each name with its corresponding index
    pairs = list(zip(name_sorted_gfs, name_sorted_gfs_idxs))

    # Sort the pairs. Each pair is a tuple, and sorted() will sort them based on the first element of the tuple by default
    sorted_pairs = sorted(pairs)

    # Unpack sorted_pairs list into two lists (name_sorted_gfs, name_sorted_gfs_idxs) with corresponding elements grouped together
    # MyPy really gets confused here. zips are just too complex for poor MyPy.
    name_sorted_gfs, name_sorted_gfs_idxs = zip(*sorted_pairs)  # type: ignore

    # Iterate over pairs of grid function names and their indices.
    # For each gridfunction, check if there are points in memory that need read.
    # If points exist, they are parsed and used to generate C code for
    # reading the grid function from memory, possibly utilizing SIMD instructions
    # for certain data types. The generated C code is then added to the
    # 'read_gf_from_memory_Ccode' string.
    for gfname, idx in zip(name_sorted_gfs, name_sorted_gfs_idxs):
        gf = gri.glb_gridfcs_dict[gfname]
        c_type_alias = gf.c_type_alias
        if enable_simd:
            if par.parval_from_str("Infrastructure") in ("NRPy", "BHaH", "BaseETK"):
                c_type_alias = "REAL_SIMD_ARRAY"
            else:
                raise ValueError("FIXME: Please specify the c_type for SIMD")
        points = sorted_list_of_points_read_from_memory[idx]
        # point is a string of integer offsets that e.g., looks like: -1,2,0
        if points:
            for point in points:
                # convert point to set of integers
                i0_offset, i1_offset, i2_offset = map(int, point.split(","))
                read_gf_from_memory_Ccode += f"const {c_type_alias} {fd_temp_variable_name(gf.name, i0_offset, i1_offset, i2_offset)} = {gf.read_gf_from_memory_Ccode_onept(i0_offset, i1_offset, i2_offset, enable_simd=enable_simd)};\n"

    return read_gf_from_memory_Ccode


def proto_FD_operators_to_sympy_expressions(
    list_of_proto_deriv_symbs: List[sp.Symbol],
    fdcoeffs: List[List[sp.Rational]],
    fdstencl: List[List[List[int]]],
    enable_simd: bool = False,
) -> Tuple[List[sp.Basic], List[str]]:
    """
    Convert finite difference (FD) operators to SymPy expressions.

    :param list_of_proto_deriv_symbs: List of prototype derivative variables (e.g., [sp.Symbol("dDD12")])
    :param fdcoeffs: List of finite difference coefficients.
    :param fdstencl: List of finite difference stencils.
    :param enable_simd: Whether to enable SIMD.

    :return: Tuple containing the list of SymPy expressions for the finite difference operators and the corresponding left-hand side variable names.

    >>> import nrpy.indexedexp as ixp
    >>> gri.glb_gridfcs_dict.clear()
    >>> par.set_parval_from_str("Infrastructure", "BHaH")
    >>> dum = gri.register_gridfunctions("FDPROTO")[0]
    >>> dum_dD   = ixp.declarerank1("FDPROTO_dD")
    >>> dum_dupD = ixp.declarerank1("FDPROTO_dupD")
    >>> dum_dKOD = ixp.declarerank1("FDPROTO_dKOD")
    >>> dum_dDD  = ixp.declarerank2("FDPROTO_dDD", symmetry="sym01")
    >>> vU = ixp.declarerank1("vU")
    >>> a0, a1, b, c = par.register_CodeParameters(c_type_alias="REAL", module=__name__, names=["a0", "a1", "b", "c"], defaultvalues=1)
    >>> exprlist = [b*dum_dDD[1][0] + c*dum - a0*dum_dupD[2], c*dum_dKOD[1] + a1*a0*vU[1]]
    >>> free_symbols_list = []
    >>> for expr in exprlist:
    ...     free_symbols_list.extend(expr.free_symbols)
    >>> list_of_proto_deriv_symbs = extract_list_of_deriv_var_strings_from_sympyexpr_list(free_symbols_list, sp.Symbol("unset"))
    >>> print(list_of_proto_deriv_symbs)
    [FDPROTO_dDD01, FDPROTO_dKOD1, FDPROTO_ddnD2, FDPROTO_dupD2]
    >>> list_of_proto_deriv_ops = [str(item).split("_")[1] for item in list_of_proto_deriv_symbs]
    >>> print(list_of_proto_deriv_ops)
    ['dDD01', 'dKOD1', 'ddnD2', 'dupD2']
    >>> fdcoeffs = [[] for _ in list_of_proto_deriv_ops]
    >>> fdstencl = [[[] for _ in range(4)] for __ in list_of_proto_deriv_ops]
    >>> for i, deriv_op in enumerate(list_of_proto_deriv_ops): fdcoeffs[i], fdstencl[i] = compute_fdcoeffs_fdstencl(deriv_op, 2)
    >>> print(fdstencl)
    [[[-1, -1, 0], [1, -1, 0], [-1, 1, 0], [1, 1, 0]], [[0, -2, 0], [0, -1, 0], [0, 0, 0], [0, 1, 0], [0, 2, 0]], [[0, 0, -2], [0, 0, -1], [0, 0, 0]], [[0, 0, 0], [0, 0, 1], [0, 0, 2]]]
    >>> FDexprs, FDlhsvarnames = proto_FD_operators_to_sympy_expressions(list_of_proto_deriv_symbs, fdcoeffs, fdstencl)
    >>> for i, lhs in enumerate(FDlhsvarnames):
    ...     print(f"{lhs} = {FDexprs[i]}")
    const REAL FDPROTO_dDD01 = invdxx0*invdxx1*(FDPROTO_i0m1_i1m1/4 - FDPROTO_i0m1_i1p1/4 - FDPROTO_i0p1_i1m1/4 + FDPROTO_i0p1_i1p1/4)
    const REAL FDPROTO_dKOD1 = invdxx1*(-3*FDPROTO/8 + FDPROTO_i1m1/4 - FDPROTO_i1m2/16 + FDPROTO_i1p1/4 - FDPROTO_i1p2/16)
    const REAL UpwindAlgInputFDPROTO_ddnD2 = invdxx2*(3*FDPROTO/2 - 2*FDPROTO_i2m1 + FDPROTO_i2m2/2)
    const REAL UpwindAlgInputFDPROTO_dupD2 = invdxx2*(-3*FDPROTO/2 + 2*FDPROTO_i2p1 - FDPROTO_i2p2/2)
    >>> FDexprs, FDlhsvarnames = proto_FD_operators_to_sympy_expressions(list_of_proto_deriv_symbs, fdcoeffs, fdstencl, enable_simd=True)
    >>> for lhs in FDlhsvarnames:
    ...     print(f"{lhs}")
    const REAL_SIMD_ARRAY FDPROTO_dDD01
    const REAL_SIMD_ARRAY FDPROTO_dKOD1
    const REAL_SIMD_ARRAY UpwindAlgInputFDPROTO_ddnD2
    const REAL_SIMD_ARRAY UpwindAlgInputFDPROTO_dupD2
    """
    # Set invdxx symbol list:
    invdxx = [sp.sympify(f"invdxx{d}") for d in range(3)]
    FDexprs = [sp.sympify(0)] * len(list_of_proto_deriv_symbs)
    FDlhsvarnames = [""] * len(list_of_proto_deriv_symbs)

    # Step 5.a.ii.A: Output finite difference expressions to Coutput string
    list_of_proto_deriv_ops = [
        str(item).split("_")[1] for item in list_of_proto_deriv_symbs
    ]
    for i, proto_deriv_var_symbol in enumerate(list_of_proto_deriv_symbs):
        # Unpack:
        operator = list_of_proto_deriv_ops[i]
        proto_deriv_var = str(proto_deriv_var_symbol)
        if "_dupD" in proto_deriv_var or "_ddnD" in proto_deriv_var:
            proto_deriv_var = f"UpwindAlgInput{proto_deriv_var}"

        # First add element to FDlhsvarnames:
        c_type_alias = gri.glb_gridfcs_dict[
            next(iter(gri.glb_gridfcs_dict))
        ].c_type_alias  # next(iter( grabs any key in the dict.
        if enable_simd:
            if c_type_alias in ("REAL", "CCTK_REAL", "double"):
                c_type_alias = "REAL_SIMD_ARRAY"
        FDlhsvarnames[i] = f"const {c_type_alias} {proto_deriv_var}"

        # Add element to FDexprs:
        for j in range(len(fdcoeffs[i])):
            varname = fd_temp_variable_name(
                "FDPROTO", fdstencl[i][j][0], fdstencl[i][j][1], fdstencl[i][j][2]
            )
            FDexprs[i] += fdcoeffs[i][j] * sp.sympify(varname)

        # Check if the last character in the op string is an integer
        if not operator[-1].isdigit():
            raise ValueError(f"Error: Expected an integer at the end of op: {operator}")

        # Extract the integer at the end of the op string
        direction = int(operator[-1])

        # Second-order derivatives:
        if operator.startswith("dDD"):
            if not operator[-2].isdigit():
                raise ValueError(
                    f"Error: Expected two integers at the end of op: {operator}"
                )
            direction1 = int(operator[-2])
            direction2 = direction
            FDexprs[i] *= invdxx[direction1] * invdxx[direction2]
        # First-order or Kreiss-Oliger derivatives:
        elif operator.startswith(("dKOD", "dD", "dupD", "ddnD")):
            FDexprs[i] *= invdxx[direction]
        else:
            raise ValueError(
                f"Error: Was unable to parse derivative operator: {operator}"
            )

    return FDexprs, FDlhsvarnames


################
# def output_finite_difference_functions_h(path=os.path.join(".")):
#     with open(os.path.join(path, "finite_difference_functions.h"), "w") as file:
#         file.write(
#             """
# #ifndef __FD_FUNCTIONS_H__
# #define __FD_FUNCTIONS_H__
# #include "math.h"
# #include "stdio.h"
# #include "stdlib.h"
# """
#         )
#         UNUSED = "__attribute__((unused))"
#         NOINLINE = "__attribute__((noinline))"
#         if par.parval_from_str("grid::GridFuncMemAccess") == "ETK":
#             UNUSED = "CCTK_ATTRIBUTE_UNUSED"
#             NOINLINE = "CCTK_ATTRIBUTE_NOINLINE"
#         file.write("#define _UNUSED   " + UNUSED + "\n")
#         file.write("#define _NOINLINE " + NOINLINE + "\n")
#
#         for key, item in outC_function_dict.items():
#             if "__FD_OPERATOR_FUNC__" in item:
#                 file.write(
#                     item.replace(
#                         "const REAL_SIMD_ARRAY _NegativeOne_ =",
#                         "const REAL_SIMD_ARRAY " + UNUSED + " _NegativeOne_ =",
#                     )
#                 )  # Many of the NegativeOne's get optimized away in the SIMD postprocessing step. No need for all the warnings
#
#         # Clear all FD functions from outC_function_dict after outputting to finite_difference_functions.h.
#         #   Otherwise c_codegen will be outputting these as separate individual C codes & attempting to build them in Makefile.
#         key_list_del = []
#         element_del = []
#         for i, func in enumerate(outC_function_master_list):
#             if "__FD_OPERATOR_FUNC__" in func.desc:
#                 if func.name not in key_list_del:
#                     key_list_del += [func.name]
#                 if func not in element_del:
#                     element_del += [func]
#         for func in element_del:
#             outC_function_master_list.remove(func)
#         for key in key_list_del:
#             outC_function_dict.pop(key)
#             if key in outC_function_prototype_dict:
#                 outC_function_prototype_dict.pop(key)
#         file.write("#endif // #ifndef __FD_FUNCTIONS_H__\n")


################


# def add_FD_func_to_outC_function_dict(fd_order, enable_simd=False):
#     # Step 5.a.ii.A: First construct a list of all the unique finite difference functions
#     c_type = "REAL"
#     if par.parval_from_str("grid::GridFuncMemAccess") == "ETK":
#         c_type = "CCTK_REAL"
#     func_prefix = "order_" + str(fd_order) + "_"
#     if enable_simd:
#         c_type = "REAL_SIMD_ARRAY"
#         func_prefix = "SIMD_" + func_prefix
#
#     # Stores the needed calls to the functions we're adding to outC_function_dict:
#     FDfunccall_list = []
#     for op in list_of_uniq_deriv_operators:
#         which_op_idx = find_which_op_idx(op, list_of_deriv_operators)
#
#         rhs_expr = sp.sympify(0)
#         for j in range(len(fdcoeffs[which_op_idx])):
#             var = sp.sympify("f" + varsuffix("f", fdstencl[which_op_idx][j], FDparams))
#             rhs_expr += fdcoeffs[which_op_idx][j] * var
#
#         # Multiply each expression by the appropriate power
#         #   of 1/dx[i]
#         used_invdxx = [False, False, False]
#         invdxx = [sp.sympify(f"invdxx{d}") for d in range(3)]
#
#         # Check if the last character in the op string is an integer
#         if not op[-1].isdigit():
#             raise ValueError(f"Error: Expected an integer at the end of op: {op}")
#
#         # Extract the integer at the end of the op string
#         direction = int(op[-1])
#
#         # Second-order derivatives:
#         if op.startswith("dDD"):
#             if not op[-2].isdigit():
#                 raise ValueError(f"Error: Expected two integers at the end of op: {op}")
#             direction1 = int(op[-2])
#             direction2 = direction
#             used_invdxx[direction1] = used_invdxx[direction2] = True
#             rhs_expr *= invdxx[direction1] * invdxx[direction2]
#         # First-order or Kreiss-Oliger derivatives:
#         elif op.startswith(("dKOD", "dD", "dupD", "ddnD")):
#             rhs_expr *= invdxx[direction]
#             used_invdxx[direction] = True
#         else:
#             raise ValueError(f"Error: Was unable to parse derivative operator: {op}")
#
#         invdxx_params = [f"const {c_type} invdxx{d}" for d in range(3) if used_invdxx[d]]
#         fdcoeffs_params = [
#             f"const {c_type} f{varsuffix('', fdstencl[which_op_idx][j], FDparams)}"
#             for j in range(len(fdcoeffs[which_op_idx]))
#         ]
#
#         outfunc_params = ", ".join(invdxx_params + fdcoeffs_params)
#
#         for i, deriv_operator in enumerate(list_of_deriv_operators):
#             if deriv_operator == op:
#                 invdxx_funccall = [f"invdxx{d}" for d in range(3) if used_invdxx[d]]
#                 gfname = list_of_base_gridfunction_names_in_derivs[i]
#                 fdcoeffs_funccall = [
#                     f"{gfname}{varsuffix(gfname, fdstencl[which_op_idx][j], FDparams)}"
#                     for j in range(len(fdcoeffs[which_op_idx]))
#                 ]
#
#                 funccall = f"{type__var(list_of_deriv_vars[i], FDparams)} = {func_prefix}f_{op}({', '.join(invdxx_funccall + fdcoeffs_funccall)});"
#                 FDfunccall_list.append(funccall)
#
#         func_name = f"{func_prefix}f_{op}"
#         if func_name not in outC_function_dict:
#             p = f"preindent=1,enable_simd={FDparams.enable_simd},outCverbose=False,cse_preprocess=True,include_braces=False"
#             outFDstr = outputC(rhs_expr, "retval", "returnstring", params=p).replace(
#                 "retval = ", "return "
#             )
#
#             op_description = (
#                 str(op)
#                 .replace("dDD", "second derivative: ")
#                 .replace("dD", "first derivative: ")
#                 .replace("dKOD", "Kreiss-Oliger derivative: ")
#                 .replace("dupD", "upwinded derivative: ")
#                 .replace("ddnD", "downwinded derivative: ")
#             )
#             desc = f" * (__FD_OPERATOR_FUNC__) Finite difference operator for {op_description} direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively."
#
#             register_CFunction(
#                 desc=desc,
#                 c_type=f"static {c_type} _NOINLINE _UNUSED",
#                 name=func_name,
#                 enableCparameters=False,
#                 params=outfunc_params,
#                 preloop="",
#                 body=outFDstr,
#             )
#
#     return FDfunccall_list


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
