"""
Provide helper functions for c_codegen to generate finite-difference C-code kernels.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import sys  # Standard Python module for multiplatform OS-level functions
from typing import Union, List, Tuple, Any, Dict
from operator import itemgetter
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import nrpy.params as par  # NRPy+: parameter interface
import nrpy.grid as gri  # NRPy+: Functions having to do with numerical grids
import nrpy.c_function as cfc
from nrpy.helpers.generic import superfast_uniq
from nrpy.helpers.cse_preprocess_postprocess import cse_preprocess

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
    Set up finite difference matrix and return its inverse.

    :param stencil_width: Width of the stencil.
    :param UPDOWNWIND_stencil_shift: Shift in the stencil for upwind or downwind.

    :return: Inverse of the finite difference matrix.

    Doctest:
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
    """
    Custom dictionary for storing FD coefficient matrices, as the inversion is expensive.

    :param key: A tuple containing stencil width and UPDOWNWIND stencil shift.
    :return: A sympy matrix associated with the provided key.
    """

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
    Set up finite difference matrix and return its inverse.
    If the inputs have already been processed through
    setup_FD_matrix__return_inverse_lowlevel() (i.e., the matrix has already
    been inverted), this function will simply return the inverse matrix stored within
    FD_Matrix_dict((stencil_width, UPDOWNWIND_stencil_shift). Otherwise, it will
    populate the FD_Matrix_dict with this matrix.

    :param stencil_width: Width of the stencil.
    :param UPDOWNWIND_stencil_shift: Shift in the stencil for upwind or downwind.
    :return: Inverse of the finite difference matrix.

    Doctest:
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

    :param derivstring: The string representation of the derivative type. Can include patterns
                        like "dKOD", "dupD", "ddnD", "dfullupD", "dfulldnD", "DDD", and "DD" to
                        dictate how the function processes and returns coefficients.
    :param fd_order: Order of the finite differencing.

    :return: A tuple containing two lists. The first list contains the finite difference coefficients
             as sympy.Rational numbers, and the second list contains the stencils as lists of integers.

    Note:
    This function computes the finite difference coefficients and stencil points for various derivative types
    specified in `derivstring`. The coefficients are determined using the inverse of the finite difference matrix,
    which is constructed and inverted in the `setup_FD_matrix__return_inverse` function.
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
    Determine the type of a given variable.

    This function checks whether a given variable is in the list of global grid functions or
    C parameters, and returns a string indicating its type: either 'gridfunction', 'Cparameter',
    or 'other'.

    :param var: The variable whose type is to be determined.

    :return: A string indicating the type of the variable. Can be 'gridfunction', 'Cparameter', or 'other'.

    :raises ValueError: If the variable is both a gridfunction and a Cparameter, or if it is neither and
                        cannot be identified.

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

    :param list_of_free_symbols: List of free symbols from SymPy expressions.
    :param upwind_control_vec: Upwind control vector.

    :returns: List of derivative variables, creating _ddnD in case upwinding is enabled with control vector.

    :raises ValueError: If a variable in the SymPy expression isn't registered as a gridfunction or Cparameter.

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
    """
    Generate a unique variable name for use in finite-difference code.

    This function takes into account the gridfunction name, and the offsets in the three dimensions to generate
    a unique name for a point of data for a grid function.

    :param gf_basename: The name of the grid function.
    :param i0_offset: The offset for i0.
    :param i1_offset: The offset for i1.
    :param i2_offset: The offset for i2.

    :return: A unique name for a point of data for the grid function.

    :raises ValueError: If the gf_basename is not provided.

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
    Extract lists of base gridfunctions and derivative operators from list_of_deriv_vars.

    This function will take a list of derivative variables and extract the base gridfunctions
    and derivative operators from them.

    :param list_of_deriv_vars: List of derivative variables.

    :return: Tuple containing lists of base gridfunctions and derivative operators.

    :raises TypeError: If the number of integers appearing in the suffix of a variable name
                       does not match the number of U's + D's in the variable name.

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
    Generate C code to read grid functions from memory for finite difference derivative calculations.

    The generated C code will take into account the grid functions and the finite difference stencil
    provided to optimize memory reads and reduce cache misses.

    :param list_of_base_gridfunction_names_in_derivs: A list of base grid function names involved in derivative computations.
    :param fdstencl: A list of lists representing the finite difference stencil to be used.
    :param free_symbols_list: A list of symbols present in all SymPy expressions passed to the parent function.
                              If there's only one symbol, it can be passed directly.
    :param mem_alloc_style: Specifies the memory allocation style. Options are '012' or '210'.
                            Default '210' means the innermost loop is the "x" direction.
    :param enable_simd: Indicates whether SIMD (Single Instruction, Multiple Data) should be enabled or not.

    :return: A string containing C code representing the reading of grid functions at the necessary points in memory.

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
        Generate a unique index ensuring the indices are ordered in memory.

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
        by transforming them into unique indices and ordering them. The sorted list
        of memory access points enables optimal access pattern during data processing.

        :param list_of_points_read_from_memory: A nested list where each sublist corresponds to a grid function and
            contains points read from memory for that grid function.
        :param mem_alloc_style: A string denoting the memory allocation style used in transforming the memory
            address to a unique index.

        :return: A nested list where each sublist corresponds to a grid function and contains the sorted points
            read from memory for that grid function.
        """

        def memory_idx_from_str(idx_str: str, mem_alloc_style: str) -> int:
            """
            Convert a string of indices into a unique memory index.

            This helper function takes a string representation of indices and returns a unique
            index based on the provided memory allocation style.

            :param idx_str: A string representation of indices, separated by commas.
            :param mem_alloc_style: The memory allocation style, either '210' or '012'.

            :return: A unique memory index integer.
            """
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
            if par.parval_from_str("Infrastructure") in ("BHaH", "CarpetX", "ETLegacy"):
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


class FDFunction:
    """
    A class to represent Finite-Difference (FD) functions in C.

    :param c_type_alias: The alias for the C data type used in the function.
    :param fd_order: The order of accuracy for the finite difference scheme.
    :param operator: The operator with respect to which the derivative is taken.
    :param symbol_to_Rational_dict: Dictionary mapping sympy symbols to their corresponding sympy Rationals.
    :param FDexpr: The sympy expression representing the finite-difference formula.
    :param enable_simd: A flag to specify if SIMD instructions should be used.
    """

    def __init__(
        self,
        c_type_alias: str,
        fd_order: int,
        operator: str,
        symbol_to_Rational_dict: Dict[sp.Basic, sp.Rational],
        FDexpr: sp.Basic,
        enable_simd: bool,
    ) -> None:
        self.c_type_alias = c_type_alias
        self.fd_order = fd_order
        self.operator = operator
        self.symbol_to_Rational_dict = symbol_to_Rational_dict
        self.FDexpr = FDexpr
        self.enable_simd = enable_simd
        self.c_function_name = "SIMD_" if enable_simd else ""
        self.c_function_name += f"fd_function_{self.operator}_fdorder{self.fd_order}"

        self.CFunction: cfc.CFunction

    def c_function_call(self, gf_name: str, deriv_var: str) -> str:
        """
        Generate the C function call for a given grid function name and derivative variable.

        :param gf_name: The name of the grid function.
        :param deriv_var: The variable that represents the derivative.

        :return: The C function call as a string.
        """
        if "_dupD" in deriv_var or "_ddnD" in deriv_var:
            deriv_var = f"UpwindAlgInput{deriv_var}"
        c_function_call = (
            f"const {self.c_type_alias} {deriv_var} = {self.c_function_name}("
        )
        c_function_call += ",".join(
            sorted(
                str(symb)
                for symb in self.FDexpr.free_symbols
                if "FDPart1_" not in str(symb)
            )
        )

        c_function_call = c_function_call.replace("FDPROTO", gf_name)
        c_function_call += ")"
        return c_function_call

    def CFunction_fd_function(self, FDexpr_c_code: str) -> cfc.CFunction:
        """
        Generate a C function based on the given finite-difference expression.

        :param FDexpr_c_code: The finite-difference expression in C code format.

        :return: A cfc.CFunction object that encapsulates the C function details.
        """
        includes: List[str] = []
        c_type_alias = self.c_type_alias
        name = self.c_function_name
        params = ""
        params += ",".join(
            sorted(
                f"const {c_type_alias} {str(symb)}"
                for symb in self.FDexpr.free_symbols
                if "FDPart1_" not in str(symb)
            )
        )

        body = f"{FDexpr_c_code}\n return FD_result;"

        return cfc.CFunction(
            includes=includes,
            desc=f"Finite difference function for operator {self.operator}, with FD accuracy order {self.fd_order}.",
            c_type=f"static {c_type_alias}",
            name=name,
            params=params,
            body=body,
        )


FDFunctions_dict: Dict[str, FDFunction] = {}


def construct_FD_functions_prefunc() -> str:
    """
    Construct the prefunc (CFunction) strings for all finite-difference functions stored in FDFunctions_dict.

    :return: The concatenated prefunc (CFunction) strings as a single string.
    """
    prefunc = ""
    for fd_func in FDFunctions_dict.values():
        prefunc += fd_func.CFunction.full_function
    return prefunc


def proto_FD_operators_to_sympy_expressions(
    list_of_proto_deriv_symbs: List[sp.Symbol],
    fd_order: int,
    fdcoeffs: List[List[sp.Rational]],
    fdstencl: List[List[List[int]]],
    enable_simd: bool = False,
) -> Tuple[List[sp.Basic], List[str], Dict[sp.Basic, sp.Rational]]:
    """
    Convert finite difference (FD) operators to SymPy expressions.

    :param list_of_proto_deriv_symbs: List of prototype derivative variables (e.g., [sp.Symbol("dDD12")])
    :param fd_order: Finite-difference accuracy order
    :param fdcoeffs: List of finite difference coefficients.
    :param fdstencl: List of finite difference stencils.
    :param enable_simd: Whether to enable SIMD.

    :return: Tuple containing the list of SymPy expressions for the finite difference operators and the corresponding left-hand side variable names.

    >>> import nrpy.indexedexp as ixp
    >>> gri.glb_gridfcs_dict.clear()
    >>> par.set_parval_from_str("Infrastructure", "BHaH")
    >>> fd_order = 2
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
    >>> for i, deriv_op in enumerate(list_of_proto_deriv_ops): fdcoeffs[i], fdstencl[i] = compute_fdcoeffs_fdstencl(deriv_op, fd_order)
    >>> print(fdstencl)
    [[[-1, -1, 0], [1, -1, 0], [-1, 1, 0], [1, 1, 0]], [[0, -2, 0], [0, -1, 0], [0, 0, 0], [0, 1, 0], [0, 2, 0]], [[0, 0, -2], [0, 0, -1], [0, 0, 0]], [[0, 0, 0], [0, 0, 1], [0, 0, 2]]]
    >>> FDexprs, FDlhsvarnames, symbol_to_Rational_dict = proto_FD_operators_to_sympy_expressions(list_of_proto_deriv_symbs, fd_order, fdcoeffs, fdstencl)
    >>> for i, lhs in enumerate(FDlhsvarnames):
    ...     print(f"{lhs} = {FDexprs[i]}")
    const REAL FDPROTO_dDD01 = FDPart1_Rational_1_4*invdxx0*invdxx1*(FDPROTO_i0m1_i1m1 - FDPROTO_i0m1_i1p1 - FDPROTO_i0p1_i1m1 + FDPROTO_i0p1_i1p1)
    const REAL FDPROTO_dKOD1 = invdxx1*(-FDPROTO*FDPart1_Rational_3_8 + FDPart1_Rational_1_16*(-FDPROTO_i1m2 - FDPROTO_i1p2) + FDPart1_Rational_1_4*(FDPROTO_i1m1 + FDPROTO_i1p1))
    const REAL UpwindAlgInputFDPROTO_ddnD2 = invdxx2*(FDPROTO*FDPart1_Rational_3_2 - FDPROTO_i2m1*FDPart1_Integer_2 + FDPROTO_i2m2*FDPart1_Rational_1_2)
    const REAL UpwindAlgInputFDPROTO_dupD2 = invdxx2*(-FDPROTO*FDPart1_Rational_3_2 + FDPROTO_i2p1*FDPart1_Integer_2 - FDPROTO_i2p2*FDPart1_Rational_1_2)
    >>> FDexprs, FDlhsvarnames, symbol_to_Rational_dict = proto_FD_operators_to_sympy_expressions(list_of_proto_deriv_symbs, fd_order, fdcoeffs, fdstencl, enable_simd=True)
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
    symbol_to_Rational_dicts: List[Dict[sp.Basic, sp.Rational]] = [{}] * len(
        list_of_proto_deriv_symbs
    )

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

        # Factorize proto_FDexprs and convert Rationals to symbols, store symbol->Rational dictionary.
        processed_FDexpr, symbol_to_Rational_dicts[i] = cse_preprocess(
            FDexprs[i],
            prefix="FDPart1",
            declare_neg1_as_symbol=enable_simd,
            negative=enable_simd,
            factor=True,
        )
        FDexprs[i] = processed_FDexpr[0]

        FDFunctions_dict[operator] = FDFunction(
            c_type_alias=c_type_alias,
            fd_order=fd_order,
            operator=operator,
            symbol_to_Rational_dict=symbol_to_Rational_dicts[i],
            FDexpr=FDexprs[i],
            enable_simd=enable_simd,
        )

    symbol_to_Rational_dict = {}

    for i, op in enumerate(FDFunctions_dict.keys()):
        FDFunctions_dict[op].FDexpr = FDexprs[i]
        symbol_to_Rational_dict.update(symbol_to_Rational_dicts[i])
    return FDexprs, FDlhsvarnames, symbol_to_Rational_dict


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
