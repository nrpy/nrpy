"""
indexedexp.py: functions related to indexed expressions,
 including e.g., tensors and pseudotensors:

Authors: Zachariah Etienne, Kenneth Sible, Steven Brandt
"""

# Step 1: Load needed modules
import string
from typing import Union, List, Optional, Any, Tuple, cast, Sequence
import sys  # Standard Python module for multiplatform OS-level functions
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import nrpy.helpers.functional as func  # NRPy+: Python toolkit for functional programming
import nrpy.params as par


# Used to specify axes of symmetry so that derivatives across these axes get set to zero
par.register_param(py_type=str, module=__name__, name="symmetry_axes", value="")

_symbol_type = sp.Expr
_rank1_type = Sequence[_symbol_type]
_rank2_type = Sequence[_rank1_type]
_rank3_type = Sequence[_rank2_type]
_rank4_type = Sequence[_rank3_type]
_recur_symbol_type = Union[sp.Expr, Sequence["_recur_symbol_type"]]


# This function used to be called _init()
def create_tensor_symbolic(
    shape: Union[int, List[int]],
    symbol: Union[str, None] = "",
    preindex: Optional[List[int]] = None,
    character_zero_index: Optional[str] = "",
) -> Sequence[_recur_symbol_type]:
    """
    Creates a multi-dimensional symbolic tensor of a specified shape.

    Args:
    shape (int or list of int): The shape of the tensor. If an int is given, a 1D tensor will be created.
    symbol (str, optional): The base symbol name for elements in the tensor. Defaults to "".
    preindex (list, optional): The index to start from when naming elements. Defaults to None.
    character_zero_index (str, optional): The starting character for zero index. Defaults to "".

    Returns:
    list: A list of sympy symbols representing the tensor.

    >>> create_tensor_symbolic(2, 'a')
    [a0, a1]

    >>> create_tensor_symbolic([2, 2], 'a')
    [[a00, a01], [a10, a11]]

    >>> create_tensor_symbolic([2, 2], 'a', [3,1,4,1])
    [[a314100, a314101], [a314110, a314111]]

    >>> create_tensor_symbolic([2, 2], 'a', character_zero_index='x')
    [[axx, axy], [ayx, ayy]]

    >>> create_tensor_symbolic([2, 2], 'a', character_zero_index='A')
    [[aAA, aAB], [aBA, aBB]]
    """
    # If shape is a single integer, convert it to a list to represent a 1D tensor
    shape = [shape] if isinstance(shape, int) else shape

    # If index is not provided, initialize it as an empty list. This represents the starting point for naming elements.
    preindex = preindex or []

    # Create an alphabet list for character indexing
    alphabet = list(string.ascii_uppercase + string.ascii_lowercase)

    # Create a mapping of integers to characters
    index_to_char = dict(
        enumerate(alphabet)
    )  # unnecessary comprehension: {i: char for i, char in enumerate(alphabet)}

    # Get the starting index for character_zero_index if provided
    start_idx = alphabet.index(character_zero_index) if character_zero_index else 0

    # Generate the tensor
    tensor = [
        sp.Symbol(
            symbol
            + "".join(
                index_to_char[start_idx + n] if character_zero_index else str(n)
                for n in preindex + [i]
            ),
        )
        if symbol
        else sp.sympify(0)
        for i in range(shape[0])
    ]

    # If the shape has more than one dimension, recursively call create_symbolic_tensor function to handle the remaining dimensions.
    if len(shape) > 1:
        tensor = [
            create_tensor_symbolic(
                shape[1:], symbol, preindex + [i], character_zero_index
            )
            for i in range(shape[0])
        ]

    return tensor


def declare_indexedexp(
    idx_expr_basename: Union[str, None], **kwargs: Any
) -> _recur_symbol_type:
    """
    Generate an indexed expression of specified rank and dimension
       Test 1: convert a symmetric rank-2 tensor to a 1D list & find the number of unique indices.
    >>> ixp = declare_indexedexp('M', rank=2, dimension=3, symmetry='sym01')
    >>> assert func.pipe(ixp, lambda x: func.repeat(func.flatten, x, 1), set, len) == 6

       Test 2: convert a symmetric rank-3 tensor to a 1D list & find the number of unique indices.
    >>> ixp = declare_indexedexp('M', rank=3, dimension=3, symmetry='sym01')
    >>> assert len(set(func.repeat(func.flatten, ixp, 2))) == 18
    >>> ixp = declare_indexedexp('M', rank=3, dimension=3, symmetry='sym02')
    >>> assert len(set(func.repeat(func.flatten, ixp, 2))) == 18
    >>> ixp = declare_indexedexp('M', rank=3, dimension=3, symmetry='sym12')
    >>> assert len(set(func.repeat(func.flatten, ixp, 2))) == 18
    >>> ixp = declare_indexedexp('M', rank=3, dimension=3, symmetry='sym012')
    >>> assert len(set(func.repeat(func.flatten, ixp, 2))) == 10

       Test 3: convert a symmetric rank-4 tensor to a 1D list & find the number of unique indices.
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='sym01')
    >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 54
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='sym02')
    >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 54
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='sym03')
    >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 54
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='sym12')
    >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 54
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='sym13')
    >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 54
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='sym23')
    >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 54
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='sym012')
    >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 30
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='sym013')
    >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 30
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='sym01_sym23')
    >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 36
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='sym02_sym13')
    >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 36
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='sym023')
    >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 30
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='sym03_sym12')
    >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 36
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='sym123')
    >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 30
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='sym0123')
    >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 15

       Test 4: convert an antisymmetric rank-4 tensor to a 1D list & find the number of unique indices.
    >>> ixp = declare_indexedexp('M', rank=2, dimension=3, symmetry='anti01')
    >>> assert len(set(map(abs, func.repeat(func.flatten, ixp, 1))).difference({0})) == 3
    >>> ixp = declare_indexedexp('M', rank=3, dimension=3, symmetry='anti012')
    >>> assert len(set(map(abs, func.repeat(func.flatten, ixp, 2))).difference({0})) == 1
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='anti0123')
    >>> assert len(set(map(abs, func.repeat(func.flatten, ixp, 3))).difference({0})) == 0
    """

    def symmetrize(
        rank: int,
        indexedexp: Sequence[_recur_symbol_type],
        symmetry: str,
        dimension: int,
    ) -> Sequence[_recur_symbol_type]:
        """
        This function symmetrizes the provided tensor of the given rank according to the specified symmetry.

        Args:
        rank (int): The rank of the tensor.
        indexedexp (List[Union[sp.Symbol, sp.Expr, int, float]]): The tensor to symmetrize.
        symmetry (str): The symmetry to apply.
        dimension (int): The dimension of the tensor.

        Returns:
        List[Union[sp.Symbol, sp.Expr, int, float]]: The symmetrized tensor.

        Raises:
        ValueError: If an unsupported rank or symmetry option is provided.
        """
        if rank == 1:
            if symmetry == "nosym":
                return indexedexp
            raise ValueError("cannot symmetrize indexed expression of rank 1")
        if rank == 2:
            indexedexp = symmetrize_rank2(
                cast(_rank2_type, indexedexp), symmetry, dimension
            )
        elif rank == 3:
            indexedexp = symmetrize_rank3(
                cast(_rank3_type, indexedexp), symmetry, dimension
            )
        elif rank == 4:
            indexedexp = symmetrize_rank4(
                cast(_rank4_type, indexedexp), symmetry, dimension
            )
        else:
            raise ValueError("unsupported rank for indexed expression")
        return indexedexp

    if idx_expr_basename is not None:

        def valid_symbol(sym: Any) -> bool:
            if not isinstance(sym, str):
                return False
            for char in sym:
                if not (char.isalpha() or char.isdigit() or char in ["_", "-", ">"]):
                    return False
            return True

        if not valid_symbol(idx_expr_basename):
            raise ValueError(
                'symbol must be an alphanumeric string, underscores and "->" allowed.'
            )
    old_dimension = kwargs.get("DIM")
    if old_dimension:
        raise ValueError("DIM= is deprecated, use dimension=")
    dimension = kwargs.get("dimension")
    if not dimension or dimension < 0:
        dimension = 3  # set default dimension to 3
    elif not isinstance(dimension, int) or dimension <= 0:
        raise ValueError(
            "dimension=N argument must be set, and N must be a positive integer"
        )
    rank = kwargs.get("rank")
    if not rank or not isinstance(rank, int) or not 0 < rank <= 4:
        raise ValueError(
            "rank=N argument must be set, and N must be between 1 and 4 inclusive"
        )
    indexedexp: Sequence[_recur_symbol_type] = create_tensor_symbolic(
        rank * [dimension], idx_expr_basename
    )
    symmetry = kwargs.get("symmetry")
    if symmetry:
        indexedexp = symmetrize(
            rank, cast(_rank2_type, indexedexp), symmetry, dimension
        )
    return zero_out_derivatives_across_symmetry_axes(indexedexp)


def symmetrize_rank2(
    indexedexp: _rank2_type,
    symmetry: str,
    dimension: int,
) -> _rank2_type:
    """
    This function takes a 2D list (matrix), a symmetry option, and a dimension number,
    and returns the matrix after applying the specified symmetry operation.

    Args:
    indexedexp (_rank2_type): The 2D list or matrix on which to perform the operation.
    symmetry (str): The symmetry operation to perform. Options include "sym01" or "nosym".
    dimension (int): The dimension of the matrix.

    Returns:
    _rank2_type: The matrix after applying the symmetry operation.

    Raises:
    ValueError: If an unsupported symmetry option is provided.

    Doctest:
    >>> a, b, c, d = sp.symbols('a b c d', real=True)
    >>> print(symmetrize_rank2([[a, b], [c, d]], "sym01", 2))
    [[a, b], [b, d]]
    """
    for sym in symmetry.split("_"):
        sign = 1 if sym[:3] == "sym" else -1
        for i, j in func.product(range(dimension), repeat=2):
            if sym[-2:] == "01":
                assert isinstance(indexedexp, list) and isinstance(indexedexp[0], list)
                if j < i:
                    indexedexp[i][j] = sign * indexedexp[j][i]
                elif i == j and sign < 0:
                    indexedexp[i][j] = 0
            elif sym == "nosym":
                pass
            else:
                raise ValueError(f"unsupported symmetry option '{str(sym)}'")
    return indexedexp


def symmetrize_rank3(
    indexedexp: _rank3_type,
    symmetry: str,
    dimension: int,
) -> _rank3_type:
    """
    This function takes a 3D list (tensor), a symmetry option, and a dimension number,
    and returns the tensor after applying the specified symmetry operation.

    Args:
    indexedexp (_rank3_type): The 3D list or tensor on which to perform the operation.
    symmetry (str): The symmetry operation to perform. Options include "sym01", "sym02", "sym12", "nosym", and others.
    dimension (int): The dimension of the tensor.

    Returns:
    _rank3_type: The tensor after applying the symmetry operation.

    Raises:
    ValueError: If an unsupported symmetry option is provided.

    Doctest:
    >>> a, b, c, d, e, f, g, h, i0 = sp.symbols('a b c d e f g h i0', real=True)
    >>> rank3 = symmetrize_rank3([[[a, b, c], [d, e, f], [g, h, i0]],
    ...                           [[a, b, c], [d, e, f], [g, h, i0]],
    ...                           [[a, b, c], [d, e, f], [g, h, i0]]], "sym01", 3)
    >>> for i in range(3):
    ...     for j in range(i, 3):
    ...         for k in range(3):
    ...             if rank3[i][j][k] - rank3[j][i][k]:
    ...                 raise ValueError("ERROR! symmetrize_rank3 didn't work.")
    """
    symmetry_list = []
    for sym in symmetry.split("_"):
        index = 3 if sym[:3] == "sym" else 4
        if len(sym[index:]) == 3:
            prefix = sym[:index]
            symmetry_list.append(prefix + sym[index : (index + 2)])
            symmetry_list.append(prefix + sym[(index + 1) : (index + 3)])
        else:
            symmetry_list.append(sym)
    for sym in (
        symmetry_list[k] for n in range(len(symmetry_list), 0, -1) for k in range(n)
    ):
        sign = 1 if sym[:3] == "sym" else -1
        for i, j, k in func.product(range(dimension), repeat=3):
            if sym[-2:] == "01":
                if j < i:
                    indexedexp[i][j][k] = sign * indexedexp[j][i][k]
                elif i == j and sign < 0:
                    indexedexp[i][j][k] = 0
            elif sym[-2:] == "02":
                if k < i:
                    indexedexp[i][j][k] = sign * indexedexp[k][j][i]
                elif i == k and sign < 0:
                    indexedexp[i][j][k] = 0
            elif sym[-2:] == "12":
                if k < j:
                    indexedexp[i][j][k] = sign * indexedexp[i][k][j]
                elif j == k and sign < 0:
                    indexedexp[i][j][k] = 0
            elif sym == "nosym":
                pass
            else:
                raise ValueError(f"unsupported symmetry option '{str(sym)}'")
    return indexedexp


def symmetrize_rank4(
    indexedexp: _rank4_type,
    symmetry: str,
    dimension: int,
) -> _rank4_type:
    """
    This function takes a 4D list (4-tensor), a symmetry option, and a dimension number,
    and returns the 4-tensor after applying the specified symmetry operation.

    Args:
    indexedexp (_rank4_type): The 4D list or 4-tensor on which to perform the operation.
    symmetry (str): The symmetry operation to perform. Options include "sym01", "sym02", "sym12", "sym03", "sym13", "sym23", "nosym", and others.
    dimension (int): The dimension of the 4-tensor.

    Returns:
    _rank4_type: The 4-tensor after applying the symmetry operation.

    Raises:
    ValueError: If an unsupported symmetry option is provided.

    Doctest:
    >>> aDDDD = declarerank4("aDDDD", dimension=3)
    >>> asymDDDD = symmetrize_rank4(aDDDD, symmetry="sym12", dimension=3)
    >>> for i in range(3):
    ...     for j in range(3):
    ...         for k in range(3):
    ...             if asymDDDD[i][j][k][2] - asymDDDD[i][k][j][2]:
    ...                 raise ValueError("ERROR! symmetrize_rank4 didn't work.")

    """
    symmetry_list = []
    for sym in symmetry.split("_"):
        index = 3 if sym[:3] == "sym" else 4
        if len(sym[index:]) in (3, 4):
            prefix = sym[:index]
            symmetry_list.append(prefix + sym[index : (index + 2)])
            symmetry_list.append(prefix + sym[(index + 1) : (index + 3)])
            if len(sym[index:]) == 4:
                symmetry_list.append(prefix + sym[(index + 2) : (index + 4)])
        else:
            symmetry_list.append(sym)
    for sym in (
        symmetry_list[k] for n in range(len(symmetry_list), 0, -1) for k in range(n)
    ):
        sign = 1 if sym[:3] == "sym" else -1
        for i, j, k, l in func.product(range(dimension), repeat=4):
            if sym[-2:] == "01":
                if j < i:
                    indexedexp[i][j][k][l] = sign * indexedexp[j][i][k][l]
                elif i == j and sign < 0:
                    indexedexp[i][j][k][l] = 0
            elif sym[-2:] == "02":
                if k < i:
                    indexedexp[i][j][k][l] = sign * indexedexp[k][j][i][l]
                elif i == k and sign < 0:
                    indexedexp[i][j][k][l] = 0
            elif sym[-2:] == "03":
                if l < i:
                    indexedexp[i][j][k][l] = sign * indexedexp[l][j][k][i]
                elif i == l and sign < 0:
                    indexedexp[i][j][k][l] = 0
            elif sym[-2:] == "12":
                if k < j:
                    indexedexp[i][j][k][l] = sign * indexedexp[i][k][j][l]
                elif j == k and sign < 0:
                    indexedexp[i][j][k][l] = 0
            elif sym[-2:] == "13":
                if l < j:
                    indexedexp[i][j][k][l] = sign * indexedexp[i][l][k][j]
                elif j == l and sign < 0:
                    indexedexp[i][j][k][l] = 0
            elif sym[-2:] == "23":
                if l < k:
                    indexedexp[i][j][k][l] = sign * indexedexp[i][j][l][k]
                elif k == l and sign < 0:
                    indexedexp[i][j][k][l] = 0
            elif sym == "nosym":
                pass
            else:
                raise ValueError(f"unsupported symmetry option '{sym}'")
    return indexedexp


# typehinting: Must allow for sp.Expr in the zerorank*()s, so that it can be modified later
def zerorank1(dimension: int = -1) -> List[Union[sp.Expr, sp.Symbol]]:
    """Initialize rank-1 indexedexp to zero."""
    return declare_indexedexp(None, rank=1, dimension=dimension)  # type:ignore


def zerorank2(dimension: int = -1) -> List[List[Union[sp.Expr, sp.Symbol]]]:
    """Initialize rank-2 indexedexp to zero."""
    return declare_indexedexp(None, rank=2, dimension=dimension)  # type:ignore


def zerorank3(dimension: int = -1) -> List[List[List[Union[sp.Expr, sp.Symbol]]]]:
    """Initialize rank-3 indexedexp to zero."""
    return declare_indexedexp(None, rank=3, dimension=dimension)  # type:ignore


def zerorank4(dimension: int = -1) -> List[List[List[List[Union[sp.Expr, sp.Symbol]]]]]:
    """Initialize rank-4 indexedexp to zero."""
    return declare_indexedexp(None, rank=4, dimension=dimension)  # type:ignore


def get_rank(IDX_EXPR: Sequence[_recur_symbol_type]) -> int:
    """
    Gets the rank of the expression based on the nested list structure.

    :param IDX_EXPR: The expression to determine the rank of.
    :return: The rank of the expression.

    >>> i = sp.sympify(1)
    >>> get_rank([i, i, i])
    1
    >>> get_rank([[i, i], [i, i]])
    2
    >>> get_rank([[[i]]])
    3
    >>> get_rank([[[[[i]]]]])
    5
    """
    rank = 0
    temp = IDX_EXPR

    while isinstance(temp, list):
        rank += 1
        if len(temp) == 0:
            break
        temp = temp[0]

    return rank


def zero_out_derivatives_across_symmetry_axes(
    IDX_EXPR: Sequence[_recur_symbol_type],
) -> _recur_symbol_type:
    """
    Checks if an index object performs a derivative across a symmetry axis.

    :param idxobj_str: The string representing the index object.
    :return: True if the index object performs a derivative across a symmetry axis, False otherwise.
    >>> zero_out_derivatives_across_symmetry_axes(declarerank1("trK_dD"))
    [trK_dD0, trK_dD1, trK_dD2]
    >>> zero_out_derivatives_across_symmetry_axes(declarerank2("trK_dDD"))
    [[trK_dDD00, trK_dDD01, trK_dDD02], [trK_dDD10, trK_dDD11, trK_dDD12], [trK_dDD20, trK_dDD21, trK_dDD22]]
    >>> zero_out_derivatives_across_symmetry_axes(declarerank2("trK_dDD", symmetry="sym01"))
    [[trK_dDD00, trK_dDD01, trK_dDD02], [trK_dDD01, trK_dDD11, trK_dDD12], [trK_dDD02, trK_dDD12, trK_dDD22]]
    >>> par.set_parval_from_str("symmetry_axes", "2")
    >>> zero_out_derivatives_across_symmetry_axes(declarerank1("trK_dD"))
    [trK_dD0, trK_dD1, 0]
    >>> zero_out_derivatives_across_symmetry_axes(declarerank2("trK_dDD"))
    [[trK_dDD00, trK_dDD01, 0], [trK_dDD10, trK_dDD11, 0], [0, 0, 0]]
    >>> zero_out_derivatives_across_symmetry_axes(declarerank2("trK_dDD", symmetry="sym01"))
    [[trK_dDD00, trK_dDD01, 0], [trK_dDD01, trK_dDD11, 0], [0, 0, 0]]
    >>> zero_out_derivatives_across_symmetry_axes(declarerank3("trK_dDDD"))
    Traceback (most recent call last):
    ...
    ValueError: Error. Derivative order > 2 not supported. Failed expression: trK_dDDD000
    >>> par.set_parval_from_str("symmetry_axes", "01")
    >>> zero_out_derivatives_across_symmetry_axes(declarerank2("trK_dDD"))
    [[0, 0, 0], [0, 0, 0], [0, 0, trK_dDD22]]
    >>> par.set_parval_from_str("symmetry_axes", "0")
    >>> zero_out_derivatives_across_symmetry_axes(declarerank4("aDD_dDD", symmetry="sym01_sym23"))
    [[[[0, 0, 0], [0, aDD_dDD0011, aDD_dDD0012], [0, aDD_dDD0012, aDD_dDD0022]], [[0, 0, 0], [0, aDD_dDD0111, aDD_dDD0112], [0, aDD_dDD0112, aDD_dDD0122]], [[0, 0, 0], [0, aDD_dDD0211, aDD_dDD0212], [0, aDD_dDD0212, aDD_dDD0222]]], [[[0, 0, 0], [0, aDD_dDD0111, aDD_dDD0112], [0, aDD_dDD0112, aDD_dDD0122]], [[0, 0, 0], [0, aDD_dDD1111, aDD_dDD1112], [0, aDD_dDD1112, aDD_dDD1122]], [[0, 0, 0], [0, aDD_dDD1211, aDD_dDD1212], [0, aDD_dDD1212, aDD_dDD1222]]], [[[0, 0, 0], [0, aDD_dDD0211, aDD_dDD0212], [0, aDD_dDD0212, aDD_dDD0222]], [[0, 0, 0], [0, aDD_dDD1211, aDD_dDD1212], [0, aDD_dDD1212, aDD_dDD1222]], [[0, 0, 0], [0, aDD_dDD2211, aDD_dDD2212], [0, aDD_dDD2212, aDD_dDD2222]]]]
    """
    symmetry_axes = par.parval_from_str("indexedexp::symmetry_axes")
    if symmetry_axes == "":
        return IDX_EXPR

    rank = get_rank(IDX_EXPR)
    if rank > 4 or rank == 0:
        raise ValueError(
            f"Error: Only ranks from 1 to 4 are supported; {IDX_EXPR} does not seem to fit this pattern."
        )

    def performs_derivative_across_symmetry_axis(idxobj_str: str) -> bool:
        if "_d" in idxobj_str:
            # First we find the order of the derivative:
            deriv_order = 1
            for i in range(len(idxobj_str) - 1):
                if idxobj_str[i : i + 2] == "_d":
                    # The order of the derivative is given by the number of D's in a row after the _d:
                    if "DDD" in idxobj_str[i + 2 :]:
                        raise ValueError(
                            f"Error. Derivative order > 2 not supported. Failed expression: {idxobj_str}"
                        )
                    if "DD" in idxobj_str[i + 2 :]:
                        deriv_order = 2
            end_idx_of_idxobj_str = len(idxobj_str) - 1
            for j in range(
                end_idx_of_idxobj_str, end_idx_of_idxobj_str - deriv_order, -1
            ):
                if idxobj_str[j] in symmetry_axes:
                    return True
        return False

    if rank == 1:
        assert isinstance(IDX_EXPR, list)
        DIM = len(IDX_EXPR)
        for i0 in range(DIM):
            if performs_derivative_across_symmetry_axis(str(IDX_EXPR[i0])):
                IDX_EXPR[i0] = sp.sympify(0)
    if rank == 2:
        assert isinstance(IDX_EXPR, list) and isinstance(IDX_EXPR[0], list)
        DIM = len(IDX_EXPR[0])
        for i0 in range(DIM):
            for i1 in range(DIM):
                if performs_derivative_across_symmetry_axis(str(IDX_EXPR[i0][i1])):
                    IDX_EXPR[i0][i1] = sp.sympify(0)
    if rank == 3:
        assert (
            isinstance(IDX_EXPR, list)
            and isinstance(IDX_EXPR[0], list)
            and isinstance(IDX_EXPR[0][0], list)
        )
        DIM = len(IDX_EXPR[0][0])
        for i0 in range(DIM):
            for i1 in range(DIM):
                for i2 in range(DIM):
                    if performs_derivative_across_symmetry_axis(
                        str(IDX_EXPR[i0][i1][i2])
                    ):
                        IDX_EXPR[i0][i1][i2] = sp.sympify(0)
    if rank == 4:
        assert (
            isinstance(IDX_EXPR, list)
            and isinstance(IDX_EXPR[0], list)
            and isinstance(IDX_EXPR[0][0], list)
            and isinstance(IDX_EXPR[0][0][0], list)
        )
        DIM = len(IDX_EXPR[0][0][0])
        for i0 in range(DIM):
            for i1 in range(DIM):
                for i2 in range(DIM):
                    for i3 in range(DIM):
                        if performs_derivative_across_symmetry_axis(
                            str(IDX_EXPR[i0][i1][i2][i3])
                        ):
                            IDX_EXPR[i0][i1][i2][i3] = sp.sympify(0)
    return IDX_EXPR


def declarerank1(idx_expr_basename: str, **kwargs: Any) -> _rank1_type:
    """Declare rank-1 indexedexp, using idx_expr_basename as base name."""
    kwargs["rank"] = 1
    return cast(_rank1_type, declare_indexedexp(idx_expr_basename, **kwargs))


def declarerank2(idx_expr_basename: str, **kwargs: Any) -> _rank2_type:
    """Declare rank-2 indexedexp, using idx_expr_basename as base name."""
    kwargs["rank"] = 2
    return cast(_rank2_type, declare_indexedexp(idx_expr_basename, **kwargs))


def declarerank3(idx_expr_basename: str, **kwargs: Any) -> _rank3_type:
    """Declare rank-3 indexedexp, using idx_expr_basename as base name."""
    kwargs["rank"] = 3
    return cast(_rank3_type, declare_indexedexp(idx_expr_basename, **kwargs))


def declarerank4(idx_expr_basename: str, **kwargs: Any) -> _rank4_type:
    """Declare rank-4 indexedexp, using idx_expr_basename as base name."""
    kwargs["rank"] = 4
    return cast(_rank4_type, declare_indexedexp(idx_expr_basename, **kwargs))


class NonInvertibleMatrixError(ZeroDivisionError):
    """Matrix Not Invertible; Division By Zero"""


# We use the following functions to evaluate 3-metric inverses
def symm_matrix_inverter2x2(
    a: List[List[Union[sp.Expr, sp.Symbol]]]
) -> Tuple[List[List[Union[sp.Expr, sp.Symbol]]], Union[sp.Expr, sp.Symbol]]:
    """
    Calculate the inverse and determinant of a symmetric 2x2 matrix.

    This function first calculates the determinant of the input matrix. If the determinant is zero,
    it raises a NonInvertibleMatrixError. It then calculates the inverse matrix.

    Args:
        a: A symmetric 2x2 matrix of sympy Expressions or Symbols.

    Returns:
        Tuple containing the inverse matrix and the determinant.

    Raises:
        NonInvertibleMatrixError: If the input matrix has determinant zero.

    Examples:
    >>> from sympy import symbols
    >>> b, c, d = symbols('b c d')
    >>> matrix = [[b, c], [c, d]]
    >>> inverse, determinant = symm_matrix_inverter2x2(matrix)
    """

    # It is far more efficient to write out the matrix determinant and inverse by hand
    #   instead of using SymPy's built-in functions, since the matrix is symmetric.
    outDET = a[0][0] * a[1][1] - a[0][1] ** 2
    if outDET == 0:
        raise NonInvertibleMatrixError("matrix has determinant zero")

    outINV = [[sp.sympify(0) for _ in range(2)] for __ in range(2)]

    # First fill in the upper-triangle of the gPhysINV matrix...
    outINV[0][0] = a[1][1] / outDET
    outINV[0][1] = -a[0][1] / outDET
    outINV[1][1] = a[0][0] / outDET
    outINV[1][0] = outINV[0][1]
    return outINV, outDET


def symm_matrix_inverter3x3(
    a: Sequence[Sequence[Union[sp.Expr, sp.Symbol]]]
) -> Tuple[List[List[Union[sp.Expr, sp.Symbol]]], Union[sp.Expr, sp.Symbol]]:
    """
    Calculate the inverse and determinant of a symmetric 3x3 matrix.

    This function first calculates the determinant of the input matrix. If the determinant is zero,
    it raises a NonInvertibleMatrixError. It then calculates the inverse matrix.

    Args:
        a: A symmetric 3x3 matrix of sympy Expressions or Symbols.

    Returns:
        Tuple containing the inverse matrix and the determinant.

    Raises:
        NonInvertibleMatrixError: If the input matrix has determinant zero.

    Examples:
    >>> from sympy import symbols
    >>> a0, b, c, d, e, f = symbols('a0 b c d e f')
    >>> matrix = [[a0, b, c], [b, d, e], [c, e, f]]
    >>> inverse, determinant = symm_matrix_inverter3x3(matrix)
    """
    # It is far more efficient to write out the matrix determinant and inverse by hand
    #   instead of using SymPy's built-in functions, since the matrix is symmetric.
    outDET = (
        -a[0][2] ** 2 * a[1][1]
        + 2 * a[0][1] * a[0][2] * a[1][2]
        - a[0][0] * a[1][2] ** 2
        - a[0][1] ** 2 * a[2][2]
        + a[0][0] * a[1][1] * a[2][2]
    )
    if outDET == 0:
        # print(a)
        raise NonInvertibleMatrixError("matrix has determinant zero")

    outINV = [[sp.sympify(0) for _ in range(3)] for __ in range(3)]

    # First fill in the upper-triangle of the gPhysINV matrix...
    outINV[0][0] = (-a[1][2] ** 2 + a[1][1] * a[2][2]) / outDET
    outINV[0][1] = (+a[0][2] * a[1][2] - a[0][1] * a[2][2]) / outDET
    outINV[0][2] = (-a[0][2] * a[1][1] + a[0][1] * a[1][2]) / outDET
    outINV[1][1] = (-a[0][2] ** 2 + a[0][0] * a[2][2]) / outDET
    outINV[1][2] = (+a[0][1] * a[0][2] - a[0][0] * a[1][2]) / outDET
    outINV[2][2] = (-a[0][1] ** 2 + a[0][0] * a[1][1]) / outDET
    outINV[1][0] = outINV[0][1]
    outINV[2][0] = outINV[0][2]
    outINV[2][1] = outINV[1][2]
    return outINV, outDET


# Validation test code for symm_matrix_inverter4x4():
# import indexedexp as ixp
# R4DD = ixp.declarerank2("R4DD", "sym01", dimension=4)
#
# # Compute R4DD's inverse:
# R4DDinv, det = ixp.symm_matrix_inverter4x4(R4DD)
#
# # Next matrix multiply: IsUnit = R^{-1} R
# IsUnit = ixp.zerorank2(dimension=4)
# for i in range(4):
#     for j in range(4):
#         for k in range(4):
#             IsUnit[i][j] += R4DDinv[i][k] * R4DD[k][j]
#             # If you'd like to check R R^{-1} instead:
#             # IsUnit[i][j] += R4DD[i][k] * R4DDinv[k][j]
#
# # Next check, is IsUnit == Unit matrix?!
# from UnitTesting.assert_equal import check_zero
# for diag in range(4):
#     print(check_zero(IsUnit[diag][diag]-1))
# for offdiag_i in range(4):
#     for offdiag_j in range(4):
#         if offdiag_i != offdiag_j:
#             print(check_zero(IsUnit[offdiag_i][offdiag_j]))
# # ^^ all should output as True.
def symm_matrix_inverter4x4(
    a: List[List[Union[sp.Expr, sp.Symbol]]]
) -> Tuple[List[List[Union[sp.Expr, sp.Symbol]]], Union[sp.Expr, sp.Symbol]]:
    """
    Function to calculate the inverse and determinant of a 4x4 symmetric matrix.

    Args:
        a (List[List[Union[sp.Expr, sp.Symbol]]]): The 4x4 symmetric matrix to be inverted.
            The matrix elements can be either symbolic expressions or symbols.

    Returns:
        Tuple[List[List[Union[sp.Expr, sp.Symbol]]], Union[sp.Expr, sp.Symbol]]:
            A tuple containing two elements.
            The first element is a 4x4 matrix (list of list) which is the inverse of the input matrix.
            The second element is the determinant of the input matrix.

    Raises:
        NonInvertibleMatrixError: Raised when the matrix has determinant zero and hence, non-invertible.

    # >>> from nrpy.helpers.assert_equal import check_zero
    # >>> R4DD = declarerank2("R4DD", "sym01", dimension=4)
    # >>> R4DDinv, det = symm_matrix_inverter4x4(R4DD)
    # >>> IsUnit = zerorank2(dimension=4)  # Next matrix multiply: IsUnit = R^{-1} R
    # >>> for i in range(4):
    # ...     for j in range(4):
    # ...         for k in range(4):
    # ...             IsUnit[i][j] += R4DDinv[i][k] * R4DD[k][j]
    # >>> for diag in range(4):
    # ...    if not check_zero(IsUnit[diag][diag]-1)):
    # ...        print("Error!")
    # >>> for offdiag_i in range(4):
    # ...    for offdiag_j in range(4):
    # ...        if offdiag_i != offdiag_j:
    # ...            if not check_zero(IsUnit[offdiag_i][offdiag_j]):
    # ...                print("Error!")
    """

    # It is far more efficient to write out the matrix determinant and inverse by hand
    #   instead of using SymPy's built-in functions, since the matrix is symmetric.
    outDET = (
        +a[0][2] * a[0][2] * a[1][3] * a[1][3]
        + a[0][3] * a[0][3] * a[1][2] * a[1][2]
        + a[0][1] * a[0][1] * a[2][3] * a[2][3]
        - a[0][0] * a[1][3] * a[1][3] * a[2][2]
        - a[0][3] * a[0][3] * a[1][1] * a[2][2]
        - a[0][0] * a[1][1] * a[2][3] * a[2][3]
        - 2
        * (
            +a[0][1] * a[0][2] * a[1][3] * a[2][3]
            - a[0][0] * a[1][2] * a[1][3] * a[2][3]
            - a[0][3]
            * (
                -a[0][2] * a[1][2] * a[1][3]
                + a[0][1] * a[1][3] * a[2][2]
                + a[0][2] * a[1][1] * a[2][3]
                - a[0][1] * a[1][2] * a[2][3]
            )
        )
        - a[3][3]
        * (
            +a[0][2] * a[0][2] * a[1][1]
            - a[0][1] * a[0][2] * a[1][2]
            - a[0][1] * a[0][2] * a[1][2]
            + a[0][0] * a[1][2] * a[1][2]
            + a[0][1] * a[0][1] * a[2][2]
            - a[0][0] * a[1][1] * a[2][2]
        )
    )
    if outDET == 0:
        raise NonInvertibleMatrixError("matrix has determinant zero")

    outINV = [[sp.sympify(0) for _ in range(4)] for __ in range(4)]

    # First fill in the upper-triangle of the gPhysINV matrix...
    outINV[0][0] = (
        -a[1][3] * a[1][3] * a[2][2]
        + 2 * a[1][2] * a[1][3] * a[2][3]
        - a[1][1] * a[2][3] * a[2][3]
        - a[1][2] * a[1][2] * a[3][3]
        + a[1][1] * a[2][2] * a[3][3]
    ) / outDET
    outINV[1][1] = (
        -a[0][3] * a[0][3] * a[2][2]
        + 2 * a[0][2] * a[0][3] * a[2][3]
        - a[0][0] * a[2][3] * a[2][3]
        - a[0][2] * a[0][2] * a[3][3]
        + a[0][0] * a[2][2] * a[3][3]
    ) / outDET
    outINV[2][2] = (
        -a[0][3] * a[0][3] * a[1][1]
        + 2 * a[0][1] * a[0][3] * a[1][3]
        - a[0][0] * a[1][3] * a[1][3]
        - a[0][1] * a[0][1] * a[3][3]
        + a[0][0] * a[1][1] * a[3][3]
    ) / outDET
    outINV[3][3] = (
        -a[0][2] * a[0][2] * a[1][1]
        + 2 * a[0][1] * a[0][2] * a[1][2]
        - a[0][0] * a[1][2] * a[1][2]
        - a[0][1] * a[0][1] * a[2][2]
        + a[0][0] * a[1][1] * a[2][2]
    ) / outDET
    outINV[0][1] = (
        +a[0][3] * a[1][3] * a[2][2]
        - a[0][3] * a[1][2] * a[2][3]
        - a[0][2] * a[1][3] * a[2][3]
        + a[0][1] * a[2][3] * a[2][3]
        + a[0][2] * a[1][2] * a[3][3]
        - a[0][1] * a[2][2] * a[3][3]
    ) / outDET
    outINV[0][2] = (
        -a[0][3] * a[1][2] * a[1][3]
        + a[0][2] * a[1][3] * a[1][3]
        + a[0][3] * a[1][1] * a[2][3]
        - a[0][1] * a[1][3] * a[2][3]
        - a[0][2] * a[1][1] * a[3][3]
        + a[0][1] * a[1][2] * a[3][3]
    ) / outDET
    outINV[0][3] = (
        -a[0][2] * a[1][2] * a[1][3]
        + a[0][1] * a[1][3] * a[2][2]
        + a[0][3] * a[1][2] * a[1][2]
        - a[0][3] * a[1][1] * a[2][2]
        + a[0][2] * a[1][1] * a[2][3]
        - a[0][1] * a[1][2] * a[2][3]
    ) / outDET
    outINV[1][2] = (
        +a[0][3] * a[0][3] * a[1][2]
        + a[0][0] * a[1][3] * a[2][3]
        - a[0][3] * a[0][2] * a[1][3]
        - a[0][3] * a[0][1] * a[2][3]
        + a[0][1] * a[0][2] * a[3][3]
        - a[0][0] * a[1][2] * a[3][3]
    ) / outDET
    outINV[1][3] = (
        +a[0][2] * a[0][2] * a[1][3]
        + a[0][1] * a[0][3] * a[2][2]
        - a[0][0] * a[1][3] * a[2][2]
        + a[0][0] * a[1][2] * a[2][3]
        - a[0][2] * a[0][3] * a[1][2]
        - a[0][2] * a[0][1] * a[2][3]
    ) / outDET
    outINV[2][3] = (
        +a[0][2] * a[0][3] * a[1][1]
        - a[0][1] * a[0][3] * a[1][2]
        - a[0][1] * a[0][2] * a[1][3]
        + a[0][0] * a[1][2] * a[1][3]
        + a[0][1] * a[0][1] * a[2][3]
        - a[0][0] * a[1][1] * a[2][3]
    ) / outDET

    # Then we fill the lower triangle of the symmetric matrix
    outINV[1][0] = outINV[0][1]
    outINV[2][0] = outINV[0][2]
    outINV[2][1] = outINV[1][2]
    outINV[3][0] = outINV[0][3]
    outINV[3][1] = outINV[1][3]
    outINV[3][2] = outINV[2][3]

    return outINV, outDET


# SymPy's generic matrix inverter takes a long time to invert 3x3 matrices, so here we have an optimized version.
# We use the following functions to evaluate 3-metric inverses
def generic_matrix_inverter2x2(
    a: List[List[Union[sp.Expr, sp.Symbol]]]
) -> Tuple[List[List[Union[sp.Expr, sp.Symbol]]], Union[sp.Expr, sp.Symbol]]:
    """
    Calculate the inverse and determinant of a general 2x2 matrix.

    This function first calculates the determinant of the input matrix. If the determinant is zero,
    it raises a NonInvertibleMatrixError. It then calculates the inverse matrix.

    Args:
        a: A 2x2 matrix of sympy Expressions or Symbols.

    Returns:
        Tuple containing the inverse matrix and the determinant.

    Raises:
        NonInvertibleMatrixError: If the input matrix has determinant zero.

    Examples:
    >>> matrix = declarerank2("gDD", dimension=2)
    >>> gUU, detg = generic_matrix_inverter2x2(matrix)
    """
    outDET = a[0][0] * a[1][1] - a[0][1] * a[1][0]
    if outDET == 0:
        raise NonInvertibleMatrixError("matrix has determinant zero")

    outINV = [[sp.sympify(0) for _ in range(2)] for __ in range(2)]

    outINV[0][0] = a[1][1] / outDET
    outINV[0][1] = -a[0][1] / outDET
    outINV[1][1] = a[0][0] / outDET
    outINV[1][0] = -a[1][0] / outDET
    return outINV, outDET


def generic_matrix_inverter3x3(
    a: List[List[Union[sp.Expr, sp.Symbol]]]
) -> Tuple[List[List[Union[sp.Expr, sp.Symbol]]], Union[sp.Expr, sp.Symbol]]:
    """
    Calculate the inverse and determinant of a general 3x3 matrix.

    This function first calculates the determinant of the input matrix. If the determinant is zero,
    it raises a NonInvertibleMatrixError. It then calculates the inverse matrix.

    Args:
        a: A 3x3 matrix of sympy Expressions or Symbols.

    Returns:
        Tuple containing the inverse matrix and the determinant.

    Raises:
        NonInvertibleMatrixError: If the input matrix has determinant zero.

    Examples:
    >>> gDD = declarerank2("gDD")
    >>> gUU, detg = generic_matrix_inverter3x3(gDD)
    """
    outDET = (
        -a[0][2] * a[1][1] * a[2][0]
        + a[0][1] * a[1][2] * a[2][0]
        + a[0][2] * a[1][0] * a[2][1]
        - a[0][0] * a[1][2] * a[2][1]
        - a[0][1] * a[1][0] * a[2][2]
        + a[0][0] * a[1][1] * a[2][2]
    )
    if outDET == 0:
        raise NonInvertibleMatrixError("matrix has determinant zero")

    outINV = [[sp.sympify(0) for _ in range(3)] for __ in range(3)]

    outINV[0][0] = -a[1][2] * a[2][1] + a[1][1] * a[2][2]
    outINV[0][1] = a[0][2] * a[2][1] - a[0][1] * a[2][2]
    outINV[0][2] = -a[0][2] * a[1][1] + a[0][1] * a[1][2]
    outINV[1][0] = a[1][2] * a[2][0] - a[1][0] * a[2][2]
    outINV[1][1] = -a[0][2] * a[2][0] + a[0][0] * a[2][2]
    outINV[1][2] = a[0][2] * a[1][0] - a[0][0] * a[1][2]
    outINV[2][0] = -a[1][1] * a[2][0] + a[1][0] * a[2][1]
    outINV[2][1] = a[0][1] * a[2][0] - a[0][0] * a[2][1]
    outINV[2][2] = -a[0][1] * a[1][0] + a[0][0] * a[1][1]

    for i in range(3):
        for j in range(3):
            outINV[i][j] /= outDET

    return outINV, outDET


def generic_matrix_inverter4x4(
    a: List[List[Union[sp.Expr, sp.Symbol]]]
) -> Tuple[List[List[Union[sp.Expr, sp.Symbol]]], Union[sp.Expr, sp.Symbol]]:
    """
    Calculate the inverse and determinant of a general 4x4 matrix.

    This function first calculates the determinant of the input matrix. If the determinant is zero,
    it raises a NonInvertibleMatrixError. It then calculates the inverse matrix.

    Args:
        a: A 3x3 matrix of sympy Expressions or Symbols.

    Returns:
        Tuple containing the inverse matrix and the determinant.

    Raises:
        NonInvertibleMatrixError: If the input matrix has determinant zero.

    Examples:
    >>> g4DD = declarerank2("g4DD", dimension=4)
    >>> g4UU, det4g = generic_matrix_inverter4x4(g4DD)
    """
    # A = {{a00, a01, a02, a03},
    #      {a10, a11, a12, a13},
    #      {a20, a21, a22, a23},
    #      {a30, a31, a32, a33}}
    # A // MatrixForm
    # CForm[FullSimplify[Det[A]]] >>> t2.txt
    # cat t2.txt | sed "s/ //g" |sed "s/ //g;s/\([0-3]\)/[\1]/g"
    outDET = (
        a[0][1] * a[1][3] * a[2][2] * a[3][0]
        - a[0][1] * a[1][2] * a[2][3] * a[3][0]
        - a[0][0] * a[1][3] * a[2][2] * a[3][1]
        + a[0][0] * a[1][2] * a[2][3] * a[3][1]
        - a[0][1] * a[1][3] * a[2][0] * a[3][2]
        + a[0][0] * a[1][3] * a[2][1] * a[3][2]
        + a[0][1] * a[1][0] * a[2][3] * a[3][2]
        - a[0][0] * a[1][1] * a[2][3] * a[3][2]
        + a[0][3]
        * (
            a[1][2] * a[2][1] * a[3][0]
            - a[1][1] * a[2][2] * a[3][0]
            - a[1][2] * a[2][0] * a[3][1]
            + a[1][0] * a[2][2] * a[3][1]
            + a[1][1] * a[2][0] * a[3][2]
            - a[1][0] * a[2][1] * a[3][2]
        )
        + (
            a[0][1] * a[1][2] * a[2][0]
            - a[0][0] * a[1][2] * a[2][1]
            - a[0][1] * a[1][0] * a[2][2]
            + a[0][0] * a[1][1] * a[2][2]
        )
        * a[3][3]
        + a[0][2]
        * (
            -(a[1][3] * a[2][1] * a[3][0])
            + a[1][1] * a[2][3] * a[3][0]
            + a[1][3] * a[2][0] * a[3][1]
            - a[1][0] * a[2][3] * a[3][1]
            - a[1][1] * a[2][0] * a[3][3]
            + a[1][0] * a[2][1] * a[3][3]
        )
    )
    if outDET == 0:
        raise NonInvertibleMatrixError("matrix has determinant zero")

    outINV = [[sp.sympify(0) for _ in range(4)] for __ in range(4)]

    # CForm[FullSimplify[Inverse[A]*Det[A]]] >>> t.txt
    # cat t.txt | sed "s/,/\n/g;s/List(//g;s/))/)/g;s/)//g;s/(//g"|grep -v ^$|sed "s/ //g;s/\([0-3]\)/[\1]/g"| awk '{line[NR]=$0}END{count=1;for(i=0;i<4;i++) { for(j=0;j<4;j++) { printf "outINV[%d][%d] = %s\n", i,j,line[count];count++; }}}'
    outINV[0][0] = (
        -a[1][3] * a[2][2] * a[3][1]
        + a[1][2] * a[2][3] * a[3][1]
        + a[1][3] * a[2][1] * a[3][2]
        - a[1][1] * a[2][3] * a[3][2]
        - a[1][2] * a[2][1] * a[3][3]
        + a[1][1] * a[2][2] * a[3][3]
    )
    outINV[0][1] = (
        a[0][3] * a[2][2] * a[3][1]
        - a[0][2] * a[2][3] * a[3][1]
        - a[0][3] * a[2][1] * a[3][2]
        + a[0][1] * a[2][3] * a[3][2]
        + a[0][2] * a[2][1] * a[3][3]
        - a[0][1] * a[2][2] * a[3][3]
    )
    outINV[0][2] = (
        -a[0][3] * a[1][2] * a[3][1]
        + a[0][2] * a[1][3] * a[3][1]
        + a[0][3] * a[1][1] * a[3][2]
        - a[0][1] * a[1][3] * a[3][2]
        - a[0][2] * a[1][1] * a[3][3]
        + a[0][1] * a[1][2] * a[3][3]
    )
    outINV[0][3] = (
        a[0][3] * a[1][2] * a[2][1]
        - a[0][2] * a[1][3] * a[2][1]
        - a[0][3] * a[1][1] * a[2][2]
        + a[0][1] * a[1][3] * a[2][2]
        + a[0][2] * a[1][1] * a[2][3]
        - a[0][1] * a[1][2] * a[2][3]
    )
    outINV[1][0] = (
        a[1][3] * a[2][2] * a[3][0]
        - a[1][2] * a[2][3] * a[3][0]
        - a[1][3] * a[2][0] * a[3][2]
        + a[1][0] * a[2][3] * a[3][2]
        + a[1][2] * a[2][0] * a[3][3]
        - a[1][0] * a[2][2] * a[3][3]
    )
    outINV[1][1] = (
        -a[0][3] * a[2][2] * a[3][0]
        + a[0][2] * a[2][3] * a[3][0]
        + a[0][3] * a[2][0] * a[3][2]
        - a[0][0] * a[2][3] * a[3][2]
        - a[0][2] * a[2][0] * a[3][3]
        + a[0][0] * a[2][2] * a[3][3]
    )
    outINV[1][2] = (
        a[0][3] * a[1][2] * a[3][0]
        - a[0][2] * a[1][3] * a[3][0]
        - a[0][3] * a[1][0] * a[3][2]
        + a[0][0] * a[1][3] * a[3][2]
        + a[0][2] * a[1][0] * a[3][3]
        - a[0][0] * a[1][2] * a[3][3]
    )
    outINV[1][3] = (
        -a[0][3] * a[1][2] * a[2][0]
        + a[0][2] * a[1][3] * a[2][0]
        + a[0][3] * a[1][0] * a[2][2]
        - a[0][0] * a[1][3] * a[2][2]
        - a[0][2] * a[1][0] * a[2][3]
        + a[0][0] * a[1][2] * a[2][3]
    )
    outINV[2][0] = (
        -a[1][3] * a[2][1] * a[3][0]
        + a[1][1] * a[2][3] * a[3][0]
        + a[1][3] * a[2][0] * a[3][1]
        - a[1][0] * a[2][3] * a[3][1]
        - a[1][1] * a[2][0] * a[3][3]
        + a[1][0] * a[2][1] * a[3][3]
    )
    outINV[2][1] = (
        a[0][3] * a[2][1] * a[3][0]
        - a[0][1] * a[2][3] * a[3][0]
        - a[0][3] * a[2][0] * a[3][1]
        + a[0][0] * a[2][3] * a[3][1]
        + a[0][1] * a[2][0] * a[3][3]
        - a[0][0] * a[2][1] * a[3][3]
    )
    outINV[2][2] = (
        -a[0][3] * a[1][1] * a[3][0]
        + a[0][1] * a[1][3] * a[3][0]
        + a[0][3] * a[1][0] * a[3][1]
        - a[0][0] * a[1][3] * a[3][1]
        - a[0][1] * a[1][0] * a[3][3]
        + a[0][0] * a[1][1] * a[3][3]
    )
    outINV[2][3] = (
        a[0][3] * a[1][1] * a[2][0]
        - a[0][1] * a[1][3] * a[2][0]
        - a[0][3] * a[1][0] * a[2][1]
        + a[0][0] * a[1][3] * a[2][1]
        + a[0][1] * a[1][0] * a[2][3]
        - a[0][0] * a[1][1] * a[2][3]
    )
    outINV[3][0] = (
        a[1][2] * a[2][1] * a[3][0]
        - a[1][1] * a[2][2] * a[3][0]
        - a[1][2] * a[2][0] * a[3][1]
        + a[1][0] * a[2][2] * a[3][1]
        + a[1][1] * a[2][0] * a[3][2]
        - a[1][0] * a[2][1] * a[3][2]
    )
    outINV[3][1] = (
        -a[0][2] * a[2][1] * a[3][0]
        + a[0][1] * a[2][2] * a[3][0]
        + a[0][2] * a[2][0] * a[3][1]
        - a[0][0] * a[2][2] * a[3][1]
        - a[0][1] * a[2][0] * a[3][2]
        + a[0][0] * a[2][1] * a[3][2]
    )
    outINV[3][2] = (
        a[0][2] * a[1][1] * a[3][0]
        - a[0][1] * a[1][2] * a[3][0]
        - a[0][2] * a[1][0] * a[3][1]
        + a[0][0] * a[1][2] * a[3][1]
        + a[0][1] * a[1][0] * a[3][2]
        - a[0][0] * a[1][1] * a[3][2]
    )
    outINV[3][3] = (
        -a[0][2] * a[1][1] * a[2][0]
        + a[0][1] * a[1][2] * a[2][0]
        + a[0][2] * a[1][0] * a[2][1]
        - a[0][0] * a[1][2] * a[2][1]
        - a[0][1] * a[1][0] * a[2][2]
        + a[0][0] * a[1][1] * a[2][2]
    )

    for mu in range(4):
        for nu in range(4):
            outINV[mu][nu] /= outDET

    return outINV, outDET


# Define the rank-3 version of the Levi-Civita symbol.
def LeviCivitaSymbol_dim3_rank3() -> List[List[List[int]]]:
    """
    Calculate the Levi-Civita symbol for 3 dimensions.

    This function creates a 3x3x3 rank-3 tensor where each element corresponds to the Levi-Civita symbol
    epsilon_ijk, calculated using the formula (i - j) * (j - k) * (k - i) * 1/2.

    Returns:
        A 3x3x3 tensor with integer elements, representing the Levi-Civita symbol.

    Example:
    >>> LeviCivitaSymbol_dim3_rank3()
    [[[0, 0, 0], [0, 0, 1], [0, -1, 0]], [[0, 0, -1], [0, 0, 0], [1, 0, 0]], [[0, 1, 0], [-1, 0, 0], [0, 0, 0]]]
    """
    # A MyPy-friendly version of zerorank3(dimension=3) for this case:
    LeviCivitaSymbol = [
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
    ]

    for i in range(3):
        for j in range(3):
            for k in range(3):
                # From https://codegolf.stackexchange.com/questions/160359/levi-civita-symbol :
                LeviCivitaSymbol[i][j][k] = (
                    (i - j) * (j - k) * (k - i) * sp.Rational(1, 2)
                )
    return LeviCivitaSymbol


# Define the UUU rank-3 version of the Levi-Civita *tensor*; UUU divides by sqrtgammaDET
def LeviCivitaTensorUUU_dim3_rank3(
    sqrtgammaDET: sp.Expr,
) -> Union[List[List[List[sp.Expr]]], List[List[List[int]]]]:
    """
    Calculate the Levi-Civita tensor with upper indices for 3 dimensions.

    This function creates a 3x3x3 rank-3 tensor where each element corresponds to the Levi-Civita symbol
    epsilon^ijk, calculated by dividing each Levi-Civita symbol by sqrtgammaDET.

    Args:
        sqrtgammaDET: A sympy expression that is used to divide each element of the Levi-Civita symbol.

    Returns:
        A 3x3x3 tensor with elements as sympy expressions, representing the Levi-Civita tensor.

    Example:
    >>> sqrtgammaDET = sp.Symbol("sqrtgammaDET")
    >>> LeviCivitaTensorUUU_dim3_rank3(sqrtgammaDET)
    [[[0, 0, 0], [0, 0, 1/sqrtgammaDET], [0, -1/sqrtgammaDET, 0]], [[0, 0, -1/sqrtgammaDET], [0, 0, 0], [1/sqrtgammaDET, 0, 0]], [[0, 1/sqrtgammaDET, 0], [-1/sqrtgammaDET, 0, 0], [0, 0, 0]]]
    """
    # Here, we import the Levi-Civita tensor and compute the tensor with upper indices
    LeviCivitaSymbolDDD = LeviCivitaSymbol_dim3_rank3()
    # A MyPy-friendly version of zerorank3(dimension=3) for this case:
    LeviCivitaTensorUUU = [
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
    ]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                LeviCivitaTensorUUU[i][j][k] = (
                    LeviCivitaSymbolDDD[i][j][k] / sqrtgammaDET
                )
    return LeviCivitaTensorUUU


# Define the DDD rank-3 version of the Levi-Civita *tensor*; DDD multiplies by sqrtgammaDET
def LeviCivitaTensorDDD_dim3_rank3(
    sqrtgammaDET: sp.Expr,
) -> Union[List[List[List[sp.Expr]]], List[List[List[int]]]]:
    """
    Calculate the Levi-Civita tensor with lower indices for 3 dimensions.

    This function creates a 3x3x3 rank-3 tensor where each element corresponds to the Levi-Civita symbol
    epsilon_ijk, calculated by multiplying each Levi-Civita symbol by sqrtgammaDET.

    Args:
        sqrtgammaDET: A sympy expression that is used to multiply each element of the Levi-Civita symbol.

    Returns:
        A 3x3x3 tensor with elements as sympy expressions, representing the Levi-Civita tensor.

    Example:
    >>> sqrtgammaDET = sp.Symbol("sqrtgammaDET")
    >>> LeviCivitaTensorDDD_dim3_rank3(sqrtgammaDET)
    [[[0, 0, 0], [0, 0, sqrtgammaDET], [0, -sqrtgammaDET, 0]], [[0, 0, -sqrtgammaDET], [0, 0, 0], [sqrtgammaDET, 0, 0]], [[0, sqrtgammaDET, 0], [-sqrtgammaDET, 0, 0], [0, 0, 0]]]
    """
    # Here, we import the Levi-Civita tensor and compute the tensor with lower indices
    LeviCivitaSymbolDDD = LeviCivitaSymbol_dim3_rank3()
    # A MyPy-friendly version of zerorank3(dimension=3) for this case:
    LeviCivitaTensorDDD = [
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
    ]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                LeviCivitaTensorDDD[i][j][k] = (
                    LeviCivitaSymbolDDD[i][j][k] * sqrtgammaDET
                )
    return LeviCivitaTensorDDD


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
