"""
indexedexp.py: Functions related to indexed expressions including e.g. tensors and pseudotensors.

Authors: Zachariah Etienne, Kenneth Sible, Steven Brandt
"""

# Step 1: Load needed modules
import string
import sys  # Standard Python module for multiplatform OS-level functions
from typing import Any, List, Optional, Sequence, Tuple, Union, cast

import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends

import nrpy.helpers.functional as func  # NRPy+: Python toolkit for functional programming
import nrpy.params as par

# Used to specify axes of symmetry so that derivatives across these axes get set to zero
par.register_param(py_type=str, module=__name__, name="symmetry_axes", value="")

# Define common type hints for tensor structures
_symbol_type = Any  # sp.Expr
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
    Create a multi-dimensional symbolic tensor of a specified shape.

    :param shape: The shape of the tensor. A list of integers, one for each dimension's size.
    :param symbol: The base symbol name for elements in the tensor. Defaults to "".
    :param preindex: A list of parent indices for recursive construction. Defaults to None.
    :param character_zero_index: The starting character for zero index (e.g., 'x' or 'A'). Defaults to "".
    :return: A list of sympy symbols representing the tensor.

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
    # Ensure shape is a list and get current dimension size and remaining shape.
    shape_list = [shape] if isinstance(shape, int) else shape
    current_dim_size = shape_list[0]
    remaining_shape = shape_list[1:]
    indices_prefix = preindex if preindex is not None else []

    # Base case: If no dimensions remain, create the final (rank-1) slice of the tensor.
    if not remaining_shape:
        # Set up for character-based indexing if requested
        if character_zero_index:
            alphabet = string.ascii_uppercase + string.ascii_lowercase
            index_to_char = dict(enumerate(alphabet))
            start_idx = alphabet.index(character_zero_index)

        tensor_slice = []
        for i in range(current_dim_size):
            final_indices = indices_prefix + [i]
            if symbol is None:
                # If no symbol name is provided, create a tensor of zeros.
                element = sp.sympify(0)
            else:
                if character_zero_index:
                    idx_str = "".join(
                        index_to_char[start_idx + n] for n in final_indices
                    )
                else:
                    idx_str = "".join(map(str, final_indices))
                element = sp.Symbol(symbol + idx_str)
            tensor_slice.append(element)
        return tensor_slice

    # Recursive step: For higher-rank tensors, build each slice of the current dimension.
    return [
        create_tensor_symbolic(
            remaining_shape, symbol, indices_prefix + [i], character_zero_index
        )
        for i in range(current_dim_size)
    ]


def _is_valid_symbol_name(sym: Any) -> bool:
    """
    Validate that the symbol name is a string with allowed characters.

    :param sym: The symbol name to validate.
    :return: True if the symbol name is valid, False otherwise.
    """
    if not isinstance(sym, str):
        return False
    # Allow letters, numbers, underscore, and "->" for derivative notation.
    return all(c.isalnum() or c in ["_", "-", ">"] for c in sym)


def _get_tensor_element(
    tensor: Sequence[_recur_symbol_type], indices: Tuple[int, ...]
) -> _symbol_type:
    """
    Get an element from a nested list tensor using a tuple of indices.

    :param tensor: The tensor (nested list) to access.
    :param indices: A tuple of indices specifying the element's position.
    :return: The tensor element at the specified indices.
    """
    element = tensor
    for index in indices:
        element = element[index]
    return element


def _set_tensor_element(
    tensor: Sequence[_recur_symbol_type], indices: Tuple[int, ...], value: _symbol_type
) -> None:
    """
    Set an element in a nested list tensor using a tuple of indices.

    :param tensor: The tensor (nested list) to modify.
    :param indices: A tuple of indices specifying the element's position.
    :param value: The new value to set at the specified position.
    """
    parent = tensor
    for index in indices[:-1]:
        parent = parent[index]
    cast(List[_recur_symbol_type], parent)[indices[-1]] = value


def _expand_symmetry_groups(symmetry: str) -> List[str]:
    """
    Expand compact symmetry notations into a list of pairwise symmetries.

    This function preserves the exact, sometimes non-obvious, behavior of the
    original code to ensure functional equivalence.
    For example: "sym01_sym23" -> ["sym01", "sym23"]
                 "sym012" -> ["sym01", "sym12"]
                 "sym0123" -> ["sym01", "sym12", "sym23"]

    :param symmetry: The compact symmetry string (e.g., "sym012").
    :return: A list of pairwise symmetry strings (e.g., ["sym01", "sym12"]).
    """
    symmetry_list = []
    for sym_group in symmetry.split("_"):
        prefix_len = 3 if sym_group.startswith("sym") else 4
        indices_str = sym_group[prefix_len:]
        # This logic for expanding groups like 'sym012' is preserved from the original.
        # It produces a chain of adjacent swaps (e.g., 0<->1, 1<->2), which, when
        # applied transitively, achieves the intended symmetry.
        if len(indices_str) in (3, 4):
            prefix = sym_group[:prefix_len]
            for i in range(len(indices_str) - 1):
                symmetry_list.append(prefix + indices_str[i : i + 2])
        else:
            symmetry_list.append(sym_group)
    return symmetry_list


def _apply_symmetrization_globally(
    indexedexp: Sequence[_recur_symbol_type],
    rank: int,
    symmetry: str,
    dimension: int,
) -> Sequence[_recur_symbol_type]:
    """
    Apply symmetry properties to a tensor of any rank. This function modifies the tensor in-place.

    :param indexedexp: The tensor to symmetrize, modified in-place.
    :param rank: The rank of the tensor.
    :param symmetry: The symmetry string to apply.
    :param dimension: The dimension of the tensor.
    :raises ValueError: If an unsupported symmetry option is provided.
    :return: The symmetrized tensor.
    """
    symmetry_list = _expand_symmetry_groups(symmetry)
    sym_map = {f"{i}{j}": (i, j) for i in range(rank) for j in range(i + 1, rank)}

    # This nested loop structure is a direct, more readable translation of the original's
    # generator expression. It acts as a fixed-point iteration hack, repeatedly applying
    # symmetries to ensure transitive relations (e.g., a=b and b=c implies a=c) are
    # fully resolved. It must be preserved for functional equivalence.
    for n in range(len(symmetry_list), 0, -1):
        for k in range(n):
            sym = symmetry_list[k]
            sign = 1 if sym.startswith("sym") else -1
            sym_indices_str = sym[-2:]

            if sym_indices_str in sym_map:
                idx1, idx2 = sym_map[sym_indices_str]
                for indices in func.product(range(dimension), repeat=rank):
                    if indices[idx1] > indices[idx2]:
                        # Swap indices to find the canonical element
                        swapped_indices = list(indices)
                        swapped_indices[idx1], swapped_indices[idx2] = (
                            swapped_indices[idx2],
                            swapped_indices[idx1],
                        )
                        val = _get_tensor_element(indexedexp, tuple(swapped_indices))
                        _set_tensor_element(indexedexp, indices, sign * val)
                    elif sign == -1 and indices[idx1] == indices[idx2]:
                        # Zero out diagonal elements for anti-symmetry
                        _set_tensor_element(indexedexp, indices, sp.sympify(0))

            elif sym != "nosym":
                raise ValueError(f"unsupported symmetry option '{sym}'")

    return indexedexp


def declare_indexedexp(
    idx_expr_basename: Union[str, None], **kwargs: Any
) -> _recur_symbol_type:
    """
    Generate an indexed expression of specified rank and dimension.

    This function creates a symbolic representation for tensors or indexed expressions
    based on provided specifications like rank, dimension, and symmetry.

    :param idx_expr_basename: Base name for the indexed expression. If None,
                              a symbolic representation is not generated.
    :param kwargs: Properties of the indexed expression such as 'rank', 'dimension', and 'symmetry'.
    :return: An indexed expression or tensor as specified.

    :raises ValueError: If 'symmetry' is not a recognized pattern, if 'dimension' or 'rank'
                        are not correctly specified, or if 'idx_expr_basename' is invalid.

    Doctest 1: convert a symmetric rank-2 tensor to a 1D list & find the number of unique indices.
    >>> ixp = declare_indexedexp('M', rank=2, dimension=3, symmetry='sym01')
    >>> assert func.pipe(ixp, lambda x: func.repeat(func.flatten, x, 1), set, len) == 6

    Doctest 2: Attempt to create a tensor with invalid symmetry.
    >>> try:
    ...     Merror = declare_indexedexp('M', rank=2, dimension=3, symmetry='01')
    ... except ValueError as e:
    ...     assert str(e) == "Unsupported symmetry '01' for indexed expression. Valid symmetry options must start with 'sym', 'anti', or 'nosym'"

    Doctest 3: convert a symmetric rank-3 tensor to a 1D list & find the number of unique indices.
    >>> ixp = declare_indexedexp('M', rank=3, dimension=3, symmetry='sym01')
    >>> assert len(set(func.repeat(func.flatten, ixp, 2))) == 18
    >>> ixp = declare_indexedexp('M', rank=3, dimension=3, symmetry='sym02')
    >>> assert len(set(func.repeat(func.flatten, ixp, 2))) == 18
    >>> ixp = declare_indexedexp('M', rank=3, dimension=3, symmetry='sym12')
    >>> assert len(set(func.repeat(func.flatten, ixp, 2))) == 18
    >>> ixp = declare_indexedexp('M', rank=3, dimension=3, symmetry='sym012')
    >>> assert len(set(func.repeat(func.flatten, ixp, 2))) == 10

    Doctest 4: convert a symmetric rank-4 tensor to a 1D list & find the number of unique indices.
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

    Doctest 5: convert an antisymmetric rank-4 tensor to a 1D list & find the number of unique indices.
    >>> ixp = declare_indexedexp('M', rank=2, dimension=3, symmetry='anti01')
    >>> assert len(set(map(abs, func.repeat(func.flatten, ixp, 1))).difference({0})) == 3
    >>> ixp = declare_indexedexp('M', rank=3, dimension=3, symmetry='anti012')
    >>> assert len(set(map(abs, func.repeat(func.flatten, ixp, 2))).difference({0})) == 1
    >>> ixp = declare_indexedexp('M', rank=4, dimension=3, symmetry='anti0123')
    >>> assert len(set(map(abs, func.repeat(func.flatten, ixp, 3))).difference({0})) == 0
    """
    if idx_expr_basename is not None and not _is_valid_symbol_name(idx_expr_basename):
        raise ValueError(
            'symbol must be an alphanumeric string, with underscores and "->" allowed.'
        )

    # Deprecation check for old 'DIM' parameter
    if "DIM" in kwargs:
        raise ValueError("DIM= is deprecated, use dimension=")

    # FIX: The following logic is restored from the original code to ensure
    # functional equivalence. The zerorank*() functions call this with
    # dimension=-1 to signal that the default dimension should be used.
    dimension = kwargs.get("dimension")
    if dimension is None or dimension < 0:
        dimension = 3  # set default dimension to 3
    elif not isinstance(dimension, int) or dimension <= 0:
        raise ValueError("dimension must be a positive integer.")

    # Validate and set rank
    rank = kwargs.get("rank")
    if not rank or not isinstance(rank, int) or not 0 < rank <= 4:
        raise ValueError("rank must be an integer between 1 and 4.")

    # Create the base tensor full of unique symbolic elements
    indexedexp = create_tensor_symbolic(rank * [dimension], idx_expr_basename)

    # Apply symmetry if specified
    symmetry = kwargs.get("symmetry")
    if symmetry:
        if not symmetry.startswith(("sym", "anti", "nosym")):
            raise ValueError(
                f"Unsupported symmetry '{symmetry}' for indexed expression. "
                "Valid symmetry options must start with 'sym', 'anti', or 'nosym'"
            )
        if rank > 1:
            indexedexp = _apply_symmetrization_globally(
                indexedexp, rank, symmetry, dimension
            )
        elif symmetry != "nosym":
            raise ValueError("Cannot symmetrize an indexed expression of rank 1.")

    # Zero out derivatives across symmetry axes, if any are defined
    return zero_out_derivatives_across_symmetry_axes(indexedexp)


def symmetrize_rank2(
    indexedexp: _rank2_type,
    symmetry: str,
    dimension: int,
) -> _rank2_type:
    """
    Symmetrize a 2D list (matrix) based on the given symmetry and dimension.

    :param indexedexp: The 2D list or matrix on which to perform the operation.
    :param symmetry: The symmetry operation to perform. Options include "sym01" or "nosym".
    :param dimension: The dimension of the matrix.
    :return: The matrix after applying the symmetry operation.

    Doctest:
    >>> a, b, c, d = sp.symbols('a b c d', real=True)
    >>> print(symmetrize_rank2([[a, b], [c, d]], "sym01", 2))
    [[a, b], [b, d]]
    """
    return cast(
        _rank2_type,
        _apply_symmetrization_globally(indexedexp, 2, symmetry, dimension),
    )


def symmetrize_rank3(
    indexedexp: _rank3_type,
    symmetry: str,
    dimension: int,
) -> _rank3_type:
    """
    Symmetrize a 3D list (tensor) based on the given symmetry and dimension.

    :param indexedexp: The 3D list or tensor on which to perform the operation.
    :param symmetry: The symmetry operation to perform. Options include "sym01", "sym02", "sym12", "nosym", and others.
    :param dimension: The dimension of the tensor.
    :return: The tensor after applying the symmetry operation.

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
    return cast(
        _rank3_type,
        _apply_symmetrization_globally(indexedexp, 3, symmetry, dimension),
    )


def symmetrize_rank4(
    indexedexp: _rank4_type,
    symmetry: str,
    dimension: int,
) -> _rank4_type:
    """
    Symmetrize a 4D list (4-tensor) based on the given symmetry and dimension.

    :param indexedexp: The 4D list or 4-tensor on which to perform the operation.
    :param symmetry: The symmetry operation to perform. Options include "sym01", "sym02", "sym12", "sym03", "sym13", "sym23", "nosym", and others.
    :param dimension: The dimension of the 4-tensor.
    :return: The 4-tensor after applying the symmetry operation.

    Doctest:
    >>> aDDDD = declarerank4("aDDDD", dimension=3)
    >>> asymDDDD = symmetrize_rank4(aDDDD, symmetry="sym12", dimension=3)
    >>> for i in range(3):
    ...     for j in range(3):
    ...         for k in range(3):
    ...             if asymDDDD[i][j][k][2] - asymDDDD[i][k][j][2]:
    ...                 raise ValueError("ERROR! symmetrize_rank4 didn't work.")

    """
    return cast(
        _rank4_type,
        _apply_symmetrization_globally(indexedexp, 4, symmetry, dimension),
    )


# typehinting: Must allow for sp.Expr in the zerorank*()s, so that it can be modified later
def zerorank1(dimension: int = -1) -> List[Union[sp.Expr, sp.Symbol]]:
    """
    Initialize rank-1 indexed expression to zero.

    :param dimension: The dimension of the rank-1 indexed expression. Default is -1.
    :return: A rank-1 indexed expression initialized to zero.
    """
    return declare_indexedexp(None, rank=1, dimension=dimension)  # type:ignore


def zerorank2(dimension: int = -1) -> List[List[Union[sp.Expr, sp.Symbol]]]:
    """
    Initialize rank-2 indexed expression to zero.

    :param dimension: The dimension of the rank-2 indexed expression. Default is -1.
    :return: A rank-2 indexed expression initialized to zero.
    """
    return declare_indexedexp(None, rank=2, dimension=dimension)  # type:ignore


def zerorank3(dimension: int = -1) -> List[List[List[Union[sp.Expr, sp.Symbol]]]]:
    """
    Initialize rank-3 indexed expression to zero.

    :param dimension: The dimension of the rank-3 indexed expression. Default is -1.
    :return: A rank-3 indexed expression initialized to zero.
    """
    return declare_indexedexp(None, rank=3, dimension=dimension)  # type:ignore


def zerorank4(dimension: int = -1) -> List[List[List[List[Union[sp.Expr, sp.Symbol]]]]]:
    """
    Initialize rank-4 indexed expression to zero.

    :param dimension: The dimension of the rank-4 indexed expression. Default is -1.
    :return: A rank-4 indexed expression initialized to zero.
    """
    return declare_indexedexp(None, rank=4, dimension=dimension)  # type:ignore


def get_rank(IDX_EXPR: Sequence[_recur_symbol_type]) -> int:
    """
    Get the rank of the expression based on the nested list structure.

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
    temp_expr = IDX_EXPR
    while isinstance(temp_expr, list):
        rank += 1
        if not temp_expr:
            break
        temp_expr = temp_expr[0]
    return rank


def _symbol_has_derivative_across_symmetry_axis(
    symbol_name: str, symmetry_axes: str
) -> bool:
    """
    Check if a symbol's name indicates it's a derivative taken across a symmetry axis.

    This function assumes derivative indices are the final characters of the symbol name.
    For example, if '2' is a symmetry axis, 'alpha_dD2' should be zeroed.

    :param symbol_name: The string name of the symbol to check (e.g., 'alpha_dD0').
    :param symmetry_axes: A string containing axes of symmetry (e.g., "01").
    :return: True if the symbol should be zeroed.
    :raises ValueError: For unsupported derivative orders (greater than 2).
    """
    if "_d" not in symbol_name:
        return False

    # Find the derivative order (number of 'D's after '_d').
    # This logic is preserved from the original code.
    deriv_order = 1
    if "DDD" in symbol_name:
        raise ValueError(
            f"Error: Derivative order > 2 not supported. Failed expression: {symbol_name}"
        )
    if "DD" in symbol_name:
        deriv_order = 2

    # The last `deriv_order` characters in the symbol name are the derivative indices.
    # Check if any of these indices are in the set of symmetry axes.
    derivative_indices = symbol_name[-deriv_order:]
    return any(idx in symmetry_axes for idx in derivative_indices)


def _recursive_zero_out(
    expression: _recur_symbol_type, symmetry_axes: str
) -> _recur_symbol_type:
    """
    Recursively traverse a tensor, zeroing elements that have derivs on symmetry axes.

    :param expression: The tensor or sub-tensor to process.
    :param symmetry_axes: A string containing axes of symmetry.
    :return: The processed tensor with derivatives on symmetry axes zeroed.
    """
    if not isinstance(expression, list):
        # Base case: we are at a SymPy expression (a scalar).
        if _symbol_has_derivative_across_symmetry_axis(str(expression), symmetry_axes):
            return sp.sympify(0)
        return expression

    # Recursive step: process each item in the list and reconstruct the list.
    return [_recursive_zero_out(item, symmetry_axes) for item in expression]


def zero_out_derivatives_across_symmetry_axes(
    IDX_EXPR: Sequence[_recur_symbol_type],
) -> _recur_symbol_type:
    """
    Zero derivatives across specified symmetry axes in an indexed expression.

    :param IDX_EXPR: Indexed expression to process.
    :return: The modified expression with derivatives across symmetry axes set to zero.
    :raises ValueError: If the expression indicates a derivative of order greater than 2,
                        as this is not supported.

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
    ValueError: Error: Derivative order > 2 not supported. Failed expression: trK_dDDD000
    >>> par.set_parval_from_str("symmetry_axes", "01")
    >>> zero_out_derivatives_across_symmetry_axes(declarerank2("trK_dDD"))
    [[0, 0, 0], [0, 0, 0], [0, 0, trK_dDD22]]
    >>> par.set_parval_from_str("symmetry_axes", "0")
    >>> zero_out_derivatives_across_symmetry_axes(declarerank4("aDD_dDD", symmetry="sym01_sym23"))
    [[[[0, 0, 0], [0, aDD_dDD0011, aDD_dDD0012], [0, aDD_dDD0012, aDD_dDD0022]], [[0, 0, 0], [0, aDD_dDD0111, aDD_dDD0112], [0, aDD_dDD0112, aDD_dDD0122]], [[0, 0, 0], [0, aDD_dDD0211, aDD_dDD0212], [0, aDD_dDD0212, aDD_dDD0222]]], [[[0, 0, 0], [0, aDD_dDD0111, aDD_dDD0112], [0, aDD_dDD0112, aDD_dDD0122]], [[0, 0, 0], [0, aDD_dDD1111, aDD_dDD1112], [0, aDD_dDD1112, aDD_dDD1122]], [[0, 0, 0], [0, aDD_dDD1211, aDD_dDD1212], [0, aDD_dDD1212, aDD_dDD1222]]], [[[0, 0, 0], [0, aDD_dDD0211, aDD_dDD0212], [0, aDD_dDD0212, aDD_dDD0222]], [[0, 0, 0], [0, aDD_dDD1211, aDD_dDD1212], [0, aDD_dDD1212, aDD_dDD1222]], [[0, 0, 0], [0, aDD_dDD2211, aDD_dDD2212], [0, aDD_dDD2212, aDD_dDD2222]]]]
    """
    symmetry_axes = par.parval_from_str("indexedexp::symmetry_axes")
    if not symmetry_axes:
        return IDX_EXPR

    # This rank check is preserved from the original implementation.
    rank = get_rank(IDX_EXPR)
    if rank > 4 or rank == 0:
        raise ValueError(
            f"Error: Only ranks from 1 to 4 are supported; {IDX_EXPR} does not seem to fit this pattern."
        )

    # Use a recursive helper to traverse the tensor structure of any rank.
    return _recursive_zero_out(IDX_EXPR, symmetry_axes)


def declarerank1(idx_expr_basename: str, **kwargs: Any) -> _rank1_type:
    """
    Declare a rank-1 indexed expression using the provided base name.

    :param idx_expr_basename: The base name for the indexed expression.
    :param kwargs: Additional arguments passed to declare_indexedexp.
    :return: The declared rank-1 indexed expression.
    """
    kwargs["rank"] = 1
    return cast(_rank1_type, declare_indexedexp(idx_expr_basename, **kwargs))


def declarerank2(idx_expr_basename: str, **kwargs: Any) -> _rank2_type:
    """
    Declare a rank-2 indexed expression using the provided base name.

    :param idx_expr_basename: The base name for the indexed expression.
    :param kwargs: Additional arguments passed to declare_indexedexp.
    :return: The declared rank-2 indexed expression.
    """
    kwargs["rank"] = 2
    return cast(_rank2_type, declare_indexedexp(idx_expr_basename, **kwargs))


def declarerank3(idx_expr_basename: str, **kwargs: Any) -> _rank3_type:
    """
    Declare a rank-3 indexed expression using the provided base name.

    :param idx_expr_basename: The base name for the indexed expression.
    :param kwargs: Additional arguments passed to declare_indexedexp.
    :return: The declared rank-3 indexed expression.
    """
    kwargs["rank"] = 3
    return cast(_rank3_type, declare_indexedexp(idx_expr_basename, **kwargs))


def declarerank4(idx_expr_basename: str, **kwargs: Any) -> _rank4_type:
    """
    Declare a rank-4 indexed expression using the provided base name.

    :param idx_expr_basename: The base name for the indexed expression.
    :param kwargs: Additional arguments passed to declare_indexedexp.
    :return: The declared rank-4 indexed expression.
    """
    kwargs["rank"] = 4
    return cast(_rank4_type, declare_indexedexp(idx_expr_basename, **kwargs))


class NonInvertibleMatrixError(ZeroDivisionError):
    """Matrix Not Invertible; Division By Zero."""


# The following matrix inverter functions are hard-coded for performance.
# They implement the analytic solution for the inverse of a matrix of a given size.


def symm_matrix_inverter2x2(
    a: List[List[Union[sp.Expr, sp.Symbol]]],
) -> Tuple[List[List[Union[sp.Expr, sp.Symbol]]], Union[sp.Expr, sp.Symbol]]:
    """
    Calculate the inverse and determinant of a symmetric 2x2 matrix.

    :param a: A symmetric 2x2 matrix of sympy Expressions or Symbols.
    :return: Tuple containing:
             - A symmetric 2x2 inverse matrix of sympy Expressions or Symbols.
             - A determinant of the input matrix as a sympy Expression or Symbol.
    :raises NonInvertibleMatrixError: If the input matrix has determinant zero.

    Doctest:
    >>> from sympy import symbols
    >>> b, c, d = symbols('b c d')
    >>> matrix = [[b, c], [c, d]]
    >>> inverse, determinant = symm_matrix_inverter2x2(matrix)
    """
    # Using the analytic solution for the determinant of a symmetric 2x2 matrix.
    outDET = a[0][0] * a[1][1] - a[0][1] ** 2
    if outDET == 0:
        raise NonInvertibleMatrixError("matrix has determinant zero")

    inv_det = 1 / outDET
    outINV = [[sp.sympify(0) for _ in range(2)] for __ in range(2)]

    outINV[0][0] = a[1][1] * inv_det
    outINV[0][1] = -a[0][1] * inv_det
    outINV[1][1] = a[0][0] * inv_det
    # Symmetric matrix property
    outINV[1][0] = outINV[0][1]
    return outINV, outDET


def symm_matrix_inverter3x3(
    a: Sequence[Sequence[Union[sp.Expr, sp.Symbol]]],
) -> Tuple[List[List[Union[sp.Expr, sp.Symbol]]], Union[sp.Expr, sp.Symbol]]:
    """
    Calculate the inverse and determinant of a symmetric 3x3 matrix.

    :param a: A symmetric 3x3 matrix of sympy Expressions or Symbols.
    :return: Tuple containing:
             - A symmetric 3x3 inverse matrix of sympy Expressions or Symbols.
             - A determinant of the input matrix as a sympy Expression or Symbol.
    :raises NonInvertibleMatrixError: If the input matrix has determinant zero.

    Doctest:
    >>> from sympy import symbols
    >>> a0, b, c, d, e, f = symbols('a0 b c d e f')
    >>> matrix = [[a0, b, c], [b, d, e], [c, e, f]]
    >>> inverse, determinant = symm_matrix_inverter3x3(matrix)
    """
    # Using the analytic solution for the determinant of a symmetric 3x3 matrix.
    outDET = (
        -a[0][2] ** 2 * a[1][1]
        + 2 * a[0][1] * a[0][2] * a[1][2]
        - a[0][0] * a[1][2] ** 2
        - a[0][1] ** 2 * a[2][2]
        + a[0][0] * a[1][1] * a[2][2]
    )
    if outDET == 0:
        raise NonInvertibleMatrixError("matrix has determinant zero")

    inv_det = 1 / outDET
    outINV = [[sp.sympify(0) for _ in range(3)] for __ in range(3)]

    # Using the analytic solution for the matrix inverse.
    outINV[0][0] = (-a[1][2] ** 2 + a[1][1] * a[2][2]) * inv_det
    outINV[0][1] = (a[0][2] * a[1][2] - a[0][1] * a[2][2]) * inv_det
    outINV[0][2] = (-a[0][2] * a[1][1] + a[0][1] * a[1][2]) * inv_det
    outINV[1][1] = (-a[0][2] ** 2 + a[0][0] * a[2][2]) * inv_det
    outINV[1][2] = (a[0][1] * a[0][2] - a[0][0] * a[1][2]) * inv_det
    outINV[2][2] = (-a[0][1] ** 2 + a[0][0] * a[1][1]) * inv_det
    # Symmetric matrix property
    outINV[1][0] = outINV[0][1]
    outINV[2][0] = outINV[0][2]
    outINV[2][1] = outINV[1][2]
    return outINV, outDET


# Validation test code for symm_matrix_inverter4x4() is available in the original file.
def symm_matrix_inverter4x4(
    a: List[List[Union[sp.Expr, sp.Symbol]]],
) -> Tuple[List[List[Union[sp.Expr, sp.Symbol]]], Union[sp.Expr, sp.Symbol]]:
    """
    Calculate the inverse and determinant of a 4x4 symmetric matrix.

    :param a: The 4x4 symmetric matrix to be inverted.
              The matrix elements can be either symbolic expressions or symbols.
    :return: Tuple containing:
             - A 4x4 inverse matrix of sympy Expressions or Symbols.
             - A determinant of the input matrix as a sympy Expression or Symbol.
    :raises NonInvertibleMatrixError: If the input matrix has determinant zero.

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
    # Hard-coded formula for the determinant of a 4x4 symmetric matrix for performance.
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

    inv_det = 1 / outDET
    outINV = [[sp.sympify(0) for _ in range(4)] for __ in range(4)]

    # Hard-coded formula for the inverse of a 4x4 symmetric matrix for performance.
    outINV[0][0] = (
        -a[1][3] * a[1][3] * a[2][2]
        + 2 * a[1][2] * a[1][3] * a[2][3]
        - a[1][1] * a[2][3] * a[2][3]
        - a[1][2] * a[1][2] * a[3][3]
        + a[1][1] * a[2][2] * a[3][3]
    ) * inv_det
    outINV[1][1] = (
        -a[0][3] * a[0][3] * a[2][2]
        + 2 * a[0][2] * a[0][3] * a[2][3]
        - a[0][0] * a[2][3] * a[2][3]
        - a[0][2] * a[0][2] * a[3][3]
        + a[0][0] * a[2][2] * a[3][3]
    ) * inv_det
    outINV[2][2] = (
        -a[0][3] * a[0][3] * a[1][1]
        + 2 * a[0][1] * a[0][3] * a[1][3]
        - a[0][0] * a[1][3] * a[1][3]
        - a[0][1] * a[0][1] * a[3][3]
        + a[0][0] * a[1][1] * a[3][3]
    ) * inv_det
    outINV[3][3] = (
        -a[0][2] * a[0][2] * a[1][1]
        + 2 * a[0][1] * a[0][2] * a[1][2]
        - a[0][0] * a[1][2] * a[1][2]
        - a[0][1] * a[0][1] * a[2][2]
        + a[0][0] * a[1][1] * a[2][2]
    ) * inv_det
    outINV[0][1] = (
        +a[0][3] * a[1][3] * a[2][2]
        - a[0][3] * a[1][2] * a[2][3]
        - a[0][2] * a[1][3] * a[2][3]
        + a[0][1] * a[2][3] * a[2][3]
        + a[0][2] * a[1][2] * a[3][3]
        - a[0][1] * a[2][2] * a[3][3]
    ) * inv_det
    outINV[0][2] = (
        -a[0][3] * a[1][2] * a[1][3]
        + a[0][2] * a[1][3] * a[1][3]
        + a[0][3] * a[1][1] * a[2][3]
        - a[0][1] * a[1][3] * a[2][3]
        - a[0][2] * a[1][1] * a[3][3]
        + a[0][1] * a[1][2] * a[3][3]
    ) * inv_det
    outINV[0][3] = (
        -a[0][2] * a[1][2] * a[1][3]
        + a[0][1] * a[1][3] * a[2][2]
        + a[0][3] * a[1][2] * a[1][2]
        - a[0][3] * a[1][1] * a[2][2]
        + a[0][2] * a[1][1] * a[2][3]
        - a[0][1] * a[1][2] * a[2][3]
    ) * inv_det
    outINV[1][2] = (
        +a[0][3] * a[0][3] * a[1][2]
        + a[0][0] * a[1][3] * a[2][3]
        - a[0][3] * a[0][2] * a[1][3]
        - a[0][3] * a[0][1] * a[2][3]
        + a[0][1] * a[0][2] * a[3][3]
        - a[0][0] * a[1][2] * a[3][3]
    ) * inv_det
    outINV[1][3] = (
        +a[0][2] * a[0][2] * a[1][3]
        + a[0][1] * a[0][3] * a[2][2]
        - a[0][0] * a[1][3] * a[2][2]
        + a[0][0] * a[1][2] * a[2][3]
        - a[0][2] * a[0][3] * a[1][2]
        - a[0][2] * a[0][1] * a[2][3]
    ) * inv_det
    outINV[2][3] = (
        +a[0][2] * a[0][3] * a[1][1]
        - a[0][1] * a[0][3] * a[1][2]
        - a[0][1] * a[0][2] * a[1][3]
        + a[0][0] * a[1][2] * a[1][3]
        + a[0][1] * a[0][1] * a[2][3]
        - a[0][0] * a[1][1] * a[2][3]
    ) * inv_det

    # Symmetric matrix property fills the lower triangle.
    outINV[1][0] = outINV[0][1]
    outINV[2][0] = outINV[0][2]
    outINV[2][1] = outINV[1][2]
    outINV[3][0] = outINV[0][3]
    outINV[3][1] = outINV[1][3]
    outINV[3][2] = outINV[2][3]

    return outINV, outDET


def generic_matrix_inverter2x2(
    a: List[List[Union[sp.Expr, sp.Symbol]]],
) -> Tuple[List[List[Union[sp.Expr, sp.Symbol]]], Union[sp.Expr, sp.Symbol]]:
    """
    Calculate the inverse and determinant of a general 2x2 matrix.

    This function first calculates the determinant of the input matrix. If the determinant is zero,
    it raises a NonInvertibleMatrixError. It then calculates the inverse matrix.

    :param a: A 2x2 matrix of sympy Expressions or Symbols.
    :return: Tuple containing the inverse matrix and the determinant.
    :raises NonInvertibleMatrixError: If the input matrix has determinant zero.

    Doctest:
    >>> matrix = declarerank2("gDD", dimension=2)
    >>> gUU, detg = generic_matrix_inverter2x2(matrix)
    """
    outDET = a[0][0] * a[1][1] - a[0][1] * a[1][0]
    if outDET == 0:
        raise NonInvertibleMatrixError("matrix has determinant zero")

    inv_det = 1 / outDET
    outINV = [[sp.sympify(0) for _ in range(2)] for __ in range(2)]

    outINV[0][0] = a[1][1] * inv_det
    outINV[0][1] = -a[0][1] * inv_det
    outINV[1][0] = -a[1][0] * inv_det
    outINV[1][1] = a[0][0] * inv_det
    return outINV, outDET


def generic_matrix_inverter3x3(
    a: List[List[Union[sp.Expr, sp.Symbol]]],
) -> Tuple[List[List[Union[sp.Expr, sp.Symbol]]], Union[sp.Expr, sp.Symbol]]:
    """
    Calculate the inverse and determinant of a general 3x3 matrix.

    This function first calculates the determinant of the input matrix. If the determinant is zero,
    it raises a NonInvertibleMatrixError. It then calculates the inverse matrix.

    :param a: A 3x3 matrix of sympy Expressions or Symbols.
    :return: Tuple containing the inverse matrix and the determinant.
    :raises NonInvertibleMatrixError: If the input matrix has determinant zero.

    Doctest:
    >>> gDD = declarerank2("gDD")
    >>> gUU, detg = generic_matrix_inverter3x3(gDD)
    """
    # Hard-coded determinant formula.
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

    inv_det = 1 / outDET
    outINV = [[sp.sympify(0) for _ in range(3)] for __ in range(3)]

    # Hard-coded matrix inverse formula (Cramer's rule).
    outINV[0][0] = (-a[1][2] * a[2][1] + a[1][1] * a[2][2]) * inv_det
    outINV[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * inv_det
    outINV[0][2] = (-a[0][2] * a[1][1] + a[0][1] * a[1][2]) * inv_det
    outINV[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * inv_det
    outINV[1][1] = (-a[0][2] * a[2][0] + a[0][0] * a[2][2]) * inv_det
    outINV[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * inv_det
    outINV[2][0] = (-a[1][1] * a[2][0] + a[1][0] * a[2][1]) * inv_det
    outINV[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * inv_det
    outINV[2][2] = (-a[0][1] * a[1][0] + a[0][0] * a[1][1]) * inv_det

    return outINV, outDET


def generic_matrix_inverter4x4(
    a: List[List[Union[sp.Expr, sp.Symbol]]],
) -> Tuple[List[List[Union[sp.Expr, sp.Symbol]]], Union[sp.Expr, sp.Symbol]]:
    """
    Calculate the inverse and determinant of a general 4x4 matrix.

    This function first calculates the determinant of the input matrix. If the determinant is zero,
    it raises a NonInvertibleMatrixError. It then calculates the inverse matrix.

    :param a: A 4x4 matrix of sympy Expressions or Symbols.
    :return: Tuple containing the inverse matrix and the determinant.
    :raises NonInvertibleMatrixError: If the input matrix has determinant zero.

    Doctest:
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

    inv_det = 1 / outDET
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
    ) * inv_det
    outINV[0][1] = (
        a[0][3] * a[2][2] * a[3][1]
        - a[0][2] * a[2][3] * a[3][1]
        - a[0][3] * a[2][1] * a[3][2]
        + a[0][1] * a[2][3] * a[3][2]
        + a[0][2] * a[2][1] * a[3][3]
        - a[0][1] * a[2][2] * a[3][3]
    ) * inv_det
    outINV[0][2] = (
        -a[0][3] * a[1][2] * a[3][1]
        + a[0][2] * a[1][3] * a[3][1]
        + a[0][3] * a[1][1] * a[3][2]
        - a[0][1] * a[1][3] * a[3][2]
        - a[0][2] * a[1][1] * a[3][3]
        + a[0][1] * a[1][2] * a[3][3]
    ) * inv_det
    outINV[0][3] = (
        a[0][3] * a[1][2] * a[2][1]
        - a[0][2] * a[1][3] * a[2][1]
        - a[0][3] * a[1][1] * a[2][2]
        + a[0][1] * a[1][3] * a[2][2]
        + a[0][2] * a[1][1] * a[2][3]
        - a[0][1] * a[1][2] * a[2][3]
    ) * inv_det
    outINV[1][0] = (
        a[1][3] * a[2][2] * a[3][0]
        - a[1][2] * a[2][3] * a[3][0]
        - a[1][3] * a[2][0] * a[3][2]
        + a[1][0] * a[2][3] * a[3][2]
        + a[1][2] * a[2][0] * a[3][3]
        - a[1][0] * a[2][2] * a[3][3]
    ) * inv_det
    outINV[1][1] = (
        -a[0][3] * a[2][2] * a[3][0]
        + a[0][2] * a[2][3] * a[3][0]
        + a[0][3] * a[2][0] * a[3][2]
        - a[0][0] * a[2][3] * a[3][2]
        - a[0][2] * a[2][0] * a[3][3]
        + a[0][0] * a[2][2] * a[3][3]
    ) * inv_det
    outINV[1][2] = (
        a[0][3] * a[1][2] * a[3][0]
        - a[0][2] * a[1][3] * a[3][0]
        - a[0][3] * a[1][0] * a[3][2]
        + a[0][0] * a[1][3] * a[3][2]
        + a[0][2] * a[1][0] * a[3][3]
        - a[0][0] * a[1][2] * a[3][3]
    ) * inv_det
    outINV[1][3] = (
        -a[0][3] * a[1][2] * a[2][0]
        + a[0][2] * a[1][3] * a[2][0]
        + a[0][3] * a[1][0] * a[2][2]
        - a[0][0] * a[1][3] * a[2][2]
        - a[0][2] * a[1][0] * a[2][3]
        + a[0][0] * a[1][2] * a[2][3]
    ) * inv_det
    outINV[2][0] = (
        -a[1][3] * a[2][1] * a[3][0]
        + a[1][1] * a[2][3] * a[3][0]
        + a[1][3] * a[2][0] * a[3][1]
        - a[1][0] * a[2][3] * a[3][1]
        - a[1][1] * a[2][0] * a[3][3]
        + a[1][0] * a[2][1] * a[3][3]
    ) * inv_det
    outINV[2][1] = (
        a[0][3] * a[2][1] * a[3][0]
        - a[0][1] * a[2][3] * a[3][0]
        - a[0][3] * a[2][0] * a[3][1]
        + a[0][0] * a[2][3] * a[3][1]
        + a[0][1] * a[2][0] * a[3][3]
        - a[0][0] * a[2][1] * a[3][3]
    ) * inv_det
    outINV[2][2] = (
        -a[0][3] * a[1][1] * a[3][0]
        + a[0][1] * a[1][3] * a[3][0]
        + a[0][3] * a[1][0] * a[3][1]
        - a[0][0] * a[1][3] * a[3][1]
        - a[0][1] * a[1][0] * a[3][3]
        + a[0][0] * a[1][1] * a[3][3]
    ) * inv_det
    outINV[2][3] = (
        a[0][3] * a[1][1] * a[2][0]
        - a[0][1] * a[1][3] * a[2][0]
        - a[0][3] * a[1][0] * a[2][1]
        + a[0][0] * a[1][3] * a[2][1]
        + a[0][1] * a[1][0] * a[2][3]
        - a[0][0] * a[1][1] * a[2][3]
    ) * inv_det
    outINV[3][0] = (
        a[1][2] * a[2][1] * a[3][0]
        - a[1][1] * a[2][2] * a[3][0]
        - a[1][2] * a[2][0] * a[3][1]
        + a[1][0] * a[2][2] * a[3][1]
        + a[1][1] * a[2][0] * a[3][2]
        - a[1][0] * a[2][1] * a[3][2]
    ) * inv_det
    outINV[3][1] = (
        -a[0][2] * a[2][1] * a[3][0]
        + a[0][1] * a[2][2] * a[3][0]
        + a[0][2] * a[2][0] * a[3][1]
        - a[0][0] * a[2][2] * a[3][1]
        - a[0][1] * a[2][0] * a[3][2]
        + a[0][0] * a[2][1] * a[3][2]
    ) * inv_det
    outINV[3][2] = (
        a[0][2] * a[1][1] * a[3][0]
        - a[0][1] * a[1][2] * a[3][0]
        - a[0][2] * a[1][0] * a[3][1]
        + a[0][0] * a[1][2] * a[3][1]
        + a[0][1] * a[1][0] * a[3][2]
        - a[0][0] * a[1][1] * a[3][2]
    ) * inv_det
    outINV[3][3] = (
        -a[0][2] * a[1][1] * a[2][0]
        + a[0][1] * a[1][2] * a[2][0]
        + a[0][2] * a[1][0] * a[2][1]
        - a[0][0] * a[1][2] * a[2][1]
        - a[0][1] * a[1][0] * a[2][2]
        + a[0][0] * a[1][1] * a[2][2]
    ) * inv_det

    return outINV, outDET


def LeviCivitaSymbol_dim3_rank3() -> List[List[List[int]]]:
    """
    Calculate the Levi-Civita symbol for 3 dimensions.

    This function creates a 3x3x3 rank-3 tensor where each element corresponds to the Levi-Civita symbol
    epsilon_ijk, calculated using the formula (i - j) * (j - k) * (k - i) / 2.

    :return: A 3x3x3 tensor with integer elements, representing the Levi-Civita symbol.

    Doctest:
    >>> LeviCivitaSymbol_dim3_rank3()
    [[[0, 0, 0], [0, 0, 1], [0, -1, 0]], [[0, 0, -1], [0, 0, 0], [1, 0, 0]], [[0, 1, 0], [-1, 0, 0], [0, 0, 0]]]
    """
    # Initialize a 3x3x3 tensor of zeros.
    levi_civita = zerorank3(dimension=3)

    for i, j, k in func.product(range(3), repeat=3):
        # This formula concisely generates the Levi-Civita symbol from indices.
        # See: https://codegolf.stackexchange.com/questions/160359/levi-civita-symbol
        levi_civita[i][j][k] = (i - j) * (j - k) * (k - i) // 2
    return cast(List[List[List[int]]], levi_civita)


def LeviCivitaTensorUUU_dim3_rank3(
    sqrtgammaDET: sp.Expr,
) -> Union[List[List[List[sp.Expr]]], List[List[List[int]]]]:
    """
    Calculate the Levi-Civita tensor with upper indices for 3 dimensions.

    This function creates a 3x3x3 rank-3 tensor where each element corresponds to the Levi-Civita symbol
    epsilon^ijk, calculated by dividing each Levi-Civita symbol by sqrtgammaDET.

    :param sqrtgammaDET: A sympy expression that is used to divide each element of the Levi-Civita symbol.
    :return: A 3x3x3 tensor with elements as sympy expressions, representing the Levi-Civita tensor.

    Doctest:
    >>> sqrtgammaDET = sp.Symbol("sqrtgammaDET")
    >>> LeviCivitaTensorUUU_dim3_rank3(sqrtgammaDET)
    [[[0, 0, 0], [0, 0, 1/sqrtgammaDET], [0, -1/sqrtgammaDET, 0]], [[0, 0, -1/sqrtgammaDET], [0, 0, 0], [1/sqrtgammaDET, 0, 0]], [[0, 1/sqrtgammaDET, 0], [-1/sqrtgammaDET, 0, 0], [0, 0, 0]]]
    """
    levi_civita_symbol = LeviCivitaSymbol_dim3_rank3()
    # To get the contravariant tensor density, divide by the sqrt of the metric determinant.
    return [
        [[elem / sqrtgammaDET for elem in row] for row in matrix]
        for matrix in levi_civita_symbol
    ]


def LeviCivitaTensorDDD_dim3_rank3(
    sqrtgammaDET: sp.Expr,
) -> Union[List[List[List[sp.Expr]]], List[List[List[int]]]]:
    """
    Calculate the Levi-Civita tensor with lower indices for 3 dimensions.

    This function creates a 3x3x3 rank-3 tensor where each element corresponds to the Levi-Civita symbol
    epsilon_ijk, calculated by multiplying each Levi-Civita symbol by sqrtgammaDET.

    :param sqrtgammaDET: A sympy expression that is used to multiply each element of the Levi-Civita symbol.
    :return: A 3x3x3 tensor with elements as sympy expressions, representing the Levi-Civita tensor.

    Doctest:
    >>> sqrtgammaDET = sp.Symbol("sqrtgammaDET")
    >>> LeviCivitaTensorDDD_dim3_rank3(sqrtgammaDET)
    [[[0, 0, 0], [0, 0, sqrtgammaDET], [0, -sqrtgammaDET, 0]], [[0, 0, -sqrtgammaDET], [0, 0, 0], [sqrtgammaDET, 0, 0]], [[0, sqrtgammaDET, 0], [-sqrtgammaDET, 0, 0], [0, 0, 0]]]
    """
    levi_civita_symbol = LeviCivitaSymbol_dim3_rank3()
    # To get the covariant tensor density, multiply by the sqrt of the metric determinant.
    return [
        [[elem * sqrtgammaDET for elem in row] for row in matrix]
        for matrix in levi_civita_symbol
    ]


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    # Validation test code for symm_matrix_inverter4x4():
    import nrpy.validate_expressions.validate_expressions as ve

    R4DD = cast(
        List[List[sp.Expr]], declarerank2("R4DD", symmetry="sym01", dimension=4)
    )

    # Compute R4DD's inverse:
    R4DDinv, det = symm_matrix_inverter4x4(R4DD)

    # Next matrix multiply: IsUnit = R^{-1} R
    IsUnit = zerorank2(dimension=4)
    for ii in range(4):
        for jj in range(4):
            for kk in range(4):
                IsUnit[ii][jj] += R4DDinv[ii][kk] * R4DD[kk][jj]
                # If you'd like to check R R^{-1} instead:
                # IsUnit[i][j] += R4DD[i][k] * R4DDinv[k][j]

    # Next check, is IsUnit == Unit matrix?!
    for diag in range(4):
        if not ve.check_zero(IsUnit[diag][diag] - 1):
            raise ValueError("4x4 matrix inversion failed, along diagonal!")
    for offdiag_i in range(4):
        for offdiag_j in range(4):
            if offdiag_i != offdiag_j:
                if not ve.check_zero(IsUnit[offdiag_i][offdiag_j]):
                    raise ValueError(
                        "4x4 matrix inversion failed, in off-diagonal element!"
                    )
    # ^^ all should output as True.

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
