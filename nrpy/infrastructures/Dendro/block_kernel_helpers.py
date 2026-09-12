# nrpy/infrastructures/Dendro/block_kernel_helpers.py
"""
Lower a registered NRPy expression set into Dendro kernel bodies.

Everything here is formulation-agnostic: it turns a mapping of right-hand-side
symbol to SymPy expression, plus an optional upwind control vector, into the
pointer bindings, point loop and parameter lists a Dendro CFunction needs.  The
physics that produces those expressions lives under ``general_relativity/``;
this module authors no field name, no physics default and no expression.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Callable, Iterable, List, Optional, Sequence, Set, Tuple, Union

import sympy as sp

import nrpy.grid as gri
import nrpy.params as par
from nrpy.finite_difference import (
    DERIVATIVE_FAMILIES,
    extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars,
    extract_list_of_deriv_var_strings_from_sympyexpr_list,
    stencil_reach_per_axis,
)
from nrpy.helpers.expression_utils import get_params_commondata_symbols_from_expr_list
from nrpy.infrastructures.Dendro import CodeParameters
from nrpy.infrastructures.Dendro import gridfunction_name_decorations as gf_names
from nrpy.infrastructures.Dendro.simple_loop import simple_loop


def cparam_declarations(names: Sequence[str]) -> str:
    """
    Build the trailing ``const <type> <cp>, ...`` parameter declarations.

    Each parameter is declared with its own registered ``cparam_type``, mapped
    by :func:`nrpy.infrastructures.Dendro.CodeParameters.c_type`, which is the
    same mapping the ``params_struct`` member uses.  Declaring everything as the
    scalar alias would coerce an ``int`` or ``bool`` parameter to floating point
    at the host call boundary.

    :param names: CodeParameter names to declare.
    :return: Comma-joined declarations, with no leading comma.

    Doctests:
    >>> import nrpy.params as par
    >>> _ = par.register_CodeParameter("REAL", __name__, "doctest_eta", 2.0)
    >>> _ = par.register_CodeParameter("int", __name__, "doctest_nsteps", 4)
    >>> cparam_declarations(("doctest_eta", "doctest_nsteps"))
    'const DendroScalar doctest_eta, const int doctest_nsteps'
    >>> cparam_declarations(())
    ''
    >>> del par.glb_code_params_dict["doctest_eta"]
    >>> del par.glb_code_params_dict["doctest_nsteps"]
    """
    return ", ".join(
        f"const {CodeParameters.c_type(par.glb_code_params_dict[name].cparam_type)}"
        f" {name}"
        for name in names
    )


def cparam_arguments(names: Sequence[str]) -> str:
    """
    Build the argument-forwarding list for the used CodeParameters.

    :param names: CodeParameter names to forward.
    :return: Comma-joined name list, with no trailing comma.

    Doctests:
    >>> cparam_arguments(("eta", "kappa1"))
    'eta, kappa1'
    """
    return ", ".join(names)


def used_codeparameters(expressions: Iterable[sp.Expr]) -> Tuple[str, ...]:
    """
    Return the registered CodeParameters the expressions actually reference.

    Reading the free symbols the equations contain is exact; scanning the
    emitted C text for names is not.  The result is sorted so the emitted
    CFunction parameter order is deterministic and the caller can forward the
    values in the same order.

    :param expressions: The expressions the kernel is generated from.
    :return: Registered CodeParameter names, sorted.

    Doctests:
    >>> import nrpy.params as par
    >>> _saved_params = dict(par.glb_code_params_dict)
    >>> try:
    ...     _z = par.register_CodeParameter("REAL", __name__, "zz_fixture_used", 1.0)
    ...     _a = par.register_CodeParameter(
    ...         "REAL", __name__, "aa_fixture_common", 2.0, commondata=True
    ...     )
    ...     used_codeparameters([_z + _a])
    ... finally:
    ...     par.glb_code_params_dict.clear()
    ...     par.glb_code_params_dict.update(_saved_params)
    ('aa_fixture_common', 'zz_fixture_used')
    >>> all(par.glb_code_params_dict[name] is value for name, value in _saved_params.items())
    True
    """
    param_symbols, commondata_symbols = get_params_commondata_symbols_from_expr_list(
        list(expressions)
    )
    return tuple(sorted(param_symbols + commondata_symbols))


def accessed_gridfunctions(expressions: Iterable[sp.Expr]) -> Set[str]:
    """
    Return the registered gridfunctions the expressions read.

    ``c_codegen`` names a derivative ``<field>_<op><component><direction>``, so
    a kernel that only differentiates a field never mentions the field's own
    symbol.  Intersecting raw free symbols with the registry therefore misses
    it, and the emitted kernel then reads a pointer nothing bound.  The
    canonical NRPy derivative extraction resolves each derivative symbol back
    to the field it differentiates, and the result is unioned with the fields
    the expressions read directly.

    The canonical extraction recognizes every family in
    :data:`nrpy.finite_difference.C_CODEGEN_DERIVATIVE_FAMILIES`, which is the
    set ``c_codegen`` can generate C for.

    :param expressions: The expressions the kernel is generated from.
    :return: Registered gridfunction names the expressions read.

    Doctests:
    >>> import nrpy.indexedexp as ixp
    >>> import nrpy.params as par
    >>> gri.glb_gridfcs_dict.clear()
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> _ = gri.register_gridfunctions_for_single_rank1("lambdaU", group="EVOL")
    >>> _ = gri.register_gridfunctions("cf", group="EVOL")
    >>> lambdaU_dD = ixp.declarerank2("lambdaU_dD")
    >>> cf_sym = sp.Symbol("cf")
    >>> sorted(accessed_gridfunctions([lambdaU_dD[0][1] + cf_sym]))
    ['cf', 'lambdaU0']
    >>> gri.glb_gridfcs_dict.clear()
    """
    free_symbols: List[sp.Basic] = []
    for expr in expressions:
        free_symbols.extend(expr.free_symbols)
    # The sentinel is the symbol form ``sp.Symbol("unset")`` that ``c_codegen``
    # itself passes.  With that form the canonical extraction also appends the
    # ``_ddnD`` twin of every ``_dupD`` derivative, which adds no base
    # gridfunction: a twin and its parent differentiate the same field.
    deriv_vars = extract_list_of_deriv_var_strings_from_sympyexpr_list(
        free_symbols, sp.Symbol("unset")
    )
    base_gridfunctions, _deriv_operators = (
        extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars(deriv_vars)
    )
    directly_read = {str(symbol) for symbol in free_symbols}
    return (set(base_gridfunctions) | directly_read) & set(gri.glb_gridfcs_dict)


def by_position(_name: str, position: int) -> str:
    """
    Return the registered registry position as the component index expression.

    The standalone-host translation units compile the kernels without the
    generated state header, so the bindings use the integer registry position
    rather than ``to_index(EvolVar::...)``; both orderings come from the same
    NRPy list.

    :param _name: Exact gridfunction name (unused; the index is positional).
    :param position: Position in the registered EVOL order.
    :return: The component index expression.
    """
    return str(position)


def output_component_bindings(
    names: Sequence[str],
    scalar_type: str,
    *,
    array: str,
    role: Callable[[str], str],
    const_pointee: bool,
    index_expression: Callable[[str, int], str],
    base_offset: Optional[str] = "geom.component_offset",
    flat_stride: Optional[str] = None,
) -> str:
    """
    Render one role-prefixed pointer binding per field, in the given order.

    :param names: Exact registered gridfunction names, in registry order.
    :param scalar_type: Generated scalar alias (e.g. ``DendroScalar``).
    :param array: Pointer-array name or flat base-pointer name.
    :param role: Function applying the pointer-role prefix.
    :param const_pointee: Whether the pointed-to scalar is const.
    :param index_expression: Maps a name and sequence position to an index.
    :param base_offset: Per-component base offset, or ``None`` for rebased data.
    :param flat_stride: Optional stride between fields in a flat layout.
    :return: C++ binding statements, one per line.

    Doctests:
    >>> output_component_bindings(
    ...     ("vU0", "vU1"), "DendroScalar", array="aux",
    ...     role=gf_names.input_pointer, const_pointee=True,
    ...     index_expression=lambda _name, position: str(position),
    ...     base_offset="offset",
    ... ).splitlines()
    ['const DendroScalar* const in_vU0 = aux[0] + offset;', 'const DendroScalar* const in_vU1 = aux[1] + offset;']
    """
    qualifier = (
        f"const {scalar_type}* const" if const_pointee else f"{scalar_type}* const"
    )
    lines: List[str] = []
    for position, name in enumerate(names):
        index = index_expression(name, position)
        offset_term = f" + {base_offset}" if base_offset is not None else ""
        if flat_stride is None:
            source = f"{array}[{index}]{offset_term}"
        else:
            source = (
                f"{array}{offset_term}"
                f" + static_cast<std::ptrdiff_t>({index}) * {flat_stride}"
            )
        lines.append(f"{qualifier} {role(name)} = {source};")
    return "\n".join(lines)


def block_pointer_bindings(evol_order: Sequence[str], scalar_type: str) -> str:
    """
    Emit the per-field input and RHS pointer bindings for the block layout.

    The bindings are rendered by the single shared emitter
    (:func:`nrpy.infrastructures.Dendro.block_kernel_helpers.output_component_bindings`)
    from the registry order, so no field name is hardcoded and the roles and
    per-component base offset stay aligned across the block-layout adapters.
    Every binding adds ``geom.component_offset``: the pointer arrays are
    allocation-relative, so a nonzero per-component base must be applied or
    multi-block layouts read the wrong component.

    :param evol_order: The EVOL names, in registry order.
    :param scalar_type: The registered Dendro scalar alias.
    :return: The binding statements.
    """
    return (
        output_component_bindings(
            evol_order,
            scalar_type,
            array="in_gfs",
            role=gf_names.input_pointer,
            const_pointee=True,
            index_expression=by_position,
        )
        + "\n"
        + output_component_bindings(
            evol_order,
            scalar_type,
            array="rhs_gfs",
            role=gf_names.rhs_pointer,
            const_pointee=False,
            index_expression=by_position,
        )
    )


def flat_block_pointer_bindings(evol_order: Sequence[str], scalar_type: str) -> str:
    """
    Emit the per-field bindings for the local-time-stepping flat-block layout.

    In this layout field ``f`` occupies ``in_gfs_flat + f * (nx * ny * nz)``.
    Extents are hoisted into a ``ptrdiff_t`` local and the same shared emitter
    renders both layouts.

    ``geom.component_offset`` is deliberately *not* applied here.  The adapter
    forwards these pointers together with the unchanged ``geom`` to the
    per-block kernel, whose own bindings apply the offset; applying it in both
    places addressed ``base + f * vol + 2 * component_offset``, so a nonzero
    offset read and wrote the wrong cells and could run past the allocation.
    The single per-block numerical body owns the offset.

    :param evol_order: The EVOL names, in registry order.
    :param scalar_type: The registered Dendro scalar alias.
    :return: The binding statements.

    Doctests:
    >>> bindings = flat_block_pointer_bindings(("aa", "bb"), "DendroScalar")
    >>> print(bindings.splitlines()[1])
    const DendroScalar* const in_aa = in_gfs_flat + static_cast<std::ptrdiff_t>(0) * vol;
    >>> "geom.component_offset" in bindings
    False
    """
    return "\n".join(
        [
            "const std::ptrdiff_t vol = static_cast<std::ptrdiff_t>(geom.nx)"
            " * geom.ny * geom.nz;",
            output_component_bindings(
                evol_order,
                scalar_type,
                array="in_gfs_flat",
                role=gf_names.input_pointer,
                const_pointee=True,
                index_expression=by_position,
                base_offset=None,
                flat_stride="vol",
            ),
            output_component_bindings(
                evol_order,
                scalar_type,
                array="rhs_gfs_flat",
                role=gf_names.rhs_pointer,
                const_pointee=False,
                index_expression=by_position,
                base_offset=None,
                flat_stride="vol",
            ),
        ]
    )


def point_loop(kernel: str, padding: str = "geom.padding") -> str:
    """
    Wrap a point kernel in the NRPy Dendro interior point loop.

    The block-geometry field names are the host contract, so they are spelled
    here once rather than at every builder.

    :param kernel: The ``c_codegen`` point kernel.
    :param padding: Ghost points to skip on each side; ``"0"`` for a pass that
        writes the padded halo as well, such as initial data.
    :return: The interior-loop-wrapped body.
    """
    return simple_loop(
        kernel,
        nx="geom.nx",
        ny="geom.ny",
        nz="geom.nz",
        padding=padding,
        pmin_padded="geom.pmin_padded",
        dx="geom.dx",
    )


def emitted_derivative_operators(
    expressions: Iterable[sp.Expr],
    upwind_control_vec: Union[List[sp.Basic], sp.Basic, str],
) -> Tuple[str, ...]:
    """
    Return the distinct derivative operators the expressions require.

    Taken from the expressions through the canonical NRPy extraction, not from
    the emitted C text: the operator names in the emitted kernel are a rendering
    of these, and reading the rendering back would be a second parser of the
    same naming scheme.

    :param expressions: The expressions the kernel is generated from.
    :param upwind_control_vec: The upwind control vector, or the string sentinel
        when upwinding is not enabled.
    :return: The operators, sorted, e.g. ``('dD0', 'dKOD1', 'dupD2')``.

    Doctests:
    >>> import nrpy.indexedexp as ixp
    >>> import nrpy.params as par
    >>> gri.glb_gridfcs_dict.clear()
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> _ = gri.register_gridfunctions("cf", group="EVOL")
    >>> cf_dD = ixp.declarerank1("cf_dD")
    >>> cf_dKOD = ixp.declarerank1("cf_dKOD")
    >>> emitted_derivative_operators([cf_dD[0] + cf_dKOD[1]], "unset")
    ('dD0', 'dKOD1')
    >>> gri.glb_gridfcs_dict.clear()
    """
    free_symbols: List[sp.Basic] = []
    for expr in expressions:
        free_symbols.extend(expr.free_symbols)
    deriv_vars = extract_list_of_deriv_var_strings_from_sympyexpr_list(
        free_symbols, upwind_control_vec, families=DERIVATIVE_FAMILIES
    )
    _base_gridfunctions, deriv_operators = (
        extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars(deriv_vars)
    )
    return tuple(sorted(set(deriv_operators)))


def padding_from_derivative_operators(
    expressions: Sequence[sp.Expr],
    upwind_control_vec: Union[List[sp.Basic], sp.Basic, str],
    fd_order: int,
) -> int:
    """
    Return the ghost points the expressions' derivatives reach.

    One number, the widest axis: the emitted point loop takes a single
    ``geom.padding`` and Dendro sizes a block's padding from its element order,
    so per-axis padding is unrepresentable in the host contract.

    This is not ``fd_order // 2``: the upwinded and Kreiss-Oliger families reach
    one point further than the centered ones (at fd_order 4, ``dupD`` reaches 3
    while ``dD`` reaches 2), so a radius-derived padding would read past the end
    of a Dendro block.  The reach itself comes from
    :func:`nrpy.finite_difference.stencil_reach_per_axis`, which reads the same
    coefficient source the kernel's C code is generated from.

    :param expressions: The expressions the kernel is generated from.
    :param upwind_control_vec: The upwind control vector, or the string sentinel
        when upwinding is not enabled.
    :param fd_order: The finite-difference order.
    :return: The widest numerical stencil reach over all axes.

    Doctests:
    >>> import nrpy.indexedexp as ixp
    >>> import nrpy.params as par
    >>> gri.glb_gridfcs_dict.clear()
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> _ = gri.register_gridfunctions("cf", group="EVOL")
    >>> cf_dD = ixp.declarerank1("cf_dD")
    >>> padding_from_derivative_operators([cf_dD[0] + cf_dD[1] + cf_dD[2]], "unset", 4)
    2

    The per-order reach is pinned in
    :func:`nrpy.finite_difference.stencil_reach_per_axis`.  This wrapper takes
    the maximum across axes, so the case below mixes a centered axis with an
    upwinded one: a uniform expression would pass just as well if this returned
    the minimum.

    >>> cf_dupD = ixp.declarerank1("cf_dupD")
    >>> mixed_axes = cf_dD[0] + cf_dupD[1] + cf_dD[2]
    >>> padding_from_derivative_operators([mixed_axes], "unset", 4)
    3

    Algebraic expressions require no numerical neighbours.  A derivative may
    use only one axis; Dendro still represents its reach with one uniform
    number.

    >>> padding_from_derivative_operators([sp.Symbol("cf")], "unset", 4)
    0
    >>> padding_from_derivative_operators([cf_dD[0]], "unset", 4)
    2
    >>> cf_dKOD = ixp.declarerank1("cf_dKOD")
    >>> padding_from_derivative_operators([cf_dKOD[2]], "unset", 4)
    3
    >>> gri.glb_gridfcs_dict.clear()
    """
    padding = stencil_reach_per_axis(expressions, upwind_control_vec, fd_order)
    return max(padding)


def upwind_control_fields_from_control_vec(
    upwind_control_vec: Iterable[sp.Expr], evol_order: Sequence[str]
) -> Tuple[str, ...]:
    """
    Return the EVOL fields that appear in the upwind control vector.

    Derived from the expressions rather than hardcoded, so the generated
    upwind-selection self-test drives exactly the fields the kernel switches
    on.

    :param upwind_control_vec: The control vector components.
    :param evol_order: The EVOL names, in registry order.
    :return: The control field names, in registry order.
    """
    # ``str`` round-trips the component so every control symbol is rebuilt
    # under SymPy's default assumptions, matching the plain ``sp.Symbol(name)``
    # built below: two symbols of the same name compare equal only when their
    # assumptions agree.  SymPy types ``free_symbols`` as a set of ``Basic``,
    # so the accumulator is annotated ``Set[sp.Basic]``.
    control_symbols: Set[sp.Basic] = set()
    for component in upwind_control_vec:
        control_symbols |= sp.sympify(str(component)).free_symbols
    return tuple(name for name in evol_order if sp.Symbol(name) in control_symbols)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
