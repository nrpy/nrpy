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

from typing import Iterable, List, Sequence, Set, Tuple, Union

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
from nrpy.infrastructures.Dendro import state_h
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


def block_pointer_bindings(evol_order: Sequence[str], scalar_type: str) -> str:
    """
    Emit the per-field input and RHS pointer bindings for the block layout.

    The bindings are rendered by the single shared emitter
    (:func:`nrpy.infrastructures.Dendro.state_h.output_component_bindings`)
    from the registry order, so no field name is hardcoded and the roles and
    per-component base offset cannot drift from the state-header renderer.
    Every binding adds ``geom.component_offset``: the pointer arrays are
    allocation-relative, so a nonzero per-component base must be applied or
    multi-block layouts read the wrong component.

    :param evol_order: The EVOL names, in registry order.
    :param scalar_type: The registered Dendro scalar alias.
    :return: The binding statements.
    """
    return (
        state_h.output_component_bindings(
            evol_order,
            scalar_type,
            array="in_gfs",
            role=gf_names.input_pointer,
            const_pointee=True,
            index_expression=by_position,
        )
        + "\n"
        + state_h.output_component_bindings(
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
            state_h.output_component_bindings(
                evol_order,
                scalar_type,
                array="in_gfs_flat",
                role=gf_names.input_pointer,
                const_pointee=True,
                index_expression=by_position,
                base_offset=None,
                flat_stride="vol",
            ),
            state_h.output_component_bindings(
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
    one point further than the centred ones (at fd_order 4, ``dupD`` reaches 3
    while ``dD`` reaches 2), so a radius-derived padding would read past the end
    of a Dendro block.  The reach itself comes from
    :func:`nrpy.finite_difference.stencil_reach_per_axis`, which reads the same
    coefficient source the kernel's C code is generated from.

    :param expressions: The expressions the kernel is generated from.
    :param upwind_control_vec: The upwind control vector, or the string sentinel
        when upwinding is not enabled.
    :param fd_order: The finite-difference order.
    :return: The ghost points required on every axis.
    :raises ValueError: If the expressions reach no neighbour on some axis.

    Doctests:
    >>> import nrpy.indexedexp as ixp
    >>> import nrpy.params as par
    >>> gri.glb_gridfcs_dict.clear()
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> _ = gri.register_gridfunctions("cf", group="EVOL")
    >>> cf_dD = ixp.declarerank1("cf_dD")
    >>> padding_from_derivative_operators([cf_dD[0] + cf_dD[1] + cf_dD[2]], "unset", 4)
    2
    >>> try:
    ...     padding_from_derivative_operators([sp.Symbol("cf")], "unset", 4)
    ... except ValueError as error:
    ...     print(str(error).split(";")[0])
    The expressions reach no ghost points ((0, 0, 0))
    >>> gri.glb_gridfcs_dict.clear()
    """
    padding = stencil_reach_per_axis(expressions, upwind_control_vec, fd_order)
    if min(padding) < 1:
        raise ValueError(
            f"The expressions reach no ghost points ({padding}); a direct-FD "
            "right-hand side must read neighbours on every axis."
        )
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
    print(f"Doctest passed: All {results.attempted} test(s) passed")
