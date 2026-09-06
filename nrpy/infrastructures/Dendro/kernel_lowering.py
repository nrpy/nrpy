# nrpy/infrastructures/Dendro/kernel_lowering.py
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

import re
from typing import Any, Dict, Iterable, List, Sequence, Set, Tuple

import sympy as sp

import nrpy.grid as gri
from nrpy.finite_difference import (
    compute_fdcoeffs_fdstencl,
    extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars,
    extract_list_of_deriv_var_strings_from_sympyexpr_list,
)
from nrpy.helpers.expression_utils import get_params_commondata_symbols_from_expr_list
from nrpy.infrastructures.Dendro import Dendro_state_h, naming
from nrpy.infrastructures.Dendro.simple_loop import simple_loop

# Finite-difference operator families emitted by c_codegen.  Field names use
# uppercase ``DD`` / ``U`` components, so these lowercase tokens identify
# derivative operators unambiguously.  The digit run that follows is
# ``<tensor component indices><derivative direction indices>`` -- e.g.
# ``hDD_dDD0112`` is the mixed second derivative in directions (1, 2) of the
# component (0, 1) -- so the run must be matched in full (``\d+``) and the
# *trailing* index/indices taken as the direction.  Matching only one or two
# digits made every derivative of a rank-1 or rank-2 field invisible and
# recorded a rank-1 first derivative under its component index instead of its
# direction.
_OPERATOR_RE = re.compile(r"_(dD|dDD|dupD|ddnD|dKOD|dfullupD|dfulldnD)(\d+)\b")

# Families carrying a single direction index; the rest carry two.  The
# classification matches nrpy/finite_difference.py, where `dfullupD` and
# `dfulldnD` are first-derivative (single-direction) operators.
_SINGLE_DIRECTION_FAMILIES = frozenset(
    {"dD", "dupD", "ddnD", "dKOD", "dfullupD", "dfulldnD"}
)


def cparam_declarations(names: Sequence[str], scalar_type: str) -> str:
    """
    Build the trailing ``const <scalar> <cp>, ...`` parameter declarations.

    :param names: CodeParameter names to declare.
    :param scalar_type: The registered Dendro scalar alias.
    :return: Comma-joined declarations, with no leading comma.

    Doctests:
    >>> cparam_declarations(("eta", "kappa1"), "DendroScalar")
    'const DendroScalar eta, const DendroScalar kappa1'
    >>> cparam_declarations((), "DendroScalar")
    ''
    """
    return ", ".join(f"const {scalar_type} {name}" for name in names)


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

    :param expressions: The lowered expressions.
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

    The canonical extraction recognizes the ``_dD``, ``_dDD``, ``_dKOD``,
    ``_dupD`` and ``_ddnD`` families.  It does not recognize ``_dfullupD`` or
    ``_dfulldnD``; no Dendro formulation emits those, and ``emitted_operators``
    below covers them when scanning emitted C text for the padding derivation.

    :param expressions: The lowered expressions.
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
    deriv_vars = extract_list_of_deriv_var_strings_from_sympyexpr_list(
        free_symbols, "unset"
    )
    base_gridfunctions, _deriv_operators = (
        extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars(deriv_vars)
    )
    directly_read = {str(symbol) for symbol in free_symbols}
    return (set(base_gridfunctions) | directly_read) & set(gri.glb_gridfcs_dict)


def _by_position(_name: str, position: int) -> str:
    """
    Return the registered registry position as the component index expression.

    The mock-vehicle translation units compile the kernels without the
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
    (:func:`nrpy.infrastructures.Dendro.Dendro_state_h.output_component_bindings`)
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
        Dendro_state_h.output_component_bindings(
            evol_order,
            scalar_type,
            array="in_gfs",
            role=naming.input_pointer,
            const_pointee=True,
            index_expression=_by_position,
        )
        + "\n"
        + Dendro_state_h.output_component_bindings(
            evol_order,
            scalar_type,
            array="rhs_gfs",
            role=naming.rhs_pointer,
            const_pointee=False,
            index_expression=_by_position,
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
            Dendro_state_h.output_component_bindings(
                evol_order,
                scalar_type,
                array="in_gfs_flat",
                role=naming.input_pointer,
                const_pointee=True,
                index_expression=_by_position,
                base_offset=None,
                flat_stride="vol",
            ),
            Dendro_state_h.output_component_bindings(
                evol_order,
                scalar_type,
                array="rhs_gfs_flat",
                role=naming.rhs_pointer,
                const_pointee=False,
                index_expression=_by_position,
                base_offset=None,
                flat_stride="vol",
            ),
        ]
    )


def point_loop(kernel: str) -> str:
    """
    Wrap a point kernel in the NRPy Dendro interior point loop.

    :param kernel: The ``c_codegen`` point kernel.
    :return: The interior-loop-wrapped body.
    """
    return simple_loop(
        kernel,
        nx="geom.nx",
        ny="geom.ny",
        nz="geom.nz",
        padding="geom.padding",
        pmin_padded="geom.pmin_padded",
        dx="geom.dx",
    )


def emitted_operators(kernel: str, fd_order: int) -> List[Dict[str, Any]]:
    """
    Derive the derivative operators the emitted kernel actually contains.

    The families are scanned from the generated code itself, and each
    operator's exact rational stencil comes from the same coefficient source
    the kernel was lowered with, so this is not a second stencil model.

    :param kernel: The emitted point kernel.
    :param fd_order: The finite-difference order.
    :return: One record per distinct operator: name, order, exact rational
        coefficient strings, signed offsets, and per-axis and total reach.
    :raises ValueError: If a derivative token in the kernel is malformed.
    """
    operators: List[Tuple[str, str, str]] = []
    for match in _OPERATOR_RE.finditer(kernel):
        base = match.group(1)
        idx = match.group(2)
        directions = idx[-1] if base in _SINGLE_DIRECTION_FAMILIES else idx[-2:]
        if len(directions) != (1 if base in _SINGLE_DIRECTION_FAMILIES else 2):
            raise ValueError(
                f"Malformed derivative token {base}{idx} in the emitted kernel."
            )
        derivstring = f"{base}{directions}"
        if derivstring not in {op[0] for op in operators}:
            operators.append((derivstring, base, directions))
    operator_records = []
    for derivstring, _base, _idx in operators:
        # Pass the raw fd_order: compute_fdcoeffs_fdstencl applies the +2 for
        # dKOD internally (exactly as c_codegen lowers the kernel), so this
        # reproduces the operator offsets from one coefficient source rather
        # than from a second stencil generator.
        coeffs, stencils = compute_fdcoeffs_fdstencl(derivstring, fd_order)
        max_offset = 0
        per_axis = [0, 0, 0]
        for stencil in stencils:
            for axis, step in enumerate(stencil):
                per_axis[axis] = max(per_axis[axis], abs(step))
                max_offset = max(max_offset, abs(step))
        operator_records.append(
            {
                "operator": derivstring,
                "fd_order": fd_order,
                "coefficients": [str(coeff) for coeff in coeffs],
                "offsets": [list(stencil) for stencil in stencils],
                "max_offset_per_axis": per_axis,
                "max_offset": max_offset,
            }
        )
    return operator_records


def padding_from_operators(
    operator_records: Sequence[Dict[str, Any]],
) -> Tuple[int, int, int]:
    """
    Return the ghost points the emitted operators reach, per axis.

    This is not ``fd_order // 2``: the upwinded and Kreiss-Oliger families
    reach one point further than the centred ones (at fd_order 4, ``dupD``
    reaches 3 while ``dD`` reaches 2), so a radius-derived padding would read
    past the end of a Dendro block.

    :param operator_records: Records from :func:`emitted_operators`.
    :return: The (px, py, pz) ghost points required.
    :raises ValueError: If the kernel reaches no neighbour on some axis.
    """
    per_axis = [0, 0, 0]
    for record in operator_records:
        for axis, reach in enumerate(record["max_offset_per_axis"]):
            per_axis[axis] = max(per_axis[axis], int(reach))
    padding = (per_axis[0], per_axis[1], per_axis[2])
    if min(padding) < 1:
        raise ValueError(
            f"The emitted kernel needs no ghost points ({padding}); a "
            "direct-FD right-hand side must read neighbours on every axis."
        )
    return padding


def upwind_control_fields(
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
    control_symbols: Set[sp.Symbol] = set()
    for component in upwind_control_vec:
        control_symbols |= set(sp.sympify(str(component)).free_symbols)
    return tuple(name for name in evol_order if sp.Symbol(name) in control_symbols)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
