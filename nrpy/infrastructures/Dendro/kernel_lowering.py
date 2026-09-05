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
import nrpy.params as par
from nrpy.finite_difference import compute_fdcoeffs_fdstencl
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

# A derivative symbol as c_codegen names it: <family>_<op><digits>, where the
# digit run carries the tensor component indices followed by the derivative
# direction indices.
_DERIVATIVE_SYMBOL_RE = re.compile(
    r"^(?P<family>.+)_(?P<op>dD|dDD|dupD|ddnD|dKOD|dfullupD|dfulldnD)"
    r"(?P<digits>[0-9]+)$"
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
    free_symbol_names = {
        str(symbol) for expr in expressions for symbol in expr.free_symbols
    }
    return tuple(sorted(free_symbol_names & set(par.glb_code_params_dict)))


def base_gridfunction_of(symbol_name: str) -> str:
    """
    Return the gridfunction a symbol reads, resolving derivative symbols.

    ``c_codegen`` names a derivative ``<family>_<op><component><direction>``,
    so a kernel that only differentiates a field never mentions the field's own
    symbol.  Intersecting raw free symbols with the registry therefore misses
    it, and the emitted kernel then reads a pointer nothing bound.

    :param symbol_name: A free-symbol name from a lowered expression.
    :return: The registered gridfunction name, or the input unchanged when the
        symbol is not a derivative.

    Doctests:
    >>> base_gridfunction_of("lambdaU_dD00")
    'lambdaU0'
    >>> base_gridfunction_of("hDD_dDD0112")
    'hDD01'
    >>> base_gridfunction_of("alpha_dupD2")
    'alpha'
    >>> base_gridfunction_of("cf")
    'cf'
    """
    match = _DERIVATIVE_SYMBOL_RE.match(symbol_name)
    if match is None:
        return symbol_name
    family, operator, digits = match.group("family", "op", "digits")
    directions = 1 if operator in _SINGLE_DIRECTION_FAMILIES else 2
    if len(digits) < directions:
        return symbol_name
    return family + digits[: len(digits) - directions]


def accessed_gridfunctions(expressions: Iterable[sp.Expr]) -> Set[str]:
    """
    Return the registered gridfunctions the expressions read.

    Derivative symbols are resolved back to the field they differentiate, so a
    kernel that only differentiates a field still reports it.

    :param expressions: The lowered expressions.
    :return: Registered gridfunction names the expressions read.
    """
    names = {
        base_gridfunction_of(str(symbol))
        for expr in expressions
        for symbol in expr.free_symbols
    }
    return names & set(gri.glb_gridfcs_dict)


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
    Extents are hoisted into a ``ptrdiff_t`` local and the per-component base
    ``geom.component_offset`` is applied exactly as in the block layout; the
    same shared emitter renders both layouts.

    :param evol_order: The EVOL names, in registry order.
    :param scalar_type: The registered Dendro scalar alias.
    :return: The binding statements.
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
                flat_stride="vol",
            ),
            Dendro_state_h.output_component_bindings(
                evol_order,
                scalar_type,
                array="rhs_gfs_flat",
                role=naming.rhs_pointer,
                const_pointee=False,
                index_expression=_by_position,
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
