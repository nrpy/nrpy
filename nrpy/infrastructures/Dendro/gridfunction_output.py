# nrpy/infrastructures/Dendro/gridfunction_output.py
"""
Pure state translator: render Dendro state artifacts from the frozen snapshot only.

This module contains no list of fCCZ4 field names and imports no equation
modules: every name, order, count, and metadata value is rendered from the
frozen NRPy gridfunction records.
"""

import json
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from nrpy.infrastructures.Dendro import naming
from nrpy.infrastructures.Dendro.freeze import FrozenNRPyDendroSnapshot

_BANNER_BODY = (
    "Generator: nrpy.infrastructures.Dendro\n"
    "Regenerate with: python -m nrpy.examples.dendro_fccz4\n"
)


def _banner(
    hash_line: str, source_object: str = "frozen NRPy gridfunction registry"
) -> str:
    """
    Build the generated-file banner.

    :param hash_line: Module ABI hash line.
    :param source_object: Frozen source the file was derived from.
    :return: Banner text with `//` comment prefixes (section 12.7).
    """
    lines = [
        "// GENERATED FILE - DO NOT EDIT",
        f"// {_BANNER_BODY.splitlines()[0]}",
        f"// {_BANNER_BODY.splitlines()[1]}",
        f"// Source object: {source_object}",
        f"// Module ABI: {hash_line}",
    ]
    return "\n".join(lines)


def state_records_by_group(snapshot: FrozenNRPyDendroSnapshot) -> Dict[str, List[str]]:
    """
    Group the frozen gridfunction names by group, in NRPy list order.

    :param snapshot: The frozen snapshot.
    :return: Mapping of group name to ordered list of exact names.
    """
    by_group: Dict[str, List[str]] = {"EVOL": [], "AUXEVOL": [], "DIAG": [], "AUX": []}
    for fg in snapshot.gridfunctions:
        by_group.setdefault(fg.group, []).append(fg.name)
    return {group: names for group, names in by_group.items() if names}


def render_fccz4_state_hpp(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the generated state header.

    It carries the EVOL enum, name array and section 5.5 metadata, and the
    section 14.7 exact-name lookup over every registered group.  Every name,
    count and metadata value is rendered from the frozen records in NRPy
    order; no field name is hardcoded here.

    :param snapshot: The frozen snapshot.
    :return: The complete C++ header text.
    """
    evol = [fg for fg in snapshot.gridfunctions if fg.group == "EVOL"]
    lines: List[str] = []
    lines.append(_banner(snapshot.hashes.module_abi_hash))
    lines.append("#pragma once")
    lines.append("")
    lines.append("#include <array>")
    lines.append("#include <cstddef>")
    lines.append("#include <optional>")
    lines.append("#include <string_view>")
    lines.append("")
    lines.append('#include "fccz4_types.hpp"')
    lines.append("")
    lines.append("namespace fccz4::generated {")
    lines.append("")
    lines.append("enum class EvolVar : unsigned {")
    for fg in evol:
        lines.append(f"    {naming.enum_member(fg.name)} = {fg.index},")
    lines.append("    count,")
    lines.append("};")
    lines.append("")
    lines.append(
        "[[nodiscard]] inline constexpr unsigned to_index(EvolVar value) noexcept {"
    )
    lines.append("    return static_cast<unsigned>(value);")
    lines.append("}")
    lines.append("")
    lines.append("inline constexpr unsigned NUM_EVOL_GFS = to_index(EvolVar::count);")
    lines.append("")
    lines.append(
        "inline constexpr std::array<std::string_view, NUM_EVOL_GFS> EVOL_GF_NAMES = {"
    )
    for fg in evol:
        lines.append(f'    "{fg.name}",')
    lines.append("};")
    lines.append("")
    lines.append("static_assert(NUM_EVOL_GFS == EVOL_GF_NAMES.size());")
    lines.append("")
    # Section 5.5: the emitter also produces rank, asymptotic value and
    # wavespeed metadata, rendered from the frozen gridfunction records.  The
    # host needs these in C++ (outer boundaries need f_infinity, a CFL step
    # needs the wavespeed); leaving them in state.json alone would force a
    # hand-maintained table, which section 5.1 exists to prevent.
    lines.append("inline constexpr std::array<unsigned, NUM_EVOL_GFS> EVOL_GF_RANK = {")
    for fg in evol:
        lines.append(f"    {int(fg.rank)},")
    lines.append("};")
    lines.append("")
    scalar_type = evol[0].gf_type
    lines.append(
        f"inline constexpr std::array<{scalar_type}, NUM_EVOL_GFS>"
        " EVOL_GF_F_INFINITY = {"
    )
    for fg in evol:
        lines.append(
            f"    {_cxx_scalar_literal(fg.f_infinity, fg.name, 'f_infinity')},"
        )
    lines.append("};")
    lines.append("")
    lines.append(
        f"inline constexpr std::array<{scalar_type}, NUM_EVOL_GFS>"
        " EVOL_GF_WAVESPEED = {"
    )
    for fg in evol:
        lines.append(f"    {_cxx_scalar_literal(fg.wavespeed, fg.name, 'wavespeed')},")
    lines.append("};")
    lines.append("")
    control_indices = _upwind_control_indices(snapshot, [fg.name for fg in evol])
    lines.append(
        "// Indices of the evolved fields that drive NRPy's upwind selection"
        " (section 8.5)."
    )
    lines.append(
        "inline constexpr unsigned NUM_UPWIND_CONTROL_GFS = " f"{len(control_indices)};"
    )
    lines.append(
        "inline constexpr std::array<unsigned, NUM_UPWIND_CONTROL_GFS>"
        " EVOL_UPWIND_CONTROL_INDICES = {"
    )
    for index in control_indices:
        lines.append(f"    {index},")
    lines.append("};")
    lines.append("")
    # Step 5: the section 14.7 strict, case-sensitive exact-name lookup.
    # Matching is case-sensitive because NRPy tensor-variance suffixes are, and
    # an unknown name resolves to an empty optional rather than to a silently
    # wrong index.  Reporting an unknown name is the caller's obligation;
    # variable_count and variable_name let it list every valid generated name
    # without a second table.
    by_group = state_records_by_group(snapshot)
    # Whitepaper section 14.7 fixes both the member names and their order:
    # evolved, diagnostic, auxevol, auxiliary.  "evolved" is rendered above.
    groups = (
        ("DIAG", "DIAG_GF_NAMES", "NUM_DIAG_GFS", "diagnostic"),
        ("AUXEVOL", "AUXEVOL_GF_NAMES", "NUM_AUXEVOL_GFS", "auxevol"),
        ("AUX", "AUX_GF_NAMES", "NUM_AUX_GFS", "auxiliary"),
    )
    for group, array_name, count_name, _member in groups:
        names = by_group.get(group, [])
        entries = "".join(f'    "{name}",\n' for name in names)
        lines.append(
            f"inline constexpr unsigned {count_name} = {len(names)};\n"
            f"inline constexpr std::array<std::string_view, {count_name}>"
            f" {array_name} = {{\n"
            f"{entries}"
            "};\n"
        )
    members = "\n".join(f"        {member}," for _g, _a, _c, member in groups)
    count_cases = "\n".join(
        f"        case VariableRef::Group::{member}: return {count_name};"
        for _g, _a, count_name, member in groups
    )
    name_cases = "\n".join(
        f"        case VariableRef::Group::{member}:"
        f" return index < {count_name} ? {array_name}[index]"
        # An empty but non-null view: a null data() reaches printf-family
        # callers, which -Wnonnull / -Wformat-overflow reject under -O2.
        ' : std::string_view{""};'
        for _g, array_name, count_name, member in groups
    )
    group_members = "\n".join(
        f"    VariableRef::Group::{member}," for _g, _a, _c, member in groups
    )
    lines.append(rf"""// Section 14.7: strict, case-sensitive exact-name selection.
struct VariableRef {{
    enum class Group : unsigned {{
        evolved,
{members}
    }} group;  // END ENUM: generated variable groups
    unsigned index;
}};  // END STRUCT: exact-name reference

[[nodiscard]] inline constexpr unsigned variable_count(VariableRef::Group group) noexcept {{
    switch (group) {{
        case VariableRef::Group::evolved: return NUM_EVOL_GFS;
{count_cases}
    }}  // END SWITCH: generated group cardinality
    return 0;
}}  // END FUNCTION: variable_count

[[nodiscard]] inline constexpr std::string_view variable_name(VariableRef::Group group, unsigned index) noexcept {{
    switch (group) {{
        case VariableRef::Group::evolved: return index < NUM_EVOL_GFS ? EVOL_GF_NAMES[index] : std::string_view{{""}};
{name_cases}
    }}  // END SWITCH: generated group name array
    return std::string_view{{""}};
}}  // END FUNCTION: variable_name

inline constexpr std::array<VariableRef::Group, {len(groups) + 1}> VARIABLE_GROUPS = {{
    VariableRef::Group::evolved,
{group_members}
}};

[[nodiscard]] inline constexpr std::optional<VariableRef> find_variable(std::string_view exact_name) noexcept {{
    for (const VariableRef::Group group : VARIABLE_GROUPS) {{
        const unsigned count = variable_count(group);
        for (unsigned index = 0; index < count; ++index) {{
            if (variable_name(group, index) == exact_name) {{
                return VariableRef{{group, index}};
            }}  // END IF: exact name matched
        }}  // END LOOP: for index over group
    }}  // END LOOP: for group over all groups
    return std::nullopt;
}}  // END FUNCTION: find_variable
""")
    lines.append("}  // namespace fccz4::generated")
    lines.append("")
    return "\n".join(lines)


def _cxx_scalar_literal(value: str, gf_name: str, field: str) -> str:
    """
    Render one frozen scalar metadata value as a C++ floating-point literal.

    :param value: The frozen value, as stored in the snapshot record.
    :param gf_name: Gridfunction name (for the error message).
    :param field: Record field name (for the error message).
    :return: A C++ literal.
    :raises ValueError: If the value is not a finite real number, which would
        otherwise emit an uncompilable symbolic expression such as ``sqrt(2)``.
    """
    try:
        number = float(value)
    except (TypeError, ValueError):
        try:
            import sympy as sp

            number = float(sp.sympify(value))
        except (TypeError, ValueError, AttributeError) as exc:
            raise ValueError(
                f"Gridfunction {gf_name!r} has non-numeric {field} {value!r}: "
                "generated metadata must be a finite real number."
            ) from exc
    if number != number or number in (float("inf"), float("-inf")):
        raise ValueError(f"Gridfunction {gf_name!r} has non-finite {field} {value!r}.")
    return repr(number)


def _upwind_control_indices(
    snapshot: FrozenNRPyDendroSnapshot, evol_names: List[str]
) -> Tuple[int, ...]:
    """
    Return the EVOL indices of the builder-recorded upwind control fields.

    The control fields come from the builder sidecar (they are derived from
    the shared expression factory's upwind control vector), so this renderer
    invents nothing; a profile that records none renders an empty array.

    :param snapshot: The frozen snapshot.
    :param evol_names: Ordered EVOL names.
    :return: The control-field indices, in EVOL order.
    :raises ValueError: If a recorded control field is not an EVOL field.
    """
    provenance = dict(snapshot.extras).get("provenance")
    if not isinstance(provenance, dict):
        return ()
    recorded = provenance.get("upwind_control_fields") or ()
    positions = {name: index for index, name in enumerate(evol_names)}
    indices = []
    for name in recorded:
        if name not in positions:
            raise ValueError(
                f"Upwind control field {name!r} is not a registered EVOL field."
            )
        indices.append(positions[name])
    return tuple(sorted(indices))


def render_component_bindings(
    names: Sequence[str],
    scalar_type: str,
    *,
    array: str,
    role: Callable[[str], str],
    const_pointee: bool,
    index_expression: Callable[[str, int], str],
    base_offset: str = "geom.component_offset",
    flat_stride: Optional[str] = None,
) -> str:
    """
    Render one role-prefixed pointer binding per field, in the given order.

    This is the single binding emitter: the frozen-snapshot renderer and the
    per-block/flat-block kernel builders all go through it, so the pointer
    role, the per-component base offset, and the field ordering cannot drift
    between them (whitepaper sections 8.2/9.11).

    :param names: Exact registered gridfunction names, in NRPy registry order.
    :param scalar_type: Generated scalar alias (e.g. ``DendroScalar``).
    :param array: Name of the caller's pointer array (component layout) or of
        the flat base pointer (flat layout).
    :param role: Role-prefix function from :mod:`naming` (e.g.
        :func:`naming.input_pointer`).
    :param const_pointee: Emit ``const <scalar>* const`` rather than
        ``<scalar>* const``.
    :param index_expression: Maps (name, position) to the component index
        expression, so a caller may use either the generated ``EvolVar`` enum
        or the frozen integer position.
    :param base_offset: Per-component base offset expression (section 8.2).
    :param flat_stride: When given, field ``f`` lives at
        ``array + base_offset + f * flat_stride`` (the LTS flat-block layout)
        instead of ``array[f] + base_offset``.
    :return: The C++ binding statements, one per line.
    """
    qualifier = (
        f"const {scalar_type}* const" if const_pointee else f"{scalar_type}* const"
    )
    lines: List[str] = []
    for position, name in enumerate(names):
        index = index_expression(name, position)
        if flat_stride is None:
            source = f"{array}[{index}] + {base_offset}"
        else:
            source = (
                f"{array} + {base_offset}"
                f" + static_cast<std::ptrdiff_t>({index}) * {flat_stride}"
            )
        lines.append(f"{qualifier} {role(name)} = {source};")
    return "\n".join(lines)


def evol_names_and_scalar_type(
    snapshot: FrozenNRPyDendroSnapshot,
) -> Tuple[Tuple[str, ...], str]:
    """
    Return the frozen EVOL names in NRPy order and their common scalar alias.

    :param snapshot: The frozen snapshot.
    :return: (ordered EVOL names, scalar type).
    :raises ValueError: If there are no EVOL fields or their types differ.
    """
    evol = [fg for fg in snapshot.gridfunctions if fg.group == "EVOL"]
    if not evol:
        raise ValueError("Cannot render pointer bindings with no EVOL fields.")
    scalar_types = {fg.gf_type for fg in evol}
    if len(scalar_types) != 1:
        raise ValueError(
            f"Refusing to render pointer bindings with mixed scalar types: "
            f"{sorted(scalar_types)}."
        )
    return tuple(fg.name for fg in evol), evol[0].gf_type


def render_state_json(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the canonical state.json manifest.

    Records are rendered in snapshot order, which is NRPy group-list order
    (whitepaper section 12.5: preserve NRPy's returned gridfunction order
    exactly), with sorted JSON keys.

    :param snapshot: The frozen snapshot.
    :return: Canonical JSON text.
    """
    records = [
        {
            "group": fg.group,
            "index": fg.index,
            "name": fg.name,
            "rank": fg.rank,
            "dimension": fg.dimension,
            "gf_type": fg.gf_type,
            "f_infinity": fg.f_infinity,
            "wavespeed": fg.wavespeed,
            "is_basename": fg.is_basename,
            "description": fg.description,
        }
        for fg in snapshot.gridfunctions
    ]
    payload = {
        "state_schema_hash": snapshot.hashes.state_schema_hash,
        "records": records,
    }
    return json.dumps(payload, sort_keys=True, indent=2) + "\n"


def state_name_index(snapshot: FrozenNRPyDendroSnapshot) -> Dict[str, Tuple[str, int]]:
    """
    Return a mapping of exact name to (group, index) for all frozen fields.

    :param snapshot: The frozen snapshot.
    :return: Mapping of exact NRPy name to (group, index).
    """
    return {fg.name: (fg.group, fg.index) for fg in snapshot.gridfunctions}


def render_fccz4_provenance_hpp(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the generated constants and hash header (fccz4_provenance.hpp).

    Carries the generated FD order, required padding, minimum element order
    under the current Dendro rule (section 8.6), and the independent semantic
    hashes (section 9.10) as C++ constants, all rendered from the frozen
    snapshot (section 9.11: state cardinality and hash constants).

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: The complete C++ header text.
    """
    gen = dict(snapshot.generation_parameters)
    padding = list(snapshot.required_padding)
    lines: List[str] = []
    lines.append(
        _banner(snapshot.hashes.module_abi_hash, "frozen snapshot constants and hashes")
    )
    lines.append("#pragma once")
    lines.append("")
    lines.append("namespace fccz4::generated {")
    lines.append("")
    lines.append(
        "inline constexpr unsigned FD_ORDER = " f"{int(str(gen.get('fd_order', 0)))};"
    )
    lines.append(f"inline constexpr unsigned REQUIRED_PADDING_X = {padding[0]};")
    lines.append(f"inline constexpr unsigned REQUIRED_PADDING_Y = {padding[1]};")
    lines.append(f"inline constexpr unsigned REQUIRED_PADDING_Z = {padding[2]};")
    lines.append(f"inline constexpr unsigned REQUIRED_PADDING = {max(padding)};")
    lines.append(
        "inline constexpr unsigned MINIMUM_ELEMENT_ORDER_UNDER_CURRENT_RULE = "
        f"{2 * max(padding)};"
    )
    lines.append("")
    lines.append(
        f'inline constexpr char const* MODULE_ABI_HASH = "{snapshot.hashes.module_abi_hash}";'
    )
    lines.append(
        f'inline constexpr char const* STATE_SCHEMA_HASH = "{snapshot.hashes.state_schema_hash}";'
    )
    lines.append(
        f'inline constexpr char const* PARAMETER_SCHEMA_HASH = "{snapshot.hashes.parameter_schema_hash}";'
    )
    lines.append(
        f'inline constexpr char const* EQUATION_HASH = "{snapshot.hashes.equation_hash}";'
    )
    lines.append(
        f'inline constexpr char const* STENCIL_HASH = "{snapshot.hashes.stencil_hash}";'
    )
    lines.append(
        f'inline constexpr char const* CFUNCTION_API_HASH = "{snapshot.hashes.cfunction_api_hash}";'
    )
    lines.append(
        f'inline constexpr char const* CFUNCTION_SOURCE_HASH = "{snapshot.hashes.cfunction_source_hash}";'
    )
    lines.append("")
    lines.append("}  // namespace fccz4::generated")
    lines.append("")
    return "\n".join(lines)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
