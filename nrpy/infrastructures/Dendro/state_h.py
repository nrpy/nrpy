# nrpy/infrastructures/Dendro/state_h.py
"""
Emit the generated state header for a Dendro solver.

Every name, order, count and metadata value is read from the NRPy gridfunction
registry, as BHaH's ``BHaH_defines_h`` reads ``par.glb_code_params_dict``.  This
module contains no field list and imports no equation module.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import math
from typing import List, Tuple, cast

import sympy as sp

import nrpy.grid as gri
from nrpy.infrastructures.Dendro import CFunction_roles as roles
from nrpy.infrastructures.Dendro import gridfunction_name_decorations as gf_names
from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.header_guards import header_guard

BANNER = generated_file_banner()


def state_records() -> List[Tuple[str, int, str]]:
    """
    Return registered gridfunctions as ``(group, index_in_group, name)``.

    Both membership and order come from
    :func:`nrpy.grid.GridFunction.gridfunction_lists`, which is NRPy's single
    ordering authority.  Every Dendro builder indexes its component arrays from
    that same list, so the emitted enum, metadata arrays and kernel pointer
    bindings cannot disagree about which component is which.

    :return: One record per registered gridfunction.

    Doctests:
    >>> import nrpy.params as par
    >>> _saved_fields = dict(gri.glb_gridfcs_dict)
    >>> try:
    ...     gri.glb_gridfcs_dict.clear()
    ...     par.set_parval_from_str("Infrastructure", "Dendro")
    ...     _ = gri.register_gridfunctions(["bXX", "aYY"], group="EVOL")
    ...     _ = gri.register_gridfunctions_for_single_rank1(
    ...         "vU", dimension=2, group="AUXEVOL"
    ...     )
    ...     _ = gri.register_gridfunctions("zAux", group="AUX")
    ...     _ = gri.register_gridfunctions("hZZ", group="DIAG")
    ...     assert state_records() == [
    ...         ('EVOL', 0, 'aYY'), ('EVOL', 1, 'bXX'),
    ...         ('AUXEVOL', 0, 'vU0'), ('AUXEVOL', 1, 'vU1'),
    ...         ('AUX', 0, 'zAux'), ('DIAG', 0, 'hZZ')]
    ...     assert group_names("EVOL") == ['aYY', 'bXX']
    ... finally:
    ...     gri.glb_gridfcs_dict.clear()
    ...     gri.glb_gridfcs_dict.update(_saved_fields)
    >>> all(gri.glb_gridfcs_dict[name] is value for name, value in _saved_fields.items())
    True
    """
    evol, auxevol, diag, aux = gri.GridFunction.gridfunction_lists()
    records: List[Tuple[str, int, str]] = []
    for group, names in (
        ("EVOL", evol),
        ("AUXEVOL", auxevol),
        ("AUX", aux),
        ("DIAG", diag),
    ):
        for index, name in enumerate(names):
            records.append((group, index, name))
    return records


def group_names(group: str) -> List[str]:
    """
    Return the registered gridfunction names in one group, in registry order.

    Read the core gridfunction registry at the point of use.

    :param group: Registry group name, one of ``EVOL``, ``AUXEVOL``, ``AUX``, ``DIAG``.
    :return: The group's gridfunction names, in NRPy registry order.
    :raises ValueError: If the group name is not a registered Dendro group.
    """
    evol, auxevol, diag, aux = gri.GridFunction.gridfunction_lists()
    names_by_group = {
        "EVOL": evol,
        "AUXEVOL": auxevol,
        "AUX": aux,
        "DIAG": diag,
    }
    if group not in names_by_group:
        raise ValueError(
            f"Unknown Dendro registry group {group!r}; expected one of "
            f"{', '.join(sorted(names_by_group))}."
        )
    return list(names_by_group[group])


def output_state_h(solver_stem: str, solver_namespace: str) -> str:
    """
    Emit the generated state header.

    It carries the EVOL enum, name array, per-field metadata, and the
    exact-name lookup over every registered group.  Every name, count and
    metadata value comes from the gridfunction registry; no field name is
    hardcoded here.

    :param solver_stem: Lowercase formulation stem for emitted header names.
    :param solver_namespace: Solver namespace, following Dendro's lowercase
        formulation habit (``namespace bssn``).
    :return: The complete C++ header text.
    :raises ValueError: If a recorded upwind control field is not an EVOL field.

    Doctests:
    >>> import nrpy.params as par
    >>> gri.glb_gridfcs_dict.clear()
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> _ = gri.register_gridfunctions(["bXX", "aYY"], group="EVOL")
    >>> par.glb_extras_dict.get("Dendro", {}).pop("upwind_control_fields", None) and None
    >>> try:
    ...     output_state_h("bssn", "bssn")
    ... except ValueError as error:
    ...     print(str(error).splitlines()[0])
    No Dendro kernel has recorded an upwind control set; register the right-hand-side CFunctions before emitting the state header.
    >>> roles.set_upwind_control_fields(("bXX",))
    >>> header = output_state_h("bssn", "bssn")
    >>> "#ifndef BSSN_STATE_H" in header
    True
    >>> header.rstrip().endswith("#endif  // BSSN_STATE_H")
    True
    >>> "inline constexpr unsigned NUM_UPWIND_CONTROL_GFS = 1;" in header
    True
    >>> "enum class EvolVar : unsigned {" in header
    True
    >>> [line.strip() for line in header.splitlines() if line.strip().endswith(("= 0,", "= 1,"))]
    ['aYY = 0,', 'bXX = 1,']
    >>> from nrpy.helpers.generic import clang_format
    >>> "}  // END NAMESPACE: bssn::generated" in clang_format(header)
    True
    >>> '#include "bssn_types.h"' in header
    True
    """
    evol = [
        (index, name, gri.glb_gridfcs_dict[name])
        for group, index, name in state_records()
        if group == "EVOL"
    ]
    lines: List[str] = []
    opening, closing = header_guard(f"{solver_stem}_state.h")
    lines.append(BANNER.rstrip("\n"))
    lines.append(opening)
    lines.append("")
    lines.append("#include <array>")
    lines.append("#include <cstddef>")
    lines.append("#include <optional>")
    lines.append("#include <string_view>")
    lines.append("")
    lines.append(f'#include "{solver_stem}_types.h"')
    lines.append("")
    lines.append(f"namespace {solver_namespace}::generated {{")
    lines.append("")
    lines.append("enum class EvolVar : unsigned {")
    for index, name, _gf in evol:
        lines.append(f"    {gf_names.enum_member(name)} = {index},")
    lines.append("    END,")
    lines.append("};  // END ENUM: EvolVar")
    lines.append("")
    lines.append(
        "[[nodiscard]] inline constexpr unsigned to_index(EvolVar value) noexcept {"
    )
    lines.append("    return static_cast<unsigned>(value);")
    lines.append("}")
    lines.append("")
    lines.append("inline constexpr unsigned NUM_EVOL_GFS = to_index(EvolVar::END);")
    lines.append("")
    lines.append(
        "inline constexpr std::array<std::string_view, NUM_EVOL_GFS> EVOL_GF_NAMES = {"
    )
    for _index, name, _gf in evol:
        lines.append(f'    "{name}",')
    lines.append("};  // END ARRAY: EVOL_GF_NAMES")
    lines.append("")
    lines.append("static_assert(NUM_EVOL_GFS == EVOL_GF_NAMES.size());")
    lines.append("")
    # The emitter also produces rank, asymptotic value and
    # wavespeed metadata, rendered from the registered gridfunction records.  The
    # host needs these in C++ (outer boundaries need f_infinity, a CFL step
    # needs the wavespeed); anything less would force a hand-maintained
    # table that could drift from the registry.
    lines.append("inline constexpr std::array<unsigned, NUM_EVOL_GFS> EVOL_GF_RANK = {")
    for _index, _name, gf in evol:
        lines.append(f"    {int(gf.rank)},")
    lines.append("};  // END ARRAY: EVOL_GF_RANK")
    lines.append("")
    scalar_type = gri.DENDRO_SCALAR_TYPE
    lines.append(
        f"inline constexpr std::array<{scalar_type}, NUM_EVOL_GFS>"
        " EVOL_GF_F_INFINITY = {"
    )
    for _index, name, gf in evol:
        lines.append(
            f"    {_cxx_scalar_literal(str(gf.f_infinity), name, 'f_infinity')},"
        )
    lines.append("};  // END ARRAY: EVOL_GF_F_INFINITY")
    lines.append("")
    lines.append(
        f"inline constexpr std::array<{scalar_type}, NUM_EVOL_GFS>"
        " EVOL_GF_WAVESPEED = {"
    )
    for _index, name, gf in evol:
        lines.append(
            f"    {_cxx_scalar_literal(str(gf.wavespeed), name, 'wavespeed')},"
        )
    lines.append("};  // END ARRAY: EVOL_GF_WAVESPEED")
    lines.append("")
    # The control fields come from the right-hand-side builder, which derives
    # them from the shared expression factory's upwind control vector, so this
    # renderer invents nothing.  Reading them through
    # roles.upwind_control_fields() means an unrecorded set raises here rather
    # than emitting an empty table that would silently disable the generated
    # upwind self-test.
    evol_positions = {name: index for index, (_i, name, _g) in enumerate(evol)}
    control_index_list: List[int] = []
    for control_name in roles.upwind_control_fields():
        if control_name not in evol_positions:
            raise ValueError(
                f"Upwind control field {control_name!r} is not a registered EVOL field."
            )
        control_index_list.append(evol_positions[control_name])
    control_indices = tuple(sorted(control_index_list))
    lines.append("// Indices of the evolved fields that drive NRPy's upwind selection.")
    lines.append(
        "inline constexpr unsigned NUM_UPWIND_CONTROL_GFS = " f"{len(control_indices)};"
    )
    lines.append(
        "inline constexpr std::array<unsigned, NUM_UPWIND_CONTROL_GFS>"
        " EVOL_UPWIND_CONTROL_INDICES = {"
    )
    for index in control_indices:
        lines.append(f"    {index},")
    lines.append("};  // END ARRAY: EVOL_UPWIND_CONTROL_INDICES")
    lines.append("")
    # The strict, case-sensitive exact-name lookup.
    # Matching is case-sensitive because NRPy tensor-variance suffixes are, and
    # an unknown name resolves to an empty optional rather than to a silently
    # wrong index.  Reporting an unknown name is the caller's obligation;
    # variable_count and variable_name let it list every valid generated name
    # without a second table.
    by_group = {g: group_names(g) for g in ("EVOL", "AUXEVOL", "AUX", "DIAG")}
    # The member names and their order are fixed:
    # evolved, diagnostic, auxevol, auxiliary.  "evolved" is rendered above.
    groups = (
        ("DIAG", "DIAG_GF_NAMES", "NUM_DIAG_GFS", "DIAG"),
        ("AUXEVOL", "AUXEVOL_GF_NAMES", "NUM_AUXEVOL_GFS", "AUXEVOL"),
        ("AUX", "AUX_GF_NAMES", "NUM_AUX_GFS", "AUX"),
    )
    for group, array_name, count_name, _member in groups:
        names = by_group.get(group, [])
        entries = "".join(f'    "{name}",\n' for name in names)
        lines.append(
            f"inline constexpr unsigned {count_name} = {len(names)};\n"
            f"inline constexpr std::array<std::string_view, {count_name}>"
            f" {array_name} = {{\n"
            f"{entries}"
            f"}};  // END ARRAY: {array_name}\n"
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
    lines.append(rf"""// Strict, case-sensitive exact-name selection.
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
}};  // END ARRAY: VARIABLE_GROUPS

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
    lines.append("// clang-format off")
    lines.append(f"}}  // END NAMESPACE: {solver_namespace}::generated")
    lines.append("// clang-format on")
    lines.append("")
    lines.append(closing)
    lines.append("")
    return "\n".join(lines)


def _cxx_scalar_literal(value: str, gf_name: str, field: str) -> str:
    """
    Render one registered scalar metadata value as a C++ floating-point literal.

    :param value: The registered metadata value.
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
            number = float(cast(sp.Expr, sp.sympify(value)))
        except (TypeError, ValueError, AttributeError) as exc:
            raise ValueError(
                f"Gridfunction {gf_name!r} has non-numeric {field} {value!r}: "
                "generated metadata must be a finite real number."
            ) from exc
    if not math.isfinite(number):
        raise ValueError(f"Gridfunction {gf_name!r} has non-finite {field} {value!r}.")
    return repr(number)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
