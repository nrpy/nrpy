# nrpy/infrastructures/Dendro/CodeParameters_output.py
"""
Pure CodeParameter translator over the frozen snapshot.

Every name, type, default, and consumer record is rendered from the frozen
``FrozenCodeParameter`` records; no physics parameter table is authored here
and no equation module is imported (whitepaper sections 6.2/6.5/9.12).
"""

import json
from typing import Any, List

from nrpy.infrastructures.Dendro.freeze import FrozenNRPyDendroSnapshot

_BANNER_BODY = (
    "Generator: nrpy.infrastructures.Dendro\n"
    "Regenerate with: python -m nrpy.examples.dendro_fccz4\n"
)


def _banner(hash_line: str, source_object: str) -> str:
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


def c_type(cparam_type: str, scalar_type: str) -> str:
    """
    Map a registered CodeParameter type to its generated C++ type.

    Whitepaper section 6.5: ``REAL`` maps to the registered Dendro scalar
    alias; ``float``/``double``/``int``/``bool`` map to the same spelled
    built-in types; ``REAL[n]``/``int[n]`` map to ``std::array``.
    ``#define`` parameters are compile-time constants and are rendered as
    ``inline constexpr``.  Unsupported types fail generation.

    :param cparam_type: Registered ``cparam_type`` string.
    :param scalar_type: Frozen Dendro scalar alias (e.g. ``DendroScalar``).
    :return: The generated C++ type spelling.
    :raises ValueError: If the type has no section 6.5 mapping.
    """
    if cparam_type == "REAL":
        return scalar_type
    if cparam_type in ("float", "double", "int", "bool"):
        return cparam_type
    if cparam_type in ("DendroScalar",):
        return scalar_type
    if cparam_type == "#define":
        raise ValueError(
            "`#define` CodeParameters are compile-time constants: emit them as "
            "`inline constexpr`, never as a `GeneratedParams` member (a "
            "placeholder `auto` non-static data member is ill-formed C++17)."
        )
    if cparam_type.startswith(("REAL[", "int[")):
        raise ValueError(
            f"CodeParameter type {cparam_type!r} has no implemented renderer: "
            "array parameters need element-wise defaults, validation and "
            "printing before they can be generated (section 6.5)."
        )
    raise ValueError(f"CodeParameter type {cparam_type!r} has no section 6.5 mapping.")


def _toml_value(cparam_type: str, default_value: Any) -> Any:
    """
    Render the TOML sample value for a CodeParameter default.

    :param cparam_type: Registered ``cparam_type`` string.
    :param default_value: Registered default value.
    :return: A JSON-serializable value for the sample TOML entry.
    :raises ValueError: If the type has no TOML mapping.
    """
    if cparam_type in ("int",):
        return int(default_value)
    if cparam_type == "bool":
        return bool(default_value)
    if cparam_type in ("REAL", "float", "double", "DendroScalar"):
        return float(default_value)
    raise ValueError(f"No TOML mapping for cparam_type {cparam_type!r}.")


def render_fccz4_parameters_hpp(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the generated parameter struct header (fccz4_parameters.hpp).

    The struct holds one member per used frozen CodeParameter in sorted
    name order, typed per section 6.5.  No default is baked into the struct:
    defaults live in the registered CFunction ``fccz4_params_set_defaults``
    and in the sample TOML.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: The complete C++ header text.
    """
    scalar_type = str(
        dict(snapshot.generation_parameters).get("Dendro_scalar_type", "DendroScalar")
    )
    lines: List[str] = []
    lines.append(
        _banner(snapshot.hashes.module_abi_hash, "frozen NRPy CodeParameter registry")
    )
    lines.append("#pragma once")
    lines.append("")
    lines.append("#include <array>")
    lines.append("")
    lines.append('#include "fccz4_types.hpp"')
    lines.append("")
    lines.append("namespace fccz4::generated {")
    lines.append("")
    lines.append("struct GeneratedParams {")
    for cp in snapshot.codeparameters:
        lines.append(f"    {c_type(cp.cparam_type, scalar_type)} {cp.name};")
    lines.append("};")
    lines.append("")
    lines.append(
        "inline constexpr unsigned NUM_CODE_PARAMETERS = "
        f"{len(snapshot.codeparameters)};"
    )
    lines.append("")
    lines.append("}  // namespace fccz4::generated")
    lines.append("")
    return "\n".join(lines)


def render_parameters_json(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the canonical parameters.json manifest.

    Records are rendered in snapshot (sorted name) order with sorted JSON
    keys (whitepaper section 12.5 canonical JSON).

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: Canonical JSON text.
    """
    records = [
        {
            "name": cp.name,
            "module": cp.module,
            "cparam_type": cp.cparam_type,
            "default_value": cp.default_value,
            "description": cp.description,
            "commondata": cp.commondata,
            "add_to_parfile": cp.add_to_parfile,
            "consumers": list(cp.consumers),
        }
        for cp in snapshot.codeparameters
    ]
    payload = {
        "parameter_schema_hash": snapshot.hashes.parameter_schema_hash,
        "records": records,
    }
    return json.dumps(payload, sort_keys=True, indent=2) + "\n"


def render_parameter_toml_sample(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the sample TOML section for the used runtime parameters.

    Only parameters with ``add_to_parfile`` appear; every value is the
    registered default (sections 6.2/6.5).  The section is embedded into the
    generated pars TOML files.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: TOML text (with trailing newline), empty if no parfile params.
    """
    lines: List[str] = ["[fccz4.params]"]
    for cp in snapshot.codeparameters:
        if not cp.add_to_parfile:
            continue
        lines.append(f"{cp.name} = {_format_toml_scalar(cp)}")
    return "\n".join(lines) + "\n"


def _format_toml_scalar(cp: Any) -> str:
    """
    Format one TOML scalar from a frozen CodeParameter record.

    :param cp: The frozen CodeParameter record.
    :return: The TOML literal.
    """
    value = _toml_value(cp.cparam_type, cp.default_value)
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    return repr(value)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
