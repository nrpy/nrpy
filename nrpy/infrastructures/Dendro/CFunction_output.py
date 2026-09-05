# nrpy/infrastructures/Dendro/CFunction_output.py
"""
Pure CFunction translator: emit one source file per registered CFunction.

Every artifact is derived from the frozen ``FrozenCFunction`` records
(whitepaper section 9.13): the registered ``full_function`` is written as the
file content (plus the deterministic generated-file banner), declarations are
emitted from the registered prototypes, and the CMake source list is derived
from the same objects.  No signature is reconstructed, no loop is wrapped,
and no registered name is changed.
"""

import os
from typing import Dict, List, Tuple

from nrpy.infrastructures.Dendro.freeze import FrozenCFunction, FrozenNRPyDendroSnapshot

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


def derived_source_path(fc: FrozenCFunction) -> str:
    """
    Derive the generated source path for one frozen CFunction.

    The registered subdirectory is already module-root-relative and includes
    the ``generated/src`` prefix (e.g. ``generated/src/rhs``), so the output
    path is exactly ``<registered-subdirectory>/<name>.cpp`` (section 9.13).

    :param fc: The frozen CFunction record.
    :return: Relative path from the module root (forward slashes).
    """
    subdirectory = fc.subdirectory.replace(os.sep, "/")
    return f"{subdirectory}/{fc.name}.cpp"


def render_CFunction_source(
    snapshot: FrozenNRPyDendroSnapshot, fc: FrozenCFunction
) -> str:
    """
    Render one generated CFunction source file.

    The registered full_function is emitted verbatim after the banner and a
    single generated preamble include that supplies the host types (the mock
    vehicle's Dendro types; frozen after the I0-1 gates), the generated
    scalar/state/parameter registries, the CFunction declarations, UPWIND_ALG,
    and the standard headers the body uses.  The preamble itself contains no
    concrete field name, no FD coefficient, and no numerical loop
    (section 9.16).

    :param snapshot: The frozen NRPy Dendro snapshot.
    :param fc: The frozen CFunction record to emit.
    :return: The complete source file text (banner + preamble include + full_function).
    """
    return (
        _banner(snapshot.hashes.module_abi_hash, f"registered CFunction {fc.name}")
        + "\n"
        + '#include "generated_project_preamble.hpp"\n'
        + fc.full_function
        + "\n"
    )


def render_CFunction_declarations(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the generated CFunction declaration header.

    One declaration per frozen CFunction, from the registered prototype
    (section 9.13: declarations from the same objects).

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: The complete C++ header text.
    """
    lines: List[str] = []
    lines.append(
        _banner(snapshot.hashes.module_abi_hash, "registered CFunction prototypes")
    )
    lines.append("#pragma once")
    lines.append("")
    lines.append('#include "fccz4_types.hpp"')
    lines.append("")
    # Host types appear in registered prototypes (the mock-vehicle names;
    # the production host types freeze after the I0-1 Dendrolib gates,
    # whitepaper sections 12.5/13.3).  Forward declarations keep this header
    # self-contained; the consuming translation unit includes the host header
    # for the complete types.
    lines.append("struct BlockGeometry;")
    lines.append("struct MockWorld;")
    lines.append("")
    for fc in snapshot.cfunctions:
        lines.append(fc.prototype)
    lines.append("")
    return "\n".join(lines)


def CFunction_artifacts(snapshot: FrozenNRPyDendroSnapshot) -> Dict[str, str]:
    """
    Render all CFunction artifacts keyed by module-relative path.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: Mapping of ``Dendro-GR/FCCZ4_GR/<path>`` to file text.
    """
    artifacts: Dict[str, str] = {}
    for fc in snapshot.cfunctions:
        rel = derived_source_path(fc)
        artifacts[f"Dendro-GR/FCCZ4_GR/{rel}"] = render_CFunction_source(snapshot, fc)
    artifacts["Dendro-GR/FCCZ4_GR/generated/include/fccz4_cfunctions.hpp"] = (
        render_CFunction_declarations(snapshot)
    )
    return artifacts


def CFunction_cmake_source_list(snapshot: FrozenNRPyDendroSnapshot) -> Tuple[str, ...]:
    """
    Derive the generated CMake source list from the frozen CFunctions.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: Module-relative source paths, in frozen (sorted name) order.
    """
    return tuple(derived_source_path(fc) for fc in snapshot.cfunctions)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
