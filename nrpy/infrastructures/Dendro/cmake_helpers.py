# nrpy/infrastructures/Dendro/cmake_helpers.py
"""
CMake emitters for a generated Dendro solver.

Every emitted name derives from the ``solver_name``, ``solver_prefix`` and
``exec_or_library_name`` arguments threaded from the calling example, matching
how BHaH threads ``project_name`` and ETLegacy threads ``thorn_name``.

NRPy provenance appears in each generated module directory and CMake project.
Names inside the module retain the formulation stem and required Dendro target
names, such as ``BSSN_ENABLE_CUDA``, ``bssn_common``, and ``bssnCtx.cpp``.

No field list, physics parameter, finite-difference coefficient or numerical
loop appears here; the generated source list is read from the CFunction
registry.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import os
from typing import Dict, List, NamedTuple, Sequence, Tuple

import nrpy.c_function as cfc
from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.header_guards import header_guard

BANNER = generated_file_banner("#").splitlines()
SRC = "${CMAKE_CURRENT_SOURCE_DIR}"


class ModuleLayout(NamedTuple):
    """
    Relative paths of the pieces of a generated Dendro-GR solver module.

    A generated solver is a module inside a Dendro-GR checkout, so its paths are
    fixed by Dendro-GR's own layout rather than by NRPy.  The examples assemble
    the project themselves, as every established infrastructure's example does,
    and take the layout from here so the two formulations cannot drift.

    :param root: The module root, ``Dendro-GR/<solver_name>/``.
    :param generated_include: Generated headers directory.
    :param generated_cmake: Generated CMake fragment directory.
    :param include: Hand-assembled headers directory.
    :param src: Hand-assembled sources directory.
    :param pars: Parameter-file directory.
    :param tests: Generated self-test directory.
    """

    root: str
    generated_include: str
    generated_cmake: str
    include: str
    src: str
    pars: str
    tests: str


def module_layout(solver_name: str) -> ModuleLayout:
    """
    Return the relative paths of one generated solver module's pieces.

    :param solver_name: NRPy module directory name, e.g. ``nrpy_bssn``.
    :return: The module layout.

    Doctests:
    >>> layout = module_layout("nrpy_fccz4")
    >>> layout.root
    'Dendro-GR/nrpy_fccz4/'
    >>> layout.generated_include
    'Dendro-GR/nrpy_fccz4/generated/include/'
    >>> layout.tests
    'Dendro-GR/nrpy_fccz4/tests/'
    """
    root = f"Dendro-GR/{solver_name}/"
    return ModuleLayout(
        root=root,
        generated_include=root + "generated/include/",
        generated_cmake=root + "generated/cmake/",
        include=root + "include/",
        src=root + "src/",
        pars=root + "pars/",
        tests=root + "tests/",
    )


def registered_CFunctions() -> List[Tuple[str, cfc.CFunction]]:
    """
    Return the registered CFunctions in sorted name order.

    :return: ``(name, CFunction)`` pairs sorted by registered name.
    """
    return sorted(cfc.CFunction_dict.items())


def derived_source_path(name: str, subdirectory: str) -> str:
    """
    Derive the generated source path for one registered CFunction.

    The registered subdirectory is already solver-root-relative and includes the
    ``generated/src`` prefix, so the output path is exactly
    ``<registered-subdirectory>/<name>.cpp``.

    :param name: Registered CFunction name.
    :param subdirectory: Registered subdirectory for the CFunction.
    :return: Path relative to the solver root, with forward slashes.

    Doctests:
    >>> derived_source_path("bssn_rhs_eval", "generated/src/rhs_eval")
    'generated/src/rhs_eval/bssn_rhs_eval.cpp'
    """
    return f"{subdirectory.replace(os.sep, '/')}/{name}.cpp"


def output_function_prototypes_h(solver_stem: str) -> str:
    """
    Emit the generated CFunction declaration header.

    One declaration per registered CFunction, from the registered prototype.

    :param solver_stem: Lowercase formulation stem for the emitted header names,
        following Dendro's own habit of naming solver files for the formulation
        (``bssnCtx.h``), never for the framework.
    :return: The complete C++ header text.
    """
    opening, closing = header_guard(f"{solver_stem}_function_prototypes.h")
    lines: List[str] = [generated_file_banner().rstrip("\n"), opening, ""]
    lines.append(f'#include "{solver_stem}_types.h"')
    lines.append("")
    # Host types appear in registered prototypes.  Forward declarations keep the
    # header self-contained; the consuming translation unit includes the host
    # header for the complete types.
    lines.append("struct block_geometry_struct;")
    lines.append("struct standalone_host_mesh_struct;")
    lines.append("")
    for _name, cfunc in registered_CFunctions():
        lines.append(cfunc.function_prototype)
    lines.append("")
    lines.append(closing)
    lines.append("")
    return "\n".join(lines)


def output_CFunctions_function_prototypes_and_construct_CMakeLists(
    solver_name: str,
    solver_stem: str,
    solver_prefix: str,
    exec_or_library_name: str,
    test_sections: Sequence[str],
    standalone_application_ctest: Sequence[str],
    real_application_ctest: Sequence[str],
    real_host_available: bool,
) -> Dict[str, str]:
    """
    Emit the CFunction sources, the prototypes header and the CMake files.

    Like
    ``BHaH/Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile``,
    this function emits one source per registered CFunction, the declaration
    header, and the files needed to build them.

    :param solver_name: NRPy module directory name, e.g. ``nrpy_bssn``.
    :param solver_stem: Lowercase formulation stem for emitted header names.
    :param solver_prefix: Bare formulation prefix for the CMake variables.
    :param exec_or_library_name: Name of the solver executable target.
    :param test_sections: Explicit application-owned self-test section names.
    :param standalone_application_ctest: Standalone application CTest lines.
    :param real_application_ctest: Real-host application CTest lines.
    :param real_host_available: Whether this formulation has qualified real-host
        integration.
    :return: Mapping of ``Dendro-GR/<solver_name>/<path>`` to file text.
    """
    layout = module_layout(solver_name)
    prefix = layout.root
    artifacts: Dict[str, str] = {}
    # One raw source file per registered CFunction; the project writer applies
    # clang-format once alongside the other generated C/C++ source files.  The umbrella
    # include of ``<stem>_defines.h`` rides on the CFunction's own ``includes``
    # field, which supplies the host types, the generated state and parameter
    # declarations, and the standard headers the body uses -- the role
    # ``BHaH_defines.h`` plays for every BHaH CFunction.
    for name, cfunc in registered_CFunctions():
        artifacts[prefix + derived_source_path(name, cfunc.subdirectory)] = (
            cfunc.raw_function
        )
    artifacts[prefix + f"generated/include/{solver_stem}_function_prototypes.h"] = (
        output_function_prototypes_h(solver_stem)
    )
    # This module also emits three CMake files.  Together with the per-CFunction
    # sources and prototypes header, they provide the same build inputs that
    # BHaH's Makefile_helpers provides for a Makefile.
    artifacts[layout.generated_cmake + "generated_sources.cmake"] = (
        output_generated_sources_cmake(solver_prefix)
    )
    artifacts[layout.tests + "CMakeLists.txt"] = output_tests_cmake(
        solver_prefix, solver_stem, test_sections
    )
    artifacts[layout.root + "CMakeLists.txt"] = output_solver_cmake(
        solver_name,
        solver_prefix,
        solver_stem,
        exec_or_library_name,
        standalone_application_ctest,
        real_application_ctest,
        real_host_available,
    )
    return artifacts


def CFunction_cmake_source_list() -> Tuple[str, ...]:
    """
    Derive the generated CMake source list from the registered CFunctions.

    :return: Solver-relative source paths, in sorted name order.
    """
    return tuple(
        derived_source_path(name, cfunc.subdirectory)
        for name, cfunc in registered_CFunctions()
    )


def output_generated_sources_cmake(solver_prefix: str) -> str:
    """
    Emit ``generated/cmake/generated_sources.cmake``.

    One entry per registered CFunction, in sorted name order.

    :param solver_prefix: Bare formulation prefix for CMake variables, e.g. ``BSSN``.
    :return: The CMake file text.
    """
    lines: List[str] = list(BANNER) + [
        f"# {solver_prefix}_MODULE_ROOT is set by the including scope so this list",
        "# is correct from either the solver or the tests directory.",
        f"set({solver_prefix}_NRPY_GENERATED_SOURCES",
    ]
    for path in CFunction_cmake_source_list():
        lines.append("  ${" + solver_prefix + "_MODULE_ROOT}/" + path)
    lines += [")", ""]
    return "\n".join(lines)


def output_solver_cmake(
    solver_name: str,
    solver_prefix: str,
    solver_stem: str,
    exec_or_library_name: str,
    standalone_application_ctest: Sequence[str],
    real_application_ctest: Sequence[str],
    real_host_available: bool,
) -> str:
    """
    Emit the generated solver's ``CMakeLists.txt``.

    :param solver_name: NRPy module directory and CMake project name, e.g. ``nrpy_bssn``.
    :param solver_prefix: Bare formulation prefix for CMake variables, e.g. ``BSSN``.
    :param solver_stem: Lowercase stem the emitters use for solver file names.
    :param exec_or_library_name: Name of the solver executable target.
    :param standalone_application_ctest: Standalone application CTest lines.
    :param real_application_ctest: Real-host application CTest lines.
    :param real_host_available: Whether to emit the selectable real-host branch.
    :return: The CMake file text.

    A formulation without qualified real-host support has no selectable host
    mode or real-host tests:

    >>> standalone_only = output_solver_cmake(
    ...     "nrpy_bssn", "BSSN", "bssn", "bssnSolver",
    ...     ("add_test(NAME bssn_standalone)",), (),
    ...     real_host_available=False,
    ... )
    >>> "BSSN_STANDALONE_HOST" in standalone_only
    False
    >>> "bssn_real_minkowski" in standalone_only
    False

    A qualified formulation retains the selectable real-host branch:

    >>> selectable = output_solver_cmake(
    ...     "fCCZ4", "FCCZ4", "fccz4", "fccz4Solver",
    ...     ("add_test(NAME fccz4_standalone)",),
    ...     ("add_test(NAME fccz4_real_minkowski)",),
    ...     real_host_available=True,
    ... )
    >>> "option(FCCZ4_STANDALONE_HOST" in selectable
    True
    >>> "fccz4_real_minkowski" in selectable
    True

    """
    # The file names come from solver_stem, which the examples also use to
    # name the emitted sources; deriving a second stem from solver_prefix
    # would emit a target whose sources do not exist whenever the two differ.
    stem = solver_stem
    lines: List[str] = list(BANNER) + [
        "cmake_minimum_required(VERSION 3.13)",
        f"project({solver_name} CXX)",
        "",
        "# The solver runs under MPI so the checks can compare all ranks.",
        "# MPI is a standard host dependency, not Dendrolib; the real Dendro-GR",
        "# build also provides it.",
        "find_package(MPI REQUIRED COMPONENTS CXX)",
        "",
    ]
    if real_host_available:
        lines += [
            f'option({solver_prefix}_STANDALONE_HOST "Build the solver standalone, against the NRPy-supplied Dendro host declarations rather than a real Dendro-GR build" ON)',
            "",
        ]
    lines += [
        f"set({solver_prefix}_MODULE_ROOT {SRC})",
        f"include({SRC}/generated/cmake/generated_sources.cmake)",
        "",
        f"add_library({stem}_common OBJECT",
        f"  src/{stem}Ctx.cpp",
        "  ${" + solver_prefix + "_NRPY_GENERATED_SOURCES}",
        ")",
        "",
        f"target_compile_features({stem}_common PUBLIC cxx_std_17)",
    ]
    standalone_host_lines = [
        "# The standalone host declarations are reached only through this",
        "# definition, so a real-host build cannot pull them in silently.",
        f"target_compile_definitions({stem}_common PUBLIC NRPY_DENDRO_STANDALONE_HOST)",
        f"target_include_directories({stem}_common PUBLIC",
        f"  {SRC}/standalone_host",
        ")",
    ]
    if real_host_available:
        lines += [f"if({solver_prefix}_STANDALONE_HOST)"]
        lines += [f"  {line}" if line else line for line in standalone_host_lines]
        lines += [
            "else()",
            "  if(NOT TARGET dendro5 OR NOT TARGET dendro_config OR NOT TARGET toml11::toml11 OR NOT TARGET bssn_common)",
            '    message(FATAL_ERROR "Real GR host requires Dendro-GR targets dendro5, dendro_config, toml11::toml11, and bssn_common")',
            "  endif()",
            f"  target_link_libraries({stem}_common PUBLIC dendro5 toml11::toml11 bssn_common)",
            "  get_target_property(NRPY_DENDRO_INCLUDE_DIRS dendro_config INTERFACE_INCLUDE_DIRECTORIES)",
            "  get_target_property(NRPY_BSSN_INCLUDE_DIRS bssn_common INTERFACE_INCLUDE_DIRECTORIES)",
            "  # Keep -Wall useful for generated sources without diagnosing the",
            "  # fixed external host's headers as though NRPy owned them.",
            f"  target_include_directories({stem}_common SYSTEM PRIVATE",
            "    ${NRPY_DENDRO_INCLUDE_DIRS}",
            "    ${NRPY_BSSN_INCLUDE_DIRS}",
            "  )",
            "endif()",
        ]
    else:
        lines += standalone_host_lines
    lines += [
        f"target_compile_options({stem}_common PRIVATE -Wall)",
        "",
        f"target_include_directories({stem}_common PUBLIC",
        f"  {SRC}/include",
        f"  {SRC}/generated/include",
        ")",
        "",
        f"if(TARGET {exec_or_library_name})",
        f'  message(FATAL_ERROR "{exec_or_library_name} target already defined; use a unique generated solver target")',
        "else()",
        f"  add_executable({exec_or_library_name} src/{stem}_main.cpp)",
        f"  target_compile_options({exec_or_library_name} PRIVATE -Wall)",
        f"  target_link_libraries({exec_or_library_name} PRIVATE {stem}_common MPI::MPI_CXX)",
    ]
    if real_host_available:
        lines += [
            f"  if(NOT {solver_prefix}_STANDALONE_HOST)",
            f"    target_link_libraries({exec_or_library_name} PRIVATE bssn_common)",
            f"    target_include_directories({exec_or_library_name} SYSTEM PRIVATE",
            "      ${NRPY_DENDRO_INCLUDE_DIRS}",
            "      ${NRPY_BSSN_INCLUDE_DIRS}",
            "    )",
            "  endif()",
        ]
    lines += [
        "endif()",
        "",
        f'option({solver_prefix}_ENABLE_CUDA "Build a qualified generated CUDA backend" OFF)',
        f"if({solver_prefix}_ENABLE_CUDA)",
        '  message(FATAL_ERROR "This generated profile is CPU-qualified only")',
        "endif()",
        "",
        "enable_testing()",
    ]
    if real_host_available:
        lines += [f"if({solver_prefix}_STANDALONE_HOST)", "add_subdirectory(tests)", ""]
        lines.extend(standalone_application_ctest)
        lines.append("else()")
        lines.extend(real_application_ctest)
        lines.extend(("endif()", ""))
    else:
        lines += ["add_subdirectory(tests)", ""]
        lines.extend(standalone_application_ctest)
        lines.append("")
    return "\n".join(lines)


def output_tests_cmake(
    solver_prefix: str, solver_stem: str, test_sections: Sequence[str]
) -> str:
    """
    Emit the generated solver's ``tests/CMakeLists.txt``.

    The generated source list is included from the emitted CMake file, so this
    target never carries a source inventory.  Paths are relative to the tests
    directory, so the target is correct both in the standalone-host
    build and when the solver is a subdirectory of a larger Dendro-GR project.

    :param solver_prefix: Bare formulation prefix for CMake variables, e.g. ``BSSN``.
    :param solver_stem: Lowercase stem the emitters use for solver file names.
    :param test_sections: Explicit application-owned self-test section names.
    :return: The CMake file text.

    """
    stem = solver_stem
    lines: List[str] = list(BANNER) + [
        f"set({solver_prefix}_MODULE_ROOT {SRC}/..)",
        f"include({SRC}/../generated/cmake/generated_sources.cmake)",
        f"add_executable({stem}_self_tests",
        f"  {SRC}/{stem}_self_tests.cpp",
        "  ${" + solver_prefix + "_NRPY_GENERATED_SOURCES}",
        ")",
        f"target_include_directories({stem}_self_tests PRIVATE",
        f"  {SRC}/../standalone_host",
        f"  {SRC}/../include",
        f"  {SRC}/../generated/include",
        f"  {SRC}",
        ")",
        f"target_compile_features({stem}_self_tests PRIVATE cxx_std_17)",
        "# This target compiles the generated sources directly rather than",
        "# linking the solver object library, so it carries the same",
        "# standalone-host definition.",
        f"target_compile_definitions({stem}_self_tests PRIVATE NRPY_DENDRO_STANDALONE_HOST)",
        f"target_compile_options({stem}_self_tests PRIVATE -Wall)",
        "",
    ]
    for section in test_sections:
        lines.append(
            f"add_test(NAME {stem}_{section} COMMAND {stem}_self_tests {section})"
        )
        lines.append(f"set_tests_properties({stem}_{section} PROPERTIES TIMEOUT 120)")
    lines.append("")
    return "\n".join(lines)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
