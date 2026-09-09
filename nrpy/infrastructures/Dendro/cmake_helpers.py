# nrpy/infrastructures/Dendro/cmake_helpers.py
"""
CMake emitters for a generated Dendro solver.

Every emitted name derives from the ``solver_name``, ``solver_prefix`` and
``exec_or_library_name`` arguments threaded from the calling example, matching
how BHaH threads ``project_name`` and ETLegacy threads ``thorn_name``.

Dendro's own conventions govern the emitted identifiers: the solver directory
carries the ``_GR`` suffix (``BSSN_GR``) while CMake variables and library
targets carry the bare formulation prefix (``BSSN_ENABLE_CUDA``,
``bssn_common``), and solver sources are named for the formulation
(``bssnCtx.cpp``), never for the framework.

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

# The Dendrolib commit the generated kernels' block-layout assumptions were
# proven against, by nrpy/infrastructures/Dendro/tests_infra.  Emitted into the
# generated project so a consumer building against a different commit is warned.
PROVEN_DENDROLIB_COMMIT = "246043709e806021fcfc011fe657b8bf964cae4c"


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

    :param solver_name: Solver directory name, e.g. Dendro-GR's own ``BSSN_GR``.
    :return: The module layout.

    Doctests:
    >>> layout = module_layout("FCCZ4_GR")
    >>> layout.root
    'Dendro-GR/FCCZ4_GR/'
    >>> layout.generated_include
    'Dendro-GR/FCCZ4_GR/generated/include/'
    >>> layout.tests
    'Dendro-GR/FCCZ4_GR/tests/'
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
    lines.append("struct BlockGeometry;")
    lines.append("struct StandaloneHostMesh;")
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
) -> Dict[str, str]:
    """
    Emit the CFunction sources, the prototypes header and the CMake files.

    The name is the whole contract, as it is for
    ``BHaH/Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile``:
    one source per registered CFunction, the declaration header, and the build
    files a consumer needs to compile them.

    :param solver_name: Solver directory name, e.g. Dendro's own ``BSSN_GR``.
    :param solver_stem: Lowercase formulation stem for emitted header names.
    :param solver_prefix: Bare formulation prefix for the CMake variables.
    :param exec_or_library_name: Name of the solver executable target.
    :param test_sections: Explicit application-owned self-test section names.
    :param standalone_application_ctest: Standalone application CTest lines.
    :param real_application_ctest: Real-host application CTest lines.
    :return: Mapping of ``Dendro-GR/<solver_name>/<path>`` to file text.
    """
    layout = module_layout(solver_name)
    prefix = layout.root
    artifacts: Dict[str, str] = {}
    # One source file per registered CFunction, written verbatim: the umbrella
    # include of ``<stem>_defines.h`` rides on the CFunction's own ``includes``
    # field, which supplies the host types, the generated state and parameter
    # declarations, and the standard headers the body uses -- the role
    # ``BHaH_defines.h`` plays for every BHaH CFunction.
    for name, cfunc in registered_CFunctions():
        artifacts[prefix + derived_source_path(name, cfunc.subdirectory)] = (
            cfunc.full_function
        )
    artifacts[prefix + f"generated/include/{solver_stem}_function_prototypes.h"] = (
        output_function_prototypes_h(solver_stem)
    )
    # The three CMake files this module also emits, so the name is the whole
    # contract: per-CFunction sources, the prototypes header, and the build
    # files, exactly as BHaH's Makefile_helpers does for a Makefile.
    artifacts[layout.generated_cmake + "nrpy_generated_sources.cmake"] = (
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
    Emit ``generated/cmake/nrpy_generated_sources.cmake``.

    One entry per registered CFunction, in sorted name order.

    :param solver_prefix: Bare formulation prefix for CMake variables, e.g. ``BSSN`` for Dendro's own ``BSSN_GR``.
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
) -> str:
    """
    Emit the generated solver's ``CMakeLists.txt``.

    :param solver_name: Solver directory and CMake project name, e.g. Dendro's own ``BSSN_GR``.
    :param solver_prefix: Bare formulation prefix for CMake variables, e.g. ``BSSN`` for Dendro's own ``BSSN_GR``.
    :param solver_stem: Lowercase stem the emitters use for solver file names.
    :param exec_or_library_name: Name of the solver executable target.
    :param standalone_application_ctest: Standalone application CTest lines.
    :param real_application_ctest: Real-host application CTest lines.
    :return: The CMake file text.

    Doctests:
    >>> text = output_solver_cmake("WAVE", "WAVE", "wave", "waveSolver", ("add_test(NAME wave_run COMMAND waveSolver)",), ())
    >>> "project(WAVE CXX)" in text
    True
    >>> "add_executable(waveSolver" in text
    True
    >>> "add_test(NAME wave_run COMMAND waveSolver)" in text
    True
    >>> "${WAVE_NRPY_GENERATED_SOURCES}" in text
    True
    >>> "src/bssnCtx.cpp" in output_solver_cmake("Z_GR", "ZORP", "zrp", "zSolver", (), ())
    False
    >>> "src/zrpCtx.cpp" in output_solver_cmake("Z_GR", "ZORP", "zrp", "zSolver", (), ())
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
        "# Dendrolib source pin.  The block geometry, layout, padding, offset,",
        "# origin and halo axes these kernels assume were proven against this",
        "# commit.  Dendro-GR declares DENDRO_dendrolib_GIT_TAG as a cache",
        "# variable before adding this subdirectory, so setting it here would be",
        "# discarded; the proven commit is recorded in a solver-scoped variable",
        "# and compared instead, which is the only form that can warn.",
        f'set({solver_prefix}_PROVEN_DENDROLIB_COMMIT "{PROVEN_DENDROLIB_COMMIT}")',
        f"if(DEFINED DENDRO_dendrolib_GIT_TAG AND NOT DENDRO_dendrolib_GIT_TAG STREQUAL {solver_prefix}_PROVEN_DENDROLIB_COMMIT)",
        "  message(WARNING",
        f'          "These kernels were proven against Dendrolib commit ${{{solver_prefix}_PROVEN_DENDROLIB_COMMIT}}, "',
        '          "but the build is configured for ${DENDRO_dendrolib_GIT_TAG}.  Re-run "',
        '          "nrpy/infrastructures/Dendro/tests_infra before relying on the block-layout claims.")',
        "endif()",
        "",
        "# The solver lifecycle runs under MPI for the rank-agreement gates.",
        "# MPI is a standard host dependency, not Dendrolib; the real Dendro-GR",
        "# build also provides it.",
        "find_package(MPI REQUIRED COMPONENTS CXX)",
        "",
        f'option({solver_prefix}_STANDALONE_HOST "Build the solver standalone, against the NRPy-supplied Dendro host declarations rather than a real Dendro-GR build" ON)',
        f"set({solver_prefix}_MODULE_ROOT {SRC})",
        f"include({SRC}/generated/cmake/nrpy_generated_sources.cmake)",
        "",
        f"add_library({stem}_common OBJECT",
        f"  src/{stem}Ctx.cpp",
        "  ${" + solver_prefix + "_NRPY_GENERATED_SOURCES}",
        ")",
        "",
        f"target_compile_features({stem}_common PUBLIC cxx_std_17)",
        "# The standalone host declarations are reached only through this",
        "# definition, so a real-host build cannot pull them in silently.",
        f"if({solver_prefix}_STANDALONE_HOST)",
        f"  target_compile_definitions({stem}_common PUBLIC NRPY_DENDRO_STANDALONE_HOST)",
        "else()",
        "  if(NOT TARGET dendro5 OR NOT TARGET toml11::toml11)",
        '    message(FATAL_ERROR "Real host requires Dendro-GR targets dendro5 and toml11::toml11")',
        "  endif()",
        f"  target_link_libraries({stem}_common PUBLIC dendro5 toml11::toml11)",
        "endif()",
        f"target_compile_options({stem}_common PRIVATE -Wall)",
        "",
        "# The standalone host header is PRIVATE: it must not leak onto the include",
        "# path of every consumer of the solver.",
        f"target_include_directories({stem}_common PUBLIC",
        f"  {SRC}/include",
        f"  {SRC}/generated/include",
        ")",
        f"target_include_directories({stem}_common PRIVATE {SRC}/standalone_host)",
        "",
        f"if(TARGET {exec_or_library_name})",
        f'  message(FATAL_ERROR "{exec_or_library_name} target already defined; use a unique generated solver target")',
        "else()",
        f"  add_executable({exec_or_library_name} src/{stem}_main.cpp)",
        f"  target_include_directories({exec_or_library_name} PRIVATE {SRC}/standalone_host)",
        f"  target_compile_options({exec_or_library_name} PRIVATE -Wall)",
        f"  target_link_libraries({exec_or_library_name} PRIVATE {stem}_common MPI::MPI_CXX)",
        "endif()",
        "",
        f'option({solver_prefix}_ENABLE_CUDA "Build a qualified generated CUDA backend" OFF)',
        f"if({solver_prefix}_ENABLE_CUDA)",
        '  message(FATAL_ERROR "This generated profile is CPU-qualified only")',
        "endif()",
        "",
        "enable_testing()",
        f"if({solver_prefix}_STANDALONE_HOST)",
        "add_subdirectory(tests)",
        "",
    ]
    lines.extend(standalone_application_ctest)
    lines.append("else()")
    lines.extend(real_application_ctest)
    lines.extend(("endif()", ""))
    return "\n".join(lines)


def output_tests_cmake(
    solver_prefix: str, solver_stem: str, test_sections: Sequence[str]
) -> str:
    """
    Emit the generated solver's ``tests/CMakeLists.txt``.

    The generated source list is included from the emitted CMake file, so this
    target never carries a source inventory.  Paths are relative to the tests
    directory, so the target is correct both in the standalone standalone-host
    build and when the solver is a subdirectory of a larger Dendro-GR project.

    :param solver_prefix: Bare formulation prefix for CMake variables, e.g. ``BSSN`` for Dendro's own ``BSSN_GR``.
    :param solver_stem: Lowercase stem the emitters use for solver file names.
    :param test_sections: Explicit application-owned self-test section names.
    :return: The CMake file text.

    Doctests:
    >>> text = output_tests_cmake("WAVE", "wave", ("state", "rhs"))
    >>> "include(${CMAKE_CURRENT_SOURCE_DIR}/../generated/cmake/nrpy_generated_sources.cmake)" in text
    True
    >>> "set(WAVE_MODULE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)" in text
    True
    >>> "${WAVE_NRPY_GENERATED_SOURCES}" in text
    True
    >>> "add_test(NAME wave_state" in text and "add_test(NAME wave_rhs" in text
    True
    """
    stem = solver_stem
    lines: List[str] = list(BANNER) + [
        f"set({solver_prefix}_MODULE_ROOT {SRC}/..)",
        f"include({SRC}/../generated/cmake/nrpy_generated_sources.cmake)",
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
