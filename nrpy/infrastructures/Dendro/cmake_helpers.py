# nrpy/infrastructures/Dendro/cmake_helpers.py
"""
CMake emitters for a generated Dendro solver.

Every emitted name derives from the solver and target names threaded from the
calling example, matching
how BHaH threads ``project_name`` and ETLegacy threads ``thorn_name``.

NRPy provenance appears in each generated module directory and CMake project.
Names inside the module retain the formulation stem and required Dendro target
names, such as ``BSSN_ENABLE_CUDA``, ``nrpy_bssn_dendro``, and ``bssnCtx.cpp``.

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
SRC = "${CMAKE_CURRENT_LIST_DIR}"


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
    production_target: str,
    qualification_target: str,
    test_sections: Sequence[str],
    standalone_application_ctest: Sequence[str],
    real_application_ctest: Sequence[str],
    module_cmake_includes: Sequence[str] = (),
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
    :param production_target: Name of the Dendro production library target.
    :param qualification_target: Name of the generated qualification driver.
    :param test_sections: Explicit application-owned self-test section names.
    :param standalone_application_ctest: Standalone application CTest lines.
    :param real_application_ctest: Real-host application CTest lines.
    :param module_cmake_includes: Solver-relative CMake files to include after
        the production library and qualification-driver definitions.
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
        production_target,
        qualification_target,
        standalone_application_ctest,
        real_application_ctest,
        module_cmake_includes,
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
    production_target: str,
    qualification_target: str,
    standalone_application_ctest: Sequence[str],
    real_application_ctest: Sequence[str],
    module_cmake_includes: Sequence[str] = (),
) -> str:
    """
    Emit the generated solver's ``CMakeLists.txt``.

    :param solver_name: NRPy module directory and CMake project name, e.g. ``nrpy_bssn``.
    :param solver_prefix: Bare formulation prefix for CMake variables, e.g. ``BSSN``.
    :param solver_stem: Lowercase stem the emitters use for solver file names.
    :param production_target: Name of the Dendro production library target.
    :param qualification_target: Name of the generated qualification driver.
    :param standalone_application_ctest: Standalone application CTest lines.
    :param real_application_ctest: Real-host application CTest lines.
    :param module_cmake_includes: Solver-relative CMake files to include after
        the standard generated targets.
    :return: The CMake file text.

    >>> text = output_solver_cmake(
    ...     "nrpy_bssn", "BSSN", "bssn", "nrpy_bssn_dendro",
    ...     "nrpy_bssn_dendro_qualify", ("add_test(NAME bssn_standalone)",),
    ...     ("add_test(NAME bssn_real_minkowski)",),
    ...     ("generated/cmake/dendro_gr_host.cmake",),
    ... )
    >>> all(name in text for name in ("nrpy_bssn_dendro", "nrpy_bssn_dendro_qualify"))
    True
    >>> "bssn_common" in text
    False
    >>> obsolete = ("BSSN_STANDALONE_HOST", "bssnSolver", "_nrpy_dendro_legacy_driver_request")
    >>> any(name in text for name in obsolete)
    False
    >>> 'include("${CMAKE_CURRENT_LIST_DIR}/generated/cmake/dendro_gr_host.cmake")' in text
    True

    """
    stem = solver_stem
    lines: List[str] = list(BANNER) + [
        "cmake_minimum_required(VERSION 3.13)",
        "set(_nrpy_dendro_top_level OFF)",
        "if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)",
        f"  project({solver_name} CXX)",
        "  set(_nrpy_dendro_top_level ON)",
        "endif()",
        "",
        "set(_nrpy_dendro_option_default OFF)",
        "if(_nrpy_dendro_top_level)",
        "  set(_nrpy_dendro_option_default ON)",
        "endif()",
        "option(NRPY_DENDRO_BUILD_DRIVERS",
        '       "Build generated Dendro qualification drivers"',
        "       ${_nrpy_dendro_option_default})",
        "option(NRPY_DENDRO_BUILD_TESTS",
        '       "Build generated Dendro numerical tests"',
        "       ${_nrpy_dendro_option_default})",
        "",
        f'set({solver_prefix}_MODULE_ROOT "{SRC}")',
        f'include("{SRC}/generated/cmake/generated_sources.cmake")',
        "",
        "if(TARGET dendro5)",
        f"  add_library({production_target}",
        f'    "{SRC}/src/{stem}Ctx.cpp"',
        "  ${" + solver_prefix + "_NRPY_GENERATED_SOURCES}",
        "  )",
        f"  target_compile_features({production_target} PUBLIC cxx_std_17)",
        f"  target_compile_options({production_target} PRIVATE -Wall)",
        f"  target_include_directories({production_target} PUBLIC",
        f'    "{SRC}/include"',
        f'    "{SRC}/generated/include"',
        "  )",
        f"  target_link_libraries({production_target} PUBLIC dendro5)",
        "elseif(NOT _nrpy_dendro_top_level)",
        f'  message(FATAL_ERROR "{production_target} requires the Dendro target dendro5")',
        "endif()",
        "",
        "if(NRPY_DENDRO_BUILD_DRIVERS)",
        "  find_package(MPI REQUIRED COMPONENTS CXX)",
        "  if(_nrpy_dendro_top_level)",
        f"    add_executable({qualification_target}",
        f'      "{SRC}/src/{stem}_main.cpp"',
        f'      "{SRC}/src/{stem}Ctx.cpp"',
        "      ${" + solver_prefix + "_NRPY_GENERATED_SOURCES}",
        "    )",
        f"    target_compile_definitions({qualification_target} PRIVATE NRPY_DENDRO_STANDALONE_HOST)",
        f"    target_include_directories({qualification_target} PRIVATE",
        f'      "{SRC}/standalone_host"',
        f'      "{SRC}/include"',
        f'      "{SRC}/generated/include"',
        "    )",
        f"    target_link_libraries({qualification_target} PRIVATE MPI::MPI_CXX)",
        "  else()",
        "    if(NOT TARGET toml11::toml11)",
        f'      message(FATAL_ERROR "{qualification_target} requires toml11::toml11")',
        "    endif()",
        f'    add_executable({qualification_target} "{SRC}/src/{stem}_main.cpp")',
        f"    target_link_libraries({qualification_target} PRIVATE {production_target} MPI::MPI_CXX toml11::toml11)",
        "  endif()",
        f"  target_compile_features({qualification_target} PRIVATE cxx_std_17)",
        f"  target_compile_options({qualification_target} PRIVATE -Wall)",
        "endif()",
        "",
        f'option({solver_prefix}_ENABLE_CUDA "Build a qualified generated CUDA backend" OFF)',
        f"if({solver_prefix}_ENABLE_CUDA)",
        '  message(FATAL_ERROR "This generated profile is CPU-qualified only")',
        "endif()",
        "",
    ]
    for cmake_include in module_cmake_includes:
        lines.extend(
            (
                f'include("${{CMAKE_CURRENT_LIST_DIR}}/{cmake_include}")',
                "",
            )
        )
    lines += [
        "if(NRPY_DENDRO_BUILD_TESTS)",
        "  if(_nrpy_dendro_top_level)",
        "    enable_testing()",
        "  endif()",
        f'  add_subdirectory("{SRC}/tests" "${{CMAKE_CURRENT_BINARY_DIR}}/{stem}_tests")',
    ]
    lines += ["  if(NRPY_DENDRO_BUILD_DRIVERS)", "    if(_nrpy_dendro_top_level)"]
    lines.extend(f"      {line}" for line in standalone_application_ctest)
    lines.append("    else()")
    lines.extend(f"      {line}" for line in real_application_ctest)
    lines.extend(("    endif()", "  endif()", "endif()", ""))
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
        f'get_filename_component({solver_prefix}_MODULE_ROOT "{SRC}/.." ABSOLUTE)',
        f'include("${{{solver_prefix}_MODULE_ROOT}}/generated/cmake/generated_sources.cmake")',
        f"add_executable({stem}_self_tests",
        f'  "{SRC}/{stem}_self_tests.cpp"',
        "  ${" + solver_prefix + "_NRPY_GENERATED_SOURCES}",
        ")",
        f"target_include_directories({stem}_self_tests PRIVATE",
        f'  "{SRC}/../standalone_host"',
        f'  "{SRC}/../include"',
        f'  "{SRC}/../generated/include"',
        f'  "{SRC}"',
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
