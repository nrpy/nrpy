# nrpy/infrastructures/Dendro/cmake.py
"""
CMake renderers for the generated Dendro project.

The module CMake (section 12.4) contains no concrete field list, no physics
parameter, no FD coefficient, and no numerical loop; the generated source
list is derived from the frozen CFunction set (section 9.13/12.1).  The
root-CMake integration marker block (section 12.3) is owned solely by the
installer (copy_tool.marker_lines), the single source that applies it.
"""

from nrpy.infrastructures.Dendro import CFunction_output
from nrpy.infrastructures.Dendro.freeze import FrozenNRPyDendroSnapshot


def render_nrpy_generated_sources_cmake(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render generated/cmake/nrpy_generated_sources.cmake.

    The source list is derived from the frozen CFunctions: one entry per
    registered CFunction, in frozen (sorted name) order (section 12.1).

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: The CMake file text.
    """
    lines = [
        "# GENERATED FILE - DO NOT EDIT",
        "# Generator: nrpy.infrastructures.Dendro",
        f"# Module ABI: {snapshot.hashes.module_abi_hash}",
        "# Source object: frozen CFunction set (one source per registered CFunction)",
        "# Source list derived from the frozen CFunction set (section 12.1).",
        "# FCCZ4_MODULE_ROOT is set by the including scope so this list is",
        "# correct whether included from the module or the tests directory.",
        "set(FCCZ4_NRPY_GENERATED_SOURCES",
    ]
    for path in CFunction_output.CFunction_cmake_source_list(snapshot):
        lines.append(f"  ${{FCCZ4_MODULE_ROOT}}/{path}")
    lines.append(")")
    lines.append("")
    return "\n".join(lines)


def render_module_cmake(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the FCCZ4_GR module CMakeLists.txt (section 12.4).

    The target declaration links the runtime context source plus every
    generated source; no fCCZ4 order switch or field list appears here.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: The CMake file text.
    """
    lines = [
        "# GENERATED FILE - DO NOT EDIT",
        "# Generator: nrpy.infrastructures.Dendro",
        f"# Module ABI: {snapshot.hashes.module_abi_hash}",
        "cmake_minimum_required(VERSION 3.13)",
        "project(FCCZ4_GR CXX)",
        "",
        "# The fccz4Solver lifecycle (section 16.8) runs under MPI for the",
        "# -n 1/-n 2 rank-agreement gates; MPI is a standard host dependency",
        "# (not Dendrolib).  The real Dendro-GR build also provides MPI.",
        "find_package(MPI REQUIRED COMPONENTS CXX)",
        "",
        "option(FCCZ4_USE_HOST_MOCK",
        '       "Build against the NRPy mock Dendro host header (required in this profile)" ON)',
        "if(NOT FCCZ4_USE_HOST_MOCK)",
        "  message(FATAL_ERROR",
        '    "This generated profile is mock-vehicle only: the real ot::Block/DVector"',
        '    " signatures freeze after the pinned-Dendrolib capability gates (section 13.2).")',
        "endif()",
        "",
        "set(FCCZ4_MODULE_ROOT ${CMAKE_CURRENT_SOURCE_DIR})",
        "include(${CMAKE_CURRENT_SOURCE_DIR}/generated/cmake/nrpy_generated_sources.cmake)",
        "",
        "add_library(fccz4_common OBJECT",
        "  src/fccz4Ctx.cpp",
        "  ${FCCZ4_NRPY_GENERATED_SOURCES}",
        ")",
        "",
        "target_compile_features(fccz4_common PUBLIC cxx_std_17)",
        "target_compile_options(fccz4_common PRIVATE -Wall -Wextra -Werror)",
        "",
        "# The mock host header is an implementation detail of this profile, so it is",
        "# PRIVATE: it must not leak onto the include path of every consumer of the",
        "# module (section 13.2 keeps the real host types unfrozen).",
        "target_include_directories(fccz4_common PUBLIC",
        "  ${CMAKE_CURRENT_SOURCE_DIR}/include",
        "  ${CMAKE_CURRENT_SOURCE_DIR}/generated/include",
        ")",
        "target_include_directories(fccz4_common PRIVATE",
        "  ${CMAKE_CURRENT_SOURCE_DIR}/host_mock",
        ")",
        "",
        "# This profile is qualified against the NRPy mock host only.  Linking the",
        "# real Dendro-GR libraries would compile the mock ot::Block/DVector doubles",
        "# against the real host: refuse rather than produce a build whose types are",
        "# not the pinned ones (sections 13.2/13.3).",
        "if(TARGET dendro5)",
        "  message(FATAL_ERROR",
        '    "A real dendro5 target is present, but this generated profile is"',
        '    " mock-vehicle only.  Regenerate after the pinned-Dendrolib capability"',
        '    " gates (section 13.2) before building against the real host.")',
        "endif()",
        "",
        "if(TARGET fccz4Solver)",
        '  message(STATUS "fccz4Solver target already defined by the host")',
        "else()",
        "  add_executable(fccz4Solver src/fccz4_main.cpp)",
        "  target_include_directories(fccz4Solver PRIVATE",
        "    ${CMAKE_CURRENT_SOURCE_DIR}/host_mock)",
        "  target_compile_options(fccz4Solver PRIVATE -Wall -Wextra -Werror)",
        "  target_link_libraries(fccz4Solver PRIVATE fccz4_common MPI::MPI_CXX)",
        "endif()",
        "",
        'option(FCCZ4_ENABLE_CUDA "Build a qualified generated CUDA backend" OFF)',
        "if(FCCZ4_ENABLE_CUDA)",
        '  message(FATAL_ERROR "This generated profile is CPU-qualified only")',
        "endif()",
        "",
        "# add_test() only registers with CTest when testing is enabled in scope.",
        "enable_testing()",
        "add_subdirectory(tests)",
        "",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
