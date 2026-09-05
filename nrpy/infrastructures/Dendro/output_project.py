# nrpy/infrastructures/Dendro/output_project.py
"""
Write the complete generated Dendro solver project to disk.

Each artifact comes from the emitter that owns it, and every emitter reads the
NRPy registries directly, as BHaH's emitters read ``par.glb_code_params_dict``.
This module owns no formulation choice and holds no state of its own: it maps
emitter output onto project-relative paths and writes it.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import os
from pathlib import Path
from typing import Dict, Union

import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.Dendro import (
    CodeParameters,
    Dendro_main_cpp,
    Dendro_parameter_file,
    Dendro_preamble_h,
    Dendro_README_md,
    Dendro_self_tests_cpp,
    Dendro_solver_context,
    Dendro_state_h,
    Dendro_types_h,
    cmake_helpers,
    output_CFunctions,
)
from nrpy.infrastructures.Dendro import registration as reg


def output_project(
    *,
    project_dir: Union[str, "os.PathLike[str]"],
    solver_name: str,
    solver_prefix: str,
    solver_stem: str,
    solver_namespace: str,
    exec_or_library_name: str,
    profile_name: str,
    generator_module: str,
) -> str:
    """
    Write the complete generated project to ``project_dir``.

    :param project_dir: Target project directory.
    :param solver_name: Solver directory and CMake project name, e.g. Dendro's own ``BSSN_GR``.
    :param solver_prefix: Bare formulation prefix for CMake variables, e.g. ``BSSN``.
    :param solver_stem: Lowercase formulation stem for emitted file names, e.g. ``bssn``.
    :param solver_namespace: Solver namespace, following Dendro's lowercase formulation habit.
    :param exec_or_library_name: Name of the solver executable target.
    :param profile_name: Name of the generation profile.
    :param generator_module: Dotted module name of the generating example.
    :return: The project directory path.
    """
    project_path = Path(project_dir)
    scalar_type = str(par.parval_from_str("Dendro_scalar_type"))
    required_padding = reg.required_padding()
    root = f"Dendro-GR/{solver_name}/"
    generated_include = root + "generated/include/"

    artifacts: Dict[str, str] = {
        generated_include
        + f"{solver_stem}_types.h": Dendro_types_h.output_Dendro_types_h(
            solver_namespace
        ),
        generated_include
        + f"{solver_stem}_constants.h": Dendro_state_h.output_Dendro_constants_h(
            solver_namespace, required_padding
        ),
        generated_include
        + f"{solver_stem}_state.h": Dendro_state_h.output_Dendro_state_h(
            solver_stem, solver_namespace
        ),
        generated_include
        + f"{solver_stem}_parameters.h": CodeParameters.output_Dendro_parameters_h(
            solver_stem, solver_namespace
        ),
        generated_include
        + "generated_project_preamble.h": Dendro_preamble_h.output_Dendro_preamble_h(
            solver_stem
        ),
        root
        + "generated/cmake/nrpy_generated_sources.cmake": (
            cmake_helpers.output_generated_sources_cmake(solver_prefix)
        ),
        root
        + f"include/{solver_stem}Ctx.h": (
            Dendro_solver_context.output_Dendro_solver_context_h(
                solver_stem, solver_namespace, scalar_type
            )
        ),
        root
        + f"src/{solver_stem}Ctx.cpp": (
            Dendro_solver_context.output_Dendro_solver_context_cpp(
                solver_stem, solver_namespace, scalar_type
            )
        ),
        root
        + f"src/{solver_stem}_main.cpp": Dendro_main_cpp.output_Dendro_main_cpp(
            solver_stem, solver_namespace, scalar_type, exec_or_library_name
        ),
        root
        + f"pars/{solver_stem}_minkowski.par": (
            Dendro_parameter_file.output_Dendro_parameter_file(
                solver_stem, profile_name, required_padding
            )
        ),
        root
        + "tests/test_generated.cpp": (
            Dendro_self_tests_cpp.output_Dendro_self_tests_cpp(
                solver_stem, solver_namespace, scalar_type
            )
        ),
        root
        + "tests/CMakeLists.txt": cmake_helpers.output_tests_cmake(
            solver_prefix, solver_stem
        ),
        root
        + "CMakeLists.txt": cmake_helpers.output_solver_cmake(
            solver_name, solver_prefix, solver_stem, exec_or_library_name
        ),
        root
        + "README.md": Dendro_README_md.output_solver_README_md(
            solver_name, generator_module
        ),
    }
    artifacts.update(output_CFunctions.CFunction_artifacts(solver_name, solver_stem))
    artifacts["README.md"] = Dendro_README_md.output_project_README_md(
        project_path.name, solver_name, generator_module
    )
    for relative_path, text in sorted(artifacts.items()):
        target = project_path / relative_path
        target.parent.mkdir(parents=True, exist_ok=True)
        with open(target, "w", encoding="utf-8", newline="\n") as file:
            file.write(text)
    # The mock host header is a fixed source asset of this package, copied
    # verbatim exactly as BHaH copies simd_intrinsics.h.  The real Dendrolib
    # header arrives once Dendrolib is pinned.
    copy_files(
        package="nrpy.infrastructures.Dendro.host_mock",
        filenames_list=["dendro_mock.h"],
        project_dir=str(project_path / "Dendro-GR" / solver_name),
        subdirectory="host_mock",
    )
    return str(project_path)
