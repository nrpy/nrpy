"""
Generate an NRPy-authored BSSN solver for Dendro-GR.

This is the second formulation lowered through the Dendro infrastructure. It
exists as much to test the infrastructure as to produce a solver: adding it
required the formulation-agnostic lowering to be extracted into
``block_kernel_helpers``, but no existing emitter changed behavior and the fCCZ4
output is unaffected.

The module-level names identify NRPy as the source: the solver directory and
CMake project are ``nrpy_bssn``, and the C++ namespace is ``nrpy::bssn``.
Names inside that module follow the formulation: CMake variables carry the
``BSSN_`` prefix, the production library is ``nrpy_bssn_dendro``, the
qualification driver is ``nrpy_bssn_dendro_qualify``, and the context source
is ``bssnCtx.cpp``.

Run as a module:

    python -m nrpy.examples.dendro_bssn --project-dir project/dendro_bssn --fd-order 6 --no-ko

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

#########################################################
# Step P1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import argparse
import os
from pathlib import Path
from typing import Dict

import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.conditional_file_updater import ConditionalFileUpdater
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.Dendro import CFunction_roles as roles
from nrpy.infrastructures.Dendro import (
    CodeParameters,
    Dendro_defines_h,
    cmake_helpers,
    constants_h,
    parfile,
    state_h,
    types_h,
)
from nrpy.infrastructures.Dendro.general_relativity import (
    bssn_host_adapter,
    constraints_eval,
    enforce_detgbar_equals_detghat_trAzero,
    initial_data,
    main_cpp,
    rhs_eval,
    self_tests_cpp,
    solver_context,
)

# Code-generation-time parameters:
# NRPy authors the module directory, CMake project, and C++ namespace.  Names
# inside that module retain the conventional BSSN stem and Dendro target names.
# These are arguments to the infrastructure, not registered CodeParameters.
solver_name = "nrpy_bssn"
solver_prefix = "BSSN"
solver_stem = "bssn"
solver_namespace = "nrpy::bssn"
production_target = "nrpy_bssn_dendro"
qualification_target = "nrpy_bssn_dendro_qualify"
profile_name = "bssn_cartesian_vacuum"

CoordSystem = "Cartesian"
LapseEvolutionOption = "OnePlusLog"
ShiftEvolutionOption = "GammaDriving2ndOrder_Covariant__Hatted"
enable_parallel_codegen = True


def parse_args() -> argparse.Namespace:
    """
    Parse the command line.

    :return: The parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description="Generate an NRPy-authored BSSN solver for Dendro-GR"
    )
    parser.add_argument("--project-dir", default=os.path.join("project", "dendro_bssn"))
    parser.add_argument(
        "--fd-order",
        type=int,
        choices=(4, 6, 8),
        default=6,
        help="centered finite-difference order; Dendro padding is 2, 3, or 4",
    )
    # argparse.BooleanOptionalAction needs Python 3.9; the supported floor is
    # 3.7, so the two flags are declared explicitly.
    parser.add_argument("--ko", dest="ko", action="store_true")
    parser.add_argument("--no-ko", dest="ko", action="store_false")
    parser.set_defaults(ko=False)
    parser.add_argument(
        "--dendro-gr-host",
        action="store_true",
        help="emit the adapter and CMake target for the Dendro-GR BSSN application",
    )
    # No --parallelization flag: the qualified CPU profile is serial point
    # loops (the kernel runs inside Dendro's own block traversal), and the
    # builders assert that.
    return parser.parse_args()


def main() -> None:
    """Generate the complete Dendro BSSN project."""
    #########################################################
    # Step 1: Parse arguments and set NRPy code-generation parameters.
    args = parse_args()

    par.set_parval_from_str("Infrastructure", "Dendro")
    par.set_parval_from_str("fp_type", "double")
    # The qualified CPU profile emits serial point loops.  The NRPy default is
    # "openmp", so pin it here; the builders assert it.
    par.set_parval_from_str("parallelization", "none")
    par.set_parval_from_str("fd_order", args.fd_order)
    par.set_parval_from_str(
        "EvolvedConformalFactor_cf", "chi" if args.dendro_gr_host else "W"
    )
    par.set_parval_from_str("detgbarOverdetghat_equals_one", True)
    par.set_parval_from_str("enable_parallel_codegen", enable_parallel_codegen)

    #########################################################
    # Step 2: Register independent right-hand-side and constraint C functions.
    rhs_eval.register_CFunctions_rhs_eval(
        solver_stem=solver_stem,
        fd_order=args.fd_order,
        enable_KreissOliger_dissipation=args.ko,
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )
    constraints_eval.register_CFunctions_constraints_eval(
        solver_stem=solver_stem, CoordSystem=CoordSystem
    )

    #########################################################
    # Step 3: Generate functions in parallel.
    #         This Python multiprocessing does not change the serial point
    #         loops generated for Dendro block traversal.
    if enable_parallel_codegen:
        pcg.do_parallel_codegen()

    rhs_by_symbol_name = rhs_eval.rhs_expressions(
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
        enable_KreissOliger_dissipation=args.ko,
    )
    diagnostics_by_name = constraints_eval.diagnostic_expressions(
        CoordSystem=CoordSystem
    )
    ko_fd_order = rhs_eval.DENDRO_FD_PROFILES[args.fd_order][0]

    #########################################################
    # Step 4: Register C functions that consume the merged EVOL registry.
    # Smooth perturbation also registers its amplitude and wavelength
    # CodeParameters. Register the parameter C functions after all other
    # registrations finish.
    initial_data.register_CFunctions_minkowski_initial_data(solver_stem=solver_stem)
    initial_data.register_CFunctions_smooth_perturbation(solver_stem=solver_stem)

    # The smooth ADM conversion, the separate connection-initialization pass,
    # and the det(gammabar)/tr(Abar) enforcement.  Both already read the registered BSSN
    # quantities, so they are shared with the fCCZ4 profile unchanged.
    initial_data.register_CFunctions_ADM_to_BSSN(
        solver_stem=solver_stem, CoordSystem=CoordSystem
    )
    enforce_detgbar_equals_detghat_trAzero.register_CFunctions_enforce_detgbar_equals_detghat_trAzero(
        solver_stem=solver_stem,
        solver_namespace=solver_namespace,
        CoordSystem=CoordSystem,
    )

    # The parameter C functions come last, after every CodeParameter the
    # scientific kernels register is in the registry.
    CodeParameters.register_CFunctions_parameters(solver_stem, solver_namespace)

    #########################################################
    # Step 5: Generate Dendro header and source files, parameter files, tests,
    #         and CMakeLists.txt.
    layout = cmake_helpers.module_layout(solver_name)
    required_padding = roles.required_padding()
    artifacts: Dict[str, str] = {
        layout.generated_include
        + f"{solver_stem}_types.h": types_h.output_types_h(
            solver_stem,
            solver_namespace,
            enforce_detgbar_equals_detghat_trAzero.status_struct_declaration(),
        ),
        layout.generated_include
        + f"{solver_stem}_constants.h": constants_h.output_constants_h(
            solver_stem,
            solver_namespace,
            args.fd_order,
            ko_fd_order,
            ko_fd_order + 2,
            required_padding,
            args.ko,
        ),
        layout.generated_include
        + f"{solver_stem}_state.h": state_h.output_state_h(
            solver_stem, solver_namespace
        ),
        layout.generated_include
        + f"{solver_stem}_parameters.h": CodeParameters.output_parameters_h(
            solver_stem, solver_namespace
        ),
        layout.generated_include
        + f"{solver_stem}_defines.h": Dendro_defines_h.output_Dendro_defines_h(
            solver_stem
        ),
        layout.include
        + f"{solver_stem}Ctx.h": solver_context.output_solver_context_h(
            solver_stem, solver_namespace
        ),
        layout.src
        + f"{solver_stem}Ctx.cpp": solver_context.output_solver_context_cpp(
            solver_stem, solver_namespace
        ),
        layout.src
        + f"{solver_stem}_main.cpp": main_cpp.output_main_cpp(
            solver_stem,
            solver_namespace,
            qualification_target,
            profile_name,
        ),
        layout.pars
        + f"{solver_stem}_minkowski.par": parfile.generate_default_parfile(
            solver_stem,
            profile_name,
            args.fd_order,
            ko_fd_order,
            ko_fd_order + 2,
            required_padding,
            args.ko,
        ),
    }
    artifacts.update(
        {
            layout.root + relative_path: text
            for relative_path, text in self_tests_cpp.output_self_test_artifacts(
                solver_stem,
                solver_namespace,
                rhs_by_symbol_name,
                diagnostics_by_name,
                fd_order=args.fd_order,
                ko_fd_order=ko_fd_order,
                enable_ko=args.ko,
            ).items()
        }
    )
    if args.dendro_gr_host:
        artifacts.update(
            {
                layout.root + relative_path: text
                for relative_path, text in bssn_host_adapter.output_bssn_host_files(
                    args.ko
                ).items()
            }
        )
    artifacts.update(
        cmake_helpers.output_CFunctions_function_prototypes_and_construct_CMakeLists(
            solver_name,
            solver_stem,
            solver_prefix,
            production_target,
            qualification_target,
            self_tests_cpp.test_sections(),
            main_cpp.standalone_ctest_statements(solver_stem, qualification_target),
            main_cpp.real_ctest_statements(solver_stem, qualification_target),
            (("generated/cmake/dendro_gr_host.cmake",) if args.dendro_gr_host else ()),
        )
    )
    for relative_path, text in sorted(artifacts.items()):
        target = Path(args.project_dir) / relative_path
        target.parent.mkdir(parents=True, exist_ok=True)
        # Generated C and C++ source files are clang-formatted, as every other NRPy
        # infrastructure formats what it emits; CMake and TOML are left as
        # written.
        with ConditionalFileUpdater(
            target,
            encoding="utf-8",
            do_format=target.suffix in (".h", ".hpp", ".c", ".cpp", ".cu"),
        ) as file:
            file.write(text)
    solver_dir = Path(args.project_dir) / layout.root
    if not args.dendro_gr_host:
        resolved_solver_dir = solver_dir.resolve()
        for relative_path in (
            "generated/cmake/dendro_gr_host.cmake",
            "src/bssn_dendro_gr_adapter.cpp",
        ):
            target = resolved_solver_dir / relative_path
            if target.exists() or target.is_symlink():
                target.unlink()

    # Copy nrpy.infrastructures.Dendro.standalone_host/dendro_standalone_host.h
    # to <project-dir>/Dendro-GR/<solver_name>/standalone_host/dendro_standalone_host.h.
    copy_files(
        package="nrpy.infrastructures.Dendro.standalone_host",
        filenames_list=["dendro_standalone_host.h"],
        project_dir=str(solver_dir),
        subdirectory="standalone_host",
    )

    copy_files(
        package="nrpy.infrastructures.Dendro",
        filenames_list=["block_geometry.h"],
        project_dir=str(solver_dir),
        subdirectory="include",
    )

    EVOL, _AUXEVOL, _DIAG, _AUX = gri.GridFunction.gridfunction_lists()
    print(f"Finished generating {solver_name} in {args.project_dir}.")
    print(f"  profile: {profile_name}")
    print(f"  evolved variables: {len(EVOL)}")
    print(f"  finite-difference order: {args.fd_order}")
    print(f"  KO finite-difference order: {ko_fd_order}")
    print(f"  effective KO difference order: {ko_fd_order + 2}")
    print(f"  Kreiss-Oliger dissipation: {'enabled' if args.ko else 'disabled'}")
    print(f"  required ghost points: {roles.required_padding()}")
    if args.dendro_gr_host:
        print("  Dendro-GR host target: nrpy_bssnSolver")
    print("Now build and run the generated self-tests with:")
    build_dir = solver_dir / "build"
    print(f"  cmake -S {solver_dir} -B {build_dir}")
    print(f"  cmake --build {build_dir} && ctest --test-dir {build_dir}")


if __name__ == "__main__":
    main()
