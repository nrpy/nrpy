"""
Generate an NRPy-authored fCCZ4 solver for Dendro-GR.

Run as a module:

    python -m nrpy.examples.dendro_fccz4 \
        --project-dir project/dendro_fccz4 --fd-order 4 --no-ko

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import argparse
import os
from pathlib import Path
from typing import Dict

import nrpy.grid as gri
import nrpy.params as par
from nrpy.helpers.conditional_file_updater import ConditionalFileUpdater
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.Dendro import CFunction_roles as roles
from nrpy.infrastructures.Dendro import (
    CodeParameters,
    Dendro_defines_h,
    cmake_helpers,
    cmdline_input_and_parfiles,
    constants_h,
    main_cpp,
    self_tests_cpp,
    solver_context,
    state_h,
    types_h,
)
from nrpy.infrastructures.Dendro.general_relativity import (
    constraints_eval,
    enforce_detgbar_equals_detghat_trAzero,
    initial_data,
    rhs_eval,
)

# Dendro names a solver directory for its formulation (its own is BSSN_GR) and
# namespaces the solver by the lowercase formulation (namespace bssn), so this
# generated solver does the same.  These are arguments to the infrastructure,
# not registered CodeParameters: NRPy threads project_name, thorn_name and
# exec_or_library_name the same way.
solver_name = "FCCZ4_GR"
solver_prefix = "FCCZ4"
solver_stem = "fccz4"
solver_namespace = "fccz4"
exec_or_library_name = "fccz4Solver"
profile_name = "fccz4_cartesian_vacuum"

CoordSystem = "Cartesian"
LapseEvolutionOption = "OnePlusLog"
ShiftEvolutionOption = "GammaDriving2ndOrder_Covariant__Hatted"


def parse_args() -> argparse.Namespace:
    """
    Parse the command line.

    :return: The parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description="Generate an NRPy-authored fCCZ4 solver for Dendro-GR"
    )
    parser.add_argument(
        "--project-dir", default=os.path.join("project", "dendro_fccz4")
    )
    # fd_order 8 reaches five ghost points.  Padding 5 is proven on the
    # pinned Dendrolib at element order 10, so the limit is this
    # generator's qualified set rather than the host.
    parser.add_argument(
        "--fd-order",
        type=int,
        choices=(2, 4, 6),
        default=4,
        help="finite-difference order; 8 is not in this generator's "
        "qualified set, though the pinned host proves padding 5",
    )
    # argparse.BooleanOptionalAction needs Python 3.9; the supported floor is
    # 3.7, so the two flags are declared explicitly.
    parser.add_argument("--ko", dest="ko", action="store_true")
    parser.add_argument("--no-ko", dest="ko", action="store_false")
    parser.set_defaults(ko=False)
    # No --parallelization flag: the qualified CPU profile is serial point
    # loops (the kernel runs inside Dendro's own block traversal), and the
    # builders assert that.  Offering a flag the generator discards would
    # silently produce an unqualified configuration.
    return parser.parse_args()


def main() -> None:
    """Generate the complete Dendro fCCZ4 project."""
    args = parse_args()

    #########################################################
    # Step 1: Set the generation profile.
    par.set_parval_from_str("Infrastructure", "Dendro")
    par.set_parval_from_str("fp_type", "double")
    # The qualified CPU profile emits serial point loops.  The NRPy default is
    # "openmp", so pin it here; the builders assert it.
    par.set_parval_from_str("parallelization", "none")
    par.set_parval_from_str("fd_order", args.fd_order)
    par.set_parval_from_str("EvolvedConformalFactor_cf", "chi")
    par.set_parval_from_str("detgbarOverdetghat_equals_one", True)

    #########################################################
    # Step 2: Register the generated C functions.  The right-hand side goes
    #         first: it registers the exact gridfunctions and physics
    #         CodeParameters through the shared fCCZ4 expression bundle, and
    #         records the ghost points the emitted operators reach.
    rhs_eval.register_CFunctions_rhs_eval(
        solver_stem=solver_stem,
        enable_fCCZ4=True,
        fd_order=args.fd_order,
        enable_KreissOliger_dissipation=args.ko,
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )

    # Minkowski initial data and the smooth analytic perturbation the
    # lifecycle gates evolve.  Both are NRPy-authored kernels.
    initial_data.register_CFunctions_minkowski_initial_data(solver_stem=solver_stem)
    initial_data.register_CFunctions_smooth_perturbation(solver_stem=solver_stem)

    # The smooth ADM conversion, the separate connection-initialization pass,
    # and the det(gammabar)/tr(Abar) enforcement.  The conversion registers the ADM source
    # fields as AUXEVOL; the enforcement is scheduled after initial data and
    # after every accepted timestep by the host context.
    initial_data.register_CFunctions_ADM_to_BSSN(
        solver_stem=solver_stem, CoordSystem=CoordSystem
    )
    enforce_detgbar_equals_detghat_trAzero.register_CFunctions_enforce_detgbar_equals_detghat_trAzero(
        solver_stem=solver_stem,
        solver_namespace=solver_namespace,
        CoordSystem=CoordSystem,
    )

    # The constraint diagnostics.  H_Z4 and the connection constraint are
    # registered as DIAG gridfunctions: they are recomputed from the evolved
    # state and are never checkpoint state.
    constraints_eval.register_CFunctions_constraints_eval(
        solver_stem=solver_stem,
        enable_fCCZ4=True,
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
        enable_KreissOliger_dissipation=args.ko,
    )

    # The parameter C functions come last, after every CodeParameter the
    # scientific kernels register is in the registry.
    CodeParameters.register_CFunctions_parameters(solver_stem, solver_namespace)

    #########################################################
    # Step 3: Assemble and write the project.  The assembly lives here, in the
    # example, exactly as it does for BHaH, ETLegacy, CarpetX and superB: an
    # example reads top to bottom as the complete recipe for one solver.
    layout = cmake_helpers.module_layout(solver_name)
    required_padding = roles.required_padding()
    artifacts: Dict[str, str] = {
        layout.generated_include
        + f"{solver_stem}_types.h": types_h.output_types_h(
            solver_stem, solver_namespace
        ),
        layout.generated_include
        + f"{solver_stem}_constants.h": constants_h.output_constants_h(
            solver_stem, solver_namespace, required_padding, args.ko
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
            solver_stem, solver_prefix
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
            solver_stem, solver_namespace, exec_or_library_name
        ),
        layout.pars
        + f"{solver_stem}_minkowski.par": cmdline_input_and_parfiles.generate_default_parfile(
            solver_stem, profile_name, required_padding, args.ko
        ),
        layout.tests
        + f"{solver_stem}_self_tests.cpp": self_tests_cpp.output_self_tests_cpp(
            solver_stem, solver_namespace
        ),
    }
    artifacts.update(
        cmake_helpers.output_CFunctions_function_prototypes_and_construct_CMakeLists(
            solver_name, solver_stem, solver_prefix, exec_or_library_name
        )
    )
    for relative_path, text in sorted(artifacts.items()):
        target = Path(args.project_dir) / relative_path
        target.parent.mkdir(parents=True, exist_ok=True)
        # C and C++ artifacts are clang-formatted, as every other NRPy
        # infrastructure formats what it emits; CMake and TOML are left as
        # written.
        with ConditionalFileUpdater(
            target,
            encoding="utf-8",
            do_format=target.suffix in (".h", ".hpp", ".c", ".cpp", ".cu"),
        ) as file:
            file.write(text)
    # The standalone-host header is a fixed source asset of this package, copied
    # verbatim exactly as BHaH copies simd_intrinsics.h in this same position.
    copy_files(
        package="nrpy.infrastructures.Dendro.standalone_host",
        filenames_list=["dendro_standalone_host.h"],
        project_dir=str(Path(args.project_dir) / layout.root),
        subdirectory="standalone_host",
    )

    EVOL, _AUXEVOL, _DIAG, _AUX = gri.GridFunction.gridfunction_lists()
    print(f"Finished generating {solver_name} in {args.project_dir}.")
    print(f"  profile: {profile_name}")
    print(f"  evolved variables: {len(EVOL)}")
    print(f"  finite-difference order: {args.fd_order}")
    print(f"  Kreiss-Oliger dissipation: {'enabled' if args.ko else 'disabled'}")
    print(f"  required ghost points: {roles.required_padding()}")
    print("Now build and run the generated self-tests with:")
    print(f"  cmake -S {args.project_dir}/Dendro-GR/{solver_name} -B build")
    print("  cmake --build build && ctest --test-dir build")


if __name__ == "__main__":
    main()
