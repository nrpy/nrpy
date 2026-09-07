"""
Generate an NRPy-authored BSSN solver for Dendro-GR.

This is the second formulation lowered through the Dendro infrastructure. It
exists as much to test the infrastructure as to produce a solver: adding it
required the formulation-agnostic lowering to be extracted into
``block_kernel_helpers``, but no existing emitter changed behaviour and the fCCZ4
output is unaffected.

The emitted names follow Dendro-GR's own BSSN solver rather than NRPy's
vocabulary: the solver directory is ``BSSN_GR``, CMake variables carry the
``BSSN_`` prefix, the object library is ``bssn_common``, the executable is
``bssnSolver``, the namespace is ``bssn``, and the context source is
``bssnCtx.cpp``.

Run as a module:

    python -m nrpy.examples.dendro_bssn \
        --project-dir project/dendro_bssn --fd-order 4 --no-ko

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
    Dendro_include_header,
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

# Dendro-GR's own BSSN solver: directory BSSN_GR, namespace bssn, sources
# bssnCtx.cpp / bssn_constraints.cpp, object library bssn_common, executable
# bssnSolver.  These are arguments to the infrastructure, not registered
# CodeParameters.
solver_name = "BSSN_GR"
solver_prefix = "BSSN"
solver_stem = "bssn"
solver_namespace = "bssn"
exec_or_library_name = "bssnSolver"
profile_name = "bssn_cartesian_vacuum"

CoordSystem = "Cartesian"
LapseEvolutionOption = "OnePlusLog"
ShiftEvolutionOption = "GammaDriving2ndOrder_Covariant__Hatted"


def parse_args() -> argparse.Namespace:
    """
    Parse the command line.

    :return: The parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description="Generate an NRPy-authored BSSN solver for Dendro-GR"
    )
    parser.add_argument("--project-dir", default=os.path.join("project", "dendro_bssn"))
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
    # builders assert that.
    return parser.parse_args()


def main() -> None:
    """Generate the complete Dendro BSSN project."""
    args = parse_args()

    #########################################################
    # Step 1: Set the generation profile.
    par.set_parval_from_str("Infrastructure", "Dendro")
    par.set_parval_from_str("fp_type", "double")
    # The qualified CPU profile emits serial point loops.  The NRPy default is
    # "openmp", so pin it here; the builders assert it.
    par.set_parval_from_str("parallelization", "none")
    par.set_parval_from_str("fd_order", args.fd_order)
    par.set_parval_from_str("EvolvedConformalFactor_cf", "W")
    par.set_parval_from_str("detgbarOverdetghat_equals_one", True)

    #########################################################
    # Step 2: Register the generated C functions.  The right-hand side goes
    #         first: it registers the exact gridfunctions and physics
    #         CodeParameters, and records the ghost points the emitted
    #         operators reach.
    rhs_eval.register_CFunctions_rhs_eval(
        solver_stem=solver_stem,
        fd_order=args.fd_order,
        enable_KreissOliger_dissipation=args.ko,
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )

    # Minkowski initial data and the smooth analytic perturbation the
    # lifecycle gates evolve.  Both are formulation-agnostic: they write every
    # registered EVOL field to its registered asymptotic value.
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

    # The constraint diagnostics: the Hamiltonian constraint and the three
    # momentum constraint components, registered as DIAG gridfunctions.
    constraints_eval.register_CFunctions_constraints_eval(
        solver_stem=solver_stem, CoordSystem=CoordSystem
    )

    # The parameter C functions come last, after every CodeParameter the
    # scientific kernels register is in the registry.
    CodeParameters.register_CFunctions_parameters(solver_namespace, solver_stem)

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
        + f"{solver_stem}_defines.h": Dendro_include_header.output_include_header(
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
