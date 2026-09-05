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

import nrpy.grid as gri
import nrpy.params as par
from nrpy.infrastructures.Dendro import registration as reg
from nrpy.infrastructures.Dendro.general_relativity import (
    diagnostics,
    initial_data,
    projection,
    rhs_eval,
)
from nrpy.infrastructures.Dendro.output_project import output_project
from nrpy.infrastructures.Dendro.runtime import parameters

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
    # fd_order 8 reaches five ghost points, above the max_proven_padding of 4
    # recorded in dendrolib_capabilities.json, so it is not offered here.
    parser.add_argument(
        "--fd-order",
        type=int,
        choices=(2, 4, 6),
        default=4,
        help="finite-difference order; 8 is capability-gated by "
        "dendrolib_capabilities.json (max_proven_padding 4)",
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
    par.set_parval_from_str("Dendro_scalar_type", "DendroScalar")
    par.set_parval_from_str("Dendro_enable_KreissOliger_dissipation", args.ko)

    #########################################################
    # Step 2: Register the generated C functions.  The right-hand side goes
    #         first: it registers the exact gridfunctions and physics
    #         CodeParameters through the shared fCCZ4 expression bundle, and
    #         records the ghost points the emitted operators reach.
    rhs_eval.register_CFunctions_rhs_eval(
        fd_order=args.fd_order,
        enable_KreissOliger_dissipation=args.ko,
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )

    # Minkowski initial data and the smooth analytic perturbation the
    # lifecycle gates evolve.  Both are NRPy-authored kernels.
    initial_data.register_CFunctions_minkowski_initial_data()
    initial_data.register_CFunctions_perturbation()

    # The smooth ADM conversion, the separate connection-initialization pass,
    # and the algebraic projection.  The conversion registers the ADM source
    # fields as AUXEVOL; the projection is scheduled after initial data and
    # after every accepted timestep by the host context.
    initial_data.register_CFunctions_initial_data_conversion(CoordSystem=CoordSystem)
    projection.register_CFunctions_projection(
        solver_namespace=solver_namespace, CoordSystem=CoordSystem
    )

    # The constraint diagnostics.  H_Z4 and the connection constraint are
    # registered as DIAG gridfunctions: they are recomputed from the evolved
    # state and are never checkpoint state.
    diagnostics.register_CFunctions_diagnostics(
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )

    # The parameter C functions come last, after every CodeParameter the
    # scientific kernels register is in the registry.
    parameters.register_CFunctions_parameters(solver_namespace, solver_stem)

    #########################################################
    # Step 3: Write the project to project_dir.
    output_project(
        project_dir=args.project_dir,
        solver_name=solver_name,
        solver_prefix=solver_prefix,
        solver_stem=solver_stem,
        solver_namespace=solver_namespace,
        exec_or_library_name=exec_or_library_name,
        profile_name=profile_name,
        generator_module="nrpy.examples.dendro_fccz4",
    )

    EVOL, _AUXEVOL, _DIAG, _AUX = gri.GridFunction.gridfunction_lists()
    print(f"Finished generating {solver_name} in {args.project_dir}.")
    print(f"  profile: {profile_name}")
    print(f"  evolved variables: {len(EVOL)}")
    print(f"  finite-difference order: {args.fd_order}")
    print(f"  Kreiss-Oliger dissipation: {'enabled' if args.ko else 'disabled'}")
    print(f"  required ghost points: {list(reg.required_padding())}")
    print("Now build and run the generated self-tests with:")
    print(f"  cmake -S {args.project_dir}/Dendro-GR/{solver_name} -B build")
    print("  cmake --build build && ctest --test-dir build")


if __name__ == "__main__":
    main()
