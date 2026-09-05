#!/usr/bin/env python3
"""
NRPy Dendro fCCZ4 example generator (whitepaper section 11).

Generates the complete NRPy-authored fCCZ4 module for Dendro-GR from a fresh
process profile.  The registration sequence (section 11.3) is fixed:

1. Import modules that register the existing shared NRPy parameters.
2. Enter a generation transaction in a fresh process.
3. Set ``Infrastructure="Dendro"`` before any fCCZ4 gridfunction registration.
4. Set ``fp_type``, ``Dendro_scalar_type``, ``parallelization``, and
   ``fd_order`` before code generation.
5. Register the host parameter table; build the shared fCCZ4 expression
   bundle, which registers the exact gridfunctions and physics CodeParameters.
6. Build and register the direct-FD RHS CFunctions (PR 5), the Minkowski
   initial data CFunctions (PR 7), the ADM conversion, the connection
   initialization and the algebraic projection (PR 8), and the constraint
   diagnostics (PR 9).
7. Register the generated parameter CFunctions after all use closure is
   known (section 6.4).
8. Freeze, hash, and seal all registries.
9. Export only from the frozen snapshot (section 12).

Run as a module:

    python -m nrpy.examples.dendro_fccz4 \
        --project-dir project/dendro_fccz4 --fd-order 4 --no-ko
"""

from __future__ import annotations

import argparse
from pathlib import Path

# pylint: disable=unused-import
import nrpy.finite_difference  # noqa: F401  (registers fd_order)
import nrpy.grid as gri
import nrpy.params as par
from nrpy.infrastructures.Dendro import (  # noqa: F401
    generation_parameters,
)
from nrpy.infrastructures.Dendro import registration as reg
from nrpy.infrastructures.Dendro.freeze import (
    assert_mutable_registries_match,
    freeze_nrpy_dendro_environment,
)
from nrpy.infrastructures.Dendro.general_relativity import (
    diagnostics,
    initial_data,
    projection,
    rhs_eval,
)
from nrpy.infrastructures.Dendro.output_project import output_project
from nrpy.infrastructures.Dendro.runtime import parameters
from nrpy.infrastructures.Dendro.transaction import (
    dendro_generation_transaction,
)


def parse_args() -> argparse.Namespace:
    """
    Parse the command line.

    :return: The parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description="Generate an NRPy-authored fCCZ4 module for Dendro-GR"
    )
    parser.add_argument(
        "--project-dir", type=Path, default=Path("project/dendro_fccz4")
    )
    parser.add_argument("--fd-order", type=int, choices=(2, 4, 6, 8), default=4)
    parser.add_argument("--ko", action=argparse.BooleanOptionalAction, default=False)
    # No --parallelization flag: whitepaper section 7.2 pins the qualified CPU
    # profile to serial point loops (the kernel runs inside Dendro's own block
    # traversal), and the builders assert that.  Offering a flag the generator
    # discards would silently produce an unqualified configuration (section 9.3).
    return parser.parse_args()


def register_profile(args: argparse.Namespace) -> None:
    """
    Register the generation profile into the NRPy registries.

    Every scientific choice is written to a registered NRPy parameter before
    any expression is built (sections 4.1/6.1/11.2/11.3); no parallel
    configuration object exists.

    :param args: The parsed command line.
    """
    par.set_parval_from_str("Infrastructure", "Dendro")
    par.set_parval_from_str("fp_type", "double")
    # Section 7.2: the qualified CPU profile emits serial point loops (the
    # kernel runs inside Dendro's own block traversal).  The NRPy default is
    # "openmp", so pin it here; the builders assert it.
    par.set_parval_from_str("parallelization", "none")
    par.set_parval_from_str("fd_order", args.fd_order)
    par.set_parval_from_str("EvolvedConformalFactor_cf", "chi")
    par.set_parval_from_str("detgbarOverdetghat_equals_one", True)

    par.set_parval_from_str("Dendro_module_name", "FCCZ4_GR")
    par.set_parval_from_str("Dendro_scalar_type", "DendroScalar")
    par.set_parval_from_str("Dendro_derivative_backend", "full_stencil")
    par.set_parval_from_str("Dendro_enable_KO", args.ko)
    par.set_parval_from_str("Dendro_fccz4_CoordSystem", "Cartesian")
    par.set_parval_from_str("Dendro_fccz4_LapseEvolutionOption", "OnePlusLog")
    par.set_parval_from_str(
        "Dendro_fccz4_ShiftEvolutionOption",
        "GammaDriving2ndOrder_Covariant__Hatted",
    )

    # PR 5: direct-FD RHS CFunctions (builds the shared expression bundle,
    # which registers the exact gridfunctions and physics CodeParameters).
    build = rhs_eval.build_fccz4_rhs(
        fd_order=args.fd_order,
        enable_KO=args.ko,
    )
    reg.register_dendro_CFunction(
        role="rhs_block",
        entry_point=True,
        name=rhs_eval.BLOCK_CFUNCTION,
        desc="Per-block direct-FD fCCZ4 RHS (25 fields).",
        subdirectory="generated/src/rhs",
        params=build.block_params,
        body=build.block_body,
    )
    reg.register_dendro_CFunction(
        role="rhs",
        entry_point=True,
        calls=(rhs_eval.BLOCK_CFUNCTION,),
        name=rhs_eval.GLOBAL_CFUNCTION,
        desc="All-block direct-FD fCCZ4 RHS (NRPy block loop).",
        subdirectory="generated/src/rhs",
        params=build.global_params,
        body=build.global_body,
    )
    reg.register_dendro_CFunction(
        role="rhs_flat_block",
        entry_point=True,
        calls=(rhs_eval.BLOCK_CFUNCTION,),
        name=rhs_eval.FLAT_BLOCK_CFUNCTION,
        desc="LTS flat-block adapter (same numerical body, flat layout).",
        subdirectory="generated/src/rhs",
        params=build.flat_block_params,
        body=build.flat_block_body,
    )

    # Record the single-derivation audit records (section 15.4 operator
    # manifest, provenance) in the builder sidecar; freeze validates and
    # carries them into the snapshot extras.
    dendro_extras = par.glb_extras_dict.setdefault("Dendro", {})
    dendro_extras["builder_extras"] = {
        "operator_manifest": build.operator_manifest,
        "rhs_canonical": dict(build.canonical_expression_digests),
        "provenance": {
            "rhs_symbols": list(build.rhs_symbols),
            "upwind_control_fields": list(build.upwind_control_fields),
            "used_codeparameters": list(build.used_codeparameters),
            "padding": list(build.padding),
        },
    }

    # PR 7: Minkowski initial data CFunctions (section 14.1) and the smooth
    # analytic perturbation the lifecycle gates evolve (section 16.8).  Both
    # are NRPy-authored kernels, so neither lives in a fixed template.
    initial_data.register_minkowski_CFunctions()
    initial_data.register_perturbation_CFunctions()

    # PR 8: the smooth ADM conversion (section 14.2), the separate connection
    # initialization pass (section 14.3), and the algebraic projection
    # (section 14.5).  The conversion registers the ADM source fields as
    # AUXEVOL; the projection is scheduled after initial data and after every
    # accepted timestep by the host context.
    initial_data.register_initial_data_conversion_CFunctions()
    projection.register_projection_CFunctions()

    # PR 9: the first diagnostic set (section 14.6).  H_Z4 and the connection
    # constraint are registered as DIAG gridfunctions: they are recomputed
    # from the evolved state and are never checkpoint state.
    diagnostics.register_diagnostics_CFunctions()

    # Section 6.4: parameter CFunctions last, after all use closure is known.
    parameters.register_parameter_CFunctions_last()


def main() -> None:
    """Generate the complete Dendro fCCZ4 project."""
    args = parse_args()
    module_abi = ""
    with dendro_generation_transaction() as tx:
        register_profile(args)
        # Section 11.4 required pre-freeze assertions.
        assert par.parval_from_str("Infrastructure") == "Dendro"
        assert par.parval_from_str("fp_type") == "double"
        assert par.parval_from_str("Dendro_scalar_type") == "DendroScalar"
        assert par.parval_from_str("Dendro_derivative_backend") == "full_stencil"
        assert par.parval_from_str("fd_order") == args.fd_order
        EVOL, _AUXEVOL, _DIAG, _AUX = gri.GridFunction.gridfunction_lists()
        assert len(EVOL) == 25
        assert "Theta_fCCZ4" in EVOL

        # Section 9.9 steps 4/9/11 are enforced inside freeze from these
        # builder records; requiring them here keeps the production path from
        # ever freezing an unchecked profile.
        builder_records = par.glb_extras_dict["Dendro"]["builder_extras"]
        assert "operator_manifest" in builder_records
        assert "rhs_canonical" in builder_records
        assert builder_records["provenance"]["rhs_symbols"]

        snapshot = freeze_nrpy_dendro_environment(
            profile_name="fccz4_cartesian_vacuum", tx=tx
        )
        assert_mutable_registries_match(snapshot)
        output_project(snapshot=snapshot, project_dir=args.project_dir)
        # Section 4.4 phase machine: the export advanced FROZEN -> EMITTED ->
        # VERIFIED only after the staged project verified against the hashes.
        tx.advance_to("EMITTED")
        tx.advance_to("VERIFIED")
        module_abi = snapshot.hashes.module_abi_hash
    print("NRPy Dendro generation complete")
    print("  module: FCCZ4_GR")
    print("  profile: fccz4_cartesian_vacuum")
    print("  evolved variables: 25")
    print(f"  finite-difference order: {args.fd_order}")
    print("  KO owner: nrpy_full_stencil" if args.ko else "  KO owner: disabled")
    print(f"  output: {args.project_dir}")
    print(f"  module ABI: {module_abi}")


if __name__ == "__main__":
    main()
