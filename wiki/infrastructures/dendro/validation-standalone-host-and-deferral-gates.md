# Production Validation And Deferred Checks

> Define checks for complete generated Dendro applications and state current limits. · Status: provisional
> Up: [Dendro](index.md)

## Summary

Production examples generate complete Dendro applications.

## Detail

Required generator checks include isolated static analysis, atomic registry
merge failure, BSSN and fCCZ4 import and registration, two byte-identical clean
generations, presence of all order-specific kernels, direct Python/CFunction/C++
name correspondence, final-state pointer indices, and rejection of native
`BSSN_GR` sources in each explicit CMake source list.

Required application checks configure and build both sibling applications with
Dendrolib. Both BSSN and fCCZ4 checks must generate W and chi variants.
MPI TwoPunctures runs must exercise finite-difference orders 4, 6,
and 8 through initialization and the first completed diagnostic output. Checks
must cover nontrivial finite fields, `alpha=W=sqrt(chi)` for both formulations,
Lambda initialization, algebraic projection, block-boundary continuity,
constraints, wave extraction, apparent horizons, forced remeshing, checkpoint,
and restart. W/chi checkpoint-formulation mismatches must be rejected.
Temporary direct numerical comparisons must cover interior stencils, block
boundaries, conformal-factor algebra, SSL, CAHD, and separate Ricci/RHS
coupling. Native-control comparisons must use the same initial data, post-remesh
state, puncture excision, momentum-component convention, and unique-node RMS;
conformal-factor volume-weighted RMS is a distinct diagnostic. Compare by physical time and check
the initial diagnostic before interpreting long evolutions. Record the evolved
conformal-factor choice, native full-psi lapse replacement, eta prescription,
KO strength, native CAKO state, and any post-merger CAKO switch before
attributing any later difference to the formulation.

Claim evidence:
- Claim: Complete BSSN and fCCZ4 qualification requires W and chi generation and rejection of incompatible restarts; BSSN native comparisons require matching initial data, excision, momentum components, and post-remesh unique-node RMS. Native comparisons must record the conformal-factor choice, full-psi lapse replacement, eta prescription, KO strength, CAKO state, and any post-merger CAKO switch; matching a parfile alone is insufficient.
- Role: normative rule
- Deciding authority: this page, `Required application checks`.
- Corroboration: `nrpy/examples/dendro_bssn.py` and `nrpy/examples/dendro_fccz4.py`, `parse_args`; `nrpy/infrastructures/Dendro/checkpoint.py`, formulation metadata; `nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py`, momentum lowering; `nrpy/infrastructures/Dendro/solver_context.py`, diagnostic scheduling and node reduction; `BSSN_GR/src/parameters.cpp`, native lapse and CAKO settings; `BSSN_GR/src/bssngr_main.cpp`, post-merger CAKO switch.

AMR parameter parity is not proof of identical remesh histories. The generated
optional Nyquist path uses the z-coordinate puncture separation, but native
Dendro-BSSN currently duplicates x separation into the z component of its
relative-position history. Compare the resulting meshes before interpreting
constraint differences as differences between evolution equations.

The `dendro-validation` workflow job runs
`nrpy/examples/tests/dendro_application_check.py` once per formulation and
configures these required checks: W and chi generation (byte-identical repeat
generation), complete builds, MPI TwoPunctures runs at FD4, FD6, and FD8 through
the first diagnostic output with the initial Hamiltonian-constraint norm ordered
by FD order (an ordering check, not a convergence test), finite diagnostics on
every run, constraint output, wave extraction with odd-m modes near zero and the
reflection relation C(l,-m) = (-1)^l conj C(l,m), apparent-horizon irreducible
masses matching the puncture ADM masses, a forced remesh whose evolved state and
node counts match a stored reference, checkpoint and byte-identical restore,
point reflection of the checkpointed puncture centers and their motion along the
puncture momenta, and rejection of both W/chi checkpoint-formulation mismatches.
The stored-reference comparison of run A's evolved constraint, ADM, and horizon
values is a regression check on the evolution, not a correctness proof. Lambda
initialization is asserted through the stored reference, which includes the
step-0 Lambda constraint (and, for fCCZ4, the Z4 diagnostics), and the algebraic
projection through checkpointing: the checkpoint writer refuses to write, and a
restore refuses to read, a checkpoint whose projection residual is nonfinite or
above tolerance. Halo exchange between ranks is covered by 1-, 3-, and 4-rank
runs whose diagnostics must agree within a relative tolerance; block boundaries
within a rank are covered only indirectly, through the numerical checks and the
stored reference. The job does not inspect field data, does not check
`alpha=W=sqrt(chi)` pointwise, and it runs none of the temporary direct
numerical comparisons or native-control comparisons above; those remain
review-time checks.

Claim evidence:
- Claim: `dendro-validation` configures the listed required application checks through `dendro_application_check.py`, covers the evolution and the forced remesh through a stored-reference regression check, covers Lambda initialization through the stored step-0 constraint row and the algebraic projection through the checkpoint write-time and restore-time residual checks, covers halo exchange through rank-count agreement and block boundaries within a rank only through the numerical checks and the stored reference, and omits field-data inspection, the pointwise initial-lapse check, the temporary direct numerical comparisons, and native-control comparisons.
- Role: CI behavior
- Deciding authority: [dendro_application_check.py](../../../nrpy/examples/tests/dendro_application_check.py), `Leg.run_variant`, `Leg.check_run_a`, `Leg.check_reference`, `Leg.check_orders`, `Leg.compare_runs`, `Leg.run_negatives`; [main.yml](../../../.github/workflows/main.yml), `dendro-validation`; [checkpoint.py](../../../nrpy/infrastructures/Dendro/checkpoint.py), `output_checkpoint_cpp`, write-time and restore-time projection-residual checks
- Corroboration: [Generated Project CI](../../validation/generated-project-ci.md), Dendro job description

No new test file, test case, doctest prompt, or stored oracle may be added
without express user permission. Runtime results belong in active review or CI
output, not as KB snapshots.

## Sources

- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - complete BSSN application generation.
- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - complete fCCZ4 application generation.
- [CMakeLists.py](../../../nrpy/infrastructures/Dendro/CMakeLists.py) - explicit generated source manifest.
- [dendro_application_check.py](../../../nrpy/examples/tests/dendro_application_check.py) - configured CI checks for both applications.
- [dendro_application_check_reference.py](../../../nrpy/examples/tests/dendro_application_check_reference.py) - stored run-A reference values.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - production evolution and service scheduling.
- [Constraints And Diagnostic Norms](constraints-and-diagnostic-norms.md) - momentum convention and comparable RMS diagnostics.
- [Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md) - remesh path and limits of native parity.
- [main.yml](../../../.github/workflows/main.yml) - currently configured CI jobs.
- [Code Test Policy](../../validation/code-test-policy.md) - permitted test changes and proof limits.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Generated Project CI](../../validation/generated-project-ci.md)
- Validates: [Project Assembly And Generating Functions](project-assembly-and-emitters.md)
