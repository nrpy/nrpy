# Production Validation And Deferred Checks

> Define checks for complete generated Dendro applications and state current limits. · Status: provisional
> Up: [Dendro](index.md)

## Summary

Production examples generate complete Dendro applications. Existing historical
self-test emitters remain preserved, but production examples do not import
them.

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

The checked-in GitHub workflow has not yet been updated for this application
layout. Therefore workflow configuration does not currently prove these Dendro
generation, build, MPI, remesh, restart, or diagnostic checks.

No new test file, test case, doctest prompt, or stored oracle may be added
without express user permission. Runtime results belong in active review or CI
output, not as KB snapshots.

## Sources

- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - complete BSSN application generation.
- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - complete fCCZ4 application generation.
- [CMakeLists.py](../../../nrpy/infrastructures/Dendro/CMakeLists.py) - explicit generated source manifest.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - production evolution and service scheduling.
- [Constraints And Diagnostic Norms](constraints-and-diagnostic-norms.md) - momentum convention and comparable RMS diagnostics.
- [Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md) - remesh path and limits of native parity.
- [main.yml](../../../.github/workflows/main.yml) - currently configured CI jobs.
- [Code Test Policy](../../validation/code-test-policy.md) - permitted test changes and proof limits.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Generated Project CI](../../validation/generated-project-ci.md)
- Validates: [Project Assembly And Generating Functions](project-assembly-and-emitters.md)
