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
Dendrolib. MPI TwoPunctures runs must exercise finite-difference orders 4, 6,
and 8 through initialization and the first completed diagnostic output. Checks
must cover nontrivial finite fields, `alpha=W`, Lambda initialization,
algebraic projection, block-boundary continuity, constraints, wave extraction,
apparent horizons, forced remeshing, checkpoint, and restart. Temporary direct
numerical comparisons must cover interior stencils, block boundaries,
conformal-factor algebra, SSL, CAHD, and separate Ricci/RHS coupling.

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
- [main.yml](../../../.github/workflows/main.yml) - currently configured CI jobs.
- [Code Test Policy](../../validation/code-test-policy.md) - permitted test changes and proof limits.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Generated Project CI](../../validation/generated-project-ci.md)
- Validates: [Project Assembly And Generating Functions](project-assembly-and-emitters.md)
