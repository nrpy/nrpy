# Octree Grid, AMR, And Time Stepping

> Explain Dendro blocks, data movement, RK stages, remeshing, and generated service scheduling. · Status: confirmed
> Up: [Dendro](index.md)

## Summary

Dendrolib represents the adaptive mesh with balanced octants and provides
padded regular `ot::Block` regions for finite-difference kernels. Generated
solver context owns evolved vectors, six-component Ricci scratch storage,
ghost exchange, RK stages, AMR transfer, checkpointing, and diagnostics.

## Detail

Zipped vectors follow octree degrees of freedom. Unzipped vectors provide
regular padded block storage. Before a Ricci/RHS traversal, solver context zips
stage state as needed, exchanges ghosts, unzips data, and fills physical
boundaries. It then visits each local block once and calls Ricci followed by
RHS. Ricci scratch requires no exchange because RHS consumes it immediately
within the same block.

`Ctx::zip()` writes the locally owned continuous-Galerkin nodes; it does not
populate ghost nodes. Any newly zipped diagnostic field used by an element
interpolator must therefore complete `readFromGhostBegin/End` before
interpolation. Wave extraction follows this rule for both Psi4 components.
Without that exchange, interpolation can read uninitialized ghost storage even
when every pointwise block value is finite.

After initial-data conversion, each RK stage before the next exchange and RHS
evaluation, and AMR transfer, solver context floors `alpha` at `CHI_FLOOR` and
W at `sqrt(CHI_FLOOR)`, then applies algebraic projection. Remeshing uses the
fixed canonical field list for transfer. Checkpoints record formulation, field
names and order, parameters, iteration, and time. Restore rejects an
incompatible formulation or field layout and preserves stored projected values.

Initial-grid convergence may reduce the active communicator while retaining
all global MPI ranks. Every global rank must still enter the Dendro remesh
decision collectives; only block-local geometry work is conditional on
`mesh->isActive()`. A rank-local early return before `isReMeshUnzip` mismatches
collectives when the active communicator later expands.

Generated scheduling includes physical boundaries, BSSN-to-ADM conversion,
constraint norms, apparent-horizon calls, wave extraction, ADM quantities,
field output, timing, checkpoint, and restart. TwoPunctures is application-local
initial data rather than a call into chi-specific `BSSN_GR` code.

## Sources

- [Dendrolib block.h](https://github.com/paralab/Dendro-5.01/blob/master/include/block.h) - `ot::Block` geometry.
- [Dendrolib mesh.h](https://github.com/paralab/Dendro-5.01/blob/master/include/mesh.h) - zip, unzip, remesh, and intergrid transfer.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - generated context and service scheduling.
- [checkpoint.py](../../../nrpy/infrastructures/Dendro/checkpoint.py) - checkpoint and restore generation.
- [physical_boundary.py](../../../nrpy/infrastructures/Dendro/general_relativity/physical_boundary.py) - boundary kernel registration.
- [BSSN_to_ADM.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN_to_ADM.py) - geometry conversion.
- [floor_the_lapse_and_conformal_factor.py](../../../nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py) - W-form lapse and conformal-factor floors.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Gridfunctions, Naming, And Loops](gridfunctions-naming-and-loops.md)
- Validated by: [Production Validation And Deferred Checks](validation-standalone-host-and-deferral-gates.md)
