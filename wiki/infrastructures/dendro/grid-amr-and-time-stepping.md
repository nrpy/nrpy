# Octree Grid, AMR, And Time Stepping

> Explain Dendro blocks, data movement, RK stages, remeshing, and generated service scheduling. · Status: provisional
> Up: [Dendro](index.md)

## Summary

Dendrolib represents the adaptive mesh with balanced octants and provides
padded regular `ot::Block` regions for finite-difference kernels. Generated
solver context owns evolved vectors, six-component Ricci scratch storage,
ghost exchange, RK stages, AMR transfer, checkpointing, and diagnostics.
Binary-puncture grids use an analytic Dendro-GR seed for initial octree
construction, followed by TwoPunctures data for the evolved state.

## Detail

Zipped vectors follow octree degrees of freedom. Unzipped vectors provide
regular padded block storage. Before a Ricci/RHS traversal, solver context zips
stage state as needed, exchanges ghosts, unzips data, and fills physical
boundaries. It then visits each local block once and calls Ricci followed by
RHS. Ricci scratch requires no exchange because RHS consumes it immediately
within the same block. The generated `physical_boundary_ghosts` pass
extrapolates only exterior physical padding, after inter-block exchange and
unzip and before centered derivatives. It uses five interior points for FD4
and six for FD6/FD8, then fills faces, edges, and corners successively; this
avoids using stale exterior ghosts in Ricci, RHS, constraints, and wave
extraction.

Claim evidence:
- Claim: Physical exterior ghosts are filled after inter-block exchange and before centered derivatives, using five interior points for FD4 or six for FD6/FD8.
- Role: public/numerical contract
- Deciding authority: `nrpy/infrastructures/Dendro/general_relativity/physical_boundary_ghosts.py`, `register_CFunction_physical_boundary_ghosts`.
- Corroboration: `nrpy/infrastructures/Dendro/solver_context.py`, `output_solver_context_cpp`, calls the generated ghost fill before derivative kernels.

`Ctx::zip()` writes the locally owned continuous-Galerkin nodes; it does not
populate ghost nodes. Any newly zipped diagnostic field used by an element
interpolator must therefore complete `readFromGhostBegin/End` before
interpolation. Wave extraction follows this rule for both Psi4 components.
Without that exchange, interpolation can read uninitialized ghost storage even
when every pointwise block value is finite.

Floors and algebraic projection act directly on owned nodes of the zipped
evolved state, without a padded-block unzip/zip. After the final RK projection
and any remesh transfer, solver context exchanges evolved-state ghosts before
puncture tracking and field output. Intermediate RK stages receive their halo
exchange during the following RHS evaluation.

Claim evidence:
- Claim: Floors and algebraic projection use owned zipped nodes, and evolved-state ghosts are refreshed after the final projection and any remesh before puncture tracking and output.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::post_timestep` and `Ctx::evolve_excision_centers` within `output_solver_context_cpp`.
- Corroboration: `nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py` and `enforce_detgbar_equals_detghat_trAzero.py`, node-range kernel signatures; `nrpy/infrastructures/Dendro/main_cpp.py`, evolution/remesh/output order.

After initial-data conversion, each RK stage before the next exchange and RHS
evaluation, and AMR transfer, solver context floors `alpha` at `CHI_FLOOR` and
W at `sqrt(CHI_FLOOR)` or chi at `CHI_FLOOR`, then applies algebraic
projection. Remeshing transfers the fixed canonical evolved-field list.
Checkpoints record formulation, field names and order, emitted `CodeParameter`
values, iteration, time, puncture-center history, merger time, and whether
a checkpoint has been written after merger. Restore retains that state so
post-merger AMR uses the same coarsening factor as uninterrupted evolution when
AMR controls in the restart TOML are unchanged. Restore rejects an incompatible
formulation or field layout and preserves stored projected values. Keeping
puncture history across restart supports history-dependent AMR decisions.

Claim evidence:
- Claim: Restart preserves the checkpoint-written-after-merger state used to select the post-merger AMR coarsening factor.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/checkpoint.py`, `output_checkpoint_cpp`; `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::write_checkpt`, `Ctx::restore_checkpt`, and `Ctx::is_remesh` within `output_solver_context_cpp`.
- Corroboration: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`, schedules remeshing and checkpoint writes around time stepping.

The generated binary-black-hole path parses native `BH_WAMR` mode 4,
selected refinement variables, constant or causal mode-6 wavelet tolerance,
wavelet coarsening factors, puncture-centered level floors, post-merger remesh
cadence, and optional wave-zone Nyquist refinement. The analytic initial-grid
seed comes from Dendro-GR's `punctureDataPhysicalCoord`, with its source and
MIT license embedded in the generated entry point. A separate single-rank
`--tpid` run computes TwoPunctures coefficients once. Fresh evolution loads
those coefficients before evolved-state conversion; after any initial-grid
remesh, the solver reconstructs that state from the same data. Grid
construction remains distinct from evolved initial data. Unsupported refinement or
tolerance modes fail at startup instead of silently substituting a constant
tolerance. These choices
target a comparable grid structure under the same parameter file; they do
not establish bitwise identity of the two remesh histories. The generated
Nyquist path computes all three components of puncture separation from their
corresponding coordinates. Native Dendro-BSSN currently duplicates x separation
into the z component of its relative-position history, so enabling Nyquist
refinement can produce different remesh decisions even with the same parameters.

Claim evidence:
- Claim: The generated binary-puncture path parses native-style wavelet/geometric AMR controls, loads separately solved TwoPunctures data for evolution, and uses the licensed native analytic octree seed; its Nyquist history uses the z-coordinate separation, unlike the native path's duplicated x separation.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`; `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::get_wtol_function` and `Ctx::is_remesh` within `output_solver_context_cpp`.
- Corroboration: `BSSN_GR/src/grUtils.cpp`, `punctureDataPhysicalCoord`; `BSSN_GR/src/dataUtils.cpp`, `calculate_relative_position_history` and `isRemeshBH`.

For compatibility with Dendro-GR BSSN_GR remeshing, the wavelet refinement
test sees the Gamma-driver auxiliary field `betU` scaled by 4/3. BSSN_GR
evolves the shift as `d_t beta^i = (3/4) B^i` with `BSSN_LAMBDA_F = (1, 0)`,
while BSSN and fCCZ4 here use `GammaDriving2ndOrder_Covariant__Hatted`, which
evolves `d_t beta^i = B^i`; with `BSSN_LAMBDA = (1, 1, 1, 1)` and the same
damping `eta` the two auxiliary fields satisfy
`B^i(BSSN_GR) = (4/3) B^i(NRPy)`. BSSN_GR's default damping is the
radius-dependent RIT profile (2.0 near the origin, about 0.25 beyond
r ≈ 60), while NRPy uses the constant `eta`, so the relation holds only where
the two agree. Because the test compares every refinement field's wavelet
coefficients with one tolerance, an unscaled `B` can coarsen earlier than
BSSN_GR where `B` determines the refinement decision. `Ctx::is_remesh`
multiplies `betU0`–`betU2` by 4/3 in the unzipped work vector, which every
other user refills with `unzip` before reading, after the physical-boundary
fill and before `isReMeshUnzip`; the evolved state is not changed.

Claim evidence:
- Claim: `Ctx::is_remesh` scales `betU` by 4/3 in its unzipped work vector before the wavelet refinement test, which puts the tested auxiliary shift-driver field in BSSN_GR's normalization (equal to BSSN_GR's `B` for `BSSN_LAMBDA_F = (1, 0)`, `BSSN_LAMBDA = (1, 1, 1, 1)` and equal `eta`); the evolved state is unchanged.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::is_remesh` within `output_solver_context_cpp`.
- Corroboration: `nrpy/equations/general_relativity/BSSN_gauge_RHSs.py`, `GammaDriving2ndOrder_Covariant__Hatted`; `BSSN_GR/src/bssneqs_SSL_HD_dxsq.cpp`, `b_rhs` and `B_rhs`; `BSSN_GR/src/rhs.cpp`, RIT `eta` profile.

Puncture-center tracking reads the evolved `vetU` shift components, not the
`betU` auxiliary shift-driver components. The tracked centers set excision
regions and puncture-centered AMR. Diagnostics are scheduled after a remesh
and transfer, if one occurred, and after puncture-center advancement.

Claim evidence:
- Claim: The generated binary-black-hole path tracks centers using `vetU` and schedules diagnostics after remeshing and puncture advancement.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`; `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::is_remesh`, `Ctx::evolve_excision_centers`, and `Ctx::diagnostic_output` within `output_solver_context_cpp`.
- Corroboration: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`, schedules remesh, puncture advancement, and diagnostic output in that order.

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
- [floor_the_lapse_and_conformal_factor.py](../../../nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py) - representation-dependent lapse and conformal-factor floors.
- [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/Dendro/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - owned-node algebraic projection.
- [physical_boundary_ghosts.py](../../../nrpy/infrastructures/Dendro/general_relativity/physical_boundary_ghosts.py) - physical exterior-padding extrapolation.
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - analytic seed, AMR parameters, and scheduling.
- [Dendro-GR dataUtils.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/dataUtils.cpp) - native black-hole refinement path.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Gridfunctions, Naming, And Loops](gridfunctions-naming-and-loops.md)
- See also: [Constraints And Diagnostic Norms](constraints-and-diagnostic-norms.md)
- Validated by: [Production Validation And Deferred Checks](validation-standalone-host-and-deferral-gates.md)
