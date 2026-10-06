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
regular padded block storage. Before a Ricci/RHS traversal, solver context
exchanges ghosts and unzips the evolved state once. It then visits each local
block once. For each block it fills the exterior physical padding, calls Ricci
followed by RHS, and finally calls `physical_boundary`, which replaces the
right-hand sides on the physical-face nodes (described below). `Ctx::rhs_blkwise`,
the entry for a requested list of blocks, repeats the same per-block sequence.
Ricci scratch requires no exchange because RHS consumes it immediately within
the same block. The generated `physical_boundary_ghosts` pass extrapolates only
exterior physical padding, after inter-block exchange and unzip and before
centered derivatives. It uses five interior points for FD4
and six for FD6/FD8, then fills faces, edges, and corners successively; this
avoids using stale exterior ghosts in Ricci, RHS, constraints, and wave
extraction.

Claim evidence:
- Claim: Physical exterior ghosts are filled after inter-block exchange and before centered derivatives, using five interior points for FD4 or six for FD6/FD8.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/general_relativity/physical_boundary_ghosts.py`, `register_CFunction_physical_boundary_ghosts`.
- Corroboration: `nrpy/infrastructures/Dendro/solver_context.py`, `output_solver_context_cpp`, calls the generated ghost fill before derivative kernels.

After each block's Ricci and RHS evaluation, `physical_boundary` replaces the
right-hand side of every evolved field at the nodes on the physical faces of the
domain by an outgoing-radiation condition:

```text
d_t f = -(x^i d_i f + n_f (f - f_inf)) / r,    r = sqrt(x^2 + y^2 + z^2).
```

The derivative `d_i f` is the centered finite-difference derivative of the
current order, which reads the extrapolated exterior padding. The exponent is
`n_f = 2` for the `lambdaU` and `aDD` components and `n_f = 1` for every other
evolved field, and the asymptotic value is `f_inf = 1` for `alpha` and `cf` and 0
for every other field. The condition has no field-dependent wave speed and acts
on face nodes only; interior nodes keep the equation's right-hand side. The
padding is a polynomial extrapolation from five (FD4) or six (FD6, FD8) interior
nodes, and the generated code states that six points give fourth-order boundary
second derivatives at FD6 and FD8.

Claim evidence:
- Claim: After each block's Ricci and RHS evaluation `Ctx::rhs` (and `Ctx::rhs_blkwise`) calls `physical_boundary`, which sets the right-hand side of every evolved field on physical-face nodes to `-(x^i d_i f + n_f (f - f_inf)) / r` with `r` the coordinate radius, `n_f = 2` for `lambdaU` and `aDD` and 1 otherwise, `f_inf = 1` for `alpha` and `cf` and 0 otherwise, and the centered finite-difference derivative of the current order; the kernel contains no field-dependent wave speed.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/general_relativity/physical_boundary.py`, `register_CFunction_physical_boundary`; `nrpy/infrastructures/Dendro/state_h.py`, `register_canonical_gridfunctions` (the `f_infinity` values); `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::rhs` and `Ctx::rhs_blkwise` within `output_solver_context_cpp`.
- Corroboration: `nrpy/infrastructures/Dendro/general_relativity/physical_boundary_ghosts.py`, `register_CFunction_physical_boundary_ghosts`, the extrapolation point counts and the stated boundary accuracy.

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
AMR controls in the restart TOML are unchanged. Restore preserves stored
projected values; the tests it applies are listed below. Keeping puncture
history across restart supports history-dependent AMR decisions.

A restore is accepted only when the checkpoint metadata passes every one of
these tests. A rejected checkpoint prints `Checkpoint metadata does not match
<formulation>`, whichever metadata test failed (a mismatch of the element order
also prints `Checkpoint element order N differs from BSSN_ELE_ORDER M`), and the
run then aborts with the generic `checkpoint restore failed` message:

- The stored formulation, field count and names, and parameter-name list equal
  the generated solver's.
- Every emitted `CodeParameter` equals its stored value exactly. For BSSN these
  are `C_CAHD`, both Kreiss-Oliger strengths (file key `KO_DISS_SIGMA`), `SSL_h`,
  `SSL_sigma`, `chi_floor`, and `eta`, plus `CFL_FACTOR`, `YBS_chi`, and
  `C_YBS_mom` when generated, and `kappa1` and `kappa2` for fCCZ4. Changing
  `ETA_CONST`, `KO_DISS_SIGMA`, `BSSN_CAHD_C`, `BSSN_SSL_H`, `BSSN_SSL_SIGMA`, or
  `CHI_FLOOR` in the restart file is therefore rejected.
- The stored domain bounds equal the current bounds.
- The stored time and time step are finite and the time step is positive, the
  stored projected algebraic residual is finite and at most 1e-10, and the
  puncture history is finite, strictly increasing in time, and ends at the stored
  time and centers.
- The stored element order is 4, 6, or 8 and equals `BSSN_ELE_ORDER`.

The residual is the global maximum over nodes of `|det(gammabar) - 1|` and
`|gammabar^ij Abar_ij|`; it is infinite when any node has a nonpositive or
nonfinite determinant or a nonfinite trace. After the metadata passes, the
residual of the restored state must again be at most 1e-10, and the launch must provide at least as many MPI ranks as the stored
active communicator, because the mesh is rebuilt on that many ranks with the
stored element order. Host keys that are not emitted `CodeParameter` values, other
than the domain bounds and `BSSN_ELE_ORDER`, are not compared (for example the
mesh, output, and extraction keys).

Claim evidence:
- Claim: A restore is accepted only when the stored formulation, field list, parameter names, every emitted `CodeParameter` value, domain bounds, time data, puncture history, element order (4, 6, or 8, and equal to `BSSN_ELE_ORDER`), and projected algebraic residual (at most 1e-10, rechecked on the restored state, and infinite for a nonpositive or nonfinite determinant or a nonfinite trace) pass, and when the launch has at least the stored number of active MPI ranks; a metadata failure prints `Checkpoint metadata does not match <formulation>` (an element-order mismatch also prints the stored and the running order) and the run aborts with `checkpoint restore failed`.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/checkpoint.py`, `output_checkpoint_cpp` (the metadata tests and the rank-count test); `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::restore_checkpt` within `output_solver_context_cpp` (the residual recheck).
- Corroboration: `nrpy/examples/tests/dendro_application_check.py`, `Leg.run_negatives`, the W/chi cross-restore rejection, and `Leg.run_variant`, the restore of an unchanged run; no CI case exercises the element-order mismatch or the refusal of a checkpoint whose trace is nonfinite.

Claim evidence:
- Claim: Restart preserves the checkpoint-written-after-merger state used to select the post-merger AMR coarsening factor.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/checkpoint.py`, `output_checkpoint_cpp`; `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::write_checkpt`, `Ctx::restore_checkpt`, and `Ctx::is_remesh` within `output_solver_context_cpp`.
- Corroboration: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`, schedules remeshing and checkpoint writes around time stepping.

New RK4 timesteps use `BSSN_CFL_FACTOR` times the smallest physical axis
spacing at the current global maximum octree level. The driver recomputes
this spacing after initial-grid convergence and after evolution remeshing;
`dx_min` reports that same minimum. Independent domain bounds and CFL must
produce finite positive widths and a finite positive CFL factor. This rule
preserves cubic-domain timesteps but reduces timesteps when Y or Z is finer
than X.

With `BSSN_RESTORE_SOLVER = 1` and no checkpoint metadata file in either slot
under `BSSN_CHKPT_FILE_PREFIX`, rank 0 prints a warning naming that prefix and
evolution starts from the initial data exactly as with `BSSN_RESTORE_SOLVER = 0`,
including loading the TwoPunctures coefficients and appending to any existing
diagnostic files. Rank 0 alone inspects the metadata files and broadcasts the
slot, so all ranks agree. A parameter file that always sets
`BSSN_RESTORE_SOLVER = 1` therefore serves the first submission and every
resubmission of a job; a wrong working directory or prefix also starts a fresh
evolution, and the warning is the only signal. Existing metadata that cannot be
restored (for example
an incompatible formulation or field layout, or a missing octree or state file)
remains a fatal error. The separate `--tpid` mode still generates initial-data
coefficients without requiring a checkpoint. A valid restore retains its stored
timestep; the driver rejects a timestep above the current minimum-axis CFL
bound before initializing RK4. Registered nonfinite floating parameters are
rejected before puncture-file access.

Claim evidence:
- Claim: New timesteps use the minimum physical axis spacing; a restore request with checkpoint metadata preserves a stored timestep only if it satisfies the current CFL bound, and a restore request without checkpoint metadata in either slot prints a warning on rank 0 and starts from the initial data as with `BSSN_RESTORE_SOLVER = 0`, while unreadable existing metadata stays fatal. The driver checks generated parameter validation before loading puncture data.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`.
- Corroboration: `nrpy/infrastructures/Dendro/CodeParameters.py`, generated parameter validation; Dendrolib `Block::computeDx`, `computeDy`, and `computeDz` define the physical axis spacings; `Leg.run_negatives` in `nrpy/examples/tests/dendro_application_check.py` runs the restore-without-metadata case; Dendro-GR `BSSN_GR/src/bssnCtx.cpp`, `BSSNCtx::initialize` and `BSSNCtx::restore_checkpt`, continue as if `BSSN_RESTORE_SOLVER` were false when no checkpoint files are found.

The generated driver restores only the metadata slot (index 0 or 1) with the
newest file modification time, and has no fallback to the other slot. It writes a
checkpoint after the step's output at every nonzero multiple of
`BSSN_CHECKPT_FREQ`, into slot `(step / BSSN_CHECKPT_FREQ) mod 2`, and the loop
ends without a final write, so a run that ends between multiples leaves no
checkpoint of its last steps. Within one write the octree and state files are
published first and the metadata file last, and the horizon search-state file
follows. Native `BSSN_GR` instead honors an explicit `BSSN_RESTORE_CHECKPT_SLOT` when that
slot exists (a missing slot falls back to the steps below),
otherwise reads the `.latest` sentinel it writes after each complete normal
checkpoint, otherwise compares the step numbers of the slots, restores the other
normal slot when a slot it chose itself is incomplete, and aborts when the
metadata, octree, or state file of the selected existing slot cannot be read. Its merger snapshot in
slot 3 is never named by `.latest` and is restorable only through an explicit
slot request.

Claim evidence:
- Claim: The generated driver restores only the metadata slot with the newest modification time and falls back to no other slot; it writes checkpoints only at nonzero multiples of `BSSN_CHECKPT_FREQ`, into slot `(step / BSSN_CHECKPT_FREQ) mod 2`, publishing the metadata after the octree and state files and then writing the horizon search-state file, with no write after the last step. Native `BSSN_GR` restores the slot named by `BSSN_RESTORE_CHECKPT_SLOT` when it exists (a missing slot falls back to auto-detection), otherwise the one the `.latest` sentinel names, otherwise the newer by step number, falls back to the other normal slot when a slot it chose itself is incomplete, aborts when an existing selected slot cannot be read, and writes its merger snapshot in slot 3 without publishing it in `.latest`.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp` (slot choice, write-then-evolve loop); `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::write_checkpt`; `nrpy/infrastructures/Dendro/checkpoint.py`, `output_checkpoint_cpp` (file publication order); Dendro-GR `BSSN_GR/src/bssnCtx.cpp`, `BSSNCtx::write_checkpt` and `BSSNCtx::restore_checkpt`.
- Corroboration: none available; no CI case checks slot choice, the write schedule, or the native restore order.

When the apparent-horizon finder is enabled (`AEH_SOLVER_FREQ > 0`), each
checkpoint also writes the finder's search state to
`<BSSN_CHKPT_FILE_PREFIX>_aeh_solver_checkpt-cp<index>.json`, with the 0/1 index
of the solver checkpoint, as Dendro-GR `BSSN_GR` names it. The search state is
the horizon count, the binary-black-hole flag, the previous three horizon
shapes, centers, times and radii, and the active, failure and fixed-radius-guess
flags; it seeds the next horizon find. Restore reads it back. A checkpoint
written without a horizon file leaves the finder in its initial state.
Dendrolib's `create_checkpoint` returns nothing: when it cannot open the file it
prints `file open failed for BAH checkpoint` to standard output and returns, and
the solver does not check a failed write. The file for each index is overwritten
in place and the index alternates, so after a failed open the file from two
checkpoints earlier remains, and a later restore reads it without a message and
seeds the finder with that older search state instead of its initial one.
`restore_checkpoint` prints `file open failed! Could not restore AH solver!` and
returns when the file is missing, and a file it cannot parse makes it throw. `BSSN_GR`
also writes a one-time merger checkpoint with index 3, which the generated
solver does not. Because a checkpoint is written after its step's output, a
restored run skips the initial output instead of repeating that step's output
rows and horizon find. The finder writes a `BHaHAHA_diagnostics` row for each
successfully found horizon only on horizon-find steps (multiples of
`AEH_SOLVER_FREQ`) that are also multiples of `BSSN_IO_OUTPUT_FREQ`; with
`BSSN_IO_OUTPUT_FREQ = 0` it writes none, although horizons are still found and
their search state is checkpointed.

Claim evidence:
- Claim: With the apparent-horizon finder enabled, each checkpoint writes the finder's search state with Dendrolib's `AEH_BHaHAHA::create_checkpoint` to `<BSSN_CHKPT_FILE_PREFIX>_aeh_solver_checkpt-cp<index>.json`, using the same 0/1 index as the solver checkpoint, and restore reads it with `AEH_BHaHAHA::restore_checkpoint`; if the file is missing, the finder keeps its initial state. A failed open at write time prints a message and leaves any earlier file at that index in place, and a horizon file that cannot be parsed makes restore throw. A restored run skips the initial output at the restored step. The finder writes a diagnostics row for each successfully found horizon only on horizon-find steps (multiples of `AEH_SOLVER_FREQ`) that are also multiples of `BSSN_IO_OUTPUT_FREQ`, and none when that frequency is 0.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::write_checkpt`, `Ctx::restore_checkpt`, and `Ctx::apparent_horizon_output` within `output_solver_context_cpp`; `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`; Dendrolib `src/aeh_bhahaha.cpp`, `AEH_BHaHAHA::create_checkpoint` and `AEH_BHaHAHA::restore_checkpoint` (file contents and missing-file return) and `AEH_BHaHAHA::find_horizons` (the `file_output_freq_` guard).
- Corroboration: Dendro-GR `BSSN_GR/src/bssnCtx.cpp`, `BSSNCtx::write_checkpt` and `BSSNCtx::restore_checkpt`, native file name and calls.

The generated binary-black-hole path parses native `BH_WAMR` mode 4,
selected refinement variables, constant or causal mode-6 wavelet tolerance,
wavelet coarsening factors, puncture-centered level floors, post-merger remesh
cadence, and optional wave-zone Nyquist refinement. The analytic initial-grid
seed comes from Dendro-GR's `punctureDataPhysicalCoord`, with its source and
MIT license embedded in the generated entry point. Its lapse and chi floors use
`std::fmax`, so at a sample exactly on a puncture, where the Dendro-GR
expressions give NaN, the lapse and chi take the floor value instead and can
still trigger refinement there. A separate single-rank
`--tpid` run computes TwoPunctures coefficients once. Fresh evolution loads
those coefficients before evolved-state conversion; after any initial-grid
remesh, the solver reconstructs that state from the same data. Grid
construction remains distinct from evolved initial data. Unsupported refinement or
tolerance modes fail at startup instead of silently substituting a constant
tolerance. These choices
target a comparable grid structure under the same parameter file; they do
not establish bitwise identity of the two remesh histories. The generated
Nyquist path and native `calculate_relative_position_history` both build the
relative-position history from all three components of the puncture separation.
The tracked puncture centers that feed that history, and the time integrator,
still differ between the two codes (see the comparison table in [BSSN
Application Wiring](bssn-application-wiring.md)), so enabling Nyquist
refinement can still produce different remesh decisions even with the same
parameters.

Claim evidence:
- Claim: The generated binary-puncture path parses native-style wavelet/geometric AMR controls, loads separately solved TwoPunctures data for evolution, and uses the licensed native analytic octree seed, whose lapse and chi floors return the floor value where the native expressions give NaN at a sample exactly on a puncture; its Nyquist history uses all three coordinate separations, as the native `calculate_relative_position_history` does.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`; `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::get_wtol_function` and `Ctx::is_remesh` within `output_solver_context_cpp`.
- Corroboration: `BSSN_GR/src/grUtils.cpp`, `punctureDataPhysicalCoord`; `BSSN_GR/src/dataUtils.cpp`, `calculate_relative_position_history` and `isRemeshBH`.

Initial refinement runs at most `BSSN_INIT_GRID_ITER - 1` passes (the full count
when that value is 1; q1's value of 10 gives nine). Each pass decides refinement
from the values transferred from the previous mesh, the first pass from the
freshly initialized data, and the loop stops when no refinement is due or a pass
leaves the global element and node counts unchanged. No lapse and chi floor or
algebraic projection runs between these transfers; after the loop, if any pass
remeshed, the solver rebuilds the initial data once on the final mesh. Evolution
remeshing instead applies the floors and projection right after each transfer.

Claim evidence:
- Claim: Initial refinement runs at most `BSSN_INIT_GRID_ITER - 1` passes (the full count when it is 1), each deciding refinement from the transferred values, and stops when no refinement is due or the global element and node counts are unchanged; no floor or projection runs between the initial transfers, the initial data are rebuilt once on the final mesh after any remesh, and evolution remeshing applies `post_timestep` right after each transfer.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp` (initial-grid loop and evolution loop).
- Corroboration: none available; no CI case checks the pass count or the floors between passes.

The generated solver moves each tracked puncture center once per step, after any
remesh of that step, by `-dt beta^i`, where `dt` is the time since the previous
update and `beta^i` is the shift interpolated from the evolved state at the
center's previous position. It uses no predictor and no stored velocity, and a
nonfinite interpolated shift stops the run. The centers set the excision regions
of the constraint norms, the puncture-centered refinement floors, the merged
test, and the puncture-center history that the Nyquist path uses.

Claim evidence:
- Claim: After each step and any remesh, the generated solver displaces each tracked puncture center once by `-dt beta^i`, with `dt` the time since the previous center update and `beta^i` interpolated from the evolved state at the previous center; there is no predictor or stored velocity, and a nonfinite shift stops the run. The centers feed the excision regions, the puncture-centered refinement floors, the merged test, and the center history.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::evolve_excision_centers` within `output_solver_context_cpp`; `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp` (the evolution loop order).
- Corroboration: `nrpy/examples/tests/dendro_application_check.py`, `Leg.run_variant`, the point reflection of the checkpointed centers and their motion along the puncture momenta; it does not check the displacement formula.

The generated refinement applies the following fixed constants and switching
rules in addition to the keys in [Runtime Parameter Keys](runtime-parameters.md).

- The punctures count as merged when the coordinate distance between the two
  tracked centers is below 0.1. The remesh cadence and the level-floor rule use
  the current separation. The coarsening factor switches to
  `BSSN_DENDRO_AMR_FAC_POST_MERGER` (when positive) only after a checkpoint has
  been written while the punctures are merged.
- Inside `BSSN_BH{1,2}_AMR_R` of a puncture the wavelet tolerance is multiplied
  by 1e12.
- Level floors apply to an element whose nearest corner lies within
  `orbital_radius = max(M1, M2) / (M1 + M2) * separation + 8` of the origin: the
  element requests at least level 9. Elements farther than `orbital_radius`
  from a puncture receive no floor from that puncture.
- While the punctures are separate, the floor around each puncture is
  `BSSN_BH{1,2}_MAX_LEV - 2` inside `BSSN_BH{1,2}_AMR_R`, and each coarser level down
  to 10 applies inside a radius larger by `BSSN_AMR_R_RATIO` than the one before
  (the loop stops when the level reaches 9).
- After the merger the floor is `ceil(2 + log2(25 Lx / (order r_lim))) - 2` inside
  `r_lim = max(BSSN_BH1_AMR_R, BSSN_BH2_AMR_R, 1.55 (M1 + M2))`, where `Lx` is the
  x-width of the domain and `order` the element order, and each coarser level
  down to 10 applies inside a radius twice as large (the loop stops when the level
  reaches 9). The puncture-radius loops never apply level 9; the origin rule above makes the level-9 request.
- Each floor above is a requested level, and the largest request is limited to
  `BSSN_MAXDEPTH - 2`. The origin floor is therefore `min(9, BSSN_MAXDEPTH - 2)`:
  level 9 for the packaged q1 depth of 14 and level 8 for `BSSN_MAXDEPTH = 10`.
- In wavelet mode 6 the radial tolerance is `BSSN_WAVELET_TOL` for `r <= 8`,
  interpolates logarithmically in `r` to `BSSN_GW_REFINE_WTOL` at the first
  extraction radius, equals `BSSN_GW_REFINE_WTOL` out to the last extraction
  radius, and is `BSSN_WAVELET_TOL_MAX` beyond it. Before the causal time
  `max(r, (r + 120) / sqrt(2))` the tolerance is `BSSN_WAVELET_TOL_MAX`, and over
  the next 100 time units it interpolates logarithmically to the radial value.
  The native keys `BSSN_WAVELET_TOL_FUNCTION_R0` and `BSSN_WAVELET_TOL_FUNCTION_R1`
  are not read.

Claim evidence:
- Claim: The generated solver treats the punctures as merged below a separation of 0.1; multiplies the wavelet tolerance by 1e12 inside `BSSN_BH{1,2}_AMR_R`; applies level floors with the constants 9, 8, 1.55, 25, 2, the levels 10 and above of the puncture loops, and the ratio `BSSN_AMR_R_RATIO` as written above, and limits every requested level to `BSSN_MAXDEPTH - 2`; uses the mode-6 radial and causal tolerance profile with the fixed radius 8 and the time constants 120 and 100; and does not read `BSSN_WAVELET_TOL_FUNCTION_R0` or `_R1`.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::is_remesh`, `Ctx::is_remesh_due`, and `Ctx::get_wtol_function` within `output_solver_context_cpp`; `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp` (the key reads).
- Corroboration: none available: no test drives the merged branch or the mode-6 profile; the constants follow Dendro-GR `BSSN_GR/src/dataUtils.cpp`, `isRemeshBH`.

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

Grid size is written initially and after each scheduled remesh check to
`dat/dgr_GridInfo.dat` with the default prefix. An actual remesh also updates
the effective VTU and gravitational-wave frequencies from the current maximum
mesh level. See [Grid size and native output cadence](constraints-and-diagnostic-norms.md#grid-size-and-native-output-cadence)
for the columns, count reduction, scaling, and restart behavior.

## Sources

- [Dendrolib block.h](https://github.com/paralab/Dendro-5.01/blob/master/include/block.h) - `ot::Block` geometry.
- [Dendrolib mesh.h](https://github.com/paralab/Dendro-5.01/blob/master/include/mesh.h) - zip, unzip, remesh, and intergrid transfer.
- [Dendrolib aeh_bhahaha.cpp](https://github.com/paralab/Dendro-5.01/blob/master/src/aeh_bhahaha.cpp) - apparent-horizon finder checkpoint write and restore.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - generated context and service scheduling.
- [checkpoint.py](../../../nrpy/infrastructures/Dendro/checkpoint.py) - checkpoint and restore generation.
- [physical_boundary.py](../../../nrpy/infrastructures/Dendro/general_relativity/physical_boundary.py) - boundary kernel registration.
- [BSSN_to_ADM.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN_to_ADM.py) - geometry conversion.
- [floor_the_lapse_and_conformal_factor.py](../../../nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py) - representation-dependent lapse and conformal-factor floors.
- [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/Dendro/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - owned-node algebraic projection.
- [physical_boundary_ghosts.py](../../../nrpy/infrastructures/Dendro/general_relativity/physical_boundary_ghosts.py) - physical exterior-padding extrapolation.
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - analytic seed, AMR parameters, and scheduling.
- [Dendro-GR dataUtils.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/dataUtils.cpp) - native black-hole refinement path.
- [Dendro-GR bssnCtx.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/bssnCtx.cpp) - native apparent-horizon checkpoint file name and restore.
- [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py) - shift-driver equation.
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - generated parameter validation and the q1 key names.
- [state_h.py](../../../nrpy/infrastructures/Dendro/state_h.py) - the `f_infinity` values registered for the evolved fields.
- [dendro_application_check.py](../../../nrpy/examples/tests/dendro_application_check.py) - CI checks of restart, reflection, and rejections.
- [Dendro-GR rhs.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/rhs.cpp) - native RIT `eta` selection.
- [Dendro-GR bssneqs_SSL_HD_dxsq.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/bssneqs_SSL_HD_dxsq.cpp) - native `B_rhs`.
- [Dendro-GR grUtils.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/grUtils.cpp) - native analytic puncture seed `punctureDataPhysicalCoord`.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Gridfunctions, Naming, And Loops](gridfunctions-naming-and-loops.md)
- See also: [Constraints And Diagnostic Norms](constraints-and-diagnostic-norms.md)
- Validated by: [Production Validation And Deferred Checks](production-validation-and-deferred-checks.md)
