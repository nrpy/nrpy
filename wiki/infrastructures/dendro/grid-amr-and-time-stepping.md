# Octree Grid, AMR, And Time Stepping

> Dendro octree blocks, ghost exchange, outgoing-radiation boundaries, checkpoint/restore checks, puncture tracking, and adaptive mesh refinement in NRPy's Dendro applications. · Status: provisional
> Up: [Dendro](index.md)

## Summary

Dendrolib represents adaptive meshes with balanced octants and supplies padded `ot::Block` arrays for finite-difference kernels. Generated `Ctx` schedules evolved-state ghost exchange, 6-component `RbarDD` scratch storage, RK4 substeps, interior-node floors and algebraic projection, checkpointing, apparent-horizon state persistence, and `BH_WAMR` wavelet/puncture-floor remeshing.

## Detail

### Block Traversal, Physical Boundaries, And Halo Exchange

Zipped vectors index continuous-Galerkin octree nodes, whereas unzipped vectors hold padded regular blocks. Before each Ricci/RHS pass (`Ctx::rhs` or block-list `Ctx::rhs_blkwise`), `Ctx` exchanges ghosts and unzips the evolved state once, then visits each local block once:
1. `physical_boundary_ghosts` extrapolates exterior physical padding successively across faces, edges, and corners using `5` interior nodes at FD4 or `6` at FD6/FD8 (giving 4th-order boundary second derivatives at FD6/FD8).
2. `Ricci_eval_order_N` populates local `RbarDD` scratch and `rhs_eval_order_N` consumes it immediately without inter-block exchange.
3. `physical_boundary` overwrites evolved RHSs on physical domain-face nodes only (using centered derivatives `d_i f` over extrapolated padding and no field-dependent wave speed):

```text
d_t f = -(x^i d_i f + n_f (f - f_inf)) / r,    r = sqrt(x^2 + y^2 + z^2),
```

with falloff `n_f = 2` for `lambdaU`/`aDD` (`n_f = 1` otherwise) and `f_inf = 1` for `alpha`/`cf` (`0` otherwise).

Because `Ctx::zip()` populates only rank-local nodes without ghost values, any newly zipped diagnostic array read by an element interpolator (such as both `Psi4` components during wave extraction) completes `readFromGhostBegin/End` before interpolation. Floors (`alpha >= CHI_FLOOR`; `W >= sqrt(CHI_FLOOR)` or `chi >= CHI_FLOOR`) and `enforce_detgbar_equals_detghat_trAzero` run directly on zipped interior nodes after initial-data conversion, each RK4 stage, and AMR transfer; after the final RK4 projection and any remesh, `Ctx` exchanges evolved-state ghosts before puncture tracking and diagnostics.

Evidence (descriptive behavior): `Dendro/general_relativity/physical_boundary_ghosts.py`; `Dendro/general_relativity/physical_boundary.py`; `Dendro/state_h.py`; `Dendro/solver_context.py::output_solver_context_cpp`.

### Checkpoints, Restore Validation, And Apparent-Horizon State

New RK4 timesteps equal `BSSN_CFL_FACTOR * min(dx, dy, dz)` at the global finest octree level (`dx_min`), recomputed after initial convergence and evolution remeshing; nonfinite registered floating parameters fail before puncture-file reads, and domain widths plus `BSSN_CFL_FACTOR` must be finite and positive.

Checkpoints write after step output at nonzero multiples of `BSSN_CHECKPT_FREQ` into alternating slot `(step / BSSN_CHECKPT_FREQ) mod 2` (publishing octree and state before metadata, then horizon state; no final write occurs if evolution stops between multiples). On `BSSN_RESTORE_SOLVER = 1`, rank `0` selects the metadata file with the newer POSIX modification time and makes no fallback to the other slot (unlike native `BSSN_GR`, which honors `BSSN_RESTORE_CHECKPT_SLOT`, `.latest`, step comparison, incomplete-slot fallback, and merger slot `3`). If neither slot 0 nor slot 1 metadata exists under `BSSN_CHKPT_FILE_PREFIX`, rank `0` logs a warning and starts fresh from TwoPunctures initial data; if metadata exists, restore aborts (`Checkpoint metadata does not match <formulation>` / `checkpoint restore failed`) unless:
- Formulation tag, field count/names, and parameter-name list match the executable.
- Every emitted `CodeParameter` (`C_CAHD`, both KO strengths via `KO_DISS_SIGMA`, `SSL_h`, `SSL_sigma`, `chi_floor`, `eta`, optional `CFL_FACTOR`/`YBS_chi`/`C_YBS_mom`, plus fCCZ4 `kappa1`/`kappa2`), domain bounds, and `BSSN_ELE_ORDER` (`4`, `6`, or `8`) match stored values exactly (non-`CodeParameter` mesh/output keys may change).
- Stored time, positive `dt <= CFL bound`, strictly increasing puncture history ending at stored time/centers, and both stored and post-restore global max algebraic residuals `max(|det(gammabar) - 1|, |gammabar^ij Abar_ij|)` are finite and `<= 1e-10` (with positive finite `det(gammabar)`), and MPI rank count is `>=` stored active communicator size.

When `AEH_SOLVER_FREQ > 0`, checkpoints also write `<BSSN_CHKPT_FILE_PREFIX>_aeh_solver_checkpt-cp<0|1>.json` (horizon count, BBH flag, last three shapes/centers/times/radii, active/failure/fixed-radius flags). Because Dendrolib's `AEH_BHaHAHA::create_checkpoint` returns `void` after printing `file open failed for BAH checkpoint`, an unchecked open failure leaves the two-checkpoints-older JSON file in place for a later restore; missing JSON on restore logs `file open failed! Could not restore AH solver!` and leaves initial finder state, while malformed JSON throws. Restored runs skip step-0 output; `BHaHAHA_diagnostics` rows write only on steps dividing both `AEH_SOLVER_FREQ` and nonzero `BSSN_IO_OUTPUT_FREQ`.

Evidence (descriptive behavior): `Dendro/checkpoint.py::output_checkpoint_cpp`; `Dendro/main_cpp.py::output_main_cpp`; `Dendro/solver_context.py`; [Dendrolib aeh_bhahaha.cpp](https://github.com/paralab/Dendro-5.01/blob/master/src/aeh_bhahaha.cpp).

### Initial Mesh Refinement, Puncture Tracking, And Fixed AMR Rules

Initial octree construction evaluates Dendro-GR's analytic `punctureDataPhysicalCoord` (`std::fmax` floors prevent `NaN` on exact puncture coordinates), loads single-rank `--tpid` spectral data, and runs at most `BSSN_INIT_GRID_ITER - 1` refinement passes (`1` when `BSSN_INIT_GRID_ITER == 1`; `9` for `q1`'s `10`), stopping early when element and node counts stabilize. No floor or projection runs between initial passes; if any pass remeshed, TwoPunctures initial data are rebuilt once on the final mesh. Unsupported refinement or wavelet-tolerance modes abort at startup.

After each step and any remesh, `Ctx::evolve_excision_centers` advances each puncture by `-dt beta^i(x_prev)` using evolved `vetU` (aborting on nonfinite shift, with no predictor or stored velocity; differences from native Heun tracking and RK schemes mean wave-zone Nyquist histories can diverge from `BSSN_GR`).

Beyond [Runtime Parameter Keys](runtime-parameters.md), `BH_WAMR` mode `4` enforces fixed switching rules:
- **Merger threshold**: coordinate separation `< 0.1` switches remesh cadence and level floors; coarsening switches to `BSSN_DENDRO_AMR_FAC_POST_MERGER > 0` only after a post-merger checkpoint is written.
- **Puncture wavelet relaxation & origin floor**: wavelet tolerance multiplies by `1e12` within `BSSN_BH{1,2}_AMR_R`; elements whose nearest corner lies within `orbital_radius = max(M1, M2)/(M1 + M2) * separation + 8` of the origin request at least `min(9, BSSN_MAXDEPTH - 2)` (elements farther than `orbital_radius` from a puncture receive no puncture floor).
- **Pre- and post-merger radial level floors** (all capped at `BSSN_MAXDEPTH - 2`): before merger, `BSSN_BH{1,2}_MAX_LEV - 2` inside `BSSN_BH{1,2}_AMR_R`, stepping down to level `10` as radius multiplies by `BSSN_AMR_R_RATIO`; after merger, `ceil(2 + log2(25 Lx / (order * r_lim))) - 2` inside `r_lim = max(BSSN_BH1_AMR_R, BSSN_BH2_AMR_R, 1.55 (M1 + M2))`, stepping down to level `10` per factor of `2` in radius.
- **Causal wavelet mode `6`** (`BSSN_WAVELET_TOL_FUNCTION_R{0,1}` unread): radial target equals `BSSN_WAVELET_TOL` for `r <= 8`, log-interpolates to `BSSN_GW_REFINE_WTOL` at the first extraction radius, stays `BSSN_GW_REFINE_WTOL` to the last extraction radius, and is `BSSN_WAVELET_TOL_MAX` beyond; before `t_causal = max(r, (r + 120) / sqrt(2))` tolerance is `BSSN_WAVELET_TOL_MAX`, log-interpolating to the radial target over `100` time units.

Evidence (descriptive behavior): `Dendro/solver_context.py::output_solver_context_cpp`; `Dendro/main_cpp.py`; [Dendro-GR dataUtils.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/dataUtils.cpp).

### Shift-Driver Wavelet Scaling And Remesh Collectives

Because native `BSSN_GR` evolves `d_t beta^i = (3/4) B^i` whereas NRPy's `GammaDriving2ndOrder_Covariant__Hatted` evolves `d_t beta^i = B^i + advection` (`B_native = (4/3) betU` wherever `eta` agrees), `Ctx::is_remesh` multiplies `betU0..2` by `4/3` inside the temporary unzipped work vector after physical-boundary fill and before `isReMeshUnzip`, preventing premature coarsening without mutating evolved state.

When initial convergence shrinks the active communicator, only block-local geometry work is gated on `mesh->isActive()`; every global MPI rank still enters `isReMeshUnzip` collectives so communicator re-expansion cannot deadlock. Grid counts write to `dat/dgr_GridInfo.dat` initially and after each remesh test, and actual remeshes rescale VTU/GW cadences ([Grid size and native output cadence](constraints-and-diagnostic-norms.md#grid-size-and-native-output-cadence)).

Evidence (descriptive behavior): `Dendro/solver_context.py::output_solver_context_cpp`; `BSSN_gauge_RHSs.py`; [Dendro-GR bssneqs_SSL_HD_dxsq.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/bssneqs_SSL_HD_dxsq.cpp).

## Sources

- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - `Ctx` block traversal, remeshing, puncture tracking, checkpointing
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - `<stem>_params_validate`, `runtime_parameter_names`
- [checkpoint.py](../../../nrpy/infrastructures/Dendro/checkpoint.py) - `output_checkpoint_cpp`
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - `output_main_cpp`, initial-grid loop, slot selection, CFL step
- [state_h.py](../../../nrpy/infrastructures/Dendro/state_h.py) - `register_canonical_gridfunctions` (`f_infinity`)
- [physical_boundary_ghosts.py](../../../nrpy/infrastructures/Dendro/general_relativity/physical_boundary_ghosts.py) - `register_CFunction_physical_boundary_ghosts`
- [physical_boundary.py](../../../nrpy/infrastructures/Dendro/general_relativity/physical_boundary.py) - `register_CFunction_physical_boundary`
- [BSSN_to_ADM.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN_to_ADM.py) - `register_CFunction_BSSN_to_ADM`
- [floor_the_lapse_and_conformal_factor.py](../../../nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py) - `CHI_FLOOR` kernel
- [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/Dendro/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - interior-node algebraic projection
- [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py) - `GammaDriving2ndOrder_Covariant__Hatted`
- [dendro_application_check.py](../../../nrpy/examples/tests/dendro_application_check.py) - `Leg.run_variant`, `Leg.run_negatives`
- [Dendrolib block.h](https://github.com/paralab/Dendro-5.01/blob/master/include/block.h) - `ot::Block`
- [Dendrolib mesh.h](https://github.com/paralab/Dendro-5.01/blob/master/include/mesh.h) - `zip`, `unzip`, `isReMeshUnzip`
- [Dendrolib aeh_bhahaha.cpp](https://github.com/paralab/Dendro-5.01/blob/master/src/aeh_bhahaha.cpp) - `AEH_BHaHAHA` checkpointing and `find_horizons`
- [Dendro-GR dataUtils.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/dataUtils.cpp) - `isRemeshBH`, `calculate_relative_position_history`
- [Dendro-GR bssnCtx.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/bssnCtx.cpp) - `BSSNCtx::write_checkpt`, `BSSNCtx::restore_checkpt`
- [Dendro-GR rhs.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/rhs.cpp) - RIT `eta` profile
- [Dendro-GR bssneqs_SSL_HD_dxsq.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/bssneqs_SSL_HD_dxsq.cpp) - native `B_rhs`
- [Dendro-GR grUtils.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/grUtils.cpp) - `punctureDataPhysicalCoord`

## See Also

- Depends on: [Gridfunctions, Naming, And Loops](gridfunctions-naming-and-loops.md)
- See also: [BSSN Application Wiring](bssn-application-wiring.md)
- Validated by: [Production Validation And Deferred Checks](production-validation-and-deferred-checks.md)
