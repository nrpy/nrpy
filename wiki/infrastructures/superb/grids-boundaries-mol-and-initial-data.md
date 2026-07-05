# Grids, Boundaries, MoL, And Initial Data

> Chare-local grid setup, boundary exchange, Method of Lines phases, initial-data staging, and NRPyElliptic integration hooks in superB. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [superB](index.md)

## Summary

superB decomposes each numerical grid across the `Nchare0`, `Nchare1`, and
`Nchare2` chare-array dimensions. Each `Timestepping` chare owns a rectangular
local grid with ghost zones, local coordinate arrays, local boundary metadata,
and communication buffers. The generated runtime synchronizes ordinary ghost
faces between neighboring chares, then handles nonlocal curvilinear inner
boundary sources separately through source-point messages.

## Detail

`numerical_grid_params_Nxx_dxx_xx_chare` derives local grid parameters from
global grid parameters and the current `chare_index`. It requires each global
`Nxx*` dimension to divide evenly by the matching `Nchare*` count, rejects
multi-chare subdomains smaller than `NGHOSTS`, copies the global spacings and
inverse spacings, computes local `xxmin*` and `xxmax*`, and allocates
cell-centered `xx[0..2]` arrays over `Nxx_plus_2NGHOSTS*` local extents.

`numerical_grids_chare` is the chare-local assembly point. For every grid, it
copies grid identity fields, calls the local-grid setup function, optionally
allocates reference-metric precompute storage, builds the chare communication
map, and, when CurviBCs are enabled, builds the chare-local boundary-condition
structs. It also initializes the per-grid diagnostic struct and invokes the
setup modes for 1D and 2D diagnostics; the diagnostics dispatcher behavior is
covered by the diagnostics leaf, not here.

Index ownership is split between a local map and header-level arithmetic.
`charecommstruct_set_up` allocates `localidx3pt_to_globalidx3pt` for every local
point including ghost zones. The global-owner and global-to-local conversions
are computed on demand through `IDX3_OF_CHARE`, `MAP_LOCAL_TO_GLOBAL_IDX*`,
`MAP_GLOBAL_TO_LOCAL_IDX*`, `globalidx3pt_to_chareidx3`, and
`globalidx3pt_to_localidx3pt` in `superB.h`.

CurviBC setup starts from the global `bc_struct` and filters it into
chare-local storage. `bcstruct_chare_set_up` counts inner boundary points whose
destination lies on the current chare, separates points whose source is local
from points whose source belongs to another chare, converts local entries from
global to local indices, records parity, and builds `nonlocalinnerbc_struct`
arrays for both incoming source chares and outgoing destination chares. The same
function filters pure outer boundary points into local `pure_outer_bc_array`
entries. `apply_bcs_inner_only_nonlocal` later consumes received source-point
buffers and applies parity-aware inner BC updates to the nonlocal entries.

Neighbor ghost exchange is face based. `send_neighbor_data` packs `NGHOSTS`
interior faces for the east-west, north-south, or top-bottom direction into
temporary face buffers and sends them to the adjacent chare entry method.
`process_ghost` writes the received values into the matching ghost slab. The
generated timestepping flow performs these exchanges for `Y_N_GFS`, for
synchronized auxiliary-evolution gridfunctions when present, and after each RK
substep for the gridfunction sets selected by the MoL generator.

Nonlocal inner BC exchange is point based. During setup, each chare sends the
global source-point indices it needs to the source-owning chares, and source
chares cache the points they must later send. During evolution,
`send_nonlocalinnerbc_data` gathers those source values by global index,
dispatches type-specific `receiv_nonlocalinnerbc_data_*` entry methods, and
`process_nonlocalinnerbc` calls `apply_bcs_inner_only_nonlocal` on the received
buffers.

The superB MoL generator splits each RK substep into explicit phases:
`MOL_PRE_RK_UPDATE` for RHS evaluation, `MOL_RK_UPDATE` for the RK update,
`MOL_POST_RK_UPDATE_APPLY_BCS` when extrapolation-style outer BCs are requested,
and `MOL_POST_RK_UPDATE` for post-RHS work. `generate_rhs_output_exprs` and
`generate_post_rhs_output_list` choose which MoL storage arrays need
synchronization for each substep. `MoL_sync_data_defines` builds the runtime
sync lists: evolved and auxiliary-evolution entries are selected from
`BHaHGridFunction` objects whose `sync_gf_in_superB` metadata is true, while
the auxiliary sync list is currently limited to the optional Psi4 diagnostic
pair, `DIAG_PSI4_RE` and `DIAG_PSI4_IM`, when `enable_psi4` is true.

Initial data has two runtime flows. Standard non-NRPyElliptic superB runs call
`initial_data` in `INITIALDATA_BIN_ONE`, apply inner-only BCs, synchronize
nonlocal inner BC data for the first `y_n` initial-data stage and auxevol data,
call `INITIALDATA_BIN_TWO`, apply outer-extrapolated plus inner BCs, then repeat
the needed nonlocal syncs. The initial-data function itself registers exact ADM
initial data when available, registers ADM-to-BSSN converters for the requested
coordinate systems, optionally reads checkpoint data, performs the two converter
bins, and dispatches the named BC application stages.

The generated `post_non_y_n_auxevol_mallocs` hook runs after memory for
non-`y_n` and auxiliary-evolution gridfunctions has been allocated. In the
standard flow it runs after `INITIALDATA_BIN_TWO` and its outer-extrapolated
plus inner BC application; in the NRPyElliptic flow it runs after the single
`initial_data` call and before `send_wavespeed_at_outer_boundary(grid)`
broadcasts the selected outer-boundary wavespeed. Current examples use this
hook for CAHD prefactor setup and for setting NRPyElliptic auxiliary-evolution
gridfunctions to constants.

NRPyElliptic projects use a narrower superB integration path. The example
registers NRPyElliptic initial-guess, auxevol, RHS, residual, and stop-condition
functions, but superB only wires their runtime hooks: `nrpyelliptic_project=True`
selects a single initial-data call, registers a per-grid
`wavespeed_at_outer_boundary` parameter, and emits
`send_wavespeed_at_outer_boundary` so the chare owning a selected outer-boundary
point broadcasts that wavespeed to every chare before evolution. The example's
RHS string then uses that value in radiation BC wavespeed inputs. The generated
volume-integration report updates `commondata.log10_current_residual` from the
`DIAG_RESIDUALGF` RMS on `sphere_R_80`, and the example passes a
`post_MoL_step_forward_in_time` hook that calls `stop_conditions_check` and
ends through `mainProxy.done()` when `commondata.stop_relaxation` is set.

## Sources

- [numerical_grids.py](../../../nrpy/infrastructures/superB/numerical_grids.py) - `register_CFunction_numerical_grid_params_Nxx_dxx_xx_chare`, `register_CFunction_numerical_grids_chare`, `register_CFunctions`
- [chare_communication_maps.py](../../../nrpy/infrastructures/superB/chare_communication_maps.py) - `register_CFunction_charecommstruct_set_up`, `chare_comm_register_C_functions`
- [CurviBoundaryConditions.py](../../../nrpy/infrastructures/superB/CurviBoundaryConditions.py) - `register_CFunction_apply_bcs_inner_only_nonlocal`, `register_CFunction_bcstruct_chare_set_up`, `CurviBoundaryConditions_register_C_functions`
- [MoL.py](../../../nrpy/infrastructures/superB/MoL.py) - `register_CFunctions`, `register_CFunction_MoL_step_forward_in_time`, `generate_rhs_output_exprs`, `generate_post_rhs_output_list`, `register_CFunction_MoL_sync_data_defines`
- [initial_data.py](../../../nrpy/infrastructures/superB/initial_data.py) - `register_CFunction_initial_data`, `register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN`
- [timestepping_chare.py](../../../nrpy/infrastructures/superB/timestepping_chare.py) - `generate_mol_step_forward_code`, `generate_send_neighbor_data_code`, `generate_process_ghost_code`, `generate_send_nonlocalinnerbc_data_code`, `generate_process_nonlocalinnerbc_code`, `output_timestepping_h_cpp_ci_register_CFunctions`
- [superB_blackhole_spectroscopy.py](../../../nrpy/examples/superB_blackhole_spectroscopy.py) - `post_non_y_n_auxevol_mallocs`, `cahdprefactor_auxevol_gridfunction`
- [superB_nrpyelliptic_conformally_flat.py](../../../nrpy/examples/superB_nrpyelliptic_conformally_flat.py) - `post_MoL_step_forward_in_time`, `rhs_string`, `register_CFunction_residual_H_compute_all_points`, `register_CFunction_stop_conditions_check`
- [superB.h](../../../nrpy/infrastructures/superB/superB/superB.h) - `IDX3_OF_CHARE`, `MAP_LOCAL_TO_GLOBAL_IDX0`, `MAP_GLOBAL_TO_LOCAL_IDX0`, `MOL_PRE_RK_UPDATE`, `INITIALDATA_BIN_ONE`, `nonlocalinnerbc_struct`

## See Also

- [superB](index.md)
- [Chare Entrypoints And Runtime](chare-entrypoints-and-runtime.md)
- [Diagnostics And Observables](diagnostics-and-observables.md)
- [BHaH Lifecycle](../bhah-lifecycle.md)
- [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
