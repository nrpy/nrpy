# Grids, Boundaries, MoL, And Initial Data

> Chare-local grid setup, boundary exchange, Method of Lines phases, initial-data staging, and NRPyElliptic integration hooks in superB. · Status: confirmed
> Up: [superB](index.md)

## Summary

superB decomposes each numerical grid across the `Nchare0`, `Nchare1`, and
`Nchare2` chare-array dimensions. Each `Timestepping` chare holds a rectangular
local grid with ghost zones, local coordinate arrays, local boundary metadata,
and communication buffers. The generated runtime synchronizes ordinary ghost
faces between neighboring chares, then handles nonlocal curvilinear inner
boundary sources separately through source-point messages.

## Detail

### Chare Grid Decomposition, Ghost Exchange, And Nonlocal CurviBCs

`numerical_grid_params_Nxx_dxx_xx_chare` derives local grid parameters from
global grid parameters and the current `chare_index`. It requires each global
`Nxx*` dimension to divide evenly by the matching `Nchare*` count, rejects
multi-chare subdomains smaller than `NGHOSTS`, copies the global spacings and
inverse spacings, computes local `xxmin*` and `xxmax*`, and allocates
cell-centered `xx[0..2]` arrays over `Nxx_plus_2NGHOSTS*` local extents.

`numerical_grids_chare` is the chare-local assembly point. For every grid, it
copies grid identity fields, calls the local-grid setup function, optionally
allocates reference-metric precompute storage, and, when CurviBCs are enabled,
builds the chare-local boundary-condition structs. It also initializes the
per-grid diagnostic struct and invokes the setup modes for 1D and 2D
diagnostics; the diagnostics dispatcher behavior is covered by the diagnostics
leaf, not here.

Subdomain ownership is arithmetic rather than a stored chare communication map.
`superB.h` provides `IDX3_OF_CHARE`, `MAP_LOCAL_TO_GLOBAL_IDX*`,
`MAP_GLOBAL_TO_LOCAL_IDX*`, `globalidx3pt_to_chareidx3`, and
`globalidx3pt_to_localidx3pt`. Nonlocal inner-boundary setup and send code maps
global source and destination indices to target chares and local indices on
demand.

For `GeneralRFM_inverse_fisheyeN*`, the map is injective with
Cartesian topology, so there are no inner coordinate boundaries.
`bcstruct_chare_set_up` directly intersects each global outer face with the
chare's ownership box, preserving ghost-layer order and face directions.
It needs neither Cartesian inversion nor a global boundary array. The
8M-separation inverse-spheroidal example disables global boundary construction,
initializes its unused global boundary struct to zero for cleanup/PUP, and
disables the global CFL scan. With `enable_distributed_cfl=True`, each
`Timestepping` chare scans its local grid, including ghosts, at the start of
`start()`. A Charm++ minimum reduction supplies the same sampled timestep to
all chares before horizon/interpolation services or evolution begin. Internal
ghost overlap repeats samples but does not change the minimum; outer ghosts
remain included. The example uses the ordinary `cfl_limited_timestep` scan,
not the optional conservative bound. Other examples retain their existing
CFL setup and global-list boundary filtering path.

The inverse-spheroidal spectroscopy example enables
`enable_xx_to_Cart_precompute=True` for grid allocation, initial data,
radiation boundaries, timestepping cleanup, and PUP registration.
`numerical_grids.py` allocates `griddata.xx_to_Cart[3]` over each chare's full
local grid, including ghosts. Initial-data setup fills these arrays while
precomputing the reference metric. Radiation boundaries read the stored
Cartesian positions and compute geometry outside the evolved-gridfunction loop.
`superB_pup.py` serializes the arrays and allocates them on unpacking;
`timestepping_chare.py` frees them during cleanup. Coordinate arrays are
separate from AUXEVOL gridfunctions and inter-chare gridfunction exchanges.

Evidence (descriptive behavior): `numerical_grids.py::register_CFunction_numerical_grids_chare`;
`initial_data.py::register_CFunction_initial_data`;
`superB_pup.py::register_CFunction_superB_pup_routines`;
`timestepping_chare.py::output_timestepping_cpp`.

The grid transport paths are ordinary Charm++ entry-method sends generated for
the `Timestepping` array:

| Flow | Routing | Message buffer | Source-code details |
| --- | --- | --- | --- |
| ghost face slabs | indexed point-to-point `Timestepping` sends | parameter-marshaled `REAL` face slabs | face area and `NGHOSTS` from superB code |
| nonlocal setup | sparse point-to-point sends | source index arrays | destination tells source which points are needed |
| nonlocal data | sparse point-to-point sends | parameter-marshaled `REAL` point data | phase-specific receive methods protect SDAG matching |
| MoL/BC ordering | SDAG local order plus `when` waits | messages buffer until wait reached | no global Charm barrier implied |

`Nchare0`, `Nchare1`, and `Nchare2` define grid and chare decomposition, not
processing-element placement. The chare/global/local index-mapping macros above
belong to grid indexing and sparse source lookup; load-balancing policy and PE
placement belong to [Chare Entrypoints And Runtime](chare-entrypoints-and-runtime.md). Current ghost and nonlocal boundary paths
do not use TRAM, multicast sections, priority sends, immediate entry methods,
expedited handling, or zero-copy transfer. Immediate entry methods are not a
valid replacement for these generated `Timestepping` chare-array receives: the
current `.ci` signatures are ordinary array entry methods carrying marshaled
buffers and matched by SDAG `when` clauses.

CurviBC setup starts from the global `bc_struct` eigencoordinate and parity maps
defined in [BHaH Curvilinear Boundaries And Stencils](../bhah/curvilinear-boundaries-and-stencils.md)
and partitions them into chare-local storage: `bcstruct_chare_set_up` separates
local-source inner boundary points from remote-source points (`nonlocalinnerbc_struct`
incoming/outgoing lists), filters local `pure_outer_bc_array` entries, and pairs
with `apply_bcs_inner_only_nonlocal` for parity-signed remote updates.

Neighbor ghost exchange is face based. `send_neighbor_data` packs `NGHOSTS`
interior faces for the east-west, north-south, or top-bottom direction into
temporary face buffers and sends them to the adjacent chare entry method.
`process_ghost` writes the received values into the matching ghost slab. The
generated timestepping flow performs these exchanges for `Y_N_GFS`, for
synchronized auxiliary-evolution gridfunctions when present, and after each RK
substep for the gridfunction sets selected by the MoL generator.

Nonlocal inner BC exchange is point based. During setup, each chare sends the
global source-point indices it needs to the source chares, and source
chares cache the points they must later send. During evolution,
`send_nonlocalinnerbc_data` gathers those source values by global index,
dispatches type-specific `receiv_nonlocalinnerbc_data_*` entry methods, and
`process_nonlocalinnerbc` calls `apply_bcs_inner_only_nonlocal` on the received
buffers.

### Method Of Lines Phases, Initial Data, And NRPyElliptic Hooks

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
bins, and dispatches the named BC application stages. The checked collision
example opts into the
same combined conformal determinant/trace projection after checkpoint recovery
and after final outer-plus-inner BC handling; its MoL post-RHS hook applies it
again after each RK update. `enable_conformal_projection` defaults to `False`,
so this projection schedule runs only for callers that opt in and supply the Method of
Lines hook. The checked collision example uses BSSN equation modules and does not wire
fCCZ4 RHSs, gauge RHSs, or `Theta_fCCZ4`.

Evidence (scientific/interface requirement): `superB/initial_data.py::register_CFunction_initial_data`; `superB_two_blackholes_collide.py`; corroboration: `superB/MoL.py::register_CFunctions`.

Both checked superB black-hole examples explicitly set
`enable_YBS_Gamma_constraint_adjustment = False` and forward it to the shared
BHaH `register_CFunction_rhs_eval`. superB maintains no separate YBS RHS fork:
changing that source constant routes through the BHaH registrar to the BSSN
RHS and gauge equation modules, while the enabled-only runtime `YBS_chi` parameter uses
the reused BHaH CodeParameters flow. Neither example exposes a command-line
switch for the code-generation choice.

Evidence (descriptive behavior): `superB_two_blackholes_collide.py`; `superB_blackhole_spectroscopy.py`; corroboration: `BHaH/general_relativity/rhs_eval.py::register_CFunction_rhs_eval`.

The generated `post_non_y_n_auxevol_mallocs` hook runs after memory for
non-`y_n` and auxiliary-evolution gridfunctions has been allocated. In the
standard flow it runs after `INITIALDATA_BIN_TWO` and its outer-extrapolated
plus inner BC application; in the NRPyElliptic flow it runs after the single
`initial_data` call and before `send_wavespeed_at_outer_boundary(grid)`
broadcasts the selected outer-boundary wavespeed. Current examples use this
hook for shared raw `dsmin` setup used by CAHD and for setting NRPyElliptic auxiliary-evolution
gridfunctions to constants.

NRPyElliptic projects use a narrower superB integration path. The example
registers NRPyElliptic initial-guess, auxevol, RHS, residual, and stop-condition
functions, but superB only wires their runtime hooks: `nrpyelliptic_project=True`
selects a single initial-data call, registers a per-grid
`wavespeed_at_outer_boundary` parameter, and emits
`send_wavespeed_at_outer_boundary` so the chare containing a selected outer-boundary
point broadcasts that wavespeed to every chare before evolution. The example's
RHS string then uses that value in radiation BC wavespeed inputs. The generated
volume-integration report updates `commondata.log10_current_residual` from the
`DIAG_RESIDUALGF` RMS on `sphere_R_80`, and the example passes a
`post_MoL_step_forward_in_time` hook that calls `stop_conditions_check` and
ends through `mainProxy.done()` when `commondata.stop_relaxation` is set.

Current configured CI builds the elliptic and spectroscopy variants without
running them, then runs only the collision variant. Source inspection establishes
emitted control flow, but configured automation does not establish elliptic
residual-stop behavior, spectroscopy/Psi4 synchronization, checkpoint/restart,
direct message-buffer checks, or distributed runtime outcomes.

## Sources

- [numerical_grids.py](../../../nrpy/infrastructures/superB/numerical_grids.py) - `register_CFunction_numerical_grid_params_Nxx_dxx_xx_chare`, `register_CFunction_numerical_grids_chare`, `register_CFunctions`
- [CurviBoundaryConditions.py](../../../nrpy/infrastructures/superB/CurviBoundaryConditions.py) - `register_CFunction_apply_bcs_inner_only_nonlocal`, `register_CFunction_bcstruct_chare_set_up`, `CurviBoundaryConditions_register_C_functions`
- [MoL.py](../../../nrpy/infrastructures/superB/MoL.py) - `register_CFunctions`, `register_CFunction_MoL_step_forward_in_time`, `generate_rhs_output_exprs`, `generate_post_rhs_output_list`, `register_CFunction_MoL_sync_data_defines`
- [initial_data.py](../../../nrpy/infrastructures/superB/initial_data.py) - `register_CFunction_initial_data`, `register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN`
- [timestepping_chare.py](../../../nrpy/infrastructures/superB/timestepping_chare.py) - `generate_mol_step_forward_code`, `generate_send_neighbor_data_code`, `generate_process_ghost_code`, `generate_send_nonlocalinnerbc_data_code`, `generate_process_nonlocalinnerbc_code`, `output_timestepping_h_cpp_ci_register_CFunctions`
- [superB_blackhole_spectroscopy.py](../../../nrpy/examples/superB_blackhole_spectroscopy.py) - `post_non_y_n_auxevol_mallocs`, `dsmin_auxevol_gridfunction`
- [superB_nrpyelliptic_conformally_flat.py](../../../nrpy/examples/superB_nrpyelliptic_conformally_flat.py) - `post_MoL_step_forward_in_time`, `rhs_string`, `register_CFunction_residual_H_compute_all_points`, `register_CFunction_stop_conditions_check`
- [superB_two_blackholes_collide.py](../../../nrpy/examples/superB_two_blackholes_collide.py) - opt-in conformal projection in initial-data and Method of Lines registration
- [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py) - shared BHaH BSSN/fCCZ4 and YBS option boundary used by superB
- [superB.h](../../../nrpy/infrastructures/superB/superB/superB.h) - `IDX3_OF_CHARE`, `MAP_LOCAL_TO_GLOBAL_IDX0`, `MAP_GLOBAL_TO_LOCAL_IDX0`, `globalidx3pt_to_chareidx3`, `globalidx3pt_to_localidx3pt`, `MOL_PRE_RK_UPDATE`, `INITIALDATA_BIN_ONE`, `nonlocalinnerbc_struct`
- [main.yml](../../../.github/workflows/main.yml) - `charmpp-validation`

## See Also

- Depends on: [BHaH Curvilinear Boundaries And Stencils](../bhah/curvilinear-boundaries-and-stencils.md)
- See also: [Chare Entrypoints And Runtime](chare-entrypoints-and-runtime.md)
- See also: [Diagnostics And Observables](diagnostics-and-observables.md)
