# Diagnostics Output And Checkpointing

> Explain BHaH diagnostics scheduling, temporary diagnostic buffers, raytracing export, progress output, and checkpoint/restart files. Status: confirmed. Last reconciled: 2026-06-29
> Up: [BHaH](index.md)

## Summary

BHaH diagnostics are generated as a scheduled `diagnostics()` driver plus
family-specific helper CFunctions. On each output step the driver allocates
temporary per-grid `diagnostic_gfs`, refreshes host data when CUDA output needs
host-side I/O, populates diagnostic channels, calls enabled output families, and
frees the temporary buffers. The same driver advances the progress indicator on
every timestep, independent of whether diagnostics were emitted.

Checkpointing is a separate generated pair, `read_checkpoint()` and
`write_checkpoint()`. Checkpoints are convergence-factor named binary files that
store `commondata`, optional BHaHAHA horizon history, per-grid metadata, and a
compact selected-point payload for all evolved gridfunctions. Restart validates
grid shape and coordinate-system hash before rebuilding `y_n_gfs`.

## Detail

`register_all_diagnostics` stages helper headers, validates optional
raytracing-data export, registers the top-level `_register_CFunction_diagnostics`
driver, and then registers per-coordinate nearest-diagnostic and volume-element
helpers. Raytracing export is deliberately narrow: it rejects CUDA generation,
requires `enable_rfm_precompute=True`, requires
`enable_RbarDD_gridfunctions=True`, allows exactly one coordinate system, and
currently supports only `Cartesian` and `Spherical`; the generated exporter also
aborts unless `commondata->NUMGRIDS == 1`.

`diagnostics()` is called once per timestep. Its cadence test is time-based:
`fabs(round(time / diagnostics_output_every) * diagnostics_output_every - time)
< 0.5 * dt`. The `diagnostics_output_every` commondata parameter must be
positive. On output steps the driver allocates
`REAL *diagnostic_gfs[MAXNUMGRIDS]`, one buffer per active grid with
`TOTAL_NUM_DIAG_GFS * Nxx_plus_2NGHOSTS_tot` values. Disabled MoL scratch
free/restore hooks are present only as commented generated C. On CUDA builds the
driver signature includes `griddata_device`; inside `#ifdef __CUDACC__` it
copies all `NUM_EVOL_GFS` current-time data from device to host, optionally
copies `T4UU` auxiliary data into `DIAG_T4UU` channels, and synchronizes the CUDA
stream before host-side diagnostics write files.

`diagnostics_gfs_h_create` runs after parallel codegen. It filters
`nrpy.grid.glb_gridfcs_dict` for gridfunctions whose group is `DIAG`, sorts them
case-insensitively, and writes `diagnostics/diagnostic_gfs.h`. The generated
header contains one enum token per diagnostic gridfunction, a
`TOTAL_NUM_DIAG_GFS` counter, and `diagnostic_gf_names[]` with C/C++
initializer handling through `DIAG_INIT`.

`diagnostic_gfs_set` is the normal GR producer for diagnostic channels. It
registers `DIAG_HAMILTONIAN`, `DIAG_MSQUARED`, `DIAG_LAPSE`, `DIAG_W`,
`DIAG_GRIDINDEX`, `DIAG_RBARDD`, optional `DIAG_T4UU`, and optional
`DIAG_PSI4_RE/IM`. At runtime it loops over grids, calls `Ricci_eval` or
`Ricci_eval_host`, calls `constraints_eval`, optionally calls `psi4`, applies
inner boundary conditions to interpolation-sensitive diagnostic channels, and
copies lapse, conformal factor, and grid index into `diagnostic_gfs`.

Nearest diagnostics are a dispatcher plus three helper samplers. Users select
`which_gfs_0d`, `which_gfs_1d`, and `which_gfs_2d` in the generated
`diagnostics_nearest()` user-edit block; the dispatcher loops over grids and
passes the caller-owned `diagnostic_gfs` buffers to the helpers. The 0D helper
`diagnostics_nearest_grid_center` appends one row to a filename that records
grid, coordinate-system name, and convergence factor, choosing the grid-center
index for Cartesian-like systems and `NGHOSTS` in the radial direction for
spherical, cylindrical, and SymTP-style systems. The 1D helper samples nearest
y- and z-axis lines, converts native coordinates through `xx_to_Cart`, sorts by
physical axis coordinate, and writes per-time files `out1d-y-*` and `out1d-z-*`.
The 2D helper samples nearest xy and yz planes, including multi-slice cases
such as opposite phi quadrants, and writes per-time `out2d-xy-*` and
`out2d-yz-*` files.

`diagnostics_nearest_common.h` supplies the shared text-output contract: time
comments use `# [time] = ...`, headers list coordinate columns plus diagnostic
names, rows use scientific notation, and `open_outfile` builds filenames with
or without embedded time. Files are opened with `"w"` only at `commondata->nn ==
0`; later calls append.

Volume diagnostics use `diagnostics_volume_integration()` and the copied
`diagnostics_volume_integration_helpers.h`. The generated routine builds a
small recipe array with user-editable spherical include/exclude rules and
integrand specs. It then calls `diags_integration_execute_recipes` with
`gridfuncs_diags`; the default examples integrate squared Hamiltonian and
momentum constraints over the whole domain and outside a radius. The
coordinate-specialized `sqrt_detgammahat_d3xx_volume_element` helper evaluates
`sqrt(detgammahat) * abs(dxx0 * dxx1 * dxx2)` by reference for ordinary
reference metrics. For `GeneralRFM` it is intentionally inert because the
active integration path reads the `DETGAMMAHATGF`-backed volume element from
the helper header.

Raytracing output is an optional diagnostics-side stage-1 export.
`output_raytracing_data` writes a time-stamped binary stage-1 payload through a
unique temporary sibling and installs it with `link()` so an existing final file
is not overwritten. Before writing payload records it refreshes same-slice
Ricci/RHS data with `Ricci_eval(...)` and `rhs_eval(...)`, then evaluates final
Cartesian coordinates, the ten unique covariant four-metric components, and the
forty unique four-Christoffel components at interior logical-grid points. The
header records binary64, little-endian, format-version, coordinate-system,
conformal-factor convention, logical-grid extents, offsets, component names,
and that ghost zones are excluded. `combine_raytracing_time_slices.py` is the
stage-2 containerizer: it strictly parses and validates stage-1 files, sorts
them by physical simulation time, validates compatibility and unique times,
optionally records coordinate-table and axisymmetry metadata, then atomically
writes a read-only combined container. It copies stage-1 point payloads; it does
not recompute metrics, Christoffels, transforms, or a true spatial index.

`progress_indicator` registers `start_wallclock_time` and
`output_progress_every`. Generated C initializes the wall-clock reference at
`nn_0`, returns early unless progress output is enabled and `nn` is on the
requested cadence, computes elapsed time and ETA when enabled, and writes the
default single-line stderr status with iteration, time, percent complete,
inverse timestep, physical time per hour, and ETA. `diagnostics()` calls it
after the diagnostics-output conditional, so progress cadence is independent of
diagnostics cadence.

`register_CFunctions` in `checkpointing.py` registers `read_checkpoint` and
`write_checkpoint`, plus the commondata parameter `checkpoint_every`.
Checkpoint filenames are `checkpoint-conv_factor%.2f.dat`. A nonpositive
`checkpoint_every` disables writes; otherwise `write_checkpoint()` uses the same
nearest-output-time test as diagnostics. Writes are fatal on open or partial
`fwrite` failure. CUDA writes first copy device params and every evolved
gridfunction to host. For each grid, the writer stores `params_struct`, a
selected point count, selected flat point indices, and a compact
`NUM_EVOL_GFS * count` payload. With multipatch enabled, selected points are
owned interiors, buffer-zone points, and outer-boundary points; otherwise every
point is selected. MoL intermediate-stage storage is freed before compacting and
reallocated after the grid payload is written.

`read_checkpoint()` returns `0` when the named checkpoint file does not exist.
If `griddata == NULL`, it reads only `commondata`, closes the file, and returns
`1`; this supports rebuilding grid metadata before reading the payload. A full
read verifies every `FREAD`, checks per-grid dimensions and `CoordSystem_hash`
against the rebuilt grid, rejects invalid selected-point counts, reallocates
host or pinned/device `y_n_gfs`, zero-fills the rebuilt evolved storage, scatters
the compact payload into `IDX4pt(gf, point)` positions, and copies restored
data back to the device when CUDA is active. At the end it sets `t_0 = time` and
`nn_0 = nn` so progress and restart bookkeeping start from the restored state.

BHaHAHA checkpointing is guarded by `enable_bhahaha`. Writes sanitize pointer
fields in a `commondata` copy before serializing it, validate that each horizon
has either all three previous-shape arrays or none, and then store a per-horizon
presence byte plus `prev_horizon_m1/m2/m3` data at the finest multigrid
resolution. Reads sanitize deserialized commondata pointers, read the presence
bytes, validate the horizon-shape dimensions, allocate the three previous-shape
arrays, and refill them from the checkpoint. This keeps raw pointer values from
the old process out of the restarted runtime while preserving horizon-history
data.

## Sources

- [diagnostics.py](../../../nrpy/infrastructures/BHaH/diagnostics/diagnostics.py) - `register_all_diagnostics`, `_register_CFunction_diagnostics`
- [diagnostic_gfs_h_create.py](../../../nrpy/infrastructures/BHaH/diagnostics/diagnostic_gfs_h_create.py) - `diagnostics_gfs_h_create`
- [diagnostic_gfs_set.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostic_gfs_set.py) - `register_CFunction_diagnostic_gfs_set`
- [diagnostics_nearest.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostics_nearest.py) - `register_CFunction_diagnostics_nearest`
- [diagnostics_nearest_grid_center.py](../../../nrpy/infrastructures/BHaH/diagnostics/diagnostics_nearest_grid_center.py) - `register_CFunction_diagnostics_nearest_grid_center`
- [diagnostics_nearest_1d_y_and_z_axes.py](../../../nrpy/infrastructures/BHaH/diagnostics/diagnostics_nearest_1d_y_and_z_axes.py) - `register_CFunction_diagnostics_nearest_1d_y_and_z_axes`, `bhah_axis_configs`
- [diagnostics_nearest_2d_xy_and_yz_planes.py](../../../nrpy/infrastructures/BHaH/diagnostics/diagnostics_nearest_2d_xy_and_yz_planes.py) - `register_CFunction_diagnostics_nearest_2d_xy_and_yz_planes`, `bhah_plane_configs`
- [diagnostics_nearest_common.h](../../../nrpy/infrastructures/BHaH/diagnostics/diagnostics_nearest_common.h) - `open_outfile`, `diag_write_header`, `diag_write_row`
- [diagnostics_volume_integration.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostics_volume_integration.py) - `register_CFunction_diagnostics_volume_integration`
- [diagnostics_volume_integration_helpers.h](../../../nrpy/infrastructures/BHaH/diagnostics/diagnostics_volume_integration_helpers.h) - `diags_integration_execute_recipes`, `DETGAMMAHATGF`
- [sqrt_detgammahat_d3xx_volume_element.py](../../../nrpy/infrastructures/BHaH/diagnostics/sqrt_detgammahat_d3xx_volume_element.py) - `register_CFunction_sqrt_detgammahat_d3xx_volume_element`
- [output_raytracing_data.py](../../../nrpy/infrastructures/BHaH/diagnostics/output_raytracing_data.py) - `register_CFunction_output_raytracing_data`
- [combine_raytracing_time_slices.py](../../../nrpy/infrastructures/BHaH/diagnostics/combine_raytracing_time_slices.py) - `parse_stage1_file`, `write_combined_file_atomically`, `main`
- [progress_indicator.py](../../../nrpy/infrastructures/BHaH/diagnostics/progress_indicator.py) - `register_CFunction_progress_indicator`
- [checkpointing.py](../../../nrpy/infrastructures/BHaH/checkpointing.py) - `register_CFunction_read_checkpoint`, `register_CFunction_write_checkpoint`, `register_CFunctions`
- [two_blackholes_collide.py](../../../nrpy/examples/two_blackholes_collide.py) - `BHaH.diagnostics.diagnostic_gfs_h_create.diagnostics_gfs_h_create`, `BHaH.diagnostics.progress_indicator.register_CFunction_progress_indicator`
- [SOURCES.md](../../../raw/SOURCES.md) - `infrastructure-modules-and-embedded-headers`

## See Also

- [BHaH](index.md)
- [GR Application Wiring](gr-application-wiring.md)
- [Black Hole Evolution](../../examples/black-hole-evolution.md)
- [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
