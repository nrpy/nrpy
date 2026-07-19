# Diagnostics Output And Checkpointing

> Explain BHaH diagnostics scheduling, temporary diagnostic buffers, raytracing export, progress output, and checkpoint/restart files. Status: confirmed. Last reconciled: 07-19-2026
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
checkpoint metadata, grid shape, coordinate-system hash, allocation sizes, and
serialized point indices before reallocating and restoring `y_n_gfs`.

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

`register_CFunction_read_checkpoint` in `read_checkpoint.py` and
`register_CFunction_write_checkpoint` in `write_checkpoint.py` independently
register the generated `read_checkpoint()` and `write_checkpoint()` CFunctions.
Within this split path, the writer registrar owns the commondata
`checkpoint_every` `CodeParameter`; checkpoint-using examples call the reader
and writer registrars separately. Checkpoint filenames remain
`checkpoint-conv_factor%.2f.dat`. A
nonpositive `checkpoint_every` disables writes; otherwise `write_checkpoint()`
uses the same nearest-output-time test as diagnostics. Writes are fatal on open
or partial `fwrite` failure. CUDA writes first copy device params and every
evolved gridfunction to host. For each grid, the writer stores `params_struct`,
a selected point count, selected flat point indices, and a compact
`NUM_EVOL_GFS * count` payload. With multipatch enabled, selected points are
owned interiors, buffer-zone points, and outer-boundary points; otherwise every
point is selected. MoL intermediate-stage storage is freed before compacting and
reallocated after the grid payload is written.

Current `BHaH` package aggregation and updated checkpoint-using examples use
the split owners. The tracked `checkpointing.py` nevertheless remains directly
importable and still defines its own older reader, writer, and combined
`register_CFunctions` wrapper; that old writer independently registers another
`checkpoint_every` owner. The module is not a forwarding or deprecation shim,
and its reader lacks the split reader's checks described below. Explicit legacy
submodule imports can therefore still reach a duplicate, less-hardened
implementation; no equivalence, compatibility, or current-canonical guarantee
applies to that direct-import path.

`read_checkpoint()` returns `0` when the named checkpoint file does not exist.
After reading `commondata`, it requires `NUMGRIDS` in `1..MAXNUMGRIDS`. With
BHaHAHA enabled, it also checks the deserialized horizon count against the
`bhahaha_params_and_data` capacity and the positive multigrid-resolution count
against both angular-dimension array capacities before those values control
pointer sanitization or loops. If `griddata == NULL`, the reader then closes the
file and returns `1`; metadata-only restart therefore retains these
commondata-level checks while deferring payload restoration until grids have
been rebuilt.

For BHaHAHA horizon-history payloads, the reader bounds the finest-resolution
index, requires positive angular dimensions, and rejects overflow in the
angular point product and `REAL` allocation size. For each rebuilt grid, it
requires positive dimensions, checks their products against `SIZE_MAX`, and
requires the point total to fit the retained `int` representation. It verifies
checkpoint dimensions and `CoordSystem_hash` against the rebuilt grid, bounds
the compact selected-point count by the rebuilt point total, rejects overflow
in evolved-gridfunction element counts and allocation byte counts, and skips
compact allocations and reads when that count is zero. Before scatter it checks
every serialized flat point index against the rebuilt grid, then uses explicit
`size_t` offsets for both destination and compact-payload indexing. A successful
full read reallocates host or pinned/device `y_n_gfs`, zero-fills it, scatters
the compact values, copies restored data to the device when CUDA is active, and
sets `t_0 = time` and `nn_0 = nn`.

BHaHAHA checkpointing remains guarded by `enable_bhahaha`. Writes sanitize
pointer fields in a `commondata` copy, require each horizon to have all three
previous-shape arrays or none, and store a presence byte plus
`prev_horizon_m1/m2/m3` data at the finest multigrid resolution. The split
writer does not have equivalent reader hardening: its BHaHAHA path indexes the
angular-dimension arrays and computes their `size_t` product before checking the
finest-resolution index and positive dimensions, and it does not bound those
indices by array capacity or check product/allocation overflow. Malformed
writer-side BHaHAHA metadata is therefore not guaranteed to fail safely. Reads
still sanitize deserialized pointers, allocate present shape arrays, and refill
them from the checkpoint, preserving horizon-history data without reusing old
process pointer values.

Claim evidence:
- Claim: Package aggregation and updated examples use separate checkpoint reader/writer registrars, with split writer owning `checkpoint_every` for that path; split reader requires valid `NUMGRIDS`, BHaHAHA capacities and finest-resolution index, positive overflow-safe dimensions/products, representable evolved-gridfunction and allocation counts, bounded compact counts and serialized point indices, zero-count-safe reads, and `size_t` scatter offsets, while split writer provides no matching BHaHAHA validation-order, capacity, or overflow guarantee; direct import of retained `checkpointing.py` still reaches an older duplicate whose writer independently registers `checkpoint_every` and which has no forwarding, deprecation, equivalence, compatibility, or current-canonical guarantee.
- Role: descriptive behavior
- Deciding authority: registered primary code `nrpy/infrastructures/BHaH/read_checkpoint.py::register_CFunction_read_checkpoint`, `nrpy/infrastructures/BHaH/write_checkpoint.py::register_CFunction_write_checkpoint`, and `nrpy/infrastructures/BHaH/checkpointing.py::register_CFunctions`
- Corroboration: `nrpy/infrastructures/BHaH/__init__.py` package import list and representative `nrpy/examples/blackhole_spectroscopy.py` split registrar calls; no independent registered test exercises malformed-input, direct-import, or restart paths
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, clang-format 22.1.8; backend=OpenMP C and CUDA source; precision=not-applicable; GPU=not-run; restart=not-run; distributed=not-applicable; error_path=not-run; options=default reader/writer source baselines plus BHaHAHA reader ordering assertions; date=07-19-2026`

Validation here is source-inspection scoped. Current Ubuntu/macOS codegen jobs
generate and build default BHaH examples without running a
write/restart/read sequence. File integrity, restored values, CUDA restart,
multipatch selection, and BHaHAHA horizon-history restart are therefore
`not-run` runtime outcomes in this KB audit. The four default reader/writer
`.c`/`.cu` baselines were regenerated in an isolated copy, matched the intended
files byte-for-byte, and passed a second fresh-process comparison. That proves
normalized emitted-source regression only; it does not establish compilation,
restart, malformed-input behavior, BHaHAHA source variants, GPU execution, or
runtime results.

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
- [BHaH package initializer](../../../nrpy/infrastructures/BHaH/__init__.py) - package import list for `read_checkpoint` and `write_checkpoint`
- [read_checkpoint.py](../../../nrpy/infrastructures/BHaH/read_checkpoint.py) - `register_CFunction_read_checkpoint`
- [write_checkpoint.py](../../../nrpy/infrastructures/BHaH/write_checkpoint.py) - `register_CFunction_write_checkpoint`
- [checkpointing.py](../../../nrpy/infrastructures/BHaH/checkpointing.py) - retained `register_CFunction_read_checkpoint`, `register_CFunction_write_checkpoint`, `register_CFunctions`
- [blackhole_spectroscopy.py](../../../nrpy/examples/blackhole_spectroscopy.py) - `BHaH.read_checkpoint.register_CFunction_read_checkpoint`, `BHaH.write_checkpoint.register_CFunction_write_checkpoint`
- [two_blackholes_collide.py](../../../nrpy/examples/two_blackholes_collide.py) - `BHaH.diagnostics.diagnostic_gfs_h_create.diagnostics_gfs_h_create`, `BHaH.diagnostics.progress_indicator.register_CFunction_progress_indicator`
- [SOURCES.md](../../../raw/SOURCES.md) - `infrastructure-modules-and-embedded-headers`
- [main.yml](../../../.github/workflows/main.yml) - `codegen-ubuntu`, `codegen-mac`

## See Also

- [BHaH](index.md)
- [GR Application Wiring](gr-application-wiring.md)
- [Black Hole Evolution](../../examples/black-hole-evolution.md)
- [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
