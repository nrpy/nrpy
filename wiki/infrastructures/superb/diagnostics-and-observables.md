# Diagnostics And Observables

> Explain superB diagnostics dispatch, nearest-output helpers, CkIO output paths, and reduction-backed observables. · Status: confirmed · Last reconciled: 07-07-2026
> Up: [superB](index.md)

## Summary

superB diagnostics are staged through a generated `diagnostics()` dispatcher,
the `Timestepping` chare's scheduled output block, nearest-point sampling
helpers, output-only CkIO sessions for 1D/2D nearest files, direct `FILE *`
center output, and reduction-backed volume-integral reporting. The diagnostic
gridfunctions sampled by this machinery are chosen by the project variant: the
GR nearest dispatcher samples Hamiltonian and `MSQUARED` diagnostics, while the
NRPyElliptic nearest dispatcher samples residual and `uu` diagnostics.

## Detail

`register_all_diagnostics` registers the common top-level diagnostics driver,
copies helper headers needed by enabled diagnostics families, registers
nearest-sampling helpers for each coordinate system, optionally registers
volume-integration support, and contributes `diagnostic_struct` to
`griddata_struct`. The generated `diagnostics()` C function is a phase
dispatcher: setup phases call nearest setup routines, write phases call nearest
write routines, and the scheduling decision is left to the Charm++ runtime path.
It also registers the `diagnostics_output_every` commondata parameter consumed
by `Timestepping`.

Nearest diagnostics are organized as 0D, 1D, and 2D outputs. The center helper
chooses a coordinate-system-specific grid index nearest the physical center:
Cartesian-like systems use midpoint indices, while spherical, cylindrical, and
`SymTP` families use `NGHOSTS` in the radial-like direction and midpoints in
the other directions. The 1D helper builds samples along the nearest physical
`y` and `z` axes, converts native coordinates to Cartesian coordinates, sorts
points by the physical axis coordinate, records chare-local indices and byte
offsets in `diagnostic_struct`, and later writes rows for the selected
diagnostic gridfunctions. The 2D helper does the analogous setup and write path
for nearest `xy` and `yz` planes, including coordinate-family-specific slice
selection and Cartesian output coordinates.

`diagnostic_struct` is the persistent per-grid, per-chare bookkeeping object for
the 1D/2D nearest-output path. It stores total point counts, chare-owned point
counts, local point indices, per-point output offsets, per-row byte sizes, total
file sizes, and filename components for `out1d-y`, `out1d-z`, `out2d-xy`, and
`out2d-yz`. The phase constants in `superB.h` define the contract between
`Timestepping`, `diagnostics()`, and the nearest helpers:
`DIAGNOSTICS_SETUP_1D`, `DIAGNOSTICS_SETUP_2D`,
`DIAGNOSTICS_WRITE_CENTER`, `DIAGNOSTICS_WRITE_Y`,
`DIAGNOSTICS_WRITE_Z`, `DIAGNOSTICS_WRITE_XY`, `DIAGNOSTICS_WRITE_YZ`,
and `DIAGNOSTICS_VOLUME`.

During timestepping, the generated control flow computes
`write_diagnostics_this_step` from `commondata.time`, `commondata.dt`, and
`commondata.diagnostics_output_every`. On output steps it allocates temporary
per-grid `diagnostic_gfs`, initializes them to `NAN`, fills them with
`diagnostic_gfs_set`, and aliases MoL's diagnostic-output pointer to that
storage. Center diagnostics are then run directly through
`diagnostics_ckio(Ck::IO::Session(), DIAGNOSTICS_WRITE_CENTER)`, but the center
helper writes `out0d` through a normal `FILE *`; no CkIO session is involved.
The 1D and 2D outputs are routed through file opens, CkIO write sessions,
callback-triggered `diagnostics_ckio` calls, completion callbacks, and file
closes before timestepping resumes. `DIAGNOSTICS_VOLUME` is also invoked through
`diagnostics_ckio`, but it bypasses the regular row-output dispatcher and enters
the volume-reduction path.

superB uses CkIO only for nearest 1D/2D output files. On each output step, the
root diagnostics branch opens each `out1d-*` and `out2d-*` file. The generated
CkIO callbacks currently target `thisProxy`, the full `Timestepping` array
proxy; Charm++ treats `CkCallback(ep, CkArrayID)` as an array-broadcast
callback. The intended root branch starts the write session with the total byte
size saved in `diagnostic_struct`, broadcasts the session token to all
`Timestepping` elements through `diagnostics_ckio`, and closes the file after
CkIO reports completion, but callback delivery itself is not encoded as a
root-only target. The nearest helpers write headers and fixed-size rows with
`Ck::IO::write`.
Charm++ defines those write offsets as relative to the whole file, not to the
session offset, while the session itself is a byte window; superB's precomputed
header and row offsets must therefore fit inside that file-global window. CkIO
uses seek/read/write at the lowest level, so this output path requires seekable
files. The generated code constructs default `Ck::IO::Options` values for
opens, but this page does not prescribe CkIO option tuning. Current superB
control flow should not be read as checkpointing an open CkIO file: write
sessions complete, files close, SDAG reaches the `closed_*` callbacks, and only
then does the timestepping continuation and checkpoint-scheduling path proceed.

| Step | superB entry or call | Charm++ message/meaning |
| --- | --- | --- |
| open | root diagnostics branch calls `Ck::IO::open` | schedules file open with a `thisProxy` array-broadcast callback |
| file ready | `ready_*` | receives `FileReadyMsg`; generated callback target is array-shaped, while the intended consuming SDAG branch is root-owned |
| start session | `Ck::IO::startSession` | intended root branch starts byte-window session |
| session ready | `start_write_*` | receives `SessionReadyMsg` and broadcasts `diagnostics_ckio` |
| write | all `Timestepping` elements call `Ck::IO::write` | multi-PE row/header writes into session window |
| write complete | `test_written_*` | receives completion `CkReductionMsg *`; delete/ignore payload as completion, not reducer data |
| close | root calls `Ck::IO::close` | closes after completion |
| closed | `closed_*` | SDAG resumes only after close callback |

Volume diagnostics use Charm++ reductions rather than CkIO row writes.
`diagnostics_ckio` checks for `DIAGNOSTICS_VOLUME` before calling the generated
`diagnostics()` row dispatcher, so volume output goes directly to
`contribute_localsums_for_diagnostic_volume_integ`. Each chare evaluates the
default integration recipes against `diagnostic_gfs`, flattens one
`std::vector<double>` per contribution in recipe order, and stores each recipe
as selected volume followed by that recipe's integrands. It contributes the
flattened vector with `CkReduction::sum_double`. The callback is constructed as
`CkIndex_Timestepping::report_sums_for_volume(NULL)` on
`Timestepping(0,0,0)`, not as a typed `[reductiontarget]` callback, so the root
entry receives a real `CkReductionMsg *` payload and reconstructs doubles from
`getSize()` and `getData()`. This is distinct from CkIO callback messages in
`test_written_*` and `closed_*`: `test_written_*` receives the empty CkIO
completion `CkReductionMsg *`, and both handlers delete callback messages
instead of parsing reducer data. The reporting entry method reconstructs the
recipe results, writes `out3d-integrals-conv_factor...` files from the root
chare, and for NRPyElliptic updates `commondata.log10_current_residual` from
the residual RMS for the named `sphere_R_80` recipe.

The GR and NRPyElliptic nearest diagnostic variants share the same orchestration
and helper calls; they differ in the diagnostic gridfunctions selected in their
`USER-EDIT` blocks. The GR dispatcher samples `DIAG_HAMILTONIANGF` and
`DIAG_MSQUAREDGF` for 0D, 1D, and 2D nearest outputs. The NRPyElliptic
dispatcher samples `DIAG_RESIDUALGF` and `DIAG_UUGF` for the same output
dimensions. This page treats those as superB output choices and leaves the
equation definitions behind those gridfunctions to the equation pages.

## Sources

- [diagnostics.py](../../../nrpy/infrastructures/superB/diagnostics/diagnostics.py) - `register_all_diagnostics`, `_register_CFunction_diagnostics`
- [diagnostics_nearest_grid_center.py](../../../nrpy/infrastructures/superB/diagnostics/diagnostics_nearest_grid_center.py) - `get_center_index_exprs_for_coordsystem`, `register_CFunction_diagnostics_nearest_grid_center`
- [diagnostics_nearest_1d_y_and_z_axes.py](../../../nrpy/infrastructures/superB/diagnostics/diagnostics_nearest_1d_y_and_z_axes.py) - `register_CFunction_diagnostics_nearest_1d_y_and_z_axes`
- [diagnostics_nearest_2d_xy_and_yz_planes.py](../../../nrpy/infrastructures/superB/diagnostics/diagnostics_nearest_2d_xy_and_yz_planes.py) - `register_CFunction_diagnostics_nearest_2d_xy_and_yz_planes`
- [general_relativity/diagnostics_nearest.py](../../../nrpy/infrastructures/superB/general_relativity/diagnostics_nearest.py) - `register_CFunction_diagnostics_nearest`
- [nrpyelliptic/diagnostics_nearest.py](../../../nrpy/infrastructures/superB/nrpyelliptic/diagnostics_nearest.py) - `register_CFunction_diagnostics_nearest`
- [timestepping_chare.py](../../../nrpy/infrastructures/superB/timestepping_chare.py) - `generate_diagnostics_code`, `diagnostics_ckio`, `contribute_localsums_for_diagnostic_volume_integ`, `report_sums_for_volume`
- [superB.h](../../../nrpy/infrastructures/superB/superB/superB.h) - `diagnostic_struct`, `DIAGNOSTICS_*`
- [Charm++ libraries manual](https://github.com/charmplusplus/charm/blob/main/doc/libraries/manual.rst) - `CkIO`, `Using CkIO`, `Parallel Output API`
- [Charm++ language manual](https://github.com/charmplusplus/charm/blob/main/doc/charm%2B%2B/manual.rst) - `Reductions on Chare Arrays`, `Built-in Reduction Types`

## See Also

- [superB](index.md)
- [Chare Entrypoints And Runtime](chare-entrypoints-and-runtime.md)
- [GR, BHaHAHA, Psi4, And Interpolation](gr-bhahaha-psi4-and-interpolation.md)
- [Lifecycle And Project Assembly](lifecycle-and-project-assembly.md)
