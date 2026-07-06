# Chare Entrypoints And Runtime

> Structural map of superB Charm++ chares, entry methods, proxies, SDAG, reductions, checkpointing, and PUP contracts. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [superB](index.md)

## Summary

superB emits a Charm++ program around one `Main` mainchare, one required
three-dimensional `Timestepping` chare array, and optional `Interpolator3d`
and `Horizon_finder` arrays when the BHaHAHA-facing service path is enabled.
The Python generators own the handwritten source text for the generated
`.h`, `.cpp`, and `.ci` files; the generated project files themselves remain
output artifacts, not KB source evidence.

## Detail

Charm++ vocabulary is used in the ordinary sense from the upstream Charm++
repository: `.ci` files declare chares, chare arrays, entry methods, messages,
and readonly variables; generated proxy classes invoke entry methods
asynchronously; SDAG expresses ordered waits across entry methods; reductions
combine contributions; and PUP routines serialize object state for migration
and checkpoint/restart. superB's local behavior is defined by the NRPy
generators that emit those Charm++ constructs.

`main_chare.py` registers `Nchare0`, `Nchare1`, and `Nchare2` as commondata
CodeParameters and emits `commondata_object.h`, `main.h`, `main.cpp`, and
`main.ci` through `output_commondata_object_h_and_main_h_cpp_ci`. The emitted
`CommondataObject` wraps `commondata_struct` and PUPs it with
`pup_commondata_struct`. The emitted `Main` class derives from `CBase_Main`,
stores wall-clock state, and conditionally PUPs that state when Charm++
checkpointing is enabled.

`Main` is the startup coordinator. Its `.ci` module declares readonly proxies
for `mainProxy` and `timesteppingArray`; with BHaHAHA enabled it also declares
readonly `horizon_finderProxy` and `interpolator3dArray`. Its constructor
parses command-line and parameter-file input into commondata, creates the
`Timestepping` array from the commondata chare counts, optionally creates the
service arrays, and invokes `timesteppingArray.start()`. In BHaHAHA builds,
`Main` waits for both `timestepping_done()` and `horizon_finder_done()` before
calling `done()`. When BHaHAHA is enabled without Charm++ checkpointing, the
same generator emits `RRMap_with_offset` as a Charm++ array-map group for
placement of the three-dimensional arrays.

`timestepping_chare.py` emits the required `Timestepping` chare array. The
header generator defines `class Timestepping : public CBase_Timestepping` with
`Timestepping_SDAG_CODE`, commondata, global and chare-local griddata pointers,
diagnostic file handles, loop counters, and optional checkpoint/BHaHAHA state.
Its `.ci` generator declares `array [3D] Timestepping`, the constructor entry
method, the long-running `start()` SDAG entry method, diagnostic file-session
entrypoints, ghost and nonlocal-inner-boundary receive entrypoints,
`continue_timestepping()`, volume-reduction reporting, optional
`startCheckpoint()` and `recvCheckPointDone()`, and optional completion
reduction target entrypoints.

The `Timestepping` SDAG entry method is the structural runtime spine. It
orders initialization, staged data synchronization, diagnostic sessions,
optional service-array messages, optional checkpoint coordination, Method of
Lines substeps, and final completion. This page deliberately stops at the
entrypoint and synchronization structure; chare-local grid setup, ghost data
contents, MoL phase details, and initial-data mechanics belong to the grids and
boundaries leaf.

Reductions are used in two distinct structural roles. Volume diagnostics build
local vectors and call `contribute(..., CkReduction::sum_double, cb)`, with
`report_sums_for_volume(CkReductionMsg *msg)` as the callback target on
`Timestepping[0,0,0]`. BHaHAHA-enabled completion uses Charm++ reduction
callbacks to notify `Main` after all timestepping elements finish and after all
horizon-finder elements finish. Checkpointing, when enabled, is coordinated by
a reduction to `startCheckpoint()` on `Timestepping[0,0,0]`, which calls
`CkStartCheckpoint("log", cb)` and then waits for `recvCheckPointDone()`.

`Interpolator3d` is an optional `array[3D]` service chare emitted by
`interpolator3d_chare.py`. The generator emits `GriddataObject`,
`InterpBufMsg`, `interpolator3d.h`, `interpolator3d.cpp`,
`interpolation_buffer_utils.cpp`, and `interpolator3d.ci`. Structurally, the
class derives from `CBase_Interpolator3d`, embeds `Interpolator3d_SDAG_CODE`,
declares entrypoints for receiving source gridfunction buffers, starting
interpolation, aggregating interpolation buffers, and reporting concatenation
completion, and implements a PUP method for migration/checkpoint state. The
details of what is interpolated for BHaHAHA or Psi4 belong to the GR service
leaf.

`Horizon_finder` is an optional `array[1D]` service chare emitted by
`horizon_finder_chare.py`. The class derives from `CBase_Horizon_finder`,
embeds `Horizon_finder_SDAG_CODE`, declares a `start()` SDAG entry method plus
`ready_for_interpolation()`, `report_interpolation_results(...)`,
`horizon_finding_complete()`, and a completion reduction target that calls
`mainProxy.horizon_finder_done()`. Like `Interpolator3d`, it owns a PUP method
for object state and has its detailed service flow documented elsewhere.

superB PUP support has two layers. Chare classes PUP their own object state:
`Main` conditionally PUPs startup/completion fields, `Timestepping` calls
`CBase_Timestepping::pup(p)`, `__sdag_pup(p)`, `pup_commondata_struct`,
`pup_griddata`, and `pup_griddata_chare`, and the optional service chares PUP
their SDAG state and owned buffers. The shared C-function registration in
`superB_pup.py` emits struct-level PUP routines for commondata, params,
boundary, MoL, chare-communication, diagnostics, nonlocal-inner-boundary,
temporary-buffer, `griddata`, and `griddata_chare` state.

## Sources

- [main_chare.py](../../../nrpy/infrastructures/superB/main_chare.py) - `Nchare0`, `Nchare1`, `Nchare2`, `output_commondata_object_h_and_main_h_cpp_ci`, `output_main_h`, `output_main_cpp`, `output_main_ci`, `RRMap_with_offset`
- [timestepping_chare.py](../../../nrpy/infrastructures/superB/timestepping_chare.py) - `output_timestepping_h`, `output_timestepping_cpp`, `output_timestepping_ci`, `output_timestepping_h_cpp_ci_register_CFunctions`, `generate_PUP_code`
- [interpolator3d_chare.py](../../../nrpy/infrastructures/superB/interpolator3d_chare.py) - `output_griddata_object_h`, `output_interp_buf_msg_h`, `output_interpolator3d_h_cpp_ci`, `Interpolator3d`
- [horizon_finder_chare.py](../../../nrpy/infrastructures/superB/horizon_finder_chare.py) - `output_horizon_finder_h_cpp_ci`, `Horizon_finder`
- [superB_pup.py](../../../nrpy/infrastructures/superB/superB/superB_pup.py) - `pup_griddata`, `pup_griddata_chare`, `register_CFunction_superB_pup_routines`
- [quickstart.rst](https://github.com/charmplusplus/charm/blob/main/doc/quickstart.rst) - `Parallel "Hello World" with Charm++`
- [manual.rst](https://github.com/charmplusplus/charm/blob/main/doc/charm%2B%2B/manual.rst) - `Execution Model`, `Proxies and the charm interface file`, `Structured Dagger`, `Reductions on Chare Arrays`, `Serialization Using the PUP Framework`, `Checkpoint/Restart-Based Fault Tolerance`

## See Also

- [superB](index.md)
- [Lifecycle And Project Assembly](../bhah/lifecycle-and-project-assembly.md)
- [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
