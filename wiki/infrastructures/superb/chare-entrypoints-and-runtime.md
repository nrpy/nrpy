# Chare Entrypoints And Runtime

> Structural map of superB Charm++ chares, entry methods, proxies, SDAG, reductions, checkpointing, and PUP contracts. · Status: confirmed · Last reconciled: 07-07-2026
> Up: [superB](index.md)

## Summary

superB emits a Charm++ program around one `Main` mainchare, one required
three-dimensional `Timestepping` chare array, and optional `Interpolator3d`
and `Horizon_finder` arrays when the BHaHAHA-facing service path is enabled.
Runtime communication is by generated Charm++ proxy handles and entry methods,
not direct object pointers. The Python generators own the handwritten source
text for generated `.h`, `.cpp`, and `.ci` files; generated project files
themselves remain output artifacts, not KB source evidence.

## Detail

Charm++ vocabulary is used in the ordinary upstream sense: `.ci` files declare
chares, chare arrays, entry methods, messages, readonly variables, and array
maps; generated proxy classes invoke entry methods asynchronously; SDAG
expresses ordered waits across entry methods; reductions combine
contributions; and PUP routines serialize object state for migration and
checkpoint/restart. superB behavior is defined by the local NRPy generators
that emit those constructs.

`main_chare.py` registers `Nchare0`, `Nchare1`, and `Nchare2` as commondata
CodeParameters and emits `commondata_object.h`, `main.h`, `main.cpp`, and
`main.ci`. `CommondataObject` wraps `commondata_struct` and PUPs it with
`pup_commondata_struct`. `Main` derives from `CBase_Main`, stores wall-clock
state, and conditionally PUPs startup and completion fields when Charm++
checkpointing is enabled.

`Main` is the startup coordinator. `main.ci` declares readonly proxy handles
for `mainProxy` and `timesteppingArray`; BHaHAHA builds also declare readonly
`horizon_finderProxy` and `interpolator3dArray`. These readonly variables are
startup-published Charm globals or handles, not shared mutable transport state.
Charm++ initializes readonly variables in the mainchare startup path and
broadcasts them; assignments outside the mainchare constructor are invalid
under the readonly contract even though the translator cannot enforce that
restriction.

`Main` assigns all runtime handles before dependent use. In non-BHaHAHA
builds, it assigns `mainProxy = thisProxy`, creates `timesteppingArray`, then
calls `timesteppingArray.start()`. In BHaHAHA builds, it assigns `mainProxy`,
creates and assigns `horizon_finderProxy`, creates and assigns
`timesteppingArray`, creates and assigns `interpolator3dArray`, then calls
`timesteppingArray.start()`. `Main` waits for both `timestepping_done()` and
`horizon_finder_done()` before `done()` when the horizon service is enabled.

Generated names have different roles. `CProxy_*` values are proxy handles.
`thisProxy` is the generated proxy handle stored in `CBase_*`; it is not a
direct object pointer. `CkIndex_*` names an entry-point index used by generated
callbacks and reduction targets, such as
`CkIndex_Timestepping::report_sums_for_volume(NULL)`. `CkArrayIndex3D` and
`CkIndex3D` are array-coordinate representations; for `Timestepping` and
`Interpolator3d`, `thisIndex.x`, `thisIndex.y`, and `thisIndex.z` are the
element coordinates.

The same proxy can broadcast or target one element depending on syntax.
Unindexed array-proxy calls broadcast to the whole array, as in
`timesteppingArray.start()`, `thisProxy.diagnostics_ckio(...)`, and
`interpolator3dArray.start_interpolation(...)`. Indexed sends target one
element, as in `thisProxy[CkArrayIndex3D(...)]`, `thisProxy(0, 0, 0)`,
`interpolator3dArray[CkArrayIndex3D(...)]`, and
`thisProxy[CkArrayIndex1D(thisIndex)]`.

`timestepping_chare.py` emits the required `Timestepping` `array [3D]`. The
header generator defines `class Timestepping : public CBase_Timestepping` with
`Timestepping_SDAG_CODE`, commondata, global and chare-local griddata pointers,
diagnostic file handles, loop counters, and optional checkpoint/BHaHAHA state.
The `.ci` generator declares the constructor, long-running `start()` SDAG
entry, diagnostic file-session entries, ghost and nonlocal-inner-boundary
receive entries, `continue_timestepping()`, volume-reduction reporting,
optional `startCheckpoint()` and `recvCheckPointDone()`, and optional
completion reduction targets.

`Interpolator3d` is an optional `array [3D]` service chare emitted by
`interpolator3d_chare.py`. It derives from `CBase_Interpolator3d`, embeds
`Interpolator3d_SDAG_CODE`, receives source gridfunction buffers, starts
interpolation, aggregates interpolation buffers, reports concatenation
completion, and PUPs migration/checkpoint state. `Horizon_finder` is an
optional `array [1D]` service chare emitted by `horizon_finder_chare.py`. It
derives from `CBase_Horizon_finder`, embeds `Horizon_finder_SDAG_CODE`, runs a
`start()` SDAG entry, receives interpolation results, signals horizon
completion, and reduces completion to `mainProxy.horizon_finder_done()`.

The `Timestepping` SDAG entry is the structural runtime spine. It orders
initialization, staged data synchronization, diagnostic sessions, optional
service-array messages, optional checkpoint coordination, Method of Lines
substeps, and final completion. Only SDAG entries containing `when` get the
generated direct-call prohibition; entries without `when` receive normal local
wrapper code. Incoming SDAG entries initialize `SDAG::Dependency` as needed,
then buffer closures or messages with `pushBuffer`. A reached `when` calls
`tryFindMessage`: if matching buffers already exist, it consumes them and
continues; otherwise it registers an `SDAG::Continuation` for a later arrival.
If multiple `when` waits can match the same incoming entry, delivery choice is
not specified unless reference numbers distinguish the waits.

SDAG reference-number behavior is split by entry shape. For marshalled SDAG
entries that need a reference number, generated code uses the first argument as
the reference number when the closure does not already have one. For message
and callback paths, reference numbers travel through `CkSetRefNum` and
`CkGetRefNum`.

superB PUP support has two layers. Chare classes PUP their own object state:
`Main` conditionally PUPs startup/completion fields; `Timestepping` PUPs
commondata, global griddata, chare-local griddata, and optional state; and the
optional service chares PUP their SDAG-visible state and owned buffers. The
shared C-function registration in `superB_pup.py` emits struct-level PUP
routines for commondata, params, boundary, MoL, chare communication,
diagnostics, nonlocal-inner-boundary, temporary-buffer, `griddata`, and
`griddata_chare` state. In current Charm++ generated code, `_sdag_pup(PUP::er&)`
is the real SDAG serializer through the generated recursive PUP path and PUPs
`__dep`; generated `__sdag_pup(PUP::er&)` is an empty compatibility stub.

Reductions appear in four distinct superB roles:

- Volume data reduction: each `Timestepping` element flattens volume and
  integrand values, then calls `contribute(outdoubles,
  CkReduction::sum_double, cb)`. The callback is
  `CkIndex_Timestepping::report_sums_for_volume(NULL)` on
  `Timestepping[0,0,0]`; `report_sums_for_volume` is not declared
  `[reductiontarget]`. It receives a `CkReductionMsg *` as a message-client
  callback and parses the summed `double` payload.
- Zero-data completion reductions: BHaHAHA completion uses zero-data
  reductions with `CkReductionTarget` callbacks so the root `Timestepping` or
  root `Horizon_finder` element can notify `Main`.
- Checkpoint coordination reduction: when checkpointing is enabled, all
  `Timestepping` elements reduce to `startCheckpoint()` on
  `Timestepping[0,0,0]`; that root calls `CkStartCheckpoint("log", cb)` and
  later waits for `recvCheckPointDone()`.
- CkIO completion callbacks: diagnostic write and close callbacks receive
  `CkReductionMsg *` as completion messages from the CkIO path. They are not
  user volume-reduction payloads.

`CkReduction::reducerType` is Charm++'s reducer selector type.
`CkReduction::sum_double` is the reducer used by superB volume diagnostics to
sum vectors of `double` values. Zero-data completion and checkpoint
coordination reductions use no user data payload.

Placement is separate from decomposition. `Nchare0`, `Nchare1`, and `Nchare2`
define the grid/chare decomposition and array extents; they are not PE
placement controls. `RRMap_with_offset` is emitted only for BHaHAHA builds
without Charm++ checkpointing. It is a `CkArrayMap` group whose `procNum` maps
flattened `(x, y, z)` coordinates plus an offset modulo `CkNumPes()`, and it is
installed through `CkArrayOptions::setMap` for the `Timestepping` and
`Interpolator3d` arrays. A map controls each element's initial/home PE; after
migration, current location is managed by Charm++ location state, not by
recomputing `RRMap_with_offset`. Active ck-ldb/AtSync load balancing is not
enabled by default in the superB generators or current CI coverage.

Checkpoint support is present but restart is not CI-proven. The current
workflow builds elliptic, spectroscopy, and collision superB projects, then
runs only the collision executable; it does not runtime-test restart. Charm++
disk checkpointing writes to the directory argument passed to
`CkStartCheckpoint`, so in the generated superB path the checkpoint directory is
`log`. `ck-cp` is not the disk checkpoint path. The `[migratable] Main`
checkpoint/restart risk remains unresolved. The generated migration constructor
also writes `mainProxy = thisProxy`; treat that as part of the unresolved
checkpoint caveat, not as a validated post-startup readonly update pattern.
`CkLocMgr::pup` also has a migration-buffer checkpoint risk: current Charm++
source aborts if the location manager is pupped while pending migration
messages are buffered.

## Sources

- [main_chare.py](../../../nrpy/infrastructures/superB/main_chare.py) - `Nchare0`, `Nchare1`, `Nchare2`, `output_commondata_object_h_and_main_h_cpp_ci`, `output_main_h`, `output_main_cpp`, `output_main_ci`, `RRMap_with_offset`
- [timestepping_chare.py](../../../nrpy/infrastructures/superB/timestepping_chare.py) - `output_timestepping_h`, `output_timestepping_cpp`, `output_timestepping_ci`, `output_timestepping_h_cpp_ci_register_CFunctions`, `generate_PUP_code`
- [interpolator3d_chare.py](../../../nrpy/infrastructures/superB/interpolator3d_chare.py) - `output_griddata_object_h`, `output_interp_buf_msg_h`, `output_interpolator3d_h_cpp_ci`, `Interpolator3d`
- [horizon_finder_chare.py](../../../nrpy/infrastructures/superB/horizon_finder_chare.py) - `output_horizon_finder_h_cpp_ci`, `Horizon_finder`
- [superB_pup.py](../../../nrpy/infrastructures/superB/superB/superB_pup.py) - `pup_griddata`, `pup_griddata_chare`, `register_CFunction_superB_pup_routines`
- [main.yml](../../../.github/workflows/main.yml) - `charmpp-validation`
- [Charm++ manual.rst](https://github.com/charmplusplus/charm/blob/main/doc/charm%2B%2B/manual.rst) - readonly variables, array maps, reductions, PUP, checkpoint/restart
- [Charm++ xi-Chare.C](https://github.com/charmplusplus/charm/blob/main/src/xlat-i/xi-Chare.C) - `CBase_*`, `thisProxy`, generated recursive `_sdag_pup`
- [Charm++ CParsedFile.C](https://github.com/charmplusplus/charm/blob/main/src/xlat-i/sdag/CParsedFile.C) - `_sdag_pup(PUP::er&)`, `__sdag_pup(PUP::er&)`
- [Charm++ CEntry.C](https://github.com/charmplusplus/charm/blob/main/src/xlat-i/sdag/CEntry.C) - SDAG direct-call guard, closure/message buffering, marshalled-entry reference numbers
- [Charm++ When.C](https://github.com/charmplusplus/charm/blob/main/src/xlat-i/sdag/constructs/When.C) - `when` buffer matching and continuation registration
- [Charm++ sdag.h](https://github.com/charmplusplus/charm/blob/main/src/ck-core/sdag.h) - `SDAG::Dependency`, `SDAG::Buffer`, `SDAG::Continuation`
- [Charm++ xi-Parameter.C](https://github.com/charmplusplus/charm/blob/main/src/xlat-i/xi-Parameter.C) - SDAG callback/message reference-number unmarshalling
- [Charm++ ckcallback.C](https://github.com/charmplusplus/charm/blob/main/src/ck-core/ckcallback.C) - `CkSetRefNum` callback delivery
- [Charm++ ckarrayindex.h](https://github.com/charmplusplus/charm/blob/main/src/ck-core/ckarrayindex.h) - `CkIndex3D`, `CkArrayIndex3D`
- [Charm++ ckarray.h](https://github.com/charmplusplus/charm/blob/main/src/ck-core/ckarray.h) - `ArrayElementT`, `thisIndex`, array reductions
- [Charm++ ckreduction.h](https://github.com/charmplusplus/charm/blob/main/src/ck-core/ckreduction.h) - `CkReduction::reducerType`, `CkReduction::sum_double`, `CkReductionMsg`
- [Charm++ cklocation.h](https://github.com/charmplusplus/charm/blob/main/src/ck-core/cklocation.h) - `CkArrayMap`, home PE, location manager roles
- [Charm++ cklocation.C](https://github.com/charmplusplus/charm/blob/main/src/ck-core/cklocation.C) - `CkLocMgr::pup` pending-migration abort

## See Also

- [superB](index.md)
- [Lifecycle And Project Assembly](lifecycle-and-project-assembly.md)
- [Diagnostics And Observables](diagnostics-and-observables.md)
- [Grids, Boundaries, MoL, And Initial Data](grids-boundaries-mol-and-initial-data.md)
- [GR, BHaHAHA, Psi4, And Interpolation](gr-bhahaha-psi4-and-interpolation.md)
- [Generated Project CI](../../validation/generated-project-ci.md)
- [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
