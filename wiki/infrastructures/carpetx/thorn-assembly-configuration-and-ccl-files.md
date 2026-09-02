# CarpetX Thorn Assembly, Configuration, And CCL Files

> Explain how CarpetX writes Cactus CCL files, `make.code.defn`, and thorn-local C++ sources from registered generator state. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [CarpetX](index.md)

## Summary

CarpetX thorn assembly writes Cactus control files, Cactus build metadata, and
registered C++ functions as generated products. Treat emitted `interface.ccl`,
`param.ccl`, `schedule.ccl`, `configuration.ccl`, `src/make.code.defn`, and
`src/*.cxx` as products of the Python generators, not source evidence for KB
claims. For CarpetX behavior, cite the CarpetX writer modules; use Cactus and
Einstein Toolkit documentation only as background for CCL, thorn, CarpetX, and
Loop terminology.

## Detail

Cactus documentation defines CCL files as the text files that describe thorn
configuration, and its configuration syntax documentation lists
`interface.ccl`, `param.ccl`, `schedule.ccl`, and optional
`configuration.ccl` as thorn-level files. CarpetX documentation adds that
CarpetX thorns may need explicit `configuration.ccl` requirements for CarpetX
APIs or schedule functions. NRPy CarpetX implements that shape through local
writers that create the thorn directory under the requested project directory.

`construct_configuration_ccl()` writes `configuration.ccl` with the generated
file warning and a required dependency line, `REQUIRES Loop CarpetX`. This is a
CarpetX-specific difference from ETLegacy's thorn-assembly page: CarpetX has a
dedicated `configuration.ccl` writer and explicitly requires both Loop and
CarpetX.

`construct_interface_ccl()` writes `interface.ccl` for one thorn. It emits the
generated-file warning, `implements: <thorn_name>`, `inherits: <inherits>`, and
the caller-provided `USES_INCLUDEs` block. When NewRadX support is enabled, it
adds `USES INCLUDE HEADER: newradx.hxx`.

The same interface writer emits public CarpetX gridfunction groups from
CarpetX gridfunction metadata. For evolution thorns, non-empty evolved
variables become `evol_variables` with `Timelevels=1` and an RHS tag pointing
to `evol_variables_rhs`; the matching RHS group also has `Timelevels=1` and
tags `InterpNumTimelevels=1`, `prolongation="none"`, and `checkpoint="no"`.
Non-empty auxiliary-evolved variables become `auxevol_variables` with
`Timelevels=1` and the same interpolation, prolongation, and checkpoint tags.
Non-empty auxiliary variables become `aux_variables` with `Timelevels=1`.

Parity tags are constructed locally in `interface_ccl.py` from CarpetX
gridfunction parity metadata. The emitted tag maps scalar and diagonal-tensor
parity values `0`, `4`, `7`, and `9` to `+1 +1 +1`; vector parity values `1`,
`2`, and `3` to sign flips in the x, y, and z directions; and off-diagonal
tensor parity values `5`, `6`, and `8` to paired sign flips for xy, xz, and yz
components.

`construct_param_ccl()` writes `param.ccl`. It scans registered `CFunction`
objects, keeps only functions whose `ET_thorn_name` matches the current thorn,
collects non-empty `ET_current_thorn_CodeParams_used` lists, deduplicates each
code parameter name, then emits the sorted parameter declarations from
`glb_code_params_dict` in a `restricted:` block. The returned list preserves
the deduplicated names registered to `param.ccl`.

`ScheduleCCL` stores one scheduled function entry: function name, schedule bin,
entry text, and whether the entry has already been output.
`construct_schedule_ccl()` writes `schedule.ccl` from caller-provided storage,
registered `ET_schedule_bins_entries`, and optional extra schedule entries. It
sorts registered functions by name, groups entries by each function's
`ET_thorn_name`, replaces `FUNC_NAME` placeholders with the registered function
name, and emits bins in this fixed CarpetX order: `STARTUP`, `BASEGRID`,
`CCTK_INITIAL`, `ODESolvers_RHS`, `ODESolvers_PostStep`, then any remaining
bins. Compared with ETLegacy, this order uses ODESolvers bins rather than the
MoL schedule bins, and it does not include ETLegacy's `Driver_BoundarySelect`
or `MoL_Register` steps.

`output_CFunctions_and_construct_make_code_defn()` performs the source-output
handoff. It filters registered `CFunction` objects by `subdirectory ==
thorn_name`, sorts them by `CFunction.name`, creates the thorn-local `src`
directory, writes each `CFunction.full_function` to `src/<name>.cxx`, and
writes `src/make.code.defn` with the matching sorted `SRCS` list. This is the
CarpetX counterpart to ETLegacy's generated `src/*.c`, but the CarpetX writer
emits `.cxx` files.

No current checked-in workflow command invokes `carpetx_wavetoy_thorns.py` or
`carpetx_baikal_thorns.py`. The `einsteintoolkit-validation` job invokes the
similarly named `carpet_*` ETLegacy generators instead. CarpetX generation,
Einstein Toolkit build, runtime, GPU, and restart outcomes are therefore
`not-run` in current CI evidence and in this KB audit; local generator code
decides only the emitted thorn shape described above.

## Sources

- [interface_ccl.py](../../../nrpy/infrastructures/CarpetX/interface_ccl.py) - `construct_interface_ccl`, `construct_parity_string`
- [param_ccl.py](../../../nrpy/infrastructures/CarpetX/param_ccl.py) - `construct_param_ccl`
- [schedule_ccl.py](../../../nrpy/infrastructures/CarpetX/schedule_ccl.py) - `ScheduleCCL`, `construct_schedule_ccl`
- [configuration_ccl.py](../../../nrpy/infrastructures/CarpetX/configuration_ccl.py) - `construct_configuration_ccl`
- [make_code_defn.py](../../../nrpy/infrastructures/CarpetX/make_code_defn.py) - `output_CFunctions_and_construct_make_code_defn`
- [main.yml](../../../.github/workflows/main.yml) - `einsteintoolkit-validation` (ETLegacy-only configured commands)
- [Cactus Users Guide chapter C1](https://www.cactuscode.org/documentation/usersguide/UsersGuidech9.html) - `C1.1.1 Thorns`, `C1.2 Anatomy of a Thorn`, `make.code.defn based thorn building`; accessed 07-12-2026
- [Cactus Users Guide chapter D2](https://www.cactuscode.org/documentation/usersguide/UsersGuidech12.html) - CCL configuration syntax; accessed 07-12-2026
- [CarpetX User Manual](https://einsteintoolkit.org/thornguide/CarpetX/CarpetX/documentation.html) - `Introduction` and CCL dependency background; accessed 07-12-2026

## See Also

- [CarpetX](index.md)
- [Code Parameters, Includes, And Loops](code-parameters-includes-and-loops.md)
- [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
- [C Function Registry](../../core/c-function-registry.md)
- [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)
- [ETLegacy Thorn Assembly And CCL Files](../etlegacy/thorn-assembly-and-ccl-files.md)
