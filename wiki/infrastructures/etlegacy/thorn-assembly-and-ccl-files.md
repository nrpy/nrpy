# Thorn Assembly And CCL Files

> Explain how ETLegacy writes Cactus CCL files, `make.code.defn`, and thorn-local C sources from registered generator state. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [ETLegacy](index.md)

## Summary

ETLegacy writes Cactus thorn control files, build metadata, and registered C
sources as generated output. Treat emitted `interface.ccl`, `param.ccl`,
`schedule.ccl`, `src/make.code.defn`, and `src/*.c` files as products of the
Python generators, not as source evidence for KB claims. Cite generator modules,
the `CFunction` registry contract, and generated-output boundary pages instead.

## Detail

Cactus thorn terminology and CCL file roles come from the Cactus documentation:
a thorn is an application module, and CCL files describe its interface,
parameters, schedule, and build inputs. ETLegacy implements that shape by
writing thorn directories under the generated project tree, so the durable facts
live in the ETLegacy Python writers and registry metadata.

`construct_interface_ccl()` writes `interface.ccl` for one thorn. It emits the
generated-file warning, an `implements:` declaration for `thorn_name`, an
`inherits:` declaration from the caller-provided `inherits` string, and the
caller-provided `USES_INCLUDEs` block for required includes or functions. When
`is_evol_thorn` is true, it appends declarations for Method of Lines,
boundary-condition, symmetry-table, and driver boundary-selection Cactus
functions. When `enable_NewRad` is true, it appends `ExtrapolateGammas` and
`NewRad_Apply` declarations.

The same writer appends public gridfunction group declarations for evolution
thorns. As observed source behavior, it calls
`gri.CarpetXGridFunction.gridfunction_lists()` and uses the returned evolved,
auxiliary-evolved, diagnostic, and auxiliary lists. Non-empty evolved lists
produce `evol_variables` with three timelevels and matching
`evol_variables_rhs` with one timelevel; when that evolved-list branch is
active, a non-empty auxiliary-evolved list also produces
`auxevol_variables` with one timelevel. A non-empty auxiliary list is checked
separately and produces `aux_variables` with three timelevels. The diagnostic
list is unpacked but not otherwise used by this writer.

`construct_param_ccl()` writes `param.ccl` and returns parameter names it
registered. It begins with the generated-file warning, optional
`shares_extends_str`, and a `restricted:` block. It scans `CFunction_dict` for
functions whose `ET_thorn_name` matches the requested thorn and whose
`ET_current_thorn_CodeParams_used` list is non-empty. It deduplicates those
names, sorts them alphabetically, looks each one up in
`par.glb_code_params_dict`, then emits type, name, default value, and an
all-values-accepted range note.

`ScheduleCCL` stores one schedule entry's function name, schedule bin, entry
text, and output flag. `construct_schedule_ccl()` writes `schedule.ccl` by
starting with generated-file text, then a storage allocation step from caller
input. It scans sorted `CFunction_dict` items, gathers each function's
`ET_schedule_bins_entries` under that function's `ET_thorn_name`, and prints a
warning when a registered function has no schedule metadata. Extra entries can
be appended for the requested thorn. Output order is the fixed bin sequence
`STARTUP`, `Driver_BoundarySelect`, `BASEGRID`, `CCTK_INITIAL`,
`MoL_Register`, `MoL_CalcRHS`, `MoL_PostStep`, and `MoL_PseudoEvolution`,
followed by any remaining bins. Entry text uses `FUNC_NAME` replacement before
being written.

`output_CFunctions_and_construct_make_code_defn()` performs the registered
C-source handoff. It filters `CFunction_dict` values by comparing each
`CFunction.subdirectory` to the requested `thorn_name`, sorts the matches by
`CFunction.name`, creates `project_dir/thorn_name/src`, writes each function's
`full_function` to `src/<name>.c`, and writes `src/make.code.defn` with the
corresponding `SRCS` list. This is the point where registry objects become
thorn-local C files and Cactus build metadata.

Validation scope is narrower than general Einstein Toolkit support. The
configured `einsteintoolkit-validation` job generates the ETLegacy Carpet
WaveToy and Baikal thorns, links them into the pinned ET environment, builds
that toolkit, and requests three testsuites. This workflow configuration proves
the job shape, not its latest success, and does not establish arbitrary thorn,
platform, GPU, or restart behavior. No generator, ET build, or testsuite was run
during this KB audit.

## Sources

- [interface_ccl.py](../../../nrpy/infrastructures/ETLegacy/interface_ccl.py) - `construct_interface_ccl`
- [grid.py](../../../nrpy/grid.py) - `GridFunction.gridfunction_lists`, `CarpetXGridFunction`
- [param_ccl.py](../../../nrpy/infrastructures/ETLegacy/param_ccl.py) - `construct_param_ccl`
- [schedule_ccl.py](../../../nrpy/infrastructures/ETLegacy/schedule_ccl.py) - `ScheduleCCL`, `construct_schedule_ccl`
- [make_code_defn.py](../../../nrpy/infrastructures/ETLegacy/make_code_defn.py) - `output_CFunctions_and_construct_make_code_defn`
- [c_function.py](../../../nrpy/c_function.py) - `CFunction`, `CFunction_dict`, `full_function`, `ET_schedule_bins_entries`, `ET_current_thorn_CodeParams_used`
- [main.yml](../../../.github/workflows/main.yml) - `einsteintoolkit-validation`
- [Cactus Users Guide chapter C1](https://www.cactuscode.org/documentation/usersguide/UsersGuidech9.html) - `C1.1.1 Thorns`, `C1.2 Anatomy of a Thorn`, `make.code.defn based thorn building`; accessed 07-12-2026

## See Also

- [ETLegacy](index.md)
- [Code Parameters, Includes, And Loops](code-parameters-includes-and-loops.md)
- [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
- [C Function Registry](../../core/c-function-registry.md)
