# Project Assembly And Emitters

> Explain how the Dendro infrastructure turns the NRPy registries into a complete generated solver directory, which module emits which artifact, and where the emitted names come from. · Status: provisional
> Up: [Dendro](index.md)

## Summary

A Dendro run has two halves. Builders register gridfunctions, CodeParameters,
and CFunctions into NRPy's registries; then `main` walks a fixed map
of project-relative paths to emitter output and writes it. Emitters read the
registries and Dendro CFunction-role metadata at the point of use. There is no
snapshot record set, manifest, installer, or generation transaction.

## Detail

### One module per emitted artifact

Modules are named for what they emit, following BHaH's `BHaH_defines_h.py` and
`main_c.py`:

| Module | Emits |
| --- | --- |
| `types_h` | the scalar contract plus declarations supplied explicitly by the application owner |
| `state_h` | the EVOL enum, name array, metadata, and exact-name lookup |
| `constants_h` | the generated finite-difference order, required padding and Kreiss-Oliger switch |
| `CodeParameters` | the generated `params_struct` header and parameter CFunctions |
| `Dendro_defines_h` | the `<stem>_defines.h` header every generated source includes, playing the role `BHaH_defines.h` plays in BHaH |
| `cmake_helpers` | one source file per registered CFunction, the `<stem>_function_prototypes.h` header, the CMake source list, and the solver and tests `CMakeLists.txt` |
| `solver_context` | generic host geometry, storage, transport, reductions, and exterior traversal |
| `general_relativity/solver_context` | GR context declarations, initialization, exterior values, projection, and diagnostics |
| `main_cpp` | the process, argument, mesh, and time-step shell |
| `general_relativity/main_cpp` | GR lifecycle registration and application rendering |
| `self_tests_cpp` | the generic test shell and isolated scalar/vector numerical fixture |
| `general_relativity/self_tests_cpp` | GR scientific sections and independent nonflat block reference |
| `parfile` | the sample parameter file for one profile |
| `block_kernel_helpers` | the formulation-agnostic pointer bindings, point loop, parameter lists, operator records and padding every builder lowers through |

The example is the visible assembly recipe. It passes GR declarations, context
policy, test sections, and lifecycle CTest statements into generic emitters.
`main` holds no hidden state: it maps
emitter output onto project-relative paths and writes it, then copies the
standalone host header and shared `block_geometry.h` through `nrpy.helpers.generic.copy_files`, exactly as
BHaH copies `simd_intrinsics.h`.

### Names come from the caller, not from a parameter registry

`solver_name`, `solver_prefix`, `solver_stem`, `solver_namespace`,
`exec_or_library_name`, and `profile_name` are function arguments threaded from
the example, as BHaH threads `project_name` and ETLegacy threads `thorn_name`.
None of them is a registered `CodeParameter`.

The spellings follow Dendro's own vocabulary rather than NRPy's: Dendro calls
the directory a solver (`BSSN_GR`), namespaces a solver by its lowercase
formulation (`namespace bssn`), and names solver files for the formulation
(`bssnCtx.cpp`). The generated fCCZ4 solver therefore emits `FCCZ4_GR`,
`namespace fccz4::generated`, and `fccz4Ctx.cpp`.

Claim evidence:
- Claim: the generated unit's names are function arguments threaded from the calling example rather than registered `CodeParameter`s, and their spellings follow Dendro's own conventions.
- Role: descriptive behavior
- Deciding authority: [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py), the emitter call arguments in `main`
- Corroboration: [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py), `module_layout`, which takes the solver name as an argument and registers no parameter for it

### Parameter ownership and host geometry

The generated parameter struct is the union of non-`#define` CodeParameters
recorded beside registered CFunctions. This retains parameters needed only by
generated standalone tests, including the smooth-perturbation amplitude and
wavelength, but excludes unrelated entries left in NRPy's process-global
registry.

The real-host TOML surface is narrower. The real context forwards parameter
members only when it calls the CFunction with role `rhs_eval_block`; therefore
the sample file, TOML bindings, and effective-parameter printout contain only
that CFunction's recorded parameters which also opt in through
`add_to_parfile`. Smooth-perturbation controls are not real-host inputs because
the real lifecycle does not call that standalone qualification kernel.
`name`, `fd_order`, `required_padding`, and `ko_enabled` remain meaningful
profile assertions: the parser compares them with the generated kernel profile
rather than forwarding them as physics parameters.

Claim evidence:
- Claim: the Dendro parameter struct is the registered-CFunction use closure, while its real-host TOML/sample/print interface is the `add_to_parfile` subset used by the block-RHS CFunction.
- Role: descriptive behavior
- Deciding authority: [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py), `emitted_parameter_names` and `runtime_parameter_names`; [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), `_codeparameter_tail` and `_REAL_SOURCE`
- Corroboration: [CFunction_roles.py](../../../nrpy/infrastructures/Dendro/CFunction_roles.py), `set_CFunction_codeparameters` and `CFunction_name_for_role`; [parfile.py](../../../nrpy/infrastructures/Dendro/parfile.py), `output_parfile_sample`

This uses the same ownership principle as [ETLegacy parameter assembly](../etlegacy/thorn-assembly-and-ccl-files.md): expose the parameters attributed to
the generated functions, not every global registration. The host adapters are
different, so their closures are different. ETLegacy constructs a thorn's
`param.ccl` from all matching thorn CFunctions; this Dendro real host currently
forwards parameters only at its block-RHS call.

Domain geometry is host-owned. `grid_physical_size` is NRPy's
reference-metric shorthand; for Cartesian reference metrics it supplies one
symmetric half-width for all three axes. It cannot express the generated real
host's independent, asymmetric axis bounds. The host sets its `Point` bounds,
and `block_geometry()` obtains spacing and padded origins from the Dendro mesh
and block. Consequently `grid_physical_size`, `grid_hole_radius`, Cartesian
origin entries, `xmin`/`xmax` axis entries, `NUMGRIDS`, `grid_rotates`, and
`CoordSystemName` are not Dendro real-host runtime controls merely because
other NRPy setup registered them; absent CFunction use keeps them out of the
generated parameter struct.

Claim evidence:
- Claim: Dendro domain bounds and block spacing are host-owned; NRPy's Cartesian `grid_physical_size` is a symmetric reference-metric shorthand and is not wired to that host geometry.
- Role: descriptive behavior
- Deciding authority: [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py), `_REAL_MAIN`; [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), `block_geometry`; [reference_metric.py](../../../nrpy/reference_metric.py), `ReferenceMetric.cartesian_like`
- Corroboration: [Dendro-GR](https://github.com/paralab/Dendro-GR), which uses independent per-axis domain bounds and Dendro block spacing

### Emitted layout

The `nrpy.examples.dendro_fccz4` and `nrpy.examples.dendro_bssn` examples drive
this layer. The BSSN example exists as the test that the layer
is generic. Adding it did require generic-layer work — the formulation-agnostic
lowering moved into `block_kernel_helpers` and `tensor_family_of` into `gridfunction_name_decorations` —
but no existing emitter changed behavior.

The solver is emitted at `Dendro-GR/<solver_name>/` inside the project
directory: `generated/include` and `generated/src` hold the registry-derived
artifacts, `include/` and `src/` the host context and entry point, `pars/` the
sample parameter file, `tests/` the generated self-tests, and `standalone_host/`
the optional standalone host header. Real builds use Dendrolib headers and the
shared geometry header in `include/`. The project carries no
generated README: an emitted prose file would restate what this page and the
generated `CMakeLists.txt` already carry.

### Host selection

`<PREFIX>_STANDALONE_HOST=ON` retains the standalone test vehicle. With it `OFF`,
the solver must be added to a host CMake tree defining `dendro5`,
`dendro_config`, `toml11::toml11`, and `bssn_common`; the generated context uses
actual `ot::Mesh`, `ot::Block`, `ot::DVector`, and `ts::Ctx` types. The external
host include directories are system includes for generated targets, keeping
generated-source warnings distinct from diagnostics owned by the host. Duplicate
executable names fail configuration.
The fCCZ4 example can be added as `FCCZ4_GR` beside upstream `BSSN_GR`.
Reproduction commands and selected-host requirements live in the
[host test README](../../../nrpy/infrastructures/Dendro/tests_infra/README.md#generated-real-host-qualification).

Claim evidence:
- Claim: disabling the standalone host selects real Dendrolib context types, requires host CMake targets `dendro5`, `dendro_config`, `toml11::toml11`, and `bssn_common`, and treats those targets' external include directories as system includes for generated compilation; duplicate executable names fail configuration.
- Role: descriptive behavior
- Deciding authority: [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py), `output_solver_cmake`; [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), `_REAL_HEADER`
- Corroboration: [Dendro_defines_h.py](../../../nrpy/infrastructures/Dendro/Dendro_defines_h.py), `output_Dendro_defines_h`

### Source list and padding

`cmake_helpers.CFunction_cmake_source_list` derives the CMake source list
from `cfc.CFunction_dict`, so the build source list is generated: every
registered CFunction maps to one emitted source file and one CMake entry.

The ghost points the emitted kernels need are recorded by the right-hand-side
builder through `CFunction_roles.set_required_padding` and read back by
`main` through `CFunction_roles.required_padding`. They are not
`fd_order // 2`: the upwinded and Kreiss-Oliger operator families reach one
point further than the centered ones. Algebraic expressions have numerical
reach zero, and a derivative restricted to one axis is accepted; the uniform
host value is the maximum canonical reach over all axes.

### Determinism

Regenerating an unchanged environment in a fresh process reproduces the tree
byte for byte. Nothing stamps a timestamp, an absolute path, or a hash into an
emitted file; the artifact map is written in sorted path order and every
registry read is order-stable because `GridFunction.gridfunction_lists()` sorts
case-insensitively.

Claim evidence:
- Claim: two fresh-process runs of the example generator with the same arguments produce byte-identical project trees.
- Role: descriptive behavior
- Deciding authority: `nrpy/examples/dendro_fccz4.py`, `main`
- Corroboration: `nrpy/grid.py`, `GridFunction.gridfunction_lists` case-insensitive sort

### Artifact boundary

The emitted solver and its binaries are generated products, not source
evidence. Cite the Python emitters and the registry symbols instead; see
[Generated Output Boundaries](../../architecture/generated-output-boundaries.md).

## Sources

- [generated_file_banner.py](../../../nrpy/infrastructures/Dendro/generated_file_banner.py) - `generated_file_banner`
- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - `main`, the inline project-assembly block and command-line profile
- [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py) - `module_layout`, `ModuleLayout`, `output_CFunctions_function_prototypes_and_construct_CMakeLists`, `CFunction_cmake_source_list`, `derived_source_path`, `output_solver_cmake`, `output_generated_sources_cmake`, `output_tests_cmake`
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - `output_parameters_h`, `emitted_parameter_names`, `runtime_parameter_names`, `output_toml_bindings`
- [CFunction_roles.py](../../../nrpy/infrastructures/Dendro/CFunction_roles.py) - CFunction roles and CodeParameter sidecars
- [parfile.py](../../../nrpy/infrastructures/Dendro/parfile.py) - `output_parfile_sample`
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - real-host profile checks and domain bounds
- [reference_metric.py](../../../nrpy/reference_metric.py) - `ReferenceMetric.cartesian_like`
- [state_h.py](../../../nrpy/infrastructures/Dendro/state_h.py) - `state_records`, `output_state_h`
- [constants_h.py](../../../nrpy/infrastructures/Dendro/constants_h.py) - `output_constants_h`
- [types_h.py](../../../nrpy/infrastructures/Dendro/types_h.py) - `output_types_h`
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - `output_solver_context_h`, `substitute_solver_identifiers`

- [Dendro_defines_h.py](../../../nrpy/infrastructures/Dendro/Dendro_defines_h.py) - explicit host header selection
- [block_geometry.h](../../../nrpy/infrastructures/Dendro/block_geometry.h) - shared generated geometry interface
- [tests_infra/README.md](../../../nrpy/infrastructures/Dendro/tests_infra/README.md) - real-host build and run procedure

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Gridfunctions, Naming, And Loops](gridfunctions-naming-and-loops.md)
- Validated by: [Validation, Standalone Host, And Deferral Gates](validation-standalone-host-and-deferral-gates.md)
- Contrasts with: [superB Lifecycle And Project Assembly](../superb/lifecycle-and-project-assembly.md)
- See also: [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
