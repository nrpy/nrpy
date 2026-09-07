# Project Assembly And Emitters

> Explain how the Dendro infrastructure turns the NRPy registries into a complete generated solver directory, which module emits which artifact, and where the emitted names come from. · Status: provisional · Last reconciled: 09-07-2026
> Up: [Dendro](index.md)

## Summary

A Dendro run has two halves. Builders register gridfunctions, CodeParameters,
and CFunctions into NRPy's registries; then `main` walks a fixed map
of project-relative paths to emitter output and writes it. Every emitter reads
`gri.glb_gridfcs_dict`, `par.glb_code_params_dict`, and `cfc.CFunction_dict`
directly at the point of use, as BHaH's emitters do. There is no snapshot
record set, no manifest, no installer, and no generation transaction.

## Detail

### One module per emitted artifact

Modules are named for what they emit, following BHaH's `BHaH_defines_h.py` and
`main_c.py`:

| Module | Emits |
| --- | --- |
| `types_h` | the scalar contract header and the det/trace enforcement status record |
| `state_h` | the EVOL enum, name array, metadata, and exact-name lookup |
| `constants_h` | the generated finite-difference order, required padding and Kreiss-Oliger switch |
| `CodeParameters` | the generated `params_struct` header, the sample parameter table, and the parameter CFunctions |
| `Dendro_defines_h` | the `<stem>_defines.h` header every generated source includes, playing the role `BHaH_defines.h` plays in BHaH |
| `cmake_helpers` | one source file per registered CFunction, the `<stem>_function_prototypes.h` header, the CMake source list, and the solver and tests `CMakeLists.txt` |
| `solver_context` | the host context header and source |
| `main_cpp` | the entry point and its lifecycle gates |
| `self_tests_cpp` | the generated CTest sources |
| `cmdline_input_and_parfiles` | the sample parameter file for one profile |
| `block_kernel_helpers` | the formulation-agnostic pointer bindings, point loop, parameter lists, operator records and padding every builder lowers through |

`main` owns no formulation choice and holds no state: it maps
emitter output onto project-relative paths and writes it, then copies the
standalone host header through `nrpy.helpers.generic.copy_files`, exactly as
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

### Emitted layout

Two examples drive this layer: `nrpy.examples.dendro_fccz4` and
`nrpy.examples.dendro_bssn`. The second one exists as the test that the layer
is generic. Adding it did require generic-layer work — the formulation-agnostic
lowering moved into `block_kernel_helpers` and `tensor_family_of` into `gridfunction_name_decorations` —
but no existing emitter changed behaviour.

The solver is emitted at `Dendro-GR/<solver_name>/` inside the project
directory: `generated/include` and `generated/src` hold the registry-derived
artifacts, `include/` and `src/` the host context and entry point, `pars/` the
sample parameter file, `tests/` the generated self-tests, and `standalone_host/`
the standalone host header the solver compiles against. The project carries no
generated README: an emitted prose file would restate what this page and the
generated `CMakeLists.txt` already carry.

### Source list and padding

`cmake_helpers.CFunction_cmake_source_list` derives the CMake source list
from `cfc.CFunction_dict`, so the build can never carry a hand-written source
inventory: one registered CFunction, one emitted source file, one CMake entry.

The ghost points the emitted kernels need are recorded by the right-hand-side
builder through `CFunction_roles.set_required_padding` and read back by
`main` through `CFunction_roles.required_padding`. They are not
`fd_order // 2`: the upwinded and Kreiss-Oliger operator families reach one
point further than the centered ones.

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
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3, GCC 13.3.0, CMake 3.28.3; backend=Dendro; precision=double; GPU=not-applicable; restart=not-applicable; distributed=1 and 2 MPI ranks; error_path=not-run; options=--fd-order 4 --no-ko; date=09-06-2026`

### Artifact boundary

The emitted solver and its binaries are generated products, not source
evidence. Cite the Python emitters and the registry symbols instead; see
[Generated Output Boundaries](../../architecture/generated-output-boundaries.md).

## Sources

- [generated_file_banner.py](../../../nrpy/infrastructures/Dendro/generated_file_banner.py) - `generated_file_banner`
- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - `main`, the inline project-assembly block and command-line profile
- [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py) - `module_layout`, `ModuleLayout`, `output_CFunctions_function_prototypes_and_construct_CMakeLists`, `CFunction_cmake_source_list`, `derived_source_path`, `output_solver_cmake`, `output_generated_sources_cmake`, `output_tests_cmake`
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - `output_parameters_h`
- [state_h.py](../../../nrpy/infrastructures/Dendro/state_h.py) - `state_records`, `output_state_h`
- [constants_h.py](../../../nrpy/infrastructures/Dendro/constants_h.py) - `output_constants_h`
- [types_h.py](../../../nrpy/infrastructures/Dendro/types_h.py) - `output_types_h`
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - `output_solver_context_h`, `substitute_solver_identifiers`

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Gridfunctions, Naming, And Loops](gridfunctions-naming-and-loops.md)
- Validated by: [Validation, Standalone Host, And Deferral Gates](validation-standalone-host-and-deferral-gates.md)
- Contrasts with: [superB Lifecycle And Project Assembly](../superb/lifecycle-and-project-assembly.md)
- See also: [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
