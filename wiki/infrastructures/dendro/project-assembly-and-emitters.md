# Project Assembly And Emitters

> Explain how the Dendro infrastructure turns the NRPy registries into a complete generated solver directory, which module emits which artifact, and where the emitted names come from. · Status: provisional · Last reconciled: 09-05-2026
> Up: [Dendro](index.md)

## Summary

A Dendro run has two halves. Builders register gridfunctions, CodeParameters,
and CFunctions into NRPy's registries; then `output_project` walks a fixed map
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
| `Dendro_types_h` | the scalar contract header and the projection status record |
| `Dendro_state_h` | the EVOL enum, name array, metadata, exact-name lookup, and the constants header |
| `CodeParameters` | the generated parameter struct header and the sample parameter table |
| `Dendro_preamble_h` | the preamble every generated source includes |
| `output_CFunctions` | one source file per registered CFunction, the declaration header, and the CMake source list |
| `Dendro_solver_context` | the host context header and source |
| `Dendro_main_cpp` | the entry point and its lifecycle gates |
| `Dendro_self_tests_cpp` | the generated CTest sources |
| `Dendro_parameter_file` | the sample parameter file for one profile |
| `cmake_helpers` | the solver and tests `CMakeLists.txt` |
| `Dendro_README_md` | the project and solver READMEs |

`output_project` owns no formulation choice and holds no state: it maps
emitter output onto project-relative paths and writes it, then copies the mock
host header through `nrpy.helpers.generic.copy_files`, exactly as BHaH copies
`simd_intrinsics.h`.

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
- Role: normative rule
- Deciding authority: `nrpy/infrastructures/Dendro/output_project.py`, `output_project` signature
- Corroboration: `nrpy/examples/dendro_fccz4.py`, the module-level `solver_name`/`solver_stem`/`solver_namespace` assignments

### Emitted layout

The solver is emitted at `Dendro-GR/<solver_name>/` inside the project
directory: `generated/include` and `generated/src` hold the registry-derived
artifacts, `include/` and `src/` the host context and entry point, `pars/` the
sample parameter file, `tests/` the generated self-tests, and `host_mock/` the
mock host header the solver compiles against. A top-level `README.md` sits
beside `Dendro-GR/`.

### Source list and padding

`output_CFunctions.CFunction_cmake_source_list` derives the CMake source list
from `cfc.CFunction_dict`, so the build can never carry a hand-written source
inventory: one registered CFunction, one emitted source file, one CMake entry.

The ghost points the emitted kernels need are recorded by the right-hand-side
builder through `registration.set_required_padding` and read back by
`output_project` through `registration.required_padding`. They are not
`fd_order // 2`: the upwinded and Kreiss-Oliger operator families reach one
point further than the centred ones.

### Determinism

Regenerating an unchanged environment in a fresh process reproduces the tree
byte for byte. Nothing stamps a timestamp, an absolute path, or a hash into an
emitted file; the artifact map is written in sorted path order and every
registry read is order-stable because `GridFunction.gridfunction_lists()` sorts
case-insensitively.

Claim evidence:
- Claim: two fresh-process runs of the example generator with the same arguments produce byte-identical project trees.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/output_project.py`, `output_project`
- Corroboration: `nrpy/grid.py`, `GridFunction.gridfunction_lists` case-insensitive sort
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3, GCC 13.3.0, CMake 3.28.3; backend=Dendro; precision=double; GPU=not-applicable; restart=not-applicable; distributed=1 and 2 MPI ranks; error_path=not-run; options=--fd-order 4 --no-ko; date=09-05-2026`

### Artifact boundary

The emitted solver and its binaries are generated products, not source
evidence. Cite the Python emitters and the registry symbols instead; see
[Generated Output Boundaries](../../architecture/generated-output-boundaries.md).

## Sources

- [output_project.py](../../../nrpy/infrastructures/Dendro/output_project.py) - `output_project`
- [output_CFunctions.py](../../../nrpy/infrastructures/Dendro/output_CFunctions.py) - `CFunction_artifacts`, `CFunction_cmake_source_list`, `derived_source_path`
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - `output_Dendro_parameters_h`, `output_parameter_file_sample`
- [Dendro_state_h.py](../../../nrpy/infrastructures/Dendro/Dendro_state_h.py) - `state_records`, `output_Dendro_state_h`, `output_Dendro_constants_h`
- [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py) - `output_solver_cmake`, `output_generated_sources_cmake`, `output_tests_cmake`
- [Dendro_solver_context.py](../../../nrpy/infrastructures/Dendro/Dendro_solver_context.py) - `output_Dendro_solver_context_h`, `substitute_solver_identifiers`
- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - generation entry point and command-line profile

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Gridfunctions, Naming, And Loops](gridfunctions-naming-and-loops.md)
- Validated by: [Validation, Host Mock, And Deferral Gates](validation-host-mock-and-deferral-gates.md)
- Contrasts with: [superB Lifecycle And Project Assembly](../superb/lifecycle-and-project-assembly.md)
- See also: [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
