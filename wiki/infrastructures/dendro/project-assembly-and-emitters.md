# Project Assembly And Generating Functions

> Explain how the Dendro infrastructure turns the NRPy registries into a complete generated solver directory, which module writes each generated file, and where the generated names come from. · Status: provisional
> Up: [Dendro](index.md)

## Summary

A Dendro run has two halves. Builders register gridfunctions, CodeParameters,
and CFunctions into NRPy's registries; then `main` walks a fixed map
of project-relative paths to generated text and writes it. Generating functions read the
registries and Dendro CFunction-role metadata at the point of use. There is no
intermediate copy of the registries or multi-step installation process.

## Detail

### One module per generated file

Modules are named for what they generate, following BHaH's `BHaH_defines_h.py` and
`main_c.py`:

| Module | Generates |
| --- | --- |
| `types_h` | scalar-type requirements plus declarations supplied explicitly by the formulation module |
| `state_h` | the EVOL enum, name array, metadata, and exact-name lookup |
| `constants_h` | the generated regular finite-difference order, KO base and effective orders, required padding, and Kreiss-Oliger switch |
| `CodeParameters` | the generated `params_struct` header and parameter CFunctions |
| `Dendro_defines_h` | the `<stem>_defines.h` header every generated source includes, playing the role `BHaH_defines.h` plays in BHaH |
| `cmake_helpers` | one source file per registered CFunction, the `<stem>_function_prototypes.h` header, the CMake source list, and the solver and tests `CMakeLists.txt` |
| `solver_context` | generic host geometry, storage, transport, reductions, and exterior traversal |
| `general_relativity/solver_context` | GR context declarations, initialization, exterior values, projection, and diagnostics |
| `main_cpp` | the process entry point, arguments, mesh construction, and time stepping |
| `general_relativity/main_cpp` | GR initialization, evolution, diagnostics, and CTest registration |
| `self_tests_cpp` | the generic test program and isolated scalar/vector numerical test |
| `general_relativity/self_tests_cpp` | GR scientific sections and independent nonflat block reference |
| `parfile` | the sample parameter file for one profile |
| `block_kernel_helpers` | the formulation-agnostic pointer bindings, point loop, parameter lists, operator records and padding every builder lowers through |

The example generator shows how the inputs are combined. It passes GR
declarations, context choices, test sections, and CTest statements into generic
generating functions.
`main` holds no hidden state: it maps
generated text onto project-relative paths and writes it, then copies the
standalone host header and shared `block_geometry.h` through `nrpy.helpers.generic.copy_files`, exactly as
BHaH copies `simd_intrinsics.h`.

### Names come from the caller, not from a parameter registry

`solver_name`, `solver_prefix`, `solver_stem`, `solver_namespace`,
`production_target`, `qualification_target`, and `profile_name` are function
arguments threaded from the example, as BHaH threads `project_name` and
ETLegacy threads `thorn_name`. None of them is a registered `CodeParameter`.

Module-level names identify NRPy as the author. The generated directories and
CMake projects are `nrpy_bssn` and `nrpy_fccz4`; their namespaces are
`nrpy::bssn` and `nrpy::fccz4`. Names inside each module retain the formulation
stem. Production libraries are `nrpy_bssn_dendro` and `nrpy_fccz4_dendro`;
qualification executables append `_qualify` to those names.

Claim evidence:
- Claim: the generated unit's names are function arguments threaded from the calling example rather than registered `CodeParameter`s; module directories, CMake projects, and namespaces identify NRPy, while child files, functions, and required Dendro targets retain their formulation names.
- Role: descriptive behavior
- Deciding authority: [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py), the generating-function arguments in `main`
- Corroboration: [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py), `module_layout`, which takes the solver name as an argument and registers no parameter for it

### Parameter selection and host geometry

The generated parameter struct is the union of non-`#define` CodeParameters
recorded beside registered CFunctions. This retains parameters needed only by
generated standalone tests, including the smooth-perturbation amplitude and
wavelength, but excludes unrelated entries left in NRPy's process-global
registry.

The real-host accepted TOML inputs are narrower. The real context forwards parameter
members only when it calls the CFunction with role `rhs_eval_block`; therefore
the sample file, TOML bindings, and effective-parameter printout contain only
that CFunction's recorded parameters which also opt in through
`add_to_parfile`. Smooth-perturbation controls are not real-host inputs because
the real-host executable does not call that standalone qualification kernel.
`name`, `fd_order`, `ko_fd_order`, `ko_effective_difference_order`,
`required_padding`, and `ko_enabled` remain meaningful profile assertions: the
parser compares them with the generated kernel profile rather than forwarding
them as physics parameters.

Claim evidence:
- Claim: the Dendro parameter struct is the registered-CFunction use closure, while its real-host TOML/sample/print interface is the `add_to_parfile` subset used by the block-RHS CFunction.
- Role: descriptive behavior
- Deciding authority: [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py), `emitted_parameter_names` and `runtime_parameter_names`; [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), `codeparameter_tail` and `_REAL_SOURCE`
- Corroboration: [CFunction_roles.py](../../../nrpy/infrastructures/Dendro/CFunction_roles.py), `set_CFunction_codeparameters` and `CFunction_name_for_role`; [parfile.py](../../../nrpy/infrastructures/Dendro/parfile.py), `output_parfile_sample`

This uses the same parameter-selection rule as [ETLegacy parameter assembly](../etlegacy/thorn-assembly-and-ccl-files.md): expose the parameters attributed to
the generated functions, not every global registration. The host adapters are
different, so their closures are different. ETLegacy constructs a thorn's
`param.ccl` from all matching thorn CFunctions; this Dendro real host currently
forwards parameters only at its block-RHS call.

The host defines domain geometry. `grid_physical_size` is NRPy's
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
this layer. The BSSN example exists as the test that the layer is generic.
Adding it did require generic-layer work — the formulation-agnostic lowering
moved into `block_kernel_helpers` and `tensor_family_of` into
`gridfunction_name_decorations`.

The solver is generated at `Dendro-GR/<solver_name>/` inside the project
directory: `generated/include` and `generated/src` hold headers and source files
derived from the registries, `include/` and `src/` the host context and entry point, `pars/` the
sample parameter file, `tests/` the generated self-tests, and `standalone_host/`
the standalone host header. Real builds use Dendrolib headers and the
shared geometry header in `include/`. The project carries no
generated README: another generated prose file would restate what this page and the
generated `CMakeLists.txt` already carry.

Each generated project contains one finite-difference profile. Generate a
separate project when an application needs another member of the 4/2, 6/4, or
8/6 regular/KO order set. This keeps every source file, constant, parameter
check, padding requirement, and compiled library in one project consistent;
the Dendro application may then select or link the generated project it needs.

### Host selection

Every generated solver uses the same two CMake modes. Configuring its directory
as the top-level project builds the standalone qualification executable and
generated numerical tests by default. Adding it to a Dendro tree that defines
`dendro5` builds the production library; drivers and tests default off in that
embedded mode. `NRPY_DENDRO_BUILD_DRIVERS` and `NRPY_DENDRO_BUILD_TESTS` can
enable them explicitly. An embedded qualification driver also requires
`toml11::toml11` and MPI. Both BSSN and fCCZ4 use actual `ot::Mesh`,
`ot::Block`, `ot::DVector`, and `ts::Ctx` types through the same generated host interface.
The two production libraries have distinct target names and can coexist in one
Dendro application. Reproduction commands and selected-host requirements live
in the [host test README](../../../nrpy/infrastructures/Dendro/tests_infra/README.md#generated-real-host-qualification).

Claim evidence:
- Claim: top-level generated builds select the standalone qualification host, while embedded builds require `dendro5` and expose separately named BSSN and fCCZ4 production libraries; optional embedded qualification drivers also require MPI and `toml11::toml11`.
- Role: descriptive behavior
- Deciding authority: [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py), `output_solver_cmake`; [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), `_REAL_HEADER`
- Corroboration: [Dendro_defines_h.py](../../../nrpy/infrastructures/Dendro/Dendro_defines_h.py), `output_Dendro_defines_h`

### Source list and padding

`cmake_helpers.CFunction_cmake_source_list` derives the CMake source list
from `cfc.CFunction_dict`, so the build source list is generated: every
registered CFunction maps to one generated source file and one CMake entry.

The ghost points the generated kernels need are recorded by the right-hand-side
builder through `CFunction_roles.set_required_padding` and read back by `main`
through `CFunction_roles.required_padding`. Dendro profiles select regular/KO
base order pairs 4/2, 6/4, and 8/6, giving exact uniform padding 2, 3, and 4.
Algebraic expressions have numerical reach zero, and a derivative restricted to
one axis is accepted; the uniform host value is the maximum canonical reach
over all axes. Dendro owns every block's interior dimensions; generated code
checks only that each padded axis can contain the required interior and halo.

### Determinism

Regenerating an unchanged environment in a fresh process reproduces the tree
byte for byte. Nothing stamps a timestamp, an absolute path, or a hash into a
generated file; the output-file map is written in sorted path order and every
registry read is order-stable because `GridFunction.gridfunction_lists()` sorts
case-insensitively.

Claim evidence:
- Claim: two fresh-process runs of the example generator with the same arguments produce byte-identical project trees.
- Role: descriptive behavior
- Deciding authority: `nrpy/examples/dendro_fccz4.py`, `main`
- Corroboration: `nrpy/grid.py`, `GridFunction.gridfunction_lists` case-insensitive sort

### Generated file boundary

The generated solver files and binaries are outputs, not source evidence. Cite the
Python generating functions and the registry symbols instead; see
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
- Validated by: [Validation, Standalone Host, And Deferred Tests](validation-standalone-host-and-deferral-gates.md)
- Contrasts with: [superB Lifecycle And Project Assembly](../superb/lifecycle-and-project-assembly.md)
- See also: [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
