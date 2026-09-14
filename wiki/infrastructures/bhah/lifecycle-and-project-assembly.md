# Lifecycle And Project Assembly

> Explain how BHaH standalone applications register runtime functions, assemble generated projects, and split executable and library entrypoints. · Status: confirmed
> Up: [BHaH](index.md)

## Summary

BHaH is NRPy's standalone generated-application infrastructure. Generators
register C functions into the global `CFunction` registry, assemble generated
source, headers, prototypes, and a dependency-aware Makefile under the generated
project directory, then build an executable lifecycle around
`register_CFunction_main_c` or a library-facing API around
`register_CFunctions_bhah_lib`.

## Detail

`register_CFunction_main_c` is the final assembler for normal BHaH
executables. Before registering generated `main`, it checks that required
runtime pieces already exist in `CFunction_dict`: parameter/default setters,
the command-line and parameter-file parser, numerical grid and timestep setup,
Method of Lines allocation/free/step functions, diagnostics, and initial data.
Detailed parameter, header, default-parfile, and CLI behavior belongs in
[Runtime Data, Parameters, Headers, And CLI](runtime-data-parameters-headers-and-cli.md).

The generated `main` owns the high-level standalone lifecycle. It initializes
`commondata`, parses inputs, optionally reads checkpoint metadata before full
grid allocation, allocates `griddata` for `MAXNUMGRIDS`, defaults grid-local
parameters, calls `numerical_grids_and_timestep`, allocates evolved and
auxiliary gridfunction storage, sets initial data, and enters
`while(commondata.time < commondata.t_final)`. Each iteration can run injected
pre-diagnostics code, calls `diagnostics`, can run injected pre-step code, calls
`MoL_step_forward_in_time`, and can run injected post-step code. Cleanup frees
griddata state, with CUDA and BHaHAHA-specific branches when those functions
are registered.

`GridCommonData` and `register_griddata_commondata` provide the project-wide
struct extension mechanism consumed by generated BHaH headers. Registration
stores declarations under `par.glb_extras_dict["griddata_struct"]` for
grid-local data or `par.glb_extras_dict["commondata_struct"]` for data common
to all grids, rejects duplicate declarations for the same module, and keeps a
description with each declaration. `register_CFunction_griddata_free` emits the
host cleanup routine for per-grid allocations; in CUDA mode it also registers
the device cleanup companion.

`register_CFunctions_bhah_lib` registers a smaller library surface:
`bhah_initialize`, `bhah_evolve`, `bhah_diagnostics`, and `bhah_finalize`.
The supplemental defines introduce `BHaH_struct`, which wraps pointers to
`commondata_struct` and `griddata_struct`. This library path initializes the
same core runtime state, advances while `commondata->time` is below
`commondata->t_final`, dispatches diagnostics, and finalizes allocated grid and
reference-metric data
through library-call entrypoints instead of a generated process `main`.

`nrpy.examples.manga_bhah_lib` selects this library boundary explicitly:
`exec_or_library_name="libbhah_lib"` and `create_lib=True` produce
`libbhah_lib.so` on Linux or `libbhah_lib.dylib` on Darwin. Its final generic
prints still say to run `./bhah_lib` and find parameters in `bhah_lib.par`.
Those messages are stale: this route has no executable entrypoint, its parser
registration is commented out, and it does not call the default-parfile writer.

Claim status: stale; contradiction: CONTR-0001.
See [CONTR-0001](../../contradictions.md#contr-0001) for authority, affected
pages, validation limits, and the executable resolution test; that contradiction
record owns the claim evidence.

`output_CFunctions_function_prototypes_and_construct_Makefile` turns the
registered `CFunction_dict` into a buildable generated project. It validates
library/executable options, writes registered C functions and
`BHaH_function_prototypes.h`, chooses compiler flags, includes additional
directories and libraries, and emits a Makefile whose target is either an
executable, a shared library, or a static archive. Executable builds may recurse
into additional generated subdirectories before the final link.

The generated Makefile represents each registered C-function source exactly
once in an explicit `ADD_SOURCE` record. That record pairs the source path with
each registered include that resolves inside the generated project at
generation time; `SOURCES`, `OBJECTS`, and `DEPFILES` are then derived from the
records without `find` or wildcard source discovery. A deliberate hand edit
therefore adds or removes a source in one record, but the Makefile remains
generated output and regeneration replaces such edits.

Each `ADD_SOURCE` record makes its resolved local headers immediate object
prerequisites. The C, C++, and CUDA compile rules also use `-MMD`, `-MP`,
`-MF`, and `-MT` to write dependency files under `.deps/` and load them on later
Make runs. After a successful compilation, those files track the active
non-system direct, transitive, and conditional includes seen by the compiler,
so changing one rebuilds the affected objects. System headers and inactive
conditional branches are not recorded.

Generated `clean` is visible and bounded. Its one root cleanup recipe is not
silenced, so Make prints one expanded `rm` command. That command removes the
exact final target and root `.deps/`; `.o` and `.d` files at depths zero through
two; and root `.txt`, `.gp`, `.dat`, `.out`, `.log`, `.avi`, `.png`, and `.bin`
files. It then invokes each registered subproject's `clean` with silent,
no-directory-printing Make flags. It preserves `.par`, unrelated suffixes,
runtime-suffix files below the root, and `.o`/`.d` files at depth three or
deeper. Cleanup uses explicit globs rather than `find`, recursive wildcards, a
filename manifest, or source-discovery logic.

Claim evidence:
- Claim: Generated BHaH Makefiles use one explicit `ADD_SOURCE` record per registered C-function source, derive object and dependency inventories without source-discovery commands, make resolved registered project-local headers immediate prerequisites, and use compiler dependency files under `.deps/`; generated `clean` prints one root `rm` command while removing the exact target, root dependency directory, `.o`/`.d` files only through depth two, and seven runtime suffix families only at the project root, delegates silent cleanup to registered subprojects, and preserves `.par`, unrelated suffixes, nested runtime files, and depth-three-or-deeper `.o`/`.d` files.
- Role: descriptive behavior
- Deciding authority: [Makefile_helpers.py](../../../nrpy/infrastructures/BHaH/Makefile_helpers.py), `_generate_c_files_and_header`, `_construct_makefile_content`, and `output_CFunctions_function_prototypes_and_construct_Makefile`
- Corroboration: emitted-Makefile assertions in the same owner module

Compiler selection replaces GNU Make's built-in `cc` only when that default is
active, preserving environment and command-line choices; CUDA Makefiles select
NVCC unless the command line overrides `CC`. Preprocessor, C, C++, CUDA,
link-driver, and library options remain in separate composed variables. On CPU,
`OPENMP=1` directly adds `-fopenmp` to C, C++, and link flags; an unsupported
compiler or linker fails normally. `OPENMP=0` explicitly omits the link flag and
uses unknown-pragma warning suppression for C and C++. No capability probe or
fallback is emitted. The linker is selected from live source records: host C++ sources
select `CXX`, ordinary CPU C sources select `CC`, and CUDA projects select
NVCC. Additional generated projects are rechecked before dependent object
builds and the final link. Static archives are removed before recreation so
deleted object members cannot survive.

Claim evidence:
- Claim: Generated BHaH Makefiles preserve CPU origin-aware compiler selection and CUDA command-line `CC` overrides, separate build-flag roles, apply direct CPU `-fopenmp` compile/link flags for `OPENMP=1` with no probe or fallback, retain explicit `OPENMP=0` opt-out behavior, choose the linker from live source records, recheck additional projects before dependent builds, and recreate static archives without stale members.
- Role: descriptive behavior
- Deciding authority: [Makefile_helpers.py](../../../nrpy/infrastructures/BHaH/Makefile_helpers.py), `_generate_c_files_and_header`, `_construct_makefile_content`, and `output_CFunctions_function_prototypes_and_construct_Makefile`
- Corroboration: none available; emitted-Makefile assertions live in the same owner module

The generated `make valgrind` target is executable-oriented. For a CPU
executable it cleans, rebuilds with the debug C flags and `OPENMP=0`, then runs
Valgrind with one OpenMP thread and a nonzero error exit. For a CUDA executable,
the same target name uses Compute Sanitizer instead: before cleaning, it checks
for `libsanitizer-collection.so` under `CUDA_SANITIZER_DIR`, whose default is
the Debian toolkit path `/usr/lib/nvidia-cuda-toolkit/compute-sanitizer`. A
missing library error supplies a `find` command and the
`make valgrind CUDA_SANITIZER_DIR=/path/to/compute-sanitizer` override. When the
library exists, the target cleans, adds `-lineinfo` to `NVCCFLAGS`, rebuilds,
and runs Compute Sanitizer memcheck with the explicit injection path.
Consequently `-lineinfo` reaches `.cu` compilation and the NVCC final link, but
mixed host `.cc`, `.cpp`, or `.cxx` compilation continues to use `CXXFLAGS`.
The pre-clean guard checks only the injection library, not the
`compute-sanitizer` executable, so a missing executable can still fail after
cleaning and rebuilding. CUDA library targets reject `make valgrind` because
they have no executable harness; CPU library targets only perform the debug
rebuild and do not invoke Valgrind.

Claim evidence:
- Claim: Generated CPU executable `make valgrind` rebuilds without OpenMP and runs Valgrind, while generated CUDA executable `make valgrind` checks only the configurable injection-library directory before cleaning, adds `-lineinfo` to NVCC `.cu` compilation and final linking, and runs Compute Sanitizer; CUDA libraries reject the target for lack of an executable harness, and CPU libraries rebuild without running Valgrind.
- Role: descriptive behavior
- Deciding authority: [Makefile_helpers.py](../../../nrpy/infrastructures/BHaH/Makefile_helpers.py), `_construct_makefile_content`
- Corroboration: none available; emitted-Makefile assertions live in the same owner module

`compile_Makefile` is the programmatic build wrapper. It autodetects a compiler
when requested, regenerates the prototype/header/Makefile assets through
`output_CFunctions_function_prototypes_and_construct_Makefile`, runs parallel
`make`, and retries once with debug flags if the expected executable does not
appear. User-facing build guidance still comes from the generated project path:
standalone BHaH examples write under the generated project directory, then
compile with `make` inside that generated directory.

Generated `project/**` files are products of these Python generators, not KB
source evidence. Documentation should cite generator modules and stable symbols
such as `register_CFunction_main_c`,
`output_CFunctions_function_prototypes_and_construct_Makefile`, and
`BHaH_struct`, unless maintainers deliberately register a generated artifact as
frozen evidence.

## Sources

- [README.md](../../../README.md) - `Standalone BHaH Generators`, `What Gets Generated?`
- [main_c.py](../../../nrpy/infrastructures/BHaH/main_c.py) - `register_CFunction_main_c`
- [bhah_lib.py](../../../nrpy/infrastructures/BHaH/bhah_lib.py) - `register_CFunctions_bhah_lib`, `BHaH_struct`
- [manga_bhah_lib.py](../../../nrpy/examples/manga_bhah_lib.py) - commented parser registration, `exec_or_library_name="libbhah_lib"`, `create_lib=True`, final build/run/parfile prints
- [Makefile_helpers.py](../../../nrpy/infrastructures/BHaH/Makefile_helpers.py) - `output_CFunctions_function_prototypes_and_construct_Makefile`, `compile_Makefile`
- [griddata_commondata.py](../../../nrpy/infrastructures/BHaH/griddata_commondata.py) - `GridCommonData`, `register_griddata_commondata`, `register_CFunction_griddata_free`

## See Also

- Parent: [BHaH](index.md)
- Depends on: [Runtime Data, Parameters, Headers, And CLI](runtime-data-parameters-headers-and-cli.md)
- Depends on: [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
- See also: [Build And Run](../../architecture/build-and-run.md)
- Depends on: [C Function Registry](../../core/c-function-registry.md)
- Depends on: [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)
