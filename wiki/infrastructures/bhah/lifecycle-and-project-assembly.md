# Lifecycle And Project Assembly

> Explain how BHaH standalone applications register runtime functions, assemble generated projects, and split executable and library entrypoints. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [BHaH](index.md)

## Summary

BHaH is NRPy's standalone generated-application infrastructure. Generators
register C functions into the global `CFunction` registry, assemble generated
source, headers, prototypes, and a Makefile under the generated project
directory, then build
an executable lifecycle around `register_CFunction_main_c` or a library-facing
API around `register_CFunctions_bhah_lib`.

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
pages, validation limits, and the executable resolution test. This is P1
descriptive contradiction record.

`output_CFunctions_function_prototypes_and_construct_Makefile` turns the
registered `CFunction_dict` into a buildable generated project. It validates
library/executable options, writes registered C functions and
`BHaH_function_prototypes.h`, chooses compiler flags, includes additional
directories and libraries, and emits a Makefile whose target is either an
executable, a shared library, or a static archive. Executable builds may recurse
into additional generated subdirectories before the final link. The helper also
emits `valgrind` and `clean` targets.

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

- [BHaH](index.md)
- [Runtime Data, Parameters, Headers, And CLI](runtime-data-parameters-headers-and-cli.md)
- [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
- [Build And Run](../../architecture/build-and-run.md)
- [C Function Registry](../../core/c-function-registry.md)
- [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)
