# Overview

> Map NRPy's repository layout, packaging metadata, and symbolic-to-generated-project flow. · Status: confirmed
> Up: [Architecture](index.md)

## Summary

NRPy is a Python/SymPy symbolic code-generation library (`nrpy`) for numerical relativity and relativistic astrophysics. Runnable scripts in `nrpy/examples/` combine symbolic equations (`nrpy/equations/`), shared symbolic/loop utilities (`nrpy/helpers/`), C/Python code generation (`c_codegen`, `py_codegen`), function registries (`CFunction_dict`, `PyFunction_dict`), and target backends (`nrpy/infrastructures/`) to emit complete solver projects or Cactus thorns under `project/`.

## Detail

### Repository Layout And Packaging Metadata

`nrpy/examples/` holds runnable project generators rather than handwritten C/C++ solver binaries:
- `BHaH` emits standalone C/CUDA executables or libraries (`project/<name>/`), `ETLegacy` and `CarpetX` emit Einstein Toolkit Cactus thorns, `superB` emits Charm++ chare projects, `Dendro` emits standalone MPI octree CMake applications ([Dendro](../infrastructures/dendro/index.md)), and `JAX` emits installable Python/JAX packages ([Generated Backend Comparison](../generated-backend-comparison.md)).
- In `setup.py`, `nrpy` reads runtime dependencies from `requirements.txt`, extracts the first `version = ...` line from `nrpy/release.txt` via `get_nrpy_version()` (raising if absent; `release.txt` itself is not added to `package_data`), assigns non-test `*.h` headers to their nearest Python subpackage via `discover_header_package_data()`, and includes `nrpy/py.typed`.
- `CITATION.md` lists the general SENR/NRPy+ paper alongside project papers for `NRPyElliptic`, GPU/CUDA relaxation, `superB`, `BHaHAHA`, `GRoovy`, and `SEOBNR`/`SEBOB` waveform generators.

### Symbolic-To-Project Generation Sequence

A standalone generator (such as `nrpy/examples/wave_equation_cartesian.py`) executes four steps:
1. **Configuration**: Set `Infrastructure` (e.g., `"BHaH"`), `project_name`, finite-difference order, and floating-point type (`REAL`), and clear `project/<project_name>/` before registering functions.
2. **Symbolic expressions to C (`ccg.c_codegen`)**: Construct SymPy scalar or indexed expressions from `nrpy/equations/`, then lower them to C assignments via `ccg.c_codegen` (`CCodeGen` controls CSE, `enable_fd_codegen`, `enable_simd`, floating-point type aliases, and gridfunction array reads, often wrapped in `BHaH.simple_loop.simple_loop`).
3. **Function registration (`cfc.register_CFunction`)**: Store each `CFunction` (subdirectory, includes, `prefunc`, Doxygen description, `cfunc_type`, `name`, `params`, `body`, `postfunc`) in `cfc.CFunction_dict`, rejecting duplicate names. `CFunction.generate_full_function` assembles includes, Doxygen header, prototype, optional `set_CodeParameters.h`, and body.
4. **Backend project assembly**: Register timestepping (MoL), boundary conditions, and diagnostics; write `CodeParameter` headers, default `.par`/`.toml` files, and `BHaH_defines.h`; register `main.c` (`BHaH.main_c.register_CFunction_main_c` verifies required kernels exist before emitting CLI parsing, grid/timestep setup, allocation, initial data, evolution/diagnostics, and cleanup); copy SIMD/CUDA headers when enabled; and emit `BHaH_function_prototypes.h` plus the `Makefile` or `CMakeLists.txt`.

## Sources

- [README.md](../../README.md) - `# NRPy`, `## Project Families and Example Generators`, `## Repository Map`, `## What Gets Generated?`
- [CITATION.md](../../CITATION.md) - general NRPy and project-specific paper citations
- [setup.py](../../setup.py) - `get_nrpy_version`, `discover_header_package_data`, `setup(...)`
- [wave_equation_cartesian.py](../../nrpy/examples/wave_equation_cartesian.py) - `project_name`, `register_CFunction_rhs_eval`, `BHaH.main_c.register_CFunction_main_c`
- [c_function.py](../../nrpy/c_function.py) - `CFunction`, `CFunction_dict`, `register_CFunction`
- [c_codegen.py](../../nrpy/c_codegen.py) - `CCodeGen`, `c_codegen`
- [main_c.py](../../nrpy/infrastructures/BHaH/main_c.py) - `register_CFunction_main_c`

## See Also

- Depends on: [Generated Output Boundaries](generated-output-boundaries.md)
- See also: [Build And Run](build-and-run.md)
- See also: [CFunction Registry](../core/c-function-registry.md)
