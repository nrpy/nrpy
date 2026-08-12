# Symbolic Codegen Lifecycle

> Explain the path from symbolic expressions to registered C functions and generated projects. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [Architecture](index.md)

## Summary

NRPy examples configure generation parameters, build symbolic expressions, turn those expressions into C with `ccg.c_codegen`, register emitted C functions through `cfc.register_CFunction`, and then ask an infrastructure layer such as BHaH to write headers, parameter files, a `main` function, and a build system.

## Detail

A generator starts by selecting infrastructure and code-generation settings. The wave-equation Cartesian example sets `Infrastructure` to `BHaH`, declares a `project_name`, configures finite-difference order and floating-point type, and clears the generated project directory before registering any C functions.

Core generated functions are registered into `cfc.CFunction_dict`. The `CFunction` class stores the subdirectory, include list, optional helper code, description, return type, name, parameters, body, and optional post-function text. `register_CFunction` rejects duplicate function names and inserts a `CFunction` instance into the registry. `CFunction.generate_full_function` combines includes, generated Doxygen-style description text, the function prototype, optional `set_CodeParameters` include, and the body into the final C text.

Symbolic-to-C conversion happens through `ccg.c_codegen`. The `CCodeGen` configuration controls options such as CSE, finite-difference codegen, SIMD, floating-point aliases, and automatic gridfunction memory reads. In the wave-equation example, an exact single-point solution is emitted with `ccg.c_codegen`, while the right-hand-side kernel uses `ccg.c_codegen(..., enable_fd_codegen=True, enable_simd=enable_simd)` inside a BHaH simple loop.

After registering the problem-specific functions, the example registers Method of Lines timestepping and diagnostics support, writes code-parameter headers and the default parameter file, emits `BHaH_defines.h`, registers the generic BHaH `main` function, registers griddata cleanup, copies SIMD support if enabled, and writes C function prototypes plus the Makefile. `register_CFunction_main_c` checks that required functions are already present, then assembles the C lifecycle: parse inputs, set up grids and timestep, allocate gridfunctions, set initial data, run diagnostics and time stepping, and free memory.

## Sources

- [README.md](../../README.md) - `## Repository Map`, `## What Gets Generated?`
- [nrpy/examples/wave_equation_cartesian.py](../../nrpy/examples/wave_equation_cartesian.py) - `project_name`, `register_CFunction_rhs_eval`, `BHaH.main_c.register_CFunction_main_c`
- [nrpy/c_function.py](../../nrpy/c_function.py) - `CFunction`, `CFunction_dict`, `register_CFunction`
- [nrpy/c_codegen.py](../../nrpy/c_codegen.py) - `CCodeGen`, `c_codegen`
- [nrpy/infrastructures/BHaH/main_c.py](../../nrpy/infrastructures/BHaH/main_c.py) - `register_CFunction_main_c`

## See Also

- [Architecture](index.md)
- [Overview](overview.md)
- [Build And Run](build-and-run.md)
- [Generated Output Boundaries](generated-output-boundaries.md)
