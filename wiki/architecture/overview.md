# Overview

> Map the NRPy project shape from symbolic source to generated artifacts. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Architecture](index.md)

## Summary

NRPy is a Python/SymPy symbolic code generation toolkit for numerical relativity and relativistic astrophysics. Its runnable example generators assemble symbolic equations, infrastructure code generation, and helper utilities into generated projects such as standalone BHaH applications, Einstein Toolkit thorns, Charm++ projects, and JAX/Python workflows.

## Detail

The repository is organized around generators rather than handwritten end-user executables. `nrpy/examples` contains runnable project entry points. Those examples pull symbolic physics from `nrpy/equations`, backend-specific project assembly from `nrpy/infrastructures`, and shared support from `nrpy/helpers`.

The top-level flow is: choose or write an example generator, assemble symbolic expressions and C functions, select an infrastructure target, then emit a project under a generated-output location such as `project/<name>/`. BHaH is the standalone application infrastructure and is the closest public path to single-patch BH@H numerical relativity workflows. ETLegacy and CarpetX produce Einstein Toolkit thorns, `superB` produces Charm++-based projects, and SEBOB/JAX workflows cover compact-object and Python/JAX project generation.

Packaging metadata presents the same identity: the installable package is named `nrpy`, describes itself as Python/SymPy symbolic code generation for numerical relativity and relativistic astrophysics, reads runtime requirements from `requirements.txt`, and includes package header data needed by generated projects.

## Sources

- [README.md](../../README.md) - `# NRPy`, `## Project Families and Example Generators`, `## Repository Map`
- [setup.py](../../setup.py) - `setup(...)`, `discover_header_package_data`

## See Also

- [Architecture](index.md)
- [Build And Run](build-and-run.md)
- [Symbolic Codegen Lifecycle](symbolic-codegen-lifecycle.md)
- [Generated Output Boundaries](generated-output-boundaries.md)
