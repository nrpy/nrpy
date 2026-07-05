# Generated Output Boundaries

> Distinguish handwritten source from generated projects, generated thorns, generated language targets, and artifacts. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [Architecture](index.md)

## Summary

NRPy source files are the handwritten generators, symbolic equations, infrastructure modules, helpers, packaging files, and CI configuration. Generated C/CUDA projects, Einstein Toolkit thorns, Charm++ projects, JAX projects, compiled binaries, archives, images, and scratch outputs are products or artifacts unless a selected generated file is deliberately registered as frozen evidence.

## Detail

Most example generators write under `project/<name>/`. Standalone BHaH examples usually create a buildable project containing C source, headers, a Makefile, a parameter file, and an executable target for Cartesian or curvilinear-coordinate workflows. Waveform, geodesic, and some physics generators also produce standalone projects, often with a GSL build dependency.

Einstein Toolkit and CarpetX generators produce thorn directories rather than standalone executables. Those thorns are generated into `project/<name>/` and then copied or linked into an Einstein Toolkit checkout for build and run. `superB` generators produce Charm++ projects, and `sebobv1_jax` generates a Python/JAX project instead of a C executable.

CI deliberately creates and builds generated outputs as validation evidence. The Ubuntu and macOS codegen jobs install NRPy, generate many projects, and run `make` where applicable. Separate Einstein Toolkit and Charm++ jobs generate thorns or Charm++ projects inside CI containers and validate that the generated products build or run.

The contribution boundary is stricter than the generated-project capability. Binary files, images, archives, compiled artifacts, generated projects, and similar non-text outputs should not be committed unless maintainers approve or a specific generated artifact is intentionally registered as source evidence. For documentation and review, cite the generator source and stable symbols rather than transient generated output whenever possible.

## Sources

- [README.md](../../README.md) - `## What Gets Generated?`, `## Project Families and Example Generators`
- [.github/workflows/main.yml](../../.github/workflows/main.yml) - `einsteintoolkit-validation`, `charmpp-validation`, `codegen-ubuntu`, `codegen-mac`
- [raw/source-docs/original-agents.md](../../raw/source-docs/original-agents.md) - `## Scope`
- [coding_style.md](../../coding_style.md) - `## Python Coding Style`, `### Formatting`

## See Also

- Parent: [Architecture](index.md)
- Depends on: [Build And Run](build-and-run.md)
- Depends on: [Contribution Style And Static Analysis](contribution-style-and-static-analysis.md)
- See also: [Overview](overview.md)
