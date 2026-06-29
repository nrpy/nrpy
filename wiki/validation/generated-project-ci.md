# Generated Project CI

> CI coverage for generated projects, external backend validation, and waveform consistency checks. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Validation](index.md)

## Summary

NRPy CI does more than lint Python. It generates and builds standalone projects
on Ubuntu and macOS, validates generated Einstein Toolkit thorns, validates
Charm++/`superB` workflows, and compares current SEOB/SEBOB waveform outputs
against a trusted commit with on-the-fly accuracy checks.

## Detail

The `codegen-ubuntu` and `codegen-mac` jobs install the package, generate many
example projects under a temporary workspace, and run `make` plus `make clean`
for buildable C projects. The example set includes wave-equation projects,
elliptic and black-hole workflows, SEOBNRv5 variants, TOV and hydro examples,
the apparent-horizon library generator, `sebobv2`, and `sebobv1_jax`. Ubuntu
installs GSL through `apt`; macOS installs it through Homebrew.

The `einsteintoolkit-validation` job uses an Apptainer image to generate Carpet
WaveToy and Baikal thorns, links the generated thorns into an Einstein Toolkit
checkout, builds the toolkit, and runs the Baikal, BaikalVacuum, and WaveToyNRPy
testsuites. The job fails if the Einstein Toolkit summary reports any failed
tests.

The `charmpp-validation` job uses an Apptainer image with Charm++ to generate
`superB` elliptic, black-hole spectroscopy, and two-black-hole workflows. It
builds each generated project and runs `superB_two_blackholes_collide` through
`charmrun +p2`.

The SEOB/SEBOB jobs checkout a trusted commit, build trusted and current
executables, then run comparison scripts. `sebob_consistency_check.py` compares
the current and trusted SEOBNRv5-family outputs across multiple approximants;
`sebobv2_consistency_check.py` does the same for `sebobv2`. Both scripts compute
an RMSE-style waveform difference and require the current-vs-trusted median
error to stay within a perturbation-derived roundoff baseline.

These jobs intentionally create generated `project/` outputs. Treat those as CI
products, not committed documentation or hand-authored source, unless a selected
generated file has been deliberately registered as frozen evidence.

## Sources

- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `codegen-ubuntu`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `codegen-mac`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `einsteintoolkit-validation`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `charmpp-validation`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `sebob-consistency-test`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `sebobv2-consistency-test`
- [../../README.md](../../README.md) - `## Project Families and Example Generators`
- [../../README.md](../../README.md) - `## What Gets Generated?`
- [../../nrpy/examples/tests/sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py) - `calculate_rmse`
- [../../nrpy/examples/tests/sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py) - `process_input_set`
- [../../nrpy/examples/tests/sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py) - `calculate_rmse`
- [../../nrpy/examples/tests/sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py) - `process_input_set`

## See Also

- [Validation](index.md)
- [Static Analysis](static-analysis.md)
- [Glossary](../glossary.md)
- [Workflows](../workflows.md)
