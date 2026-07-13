# Generated Project CI

> CI coverage for generated projects, external backend validation, and waveform consistency checks. · Status: confirmed · Last reconciled: 07-13-2026
> Up: [Validation](index.md)

## Summary

Workflow YAML configures seven jobs. Two generate/build standalone projects on
Ubuntu and macOS; one builds and regression-tests ETLegacy thorns; one builds
three Charm++/superB projects and runs one; two build/run trusted/current
waveforms and compare output; one performs Python static analysis. These are
configured cells, not latest-pass claims.

## Detail

Configured GitHub job map:

| Job | Matrix/environment | Generate/build scope | Run/result-check scope |
| --- | --- | --- | --- |
| `static-analysis` | Ubuntu 22.04/24.04 and selected Python 3.7.13, 3.8.12, 3.9.19, `3.x` cells | No generated project coverage | Python file execution plus version-dependent static checks; see [Static Analysis](static-analysis.md) |
| `codegen-ubuntu` | Ubuntu 22.04/24.04; three matrix cells after one exclusion | Installs NRPy, generates in `tmp/`, and builds 21 default C/library projects: standalone elliptic; three wave projects; collision, spectroscopy, and spinning BH; PN momenta; all nine SEOBNRv5 approximant/calibration variants; TOVola; hydro-without-hydro; BHaHAHA; `sebobv2`. It generates `sebobv1_jax` without package install/build. | No generated executable, test, or numerical result is run; `make clean` follows each C/library build. MANGA commands are commented out. |
| `codegen-mac` | macOS 14/26 with Python 3.9, 3.10, 3.11, `3.x` | Same 21 default C/library builds and JAX generation as Ubuntu; GSL installed with Homebrew | No generated executable, test, or numerical result is run. |
| `einsteintoolkit-validation` | Ubuntu 24.04, Apptainer 1.3.2, ET 2024-06 beta image | Generates only `carpet_wavetoy_thorns.py` and `carpet_baikal_thorns.py`, links ETLegacy thorns/fixtures into ET, then builds ET | Runs `Baikal`, `BaikalVacuum`, and `WaveToyNRPy` Cactus testsuites and fails on nonzero reported failures. No `carpetx_*` generation/build/run. |
| `charmpp-validation` | Ubuntu 24.04, Apptainer image, paths pinned to Charm++ 8.0.0 | Generates and builds `superB_nrpyelliptic_conformally_flat`, `superB_blackhole_spectroscopy`, and `superB_two_blackholes_collide` with `make -j2` | Runs only `./charmrun +p2 ./superB_two_blackholes_collide`; no explicit scientific-output assertion beyond process success. |
| `sebob-consistency-test` | Ubuntu 22.04/24.04; three matrix cells after one exclusion | Checks out trusted commit `785467615d63669a98fe85c6686c2388a324139e`; generates/builds trusted and current copies of all nine SEOBNRv5 variants | Runs nine helper invocations. Each rebuilds both executables, runs ten deterministic inputs, and requires median current/trusted amplitude-plus-phase error not exceed its perturbation-derived baseline. |
| `sebobv2-consistency-test` | Same Ubuntu matrix shape | Generates/builds trusted and current `sebobv2` at the same trusted commit | Runs one helper invocation with ten deterministic inputs and the same median-error criterion. |

Einstein Toolkit documents Cactus testsuites as regression comparisons, not
convergence or physics-correctness tests; see official
[Adding a test case](https://docs.einsteintoolkit.org/et-docs/Adding_a_test_case),
heading `A test case is...`. Therefore the ET job establishes configured
ETLegacy regression checks under its image, not broad numerical validation.

Official Charm++ [Quickstart](https://charm.readthedocs.io/en/v8.0.0/quickstart.html),
headings `Compiling the Example` and `Running the Example`, supports the
external `charmc`/`charmrun +pN` toolchain shape. NRPy's workflow alone decides
its pinned Charm++ path and exact tested projects. It does not runtime-test
restart, Psi4 output, elliptic residual-stop behavior, active load balancing,
TRAM, sections, priority, immediate, expedited, or zero-copy paths.

Waveform comparison helpers directly parse generated executable stdout and
compute their own amplitude-plus-phase regression metric. The Cactus testsuite
route instead delegates numerical comparison to its fixture and tolerance
configuration; the checked-in `WaveToyNRPy` test sets `RELTOL 1e-11`. NRPy's
workflow parses the testsuite summary and fails on a nonzero failure count.
Neither regression route proves physical accuracy beyond its stated fixtures or
inputs, and workflow configuration does not prove the jobs most recently passed.

The local `.github/full_nrpy_local_ci.sh` helper is separate from GitHub job
coverage. It installs dependencies, performs broad static analysis, invokes 28
configured generator commands, and builds 21 non-Carpet/non-superB C/library
projects. It generates two superB projects without building them and all four
Carpet/CarpetX families without ET compilation; it omits the superB elliptic,
Kasner, GRoovy, MANGA, and all geodesic generators. It then configures builds
with `--cuda` for curvilinear wave, multicoordinate wave, standalone elliptic,
three black-hole examples, hydro-without-hydro, and TOVola. TOVola has no
argument parser or CUDA branch, so its extra `--cuda` token is ignored and that
cell is an ordinary C build. The helper installs no CUDA toolkit, declares no
GPU runner, runs no generated executable, and checks no GPU result. Treat it as
a local command recipe requiring a prepared environment, not CI pass evidence.

Explicitly unsupported or unverified by these configurations: CarpetX build or
runtime; JAX generated-package install/import/basic test or accelerator runtime;
any CUDA executable/GPU result; standalone evolution runtime; restart behavior;
geodesic/raytracing projects; GRoovy; active MANGA build; Kasner; and scientific
correctness beyond stated ET regression and waveform-comparison assertions.

These jobs intentionally create generated `project/` outputs. Treat those as CI
products, not committed documentation or hand-authored source, unless a selected
generated file has been deliberately registered as frozen evidence.

## Sources

- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `codegen-ubuntu`; [.github/full_nrpy_local_ci.sh](../../.github/full_nrpy_local_ci.sh) - `example_scripts`, `cuda_example_scripts`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `codegen-mac`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `einsteintoolkit-validation`; official Einstein Toolkit [Adding a test case](https://docs.einsteintoolkit.org/et-docs/Adding_a_test_case) - `A test case is...`
- [test.ccl](../../nrpy/examples/et_WaveToyfiles/test/test.ccl) - `TEST WaveToyNRPy_test`, `RELTOL 1e-11`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `charmpp-validation`; official Charm++ [Quickstart](https://charm.readthedocs.io/en/v8.0.0/quickstart.html) - `Compiling the Example`, `Running the Example`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `sebob-consistency-test`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `sebobv2-consistency-test`
- [../../README.md](../../README.md) - `## Project Families and Example Generators`
- [../../README.md](../../README.md) - `## What Gets Generated?`
- [../../nrpy/examples/tests/sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py) - `calculate_rmse`
- [../../nrpy/examples/tests/sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py) - `process_input_set`
- [../../nrpy/examples/tests/sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py) - `calculate_rmse`
- [../../nrpy/examples/tests/sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py) - `process_input_set`

## See Also

- Parent: [Validation](index.md)
- Contrasts with: [Static Analysis](static-analysis.md)
- See also: [Glossary](../glossary.md)
- See also: [Workflows](../workflows.md)
