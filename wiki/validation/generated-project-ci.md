# Generated Project CI

> CI coverage for generated projects, external backend validation, and waveform consistency checks. · Status: confirmed
> Up: [Validation](index.md)

## Summary

Workflow YAML separates static analysis, Ubuntu/macOS code generation,
ETLegacy regression, Dendro qualification, Charm++/superB, and trusted/current
waveform consistency routes. These are configured routes, not execution-result
snapshots.

## Detail

Configured GitHub job map:

| Job | Configured context | Generate/build scope | Run/result-check scope |
| --- | --- | --- | --- |
| `static-analysis` | Configured Linux/Python matrix | No generated project coverage | Python file execution plus version-dependent static checks; see [Static Analysis](static-analysis.md) |
| `codegen-ubuntu` | Configured Ubuntu/Python matrix | Installs NRPy, generates in `tmp/`, and builds the selected default C/library projects with `make`, spanning elliptic, wave, black-hole, PN, SEOBNR, TOV, hydro, BHaHAHA, and `sebobv2` routes. It generates `sebobv1_jax` without package install/build. | The `make` builds run no generated executable, and `make clean` follows each; MANGA commands are commented out. |
| `codegen-mac` | Configured macOS/Python matrix | Same selected default C/library builds and JAX generation as Ubuntu; no Dendro generation or build; GSL installed with Homebrew | No generated executable, test, or numerical result is run. |
| `einsteintoolkit-validation` | Configured Ubuntu/Apptainer Einstein Toolkit image | Generates `carpet_wavetoy_thorns.py` and `carpet_baikal_thorns.py`, links ETLegacy thorns/fixtures into ET, then builds ET | Runs the configured Baikal, BaikalVacuum, and WaveToyNRPy Cactus testsuites and fails on reported failures. No `carpetx_*` generation/build/run. |
| `dendro-validation` | Configured Ubuntu/Apptainer image, digest-verified before use; 60-minute job limit, five-minute limit per MPI invocation, and one OpenMP thread per rank | Generates BSSN and fCCZ4 at finite-difference orders 4, 6, and 8 with KO dissipation enabled and disabled; builds and runs every generated standalone CTest check; runs focused AddressSanitizer and UndefinedBehaviorSanitizer storage checks for every profile; embeds both formulations in Dendro-GR for FD4/6/8 with KO and FD6 without KO | Runs each selected real-host profile on one and two MPI ranks, runs one two-rank Minkowski step, checks injected transport, callback, element-order, and TOML-profile failures, and prints transient FD6 compile, size, block-RHS time, and allocation measurements without treating them as fixed thresholds. |
| `charmpp-validation` | Configured Ubuntu/Apptainer Charm++ context | Generates and builds the configured superB elliptic, spectroscopy, and collision projects | Runs the configured collision executable through `charmrun`; no explicit scientific-output assertion beyond process success. |
| `sebob-consistency-test` | Configured Ubuntu matrix | Checks out the workflow-selected trusted revision; generates/builds trusted and current SEOBNRv5 variants | Each helper invocation rebuilds both executables, uses exactly ten deterministic inputs, and requires median current/trusted amplitude-plus-phase error not exceed the perturbation-derived baseline. |
| `sebobv2-consistency-test` | Same Ubuntu matrix shape | Generates/builds trusted and current `sebobv2` at the workflow-selected trusted revision | Uses the same ten-input and median-error criterion. |

Claim evidence:
- Claim: Each `sebob-consistency-test` helper invocation uses exactly ten deterministic input sets.
- Role: CI behavior
- Deciding authority: [`sebob_consistency_check.py`](../../nrpy/examples/tests/sebob_consistency_check.py), module `__main__` entry point, `num_sets`
- Corroboration: `none available`; the workflow invokes the helper but does not independently restate its input count

Claim evidence:
- Claim: The `sebobv2-consistency-test` helper invocation uses exactly ten deterministic input sets.
- Role: CI behavior
- Deciding authority: [`sebobv2_consistency_check.py`](../../nrpy/examples/tests/sebobv2_consistency_check.py), module `__main__` entry point, `num_sets`
- Corroboration: `none available`; the workflow invokes the helper but does not independently restate its input count

The `codegen-ubuntu` and `codegen-mac` rows' generate/build-scope cells state
that both jobs build all twelve SEOBNRv5 approximant/calibration variants (3
approximants × {production, `-calibration_no_spin`, `-calibration_spin`,
`-nrpy_calibrated`}) alongside the other named default projects.

Claim evidence:
- Claim: `codegen-ubuntu` and `codegen-mac` each generate and build 24 default C/library projects, including all twelve SEOBNRv5 approximant/calibration variants (3 approximants × {production, `-calibration_no_spin`, `-calibration_spin`, `-nrpy_calibrated`}), with no generated executable, test, or numerical result run for any of them.
- Role: CI behavior
- Deciding authority: [main.yml](../../.github/workflows/main.yml), jobs `codegen-ubuntu`, `codegen-mac`

A successful named build can establish only named generation plus toolchain
compile/link compatibility. Generation completion or file existence is not a
semantic result; build does not prove runtime, and process completion does not
prove semantics. [Code Test Policy](code-test-policy.md) owns validation-layer
names, proof limits, the status-only gate, and the narrow examples exception.

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
inputs, and workflow configuration does not establish execution outcomes.
The Charm++ process-success cell is retained descriptive legacy, not precedent
for adding another generic status-only build or runtime cell.

The local `.github/full_nrpy_local_ci.sh` helper is separate from GitHub job
coverage. It installs dependencies, performs broad static analysis, invokes
its configured generators, and builds selected non-Carpet/non-superB C/library
projects. It generates superB and Carpet/CarpetX families without their
external-host builds; it omits other example families documented in the
helper's source. It then configures selected builds with `--cuda`, including
curvilinear and multicoordinate wave, standalone elliptic, black-hole,
hydro-without-hydro, and TOVola routes. TOVola has no
argument parser or CUDA branch, so its extra `--cuda` token is ignored and that
cell is an ordinary C build. The helper installs no CUDA toolkit, declares no
GPU runner, runs no generated executable, and checks no GPU result. Treat it as
a local command recipe requiring a prepared environment, not CI pass evidence.

Dendro has a dedicated configured job on every workflow invocation. Separate
generated projects cover the complete BSSN/fCCZ4, finite-difference-order
4/6/8, and KO-on/off matrix. Each project runs its full generated standalone
CTest set. Focused sanitizer builds check nonzero offsets and flat storage for
every profile. The job then reconfigures one Dendro-GR build for both
formulations at FD4/6/8 with KO and FD6 without KO. The real-host test programs
exercise block geometry, data transfer, callbacks, parameter forwarding, and
failure termination on one and two MPI ranks. Each qualification executable
also evolves Minkowski data for one RK4 step on two ranks. Host element-order
and TOML-profile mismatches must emit their expected diagnostics. Transient
FD6 timing, allocation, generated-size, and compile measurements are job output,
not stored thresholds. These commands establish configured checks, not a stored
run result. Reproduction details and
proof limits live in [Validation, Standalone Host, And Deferred
Tests](../infrastructures/dendro/validation-standalone-host-and-deferral-gates.md).

Claim evidence:
- Claim: the configured Dendro job runs every standalone BSSN/fCCZ4 order-and-KO profile, focused sanitizer storage checks, and the selected real-host order/KO matrix on one and two MPI ranks; this configuration does not establish a latest successful run.
- Role: CI behavior
- Deciding authority: [main.yml](../../.github/workflows/main.yml), `dendro-validation`
- Corroboration: [general_relativity/self_tests_cpp.py](../../nrpy/infrastructures/Dendro/general_relativity/self_tests_cpp.py), generated numerical checks; [runtime_integration_test.cpp](../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp), real-host block and failure checks; [general_relativity/main_cpp.py](../../nrpy/infrastructures/Dendro/general_relativity/main_cpp.py), Minkowski acceptance checks

Explicitly unsupported or unverified by these configurations: CarpetX build or
runtime; JAX generated-package install/import/basic test or accelerator runtime;
any CUDA executable/GPU result; Dendro general boundaries, distributed
remeshing, local time stepping, restart, output, GPU execution, or threaded
kernels;
long-time or nonlinear Dendro evolution;
geodesic/raytracing projects; GRoovy; active MANGA build; Kasner; and scientific
correctness beyond the stated regression, property, and waveform assertions.

These jobs intentionally create generated `project/` outputs. Treat those as CI
products, not committed documentation or hand-authored source, unless a selected
generated file has been deliberately registered as frozen evidence.

## Sources

- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `codegen-ubuntu`; [.github/full_nrpy_local_ci.sh](../../.github/full_nrpy_local_ci.sh) - `example_scripts`, `cuda_example_scripts`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `codegen-mac`
- [cmake_helpers.py](../../nrpy/infrastructures/Dendro/cmake_helpers.py) - `output_solver_cmake`, `output_tests_cmake`, their `add_test` registrations
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `einsteintoolkit-validation`; official Einstein Toolkit [Adding a test case](https://docs.einsteintoolkit.org/et-docs/Adding_a_test_case) - `A test case is...`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `dendro-validation`
- [general_relativity/self_tests_cpp.py](../../nrpy/infrastructures/Dendro/general_relativity/self_tests_cpp.py) - nonflat RHS, diagnostics, and KO multiprecision oracle
- [general_relativity/main_cpp.py](../../nrpy/infrastructures/Dendro/general_relativity/main_cpp.py) - real-host Minkowski acceptance checks
- [test.ccl](../../nrpy/examples/et_WaveToyfiles/test/test.ccl) - `TEST WaveToyNRPy_test`, `RELTOL 1e-11`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `charmpp-validation`; official Charm++ [Quickstart](https://charm.readthedocs.io/en/v8.0.0/quickstart.html) - `Compiling the Example`, `Running the Example`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `sebob-consistency-test`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `sebobv2-consistency-test`
- [../../README.md](../../README.md) - `## Project Families and Example Generators`
- [../../README.md](../../README.md) - `## What Gets Generated?`
- [../../nrpy/examples/tests/sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py) - `calculate_rmse`
- [../../nrpy/examples/tests/sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py) - `process_input_set`
- [../../nrpy/examples/tests/sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py) - module `__main__` entry point, `num_sets`
- [../../nrpy/examples/tests/sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py) - `calculate_rmse`
- [../../nrpy/examples/tests/sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py) - `process_input_set`
- [../../nrpy/examples/tests/sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py) - module `__main__` entry point, `num_sets`

## See Also

- Parent: [Validation](index.md)
- Depends on: [Code Test Policy](code-test-policy.md)
- Contrasts with: [Static Analysis](static-analysis.md)
- See also: [Glossary](../glossary.md)
- See also: [Workflows](../workflows.md)
