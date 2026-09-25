# Generated Project CI

> CI coverage for generated projects, external backend validation, and waveform consistency checks. · Status: confirmed
> Up: [Validation](index.md)

## Summary

Workflow YAML separates static analysis, Ubuntu/macOS code generation, ETLegacy
regression, Dendro application validation plus a weekly run of the same checks
against the head of Dendrolib's master branch, Charm++/superB, and
trusted/current waveform consistency routes. These are configured routes, not
execution-result snapshots.

## Detail

Configured GitHub job map:

| Job | Configured context | Generate/build scope | Run/result-check scope |
| --- | --- | --- | --- |
| `static-analysis` | Configured Linux/Python matrix | No generated project coverage | Python file execution plus version-dependent static checks; see [Static Analysis](static-analysis.md) |
| `codegen-ubuntu` | Configured Ubuntu/Python matrix | Installs NRPy, generates in `tmp/`, and builds the selected default C/library projects with `make`, spanning elliptic, wave, black-hole, PN, SEOBNR, TOV, hydro, BHaHAHA, and `sebobv2` routes. It generates `sebobv1_jax` without package install/build. | The `make` builds run no generated executable, and `make clean` follows each; MANGA commands are commented out. |
| `codegen-mac` | Configured macOS/Python matrix | Same selected default C/library builds and JAX generation as Ubuntu; no Dendro generation or build; GSL installed with Homebrew | No generated executable, test, or numerical result is run. |
| `einsteintoolkit-validation` | Configured Ubuntu/Apptainer Einstein Toolkit image | Generates `carpet_wavetoy_thorns.py` and `carpet_baikal_thorns.py`, links ETLegacy thorns/fixtures into ET, then builds ET | Runs the configured Baikal, BaikalVacuum, and WaveToyNRPy Cactus testsuites and fails on reported failures. No `carpetx_*` generation/build/run. |
| `dendro-validation` | Ubuntu 24.04 runner with apt-installed Open MPI, GSL, BLAS/LAPACK, and gfortran; matrix `formulation: [bssn, fccz4]`; 75-minute job timeout | Each leg runs `nrpy/examples/tests/dendro_application_check.py`, which generates the W and chi projects twice, configures and builds them against the pinned Dendrolib, solves TwoPunctures once per variant, and runs short MPI evolutions with forced remeshing, horizon finds, wave extraction, checkpoint and restore, 1/3/4-rank repeats, FD4/6/8 initialization, a stored-reference comparison of the evolved diagnostics, and process-boundary rejections. | Proves only the named layers for the helper's two CI profiles: generation determinism, compile/link compatibility, closed-form TwoPunctures ADM and horizon-mass agreement, symmetry properties, restart identity, rank-count agreement, agreement with the stored reference for the evolved diagnostics (a regression check, not a correctness proof), ordering of the initial Hamiltonian-constraint norm with FD order (not a convergence test), and the listed rejections. It is not long-time, merger, or production-resolution evidence. |
| `dendro-validation-dendrolib-master` (`dendrolib-canary.yml`) | Weekly schedule and manual dispatch only; same runner, packages, matrix, and timeout as `dendro-validation`. | Runs the same helper with `--dendrolib-ref master`: Dendrolib is cloned at the head of its master branch instead of the pinned commit, and the resolved commit is printed. | Same checks as `dendro-validation`; a pass or failure is evidence only for the printed Dendrolib commit, and pull requests are unaffected. |
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

Dendro's job runs one helper per formulation. The helper generates the W and chi
sibling applications twice with separate caches and requires byte-identical
trees, then configures and builds each with `CPU_ARCH=x86-64-v3` against the
Dendrolib and toml11 revisions pinned by the generated `CMakeLists.txt`. One
TwoPunctures solve is shared by all runs of each conformal-factor variant, and
each check made for both variants is reported as one W-and-chi result. Profile P
(maximum depth 13, FD6, eight steps, remesh, horizon, wave, and checkpoint
cadence 4) checks the solved ADM energy and angular momentum against the closed
forms for the conformal-factor rescaling, horizon irreducible masses against the
puncture ADM masses, equal-mass reflection symmetry of the horizons, vanishing
ADM momentum and in-plane angular momentum, vanishing odd-m wave modes and the
reflection relation C(l,-m) = (-1)^l conj C(l,m) over all 21 modes with l = 2-4
at two radii, point reflection of the checkpointed puncture centers at step 4
and their motion along the puncture momenta, byte-identical `dat/`, `bah/`, and
`vtu/` outputs after a stop at step 4 and restore, and agreement of a 3-rank
run, which writes a horizon checkpoint, with the diagnostics of the main runs
within a relative tolerance (horizon observables through the irreducible mass;
the finder's convergence residuals are not compared). Run A's constraint and ADM
rows and horizon observables at steps 0, 4, and 8 must also match a stored
reference, `nrpy/examples/tests/dendro_application_check_reference.py`, within a
fixed relative and absolute tolerance; because it covers the evolved state after
the forced remesh, including node counts, it is the job's regression check on
the evolution equations, gauge, and grid transfer. Profile O (maximum depth 10,
four steps) checks that the initial Hamiltonian-constraint norm falls by at
least half per order increase from FD4 to FD8 on one shared mesh, an ordering
check rather than a convergence test, plus a 1-rank repeat. The W and chi
variants must agree at the shared initial diagnostic. Every evolution run the
helper checks (runs A, B, B restored, C, and the four profile-O runs) must pass
one sanity check: finite output, constraint rows at the configured cadence, a
node ceiling, and no unread-parameter warning. Negative cases require exit
status 1-123 and a named diagnostic for a W/chi cross-restore, `--tpid` on more
than one rank, a TwoPunctures parameter mismatch, a missing TwoPunctures file,
an unsupported element order, an unsupported refinement mode, an excessive CFL
factor, `TPID_REPLACE_LAPSE_WITH_SQRT_CHI = false`, an integer given for a
real-valued parameter, and a lapse blow-up with constraint output off; a run
given an unread key must warn about it and still succeed.

The stored reference is a generated `trusted_dict` keyed by formulation and
conformal factor. It changes only through the helper's `--update-reference`
mode, which writes run A's values as a candidate instead of comparing; the
candidate is reviewed as a diff and then compared in a second run without the
flag. An intentional change to the evolution, the gauge, the remesh, or the
Dendrolib revision therefore requires regenerating the reference in the same
change, and the weekly run against Dendrolib master fails when an upstream
change alters these values beyond the tolerance. The constraint and ADM rows are
printed with up to ten significant digits; the horizon radii and circumferences
with ten, area and irreducible mass with sixteen, and time and centroid in fixed
point. The last printed digits can vary between math-library paths and rank
counts; a one-unit difference in the tenth digit lies within the tolerance. The
ADM linear-momentum components and the angular-momentum components J_x and J_y
vanish for this configuration, so their printed digits are round-off and the
absolute tolerance, not the stored value, bounds them. The check reports the
worst entry.

Claim evidence:
- Claim: the Dendro helper compares run A's constraint and ADM rows and horizon observables at steps 0, 4, and 8 with the stored `trusted_dict` reference within a fixed relative and absolute tolerance, so the weekly Dendrolib-master run fails when an upstream change alters those values beyond the tolerance; the reference changes only through `--update-reference`, which writes a candidate instead of comparing and refuses to write when any check failed or when a Dendrolib branch or tag replaces the pin.
- Role: CI behavior
- Deciding authority: [dendro_application_check.py](../../nrpy/examples/tests/dendro_application_check.py), `REFERENCE_STEPS`, `REFERENCE_RTOL`, `REFERENCE_ATOL`, `read_reference`, `Leg.check_reference`, `Leg.run`, `main`; [dendrolib-canary.yml](../../.github/workflows/dendrolib-canary.yml), `dendro-validation-dendrolib-master`
- Corroboration: [dendro_application_check_reference.py](../../nrpy/examples/tests/dendro_application_check_reference.py), `trusted_dict`; [Test Oracles And Safe Updates](test-oracles-and-safe-updates.md), `Two-Process Oracle Update`

A separate workflow, `dendrolib-canary.yml`, runs the same two helper legs
weekly and on manual dispatch with `--dendrolib-ref master`. The helper then
shallow-clones that Dendrolib branch from the repository declared in the
generated `CMakeLists.txt`, configures both variants with
`FETCHCONTENT_SOURCE_DIR_DENDROLIB` pointing at the clone, and prints the
resolved commit. Pull requests and pushes keep building against the pinned
commit, so an upstream change that breaks the generated applications, or changes
their stored-reference values beyond the tolerance, shows up in this weekly run
without failing unrelated changes; a numerical change there calls for moving the
pin and regenerating the reference in the same change; its results are evidence
only for the printed commit.

Claim evidence:
- Claim: `dendrolib-canary.yml` runs the Dendro helper legs weekly and on manual dispatch against the head of Dendrolib's master branch, printing the resolved commit, while `dendro-validation` in `main.yml` keeps the commit pinned by the generator; a pass or failure is evidence only for the printed Dendrolib commit, and the configuration does not establish any run outcome.
- Role: CI behavior
- Deciding authority: [dendrolib-canary.yml](../../.github/workflows/dendrolib-canary.yml), `on`, `dendro-validation-dendrolib-master`; [dendro_application_check.py](../../nrpy/examples/tests/dendro_application_check.py), `Leg.generate_and_build`, `--dendrolib-ref`
- Corroboration: [main.yml](../../.github/workflows/main.yml), `dendro-validation`; [CMakeLists.py](../../nrpy/infrastructures/Dendro/CMakeLists.py), the emitted Dendrolib `FetchContent_Declare`

The helper generates only with Kreiss-Oliger dissipation enabled. The apt
packages, compilers, and `requirements.txt` packages are not pinned; the helper
prints their resolved versions, and a pass is evidence only for those versions.
Every subprocess has an argument vector, a timeout, and a bounded log tail on
failure, and the work directory is removed unconditionally.

Claim evidence:
- Claim: the `dendro-validation` job generates, builds, runs, and checks both complete sibling Dendro applications in W and chi variants through `dendro_application_check.py`, with the checks and rejections listed in this section.
- Role: CI behavior
- Deciding authority: [main.yml](../../.github/workflows/main.yml), `dendro-validation`; [dendro_application_check.py](../../nrpy/examples/tests/dendro_application_check.py), `Leg.run`, `Leg.run_variant`, `Leg.check_run_a`, `Leg.check_reference`, `Leg.compare_runs`, `Leg.run_negatives`
- Corroboration: [dendro_bssn.py](../../nrpy/examples/dendro_bssn.py) and [dendro_fccz4.py](../../nrpy/examples/dendro_fccz4.py), current command-line interface; [CMakeLists.py](../../nrpy/infrastructures/Dendro/CMakeLists.py), pinned dependency revisions

Claim evidence:
- Claim: the helper's pass results cover only its two CI profiles, eight-step or shorter evolutions, and Kreiss-Oliger-enabled generation; they are not evidence for long-time, merger, production-resolution, or KO-off (`--no-ko`) Dendro behavior.
- Role: CI behavior
- Deciding authority: [dendro_application_check.py](../../nrpy/examples/tests/dendro_application_check.py), `PROFILE_P`, `PROFILE_O`, `COMMON_OVERRIDES`, `Leg.generate_and_build`
- Corroboration: [Production Validation And Deferred Checks](../infrastructures/dendro/validation-standalone-host-and-deferral-gates.md), `Required application checks`

Explicitly unsupported or unverified by these configurations: CarpetX build or
runtime; JAX generated-package install/import/basic test or accelerator runtime;
any CUDA executable/GPU result; Dendro general boundaries, local time stepping,
GPU execution, threaded kernels, three-rank horizon-checkpoint contents and
restore, or KO-off (`--no-ko`) generation and build; long-time, merger, or
production-resolution Dendro evolution; geodesic/raytracing projects; GRoovy;
active MANGA build; Kasner; and scientific correctness beyond the stated
regression, property, and waveform assertions.

These jobs intentionally create generated `project/` outputs. Treat those as CI
products, not committed documentation or hand-authored source, unless a selected
generated file has been deliberately registered as frozen evidence.

## Sources

- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `codegen-ubuntu`; [.github/full_nrpy_local_ci.sh](../../.github/full_nrpy_local_ci.sh) - `example_scripts`, `cuda_example_scripts`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `codegen-mac`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `einsteintoolkit-validation`; official Einstein Toolkit [Adding a test case](https://docs.einsteintoolkit.org/et-docs/Adding_a_test_case) - `A test case is...`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `dendro-validation`
- [dendro_application_check.py](../../nrpy/examples/tests/dendro_application_check.py) - `Leg.run`, profiles, checks, and negative cases
- [dendro_application_check_reference.py](../../nrpy/examples/tests/dendro_application_check_reference.py) - stored run-A reference values
- [solver_context.py](../../nrpy/infrastructures/Dendro/solver_context.py) - `Ctx::diagnostic_output`, `Ctx::adm_output` output precision
- [diagnostics_file_output.py](../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_file_output.py) - horizon diagnostics column formats
- [dendrolib-canary.yml](../../.github/workflows/dendrolib-canary.yml) - weekly `dendro-validation-dendrolib-master`
- [dendro_bssn.py](../../nrpy/examples/dendro_bssn.py) - current BSSN generation interface
- [dendro_fccz4.py](../../nrpy/examples/dendro_fccz4.py) - current fCCZ4 generation interface
- [CMakeLists.py](../../nrpy/infrastructures/Dendro/CMakeLists.py) - current generated source and target emission
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
