# Generated Project CI

> GitHub Actions workflows, local build scripts, waveform regression checks, and proof limits for generated NRPy projects. · Status: confirmed
> Up: [Validation](index.md)

## Summary

`.github/workflows/main.yml`, `.github/workflows/dendrolib-canary.yml`, `.github/workflows/sebobv2-projected-attachment.yml`, and `.github/workflows/spheroidal_fisheye_roundtrip.yml` configure Linux/macOS C and library builds, `sebobv1_jax` package installation, Einstein Toolkit Carpet regression testsuites, MPI Dendro BSSN/fCCZ4 checks, Charm++/superB builds and `+p2` execution, 10-input SEOBNRv5/`sebobv2` waveform comparisons, `sebobv2` projected-attachment reference/matching checks, and generated BHaH spheroidal coordinate round trips. Workflow definitions establish configured checks rather than recorded run outcomes.

## Detail

### Configured Workflow Matrix

| Workflow step | Runner environment | Generation and build scope | Execution and numerical check scope |
| --- | --- | --- | --- |
| `static-analysis` | Linux / Python matrix | No generated `project/` builds | Module execution + static linters ([Static Analysis](static-analysis.md)) |
| `codegen-ubuntu` | Ubuntu / Python matrix | Installs NRPy; generates in `tmp/` and runs `make` + `make clean` across elliptic, wave, BBH, PN, all SEOBNRv5 approximant x `{production, -calibration_no_spin, -calibration_spin, -nrpy_calibrated}` variants, TOV, hydro, BHaHAHA, and `sebobv2` (MANGA commented out); installs `sebobv1_jax` | Runs no compiled C executable; calls `SEOBNRv5_aligned_spin_coefficients` once without checking output values |
| `codegen-mac` | macOS / Python matrix (Homebrew GSL) | Same C/library `make` builds and `sebobv1_jax` install as Ubuntu; omits Dendro | Runs no compiled C executable; calls `SEOBNRv5_aligned_spin_coefficients` once without value checks |
| `einsteintoolkit-validation` | Ubuntu / Apptainer ET image | Generates `carpet_wavetoy_thorns.py` and `carpet_baikal_thorns.py`, symlinks ETLegacy thorns and test files, builds ET (no `carpetx_*`) | Runs `Baikal`, `BaikalVacuum`, and `WaveToyNRPy` (`RELTOL 1e-11`) Cactus testsuites and fails on nonzero test failures |
| `dendro-validation` | Ubuntu 24.04 + Open MPI, GSL, BLAS/LAPACK, `gfortran`; `formulation: [bssn, fccz4]`; 75-min timeout | Runs `dendro_application_check.py`: double-generates and CMake-builds `W` and `chi` variants (`-DCPU_ARCH=x86-64-v3`) against Dendrolib `master` | Solves `--tpid` and runs short MPI Profile P/O checks, restart identity, 1/3/4-rank agreement, FD4/6/8 initial `H` ordering, stored-reference comparison, and negative exit tests ([Dendro Production Validation](../infrastructures/dendro/production-validation-and-deferred-checks.md)) |
| `dendro-validation-dendrolib-master` (`dendrolib-canary.yml`) | Weekly cron + manual dispatch; same Ubuntu runner/matrix | Shallow-clones Dendrolib `master` via `--dendrolib-ref master` (`FETCHCONTENT_SOURCE_DIR_DENDROLIB`) and prints commit | Same `dendro_application_check.py` suite as `dendro-validation` |
| `charmpp-validation` | Ubuntu / Apptainer Charm++ | Generates and builds `superB` elliptic, spectroscopy, and `two_blackholes_collide` | Runs `superB_two_blackholes_collide` via `charmrun` (exit-status check only; no output assertion) |
| `sebob-consistency-test` | Ubuntu matrix | Checks out trusted git revision; generates and compiles trusted vs current SEOBNRv5 variants | Runs 10 deterministic binary inputs per variant (`sebob_consistency_check.py`), requiring median amplitude-plus-phase RMSE `<=` perturbation baseline |
| `sebobv2-consistency-test` | Ubuntu matrix | Generates and compiles trusted vs current `sebobv2` | Same 10-input median amplitude-plus-phase RMSE check (`sebobv2_consistency_check.py`) |
| `projected-reference-and-matching` (`sebobv2-projected-attachment.yml`) | Ubuntu 24.04 / Python 3.12 | Static analysis + builds scalar reference and `--projected-attachment-diagnostics` `sebobv2` | Runs `sebobv2_projected_attachment_check.py` (RKF45 `10M` orbital/spin reference, `t_attach`/QNM handoff, NQC matching, phase overlap, `--baseline` scalar `stdout` byte-identity) |
| `generated-bhah-roundtrip` (`spheroidal_fisheye_roundtrip.yml`) | Ubuntu 24.04 / Python 3.12, `float` and `double` matrix | Static analysis and BHaH N1/N8 spheroidal converter builds | Runs `spheroidal_fisheye_roundtrip_check.py` with numerical round trips, unequal xy/z maps, near-origin and difficult inversion points, and a shell of N8 Cartesian points in double precision |

Evidence (CI behavior): `main.yml`; `dendrolib-canary.yml`; `sebobv2-projected-attachment.yml`; `spheroidal_fisheye_roundtrip.yml`; `sebob_consistency_check.py`; `sebobv2_consistency_check.py`; `sebobv2_projected_attachment_check.py`; `spheroidal_fisheye_roundtrip_check.py`.

### `sebobv1_jax` install-and-call step

Both `codegen-ubuntu` and `codegen-mac` generate `sebobv1_jax`, run `pip install .` inside the generated package, import `SEOBNRv5_aligned_spin_coefficients`, and evaluate it eagerly (without `jax.jit`) on Python scalars `(1.5, 0.3, -0.2, 0.02, 2.4627455127717882e-05, 50.0)`, printing the returned `Commondata`. The step fails only on generation, install, import, or call exceptions; it asserts no numerical field values and runs no `pytest`, `float32`, traced, or GPU/TPU execution.

Evidence (CI behavior): `.github/workflows/main.yml` (`codegen-ubuntu`, `codegen-mac`).

### Proof Limits, Local CI Script, And Unverified Paths

Compilation proves only code generation plus C/C++/Charm++ compile/link compatibility, and process exit `0` (such as `charmpp-validation`, retained as historical coverage rather than precedent for new status-only checks) proves no numerical property; see [Code Test Policy](code-test-policy.md) and [Narrow Examples Product Exception](examples-product-exception.md). Cactus testsuites ([Adding a test case](https://docs.einsteintoolkit.org/et-docs/Adding_a_test_case)) and 10-input SEOBNR RMSE scripts verify regression tolerances on tested inputs rather than global convergence, while Charm++ [Quickstart](https://charm.readthedocs.io/en/v8.0.0/quickstart.html) documents only `charmc`/`charmrun +pN` invocation (leaving superB restart, `Psi4`, elliptic residual stopping, load balancing, TRAM, sections, and zero-copy paths untested in CI).

`.github/full_nrpy_local_ci.sh` is a manual developer script outside GitHub Actions: it runs static analysis, generates superB and Carpet/CarpetX thorns without compiling them, compiles selected standalone C/library examples, and configures `--cuda` builds for curvilinear/multicoordinate wave, standalone elliptic, black-hole, `hydro_without_hydro`, and `tovola_neutron_star` (`tovola_neutron_star.py` has no CLI parser, so `--cuda` is ignored and compiles plain C) without installing `nvcc` or running any binary.

For full Dendro CLI flags (`--formulation`, `--work-dir`, `--update-reference`, `--dendrolib-ref`), `dendro_application_check_reference.py` digit-rounding rules, and unexercised Dendro modes (`pars/q1.par.lowres.toml`, `--fd-order {4,8}` generation, `KO` off, `--ybs-*`, 3-rank horizon-JSON restore, ASan/UBSan, long-time runs), see [Dendro Production Validation](../infrastructures/dendro/production-validation-and-deferred-checks.md). Also unverified in GitHub Actions: CarpetX compilation/execution, CUDA kernel execution, geodesic/raytracing binaries, `groovy_TOV_BSSN.py`, `manga_bhah_lib.py`, and `kasner_exact_evolution.py`. Generated `project/` directories remain disposable build outputs ([Generated Output Boundaries](../architecture/generated-output-boundaries.md)).

Evidence (CI behavior): `.github/workflows/main.yml`; `.github/full_nrpy_local_ci.sh`; `nrpy/examples/et_WaveToyfiles/test/test.ccl`; `dendro_application_check.py`.

## Sources

- [main.yml](../../.github/workflows/main.yml) - `static-analysis`, `codegen-ubuntu`, `codegen-mac`, `einsteintoolkit-validation`, `dendro-validation`, `charmpp-validation`, `sebob*-consistency-test`
- [dendrolib-canary.yml](../../.github/workflows/dendrolib-canary.yml) - `dendro-validation-dendrolib-master`
- [sebobv2-projected-attachment.yml](../../.github/workflows/sebobv2-projected-attachment.yml) - `projected-reference-and-matching`
- [spheroidal_fisheye_roundtrip.yml](../../.github/workflows/spheroidal_fisheye_roundtrip.yml) - `generated-bhah-roundtrip`
- [full_nrpy_local_ci.sh](../../.github/full_nrpy_local_ci.sh) - `example_scripts`, `cuda_example_scripts`
- [sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py) - `calculate_rmse`, `process_input_set`, `__main__`
- [sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py) - `calculate_rmse`, `process_input_set`, `__main__`
- [sebobv2_projected_attachment_check.py](../../nrpy/examples/tests/sebobv2_projected_attachment_check.py) - `check_handoff`, `check_matching`, `check_phase_overlap`, `main`
- [spheroidal_fisheye_roundtrip_check.py](../../nrpy/examples/tests/spheroidal_fisheye_roundtrip_check.py) - `main`, generated BHaH N1/N8 coordinate checks
- [dendro_application_check.py](../../nrpy/examples/tests/dendro_application_check.py) - `Leg.run`, `PROFILE_P`, `PROFILE_O`, `Leg.run_negatives`
- [dendro_application_check_reference.py](../../nrpy/examples/tests/dendro_application_check_reference.py) - `trusted_dict`
- [CMakeLists.py](../../nrpy/infrastructures/Dendro/CMakeLists.py) - `FetchContent_Declare` for Dendrolib `master`
- [test.ccl](../../nrpy/examples/et_WaveToyfiles/test/test.ccl) - `TEST WaveToyNRPy_test`, `RELTOL 1e-11`
- [README.md](../../README.md) - `## Project Families and Example Generators`, `## What Gets Generated?`
- [Adding a test case](https://docs.einsteintoolkit.org/et-docs/Adding_a_test_case) - `A test case is...`
- [Charm++ Quickstart](https://charm.readthedocs.io/en/v8.0.0/quickstart.html) - `Compiling the Example`, `Running the Example`

## See Also

- Depends on: [Code Test Policy](code-test-policy.md)
- Contrasts with: [Static Analysis](static-analysis.md)
- See also: [Test Oracles And Safe Updates](test-oracles-and-safe-updates.md)
