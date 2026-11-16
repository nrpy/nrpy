# Generated Backend Comparison

> Compare generated NRPy backend families by emitted files, build commands, and CI validation scope. · Status: provisional
> Up: [NRPy Knowledge Base](../AGENTS.md)

## Summary

NRPy generates solver code for six host targets: standalone C/CUDA executables or libraries (`BHaH`), classic Cactus thorns (`ETLegacy`), CarpetX/AMReX Cactus thorns (`CarpetX`), Charm++ distributed chare applications (`superB`), Python/JAX packages (`JAX`), and standalone MPI octree CMake applications (`Dendro`).

## Detail

### Shared Generated-Output Rule And Backend Matrix

Across all six backends, files emitted under `project/**` (sources, headers, `.ccl`/`.ci` files, packaging files, binaries, logs, and checkpoints) are generated build or runtime outputs; KB pages cite the Python generator modules, registries, and static headers rather than generated `project/**` trees unless a specific output file is frozen as registered evidence ([Generated Output Boundaries](architecture/generated-output-boundaries.md)).

| Backend | Source family | Generated output | Build and runtime path | Validation route |
| --- | --- | --- | --- | --- |
| BHaH | `nrpy/infrastructures/BHaH/` over `CFunction_dict` (GR, GRHD, elliptic, wave, geodesics, SEOBNR/BOB, BHaHAHA). | `project/<name>/`: `.c`/`.cu`, `BHaH_defines.h`, `BHaH_function_prototypes.h`, `<name>.par`, `Makefile`, executable or `lib*.{a,so,dylib}`. | Run Python generator, `cd project/<name>`, `make`, and run `./<name>` or link `bhah_*` / `bah_*` library entrypoints (GSL required for waveform/spectroscopy examples). | `codegen-ubuntu` and `codegen-mac` generate and compile standalone OpenMP projects; SEOB/SEBOB jobs run numerical waveform checks; equation modules carry `tests/*.py` trusted dictionaries. |
| ETLegacy | `nrpy/infrastructures/ETLegacy/` (`CFunction` registry, MoL schedule bins, Cactus gridfunction groups, CCL writers). | Cactus thorn directory with `interface.ccl`, `param.ccl`, `schedule.ccl`, `src/make.code.defn`, and `src/*.c`. | Generate thorns via `carpet_wavetoy_thorns.py` or `carpet_baikal_thorns.py`, link into an Einstein Toolkit tree, and build/run via Cactus (`simfactory`/Cactus exe). | `einsteintoolkit-validation` builds Einstein Toolkit and runs `WaveToyNRPy`, `Baikal`, and `BaikalVacuum` testsuites; GR RHS modules also check trusted dictionaries. |
| CarpetX | `nrpy/infrastructures/CarpetX/` (`Loop`/`CarpetX` C++ kernels, `ODESolvers` bins, CarpetX gridfunction tags). | Cactus thorn directory with `interface.ccl`, `param.ccl`, `schedule.ccl`, `configuration.ccl`, `src/make.code.defn`, and `src/*.cxx`. | Generate thorns (`carpetx_wavetoy_thorns.py`, `carpetx_baikal_thorns.py`), place in a CarpetX-capable Einstein Toolkit tree, and build with `Loop` and `CarpetX`. | Equation/RHS trusted-expression dictionaries only; checked-in `.github/workflows/main.yml` configures no `carpetx_*` toolkit build or runtime job. |
| superB | `nrpy/infrastructures/superB/` (`Main`, `Timestepping`, horizon/diagnostics chares, `.ci` SDAG, PUP serialization). | `project/<name>/`: `.h`, `.cpp`, `.ci`, PUP files, `BHaH_defines.h`, `<name>.par`, and `charmc` `Makefile`. | Run generator, `cd project/<name>`, `make` (`charmc`), and launch `./charmrun +pN ./<name>`. | `charmpp-validation` generates and builds superB projects in a Charm++ 8.0.0 container and runs `superB_two_blackholes_collide` via `charmrun +p2`. |
| JAX | `nrpy/infrastructures/JAX/` (`PyFunction_dict` and `commondata_params_dict`, exercised by `sebobv1_jax.py`). | `project/<name>/src/<name>/`: one `.py` module per `PyFunction`, `Commondata.py`, `__init__.py`, `pyproject.toml`, `setup.cfg`, and import smoke test. | Run `python -m nrpy.examples.sebobv1_jax` and `pip install` the generated package. | `codegen-ubuntu` and `codegen-mac` install `sebobv1_jax` and call `SEOBNRv5_aligned_spin_coefficients` once without checking returned values ([Generated Project CI](validation/generated-project-ci.md#sebobv1_jax-install-and-call-step)); full waveform runtime remains provisional. |
| Dendro | `nrpy/infrastructures/Dendro/` (`dendro_bssn.py` and `dendro_fccz4.py` standalone octree AMR solvers, distinct from upstream `BSSN_GR`). | `project/Dendro_NRPy_{BSSN,fCCZ4}/`: explicit `CMakeLists.txt` C++ sources, solver context, TOML parser, TwoPunctures, RHS/constraint/Psi4 kernels, and AH wrapper. | Run generator, configure/build with CMake (fetching `Dendro-5.01` and `toml11`), solve initial data with `--tpid`, and evolve under MPI (`mpirun`). | `dendro-validation` runs `dendro_application_check.py` across deterministic generation, CMake builds, FD4/6/8 TwoPunctures runs, short MPI evolutions, remesh, horizon, wave, and checkpoint restart checks. |

Evidence (CI behavior): `.github/workflows/main.yml` (`codegen-ubuntu`, `codegen-mac`, `einsteintoolkit-validation`, `charmpp-validation`, `dendro-validation`); corroboration: `dendro_application_check.py::Leg.generate_and_build` and `Dendro/CMakeLists.py`.

## Sources

- [README.md](../README.md) - `## Project Families and Example Generators`, `## What Gets Generated?`
- [main.yml](../.github/workflows/main.yml) - `codegen-ubuntu`, `codegen-mac`, `einsteintoolkit-validation`, `charmpp-validation`, `dendro-validation`, SEOB/SEBOB jobs
- [dendro_application_check.py](../nrpy/examples/tests/dendro_application_check.py) - `Leg.generate_and_build`, `Leg.run`
- [dendro_bssn.py](../nrpy/examples/dendro_bssn.py) - `main`
- [dendro_fccz4.py](../nrpy/examples/dendro_fccz4.py) - `main`
- [CMakeLists.py](../nrpy/infrastructures/Dendro/CMakeLists.py) - `output_CFunctions_function_prototypes_and_construct_CMakeLists`
- [Karpathy LLM Wiki note](https://gist.githubusercontent.com/karpathy/442a6bf555914893e9891c11519de94f/raw/ac46de1ad27f92b28ac95459c782c07f6b8c964a/llm-wiki.md) - cross-branch synthesis filing motivation

## See Also

- Depends on: [Infrastructures](infrastructures/index.md)
- Depends on: [Generated Output Boundaries](architecture/generated-output-boundaries.md)
- Depends on: [Build And Run](architecture/build-and-run.md)
