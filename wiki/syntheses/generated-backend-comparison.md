# Generated Backend Comparison

> Compare generated backend families by output shape, build path, validation route, and artifact boundary. · Status: provisional · Last reconciled: 09-04-2026
> Up: [Syntheses](index.md)

## Summary

NRPy has six routed generated-backend families with different runtime hosts:
BHaH emits standalone C/CUDA-style projects or libraries, ETLegacy and CarpetX
emit Einstein Toolkit thorns, superB emits Charm++ projects, JAX emits a
Python/JAX package, and Dendro emits an fCCZ4 module for a Dendro-GR checkout.
Across all six, generated `project/**` trees are products, not source evidence,
unless maintainers deliberately freeze and register a selected generated
artifact.

This page is a filed query synthesis in the Karpathy LLM Wiki sense: a
cross-branch answer that should compound in the KB instead of disappearing into
chat history.

## Detail

| Backend | Source family | Generated output | Build/runtime path | Validation route | Artifact boundary |
| --- | --- | --- | --- | --- | --- |
| BHaH | Standalone BHaH infrastructure over registered `CFunction` objects, with GR, GRHD, elliptic, geodesic, waveform, and library-oriented examples. | Generated C/CUDA-style project assets under `project/<name>/`: C sources, headers, prototypes, parameter files, Makefile, executable target, or shared/static library surface. | Run the Python generator, enter the generated project directory, run `make`, then run the executable or link/use the generated library. Some workflows need GSL; library paths expose `bhah_initialize`, `bhah_evolve`, `bhah_diagnostics`, and `bhah_finalize`. | Generated-project CI builds many standalone projects on Ubuntu and macOS. Physics-specific validation also appears in equation trusted-expression pages and example-specific checks. | Generated project files, binaries, diagnostics, checkpoints, images, and archives are products. Cite Python generators, registry symbols, and selected frozen evidence instead. |
| ETLegacy | Classic Einstein Toolkit/Cactus thorn infrastructure using `CFunction` registry metadata, CCL writers, MoL schedule bins, and ET-style gridfunction groups. | Generated Cactus thorn directory with `interface.ccl`, `param.ccl`, `schedule.ccl`, `src/make.code.defn`, and thorn-local `src/*.c`. | Generate thorns from NRPy, copy or link them into an Einstein Toolkit checkout, then build and run through the Cactus/Einstein Toolkit environment. | ET validation CI generates Carpet WaveToy and Baikal-style thorns, builds Einstein Toolkit, and runs Baikal, BaikalVacuum, and WaveToyNRPy testsuites. ETLegacy GR RHS also has trusted-expression dictionaries. | Emitted CCL files and `src/*.c` are generated thorn products. Cite ETLegacy writer modules, registry metadata, and trusted-expression evidence. |
| CarpetX | CarpetX/Cactus thorn infrastructure using Loop/CarpetX dependencies, C++ source emission, ODESolvers schedule bins, and CarpetX gridfunction metadata. | Generated thorn directory with `interface.ccl`, `param.ccl`, `schedule.ccl`, `configuration.ccl`, `src/make.code.defn`, and thorn-local `src/*.cxx`. | Generate thorns from NRPy, place them in an Einstein Toolkit/CarpetX-capable checkout, then build/run in that host environment. CarpetX thorns require `Loop` and `CarpetX`; schedules use ODESolvers bins. | Current cited CI/validation pages do not establish `carpetx_*` Einstein Toolkit build/test coverage. CarpetX GR RHS has trusted-expression dictionaries; several SIMD/CAHD details are source-observed caveats rather than proven runtime guarantees. | Emitted CCL, configuration, C++ source, and built toolkit outputs are generated products. Cite CarpetX writer modules, Cactus/CarpetX background docs only for terminology, and local validation pages for NRPy facts. |
| superB | Charm++-based superB infrastructure for distributed-memory generated applications, with `Main`, `Timestepping`, optional interpolation/horizon chares, and PUP support. | Generated Charm++ project under `project/<project_name>/`: `.h`, `.cpp`, `.ci`, PUP routines, copied static headers, parameter/default files, BHaH defines/prototypes, and a Makefile using `charmc`. | Run the Python generator, enter the generated project directory, run `make`, then launch with `./charmrun +pN ./<project_name>`. Optional BHaHAHA and checkpoint paths add service chares and link inputs. | `charmpp-validation` CI generates several superB workflows, builds them in a Charm++ Apptainer image, and runs `superB_two_blackholes_collide` through `charmrun +p2`. | Generated Charm++ projects, translated `.decl.h`/`.def.h`, logs, checkpoints, binaries, and linked service outputs are products. Cite superB generator modules and static source headers. |
| JAX | Python/JAX infrastructure driven by `PyFunction_dict` and `commondata_params_dict`, currently surfaced through `sebobv1_jax`. | Generated Python package under `project/<name>/src/<name>/`, with one module per registered `PyFunction`, `Commondata.py`, package `__init__.py`, `pyproject.toml`, `setup.cfg`, requirements, README, `.gitignore`, and a minimal import smoke test. | Run `python -m nrpy.examples.sebobv1_jax` or another JAX generator. Current CI generation route does not run a following generated `make` step; runtime use of generated package behavior is narrower than C backend build validation. | Ubuntu and macOS codegen CI run the JAX generator. Current `sebobv1_jax` route is generation-only and has a documented `a_f` Commondata mismatch, so end-to-end waveform runtime validation is provisional. | Generated Python package files and packaging metadata are products. Cite JAX project generator, `PyFunction`/Commondata registry code, example source, and CI workflow rather than generated package files. |
| Dendro | Dendro infrastructure over registered `CFunction` objects, frozen into an immutable snapshot before any target file is rendered; fCCZ4 is currently the only application. | Generated `FCCZ4_GR` module under `project/<name>/Dendro-GR/`: state and parameter headers, one source per registered CFunction, module CMake, JSON manifests, a parameter profile, generated self-tests, and the mock host header; rendered verifier and installer tools are written at the project root beside the module. | Run the Python generator, configure the module with CMake, build it, then run the generated executable; the module compiles against a mock host header and needs only MPI and a C++17 compiler. Installing into a real Dendro-GR checkout is a separate explicit tool run that writes a receipt; the module refuses to configure where a real `dendro5` target is defined, so that path is not yet a working host integration. | The `dendro-fccz4` CI job generates, builds under `-Wall -Wextra -Werror`, runs ten generated self-tests, and runs Minkowski lifecycle gates under two MPI ranks. Symbolic and emitted-source contracts run as owner doctests in static analysis. | Generated module files, manifests, receipts, and binaries are products. Cite the Dendro generator modules, the frozen registry records, and the pin and capability records. |

The main backend split is not language alone. BHaH and superB both generate
standalone project trees, but BHaH's lifecycle is process/library oriented
while superB's lifecycle is Charm++ chare oriented. ETLegacy and CarpetX both
generate Cactus thorns, but ETLegacy emits C sources and MoL-oriented schedules
while CarpetX emits C++ sources, `configuration.ccl`, `Loop CarpetX`
requirements, and ODESolvers-oriented schedules. JAX is the outlier: it
consumes Python-function registries and writes a Python package rather than a
C/C++ build product. Dendro differs from all five on a different axis: its
target translators read a frozen snapshot and nothing else, so apart from the
designated self-test template no generated file introduces a state-field name
or a formulation expression that was not registered first.

Claim evidence:
- Claim: Dendro's target translators read only a frozen snapshot, so apart from the designated generated-project self-test template no generated file introduces a state-field name or a formulation expression that was not registered first. Nothing cited here decides whether another backend family also freezes its registries.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/freeze.py`, `freeze_nrpy_dendro_environment`
- Corroboration: `nrpy/infrastructures/Dendro/output_project.py`, module-docstring template policy and its designated self-test exception
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3; backend=Dendro; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=--fd-order 4 --no-ko; date=09-04-2026`

Validation is uneven by backend. Selected standalone C examples and superB have
configured build or run coverage in generated-project CI; ETLegacy has configured Einstein Toolkit build/test
coverage, while CarpetX currently has trusted-expression checks and source-level
assembly documentation but no cited `carpetx_*` CI build/test route in this page;
JAX currently has generation coverage and a narrow implemented SEOBNRv5
coefficient surface. Treat broad runtime claims for JAX and any backend option
not covered by the cited validation pages as provisional. Workflow configuration
proves configured job shape, never latest successful execution.

## Sources

- Karpathy LLM Wiki approach - query-output filing principle:
  `https://gist.githubusercontent.com/karpathy/442a6bf555914893e9891c11519de94f/raw/ac46de1ad27f92b28ac95459c782c07f6b8c964a/llm-wiki.md`.
- [README.md](../../README.md) - `## Project Families and Example Generators`, `## What Gets Generated?`
- [main.yml](../../.github/workflows/main.yml) - `codegen-ubuntu`, `einsteintoolkit-validation`, `charmpp-validation`, SEOB/SEBOB consistency jobs
- [dendro-fccz4-validation.yml](../../.github/workflows/dendro-fccz4-validation.yml) - `dendro-fccz4`
- [output_project.py](../../nrpy/infrastructures/Dendro/output_project.py) - `output_project`, template policy
- [freeze.py](../../nrpy/infrastructures/Dendro/freeze.py) - `freeze_nrpy_dendro_environment`

## See Also

- Parent: [Syntheses](index.md)
- Depends on: [Generated Output Boundaries](../architecture/generated-output-boundaries.md)
- Depends on: [Build And Run](../architecture/build-and-run.md)
- Validated by: [Generated Project CI](../validation/generated-project-ci.md)
- See also: [Infrastructures](../infrastructures/index.md)
- Depends on: [BHaH](../infrastructures/bhah/index.md) - standalone BHaH branch router.
- Depends on: [BHaH Lifecycle And Project Assembly](../infrastructures/bhah/lifecycle-and-project-assembly.md) - BHaH generated project and executable/library lifecycle.
- Depends on: [BHaH GR Application Wiring](../infrastructures/bhah/gr-application-wiring.md) - BHaH GR kernel registration and diagnostics wiring.
- Depends on: [ETLegacy](../infrastructures/etlegacy/index.md) - ETLegacy branch router.
- Depends on: [ETLegacy Thorn Assembly And CCL Files](../infrastructures/etlegacy/thorn-assembly-and-ccl-files.md) - generated Cactus CCL and C source handoff.
- Depends on: [ETLegacy GR BSSN RHS, Ricci, Constraints, And Validation](../infrastructures/etlegacy/gr-bssn-rhs-ricci-constraints-and-validation.md) - ETLegacy trusted-expression validation route.
- Depends on: [CarpetX](../infrastructures/carpetx/index.md) - CarpetX branch router.
- Depends on: [CarpetX Thorn Assembly, Configuration, And CCL Files](../infrastructures/carpetx/thorn-assembly-configuration-and-ccl-files.md) - generated Cactus/CarpetX CCL, configuration, and C++ source handoff.
- Depends on: [CarpetX GR BSSN RHS, Ricci, Constraints, And Validation](../infrastructures/carpetx/gr-bssn-rhs-ricci-constraints-and-validation.md) - CarpetX trusted-expression validation route and caveats.
- Depends on: [superB](../infrastructures/superb/index.md) - superB branch router.
- Depends on: [superB Lifecycle And Project Assembly](../infrastructures/superb/lifecycle-and-project-assembly.md) - Charm++ project assembly, `charmc`, `charmrun`, PUP, and generated assets.
- Depends on: [superB Chare Entrypoints And Runtime](../infrastructures/superb/chare-entrypoints-and-runtime.md) - Charm++ chare runtime structure.
- Depends on: [JAX](../infrastructures/jax/index.md) - JAX branch router.
- Depends on: [JAX Project Generation Lifecycle](../infrastructures/jax/project-generation-lifecycle.md) - generated Python package lifecycle and artifact boundary.
- Depends on: [JAX Commondata And PyFunction Registry](../infrastructures/jax/commondata-and-pyfunction-registry.md) - generation-time registries consumed by JAX project generation.
- Depends on: [SEBOBv1 JAX Workflow](../infrastructures/jax/sebobv1-jax-workflow.md) - current JAX example surface, CI generation route, and provisional runtime caveat.
- Depends on: [Dendro](../infrastructures/dendro/index.md) - Dendro branch router.
- Depends on: [Dendro Lifecycle And Project Assembly](../infrastructures/dendro/lifecycle-and-project-assembly.md) - staged module export, installer receipt, and generated manifests.
- Depends on: [Dendro Frozen Snapshot And Pure Translators](../infrastructures/dendro/frozen-snapshot-and-pure-translators.md) - the freeze boundary that distinguishes this backend from the other five.
