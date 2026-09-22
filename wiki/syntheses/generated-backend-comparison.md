# Generated Backend Comparison

> Compare generated backend families by output files, build path, validation route, and generated-file boundary. · Status: provisional
> Up: [Syntheses](index.md)

## Summary

NRPy routes generated-backend families with different runtime hosts:
BHaH emits standalone C/CUDA-style projects or libraries, ETLegacy and CarpetX
emit Einstein Toolkit thorns, superB emits Charm++ projects, JAX emits a
Python/JAX package, and Dendro emits complete BSSN and fCCZ4 applications beside
`BSSN_GR` in one Dendro-GR checkout.
Across them, files under generated `project/**` trees are build or runtime output,
not source evidence,
unless maintainers deliberately freeze and register selected generated files.

This page is a filed query synthesis in the Karpathy LLM Wiki sense: a
cross-branch answer that should compound in the KB instead of disappearing into
chat history.

## Detail

| Backend | Source family | Generated output | Build/runtime path | Validation route | Generated-file boundary |
| --- | --- | --- | --- | --- | --- |
| BHaH | Standalone BHaH infrastructure over registered `CFunction` objects, with GR, GRHD, elliptic, geodesic, waveform, and library-oriented examples. | Generated C/CUDA-style files under `project/<name>/`: C sources, headers, prototypes, parameter files, Makefile, executable target, or shared/static library functions. | Run the Python generator, enter the generated project directory, run `make`, then run the executable or link/use the generated library. Some generated BHaH examples require GSL; library paths expose `bhah_initialize`, `bhah_evolve`, `bhah_diagnostics`, and `bhah_finalize`. | Generated-project CI builds many standalone projects on Ubuntu and macOS. Physics-specific validation also appears in equation trusted-expression pages and example-specific checks. | Cite Python generators and registry symbols instead of generated project files, binaries, diagnostics, checkpoints, images, and archives. |
| ETLegacy | Classic Einstein Toolkit/Cactus thorn infrastructure using `CFunction` registry metadata, CCL writers, MoL schedule bins, and ET-style gridfunction groups. | Generated Cactus thorn directory with `interface.ccl`, `param.ccl`, `schedule.ccl`, `src/make.code.defn`, and thorn-local `src/*.c`. | Generate thorns from NRPy, copy or link them into an Einstein Toolkit checkout, then build and run through the Cactus/Einstein Toolkit environment. | ET validation CI generates Carpet WaveToy and Baikal-style thorns, builds Einstein Toolkit, and runs Baikal, BaikalVacuum, and WaveToyNRPy testsuites. ETLegacy GR RHS also has trusted-expression dictionaries. | Cite ETLegacy writer modules, registry metadata, and trusted-expression evidence instead of generated CCL files and `src/*.c`. |
| CarpetX | CarpetX/Cactus thorn infrastructure using Loop/CarpetX dependencies, C++ source emission, ODESolvers schedule bins, and CarpetX gridfunction metadata. | Generated thorn directory with `interface.ccl`, `param.ccl`, `schedule.ccl`, `configuration.ccl`, `src/make.code.defn`, and thorn-local `src/*.cxx`. | Generate thorns from NRPy, place them in an Einstein Toolkit/CarpetX-capable checkout, then build/run in that host environment. CarpetX thorns require `Loop` and `CarpetX`; schedules use ODESolvers bins. | Current cited CI/validation pages do not establish `carpetx_*` Einstein Toolkit build/test coverage. CarpetX GR RHS has trusted-expression dictionaries; several SIMD/CAHD details are source-observed caveats rather than proven runtime guarantees. | Cite CarpetX writer modules, Cactus/CarpetX background docs only for terminology, and local validation pages for NRPy facts instead of emitted CCL files, configuration, C++ source, and built toolkit output. |
| superB | Charm++-based superB infrastructure for distributed-memory generated applications, with `Main`, `Timestepping`, optional interpolation/horizon chares, and PUP support. | Generated Charm++ project under `project/<project_name>/`: `.h`, `.cpp`, `.ci`, PUP routines, copied static headers, parameter/default files, BHaH defines/prototypes, and a Makefile using `charmc`. | Run the Python generator, enter the generated project directory, run `make`, then launch with `./charmrun +pN ./<project_name>`. Optional BHaHAHA and checkpoint paths add service chares and link inputs. | `charmpp-validation` CI generates several superB workflows, builds them in a Charm++ Apptainer image, and runs `superB_two_blackholes_collide` through `charmrun +p2`. | Cite superB generator modules and static source headers instead of generated Charm++ projects, translated `.decl.h`/`.def.h`, logs, checkpoints, binaries, and linked service output. |
| JAX | Python/JAX infrastructure driven by `PyFunction_dict` and `commondata_params_dict`, currently surfaced through `sebobv1_jax`. | Generated Python package under `project/<name>/src/<name>/`, with one module per registered `PyFunction`, `Commondata.py`, package `__init__.py`, `pyproject.toml`, `setup.cfg`, requirements, README, `.gitignore`, and a minimal import smoke test. | Run `python -m nrpy.examples.sebobv1_jax` or another JAX generator. Current CI generation route does not run a following generated `make` step; runtime use of generated package behavior is narrower than C backend build validation. | Ubuntu and macOS codegen CI run the JAX generator. Current `sebobv1_jax` route is generation-only and has a documented `a_f` Commondata mismatch, so end-to-end waveform runtime validation is provisional. | Cite JAX project generator, `PyFunction`/Commondata registry code, example source, and CI workflow instead of generated Python package files and packaging metadata. |
| Dendro | Dendro infrastructure over NRPy gridfunction, CodeParameter, and `CFunction` registries; BSSN and fCCZ4 are complete applications. | `NRPy_BSSN_GR/` or `NRPy_fCCZ4_GR/` beside `BSSN_GR/`, with application-owned state, parameters, executable, context, TwoPunctures, kernels, runtime services, and explicit CMake sources. | Run either Python generator from the Dendro-GR tree, configure the root or generated application, build against Dendrolib, then run its solver with a Dendro parameter file. | Required validation uses deterministic generation, complete builds, and MPI TwoPunctures runs at FD4/6/8. Current checked-in workflow does not yet configure this application layout. | Cite Dendro generator modules and registrars instead of generated source, binaries, diagnostics, or checkpoints. |

The main backend split is not language alone. BHaH and superB both generate
standalone project trees, but BHaH initializes and evolves through a process or
library while superB uses Charm++ chares. ETLegacy and CarpetX both
generate Cactus thorns, but ETLegacy emits C sources and MoL-oriented schedules
while CarpetX emits C++ sources, `configuration.ccl`, `Loop CarpetX`
requirements, and ODESolvers-oriented schedules. JAX is the outlier: it
consumes Python-function registries and writes a Python package rather than a
C/C++ build. Dendro instead targets Dendrolib and emits complete sibling
applications. Each sibling owns evolution, AMR, checkpoint, boundary,
constraint, wave, horizon, and TwoPunctures paths. It neither wraps nor compiles
the chi-specific `BSSN_GR` implementation.

Claim evidence:
- Claim: the Dendro infrastructure supports fCCZ4 and BSSN applications, each generated by its own top-level example module through the same emitters. Nothing cited here decides the application inventory of another backend family.
- Role: descriptive behavior
- Deciding authority: [dendro_fccz4.py](../../nrpy/examples/dendro_fccz4.py), `main`; [dendro_bssn.py](../../nrpy/examples/dendro_bssn.py), `main`
- Corroboration: [CMakeLists.py](../../nrpy/infrastructures/Dendro/CMakeLists.py), explicit registered-CFunction source emission

Claim evidence:
- Claim: each Dendro example emits one complete sibling application with an explicit source list; current checked-in CI does not yet validate this application layout.
- Role: descriptive behavior
- Deciding authority: [CMakeLists.py](../../nrpy/infrastructures/Dendro/CMakeLists.py), `output_CFunctions_function_prototypes_and_construct_CMakeLists`; [solver_context.py](../../nrpy/infrastructures/Dendro/solver_context.py), generated context
- Corroboration: [main.yml](../../.github/workflows/main.yml), which contains no current complete-sibling application route

Validation is uneven by backend. Selected standalone C examples and superB have
configured build or run coverage in generated-project CI; ETLegacy has configured Einstein Toolkit build/test
coverage, while CarpetX currently has trusted-expression checks and source-level
assembly documentation but no cited `carpetx_*` CI build/test route in this page;
JAX currently has generation coverage and a narrow set of implemented SEOBNRv5
coefficient functions. Treat broad runtime claims for JAX and any backend option
not covered by the cited validation pages as provisional. Workflow configuration
proves configured job shape, never latest successful execution.

## Sources

- Karpathy LLM Wiki approach - query-output filing principle, raw gist:
  `https://gist.githubusercontent.com/karpathy/442a6bf555914893e9891c11519de94f/raw/ac46de1ad27f92b28ac95459c782c07f6b8c964a/llm-wiki.md`.
- [README.md](../../README.md) - `## Project Families and Example Generators`, `## What Gets Generated?`
- [main.yml](../../.github/workflows/main.yml) - `codegen-ubuntu`, `einsteintoolkit-validation`, `charmpp-validation`, SEOB/SEBOB consistency jobs
- [dendro_fccz4.py](../../nrpy/examples/dendro_fccz4.py) - complete fCCZ4 application generation
- [dendro_bssn.py](../../nrpy/examples/dendro_bssn.py) - `main`
- [CMakeLists.py](../../nrpy/infrastructures/Dendro/CMakeLists.py) - explicit application source emission
- [solver_context.py](../../nrpy/infrastructures/Dendro/solver_context.py) - generated evolution and service scheduling

## See Also

- Parent: [Syntheses](index.md)
- Depends on: [Generated Output Boundaries](../architecture/generated-output-boundaries.md)
- Depends on: [Build And Run](../architecture/build-and-run.md)
- Validated by: [Generated Project CI](../validation/generated-project-ci.md)
- See also: [Infrastructures](../infrastructures/index.md)
- Depends on: [BHaH](../infrastructures/bhah/index.md) - standalone BHaH branch router.
- Depends on: [BHaH Lifecycle And Project Assembly](../infrastructures/bhah/lifecycle-and-project-assembly.md) - BHaH project generation, initialization, evolution, diagnostics, and cleanup.
- Depends on: [BHaH GR Application Wiring](../infrastructures/bhah/gr-application-wiring.md) - BHaH GR kernel registration and diagnostics wiring.
- Depends on: [ETLegacy](../infrastructures/etlegacy/index.md) - ETLegacy branch router.
- Depends on: [ETLegacy Thorn Assembly And CCL Files](../infrastructures/etlegacy/thorn-assembly-and-ccl-files.md) - generated Cactus CCL and C source handoff.
- Depends on: [ETLegacy GR BSSN RHS, Ricci, Constraints, And Validation](../infrastructures/etlegacy/gr-bssn-rhs-ricci-constraints-and-validation.md) - ETLegacy trusted-expression validation route.
- Depends on: [CarpetX](../infrastructures/carpetx/index.md) - CarpetX branch router.
- Depends on: [CarpetX Thorn Assembly, Configuration, And CCL Files](../infrastructures/carpetx/thorn-assembly-configuration-and-ccl-files.md) - generated Cactus/CarpetX CCL, configuration, and C++ source handoff.
- Depends on: [CarpetX GR BSSN RHS, Ricci, Constraints, And Validation](../infrastructures/carpetx/gr-bssn-rhs-ricci-constraints-and-validation.md) - CarpetX trusted-expression validation route and caveats.
- Depends on: [superB](../infrastructures/superb/index.md) - superB branch router.
- Depends on: [superB Lifecycle And Project Assembly](../infrastructures/superb/lifecycle-and-project-assembly.md) - Charm++ project assembly, `charmc`, `charmrun`, PUP, and generated files.
- Depends on: [superB Chare Entrypoints And Runtime](../infrastructures/superb/chare-entrypoints-and-runtime.md) - Charm++ chare runtime structure.
- Depends on: [JAX](../infrastructures/jax/index.md) - JAX branch router.
- Depends on: [JAX Project Generation Lifecycle](../infrastructures/jax/project-generation-lifecycle.md) - generated Python package creation and output-file boundary.
- Depends on: [JAX Commondata And PyFunction Registry](../infrastructures/jax/commondata-and-pyfunction-registry.md) - generation-time registries consumed by JAX project generation.
- Depends on: [SEBOBv1 JAX Workflow](../infrastructures/jax/sebobv1-jax-workflow.md) - current JAX example functions, CI generation route, and provisional runtime caveat.
- Depends on: [Dendro](../infrastructures/dendro/index.md) - Dendro branch router.
- Depends on: [Project Assembly And Generating Functions](../infrastructures/dendro/project-assembly-and-emitters.md) - one generating function per output file, reading the registries directly.
- Depends on: [Production Validation And Deferred Checks](../infrastructures/dendro/validation-standalone-host-and-deferral-gates.md) - generation, build, and MPI application checks.
