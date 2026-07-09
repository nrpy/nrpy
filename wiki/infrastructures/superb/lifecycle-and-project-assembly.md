# Lifecycle And Project Assembly

> Explain how superB generators assemble Charm++ projects, emitted assets, build/run guidance, and top-level header/PUP registration. · Status: confirmed · Last reconciled: 07-07-2026
> Up: [superB](index.md)

## Summary

superB examples generate Charm++ projects under `project/<project_name>/`.
They first register physics and infrastructure C functions, then assemble
headers, parameter defaults, Charm++ chare sources and interfaces, PUP support,
static superB headers, and a Makefile that uses `charmc`. The generated project
is a product of the local NRPy generators, so KB citations should point to the
Python generators and static source headers rather than generated files under
`project/`.

This page owns superB build and interface-generation facts: `.ci` files,
`mainmodule` and `extern module` declarations, generated `.decl.h`/`.def.h`
headers, `charmc`/`charmxi` translation context, generated Makefile rules,
`-language charm++`, `-module CkIO`, and `charmrun` launch shape. CkIO is only a
build/link fact here; CkIO sessions and callbacks are owned by
[Diagnostics And Observables](diagnostics-and-observables.md).

## Detail

The public README classifies superB as NRPy's Charm++-based infrastructure for
distributed-memory workflows and lists three entry examples:
`superB_two_blackholes_collide`, `superB_blackhole_spectroscopy`, and
`superB_nrpyelliptic_conformally_flat`. It also distinguishes superB output
from standalone BHaH executables: superB generators produce Charm++ projects,
with GSL required by some workflows.

Each example chooses a `project_name`, derives `project_dir` as
`project/<project_name>`, registers its physics-specific C functions, and then
enters a final assembly phase. The common assembly steps are to write
`CodeParameters` headers, register the commondata default setter, generate a
default parameter file, register the command-line/parameter-file parser,
register superB PUP routines, copy the static `superB.h` and
`superB_pup_function_prototypes.h` headers, emit the `Main` and `Timestepping`
Charm++ files, register griddata cleanup, emit `BHaH_defines.h` with
`superB/superB.h` in its additional includes, and construct the Makefile.

The GR black-hole collision example adds optional BHaHAHA project assembly:
when enabled, it emits `Interpolator3d` and `Horizon_finder` chare files,
includes the BHaHAHA header, adds the BHaHAHA subdirectory to recursive make,
and links against that generated library. The spectroscopy example follows the
same general assembly pattern but also copies TwoPunctures headers, adds GSL
compiler and linker flags, and prints a checkpoint restart command when Charm++
checkpointing is enabled. The NRPyElliptic example uses the same superB
assembly route for an elliptic workflow and passes `nrpyelliptic_project=True`
into the timestepping generator.

`output_commondata_object_h_and_main_h_cpp_ci` emits `commondata_object.h`,
`main.h`, `main.cpp`, and `main.ci`. The generated `main.ci` is a Charm++
`mainmodule main` interface. It declares the `Main` mainchare, publishes
readonly proxies, always reaches `extern module timestepping`, and also reaches
`extern module horizon_finder` and `extern module interpolator3d` when BHaHAHA
support is emitted. The corresponding timestepping wrapper emits
`timestepping.h`, `timestepping.cpp`, and `timestepping.ci`, then registers the
C functions needed by that generated runtime. Optional service emitters own
`interpolator3d.ci` and `horizon_finder.ci`.

In Charm++ terminology, `.ci` files are interface files, and compiling them
generates module-specific `.decl.h` and `.def.h` files used by the C++
implementation. superB header names follow the module names: `main.ci` produces
`main.decl.h`/`main.def.h`, `timestepping.ci` produces
`timestepping.decl.h`/`timestepping.def.h`, and the optional service modules
produce `interpolator3d.*` and `horizon_finder.*` generated headers. The
generic Charm++ fact comes from upstream Charm++ documentation; the NRPy fact
is that superB emits those `.ci`/C++ pairs from its local generators. The
Charm++ interface translator is `charmxi`; generated superB Makefiles normally
reach that translator through `charmc` by compiling the `.ci` targets with
`$(CC)`.

`output_CFunctions_function_prototypes_and_construct_Makefile` writes every
registered `CFunction` from `CFunction_dict` as `.cpp`, writes
`BHaH_function_prototypes.h`, and creates a project Makefile. For superB, the
examples call it with `CC="charmc"` and link `-module CkIO`; optional settings
add GSL flags and BHaHAHA link inputs. The generated Makefile encodes expected
`.ci` outputs directly rather than relying on generated dependency output: it
has explicit rules such as `timestepping.decl.h timestepping.def.h:
timestepping.ci`, `main.decl.h main.def.h: main.ci`, and, when BHaHAHA is
enabled, matching rules for `interpolator3d.ci` and `horizon_finder.ci`. The
same Makefile handles C++ object compilation, the final `$(CC) -language
charm++` link step, optional service chares, and a `clean` target that removes
object files, generated `.decl.h`/`.def.h` headers, `charmrun`, and checkpoint
logs. The helper also offers `compile_Makefile()`, which can regenerate the
Makefile, run parallel `make`, and retry with debug flags, but the current
examples finish by printing manual build/run commands.

The examples' printed post-generation path is deliberately simple: change into
`project/<project_name>`, run `make`, then run the Charm++ executable through
`./charmrun +p4 ./<project_name>` as an example four-processor launch. The
Charm++ quickstart describes `charmc` as the wrapper used to compile Charm++
applications and `charmrun` as the launcher for Charm++ applications; superB's
specific command strings come from the NRPy examples.

Generated outputs are not durable KB sources by default. Do not cite generated
`project/**` files, generated `.decl.h`/`.def.h` headers, generated
executables, copied launchers, checkpoint logs, or runtime outputs for this
page unless a selected generated artifact has been deliberately registered as
frozen evidence. Cite local generators, static checked-in headers, examples,
and workflow files instead.

Charm++ version statements on this page are validation-scoped. The current CI
job exports `/opt/charm-8.0.0` paths before building and running generated
superB projects, so `8.0.0` is a workflow-environment fact, not a claim that
local generators target the latest upstream Charm++ source.

Top-level PUP support is registered before the final Makefile is emitted.
`register_CFunction_superB_pup_routines` registers
`superB_pup_routines.cpp` under the `superB` subdirectory with includes for
`BHaH_defines.h` and `BHaH_function_prototypes.h`. Its generated routines cover
the shared commondata, parameter, boundary-condition, Method of Lines,
communication, diagnostic, temporary-buffer, nonlocal-inner-boundary, and
griddata structures needed by superB. The copied
`superB_pup_function_prototypes.h` header declares those PUP entry points, and
`superB.h` supplies shared superB macros, Charm++/CkIO/PUP includes, and common
struct declarations consumed by generated code.

## Sources

- [README.md](../../../README.md) - `superB / Charm++ Generators`, `What Gets Generated?`
- [superB_two_blackholes_collide.py](../../../nrpy/examples/superB_two_blackholes_collide.py) - `project_name`, `STEP 3`, `output_CFunctions_function_prototypes_and_construct_Makefile`
- [superB_blackhole_spectroscopy.py](../../../nrpy/examples/superB_blackhole_spectroscopy.py) - `enable_charm_checkpointing`, `STEP 3`, checkpoint restart print
- [superB_nrpyelliptic_conformally_flat.py](../../../nrpy/examples/superB_nrpyelliptic_conformally_flat.py) - `project_name`, `nrpyelliptic_project=True`, `STEP 3`
- [Makefile_helpers.py](../../../nrpy/infrastructures/superB/Makefile_helpers.py) - `output_CFunctions_function_prototypes_and_construct_Makefile`, `compile_Makefile`
- [main_chare.py](../../../nrpy/infrastructures/superB/main_chare.py) - `output_commondata_object_h_and_main_h_cpp_ci`, `output_main_ci`
- [timestepping_chare.py](../../../nrpy/infrastructures/superB/timestepping_chare.py) - `output_timestepping_h_cpp_ci_register_CFunctions`, `output_timestepping_ci`
- [interpolator3d_chare.py](../../../nrpy/infrastructures/superB/interpolator3d_chare.py) - `output_interpolator3d_ci`
- [horizon_finder_chare.py](../../../nrpy/infrastructures/superB/horizon_finder_chare.py) - `output_horizon_finder_ci`
- [superB_pup.py](../../../nrpy/infrastructures/superB/superB/superB_pup.py) - `register_CFunction_superB_pup_routines`
- [superB.h](../../../nrpy/infrastructures/superB/superB/superB.h) - `__SUPERB_H__`, `ckio.h`, `pup.h`
- [superB_pup_function_prototypes.h](../../../nrpy/infrastructures/superB/superB/superB_pup_function_prototypes.h) - `pup_commondata_struct`, `pup_griddata_chare`
- [.github/workflows/main.yml](../../../.github/workflows/main.yml) - `charmpp-validation`
- [Charm++ Quickstart](https://github.com/charmplusplus/charm/blob/main/doc/quickstart.rst) - `Parallel "Hello World" with Charm++`, `Compiling the Example`, `Running the Example`
- [Charm++ Manual](https://github.com/charmplusplus/charm/blob/main/doc/charm%2B%2B/manual.rst) - `Charm++ Interface (.ci) Files`, `Generated Files`
- [Charm++ charmc](https://github.com/charmplusplus/charm/blob/main/src/scripts/charmc) - `charmxi`

## See Also

- [superB](index.md)
- [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
- [Build And Run](../../architecture/build-and-run.md)
- [C Function Registry](../../core/c-function-registry.md)
- [Diagnostics And Observables](diagnostics-and-observables.md)
- [Generated Project CI](../../validation/generated-project-ci.md)
