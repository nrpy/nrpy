# Lifecycle And Project Assembly

> Explain how superB generators assemble Charm++ projects, emitted assets, build/run guidance, and top-level header/PUP registration. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [superB](index.md)

## Summary

superB examples generate Charm++ projects under `project/<project_name>/`.
They first register physics and infrastructure C functions, then assemble
headers, parameter defaults, Charm++ chare sources and interfaces, PUP support,
static superB headers, and a Makefile that uses `charmc`. The generated project
is a product of the local NRPy generators, so KB citations should point to the
Python generators and static source headers rather than generated files under
`project/`.

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
`main.h`, `main.cpp`, and `main.ci`. The corresponding timestepping wrapper
emits `timestepping.h`, `timestepping.cpp`, and `timestepping.ci`, then
registers the C functions needed by that generated runtime. In Charm++
terminology, the `.ci` files are interface files, and compiling them generates
module-specific `.decl.h` and `.def.h` files used by the C++ implementation.
That generic Charm++ fact comes from upstream Charm++ documentation; the NRPy
fact is that superB emits those `.ci`/C++ pairs from its local generators.

`output_CFunctions_function_prototypes_and_construct_Makefile` writes every
registered `CFunction` from `CFunction_dict` as `.cpp`, writes
`BHaH_function_prototypes.h`, and creates a project Makefile. For superB, the
examples call it with `CC="charmc"` and link `-module CkIO`; optional settings
add GSL flags and BHaHAHA link inputs. The generated Makefile has rules for
Charm++ interface translation, C++ object compilation, the final
`$(CC) -language charm++` link step, optional service chares, and a `clean`
target that removes object files, generated `.decl.h`/`.def.h` headers,
`charmrun`, and checkpoint logs. The helper also offers `compile_Makefile()`,
which can regenerate the Makefile, run parallel `make`, and retry with debug
flags, but the current examples finish by printing manual build/run commands.

The examples' printed post-generation path is deliberately simple: change into
`project/<project_name>`, run `make`, then run the Charm++ executable through
`./charmrun +p4 ./<project_name>` as an example four-processor launch. The
Charm++ quickstart describes `charmc` as the wrapper used to compile Charm++
applications and `charmrun` as the launcher for Charm++ applications; superB's
specific command strings come from the NRPy examples.

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
- [superB_pup.py](../../../nrpy/infrastructures/superB/superB/superB_pup.py) - `register_CFunction_superB_pup_routines`
- [superB.h](../../../nrpy/infrastructures/superB/superB/superB.h) - `__SUPERB_H__`, `ckio.h`, `pup.h`
- [superB_pup_function_prototypes.h](../../../nrpy/infrastructures/superB/superB/superB_pup_function_prototypes.h) - `pup_commondata_struct`, `pup_griddata_chare`
- [Charm++ Quickstart](https://github.com/charmplusplus/charm/blob/main/doc/quickstart.rst) - `Parallel "Hello World" with Charm++`, `Compiling the Example`, `Running the Example`
- [Charm++ Manual](https://github.com/charmplusplus/charm/blob/main/doc/charm%2B%2B/manual.rst) - `Charm++ Interface (.ci) Files`, `Generated Files`

## See Also

- [superB](index.md)
- [BHaH Lifecycle](../bhah-lifecycle.md)
- [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
- [Build And Run](../../architecture/build-and-run.md)
- [C Function Registry](../../core/c-function-registry.md)
