# superB Charm++ Workflows

> Route the three superB example generators to their Charm++ project, build, run, and validation shape. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [Examples](index.md)

## Summary

The `superB_*` examples generate Charm++/superB projects, not plain standalone
executables. They target distributed-memory Charm++ runs through chares,
chare-local grids, PUP support, CkIO-linked builds, and `charmrun` launches.

The three checked examples are `superB_two_blackholes_collide`,
`superB_blackhole_spectroscopy`, and
`superB_nrpyelliptic_conformally_flat`. Each writes `project/<project_name>/`,
prints `make` as the build step, and runs as `./charmrun +pN ./<project_name>`
with `+p4` shown by the examples and `+p2` used by CI for the collision run.
Upstream Charm++ examples and tests can inform pattern matching, but they are
not NRPy CI validation evidence. Each generator deletes and recreates its fixed
project directory; preserve prior output before rerunning.

## Detail

| Example | Generate | Project | Prerequisites | Build/run shape | Validation status |
| --- | --- | --- | --- | --- | --- |
| Two black holes collide | `python -m nrpy.examples.superB_two_blackholes_collide` | `project/superB_two_blackholes_collide/` | Python, Charm++ toolchain, `make`; BHaHAHA library is generated under the project because `enable_BHaHAHA = True` | `make`, then `./charmrun +p4 ./superB_two_blackholes_collide` | CI generates, builds with `make -j2`, and runs `./charmrun +p2 ./superB_two_blackholes_collide` |
| Black-hole spectroscopy | `python -m nrpy.examples.superB_blackhole_spectroscopy` | `project/superB_blackhole_spectroscopy/` | Python, Charm++ toolchain, `make`, GSL, TwoPunctures headers copied by the generator; BHaHAHA and Psi4 services enabled | `make`, then `./charmrun +p4 ./superB_blackhole_spectroscopy`; the generator also prints `./charmrun +p4 ./superB_blackhole_spectroscopy +restart log` as the checkpoint restart command | CI generates and builds with `make -j2`; it does not run the executable or restart command in the current workflow |
| Elliptic conformally flat | `python -m nrpy.examples.superB_nrpyelliptic_conformally_flat` | `project/superB_nrpyelliptic_conformally_flat/` | Python, Charm++ toolchain, `make`; no GSL flag in the example; `--floating_point_precision` accepts `float` or `double` tolerance branches | `make`, then `./charmrun +p4 ./superB_nrpyelliptic_conformally_flat` | CI generates and builds with `make -j2`; it does not run the executable in the current workflow |

These are configured workflow facts, not latest-pass claims. The workflow pins
its container environment to `/opt/charm-8.0.0`; this says nothing about newer
Charm++ releases. Official Charm++ [Quickstart](https://charm.readthedocs.io/en/v8.0.0/quickstart.html)
headings `Compiling the Example` and `Running the Example` corroborate the
external toolchain shape: `charmc` compiles Charm++ applications and
`charmrun +pN` launches them. NRPy's exact generated assets and CI scope still
come from local generators and workflow YAML.

All three examples set the generated infrastructure parameter to BHaH but
assemble through superB helpers. That means they reuse BHaH equation,
diagnostic, reference-metric, boundary-condition, and CodeParameters machinery,
then emit Charm++ project assets: chare interface and C++ files for the `Main`
and `Timestepping` flow, static `superB` headers, PUP routines, and a Makefile
compiled with `charmc`.

Chare decomposition is selected per example through `Nchare0`, `Nchare1`, and
`Nchare2` defaults. The generated grid setup requires each global grid
dimension divided by its matching chare count to be an integer; multi-chare
dimensions are rejected only when the per-chare interior size is smaller than
`NGHOSTS`. The collision
example defaults to spherical coordinates with chares `18 x 2 x 1`. The
spectroscopy example defaults to `SinhCylindrical` with chares `4 x 1 x 4`.
The elliptic example defaults to `SinhSpherical` with chares `16 x 2 x 2`.

Communication-map and runtime internals are owned by the superB infrastructure
pages. At example level, the important facts are that each generator registers
`chare_comm_register_C_functions`, CurviBC chare routines, superB MoL routines,
and `register_CFunction_superB_pup_routines`; copies `superB.h` and
`superB_pup_function_prototypes.h`; links `-module CkIO`; and emits the
project Makefile with `CC="charmc"`. These are Charm++ projects, so run them
through `charmrun` rather than executing a plain standalone binary directly.
Low-level ownership stays in the infrastructure leaves: build and interface
rules belong to [superB Lifecycle And Project Assembly](../infrastructures/superb/lifecycle-and-project-assembly.md),
runtime chare/proxy details belong to
[superB Chare Entry Points And Runtime](../infrastructures/superb/chare-entrypoints-and-runtime.md),
CkIO and volume-output behavior belongs to
[superB Diagnostics And Observables](../infrastructures/superb/diagnostics-and-observables.md),
and CI coverage belongs to [Generated Project CI](../validation/generated-project-ci.md).

`superB_two_blackholes_collide` is the Charm++ analog of the lightweight
standalone black-hole evolution example. It evolves Brill-Lindquist initial
data in `Spherical` coordinates, uses RK4, fourth-order finite differences,
outgoing-radiation boundaries, BSSN RHS and Ricci paths, constraint evaluation,
and BHaHAHA horizon finding. Its BHaHAHA path generates a `BHaHAHA-4o`
subdirectory because `fd_order = 4`, emits `Interpolator3d` and
`Horizon_finder` service chares, adds that directory to recursive make, and
links `-lBHaHAHA-4o`.

`superB_blackhole_spectroscopy` is the heavier Charm++ spectroscopy workflow.
It evolves TwoPunctures-derived data in `SinhCylindrical` coordinates,
defaults to `initial_sep = 0.5` unless `--paper` is supplied, uses eighth-order
finite differences, RK4 by default and SSPRK33 for the paper branch, enables
Kreiss-Oliger dissipation, CAKO, SSL, BHaHAHA, and Psi4, and enables Charm++
checkpointing with `checkpoint_every = 20.0`. Because it registers
TwoPunctures support, its generated Makefile uses `gsl-config --cflags` and
`gsl-config --libs`; this is the GSL-backed superB example.

`superB_nrpyelliptic_conformally_flat` is the Charm++ analog of the standalone
NRPyElliptic conformally flat workflow. It uses `SinhSpherical` coordinates,
RK4, tenth-order finite differences, sixth-order radiation boundary stencils,
`nn_max = 10000`, `eta_damping = 11.0`, and precision-dependent residual
tolerances of `-11.2` for double and `-6.5` for float. It passes
`nrpyelliptic_project=True` to the timestepping generator, uses
`convergence_factor` in the generated command-line parser, and stops relaxation
by calling `stop_conditions_check` before `mainProxy.done()`.

Precision and GSL support are source-specific. The elliptic superB example is
the only one of these three with a command-line precision option. The collision
and spectroscopy examples do not define a precision flag in their argument
parsers. Spectroscopy is the GSL-backed workflow because its Makefile helper
adds `gsl-config` compiler and linker flags; the collision and elliptic
examples do not add GSL flags.

The relationship to standalone BHaH is reuse plus distribution. Standalone BHaH
examples generate C or CUDA-ready projects with a direct executable. superB
examples keep the same equation families and many BHaH support generators, but
replace the top-level runtime with Charm++ chares, PUP serialization,
communication maps, CkIO-aware linking, and `charmrun` launch commands. The
elliptic and black-hole examples should be read as analogs of existing
standalone workflows, while detailed superB internals remain in the
infrastructure leaves.

## Sources

- [superB_two_blackholes_collide.py](../../nrpy/examples/superB_two_blackholes_collide.py) - `project_name`, `CoordSystem`, `Nchare0`, `enable_BHaHAHA`
- [superB_two_blackholes_collide.py](../../nrpy/examples/superB_two_blackholes_collide.py) - `bhahaha.py`, `register_CFunction_superB_pup_routines`, `output_CFunctions_function_prototypes_and_construct_Makefile`
- [superB_blackhole_spectroscopy.py](../../nrpy/examples/superB_blackhole_spectroscopy.py) - `--paper`, `project_name`, `enable_charm_checkpointing`, `enable_psi4`, `addl_CFLAGS`
- [superB_blackhole_spectroscopy.py](../../nrpy/examples/superB_blackhole_spectroscopy.py) - `TwoPunctures`, `BHaHAHA`, `BHaH.checkpointing.register_CFunctions`, `output_CFunctions_function_prototypes_and_construct_Makefile`
- [superB_nrpyelliptic_conformally_flat.py](../../nrpy/examples/superB_nrpyelliptic_conformally_flat.py) - `project_name`, `--floating_point_precision`, `Nchare0`, `nrpyelliptic_project=True`
- [superB_nrpyelliptic_conformally_flat.py](../../nrpy/examples/superB_nrpyelliptic_conformally_flat.py) - `get_log10_residual_tolerance`, `post_MoL_step_forward_in_time`, `addl_libraries=["-module CkIO"]`
- [README.md](../../README.md) - `Prerequisites by Workflow`, `superB / Charm++ Generators`, `What Gets Generated?`; official [Charm++ Quickstart](https://charm.readthedocs.io/en/v8.0.0/quickstart.html) - `Compiling the Example`, `Running the Example`
- [.github/workflows/main.yml](../../.github/workflows/main.yml) - `charmpp-validation`

## See Also

- [Examples](index.md)
- [Elliptic Initial Data](elliptic-initial-data.md)
- [Black Hole Evolution](black-hole-evolution.md)
- [superB](../infrastructures/superb/index.md)
- [Generated Output Boundaries](../architecture/generated-output-boundaries.md)
- [Generated Project CI](../validation/generated-project-ci.md)
- Depends on: [superB Lifecycle And Project Assembly](../infrastructures/superb/lifecycle-and-project-assembly.md) - compiled superB assembly summary.
- Depends on: [superB Chare Entry Points And Runtime](../infrastructures/superb/chare-entrypoints-and-runtime.md) - compiled chare and runtime handle summary.
- Depends on: [superB Diagnostics And Observables](../infrastructures/superb/diagnostics-and-observables.md) - compiled CkIO and volume-output summary.
- Depends on: [superB Grids, Boundaries, MoL, And Initial Data](../infrastructures/superb/grids-boundaries-mol-and-initial-data.md) - compiled chare-grid and MoL summary.
- Depends on: [superB GR, BHaHAHA, Psi4, And Interpolation](../infrastructures/superb/gr-bhahaha-psi4-and-interpolation.md) - compiled service-chare summary.
