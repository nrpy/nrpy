# BHaH Lifecycle

> Summarize how BHaH standalone applications are registered, assembled, run, and cleaned up. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Infrastructures](index.md)

## Summary

BHaH is NRPy's standalone generated-application infrastructure. Example
generators register C functions into the global C-function registry, emit
headers and a Makefile under `project/<name>/`, and use `register_CFunction_main_c`
as the final assembler for a standard initialize, evolve, diagnose, and cleanup
lifecycle.

## Detail

`register_CFunction_main_c` emits the top-level C `main` function after checking
that required pieces have already been registered: commondata and params default
setters, the command-line and parameter-file parser, numerical grid/timestep
setup, Method of Lines allocation/free/step functions, diagnostics, and initial
data.

The generated `main` lifecycle is ordered around common BHaH runtime state. It
initializes `commondata`, parses parameter-file and command-line inputs,
allocates `griddata`, initializes grid-local parameters, builds numerical grids
and the timestep, allocates evolved and auxiliary gridfunctions, sets initial
data, enters a `while(commondata.time < commondata.t_final)` loop, runs
diagnostics, advances with Method of Lines, then frees allocated memory.

`register_CFunction_numerical_grids_and_timestep` and the wrapper
`register_CFunctions` own grid geometry and timestep setup. They register grid
parameters such as `Nxx`, `Nxx_plus_2NGHOSTS`, `dxx`, `invdxx`, `gridname`,
`CoordSystem_hash`, `grid_idx`, and `convergence_factor`, then emit code for
independent-grid or multipatch setup, optional reference-metric precompute,
optional curvilinear boundary-condition setup, CFL timestep selection, and
initial time bookkeeping.

`griddata_commondata` records extra struct declarations that belong in
`griddata_struct` or `commondata_struct`, and registers host/device cleanup
functions. `simple_loop` supplies the standard BHaH loop skeleton for all points,
interior points, or interior-plus-one-upper regions, with OpenMP pragmas, CUDA
index ranges, optional SIMD increments, and optional reference-metric reads.
The generic `loop()`, `GPUKernel`, and host/device utility helpers are owned by
the Core helper pages; this BHaH page owns the BHaH loop-region and lifecycle
choices that consume them.

`bhah_lib` is a related library interface. Its aggregate registration function
creates `bhah_initialize`, `bhah_evolve`, `bhah_diagnostics`, and
`bhah_finalize`, plus a `BHaH_struct` wrapper containing `commondata` and
`griddata` pointers.

The black-hole example shows the full GR-facing use of this lifecycle: it
registers initial data, numerical grids, diagnostics, BSSN RHS/Ricci/constraint
functions, curvilinear boundary conditions, Method of Lines timestepping,
coordinate transforms, headers, the generated `main`, grid cleanup, and the
project Makefile.

## Sources

- [README.md](../../README.md) - `Standalone BHaH Generators`, `What Gets Generated?`
- [main_c.py](../../nrpy/infrastructures/BHaH/main_c.py) - `register_CFunction_main_c`
- [bhah_lib.py](../../nrpy/infrastructures/BHaH/bhah_lib.py) - `register_CFunctions_bhah_lib`, `BHaH_struct`
- [simple_loop.py](../../nrpy/infrastructures/BHaH/simple_loop.py) - `simple_loop`
- [griddata_commondata.py](../../nrpy/infrastructures/BHaH/griddata_commondata.py) - `GridCommonData`, `register_griddata_commondata`
- [griddata_commondata.py](../../nrpy/infrastructures/BHaH/griddata_commondata.py) - `register_CFunction_griddata_free`
- [numerical_grids_and_timestep.py](../../nrpy/infrastructures/BHaH/numerical_grids_and_timestep.py) - `register_CFunctions`, `register_CFunction_numerical_grids_and_timestep`
- [two_blackholes_collide.py](../../nrpy/examples/two_blackholes_collide.py) - `project_name`, `CoordSystem`, `IDtype`

## See Also

- [Infrastructures](index.md)
- [First Wave Equation Run](../examples/first-wave-equation-run.md)
- [Black Hole Evolution](../examples/black-hole-evolution.md)
- [Loop Kernel And Device Helpers](../core/helpers/loop-kernel-and-device-helpers.md)
- [Parallel Codegen Orchestration](../core/helpers/parallel-codegen-orchestration.md)
