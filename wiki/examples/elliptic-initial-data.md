# Elliptic Initial Data

> Compare standalone NRPyElliptic and Charm++/superB conformally flat initial-data workflows. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Examples](index.md)

## Summary

NRPy has two conformally flat elliptic initial-data examples with the same
physical goal but different runtime targets. `python -m
nrpy.examples.nrpyelliptic_conformally_flat` generates a standalone BHaH
project under `project/nrpyelliptic_conformally_flat/`. `python -m
nrpy.examples.superB_nrpyelliptic_conformally_flat` generates a Charm++/superB
project under `project/superB_nrpyelliptic_conformally_flat/`.

Both examples solve the NRPyElliptic hyperbolic-relaxation system for
conformally flat binary-black-hole initial data. This page is the example leaf
for the workflow comparison; it does not create a separate non-example
NRPyElliptic infrastructure page.

## Detail

Use these commands side by side:

| Target | Generate | Project | Build | Run |
| --- | --- | --- | --- | --- |
| Standalone BHaH | `python -m nrpy.examples.nrpyelliptic_conformally_flat` | `project/nrpyelliptic_conformally_flat/` | `cd project/nrpyelliptic_conformally_flat && make` | `./nrpyelliptic_conformally_flat` |
| Charm++/superB | `python -m nrpy.examples.superB_nrpyelliptic_conformally_flat` | `project/superB_nrpyelliptic_conformally_flat/` | `cd project/superB_nrpyelliptic_conformally_flat && make` | `./charmrun +p4 ./superB_nrpyelliptic_conformally_flat` |

The standalone generator targets the BHaH infrastructure and can choose
OpenMP-style C output or CUDA output. It accepts `--cuda` and
`--floating_point_precision`, defaults to double precision, uses OpenMP when
CUDA is not requested, and writes a C or CUDA-ready standalone project with a
Makefile and parameter file. The superB elliptic generator accepts
`--floating_point_precision`, defaults to double precision, and produces a
Charm++ project that builds with `charmc`; it does not expose the standalone
`--cuda` switch.

Both examples use RK4 Method of Lines timestepping as the relaxation driver,
treat `t_final = grid_physical_size` as effectively unused by NRPyElliptic,
cap relaxation at `nn_max = 10000`, set `eta_damping = 11.0`,
`MINIMUM_GLOBAL_WAVESPEED = 0.7`, `CFL_FACTOR = 1.0`, tenth-order finite
differences, sixth-order radiation boundary stencils, outgoing-radiation
boundaries, SIMD intrinsics where supported, and a `gw150914` default data
preset with alternate `axisymmetric` and `single_puncture` branches in the
source.

The residual stop rule is precision dependent in both examples:
`log10_residual_tolerance = -11.2` for double precision and `-6.5` for float.
Progress output reports `nn / nn_max` and current-versus-target
`log10(residual)`. Both generated command-line parsers accept
`convergence_factor`. The standalone project also checkpoints during the
relaxation run, writes a forced checkpoint when the stop condition is reached,
and defaults `checkpoint_every` to `50.0`. The superB elliptic example sets its
initial-data checkpoint hook off, then stops by calling `stop_conditions_check`
inside a Charm++ serial block and finishing through `mainProxy.done()`.

The coordinate defaults differ. The standalone example currently selects
`SinhSymTP`, with `Nxx = [128, 128, 16]`, `AMAX = 1.0e6`, `bScale = 5.0`, and
`SINHWAA = 0.07`. Its source also carries inactive alternatives for
`SymTP`, `SinhCylindricalv2`, and `SinhSpherical`. The superB elliptic example
currently selects `SinhSpherical`, with `Nxx = [128, 64, 16]`, `AMPL = 1.0e6`,
`SINHW = 0.06`, and chare defaults `Nchare0 = 16`, `Nchare1 = 2`,
`Nchare2 = 2` for spherical coordinates.

The assembly path is the main practical difference. The standalone example
registers BHaH NRPyElliptic initial guess, auxiliary-evolution setup, grids,
diagnostics, RHS, residual, stop conditions, curvilinear boundary conditions,
MoL timestepping, checkpointing, headers, generated `main`, cleanup, copied
intrinsics, and a standalone Makefile. The superB example reuses the same
equation-side NRPyElliptic functions, but wraps them in superB numerical grids,
diagnostics, chare communication maps, CurviBCs, superB MoL, PUP routines,
`Main` and `Timestepping` chare files, copied `superB` headers, and a Charm++
Makefile linked with `-module CkIO`.

Validation status is split the same way. The Ubuntu and macOS codegen jobs run
`python -m nrpy.examples.nrpyelliptic_conformally_flat`, build
`project/nrpyelliptic_conformally_flat`, and run `make clean`. The dedicated
`charmpp-validation` workflow generates and builds
`superB_nrpyelliptic_conformally_flat` inside a Charm++ Apptainer image.
Generated `project/**` files are validation products, not documentation sources.

## Sources

- [nrpyelliptic_conformally_flat.py](../../nrpy/examples/nrpyelliptic_conformally_flat.py) - `project_name`, `get_log10_residual_tolerance`, `CoordSystem`, `BHaH.main_c.register_CFunction_main_c`
- [nrpyelliptic_conformally_flat.py](../../nrpy/examples/nrpyelliptic_conformally_flat.py) - `--cuda`, `--floating_point_precision`, `enable_checkpointing`, `cmdline_inputs=["convergence_factor"]`
- [superB_nrpyelliptic_conformally_flat.py](../../nrpy/examples/superB_nrpyelliptic_conformally_flat.py) - `project_name`, `get_log10_residual_tolerance`, `CoordSystem`, `Nchare0`
- [superB_nrpyelliptic_conformally_flat.py](../../nrpy/examples/superB_nrpyelliptic_conformally_flat.py) - `nrpyelliptic_project=True`, `register_CFunction_superB_pup_routines`, `output_CFunctions_function_prototypes_and_construct_Makefile`
- [README.md](../../README.md) - `Good Next Examples`, `Standalone BHaH Generators`, `superB / Charm++ Generators`
- [.github/workflows/main.yml](../../.github/workflows/main.yml) - `codegen-ubuntu`, `codegen-mac`, `charmpp-validation`
- [Generated Output Boundaries](../architecture/generated-output-boundaries.md) - compiled generated-output boundary summary
- [Generated Project CI](../validation/generated-project-ci.md) - compiled CI validation summary

## See Also

- [Examples](index.md)
- [Black Hole Evolution](black-hole-evolution.md)
- [superB Lifecycle And Project Assembly](../infrastructures/superb/lifecycle-and-project-assembly.md)
- [superB Grids, Boundaries, MoL, And Initial Data](../infrastructures/superb/grids-boundaries-mol-and-initial-data.md)
- [Generated Output Boundaries](../architecture/generated-output-boundaries.md)
