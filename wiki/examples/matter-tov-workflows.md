# Matter TOV Workflows

> Route neutron-star, static-fluid, GRoovy, and MANGA-adjacent TOV workflows. · Status: confirmed · Last reconciled: 06-30-2026
> Up: [Examples](index.md)

## Summary

The matter-facing TOV examples split into three runnable standalone BHaH
projects and one library generator. `tovola_neutron_star` solves TOV initial
data with TOVola, interpolates it through `TOVola_interp`, and validates
constraints without time evolution. `hydro_without_hydro` keeps the fluid
static while evolving BSSN spacetime variables with matter coupling.
`groovy_TOV_BSSN__Spherical` evolves a GRHD TOV star through GRoovy and
requires a local GRHayL build under the generated project. `manga_bhah_lib`
generates `project/bhah_lib/` as a library target for MANGA-facing integration,
not a normal runnable executable workflow.

## Detail

Run repository-checkout examples as modules after installing the package or
after putting `.` on `PYTHONPATH` from the repository root. The command shapes
are:

```bash
python -m nrpy.examples.tovola_neutron_star
python -m nrpy.examples.hydro_without_hydro
python -m nrpy.examples.groovy_TOV_BSSN
python -m nrpy.examples.manga_bhah_lib
```

`tovola_neutron_star` writes `project/tovola_neutron_star/`. It sets
`CoordSystem = "SinhSpherical"`, `IDtype = "TOVola_interp"`, and
`IDCoordSystem = "Spherical"`. The generator registers both `TOVola_solve` and
`TOVola_interp`, stores TOV arrays in the initial-data persist struct, calls
`TOVola_solve(commondata, &ID_persist)` before initial-data interpolation, and
frees the persisted radius, density, pressure, mass, metric, and isotropic
radius arrays afterward. It enables checkpoint-aware initial data, matter
source support through `enable_T4munu=True`, Ricci and constraint evaluation,
nearest and volume diagnostics, RK4 Method of Lines scaffolding with an empty
RHS, curvilinear boundary infrastructure, coordinate conversion helpers, and
default checkpoint registration with `default_checkpoint_every = 2.0`.

The generated `tovola_neutron_star` target is buildable with GSL flags from
`gsl-config`:

```bash
cd project/tovola_neutron_star
make
./tovola_neutron_star
```

Its expected validation milestone is generation plus `make` and `make clean`;
Ubuntu and macOS CI both install GSL and exercise that path. Runtime defaults
include `t_final = 1.0e-10`, so this is an initial-data and constraints-oriented
workflow rather than a long evolution.

`hydro_without_hydro` writes `project/hydro_without_hydro/`. It defaults to
OpenMP and `double`, accepts `--cuda` to generate CUDA source and use `nvcc`,
and accepts `--floating_point_precision`. Its default coordinate and ID choices
are `CoordSystem = "Spherical"` and `IDtype = "TOVola_interp"` with spherical
TOV initial-data coordinates. It uses `HarmonicSlicing`, frozen shift, RK4,
fourth-order finite differences, outgoing radiation boundaries, checkpointing
every `2.0`, and `t_final = 7.5` from `grid_physical_size`.

This example is "hydro without hydro" because TOVola supplies matter initial
data and `enable_T4munu=True`, while the RHS evolves the spacetime through
BSSN-style `rhs_eval` and constraint paths coupled to `T4munu`; there is no
GRHD flux update in its RHS. With separate Ricci enabled, the RHS first calls
`Ricci_eval`, then `rhs_eval`, then radiation boundary conditions. The
generated target also uses `gsl-config` for C flags and libraries:

```bash
cd project/hydro_without_hydro
make
./hydro_without_hydro
```

CI validates the default OpenMP path on Ubuntu and macOS by generating the
project, building it, and running `make clean`. CUDA is a source-supported
option from the example, but the checked-in CI evidence covers the default path.

`groovy_TOV_BSSN.py` writes `project/groovy_TOV_BSSN__Spherical/` by default.
It sets `CoordSystem = "Spherical"`, `IDType = "TOV"`, evolving spacetime on,
evolving temperature off, `OnePlusLog` lapse, covariant second-order gamma
driver shift, RK4, fourth-order finite differences, Kreiss-Oliger dissipation,
and `t_final = 2300.0`. TOV setup uses a hybrid EOS configuration with central
baryon density `1.28e-3`, `neos = 1`, `Gamma_poly_tab = 2`, and
`K_poly_tab0 = 100`, then adjusts TOVola code-parameter defaults before
registering hybrid-EOS TOV initial data and an ADM-to-BSSN reader/converter.

The GRoovy workflow is a real GRHD evolution path. It registers GRHD
gridfunctions, reconstruction, HLL fluxes, flux divergences, source terms,
primitive/conservative conversion, stress-energy tensor computation, copy and
outflow boundary conditions, basis transforms, and GRoovy diagnostics. Its
spacetime RHS calls `Ricci_eval` and `rhs_eval`; its GRHD RHS calls
`grhd_rhs_eval` with GRHayL parameter and EOS structs; after each RHS stage it
enforces the BSSN determinant condition, converts conservatives back to
primitives, applies hydrodynamic boundary handling, and recomputes
stress-energy.

GRHayL is a hard generated-project build dependency for the GRoovy TOV example.
The generator adds `ghl.h`, reconstruction, atmosphere, and con2prim headers,
passes `GRHayL/include/ghl/` as an include directory, links `-lghl` from
`GRHayL/lib`, and embeds an rpath to that generated-project-local library
directory. When run as `__main__`, it clones `https://github.com/GRHayL/GRHayL.git`
into `project/groovy_TOV_BSSN__Spherical/GRHayL`, configures it with
`--prefix=./ --buildtype=opt --disable-hdf5`, runs parallel `make`, and runs
`make install` before instructing the user to build the BHaH executable:

```bash
cd project/groovy_TOV_BSSN__Spherical
make
./groovy_TOV_BSSN__Spherical
```

This source-managed GRHayL clone means the generator itself needs network access
unless the dependency is already provided in the expected generated-project
layout. Unlike `tovola_neutron_star` and `hydro_without_hydro`, no checked-in CI
step in the inspected workflow builds `groovy_TOV_BSSN.py`.

`manga_bhah_lib.py` writes `project/bhah_lib/` and must be read as library
generation. It sets `project_name = "bhah_lib"`, uses spherical coordinates,
`IDtype = "TOVola_interp"`, spherical initial-data coordinates, BSSN evolution
registration, `enable_T4munu=True`, checkpointing every `2.0`, and a custom
diagnostic set that includes Hamiltonian-constraint and metric quantities plus
`T4UU00`. It registers `TOVola_solve`, `TOVola_interp`, checkpoint/progress
helpers, grids, diagnostics, BSSN initial data, RHS, Ricci, determinant
enforcement, constraints, curvilinear boundary conditions, RK4 Method of Lines
functions, coordinate transforms, and coordinate-system wrappers.

The library split happens at Makefile generation: `exec_or_library_name` is
`libbhah_lib`, `create_lib=True`, and `register_CFunctions_bhah_lib()` provides
the library-facing surface. The generated project still uses GSL build flags,
but the expected product is a library, not a standalone executable run:

```bash
cd project/bhah_lib
make
```

The source's final print still says to run `./bhah_lib`, but the Makefile call
is explicitly a library build. Treat that print as stale generic messaging for
this example. The checked-in CI workflows contain commented-out Ubuntu and macOS
commands for `python -m nrpy.examples.manga_bhah_lib && (cd project/bhah_lib &&
make && make clean)`, so source and workflow evidence support generation intent
and a disabled/commented CI status, not active CI coverage.

Across these examples, GSL is the common external C build dependency because
their generated Makefiles add `$(shell gsl-config --cflags)` and
`$(shell gsl-config --libs)`. GRoovy adds GRHayL on top of GSL. The runnable
examples produce executable targets named after their project directories;
`manga_bhah_lib` produces a `libbhah_lib` library target and should be kept
separate from executable run instructions.

## Sources

- [tovola_neutron_star.py](../../nrpy/examples/tovola_neutron_star.py) - `project_name`, `CoordSystem`, `IDtype`, `IDCoordSystem`
- [tovola_neutron_star.py](../../nrpy/examples/tovola_neutron_star.py) - `register_CFunction_TOVola_interp`, `register_CFunction_TOVola_solve`, `register_CFunction_initial_data`
- [tovola_neutron_star.py](../../nrpy/examples/tovola_neutron_star.py) - `register_CFunctions`, `register_CFunctions_function_prototypes_and_construct_Makefile`
- [hydro_without_hydro.py](../../nrpy/examples/hydro_without_hydro.py) - `parser`, `parallelization`, `project_name`, `CoordSystem`, `IDtype`
- [hydro_without_hydro.py](../../nrpy/examples/hydro_without_hydro.py) - `enable_T4munu`, `rhs_string`, `output_CFunctions_function_prototypes_and_construct_Makefile`
- [groovy_TOV_BSSN.py](../../nrpy/examples/groovy_TOV_BSSN.py) - `project_name`, `IDType`, `grhayl_setup_str`
- [groovy_TOV_BSSN.py](../../nrpy/examples/groovy_TOV_BSSN.py) - `register_all_grhd_gridfunctions`, `register_CFunction_grhd_rhs_eval`, `post_RHS_string`
- [groovy_TOV_BSSN.py](../../nrpy/examples/groovy_TOV_BSSN.py) - `additional_includes`, `ghl_INC_FLAG`, `repo_url`
- [manga_bhah_lib.py](../../nrpy/examples/manga_bhah_lib.py) - `project_name`, `set_of_CoordSystems`, `IDtype`
- [manga_bhah_lib.py](../../nrpy/examples/manga_bhah_lib.py) - `register_CFunctions_bhah_lib`, `create_lib=True`
- [.github/workflows/main.yml](../../.github/workflows/main.yml) - `codegen-ubuntu`, `codegen-mac`
- [README.md](../../README.md) - `Prerequisites by Workflow`

## See Also

- [Examples](index.md)
- [First Wave Equation Run](first-wave-equation-run.md)
- [Black Hole Evolution](black-hole-evolution.md)
- [Generated Project CI](../validation/generated-project-ci.md)
- [GRoovy GRHD Runtime](../infrastructures/bhah/groovy-grhd-runtime.md)
- [Lifecycle And Project Assembly](../infrastructures/bhah/lifecycle-and-project-assembly.md)
