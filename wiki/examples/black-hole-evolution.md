# Black Hole Evolution

> Route the lightweight single-patch BH@H black-hole evolution example to its generated-application context. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [Examples](index.md)

## Summary

`python -m nrpy.examples.two_blackholes_collide` is the compact public
single-patch BH@H numerical-relativity example. It evolves Brill-Lindquist
binary black-hole initial data in a spherical curvilinear setup and writes a
buildable project under `project/two_blackholes_collide/`.

## Detail

The README presents this as the closest public example of the single-patch BH@H
workflow and notes that it is intentionally lightweight. The generator itself
sets the BHaH infrastructure, defaults to OpenMP unless `--cuda` is supplied,
uses `double` precision by default, and rejects raytracing output for CUDA.

The source configuration uses `CoordSystem = "Spherical"`,
`IDtype = "BrillLindquist"`, Cartesian initial-data coordinates, `OnePlusLog`
lapse evolution, `GammaDriving2ndOrder_Covariant` shift evolution, RK4 Method
of Lines timestepping, fourth-order finite differences, reference-metric
precomputation, intrinsic support, curvilinear boundary conditions, and outgoing
radiation boundaries. The default masses are equal, with the two black holes
placed on opposite sides of the origin along the z direction.

Its registration flow is the GR-facing BHaH lifecycle in one file. It registers
Brill-Lindquist initial data, grid and timestep setup, diagnostic gridfunctions,
nearest and volume diagnostics, reference-metric precompute, BSSN RHS
evaluation, a separate Ricci evaluation path when enabled, determinant
enforcement, constraint evaluation, curvilinear boundary conditions, Method of
Lines timestepping, Cartesian/coordinate transforms, basis transforms, coordinate
wrapper functions, headers, the parameter parser, the generated `main`, grid
cleanup, copied intrinsics, and a Makefile.

When OpenMP and double precision are used, the generator also enables BHaHAHA
integration and sets apparent-horizon defaults. The generated command-line
parser accepts `convergence_factor`, matching the broader standalone BHaH
example pattern.

For the rest of the standalone BHaH numerical-relativity family, use
[Standalone GR/BHaH](standalone-gr-bhah.md). For raytracing output and
standalone geodesic generators, use [Geodesic Raytracing](geodesic-raytracing.md).
For BHaHAHA library generation behind horizon-enabled examples, use
[Apparent Horizon Library](apparent-horizon-library.md).

## Sources

- [README.md](../../README.md) - `Lightweight Single-Patch BH@H Example`, `Standalone BHaH Generators`
- [two_blackholes_collide.py](../../nrpy/examples/two_blackholes_collide.py) - `project_name`, `CoordSystem`, `IDtype`
- [two_blackholes_collide.py](../../nrpy/examples/two_blackholes_collide.py) - `BHaH.general_relativity.rhs_eval.register_CFunction_rhs_eval`
- [two_blackholes_collide.py](../../nrpy/examples/two_blackholes_collide.py) - `BHaH.main_c.register_CFunction_main_c`
- [two_blackholes_collide.py](../../nrpy/examples/two_blackholes_collide.py) - `BHaH.Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile`

## See Also

- [Examples](index.md)
- [First Wave Equation Run](first-wave-equation-run.md)
- [Standalone GR/BHaH](standalone-gr-bhah.md)
- [Geodesic Raytracing](geodesic-raytracing.md)
- [Apparent Horizon Library](apparent-horizon-library.md)
- [BHaH Lifecycle](../infrastructures/bhah-lifecycle.md)
