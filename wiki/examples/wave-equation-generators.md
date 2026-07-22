# Wave Equation Generators

> Compare the Cartesian, curvilinear, and multicoordinate wave-equation example generators. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [Examples](index.md)

## Summary

Use `python -m nrpy.examples.wave_equation_cartesian` for the minimal first
standalone run. It writes `project/wave_equation_cartesian/`, generates C
source, builds with `make`, and runs as `./wave_equation_cartesian`.

Use `python -m nrpy.examples.wave_equation_curvilinear` when you need a
single reference-metric coordinate system, curvilinear boundary conditions,
nearest and volume diagnostics, checkpointing, optional CUDA generation, or
reference-metric precompute toggles. It writes
`project/wave_equation_curvilinear/` and runs as
`./wave_equation_curvilinear` after `make`.

Use `python -m nrpy.examples.wave_equation_multicoordinates` when you need one
generated project with several registered coordinate systems. It writes
`project/wave_equation_multicoordinates/` and runs as
`./wave_equation_multicoordinates` after `make`.

## Detail

All three examples are standalone BHaH generators. Each source file removes its
own `project/<project_name>/` directory before generation, writes a default
`<project_name>.par`, registers `convergence_factor` as a generated executable
command-line input, emits `BHaH_defines.h`, generated C or CUDA sources,
function prototypes, and a Makefile, then prints the same build shape:

```bash
python -m nrpy.examples.wave_equation_cartesian
cd project/wave_equation_cartesian
make
./wave_equation_cartesian
```

Preserve any wanted prior project output before generation. Configured
Ubuntu/macOS CI generates and builds all three default OpenMP projects, but
does not run their executables or compare diagnostic errors. The local helper
also configures CUDA generation/build for curvilinear and multicoordinate
projects only; it performs no GPU runtime or numerical check.

For the other two workflows, replace `wave_equation_cartesian` with
`wave_equation_curvilinear` or `wave_equation_multicoordinates`. From a
checkout without editable install, follow [Build And Run](../architecture/build-and-run.md)
and set `PYTHONPATH` from the repository root before invoking the modules.

The Cartesian generator has no generator-time CLI parser. Its source fixes
`project_name = "wave_equation_cartesian"`, `fp_type = "double"`,
`WaveType = "SphericalGaussian"`, `MoL_method = "RK4"`, `fd_order = 4`, and a
cell-centered Cartesian grid with `Nxx0`, `Nxx1`, and `Nxx2` defaults of 64.
The generated executable accepts `convergence_factor`; the numerical-grid
setup multiplies each `Nxx` by this value before computing spacings, and
diagnostics write center-point relative errors to
`out0d-conv_factor%.2f.txt`.

The Cartesian workflow is the first-run path because the generator owns its
grid setup and manually defined boundary condition path. Initial data calls the
Cartesian exact solution at each point. RHS code uses finite-difference
codegen over the interior loop. Boundary conditions are quadratic polynomial
extrapolation on all six cube faces through `apply_bcs`, not curvilinear
radiation boundary conditions.

The curvilinear generator adds generator-time flags:

```bash
python -m nrpy.examples.wave_equation_curvilinear \
  --floating_point_precision double \
  --cuda \
  --disable_intrinsics \
  --disable_rfm_precompute
```

`--floating_point_precision` lowercases into `fp_type`. `--cuda` changes the
BHaH `parallelization` parameter from default `openmp` to `cuda`.
`--disable_intrinsics` turns off intrinsic-enabled RHS generation.
`--disable_rfm_precompute` turns off reference-metric precompute unless the
chosen coordinate system starts with `GeneralRFM`, where source raises an
error because that family requires precompute. The default source coordinate
choice is `CoordSystem = "SinhCylindrical"`, with `set_of_CoordSystems`
containing that one value and `Nxx_dict["SinhCylindrical"] = [64, 2, 64]`.
For a spherical Gaussian in this default cylindrical case, the source sets
`symmetry_axes` to `"1"`; the same file also contains a spherical branch that
would set `"12"` when a spherical coordinate choice has both angular grid
sizes equal to 2.

The curvilinear workflow registers reference-metric numerical grids with
curvilinear boundary support, per-coordinate RHS functions, Cartesian-to-grid
and grid-to-Cartesian transforms, nearest diagnostics, volume-integration
diagnostics, diagnostic gridfunction setup, checkpointing every 50.0 by
default, and a progress indicator. It registers curvilinear boundary
conditions with second-order radiation stencils. The generated RHS wrapper
applies `apply_bcs_outerradiation_and_inner` when `outer_bc_type` is
`radiation`, and the post-RHS hook applies `apply_bcs_outerextrap_and_inner`
when `outer_bc_type` is `extrapolation`.

The multicoordinate generator has the same generator-time CLI flag shapes as
the curvilinear generator:

```bash
python -m nrpy.examples.wave_equation_multicoordinates \
  --cuda \
  --floating_point_precision double \
  --disable_intrinsics \
  --disable_rfm_precompute
```

Its coordinate source differs: it sets
`set_of_CoordSystems = {"Spherical", "SinhSpherical", "Cartesian", "SinhCartesian"}`
and sets `CoordSystem_to_register_CodeParameters` to `"All"` so code
parameters are registered for all coordinate systems instead of a single
default. Its `Nxx_dict` gives `[64, 2, 2]` for `Spherical` and
`SinhSpherical`, and `[64, 64, 64]` for `Cartesian` and `SinhCartesian`.
Unlike the curvilinear file, this source does not set `symmetry_axes`.

The multicoordinate workflow registers reference-metric numerical grids,
curvilinear boundary conditions, one RHS per coordinate system, one
`xx_to_Cart` wrapper per coordinate system, nearest diagnostics, diagnostic
gridfunction setup, checkpointing every 2.0 by default, and a progress
indicator. Volume-integration diagnostics are disabled in this example. Its
radiation and extrapolation boundary hooks mirror the curvilinear generator.

Reference-metric precompute belongs to the curvilinear and multicoordinate
workflows, not the Cartesian first-run generator. Both reference-metric
generators default to precompute enabled and call
`rfm_precompute.register_CFunctions_rfm_precompute` when it remains enabled.
When precompute is enabled, generated RHS calls receive `rfmstruct` and
`auxevol_gfs`; when disabled, generated RHS calls receive raw coordinate
arrays `xx`. The Cartesian generator emits `BHaH_defines.h` with
`enable_rfm_precompute=False`.

CUDA is also limited to the curvilinear and multicoordinate generators. With
`--cuda`, source sets `parallelization = "cuda"`, registers CUDA utility
functions, emits device headers through `output_device_headers`, copies CUDA
intrinsics where the source selects them, asks the Makefile helper for `nvcc`,
and writes generated source with `.cu` extension. Without `--cuda`, those
generators use OpenMP/default C generation and `.c` source files. The Cartesian
generator always follows the C path and copies SIMD intrinsics only when its
source-level `enable_simd` constant is true.

Generated files under `project/**` are build products. Use the Python
generator files as source evidence for behavior, and use
[Generated Output Boundaries](../architecture/generated-output-boundaries.md)
and [Lifecycle And Project Assembly](../infrastructures/bhah/lifecycle-and-project-assembly.md)
for project ownership, build, and generated-source boundaries.

## Sources

- [wave_equation_cartesian.py](../../nrpy/examples/wave_equation_cartesian.py) - `project_name`, `register_CFunction_numerical_grids_and_timestep_setup`, `register_CFunction_diagnostics`, `register_CFunction_rhs_eval`, `register_CFunction_apply_bcs`
- [wave_equation_cartesian.py](../../nrpy/examples/wave_equation_cartesian.py) - `BHaH.cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser`, `BHaH.BHaH_defines_h.output_BHaH_defines_h`, `BHaH.Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile`
- [wave_equation_curvilinear.py](../../nrpy/examples/wave_equation_curvilinear.py) - `parser`, `project_name`, `CoordSystem`, `Nxx_dict`, `set_of_CoordSystems`, `symmetry_axes`
- [wave_equation_curvilinear.py](../../nrpy/examples/wave_equation_curvilinear.py) - `BHaH.numerical_grids_and_timestep.register_CFunctions`, `BHaH.rfm_precompute.register_CFunctions_rfm_precompute`, `BHaH.CurviBoundaryConditions.register_all.register_C_functions`, `BHaH.Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile`
- [wave_equation_multicoordinates.py](../../nrpy/examples/wave_equation_multicoordinates.py) - `parser`, `project_name`, `set_of_CoordSystems`, `Nxx_dict`, `CoordSystem_to_register_CodeParameters`
- [wave_equation_multicoordinates.py](../../nrpy/examples/wave_equation_multicoordinates.py) - `BHaH.numerical_grids_and_timestep.register_CFunctions`, `BHaH.rfm_precompute.register_CFunctions_rfm_precompute`, `BHaH.CurviBoundaryConditions.register_all.register_C_functions`, `BHaH.Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile`

## See Also

- [Examples](index.md)
- [First Wave Equation Run](first-wave-equation-run.md)
- [Build And Run](../architecture/build-and-run.md)
- [Generated Output Boundaries](../architecture/generated-output-boundaries.md)
- [Lifecycle And Project Assembly](../infrastructures/bhah/lifecycle-and-project-assembly.md)
