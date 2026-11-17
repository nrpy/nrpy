# Wave Equation Generators

> Compare the Cartesian, curvilinear, and multicoordinate BHaH wave-equation example generators. · Status: confirmed
> Up: [Examples](index.md)

## Summary

`wave_equation_cartesian.py` provides the minimal zero-flag standalone BHaH build and run, `wave_equation_curvilinear.py` adds single-coordinate reference-metric precomputation, curvilinear radiation/extrapolation boundaries, volume diagnostics, checkpointing, and `--cuda`, and `wave_equation_multicoordinates.py` emits multiple coordinate systems inside one executable. Each generator recreates `project/<project_name>/`, emits `<project_name>.par`, and accepts runtime `convergence_factor` on the executable CLI.

## Detail

### First Standalone Run (`wave_equation_cartesian.py`)

From the repository root (see [Build And Run](../architecture/build-and-run.md) for `PYTHONPATH=.` and `Linux`/`Darwin` executable prerequisites):

```bash
PYTHONPATH=. python3 nrpy/examples/wave_equation_cartesian.py
cd project/wave_equation_cartesian
make
./wave_equation_cartesian 2.0
```

`wave_equation_cartesian.py` takes no generator-time CLI flags: it fixes `fp_type = "double"`, `WaveType = "SphericalGaussian"`, `MoL_method = "RK4"`, `fd_order = 4`, and a cell-centered `64^3` Cartesian grid (`enable_rfm_precompute=False`, C/OpenMP only, copying SIMD headers when `enable_simd` is true). It registers exact-solution, initial-data, grid-setup, diagnostic, interior `simple_loop` FD RHS, and quadratic face-extrapolation `apply_bcs` functions; passing `convergence_factor` (`2.0` above) scales `Nxx` and writes center-point relative error to `out0d-conv_factor2.00.txt`.

Evidence (descriptive behavior): `wave_equation_cartesian.py`; `wave_equation_curvilinear.py`; `wave_equation_multicoordinates.py`.

### Curvilinear And Multicoordinate Generators (`--cuda`, `rfm_precompute`)

Both curvilinear generators accept identical generator-time flags:

```bash
PYTHONPATH=. python3 nrpy/examples/wave_equation_curvilinear.py \
  --floating_point_precision double --cuda --disable_intrinsics --disable_rfm_precompute
```

- **`wave_equation_curvilinear.py`**: Defaults to `CoordSystem = "SinhCylindrical"` (`Nxx_dict["SinhCylindrical"] = [64, 2, 64]`, setting `symmetry_axes = "1"` for `SphericalGaussian`, or `"12"` on `Spherical` with `[*, 2, 2]`). Registers `rfm` grids, `CurviBoundaryConditions` (`fd_order = 2` Sommerfeld `apply_bcs_outerradiation_and_inner` inside RHS or `apply_bcs_outerextrap_and_inner` post-RHS), coordinate transforms, nearest + volume diagnostics, and checkpointing (`50.0` interval).
- **`wave_equation_multicoordinates.py`**: Sets `set_of_CoordSystems = {"Spherical", "SinhSpherical", "Cartesian", "SinhCartesian"}` (`[64, 2, 2]` on spherical, `[64, 64, 64]` on Cartesian, no `symmetry_axes`) and `CoordSystem_to_register_CodeParameters = "All"`, emitting one RHS and `xx_to_Cart` per coordinate system, nearest diagnostics (volume integration disabled), and checkpointing (`2.0` interval).
- **Precompute and CUDA**: Both default to `enable_rfm_precompute=True` (`rfmstruct` + `auxevol_gfs` passed to RHS; `--disable_rfm_precompute` passes raw `xx` and raises on `GeneralRFM*`). `--cuda` sets `parallelization = "cuda"`, emits device headers, `.cu` kernels, `cuda_intrinsics.h`, and an `nvcc` `Makefile`. Ubuntu/macOS CI compiles default OpenMP builds for all three generators (and `.github/full_nrpy_local_ci.sh` compiles CUDA for curvilinear/multicoordinates) without running the binaries.

## Sources

- [wave_equation_cartesian.py](../../nrpy/examples/wave_equation_cartesian.py) - `project_name`, `register_CFunction_numerical_grids_and_timestep_setup`, `register_CFunction_diagnostics`, `register_CFunction_rhs_eval`, `register_CFunction_apply_bcs`
- [wave_equation_curvilinear.py](../../nrpy/examples/wave_equation_curvilinear.py) - `parser`, `CoordSystem`, `Nxx_dict`, `symmetry_axes`, `rfm_precompute`, `CurviBoundaryConditions`
- [wave_equation_multicoordinates.py](../../nrpy/examples/wave_equation_multicoordinates.py) - `parser`, `set_of_CoordSystems`, `Nxx_dict`, `CoordSystem_to_register_CodeParameters`

## See Also

- Depends on: [Wave Equation](../equations/wave-equation.md)
- Depends on: [Build And Run](../architecture/build-and-run.md)
- Depends on: [BHaH Project Assembly](../infrastructures/bhah/lifecycle-and-project-assembly.md)
