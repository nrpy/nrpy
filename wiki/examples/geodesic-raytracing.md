# Geodesic Raytracing

> Explain standalone massive and photon geodesic examples plus batch photon raytracing visualization artifacts. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [Examples](index.md)

## Summary

NRPy has three checked-in geodesic example generators. `massive_single_geodesic_integrator_analytical`
builds a single massive-particle Kerr-Schild Cartesian trajectory and uses GSL's
RKF45 ODE path. `photon_single_geodesic_integrator_analytical` builds a single photon trajectory
with the same analytic spacetime target but uses the split-pipeline photon RKF45
kernels directly. `photon_batch_geodesic_integrator_analytical` builds a tiled photon
raytracing project, defaults to OpenMP, can generate CUDA code with `--cuda`,
honors `--outdir`, and writes per-tile light-blueprint ZIP artifacts for the
lensed-image renderer and diagnostic scripts.

All generation, build, executable, trajectory, and rendering commands on this
page are manual/source-supported. Neither GitHub workflow nor the local full-CI
helper invokes these three generators. No runtime or numerical result was
reproduced during this KB audit.

## Detail

From a repository checkout, run generators from the repository root. If NRPy is
not installed editable, first set `PYTHONPATH` as described by
[Build And Run](../architecture/build-and-run.md):

```bash
export PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}."
```

For a massive single-ray run:

```bash
python -m nrpy.examples.massive_single_geodesic_integrator_analytical
cd project/massive_single_geodesic_integrator_analytical
make
./massive_single_geodesic_integrator_analytical
python3 visualize_trajectory.py --particle_type Massive
```

The generated C program initializes one massive test particle in
`KerrSchild_Cartesian`, computes `u^t` from the massive normalization
constraint, integrates an eight-component state through `gsl_odeiv2_step_rkf45`,
writes `trajectory.txt`, then reports final normalization and conserved-quantity
errors. Its generated Makefile uses `$(shell gsl-config --cflags)` and
`$(shell gsl-config --libs)`, so GSL is a build dependency sourced from the
generator. Official GSL [Using the Library](https://www.gnu.org/software/gsl/doc/html/usage.html)
headings `Compiling and Linking` and `Linking programs with the library`
describe the external header/linker prerequisite; exact `gsl-config` use is
NRPy generator behavior.

For a photon single-ray run:

```bash
python -m nrpy.examples.photon_single_geodesic_integrator_analytical
cd project/photon_single_geodesic_integrator_analytical
make
./photon_single_geodesic_integrator_analytical
python3 visualize_trajectory.py --particle_type Photon
```

The single-photon path also targets `KerrSchild_Cartesian`, but it does not use
GSL. It allocates one Structure-of-Arrays photon state, computes photon `p_t`
with `p0_reverse_kernel`, runs the split RKF45 pipeline through
`interpolation_kernel`, `calculate_ode_rhs_kernel`, `rkf45_stage_update`, and
`rkf45_finalize_and_control`, writes `trajectory.txt`, and reports null
normalization and conserved-quantity errors. Both single-ray generators copy
`visualize_trajectory.py` into the generated project and print
`pip install matplotlib numpy`; those Python visualization dependencies are
source-limited to the checked-in script imports and generator message.

`trajectory.txt` is the handoff artifact for the single-ray visualization. The
massive file header is `# tau t x y z u^t u^x u^y u^z`; the photon file header
is `# lambda t x y z p_t p_x p_y p_z aux`. `visualize_trajectory.py` loads that
text file with NumPy, validates it is present and non-empty, prints initial and
final spatial positions plus the accumulated affine or proper parameter, then
saves `trajectory_plot.png` and displays a Matplotlib 3D path with an approximate
`r=2M` horizon.

For a batch lensed-image run with default OpenMP output under `project/`:

```bash
python -m nrpy.examples.photon_batch_geodesic_integrator_analytical
cd project/photon_batch_geodesic_integrator_analytical
make
./photon_batch_geodesic_integrator_analytical
python3 visualize_lensed_image.py --source_r_min 6.0 --source_r_max 20.0 --window_width 1.0 --window_height 1.0 --window_tiles_width 2 --window_tiles_height 2 --pixel_width 600
python3 blueprint_analysis.py --window_tiles_width 2 --window_tiles_height 2 --window_width 1.0 --window_height 1.0
```

For a CUDA batch project or a different parent output directory:

```bash
python -m nrpy.examples.photon_batch_geodesic_integrator_analytical --cuda --outdir /tmp/nrpy-geodesics
cd /tmp/nrpy-geodesics/photon_batch_geodesic_integrator_analytical
make
./photon_batch_geodesic_integrator_analytical
```

The batch generator derives `project_dir` from `--outdir` plus
`photon_batch_geodesic_integrator_analytical`. It sets `parallelization` to `openmp` unless
`--cuda` is present. OpenMP uses `gcc`, `-fopenmp`, C sources, and a default
`2x2` tile grid with `scan_density = 500`; CUDA uses `nvcc`, `-lcudart`,
`-DUSE_GPU`, `.cu` sources, copied `cuda_intrinsics.h`, a default `1x1` tile
grid, and `scan_density = 1000`. The CUDA runtime/toolchain facts are limited to
the generator's compiler, flags, and copied helper declarations. NVIDIA's
official [NVCC guide](https://docs.nvidia.com/cuda/cuda-programming-guide/02-basics/nvcc.html)
heading `NVCC: The NVIDIA CUDA Compiler` establishes that `nvcc` belongs to the
CUDA Toolkit; this page does not claim any GPU model/toolkit version was tested.

The batch executable produces tiled `light_blueprint_XX_YY.zip` files in its
project directory. Those ZIPs are generated artifacts, not KB sources. The
Python schema in `blueprint_config_and_schema.py` defines `BLUEPRINT_DTYPE` for
the binary records: termination type, camera-window coordinates, source-plane
coordinates, final sphere angles, and affine/time values at window and source
intersections. The schema file explicitly says the dtype must match the C
`blueprint_data_t` layout and that termination enum constants must stay
synchronized with generated C headers; this is the main schema synchronization
risk.

`visualize_lensed_image.py` expects those per-tile ZIPs next to the script,
constructs their names from `window_tiles_width` and `window_tiles_height`,
downloads `noirlab2430b.tif` at runtime from the URL embedded in the script if
the texture is missing, generates an accretion-disk texture, and calls
`render_lensed_image.generate_static_lensed_image`. The renderer streams binary
records from the ZIPs with `BLUEPRINT_DTYPE`, maps source-plane hits to the disk
texture and sphere escapes to the starmap, processes tiles in a bounded process
pool, and writes the final PNG. `render_lensed_image.py` has a no-op `numba`
fallback if import fails, but the batch generator still prints
`pip install matplotlib numpy numba Pillow`; that dependency statement is
source-limited to the checked-in files.

`blueprint_analysis.py` is the diagnostic script for the same ZIP artifacts. It
streams each tile, counts raw termination enums, compares them to
`blueprint_config_and_schema.py`, reports camera-window in-view statistics,
prints early records, and displays heatmaps for camera-window, source-plane, and
celestial-sphere coordinates. Its warning text directs maintainers to update the
schema file when raw enum values do not match current Python constants.

This standalone batch raytracer is related to, but distinct from, the
evolution-time raytracing export enabled by
`python -m nrpy.examples.two_blackholes_collide --raytracing-outputs`. The
black-hole evolution option writes Cartesian metric and Christoffel data on
diagnostic output steps and is currently rejected for CUDA builds. See
[Black Hole Evolution](black-hole-evolution.md) and
[Geodesics And Raytracing Runtime](../infrastructures/bhah/geodesics-and-raytracing-runtime.md)
for that context.

## Sources

- [massive_single_geodesic_integrator_analytical.py](../../nrpy/examples/massive_single_geodesic_integrator_analytical.py) - `project_name`, `main_single`, `single_integrator_analytical`, `gsl-config`; official GSL [Using the Library](https://www.gnu.org/software/gsl/doc/html/usage.html) - `Compiling and Linking`
- [photon_single_geodesic_integrator_analytical.py](../../nrpy/examples/photon_single_geodesic_integrator_analytical.py) - `project_name`, `main_single`, `p0_reverse_kernel`, `rkf45_stage_update`, `trajectory.txt`
- [photon_batch_geodesic_integrator_analytical.py](../../nrpy/examples/photon_batch_geodesic_integrator_analytical.py) - `--outdir`, `--cuda`, `parallelization_mode`, `vis_command`, `blueprint_command`; official NVIDIA [NVCC guide](https://docs.nvidia.com/cuda/cuda-programming-guide/02-basics/nvcc.html) - `NVCC: The NVIDIA CUDA Compiler`
- [visualize_trajectory.py](../../nrpy/examples/geodesic_visualizations/visualize_trajectory.py) - `visualize_trajectory`, `plot_trajectory`
- [blueprint_config_and_schema.py](../../nrpy/examples/geodesic_visualizations/blueprint_config_and_schema.py) - `BLUEPRINT_DTYPE`, `TERM_SPHERE`, `TERM_SOURCE_PLANE`
- [render_lensed_image.py](../../nrpy/examples/geodesic_visualizations/render_lensed_image.py) - `generate_static_lensed_image`, `_process_blueprint_tile`, `_load_texture`
- [visualize_lensed_image.py](../../nrpy/examples/geodesic_visualizations/visualize_lensed_image.py) - `main`, `light_blueprint_{i:02d}_{j:02d}.zip`, `urlretrieve`
- [blueprint_analysis.py](../../nrpy/examples/geodesic_visualizations/blueprint_analysis.py) - `diagnose_blueprint`, `plot_heatmaps`
- [two_blackholes_collide.py](../../nrpy/examples/two_blackholes_collide.py) - `--raytracing-outputs`

## See Also

- [Examples](index.md)
- [Black Hole Evolution](black-hole-evolution.md)
- [Geodesics](../equations/general-relativity/geodesics.md)
- [Geodesics And Raytracing Runtime](../infrastructures/bhah/geodesics-and-raytracing-runtime.md)
