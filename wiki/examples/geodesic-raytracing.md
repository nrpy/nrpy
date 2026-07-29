# Geodesic Raytracing

> Explain standalone massive and photon geodesic examples plus batch photon raytracing visualization artifacts. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [Examples](index.md)

## Summary

NRPy has five checked-in geodesic example generators. `massive_single_geodesic_integrator_analytical`
builds a single massive-particle Kerr-Schild Cartesian trajectory and uses GSL's
RKF45 ODE path. `photon_single_geodesic_integrator_analytical` builds a single photon trajectory
with the same analytic spacetime target but uses the split-pipeline photon RKF45
kernels directly. `photon_batch_geodesic_integrator_analytical` builds a tiled photon
raytracing project, defaults to OpenMP, can generate CUDA code with `--cuda`,
honors `--outdir`, and writes per-tile light-blueprint binary artifacts for the
lensed-image renderer and diagnostic scripts. The numerical single- and batch-
photon generators use the same photon runtime with a combined numerical-
spacetime `.bin` dataset and currently support CPU/OpenMP
`SinhCylindricalv2n2` interpolation.

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
python -m nrpy.examples.photon_single_geodesic_integrator_analytical \
  --observer-position 0.0 -10.0 0.0 \
  --observer-look-forward 0.0 1.0 0.0 \
  --observer-up 0.0 0.0 1.0 \
  --escape-radius 25.0
cd project/photon_single_geodesic_integrator_analytical
make
./photon_single_geodesic_integrator_analytical
python3 visualize_trajectory.py --particle_type Photon
```

The single-photon path also targets `KerrSchild_Cartesian`, but it does not use
GSL. It allocates one Structure-of-Arrays photon state, constructs the complete
initial momentum from the observer tetrad, and runs the split RKF45 pipeline through
`interpolation_kernel`, `calculate_ode_rhs_kernel`, `rkf45_stage_update`, and
`rkf45_finalize_and_control`, writes `trajectory.txt`, and reports null
normalization and conserved-quantity errors. Both single-ray generators copy
`visualize_trajectory.py` into the generated project and print
`pip install matplotlib numpy`; those Python visualization dependencies are
source-limited to the checked-in script imports and generator message.

Claim evidence:
- Claim: The single-photon analytical generator constructs initial momentum from the observer tetrad, runs the split RKF45 pipeline, writes `trajectory.txt`, and reports normalization and conserved-quantity diagnostics.
- Role: generated-output boundary
- Deciding authority: `nrpy/examples/photon_single_geodesic_integrator_analytical.py` — generator registration and `single_integrator_analytical`
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/set_initial_conditions_kernel.py` — observer-tetrad initialization
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

`trajectory.txt` is the handoff artifact for the single-ray visualization. The
massive file header is `# tau t x y z u^t u^x u^y u^z`; the photon file header
is `# lambda t x y z p^t p^x p^y p^z aux`. `visualize_trajectory.py` loads that
text file with NumPy, validates it is present and non-empty, prints initial and
final spatial positions plus the accumulated affine or proper parameter, then
saves `trajectory_plot.png` and displays a Matplotlib 3D path with an approximate
`r=2M` horizon.

Claim evidence:
- Claim: `trajectory.txt` is the single-ray visualization handoff with the documented state headers, and `visualize_trajectory.py` validates the input before plotting it.
- Role: generated-output boundary
- Deciding authority: `nrpy/examples/geodesic_visualizations/visualize_trajectory.py` — `visualize_trajectory`; `nrpy/examples/photon_single_geodesic_integrator_analytical.py` — trajectory header
- Corroboration: `nrpy/examples/photon_single_geodesic_integrator_numerical.py` — numerical trajectory header
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

For a batch lensed-image run with default OpenMP output under `project/`:

```bash
python -m nrpy.examples.photon_batch_geodesic_integrator_analytical \
  --observer-position 0.0 -10.0 0.0 \
  --observer-look-forward 0.0 1.0 0.0 \
  --observer-up 0.0 0.0 1.0 \
  --observer-fov 1.0 1.0 \
  --scan-density 500 \
  --escape-radius 25.0
cd project/photon_batch_geodesic_integrator_analytical
make
./photon_batch_geodesic_integrator_analytical
python3 visualize_lensed_image.py --terminal-plane-radius 6.0 20.0 --pixel-width 600
python3 blueprint_analysis.py
```

For a CUDA batch project or a different parent output directory:

```bash
python -m nrpy.examples.photon_batch_geodesic_integrator_analytical \
  --cuda --outdir /tmp/nrpy-geodesics \
  --observer-position 0.0 -10.0 0.0 \
  --observer-look-forward 0.0 1.0 0.0 \
  --observer-up 0.0 0.0 1.0 \
  --observer-fov 1.0 1.0 \
  --scan-density 500 \
  --escape-radius 25.0
cd /tmp/nrpy-geodesics/photon_batch_geodesic_integrator_analytical
make
./photon_batch_geodesic_integrator_analytical
```

The batch generator derives `project_dir` from `--outdir` plus
`photon_batch_geodesic_integrator_analytical`. It sets `parallelization` to `openmp` unless
`--cuda` is present. OpenMP uses `gcc`, `-fopenmp`, C sources, and a default
`2x2` tile grid with a width-side `scan_density` of 500; CUDA uses `nvcc`,
`-lcudart`, `-DUSE_GPU`, `.cu` sources, copied `cuda_intrinsics.h`, a default
`1x1` tile grid, and a width-side `scan_density` of 1000. Height-side sampling
is derived internally from the fields of view and tile-grid aspect ratio.
Generated programs expose `tiles_width`, `tiles_height`, and `scan_density` as
angular ray-sampling controls; `--pixel-width` on `visualize_lensed_image.py`
controls final PNG resampling only. The visualization reads tile counts and
`alpha_w`/`alpha_h` from blueprint headers rather than accepting them as CLI
overrides. The CUDA runtime/toolchain facts are limited to
the generator's compiler, flags, and copied helper declarations. NVIDIA's
official [NVCC guide](https://docs.nvidia.com/cuda/cuda-programming-guide/02-basics/nvcc.html)
heading `NVCC: The NVIDIA CUDA Compiler` establishes that `nvcc` belongs to the
CUDA Toolkit; this page does not claim any GPU model/toolkit version was tested.

Claim evidence:
- Claim: The analytical batch generator defaults to OpenMP, selects CUDA with `--cuda`, derives its project directory from `--outdir`, and exposes tile counts and scan density as angular sampling controls while final PNG pixel width remains a visualization setting.
- Role: user-facing commands and interfaces
- Deciding authority: `nrpy/examples/photon_batch_geodesic_integrator_analytical.py` — argument parser and project generation
- Corroboration: `nrpy/examples/geodesic_visualizations/visualize_lensed_image.py` — header-driven tile and aspect-ratio handling
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

The batch executable produces tiled `light_blueprint_XX_YY.bin` files in its
project directory. Those binary files are generated artifacts, not KB sources.
The native same-build blueprint schema is version 6 with 100-byte records. The
header stores tile identity/counts and `alpha_w`/`alpha_h`. Each record stores
nonterminal/terminal-plane diagnostics, final sphere angles, affine/time values,
and normalized `image_width_fraction`/`image_height_fraction` sample coordinates.
The Python
schema in `blueprint_config_and_schema.py` explicitly says the dtype must match the C
`blueprint_data_t` layout and that termination enum constants must stay
synchronized with generated C headers; this is the main schema synchronization
risk.

Claim evidence:
- Claim: Analytical batch artifacts use native same-build blueprint schema version 6 with 100-byte records, and Python `BLUEPRINT_DTYPE` must match the generated C `blueprint_data_t` layout and termination enums.
- Role: generated-output boundary
- Deciding authority: `nrpy/examples/geodesic_visualizations/blueprint_config_and_schema.py` — `BLUEPRINT_SCHEMA_VERSION`, `BLUEPRINT_DTYPE`
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/calculate_and_fill_blueprint_data_universal.py` — generated record layout
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

`visualize_lensed_image.py` expects those per-tile `.bin` files next to the
script and discovers the tile grid and angular aspect ratio from headers. It
does not use plane coordinates to place pixels: the renderer maps each record's
normalized image fractions with floor-to-pixel placement, preserves the
observer's horizontal orientation, and applies the required vertical
raster-row flip. It
downloads `noirlab2430b.tif` at runtime from the URL embedded in the script if
the texture is missing, generates an accretion-disk texture, and calls
`render_lensed_image.generate_static_lensed_image`. The renderer streams binary
records with `BLUEPRINT_DTYPE`, maps terminal-plane hits to the disk
texture and sphere escapes to the starmap, processes tiles in a bounded process
pool, and writes the final PNG. `render_lensed_image.py` has a no-op `numba`
fallback if import fails, but the batch generator still prints
`pip install matplotlib numpy numba Pillow`; that dependency statement is
source-limited to the checked-in files.

Claim evidence:
- Claim: `visualize_lensed_image.py` consumes per-tile binary artifacts, discovers tile geometry from headers, maps normalized image fractions to pixels, and delegates final rendering to `render_lensed_image.py`.
- Role: generated-output boundary
- Deciding authority: `nrpy/examples/geodesic_visualizations/visualize_lensed_image.py` — `main`; `render_lensed_image.py` — `generate_static_lensed_image`
- Corroboration: `nrpy/examples/geodesic_visualizations/blueprint_io.py` — binary header and record readers
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

`blueprint_analysis.py` is the diagnostic script for the same binary artifacts.
It streams each tile, counts raw termination enums, compares them to
`blueprint_config_and_schema.py`, reports nonterminal-plane diagnostic statistics,
prints early records, and displays heatmaps for nonterminal-plane, terminal-plane,
and celestial-sphere coordinates. Its warning text directs maintainers to
update the schema file when raw enum values do not match current Python
constants.

Claim evidence:
- Claim: `blueprint_analysis.py` streams the same binary artifacts, reports termination diagnostics against the configured enum names, and provides nonterminal-plane, terminal-plane, and celestial-sphere heatmaps.
- Role: descriptive behavior
- Deciding authority: `nrpy/examples/geodesic_visualizations/blueprint_analysis.py` — `diagnose_blueprint`, `plot_heatmaps`
- Corroboration: `nrpy/examples/geodesic_visualizations/blueprint_config_and_schema.py` — termination constants and record fields
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

This standalone batch raytracer is related to, but distinct from, the
evolution-time raytracing export enabled by
`python -m nrpy.examples.two_blackholes_collide --raytracing-time`. The
black-hole evolution option writes Cartesian metric and Christoffel data on
diagnostic output steps and is currently rejected for CUDA builds. The optional
`--raytracing-static-christoffels` flag changes only the GammaUDD payload on the
last scheduled output; the default remains nonstatic. See
[Black Hole Evolution](black-hole-evolution.md) and
[Geodesics And Raytracing Runtime](../infrastructures/bhah/geodesics-and-raytracing-runtime.md)
for that context.

Claim evidence:
- Claim: Evolution-time raytracing export is distinct from the standalone batch raytracer; `--raytracing-static-christoffels` is optional, default-off, and changes only the GammaUDD payload on the last scheduled output.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/BHaH/diagnostics/diagnostics.py` — `enable_static_christoffels` scheduling
- Corroboration: `nrpy/infrastructures/BHaH/diagnostics/output_raytracing_data.py` — static GammaUDD selection
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

For a numerical-spacetime dataset, first generate the evolution project and
its slice-combining workflow:

```bash
python -m nrpy.examples.two_blackholes_collide \
  --raytracing-time 30.0 0.4 \
  --raytracing-data-mode g4DD \
  --raytracing-coord-system SinhCylindricalv2n2
cd project/two_blackholes_collide
./run_raytracing_data_pipeline.sh
```

The generated pipeline copies `combine_raytracing_time_slices.py`, writes
`raytracing_run_metadata.json`, runs the evolution, and combines the full
logical-grid slice payload into `project/raytracing_data/*.bin`. The numerical
single and batch generators then accept that combined filename, the
`SinhCylindricalv2n2` domain tuple, the final stored-time bound, and the nominal
`--dt-spacetime-data` used for synthetic temporal-stencil edge times. Actual
slice-to-slice time gaps are approximately that nominal spacing, but the first
stored slice time is read from the `.bin` and printed at runtime; it is not a
CLI or `.par` input. Stored slice-table times remain authoritative and may be
near, rather than exactly equal to, nominal evolution intervals.

Claim evidence:
- Claim: The numerical pipeline combines full logical-grid slice payloads into a `.bin` container; nominal `--dt-spacetime-data` describes approximate gaps and synthetic edge times, while stored slice-table times remain authoritative at runtime.
- Role: user-facing commands and interfaces
- Deciding authority: `nrpy/examples/photon_single_geodesic_integrator_numerical.py` and `photon_batch_geodesic_integrator_numerical.py` — CLI help and runtime parameters
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/time_window_manager_numerical.py` — slice-table loading and time-window construction
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

```bash
python -m nrpy.examples.photon_single_geodesic_integrator_numerical \
  --bin-name DATASET.bin \
  --coord-system-numerical SinhCylindricalv2n2 \
  --domain 30.0 0.075 0.05 1.0 4.0 \
  --t-numerical-end 30.0 \
  --dt-spacetime-data 0.4 --t-start 10.0 \
  --observer-position 0.0 -10.0 0.0 \
  --observer-look-forward 0.0 1.0 0.0 \
  --observer-up 0.0 0.0 1.0 \
  --escape-radius 25.0 \
  --eom geodesic --interpolation-method g4DD
```

The numerical batch generator has the same numerical-data arguments and adds
the batch visualization/blueprint workflow. Both generators currently require
the numerical `.bin` file to be available under `project/raytracing_data/`.

## Sources

- [massive_single_geodesic_integrator_analytical.py](../../nrpy/examples/massive_single_geodesic_integrator_analytical.py) - `project_name`, `main_single`, `single_integrator_analytical`, `gsl-config`; official GSL [Using the Library](https://www.gnu.org/software/gsl/doc/html/usage.html) - `Compiling and Linking`
- [photon_single_geodesic_integrator_analytical.py](../../nrpy/examples/photon_single_geodesic_integrator_analytical.py) - `project_name`, `main_single`, observer-tetrad initialization, `rkf45_stage_update`, `trajectory.txt`
- [photon_batch_geodesic_integrator_analytical.py](../../nrpy/examples/photon_batch_geodesic_integrator_analytical.py) - `--outdir`, `--cuda`, `parallelization_mode`, `vis_command`, `blueprint_command`; official NVIDIA [NVCC guide](https://docs.nvidia.com/cuda/cuda-programming-guide/02-basics/nvcc.html) - `NVCC: The NVIDIA CUDA Compiler`
- [photon_single_geodesic_integrator_numerical.py](../../nrpy/examples/photon_single_geodesic_integrator_numerical.py) - numerical dataset CLI, `register_CFunction_numerical_interpolation`, `--t-start`
- [photon_batch_geodesic_integrator_numerical.py](../../nrpy/examples/photon_batch_geodesic_integrator_numerical.py) - numerical dataset CLI, `combine_raytracing_time_slices.py`, batch visualization workflow
- [combine_raytracing_time_slices.py](../../nrpy/infrastructures/BHaH/diagnostics/combine_raytracing_time_slices.py) - `InputSliceInfo`, `parse_args`, `--run-metadata`
- [visualize_trajectory.py](../../nrpy/examples/geodesic_visualizations/visualize_trajectory.py) - `visualize_trajectory`, `plot_trajectory`
- [blueprint_config_and_schema.py](../../nrpy/examples/geodesic_visualizations/blueprint_config_and_schema.py) - `BLUEPRINT_DTYPE`, `STOP_CONDITION_TERMINAL_PLANE`
- [render_lensed_image.py](../../nrpy/examples/geodesic_visualizations/render_lensed_image.py) - `generate_static_lensed_image`, `_process_blueprint_tile`, `_load_texture`
- [visualize_lensed_image.py](../../nrpy/examples/geodesic_visualizations/visualize_lensed_image.py) - `main`, `light_blueprint_{i:02d}_{j:02d}.bin`, `urlretrieve`
- [blueprint_analysis.py](../../nrpy/examples/geodesic_visualizations/blueprint_analysis.py) - `diagnose_blueprint`, `plot_heatmaps`
- [two_blackholes_collide.py](../../nrpy/examples/two_blackholes_collide.py) - `--raytracing-time`, `--raytracing-data-mode`, `--raytracing-static-christoffels`

## See Also

- [Examples](index.md)
- [Black Hole Evolution](black-hole-evolution.md)
- [Geodesics](../equations/general-relativity/geodesics.md)
- [Geodesics And Raytracing Runtime](../infrastructures/bhah/geodesics-and-raytracing-runtime.md)
