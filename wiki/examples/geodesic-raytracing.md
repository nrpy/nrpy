# Geodesic Raytracing

> Generate and run standalone geodesic integrators and render the batch photon results. · Status: confirmed
> Up: [Examples](index.md)

## Summary

The standalone examples evolve one massive particle, one photon, or a tiled camera grid of photons. The numerical photon generators read combined spacetime slices on a `SinhCylindricalv2n2` grid and generate CPU/OpenMP projects. See [Geodesics And Raytracing Runtime](../infrastructures/bhah/geodesics-and-raytracing-runtime.md) for interpolation, RKF45 control, and endpoint selection.

## Detail

### Analytic Geodesic Examples

Run the generators from the repository root with NRPy importable. The massive
example generates a Kerr-Schild C integrator with a GSL build dependency.

```bash
python -m nrpy.examples.massive_single_geodesic_integrator_analytical
cd project/massive_single_geodesic_integrator_analytical
make
./massive_single_geodesic_integrator_analytical
python3 visualize_trajectory.py --particle_type Massive
```

The analytical photon generators accept observer position, look direction,
up direction, and escape radius. The single-photon example defaults to
`KerrSchild_Cartesian` and also accepts
`BrillLindquist_InitialData_Static_Cartesian` through `--spacetime`.
The analytic batch generator also requires `--observer-fov W H` and
`--scan-density N`. It defaults to OpenMP. Use `--cuda` to generate a
CUDA project and `--outdir DIR` to select the parent of its build directory.

Evidence (descriptive behavior): `massive_single_geodesic_integrator_analytical.py`, `photon_single_geodesic_integrator_analytical.py`, and `photon_batch_geodesic_integrator_analytical.py`.

### Numerical Photon Examples

Generate a combined spacetime file as described in
[Standalone GR/BHaH](standalone-gr-bhah.md#numerical-spacetime-export).
The following example uses a `g4DD` file from a run with the shown grid
and times. Replace `DATASET.bin` with the filename under
`project/raytracing_data/` and match the geometry parameters to that file.

```bash
python -m nrpy.examples.photon_single_geodesic_integrator_numerical \
  --bin-name DATASET.bin \
  --coord-system-numerical SinhCylindricalv2n2 \
  --domain 30.0 0.075 0.05 1.0 4.0 \
  --t-numerical-end 30.0 --dt-spacetime-data 0.4 --t-start 10.0 \
  --observer-position 0.0 -10.0 0.0 --observer-look-forward 0.0 1.0 0.0 \
  --observer-up 0.0 0.0 1.0 --observer-fov 1.0 1.0 --escape-radius 25.0 \
  --eom normalized --interpolation-method g4DD
cd project/photon_single_geodesic_integrator_numerical_g4DD
make
./photon_single_geodesic_integrator_numerical_g4DD
```

For a batch, substitute `photon_batch_geodesic_integrator_numerical` in the
module and generated project names and add `--scan-density N`. The batch
option `--spacetime-project blackhole_spectroscopy` selects the source project
in the printed slice-combining command. The executable reads the dataset
selected by `--bin-name`. Set `--eom geodesic` for direct affine-parameter
evolution or `--eom normalized` for coordinate-time evolution. Match
`--interpolation-method` to the stored file type, `g4DD`, `g4DD_d0`, or
`GammaUDD`. The first stored time is read from the file. The nominal
`--dt-spacetime-data` may differ from individual gaps in the slice table.

Evidence (descriptive behavior): `photon_single_geodesic_integrator_numerical.py` and `photon_batch_geodesic_integrator_numerical.py`.

### Numerical Interpolation And RKF45 Options

Both numerical generators accept `--spatial-mask-threshold FACTOR`.
It sets `numerical_spacetime_spatial_mask_threshold` in the generated
`.par`, defaults to 10, and requires a finite value greater than one.
A sufficiently large finite value suppresses exclusions. The
[spatial mask calculation](../infrastructures/bhah/geodesics-and-raytracing-runtime.md#spatial-node-exclusion-for-numerical-photons)
describes how each attempt selects and reuses the mask.

`--interpolation-half-widths SPATIAL TEMPORAL` sets the two stencil half
widths. `--time-slice-stride N` selects the retained slices. Use
`--fast-light` at generation time to sample the geometry at `t_start`
throughout evolution. To compare this calculation with the dynamical one,
generate separate projects with `--outdir` and keep the dataset, observer
and tolerances the same.

`--l-sense LENGTH` and `--tol-scale-min SCALE` enable metric-gradient
scaling for `g4DD` and `g4DD_d0`. The sensing length must be finite and
nonnegative, and the minimum scale must lie in `(0, 1]`. Both default to one.
`--rkf45-step-range H_MIN H_MAX` configures the step range.
`--rkf45-max-consecutive-min-steps N` limits consecutive attempts at the
minimum step, with zero disabling that limit. `--evolution-measure-max`
sets the upper limit on the normal-observer log-energy measure.

Evidence (descriptive behavior): `photon_single_geodesic_integrator_numerical.py` and `photon_batch_geodesic_integrator_numerical.py`.

### Single-Photon Trajectory And Plane Files

Single-photon executables write `trajectory.txt`. Numerical trajectories
include the signed null residual, and `initial_state.txt` records the
initial direct four-momentum before any normalized conversion. Both
single-photon generators enable `rkf45_trials.txt` and `rkf45_stages.txt`
by default. Use `--disable-rkf45-trial-debug` to omit them. Numerical stage
rows describe successful geometry evaluations. The numerical debug files
include the selected temporal center and stencil times, the step before
and after the temporal cap, and the mask flag, excluded-node count, plane,
and spatial centers. With metric-gradient diagnostics, the trial file also
records `trial_gradient_max`, `q`, `tol_scale`, and the effective absolute
and relative tolerances.

For a terminal plane, supply `--terminal-plane-center`,
`--terminal-plane-normal`, `--terminal-plane-up`, and
`--terminal-plane-radius MIN MAX`. For a nonterminal plane, supply
`--non-terminal-plane-center`, `--non-terminal-plane-normal`, and
`--non-terminal-plane-up`. The plane distance and integration-parameter
refinement tolerances can be set with
`--plane-crossing-residual-tolerance` and
`--plane-crossing-parameter-tolerance`. The single-photon
`plane_crossings.txt` file reports the refined event times, local plane
coordinates, state, and refinement diagnostics. The batch executable uses
per-tile crossing files.

Evidence (descriptive behavior): `single_integrator_analytical.py::single_integrator_analytical`, `single_integrator_numerical.py::single_integrator_numerical`, and `batch_integrator_numerical.py::batch_integrator_numerical`.

### Camera Sampling And Batch Images

The photon generators accept `--a-range A_MIN A_MAX` and
`--b-range B_MIN B_MAX` in the original dimensionless camera coordinates
`[0,1]`. Use one tile for a cropped batch. With
`Delta a = a_max - a_min` and `Delta b = b_max - b_min`, a cropped batch
samples `scan_density` columns and
`floor(scan_density * Delta b / Delta a * tan(alpha_h/2) / tan(alpha_w/2))`
rows, where `alpha_w` and `alpha_h` are the full camera fields of view.
Samples are at cell centers. A single-photon run selects the center of the
requested rectangle. To match a batch ray, choose ranges centered at that
ray's sample coordinates and use the same camera and spacetime.

Evidence (descriptive behavior): `main_batch.py::main` and `set_initial_conditions_kernel.py::set_initial_conditions_kernel`.

Batch rays are written to `light_blueprint_XX_YY.bin`. The native binary
header is 92 bytes and gives the tile indices, camera fields of view, and
exact scan bounds. The 124-byte records include plane coordinates or sphere
angles, termination status, and sample coordinates. Version 10 indicates a
tile with apparent-horizon stops. For those rays, `y_t` stores the containing
horizon ID and `z_t` is `NaN`. Other tiles use version 9.

In the generated batch directory, run
`python3 visualize_lensed_image.py --pixel-width 450` to render the stored
scan rectangle. Pixel width controls rendering resolution independently of
the number of traced photons. Use `--sphere-image` for a custom sky
texture or `--terminal-plane-image` for a terminal-plane texture.
`python3 blueprint_analysis.py --no-plots` reports termination counts and
writes `blueprint_domains.txt`, grouping adjacent samples by termination
status. The interactive view can display these groups and the plane and
sphere diagnostics.

Evidence (descriptive behavior): `blueprint_config_and_schema.py::BLUEPRINT_DTYPE` and `BLUEPRINT_HORIZON_VERSION`; corroboration: `main_batch.py::main`, `visualize_lensed_image.py::main`, and `blueprint_analysis.py::label_domains`.

### Puncture Association And Apparent-Horizon Containment

The generated export script invokes a NumPy-based combiner for the apparent-horizon
text outputs. It uses the evolution project directory and its matching spacetime
files. To combine the horizons of an existing spectroscopy run, use the following
command.

```bash
python -m nrpy.infrastructures.BHaH.general_relativity.geodesics.combine_apparent_horizon_time_slices \
  --project-dir project/blackhole_spectroscopy \
  --example-name blackhole_spectroscopy
```

Pass the resulting `*_apparent_horizons.bin` path to either numerical
photon generator with `--apparent-horizon-bin-path FILE`. This enables
containment checks during integration. The alternative
`--puncture-diagnostics-dir DIR` uses the saved BHaHAHA text files to
classify the single photon's final accepted state or batch rays that stop
at the evolution-measure cutoff. The single-photon executable
writes `horizon_diagnostics.txt`. The batch executable writes
`light_blueprint_puncture_XX_YY.bin` with horizon classification and
signed radial gap. With `visualize_lensed_image.py --color-by-puncture`, the
renderer takes horizon-stop IDs directly from the main blueprint and takes
energy-stop labels from the matching puncture file when present.

Evidence (descriptive behavior): `raytracing_example_helper.py::shell_generator`, `photon_single_geodesic_integrator_numerical.py`, and `photon_batch_geodesic_integrator_numerical.py`.

Evidence (descriptive behavior): `single_integrator_numerical.py::single_integrator_numerical` and `batch_integrator_numerical.py::batch_integrator_numerical`.

Evidence (descriptive behavior): `render_lensed_image.py::generate_static_lensed_image` and `_process_blueprint_tile`.

## Sources

- [massive_single_geodesic_integrator_analytical.py](../../nrpy/examples/massive_single_geodesic_integrator_analytical.py) - generator arguments and project assembly
- [photon_single_geodesic_integrator_analytical.py](../../nrpy/examples/photon_single_geodesic_integrator_analytical.py) - generator arguments and project assembly
- [photon_batch_geodesic_integrator_analytical.py](../../nrpy/examples/photon_batch_geodesic_integrator_analytical.py) - generator arguments and project assembly
- [photon_single_geodesic_integrator_numerical.py](../../nrpy/examples/photon_single_geodesic_integrator_numerical.py) - generator arguments and project assembly
- [photon_batch_geodesic_integrator_numerical.py](../../nrpy/examples/photon_batch_geodesic_integrator_numerical.py) - generator arguments and project assembly
- [single_integrator_analytical.py](../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/single_integrator_analytical.py) - single_integrator_analytical
- [single_integrator_numerical.py](../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/single_integrator_numerical.py) - single_integrator_numerical
- [batch_integrator_numerical.py](../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_numerical.py) - batch_integrator_numerical
- [main_batch.py](../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/main_batch.py) - main_batch
- [set_initial_conditions_kernel.py](../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/set_initial_conditions_kernel.py) - set_initial_conditions_kernel
- [visualize_trajectory.py](../../nrpy/examples/geodesic_visualizations/visualize_trajectory.py) - trajectory plotting
- [blueprint_config_and_schema.py](../../nrpy/examples/geodesic_visualizations/blueprint_config_and_schema.py) - blueprint headers and fields
- [visualize_lensed_image.py](../../nrpy/examples/geodesic_visualizations/visualize_lensed_image.py) - image options
- [render_lensed_image.py](../../nrpy/examples/geodesic_visualizations/render_lensed_image.py) - horizon and puncture colors
- [blueprint_analysis.py](../../nrpy/examples/geodesic_visualizations/blueprint_analysis.py) - termination diagnostics
- [raytracing_example_helper.py](../../nrpy/infrastructures/BHaH/general_relativity/geodesics/raytracing_example_helper.py) - `shell_generator`

## See Also

- Depends on: [Geodesics And Raytracing Runtime](../infrastructures/bhah/geodesics-and-raytracing-runtime.md)
- See also: [Standalone GR/BHaH](standalone-gr-bhah.md)
- See also: [Build And Run](../architecture/build-and-run.md)
