# Standalone GR/BHaH

> Compare standalone BHaH numerical-relativity generators by initial data, coordinates, formulations, diagnostics, and OpenMP/CUDA build modes. · Status: confirmed
> Up: [Examples](index.md)

## Summary

Four generators in `nrpy/examples/` emit standalone single-patch BHaH general-relativity projects under `project/<project_name>/`: `two_blackholes_collide.py` (lightweight Brill-Lindquist head-on binary), `spinning_blackhole.py` (UIUC spinning puncture), `blackhole_spectroscopy.py` (quasicircular TwoPunctures binary with optional `--fccz4`, checkpointing, and `Psi4`), and `kasner_exact_evolution.py` (cosmological `GeneralRFM_fisheyeN1` benchmark). Each recreates its `project/<project_name>/` directory before emitting C/CUDA sources, `<project_name>.par`, and a `Makefile`.

## Detail

### `two_blackholes_collide.py` And `spinning_blackhole.py`

- **`two_blackholes_collide.py`** (`README.md`'s lightweight single-patch BH@H example): Evolves equal-mass opposite-`z` `IDtype = "BrillLindquist"` data (Cartesian ID coordinates) on a `Spherical` grid with `OnePlusLog` lapse, `GammaDriving2ndOrder_Covariant` shift, RK4 MoL, `fd_order = 4`, `rfm_precompute`, Sommerfeld radiation boundaries, separate Ricci evaluation, algebraic determinant/trace projection, constraint/nearest/volume diagnostics, and `convergence_factor` CLI scaling. Accepts `--cuda`, `--floating_point_precision`, and `--raytracing-time` for numerical raytracing output.
- **`spinning_blackhole.py`**: Evolves `IDtype = "UIUCBlackHole"` (`M = 1.0`, `chiU = (0, 0, +0.8)`, spherical ID coordinates; `OffsetKerrSchild` is present as an inactive branch) on `SinhCylindrical` with the same `OnePlusLog`/`GammaDriving2ndOrder_Covariant`, RK4, `fd_order = 4`, separate Ricci, and radiation-boundary configuration. Accepts `--cuda` and `--floating_point_precision`.

### `blackhole_spectroscopy.py` (`--fccz4`, TwoPunctures, And `Psi4`)

`blackhole_spectroscopy.py` evolves `IDtype = "TP_Interp"` on `SinhCylindrical` using `NRPyPN_quasicircular_momenta` and `TwoPunctures_lib` (freeing TwoPunctures derivative buffers after interpolation), overriding the `psi^n` initial lapse with `W` for slow-start lapse (SSL; see [BHaH GR Initial Data And Psi4](../infrastructures/bhah/gr-initial-data-and-psi4.md)). It configures `fd_order = 8`, separate Ricci, radiation boundaries, checkpointing every `2.0`, `Psi4` `s = -2` spherical-harmonic extraction, and `gsl-config` `Makefile` flags.

BSSN is the default. The `--fccz4` flag switches to fCCZ4 RHS/gauge kernels and zeroes `Theta_fCCZ4` only on fresh ADM-to-BSSN/fCCZ4 conversion (checkpoint restore returns earlier without zeroing `Theta_fCCZ4`). Both formulations run `enforce_detgbar_equals_detghat_trAzero` after initial-data boundaries and in the MoL post-RHS step.

Evidence (scientific/interface requirement): `blackhole_spectroscopy.py` (`--fccz4`, MoL/ID calls); corroboration: `BHaH/general_relativity/initial_data.py`; `ADM_Initial_Data_Reader__BSSN_Converter.py`; `BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py`.

### Numerical Spacetime Export

`two_blackholes_collide.py` and `blackhole_spectroscopy.py` accept
`--raytracing-time [T_FINAL OUTPUT_EVERY]` for OpenMP, double-precision
evolution. `--raytracing-data-mode` selects the data to write. The
`--raytracing-coord-system SinhCylindricalv2n2` option selects the grid
coordinate system used by the numerical photon integrators. Set its
parameters with `--raytracing-domain` and its resolution with `--raytracing-Nxx`.
Spectroscopy accepts `--initial-sep` and `--initial-p-r` for the
TwoPunctures separation and radial momentum.

From the generated evolution directory, run
`./run_raytracing_data_pipeline.sh` to build, evolve, and combine the
numerical spacetime and apparent-horizon data under `project/raytracing_data/`.
The apparent-horizon combination step requires NumPy. The generated
`raytracing_run_metadata.json` records generation settings. Run settings
changed by editing the `.par` file are not copied back to that JSON.

Evidence (descriptive behavior): `raytracing_example_helper.py::add_raytracing_cli`, `check_raytracing_cli`, and `shell_generator`.

### `kasner_exact_evolution.py` And OpenMP/CUDA/BHaHAHA Rules

`kasner_exact_evolution.py` evolves `IDtype = "Kasner"` on `GeneralRFM_fisheyeN1` with frozen lapse/shift, RK4, `fd_order = 4`, extrapolation outer boundaries, Kasner exponent constraint checks, and nearest diagnostics. Although its CLI parser defines `--cuda`, CUDA generation fails because `initial_data.py` requests `generalrfm_precompute` support for `GeneralRFM_fisheyeN1`, which raises on `parallelization == "cuda"` before `Ricci_eval` runs; only OpenMP (`--floating_point_precision`) is usable.

Evidence (descriptive behavior): `kasner_exact_evolution.py` (`CoordSystem = "GeneralRFM_fisheyeN1"`, `--cuda`); corroboration: `BHaH/general_relativity/initial_data.py` and `BHaH/generalrfm_precompute.py::register_CFunction_generalrfm_precompute`.

On the three black-hole generators, `--cuda` emits `.cu` files, `cuda_intrinsics.h`, `nvcc` build rules, and non-`restrict` device pointers, whereas OpenMP + `double` builds link a BHaHAHA apparent-horizon subdirectory (non-`double` OpenMP BHaHAHA raises an error; Kasner omits BHaHAHA). GitHub `codegen-{ubuntu,mac}` compiles default OpenMP builds (and `.github/full_nrpy_local_ci.sh` adds CUDA builds) for the three black-hole generators without executing them.

## Sources

- [README.md](../../README.md) - `## Lightweight Single-Patch BH@H Example`, `### Standalone BHaH Generators`
- [two_blackholes_collide.py](../../nrpy/examples/two_blackholes_collide.py) - `project_name`, `CoordSystem`, `IDtype`, `--raytracing-time`, `enable_bhahaha`
- [blackhole_spectroscopy.py](../../nrpy/examples/blackhole_spectroscopy.py) - `parser`, `enable_fCCZ4`, `IDtype`, `TwoPunctures_lib`, `enable_psi4_diagnostics`, `read_checkpoint`, `write_checkpoint`
- [spinning_blackhole.py](../../nrpy/examples/spinning_blackhole.py) - `project_name`, `CoordSystem`, `IDtype`, `default_BH_spin_chiU`
- [kasner_exact_evolution.py](../../nrpy/examples/kasner_exact_evolution.py) - `parser`, `CoordSystem`, `IDtype`, `LapseEvolutionOption`, `ShiftEvolutionOption`
- [initial_data.py](../../nrpy/infrastructures/BHaH/general_relativity/initial_data.py) - `register_CFunction_initial_data`
- [ADM_Initial_Data_Reader__BSSN_Converter.py](../../nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py) - `register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN`
- [enforce_detgbar_equals_detghat_trAzero.py](../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- [generalrfm_precompute.py](../../nrpy/infrastructures/BHaH/generalrfm_precompute.py) - `register_CFunctions_generalrfm_support`, `register_CFunction_generalrfm_precompute`

- [raytracing_example_helper.py](../../nrpy/infrastructures/BHaH/general_relativity/geodesics/raytracing_example_helper.py) - generation options, run metadata, and build and export commands

## See Also

- Depends on: [BHaH Project Assembly](../infrastructures/bhah/lifecycle-and-project-assembly.md)
- Depends on: [BHaH GR Evolution And Constraints](../infrastructures/bhah/gr-evolution-and-constraints.md)
- Depends on: [BHaH GR Initial Data And Psi4](../infrastructures/bhah/gr-initial-data-and-psi4.md)
