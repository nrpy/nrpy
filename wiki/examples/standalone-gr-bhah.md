# Standalone GR/BHaH

> Route standalone BHaH numerical-relativity generators by initial data, coordinates, diagnostics, and build mode. · Status: confirmed · Last reconciled: 07-19-2026
> Up: [Examples](index.md)

## Summary

Standalone GR/BHaH examples are Python generators that write complete BHaH
projects under `project/<project_name>/`, then expect `make` and the generated
executable to run inside that project directory. The four checked-in owners are
`two_blackholes_collide.py` for compact Brill-Lindquist binary black-hole
evolution, `blackhole_spectroscopy.py` for a TwoPunctures-backed binary with
checkpointing and Psi4 diagnostics, `spinning_blackhole.py` for a single
spinning UIUC black hole, and `kasner_exact_evolution.py` for a vacuum Kasner
benchmark.

Each generator deletes and recreates its fixed project directory. Configured
GitHub codegen generates and builds the default OpenMP path for the first three,
but omits Kasner and runs no executable. The local helper configures OpenMP and
CUDA builds for the first three, omits Kasner, and also runs no executable or
scientific check. Preserve output before reruns; runtime remains manual here.

## Detail

`python -m nrpy.examples.two_blackholes_collide` generates
`project/two_blackholes_collide/`. It uses `IDtype = "BrillLindquist"` with
Cartesian initial-data coordinates on a `Spherical` evolution coordinate
system, `OnePlusLog` lapse, `GammaDriving2ndOrder_Covariant` shift, RK4 Method
of Lines timestepping, fourth-order finite differences, reference-metric
precomputation, outgoing radiation boundaries, and a separate Ricci evaluation
before BSSN RHS evaluation. Its defaults place equal masses at opposite
z-axis positions. The generator accepts `--cuda`,
`--floating_point_precision`, and `--raytracing-outputs`; raytracing output is
explicitly OpenMP-only and is rejected with `--cuda`.

`python -m nrpy.examples.blackhole_spectroscopy` generates
`project/blackhole_spectroscopy/`. It uses `IDtype = "TP_Interp"` with
Cartesian initial-data coordinates on `SinhCylindrical`, registers
`NRPyPN_quasicircular_momenta`, registers the TwoPunctures library, solves the
TwoPunctures initial-data persist structure before filling BHaH initial data,
and frees the TwoPunctures derivative storage afterward. It uses
`OnePlusLog`/`GammaDriving2ndOrder_Covariant`, eighth-order finite differences,
separate Ricci, outgoing radiation boundaries, checkpointing every `2.0` by
default, Psi4 and spin-weight minus-two spherical-harmonic diagnostics, and GSL
Makefile flags through `gsl-config`. Its source-backed generation flags are
`--cuda`, `--fccz4`, and `--floating_point_precision`. BSSN remains the
default. `--fccz4` selects fCCZ4 RHS/gauge registration and zeroes
`Theta_fCCZ4` only while fresh ADM data are converted; checkpoint loading is
attempted first and returns without that zeroing. Both choices use the same
combined determinant/trace projection after initial-data boundary handling and
in the Method of Lines post-RHS hook.

Claim evidence:
- Claim: `blackhole_spectroscopy.py` defaults to BSSN and accepts `--fccz4`; that flag selects fCCZ4 RHS/gauge registration and fresh-data-only Theta zeroing after the checkpoint branch, while both formulations use the same initial-data and post-RHS conformal projection; this source ordering does not prove restart correctness or scientific validity.
- Role: public/scientific contract
- Deciding authority: [blackhole_spectroscopy.py](../../nrpy/examples/blackhole_spectroscopy.py), `parser`, `enable_fCCZ4`, initial-data/RHS/Method of Lines registrations
- Corroboration: [initial_data.py](../../nrpy/infrastructures/BHaH/general_relativity/initial_data.py), `register_CFunction_initial_data`; [ADM_Initial_Data_Reader__BSSN_Converter.py](../../nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py), `register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN`; [enforce_detgbar_equals_detghat_trAzero.py](../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), combined projector
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=Python 3.12.3, GCC 13.3.0, GNU Make 4.3; backend=BHaH OpenMP; precision=double; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=not-run; options=default BSSN generation/build, opt-in fCCZ4 generation/build and t_final=0.1 startup through iteration 1; date=08-28-2026`

`python -m nrpy.examples.spinning_blackhole` generates
`project/spinning_blackhole/`. It evolves `IDtype = "UIUCBlackHole"` in
`SinhCylindrical` coordinates from spherical initial-data coordinates, with an
inactive source option for `OffsetKerrSchild` in the script. It uses the same
`OnePlusLog` lapse, `GammaDriving2ndOrder_Covariant` shift, RK4, fourth-order
finite differences, separate Ricci, and outgoing radiation pattern as the
compact binary example, but owns the single spinning black-hole defaults:
`M = 1.0` and spin vector `(chi_x, chi_y, chi_z) = (0.0, 0.0, +0.8)`.
Its source-backed generation flags are `--cuda` and `--floating_point_precision`.

`python -m nrpy.examples.kasner_exact_evolution` generates
`project/kasner_exact_evolution/`. It is a benchmark rather than a black-hole
horizon workflow: `IDtype = "Kasner"`, Cartesian initial-data coordinates,
`GeneralRFM_fisheyeN1` evolution coordinates, frozen lapse and shift, RK4,
fourth-order finite differences, extrapolation outer boundaries, and a guard
that requires the Kasner exponents to satisfy both Kasner constraints. It uses
Kasner-specific diagnostic gridfunction registration and nearest diagnostics.
The source keeps separate Ricci for supported paths, but disables device-side
separate Ricci when CUDA and GeneralRFM are combined, then still registers a
host-only Ricci path for CUDA. Its source-backed generation flags are `--cuda`
and `--floating_point_precision`.

All four generators default to OpenMP and switch to CUDA only when `--cuda` is
present. In CUDA mode they register CUDA host/device helpers, use `nvcc`,
choose `.cu` source output, copy `cuda_intrinsics.h`, and relax generated
pointer qualifiers from `*restrict` to `*`. In OpenMP mode, the three
black-hole examples generate or link a BHaHAHA static library subdirectory and
require double precision for that integration; non-double OpenMP BHaHAHA
generation raises an error. The Kasner benchmark does not generate BHaHAHA.

The shared BHaH runtime pieces include nearest and volume diagnostics,
diagnostic gridfunction header generation, progress output, constraint
evaluation, determinant enforcement, curvilinear boundary setup,
Cartesian/native coordinate transforms, basis transforms, generated parameter
files, `convergence_factor` command-line parsing, copied intrinsics, generated
headers, and a Makefile. Details of generated project assembly belong in
[Lifecycle And Project Assembly](../infrastructures/bhah/lifecycle-and-project-assembly.md);
details of diagnostic file formats, checkpoint payloads, and raytracing stage
files belong in
[Diagnostics Output And Checkpointing](../infrastructures/bhah/diagnostics-output-and-checkpointing.md).

For related workflows, use
[Geodesics And Raytracing Runtime](../infrastructures/bhah/geodesics-and-raytracing-runtime.md)
for standalone geodesics and evolution-time raytracing export, and
[Apparent Horizon Library](apparent-horizon-library.md) for the BHaHAHA library
generator that horizon-enabled black-hole examples call.

## Sources

- [two_blackholes_collide.py](../../nrpy/examples/two_blackholes_collide.py) - `project_name`, `CoordSystem`, `IDtype`, `--raytracing-outputs`, `enable_bhahaha`
- [blackhole_spectroscopy.py](../../nrpy/examples/blackhole_spectroscopy.py) - `project_name`, `IDtype`, `BHaH.general_relativity.TwoPunctures.TwoPunctures_lib.register_C_functions`, `enable_psi4_diagnostics`, `BHaH.read_checkpoint.register_CFunction_read_checkpoint`, `BHaH.write_checkpoint.register_CFunction_write_checkpoint`
- [initial_data.py](../../nrpy/infrastructures/BHaH/general_relativity/initial_data.py) - checkpoint-first `register_CFunction_initial_data`
- [ADM_Initial_Data_Reader__BSSN_Converter.py](../../nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py) - fresh-data `Theta_fCCZ4` initialization
- [enforce_detgbar_equals_detghat_trAzero.py](../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - shared determinant/trace projector
- [spinning_blackhole.py](../../nrpy/examples/spinning_blackhole.py) - `project_name`, `CoordSystem`, `IDtype`, `spin_alignment_vector_params`, `default_BH_spin_chiU`
- [kasner_exact_evolution.py](../../nrpy/examples/kasner_exact_evolution.py) - `project_name`, `IDtype`, `LapseEvolutionOption`, `ShiftEvolutionOption`, `use_separate_ricci`

## See Also

- [Examples](index.md)
- [Black Hole Evolution](black-hole-evolution.md)
- [Apparent Horizon Library](apparent-horizon-library.md)
- [Lifecycle And Project Assembly](../infrastructures/bhah/lifecycle-and-project-assembly.md)
- [Diagnostics Output And Checkpointing](../infrastructures/bhah/diagnostics-output-and-checkpointing.md)
- [Geodesics And Raytracing Runtime](../infrastructures/bhah/geodesics-and-raytracing-runtime.md)
