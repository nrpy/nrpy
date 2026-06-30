# SEOBNR BOB Generated Library

> BHaH generated C library flow for SEOBNRv5, BOB, and SEBOBv2 waveforms. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [BHaH](index.md)

## Summary

The BHaH SEOBNR/BOB infrastructure turns symbolic SEOBNR and BOB equation
modules into a standalone C executable or library-like generated project. The
assembly is not a grid evolution: it builds a `commondata_struct`, registers
GSL-backed root finding, ODE integration, spline interpolation, dynamics,
inspiral waveform, merger waveform, NQC, rotation, and output helpers, emits
headers and a Makefile with GSL flags, and fills `commondata` arrays for
dynamics and waveform samples.

## Detail

The aligned-spin SEOBNRv5 example registers a compact waveform pipeline. Its
generated `main` initializes `commondata`, parses parfile and command-line
inputs, computes coefficients, solves conservative and dissipative initial
conditions, integrates dynamics, generates inspiral waveform samples, applies
NQC corrections, builds the IMR waveform, optionally prints waveform rows, and
frees the `commondata` arrays. The SEBOBv2 example extends that route with
quasi-precessing spin coefficients and spin dynamics, post-adiabatic
integration, higher-mode inspiral storage, optional coprecessing-rotation
sandbox output, SEBOBv2 NQC corrections, BOBv2 merger waveform generation, and
explicit cleanup of many GSL spline accelerators stored in `commondata`.

Initial-condition and root-finding helpers are ordinary generated C functions
around GSL. `register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_conservative`
builds a two-variable `gsl_multiroot_function_fdf`, seeds `r` and `pphi` from
`initial_omega`, calls the multidimensional root wrapper, and writes the result
to `commondata`. The nodf and dissipative initial-condition routes use the same
support family, while `root_finding_1d` and `root_finding_multidimensional`
provide reusable GSL wrappers used later by post-adiabatic and NQC stages.

Dynamics generation is split into equation kernels and integration managers.
`register_CFunction_SEOBNRv5_aligned_spin_pa_integration` builds a radial grid,
solves alternating post-adiabatic `pphi` and `prstar` updates with one-dimensional
root finding, differentiates radial arrays, integrates time and phase, calls
`SEOBNRv5_aligned_spin_ode_integration`, merges post-adiabatic and ODE dynamics,
and stores augmented quantities in `commondata->dynamics_low`,
`dynamics_fine`, and `dynamics_raw`. `register_CFunction_SEOBNRv5_aligned_spin_ode_integration`
uses `gsl_odeiv2_step_rk8pd`, adaptive control, and
`SEOBNRv5_aligned_spin_right_hand_sides`; it tracks merger-side stop conditions,
splits coarse and fine dynamics, and refines peak times through spline helpers.

The inspiral waveform layer consumes stored dynamics rather than re-integrating.
`register_CFunction_SEOBNRv5_aligned_spin_waveform_from_dynamics` allocates
`waveform_low` and `waveform_fine`, loads each dynamics row into a local
`NUMVARS` array, calls `SEOBNRv5_aligned_spin_waveform`, and stores complex
strain modes by `IDX_WF`. SEBOBv2 uses the higher-mode variants in the same
subtree and defines mode indices such as `STRAIN22`, `STRAIN21`, `STRAIN33`,
`STRAIN32`, `STRAIN44`, `STRAIN43`, and `STRAIN55` in generated
`BHaH_defines.h`.

Merger and IMR assembly sit above the inspiral arrays. The aligned-spin
`SEOBNRv5_aligned_spin_IMR_waveform` interpolates inspiral modes at the
attachment time, allocates ringdown arrays, then selects either native
SEOBNRv5 merger-ringdown or BOB ringdown based on the registration flag before
concatenating inspiral and ringdown into `commondata->waveform_IMR`.
`SEBOBv2_IMR_waveform` follows the same assembly shape but calls
`BOB_v2_waveform_from_times` and stores the current `(2,2)` SEBOBv2 output.
The BOBv2 helper evaluates amplitude and wrapped phase at requested times,
unwraps phase with `SEOBNRv5_aligned_spin_unwrap`, and zero-references the
phase at the first ringdown sample.

NQC correction functions solve attachment systems and mutate the inspiral
waveform arrays before IMR assembly. `SEBOBv2_NQC_corrections` allocates basis,
radius, frequency, amplitude, and phase arrays; locates the ISCO/attachment
time; configures BOB peak data; builds cropped spline systems around the
attachment point; solves amplitude and phase correction systems with GSL
matrices; and applies the corrections to low and fine waveform samples.
The aligned-spin and higher-mode BOB NQC modules follow the same infrastructure
role for other approximant choices.

Precessing support is present as generated library helpers, not a complete
native precessing IMR replacement in these examples. `SEOBNRv5_coprecessing_angles`
builds physical Euler-angle arrays from spin dynamics, and
`SEOBNRv5_coprecessing_rotations` rotates coprecessing inspiral modes into
observer polarizations using the registered mode arrays. Shared rotation and
spin-weighted spherical-harmonic support live in BHaH helper modules and are
cited when those generated functions are assembled.

Project assembly is explicit. `sebobv2.py` writes commondata-only
`set_CodeParameters*.h` parameter-access headers, registers
`commondata_struct_set_to_default`, emits the parfile/parser with waveform
command-line inputs, copies `spline_struct.h`, adds GSL and complex-number
includes to `BHaH_defines.h`, defines SEOBNR array index macros and mode
constants, registers the generated `main`, and constructs a Makefile with
`gsl-config --cflags` and `gsl-config --libs`. CI then builds current and
trusted SEOB/SEBOB executables and compares waveform outputs through the SEOB
consistency scripts.

## Sources

- [seobnrv5_aligned_spin_inspiral.py](../../../nrpy/examples/seobnrv5_aligned_spin_inspiral.py) - `register_CFunction_main_c`, `SEOBNRv5_aligned_spin_IMR_waveform`
- [sebobv2.py](../../../nrpy/examples/sebobv2.py) - `register_CFunction_main_c`, `SEBOBv2_IMR_waveform`, `supplemental_defines_dict`
- [SEOBNRv5_aligned_spin_initial_conditions_conservative.py](../../../nrpy/infrastructures/BHaH/seobnr/initial_conditions/SEOBNRv5_aligned_spin_initial_conditions_conservative.py) - `register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_conservative`
- [root_finding_1d.py](../../../nrpy/infrastructures/BHaH/seobnr/utils/root_finding_1d.py) - `register_CFunction_root_finding_1d`
- [root_finding_multidimensional.py](../../../nrpy/infrastructures/BHaH/seobnr/utils/root_finding_multidimensional.py) - `register_CFunction_root_finding_multidimensional`
- [SEOBNRv5_aligned_spin_pa_integration.py](../../../nrpy/infrastructures/BHaH/seobnr/dynamics/SEOBNRv5_aligned_spin_pa_integration.py) - `register_CFunction_SEOBNRv5_aligned_spin_pa_integration`
- [SEOBNRv5_aligned_spin_ode_integration.py](../../../nrpy/infrastructures/BHaH/seobnr/dynamics/SEOBNRv5_aligned_spin_ode_integration.py) - `register_CFunction_SEOBNRv5_aligned_spin_ode_integration`
- [SEOBNRv5_aligned_spin_waveform_from_dynamics.py](../../../nrpy/infrastructures/BHaH/seobnr/inspiral_waveform/SEOBNRv5_aligned_spin_waveform_from_dynamics.py) - `register_CFunction_SEOBNRv5_aligned_spin_waveform_from_dynamics`
- [SEOBNRv5_aligned_spin_IMR_waveform.py](../../../nrpy/infrastructures/BHaH/seobnr/SEOBNRv5_aligned_spin_IMR_waveform.py) - `register_CFunction_SEOBNRv5_aligned_spin_IMR_waveform`
- [SEBOBv2_IMR_waveform.py](../../../nrpy/infrastructures/BHaH/seobnr/SEBOBv2_IMR_waveform.py) - `register_CFunction_SEBOBv2_IMR_waveform`
- [BOB_v2_waveform_from_times.py](../../../nrpy/infrastructures/BHaH/seobnr/merger_waveform/BOB_v2_waveform_from_times.py) - `register_CFunction_BOB_v2_waveform_from_times`
- [SEBOBv2_NQC_corrections.py](../../../nrpy/infrastructures/BHaH/seobnr/nqc_corrections/SEBOBv2_NQC_corrections.py) - `register_CFunction_SEBOBv2_NQC_corrections`
- [SEOBNRv5_coprecessing_angles.py](../../../nrpy/infrastructures/BHaH/seobnr/SEOBNRv5_coprecessing_angles.py) - `register_CFunction_SEOBNRv5_coprecessing_angles`
- [SEOBNRv5_coprecessing_rotations.py](../../../nrpy/infrastructures/BHaH/seobnr/SEOBNRv5_coprecessing_rotations.py) - `register_CFunction_SEOBNRv5_coprecessing_rotations`
- [spin_weight_minus2_spherical_harmonics.py](../../../nrpy/infrastructures/BHaH/special_functions/spin_weight_minus2_spherical_harmonics.py) - `register_CFunction_spin_weight_minus2_sph_harmonics`, `spin_weight_minus2_sph_harmonics`
- [spline_struct.h](../../../nrpy/infrastructures/BHaH/seobnr/spline_struct.h) - `spline_data`

## See Also

- [BHaH](index.md)
- [Lifecycle And Project Assembly](lifecycle-and-project-assembly.md)
- [SEOBNR And BOB](../../equations/seobnr/index.md)
- [Generated Project CI](../../validation/generated-project-ci.md)
