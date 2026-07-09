# Waveform JAX PN Generators

> Map SEOBNR/SEBOB waveform, JAX, and PN momentum example generators to generated project names, dependency classes, and example-owned consistency checks. · Status: confirmed · Last reconciled: 06-30-2026
> Up: [Examples](index.md)

## Summary

This example leaf covers four related generator families. The
`seobnrv5_aligned_spin_inspiral` route creates GSL-backed BHaH C waveform
projects whose project name is selected from the SEOBNRv5 approximant flag and
optional calibration suffix. `sebobv2` creates one GSL-backed BHaH C waveform
project named `sebobv2`. `sebobv1_jax` creates a Python/JAX project named
`sebobv1_jax`, but the current source-backed surface only initializes
SEOBNRv5 aligned-spin coefficients. `nrpypn_quasicircular_momenta` creates a
BHaH C utility project named `nrpypn_quasicircular_momenta` for 3.5PN
quasicircular binary-black-hole momenta.

The SEOB/SEBOB consistency scripts are example-owned validation helpers. They
build trusted and current generated executable directories, run matching
waveform outputs, compute amplitude-plus-phase RMSE-style differences, and
require the current-vs-trusted median error to stay within a perturbation-based
roundoff baseline.

## Detail

`nrpy.examples.seobnrv5_aligned_spin_inspiral` is the flag-driven SEOBNRv5
aligned-spin generator. If no approximant flag is supplied, it defaults to
`-seobnrv5_bob`. Normal use chooses one approximant flag, then optionally adds
one calibration flag:

| Flag family | Generated project name | NQC choice | Merger-ringdown choice | Dependency class |
| --- | --- | --- | --- | --- |
| `-seobnrv5_bob` or no approximant flag | `seobnrv5_bob` | BOB-informed NQC | BOB-informed merger-ringdown | BHaH C project with GSL |
| `-seobnrv5_nrnqc_bob` | `seobnrv5_nrnqc_bob` | native numerical-relativity NQC | BOB-informed merger-ringdown | BHaH C project with GSL |
| `-seobnrv5_nrpy` | `seobnrv5_nrpy` | native numerical-relativity NQC | native SEOBNRv5 merger-ringdown | BHaH C project with GSL |

Project-name expansion is string-based. The script starts from the approximant
base name, then appends `_calibration_no_spin` for `-calibration_no_spin` or
`_calibration_spin` for `-calibration_spin`. It rejects using both calibration
flags at once. Thus `-seobnrv5_bob -calibration_no_spin` writes
`project/seobnrv5_bob_calibration_no_spin/`, while
`-seobnrv5_nrpy -calibration_spin` writes
`project/seobnrv5_nrpy_calibration_spin/`. The calibration flags also flow into
SEOBNRv5 aligned-spin coefficient registration, so they are not only naming
suffixes.

The SEOBNRv5 waveform projects are GSL-backed C projects. The generator sets
`Infrastructure` to `BHaH`, registers C functions for commondata I/O,
root-finding, initial conditions, ODE integration, spline/interpolation-backed
waveform processing, NQC corrections, merger-ringdown, and IMR assembly, then
emits GSL headers and a Makefile using `gsl-config --cflags` and
`gsl-config --libs`. Its command-line/parfile inputs are `mass_ratio`, `chi1`,
`chi2`, `initial_omega`, `total_mass`, and `dt`.

`nrpy.examples.sebobv2` is a separate BHaH C waveform generator with fixed
`project_name = "sebobv2"`. It is also GSL-backed and uses the same
waveform-input list: `mass_ratio`, `chi1`, `chi2`, `initial_omega`,
`total_mass`, and `dt`. The source describes it as work in progress: current
behavior computes aligned-spin `(2,2)` IMR modes using SEOBNRv5 and BOBv2,
while higher-mode and broader precessing improvements are listed as incomplete
or partial. The source registers higher-mode inspiral storage and optional
coprecessing-rotation sandbox code, but its generated stdout path still prints
the `(2,2)` strain.

`nrpy.examples.sebobv1_jax` is not a GSL-backed waveform executable. It sets
`Infrastructure` to `JAX`, writes `project/sebobv1_jax/`, registers Commondata
metadata, registers the `SEOBNRv5_aligned_spin_coefficients` Python/JAX
function, runs parallel codegen, and emits a generated Python package. Its own
user-facing text says the current project initializes SEOBNRv5 aligned-spin
coefficients. Class generation, the full SEBOBv1 JAX port, waveform-generation
tests, documentation, mismatch/calibration utilities, and possible SEBOBv2
JAX work are roadmap items in the example source, not completed workflow
features. Existing JAX infrastructure pages own the generated-package lifecycle
and the narrower Commondata/PyFunction details.

`nrpy.examples.nrpypn_quasicircular_momenta` is the PN momentum utility
generator. It writes `project/nrpypn_quasicircular_momenta/`, sets
`Infrastructure` to `BHaH`, registers
`NRPyPN_quasicircular_momenta`, emits a default parfile and parser, then calls
the generated C function from `main`. Its command-line/parfile inputs are
`initial_sep`, `mass_ratio`, `bbhxy_BH_M_chix`, `bbhxy_BH_M_chiy`,
`bbhxy_BH_M_chiz`, `bbhxy_BH_m_chix`, `bbhxy_BH_m_chiy`, and
`bbhxy_BH_m_chiz`. Unlike the SEOBNR/SEBOB waveform examples above, this
source does not add GSL headers, GSL flags, or GSL libraries.

`nrpy.examples.tests.sebob_consistency_check` and
`nrpy.examples.tests.sebobv2_consistency_check` compare generated waveform
executables at workflow level. Each script takes `--trusted-exec` and
`--current-exec` as directories; the executable and parfile are expected to be
named after the directory basename. The SEOBNRv5-family script also checks
that trusted and current approximant directory names match, and it forces zero
spins for `no_spin` calibration approximants. Both scripts rebuild trusted and
current projects with `make clean` and `make`, run ten deterministic input
sets, form complex `h22` from stdout columns, unwrap phase, interpolate both
waveforms over the shared time interval, and compute normalized RMSE-style
amplitude and phase errors. A run passes only when the median current-vs-trusted
error is no larger than the median trusted-vs-perturbed baseline error.

These consistency scripts document the example workflow, not a new validation
subsystem. The broader CI page is context for where generated projects are
built and compared; this page stays scoped to what the example generators and
their local helper scripts do.

## Sources

- [seobnrv5_aligned_spin_inspiral.py](../../nrpy/examples/seobnrv5_aligned_spin_inspiral.py) - `argparse` flags, `project_name`, `register_CFunction_main_c`, `register_CFunction_cmdline_input_and_parfile_parser`, `output_CFunctions_function_prototypes_and_construct_Makefile`
- [sebobv2.py](../../nrpy/examples/sebobv2.py) - `project_name`, workflow docstring, `register_CFunction_main_c`, `register_CFunction_cmdline_input_and_parfile_parser`, `output_CFunctions_function_prototypes_and_construct_Makefile`
- [sebobv1_jax.py](../../nrpy/examples/sebobv1_jax.py) - `project_name`, `Infrastructure`, `register_commondata_params`, `register_PyFunction_SEOBNRv5_aligned_spin_coefficients`, `output_PyFunction_files_and_construct_project`
- [nrpypn_quasicircular_momenta.py](../../nrpy/examples/nrpypn_quasicircular_momenta.py) - `project_name`, `register_CFunction_main_c`, `register_CFunction_NRPyPN_quasicircular_momenta`, `register_CFunction_cmdline_input_and_parfile_parser`
- [sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py) - `run_sebob`, `calculate_rmse`, `process_input_set`, `__main__`
- [sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py) - `run_sebobv2`, `calculate_rmse`, `process_input_set`, `__main__`

## See Also

- [Examples](index.md)
- [First Wave Equation Run](first-wave-equation-run.md)
- [SEOBNR BOB Generated Library](../infrastructures/bhah/seobnr-bob-generated-library.md)
- [SEBOBv1 JAX Workflow](../infrastructures/jax/sebobv1-jax-workflow.md)
- [Python Codegen](../core/python-codegen.md)
- [Generated Project CI](../validation/generated-project-ci.md)
