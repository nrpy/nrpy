# Waveform JAX PN Generators

> Map SEOBNR/SEBOB waveform, JAX, and PN momentum example generators to generated project names, dependency classes, and example-tree consistency checks. · Status: confirmed
> Up: [Examples](index.md)

## Summary

The `seobnrv5_aligned_spin_inspiral` route creates GSL-backed BHaH C waveform
projects whose project name is selected from the SEOBNRv5 approximant flag and
optional calibration suffix. `sebobv2` creates one GSL-backed BHaH C waveform
project named `sebobv2` (with optional `--projected-attachment-diagnostics`
generation and parfile `use_projected_attachment`). `sebobv1_jax` creates a Python/JAX project named
`sebobv1_jax`, but the current source-backed code only initializes
SEOBNRv5 aligned-spin coefficients. `nrpypn_quasicircular_momenta` creates a
BHaH C utility project named `nrpypn_quasicircular_momenta` for 3.5PN
quasicircular binary-black-hole momenta.

The SEOB/SEBOB consistency scripts are `nrpy/examples/tests/` validation helpers. They
build trusted and current generated executable directories, run matching
waveform outputs, compute amplitude-plus-phase RMSE-style differences, and
require the current-vs-trusted median error to stay within a perturbation-based
roundoff baseline. For calibration-mode SEOBNRv5 approximants, the
SEOBNRv5-family script additionally requires the trusted and current
executables to succeed or fail consistently on each input; any asymmetry
(one succeeds, the other crashes) fails the run independently of the median
check.

## Detail

### SEOBNRv5 And SEBOBv2 C Waveform Generators

`nrpy.examples.seobnrv5_aligned_spin_inspiral` is the flag-driven SEOBNRv5
aligned-spin generator. If no approximant flag is supplied, it defaults to
`-seobnrv5_bob`. Normal use chooses one approximant flag, then optionally adds
one of three mutually exclusive calibration/selector flags:
`-calibration_no_spin`, `-calibration_spin`, or `-nrpy_calibrated`
(`-nrpy_calibrated` selects the NRPy-calibrated `a6`/`Delta_t_NS` fit in
production mode and cannot be combined with either calibration flag):

| Flag family | Generated project name | NQC choice | Merger-ringdown choice | Dependency class |
| --- | --- | --- | --- | --- |
| `-seobnrv5_bob` or no approximant flag | `seobnrv5_bob` | BOB-informed NQC | BOB-informed merger-ringdown | BHaH C project with GSL |
| `-seobnrv5_nrnqc_bob` | `seobnrv5_nrnqc_bob` | native numerical-relativity NQC | BOB-informed merger-ringdown | BHaH C project with GSL |
| `-seobnrv5_nrpy` | `seobnrv5_nrpy` | native numerical-relativity NQC | native SEOBNRv5 merger-ringdown | BHaH C project with GSL |

Evidence (descriptive behavior): `seobnrv5_aligned_spin_inspiral.py::argparse` flag block and the `project_name` base-name assignment in the module `__main__` body.

Project-name expansion is string-based. The script starts from the approximant
base name, then appends `_calibration_no_spin` for `-calibration_no_spin`,
`_calibration_spin` for `-calibration_spin`, or `_nrpy_calibrated` for
`-nrpy_calibrated`. It rejects combining the two calibration flags, and
rejects combining `-nrpy_calibrated` with either calibration flag. Thus
`-seobnrv5_bob -calibration_no_spin` writes
`project/seobnrv5_bob_calibration_no_spin/`, `-seobnrv5_nrpy -calibration_spin`
writes `project/seobnrv5_nrpy_calibration_spin/`, and
`-seobnrv5_bob -nrpy_calibrated` writes
`project/seobnrv5_bob_nrpy_calibrated/`. These flags also flow into SEOBNRv5
aligned-spin coefficient registration, so they are not only naming suffixes.

Evidence (descriptive behavior): `seobnrv5_aligned_spin_inspiral.py::project_name` string-assignment block and the `calibration_no_spin`/`calibration_spin`/`nrpy_calibrated` `ValueError` guards.

The SEOBNRv5 waveform projects are GSL-backed C projects. The generator sets
`Infrastructure` to `BHaH`, registers C functions for commondata I/O,
root-finding, initial conditions, ODE integration, spline/interpolation-backed
waveform processing, NQC corrections, merger-ringdown, and IMR assembly, then
emits GSL headers and a Makefile using `gsl-config --cflags` and
`gsl-config --libs`. Its command-line/parfile inputs are `mass_ratio`, `chi1`,
`chi2`, `initial_omega`, `total_mass`, and `dt`.

`nrpy.examples.sebobv2` is a separate GSL-backed BHaH C waveform generator with
fixed `project_name = "sebobv2"` and the same 6-argument scalar CLI (`mass_ratio`,
`chi1`, `chi2`, `initial_omega`, `total_mass`, `dt`). Passing
`--projected-attachment-diagnostics` at generation time sets
`enable_projected_reference_diagnostics = True` so `main` prints
`NRPY_REFERENCE`, `NRPY_SPINS`, `NRPY_ATTACH`, `NRPY_FINE`, `NRPY_NQC`,
`NRPY_CORRECTED`, and `NRPY_IMR` diagnostic records to `stderr`.
`SEOBNRv5_quasi_precessing_spin_coefficients`
also registers parfile-only `bool` `use_projected_attachment` (default `False`,
omitted from `cmdline_inputs` to preserve the 6-argument CLI): when true in `.par`,
`main` rejects CLI overrides (`argc > 2`), sets `chi1 = chi1_z` and `chi2 = chi2_z`,
and calls `SEOBNRv5_projected_attachment_reference` before special-amplitude and
`SEBOBv2_NQC_corrections` stages. Current stdout prints the `(2,2)` strain while
higher-mode inspiral storage and optional coprecessing-rotation sandbox helpers
remain partial.

Evidence (descriptive behavior): `sebobv2.py::register_CFunction_main_c`; `BHaH/seobnr/SEOBNRv5_quasi_precessing_spin_coefficients.py`; `BHaH/seobnr/SEOBNRv5_projected_attachment_reference.py`.

### JAX Package, PN Momentum Utility, And Consistency Scripts

`nrpy.examples.sebobv1_jax` is not a GSL-backed waveform executable. It sets
`Infrastructure` to `JAX`, writes `project/sebobv1_jax/`, registers Commondata
metadata, registers the `SEOBNRv5_aligned_spin_coefficients` Python/JAX
function, runs parallel codegen, and emits a generated Python package. Its own
user-facing text says the current project initializes SEOBNRv5 aligned-spin
coefficients. Class generation, the full SEBOBv1 JAX port, waveform-generation
tests, documentation, mismatch/calibration utilities, and possible SEBOBv2
JAX work are roadmap items in the example source, not completed workflow
features. Existing JAX infrastructure pages document generated-package emission
and the narrower Commondata/PyFunction details.

The generator supplies 14 Commondata names, dtypes, defaults, and
descriptions, and batch registration rejects lists of unequal length, so
generated `Commondata.py` includes the `a_f` field that the emitted coefficient
function passes to `Commondata(...)`. GitHub codegen generates and installs
the package, imports the coefficient function, and calls it once with
Python-scalar inputs; it does not run the generated basic test or check any
returned value, JAX accelerator, or numerical result.

[Generated Project CI](../validation/generated-project-ci.md#sebobv1_jax-install-and-call-step)
gives the exact configured `sebobv1_jax` CI step and its evidence.

The generated metadata declares `jax`, `jaxlib>=0.4.0`, and
`numpy>=1.21.0`; official JAX [Installation](https://docs.jax.dev/en/latest/installation.html)
under `Supported platforms` remains authority for current platform-specific JAX
installation, not NRPy's workflow.

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

For calibration-mode approximants (`_calibration_no_spin`/`_calibration_spin`
in the directory name), the SEOBNRv5-family script treats each input set's
trusted/current executable run as its own success/failure outcome instead of
letting a crash propagate: an input set where both executables succeed or
both fail contributes no regression signal (a consistent failure is excluded
from the median rather than aborting the script), but an input set where one
succeeds and the other crashes is reported as a regression and fails the run
regardless of the median comparison. This accounts for calibration-mode
approximants depending on externally supplied calibration coefficients that
this consistency script runs with placeholder parfile defaults. `sebobv2_consistency_check`
has no calibration-mode variants and is unaffected by this behavior.

Evidence (descriptive behavior): `sebob_consistency_check.py::process_input_set`; `sebob_consistency_check.py::__main__`.

`nrpy.examples.tests.sebobv2_projected_attachment_check` (`--project`, `--baseline`,
configured in `.github/workflows/sebobv2-projected-attachment.yml`) compiles an
embedded standalone RKF45 C integrator (`independent_reference.c`) in `main` to
verify `10M` orbital and spin reference recovery (`NRPY_REFERENCE`, `NRPY_SPINS`,
including starts inside `10M`) and `argc > 2` CLI-override rejection when
`use_projected_attachment = True`, runs `check_handoff` (attachment-time
retention between `NRPY_ATTACH` and `NRPY_NQC` plus `(2,2)` QNM `omega`/`tau`
equality), `check_matching` (NQC attachment), and `check_phase_overlap`
(merger-ringdown phase agreement), and with `--baseline` runs
`use_projected_attachment = False` to verify byte-identical aligned-spin scalar
`stdout`. Together with `.github/workflows/main.yml`
(which compares production, `-calibration_no_spin`, and `-calibration_spin`
SEOBNRv5 variants plus `sebobv2`, while building `-nrpy_calibrated` only), these
checks cover the tested inputs without claiming general waveform accuracy.

Evidence (CI behavior): `sebobv2_projected_attachment_check.py::main`; `sebobv2_projected_attachment_check.py::check_handoff`; `.github/workflows/sebobv2-projected-attachment.yml`; `.github/workflows/main.yml`.

## Sources

- [seobnrv5_aligned_spin_inspiral.py](../../nrpy/examples/seobnrv5_aligned_spin_inspiral.py) - `argparse` flags, `project_name`, `register_CFunction_main_c`, `register_CFunction_cmdline_input_and_parfile_parser`, `output_CFunctions_function_prototypes_and_construct_Makefile`
- [sebobv2.py](../../nrpy/examples/sebobv2.py) - `project_name`, `--projected-attachment-diagnostics`, `register_CFunction_main_c`, `register_CFunction_cmdline_input_and_parfile_parser`, `output_CFunctions_function_prototypes_and_construct_Makefile`
- [SEOBNRv5_quasi_precessing_spin_coefficients.py](../../nrpy/infrastructures/BHaH/seobnr/SEOBNRv5_quasi_precessing_spin_coefficients.py) - `register_CFunction_SEOBNRv5_quasi_precessing_spin_coefficients` (`use_projected_attachment`)
- [SEOBNRv5_projected_attachment_reference.py](../../nrpy/infrastructures/BHaH/seobnr/SEOBNRv5_projected_attachment_reference.py) - `register_CFunction_SEOBNRv5_projected_attachment_reference`
- [sebobv1_jax.py](../../nrpy/examples/sebobv1_jax.py) - `project_name`, `Infrastructure`, `register_commondata_params`, `register_PyFunction_SEOBNRv5_aligned_spin_coefficients`, `output_PyFunction_files_and_construct_project`
- [JAX Installation](https://docs.jax.dev/en/latest/installation.html) - `Supported platforms`
- [nrpypn_quasicircular_momenta.py](../../nrpy/examples/nrpypn_quasicircular_momenta.py) - `project_name`, `register_CFunction_main_c`, `register_CFunction_NRPyPN_quasicircular_momenta`, `register_CFunction_cmdline_input_and_parfile_parser`
- [sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py) - `run_sebob`, `calculate_rmse`, `process_input_set`, `__main__`
- [sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py) - `run_sebobv2`, `calculate_rmse`, `process_input_set`, `__main__`
- [sebobv2_projected_attachment_check.py](../../nrpy/examples/tests/sebobv2_projected_attachment_check.py) - `check_handoff`, `check_matching`, `check_phase_overlap`, `main`
- [main.yml](../../.github/workflows/main.yml) - `sebob-consistency-test`, `sebobv2-consistency-test`
- [sebobv2-projected-attachment.yml](../../.github/workflows/sebobv2-projected-attachment.yml) - `projected-reference-and-matching`

## See Also

- Depends on: [SEOBNR BOB Generated Library](../infrastructures/bhah/seobnr-bob-generated-library.md)
- Depends on: [SEBOBv1 JAX Workflow](../infrastructures/jax/sebobv1-jax-workflow.md)
- Depends on: [Python Codegen](../core/python-codegen.md)
