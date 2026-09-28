# SEBOBv1 JAX Workflow

> Current SEBOBv1 JAX entry workflow and the implemented SEOBNRv5 aligned-spin coefficient surface. · Status: confirmed
> Up: [JAX](index.md)

## Summary

`nrpy.examples.sebobv1_jax` is a JAX project-generation entrypoint, not a
complete waveform generator. Its implemented surface sets up generation for a
single registered `PyFunction`:
`SEOBNRv5_aligned_spin_coefficients`, which emits generated function text for
SEOBNRv5 aligned-spin calibration/remnant expressions and QNM interpolation
data, and returns a `Commondata` object whose fields, including the final spin
`a_f`, are all registered by the example. The user-facing example text says the current
project initializes SEOBNRv5 aligned-spin coefficients; its class-generation,
full SEBOBv1 port, waveform tests, documentation, calibration, and SEBOBv2
bullets are roadmap items rather than implemented behavior in this route.

## Detail

The entry script first configures generation state. It sets the global
`Infrastructure` parameter to `JAX`, chooses `project_name = "sebobv1_jax"`,
sets `project_dir` to `project/sebobv1_jax`, removes that generated project
directory with `shutil.rmtree(..., ignore_errors=True)`, enables parallel
codegen through the `enable_parallel_codegen` parameter, submits Commondata
metadata through `JAX.commondata.register_commondata_params()`, and registers
the SEOBNRv5 aligned-spin coefficient `PyFunction` through
`JAX.sebob.SEOBNRv5_aligned_spin_coefficients.register_PyFunction_SEOBNRv5_aligned_spin_coefficients()`.
When run as `python -m nrpy.examples.sebobv1_jax`, it then calls
`pcg.do_parallel_codegen()` and
`JAX.jax_project_generator.output_PyFunction_files_and_construct_project()` to
write the generated Python/JAX package.

The Commondata registration call lists mass, spin, frequency, timestep,
coefficient, stop-radius, QNM, final-mass, and final-spin field names. Exact
batch-registration semantics belong to
[Commondata And PyFunction Registry](commondata-and-pyfunction-registry.md);
this workflow page treats the entrypoint's call as generation metadata setup,
not as evidence of generated-waveform runtime behavior.

`register_PyFunction_SEOBNRv5_aligned_spin_coefficients()` is written for the
parallel-codegen registration model. During `pcg_registration_phase()` it
records its own module-qualified call through `pcg.register_func_call()` and
returns without creating the `PyFunction`. During the codegen phase it declares
generated imports for `jax`, `jax.numpy as jnp`, and `.Commondata.Commondata`;
sets the generated function name to
`SEOBNRv5_aligned_spin_coefficients`; takes
`mass_ratio, chi1, chi2, initial_omega, dt, total_mass`; computes `q`, `eta`,
`m1`, `m2`, and `dT`; then instantiates
`SEOBNR_aligned_spin_constants()` to obtain symbolic SEOBNRv5 aligned-spin
expressions.

The coefficient function emits six symbolic assignments through `py_codegen()`:
`a6`, `dSO`, `Delta_t`, `M_f`, `a_f`, and `rISCO`. That path uses the
JAX-oriented Python codegen interface, so the generated expression text is
printer output rather than hand-written Python for those symbolic quantities.
After those assignments, the registration body appends handwritten JAX array
tables for final-spin values and the real/imaginary `(2,2)` QNM data, computes
`omega_qnm` and `tau_qnm` with `jnp.interp`, and emits a generated
`Commondata(...)` return expression populated with masses, spins, initial
frequency, rescaled timestep, coefficients, `rISCO`, `rstop`, QNM values,
`M_f`, and `a_f`. Every keyword in that return is a field that the example
registers in `Commondata`. The printed symbolic expressions contain `Max` and
`Min` as `jnp.maximum` and `jnp.minimum` calls, and integers or rationals that
do not fit in a signed 32-bit integer as floats (see
[CSE And Printer Support](../../core/helpers/cse-and-printer-support.md)), so
the generated module needs no imports beyond `jax`, `jax.numpy`, and
`Commondata`, and its expressions can be evaluated with JAX arrays. The earlier
mismatch between the return and the registered fields is recorded as resolved in
[CONTR-0002](../../contradictions.md#contr-0002).

Claim evidence:
- Claim: The generated `SEOBNRv5_aligned_spin_coefficients` returns `Commondata(...)` using only keywords that `nrpy.examples.sebobv1_jax` registers as `Commondata` fields, including `a_f`. This does not claim that CI imports or calls the generated function.
- Role: descriptive behavior
- Deciding authority: [SEOBNRv5_aligned_spin_coefficients.py](../../../nrpy/infrastructures/JAX/sebob/SEOBNRv5_aligned_spin_coefficients.py), `register_PyFunction_SEOBNRv5_aligned_spin_coefficients` emitted return; [sebobv1_jax.py](../../../nrpy/examples/sebobv1_jax.py), `register_commondata_params` call
- Corroboration: [commondata.py](../../../nrpy/infrastructures/JAX/commondata.py), `register_commondata_params`, rejects unequal list lengths, so a dropped field fails generation; the CONTR-0002 resolution test in [Contradictions](../../contradictions.md#contr-0002) exercises the generated call.

The symbolic source behind this narrow surface is
`SEOBNR_aligned_spin_constants`. Its constructor defines symbolic masses and
spins, computes aligned-spin calibration parameters unless calibration modes
override them, forms `Delta_t`, computes final spin and final mass helpers, and
derives `rISCO` and `rstop`. The JAX workflow currently consumes only the
coefficient/remnant subset needed by
`SEOBNRv5_aligned_spin_coefficients`; it does not assemble SEOBNR dynamics,
inspiral modes, NQC corrections, merger-ringdown waveforms, or mismatch and
calibration utilities.

CI coverage for this route is generation-only. In both `codegen-ubuntu` and
`codegen-mac`, the workflow installs NRPy, generates and `make`-builds many C
example projects, then runs `python -m nrpy.examples.sebobv1_jax` without a
following `make` step for the generated Python/JAX project. Do not cite
generated `project/sebobv1_jax/` files as source evidence unless a maintainer
deliberately freezes and registers such output.

No generated-package install, generated import, coefficient-function call, or
returned-field assertion appears in those configured jobs. Workflow
configuration proves job shape, not generated-package runtime behavior.

## Sources

- [sebobv1_jax.py](../../../nrpy/examples/sebobv1_jax.py) - `project_name`, `Infrastructure`, `enable_parallel_codegen`, `register_commondata_params`, `register_PyFunction_SEOBNRv5_aligned_spin_coefficients`, `output_PyFunction_files_and_construct_project`
- [commondata.py](../../../nrpy/infrastructures/JAX/commondata.py) - `register_commondata_params`, `generate_commondata_dataclass`
- [jax_printer.py](../../../nrpy/helpers/jax_printer.py) - `NRPyJaxPrinter`, `_print_Integer`, `_print_Rational`, `_print_Max`, `_print_Min`
- [SEOBNRv5_aligned_spin_coefficients.py](../../../nrpy/infrastructures/JAX/sebob/SEOBNRv5_aligned_spin_coefficients.py) - `register_PyFunction_SEOBNRv5_aligned_spin_coefficients`
- [SEOBNRv5_aligned_spin_constants.py](../../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_constants.py) - `SEOBNR_aligned_spin_constants`
- [py_codegen.py](../../../nrpy/py_codegen.py) - `PyCodeGen`, `py_codegen`
- [parallel_codegen.py](../../../nrpy/helpers/parallel_codegen.py) - `pcg_registration_phase`, `register_func_call`, `do_parallel_codegen`
- [main.yml](../../../.github/workflows/main.yml) - `codegen-ubuntu`, `codegen-mac`

## See Also

- [JAX](index.md)
- [Project Generation Lifecycle](project-generation-lifecycle.md)
- [Commondata And PyFunction Registry](commondata-and-pyfunction-registry.md)
- [SEOBNR And BOB](../../equations/seobnr/index.md)
- [CSE And Printer Support](../../core/helpers/cse-and-printer-support.md)
- [Generated Project CI](../../validation/generated-project-ci.md)
