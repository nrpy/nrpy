# SEBOBv1 JAX Workflow

> JAX project-generation entry script `nrpy.examples.sebobv1_jax` and its generated `SEOBNRv5_aligned_spin_coefficients` `PyFunction`. · Status: confirmed
> Up: [JAX](index.md)

## Summary

`nrpy.examples.sebobv1_jax` generates a Python/JAX package containing one registered `PyFunction`, `SEOBNRv5_aligned_spin_coefficients`, which evaluates SEOBNRv5 aligned-spin calibration/remnant expressions and `(2, 2)` QNM interpolation tables and returns a `Commondata` dataclass (including final spin `a_f`). Full SEBOBv1 dynamics, waveform modes, calibration, and SEBOBv2 porting listed in the script's banner remain unimplemented roadmap items.

## Detail

### Project Setup, `PyFunction` Emission, And `Commondata` Fields

`sebobv1_jax.py` sets global `Infrastructure = "JAX"` and `project_name = "sebobv1_jax"`, deletes `project/sebobv1_jax` via `shutil.rmtree(..., ignore_errors=True)`, enables `enable_parallel_codegen`, registers `Commondata` fields through `JAX.commondata.register_commondata_params()` ([Commondata And PyFunction Registry](commondata-and-pyfunction-registry.md)), and calls `register_PyFunction_SEOBNRv5_aligned_spin_coefficients()`. Running `python -m nrpy.examples.sebobv1_jax` executes `pcg.do_parallel_codegen()` and `JAX.jax_project_generator.output_PyFunction_files_and_construct_project()`.

During `pcg_registration_phase()`, `register_PyFunction_SEOBNRv5_aligned_spin_coefficients()` records its module-qualified call via `pcg.register_func_call()` and returns `None`. During codegen it declares `import jax`, `import jax.numpy as jnp`, and `from .Commondata import Commondata`, takes `(mass_ratio, chi1, chi2, initial_omega, dt, total_mass)`, computes `q`, `eta`, `m1`, `m2`, and `dT`, and instantiates `SEOBNR_aligned_spin_constants()` (using only its coefficient/remnant expressions—not inspiral ODEs, NQC terms, or merger-ringdown waveforms).

It then emits six unexpanded symbolic assignments via `py_codegen()` (`a6`, `dSO`, `Delta_t`, `M_f`, `a_f`, `rISCO`, where clamping `1 - 4 nu` in `M_f` prints `jnp.maximum` via `NRPyJaxPrinter` and avoiding `sp.expand()` prevents `float32` coefficient overflow; see [CSE And Printer Support](../../core/helpers/cse-and-printer-support.md) and [Python Codegen](../../core/python-codegen.md)), appends JAX lookup arrays for `(2, 2)` QNM frequencies/damping times interpolated with `jnp.interp` (`omega_qnm`, `tau_qnm`), and returns `Commondata(...)` whose keywords match every field registered in `sebobv1_jax.py`.

Evidence (descriptive behavior): `JAX/sebob/SEOBNRv5_aligned_spin_coefficients.py`; `sebobv1_jax.py`; `JAX/commondata.py::register_commondata_params`.

### Configured CI Execution And Untested Paths

In `codegen-ubuntu` and `codegen-mac`, `.github/workflows/main.yml` runs `python -m nrpy.examples.sebobv1_jax`, installs `project/sebobv1_jax` with `pip install .`, imports `SEOBNRv5_aligned_spin_coefficients`, evaluates it eagerly once on Python scalars, and prints `Commondata` ([`sebobv1_jax` install-and-call step](../../validation/generated-project-ci.md#sebobv1_jax-install-and-call-step)). Because the generated function branches in Python on `eta` and CI applies no `jax.jit`, `float32` array, `pytest` run, GPU/TPU execution, or numerical value assertion, CI verifies only package generation, installation, import, and eager scalar execution; generated `project/sebobv1_jax/` files remain disposable outputs.

Evidence (CI behavior): `.github/workflows/main.yml` (`codegen-ubuntu`, `codegen-mac`).

## Sources

- [sebobv1_jax.py](../../../nrpy/examples/sebobv1_jax.py) - `project_name`, `Infrastructure`, `register_commondata_params`, `output_PyFunction_files_and_construct_project`
- [commondata.py](../../../nrpy/infrastructures/JAX/commondata.py) - `register_commondata_params`, `generate_commondata_dataclass`
- [jax_printer.py](../../../nrpy/helpers/jax_printer.py) - `NRPyJaxPrinter`
- [SEOBNRv5_aligned_spin_coefficients.py](../../../nrpy/infrastructures/JAX/sebob/SEOBNRv5_aligned_spin_coefficients.py) - `register_PyFunction_SEOBNRv5_aligned_spin_coefficients`
- [SEOBNRv5_aligned_spin_constants.py](../../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_constants.py) - `SEOBNR_aligned_spin_constants`
- [py_codegen.py](../../../nrpy/py_codegen.py) - `PyCodeGen`, `py_codegen`
- [parallel_codegen.py](../../../nrpy/helpers/parallel_codegen.py) - `pcg_registration_phase`, `register_func_call`, `do_parallel_codegen`
- [main.yml](../../../.github/workflows/main.yml) - `codegen-ubuntu`, `codegen-mac`

## See Also

- Depends on: [JAX Project Layout And Module Emission](project-generation-lifecycle.md)
- Depends on: [Commondata And PyFunction Registry](commondata-and-pyfunction-registry.md)
- Depends on: [CSE And Printer Support](../../core/helpers/cse-and-printer-support.md)
