# Project Generation Lifecycle

> Explain how the JAX project generator turns registered `PyFunction` objects into a generated Python package. · Status: confirmed · Last reconciled: 2026-06-30
> Up: [JAX](index.md)

## Summary

`output_PyFunction_files_and_construct_project()` is the JAX infrastructure's
package emitter. It validates the requested project path and Python package
name, creates the project root, writes registered `PyFunction` modules and
`Commondata.py` under `src/<project_path.name>/`, then emits packaging
metadata and a minimal import smoke test. In the current `sebobv1_jax` route,
the resolved `project_dir` basename and package name both equal
`sebobv1_jax`. The generated `project/<name>/` tree is an output product, not
source evidence for KB claims, unless maintainers deliberately freeze and
register a generated artifact.

## Detail

The input side is the `PyFunction` registry. `register_PyFunction()` stores
each generated Python function in `PyFunction_dict` by name, and each
`PyFunction` carries the target subdirectory, imports, optional pre-function
text, decorators, function name, parameters, body, optional post-function text,
and the generated `full_function` string. The JAX project generator consumes
that registry rather than discovering Python functions from generated files.

`output_PyFunction_files_and_construct_project()` runs the lifecycle in a fixed
order. It resolves `project_dir` to a `Path`, calls `_validate_inputs()`,
creates the project directory, delegates source emission to
`_generate_python_files_and_init()`, and finally calls
`_generate_project_metadata()`. Failures from those stages are logged and
wrapped as `RuntimeError`.

`_validate_inputs()` checks the package name before generation. The current
checks require `project_name` to be a string, allow only alphanumeric
characters and underscores, reject names that start with a digit, reject an
existing non-directory project path, and create a temporary `.write_test` file
to verify the resolved output directory can be created and written.

`_generate_python_files_and_init()` writes the generated package under
`src/<project_path.name>/`, where `project_path` is the resolved
`project_dir`. For every registered `PyFunction`, it optionally applies
`lib_function_prefix`, mutates the function name and regenerated
`full_function`, creates the requested subdirectory, creates package
`__init__.py` files for generated directories, and writes one
`<function_name>.py` module containing a module docstring plus the generated
function text. It always emits `Commondata.py` from
`generate_commondata_dataclass()`. The top-level package `__init__.py` imports
subdirectory functions, direct package functions, and `Commondata`; the
generated `__all__` is currently written only when direct imports exist and then
contains those direct imports plus `Commondata`.

`_generate_project_metadata()` emits the rest of the generated Python project:
`pyproject.toml`, `setup.cfg`, `requirements.txt`, `README.md`, `.gitignore`,
`tests/__init__.py`, and `tests/test_basic.py`. Treat this metadata as current
generator output, not production packaging advice: `setup.cfg` and `README.md`
contain placeholder author, URL, and usage text, and the generated test only
imports the package and checks that `__version__` exists.

`nrpy.examples.sebobv1_jax` is the current user-facing JAX generation route. It
sets `Infrastructure` to `JAX`, chooses `project_name = "sebobv1_jax"`, clears
`project/sebobv1_jax`, enables parallel codegen, registers `Commondata`
parameters and the SEOBNRv5 aligned-spin coefficient `PyFunction`, then calls
`pcg.do_parallel_codegen()` and
`JAX.jax_project_generator.output_PyFunction_files_and_construct_project()`
when run as a module. README guidance lists this as JAX project generation, and
the generated-project CI jobs run `python -m nrpy.examples.sebobv1_jax`; unlike
the C examples in the same jobs, that JAX route is not followed by a generated
`make` build.

Generated files under `project/<name>/`, including generated JAX package files,
are products of handwritten generators. Cite the source generator, registry,
example, README, and CI workflow for behavior rather than citing transient
generated project output. Printer-level JAX expression emission belongs to
[CSE And Printer Support](../../core/helpers/cse-and-printer-support.md), which
owns `py_codegen()` and `NRPyJaxPrinter`; this page only documents how the
already-registered Python functions are assembled into a package.

## Sources

- [jax_project_generator.py](../../../nrpy/infrastructures/JAX/jax_project_generator.py) - `_validate_inputs`, `_generate_python_files_and_init`, `_generate_project_metadata`, `output_PyFunction_files_and_construct_project`
- [py_function.py](../../../nrpy/py_function.py) - `PyFunction`, `PyFunction_dict`, `register_PyFunction`
- [sebobv1_jax.py](../../../nrpy/examples/sebobv1_jax.py) - `project_name`, `project_dir`, `output_PyFunction_files_and_construct_project`
- [README.md](../../../README.md) - `Project Families and Example Generators`, `What Gets Generated?`, `Contributor Setup`
- [main.yml](../../../.github/workflows/main.yml) - `codegen-ubuntu`, `codegen-mac`
- [Generated Output Boundaries](../../architecture/generated-output-boundaries.md) - generated project boundary
- [Generated Project CI](../../validation/generated-project-ci.md) - `codegen-ubuntu`, `codegen-mac`
- [CSE And Printer Support](../../core/helpers/cse-and-printer-support.md) - `py_codegen`, `NRPyJaxPrinter`

## See Also

- [JAX](index.md)
- [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
- [Generated Project CI](../../validation/generated-project-ci.md)
- [CSE And Printer Support](../../core/helpers/cse-and-printer-support.md)
