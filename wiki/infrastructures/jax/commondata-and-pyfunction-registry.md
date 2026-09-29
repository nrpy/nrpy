# Commondata And PyFunction Registry

> Consumer view of JAX shared-data and generated-function registries. · Status: confirmed
> Up: [JAX](index.md)

## Summary

JAX project generation is driven by two generation-time registries:
`commondata_params_dict` records fields for the generated `Commondata`
dataclass, and `PyFunction_dict` records generated Python functions that will be
written into the package source tree. `Commondata` registration defines shared
data fields by name, type, default value, and description; `PyFunction`
registration supplies function text and placement metadata. Parallel codegen is
relevant here only when registration is deferred into worker calls and merged
back before the JAX project writer consumes `PyFunction_dict`.

## Detail

`nrpy.infrastructures.JAX.commondata` owns the JAX shared-data registry.
`commondata_params_dict` is keyed by parameter name, and each entry stores
`dtype`, `default`, and `description` metadata for one generated dataclass
field. `register_commondata_param()` rejects duplicate names before inserting a
new entry, so direct single-field registration is duplicate-checked at the
registry boundary.

`register_commondata_params()` is a batch convenience wrapper. It first checks
that `names`, `dtypes`, `defaults`, and `descriptions` have the same length and
raises `ValueError` before registering any field if they do not. It then
iterates over `zip(names, dtypes, defaults, descriptions)` and calls
`register_commondata_param()` for each tuple, so no listed field is silently
dropped. `nrpy.examples.sebobv1_jax` supplies fourteen entries in each list,
including the final-spin field `a_f` that its emitted coefficient function
passes to `Commondata(...)`.

Claim evidence:
- Claim: `register_commondata_params()` raises `ValueError` and registers no field when its four input lists differ in length; otherwise it registers every listed field. This says nothing about duplicate names beyond the per-field check in `register_commondata_param()`.
- Role: descriptive behavior
- Deciding authority: [nrpy/infrastructures/JAX/commondata.py](../../../nrpy/infrastructures/JAX/commondata.py), `register_commondata_params`
- Corroboration: `none available`; the only targeted check is the `register_commondata_params` doctest in the deciding module, which exercises a length mismatch but is not a separate source.

`generate_commondata_dataclass()` turns the current registry contents into the
generated `Commondata.py` module text. It emits `from dataclasses import
dataclass`, an `@dataclass class Commondata`, and one field line per registered
entry using the stored name, dtype, default, and optional description comment,
ending the module with a newline. When no fields are registered, the generated
class body is `pass`.

`nrpy.py_function` owns generated Python function registration. A `PyFunction`
stores the generated function's registered name, import lines, optional
prefunction text, description, optional decorators, parameter list, body,
optional postfunction text, and optional subdirectory. On construction it also
builds `full_function`, the complete generated function text. For the JAX
consumer path, the most important fields are the registered `name`, `imports`,
`pyfunc_decorators`, `body`, `subdirectory`, and `full_function`; deeper
formatting and validation behavior belongs to the core helper pages.
`register_PyFunction()` rejects duplicate names in `PyFunction_dict`, then
stores a new `PyFunction` under that name.

The JAX project writer consumes both registries while constructing a generated
package. `_generate_python_files_and_init()` iterates over `PyFunction_dict`,
creates any requested subdirectories under `src/<project_name>/`, writes one
module per registered function using `pyfunc.full_function`, emits package
`__init__.py` files, then writes `Commondata.py` from
`generate_commondata_dataclass()`. It imports `commondata_params_dict` in the
same module and clears it in the local test-mode path, but the write path
itself obtains generated dataclass text through `generate_commondata_dataclass()`.

Parallel codegen does not change the JAX package layout or the meaning of
`PyFunction` entries. Its role is the deferred-registration path:
`pcg_registration_phase()` lets a registration function choose to record a
deferred call instead of doing generation immediately, `do_parallel_codegen()`
runs registered calls in workers, and `unpack_NRPy_environment_dict()` merges
the returned environment dictionaries back into parent globals, including
`PyFunction_dict`. After that merge, JAX project generation sees the same
`PyFunction_dict` shape as it would after serial registration.

## Sources

- [nrpy/infrastructures/JAX/commondata.py](../../../nrpy/infrastructures/JAX/commondata.py) - `commondata_params_dict`, `register_commondata_param`, `register_commondata_params`, `generate_commondata_dataclass`
- [nrpy/infrastructures/JAX/jax_project_generator.py](../../../nrpy/infrastructures/JAX/jax_project_generator.py) - `_generate_python_files_and_init`, `output_PyFunction_files_and_construct_project`
- [nrpy/py_function.py](../../../nrpy/py_function.py) - `PyFunction`, `PyFunction_dict`, `register_PyFunction`
- [nrpy/helpers/parallel_codegen.py](../../../nrpy/helpers/parallel_codegen.py) - `NRPyEnv`, `pcg_registration_phase`, `do_parallel_codegen`, `unpack_NRPy_environment_dict`
- [nrpy/examples/sebobv1_jax.py](../../../nrpy/examples/sebobv1_jax.py) - `JAX.commondata.register_commondata_params`

## See Also

- [JAX](index.md)
- [Project Generation Lifecycle](project-generation-lifecycle.md)
- [Python Function Registry](../../core/python-function-registry.md)
- [Python Codegen](../../core/python-codegen.md)
- [CSE And Printer Support](../../core/helpers/cse-and-printer-support.md)
- [Parallel Codegen Orchestration](../../core/helpers/parallel-codegen-orchestration.md)
