# Python Function Registry

> Core route for generated Python/JAX-compatible function objects and registration. · Status: confirmed · Last reconciled: 06-30-2026
> Up: [Core APIs](index.md)

## Summary

`nrpy.py_function` stores generated Python/JAX-compatible function text and metadata in a global registry. It is analogous to the C function registry, but the stored object emits Python function text instead of C function text. A `PyFunction` owns one assembled `full_function` string, and `register_PyFunction()` inserts the object into `PyFunction_dict` under its registered name.

## Detail

`PyFunction` requires `name`, `desc`, and `body` at construction. If any of those fields are missing or empty, construction raises `ValueError` naming the missing attribute. Optional constructor fields are `subdirectory`, `imports`, `prefunc`, `pyfunc_decorators`, `params`, and `postfunc`.

Import validation happens in two stages. Construction rejects a non-list `imports` value when imports are supplied. During `generate_full_function()`, every item in the imports list must be a string; a non-string import raises `TypeError`.

`remove_hashes()` normalizes description text into an indented triple-quoted Python docstring. It dedents the input, strips leading `#` characters and surrounding whitespace from each line, wraps the result in triple quotes, and indents the docstring by four spaces. This is how `desc` becomes the function docstring inside `full_function`.

`indent_body()` adds four spaces to each nonblank body line while preserving that line's existing relative indentation. Blank lines are kept as blank lines, so the generated function body keeps internal spacing without stripping nested indentation.

`subdirectory_depth()` reports the normalized depth of a target subdirectory. It uses path normalization, treats `.` as depth zero, and handles forms such as leading `./`, repeated slashes, and trailing slashes before counting path components.

`generate_full_function()` assembles `full_function` in a fixed order: imports, optional `prefunc`, decorators plus `def`, normalized docstring, indented body, and optional `postfunc`. Imports are emitted one per line followed by a blank line. `prefunc` is stripped of leading and trailing newlines and placed before the function definition. Decorators are emitted immediately before `def`. `postfunc` is stripped of leading and trailing newlines and appended after the function body when present.

`register_PyFunction()` is the global registry boundary. It rejects duplicate names already present in `PyFunction_dict` with `ValueError`; otherwise it constructs a `PyFunction` from the supplied fields and stores it in `PyFunction_dict` under `name`.

Ownership is split between core and infrastructure pages. This page owns the `PyFunction` object contract, registry duplicate behavior, and text assembly rules. JAX infrastructure pages own how registered `PyFunction.full_function` values are written into package files.

## Sources

- [nrpy/py_function.py](../../nrpy/py_function.py) - `PyFunction`, `PyFunction_dict`, `register_PyFunction`

## See Also

- [Core APIs](index.md)
- [C Function Registry](c-function-registry.md)
- [Parallel Codegen Orchestration](helpers/parallel-codegen-orchestration.md)
- [Commondata And PyFunction Registry](../infrastructures/jax/commondata-and-pyfunction-registry.md)
- [Project Generation Lifecycle](../infrastructures/jax/project-generation-lifecycle.md)
