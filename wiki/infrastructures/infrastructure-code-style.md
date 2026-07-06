# Infrastructure Code Style

> Infrastructure module structure, C-function registration, generated-code validation, and BHaH generator style rules. · Status: provisional · Last reconciled: 07-06-2026
> Up: [Infrastructures](index.md)

## Summary

Infrastructure code centers on explicit C-function registration routines that
read in generated execution order. Registration functions own their
`CodeParameter` registration, avoid import-time registry mutation, use standard
doctest and trusted-string patterns only when they assert meaningful behavior,
and prefer shared BHaH codegen helpers over hand-written emitted C for ordinary
per-grid or per-point kernels.

## Detail

### Module Organization

An infrastructure Python file usually registers one primary C function through
`register_CFunction_<name>()`. Module filenames use snake_case and directly
describe purpose. `__init__.py` files stay as explicit relative import
aggregators with a flat namespace.

The standard module shape is: module docstring with author info, imports, the
main registration function, optional reusable helper functions, and an optional
`if __name__ == "__main__":` doctest-runner prefix. For modules under
`nrpy/equations/`, that doctest-runner prefix is required as the opening main
block pattern; for infrastructure modules it is strongly encouraged when the
module is runnable.

Registration functions should imitate the C function they register. Declare
`desc`, `cfunc_type`, `name`, `params`, `body`, and related values on separate
lines before calling `cfc.register_CFunction()`. Keep one-off symbolic setup
linear inside the registration routine instead of pushing it into single-use
private helpers.

### CodeParameter Registration Scope

`CodeParameter` registration mutates the global code-parameter registry used
for generated headers, parameter files, and trusted-string tests. Do not
register `CodeParameter`s at module import time. Register each code parameter
inside exactly one core registration or setup function that owns the generated
code needing that parameter.

Do not register the same `CodeParameter` redundantly in multiple modules. If
several C functions need the same parameter, choose one central registration
function and document that downstream registration depends on it. Do not rely
on package imports or `__init__.py` aggregation imports to register parameters
as side effects.

### Doctest And Trusted String Patterns

Infrastructure registration functions use `Doctests:` as the final docstring
section label immediately before `>>>` lines. Older variants such as
`Doctest:` and `DocTests:` exist; use `Doctests:` in new code.

For stable generated C text, the standard doctest pattern imports
`validate_strings` and `clang_format` inside the doctest, clears
`cfc.CFunction_dict`, registers the target function, formats
`full_function`, and compares against a trusted file with `validate_strings`.
Use `file_ext="cu"` for CUDA and `"c"` otherwise.

Golden-output doctests are most appropriate in `nrpy/infrastructures/*/*.py`,
where registration functions are a core public interface. Avoid generated C
trusted files for C functions whose bodies are dominated by large SymPy-derived
kernels; SymPy and codegen details make those exact strings brittle. Prefer
symbolic validation or cheaper structural checks for such functions.

Do not add trivial doctests whose main assertion is that Python-based C/C++
generation ran, a `CFunction` registered, a registration function returned, or
an object could be called. Missing doctests are acceptable when the only
practical doctest would be such a low-value generated-text spot check. If a
`Doctests:` block contains `# FIXME`, complete it in a follow-up or remove it
when the real doctest is ready.

### Parallel Codegen Registration

Registration functions that participate in `nrpy.helpers.parallel_codegen`
discovery use the standard early-return guard. During the registration phase,
they call `pcg.register_func_call(...)` with the fully qualified function name
and locals, return `None`, and only perform actual C registration outside that
guard. These functions return `Union[None, pcg.NRPyEnv_type]`; registration
functions that do not support parallel/discovery codegen omit the pattern and
return `None`.

### Black Suppression

Use `# fmt: off` and `# fmt: on` sparingly, only when Black would destroy
intentional alignment. One accepted use is a compact group of `CodeParameter`
registrations inside the single registration function that owns those
parameters.

### C Function Registration And Helper Functions

Register C functions through `cfc.register_CFunction()` with explicit standard
fields: `subdirectory`, `includes`, optional `prefunc`, `desc`, `cfunc_type`,
`name`, `params`, `include_CodeParameters_h`, `body`, and optional `postfunc`.
Helper C functions emitted before the main function are generated as strings
and concatenated into `prefunc`.

Inline C and Python functions that do not save real code by existing
separately. Private/module-local Python helpers, including leading-underscore
helpers, need at least two real call sites unless an external API, callback, or
test harness requires a named function. Separate C helper functions only when
they are called from multiple locations or have more than roughly 10-15 lines
of actual logic. Do not keep functions whose docstring is longer than the body
or that simply return one expression.

### BHaH Symbolic Codegen Rules

For new BHaH infrastructure generators that emit ordinary per-grid or
per-point kernels, prefer infrastructure helpers over handwritten emitted C.
Use `BHaH.simple_loop.simple_loop()` for standard grid loops. When symbolic
expressions depend on registered gridfunctions, prefer `ccg.c_codegen(...,
automatically_read_gf_data_from_memory=True)` with expected array aliases such
as `in_gfs` instead of writing repetitive memory loads by hand.

Keep transformations symbolic until `c_codegen()` whenever practical. Avoid
string-based replacement of symbolic expressions or generated C when the same
result can be represented symbolically. Do not introduce new
`#include "set_CodeParameters.h"` lines inside generated function bodies for
this class of infrastructure code; derive needed `params` and `commondata`
locals from symbolic expressions using
`get_params_commondata_symbols_from_expr_list()` and
`generate_definition_header()`.

Avoid routine post-registration mutation of `cfc.CFunction_dict[...]` bodies.
Prefer explicit extension hooks or parameters in the shared registration helper
that owns the emitted C. Build symbolic expressions, `expr_list`, `lhs_list`,
tensor declarations, and registration-time metadata near the `c_codegen()` or
`simple_loop()` call that consumes them. Append emitted C bodies top-to-bottom
with `body += ...` so Python assembly order matches generated execution order.

Add short comments at jarring transitions, such as registering diagnostic
gridfunctions before parity tables or constructing a symbolic kernel
immediately before `c_codegen()`.

### Runtime Data, Gridfunctions, And Generated C Idioms

Common generated C parameter structs are `commondata_struct *restrict
commondata`, `griddata_struct *restrict griddata`, `params_struct *restrict
params`, and `bc_struct *restrict bcstruct`.

BHaH gridfunction names encode tensor type and indices. Examples include scalar
names such as `HHGF`, `VVGF`, `WWGF`, `TRKGF`, and `CFGF`; symmetric rank-2
names such as `HDD00GF`, `HDD01GF`, and `HDD11GF`; traceless extrinsic
curvature names such as `ADD00GF`; and derivative names such as
`SRC_PARTIAL_D_HDD000GF` and `SRC_PARTIAL_D_WW0GF`.

Standard groups are `EVOL` for time-evolved quantities, `AUXEVOL` for
auxiliary evolution quantities, and `AUX` for diagnostics. BHaH also uses the
core `DIAG` group in generated gridfunction registration.

Indexing macros include `IDX3(i0, i1, i2)`, `IDX4(gf, i0, i1, i2)`,
`IDX4pt(gf, idx3)`, and `IDX2(i1, i2)`. Function-local custom macros such as
`SRC_IDX4`, `DST_IDX4`, and `EX_IDX4` adapt indexing to source, destination,
or external data layouts.

SymPy integration for generated code uses `ccg.c_codegen()` to turn symbolic
expressions into C assignments, including finite-difference codegen when
enabled. Tensor derivatives are declared with `ixp.declarerankN()` before being
passed to codegen.

BHaHAHA-style error handling uses enum error codes such as `BHAHAHA_SUCCESS`
and module-specific names, exposes messages through `bah_error_message()`,
checks integer return codes immediately, and returns early when
`commondata->error_flag != BHAHAHA_SUCCESS`. Parallel regions use
`#pragma omp critical` to set shared error flags.

## Sources

- [original-agents.md](../../raw/source-docs/original-agents.md) - `## Infrastructure Code Rules`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `### Doctests`, `### Parallel Codegen Pattern`, `### Black Suppression`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `### C Function Registration from Python`, `### BHaH Symbolic Codegen Rules`, `### Inlining Rules`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `### Standard Struct Pointer Params`, `### Gridfunction Naming / Grouping`, `### Memory / Error Handling`

## See Also

- Parent: [Infrastructures](index.md)
- Depends on: [C Function Registry](../core/c-function-registry.md)
- Depends on: [Gridfunctions And Parameters](../core/gridfunctions-and-parameters.md)
- Depends on: [Parallel Codegen Orchestration](../core/helpers/parallel-codegen-orchestration.md)
- See also: [C And Embedded C Style](../architecture/c-and-embedded-c-style.md)
