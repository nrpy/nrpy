# NRPy Agent Instructions

Repo agent rules, compressed. CI enforce many.

## Scope

- NRPy: Python/SymPy -> generated C.
- No binary files in PRs unless maintainers approve.
- No images, archives, compiled artifacts, or other non-text assets unless maintainers approve.

## Required Checks

- Python changes: `black .` before commit.
- Every modified Python file: `./.github/single_file_static_analysis.sh <path-to-file.py>` before commit.
- Exception: trusted-value files in `.../tests/*.py`; generated reference values only, skip single-file static analysis.
- No hand-picked subset only. Contributions pass project static analysis before merge.
- Run Python examples directly from this directory: append `.` to `PYTHONPATH`: `export PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}."`

Static-analysis script checks:

- `black --check`
- `isort --check-only`
- `mypy --strict --allow-untyped-calls`
- `pylint` threshold `>= 9.91/10`
- `pydocstyle`
- `darglint -v 2`
- doctests via `python3 <file>`

## Python Style

### Formatting

- Black; 4 spaces.

### Naming

- Classes: `PascalCase`; constants: `UPPER_SNAKE_CASE` or leading `_`; private helpers: leading `_`.
- Boolean-like params positive.

Prefer:

- `enable_feature`
- `include_header`
- `allow_resize`
- `has_ghost_zones`

Avoid:

- `disable_feature`
- `omit_header`
- `forbid_resize`
- `no_ghost_zones`

### Imports

- Order:
  1. standard library
  2. third-party
  3. local NRPy
- Blank line between groups; follow `isort`.

Canonical aliases mandatory:

- `import sympy as sp`
- `import nrpy.indexedexp as ixp`
- `import nrpy.params as par`
- `import nrpy.grid as gri`
- `import nrpy.reference_metric as refmetric`
- `import nrpy.c_function as cfc`
- `import nrpy.c_codegen as ccg`
- `import nrpy.helpers.parallel_codegen as pcg`

- No invented aliases.

### `__init__.py`

- No module docstring, comments, or executable code beyond explicit relative import aggregation.
- Flat namespace with explicit relative imports.

### Python Docstrings

- Triple double quotes.
- Sphinx/reST style:
  - `:param name:`
  - `:return:`
  - `:raises ExceptionType:`
- `:return:`, not `:returns:`.
- No Google-style `Args:` / `Returns:` / `Raises:`.
- No type info in `:param` text.

### Python String Literals

- Prefer triple double-quoted multiline strings for multiline literals: `desc`, `body`, `prefunc`, `postfunc`, generated-code snippets, long messages, validation strings.
- `r"""..."""` when backslashes stay literal; `rf"""..."""` only when interpolation needed.
- Keep interpolated expressions simple and Python-3.7-friendly; move complex expressions to locals first.
- No replacing clear static multiline literals with adjacent string fragments or `"\n".join(...)` unless readability clearly improves.
- Single-line ordinary strings stay normal double-quoted strings unless triple quotes help.

### Module Docstrings

Every non-`__init__.py` Python file top docstring:

```python
"""
<Description paragraphs.>

Author: Name
        email obfuscated
"""
```

Multiple authors:

```python
"""
<Description.>

Authors: Name One
         email
         Name Two
         email
"""
```

Rules:

- `Author:` for one author; `Authors:` for multiple.
- No authors not already in file.
- Optional leading filename comment only exact form `# <relative path from nrpy root>.py`
- `__init__.py` and trusted test-vector files never get module docstrings.
- Only `Author:` / `Authors:` keys.

Avoid:

- singular `Author:` with multiple names
- `Email:` or `Contributor:`
- mixed metadata styles in one file

### Type Hints

- Many type hints; always include return annotations, including `-> None`.
- `typing` forms: `Dict`, `List`, `Optional`, `Tuple`, `Union`, `cast`; `typing_extensions.Literal` for constrained strings.
- `# type: ignore` only when needed, mainly third-party typing gaps.
- Avoid `Any`; if unavoidable, justify inline.
- No builtin generics like `list[X]`, `dict[X, Y]`, `tuple[X, ...]`; no `X | None`.
- No `from __future__ import annotations` unless real need.

### Comments

In docstrings:

- Inline `Note:` okay.
- reST `.. note::` okay.
- No mixing both styles in one docstring.

Procedural comment style:

- `# Step N:`
- substeps: `# Step 1.a:`
- module preamble: `# Step P1:`
- use `Step`, not `STEP`

### `if __name__ == "__main__":`

- Runnable non-test, non-`__init__.py` modules under `nrpy/equations/` and subdirs start `__main__` block with exact doctest-runner below.
- Runnable modules under `infrastructures/*/*.py`: strongly encouraged.
- Other runnable modules: recommended when useful.
- Outside `nrpy/infrastructures/*/*.py`, Python-based C/C++ codegen doctests discouraged unless signal beats cheaper checks.

Use this exact block:

```python
if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
```

- Imports used only here stay inside block.
- In `nrpy/equations/**`, symbolic validation/trusted-results generation may follow this prefix in same `__main__` block.
- Allowed follow-on steps there: `ve.process_dictionary_of_expressions(...)`, `ve.compare_or_generate_trusted_results(...)`, coordinate-system loops, validation-key loops.
- Allowance for established equation-module validation workflow, not arbitrary script logic.

### General Python Organization

- Classes group related functionality; module functions handle procedural codegen.
- Common pattern: `register_CFunction_*()`.
- Doctests live in docstrings.
- Raise `ValueError` for invalid input; use `warnings.warn()` for non-critical issues.
- Validate explicitly in constructors.

## Equation Setup Rules

### SymPy

- Avoid `sp.simplify()` in equation-building code; allowed only in tests, validation, or explicit identity checks.
- `nrpy.helpers.cached_functions.cached_simplify()` narrow exception: sparingly for small local scalar subexpressions when it helps symbolic tractability, codegen stability, or readability, especially in `nrpy/equations/general_relativity`.
- No cached simplification on whole tensors or large assembled RHS expressions when plain symbolic construction works.
- If need not obvious nearby, add short comment explaining why cached simplification warranted.
- No `sp.subs()` / `sp.replace()` for pattern-based transformation.

Allowed `.subs()` cases:

- coordinate substitution in elliptic source terms
- face-value substitution in GRHD fluxes
- systematic `nrpyABS` to `sp.Abs`
- surface radius substitution in horizon modules
- evaluation at specific parameter value like `.subs(t, t_attach)`

Preferred:

- init accumulators with `sp.sympify(0)` / `sp.sympify(1)`
- `sp.Rational()` for exact fractions
- `sp.symbols(..., real=True)` for scalars

### Indexed Expressions

- `nrpy.indexedexp` as `ixp` for tensors.

Common helpers:

- `ixp.zerorank1()` through `ixp.zerorank4()`
- `ixp.declarerank1()` through `ixp.declarerank4()`
- `ixp.symm_matrix_inverter2x2/3x3/4x4()`
- `ixp.LeviCivitaSymbol_dim3_rank3()`
- `ixp.LeviCivitaTensorUUU_dim3_rank3()`

Symmetry strings:

- `"sym01"` for metrics
- `"sym12"` for derivatives of symmetric tensors
- `"sym01_sym23"` for Riemann-like tensors
- `"nosym"` for non-symmetric arrays

### Expression Construction

- Explicit nested loops.
- No Einstein summation notation.
- No matrix-multiply operators.

### Symbol Naming

Standard suffixes:

- `U`: contravariant
- `D`: covariant
- `DD`, `UU`: rank-2
- `dD`, `dDD`: first or second partial derivatives
- `dupD`: upwinded derivative
- `dBarD`, `dHatD`: conformal or reference covariant derivative
- `rhs`: evolution RHS

Derivative naming:

- first partials use `*_dD`
- second partials use `*_dDD`
- FD-like names include `dD` or `dDD`
- upwinded use `dupD`
- declare derivatives with `ixp` helpers

### Expression Validation

- Every equation module validates symbolic expressions against trusted numerical values.

Pipeline:

1. build `results_dict`
2. run `ve.process_dictionary_of_expressions(...)`
3. run `ve.compare_or_generate_trusted_results(...)`

Use `nrpy.validate_expressions.validate_expressions as ve`.

Key APIs:

- `assert_equal(vardict_1, vardict_2)`
- `check_zero(expression)`
- `process_dictionary_of_expressions(dict)`
- `compare_against_trusted(...)`
- `output_trusted(...)`
- `compare_or_generate_trusted_results(...)`

Trusted file rules:

- Trusted vectors live under `*/tests/`; generated artifacts, minimal stable structure.
- No module docstrings, functions, or classes.
- Only `mpf` or `mpc` imports, with `# type: ignore`.
- Skip single-file static analysis for trusted-value files.
- No hand-editing trusted values; regenerate from owning module.
- If mismatch legitimate, delete stale file and rerun module.
- Commit message explains trusted-file change.

Equation output contract:

- Assign key outputs to `self.<expr_name>`.
- Validate with same `process_dictionary_of_expressions` / `compare_or_generate_trusted_results` pattern.
- Result names match `trusted_dict` keys.
- Existing naming like `*_rhs`, `*_expr_list_*`; no invented naming unless trusted files updated to match.

### Prohibited / Restricted Dependencies

- Core NRPy code must not depend on `numpy`.
- NRPy workflow: symbolic SymPy to generated C.
- Exception: visualization and post-processing scripts may use `numpy`.

Regex rule:

- No `import re` when `.replace()` or plain string methods enough.
- Regex only for real pattern-matching need.
- If regex justified, add comment why `.replace()` not enough.

## Infrastructure Code Rules

### Module Organization

- Usually one Python infrastructure file registers one main C function via `register_CFunction_*()`.
- Module names descriptive `snake_case`.
- `__init__.py` flat, explicit relative imports only.
- Standard order:
  1. module docstring with author info
  2. imports
  3. main `register_CFunction_*()` function
  4. helper functions only if truly useful
  5. doctest `__main__` block

### Doctests

- Registration-function docstrings end with `Doctests:` before first `>>>`.
- Older files may use `Doctest:` or `DocTests:`. New code uses `Doctests:`.

Preferred doctest pattern for stable generated C:

```python
Doctests:
>>> from nrpy.helpers.generic import validate_strings, clang_format
>>> import nrpy.c_function as cfc
>>> import nrpy.params as par
>>> cfc.CFunction_dict.clear()
>>> _ = register_CFunction_foo("Spherical")
>>> generated_str = clang_format(cfc.CFunction_dict["foo__rfm__Spherical"].full_function)
>>> _ = validate_strings(generated_str, "foo__openmp__Spherical", file_ext="c")
```

Rules:

- Import `validate_strings` / `clang_format` inside doctest.
- Clear `cfc.CFunction_dict` first.
- `file_ext="c"` usually, `"cu"` for CUDA.
- Python-based C/C++-generation doctests fit best in `nrpy/infrastructures/*/*.py`; outside, use only when meaningful behavior cheaper checks miss.
- No trusted output generation/checks for C functions dominated by large generated kernels; prefer symbolic validation or small structural checks.
- No trivial doctests whose main value is "registration succeeded".
- No doctest required when only realistic option is brittle generated-text spot check.
- Keep doctests lightweight.
- If `# FIXME` placeholder exists, finish or remove later.

### Parallel Codegen Pattern

If registration uses `nrpy.helpers.parallel_codegen`, use discovery guard:

```python
def register_CFunction_my_function() -> Union[None, pcg.NRPyEnv_type]:
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    return pcg.NRPyEnv()
```

Rules:

- Only for functions participating in parallel/discovery codegen.
- Register function after guard.
- Non-parallel registration functions return `-> None`.

### Black Suppression

- `# fmt: off` / `# fmt: on` rarely; only when Black destroys intentional alignment.

### C Function Registration from Python

Use `cfc.register_CFunction()` with named params like:

- `subdirectory`
- `includes`
- `prefunc`
- `desc`
- `cfunc_type`
- `name`
- `params`
- `include_CodeParameters_h`
- `body`
- `postfunc`

Structure:

- Declare `desc`, `cfunc_type`, `name`, `params`, `body` on separate lines.
- Keep function self-contained, readable in execution order.
- No splitting one-off setup into tiny helpers.

### Embedded C in Python Strings

- Raw strings `r"""..."""` for C bodies; `rf"""..."""` if interpolation needed.
- Double braces in f-strings: `{{` and `}}`.
- Embedded C indentation: 2 spaces.
- Keep interpolated expressions simple and Python-3.7-friendly.
- No churn of old valid indentation only for style.
- No mechanical re-indent of historical raw strings unless fixing readability or behavior.
- If `include_CodeParameters_h=True`, body uses `#include "set_CodeParameters.h"`.
- String replacement may adapt generated C when needed.

### BHaH Symbolic Codegen Rules

For new ordinary per-grid or per-point kernels:

- Prefer `BHaH.simple_loop.simple_loop()` over handwritten `LOOP_OMP(...)` or raw nested loops.
- Prefer `ccg.c_codegen(..., automatically_read_gf_data_from_memory=True)` for registered GF reads.
- Keep transforms symbolic until `c_codegen()`.
- No string-based replacement when symbolic expression can express same thing.
- No new `#include "set_CodeParameters.h"` inside generated bodies for this infrastructure class.
- Instead derive needed locals with `get_params_commondata_symbols_from_expr_list()` and `generate_definition_header()`.
- Avoid post-registration mutation of `cfc.CFunction_dict[...]`.
- Prefer explicit hooks or params in shared registration logic.
- Keep registration linear.
- No hiding one-off symbolic setup in tiny private helpers.
- Build expressions right before consume in `c_codegen()` / `simple_loop()`.
- Register gridfunctions and parity tables near point of need.
- No front-loading unrelated setup before core metadata.
- Build emitted C body top-to-bottom with `body += ...`.
- Add short comments at abrupt transitions.

### Inlining Rules

- Inline C/Python functions that do not save lines by being separate.

Separate function only if:

- called from multiple locations
- or contains more than about 10 to 15 lines real logic

Avoid:

- single-use 2 to 5 line helpers
- functions with docstring longer than body
- functions that only return one expression

### `prefunc`

- Helper C functions emitted into `prefunc`; concatenate with `prefunc += ...`.

### Standard Struct Pointer Params

Common patterns:

- `commondata_struct *restrict commondata`
- `griddata_struct *restrict griddata`
- `params_struct *restrict params`
- `bc_struct *restrict bcstruct`

### Gridfunction Naming / Grouping

Naming patterns:

- Scalars: `HHGF`, `VVGF`, `WWGF`, `TRKGF`, `CFGF`
- Symmetric rank-2: `HDD00GF`, `HDD01GF`, `HDD11GF`
- Traceless extrinsic curvature: `ADD00GF`
- Derivatives: `SRC_PARTIAL_D_HDD000GF`

Groups:

- `EVOL`
- `AUXEVOL`
- `AUX`

### Indexing Macros

Common macros:

- `IDX3(i0, i1, i2)`
- `IDX4(gf, i0, i1, i2)`
- `IDX4pt(gf, idx3)`
- `IDX2(i1, i2)`

- Custom variants like `SRC_IDX4`, `DST_IDX4`, `EX_IDX4` allowed where needed.

### Performance / Parallelism

SIMD:

- `sum_lagrange_x0_simd()` from `interpolation_lagrange_uniform.h` for vectorized inner loops.
- `#pragma omp simd` for GF-level vectorization.
- `simd_intrinsics.h` may be copied into project when SIMD disabled elsewhere.

OpenMP:

- `#pragma omp parallel for`
- `#pragma omp parallel for reduction(+:var)`
- `#pragma omp parallel for collapse(2)`
- `#pragma omp critical`
- `LOOP_OMP("omp parallel for", ...)`

### Memory / Error Handling

Memory:

- Prefer `BHAH_MALLOC` / `BHAH_FREE`.
- Plain `malloc()` / `free()` okay for simple cases.
- Check allocations for `NULL`.
- Prefer one large heap allocation over many small ones.
- VLAs okay for stack temporaries.

Errors:

- Enum-based error codes.
- Check error codes immediately.
- Common pattern: `if (commondata->error_flag != BHAHAHA_SUCCESS) return;`
- In parallel regions use `#pragma omp critical` when setting shared error flags.

### SymPy to C

- `ccg.c_codegen()` generates C from SymPy.
- Often pass expression lists and target names.
- Declare derivatives with `ixp.declarerankN()` before codegen.

### Preprocessor / Comment Patterns

Preprocessor:

- `#ifndef REAL` / `#define REAL double` / `#endif`.
- `#ifndef M_PI` fallback if needed.
- `#ifdef STANDALONE` for standalone tests.
- GCC optimization pragmas only when appropriate.
- `MAYBE_UNUSED` for intentionally unused variables.

C comments:

- `//`, not block comments, in function bodies and general C.
- Function docs come from Doxygen or `desc=`.
- Struct fields grouped with clear separator comments.

## C/H Style

### Formatting

- 2 spaces; no tabs; aim about 100-char lines; K&R braces.

### Naming

- Structs: `snake_case`, often `_struct` suffix.
- Functions/variables: `snake_case`.
- Macros/enums: `UPPER_SNAKE_CASE`.

### Header Guards / Includes

- `#ifndef` / `#define` / `#endif`; guard names `UPPER_SNAKE_CASE`.
- Follow surrounding file style if simple vs double-underscore variants differ.
- Standard includes first; project headers in quotes; conditional includes for platform-specific headers.

### Macros / Declarations

- Wrap macro args in parentheses.
- `#if defined(...)`, `#elif defined(...)`, `#endif`.
- Header helpers should be `static inline`.
- Keep return type on same line as function name.
- `const` and `restrict` where appropriate.
- Declare vars at first use.

### End-Curly-Brace Comments

- Every closing brace ending non-trivial new C block needs `// END ...` comment unless block under 5 lines and opening brace still visible.

Formats:

- Function: `} // END FUNCTION: name`
- `for`: `} // END LOOP: for <var> over <range/purpose>`
- `while`: `} // END WHILE: brief description`
- `if`: `} // END IF: brief condition`
- `else if`: `} // END ELSE IF: brief condition`
- `else`: `} // END ELSE: brief description`
- `switch`: `} // END SWITCH: brief description`
- OpenMP parallel: `} // END OMP PARALLEL: brief description`
- OpenMP critical: `} // END OMP CRITICAL: brief description`
- OpenMP `for` loop: `} // END LOOP: for <var> over <range/purpose>`
- Anonymous block: `} // END BLOCK: description`
- `do ... while`: `} while (...); // END DO-WHILE: brief description`

Rules:

- keyword after `//` all caps
- use colon in all forms
- every end-curly-brace comment includes brief description
- goal: reader see immediately what brace closes
- keep highest-signal context: step number, algorithm phase, special case, cleanup purpose
- do not replace informative comment with generic `END BLOCK` unless brace really closes anonymous scoping block
- describe what block does or handles, not just that code checked something
- for index loops, prefer semantic domain like theta points, phi points, radial points, resolutions, horizons, ghost-zone layers
- for anonymous scoped blocks tied to numbered comments, keep step number when useful
- OpenMP-parallelized loops still use `END LOOP`
- no trailing period
- no parentheses after function names

### Doxygen for C

- Every C function with body longer than about 10 lines needs Doxygen above declaration or definition.

Placement:

- If declared in header and defined in `.c`, put comment on declaration.
- Else put comment on definition.

Structure:

- opening `/**` on own line
- closing `*/` on own line
- first content line = brief description, no `@brief`
- always blank ` *` line after brief
- for complex functions prefer numbered steps
- tag order: `@param` block, `@return`, blank ` *`, then `@note` / `@warning` / `@pre`
- do not hide warning inside `@param`
- keep `@param` one line
- use `@param[in]` for `const` pointers
- use `@param[out]` or `@param[in,out]` for mutable pointers
- omit `@return` for `void`
- integer error returns should name outcomes or enum values

- Same conventions apply to `desc=` strings, except no literal `/**` / `*/`.

## Additional Project Rules

- If conditional Python logic builds C body text, prefer incremental `body += ...`.

## Quick Reference

Python: 4 spaces; Black/isort; Sphinx docstrings; triple-quoted multiline literals; old-style `typing`; no `__init__.py` docstrings; doctest runner required in runnable `nrpy/equations/**`; strongly encouraged in `infrastructures/*/*.py`; trusted `*/tests/*.py` files are generated values only, no module docstrings, no single-file static analysis.

C: 2 spaces; K&R braces; Doxygen comments; `// END ...` comments on non-trivial closing braces; header guards with `#ifndef`.

Testing/validation: run single-file static analysis on each modified Python file except trusted-value `*/tests/*.py`; trusted vectors under `tests/`; regenerate trusted outputs, do not hand-edit; favor symbolic validation over brittle golden generated-C output for large kernels.
