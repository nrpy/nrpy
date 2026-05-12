# NRPy Agent Instructions

Use [coding_style.md](coding_style.md) as ground truth. If `AGENTS.md` and `coding_style.md` differ, `coding_style.md` wins. This file compress repo rules for agents. CI enforce many rules.

## Scope

- NRPy generate C from Python/SymPy.
- No binary files in PRs unless maintainers approve exception.
- No images, archives, compiled artifacts, or other non-text assets unless maintainers approve exception.

## Required Checks

- Run `black .` before commit for Python changes.
- Run `./.github/single_file_static_analysis.sh <path-to-file.py>` on every modified Python file before commit.
- Exception: trusted-value files in `.../tests/*.py` skip single-file static analysis. They store generated reference values only.
- Do not analyze hand-picked subset only.
- Contributions must pass project static analysis before merge.

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

- Follow Black.
- Use 4 spaces.

### Naming

- Classes: `PascalCase`
- Constants: `UPPER_SNAKE_CASE` or leading `_`
- Private helpers: leading `_`
- Boolean-like params use positive names.

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
- Separate groups with blank line.
- Follow `isort`.

Canonical aliases mandatory:

- `import sympy as sp`
- `import nrpy.indexedexp as ixp`
- `import nrpy.params as par`
- `import nrpy.grid as gri`
- `import nrpy.reference_metric as refmetric`
- `import nrpy.c_function as cfc`
- `import nrpy.c_codegen as ccg`
- `import nrpy.helpers.parallel_codegen as pcg`

- Do not invent aliases.

### `__init__.py`

- No module docstring.
- No comments.
- No executable code beyond explicit relative import aggregation.
- Keep namespace flat with explicit relative imports.

### Python Docstrings

- Use triple double quotes.
- Use Sphinx/reST style:
  - `:param name:`
  - `:return:`
  - `:raises ExceptionType:`
- Use `:return:`, not `:returns:`.
- Never use Google-style `Args:` / `Returns:` / `Raises:`.
- Do not put type info in `:param` text.

### Python String Literals

- Prefer triple double-quoted multiline strings for multiline literals.
- Applies to `desc`, `body`, `prefunc`, `postfunc`, generated-code snippets, long messages, validation strings.
- Use `r"""..."""` when backslashes must stay literal.
- Use `rf"""..."""` only when interpolation needed.
- Keep interpolated expressions simple and Python-3.7-friendly. Move complex expressions to locals first.
- Do not replace clear static multiline literals with adjacent string fragments or `"\n".join(...)` unless readability clearly improves.
- Keep single-line ordinary strings as normal double-quoted strings unless triple quotes help.

### Module Docstrings

For every non-`__init__.py` Python file, top docstring format:

```python
"""
<Description paragraphs.>

Author: Name
        email obfuscated
"""
```

For multiple authors:

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

- Use `Author:` for one author, `Authors:` for multiple.
- Do not add authors not already authors of file.
- Optional leading filename comment only if exact form `# <relative path from nrpy root>.py`
- `__init__.py` never get module docstrings.
- Trusted test-vector files never get module docstrings.
- Use only `Author:` / `Authors:` keys.

Avoid:

- singular `Author:` with multiple names
- `Email:` or `Contributor:`
- mixed metadata styles in one file

### Type Hints

- Use type hints heavily.
- Always include return annotations, including `-> None`.
- Use `typing` forms: `Dict`, `List`, `Optional`, `Tuple`, `Union`, `cast`.
- Use `typing_extensions.Literal` for constrained strings.
- Use `# type: ignore` only when needed, mainly third-party typing gaps.
- Avoid `Any` if more precise type possible.
- If `Any` unavoidable, justify inline.
- Do not use builtin generics like `list[X]`, `dict[X, Y]`, `tuple[X, ...]`.
- Do not use `X | None`.
- Do not add `from __future__ import annotations` unless real need.

### Comments

In docstrings:

- Inline `Note:` okay.
- reST `.. note::` okay.
- Do not mix both styles in one docstring.

Procedural comment style:

- `# Step N:`
- substeps: `# Step 1.a:`
- module preamble: `# Step P1:`
- use `Step`, not `STEP`

### `if __name__ == "__main__":`

- Runnable non-test, non-`__init__.py` modules under `nrpy/equations/` and subdirs should start `__main__` block with exact doctest-runner below.
- Runnable modules under `infrastructures/*/*.py`: strongly encouraged.
- Other runnable modules: recommended when useful.
- Outside `nrpy/infrastructures/*/*.py`, Python-based C/C++ codegen doctests discouraged unless they add signal cheaper checks cannot provide.

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

- Imports used only here belong inside block.
- In `nrpy/equations/**`, symbolic validation and trusted-results generation may follow this prefix in same `__main__` block.
- Allowed follow-on steps there include `ve.process_dictionary_of_expressions(...)`, `ve.compare_or_generate_trusted_results(...)`, loops over coordinate systems, loops over validation keys.
- This allowance is for established equation-module validation workflow, not arbitrary script logic.

### General Python Organization

- Classes group related functionality.
- Module functions handle procedural codegen.
- `register_CFunction_*()` pattern common.
- Doctests live in docstrings.
- Raise `ValueError` for invalid input.
- Use `warnings.warn()` for non-critical issues.
- Validate explicitly in constructors.

## Equation Setup Rules

### SymPy

- Avoid `sp.simplify()` in equation-building code.
- Allowed only in tests, validation, or explicit identity checks.
- `nrpy.helpers.cached_functions.cached_simplify()` is narrow exception.
- Use cached simplification sparingly for small local scalar subexpressions when it materially helps symbolic tractability, codegen stability, or readability, especially in `nrpy/equations/general_relativity`.
- Do not use it on whole tensors or large assembled RHS expressions when plain symbolic construction is practical.
- If need not obvious nearby, add short comment explaining why cached simplification warranted.
- Do not use `sp.subs()` / `sp.replace()` for pattern-based transformation.

Allowed `.subs()` cases:

- coordinate substitution in elliptic source terms
- face-value substitution in GRHD fluxes
- systematic `nrpyABS` to `sp.Abs`
- surface radius substitution in horizon modules
- evaluation at specific parameter value like `.subs(t, t_attach)`

Preferred:

- init accumulators with `sp.sympify(0)` / `sp.sympify(1)`
- use `sp.Rational()` for exact fractions
- use `sp.symbols(..., real=True)` for scalars

### Indexed Expressions

- Use `nrpy.indexedexp` as `ixp` for tensors.

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

- Use explicit nested loops.
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

- Every equation module validate symbolic expressions against trusted numerical values.

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

- Trusted vectors live under `*/tests/`.
- Treat trusted vectors as generated artifacts with minimal stable structure.
- No module docstrings.
- No functions or classes.
- Only `mpf` or `mpc` imports, with `# type: ignore`.
- Skip single-file static analysis for these trusted-value files.
- Do not hand-edit trusted values.
- Regenerate from owning module.
- If mismatch legitimate, delete stale file and rerun module.
- Commit message must explain trusted-file change.

Equation output contract:

- Assign key outputs to `self.<expr_name>`.
- Validate with same `process_dictionary_of_expressions` / `compare_or_generate_trusted_results` pattern.
- Result names must match `trusted_dict` keys.
- Use existing naming like `*_rhs`, `*_expr_list_*`.
- Do not invent naming unless trusted files updated to match.

### Prohibited / Restricted Dependencies

- Core NRPy code must not depend on `numpy`.
- NRPy workflow: symbolic SymPy to generated C.
- Exception: visualization and post-processing scripts may use `numpy`.

Regex rule:

- Do not `import re` when `.replace()` or plain string methods enough.
- Use regex only for real pattern-matching need.
- If regex justified, add comment why `.replace()` not enough.

## Infrastructure Code Rules

### Module Organization

- Usually one Python infrastructure file register one main C function via `register_CFunction_*()`.
- Module names descriptive `snake_case`.
- `__init__.py` flat, explicit relative imports only.
- Standard order:
  1. module docstring with author info
  2. imports
  3. main `register_CFunction_*()` function
  4. helper functions only if truly useful
  5. doctest `__main__` block

### Doctests

- Registration-function docstrings should end with `Doctests:` before first `>>>`.
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
- Use `file_ext="c"` usually, `"cu"` for CUDA.
- Python-based C/C++-generation doctests are most appropriate in `nrpy/infrastructures/*/*.py`.
- Outside `nrpy/infrastructures/*/*.py`, use only when they verify meaningful behavior cheaper checks would miss.
- Do not generate or check trusted output files for C functions dominated by large generated kernels. Prefer symbolic validation or small structural checks there.
- Forbid trivial doctests that only check box.
- Do not write doctests whose main value is "registration succeeded".
- Do not require a doctest when only realistic option is brittle generated-text spot check.
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

- Use only for functions participating in parallel/discovery codegen.
- Register function after guard.
- Non-parallel registration functions return `-> None`.

### Black Suppression

- Use `# fmt: off` / `# fmt: on` rarely.
- Only when Black destroys intentional alignment.

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
- Do not split one-off setup into tiny helpers.

### Embedded C in Python Strings

- Use raw strings `r"""..."""` for C bodies.
- Use `rf"""..."""` if interpolation needed.
- Double braces in f-strings: `{{` and `}}`.
- Embedded C indentation uses 2 spaces.
- Keep interpolated expressions simple and Python-3.7-friendly.
- Do not churn old valid indentation only for style.
- Do not mechanically re-indent historical raw strings unless fixing readability or behavior.
- If `include_CodeParameters_h=True`, body uses `#include "set_CodeParameters.h"`.
- String replacement may adapt generated C when needed.

### BHaH Symbolic Codegen Rules

For new ordinary per-grid or per-point kernels:

- Prefer `BHaH.simple_loop.simple_loop()` over handwritten `LOOP_OMP(...)` or raw nested loops.
- Prefer `ccg.c_codegen(..., automatically_read_gf_data_from_memory=True)` for registered GF reads.
- Keep transforms symbolic until `c_codegen()`.
- Do not use string-based replacement when symbolic expression can express same thing.
- Do not add new `#include "set_CodeParameters.h"` inside generated bodies for this infrastructure class.
- Instead derive needed locals with `get_params_commondata_symbols_from_expr_list()` and `generate_definition_header()`.
- Avoid post-registration mutation of `cfc.CFunction_dict[...]`.
- Prefer explicit hooks or params in shared registration logic.
- Keep registration linear.
- Do not hide one-off symbolic setup in tiny private helpers.
- Build expressions right before consume in `c_codegen()` / `simple_loop()`.
- Register gridfunctions and parity tables near point of need.
- Do not front-load unrelated setup before core metadata.
- Build emitted C body top-to-bottom with `body += ...`.
- Add short comments at abrupt transitions.

### Inlining Rules

- Inline C and Python functions that do not save lines by being separate.

Separate function only if:

- called from multiple locations
- or contains more than about 10 to 15 lines real logic

Avoid:

- single-use 2 to 5 line helpers
- functions with docstring longer than body
- functions that only return one expression

### `prefunc`

- Helper C functions emitted into `prefunc`.
- Concatenate with `prefunc += ...`.

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

- Use `sum_lagrange_x0_simd()` from `interpolation_lagrange_uniform.h` for vectorized inner loops.
- Use `#pragma omp simd` for GF-level vectorization.
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

- Use enum-based error codes.
- Check error codes immediately.
- Common pattern: `if (commondata->error_flag != BHAHAHA_SUCCESS) return;`
- In parallel regions use `#pragma omp critical` when setting shared error flags.

### SymPy to C

- Use `ccg.c_codegen()` to generate C from SymPy.
- Often pass expression lists and target names.
- Declare derivatives with `ixp.declarerankN()` before codegen.

### Preprocessor / Comment Patterns

Preprocessor:

- Use `#ifndef REAL` / `#define REAL double` / `#endif`.
- Use `#ifndef M_PI` fallback if needed.
- Use `#ifdef STANDALONE` for standalone tests.
- Use GCC optimization pragmas only when appropriate.
- Use `MAYBE_UNUSED` for intentionally unused variables.

C comments:

- Use `//`, not block comments, in function bodies and general C.
- Function docs come from Doxygen or `desc=`.
- Struct fields grouped with clear separator comments.

## C/H Style

### Formatting

- Use 2 spaces.
- No tabs.
- Aim about 100-char lines.
- Use K&R braces.

### Naming

- Structs: `snake_case`, often `_struct` suffix
- Functions: `snake_case`
- Macros: `UPPER_SNAKE_CASE`
- Variables: `snake_case`
- Enums: `UPPER_SNAKE_CASE`

### Header Guards / Includes

- Use `#ifndef` / `#define` / `#endif`.
- Guard names `UPPER_SNAKE_CASE`.
- Follow surrounding file style if simple vs double-underscore variants differ.
- Standard includes first.
- Project headers in quotes.
- Conditional includes for platform-specific headers.

### Macros / Declarations

- Wrap macro args in parentheses.
- Use `#if defined(...)`, `#elif defined(...)`, `#endif`.
- Header helpers should be `static inline`.
- Keep return type on same line as function name.
- Use `const` and `restrict` where appropriate.
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

Python:

- 4 spaces
- Black + isort
- Sphinx docstrings
- triple-quoted multiline strings for multiline literals
- old-style `typing` annotations
- no docstrings in `__init__.py`
- doctest runner block required at end of runnable files in `nrpy/equations/**`
- doctest runner block strongly encouraged in `infrastructures/*/*.py` and recommended elsewhere
- trusted-value files in `*/tests/*.py` store reference values only; no module docstrings, no single-file static analysis

C:

- 2 spaces
- K&R braces
- Doxygen comments
- `// END ...` comments on non-trivial closing braces
- header guards with `#ifndef`

Testing / validation:

- single-file static analysis on each modified Python file except trusted-value `*/tests/*.py`
- trusted vectors under `tests/`
- regenerate trusted outputs, do not hand-edit
- favor symbolic validation over brittle golden generated-C output for large kernels
