# NRPy Project Coding Style Guide

## Overview

This document provides a comprehensive style guide for the NRPy project, a numerical relativity framework that generates C code from Python/SymPy expressions. The coding standards described herein are enforced by CI tooling and must be followed by all contributors.

---

## Table of Contents

1. [Python Coding Style](#python-coding-style)
2. [Equation Setup Patterns](#equation-setup-patterns)
3. [Infrastructure Code Patterns](#infrastructure-code-patterns)
4. [C/H Coding Style](#ch-coding-style)
5. [Static Analysis Configuration](#static-analysis-configuration)
6. [Style Comparison Summary](#style-comparison-summary)

---

## Python Coding Style

### Formatting

The project uses **Black** for automatic code formatting with the following configuration:
- Line length: 88 characters
- Indentation: 4 spaces
- String quotes: Double quotes preferred

Run `black .` before committing to ensure consistent formatting.
For any modified Python file, also run `.github/single_file_static_analysis.sh <path-to-file.py>` before committing.

### Naming Conventions

| Element | Convention | Examples |
|---------|------------|----------|
| Classes | PascalCase | `CCodeGen`, `NRPyParameter`, `CodeParameter`, `ReferenceMetric`, `BSSNRHSs` |
| Functions/Methods | snake_case | `setup_FD_matrix__return_inverse_lowlevel()`, `register_CFunction_diagnostics()` |
| Variables | snake_case | `cparam_type`, `glb_params_dict`, `enable_simd` |
| Constants | UPPER_SNAKE_CASE or leading underscore | `_UNSET_DEFAULT`, `_VALID_NRPY_PARAM_TYPES` |
| Private helpers | Leading underscore | `_format_c_offset_str()`, `_flatten_and_unique_str()`, `_parse_array_spec()` |

### Import Organization

Imports follow `isort` conventions, organized in the following order:

1. **Standard library** imports
2. **Third-party** imports
3. **Local NRPy** imports

Each group should be separated by a blank line.

```python
# Standard library
import os
import sys
from typing import Any, Dict, List, Optional, Tuple, Union, cast

# Third-party
import sympy as sp
import sympy.codegen.ast as sp_ast
from typing_extensions import Literal

# Local NRPy modules
import nrpy.finite_difference as fin
import nrpy.grid as gri
import nrpy.params as par
from nrpy.helpers.generic import superfast_uniq
```

**Canonical module aliases** — the following short aliases are used 100% consistently across the codebase. Always use exactly these aliases; do not invent alternatives:

| Import | Canonical alias |
|--------|----------------|
| `import sympy` | `import sympy as sp` |
| `import nrpy.indexedexp` | `import nrpy.indexedexp as ixp` |
| `import nrpy.params` | `import nrpy.params as par` |
| `import nrpy.grid` | `import nrpy.grid as gri` |
| `import nrpy.reference_metric` | `import nrpy.reference_metric as refmetric` |
| `import nrpy.c_function` | `import nrpy.c_function as cfc` |
| `import nrpy.c_codegen` | `import nrpy.c_codegen as ccg` |
| `import nrpy.helpers.parallel_codegen` | `import nrpy.helpers.parallel_codegen as pcg` |

### `__init__.py` Files

- **Always lack module-level docstrings.** Every `__init__.py` in the codebase is a bare import aggregation file with no docstring, no comments, and no executable code beyond `from . import ...` or `from .module import ...` statements.
- Use explicit relative imports for all submodules, maintaining a flat namespace.

```python
from . import (
    ADM_Initial_Data_Reader__BSSN_Converter,
    NRPyPN_quasicircular_momenta,
    Ricci_eval,
)
```

### Docstring Style

- Use **triple double-quotes** (`"""`) consistently.
- Follow **Sphinx/reStructuredText-style** docstrings using `:param name:`, `:return:`, `:raises ExceptionType:` fields. This is **not** Google-style (which uses `Args:` / `Returns:` / `Raises:` sections) — never use the Google-style section headers.
- Use `:return:` (**singular**, not `:returns:`).
- Module docstrings at the top of each file should describe purpose and authors.
- **Exception**: `__init__.py` files never have module docstrings (see above).
- **Exception**: Test data files (trusted dict files) never have module docstrings.

**Do not embed type information inside `:param` descriptions.** The type is already in the function signature. The following is an anti-pattern seen in a few older files:
```python
# BAD — redundant, the type is already in the signature
:param enable_rfm_precompute: (bool) Whether to enable precomputation...

# GOOD — describe meaning only
:param enable_rfm_precompute: Whether to enable precomputation...
```

```python
def _format_c_offset_str(var: str, offset: int) -> str:
    """
    Format a C-style variable with an offset, matching original output.

    :param var: The variable name string (e.g., "i0").
    :param offset: The integer offset.
    :return: A string like "i0", "i0+1", or "i0-1".
    """
```

### Module Docstring Format

Every non-`__init__.py` file must have a module docstring at the top following this exact structure:

```python
"""
<One or more paragraphs describing the module's purpose.>

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
```

For multiple authors, use the plural `Authors:` key:

```python
"""
<Description.>

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Ken Sible
         ksible **at** outlook **dot* com
"""
```

Rules:
- Use singular `Author:` for exactly one author; plural `Authors:` for two or more.
- In current NRPy Python files, single-author docstrings most often put the email on the next indented line, and multi-author docstrings often list contributors as stacked name/email lines.
- The names shown in the examples above are illustrative, not prescriptive. Do not add Zachariah B. Etienne (or any other person) to a file's `Author:`/`Authors:` metadata unless that person is already an author of the file or is being intentionally credited for that file's content.
- Author contact information may be included in whatever source-level format is most practical for the file; this style guide does not enforce a specific email layout or obfuscation pattern.
- Email obfuscation is encouraged when publishing addresses in source code.

**Common anti-patterns to avoid** (all widespread in older files — do not reproduce):

- **`Author:` (singular) with multiple names** — if there is more than one contributor, the key must be `Authors:` (plural), not `Author:`.
- **`Email:` or `Contributor:` as alternative metadata keys** — non-standard. Use `Author:` / `Authors:` only.
- **Inconsistent author metadata keys within the same file** — keep the docstring metadata readable and internally consistent, but the exact email formatting is not style-enforced.

### Type Hints

- **Extensive use** of type hints throughout the codebase.
- Uses the `typing` module: `Any`, `Dict`, `List`, `Optional`, `Tuple`, `Union`, `cast`.
- Uses `typing_extensions.Literal` for constrained string values.
- Return type annotations are always present, including `-> None`.
- `# type: ignore` comments are used selectively for third-party imports lacking stubs (e.g., `from mpmath import mpf  # type: ignore`).
- **Do not use Python 3.9+ builtin generics** (`list[X]`, `dict[X, Y]`, `tuple[X, ...]`) or the `X | None` union shorthand. Always use `List[X]`, `Dict[X, Y]`, `Optional[X]`, `Union[X, Y]` from `typing` — zero occurrences of the newer syntax exist in the codebase.
- **`from __future__ import annotations`** is used in ~14 files (mostly the `BHaH/rotation` submodule). It is not the codebase-wide standard; do not add it to files that do not already use it unless there is a specific reason (e.g., forward references in the same file).

```python
def process_data(
    data: List[float],
    threshold: Optional[float] = None
) -> Dict[str, Union[int, float]]:
    ...
```

### Comment Style

- Inline comments are used sparingly, primarily for complex algorithmic explanations.
- Block comments with `#` prefix are used for section headers.
- Comments should be clear and concise.
- **`Note:` inside docstrings** — inline `Note:` prose is acceptable for clarifying subtle behavior directly inside a docstring paragraph. A small number of files (primarily in `equations/`) use the reST `.. note::` directive block instead; either form is acceptable, but do not mix both within the same docstring.

```python
# Inline Note: (common, acceptable)
:param foo: The value to process.  Note: must be positive.

# reST block (used in a few equations/ files)
.. note::
    This function assumes foo is positive.
```

**`# Step N:` convention** — Procedural code (both at module level and inside functions) uses numbered step comments to structure logic. The format is `# Step N:` (lowercase "Step", colon at end). Sub-steps use dotted notation: `# Step 1.a:`, `# Step 1.b:`. Module-level preamble steps that precede the main logic use `# Step P1:`, `# Step P2:`, etc.

```python
# Step P1: Import needed modules.
# Step 1: Set up the reference metric.
# Step 1.a: Declare tensor components.
# Step 1.b: Apply symmetry conditions.
# Step 2: Compute Christoffel symbols.
```

The uppercase variant `# STEP N:` appears in a small number of older files and is **non-standard** — use lowercase `# Step N:` in all new code.

### `if __name__ == "__main__":` Block

Every runnable non-test, non-`__init__.py` file ends with this exact block to run its embedded doctests:

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

Additional imports needed only for the `__main__` block (e.g., `os`, `nrpy.validate_expressions.validate_expressions as ve`) are placed inside the block, not at module level.

### Class and Function Organization

- Classes group related functionality (e.g., `CCodeGen`, `ReferenceMetric`).
- Module-level functions are used for procedural code generation.
- Registration pattern is common: `register_CFunction_*()` functions.
- Doctests are embedded in docstrings for testing.

### Error Handling

- `ValueError` is raised for invalid inputs.
- `warnings.warn()` is used for non-critical issues.
- Explicit validation is performed in constructors.

---

## Equation Setup Patterns

This section documents the patterns used when building symbolic equations in the `nrpy/equations/` modules.

### SymPy Usage Guidelines

- **`sp.simplify()` — MINIMAL USE**: Avoid calling `sp.simplify()` in equation-building code. It is only acceptable in test/validation code or explicit identity checks. Some modules explicitly document "No SymPy .subs() or simplify() calls" as an enforcement rule.

- **`sp.subs()`, `sp.replace()` — NEVER USE for pattern-based expression transformation**: These methods are forbidden in core equation-building. The only acceptable uses of `.subs()` are:
  - Coordinate substitution in elliptic source terms
  - Face-value substitution in GRHD fluxes
  - Systematic `nrpyABS` → `sp.Abs` conversion
  - Surface radius substitution in horizon modules

- **Preferred SymPy patterns**:
  - Use `sp.sympify(0)` / `sp.sympify(1)` for accumulator initialization
  - Use `sp.Rational()` for exact fractions (e.g., `sp.Rational(2, 3)`)
  - Use `sp.symbols()` with `real=True` for scalar quantities

### Indexed Expressions (`indexedexp` / `ixp`)

The `nrpy.indexedexp` module (imported as `ixp`) is the primary interface for tensor operations:

| Function | Purpose |
|----------|---------|
| `ixp.zerorank1()` through `ixp.zerorank4()` | Initialize tensors to zero |
| `ixp.declarerank1()` through `ixp.declarerank4()` | Declare symbolic tensors with optional symmetry |
| `ixp.symm_matrix_inverter2x2/3x3/4x4()` | Symmetric matrix inversion |
| `ixp.LeviCivitaSymbol_dim3_rank3()` / `ixp.LeviCivitaTensorUUU_dim3_rank3()` | Levi-Civita tensors |

**Symmetry specifications**: `"sym01"` (metric tensors), `"sym12"` (derivatives of symmetric tensors), `"sym01_sym23"` (Riemann-like tensors), `"nosym"` (non-symmetric arrays)

### Expression Construction Pattern

Use explicit nested loops for element-by-element tensor accumulation. No Einstein summation convention, no matrix multiplication operators:

```python
resultDD = ixp.zerorank2()
for i in range(3):
    for j in range(3):
        for k in range(3):
            resultDD[i][j] += tensor1[i][k] * tensor2[k][j]
```

### Standardized Variable Naming Conventions

| Suffix | Meaning | Examples |
|--------|---------|----------|
| `U` | Contravariant | `betaU`, `sU`, `LambdabarU` |
| `D` | Covariant | `alpha_dD`, `sD`, `SD` |
| `DD`/`UU` | Rank-2 covariant/contravariant | `gammaDD`, `AbarUU` |
| `dD`/`dDD` | First/second partial derivative | `cf_dD`, `cf_dDD` |
| `dupD` | Upwinded derivative | `trK_dupD` |
| `dBarD`/`dHatD` | Conformal/reference-metric covariant derivative | `phi_dBarD` |
| `rhs` | Evolution equation RHS | `alpha_rhs`, `gammabar_rhsDD` |

### Derivative Variable Naming Conventions

Derivative-like symbols follow strict suffix conventions to ensure codegen/validation keys stay consistent across modules:
- First partial derivatives use `*_dD` suffixes (e.g. `alpha_dD`, `cf_dD`).
- Second partial derivatives use `*_dDD` suffixes.
- Any derivative symbol intended to represent a *finite-difference-like* derivative must include the `dD`/`dDD` components (avoid derivative names that omit these markers).
- Upwinded derivatives use suffix `dupD`.

When declaring derivative objects in SymPy, use the project’s `nrpy.indexedexp` declaration helpers so the derivative symbols carry the expected suffix/name encoding.

### Expression Validation via Trusted Dictionaries

Every equation module validates its symbolic expressions against pre-computed trusted numerical values. The validation pipeline lives in `nrpy/validate_expressions/validate_expressions.py` and follows this pattern:

**Step 1 — Build a results dictionary**: Equation modules construct a dictionary mapping symbolic expression names to their SymPy expressions:

```python
results_dict = {
    "Ricci_exprs": Ricci_exprs,
    "RbarDD": RbarDD,
    "GammabarUDD": GammabarUDD,
    # ... more expressions
}
```

**Step 2 — Process and compare/generate**: The module calls `ve.compare_or_generate_trusted_results()`:

```python
import nrpy.validate_expressions.validate_expressions as ve

results_dict = ve.process_dictionary_of_expressions(
    bq.__dict__, fixed_mpfs_for_free_symbols=True
)
ve.compare_or_generate_trusted_results(
    os.path.abspath(__file__),
    os.getcwd(),
    f"{os.path.splitext(os.path.basename(__file__))[0]}_{Coord}",
    results_dict,
)
```

**How validation works internally**:

1. **CSE optimization**: Each SymPy expression is passed through SymPy's `cse()` (common subexpression elimination) to speed up numerical evaluation.
2. **Symbol substitution**: Free symbols are mapped to deterministic `mpf` values using MD5-hashed seeds (`fixed_mpfs_for_free_symbols=True` ensures reproducibility). Special constants like `PI` and `M_SQRT1_2` are set to their correct mathematical values.
3. **High-precision evaluation**: Expressions are evaluated at 30-digit precision using `mpmath`. If a result is suspiciously close to zero, precision is doubled to 60 digits to confirm whether it should be exactly zero.
4. **Trusted file comparison**: The computed `mpf`/`mpc` values are compared against a `trusted_dict` stored in a sibling `tests/` directory. Relative error tolerance is `10^(-4/5 * precision)` (~24 significant digits).
5. **Auto-generation**: If no trusted file exists, one is automatically generated and formatted with Black.

**Trusted file format**: Test data files in `tests/` directories follow the **Trusted Vector File Contract** described below:

```python
from mpmath import mpf  # type: ignore

trusted_dict = {
    "alpha_rhs": mpf("1.72309290360185808131938983598756"),
    "bet_rhsU0": mpf("0.0"),
    "GammahatUDD_0": mpf("2.74927283879025683007599345524926"),
    # ... 35+ significant digits for all values
}
```

**Key API functions**:

| Function | Purpose |
|----------|---------|
| `assert_equal(vardict_1, vardict_2)` | Assert two expression dictionaries are numerically equal |
| `check_zero(expression)` | Check if a single expression evaluates to zero |
| `process_dictionary_of_expressions(dict)` | Convert a dict of SymPy expressions to numerical `mpf`/`mpc` values |
| `compare_against_trusted(...)` | Compare results against an existing trusted file (raises `ValueError` on mismatch) |
| `output_trusted(...)` | Write a new trusted file (formatted with Black) |
| `compare_or_generate_trusted_results(...)` | Auto-selects compare or generate based on file existence |

**Updating trusted files**: If a trusted file mismatch is legitimate (e.g., after a correct algorithm change), delete the old `tests/<basename>.py` file and rerun the module to regenerate. The commit message must indicate the reason for updating the trusted file.

### Equation Module Output Contract

Equation modules must expose the symbolic outputs that will be validated and code-generated via a predictable structure. The typical pattern is:

1. During symbolic construction, assign key output expressions to `self.<expr_name>`.
2. Build a results dictionary from the instance namespace (commonly `bq.__dict__` or an explicit dict) and pass it through `nrpy.validate_expressions.validate_expressions`.

A common usage looks like:

```python
import nrpy.validate_expressions.validate_expressions as ve

results_dict = ve.process_dictionary_of_expressions(
    bq.__dict__, fixed_mpfs_for_free_symbols=True
)
ve.compare_or_generate_trusted_results(
    os.path.abspath(__file__),
    os.getcwd(),
    f"{os.path.splitext(os.path.basename(__file__))[0]}_{CoordSystem}",
    results_dict,
)
```

Naming expectations:
- Expression names in the results dictionary must match the `trusted_dict` keys used by the corresponding `tests/<module>*.py` file.
- Use `*_rhs`, `*_expr_list_*`, or the established naming already used by the module family; do not invent new conventions without updating/aligning the trusted files.

### Trusted Vector File Contract

Trusted (golden) test-vector files under `*/tests/` are treated as generated artifacts. They must follow a minimal, stable structure:

```python
from mpmath import mpf  # type: ignore

trusted_dict = {
    "some_key": mpf("123.456..."),
}
```

Style rules:
- No module-level docstrings.
- No functions/classes.
- Only `mpf`/`mpc` imports (with `# type: ignore` as used elsewhere in the repo).
- Do not hand-edit trusted values; regenerate via the owning equation/infrastructure module.

### Prohibited Dependencies

- **`import re` — USE SPARINGLY**: Regex is used only when genuine pattern matching is required (e.g., detecting coordinate-system variants by name, transforming loop-body strings). It is not acceptable for simple string manipulation that `.replace()` can handle. Five core files (`reference_metric.py`, `rfm_wrapper_functions.py`, `CarpetX/general_relativity/rhs_eval.py`, `BHaH/general_relativity/constraints_eval.py`, `ETLegacy/general_relativity/rhs_eval.py`) use `re` for legitimate reasons; follow their pattern and add a comment explaining why `.replace()` is insufficient.

- **`numpy` — NEVER DEPEND ON NUMPY**: NRPy codes cannot ever depend on numpy. The workflow is: symbolic expression → C code. All computation is done symbolically with SymPy, then code-generated to C.

---

## Infrastructure Code Patterns

Based on analysis of the BHaHAHA infrastructure module, the following patterns are standard for NRPy infrastructure code:

### Module Organization

- Each Python file typically registers **one primary C function** via a `register_CFunction_*` function
- Python modules use **snake_case** naming that directly describes their purpose
- The `__init__.py` uses explicit relative imports for all submodules, maintaining a flat namespace
- Standard module structure:
  ```python
  """Module docstring with author info."""
  # Imports
  # Main registration function: register_CFunction_<name>()
  # Optional helper functions (Python-side symbolic computation)
  # if __name__ == "__main__": doctest block
  ```

### Doctest Conventions

#### `Doctests:` section label

Infrastructure registration functions consistently use a `Doctests:` (plural) label as the last item in the docstring, immediately before the first `>>>` line:

```python
def register_CFunction_foo(...) -> ...:
    """
    Brief description.

    :param CoordSystem: ...
    :return: ...

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> ...
    """
```

The singular `Doctest:` label also appears in a small number of older files; prefer the plural `Doctests:` in new code.

#### `validate_strings` pattern

The standard doctest idiom for verifying generated C code, when the emitted C text is stable enough for exact output comparison, is:

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

Key points:
- `clang_format` normalizes whitespace before comparison; always apply it to `full_function` before passing to `validate_strings`.
- `validate_strings` compares against a trusted file in `tests/` (auto-generated on first run).
- Set `file_ext="cu"` when the generated code is CUDA, `"c"` otherwise.
- Import `validate_strings` (and `clang_format` if needed) inside the doctest, not at module level.
- **Exception — generated-kernel-dominated C functions**: Do **not** generate or check trusted output files for C functions whose bodies primarily consist of generated kernels, especially large kernels emitted from SymPy expressions. Such output is too sensitive to SymPy version and codegen details for exact string comparison to be a reliable unit-test signal.
- For these generated-kernel-heavy functions, prefer validation at the symbolic-expression level or with cheaper structural/sanity checks instead of outputting a golden C file under `tests/`.

#### Doctest placeholders

Some functions include `Doctests:` blocks with `# FIXME` placeholders. Treat these as temporary scaffolding:
- If you add doctests, keep them lightweight (no heavy code generation, no long numeric computations).
- Prefer doctests that validate generated strings (e.g., via `validate_strings`) or cheap invariants, but skip golden-output doctests for large generated kernels as described above.
- If you see `# FIXME`, either complete it in a follow-up PR or remove the placeholder when the doctest is ready.

**Doctest prohibition:** Do not write doctests whose primary purpose is to assert that a `CFunction` successfully registers (e.g., checking membership in `cfc.CFunction_dict` or that `register_CFunction_*()` returns without error). Such doctests provide low value and tend to be brittle.

### Parallel Codegen Phase Detection

Registration functions that support parallel code generation use a standard early-return pattern:

```python
from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.helpers.parallel_codegen as pcg

def register_CFunction_my_function() -> Union[None, pcg.NRPyEnv_type]:
    """Register my C function."""
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    # ... actual registration logic
    return pcg.NRPyEnv()
```

This pattern allows the codegen system to first collect all function calls during a discovery phase, then replay them in parallel workers. Functions that do not need parallel codegen omit this pattern and return `-> None` directly.

**Rule (when using `nrpy.helpers.parallel_codegen`)**: Any C-registration function that participates in parallel/discovery codegen must implement the early-return guard shown above, and only register the C function after this guard.

### Black Format Suppression

Use `# fmt: off` and `# fmt: on` sparingly, only when Black formatting would destroy intentional alignment. The most common case is `par.CodeParameter` registration lines where positional arguments must align:

```python
# fmt: off
_ = par.CodeParameter("int", __name__, "nn_0", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("int", __name__, "nn", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("REAL", __name__, "CFL_FACTOR", 0.5, commondata=True)
# fmt: on
```

### C Function Registration Pattern

All C functions are registered using `cfc.register_CFunction()` with these standard parameters:
```python
cfc.register_CFunction(
    subdirectory="",           # Empty string for root output
    includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
    prefunc=prefunc,           # Optional: helper functions emitted before main function
    desc=desc,                 # Multi-line description
    cfunc_type=cfunc_type,     # Return type: "void", "int", "REAL", "const char *"
    name=name,                 # Function name (snake_case)
    params=params,             # Parameter string
    include_CodeParameters_h=False,
    body=body,                 # Main function body as raw string
    postfunc=postfunc,         # Optional: standalone test code after main function
)
```

#### Embedded C Code String Conventions

C code embedded in Python registration functions follows these conventions:

- **Raw strings** (`r"""..."""`) are used for C code bodies to avoid escape character issues.
- **f-strings with doubled braces** (`rf"""...{var}..."""`) are used when Python variable interpolation is needed. Braces in C code (e.g., struct initializers, loop bodies) must be doubled: `{{`, `}}`.
- **C code indentation** inside raw strings uses 2 spaces (matching the C style guide), regardless of the Python indentation level.

For readability, embedded C inside raw strings should use 2-space indentation when practical, but do not churn indentation across modules that are already valid.
- **`#include` directives** inside C bodies use `#include "set_CodeParameters.h"` when `include_CodeParameters_h=True`.
- **String replacement** is used to adapt generated C code: `.replace("auxevol_gfs[IDX4(", "commondata->interp_src_gfs[IDX4(SRC_")`.

**Raw-string indentation note**: Embedded C inside raw strings may appear with minor indentation variance across modules due to historical edits. The requirement is that the C is valid and consistently readable; do not mechanically re-indent unless you are fixing a functional or formatting issue.

```python
body = rf"""
  const int Nxx0 = params->Nxx0;
  for (int i0 = NGHOSTS; i0 < Nxx0 + NGHOSTS; i0++) {{
    {ccg.c_codegen(expr_list, lhs_list)}
  }} // END LOOP over i0
"""
```

#### Python Function Structure for C Function Registration

The Python registration function must closely imitate the C function it registers. Follow these guidelines:

- Use separate lines for declaring `desc`, `cfunc_type`, `name`, `params`, `body`, etc. variables before passing them to `register_CFunction()`
- Avoid helper functions when possible - the registration function should be self-contained

Example pattern:
```python
def register_CFunction_my_function() -> None:
    """Register my C function."""
    desc = """Description of what the function does."""
    cfunc_type = "void"
    name = "bah_my_function"
    params = "const commondata_struct *restrict commondata, ..."
    body = r"""
    // C code body
    """
    cfc.register_CFunction(
        subdirectory="",
        includes=["BHaH_defines.h"],
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )
```

### Function Inlining Guidelines

**Inline all C and Python functions that do not save lines of code by existing separately.** The worst offenders are functions that:
- Are used only once
- Are only a few lines long
- When considering Doxygen-style docstrings (C) or proper docstrings (Python), would need to have several lines and/or be called many times for separation to make sense

**Rule of thumb**: A function should only be separated if it:
- Is called from multiple locations, OR
- Is complex enough (>10-15 lines of actual logic) that separation improves readability

**Anti-patterns to avoid**:
- Single-use helper functions that are only 2-5 lines
- Functions where the docstring is longer than the function body
- Functions that simply return a single expression or trivial computation
- Helper functions that only exist to give a name to a small code fragment

**Example of what NOT to do** (Python):
```python
# BAD: This function is only 3 lines + docstring overhead
def compute_face_value(val_left: REAL, val_right: REAL) -> REAL:
    """Compute face value as average."""
    return (val_left + val_right) / 2
```

**Example of what NOT to do** (C):
```c
// BAD: Single-use, trivial function with Doxygen overhead
/**
 * @brief Compute face value.
 * @param left Left value.
 * @param right Right value.
 * @return Average of left and right.
 */
static inline REAL compute_face_value(const REAL left, const REAL right) {
    return (left + right) / 2;
}
```

**Preferred approach**: Inline the computation directly where it's used:
```python
# GOOD: Direct, self-explanatory
face_value = (val_left + val_right) / 2
```

```c
// GOOD: Inline trivial computation
const REAL face_value = (val_left + val_right) / 2;
```

### Helper Function Pattern (prefunc)

Helper C functions are generated as strings and concatenated into `prefunc`:
```python
prefunc = helper_function_1()
prefunc += helper_function_2()
```

### Standard Struct Parameters

C functions consistently use these struct pointer patterns:
- `commondata_struct *restrict commondata` - Global simulation parameters and shared data
- `griddata_struct *restrict griddata` - Per-grid data including gridfunctions and coordinates
- `params_struct *restrict params` - Grid-specific parameters (Nxx, dxx, etc.)
- `bc_struct *restrict bcstruct` - Boundary condition data

### Grid Function (GF) Naming Conventions

Grid functions follow a strict naming convention encoding tensor type and indices:
- **Scalars**: `HHGF`, `VVGF`, `WWGF`, `TRKGF`, `CFGF`
- **Rank-2 symmetric tensors**: `HDD00GF`, `HDD01GF`, `HDD11GF`, etc. (DD = covariant)
- **Traceless extrinsic curvature**: `ADD00GF`, `ADD01GF`, etc.
- **Derivatives**: `SRC_PARTIAL_D_HDD000GF`, `SRC_PARTIAL_D_WW0GF`

### GF Group Organization

- **EVOL group**: Time-evolved quantities (`hh`, `vv`)
- **AUXEVOL group**: Auxiliary quantities needed for evolution (`hDD`, `WW`, `partial_D_hDD`, `trK`, `aDD`)
- **AUX group**: Diagnostic quantities (`Theta`)

### Indexing Macros

- `IDX3(i0, i1, i2)` - 3D flat indexing
- `IDX4(gf, i0, i1, i2)` - 4D (GF + 3D) indexing
- `IDX4pt(gf, idx3)` - 4D with precomputed 3D index
- `IDX2(i1, i2)` - 2D indexing for surface data
- Custom macros defined per function: `SRC_IDX4`, `DST_IDX4`, `EX_IDX4`

### SIMD/Vectorization Patterns

- Uses `sum_lagrange_x0_simd()` from `interpolation_lagrange_uniform.h` for vectorized inner loops
- `#pragma omp simd` used for GF-level vectorization
- `simd_intrinsics.h` is copied to project when SIMD is disabled elsewhere

### OpenMP Patterns

- `#pragma omp parallel for` for outer loops
- `#pragma omp parallel for reduction(+:var)` for accumulations
- `#pragma omp parallel for collapse(2)` for nested loops
- `#pragma omp critical` for shared state updates
- `LOOP_OMP("omp parallel for", i0, start, end, ...)` macro for multi-dimensional loops

### Memory Management

- Uses `BHAH_MALLOC(ptr, size)` macro for tracked allocations
- Uses `BHAH_FREE(ptr)` macro for tracked deallocations
- Standard `malloc()`/`free()` for simple cases
- All allocations are checked for NULL with error returns
- Single large heap allocations preferred over many small ones
- VLAs (Variable Length Arrays) used for stack-allocated temporary buffers

### Error Handling

- Enum-based error codes: `bhahaha_error_codes`
- Pattern: `BHAHAHA_SUCCESS` (0), `MODULE_SPECIFIC_ERROR_NAME`
- Error messages provided via `bah_error_message()` function
- Functions return error codes (int) that are checked immediately
- Pattern: `if (commondata->error_flag != BHAHAHA_SUCCESS) return;`
- `#pragma omp critical` used to set error flags in parallel regions

### SymPy Integration for Code Generation

- Uses `ccg.c_codegen()` to convert SymPy expressions to C code
- Pattern: `ccg.c_codegen([expr1, expr2], ["var1", "var2"], enable_fd_codegen=True)`
- Tensor derivatives declared via `ixp.declarerankN()` then passed to codegen

### Preprocessor Directive Patterns

- `#ifndef REAL / #define REAL double / #endif` for type flexibility
- `#ifndef M_PI / #define M_PI ... / #endif` for portability
- `#ifdef STANDALONE` blocks for self-contained test programs
- `#pragma GCC optimize("unroll-loops")` before performance-critical functions
- `#pragma GCC reset_options` after function body
- `MAYBE_UNUSED` attribute for unused variables in certain code paths

### C Code Comment Patterns

- Block comments with `//` style, not `/* */`
- End-of-loop comments: `} // END LOOP over theta`
- End-of-block comments: `} // END IF: condition description`
- Function documentation in `desc=` parameter of `register_CFunction()`

### Struct Field Organization

Fields are grouped by purpose with clear comment separators:
```c
//==========================
// Metric and grid setup
//==========================
//==========================
// External Input Numerical Grid: Radial parameters
//==========================
```

---

## C/H Coding Style

### 1. Indentation

- Use **2 spaces** consistently.
- **No tabs** are permitted.

```c
typedef struct {
  //==========================
  // Metric and grid setup
  //==========================
  REAL *input_metric_data;
  int num_points;
} bhahaha_params_and_data_struct;
```

### 2. Naming Conventions

| Element | Convention | Examples |
|---------|------------|----------|
| Structs | snake_case with `_struct` suffix | `bhahaha_params_and_data_struct`, `diags_integration_recipe_t` |
| Functions | snake_case | `diag_write_header()`, `bah_radial_grid_cell_centered_set_up()` |
| Macros | UPPER_SNAKE_CASE | `IDX2`, `MAX_RESOLUTIONS`, `DIAG_TIME_COMMENT_FMT` |
| Variables | snake_case | `num_spaces_in_coord_names`, `interp_src_gfs` |
| Enums | UPPER_SNAKE_CASE | `INTERP_GAMMADDXXGF`, `DIAGS_EXTRACT_VOLUME` |

### 3. Brace Placement

- **K&R style**: Opening brace on the same line for functions and control structures.
- Struct definitions use opening brace on the same line.

```c
static inline void diag_write_time_comment(FILE *file_ptr, const REAL time) {
  fprintf(file_ptr, DIAG_TIME_COMMENT_FMT, time);
}
```

### 4. Header Guard Style

- Use **traditional ifndef/define pattern** with UPPER_SNAKE_CASE.
- Two variants are observed:
  - Simple: `#ifndef BHAHAHA_HEADER_H` / `#define BHAHAHA_HEADER_H`
  - Double-underscore: `#ifndef __SIMD_INTRINSICS_H__` / `#define __SIMD_INTRINSICS_H__`

```c
#ifndef BHAHAHA_HEADER_H
#define BHAHAHA_HEADER_H

// Header content

#endif // BHAHAHA_HEADER_H
```

### 5. Include Organization

- Standard library includes first (`<math.h>`, `<string.h>`, `<stdio.h>`).
- Project headers in quotes: `"BHaH_defines.h"`, `"diagnostic_gfs.h"`.
- Conditional includes for platform-specific code:

```c
#if defined(__AVX512F__)
#include <immintrin.h>
#endif
```

### 6. Comment Style

- **Doxygen-style** `/** ... */` for function documentation.
- `//` for inline comments.
- Use `@param`, `@return`, `@brief`, `@details`, `@note`, `@code` tags.
- Section headers with `// ========================` patterns.

```c
/**
 * @file diagnostics_volume_integration_helpers.h
 * @brief Provides a recipe-based interface to compute domain integrals...
 *
 * @section usage Usage
 * ...
 */
```

### 7. Macro Conventions

- Function-like macros use parenthesized arguments.
- All macro parameters are wrapped in parentheses for safety.
- Conditional compilation uses `#if defined(...)`, `#elif defined(...)`, `#endif`.

```c
#define IDX2(itheta, iphi) ((itheta) + NUM_THETA * (iphi))
```

### 8. Function Declaration Style

- `static inline` for header-defined helper functions.
- Return type on the same line as function name.
- Parameters aligned, with type qualifiers (`const`, `restrict`).

```c
static inline void diag_write_header(FILE *file_ptr, const char *coord_names, const int NUM_GFS,
                                     const int which_gfs[], const char **diagnostic_gf_names) {
  // Function body
}
```

### 9. Variable Declaration Patterns

- Variables declared at point of first use (C99 style).
- `const` used extensively for immutable parameters.
- `restrict` pointer qualifier for performance.
- Grouped declarations with aligned comments.

---

## Static Analysis Configuration

The `.github/single_file_static_analysis.sh` script enforces the following checks:

Run this script on every modified Python file before committing. This is the required pre-commit check for Python changes:

```bash
./.github/single_file_static_analysis.sh path/to/modified_file.py
```

| Tool | Purpose | Configuration |
|------|---------|---------------|
| **black** | Code formatting | `--check` mode |
| **isort** | Import sorting | `--check-only` |
| **mypy** | Type checking | `--strict --allow-untyped-calls` |
| **pylint** | Code quality | `.pylintrc`, threshold ≥ 9.91/10 |
| **pydocstyle** | Docstring style | `.pydocstyle` config |
| **darglint** | Docstring argument checking | `-v 2` (verbose) |
| **doctests** | Embedded tests | `python3 <file>` |

This indicates the project enforces **very strict** coding standards with a near-perfect pylint score requirement.

---

## Style Comparison Summary

| Aspect | Python | C |
|--------|--------|---|
| Indentation | 4 spaces | 2 spaces |
| Line length | 88 characters (black) | ~100 characters |
| Naming | snake_case for functions, PascalCase for classes | snake_case for functions/variables, UPPER_CASE for macros |
| Docstrings | Sphinx/reStructuredText-style with `:param:` / `:return:` | Doxygen-style with `@param` |
| Type hints | Extensive | N/A (C language) |
| Braces | N/A | K&R style |
| Header guards | N/A | `#ifndef FILE_H` / `#define FILE_H` |
| `__init__.py` | No docstrings, bare imports only | N/A |
| Test data files | No docstrings, `trusted_dict` only | N/A |

---

## Additional Notes

- All code contributions must pass the static analysis checks before being merged.
- For Python changes, run `.github/single_file_static_analysis.sh` on each modified Python file, not just on a hand-picked subset.
- When in doubt, follow the existing patterns in the codebase.
- This style guide is a living document and may be updated as the project evolves.
- **Author email formatting**: Source files may use any readable email formatting or obfuscation scheme. Obfuscation is encouraged, but this guide does not enforce one exact representation.
- **Doctest placeholders**: Some registration functions have `Doctests:\n    # FIXME` in their docstrings, indicating tests that need to be written.
- **`body +=` for conditional C code**: When a C function body has sections that are conditionally included based on Python parameters, build the body string incrementally with `+=`:
  ```python
  body = ""
  if enable_feature:
      body += r"""  // optional section"""
  body += r"""  // always-present section"""
  ```
  This is the accepted pattern when a single `r"""..."""` literal cannot capture all cases.
