# Python Coding Style

> Python formatting, naming, imports, docstrings, type hints, comments, and module-shape rules. · Status: provisional · Last reconciled: 07-06-2026
> Up: [Architecture](index.md)

## Summary

Python source follows Black formatting, isort import grouping, canonical NRPy
module aliases, Sphinx/reStructuredText docstrings, explicit return type
annotations, and conservative helper-function use. `__init__.py` files and
generated trusted-value files are special cases: both omit module docstrings,
and `__init__.py` files stay as bare explicit import aggregators.

## Detail

### Formatting And Artifacts

Python uses 4-space indentation and lets Black decide line wrapping. Run Black before committing Python changes, then run
`./.github/single_file_static_analysis.sh <path-to-file.py>` for each modified
Python file. Do not add binary files, images, archives, compiled artifacts, or
other non-text assets in ordinary pull requests; redesign the change as text or
discuss a maintainer exception.

### Naming And Imports

Classes use PascalCase, functions and variables use snake_case, constants use
UPPER_SNAKE_CASE or a leading underscore, and private helpers use a leading
underscore. Boolean-like options should use positive names such as
`enable_feature`, `include_header`, `allow_resize`, or `has_ghost_zones`;
avoid double-negative flags such as `disable_feature=False`.

Imports follow isort order: standard library, third-party modules, then local
NRPy imports, with a blank line between groups. Use the canonical aliases:
`sympy as sp`, `nrpy.indexedexp as ixp`, `nrpy.params as par`, `nrpy.grid as
gri`, `nrpy.reference_metric as refmetric`, `nrpy.c_function as cfc`,
`nrpy.c_codegen as ccg`, and `nrpy.helpers.parallel_codegen as pcg`.

### `__init__.py` Files

Every `__init__.py` is a bare import aggregation file. It has no module-level
docstring, no comments, and no executable code beyond explicit relative imports
such as `from . import module` or `from .module import symbol`. Keep the
namespace flat and explicit.

### Docstrings And Module Headers

Use triple double-quoted strings for docstrings. Function docstrings use
Sphinx/reStructuredText fields: `:param name:`, `:return:`, and
`:raises ExceptionType:`. Do not use Google-style `Args:`, `Returns:`, or
`Raises:` headers. Use singular `:return:`, not `:returns:`.

Do not embed redundant type information inside `:param` descriptions because
types already live in function signatures. Describe meaning and behavior.

Every non-`__init__.py` source file has a module docstring at the top that
describes purpose and authorship. A leading filename comment may appear before
the module docstring only when it exactly matches the relative `nrpy` path form
`# <relative path from nrpy root>.py`. Test data files that only carry trusted
dictionaries omit module docstrings.

Use `Author:` for one author and `Authors:` for more than one; those are the
only allowed metadata keys. Avoid nonstandard keys such as `Email:` or
`Contributor:`, a singular `Author:` with multiple names, and mixed metadata
styles in one file. Do not add authors not already present in the file.

### Python String Literals

Multiline Python string literals should use triple double quotes whenever
practical, including `desc`, `body`, `prefunc`, `postfunc`, generated-code
snippets, long messages, and validation strings. Use raw triple-quoted strings
for embedded C or other text where backslashes must stay literal. Use
`rf"""...{name}..."""` only when interpolation is needed; keep interpolation
Python-3.7-compatible and move complex expressions to local variables first.

Do not replace a readable static multiline literal with adjacent string
fragments or `"\n".join(...)`. Single-line ordinary strings stay normal
double-quoted strings unless triple quotes materially improve readability.

### Type Hints

Use type hints broadly and include return annotations, including `-> None`.
Prefer `Dict`, `List`, `Optional`, `Tuple`, `Union`, and `cast` from `typing`;
use `typing_extensions.Literal` for constrained string values. Avoid `Any` when
a more specific type, union, protocol, `object`, or helper alias can describe
the value. If `Any` is unavoidable for a third-party typing gap, document why
inline.

Do not introduce Python 3.9 builtin generics such as `list[X]` or `dict[X, Y]`,
and do not use `X | None` union shorthand. Use `List[X]`, `Dict[X, Y]`,
`Optional[X]`, and `Union[X, Y]`. A few older directories already use newer
syntax; do not spread it to new files. `from __future__ import annotations` is
not a repo-wide standard; add it only to files that already use it or when a
specific forward-reference need justifies it.

Use `# type: ignore` selectively for imports or third-party APIs that lack
usable stubs, such as the trusted-value `mpmath` import pattern.

### Comments And Procedural Structure

Use inline `Note:` prose inside docstrings for subtle behavior. Some equation
files use reST `.. note::` blocks; either form is acceptable, but do not mix
both forms in one docstring.

Procedural code uses numbered comments in the form `# Step N:`. Substeps use
`# Step 1.a:` and module preamble steps use `# Step P1:`. The uppercase
`# STEP N:` variant appears in older examples and infrastructure code; use the
lowercase form in new code.

### Main Blocks, Doctests, And Organization

Runnable non-test, non-`__init__.py` modules under `nrpy/equations/` and its
subdirectories should start their `if __name__ == "__main__":` block with the
standard doctest runner: import `doctest` and `sys`, run `doctest.testmod()`,
print a pass/fail message, and exit with status `1` on failure. The same
pattern is strongly encouraged for runnable `nrpy/infrastructures/*/*.py`
modules and recommended elsewhere when useful.

In equation modules, symbolic validation and trusted-results generation may
follow the doctest-runner prefix inside the same main block. Imports used only
by that block belong inside it. Outside `nrpy/infrastructures/*/*.py`, doctests
that invoke Python-driven C/C++ generation are discouraged unless they provide
signal that cheaper symbolic or structural checks cannot provide.

Classes should group related functionality. Module-level functions are normal
for procedural code generation, and `register_CFunction_*()` functions are the
standard registration pattern. Raise `ValueError` for invalid inputs, use
`warnings.warn()` for non-critical issues, and perform explicit constructor
validation.

Private/module-local helpers, including leading-underscore helpers, need at
least two real call sites unless an external API, callback protocol, or test
harness requires a named function. Do not keep single-use helpers for cosmetic
string manipulation or local tidiness; inline them at the point of use.

## Sources

- [original-agents.md](../../raw/source-docs/original-agents.md) - `## Python Style`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `### Python String Literals`, `### Module Docstrings`, `### Type Hints`, `### Comments`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `### if __name__ == "__main__":`, `### General Python Organization`, `## Quick Reference`

## See Also

- Parent: [Architecture](index.md)
- Validated by: [Static Analysis](../validation/static-analysis.md)
- See also: [Contribution Style And Static Analysis](contribution-style-and-static-analysis.md)
- See also: [C And Embedded C Style](c-and-embedded-c-style.md)
