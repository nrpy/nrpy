# Python Coding Style

> Python formatting, naming, canonical imports, Sphinx docstrings, Python-3.7 typing, `__main__` doctest runners, and prohibited `re` imports in `nrpy/`. · Status: provisional
> Up: [Architecture](index.md)

## Summary

Handwritten Python in NRPy uses 4-space Black formatting, `isort` import grouping with canonical aliases, bare explicit `__init__.py` aggregators, Sphinx `:param:`/`:return:`/`:raises:` docstrings (omitting `:return:` on `-> None`), Python-3.7 `typing` generics (`List`, `Dict`, `Optional`, `Union`), `# Step N:` procedural comments, and `doctest.testmod()` `__main__` blocks, while forbidding `import re` under `nrpy/`.

## Detail

### Formatting, Naming, Imports, And `__init__.py`

Black controls 4-space wrapping ([Static Analysis](../validation/static-analysis.md)), and non-text binaries/images/archives are excluded from ordinary commits ([Generated Output Boundaries](generated-output-boundaries.md)). Use `PascalCase` classes, `snake_case` functions/variables, `UPPER_SNAKE_CASE` or `_`-prefixed constants, `_`-prefixed private helpers (requiring `>= 2` real call sites unless required by an external callback; inline single-use helpers), and positive boolean names (`enable_*`, `include_*`, `has_*`).

Group imports with `isort` (stdlib, third-party, local `nrpy` separated by blank lines) using canonical aliases: `sympy as sp`, `nrpy.indexedexp as ixp`, `nrpy.params as par`, `nrpy.grid as gri`, `nrpy.reference_metric as refmetric`, `nrpy.c_function as cfc`, `nrpy.c_codegen as ccg`, and `nrpy.helpers.parallel_codegen as pcg`. Keep `__init__.py` files free of docstrings, comments, and logic beyond explicit `from . import module` / `from .module import sym` (`BHaH/fisheye/__init__.py` `__all__` is a legacy exception).

### Docstrings, String Literals, And Python-3.7 Type Hints

Every non-`__init__.py`, non-`trusted_dict` Python module starts with optional `# <relative/path/from/nrpy>.py` and a triple-double-quoted module docstring carrying `Author:` (one name) or `Authors:` (multiple names; no `Email:` or `Contributor:`). Function docstrings use Sphinx `:param name:` (describe semantics without repeating signature types), singular `:return:`, and `:raises Exc:`—never Google `Args:`/`Returns:`. Functions returning `-> None` omit `:return:` completely (no `:return None:` or `# noqa: DAR202`), whereas `Union[None, X]` / `Optional[X]` document both cases. Do not mix inline `Note:` and `.. note::` in one docstring, and use lowercase `# Step N:`, `# Step 1.a:`, and `# Step P1:` in procedural code.

Prefer `"""..."""` for multiline text (`r"""..."""` for literal backslashes, `rf"""...{var}..."""` with Python-3.7-safe locals when interpolating; avoid `"...\n".join(...)` for static blocks). Annotate all parameters and return types (`-> None` included) using `typing` (`Dict`, `List`, `Optional`, `Tuple`, `Union`, `cast`, `typing_extensions.Literal`); never introduce Python 3.9+ `list[X]` or 3.10+ `X | None`, avoid `from __future__ import annotations` unless already present or needed for forward references, avoid `Any` when a concrete type/union fits, and scope `# type: ignore` to untyped third-party imports such as `mpmath`.

### `__main__` Runners And Prohibited `import re`

Runnable non-`__init__.py` modules start `if __name__ == "__main__":` with `doctest.testmod()` exiting `1` on failure (followed in equation modules by optional `trusted_dict` checks, keeping test-only imports inside doctests or `__main__`; see [Code Test Policy](../validation/code-test-policy.md)). Raise `ValueError` on invalid inputs and `warnings.warn()` on non-fatal diagnostics.

`import re` is prohibited anywhere under `nrpy/`: use plain `str` methods instead (`coding_style.md` grandfathers named legacy sites until their next substantive edit; repository tools outside `nrpy/` such as `tools/kb_lint.py` are exempt).

Evidence (normative rule): `coding_style.md` (`## Python Coding Style`, `### Docstring Style`, `### Prohibited Dependencies`).

## Sources

- [coding_style.md](../../coding_style.md) - `## Python Coding Style`, `### Docstring Style`, `### Module Docstring Format`, `### Prohibited Dependencies`

## See Also

- Depends on: [Contribution Style And Static Analysis](contribution-style-and-static-analysis.md)
- Depends on: [Code Test Policy](../validation/code-test-policy.md)
- Depends on: [Test Oracles And Safe Updates](../validation/test-oracles-and-safe-updates.md)
