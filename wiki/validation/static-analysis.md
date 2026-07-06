# Static Analysis

> Local and CI static-analysis behavior for NRPy Python changes. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [Validation](index.md)

## Summary

For Python source changes, run `black .` before commit and run
`./.github/single_file_static_analysis.sh <path.py>` for every modified Python
file, except generated trusted-value files under `*/tests/*.py`. The
single-file script checks formatting, imports, typing, linting, docstrings,
docstring arguments, and doctests.

## Detail

The preserved agent instructions require the single-file static-analysis script
on every modified Python file and explicitly forbid using only a hand-picked
subset. They also record the trusted-value exception: generated reference-value
files under `*/tests/*.py` are not hand-authored source modules and skip the
single-file static-analysis script.

The single-file script expects one Python file argument. It verifies the file
exists, prints the Python version, then runs `black --check`, `isort
--check-only`, `mypy --strict --pretty --allow-untyped-calls`, `pylint` with
`.pylintrc`, `pydocstyle`, `darglint -v 2`, and `python3 <file>` for doctests.
Its Pylint gate fails below `9.91`.

The GitHub `static-analysis` job runs across Ubuntu and Python-version matrices.
It clears the NRPy cache, installs runtime and development dependencies, skips
`__init__.py`, `project/`, `build/`, `*/tests/*`, `*manga*`, and visualization
scripts, then runs doctests for non-example modules, Black on newer Python,
isort and mypy on Python versions that support those checks, Pylint, pydocstyle,
and darglint. That CI job has its own inline Pylint threshold, while the
required local single-file script uses the stricter `9.91` threshold.

The config files supply the tool policy: `pyproject.toml` sets isort to the
Black profile; `.mypy.ini` enables strict typed definitions, no implicit
optional, redundant-cast and unused-ignore warnings, return-`Any` warnings, and
strict equality; `.pydocstyle` lists the accepted ignore set; `.darglint` uses
Sphinx docstring style; `.pylintrc` and `.pylintrc_python36` set naming,
formatting, disabled warnings, and version-specific Pylint behavior.

## Sources

- [../../raw/source-docs/original-agents.md](../../raw/source-docs/original-agents.md) - `## Required Checks`
- [../../.github/single_file_static_analysis.sh](../../.github/single_file_static_analysis.sh) - `run_test_step`
- [../../.github/workflows/main.yml](../../.github/workflows/main.yml) - `static-analysis`
- [../../pyproject.toml](../../pyproject.toml) - `[tool.isort]`
- [../../.mypy.ini](../../.mypy.ini) - `[mypy]`
- [../../.pydocstyle](../../.pydocstyle) - `[pydocstyle]`
- [../../.pylintrc](../../.pylintrc) - `[MESSAGES CONTROL]`
- [../../.pylintrc_python36](../../.pylintrc_python36) - `[MESSAGES CONTROL]`
- [../../.darglint](../../.darglint) - `[darglint]`

## See Also

- [Validation](index.md)
- [Generated Project CI](generated-project-ci.md)
- [Lint Checks](../lint/CHECKS.md)
- [Workflows](../workflows.md)
