# Static Analysis

> Local and CI static-analysis behavior for NRPy Python changes. · Status: confirmed · Last reconciled: 07-13-2026
> Up: [Validation](index.md)

## Summary

For Python source changes, run `black .` only in an owned clean detached
worktree or local copy outside shared `/work`. Run
`./.github/single_file_static_analysis.sh <path.py>` for every modified Python
file. The single-file script checks formatting, imports, typing, linting, and
docstrings, then executes the target file. That final step runs doctests only
when the target's `__main__` path invokes a doctest runner.

## Detail

Current `coding_style.md` owns the contributor commands: `black .` before
commit and the single-file script for each modified Python file. KB safe-
reproduction governance further requires blanket formatters to run only in an
owned clean detached worktree or local copy outside shared `/work`. Generated
reference-value files under `*/tests/*.py` are regenerated evidence rather than
hand-authored modules, but current governance declares no exception from the
single-file script when such a Python file is modified.

The single-file script expects exactly one Python file argument and must be run
from repository root because it resolves `.pylintrc` there. It verifies the file
exists, prints the Python version, then schedules `black --check`, `isort
--check-only`, `mypy --strict --pretty --allow-untyped-calls`, `pylint` with
`.pylintrc`, `pydocstyle`, `darglint -v 2`, and `python3 <file>`. Its Pylint
gate fails below `9.91`. Executing a file can run module-specific validation,
generation, or other `__main__` side effects; inspect unfamiliar targets first.

The configured GitHub `static-analysis` matrix has Ubuntu 22.04 and 24.04 with
Python `3.7.13`, `3.8.12`, `3.9.19`, and `3.x`, excluding only Ubuntu 24.04 with
3.7.13. Its file discovery skips `__init__.py`, `project/`, `build/`, every
`*/tests/*` path, and every path containing `manga`. It also names
`./nrpy/examples/visualization_scripts/*`, but that path does not match the
current `nrpy/examples/geodesic_visualizations/` directory; those current
companion scripts are therefore not excluded by that rule.

For discovered files, CI executes non-root, non-example modules as its doctest
step; applies Black only on Python minor version 10 or newer; skips isort and
mypy on 3.7.13 and 3.8.12; and applies Pylint, pydocstyle, and darglint on all
matrix cells. Its inline Pylint threshold is `9.5`, not the local script's
`9.91`. Workflow YAML proves this configured matrix and command shape, not a
latest successful run.

The config files supply the tool policy: `pyproject.toml` sets isort to the
Black profile; `.mypy.ini` enables strict typed definitions, no implicit
optional, redundant-cast and unused-ignore warnings, return-`Any` warnings, and
strict equality; `.pydocstyle` lists the accepted ignore set; `.darglint` uses
Sphinx docstring style; `.pylintrc` and `.pylintrc_python36` set naming,
formatting, disabled warnings, and version-specific Pylint behavior.

## Sources

- [../../coding_style.md](../../coding_style.md) - `## Python Coding Style`, `### Formatting`
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
