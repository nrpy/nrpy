# Static Analysis

> Local and CI static-analysis behavior for NRPy Python changes. · Status: confirmed · Last reconciled: 07-13-2026
> Up: [Validation](index.md)

## Summary

For modified handwritten Python, run `black .` from repository root only in an
isolated, user-owned intended-change worktree or copy with no unrelated
modifications, then inspect its diff. Run
`./.github/single_file_static_analysis.sh <path.py>` for each file individually;
the full gate must pass and report Pylint **10.00/10.00**. Current wrapper and
CI thresholds do not by themselves guarantee that policy result.

## Detail

Current `coding_style.md` supplies the contributor commands, while [Code Test
Policy](code-test-policy.md) owns their prospective test-policy requirements.
Generated trusted `*/tests/*.py` data is exempt from per-file static analysis;
its handwritten owner is not. A directory path never exempts legacy or
exceptional handwritten Python, and this data exception does not admit new
executable Python under ordinary `tests/` stores. See [Test Oracles And Safe
Updates](test-oracles-and-safe-updates.md) for the store and update rules.

The single-file script expects exactly one Python file argument and must be run
from repository root because it resolves `.pylintrc` there. It verifies the file
exists, prints the Python version, then schedules `black --check`, `isort
--check-only`, `mypy --strict --pretty --allow-untyped-calls`, `pylint` with
`.pylintrc`, `pydocstyle`, `darglint -v 2`, and `python3 <file>`. Its Pylint
gate fails below `9.91`. Executing a file can run module-specific validation,
generation, or other `__main__` side effects, and analysis tools can create or
update caches or other filesystem output. The wrapper's mypy command supplies no
`--cache-dir`, and `.mypy.ini` sets no `cache_dir`, so the invocation must not be
assumed write-free. Inspect the target, the wrapper's command construction, and
every invoked tool. The wrapper interpolates the target path into command strings
dispatched through `eval`. Before invocation, require a repository-relative
argument whose resolved target stays inside the repository, is a regular
non-symlink file, has no `..` component, does not begin with `-`, and matches
`^[A-Za-z0-9_./-]+$`. If that gate fails, do not invoke the wrapper; hardening it
to argument-vector execution requires separately authorized source work. For an
accepted path, run the wrapper only in an isolated, user-owned intended-change
worktree or copy with no unrelated modifications. Disable or redirect every
writable tool cache and filesystem output to an owned disposable location, and
inspect repository status afterward.

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
`9.91`. Both `.pylintrc` and `.pylintrc_python36` set `fail-under=10`, but the
wrapper and workflow suppress Pylint's own nonzero result while parsing ratings
against their weaker floors. Therefore script or job success under these paths
does not guarantee the required Pylint 10.00. Workflow YAML proves configured
matrix and command shape, not a latest successful run.

Claim status: stale; contradiction: CONTR-0003.
See [CONTR-0003](../contradictions.md#contr-0003) for the enforcement gap.

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
- [../../.pylintrc](../../.pylintrc) - `[MASTER]`
- [../../.pylintrc_python36](../../.pylintrc_python36) - `[MASTER]`
- [../../.darglint](../../.darglint) - `[darglint]`

## See Also

- Parent: [Validation](index.md)
- Depends on: [Code Test Policy](code-test-policy.md)
- Depends on: [Test Oracles And Safe Updates](test-oracles-and-safe-updates.md)
- See also: [Generated Project CI](generated-project-ci.md)
- See also: [Lint Checks](../lint/CHECKS.md)
- See also: [Workflows](../workflows.md)
