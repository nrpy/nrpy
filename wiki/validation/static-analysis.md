# Static Analysis

> Local and CI static-analysis tools, Pylint score rules, cache isolation, and Python 3.7 compatibility for NRPy Python modules. · Status: confirmed
> Up: [Validation](index.md)

## Summary

Modified handwritten Python modules require `black .` (run only in an isolated worktree or copy with no unrelated modifications) and `.github/single_file_static_analysis.sh <file>`. Newly added handwritten Python files must score Pylint **10.00/10.00**, while existing tracked handwritten files must not regress below their pre-change Pylint score (even when at or below `9.5`); neither `.github/single_file_static_analysis.sh` (`9.91`) nor CI (`9.5`) enforces that distinction automatically.

## Detail

### Pylint Policy, Trusted-Data Exemption, And Wrapper Execution

`coding_style.md` defines the static-analysis commands, while [Code Test Policy](code-test-policy.md) and [Test Oracles And Safe Updates](test-oracles-and-safe-updates.md) govern test placement and `*/tests/*.py` reference files. Generated `trusted_dict` modules under `*/tests/*.py` are exempt from per-file static analysis, whereas their handwritten defining modules (and any handwritten Python inside `tests/`) are not. Pre-change git tracking determines whether a file is new (`10.00/10.00` required) or existing (measure baseline and candidate with identical Pylint version/config and reject any decrease).

`.github/single_file_static_analysis.sh <file>` must run from repository root (where it resolves `.pylintrc`) on one validated repo-relative regular non-symlink path matching `^[A-Za-z0-9_./-]+$` (no `..` or leading `-`, because the script dispatches command strings via `eval`). It prints the Python version and runs `black --check`, `isort --check-only`, `mypy --strict --pretty --allow-untyped-calls`, `pylint --rcfile=.pylintrc`, `pydocstyle`, `darglint -v 2`, and `python3 <file>`. Because `python3 <file>` executes `__main__` and `mypy`/`pylint` write caches when unconfigured, run the wrapper only in an isolated copy with `HOME`, `XDG_CACHE_HOME`, `PYLINTHOME`, and `MYPY_CACHE_DIR` redirected to disposable directories.

Evidence (normative rule): `coding_style.md` (`## Static Analysis Configuration`); `.github/single_file_static_analysis.sh::run_test_step`.

### CI Matrix, Unchecked Pylint Thresholds, And Python 3.7 Floor

GitHub Actions `static-analysis` runs on Ubuntu 22.04 and 24.04 across Python `3.7.13`, `3.8.12`, `3.9.19`, and `3.x` (excluding Ubuntu 24.04 + `3.7.13`). File discovery skips `__init__.py`, `project/`, `build/`, `*/tests/*`, `*manga*`, and stale `./nrpy/examples/visualization_scripts/*` (leaving `nrpy/examples/geodesic_visualizations/` included). CI executes non-root, non-example modules for doctests, runs `black` only on Python `>= 3.10`, skips `isort`/`mypy` on `3.7.13` and `3.8.12`, and runs `pylint`, `pydocstyle`, and `darglint` on all cells.

Although `.pylintrc` and `.pylintrc_python36` set `fail-under=10`, both `.github/single_file_static_analysis.sh` (`9.91` floor) and `.github/workflows/main.yml` (`9.5` floor) suppress Pylint's exit code and compare parsed scores against fixed floors without checking pre-change baselines or new-file status—so a new file at `9.94` or a `10.00 -> 9.93` regression passes automation while a legacy `< 9.5` edit fails CI; reviewers must verify per-file scores manually.

Claim status: stale; contradiction: CONTR-0003. See [CONTR-0003](../contradictions.md#contr-0003).

Because `static-analysis` executes on Python `3.7.13`, `3.7` is the active syntax and runtime floor: newer syntax fails at parse time, while post-`3.7` runtime operators that parse on `3.8+` (such as PEP 584 `dict |= other`) raise `TypeError` on `3.7`/`3.8` despite passing `py_compile` on modern Python. Tool configs live in `pyproject.toml` (`isort` Black profile), `.mypy.ini` (strict typing, no implicit `Optional`, `warn_return_any`, `strict_equality`), `.pydocstyle`, `.darglint` (`sphinx`), `.pylintrc`, and `.pylintrc_python36`.

Evidence (CI behavior): `.github/workflows/main.yml` (`static-analysis`); `.github/single_file_static_analysis.sh`; `.pylintrc`; `setup.py`.

## Sources

- [coding_style.md](../../coding_style.md) - `## Python Coding Style`, `## Static Analysis Configuration`
- [single_file_static_analysis.sh](../../.github/single_file_static_analysis.sh) - `run_test_step`
- [main.yml](../../.github/workflows/main.yml) - `static-analysis`
- [setup.py](../../setup.py) - `python_requires` and classifiers
- [pyproject.toml](../../pyproject.toml) - `[tool.isort]`
- [.mypy.ini](../../.mypy.ini) - `[mypy]`
- [.pydocstyle](../../.pydocstyle) - `[pydocstyle]`
- [.pylintrc](../../.pylintrc) - `[MASTER]`
- [.pylintrc_python36](../../.pylintrc_python36) - `[MASTER]`
- [.darglint](../../.darglint) - `[darglint]`

## See Also

- Depends on: [Code Test Policy](code-test-policy.md)
- Depends on: [Test Oracles And Safe Updates](test-oracles-and-safe-updates.md)
- See also: [Generated Project CI](generated-project-ci.md)
