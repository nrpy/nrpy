# Build And Run

> Package installation, contributor environment setup, `PYTHONPATH` execution, backend toolchain prerequisites, and the `nrpyinline` CLI entry point. · Status: confirmed
> Up: [Architecture](index.md)

## Summary

Install published NRPy with `python -m pip install nrpy` (`python_requires >= 3.6`), or install pinned contributor dependencies on Python `3.10+` (`black==26.5.1` floor) followed by `python -m pip install -e .`. Without an editable install, append `.` to `PYTHONPATH` before running `python -m nrpy.examples.<name>` from repository root.

## Detail

### Installation, Contributor Setup, And Backend Toolchains

Published and pinned contributor installation commands:

```bash
# End-user package install (Python >= 3.6; appends clang-format if missing)
python -m pip install nrpy

# Pinned contributor setup (requires Python >= 3.10 due to black==26.5.1)
python -m pip install -U -r requirements-dev.txt
python -m pip install -e .

# Direct checkout execution without editable install
export PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}."
```

Older Python cells in `.github/workflows/main.yml` strip version pins from `requirements-dev.txt`; do not use unpinned fallback installs for local contributor work.

Toolchain requirements depend on the target backend:
- **Standalone BHaH C/CUDA**: C compiler (plus `nvcc` for `--cuda`), `make`, Linux or macOS for executable `Makefile` generation (`create_lib=False` rejects other `platform.system()` values; see [Wave Equation Generators](../examples/wave-equation-generators.md)), and GSL for TwoPunctures, TOVola, SEOBNR, PN, and massive-geodesic examples.
- **Einstein Toolkit & CarpetX / superB / JAX**: Python generates thorn, chare, or `PyFunction` trees; compiling/running requires an Einstein Toolkit checkout, Charm++ (`charmc`/`charmrun`), or `jax`/`jaxlib`/`numpy`.
- **Dendro**: CMake `>= 3.18`, GNU C/C++17/Fortran, MPI (C/C++/Fortran), OpenMP, GSL, BLAS/LAPACK, Git, and network access on first configure for `Dendrolib`/`toml11`/`spdlog` ([Dendro Project Assembly](../infrastructures/dendro/project-assembly-and-emitters.md)).

Evidence (descriptive behavior): `setup.py`; `requirements-dev.txt`; `BHaH/Makefile_helpers.py::_validate_inputs`; `Dendro/CMakeLists.py::build_and_run_instructions`.

### Embedded-Snippet Command: `nrpyinline`

`setup.py` registers console script `nrpyinline=bin.nrpyinline:main`:

```bash
nrpyinline <text-or-source-file>
```

`run_script_from_file()` strips comment prefixes, captures blocks delimited by exact `NRPYSTART` and `NRPYEND` lines (raising `ValueError` on nested starts, unmatched ends, or unclosed blocks), and executes dedented Python blocks in one shared namespace with current-process permissions—run it only on trusted files (`main()` returns exit `1` on missing files, marker errors, or `SyntaxError`, while runtime exceptions propagate with tracebacks). Do not invoke `python bin/nrpyinline.py <file>` directly: its `__main__` block runs doctests and ignores CLI file arguments.

Evidence (descriptive behavior): `setup.py`; `bin/nrpyinline.py::run_script_from_file`; `bin/nrpyinline.py::main`.

## Sources

- [README.md](../../README.md) - `## Installation`, `## Prerequisites by Workflow`, `## First Successful Run`, `## Contributor Setup`
- [setup.py](../../setup.py) - `check_python_version`, `read_requirements_file`, `entry_points`
- [requirements.txt](../../requirements.txt) - runtime dependencies
- [requirements-dev.txt](../../requirements-dev.txt) - pinned contributor dependencies
- [Black 26.5.1 package metadata](https://pypi.org/pypi/black/26.5.1/json) - `requires_python >= 3.10`
- [bin/nrpyinline.py](../../bin/nrpyinline.py) - `strip_comment_markers`, `run_script_from_file`, `main`
- [Makefile_helpers.py](../../nrpy/infrastructures/BHaH/Makefile_helpers.py) - `_validate_inputs`
- [CMakeLists.py](../../nrpy/infrastructures/Dendro/CMakeLists.py) - `build_and_run_instructions`
- [dendro_bssn.py](../../nrpy/examples/dendro_bssn.py) - `main`
- [main.yml](../../.github/workflows/main.yml) - `static-analysis`, `codegen-ubuntu`, `codegen-mac`

## See Also

- See also: [Overview](overview.md)
- See also: [Contribution Style And Static Analysis](contribution-style-and-static-analysis.md)
- Example: [Wave Equation Generators](../examples/wave-equation-generators.md)
