# Build And Run

> Compile the supported install paths, command entry points, prerequisites, and first runnable workflow. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [Architecture](index.md)

## Summary

End users can install the published package with `python -m pip install nrpy`. Contributors using the exact current pinned development-requirements command need at least Python 3.10, then install the repository in editable mode. This is a minimum derived from one pinned tool, not a guarantee that every future dependency resolves on every Python 3.10+ environment. Without an editable install, direct example execution from the repository root requires appending `.` to `PYTHONPATH`.

## Detail

The published-package path is:

```bash
python -m pip install nrpy
```

The current pinned contributor path is:

```bash
python -m pip install -U -r requirements-dev.txt
python -m pip install -e .
```

Run that exact `requirements-dev.txt` command with Python 3.10 or newer. The file pins `black==26.5.1`, whose official package metadata requires Python 3.10 or newer. This is narrower than `setup.py`'s Python 3.6 minimum for the installable runtime package. The CI workflow handles selected older Python jobs differently: it strips pinned versions from `requirements-dev.txt` and lets pip select compatible releases. That CI fallback is not a reproducible pinned local tool set, so this page does not present it as the normal contributor install command.

After editable install, generators can run as modules from the repository root, for example `python -m nrpy.examples.wave_equation_cartesian`. If running Python examples directly from the checkout without relying on editable install, set:

```bash
export PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}."
```

The first recommended standalone workflow is `python -m nrpy.examples.wave_equation_cartesian`, which writes `project/wave_equation_cartesian/`. Build with `make` inside that directory and run `./wave_equation_cartesian`. The project also includes `wave_equation_cartesian.par`, and the example accepts `convergence_factor` as a command-line override.

Prerequisites depend on the generated project family. Standalone BHaH workflows need a C compiler and `make`; some standalone, waveform, geodesic, and additional physics workflows also need GSL. Einstein Toolkit and CarpetX generators can generate thorns from the Python package, but compiling or running those thorns requires an Einstein Toolkit environment. `superB` workflows require a Charm++ toolchain, and generated JAX workflows need suitable Python/JAX tooling in the generated environment.

The package setup requires Python 3.6 or newer, reads `requirements.txt`, and appends `clang-format` if the executable is not already present. The GitHub workflow configures install-and-generate jobs that install the package, generate multiple projects on Ubuntu and macOS, and build generated C projects where the job lists a `make` step. Workflow configuration establishes intended CI coverage, not the result of any particular run.

### Embedded-Snippet Command: `nrpyinline`

`setup.py` declares one console-script entry point, `nrpyinline=bin.nrpyinline:main`. A normal or editable package install therefore exposes:

```bash
nrpyinline <text-or-source-file>
```

`run_script_from_file()` scans that file for lines which reduce to the exact markers `NRPYSTART` and `NRPYEND` after its comment-marker cleanup. It dedents and executes each captured Python block in one shared namespace, so definitions from an earlier block remain available to later blocks. Nested starts, an end without a start, and an unclosed final block raise `ValueError`.

This command executes embedded Python with the permissions of the current process; it is not a parser, validator, or sandbox. Use it only on trusted files. `main()` converts marker errors, missing files, and syntax errors to exit status 1, while other exceptions from embedded code propagate after `run_script_from_file()` prints a traceback.

Do not substitute `python bin/nrpyinline.py <file>` for the installed command. In the current source, the module's `__main__` block runs its doctests and does not call `main()`, so that direct form ignores the file argument as a command input.

## Sources

- [README.md](../../README.md) - `## Installation`, `## Prerequisites by Workflow`, `## First Successful Run`, `## Contributor Setup`
- [requirements.txt](../../requirements.txt) - runtime dependency list
- [requirements-dev.txt](../../requirements-dev.txt) - development dependency list
- [setup.py](../../setup.py) - `check_python_version`, `read_requirements_file`, `setup(...)`, `entry_points["console_scripts"]`
- [bin/nrpyinline.py](../../bin/nrpyinline.py) - `strip_comment_markers`, `run_script_from_file`, `main`, `if __name__ == "__main__"`
- [raw/source-docs/original-agents.md](../../raw/source-docs/original-agents.md) - `## Required Checks`
- [.github/workflows/main.yml](../../.github/workflows/main.yml) - `static-analysis`, `codegen-ubuntu`, `codegen-mac`
- [Black 26.5.1 package metadata](https://pypi.org/pypi/black/26.5.1/json) - `requires_python` and installation requirement; accessed 07-12-2026

## See Also

- Parent: [Architecture](index.md)
- See also: [Overview](overview.md)
- See also: [Symbolic Codegen Lifecycle](symbolic-codegen-lifecycle.md)
- See also: [Contribution Style And Static Analysis](contribution-style-and-static-analysis.md)
