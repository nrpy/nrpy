# Build And Run

> Compile the supported install paths, prerequisites, and first runnable workflow. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [Architecture](index.md)

## Summary

End users can install the published package with `python -m pip install nrpy`. Contributors normally install development requirements and then install the repository in editable mode. Without an editable install, direct example execution from the repository root requires appending `.` to `PYTHONPATH`.

## Detail

The published-package path is:

```bash
python -m pip install nrpy
```

The contributor path is:

```bash
python -m pip install -U -r requirements-dev.txt
python -m pip install -e .
```

After editable install, generators can run as modules from the repository root, for example `python -m nrpy.examples.wave_equation_cartesian`. If running Python examples directly from the checkout without relying on editable install, set:

```bash
export PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}."
```

The first recommended standalone workflow is `python -m nrpy.examples.wave_equation_cartesian`, which writes `project/wave_equation_cartesian/`. Build with `make` inside that directory and run `./wave_equation_cartesian`. The project also includes `wave_equation_cartesian.par`, and the example accepts `convergence_factor` as a command-line override.

Prerequisites depend on the generated project family. Standalone BHaH workflows need a C compiler and `make`; some standalone, waveform, geodesic, and additional physics workflows also need GSL. Einstein Toolkit and CarpetX generators can generate thorns from the Python package, but compiling or running those thorns requires an Einstein Toolkit environment. `superB` workflows require a Charm++ toolchain, and generated JAX workflows need suitable Python/JAX tooling in the generated environment.

The package setup requires Python 3.6 or newer, reads `requirements.txt`, and appends `clang-format` if the executable is not already present. The GitHub workflow validates install-and-generate paths by installing the package, generating multiple projects on Ubuntu and macOS, and building generated C projects where applicable.

## Sources

- [README.md](../../README.md) - `## Installation`, `## Prerequisites by Workflow`, `## First Successful Run`, `## Contributor Setup`
- [requirements.txt](../../requirements.txt) - runtime dependency list
- [requirements-dev.txt](../../requirements-dev.txt) - development dependency list
- [setup.py](../../setup.py) - `check_python_version`, `read_requirements_file`, `setup(...)`
- [raw/source-docs/original-agents.md](../../raw/source-docs/original-agents.md) - `## Required Checks`
- [.github/workflows/main.yml](../../.github/workflows/main.yml) - `static-analysis`, `codegen-ubuntu`, `codegen-mac`

## See Also

- [Architecture](index.md)
- [Overview](overview.md)
- [Symbolic Codegen Lifecycle](symbolic-codegen-lifecycle.md)
- [Contribution Style And Static Analysis](contribution-style-and-static-analysis.md)
