[![Python CI](https://github.com/nrpy/nrpy/actions/workflows/main.yml/badge.svg)](https://github.com/nrpy/nrpy/actions/workflows/main.yml)

# NRPy

NRPy is a Python/SymPy-based symbolic code generation toolkit for numerical relativity and relativistic astrophysics. It uses symbolic equation modules to generate standalone C/CUDA applications, Einstein Toolkit thorns, Charm++ projects, and JAX/Python workflows. In its current public form, NRPy provides the core infrastructure for single-patch BlackHoles@Home (BH@H) numerical relativity simulations in Cartesian, curvilinear, and singular coordinate systems.

At a high level:

- `nrpy/examples` contains runnable project generators.
- `nrpy/equations`, `nrpy/infrastructures`, and `nrpy/helpers` contain the symbolic physics, backend-specific code generation, and shared support code that those generators use.

## Citation

NRPy is the product of many scientists' work. [Please cite the relevant NRPy papers in your publications](CITATION.md); these citations go a long way toward proper acknowledgement of these efforts and support many junior scientists' careers.

## Installation

For end users who want to generate projects from the published package:

```bash
python -m pip install nrpy
```

For contributors working from a clone of this repository:

```bash
git clone https://github.com/nrpy/nrpy.git
cd nrpy
python -m pip install -U -r requirements-dev.txt
python -m pip install -e .
```

After an editable install, you can run generators from the repository root with `python -m nrpy.examples.<example_name>`.

## Prerequisites by Workflow

`python -m pip install nrpy` installs the Python package, but many generated projects also need external tools:

- Standalone BHaH examples: a C compiler and `make`; some also link against GSL
- Waveform and geodesic examples: a C compiler, `make`, and GSL
- Einstein Toolkit / Carpet / CarpetX generators: the Python package install is enough to generate thorns; an existing Einstein Toolkit environment is needed to build and run the generated thorns
- `superB` examples: a Charm++ toolchain; some also link against GSL
- JAX example generation: Python tooling for the generated project, plus JAX in that generated environment as needed

If you are new to NRPy, start with one of the standalone BHaH examples below that does not require GSL. Those examples are the closest path to the current public single-patch BH@H workflow.

## First Successful Run

This is the safest first workflow because it only requires Python, a C compiler, and `make`.

### 1. Generate a standalone project

```bash
python -m nrpy.examples.wave_equation_cartesian
```

This writes a buildable project to `project/wave_equation_cartesian/`.

### 2. Build it

```bash
cd project/wave_equation_cartesian
make
```

### 3. Run it

```bash
./wave_equation_cartesian
```

The generator also creates a parameter file named `wave_equation_cartesian.par`. Many standalone BHaH examples accept command-line overrides as well; for this example, `convergence_factor` is one supported input.

### 4. Inspect the outputs

Expect generated diagnostics and data products in the project directory. The exact filenames vary by example, but the important milestone is:

- generation succeeds
- `make` succeeds
- the executable runs without further manual project editing

For this example, one simple command-line override is `convergence_factor`. For example:

```bash
./wave_equation_cartesian 2.0
```

This reruns the example at higher resolution and writes output files such as `out0d-conv_factor2.00.txt`.

## Lightweight Single-Patch BH@H Example

If you want the closest public example of the single-patch BH@H numerical relativity workflow, use:

```bash
python -m nrpy.examples.two_blackholes_collide
```

This example:

- evolves Brill-Lindquist binary black hole initial data
- uses a single-patch curvilinear spherical-coordinate setup
- serves as a compact public example of the single-patch BH@H workflow
- is intentionally lightweight: it is cheap enough to run in seconds even on very modest hardware

## Good Next Examples

Once the first standalone workflow works, these are reasonable next steps:

- Single-patch BH@H NR example: `python -m nrpy.examples.two_blackholes_collide`
- Lightweight elliptic initial-data workflow: `python -m nrpy.examples.nrpyelliptic_conformally_flat`
- Waveform generation, requires GSL: `python -m nrpy.examples.seobnrv5_aligned_spin_inspiral -seobnrv5_bob`
- Einstein Toolkit / CarpetX thorn generation; building the generated thorns requires an ET environment: `python -m nrpy.examples.carpetx_wavetoy_thorns`

## Project Families and Example Generators

### Standalone BHaH Generators

BHaH is NRPy's standalone application infrastructure. In this public repository, it provides the core capabilities for single-patch numerical relativity evolutions in Cartesian, curvilinear, and singular coordinate systems, and these generators typically produce buildable C or CUDA-ready projects under `project/<name>/`.

- Wave equation solvers: `python -m nrpy.examples.wave_equation_cartesian`, `python -m nrpy.examples.wave_equation_curvilinear`, `python -m nrpy.examples.wave_equation_multicoordinates`
- Black hole evolution and diagnostics: `python -m nrpy.examples.two_blackholes_collide`, `python -m nrpy.examples.blackhole_spectroscopy` (requires GSL), `python -m nrpy.examples.spinning_blackhole`
- Elliptic / initial-data workflows: `python -m nrpy.examples.nrpyelliptic_conformally_flat`

### Einstein Toolkit and CarpetX Generators

These generators target the Einstein Toolkit rather than producing standalone executables. ETLegacy refers to the classic Einstein Toolkit / Carpet-style target; CarpetX targets the newer CarpetX driver stack.

- ETLegacy / Carpet workflows: `python -m nrpy.examples.carpet_wavetoy_thorns`, `python -m nrpy.examples.carpet_baikal_thorns`
- CarpetX workflows: `python -m nrpy.examples.carpetx_wavetoy_thorns`, `python -m nrpy.examples.carpetx_baikal_thorns`

These generators write thorn directories under `project/<name>/`. You only need an Einstein Toolkit checkout when you want to compile or run those generated thorns.

### superB / Charm++ Generators

These generators target `superB`, NRPy's Charm++-based infrastructure for distributed-memory workflows and a scalable extension of single-patch BH@H-generated code.

- GR evolution and spectroscopy: `python -m nrpy.examples.superB_two_blackholes_collide`, `python -m nrpy.examples.superB_blackhole_spectroscopy`
- Elliptic initial data: `python -m nrpy.examples.superB_nrpyelliptic_conformally_flat`

### Compact-Object and Waveform Generators

These generators focus on compact-binary dynamics, inspiral, and waveform modeling.

- SEOBNRv5 / BOB waveform projects: `python -m nrpy.examples.seobnrv5_aligned_spin_inspiral -seobnrv5_bob`, `python -m nrpy.examples.seobnrv5_aligned_spin_inspiral -seobnrv5_nrpy`
- SEBOBv2 standalone project, requires GSL: `python -m nrpy.examples.sebobv2`
- JAX project generation: `python -m nrpy.examples.sebobv1_jax`
- Post-Newtonian utility: `python -m nrpy.examples.nrpypn_quasicircular_momenta`

### Additional Physics Generators

These generators cover additional relativistic physics applications beyond the main PDE and thorn-generation workflows.

- Neutron-star and matter-related workflows: `python -m nrpy.examples.tovola_neutron_star` (requires GSL), `python -m nrpy.examples.hydro_without_hydro` (requires GSL)
- Geodesic integration: `python -m nrpy.examples.mass_geodesic_integrator`, `python -m nrpy.examples.photon_geodesic_integrator`

### Specialized Utilities

- Apparent horizon library generation: `python -m nrpy.examples.bhahaha`

## What Gets Generated?

- Standalone BHaH examples usually generate a project directory with C source, headers, a Makefile, a parameter file, and a runnable executable target for single-patch Cartesian or curvilinear-coordinate workflows.
- Waveform and geodesic generators also produce standalone projects, but they typically require GSL at build time.
- Some additional standalone physics examples, including `blackhole_spectroscopy`, `hydro_without_hydro`, `tovola_neutron_star`, and `sebobv2`, also require GSL at build time.
- Einstein Toolkit and CarpetX generators produce thorn directories that you copy or link into an Einstein Toolkit checkout; an ET environment is needed for the subsequent build/run step, not for thorn generation itself.
- `superB` generators produce Charm++ projects rather than plain standalone executables, and some of them also require GSL.
- `sebobv1_jax` generates a Python/JAX project instead of a C executable.

Most generators print project-specific next steps when they finish, but the exact instructions vary by example.

## Repository Map

- `nrpy/examples`: runnable generators and project entry points
- `nrpy/infrastructures`: backend-specific code generation layers, including BHaH, ETLegacy, CarpetX, `superB`, and JAX
- `nrpy/equations`: symbolic physics and equation modules
- `nrpy/helpers`: shared code generation and utility helpers
- `nrpy/tests`: reference-metric validation modules
- `bin/nrpyinline.py`: advanced utility for extracting and running embedded NRPy snippets from text or source files

## Contributor Setup

The recommended local workflow is:

```bash
python -m pip install -U -r requirements-dev.txt
python -m pip install -e .
```

You do not need to modify `PYTHONPATH` if you install the repository in editable mode.

Useful local sanity checks include:

```bash
python -m nrpy.examples.wave_equation_cartesian
python -m nrpy.examples.sebobv1_jax
```

For broader validation, see `.github/workflows/main.yml`, which exercises static analysis, standalone project generation/builds, Einstein Toolkit validation, Charm++ / `superB` workflows, and SEOB consistency checks.

## Legacy Note

Older notebook-centered NRPy+ (NRPy 1.0) material lives in the separate [`nrpytutorial`](https://github.com/zachetienne/nrpytutorial) repository. This repository focuses on the current package and generator-based workflow.
