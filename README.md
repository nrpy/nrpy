[![Python CI](https://github.com/nrpy/nrpy/actions/workflows/main.yml/badge.svg)](https://github.com/nrpy/nrpy/actions/workflows/main.yml)

# NRPy 2: Python/SymPy-Based Code Generation for Numerical Relativity... and Beyond!

## Quick start, Step 1:

```
pip install nrpy
```

## Quick start, Step 2: Choose a project to build.

### BlackHoles@Home infrastructure (standalone): Choose a project, run the provided command, then follow the instructions for compiling & running the generated C code.

* **Wave equation solver**
  - **Cartesian coordinates**:
    ```
    python3 -m nrpy.examples.wave_equation_cartesian
    ```
  - **Curvilinear coordinates**:
    ```
    python3 -m nrpy.examples.wave_equation_curvilinear
    ```
* **General relativity**
  - **Two black holes collide**:
    ```
    python3 -m nrpy.examples.two_blackholes_collide
    ```
  - **Black hole spectroscopy**:
    ```
    python3 -m nrpy.examples.blackhole_spectroscopy
    ```
  - **Spinning black hole**:
    ```
    python3 -m nrpy.examples.spinning_blackhole
    ```
  - **Binary black hole initial data, courtesy NRPyElliptic**:
    ```
    python3 -m nrpy.examples.nrpyelliptic_conformally_flat
    ```

### Einstein Toolkit infrastructure: Choose a project to build, run the provided command. Check the `examples/et_*` directory for a ThornList and parameter file. Thorns will be output to `project/`

* **Wave equation solver**
  - **Cartesian coordinates, with Carpet AMR infrastructure**:
    ```
    python3 -m nrpy.examples.carpet_wavetoy_thorns
    ```
* **General relativity: Generate Baikal and BaikalVacuum thorns**
  - **Cartesian coordinates, with Carpet AMR infrastructure**:
    ```
    python3 -m nrpy.examples.carpet_baikal_thorns
    ```
### superB infrastructure: Choose a project, run the provided command, then follow the instructions for installing Charm++, compiling and running the generated C++ code.

* **General relativity**
  - **Two black holes collide**:
    ```
    python3 -m nrpy.examples.superB_two_blackholes_collide
    ```
  - **Black hole spectroscopy**:
    ```
    python3 -m nrpy.examples.superB_blackhole_spectroscopy
    ```
  - **Binary black hole initial data, courtesy NRPyElliptic**:
    ```
    python3 -m nrpy.examples.superB_nrpyelliptic_conformally_flat
    ```

## Quick Start, Step 3

1. If working with a BlackHoles@Home project: follow the directions at the end of the code generation, starting with "Now go into `[directory name]` and type `make` to build, then `[executable]` to run. Parameter file can be found in `[parameter filename]`."
  1. Parameter files are text files, making it easy to adjust simulation parameters.
  1. In addition, parameters can be set at the command line. For example, `wavetoy` has a parameter `convergence_factor` that increases the resolution (technically `Nx=Ny=Nz`) by this factor. To output at twice the resolution, simply run `./wavetoy 2.0`, and a new file will be output `out0d-conv_factor2.00.txt`, which contains data at 2x the resolution.
  1. Analyze the output from `out0d-conv_factor1.00.txt` and `out0d-conv_factor2.00.txt` in e.g., `gnuplot`.
1. If working with an Einstein Toolkit project, the output will be Einstein Toolkit modules (thorns). You'll want to either copy or link them to an arrangement in `arrangements/[subdirectory]/`, then add the thorns to your `ThornList`, and compile.
1. If working with a superB project: install Charm++ following the instructions in [Charm++ documentation](https://charm.readthedocs.io/en/latest/charm++/manual.html#installing-charm). Then, follow the directions at the end of the code generation, starting with "Now go into `[directory name]` and type `make` to build, then `[executable]` to run. Parameter file can be found in `[parameter filename]`." As in a BlackHoles@Home project parameters can be changed and output from `out0d-conv_factor1.00.txt` and `out0d-conv_factor2.00.txt` can be analyzed using e.g., `gnuplot`.

# Contributing to NRPy 2 and running locally

Want to contribute to NRPy 2? Great! First clone the NRPy 2.0 repo:
```
git clone https://github.com/nrpy/nrpy.git
```

Next, you'll want to make sure your development environment is consistent with what GitHub Actions expects:
```
cd nrpy
pip install -U -r requirements-dev.txt
```

Finally, to run anything in the NRPy repo, you'll need to set your `PYTHONPATH` appropriately. If you're using bash, attach the following line to the bottom of your `.bashrc` file:
```
export PYTHONPATH=$PYTHONPATH:.
```

Once this is set up, you can run any Python script in the NRPy 2 repo from the repository's root directory. For example,
```
python3 nrpy/helpers/cse_preprocess_postprocess.py
```
will run all the doctests in that file.

# Key Improvements over `nrpytutorial` version of NRPy (NRPy 1):

## Easy Installation
- NRPy has been transformed into a proper Python project and can now be installed via pip! Use the command `pip install nrpy` to get started.
- Visit our [PyPI page](https://pypi.org/project/nrpy) for more information.
  - With pip, it's easier than ever to build your own projects based on NRPy.
  - You can now generate a complete C project from start to finish without the need for running a Jupyter notebook.
    - For instance, running `pip install nrpy && python3 -m nrpy.examples.two_blackholes_collide` will generate a C code project that evolves Brill-Lindquist forward in time using the BSSN formulation.
  - Check out [GitHub README](https://github.com/nrpy/nrpy/blob/main/README.md) for instructions on generating other pre-baked example codes... or use them to generate your own codes!

## Python 3.6+ Features
- NRPy now makes use of Python features introduced in version 3.6 and above, such as f-strings.
- The code is now more Pythonic and incorporates objects where useful. Global variables are now a thing of the past!

## User-friendly
- It's much simpler to work with NRPy now; you no longer have to read the source code of each function you call.
  - Facilitating this, you'll find:
    - Docstrings for all functions, classes, and modules.
    - Type hints across all modules; `mypy --strict` passes.
    - Numerous doctests.
    - Code formatted with Black.
    - Stricter linting.

## Improved Continuous Integration
- GitHub Actions now checks all files within the repo and will fail if any of the following conditions are met:
  - Doctest failure
  - pylint score less than 9.5
  - Black needs to reformat any `.py` file
  - `mypy --strict` fails on any `.py` file
  - Generating and compiling all examples from the pip-installed NRPy fresh from the latest git commit fails.

## More Extensible
- The "SENR" infrastructure has been replaced with "BHaH" (BlackHoles@Home). All BHaH-specific functionality is located in `nrpy/infrastructures/BHaH/`.
  - While BHaH currently only supports single grids, multi-patch support will be added soon.
  - You'll notice the old paramstruct has been broken into commondata_struct (data shared by all grids) and griddata (contains data specific to a particular grid).
    - Adding multi-patch support is a matter of setting commondata.NUMGRIDS > 1 and all the special algorithms.
  - There is a common but slightly customizable `main.c` file used by all BHaH codes, see `nrpy/infrastructures/BHaH/main_c.py`. This should greatly minimize code duplication in BHaH codes.
  - Parameter files are now supported, as well as advanced command-line input.
  - The `infrastructures/` directory includes helper functions for specific infrastructures. It currently contains BHaH and CarpetX subdirectories, with more to come.
  - `Cparameters` has been renamed to `CodeParameters`, allowing for future extensions of NRPy to output kernels in Python, Fortran, etc.
  - Rewritten expression validation infrastructure, to make it easier to validate newly added sympy expressions -- regardless of how complex they are.

## Plans for Old nrpytutorial Code
- We'll migrate the Jupyter notebooks to [a new nrpytutorial GitHub repo](https://github.com/nrpy/nrpytutorial) as they are updated to NRPy 2.
- All the core `.py` files from nrpytutorial have been modernized & ported to NRPy 2.
- What `.py` files remain in nrpytutorial will be ported to NRPy 2.
