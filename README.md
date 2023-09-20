[![Python CI](https://github.com/nrpy/nrpy/actions/workflows/main.yml/badge.svg)](https://github.com/nrpy/nrpy/actions/workflows/main.yml)

# NRPy+ 2.0: Next-generation Python/SymPy-Based Code Generation for Numerical Relativity... and Beyond!

## Quick start:
1. `pip install nrpy`
2. Choose a project to build, run the provided command, then follow the instructions for compiling & running the generated C code.
   3. Wave equation solver
      3. **Cartesian coordinates**: `python3 -m nrpy.examples.wave_equation_cartesian`
      3. **Curvilinear coordinates**: `python3 -m nrpy.examples.wave_equation_curvilinear`
   4. General relativity
      5. **Two black holes collide**: `python3 -m nrpy.examples.two_blackholes_collide`

## Slow start, in case you want to develop NRPy+ directly:
1. Clone this repository.
2. `cd nrpy` , then e.g., run an example code: `python3 nrpy/examples/wave_equation_cartesian` , then follow the instructions for compiling & running the generated C code.
1. Run `./wavetoy` . It'll output at various times the solution from the point closest to r=0, as well as the relative error compared to the exact solution at that point to a file `out0d-conv_factor1.00.txt` .
1. All BlackHoles@Home (BHaH) infrastructure codes now support parameter files, and each code generation automatically outputs a parfile named `[project name].par`. So you'll find an editable `wavetoy.par`.
1. In addition, users can override certain parameter file parameters at the command line. E.g., `wavetoy` has a parameter `convergence_factor` that increases the resolution (technically `Nx=Ny=Nz`) by this factor. To output at twice the resolution, simply run `./wavetoy 2.0`, and a new file will be output `out0d-conv_factor2.00.txt`, which contains data at 2x the resolution.
1. Analyze the output from `out0d-conv_factor1.00.txt` and `out0d-conv_factor2.00.txt` in e.g., `gnuplot`.


## Key Improvements over NRPy+ 1.0:

### Easy Installation
- NRPy+ has been transformed into a proper Python project and can now be installed via pip! Use the command `pip install nrpy` to get started.
- Visit our [PyPI page](https://pypi.org/project/nrpy) for more information.
  - With pip, it's easier than ever to build your own projects based on NRPy+.
  - You can now generate a complete C project from start to finish without the need for running a Jupyter notebook.
    - For instance, running `pip install nrpy && python3 -m nrpy.examples.two_blackholes_collide` will generate a C code project that evolves Brill-Lindquist forward in time using the BSSN formulation.
  - Check out [GitHub README](https://github.com/nrpy/nrpy/blob/main/README.md) for instructions on generating other pre-baked example codes... or use them to generate your own codes!

### Python 3.6+ Features
- NRPy+ now makes use of Python features introduced in version 3.6 and above, such as f-strings.
- The code is now more Pythonic and incorporates objects where useful. Global variables are now a thing of the past!

### User-friendly
- It's much simpler to work with NRPy+ now; you no longer have to read the source code of each function you call.
  - Facilitating this, you'll find:
    - Docstrings for all functions, classes, and modules.
    - Type hints across all modules; `mypy --strict` passes.
    - Numerous doctests.
    - Code formatted with Black.
    - Stricter linting.

### Improved Continuous Integration
- GitHub Actions now checks all files within the repo and will fail if any of the following conditions are met:
  - Doctest failure
  - pylint score less than 9.5
  - Black needs to reformat any `.py` file
  - `mypy --strict` fails on any `.py` file
  - Generating and compiling all examples from the pip-installed NRPy+ fresh from the latest git commit fails.

### More Extensible
- The "SENR" infrastructure has been replaced with "BHaH" (BlackHoles@Home). All BHaH-specific functionality is located in `nrpy/infrastructures/BHaH/`.
  - While BHaH currently only supports single grids, multi-patch support will be added soon.
  - You'll notice the old paramstruct has been broken into commondata_struct (data shared by all grids) and griddata (contains data specific to a particular grid).
    - Adding multi-patch support is a matter of setting commondata.NUMGRIDS > 1 and all the special algorithms.
  - There is a common but slightly customizable `main.c` file used by all BHaH codes, see `nrpy/infrastructures/BHaH/main_c.py`. This should greatly minimize code duplication in BHaH codes.
  - Parameter files are now supported, as well as advanced command-line input.
  - The `infrastructures/` directory includes helper functions for specific infrastructures. It currently contains BHaH and CarpetX subdirectories, with more to come.
  - `Cparameters` has been renamed to `CodeParameters`, allowing for future extensions of NRPy+ to output kernels in Python, Fortran, etc.
  - Rewritten expression validation infrastructure, to make it easier to validate newly added sympy expressions -- regardless of how complex they are.

## Plans for Old nrpytutorial Code
- We'll migrate the Jupyter notebooks to [a new nrpytutorial GitHub repo](https://github.com/nrpy/nrpytutorial) as they are updated to NRPy+ 2.0.
- All the core `.py` files from nrpytutorial have been modernized & ported to NRPy+ 2.0.
- What `.py` files remain in nrpytutorial will be ported to NRPy+ 2.0.