[![Python CI](https://github.com/nrpy/nrpy/actions/workflows/master.yml/badge.svg)](https://github.com/nrpy/nrpy/actions/workflows/master.yml)

# NRPy+: Python-Based Code Generation for Numerical Relativity... and Beyond!

## Quick start:
1. `pip install nrpy`
2. Solving the wave equation in Cartesian coordinates: `python3 -m nrpy.examples.wave_equation_cartesian` , then follow the instructions for compiling & running the generated C code.
3. Solving the wave equation in curvilinear coordinates: `python3 -m nrpy.examples.wave_equation_curvilinear` , then follow the instructions for compiling & running the generated C code.
4. **Two black holes collide**: Solving the general relativistic field equations in curvilinear coordinates: `python3 -m nrpy.examples.two_blackholes_collide` , then follow the instructions for compiling & running the generated C code.

## Slow start, in case you want to develop NRPy+ directly:
1. Clone this repository.
2. `cd nrpy` , then e.g., run an example code: `python3 nrpy/examples/wave_equation_cartesian` , then follow the instructions for compiling & running the generated C code.
1. Run `./wavetoy` . It'll output at various times the solution from the point closest to r=0, as well as the relative error compared to the exact solution at that point to a file `out0d-conv_factor1.00.txt` .
1. All BlackHoles@Home (BHaH) infrastructure codes now support parameter files, and each code generation automatically outputs a parfile named `[project name].par`. So you'll find an editable `wavetoy.par`.
1. In addition, users can override certain parameter file parameters at the command line. E.g., `wavetoy` has a parameter `convergence_factor` that increases the resolution (technically `Nx=Ny=Nz`) by this factor. To output at twice the resolution, simply run `./wavetoy 2.0`, and a new file will be output `out0d-conv_factor2.00.txt`, which contains data at 2x the resolution.
1. Analyze the output from `out0d-conv_factor1.00.txt` and `out0d-conv_factor2.00.txt` in e.g., `gnuplot`.
