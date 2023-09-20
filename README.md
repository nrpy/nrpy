[![Python CI](https://github.com/zachetienne/nrpynextgen/actions/workflows/master.yml/badge.svg)](https://github.com/zachetienne/nrpynextgen/actions/workflows/master.yml)

# Next-generation NRPy+

Instructions for generating the wavetoy example:

1. Set up a Python virtual environment: `python3 -m venv ~/virt`
2. Activate the virtual environment: `source ~/virt/bin/activate`
1. Navigate to the `nrpynextgen` root directory.
1. Install requirements: `pip install -r requirements.txt`
1. Set up your environment paths: `export PYTHONPATH=.:$PYTHONPATH`
1. Generate the wave equation C code: `python3 s2f/wave_equation_cartesian.py` . It will output a full C code to `projects/wavetoy`
1. Navigate to `projects/wavetoy` and run `make` . It will compile the wavetoy C code.
1. Run `./wavetoy` . It'll output at various times the solution from the point closest to r=0, as well as the relative error compared to the exact solution at that point to a file `out0d-conv_factor1.00.txt` .
1. All BlackHoles@Home (BHaH) infrastructure codes now support parameter files, and each code generation automatically outputs a parfile named `[project name].par`. So you'll find an editable `wavetoy.par`.
1. In addition, users can override certain parameter file parameters at the command line. E.g., `wavetoy` has a parameter `convergence_factor` that increases the resolution (technically `Nx=Ny=Nz`) by this factor. To output at twice the resolution, simply run `./wavetoy 2.0`, and a new file will be output `out0d-conv_factor2.00.txt`, which contains data at 2x the resolution.
1. Analyze the output from `out0d-conv_factor1.00.txt` and `out0d-conv_factor2.00.txt` in e.g., `gnuplot`.
