# First Wave Equation Run

> Describe the first standalone project generation, build, run, and output milestone. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [Examples](index.md)

## Summary

The recommended first run is `python -m nrpy.examples.wave_equation_cartesian`.
It needs Python, a C compiler, and `make`; it writes
`project/wave_equation_cartesian/`, which can then be built with `make` and run
as `./wave_equation_cartesian`. The generator deletes and recreates that fixed
project directory, so preserve any wanted prior output first.

## Detail

For a package install, run the generator with:

```bash
python -m nrpy.examples.wave_equation_cartesian
```

From a repository checkout without an editable install, first append `.` to
`PYTHONPATH` from the repository root. After generation, build and run:

```bash
cd project/wave_equation_cartesian
make
./wave_equation_cartesian
```

The generator creates a parameter file named `wave_equation_cartesian.par`.
The important first milestone is that project generation succeeds, `make`
succeeds, and the executable runs without manual edits inside the generated
project. This is a user-run milestone, not a claim reproduced during this KB
audit. Configured Ubuntu/macOS CI generates and builds the default project but
does not run the executable or inspect its numerical output. The README also
documents a simple command-line override:
`./wave_equation_cartesian 2.0`, which uses `convergence_factor` and writes
output such as `out0d-conv_factor2.00.txt`.

The example source sets `project_name = "wave_equation_cartesian"`, uses the
BHaH infrastructure, chooses double precision, a `SphericalGaussian` wave, RK4
Method of Lines timestepping, fourth-order finite differences, and a Cartesian
cell-centered grid. It registers exact-solution, initial-data, numerical-grid,
diagnostic, RHS, and boundary-condition C functions, then writes code-parameter
headers, a default parfile, `BHaH_defines.h`, the generated `main`, cleanup
functions, copied SIMD intrinsics when enabled, function prototypes, and the
Makefile.

The RHS path is the compact symbolic-codegen pattern: the example builds wave
equation SymPy RHS expressions, passes them through `ccg.c_codegen`, and wraps
the emitted assignments in `BHaH.simple_loop.simple_loop` over the interior.
For the broader Cartesian, curvilinear, and multicoordinate wave-equation
family, use [Wave Equation Generators](wave-equation-generators.md).

## Sources

- [README.md](../../README.md) - `Installation`, `First Successful Run`
- [README.md](../../README.md) - `Prerequisites by Workflow`
- [wave_equation_cartesian.py](../../nrpy/examples/wave_equation_cartesian.py) - `project_name`, `register_CFunction_rhs_eval`
- [wave_equation_cartesian.py](../../nrpy/examples/wave_equation_cartesian.py) - `register_CFunction_initial_data`, `register_CFunction_apply_bcs`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `Required Checks`

## See Also

- [Examples](index.md)
- [Wave Equation Generators](wave-equation-generators.md)
- [Black Hole Evolution](black-hole-evolution.md)
- [Lifecycle And Project Assembly](../infrastructures/bhah/lifecycle-and-project-assembly.md)
