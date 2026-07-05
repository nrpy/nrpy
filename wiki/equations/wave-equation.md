# Wave Equation

> Map the scalar wave-equation RHS and exact-solution modules. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [Equations](index.md)

## Summary

The wave-equation modules provide a compact first-order-in-time scalar wave
system for code generation. Cartesian and reference-metric curvilinear RHS
builders own the evolution expressions, while the initial-data module owns
Cartesian exact solutions for plane-wave and spherical-Gaussian tests.

## Detail

`WaveEquation_RHSs` constructs the Cartesian system using the evolution
gridfunctions `uu` and `vv`. It registers the common `wavespeed` code
parameter, treats `vv` as `uu_rhs`, and sets `vv_rhs` to the Cartesian Laplacian
of `uu` multiplied by `wavespeed**2`.

`WaveEquationCurvilinear_RHSs` generalizes the RHS path to any named NRPy
reference metric accepted by the constructor. It selects either the ordinary
reference metric or its `_rfm_precompute` variant, contracts
`GammahatUDD` with `ghatUU`, declares `uu_dD` and symmetric `uu_dDD`, and forms
`vv_rhs` from the reference-metric Laplacian and contracted-Christoffel term
before multiplying by `wavespeed**2`.

Both RHS classes register `uu` and `vv` as `EVOL` gridfunctions when they are
not already present. The registered asymptotic values are `f_infinity=[2.0,
0.0]`, matching the exact-solution convention that `uu` is shifted away from
zero so relative-error diagnostics are well defined.

`WaveEquation_solution_Cartesian` owns the exact-solution and initial-data
expressions. It registers `wavespeed`, registers `time` either in the local
module or in the BHaH Method-of-Lines module when `Infrastructure` is `BHaH`,
and dispatches by `WaveType`. `PlaneWave` registers `kk0`, `kk1`, and `kk2`,
normalizes that wave vector, and stores `uu_exactsoln` and `vv_exactsoln`.
`SphericalGaussian` registers `sigma`, builds the outgoing-plus-ingoing
spherical Gaussian, differentiates it for `vv_exactsoln`, and also stores
origin-safe `uu_exactsoln_r0` and `vv_exactsoln_r0` expressions.

Validation follows the trusted-expression pipeline in each module's
`__main__` block. The Cartesian RHS has a single trusted variant; the
curvilinear RHS checks Cartesian, Spherical, SinhSpherical, SinhCylindrical,
and SinhSymTP variants; the initial-data module checks both `PlaneWave` and
`SphericalGaussian` exact solutions.

## Sources

- [WaveEquation_RHSs.py](../../nrpy/equations/wave_equation/WaveEquation_RHSs.py) - `WaveEquation_RHSs`, `uu_rhs`, `vv_rhs`
- [WaveEquationCurvilinear_RHSs.py](../../nrpy/equations/wave_equation/WaveEquationCurvilinear_RHSs.py) - `WaveEquationCurvilinear_RHSs`, `uu_rhs`, `vv_rhs`
- [WaveEquation_Solutions_InitialData.py](../../nrpy/equations/wave_equation/WaveEquation_Solutions_InitialData.py) - `WaveEquation_solution_Cartesian`, `PlaneWave`, `SphericalGaussian`
- [WaveEquation_RHSs_WaveEquation.py](../../nrpy/equations/wave_equation/tests/WaveEquation_RHSs_WaveEquation.py) - `trusted_dict`
- [WaveEquationCurvilinear_RHSs_Spherical.py](../../nrpy/equations/wave_equation/tests/WaveEquationCurvilinear_RHSs_Spherical.py) - `trusted_dict`
- [WaveEquation_Solutions_InitialData_PlaneWave.py](../../nrpy/equations/wave_equation/tests/WaveEquation_Solutions_InitialData_PlaneWave.py) - `trusted_dict`
- [WaveEquation_Solutions_InitialData_SphericalGaussian.py](../../nrpy/equations/wave_equation/tests/WaveEquation_Solutions_InitialData_SphericalGaussian.py) - `trusted_dict`

## See Also

- [Equations](index.md)
- [Conformally Flat Elliptic](conformally-flat-elliptic.md)
- [Geometry And Special-Function Support](geometry-and-special-function-support.md)
- [Trusted Expression Pipeline](trusted-expression-pipeline.md)
- [First Wave Equation Run](../examples/first-wave-equation-run.md)
- [Reference Metrics](../core/reference-metrics.md)
