# Wave Equation

> Map the scalar wave-equation RHS and exact-solution modules. · Status: contested
> Up: [Equations](index.md)

## Summary

The wave-equation modules provide a compact first-order-in-time scalar wave
system for code generation. Cartesian and reference-metric curvilinear RHS
builders construct the evolution expressions, while the initial-data module constructs
Cartesian exact solutions for plane-wave and spherical-Gaussian tests.

## Detail

`WaveEquation_RHSs` constructs the Cartesian system using the evolution
gridfunctions `uu` and `vv`. It requests the common `wavespeed` code
parameter subject to the registry limitation below, treats `vv` as `uu_rhs`,
and sets `vv_rhs` to the Cartesian Laplacian of `uu` multiplied by
`wavespeed**2`.

`WaveEquationCurvilinear_RHSs` generalizes the RHS path to any named NRPy
reference metric accepted by the constructor. It selects either the ordinary
reference metric or its `_rfm_precompute` variant, forms the tensor contraction of
`GammahatUDD` with `ghatUU`, declares `uu_dD` and symmetric `uu_dDD`, and forms
`vv_rhs` from the reference-metric Laplacian and contracted-Christoffel term
before multiplying by `wavespeed**2`.

Both RHS classes and the initial-data constructor test only whether `uu` is
registered. When it is absent, they register `uu` and `vv` together as `EVOL`
with `f_infinity=[2.0, 0.0]`. When `uu` is already present, they neither ensure
that `vv` exists nor check the existing groups and asymptotic values. Callers
reusing registrations must ensure that both compatible `EVOL` entries exist.
The new-registration values match the exact solutions' offset of `uu` from zero.

The Cartesian RHS and initial-data constructors test `par.glb_params_dict`
before registering `wavespeed`; initial data uses that dictionary for `time`
as well. CodeParameters are stored in `par.glb_code_params_dict`, so a same-name
NRPy parameter can suppress required CodeParameter registration. The
curvilinear RHS calls `register_CodeParameter` for `wavespeed` unconditionally.

Claim status: contested; contradiction: CONTR-0021.
See [CONTR-0021](../contradictions.md#contr-0021).

Evidence (descriptive behavior): `wave_equation/WaveEquation_RHSs.py::WaveEquation_RHSs.__init__`; `wave_equation/WaveEquationCurvilinear_RHSs.py::WaveEquationCurvilinear_RHSs.__init__`; `wave_equation/WaveEquation_Solutions_InitialData.py::WaveEquation_solution_Cartesian.__init__`.

`WaveEquation_solution_Cartesian` constructs the exact-solution and initial-data
expressions. Subject to those guards, it registers `wavespeed` and `time`;
`time` uses the local module or the BHaH Method-of-Lines module when
`Infrastructure` is `BHaH`. The constructor dispatches by `WaveType`. `PlaneWave` registers `kk0`, `kk1`, and `kk2`,
normalizes that wave vector, and stores `uu_exactsoln` and `vv_exactsoln`.
`SphericalGaussian` registers `sigma`, builds the outgoing-plus-ingoing
spherical Gaussian, differentiates it for `vv_exactsoln`, and also stores
separate `uu_exactsoln_r0` and `vv_exactsoln_r0` expressions for evaluation at
the origin. The plane-wave vector must have nonzero norm, and the Gaussian
width `sigma` must be nonzero. The primary Gaussian expressions divide by
`r`; the caller must explicitly select the separate origin expressions at
`r=0`. Those expressions are not equivalent replacements for the primary
solution away from the origin.

Evidence (descriptive behavior): `wave_equation/WaveEquation_Solutions_InitialData.py::WaveEquation_solution_Cartesian`.

Validation compares trusted-expression dictionaries in each module's
`__main__` block. The Cartesian RHS has a single trusted variant; the
curvilinear RHS checks Cartesian, Spherical, SinhSpherical, SinhCylindrical,
and SinhSymTP variants, all with `enable_rfm_precompute=False`; the
initial-data module checks both `PlaneWave` and `SphericalGaussian` exact
solutions. Thus the implemented precompute branch has no wave-module trusted
variant in this entry point; for standalone C/CUDA build and execution see
[Wave Equation Generators](../examples/wave-equation-generators.md).

## Sources

- [params.py](../../nrpy/params.py) - `CodeParameter._register_globally`
- [WaveEquation_RHSs.py](../../nrpy/equations/wave_equation/WaveEquation_RHSs.py) - `WaveEquation_RHSs`, `uu_rhs`, `vv_rhs`
- [WaveEquationCurvilinear_RHSs.py](../../nrpy/equations/wave_equation/WaveEquationCurvilinear_RHSs.py) - `WaveEquationCurvilinear_RHSs`, `uu_rhs`, `vv_rhs`
- [WaveEquation_Solutions_InitialData.py](../../nrpy/equations/wave_equation/WaveEquation_Solutions_InitialData.py) - `WaveEquation_solution_Cartesian`, `PlaneWave`, `SphericalGaussian`
- [WaveEquation_RHSs_WaveEquation.py](../../nrpy/equations/wave_equation/tests/WaveEquation_RHSs_WaveEquation.py) - `trusted_dict`
- [WaveEquationCurvilinear_RHSs_Spherical.py](../../nrpy/equations/wave_equation/tests/WaveEquationCurvilinear_RHSs_Spherical.py) - `trusted_dict`
- [WaveEquation_Solutions_InitialData_PlaneWave.py](../../nrpy/equations/wave_equation/tests/WaveEquation_Solutions_InitialData_PlaneWave.py) - `trusted_dict`
- [WaveEquation_Solutions_InitialData_SphericalGaussian.py](../../nrpy/equations/wave_equation/tests/WaveEquation_Solutions_InitialData_SphericalGaussian.py) - `trusted_dict`

## See Also

- See also: [Conformally Flat Elliptic](conformally-flat-elliptic.md)
- See also: [Geometry And Special-Function Support](geometry-and-special-function-support.md)
- See also: [Expression Validation Helpers](../validation/expression-validation-helpers.md)
