# Wave Equation

> Map the scalar wave-equation RHS and exact-solution modules. · Status: contested
> Up: [Equations](index.md)

## Summary

The wave-equation modules provide a compact first-order-in-time scalar wave
system for code generation. Cartesian and reference-metric curvilinear RHS
builders own the evolution expressions, while the initial-data module owns
Cartesian exact solutions for plane-wave and spherical-Gaussian tests.

## Detail

`WaveEquation_RHSs` constructs the Cartesian system using the evolution
gridfunctions `uu` and `vv`. It requests the common `wavespeed` code
parameter subject to the registry limitation below, treats `vv` as `uu_rhs`,
and sets `vv_rhs` to the Cartesian Laplacian of `uu` multiplied by
`wavespeed**2`.

`WaveEquationCurvilinear_RHSs` generalizes the RHS path to any named NRPy
reference metric accepted by the constructor. It selects either the ordinary
reference metric or its `_rfm_precompute` variant, contracts
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

Claim evidence:
- Claim: Wave constructors guard paired gridfunction registration with `uu` alone; Cartesian RHS and initial-data parameter guards inspect the NRPy-parameter dictionary instead of the CodeParameter dictionary.
- Role: descriptive behavior
- Deciding authority: [WaveEquation_RHSs.py](../../nrpy/equations/wave_equation/WaveEquation_RHSs.py), `WaveEquation_RHSs.__init__`; [WaveEquationCurvilinear_RHSs.py](../../nrpy/equations/wave_equation/WaveEquationCurvilinear_RHSs.py), `WaveEquationCurvilinear_RHSs.__init__`; [WaveEquation_Solutions_InitialData.py](../../nrpy/equations/wave_equation/WaveEquation_Solutions_InitialData.py), `WaveEquation_solution_Cartesian.__init__`
- Corroboration: [params.py](../../nrpy/params.py), `CodeParameter._register_globally`, stores CodeParameters in `glb_code_params_dict`.

`WaveEquation_solution_Cartesian` owns the exact-solution and initial-data
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

Claim evidence:
- Claim: Exact wave solutions require nonzero wave-vector norm or Gaussian width; the Gaussian origin requires explicit use of the separate origin expressions.
- Role: descriptive behavior
- Deciding authority: [WaveEquation_Solutions_InitialData.py](../../nrpy/equations/wave_equation/WaveEquation_Solutions_InitialData.py), `PlaneWave, SphericalGaussian`
- Corroboration: `none available`; the named implementation decides the stated behavior.

Validation follows the trusted-expression pipeline in each module's
`__main__` block. The Cartesian RHS has a single trusted variant; the
curvilinear RHS checks Cartesian, Spherical, SinhSpherical, SinhCylindrical,
and SinhSymTP variants, all with `enable_rfm_precompute=False`; the
initial-data module checks both `PlaneWave` and `SphericalGaussian` exact
solutions. Thus the implemented precompute branch has no wave-module trusted
variant in this entry point.

These trusted files sample symbolic outputs at one deterministic free-symbol
assignment. They do not prove that the continuous expressions satisfy the PDE
for all coordinates, that finite-difference generated code builds or runs, or
that a numerical evolution converges. Build/run and numerical-result evidence
belongs to [First Wave Equation Run](../examples/first-wave-equation-run.md) and
the relevant infrastructure and validation owners.

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

- Parent: [Equations](index.md)
- See also: [Conformally Flat Elliptic](conformally-flat-elliptic.md)
- See also: [Geometry And Special-Function Support](geometry-and-special-function-support.md)
- See also: [Trusted Expression Pipeline](trusted-expression-pipeline.md)
- See also: [First Wave Equation Run](../examples/first-wave-equation-run.md)
- See also: [Reference Metrics](../core/reference-metrics.md)
