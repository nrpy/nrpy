# Horizon Diagnostics

> Map equation-side BHaHAHA apparent-horizon geometry and spin diagnostics. · Status: confirmed
> Up: [General Relativity](index.md)

## Summary

The `nrpy.equations.general_relativity.bhahaha` modules define symbolic
apparent-horizon geometry and spin diagnostics on horizon surfaces
`r = h(theta, phi)` in a spherical reference-metric chart: two-surface area and
circumference integrands, the expansion function, and SpECTRE-style Omega-based
spin estimators.

## Detail

`ExpansionFunctionThetaClass` computes the apparent-horizon expansion
`Theta = D_i s^i - K + s^i s^j K_ij`. It accepts only `CoordSystem="Spherical"`
and supports the cached keys `Spherical` and `Spherical_rfm_precompute`. The
construction declares horizon-shape derivatives, BSSN-style metric and
extrinsic-curvature inputs, builds the physical metric `gammaDD`, inverse
metric, determinant derivatives, unit normal `sU`, `KDD`, and stores `Theta`.

`area.py` reuses the spherical expansion geometry to provide `area`, `area2`,
`area3`, `compute_q2DD`, `circumferential_arclength`, and
`circumference_metric_roots`. It also provides `spin_NewtonRaphson` and
`spin_HalleysMethod` iteration expressions for circumference-ratio spin
estimation. Its script validation checks that the three area expressions agree
before exporting representative two-surface geometry expressions.

`SpECTRESpinEstimateClass` constructs an intrinsic two-surface estimator using
the induced metric, `X_B = e_B^i K_ij s^j`, the two-surface Ricci scalar `R`, and
the spin function `Omega`. The public API returns integrands through
`get_public_integrands`, gridfunction assignment recipes through
`get_gridfunction_assignments`, and algebraic post-integration helpers for
centroids, spin magnitude, and nominal/fallback spin-vector components. It is
spherical-only, provides the same two cached coordinate keys, and leaves
near-zero branch selection to generated C code. Its scalar public integrands
require `dA`; `area_density` supplies `sqrtq` separately.

Evidence (descriptive behavior): `general_relativity/bhahaha/SpECTRESpinEstimate.py::get_public_integrands`.

The SpECTRE class obtains the physical metric, extrinsic curvature, and unit
normal in the spherical basis from `ExpansionFunctionTheta` and evaluates them
on the horizon by substituting `r = h(theta, phi)`. It uses `hh_dD[1]` for the
theta derivative and `hh_dD[2]` for the phi derivative, as in the expansion
function. The generated spin integrator assigns the induced metric and `X_B`,
computes their angular derivatives and the two-surface curvature, then assigns
and differentiates the curvature-dependent flux.

`Salpha_and_magnitude_from_zOmega` requires area-weighted integrals of externally
supplied, appropriately normalized modes. It returns
`S_alpha=integral(z_alpha*Omega*dA)/(8*pi)` and their Euclidean norm; this
module neither solves for nor normalizes the modes. Rescaling them rescales
the reported spin.

Evidence (descriptive behavior): `general_relativity/bhahaha/SpECTRESpinEstimate.py::SpECTRESpinEstimateClass.__init__`, `get_gridfunction_assignments`, and `Salpha_and_magnitude_from_zOmega`; none available for a complete numerical spin calculation.

The cited representative trusted dictionaries cover the ordinary spherical
variants for each assigned equation-side family. The source modules also expose
rfm-precompute keys where noted above, but those broader variants are not
enumerated as source evidence here. Trusted files validate symbolic
expressions; they do not document infrastructure generation.

### Lowest-shear AKV poles

The generated BHaHAHA spin integrator solves the squared Killing-residual
eigenproblem `K z = lambda M z` with the constant mode and centered-phi Nyquist
mode removed. It orders the three eigenpairs by increasing eigenvalue and
normalizes each potential by `integral(z_alpha^2 dA) = A^3/(48*pi^2)`. It uses
the lowest-shear potential `zU0` to compute the signed angular momentum
`J0 = integral(zU0*Omega dA)/(8*pi)` and refines the maximum and minimum with
quadratic fits to the potential and horizon radius on reflected
`5 × 5` stencils. The fits use local Cartesian tangent coordinates and
require the appropriate definite Hessian and a displacement within the
local trust region. Procrustes-aligned modes seed later solves. The measured
poles come from the lowest eigenmode of the current solve.
The direction from the minimum to the maximum is measured in background
Cartesian coordinates. The dimensionless vector is
`J0` times this unit direction divided by
`M_irr^2 + J0^2/(4*M_irr^2)`, where `M_irr^2 = A/(16*pi)`.
Changing the sign of `zU0` reverses `J0` and swaps the poles, preserving the
vector. At zero spin the vector vanishes and there is no physical spin axis.

Evidence (descriptive behavior): [spectre_spin_integrator.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/spectre_spin_integrator.py) `bah_compute_spectre_spin_potentials` and `diagnostics_spectre_spin`.

## Sources

- [ExpansionFunctionTheta.py](../../../nrpy/equations/general_relativity/bhahaha/ExpansionFunctionTheta.py) - `ExpansionFunctionThetaClass`, `ExpansionFunctionTheta`
- [area.py](../../../nrpy/equations/general_relativity/bhahaha/area.py) - `area`, `area2`, `area3`, `compute_q2DD`, `circumference_metric_roots`
- [SpECTRESpinEstimate.py](../../../nrpy/equations/general_relativity/bhahaha/SpECTRESpinEstimate.py) - `SpECTRESpinEstimateClass`, `get_public_integrands`, `compute_spin_vectors_for_c`
- [area_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/area_Spherical.py) - `trusted_dict`
- [ExpansionFunctionTheta_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/ExpansionFunctionTheta_Spherical.py) - `trusted_dict`
- [SpECTRESpinEstimate_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/SpECTRESpinEstimate_Spherical.py) - `trusted_dict`

- [spectre_spin_integrator.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/spectre_spin_integrator.py) - `bah_compute_spectre_spin_potentials`, `diagnostics_spectre_spin`

## See Also

- See also: [BSSN Family](bssn-family.md)
- Depends on: [Reference Metrics](../../core/reference-metrics.md)
- Validated by: [Expression Validation Helpers](../../validation/expression-validation-helpers.md)
