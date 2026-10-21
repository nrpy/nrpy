# Horizon Diagnostics

> Map equation-side BHaHAHA apparent-horizon geometry and spin diagnostics. · Status: confirmed
> Up: [General Relativity](index.md)

## Summary

The `nrpy.equations.general_relativity.bhahaha` modules define symbolic
apparent-horizon geometry and spin diagnostics on surfaces
`r = h(theta, phi)` in a spherical reference-metric chart. This page covers only
equation-side modules: surface area and circumference ingredients, the
expansion function, and SpECTRE-style Omega-based spin estimators.

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
before exporting representative surface-geometry expressions.

`SpECTRESpinEstimateClass` constructs an intrinsic two-surface estimator using
the induced metric, `X_B = e_B^i K_ij s^j`, the surface Ricci scalar `R`, and
the spin function `Omega`. The public API returns integrands through
`get_public_integrands`, gridfunction assignment recipes through
`get_gridfunction_assignments`, and algebraic post-integration helpers for
centroids, spin magnitude, and nominal/fallback spin-vector components. It is
spherical-only, provides the same two cached coordinate keys, and leaves
near-zero branch selection to generated C code. Its scalar public integrands
require `dA`; `area_density` supplies `sqrtq` separately.

Claim evidence:
- Claim: SpECTRE scalar public integrands require an external area measure; `area_density` supplies `sqrtq` separately.
- Role: descriptive behavior
- Deciding authority: [SpECTRESpinEstimate.py](../../../nrpy/equations/general_relativity/bhahaha/SpECTRESpinEstimate.py), `get_public_integrands`
- Corroboration: none available; the expressions' explicit factors decide their integration measure.

The SpECTRE class expects physical metric, extrinsic curvature, and normal
components already evaluated on the surface in the spherical basis. Its
`hh_dD[0]` and `hh_dD[1]` mean theta and phi derivatives; the other horizon
modules use slots 1 and 2. Callers must map these indices explicitly.
`get_gridfunction_assignments` does not generate the required derivatives:
assign the induced metric and `X_B`, compute their angular derivatives and
surface curvature, then construct and differentiate the curvature-dependent
flux where needed, using consistent finite-difference stencils.

`Salpha_and_magnitude_from_zOmega` requires area-weighted integrals of externally
supplied, appropriately normalized modes. It returns
`S_alpha=integral(z_alpha*Omega*dA)/(8*pi)` and their Euclidean norm; this
module neither solves for nor normalizes the modes. Rescaling them rescales
the reported spin.

Claim evidence:
- Claim: SpECTRE spin expressions require on-surface physical tensors, their own two-dimensional angular-index convention, caller-computed derivatives in dependency order, and normalized modes with area-weighted integrals.
- Role: descriptive behavior
- Deciding authority: [SpECTRESpinEstimate.py](../../../nrpy/equations/general_relativity/bhahaha/SpECTRESpinEstimate.py), `SpECTRESpinEstimateClass.__init__`, `get_gridfunction_assignments`, and `Salpha_and_magnitude_from_zOmega`
- Corroboration: [ExpansionFunctionTheta.py](../../../nrpy/equations/general_relativity/bhahaha/ExpansionFunctionTheta.py), `ExpansionFunctionThetaClass`, for the differing three-dimensional derivative slots; none available for a complete numerical spin calculation.

The cited representative trusted dictionaries cover the ordinary spherical
variants for each assigned equation-side family. The source modules also expose
rfm-precompute keys where noted above, but those broader variants are not
enumerated as source evidence here. Trusted files validate symbolic
expressions; they do not document infrastructure generation.

## Sources

- [ExpansionFunctionTheta.py](../../../nrpy/equations/general_relativity/bhahaha/ExpansionFunctionTheta.py) - `ExpansionFunctionThetaClass`, `ExpansionFunctionTheta`
- [area.py](../../../nrpy/equations/general_relativity/bhahaha/area.py) - `area`, `area2`, `area3`, `compute_q2DD`, `circumference_metric_roots`
- [SpECTRESpinEstimate.py](../../../nrpy/equations/general_relativity/bhahaha/SpECTRESpinEstimate.py) - `SpECTRESpinEstimateClass`, `get_public_integrands`, `compute_spin_vectors_for_c`
- [area_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/area_Spherical.py) - `trusted_dict`
- [ExpansionFunctionTheta_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/ExpansionFunctionTheta_Spherical.py) - `trusted_dict`
- [SpECTRESpinEstimate_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/SpECTRESpinEstimate_Spherical.py) - `trusted_dict`

## See Also

- [General Relativity](index.md)
- [BSSN Family](bssn-family.md)
- [Metric Conversions And Matter](metric-conversions-and-matter.md)
- [Initial Data](initial-data.md)
- [Reference Metrics](../../core/reference-metrics.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
