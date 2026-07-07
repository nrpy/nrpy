# Horizon Diagnostics

> Map equation-side BHaHAHA apparent-horizon geometry and spin diagnostics. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [General Relativity](index.md)

## Summary

The `nrpy.equations.general_relativity.bhahaha` modules define symbolic
apparent-horizon geometry and spin diagnostics on surfaces
`r = h(theta, phi)` in a spherical reference-metric chart. This page covers only
equation-side modules: surface area and circumference ingredients, the
expansion function, approximate-Killing-vector spin integrands,
vorticity-dipole spin integrands, and SpECTRE-style Omega-based spin
estimators.

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

`ApproxKillingSpinClass` builds single-pass apparent-horizon integrands for
approximate-Killing-vector spin. Public outputs are `sqrtq`,
`Hmn_integrand`, `Nmn_integrand`, and `Jm_integrand`. The class accepts
`Spherical` and `Spherical_rfm_precompute`, constructs the induced 2-metric and
surface connection, uses the standard `l=1` rotation basis on the sphere, and
leaves the global `1/(8*pi)` angular-momentum prefactor to the caller.

`HorizonSpinVorticityDipoleClass` builds quasilocal spin-vector densities from
the horizon vorticity one-form. It supports `K_input="BSSN"` or `"external"`,
`normal_orientation` values `outward` and `inward`, and `K_sign` values `+1`
and `-1`, while the dictionary accessor exposes the usual spherical and
rfm-precompute keys. Its public `JCart_densU` components include the `1/(8*pi)`
factor; helper methods construct scalar spin density about an axis and
Omega-moment densities from a supplied Omega expression.

`SpECTRESpinEstimateClass` constructs an intrinsic two-surface estimator using
the induced metric, `X_B = e_B^i K_ij s^j`, the surface Ricci scalar `R`, and
the spin function `Omega`. The public API returns integrands through
`get_public_integrands`, gridfunction assignment recipes through
`get_gridfunction_assignments`, and algebraic post-integration helpers for
centroids, spin magnitude, and nominal/fallback spin-vector components. It is
spherical-only, provides the same two cached coordinate keys, and leaves
near-zero branch selection to generated C code.

The cited representative trusted dictionaries cover the ordinary spherical
variants for each assigned equation-side family. The source modules also expose
rfm-precompute keys where noted above, but those broader variants are not
enumerated as source evidence here. Trusted files validate symbolic
expressions; they do not document infrastructure generation.

## Sources

- [ExpansionFunctionTheta.py](../../../nrpy/equations/general_relativity/bhahaha/ExpansionFunctionTheta.py) - `ExpansionFunctionThetaClass`, `ExpansionFunctionTheta`
- [area.py](../../../nrpy/equations/general_relativity/bhahaha/area.py) - `area`, `area2`, `area3`, `compute_q2DD`, `circumference_metric_roots`
- [approx_killing_vector_spin.py](../../../nrpy/equations/general_relativity/bhahaha/approx_killing_vector_spin.py) - `ApproxKillingSpinClass`, `ApproxKillingSpin`
- [HorizonSpinVorticityDipole.py](../../../nrpy/equations/general_relativity/bhahaha/HorizonSpinVorticityDipole.py) - `HorizonSpinVorticityDipoleClass`, `JCart_densU`
- [SpECTRESpinEstimate.py](../../../nrpy/equations/general_relativity/bhahaha/SpECTRESpinEstimate.py) - `SpECTRESpinEstimateClass`, `get_public_integrands`, `compute_spin_vectors_for_c`
- [area_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/area_Spherical.py) - `trusted_dict`
- [ExpansionFunctionTheta_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/ExpansionFunctionTheta_Spherical.py) - `trusted_dict`
- [approx_killing_vector_spin_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/approx_killing_vector_spin_Spherical.py) - `trusted_dict`
- [HorizonSpinVorticityDipole_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/HorizonSpinVorticityDipole_Spherical.py) - `trusted_dict`
- [SpECTRESpinEstimate_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/SpECTRESpinEstimate_Spherical.py) - `trusted_dict`

## See Also

- [General Relativity](index.md)
- [BSSN Family](bssn-family.md)
- [Metric Conversions And Matter](metric-conversions-and-matter.md)
- [Initial Data](initial-data.md)
- [Reference Metrics](../../core/reference-metrics.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
