# Horizon Diagnostics

> Map equation-side BHaHAHA apparent-horizon geometry and spin diagnostics. · Status: confirmed · Last reconciled: 08-10-2026
> Up: [General Relativity](index.md)

## Summary

The `nrpy.equations.general_relativity.bhahaha` modules define symbolic
apparent-horizon geometry and spin diagnostics on surfaces
`r = h(theta, phi)` in a spherical reference-metric chart. This page covers the
equation-side modules and the comparison branch's runtime AKV-pole contract:
surface area and circumference ingredients, the
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

### Lowest-shear AKV-pole spin

The comparison branch's AKV-pole diagnostic uses the positive-definite
squared-Killing-residual generalized eigenproblem `K z = lambda M z`, with
constants and the centered-phi Nyquist checkerboard removed in the reduced
space. This is the Owen/SpEC-line variation of the AKV construction used by
the runtime branch, not the original curvature-weighted Cook--Whiting
eigenproblem. The three returned eigenpairs are put in ascending algebraic
eigenvalue order before expansion and normalization, so `zU0` is the strict
lowest-eigenvalue, lowest-shear potential rather than the mode with largest
angular momentum or a Procrustes-rotated combination. Near degeneracy is not
an error; Procrustes-aligned combinations remain seeds only for later solves.

Every eigenpotential retains the existing integral normalization

\[
\oint z_\alpha^2\,dA=\frac{A^3}{48\pi^2}.
\]

In particular, the lowest-mode signed angular momentum and its magnitude are

\[
J_0=\frac{1}{8\pi}\oint z_0\Omega\,dA
   =\frac{\mathrm{ZOU}[0]}{8\pi},
\qquad |J_0|=\lVert S_{\mathrm{AKV}}\rVert.
\]

The extrema determine only a coordinate direction and do not renormalize
`zU0`. Each coarse global extremum is refined from unique points in a reflected
`5 x 5` scalar stencil. The runtime maps that neighborhood to a nonsingular
Cartesian tangent basis, fits quadratics to both `zU0` and the horizon radius,
requires the appropriate definite Hessian and a local trust region, and
evaluates both refined horizon points at the same fitted angular locations.
If `x_max` and `x_min` denote those background-Cartesian points, then

\[
\hat n_{\mathrm{AKV}}^i=
\frac{x_{\max}^i-x_{\min}^i}{\lVert x_{\max}-x_{\min}\rVert},
\qquad
S_{\mathrm{AKV}}^i=J_0\hat n_{\mathrm{AKV}}^i.
\]

The diagnostic uses its own Christodoulou mass rather than either three-mode
magnitude:

\[
M_{\mathrm{irr}}^2=\frac{A}{16\pi},\qquad
M_{\mathrm{AKV}}^2=M_{\mathrm{irr}}^2+
\frac{J_0^2}{4M_{\mathrm{irr}}^2},\qquad
\chi_{\mathrm{AKV}}^i=\frac{S_{\mathrm{AKV}}^i}{M_{\mathrm{AKV}}^2}.
\]

An eigensolver sign change `z0 -> -z0` reverses `ZOU[0]` and exchanges the
maximum and minimum, so both signs cancel in the full spin vector. The AKV
operator and potential are intrinsic to the horizon, but the reported
three-dimensional pole direction depends on the background Cartesian
coordinates. At exactly zero spin the vector is zero and a normalized physical
direction is undefined. Failure of either pole fit invalidates only this
AKV-pole vector; it does not invalidate the existing SpECTRE-style or direct
Gram-matrix comparison vectors.

Claim evidence:
- Claim: The comparison branch defines the AKV-pole vector from the signed integral of the strict lowest-eigenvalue normalized AKV potential and its refined background-Cartesian extrema, uses an independent Christodoulou mass, preserves eigenvector-sign invariance, and makes pole-fit failure local to this vector.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/BHaH/BHaHAHA/spectre_spin_integrator.py` `bah_compute_spectre_spin_potentials`, `diagnostics_spectre_spin`
- Corroboration: `nrpy/equations/general_relativity/bhahaha/SpECTRESpinEstimate.py` `get_public_integrands`; `nrpy/infrastructures/BHaH/BHaHAHA/BHaHAHA_header.h` `bhahaha_diagnostics_struct`
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-run; tool_version=not-run; backend=not-run; precision=not-run; GPU=not-run; restart=not-run; distributed=not-run; error_path=not-run; options=not-run; date=08-10-2026`

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
- [spectre_spin_integrator.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/spectre_spin_integrator.py) - `bah_compute_spectre_spin_potentials`, `diagnostics_spectre_spin`
- [BHaHAHA_header.h](../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaHAHA_header.h) - `bhahaha_diagnostics_struct`
- [area_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/area_Spherical.py) - `trusted_dict`
- [ExpansionFunctionTheta_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/ExpansionFunctionTheta_Spherical.py) - `trusted_dict`
- [approx_killing_vector_spin_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/approx_killing_vector_spin_Spherical.py) - `trusted_dict`
- [HorizonSpinVorticityDipole_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/HorizonSpinVorticityDipole_Spherical.py) - `trusted_dict`
- [SpECTRESpinEstimate_Spherical.py](../../../nrpy/equations/general_relativity/bhahaha/tests/SpECTRESpinEstimate_Spherical.py) - `trusted_dict`

## See Also

- Parent: [General Relativity](index.md)
- See also: [BSSN Family](bssn-family.md)
- See also: [Metric Conversions And Matter](metric-conversions-and-matter.md)
- See also: [Initial Data](initial-data.md)
- Depends on: [Reference Metrics](../../core/reference-metrics.md)
- Validated by: [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
