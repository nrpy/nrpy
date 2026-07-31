# Geodesics

> Map analytic spacetime metrics, geodesic RHS construction, and conserved diagnostics. · Status: confirmed · Last reconciled: 07-30-2026
> Up: [General Relativity](index.md)

## Summary

The geodesics package builds symbolic metric, Christoffel, equation-of-motion,
constraint, and diagnostic expressions for particle paths. Its analytic
registry currently implements `KerrSchild_Cartesian`; `GeodesicEquations`
supports massive and photon RHS variants; `GeodesicDiagnostics` provides
energy, Cartesian angular momentum, and Kerr-Schild Carter-constant diagnostics
where the required symmetries are available.

## Detail

`AnalyticSpacetimes` dispatches by spacetime name. The implemented branch is
`KerrSchild_Cartesian`, which registers `M_scale` and `a_spin`, defines
coordinates `[t, x, y, z]`, and stores the covariant four-metric `g4DD`.
Unsupported analytic spacetime names raise `ValueError`; the lazy
`Analytic_Spacetimes` dictionary caches constructed instances.

`GeodesicEquations` obtains `g4DD` and `xx` from the analytic registry, computes
`g4DD_dD`, builds `Gamma4UDD`, and then chooses the RHS and initialization
constraint by `particle_type`. The massive path returns eight ODE right-hand
sides from position and four-velocity evolution, plus `u0_massive` from the
timelike normalization condition. The photon path returns nine ODE right-hand
sides from position and four-momentum evolution; the ninth component is the
normal-observer path-length diagnostic, and `p0_photon` comes from the null constraint.
For photon evolution, the module provides four explicit RHS forms: ordinary
`geodesic_eom_rhs_photon` and
`geodesic_eom_rhs_photon_christoffel`, plus normalized
`geodesic_eom_rhs_photon_normalized` and
`geodesic_eom_rhs_photon_normalized_christoffel`. The first member of each pair
consumes metric derivatives; the second consumes Christoffels. EOM-family
selection and geometry-payload selection are separate generation choices.
Unsupported particle types raise `ValueError`.

The photon normalization helpers make the direct and normalized state
conventions explicit. They define `u = ln|alpha p^0|` and
`Pi_i = p_i/(alpha p^0)`, so the normalized null constraint is
`gamma^ij Pi_i Pi_j = 1`; the direct constraint is
`g_mu_nu p^mu p^nu = 0`. `normal_observer_log_energy` evaluates the same
normal-observer log-energy measure for direct states by constructing the needed
inverse-metric component from the covariant `g4DD`. The measure is an
upper-only evolution cutoff, not a claim that energy is conserved in a general
spacetime.

Claim evidence:
- Claim: Photon direct and normalized states use the documented null constraints, and direct states can be mapped to the common normal-observer log-energy measure `ln|alpha p^0|`.
- Role: public/scientific contract
- Deciding authority: `nrpy/equations/general_relativity/geodesics/geodesics.py` — `photon_momentum_to_normalized_quantities`, `normalization_constraint_photon`, `normalization_constraint_photon_normalized`
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/normal_observer_log_energy.py` — `normal_observer_log_energy`
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-30-2026`

Claim evidence:
- Claim: `GeodesicEquations` exposes direct and normalized photon RHS forms for metric-derivative and Christoffel payloads, with EOM-family and geometry-payload selection independent.
- Role: public/scientific contract
- Deciding authority: `nrpy/equations/general_relativity/geodesics/geodesics.py` — `GeodesicEquations`
- Corroboration: `nrpy/equations/general_relativity/geodesics/geodesics.py` — photon RHS construction and payload selection branches
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

The same module also exposes reusable recipes for numerical spacetime data:
`symbolic_numerical_christoffel_recipe`,
`symbolic_g4DD_recipe_from_bssn_grid_basis`,
`symbolic_g4DD_dt_recipe_from_bssn_grid_basis`, and
`symbolic_christoffel_recipe_from_bssn_grid_basis`. The BSSN-grid helpers use
reference-metric Jacobians and assume a time-independent spatial map lifted
into four dimensions; the Christoffel recipe can optionally use static time
derivatives for an explicitly requested approximate endpoint payload.

Claim evidence:
- Claim: The numerical-spacetime recipes use the BSSN-grid reference-metric transformation and can request static time derivatives only for an explicitly selected approximate endpoint payload.
- Role: public/scientific contract
- Deciding authority: `nrpy/equations/general_relativity/geodesics/geodesics.py` — `symbolic_g4DD_dt_recipe_from_bssn_grid_basis`, `symbolic_christoffel_recipe_from_bssn_grid_basis`
- Corroboration: `nrpy/infrastructures/BHaH/diagnostics/output_raytracing_data.py` — mode-selected exporter payload
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

`GeodesicDiagnostics` can use analytic metrics or a generic `"Numerical"`
metric placeholder. It always constructs `E_expr = -p_0`, builds scalar `Lz_expr` only
for spacetime names ending in `Cartesian`, and builds `Q_expr` only for
`KerrSchild_Cartesian`. The Carter-constant path distinguishes massive and
photon particles and relies on the module assumptions documented in
`conserved_quantities.py`.

Claim evidence:
- Claim: `GeodesicDiagnostics` exposes scalar `Lz_expr` and the axial method `compute_angular_momentum_z_cartesian()` for Cartesian spacetime names; the removed vector `L_exprs` and old Cartesian method are not part of the current API.
- Role: public/scientific contract
- Deciding authority: `nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/conserved_quantities.py` — `GeodesicDiagnostics`
- Corroboration: `nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/tests/conserved_quantities_KerrSchild_Cartesian_massive.py` and `conserved_quantities_KerrSchild_Cartesian_photon.py` — trusted diagnostics
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

Representative trusted dictionaries cover the analytic metric, massive and
photon geodesic equations, and massive and photon conserved diagnostics for
`KerrSchild_Cartesian`. The module-level validation also checks metric and
Christoffel symmetry plus identity-map and spherical-to-Cartesian BSSN transform
cases before emitting trusted expressions.

## Sources

- [analytic_spacetimes.py](../../../nrpy/equations/general_relativity/geodesics/analytic_spacetimes.py) - `AnalyticSpacetimes`, `Analytic_Spacetimes`
- [geodesics.py](../../../nrpy/equations/general_relativity/geodesics/geodesics.py) - `GeodesicEquations`, `geodesic_eom_rhs_massive`, `geodesic_eom_rhs_photon`, `geodesic_eom_rhs_photon_christoffel`, `geodesic_eom_rhs_photon_normalized`, `geodesic_eom_rhs_photon_normalized_christoffel`, `symbolic_g4DD_dt_recipe_from_bssn_grid_basis`, `symbolic_christoffel_recipe_from_bssn_grid_basis`
- [conserved_quantities.py](../../../nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/conserved_quantities.py) - `GeodesicDiagnostics`, `compute_energy`, `compute_angular_momentum_z_cartesian`, `compute_carter_constant_KerrSchild_Cartesian`
- [analytic_spacetimes_KerrSchild_Cartesian.py](../../../nrpy/equations/general_relativity/geodesics/tests/analytic_spacetimes_KerrSchild_Cartesian.py) - `trusted_dict`
- [geodesics_KerrSchild_Cartesian_massive.py](../../../nrpy/equations/general_relativity/geodesics/tests/geodesics_KerrSchild_Cartesian_massive.py) - `trusted_dict`
- [geodesics_KerrSchild_Cartesian_photon.py](../../../nrpy/equations/general_relativity/geodesics/tests/geodesics_KerrSchild_Cartesian_photon.py) - `trusted_dict`
- [conserved_quantities_KerrSchild_Cartesian_massive.py](../../../nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/tests/conserved_quantities_KerrSchild_Cartesian_massive.py) - `trusted_dict`
- [conserved_quantities_KerrSchild_Cartesian_photon.py](../../../nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/tests/conserved_quantities_KerrSchild_Cartesian_photon.py) - `trusted_dict`

## See Also

- [General Relativity](index.md)
- [BSSN Family](bssn-family.md)
- [Metric Conversions And Matter](metric-conversions-and-matter.md)
- [Psi4 And Tetrads](psi4-and-tetrads.md)
- [Reference Metrics](../../core/reference-metrics.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
