# Geodesics

> Map analytic spacetime metrics, geodesic RHS construction, and conserved diagnostics. · Status: confirmed · Last reconciled: 2026-06-29
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
Eulerian path-length diagnostic, and `p0_photon` comes from the null constraint.
Unsupported particle types raise `ValueError`.

The same module also exposes reusable recipes for numerical spacetime data:
`symbolic_numerical_christoffel_recipe`,
`symbolic_g4DD_recipe_from_bssn_grid_basis`,
`symbolic_christoffel_recipe_from_grid_basis`, and
`symbolic_christoffel_recipe_from_bssn_grid_basis`. The BSSN-grid helpers accept
target bases `Cartesian` and `Spherical`, use reference-metric Jacobians, and
assume a time-independent spatial map lifted into four dimensions.

`GeodesicDiagnostics` can use analytic metrics or a generic `"Numerical"`
metric placeholder. It always constructs `E_expr = -p_0`, builds `L_exprs` only
for spacetime names ending in `Cartesian`, and builds `Q_expr` only for
`KerrSchild_Cartesian`. The Carter-constant path distinguishes massive and
photon particles and relies on the module assumptions documented in
`conserved_quantities.py`.

Representative trusted dictionaries cover the analytic metric, massive and
photon geodesic equations, and massive and photon conserved diagnostics for
`KerrSchild_Cartesian`. The module-level validation also checks metric and
Christoffel symmetry plus identity-map and spherical-to-Cartesian BSSN transform
cases before emitting trusted expressions.

## Sources

- [analytic_spacetimes.py](../../../nrpy/equations/general_relativity/geodesics/analytic_spacetimes.py) - `AnalyticSpacetimes`, `Analytic_Spacetimes`
- [geodesics.py](../../../nrpy/equations/general_relativity/geodesics/geodesics.py) - `GeodesicEquations`, `geodesic_eom_rhs_massive`, `geodesic_eom_rhs_photon`, `symbolic_christoffel_recipe_from_bssn_grid_basis`
- [conserved_quantities.py](../../../nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/conserved_quantities.py) - `GeodesicDiagnostics`, `compute_energy`, `compute_angular_momentum_cartesian`, `compute_carter_constant_KerrSchild_Cartesian`
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
