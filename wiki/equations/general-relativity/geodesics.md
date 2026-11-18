# Geodesics

> Map analytic spacetime metrics, geodesic RHS construction, and conserved diagnostics. · Status: confirmed
> Up: [General Relativity](index.md)

## Summary

The geodesics package builds symbolic metric, Christoffel, equation-of-motion,
constraint, and diagnostic expressions for particle paths. Its analytic
registry implements `KerrSchild_Cartesian`, an analytic-like static
continuation of coincident Brill-Lindquist initial data, and a time-dependent
Cartesian representation of flat spacetime; `GeodesicEquations`
supports massive and photon RHS variants; `GeodesicDiagnostics` provides
energy, axial angular momentum, and Kerr-Schild Carter-constant diagnostics
where the required symmetries are available.

## Detail

### Analytic Metrics

`AnalyticSpacetimes` dispatches by spacetime name. `KerrSchild_Cartesian`
registers `M_scale` and `a_spin`. The
`BrillLindquist_InitialData_Static_Cartesian` recipe registers `M_total` and
uses coincident Brill-Lindquist spatial data with
`gamma_ij = psi^4 delta_ij`, `alpha = psi^-2`, and zero shift. This second
branch defines a static four-metric for controlled tests.
`Minkowski_TimeDependentAxisymmetric_Cartesian` defines inertial coordinates
`T=t`, `(X,Y)=a(t) Rot(Omega(t)+kappa(t)z)(x,y)`, and `Z=b(t)z`. It registers
the four amplitudes, four frequencies, and twist length `flat_L` as CodeParameters.
The metric is formed from the Cartesian coordinate Jacobian. This branch also
stores the inverse coordinate map for exact photon positions.
`GeodesicEquations` differentiates `g4DD` to form `g4DD_dD`. For the
time-dependent flat metric it forms `Gamma4UDD` from the Cartesian coordinate
map's inverse Jacobian and second derivatives; the generic Christoffel formula
remains available for other metrics.
The metric remains flat and rotationally symmetric about the `z`
axis when `a(t)` and `b(t)` are nonzero. All three branches define coordinates
`[t, x, y, z]` and store the covariant four-metric `g4DD`.
Unsupported analytic spacetime names raise `ValueError`. The lazy
`Analytic_Spacetimes` dictionary caches constructed instances.

Evidence (descriptive behavior): `analytic_spacetimes.py::AnalyticSpacetimes` and `Analytic_Spacetimes`.

### Geodesic Equations And Normalization

`GeodesicEquations` obtains `g4DD` and `xx` from the analytic registry, computes
`g4DD_dD`, builds `Gamma4UDD`, and then chooses the RHS and initialization
constraint by `particle_type`. The massive path returns eight ODE right-hand
sides from position and four-velocity evolution, plus `u0_massive` from the
timelike normalization condition. The photon path returns nine ODE right-hand
sides for position, momentum, and signed normal-observer path length
`dL_normal/dlambda = alpha*p^0`, with `alpha = 1/sqrt(-g^{00})` when
`g^{00} < 0`. The normalization quadratics assume signature `(-,+,+,+)`
and `g_00 < 0`. They select the positive massive root and negative photon root.
In Minkowski coordinates these give `u^0 = sqrt(1 + sum_i (u^i)^2)` and
`p^0 = -sqrt(sum_i (p^i)^2)`.

For photon evolution, the module provides four explicit RHS forms: ordinary
`geodesic_eom_rhs_photon` and
`geodesic_eom_rhs_photon_christoffel`, plus normalized
`geodesic_eom_rhs_photon_normalized` and
`geodesic_eom_rhs_photon_normalized_christoffel`. The first member of each pair
consumes metric derivatives; the second consumes Christoffels. EOM-family
selection and selection of metric derivatives or Christoffels are separate generation choices.

The photon normalization helpers make the direct and normalized state
conventions explicit. They define `u = ln|alpha p^0|` and
`Pi_i = p_i/(alpha p^0)`, so the normalized null constraint is
`gamma^ij Pi_i Pi_j = 1`; the direct constraint is
`g_mu_nu p^mu p^nu = 0`. `normal_observer_log_energy` evaluates the same
normal-observer log-energy measure for direct states by constructing the needed
inverse-metric component from covariant `g4DD`. The photon integrators use
this measure for their upper evolution cutoff.

The analytic geometry stored by `GeodesicEquations` and its returned generic
RHS expressions are separate: the RHS and normalization expressions use
`conn_Gamma4UDD`, `metric_g4DD`, and tangent-component symbols. Callers supply
the metric, connection, and tangent values needed to evaluate them.
Unsupported particle types raise `ValueError`.

Evidence (descriptive behavior): `geodesics.py::GeodesicEquations`, `photon_momentum_to_normalized_quantities`, `normalization_constraint_photon_normalized`, and `geodesic_eom_rhs_photon`.

### Numerical Spacetime And Conserved Quantities

The same module also exposes reusable recipes for numerical spacetime data:
`symbolic_numerical_christoffel_recipe`,
`symbolic_g4DD_recipe_from_bssn_grid_basis`,
`symbolic_g4DD_dt_recipe_from_bssn_grid_basis`, and
`symbolic_christoffel_recipe_from_bssn_grid_basis`. The BSSN-grid helpers use
reference-metric Jacobians and assume a time-independent spatial map lifted
into four dimensions; the Christoffel recipe can optionally use static time
derivatives for an explicitly requested approximate endpoint.

Evidence (descriptive behavior): `geodesics.py::symbolic_g4DD_dt_recipe_from_bssn_grid_basis` and `symbolic_christoffel_recipe_from_bssn_grid_basis`.

`GeodesicDiagnostics` can use analytic metrics or a generic `"Numerical"`
metric placeholder. It always constructs `E_expr = -p_0`, builds scalar `Lz_expr` only
for spacetime names ending in `Cartesian`, and builds `Q_expr` only for
`KerrSchild_Cartesian`. The Carter-constant path distinguishes massive and
photon particles and relies on the module assumptions documented in
`conserved_quantities.py`.

Evidence (descriptive behavior): `conserved_quantities.py::GeodesicDiagnostics` and `compute_angular_momentum_z_cartesian`.

Representative trusted dictionaries cover Kerr-Schild and Brill-Lindquist metric
recipes, plus massive and photon geodesic equations and conserved diagnostics
for `KerrSchild_Cartesian`. The module-level validation also checks metric and
Christoffel symmetry plus identity-map and spherical-to-Cartesian BSSN transform
cases before emitting trusted expressions.

## Sources

- [analytic_spacetimes.py](../../../nrpy/equations/general_relativity/geodesics/analytic_spacetimes.py) - `AnalyticSpacetimes`, `Analytic_Spacetimes`
- [geodesics.py](../../../nrpy/equations/general_relativity/geodesics/geodesics.py) - `GeodesicEquations`, `geodesic_eom_rhs_massive`, `geodesic_eom_rhs_photon`, `geodesic_eom_rhs_photon_christoffel`, `geodesic_eom_rhs_photon_normalized`, `geodesic_eom_rhs_photon_normalized_christoffel`, `symbolic_g4DD_dt_recipe_from_bssn_grid_basis`, `symbolic_christoffel_recipe_from_bssn_grid_basis`
- [normal_observer_log_energy.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/normal_observer_log_energy.py) - `normal_observer_log_energy`
- [conserved_quantities.py](../../../nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/conserved_quantities.py) - `GeodesicDiagnostics`, `compute_energy`, `compute_angular_momentum_z_cartesian`, `compute_carter_constant_KerrSchild_Cartesian`
- [analytic_spacetimes_KerrSchild_Cartesian.py](../../../nrpy/equations/general_relativity/geodesics/tests/analytic_spacetimes_KerrSchild_Cartesian.py) - `trusted_dict`
- [analytic_spacetimes_BrillLindquist_InitialData_Static_Cartesian.py](../../../nrpy/equations/general_relativity/geodesics/tests/analytic_spacetimes_BrillLindquist_InitialData_Static_Cartesian.py) - `trusted_dict`
- [geodesics_KerrSchild_Cartesian_massive.py](../../../nrpy/equations/general_relativity/geodesics/tests/geodesics_KerrSchild_Cartesian_massive.py) - `trusted_dict`
- [geodesics_KerrSchild_Cartesian_photon.py](../../../nrpy/equations/general_relativity/geodesics/tests/geodesics_KerrSchild_Cartesian_photon.py) - `trusted_dict`
- [conserved_quantities_KerrSchild_Cartesian_massive.py](../../../nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/tests/conserved_quantities_KerrSchild_Cartesian_massive.py) - `trusted_dict`
- [conserved_quantities_KerrSchild_Cartesian_photon.py](../../../nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/tests/conserved_quantities_KerrSchild_Cartesian_photon.py) - `trusted_dict`

## See Also

- See also: [BSSN Family](bssn-family.md)
- Depends on: [Reference Metrics](../../core/reference-metrics.md)
- Validated by: [Expression Validation Helpers](../../validation/expression-validation-helpers.md)
