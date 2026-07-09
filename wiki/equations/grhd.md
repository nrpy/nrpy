# GRHD

> Map NRPy's general relativistic hydrodynamics equation builders, flux helpers, and validation coverage. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [Equations](index.md)

## Summary

The GRHD modules build Valencia-style symbolic hydrodynamics expressions on top
of BSSN-derived ADM metric quantities. `GRHD_Equations` owns the conserved
variables, stress-energy tensors, fluxes, source terms, connection terms, and
reference-metric rescalings; companion modules compute characteristic speeds,
HLL interface fluxes, and no-branch min/max expressions used by generated C
code.

## Detail

`GRHD_Equations` initializes a coordinate-system reference metric, imports the
matching `BSSN_quantities` entry, converts BSSN quantities to ADM spatial metric
and extrinsic curvature through `BSSN_to_ADM`, and declares primitive fluid
symbols `u4Ut`, `rhob`, `P`, `h`, `Ye`, `S`, and `rescaledvU`. The physical
three-velocity `VU` and spatial four-velocity components are reconstructed from
`rescaledvU` with the reference-metric rescaling factors `ReU`.

The conserved variables are stored as object attributes. `compute_rho_star`,
`compute_Ye_star`, `compute_S_star`, `compute_tau_tilde`, and
`compute_S_tildeD` build densitized density, electron-fraction, entropy, energy,
and momentum variables from the lapse, conformal volume factor, primitive fluid
state, four-velocity, and stress-energy contractions.

Stress-energy support is split between `compute_T4UU` and `compute_T4UD`.
`compute_T4UU` uses `ADM_to_g4UU` and the perfect-fluid form of `T^{mu nu}`;
`compute_T4UD` lowers one index with `ADM_to_g4DD`. Both methods also store
rescaled tensor forms so curvilinear generated code can use reference-metric
variables without changing the Cartesian-equivalent symbolic form.

Flux terms follow the conserved-variable split. `compute_rho_star_fluxU`,
`compute_Ye_star_fluxU`, and `compute_S_star_fluxU` advect the corresponding
conserved scalars with `VU`; `compute_tau_tilde_fluxU` combines the energy
stress tensor flux with the density advection subtraction; and
`compute_S_tilde_fluxUD` stores the momentum-flux matrix. Each flux method also
stores a rescaled form divided or transformed by `ReU` where needed.

Source terms are separated from flux construction. `compute_tau_source_term`
uses extrinsic curvature, shift, lapse derivatives, and `T4UU`.
`compute_S_tilde_source_termD` combines lapse-gradient, shift-gradient, and
covariant spatial-metric derivative terms. `compute_all_connection_terms` and
`compute_S_tilde_connection_termsD` add reference-Christoffel contributions
from `GammahatUDD` to the density, electron-fraction, entropy, energy, and
momentum equations. `construct_all_equations` runs the full setup order used by
the trusted GRHD equation tests.

`find_cp_cm` computes the two characteristic speeds in one flux direction from
the contravariant four-metric, four-velocity, and sound speed squared. It uses
the no-branch maximum helper to clamp the quadratic discriminant before taking
the square root, then orders the two speeds with no-branch min/max. `find_cmax_cmin`
builds the face-centered four-metric from ADM face data, evaluates the right
and left speeds, and returns the nonnegative `cmin` and `cmax` values required
by HLL fluxes.

`calculate_GRHD_Tmunu_and_contractions` evaluates the conserved variables and
physical fluxes for one reconstructed side of a cell face. `calculate_HLL_fluxes`
does that for right and left states, obtains `cmin` and `cmax`, and applies
`HLL_solver` to `rho_star`, `Ye_star`, `S_star`, `tau_tilde`, and each
component of `S_tildeD`.

`Min_Max_and_Piecewise_Expressions.py` provides symbolic branch-avoidance
helpers. `min_noif` and `max_noif` express extrema through `nrpyAbs`, which
later becomes a generated C absolute-value call. The coordinate-bound helpers
register `TINYDOUBLE` on demand and return symbolic 0-or-1 masks for
less-than, less-or-equal, greater-than, and greater-or-equal comparisons.

Representative trusted files cover both Cartesian-equivalent and curvilinear
paths. `GRHD_equations_Cartesian.py`, `GRHD_equations_Spherical.py`, and
`GRHD_equations_SinhSpherical_rfm_precompute.py` validate conserved variables,
fluxes, source terms, connection terms, tensor rescalings, and coordinate
variants. Separate trusted dictionaries validate characteristic speeds, HLL
flux assembly, and the no-branch min/max helpers.

## Sources

- [GRHD_equations.py](../../nrpy/equations/grhd/GRHD_equations.py) - `GRHD_Equations`, `construct_all_equations`
- [GRHD_equations.py](../../nrpy/equations/grhd/GRHD_equations.py) - `compute_rho_star`, `compute_Ye_star`, `compute_S_star`, `compute_tau_tilde`, `compute_S_tildeD`
- [GRHD_equations.py](../../nrpy/equations/grhd/GRHD_equations.py) - `compute_T4UU`, `compute_T4UD`, `compute_rho_star_fluxU`, `compute_tau_tilde_fluxU`, `compute_S_tilde_fluxUD`
- [GRHD_equations.py](../../nrpy/equations/grhd/GRHD_equations.py) - `compute_tau_source_term`, `compute_all_connection_terms`, `compute_S_tilde_source_termD`, `compute_S_tilde_connection_termsD`
- [characteristic_speeds.py](../../nrpy/equations/grhd/characteristic_speeds.py) - `find_cp_cm`, `find_cmax_cmin`
- [HLL_fluxes.py](../../nrpy/equations/grhd/HLL_fluxes.py) - `calculate_GRHD_Tmunu_and_contractions`, `HLL_solver`, `calculate_HLL_fluxes`
- [Min_Max_and_Piecewise_Expressions.py](../../nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py) - `min_noif`, `max_noif`, `coord_leq_bound`, `coord_geq_bound`, `coord_less_bound`, `coord_greater_bound`
- [GRHD_equations_Cartesian.py](../../nrpy/equations/grhd/tests/GRHD_equations_Cartesian.py) - `trusted_dict`
- [GRHD_equations_Spherical.py](../../nrpy/equations/grhd/tests/GRHD_equations_Spherical.py) - `trusted_dict`
- [GRHD_equations_SinhSpherical_rfm_precompute.py](../../nrpy/equations/grhd/tests/GRHD_equations_SinhSpherical_rfm_precompute.py) - `trusted_dict`
- [characteristic_speeds.py](../../nrpy/equations/grhd/tests/characteristic_speeds.py) - `trusted_dict`
- [HLL_fluxes.py](../../nrpy/equations/grhd/tests/HLL_fluxes.py) - `trusted_dict`
- [Min_Max_and_Piecewise_Expressions.py](../../nrpy/equations/grhd/tests/Min_Max_and_Piecewise_Expressions.py) - `trusted_dict`

## See Also

- [Equations](index.md)
- [BSSN Family](general-relativity/bssn-family.md)
- [Metric Conversions And Matter](general-relativity/metric-conversions-and-matter.md)
- [Fishbone-Moncrief](general-relativity/fishbone-moncrief.md)
- [Reference Metrics](../core/reference-metrics.md)
- [Trusted Expression Pipeline](trusted-expression-pipeline.md)
