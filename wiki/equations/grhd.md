# GRHD

> Map NRPy's general relativistic hydrodynamics equation builders, flux helpers, and validation coverage. · Status: confirmed
> Up: [Equations](index.md)

## Summary

The GRHD modules build perfect-fluid, Valencia-style symbolic hydrodynamics
expressions on top of BSSN-derived ADM metric quantities. `GRHD_Equations` owns the conserved
variables, stress-energy tensors, fluxes, source terms, connection terms, and
reference-metric rescalings; companion modules compute characteristic speeds,
HLL interface fluxes, and no-branch min/max expressions used by generated C
code. This equation slice does not add magnetic-field terms or primitive-state
recovery; `GRMHDEquations` in [GRMHD](grmhd.md) subclasses `GRHD_Equations` to
add magnetic stress-energy.

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
and momentum variables from the lapse, reference-metric volume factor
`e6phi` (equal to `sqrt(gamma/gammahat)` when `det(gammabar) = det(gammahat)`),
primitive fluid state, four-velocity, and stress-energy contractions.

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
covariant spatial-metric derivative terms; the reference covariant derivative
of `gamma_ij` follows the metric split of
[Jacques et al., Eqs. (2)-(3)](https://arxiv.org/pdf/2412.03659v2). `compute_all_connection_terms` and
`compute_S_tilde_connection_termsD` add reference-Christoffel contributions
from `GammahatUDD` to the density, electron-fraction, entropy, energy, and
momentum equations. `construct_all_equations` runs the full setup order used by
the trusted GRHD equation tests.

The perfect-fluid tensor is the fluid part of
[Duez et al., Eq. (16)](https://arxiv.org/pdf/astro-ph/0503420v2).
Their Eqs. (34)-(39) give mass, momentum, and energy variables, fluxes, and the
energy source before reference-metric densitization.
[Jacques et al., Eqs. (4)-(5) and (13)-(21)](https://arxiv.org/pdf/2412.03659v2)
supply the reference volume factor, electron fraction, conserved system,
geometric sources, and connection terms. Their Eqs. (22)-(24) define component
rescaling. `compute_S_star` sets `S_star = alpha * e6phi * S * u^0`, and
`compute_S_star_fluxU` transports `S_star` with `VU` and no source term. The
caller supplies the primitive entropy variable `S` in the equation-of-state
convention of the calling code. Neither the `S_star` definition nor its
transport law comes from Duez et al. or Jacques et al.; Jacques et al.,
Eq. (11), only rewrites the chosen current's divergence.

Claim evidence:
- Claim: The perfect-fluid tensor, conserved mass, momentum, and energy variables, fluxes, and energy source follow Duez et al. Eqs. (16) and (34)-(39), with reference-metric volume factor, electron fraction, conserved system, sources, connection terms, and component rescaling from Jacques et al. Eqs. (4)-(5) and (13)-(24). Neither paper defines `S_star` or its transport law.
- Role: public/scientific contract
- Deciding authority: [Duez et al. (2005)](https://arxiv.org/pdf/astro-ph/0503420v2), Eqs. (16) and (34)-(39); [Jacques et al.](https://arxiv.org/pdf/2412.03659v2), Eqs. (4)-(5), (11), and (13)-(24)
- Corroboration: [GRHD_equations.py](../../nrpy/equations/grhd/GRHD_equations.py), method docstrings including `compute_S_star`; [GRHD_equations_Cartesian.py](../../nrpy/equations/grhd/tests/GRHD_equations_Cartesian.py), `trusted_dict`

`find_cp_cm` computes the smaller and larger roots of the characteristic
quadratic in one flux direction from the contravariant four-metric,
four-velocity, and squared fluid-frame signal speed `v02`. GRHD supplies the
sound speed squared `c_s^2`; GRMHD supplies its estimate `v_0^2`. It uses the
no-branch maximum helper to clamp the quadratic discriminant before taking
the square root, then orders the roots with no-branch min/max.
The quadratic is derived from
[Duez et al., Eqs. (49)-(50)](https://arxiv.org/pdf/astro-ph/0503420v2);
the paper does not print its expanded coefficients.
`find_cmax_cmin` builds the face-centered four-metric from ADM face data,
evaluates the right and left reconstructed states with their supplied `v02_r`
and `v02_l`, and returns nonnegative `(cmin, cmax)` HLL bounds. See
[GRMHD](grmhd.md) for its magnetic speed estimate.

Claim evidence:
- Claim: `find_cp_cm` uses one squared fluid-frame speed `v02` to return ordered roots; `find_cmax_cmin` uses the two face-state values `v02_r` and `v02_l` to return nonnegative HLL bounds. GRHD supplies sound speeds and GRMHD supplies magnetic speed estimates.
- Role: descriptive behavior
- Deciding authority: [characteristic_speeds.py](../../nrpy/equations/grhd/characteristic_speeds.py), `find_cp_cm` and `find_cmax_cmin`
- Corroboration: [GRMHD characteristic_speeds.py](../../nrpy/equations/grmhd/characteristic_speeds.py), `find_cmax_cmin` caller

These helpers assume caller-supplied physical states and usable denominators.
`flux_dirn` is used as a spatial index and is not range-checked; intended values
are `0`, `1`, or `2`. `find_cp_cm` divides by its quadratic coefficient `a`, and
`HLL_solver` divides by `cmax + cmin`; neither function supplies a zero-
denominator fallback. `find_cmax_cmin` also uses one face metric for both
reconstructed states, as stated in its source docstring.
The HLL flux and its nonnegative speed bounds are
[Duez et al., Eq. (48)](https://arxiv.org/pdf/astro-ph/0503420v2)
and the definitions immediately before it.

`calculate_Tmunu_and_contractions_from_equations` overwrites the supplied
Cartesian equation object's metric and fluid attributes, then computes its
conserved variables and physical fluxes for one reconstructed face state.
`grhd.HLL_fluxes.calculate_HLL_fluxes` passes a Cartesian `GRHD_Equations`
object for each state. `grmhd.HLL_fluxes.calculate_HLL_fluxes` passes a
Cartesian `GRMHDEquations` object with its magnetic field assigned.
The GRHD `calculate_HLL_fluxes` obtains `cmin` and `cmax` and applies
`HLL_solver` to `rho_star`, `Ye_star`, `S_star`, `tau_tilde`, and each
component of `S_tildeD`.

Claim evidence:
- Claim: `calculate_Tmunu_and_contractions_from_equations` replaces the Cartesian equation object's metric and fluid state and computes conserved variables and fluxes. The GRHD and GRMHD HLL functions pass their respective equation classes.
- Role: descriptive behavior
- Deciding authority: [HLL_fluxes.py](../../nrpy/equations/grhd/HLL_fluxes.py), `calculate_Tmunu_and_contractions_from_equations` and `calculate_HLL_fluxes`
- Corroboration: [GRMHD HLL_fluxes.py](../../nrpy/equations/grmhd/HLL_fluxes.py), `calculate_HLL_fluxes` caller

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

That last evidence has an important boundary: the
`Min_Max_and_Piecewise_Expressions.py` script replaces symbolic `nrpyAbs` with
`sin`, `cos`, or `exp` before producing trusted values. It fingerprints the
resulting sampled SymPy expressions but does not test the later `nrpyAbs`-to-C
absolute-value lowering or prove the 0/1 mask semantics over coordinate
domains. All cited trusted files are sampled expression checks, not generated
build, runtime, shock-tube, conservation, convergence, or numerical-accuracy
tests.

## Sources

- [Duez et al. (2005)](https://arxiv.org/pdf/astro-ph/0503420v2) - Eqs. (16), (34)-(39), (48)-(50)
- [Jacques et al.](https://arxiv.org/pdf/2412.03659v2) - Eqs. (2)-(5), (11), (13)-(24)
- [GRHD_equations.py](../../nrpy/equations/grhd/GRHD_equations.py) - `GRHD_Equations`, `construct_all_equations`
- [GRHD_equations.py](../../nrpy/equations/grhd/GRHD_equations.py) - `compute_rho_star`, `compute_Ye_star`, `compute_S_star`, `compute_tau_tilde`, `compute_S_tildeD`
- [GRHD_equations.py](../../nrpy/equations/grhd/GRHD_equations.py) - `compute_T4UU`, `compute_T4UD`, `compute_rho_star_fluxU`, `compute_tau_tilde_fluxU`, `compute_S_tilde_fluxUD`
- [GRHD_equations.py](../../nrpy/equations/grhd/GRHD_equations.py) - `compute_tau_source_term`, `compute_all_connection_terms`, `compute_S_tilde_source_termD`, `compute_S_tilde_connection_termsD`
- [characteristic_speeds.py](../../nrpy/equations/grhd/characteristic_speeds.py) - `find_cp_cm`, `find_cmax_cmin`
- [HLL_fluxes.py](../../nrpy/equations/grhd/HLL_fluxes.py) - `calculate_Tmunu_and_contractions_from_equations`, `HLL_solver`, `calculate_HLL_fluxes`
- [Min_Max_and_Piecewise_Expressions.py](../../nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py) - `min_noif`, `max_noif`, `coord_leq_bound`, `coord_geq_bound`, `coord_less_bound`, `coord_greater_bound`
- [GRHD_equations_Cartesian.py](../../nrpy/equations/grhd/tests/GRHD_equations_Cartesian.py) - `trusted_dict`
- [GRHD_equations_Spherical.py](../../nrpy/equations/grhd/tests/GRHD_equations_Spherical.py) - `trusted_dict`
- [GRHD_equations_SinhSpherical_rfm_precompute.py](../../nrpy/equations/grhd/tests/GRHD_equations_SinhSpherical_rfm_precompute.py) - `trusted_dict`
- [characteristic_speeds.py](../../nrpy/equations/grhd/tests/characteristic_speeds.py) - `trusted_dict`
- [HLL_fluxes.py](../../nrpy/equations/grhd/tests/HLL_fluxes.py) - `trusted_dict`
- [Min_Max_and_Piecewise_Expressions.py](../../nrpy/equations/grhd/tests/Min_Max_and_Piecewise_Expressions.py) - `trusted_dict`

## See Also

- Parent: [Equations](index.md)
- Depends on: [BSSN Family](general-relativity/bssn-family.md)
- Depends on: [Metric Conversions And Matter](general-relativity/metric-conversions-and-matter.md)
- See also: [Fishbone-Moncrief](general-relativity/fishbone-moncrief.md)
- Depends on: [Reference Metrics](../core/reference-metrics.md)
- Validated by: [Trusted Expression Pipeline](trusted-expression-pipeline.md)
- See also: [GRMHD](grmhd.md)
