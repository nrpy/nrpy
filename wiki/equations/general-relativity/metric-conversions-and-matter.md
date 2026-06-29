# Metric Conversions And Matter

> Track ADM, BSSN, four-metric, four-Christoffel, matter-source, and Lorentz-boost helpers. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [General Relativity](index.md)

## Summary

The metric-conversion modules form the boundary between ADM quantities, BSSN
quantities, and four-dimensional spacetime objects. They do not own BSSN
evolution; instead they reconstruct the tensors that neighboring equation
modules need for matter coupling, diagnostics, wave extraction, geodesics, and
generated-code setup.

Matter support is similarly factored: `T4munu.py` projects a symbolic
contravariant stress-energy tensor into BSSN source variables and returns source
terms for BSSN RHSs and constraints. `LorentzBoost.py` is a separate Cartesian
four-dimensional tensor transform helper.

## Detail

`ADM_to_BSSN` converts ADM inputs `gammaDD`, `KDD`, `betaU`, and `BU` into
BSSN-facing fields. Its object state includes `gammabarDD`, `hDD`, `cf`,
`trK`, `AbarDD`, `aDD`, `vetU`, and `betU`. It uses the requested reference
metric, honors `enable_rfm_precompute`, and can stop after the conformal factor
with `compute_cf_only`, which is used by initial-data providers when they need
the standard lapse but not the full BSSN state.

`BSSN_to_ADM` takes the opposite direction from the cached
`BSSN_quantities[...]` object for a coordinate system. It reconstructs
`betaU`, `gammaDD`, `gammaDDdD`, `gammaDDdDD`, `gammaUU`, `detgamma`,
`GammaUDD`, `KDD`, and `KDDdD`. The conformal-factor branch follows the global
`EvolvedConformalFactor_cf` parameter for `phi`, `chi`, or `W`.

`g4munu_conversions.py` exposes function-level conversions rather than a state
class. `ADM_to_g4DD` and `ADM_to_g4UU` build the covariant and contravariant
four-metric from ADM spatial metric, shift, and lapse. `BSSN_to_g4DD` and
`BSSN_to_g4UU` route through `BSSN_to_ADM`. `g4DD_to_ADM` extracts
`gammaDD`, `alpha`, and `betaU` from a four-metric, and `g4DD_to_BSSN` returns
`hDD`, `cf`, `vetU`, and `alpha` after an ADM-to-BSSN pass with dummy
extrinsic curvature and shift-driver variables.

`BSSN_to_g4Christoffel` reconstructs enough ADM and gauge derivative data to
build physical spacetime Christoffels from BSSN quantities. Its public state
includes the four-metric `g4DD`, inverse four-metric `g4UU`, metric derivatives
`g4DD_dD`, spatial-metric time derivative `gammaDDd0`, and `Gamma4UDD`. The
class validates coordinate-system and precompute inputs before using reference
metric caches.

`T4munu.py` declares symbolic `T4UU` and projects it into `SDD`, `SD`, `S`,
and `rho` using either explicit ADM fields or BSSN fields converted through
`BSSN_to_ADM`. It also returns source terms for `trK_rhs`, `Abar_rhsDD`, and
`Lambdabar_rhsU`, plus Hamiltonian and momentum-constraint source terms. The
module registers `PI` as a code parameter for those source terms.

`LorentzBoost` builds `LorentzMatrix` and `InverseLorentzMatrix` from a
Cartesian boost velocity, then applies them to contravariant vectors, covariant
vectors, and covariant rank-2 through rank-4 tensors. This helper is
four-dimensional but deliberately Cartesian; it does not use NRPy reference
metrics.

The representative trusted files show the validation coverage shape: direct
ADM-to-BSSN StaticTrumpet conversion, BSSN-to-ADM coordinate variants,
four-metric round trips and BSSN four-metric variants, four-Christoffel
variants, matter projections/source terms, and Lorentz-boost tensor outputs.

## Sources

- [ADM_to_BSSN.py](../../../nrpy/equations/general_relativity/ADM_to_BSSN.py) - `ADM_to_BSSN`, `compute_cf_only`
- [BSSN_to_ADM.py](../../../nrpy/equations/general_relativity/BSSN_to_ADM.py) - `BSSN_to_ADM`
- [g4munu_conversions.py](../../../nrpy/equations/general_relativity/g4munu_conversions.py) - `ADM_to_g4DD`, `BSSN_to_g4DD`, `ADM_to_g4UU`, `BSSN_to_g4UU`, `g4DD_to_ADM`, `g4DD_to_BSSN`
- [BSSN_to_g4Christoffel.py](../../../nrpy/equations/general_relativity/BSSN_to_g4Christoffel.py) - `BSSN_to_g4Christoffel`, `Gamma4UDD`
- [T4munu.py](../../../nrpy/equations/general_relativity/T4munu.py) - `T4UU_and_ADM_to_SDD_SD_S_rho`, `T4UU_and_BSSN_to_SDD_SD_S_rho`, `BSSN_RHSs_T4UU_source_terms`, `BSSN_constraints_T4UU_source_terms`
- [LorentzBoost.py](../../../nrpy/equations/general_relativity/LorentzBoost.py) - `LorentzBoost`
- [ADM_to_BSSN_StaticTrumpet.py](../../../nrpy/equations/general_relativity/tests/ADM_to_BSSN_StaticTrumpet.py) - `trusted_dict`
- [BSSN_to_ADM_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_to_ADM_Cartesian.py) - `trusted_dict`
- [g4munu_conversions.py](../../../nrpy/equations/general_relativity/tests/g4munu_conversions.py) - `trusted_dict`
- [BSSN_to_g4Christoffel_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_to_g4Christoffel_Cartesian.py) - `trusted_dict`
- [T4munu.py](../../../nrpy/equations/general_relativity/tests/T4munu.py) - `trusted_dict`
- [LorentzBoost.py](../../../nrpy/equations/general_relativity/tests/LorentzBoost.py) - `trusted_dict`

## See Also

- [General Relativity](index.md)
- [BSSN Family](bssn-family.md)
- [Initial Data](initial-data.md)
- [Psi4 And Tetrads](psi4-and-tetrads.md)
- [Geodesics](geodesics.md)
- [GRHD](../grhd.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
