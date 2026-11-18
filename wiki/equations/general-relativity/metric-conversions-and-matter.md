# Metric Conversions And Matter

> Symbolic conversions between ADM, BSSN, four-metric, and four-Christoffel tensors, `T4munu` stress-energy projections, and Cartesian `LorentzBoost` transforms. · Status: confirmed
> Up: [General Relativity](index.md)

## Summary

NRPy's metric-conversion modules translate between ADM spatial quantities, reference-metric BSSN variables, spacetime four-metrics `g_(mu nu)`, and spacetime Christoffel symbols `Gamma^mu_(nu rho)` for initial data, diagnostics, wave extraction, geodesics, and matter coupling. `T4munu.py` projects contravariant `T^(mu nu)` into ADM 3+1 matter terms, and `LorentzBoost.py` transforms Cartesian four-tensors.

## Detail

### ADM, BSSN, Four-Metric, And Four-Christoffel Conversions

`ADM_to_BSSN` converts `gammaDD`, `KDD`, `betaU`, and `BU` into `gammabarDD`, rescaled `hDD`, `cf`, `trK`, `AbarDD`, rescaled `aDD`, `vetU`, and `betU` using the requested `ReferenceMetric` and `enable_rfm_precompute` (`compute_cf_only=True` stops after `cf` for initial-lapse helpers). `BSSN_to_ADM` reads cached `BSSN_quantities[...]` and global `EvolvedConformalFactor_cf` (`phi`, `chi = exp(-4 phi)`, or `W = exp(-2 phi)`) to reconstruct `betaU`, `gammaDD`, `gammaDDdD`, `gammaDDdDD`, `gammaUU`, `detgamma`, `GammaUDD`, `KDD`, and `KDDdD` by rearranging Mewes et al. Eqs. (7) and (10) with chain/product differentiation.

`g4munu_conversions.py` exposes stateless functions: `ADM_to_g4DD` and `ADM_to_g4UU` form `g_(mu nu)` and `g^(mu nu)` from `gammaDD`, `betaU`, and `alpha` (Mewes et al. Eqs. (5)-(6) and matrix inversion); `BSSN_to_g4DD` and `BSSN_to_g4UU` delegate through `BSSN_to_ADM`; `g4DD_to_ADM` extracts `(gammaDD, alpha, betaU)` (round-trip checked against `ADM_to_g4DD` in its doctest, while `tests/g4munu_conversions.py` samples forward and inverse formulas separately); and `g4DD_to_BSSN` returns `(hDD, cf, vetU, alpha)` via `ADM_to_BSSN` with dummy `KDD`/`BU`.

`BSSN_to_g4Christoffel` builds the four-metric, its inverse and derivatives,
and the four-Christoffel symbols. It subtracts the residual trace of
`AbarDD` before reconstructing `KDD` and obtains `gammaDDd0` from ADM
evolution. Callers supply the coordinate-time derivatives `alpha_d0` and
`vetU_d0` of the lapse and rescaled shift. The static variants
`g4DD_dD_static` and `Gamma4UDD_static` keep spatial metric derivatives
and set coordinate-time metric derivatives to zero.

Evidence (scientific/interface requirement): `ADM_to_BSSN.py`, `BSSN_to_ADM.py`, and `general_relativity/g4munu_conversions.py`. The metric conventions follow Mewes et al., Eqs. (5)–(7) and (10).

Evidence (descriptive behavior): `BSSN_to_g4Christoffel.py::BSSN_to_g4Christoffel`.

### Stress-Energy Projections (`T4munu`) And `LorentzBoost`

`T4munu.py` declares symmetric `T4UU`, projects it into `SDD`, `SD`, `S`, and `rho` from explicit ADM tensors or via `BSSN_to_ADM`, registers `CodeParameter` `PI`, and returns matter source terms for `trK_rhs`, `Abar_rhsDD`, `Lambdabar_rhsU`, `H`, and `MU`.

`LorentzBoost` constructs `LorentzMatrix` and `InverseLorentzMatrix` from a 3-component Cartesian velocity `vU` (without `ReferenceMetric` rescaling) and transforms contravariant/covariant 4-vectors and covariant rank-2..4 tensors. It does not check `sum_i v_i^2 < 1`, and passing exact `vU = (0, 0, 0)` divides `0 / 0` in the spatial projector, yielding `nan` rather than the identity matrix; callers must handle zero boost velocity explicitly ([CONTR-0022](../../contradictions.md#contr-0022)).

Claim status: contested; contradiction: CONTR-0022.

Evidence (descriptive behavior): `general_relativity/T4munu.py`; `general_relativity/LorentzBoost.py::LorentzBoost._compute_LorentzMatrix`; `tests/g4munu_conversions.py`.

## Sources

- [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2) - Eqs. (5)-(7), (10)
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

- Depends on: [BSSN Family](bssn-family.md)
- See also: [Initial Data](initial-data.md)
- Validated by: [Expression Validation Helpers](../../validation/expression-validation-helpers.md)
