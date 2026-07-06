# BSSN Family

> Map the main BSSN equation modules and their validation expectations. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [General Relativity](index.md)

## Summary

The BSSN family is split into reusable symbolic quantities, evolution RHSs,
gauge RHSs, and constraints. The modules build SymPy expressions with explicit
indexed loops, store key outputs on objects or return values, and validate those
outputs through the trusted-expression pipeline.

## Detail

`BSSNQuantities` is the shared setup layer. It registers or declares the BSSN
gridfunctions, builds conformal metric and extrinsic-curvature quantities,
tracks reference-metric rescalings, computes inverse/derivative/Ricci-related
objects, and exposes sorted Ricci names and expressions for code generation.

`BSSNRHSs` consumes `BSSNQuantities` and constructs the non-gauge evolution
right-hand sides. Its public object state includes `cf_rhs`, `trK_rhs`,
`Lambdabar_rhsU`, `h_rhsDD`, `a_rhsDD`, and `lambda_rhsU`; it also assembles
`BSSN_RHSs_varname_to_expr_dict` so generated-code consumers can use stable
names such as `cf_rhs`, `trK_rhs`, `lambda_rhsU0`, `a_rhsDD00`, and `h_rhsDD00`.

`BSSN_gauge_RHSs` handles lapse and shift choices separately from the main RHS
class. It returns `alpha_rhs`, `vet_rhsU`, and `bet_rhsU`, validates supported
lapse options such as `OnePlusLog`, `HarmonicSlicing`, `Frozen`, and
`OnePlusLogAlt`, and validates supported shift options such as frozen and
Gamma-driver variants.

`BSSNconstraints` constructs Hamiltonian and momentum constraint expressions.
It registers diagnostic gridfunctions for `H` and `MSQUARED`, optionally
registers `MU`, and stores `H`, `MU`, `Msquared`, and rescaled `mU` outputs.

Representative trusted files pin the core RHS, quantity, and constraint
dictionaries. Gauge validation is driven by the supported lapse and shift option
names in `BSSN_gauge_RHSs`, while coordinate and reference-metric variants stay
validation evidence rather than new page scope.

The style contract for these modules is part of their interface: tensor
construction uses explicit loops, established suffixes such as `U`, `D`, `DD`,
`dD`, `dDD`, `dupD`, and `rhs`, and validation keys must stay aligned with the
corresponding trusted files.

## Sources

- [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py) - `BSSNRHSs`, `BSSN_RHSs_varname_to_expr_dict`
- [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py) - `BSSNQuantities`, `BSSN_quantities`
- [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py) - `BSSN_gauge_RHSs`
- [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py) - `BSSNconstraints`, `BSSN_constraints`
- [BSSN_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_RHSs_Cartesian.py) - `trusted_dict`
- [BSSN_quantities_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_quantities_Cartesian.py) - `trusted_dict`
- [BSSN_constraints_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_constraints_Cartesian.py) - `trusted_dict`
- [original-agents.md](../../../raw/source-docs/original-agents.md) - `## Equation Setup Rules`

## See Also

- [General Relativity](index.md)
- [Equations](../index.md)
- [Metric Conversions And Matter](metric-conversions-and-matter.md)
- [Initial Data](initial-data.md)
- [GRHD](../grhd.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
