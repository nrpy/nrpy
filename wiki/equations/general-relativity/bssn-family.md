# BSSN Family

> Map the main BSSN equation modules and their validation expectations. · Status: confirmed · Last reconciled: 07-13-2026
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

The `general_relativity/nrpylatex` subpackage is aggregate-covered here rather
than a separate equation family. `test_parse_BSSN.py::test_example_BSSN` parses
one embedded BSSN equation block with the external `nrpylatex.parse_latex`, then
uses `validate_expressions.assert_equal` to compare parsed quantities with the
handwritten Cartesian `RbarDD_gridfunctions` BSSN RHS, default gauge RHS,
constraints, and Ricci tensor. For a tuple parser result, the harness reads
parsed values from module globals; otherwise it reads the returned object's
namespace. This is local compatibility handling, not a claim about every
upstream parser version. The harness declares the conformal metric and
perturbation separately, aliases the barred metric for unadorned index
operations, and covers both braced and unbraced partial-derivative replacement
forms. Because
`assert_equal` evaluates both sides at a fixed sampled point, this is a
deterministic sampled-numerical cross-check, not a formal symbolic identity
proof and not a build, runtime, or accuracy test.

Claim evidence:
- Claim: `test_example_BSSN()` reads module globals for a tuple parser result and otherwise reads the returned object's namespace, uses separate conformal-metric and perturbation declarations plus both partial-derivative replacement forms, and performs only a deterministic sampled-numerical cross-check; this local branch handling does not establish behavior for every upstream parser version.
- Role: descriptive behavior
- Deciding authority: [test_parse_BSSN.py](../../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py), `test_example_BSSN`
- Corroboration: none available; compatibility handling is local harness behavior, and the validation result below exercises only the installed tuple-result path
- Validation: `inspected=pass; generated=not-run; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0, NRPyLaTeX 1.4.0; backend=not-applicable; precision=30 decimal digits; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=tuple-result module-globals path run, returned-object namespace branch inspected only; date=07-13-2026`

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
- [test_parse_BSSN.py](../../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py) - `test_example_BSSN`
- [validate_expressions.py](../../../nrpy/validate_expressions/validate_expressions.py) - `assert_equal`
- [original-agents.md](../../../raw/source-docs/original-agents.md) - `## Equation Setup Rules`

## See Also

- [General Relativity](index.md)
- [Equations](../index.md)
- [Metric Conversions And Matter](metric-conversions-and-matter.md)
- [Initial Data](initial-data.md)
- [GRHD](../grhd.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
