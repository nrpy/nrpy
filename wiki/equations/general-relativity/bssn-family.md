# BSSN Family

> Map the main BSSN equation modules, analytic algebraic constraints, and validation expectations. · Status: confirmed · Last reconciled: 08-26-2026
> Up: [General Relativity](index.md)

## Summary

The BSSN family is split into reusable symbolic quantities, evolution RHSs,
gauge RHSs, and constraints. The modules build SymPy expressions with explicit
indexed loops, store key outputs on objects or return values, and validate those
outputs through the trusted-expression pipeline. The symbolic evolution
contract assumes the conformal determinant and trace-free-curvature constraints;
applications must project them at each state-repair point they wire.

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

### Analytic determinant and trace constraints

`BSSNQuantities` no longer computes or exposes `trAbar`. `BSSNRHSs` therefore
uses the constrained conformal-metric equation with
`gammabar^ij Abar_ij=0` analytically and contains no off-constraint
`alpha tr(Abar)` correction. `BSSN_to_g4Christoffel` reuses
`BSSN_to_ADM.KDD`, whose owner consumes the already constrained `AbarDD`
without projecting a second private copy.

Three backend kernels provide one combined all-points projection for BSSN
storage. Each first rescales `gammabarDD` to satisfy
`det(gammabar)=det(gammahat)`, inverts that post-enforced metric, and then uses
the same metric to remove the trace from `AbarDD`. Both corrected tensors are
written from one loop. This is algebraic projection, not exponential damping,
and is compatible with BSSN or fCCZ4 states using that storage. Scheduling is
backend- and application-owned: BHaH and superB callers opt in, while ETLegacy
and CarpetX emit schedule entries. Reviewed owner paths and collision examples
do not establish an end-to-end lifecycle for the new fCCZ4 RHS, gauge, and
`Theta_fCCZ4` storage.

Claim evidence:
- Claim: BSSN and fCCZ4 share a constrained-state contract in which `trAbar` is absent from the symbolic quantity/RHS API; `BSSN_to_g4Christoffel` reuses `BSSN_to_ADM.KDD`, whose owner consumes constrained `AbarDD` without a private projection; BHaH, ETLegacy, and CarpetX provide combined determinant/trace projection kernels whose trace uses the post-enforced metric, but application lifecycle wiring is separate and the reviewed owner paths and collision examples do not establish end-to-end fCCZ4 evolution.
- Role: public/scientific contract
- Deciding authority: [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities`; [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs`; [BSSN_to_ADM.py](../../../nrpy/equations/general_relativity/BSSN_to_ADM.py), `BSSN_to_ADM`; [BSSN_to_g4Christoffel.py](../../../nrpy/equations/general_relativity/BSSN_to_g4Christoffel.py), `BSSN_to_g4Christoffel`; backend `register_CFunction_enforce_detgbar_equals_detghat_trAzero` implementations listed in Sources
- Corroboration: none available; the cited symbolic owners and backend kernels are the deciding implementation evidence, while reviewed application paths establish only their own limited wiring
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=BHaH, ETLegacy, CarpetX generated C; precision=exact symbolic off-diagonal determinant/trace identities and generated-source structure; GPU=not-run; restart=not-run; distributed=not-run; error_path=not-run; options=Cartesian and SinhSpherical reference-metric precompute in all three backends, application scheduling inspected only; date=08-26-2026`

`BSSN_gauge_RHSs` handles lapse and shift choices separately from the main RHS
class. It returns `alpha_rhs`, `vet_rhsU`, and `bet_rhsU`, validates supported
lapse options such as `OnePlusLog`, `BHSHarmonicSlicing`, `Frozen`, and
`OnePlusLogAlt`, and validates supported shift options such as frozen and
Gamma-driver variants.

`BHSHarmonicSlicing` names NRPy's Baumgarte-Hughes-Shapiro
`exp(6*phi)`-tracking prescription:
`partial_t(alpha) = partial_t(exp(6*phi))`. The implementation applies that
chain rule to the already constructed total coordinate-time conformal-factor
RHS: `alpha_rhs = -3*W**(-4)*W_rhs` when `W = exp(-2*phi)`, and
`alpha_rhs = 6*exp(6*phi)*phi_rhs` when `phi` is evolved. The option supports
`W` and `phi`, but not `chi`. Its evolution preserves the difference
`alpha-exp(6*phi)`, so it preserves `alpha=exp(6*phi)` only when that relation
is imposed initially. `OnePlusLogAlt` separately implements
`partial_0(alpha)=-alpha*(1-alpha)*K`. These distinct option formulas do not
imply that the `exp(6*phi)`-tracking branch must have a nonzero RHS on exact
stationary data. Descriptively, the BHS branch is an `exp(6*phi)`-tracking
lapse. The former `HarmonicSlicing` option string is unsupported; no alias or
compatibility path is retained.
[Initial Data](initial-data.md) owns the cited StaticTrumpet scientific pairing.

Claim evidence:
- Claim: `BSSN_gauge_RHSs(..., LapseEvolutionOption="BHSHarmonicSlicing")` sets the lapse RHS to the total coordinate-time derivative of `exp(6*phi)` by applying the displayed chain rule to `cf_rhs`; it supports `W` and `phi`, rejects `chi`, preserves `alpha-exp(6*phi)`, and therefore preserves `alpha=exp(6*phi)` only from initial data satisfying that relation. The former `HarmonicSlicing` string is rejected. `OnePlusLogAlt` separately implements `partial_0(alpha)=-alpha*(1-alpha)*K`; neither separate formula implies that the BHS branch must be nonzero on exact stationary data.
- Role: descriptive behavior
- Deciding authority: [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py), `BSSN_gauge_RHSs`; [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs`
- Corroboration: [legacy NRPy gauge tutorial](https://github.com/zachetienne/nrpytutorial/blob/a32e120f5642bee00e32e9e04dd8cb4c58ae661c/Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.ipynb), Steps 2.b and 2.d
- Validation: `inspected=pass; generated=not-run; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=exact symbolic equality and trusted-expression comparison; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=pass for chi and the removed HarmonicSlicing string; options=BHSHarmonicSlicing and OnePlusLogAlt across the module validation matrix, with W, phi, and chi checked for the BHS branch; date=08-25-2026`

General Bona-Masso harmonic slicing is
`(partial_t-beta^i*partial_i)alpha=-alpha**2*K`. At zero shift, Baumgarte and
Shapiro integrate it as `alpha=C(x)*exp(6*phi)`, equivalently preserving
`alpha*exp(-6*phi)`, and choose the special normalization `C(x)=1`. Only that
special case yields `alpha=exp(6*phi)` and hence the derivative equality used
by the NRPy branch; the branch itself instead preserves an additive difference
and is not a general harmonic-slicing implementation. Baumgarte, Hughes, and
Shapiro later used the derivative equality as their zero-shift harmonic
prescription, which explains NRPy's historical option name.

Claim evidence:
- Claim: Bona-Masso harmonic slicing is `(partial_t-beta^i*partial_i)alpha=-alpha**2*K`; at zero shift it preserves `alpha*exp(-6*phi)` and has solution `alpha=C(x)*exp(6*phi)`, with `C(x)=1` giving the derivative equality adopted by the historical NRPy branch.
- Role: public/scientific contract
- Deciding authority: [Baumgarte and de Oliveira, arXiv:2201.08857v1](https://arxiv.org/pdf/2201.08857v1), Eq. (1) and `f(alpha)=1`; [Baumgarte and Shapiro, arXiv:gr-qc/9810065v1](https://arxiv.org/pdf/gr-qc/9810065v1), Eqs. (30)-(32); [Baumgarte, Hughes, and Shapiro, arXiv:gr-qc/9902024v1](https://arxiv.org/pdf/gr-qc/9902024v1), p. 2
- Corroboration: none available; the papers define the scientific relations, while current NRPy option mapping is covered by the descriptive claim above

`BSSNconstraints` constructs Hamiltonian, momentum, and covariant conformal
connection constraint expressions. It registers diagnostic gridfunctions for
`H` and `MSQUARED`, optionally registers `MU`, and stores `H`, `MU`,
`Msquared`, and rescaled `mU` outputs. Its covariant conformal connection
constraint vector is `LambdaConstraintU = LambdabarU - DGammaU`.

Claim evidence:
- Claim: `BSSNconstraints.LambdaConstraintU` represents the covariant conformal connection constraint vector: evolved `LambdabarU` minus `DGammaU`, where `DGammaU` contracts the conformal/reference connection difference with the inverse conformal metric.
- Role: public/scientific contract
- Deciding authority: [Brown, *Covariant formulations of BSSN and the standard gauge*, arXiv:0902.3652v2](https://arxiv.org/pdf/0902.3652v2), Eqs. (12a), (12b), and (15)
- Corroboration: [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints_dict.__getitem__` doctest

For local diagnostics, `BSSNconstraints` also stores the direct conformal-metric
contraction `LambdaConstraintSquared` and its plain square root
`LambdaConstraintMagnitude`.

Claim evidence:
- Claim: `BSSNconstraints` computes `LambdaConstraintSquared` by contracting `LambdaConstraintU` twice with `gammabarDD`, then defines `LambdaConstraintMagnitude` as its plain square root; no physical-metric conformal-factor conversion is applied.
- Role: descriptive behavior
- Deciding authority: [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints.__init__`
- Corroboration: [BSSN_constraints_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_constraints_Cartesian.py), `trusted_dict`; [BSSN_constraints_Spherical.py](../../../nrpy/equations/general_relativity/tests/BSSN_constraints_Spherical.py), `trusted_dict`
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction and BHaH OpenMP C generation; precision=symbolic construction and 30-decimal-digit trusted sampling; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=Cartesian and Spherical symbolic checks through the owner doctest plus six trusted-expression coordinate variants; date=07-28-2026`

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

- [Brown, arXiv:0902.3652v2](https://arxiv.org/pdf/0902.3652v2) - Eqs. (12a), (12b), and (15); published as [Phys. Rev. D 79, 104029](https://doi.org/10.1103/PhysRevD.79.104029) (secondary metadata)
- [Baumgarte and de Oliveira, arXiv:2201.08857v1](https://arxiv.org/pdf/2201.08857v1) - Eq. (1), Bona-Masso slicing and the `f(alpha)=1` harmonic specialization
- [Baumgarte and Shapiro, arXiv:gr-qc/9810065v1](https://arxiv.org/pdf/gr-qc/9810065v1) - Eqs. (30)-(32), zero-shift harmonic slicing and its integrated lapse relation
- [Baumgarte, Hughes, and Shapiro, arXiv:gr-qc/9902024v1](https://arxiv.org/pdf/gr-qc/9902024v1) - p. 2, zero-shift harmonic-slicing relation
- [legacy NRPy gauge tutorial](https://github.com/zachetienne/nrpytutorial/blob/a32e120f5642bee00e32e9e04dd8cb4c58ae661c/Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.ipynb) - Step 2.b, `HarmonicSlicing` chain-rule derivation
- [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py) - `BSSNRHSs`, `BSSN_RHSs_varname_to_expr_dict`
- [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py) - `BSSNQuantities`, `BSSN_quantities`
- [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py) - `BSSN_gauge_RHSs`
- [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py) - `BSSNconstraints`, `BSSN_constraints`
- [BSSN_to_ADM.py](../../../nrpy/equations/general_relativity/BSSN_to_ADM.py) - `BSSN_to_ADM`, canonical constrained-state `KDD` reconstruction
- [BSSN_to_g4Christoffel.py](../../../nrpy/equations/general_relativity/BSSN_to_g4Christoffel.py) - `BSSN_to_g4Christoffel`
- [BHaH algebraic constraint projection](../../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- [ETLegacy algebraic constraint projection](../../../nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- [CarpetX algebraic constraint projection](../../../nrpy/infrastructures/CarpetX/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- [BSSN_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_RHSs_Cartesian.py) - `trusted_dict`
- [BSSN_quantities_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_quantities_Cartesian.py) - `trusted_dict`
- [BSSN_constraints_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_constraints_Cartesian.py) - `trusted_dict`
- [test_parse_BSSN.py](../../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py) - `test_example_BSSN`
- [validate_expressions.py](../../../nrpy/validate_expressions/validate_expressions.py) - `assert_equal`
- [original-agents.md](../../../raw/source-docs/original-agents.md) - `## Equation Setup Rules`

## See Also

- [General Relativity](index.md)
- [Equations](../index.md)
- [Fully Covariant Conformal Z4](fccz4.md)
- [Metric Conversions And Matter](metric-conversions-and-matter.md)
- [Initial Data](initial-data.md)
- [GRHD](../grhd.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
