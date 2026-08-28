# BSSN Family

> Map the main BSSN equation modules, analytic algebraic constraints, and validation expectations. · Status: confirmed · Last reconciled: 08-28-2026
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
`alpha tr(Abar)` correction. These symbolic classes perform no evolved-state
projection; backend and application repair code owns any required projection.

Claim evidence:
- Claim: `BSSNQuantities` does not expose `trAbar`, and `BSSNRHSs` contains no constraint-proportional `alpha tr(Abar)` term; callers own algebraic projection of the evolved state.
- Role: descriptive behavior
- Deciding authority: [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities`; [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs`
- Corroboration: [test_parse_BSSN.py](../../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py), `test_example_BSSN`
- Validation: `inspected=pass; generated=not-applicable; built=not-applicable; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=fixed sampled trusted-expression and cross-representation comparisons; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=six BSSN RHS coordinate variants and nine BSSN quantity variants; date=08-26-2026`

BHaH and ETLegacy provide one combined all-points projector for the determinant
and trace constraints. The GeneralRFM branch explicitly forms the
determinant-corrected metric, inverts it, and removes the resulting trace from
`AbarDD`. Standard orthogonal reference metrics instead invert the reduced
`H=I+h` matrix and use algebraic cancellation of the determinant-rescaling
factor; the result is the equivalent trace removal relative to the corrected
metric. The generated block loads all six `hDD` and six `aDD` components before
its first store, computes six corrected metric and six corrected curvature
outputs, then stores all twelve. For finite, nonsingular states where the
rescaling and inverse are defined, this is algebraic projection to roundoff,
not exponential damping.

BHaH callers opt in. Its `two_blackholes_collide` example enables projection
for initial data and calls the same combined projector after boundary handling
in its Method of Lines hook. ETLegacy schedules its combined projector in
`CCTK_INITIAL`, after evolved-variable boundary conditions in `MoL_PostStep`,
and before constraints in `MoL_PseudoEvolution`. The inspected ETLegacy
projector registration emits no determinant-only repair entry.

Claim evidence:
- Claim: BHaH and ETLegacy use a combined determinant/trace projector; the inspected BHaH collision caller applies it after boundary handling, while the inspected ETLegacy projector registration assigns the same combined projector to all three repair schedules and emits no determinant-only repair entry.
- Role: descriptive behavior
- Deciding authority: [BHaH algebraic constraint projection](../../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`; [BHaH initial data](../../../nrpy/infrastructures/BHaH/general_relativity/initial_data.py), `register_CFunction_initial_data`; [two_blackholes_collide.py](../../../nrpy/examples/two_blackholes_collide.py), Method of Lines registration; [ETLegacy algebraic constraint projection](../../../nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- Corroboration: none available; BHaH and ETLegacy own separate backend behavior
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=ETLegacy generated C plus BHaH source inspection; precision=exact determinant/trace identities and generated-source structure; GPU=not-run; restart=not-run; distributed=not-run; error_path=not-run; options=ETLegacy Cartesian, SinhSpherical reference-metric precompute, and GeneralRFM projector variants; BHaH Spherical collision wiring inspected only; date=08-26-2026`

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
connection constraint expressions. It unconditionally registers scalar
diagnostic gridfunctions `H`, `M`, and `LAMBDA_CONSTRAINT`, optionally
registers the `MU` components, and stores `H`, `MU`, the internal contraction
`Msquared`, and rescaled `mU`. The internal contraction is the physical norm
square `gamma_ij M^i M^j`, implemented as the conformal-metric contraction
divided by `exp_m4phi`; infrastructure producers export `M` as
`sqrt(Msquared)`.

Claim evidence:
- Claim: `BSSNconstraints` registers scalar diagnostic storage named `H`, `M`, and `LAMBDA_CONSTRAINT`; its internal `Msquared` is `gamma_ij M^i M^j` with `gamma_ij = gammabar_ij / exp_m4phi`, and the BHaH, ETLegacy, and CarpetX constraint producers export `M = sqrt(Msquared)` rather than the internal square.
- Role: descriptive behavior
- Deciding authority: [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints.__init__`; [BHaH constraints_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/constraints_eval.py), `register_CFunction_constraints_eval`; [ETLegacy BSSN_constraints.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/BSSN_constraints.py), `register_CFunction_BSSN_constraints`; [CarpetX BSSN_constraints.py](../../../nrpy/infrastructures/CarpetX/general_relativity/BSSN_constraints.py), `register_CFunction_BSSN_constraints`
- Corroboration: [BSSN_constraints_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_constraints_Cartesian.py), `trusted_dict`; [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints_dict.__getitem__` physical-ADM-norm doctest; [BHaH constraints_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/constraints_eval.py), existing code-generation doctest
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction plus BHaH, ETLegacy, and CarpetX registration; precision=symbolic comparison and 30-significant-digit deterministic trusted sampling; GPU=not-run; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=W, chi, and phi physical-norm checks plus six trusted coordinate variants; date=08-28-2026`

Its covariant conformal connection constraint vector is
`LambdaConstraintU = LambdabarU - DGammaU`.

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
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction and BHaH OpenMP C generation; precision=symbolic construction and 30-decimal-digit trusted sampling; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=Cartesian and Spherical symbolic checks through the owner doctest plus six trusted-expression coordinate variants; date=08-28-2026`

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
- [BHaH constraints_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/constraints_eval.py) - `register_CFunction_constraints_eval`
- [ETLegacy BSSN_constraints.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/BSSN_constraints.py) - `register_CFunction_BSSN_constraints`
- [CarpetX BSSN_constraints.py](../../../nrpy/infrastructures/CarpetX/general_relativity/BSSN_constraints.py) - `register_CFunction_BSSN_constraints`
- [BHaH algebraic constraint projection](../../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- [BHaH initial data](../../../nrpy/infrastructures/BHaH/general_relativity/initial_data.py) - `register_CFunction_initial_data`
- [two_blackholes_collide.py](../../../nrpy/examples/two_blackholes_collide.py) - combined-projector registration and Method of Lines hook
- [ETLegacy algebraic constraint projection](../../../nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
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
