# BSSN Family

> Map the main BSSN equation modules, analytic algebraic constraints, and validation expectations. · Status: confirmed
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

The first `BSSNQuantities` construction registers the BSSN gridfunctions; a
later construction for another cache key, made after they are registered,
declares the same names. Both paths give `real=True` SymPy symbols for the
scalars `trK`, `cf`, and `alpha`, and symbols without assumptions for the
tensor components, so expressions built from different `BSSN_quantities` cache
entries share one symbol per gridfunction and can be substituted,
differentiated, or simplified together.

Evidence (descriptive behavior): [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities.__init__`; [grid.py](../../../nrpy/grid.py), `register_gridfunctions` (no trusted file distinguishes symbol assumptions)

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

Evidence (descriptive behavior): [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities`; [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs`; corroboration: [test_parse_BSSN.py](../../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py), `test_example_BSSN`

`BSSN_algebraic_constraints(CoordSystem, enable_rfm_precompute=False)`
returns the corrected rescaled tensors `(hprimeDD, aprimeDD)`. This equation
function selects `_project_general_rfm` for names beginning with `GeneralRFM`
and `_project_orthogonal_rfm` otherwise. BHaH and ETLegacy call it to generate
one combined all-points determinant/trace projector; the backends supply
storage, loops, and schedules.

Evidence (descriptive behavior): [BSSN_algebraic_constraints.py](../../../nrpy/equations/general_relativity/BSSN_algebraic_constraints.py), `BSSN_algebraic_constraints`, `_project_general_rfm`, and `_project_orthogonal_rfm`; corroboration: [BHaH projector](../../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) and [ETLegacy projector](../../../nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`

The GeneralRFM branch explicitly forms the
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
and before constraints in `MoL_PseudoEvolution`. The ETLegacy
projector registration emits no determinant-only repair entry.

Evidence (descriptive behavior): [BHaH algebraic constraint projection](../../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`; [BHaH initial data](../../../nrpy/infrastructures/BHaH/general_relativity/initial_data.py), `register_CFunction_initial_data`; [two_blackholes_collide.py](../../../nrpy/examples/two_blackholes_collide.py), Method of Lines registration; [ETLegacy algebraic constraint projection](../../../nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`

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

Evidence (descriptive behavior): [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py), `BSSN_gauge_RHSs`; [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs`; corroboration: [legacy NRPy tutorial repository](https://github.com/zachetienne/nrpytutorial), `Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.ipynb`, Steps 2.b and 2.d

General Bona-Masso harmonic slicing is
`(partial_t-beta^i*partial_i)alpha=-alpha**2*K`. At zero shift, Baumgarte and
Shapiro integrate it as `alpha=C(x)*exp(6*phi)`, equivalently preserving
`alpha*exp(-6*phi)`, and choose the special normalization `C(x)=1`. Only that
special case yields `alpha=exp(6*phi)` and hence the derivative equality used
by the NRPy branch; the branch itself instead preserves an additive difference
and is not a general harmonic-slicing implementation. Baumgarte, Hughes, and
Shapiro later used the derivative equality as their zero-shift harmonic
prescription, which explains NRPy's historical option name.

Evidence (public/scientific contract): [Baumgarte and de Oliveira, arXiv:2201.08857v1](https://arxiv.org/pdf/2201.08857v1), Eq. (1) and `f(alpha)=1`; [Baumgarte and Shapiro, arXiv:gr-qc/9810065v1](https://arxiv.org/pdf/gr-qc/9810065v1), Eqs. (30)-(32); [Baumgarte, Hughes, and Shapiro, arXiv:gr-qc/9902024v1](https://arxiv.org/pdf/gr-qc/9902024v1), p. 2

`BSSNconstraints` constructs Hamiltonian, momentum, and covariant conformal
connection constraint expressions. It registers scalar diagnostic `H`. The
Boolean NRPy parameter `register_M_and_LAMBDA_CONSTRAINT_gridfunctions`,
registered by `BSSN_constraints.py` with default `True`, also registers `M` and
`LAMBDA_CONSTRAINT`; setting it to `False` with `par.set_parval_from_str`
before construction suppresses their storage registration while retaining the
expressions. It is an NRPy parameter, not a constructor argument. The separate
NRPy parameter `register_MU_gridfunctions` registers the `MU` components.
`BSSN_constraints` rebuilds a cached entry when either parameter differs from
its value at that entry's construction. The object stores `H`, `MU`, the
internal contraction `Msquared`, and rescaled `mU`. The internal contraction is
the physical norm square `gamma_ij M^i M^j`, implemented as the conformal-metric
contraction divided by `exp_m4phi`; the BHaH, ETLegacy, and CarpetX constraint
producers export `M` as `sqrt(Msquared)`.

Evidence (descriptive behavior): [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints.__init__`; [BHaH constraints_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/constraints_eval.py), `register_CFunction_constraints_eval`; [ETLegacy BSSN_constraints.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/BSSN_constraints.py), `register_CFunction_BSSN_constraints`; [CarpetX BSSN_constraints.py](../../../nrpy/infrastructures/CarpetX/general_relativity/BSSN_constraints.py), `register_CFunction_BSSN_constraints`; corroboration: [BSSN_constraints_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_constraints_Cartesian.py), `trusted_dict`; [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints_dict.__getitem__` physical-ADM-norm doctest; [BHaH constraints_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/constraints_eval.py), existing code-generation doctest

Its covariant conformal connection constraint vector is
`LambdaConstraintU = LambdabarU - DGammaU`.
Brown's connection equation adds `C^k Dhat_k beta^i` to the coordinate-time
RHS. Combined with the Lie stretch, its shift-gradient part is
`-DGammaU[k] Dhat_k betaU[i]`; the change vanishes when `C^i=0`.

Evidence (public/scientific contract): [Brown, *Covariant formulations of BSSN and the standard gauge*, arXiv:0902.3652v2](https://arxiv.org/pdf/0902.3652v2), Eqs. (12a), (12b), (15), (21e), and (22e); corroboration: [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs`; [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints_dict.__getitem__`

### Opt-in YBS Gamma-constraint adjustment

`BSSNRHSs(..., enable_YBS_Gamma_constraint_adjustment=True)` adds the
lower-order term `-YBS_chi * (LambdabarU[i] - DGammaU[i]) * Dbarbetacontraction`
to the unrescaled conformal-connection RHS, and the runtime default
`YBS_chi = 0.0` removes it exactly. [YBS Gamma-Constraint
Adjustment](ybs-gamma-constraint-adjustment.md) gives its mapping to the
published coefficients, timestep effect, parameter registration, gauge
forwarding, and constraint propagation.

### Constraint diagnostics and validation

For local diagnostics, `BSSNconstraints` also stores the direct conformal-metric
contraction `LambdaConstraintSquared` and its plain square root
`LambdaConstraintMagnitude`; no physical-metric conformal-factor conversion is
applied.

Evidence (descriptive behavior): [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints.__init__`; corroboration: [BSSN_constraints_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_constraints_Cartesian.py), `trusted_dict`; [BSSN_constraints_Spherical.py](../../../nrpy/equations/general_relativity/tests/BSSN_constraints_Spherical.py), `trusted_dict`

Representative trusted files pin the core RHS, quantity, and constraint
dictionaries. All BSSN RHS trusted cases jointly enable YBS Gamma and
YBS-MOM, so they pin the changed `a_rhsDD` expressions and every other output
already present in those dictionaries without separate option-combination files.
Gauge validation is driven by the supported lapse and shift option
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

Evidence (descriptive behavior): [test_parse_BSSN.py](../../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py), `test_example_BSSN`

The style contract for these modules is part of their interface: tensor
construction uses explicit loops, established suffixes such as `U`, `D`, `DD`,
`dD`, `dDD`, `dupD`, and `rhs`, and validation keys must stay aligned with the
corresponding trusted files.

### Kreiss-Oliger dissipation terms

`add_KreissOliger_dissipation_terms` in `kreiss_oliger_terms.py` adds
Kreiss-Oliger dissipation to an RHS dictionary in place. BHaH and Dendro
`register_CFunction_rhs_eval` call it, and `build_fccz4_expression_bundle`
calls it for fCCZ4 with `include_Theta_fCCZ4=True`. It registers the
CodeParameters `KreissOliger_strength_gauge` and
`KreissOliger_strength_nongauge` under the caller-supplied
`registering_module`; the parameter names are the same for every caller, and
the module name appears in the generated parameter comments and parameter-file
sections.

For each direction `k`, the helper adds a strength times `Q_dKOD[k]` times the
reference-metric factor `ReU[k]` to the right-hand side of each dissipated
field `Q`. The gauge strength multiplies the `alpha`, `vetU`, and `betU` terms,
and the non-gauge strength multiplies the `cf`, `trK`, `lambdaU`, `aDD`, and `hDD`
terms, plus `Theta_fCCZ4` when `include_Theta_fCCZ4` is set. With
`enable_CAKO`, both strengths are multiplied by the caller's conformal-factor
expression `W`; [YBS-MOM Timestep-Scaled Momentum
Adjustment](ybs-momentum-damping.md) records that weight for each
conformal-factor choice.

A gauge field receives dissipation only when its RHS entry is not exactly the
SymPy zero when the helper is called. `BSSN_gauge_RHSs` returns that zero for a
field that does not evolve, so the helper adds nothing to a `Frozen` lapse, a
`Frozen` shift and its driver, or the `betU` driver of the first-order shift
options, and it dissipates the evolved `B^i` of `NonAdvectingGammaDriving` and
of the second-order drivers. Callers insert the gauge RHSs before calling the
helper, as BHaH, Dendro, and `build_fccz4_expression_bundle` do; terms a caller
adds afterwards, such as the slow-start lapse term, do not affect the choice.

Evidence (descriptive behavior): [kreiss_oliger_terms.py](../../../nrpy/equations/general_relativity/kreiss_oliger_terms.py), `add_KreissOliger_dissipation_terms`; [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py), `BSSN_gauge_RHSs`; corroboration: [kreiss_oliger_terms_Cartesian_CAKO.py](../../../nrpy/equations/general_relativity/tests/kreiss_oliger_terms_Cartesian_CAKO.py), `trusted_dict`, which samples the CAKO and `Theta_fCCZ4` branches with every gauge RHS nonzero; no trusted file samples the helper with a zero gauge RHS entry

## Sources

- [BSSN_algebraic_constraints.py](../../../nrpy/equations/general_relativity/BSSN_algebraic_constraints.py) - `BSSN_algebraic_constraints`

- [Brown, arXiv:0902.3652v2](https://arxiv.org/pdf/0902.3652v2) - Eqs. (12a), (12b), (15), (21e), and (22e); published as [Phys. Rev. D 79, 104029](https://doi.org/10.1103/PhysRevD.79.104029) (secondary metadata)
- [Baumgarte and de Oliveira, arXiv:2201.08857v1](https://arxiv.org/pdf/2201.08857v1) - Eq. (1), Bona-Masso slicing and the `f(alpha)=1` harmonic specialization
- [Baumgarte and Shapiro, arXiv:gr-qc/9810065v1](https://arxiv.org/pdf/gr-qc/9810065v1) - Eqs. (30)-(32), zero-shift harmonic slicing and its integrated lapse relation
- [Baumgarte, Hughes, and Shapiro, arXiv:gr-qc/9902024v1](https://arxiv.org/pdf/gr-qc/9902024v1) - p. 2, zero-shift harmonic-slicing relation
- [legacy NRPy tutorial repository](https://github.com/zachetienne/nrpytutorial) - `Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.ipynb`, Step 2.b, `HarmonicSlicing` chain-rule derivation
- [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py) - `BSSNRHSs`, `BSSN_RHSs_varname_to_expr_dict`
- [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py) - `BSSNQuantities`, `BSSN_quantities`
- [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py) - `BSSN_gauge_RHSs`
- [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py) - `BSSNconstraints`, `BSSN_constraints`
- [kreiss_oliger_terms.py](../../../nrpy/equations/general_relativity/kreiss_oliger_terms.py) - `add_KreissOliger_dissipation_terms`
- [kreiss_oliger_terms_Cartesian_CAKO.py](../../../nrpy/equations/general_relativity/tests/kreiss_oliger_terms_Cartesian_CAKO.py) - `trusted_dict`
- [grid.py](../../../nrpy/grid.py) - `register_gridfunctions`
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
- [YBS Gamma-Constraint Adjustment](ybs-gamma-constraint-adjustment.md)
- [YBS-MOM Timestep-Scaled Momentum Adjustment](ybs-momentum-damping.md)
- [Fully Covariant Conformal Z4](fccz4.md)
- [Metric Conversions And Matter](metric-conversions-and-matter.md)
- [Initial Data](initial-data.md)
- [GRHD](../grhd.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
