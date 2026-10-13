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

`BSSN_algebraic_constraints(CoordSystem, enable_rfm_precompute=False)`
returns the corrected rescaled tensors `(hprimeDD, aprimeDD)`. This equation
function selects `_project_general_rfm` for names beginning with `GeneralRFM`
and `_project_orthogonal_rfm` otherwise. BHaH and ETLegacy call it to generate
one combined all-points determinant/trace projector; the backends supply
storage, loops, and schedules.

Claim evidence:
- Claim: `BSSN_algebraic_constraints` defines the symbolic determinant/trace projection and returns corrected rescaled metric and curvature tensors; BHaH and ETLegacy delegate the expressions to it.
- Role: descriptive behavior
- Deciding authority: [BSSN_algebraic_constraints.py](../../../nrpy/equations/general_relativity/BSSN_algebraic_constraints.py), `BSSN_algebraic_constraints`, `_project_general_rfm`, and `_project_orthogonal_rfm`
- Corroboration: [BHaH projector](../../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) and [ETLegacy projector](../../../nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`

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

Claim evidence:
- Claim: BHaH and ETLegacy use a combined determinant/trace projector; the BHaH collision caller applies it after boundary handling, while the ETLegacy projector registration assigns the same combined projector to all three repair schedules and emits no determinant-only repair entry.
- Role: descriptive behavior
- Deciding authority: [BHaH algebraic constraint projection](../../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`; [BHaH initial data](../../../nrpy/infrastructures/BHaH/general_relativity/initial_data.py), `register_CFunction_initial_data`; [two_blackholes_collide.py](../../../nrpy/examples/two_blackholes_collide.py), Method of Lines registration; [ETLegacy algebraic constraint projection](../../../nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- Corroboration: none available; BHaH and ETLegacy own separate backend behavior

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
- Corroboration: [legacy NRPy tutorial repository](https://github.com/zachetienne/nrpytutorial), `Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.ipynb`, Steps 2.b and 2.d

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
connection constraint expressions. It registers scalar diagnostic `H`;
`register_M_and_LAMBDA_CONSTRAINT_gridfunctions=True` also registers `M` and
`LAMBDA_CONSTRAINT` by default. Disabling that option suppresses their storage
registration while retaining the expressions. A separate option registers the
`MU` components. The object stores `H`, `MU`, the internal contraction
`Msquared`, and rescaled `mU`. The internal contraction is the physical norm
square `gamma_ij M^i M^j`, implemented as the conformal-metric contraction
divided by `exp_m4phi`; infrastructure producers export `M` as
`sqrt(Msquared)`.

Claim evidence:
- Claim: `BSSNconstraints` registers `H` and, by default, `M` and `LAMBDA_CONSTRAINT`; the latter two registrations can be disabled with `register_M_and_LAMBDA_CONSTRAINT_gridfunctions=False` while retaining their expressions; its internal `Msquared` is `gamma_ij M^i M^j` with `gamma_ij = gammabar_ij / exp_m4phi`, and the BHaH, ETLegacy, and CarpetX constraint producers export `M = sqrt(Msquared)` rather than the internal square.
- Role: descriptive behavior
- Deciding authority: [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints.__init__`; [BHaH constraints_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/constraints_eval.py), `register_CFunction_constraints_eval`; [ETLegacy BSSN_constraints.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/BSSN_constraints.py), `register_CFunction_BSSN_constraints`; [CarpetX BSSN_constraints.py](../../../nrpy/infrastructures/CarpetX/general_relativity/BSSN_constraints.py), `register_CFunction_BSSN_constraints`
- Corroboration: [BSSN_constraints_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_constraints_Cartesian.py), `trusted_dict`; [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints_dict.__getitem__` physical-ADM-norm doctest; [BHaH constraints_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/constraints_eval.py), existing code-generation doctest

Its covariant conformal connection constraint vector is
`LambdaConstraintU = LambdabarU - DGammaU`.
Brown's connection equation adds `C^k Dhat_k beta^i` to the coordinate-time
RHS. Combined with the Lie stretch, its shift-gradient part is
`-DGammaU[k] Dhat_k betaU[i]`; the change vanishes when `C^i=0`.

Claim evidence:
- Claim: `BSSNconstraints.LambdaConstraintU` is the covariant conformal connection constraint, and `BSSNRHSs` includes Brown's `C^k Dhat_k beta^i` adjustment.
- Role: public scientific equation definition
- Deciding authority: [Brown, *Covariant formulations of BSSN and the standard gauge*, arXiv:0902.3652v2](https://arxiv.org/pdf/0902.3652v2), Eqs. (12a), (12b), (15), (21e), and (22e)
- Corroboration: [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs`; [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints_dict.__getitem__`

### Opt-in YBS Gamma-constraint adjustment

`BSSNRHSs(..., enable_YBS_Gamma_constraint_adjustment=True)` adds

```text
-YBS_chi * (LambdabarU[i] - DGammaU[i]) * Dbarbetacontraction
```

to the unrescaled conformal-connection RHS, where
`Dbarbetacontraction = Dbar_j beta^j` and the parenthesized factor is the
Lambda constraint `C^i`. This is the incremental term for NRPy's existing
equation, not the full paper coefficient: the existing RHS already contains
`+(2/3) DGammaU[i] * Dbarbetacontraction`, which supplies the paper's baseline
`-(2/3) C^i Dbar_j beta^j` contribution when the connection equation is
rewritten in terms of `C^i`. The full Yo-Baumgarte-Shapiro coefficient is
`-(YBS_chi + 2/3) C^i Dbar_j beta^j`; adding that full expression again would
double-count the baseline. Yo-Lin-Cao write the same coefficient as
`-(2/3)(xi + 1)`, so NRPy's incremental parameter satisfies
`YBS_chi = 2*xi/3`. Their usual `xi=1` choice therefore maps to
`YBS_chi = 2/3`. The runtime default is `0.0`, which removes the term exactly;
a user who enables the option sets `YBS_chi` above zero in the parameter file.
The recommended maximum is `2/3` for BSSN. It is the `xi=1` value that Yo, Lin,
and Cao state is usually chosen in Eq. (47), a customary choice and not a
stability bound. The term is lower differential order and changes the local
growth rate of `C^i` by `-YBS_chi*Dbar_j beta^j`. It leaves the principal
characteristic speeds unchanged, but sufficiently large decay rates can
restrict explicit timesteps. For example, classical RK4 applied to
`dC/dt=-kappa*C` has amplification `R(-3)=1.375` when `kappa*Delta t=3`,
despite continuum decay. Larger parameter values have no supporting bound
here, and negative values amplify `C^i` wherever `Dbar_j beta^j > 0`.

Claim evidence:
- Claim: The lower-order YBS Gamma addition leaves principal characteristic speeds unchanged but can impose an explicit timestep restriction through its local decay rate.
- Role: descriptive behavior
- Deciding authority: [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS Gamma addition; [rk_butcher_table_dictionary.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/rk_butcher_table_dictionary.py), `generate_Butcher_tables`, classical RK4 coefficients
- Corroboration: the classical RK4 polynomial `R(z)=1+z+z^2/2+z^3/6+z^4/24` gives the displayed amplification directly.

The equation layer uses a plain real SymPy symbol for `YBS_chi`; constructing
enabled symbolic RHSs does not register a global `CodeParameter`. The BHaH
`register_CFunction_rhs_eval` owner registers `YBS_chi` only when the option is
enabled, for either formulation, as a runtime `REAL`, using the
default `0.0` (term disabled at run time), with `commondata=True` and
`add_to_parfile=True`; there is no Python coefficient argument.
`BSSNRHSs_dict.get_rhs` receives the Boolean
separately from the unchanged coordinate/options string and keeps enabled
expressions in a separate internal cache. With the option false, neither the
term nor BHaH parameter registration is present, and the existing base cache
is used.

`BSSN_gauge_RHSs` forwards the Boolean to the shared BSSN RHS cache. The
derivative-based shift choices `GammaDriving2ndOrder_NoCovariant`,
`GammaDriving2ndOrder_Covariant`,
`GammaDriving2ndOrder_Covariant__Hatted`, and `NonAdvectingGammaDriving`
consume the adjusted Lambda RHS by default. Supplying
`evolved_connection_rhsU` replaces that input with the caller's unrescaled
connection RHS; forwarding the Boolean alone does not adjust a replacement
vector. `Frozen` and both first-order covariant variants do not consume that
RHS and remain unchanged.

Claim evidence:
- Claim: Derivative-based shift equations use the adjusted BSSN connection RHS only when `evolved_connection_rhsU` is absent; an explicit vector supplies the unrescaled connection evolution instead.
- Role: descriptive behavior
- Deciding authority: [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py), `BSSN_gauge_RHSs`
- Corroboration: [fCCZ4_gauge_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_gauge_RHSs.py), `fCCZ4_gauge_RHSs`, supplies the fCCZ4 connection RHS.

With `Theta = Dbar_j beta^j`, the corresponding constraint propagation has
the lower-order addition `-YBS_chi*Theta*C^i`:

```text
partial_t C^i = beta^j Dhat_j C^i
                + 2 alpha exp(4 phi) M^i - YBS_chi Theta C^i.
```

This term attenuates its isolated local factor only where
`YBS_chi*Theta > 0`; it amplifies where that product is negative and is neutral
where it vanishes. Advection and the momentum-constraint source can dominate,
so a globally positive `YBS_chi` is not a uniform damping guarantee when
`Theta` changes sign. Because the addition is lower order, it does not change
the current second-order system's principal symbol or high-frequency
characteristic speeds; it can still change finite-frequency growth.

Claim evidence:
- Claim: the opt-in implementation adds only `-YBS_chi*C^i*Theta`, maps the full-paper parameters by `YBS_chi=2*xi/3`, and provides conditional rather than uniform attenuation with the displayed lower-order propagation term.
- Role: public/scientific contract
- Deciding authority: [Yo, Baumgarte, and Shapiro, arXiv:gr-qc/0209066v2](https://arxiv.org/pdf/gr-qc/0209066v2), Eq. (45); [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2), Eq. (47)
- Corroboration: [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs`; [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py), `BSSN_gauge_RHSs`; [BSSN_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_RHSs_Cartesian.py), representative jointly enabled `trusted_dict`

Claim evidence:
- Claim: the runtime default `YBS_chi = 0.0` removes the term exactly, and the recommended BSSN maximum `2/3` is the usual `xi=1` value stated in the source, a customary choice and not a derived stability bound.
- Role: descriptive behavior
- Deciding authority: [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval`; [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2), Eq. (47)
- Corroboration: [Dendro rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval`; [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS Gamma branch

For local diagnostics, `BSSNconstraints` also stores the direct conformal-metric
contraction `LambdaConstraintSquared` and its plain square root
`LambdaConstraintMagnitude`.

Claim evidence:
- Claim: `BSSNconstraints` computes `LambdaConstraintSquared` by contracting `LambdaConstraintU` twice with `gammabarDD`, then defines `LambdaConstraintMagnitude` as its plain square root; no physical-metric conformal-factor conversion is applied.
- Role: descriptive behavior
- Deciding authority: [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints.__init__`
- Corroboration: [BSSN_constraints_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_constraints_Cartesian.py), `trusted_dict`; [BSSN_constraints_Spherical.py](../../../nrpy/equations/general_relativity/tests/BSSN_constraints_Spherical.py), `trusted_dict`

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

Claim evidence:
- Claim: `test_example_BSSN()` reads module globals for a tuple parser result and otherwise reads the returned object's namespace, uses separate conformal-metric and perturbation declarations plus both partial-derivative replacement forms, and performs only a deterministic sampled-numerical cross-check; this local branch handling does not establish behavior for every upstream parser version.
- Role: descriptive behavior
- Deciding authority: [test_parse_BSSN.py](../../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py), `test_example_BSSN`
- Corroboration: none available; compatibility handling is local harness behavior, and the validation result below exercises only the installed tuple-result path

The style contract for these modules is part of their interface: tensor
construction uses explicit loops, established suffixes such as `U`, `D`, `DD`,
`dD`, `dDD`, `dupD`, and `rhs`, and validation keys must stay aligned with the
corresponding trusted files.

## Sources

- [BSSN_algebraic_constraints.py](../../../nrpy/equations/general_relativity/BSSN_algebraic_constraints.py) - `BSSN_algebraic_constraints`
- [rk_butcher_table_dictionary.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/rk_butcher_table_dictionary.py) - `generate_Butcher_tables`

- [Brown, arXiv:0902.3652v2](https://arxiv.org/pdf/0902.3652v2) - Eqs. (12a), (12b), (15), (21e), and (22e); published as [Phys. Rev. D 79, 104029](https://doi.org/10.1103/PhysRevD.79.104029) (secondary metadata)
- [Yo, Baumgarte, and Shapiro, arXiv:gr-qc/0209066v2](https://arxiv.org/pdf/gr-qc/0209066v2) - Eq. (45), adjusted Gamma-constraint coefficient
- [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2) - Eq. (47), `xi` parameterization of the adjusted coefficient
- [Baumgarte and de Oliveira, arXiv:2201.08857v1](https://arxiv.org/pdf/2201.08857v1) - Eq. (1), Bona-Masso slicing and the `f(alpha)=1` harmonic specialization
- [Baumgarte and Shapiro, arXiv:gr-qc/9810065v1](https://arxiv.org/pdf/gr-qc/9810065v1) - Eqs. (30)-(32), zero-shift harmonic slicing and its integrated lapse relation
- [Baumgarte, Hughes, and Shapiro, arXiv:gr-qc/9902024v1](https://arxiv.org/pdf/gr-qc/9902024v1) - p. 2, zero-shift harmonic-slicing relation
- [legacy NRPy tutorial repository](https://github.com/zachetienne/nrpytutorial) - `Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.ipynb`, Step 2.b, `HarmonicSlicing` chain-rule derivation
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
- [YBS-MOM Timestep-Scaled Momentum Adjustment](ybs-momentum-damping.md)
- [Fully Covariant Conformal Z4](fccz4.md)
- [Metric Conversions And Matter](metric-conversions-and-matter.md)
- [Initial Data](initial-data.md)
- [GRHD](../grhd.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
