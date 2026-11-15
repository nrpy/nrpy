# BSSN Family

> Reference-metric BSSN quantities, evolution and gauge RHSs, algebraic determinant/trace projection, constraints, Kreiss-Oliger dissipation, and `nrpylatex` cross-checks. · Status: confirmed
> Up: [General Relativity](index.md)

## Summary

NRPy's reference-metric BSSN modules split into cached quantities (`BSSN_quantities.py`), non-gauge evolution equations (`BSSN_RHSs.py`), lapse and shift drivers (`BSSN_gauge_RHSs.py`), algebraic projection (`BSSN_algebraic_constraints.py`), constraint diagnostics (`BSSN_constraints.py`), and Kreiss-Oliger dissipation (`kreiss_oliger_terms.py`). Continuum RHSs impose `det(gammabar) = det(gammahat)` and `gammabar^ij Abar_ij = 0` analytically, leaving numerical projection to backend callers.

## Detail

### `BSSNQuantities`, `BSSNRHSs`, And Algebraic Projection

`BSSNQuantities` registers BSSN gridfunctions on first construction and re-declares identical symbols on later cache keys (`real=True` for scalars `trK`, `cf`, `alpha`; assumption-free symbols for tensor components so expressions across `BSSN_quantities` keys share symbols), building rescaled conformal metrics, connections, derivatives, and sorted `Ricci_exprs`. `BSSNRHSs` reads `BSSNQuantities` and exports `cf_rhs`, `trK_rhs`, `Lambdabar_rhsU`, `h_rhsDD`, `a_rhsDD`, `lambda_rhsU`, and `BSSN_RHSs_varname_to_expr_dict` using explicit indexed loops (`U`/`D`/`DD`/`dD`/`dDD`/`dupD`/`rhs` suffixes). Brown's covariant connection adjustment adds `C^k Dhat_k beta^i` (`C^i = Lambdabar^i - DeltaGamma^i`), reducing the shift-gradient term to `-DeltaGamma^k Dhat_k beta^i`; optional `enable_YBS_Gamma_constraint_adjustment=True` adds `-YBS_chi * (LambdabarU[i] - DGammaU[i]) * Dbarbetacontraction` ([YBS Gamma-Constraint Adjustment](ybs-gamma-constraint-adjustment.md)).

Because `BSSNQuantities` omits `trAbar` and `BSSNRHSs` omits off-constraint `alpha tr(Abar)` terms, `BSSN_algebraic_constraints(CoordSystem, enable_rfm_precompute=False)` returns projected `(hprimeDD, aprimeDD)` via `_project_general_rfm` (`GeneralRFM*`: explicit determinant-rescaled metric inversion and trace subtraction) or `_project_orthogonal_rfm` (orthogonal reference metrics: `H = I + h` inversion with algebraic determinant cancellation). Generated BHaH (`opt-in` in `initial_data` and MoL `post_rhs_string`, e.g., `two_blackholes_collide.py`) and ETLegacy (`CCTK_INITIAL`, `MoL_PostStep`, `MoL_PseudoEvolution`) kernels load all 12 `hDD`/`aDD` components before storing all 12 projected values to roundoff.

Evidence (scientific/interface requirement): `BSSN_quantities.py`; `BSSN_RHSs.py`; `BSSN_algebraic_constraints.py`; [Brown](https://arxiv.org/pdf/0902.3652v2) Eqs. (12a)-(12b), (15), (21e), (22e).

### Lapse And Shift Gauge Options (`BSSN_gauge_RHSs`)

`BSSN_gauge_RHSs` returns `(alpha_rhs, vet_rhsU, bet_rhsU)` across shift options (`Frozen`, first/second-order covariant/hatted/non-covariant Gamma drivers, `NonAdvectingGammaDriving`) and four lapse options (`OnePlusLog`, `BHSHarmonicSlicing`, `Frozen`, `OnePlusLogAlt`; obsolete `"HarmonicSlicing"` raises `ValueError`).

`OnePlusLogAlt` implements `partial_0 alpha = -alpha (1 - alpha) K`. `BHSHarmonicSlicing` implements Baumgarte-Hughes-Shapiro `exp(6 phi)`-tracking (`partial_t alpha = partial_t exp(6 phi)`, i.e., `alpha_rhs = -3 W^(-4) W_rhs` for `W` or `6 exp(6 phi) phi_rhs` for `phi`; `chi` raises `ValueError`). Whereas general Bona-Masso harmonic slicing `(partial_t - beta^i partial_i) alpha = -alpha^2 K` integrates at zero shift to `alpha = C(x) exp(6 phi)` and matches `partial_t alpha = partial_t exp(6 phi)` only for `C(x) = 1`, NRPy's `BHSHarmonicSlicing` preserves the additive difference `alpha - exp(6 phi)` ([Initial Data](initial-data.md) documents its StaticTrumpet pairing).

Evidence (scientific/interface requirement): `BSSN_gauge_RHSs.py`; [Baumgarte and de Oliveira](https://arxiv.org/pdf/2201.08857v1) Eq. (1); [Baumgarte and Shapiro](https://arxiv.org/pdf/gr-qc/9810065v1) Eqs. (30)-(32); [Baumgarte, Hughes, and Shapiro](https://arxiv.org/pdf/gr-qc/9902024v1) p. 2.

### Constraints, Kreiss-Oliger Dissipation, And Validation

`BSSNconstraints` (`BSSN_constraints` cache, rebuilt whenever NRPy parameters `register_M_and_LAMBDA_CONSTRAINT_gridfunctions`—default `True`—or `register_MU_gridfunctions` change) registers diagnostic `H` (and optionally `M`, `LAMBDA_CONSTRAINT`, `MU`) and exposes `H`, `MU`, rescaled `mU`, physical norm square `Msquared = gamma_ij M^i M^j` (conformal contraction divided by `exp(-4 phi)`, exported as `M = sqrt(Msquared)` in BHaH/ETLegacy/CarpetX), and unrescaled conformal `LambdaConstraintU = LambdabarU - DGammaU`, `LambdaConstraintSquared`, and `LambdaConstraintMagnitude`.

`kreiss_oliger_terms.py::add_KreissOliger_dissipation_terms` mutates an RHS dict in place (called after gauge RHS insertion by BHaH, Dendro, and `fCCZ4_system.py` with `include_Theta_fCCZ4=True`), registering `KreissOliger_strength_gauge` (`alpha`, `vetU`, `betU`, skipping any gauge entry equal to exact SymPy `0` such as `Frozen` or first-order `betU`) and `KreissOliger_strength_nongauge` (`cf`, `trK`, `lambdaU`, `aDD`, `hDD`, optional `Theta_fCCZ4`), adding `strength * Q_dKOD[k] * ReU[k]` (multiplied by `W` when `enable_CAKO=True`).

`trusted_dict` suites pin `BSSN_quantities`, `BSSN_constraints`, `BSSN_gauge_RHSs`, `kreiss_oliger_terms` (nonzero gauge RHSs), and `BSSN_RHSs` (with YBS-Gamma and YBS-MOM jointly enabled); `nrpylatex/test_parse_BSSN.py::test_example_BSSN` cross-checks LaTeX-parsed Cartesian BSSN equations against handwritten expressions at fixed `mpf` sample points via `validate_expressions.assert_equal`.

Evidence (descriptive behavior): `equations/general_relativity/BSSN_constraints.py`; `kreiss_oliger_terms.py`; `nrpylatex/test_parse_BSSN.py::test_example_BSSN`.

## Sources

- [Brown, arXiv:0902.3652v2](https://arxiv.org/pdf/0902.3652v2) - Eqs. (12a), (12b), (15), (21e), (22e)
- [Phys. Rev. D 79, 104029](https://doi.org/10.1103/PhysRevD.79.104029) - published Brown reference
- [Baumgarte and de Oliveira, arXiv:2201.08857v1](https://arxiv.org/pdf/2201.08857v1) - Eq. (1)
- [Baumgarte and Shapiro, arXiv:gr-qc/9810065v1](https://arxiv.org/pdf/gr-qc/9810065v1) - Eqs. (30)-(32)
- [Baumgarte, Hughes, and Shapiro, arXiv:gr-qc/9902024v1](https://arxiv.org/pdf/gr-qc/9902024v1) - p. 2
- [legacy NRPy tutorial repository](https://github.com/zachetienne/nrpytutorial) - `Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.ipynb`
- [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py) - `BSSNQuantities`, `BSSN_quantities`
- [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py) - `BSSNRHSs`, `BSSN_RHSs_varname_to_expr_dict`
- [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py) - `BSSN_gauge_RHSs`
- [BSSN_algebraic_constraints.py](../../../nrpy/equations/general_relativity/BSSN_algebraic_constraints.py) - `BSSN_algebraic_constraints`
- [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py) - `BSSNconstraints`, `BSSN_constraints`
- [kreiss_oliger_terms.py](../../../nrpy/equations/general_relativity/kreiss_oliger_terms.py) - `add_KreissOliger_dissipation_terms`
- [grid.py](../../../nrpy/grid.py) - `register_gridfunctions`
- [ETLegacy BSSN_constraints.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/BSSN_constraints.py) - `register_CFunction_BSSN_constraints`
- [CarpetX BSSN_constraints.py](../../../nrpy/infrastructures/CarpetX/general_relativity/BSSN_constraints.py) - `register_CFunction_BSSN_constraints`
- [BHaH constraints_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/constraints_eval.py) - `register_CFunction_constraints_eval`
- [BHaH enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- [ETLegacy enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- [BHaH initial_data.py](../../../nrpy/infrastructures/BHaH/general_relativity/initial_data.py) - `register_CFunction_initial_data`
- [two_blackholes_collide.py](../../../nrpy/examples/two_blackholes_collide.py) - projector registration
- [BSSN_quantities_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_quantities_Cartesian.py) - `trusted_dict`
- [BSSN_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_RHSs_Cartesian.py) - `trusted_dict`
- [BSSN_constraints_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_constraints_Cartesian.py) - `trusted_dict`
- [BSSN_constraints_Spherical.py](../../../nrpy/equations/general_relativity/tests/BSSN_constraints_Spherical.py) - `trusted_dict`
- [kreiss_oliger_terms_Cartesian_CAKO.py](../../../nrpy/equations/general_relativity/tests/kreiss_oliger_terms_Cartesian_CAKO.py) - `trusted_dict`
- [test_parse_BSSN.py](../../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py) - `test_example_BSSN`
- [validate_expressions.py](../../../nrpy/validate_expressions/validate_expressions.py) - `assert_equal`

## See Also

- See also: [fCCZ4 Formulation And Gauge Equations](fccz4-formulation.md)
- See also: [YBS Gamma-Constraint Adjustment](ybs-gamma-constraint-adjustment.md)
- Validated by: [Expression Validation Helpers](../../validation/expression-validation-helpers.md)
