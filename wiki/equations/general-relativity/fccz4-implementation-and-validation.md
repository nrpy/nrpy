# fCCZ4 Implementation, Assembly, And Validation

> NRPy symbolic caching, rescaled tensor storage, YBS adjustments, `fCCZ4_system` RHS assembly, and trusted-expression coverage for fCCZ4. · Status: confirmed
> Up: [General Relativity](index.md)

## Summary

NRPy constructs fCCZ4 symbolic expressions as non-mutating deltas on option-matched cached `BSSNQuantities` and `BSSNRHSs` objects. `fCCZ4_constraints.py` caches `C^i`, spatial Z4 vectors, `RbarZ4_ij`, and `H_Z4`; `fCCZ4_RHSs.py` copies that cache and exports 18 rescaled non-gauge RHS expressions using evolved gridfunction basename `Theta_fCCZ4` while reusing `lambdaU` for `Lambdatilde^i`. `fCCZ4_system.py::build_fccz4_expression_bundle` combines non-gauge RHSs, gauge RHSs, Kreiss-Oliger dissipation, CAHD, and slow-start lapse for BHaH, whereas Dendro wires fCCZ4 directly.

## Detail

### Cached BSSN Reuse And YBS Adjustments

`FCCZ4Constraints` (`FCCZ4ConstraintsDict` / `fCCZ4_constraints`) caches the fCCZ4 constraint, spatial-Z4, Z4-Ricci, and `H_Z4` expressions, adding `T4munu` matter terms only when `enable_T4munu=True`. `FCCZ4RHSs` (`FCCZ4RHSsDict` / `fCCZ4_RHSs`) copies those tensor containers without mutating `BSSNQuantities` or `BSSNRHSs` and forms rescaled RHSs via componentwise division by `ReDD_ij` (not a contraction):

```text
(partial_t h_ij)_fCCZ4 = (partial_t h_ij)_BSSN   [with tr(Abar) = 0 imposed analytically],
(partial_t a_ij)_fCCZ4 = (partial_t a_ij)_BSSN
                         + [-2 alpha Theta Abar_ij + alpha exp(-4 phi) (deltaR_ij)^TF] / ReDD_ij,
(partial_t K)_fCCZ4    = (partial_t K)_BSSN + alpha(H_Z4 - 2 Theta K) - 3 kappa1(1 + kappa2) Theta.
```

The conformal-factor RHS is reused unchanged, `deltaLambda^i` is added before `ReU_i` rescaling, and matter terms reuse `T4munu.py` (`BSSN_RHSs_T4UU_source_terms`, `BSSN_constraints_T4UU_source_terms`).

With `enable_YBS_Gamma_constraint_adjustment=True`, `FCCZ4RHSs` requests the YBS-adjusted `BSSNRHSs` base in a separate internal cache (leaving the coordinate/options key string unchanged). Because `lambdaU` stores `LambdatildeU`, the shared BSSN term `-YBS_chi * (LambdatildeU[i] - DGammaU[i]) * Dbar_j beta^j` evaluates to `-YBS_chi * C^i * Dbar_j beta^j` (`Dbar_j beta^j = Dhat_j beta^j` when `det(gammabar) = det(gammahat)`). Combined with `+(2/3) C^i Dbar_j beta^j` in `deltaLambda^i`, the net propagation coefficient is `2/3 - YBS_chi` (zero at `YBS_chi = 2/3`; runtime default `YBS_chi = 0`; recommended maximum `YBS_chi = 4/3` to match BSSN `xi = 1` coefficient `-2/3`, not a stability bound). Enabling `enable_YBS_momentum_constraint_adjustment=True` modifies `a_rhsDD` directly using the lower-index matter-complete momentum residual without exposing `MU` or adding evolved gridfunctions.

Evidence (descriptive behavior): `fCCZ4_constraints.py::FCCZ4Constraints.__init__`; `fCCZ4_RHSs.py::FCCZ4RHSs.__init__`; `BSSN_RHSs.py::BSSNRHSs.__init__`; [Yo, Lin, and Cao](https://arxiv.org/pdf/1205.5111v2) Eq. (47).

### Reference-Rescaled Storage, Gridfunctions, And External Ricci

NRPy evolves reference-rescaled tensor and vector components:

```text
gammabar_ij   = gammahat_ij + ReDD_ij hDD_ij,   Abar_ij = ReDD_ij aDD_ij,
Lambdatilde^i = ReU_i lambdaU_i,   beta^i = ReU_i vetU_i,   B^i = ReU_i betU_i.
```

`FCCZ4RHSs` exposes Python attributes `Theta`, `Theta_dD` (centered), `Theta_dupD` (upwinded), `Theta_rhs`, `Z4constraintU`, `ZbarU`, `ZU`, `ZD`, `LambdatildeU`, `RbarZ4DD`, `RbarZ4`, `H_Z4`, and rescaled/unrescaled RHS arrays, while gridfunction storage and derivative basenames use `Theta_fCCZ4`, `Theta_fCCZ4_dD`, `Theta_fCCZ4_dupD`, and `Theta_fCCZ4_rhs` with no generic `Theta` gridfunction alias. `fCCZ4_RHSs_varname_to_expr_dict` always holds 18 sorted entries (`hDD` x6, `aDD` x6, `lambdaU` x3, conformal factor, `trK`, `Theta_fCCZ4`) whether YBS-MOM is off or on.

Enabling `enable_T4munu=True` incrementally registers only symmetric `T4UU` `AUXEVOL` gridfunctions plus the `PI` `CodeParameter`, without registering `H`, `M`, or `LAMBDA_CONSTRAINT`. When `enable_RbarDD_gridfunctions=True`, `BSSNQuantities_dict` registers symmetric `RbarDD` in `AUXEVOL` (`auxevol_gfs`) if absent, and `FCCZ4RHSs` treats `RbarDD` as `RbarNRPy_ij(Lambdatilde)` (evolved `lambdaU` in the linear derivative slot, geometric `DeltaGamma^i_jk` in nonlinear products) and adds `deltaR_ij`; callers must still generate the kernel that populates `RbarDD`.

Evidence (descriptive behavior): `BSSN_quantities.py::BSSNQuantities.__init__`; `fCCZ4_RHSs.py::FCCZ4RHSs`; `fCCZ4_constraints.py::FCCZ4Constraints`; `T4munu.py::BSSN_constraints_T4UU_source_terms`.

### Complete RHS Assembly And Application Wiring

`fCCZ4_system.py::build_fccz4_expression_bundle` copies `fCCZ4_RHSs_varname_to_expr_dict`, appends lapse and shift RHSs from `fCCZ4_gauge_RHSs`, and optionally adds Kreiss-Oliger dissipation (curvature-aware KO requires KO enabled; see [BSSN Family](bssn-family.md)), CAHD, and slow-start lapse. It returns frozen dataclass `FCCZ4ExpressionBundle` (`rhs_by_symbol_name`, shift `upwind_control_vec`, opt-in `diagnostics_by_name` holding `H_Z4` and `Z4constraintU0..2`, and `provenance_files`; internal dicts remain mutable).

CAHD registers `CodeParameter` `C_CAHD` (default `0.15`) and adds `-C_CAHD * CFL_FACTOR * dsmin * H` to `phi_rhs` (multiplied by `-2 W` for `W_rhs` or `-4 chi` for `chi_rhs`), where `H` is `BSSN_constraints` Hamiltonian built with evolved `lambdaU` rather than `H_Z4`. Instantiating `BSSN_constraints` registers diagnostic gridfunction `H` and, under default NRPy parameter `register_M_and_LAMBDA_CONSTRAINT_gridfunctions=True`, `M` and `LAMBDA_CONSTRAINT` even when `diagnostics_by_name` is disabled.

BHaH `rhs_eval.py::register_CFunction_rhs_eval` calls `build_fccz4_expression_bundle` when `enable_fCCZ4=True`, whereas Dendro `rhs_eval.py` wires fCCZ4 directly; for `blackhole_spectroscopy.py --fccz4` checkpoint-first `Theta_fCCZ4 = 0` initialization see [BHaH GR Initial Data, Diagnostics, And Psi4](../../infrastructures/bhah/gr-initial-data-and-psi4.md), and for shared algebraic projection (`enforce_detgbar_equals_detghat_trAzero`) see [BHaH GR Evolution, Ricci, And Constraints](../../infrastructures/bhah/gr-evolution-and-constraints.md).

Evidence (descriptive behavior): `fCCZ4_system.py::build_fccz4_expression_bundle`; `BSSN_constraints.py::BSSNconstraints.__init__`; `BHaH/general_relativity/rhs_eval.py`; `Dendro/general_relativity/rhs_eval.py`.

### Trusted-Expression Validation Coverage And Gaps

`fCCZ4_constraints.py`, `fCCZ4_RHSs.py`, and `fCCZ4_gauge_RHSs.py` run `doctest.testmod()` (zero inline doctests) followed by `trusted_dict` checks via `validate_expressions.py`; `fCCZ4_system.py` runs assembly doctests plus Cartesian `trusted_dict` comparisons with KO off and on (`CAHD` and both YBS flags off):

- `fCCZ4_constraints`: `W`, `phi`, and `chi` in `SinhCartesian` and `SinhSpherical`, `RbarDD` gridfunctions in both, `SinhSpherical` `rfm_precompute` + `T4munu`, and `Cartesian`.
- `fCCZ4_RHSs` (both YBS flags jointly enabled): `SinhCartesian` baseline and `RbarDD`, `SinhSpherical` baseline, `RbarDD`, and `rfm_precompute` + `T4munu`, plus `Cartesian`.
- `fCCZ4_gauge_RHSs` (both YBS flags disabled): full `4 x 7` lapse/shift matrices in `SinhCartesian` and `SinhSpherical`, plus `Cartesian` default and `SinhSpherical` `rfm_precompute` + `T4munu`.

Because `BHaH/general_relativity/rhs_eval.py` `__main__` never sets `enable_fCCZ4=True`, its joint-YBS `trusted_dict` modules pin BSSN only. Consequently, no `trusted_dict` pins the YBS-Gamma shift-driver change in `fCCZ4_gauge_RHSs` or `build_fccz4_expression_bundle` with either YBS flag enabled; fCCZ4 YBS-MOM on `a_rhsDD` is pinned only in `fCCZ4_RHSs` tests.

Evidence (descriptive behavior): `fCCZ4_constraints.py`; `fCCZ4_RHSs.py`; `fCCZ4_gauge_RHSs.py`; `fCCZ4_system.py`; `BHaH/general_relativity/rhs_eval.py`; `validate_expressions.py::compare_or_generate_trusted_results`.

## Sources

- [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py) - `FCCZ4Constraints`, `FCCZ4ConstraintsDict`, `fCCZ4_constraints`, `__main__`
- [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py) - `FCCZ4RHSs`, `FCCZ4RHSsDict`, `fCCZ4_RHSs`, `__main__`
- [fCCZ4_gauge_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_gauge_RHSs.py) - `fCCZ4_gauge_RHSs`, `__main__`
- [fCCZ4_system.py](../../../nrpy/equations/general_relativity/fCCZ4_system.py) - `FCCZ4ExpressionBundle`, `build_fccz4_expression_bundle`, `__main__`
- [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py) - `BSSNQuantities`, `BSSNQuantities_dict`
- [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py) - `BSSNRHSs`, `BSSNRHSs_dict`
- [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py) - `BSSNconstraints`
- [T4munu.py](../../../nrpy/equations/general_relativity/T4munu.py) - `BSSN_RHSs_T4UU_source_terms`, `BSSN_constraints_T4UU_source_terms`
- [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2) - Eq. (47)
- [BHaH rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py) - `register_CFunction_rhs_eval`, BSSN-only `__main__`
- [Dendro rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - `register_CFunction_rhs_eval`
- [validate_expressions.py](../../../nrpy/validate_expressions/validate_expressions.py) - `process_dictionary_of_expressions`, `compare_or_generate_trusted_results`
- [fCCZ4_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_Cartesian.py) - `trusted_dict`
- [fCCZ4_RHSs_SinhCartesian_RbarDD_gridfunctions.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhCartesian_RbarDD_gridfunctions.py) - `trusted_dict`
- [fCCZ4_RHSs_SinhSpherical_RbarDD_gridfunctions.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhSpherical_RbarDD_gridfunctions.py) - `trusted_dict`
- [fCCZ4_RHSs_SinhSpherical_rfm_precompute_T4munu.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhSpherical_rfm_precompute_T4munu.py) - `trusted_dict`
- [fCCZ4_constraints_SinhSpherical_phi.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_constraints_SinhSpherical_phi.py) - `trusted_dict`
- [fCCZ4_gauge_RHSs_Frozen_Frozen_SinhSpherical.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_gauge_RHSs_Frozen_Frozen_SinhSpherical.py) - `trusted_dict`
- [fCCZ4_gauge_RHSs_OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted_Cartesian.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_gauge_RHSs_OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted_Cartesian.py) - `trusted_dict`
- [fCCZ4_gauge_RHSs_OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted_SinhSpherical_rfm_precompute_T4munu.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_gauge_RHSs_OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted_SinhSpherical_rfm_precompute_T4munu.py) - `trusted_dict`
- [fCCZ4_system_Cartesian.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_system_Cartesian.py) - `trusted_dict`
- [fCCZ4_system_Cartesian_KO.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_system_Cartesian_KO.py) - `trusted_dict`
- [rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_SinhSpherical_RbargfsTrue_T4munuTrue_ImprovementsTrue.py](../../../nrpy/infrastructures/BHaH/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_SinhSpherical_RbargfsTrue_T4munuTrue_ImprovementsTrue.py) - `trusted_dict`

## See Also

- Depends on: [fCCZ4 Formulation And Gauge Equations](fccz4-formulation.md)
- Depends on: [BSSN Family](bssn-family.md)
- Depends on: [YBS-MOM Timestep-Scaled Momentum Adjustment](ybs-momentum-damping.md)
