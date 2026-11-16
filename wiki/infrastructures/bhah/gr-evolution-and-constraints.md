# GR Evolution, Ricci, And Constraints

> BHaH CFunction registration for BSSN and fCCZ4 RHS kernels, YBS adjustments, Ricci and stored first-derivative helpers, CPU tiling, constraint evaluation, and algebraic projection. · Status: confirmed
> Up: [BHaH](index.md)

## Summary

BHaH's general-relativity evolution modules wrap symbolic expressions from `nrpy/equations/general_relativity/` into interior OpenMP or CUDA loops without rederiving continuum equations. They register `rhs_eval`, `Ricci_eval`, optional stored first-derivative kernels (`hDDdD_eval` and `cfdD_alphadD_vetUdD_eval`), OpenMP tiled `rhs_eval_with_Ricci`, host-side `constraints_eval`, and all-points algebraic constraint projection `enforce_detgbar_equals_detghat_trAzero`.

## Detail

### RHS Registration And YBS Options

`rhs_eval.py::register_CFunction_rhs_eval(..., enable_fCCZ4=False, enable_YBS_Gamma_constraint_adjustment=False, enable_YBS_momentum_constraint_adjustment=False)` assembles BSSN in-module (`BSSN_RHSs` + `BSSN_gauge_RHSs` + Kreiss-Oliger / CAKO, CAHD, and slow-start-lapse terms on copied cached dictionaries) by default, and delegates symbolic RHS/gauge/dissipation/CAHD/SSL assembly to `fCCZ4_system.py::build_fccz4_expression_bundle` when `enable_fCCZ4=True`.
Both branches map sorted RHS names to `rhs_gfs` via `BHaHGridFunction.access_gf` and emit an interior `simple_loop` with `vetU`-controlled upwinding, optional SIMD/CUDA intrinsics, and signature `(const commondata_struct *restrict commondata, const params_struct *restrict params, const rfm_struct *restrict rfmstruct, ...)` when `enable_rfm_precompute=True` versus `(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL *restrict x0, const REAL *restrict x1, const REAL *restrict x2, ...)` when `False`.

When `enable_YBS_Gamma_constraint_adjustment=True`, `register_CFunction_rhs_eval` leaves the coordinate/options string unchanged, registers runtime `CodeParameter` `YBS_chi` (default `0`) after the parallel-registration guard, and forwards the boolean to the symbolic caches (separate internal cache per formulation). `two_blackholes_collide.py`, `blackhole_spectroscopy.py`, `superB_two_blackholes_collide.py`, and `superB_blackhole_spectroscopy.py` pass an explicit `False` constant with no CLI switch.

Whenever `enable_CAHD` or `enable_YBS_momentum_constraint_adjustment` is enabled, `register_CFunction_rhs_eval` registers `AUXEVOL` gridfunction `"dsmin"` (emitted C index macro `DSMINGF`); when `enable_YBS_momentum_constraint_adjustment=True`, it also registers runtime `CodeParameter` `C_YBS_mom` (default `0`, recommended `1.75`), modifying `a_rhsDD` in place without adding evolved fields, initial-data steps, boundary conditions, or KO terms. `blackhole_spectroscopy.py` keeps `enable_YBS_momentum_constraint_adjustment=False` and schedules `dsmin_gf.py::register_CFunction_dsmin_auxevol_gridfunction` whenever either CAHD or YBS-MOM is enabled.

Evidence (descriptive behavior): `BHaH/general_relativity/rhs_eval.py::register_CFunction_rhs_eval`; `fCCZ4_system.py::build_fccz4_expression_bundle`; `BHaH/general_relativity/dsmin_gf.py`; `blackhole_spectroscopy.py`; `two_blackholes_collide.py`.

### Ricci Evaluation, Stored First Derivatives, And CPU Tiling

`Ricci_eval.py::register_CFunction_Ricci_eval` evaluates `BSSN_quantities[CoordSystem + "_rfm_precompute"].Ricci_exprs`, writing `RBARDD*GF` in `auxevol_gfs` or `DIAG_RBARDD*GF` when `host_only_version=True` (which forces OpenMP to emit `Ricci_eval_host` for CUDA builds; `GeneralRFM` rejects CUDA).

`hDDdD_eval.py::register_CFunction_hDDdD_eval` computes centered first derivatives of `hDD` into `SCRATCH` gridfunctions `hDDdD` over the interior grown by `fd_order/2` in directions transverse to each derivative's stencil. Passing `enable_hDDdD_gridfunctions=True` to `Ricci_eval` supplies `stored_first_derivatives` to `c_codegen()`, computing mixed second derivatives as single first-derivative stencils on `hDDdD` while leaving diagonal, upwind, and KO operators and symbolic caches unchanged. `blackhole_spectroscopy.py` reuses `RK_OUTPUT_GFS` as scratch (`hDDdD_eval` -> `Ricci_eval` -> `rhs_eval` overwriting `RK_OUTPUT_GFS`), enabling it only for CUDA double precision and excluding `RK3 Ralston` and `RK3 Heun` (whose second stage keeps all three RK buffers live).

`cfdD_alphadD_vetUdD_eval.py::register_CFunction_cfdD_alphadD_vetUdD_eval` stores the `partial_0` and `partial_1` first derivatives (`j < k` mixed-derivative canonicalization, 10 gridfunctions grown by `fd_order/2` along differentiating axes) of `cf`, `alpha`, and `vetU` in `AUXEVOL` (`NUM_AUXEVOL_GFS` rises from `6` to `16`, ~`61` MB / ~`10%` on standard device + host grids, since stencil reads during `rhs_eval` forbid pointwise `SCRATCH` reuse). CUDA-capable BHaH BSSN generators (`blackhole_spectroscopy.py`, `two_blackholes_collide.py`, `spinning_blackhole.py`, `hydro_without_hydro.py`, `kasner_exact_evolution.py`) expose default-off `enable_cfdD_alphadD_vetUdD_gridfunctions_for_GPU` under `--cuda`.

Both `Ricci_eval` and `rhs_eval` accept `enable_cpu_tiling=False`. In `blackhole_spectroscopy.py`, a single OpenMP predicate (requiring separate Ricci/RHS, `rfm_precompute`, BSSN, and neither stored-derivative option) enables tiled producers and replaces full-grid substep calls with `rhs_eval_with_Ricci` while keeping full-grid kernels for diagnostics. For OpenMP `diagnostic_gfs_set(..., enable_hDDdD_gridfunctions=True)`, temporary scratch is allocated from `y_n_gfs` and freed before `constraints_eval`; CUDA diagnostics use unstored `Ricci_eval_host`.

Evidence (descriptive behavior): `Ricci_eval.py`; `hDDdD_eval.py`; `cfdD_alphadD_vetUdD_eval.py`; `finite_difference.py::select_stored_first_derivatives`; `c_codegen.py::gridfunction_management_and_FD_codegen`.

### Constraint Diagnostics And Algebraic Projection

`constraints_eval.py::register_CFunction_constraints_eval` forces OpenMP, reads `BSSN_constraints[CoordSystem + "_rfm_precompute_RbarDD_gridfunctions" + optional "_T4munu"]`, and writes `H`, `sqrt(gamma_ij M^i M^j)`, and `LAMBDA_CONSTRAINT` into `DIAG_HAMILTONIANGF`, `DIAG_MGF`, and `DIAG_LAMBDA_CONSTRAINTGF` under `diagnostics/`. Under CUDA parallelization, references to `RBARDD*` and `T4UU*` in `auxevol_gfs` are rewritten to `DIAG_*` channels populated on the host.

`enforce_detgbar_equals_detghat_trAzero.py::register_CFunction_enforce_detgbar_equals_detghat_trAzero` emits an all-points loop (taking `(const params_struct *restrict params, const rfm_struct *restrict rfmstruct, ...)` when `enable_rfm_precompute=True` or `(const params_struct *restrict params, const REAL *restrict x0, const REAL *restrict x1, const REAL *restrict x2, ...)` when `False`, plus read-only `auxevol_gfs`) that loads all 12 independent `hDD` and `aDD` components and evaluates `BSSN_algebraic_constraints` before storing: on standard orthogonal coordinate systems (`_project_orthogonal_rfm`) it inverts unrescaled `H = I + h` once (`h' = det(H)^(-1/3) H - I`, `a' = a - H (H^(-1):a) / 3`), whereas on `GeneralRFM*` (`_project_general_rfm`) it rescales `gammabar' = (abs(det(gammahat)) / det(gammabar))^(1/3) gammabar`, inverts `gammabar'`, and subtracts its trace from `Abar`.
Projection is opt-in: `blackhole_spectroscopy.py` (for both BSSN and `--fccz4`) and `two_blackholes_collide.py` (BSSN) pass `enable_conformal_projection=True` to `register_CFunction_initial_data` and append the projector to the Method of Lines `post_rhs_string`.

Evidence (descriptive behavior): `constraints_eval.py`; `enforce_detgbar_equals_detghat_trAzero.py::register_CFunction_enforce_detgbar_equals_detghat_trAzero`; `BSSN_algebraic_constraints.py::BSSN_algebraic_constraints`; `blackhole_spectroscopy.py`; `two_blackholes_collide.py`.

## Sources

- [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py) - `register_CFunction_rhs_eval`, `register_CFunction_rhs_eval_with_Ricci`
- [dsmin_gf.py](../../../nrpy/infrastructures/BHaH/general_relativity/dsmin_gf.py) - `register_CFunction_dsmin_auxevol_gridfunction`
- [Ricci_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/Ricci_eval.py) - `register_CFunction_Ricci_eval`
- [hDDdD_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/hDDdD_eval.py) - `register_CFunction_hDDdD_eval`, `register_hDDdD_gridfunctions`
- [cfdD_alphadD_vetUdD_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/cfdD_alphadD_vetUdD_eval.py) - `register_CFunction_cfdD_alphadD_vetUdD_eval`, `cfdD_alphadD_vetUdD_gridfunction_expressions`
- [constraints_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/constraints_eval.py) - `register_CFunction_constraints_eval`
- [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- [diagnostic_gfs_set.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostic_gfs_set.py) - `register_CFunction_diagnostic_gfs_set`
- [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py) - `BSSNQuantities`
- [BSSN_algebraic_constraints.py](../../../nrpy/equations/general_relativity/BSSN_algebraic_constraints.py) - `BSSN_algebraic_constraints`
- [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py) - `BSSNconstraints`
- [finite_difference.py](../../../nrpy/finite_difference.py) - `select_stored_first_derivatives`
- [c_codegen.py](../../../nrpy/c_codegen.py) - `gridfunction_management_and_FD_codegen`
- [blackhole_spectroscopy.py](../../../nrpy/examples/blackhole_spectroscopy.py) - `--fccz4`, CPU tiling, stored derivatives, projection hooks
- [two_blackholes_collide.py](../../../nrpy/examples/two_blackholes_collide.py) - BSSN RHS, stored derivatives, projection hooks
- [superB_two_blackholes_collide.py](../../../nrpy/examples/superB_two_blackholes_collide.py) - forwarded `enable_YBS_Gamma_constraint_adjustment=False`
- [superB_blackhole_spectroscopy.py](../../../nrpy/examples/superB_blackhole_spectroscopy.py) - forwarded `enable_YBS_Gamma_constraint_adjustment=False`
- [rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_SinhSpherical_RbargfsFalse_T4munuFalse_ImprovementsFalse.py](../../../nrpy/infrastructures/BHaH/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_SinhSpherical_RbargfsFalse_T4munuFalse_ImprovementsFalse.py) - `trusted_dict`
- [fCCZ4_system.py](../../../nrpy/equations/general_relativity/fCCZ4_system.py) - `build_fccz4_expression_bundle`

## See Also

- Depends on: [BSSN Family](../../equations/general-relativity/bssn-family.md)
- Depends on: [fCCZ4 Implementation, Assembly, And Validation](../../equations/general-relativity/fccz4-implementation-and-validation.md)
- Depends on: [YBS-MOM Timestep-Scaled Momentum Adjustment](../../equations/general-relativity/ybs-momentum-damping.md)
