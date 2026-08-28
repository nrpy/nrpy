# ETLegacy GR BSSN RHS, Ricci, Constraints, And Validation

> ETLegacy registration path for generated BSSN Ricci, RHS, constraints, and RHS trusted-expression evidence. · Status: confirmed · Last reconciled: 08-28-2026
> Up: [ETLegacy](index.md)

## Summary

This page owns the ETLegacy CFunction registration layer for the heavy BSSN
runtime kernels. Symbolic BSSN equation construction stays with
[BSSN Family](../../equations/general-relativity/bssn-family.md), and common
trusted-value mechanics stay with
[Trusted Expression Pipeline](../../equations/trusted-expression-pipeline.md).

`register_CFunction_Ricci_eval()` emits the Ricci tensor auxiliary kernel and
schedules it before the BSSN RHS kernel. `register_CFunction_rhs_eval()` emits
the evolved RHS kernel, including gauge RHSs and optional improvement terms.
`register_CFunction_BSSN_constraints()` emits Hamiltonian, momentum, and
conformal connection diagnostics in `MoL_PseudoEvolution`, including scalar
momentum and Lambda-constraint magnitudes.

## Detail

`register_CFunction_Ricci_eval()` preserves the caller's finite-difference
order, sets the local `fd_order` for parallel codegen, selects
`BSSN_quantities` with the optional `_rfm_precompute` key suffix, and writes
`RbarDD` auxiliary gridfunctions from `Bq.Ricci_exprs`. Its generated body uses
ETLegacy includes, optional SIMD intrinsics and SIMD inverse grid spacings,
finite-difference codegen, finite-difference helper functions, Golden Kernels,
and an interior `simple_loop`. The registered schedule is guarded by
`if(fd_order == <order>)` and places the Ricci function in `MoL_CalcRHS` as
`<thorn>_Ricci before <thorn>_RHS`, reading `hDD` and `lambdaU` gridfunctions
and writing the six `RbarDD` gridfunctions.

`register_CFunction_rhs_eval()` builds the generated BSSN RHS function for one
thorn, coordinate system, finite-difference order, gauge choice, and option
set. It selects `BSSN_RHSs` with key suffixes for reference-metric precompute,
`RbarDD` gridfunction use, and optional T4munu source terms, then calls
`BSSN_gauge_RHSs()` and appends `alpha_rhs`, `vet_rhsU*`, and `bet_rhsU*` to
the local RHS dictionary. That local dictionary is sorted before later
modification so validation and generated output use deterministic expression
names.

The RHS registration wires the major ETLegacy options directly into generated
code. SIMD mode adds the SIMD header, SIMD `invdxx*` variables, and a SIMD read
for `eta`. After option handling, the RHS and constraint registrations discover
registered non-commondata CodeParameters from their final expression lists.
SIMD functions read those discovered parameters through
`read_CodeParameters()`, while scalar functions use
`DECLARE_CCTK_PARAMETERS`. The discovered list also drives current-thorn
CFunction parameter metadata. Non-SIMD RHS code additionally defines scalar
inverse grid spacings and `UPWIND_ALG`.

This pipeline does not branch on `PI`, thorn name, or T4munu state. Matter
expressions incidentally contribute `PI` because their source terms reference
the registered non-commondata `PI` CodeParameter; vacuum expressions do not.

Claim evidence:
- Claim: ETLegacy BSSN RHS and constraint registrations discover non-commondata CodeParameters from final expressions, use that list for current-thorn metadata and SIMD parameter declarations, use Cactus parameter macros for scalar declarations, and contain no `PI`-, thorn-name-, or T4munu-specific parameter-selection branch.
- Role: descriptive behavior
- Deciding authority: [rhs_eval.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval`; [BSSN_constraints.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/BSSN_constraints.py), `register_CFunction_BSSN_constraints`
- Corroboration: [expression_utils.py](../../../nrpy/helpers/expression_utils.py), `get_params_commondata_symbols_from_expr_list`; [CodeParameters.py](../../../nrpy/infrastructures/ETLegacy/CodeParameters.py), `read_CodeParameters`
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3; backend=ETLegacy; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=Cartesian, fd_order 4, RHS KO enabled, reference-metric precompute disabled, T4munu False/True, SIMD False/True; date=08-18-2026`

Reference-metric precompute is threaded through the BSSN and reference-metric
lookup keys. Kreiss-Oliger
dissipation adds `_dKOD` derivative terms to gauge and non-gauge RHSs; CAKO
uses separate gauge and non-gauge strength parameters multiplied by the
conformal-factor-derived `W`, while the non-CAKO path uses the shared
`diss_strength` parameter name. CAHD builds the matching `BSSN_constraints`
object, registers `C_CAHD` and
`CFL_FACTOR__ignore_repeats_Carpet_timeref_factors`, computes a `dsmin`
prefactor, and adds a conformal-factor-adjusted Hamiltonian-damping term to
`cf_rhs`. SSL registers `SSL_h` and `SSL_sigma`, computes a time-dependent
Gaussian prefactor, and modifies `alpha_rhs`.

After option handling, `register_CFunction_rhs_eval()` maps RHS expression
names onto ETLegacy RHS gridfunction names and constructs the upwind control
vector from `vetU[i] * rfm.ReU[i]`. The final `c_codegen()` call enables
finite-difference codegen, optional SIMD, finite-difference helper functions,
Golden Kernels, and `upwind_control_vec=betaU`. The schedule is also guarded by
the chosen `fd_order`; it places `<thorn>_RHS` in `MoL_CalcRHS after
<thorn>_Ricci`, reads `evol_variables(everywhere)` and
`auxevol_variables(interior)`, and writes `evol_variables_rhs(interior)`.
Parameter metadata is passed through `ET_current_thorn_CodeParams_used`, with
`eta` and `fd_order` always listed. The current source then lists
`diss_strength_gauge` and `diss_strength_nongauge` when `enable_CAKO` is true,
and lists `diss_strength` otherwise, even if
`enable_KreissOliger_dissipation` is false. SSL and CAHD parameters are added
only when enabled. The expression-derived non-commondata parameters described
above are appended to that manual metadata list.

`register_CFunction_BSSN_constraints()` emits the diagnostic constraints
kernel. It selects `BSSN_constraints` with optional reference-metric precompute
and T4munu suffixes, outputs `H`, `MU0`, `MU1`, `MU2`,
`M = sqrt(gamma_ij M^i M^j)`, and
`LAMBDA_CONSTRAINT = sqrt(gammabar_ij C^i C^j)`, and runs finite-difference
codegen with finite-difference helper functions and Golden Kernels in an
interior `simple_loop`. Its schedule is guarded by the requested
finite-difference order and places `<thorn>_BSSN_constraints` in
`MoL_PseudoEvolution`, reading BSSN state and optional T4munu gridfunctions and
writing `aux_variables`.

Claim evidence:
- Claim: ETLegacy `register_CFunction_BSSN_constraints` writes `H`, `MU0` through `MU2`, `M = sqrt(BSSNconstraints.Msquared)`, and `LAMBDA_CONSTRAINT = BSSNconstraints.LambdaConstraintMagnitude` to auxiliary gridfunctions.
- Role: descriptive behavior
- Deciding authority: [BSSN_constraints.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/BSSN_constraints.py), `register_CFunction_BSSN_constraints`
- Corroboration: [core BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints.__init__`; [interface_ccl.py](../../../nrpy/infrastructures/ETLegacy/interface_ccl.py), `construct_interface_ccl`
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=ETLegacy registration; precision=symbolic code generation; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=Cartesian, fd_order 4, T4munu False/True, SIMD registration exercised; date=08-28-2026`

All three generated kernels replace the finite-difference helper prefunc text
`NO_INLINE` with `CCTK_ATTRIBUTE_NOINLINE` before registration. The local
comment says this avoids a higher-order finite-difference compile hang with
some GCC versions without changing the shared finite-difference helper for
other infrastructures.

RHS validation remains part of this page because `rhs_eval.py` validates the
ETLegacy-specific assembled RHS dictionary after ETLegacy option handling and
before final C code emission. With `validate_expressions=True`, the module
passes `local_BSSN_RHSs_varname_to_expr_dict` through
`process_dictionary_of_expressions(..., fixed_mpfs_for_free_symbols=True)`.
With `return_validation_dict=True`, it returns that processed dictionary.
Otherwise, the in-function validation path compares or generates a trusted file
whose basename includes the lapse option, shift option, coordinate system,
`T4munu{enable_T4munu}`, `KO{enable_KreissOliger_dissipation}`, and
`improvements{enable_SSL}`. The module `__main__` path drives the covariant
Cartesian OnePlusLog/Gamma-driver cases for both T4munu states and both
all-improvements states; because it requests `return_validation_dict=True` and
then calls `compare_or_generate_trusted_results()` itself, its generated
basenames omit the `KO...` segment and use the local `enable_improvements`
loop variable.

ETLegacy has six backend-local RHS trusted baselines: four covariant cases span
both `T4munu` states and both improvements states, while two noncovariant
KO-enabled cases span both `T4munu` states with improvements disabled.

## Sources

- [rhs_eval.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/rhs_eval.py) - `register_CFunction_rhs_eval`, `validate_expressions`, `__main__`
- [Ricci_eval.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/Ricci_eval.py) - `register_CFunction_Ricci_eval`
- [BSSN_constraints.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/BSSN_constraints.py) - `register_CFunction_BSSN_constraints`
- [CodeParameters.py](../../../nrpy/infrastructures/ETLegacy/CodeParameters.py) - `read_CodeParameters`
- [expression_utils.py](../../../nrpy/helpers/expression_utils.py) - `get_params_commondata_symbols_from_expr_list`
- [rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuFalse_improvementsFalse.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuFalse_improvementsFalse.py) - `trusted_dict`
- [rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuFalse_improvementsTrue.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuFalse_improvementsTrue.py) - `trusted_dict`
- [rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuTrue_improvementsFalse.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuTrue_improvementsFalse.py) - `trusted_dict`
- [rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuTrue_improvementsTrue.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuTrue_improvementsTrue.py) - `trusted_dict`
- [rhs_eval_OnePlusLog_GammaDriving2ndOrder_NoCovariant_Cartesian_T4munuFalse_KOTrue_improvementsFalse.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_NoCovariant_Cartesian_T4munuFalse_KOTrue_improvementsFalse.py) - `trusted_dict`
- [rhs_eval_OnePlusLog_GammaDriving2ndOrder_NoCovariant_Cartesian_T4munuTrue_KOTrue_improvementsFalse.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_NoCovariant_Cartesian_T4munuTrue_KOTrue_improvementsFalse.py) - `trusted_dict`

## See Also

- [ETLegacy](index.md)
- [MoL, Boundaries, Symmetry, And RHS Initialization](mol-boundaries-symmetry-and-rhs-initialization.md)
- [GR ADM/BSSN, Slicing, And Matter Coupling](gr-adm-bssn-slicing-and-matter-coupling.md)
- Depends on: [BSSN Family](../../equations/general-relativity/bssn-family.md)
- Depends on: [Symbolic Expression Utilities](../../core/helpers/symbolic-expression-utilities.md)
- [Trusted Expression Pipeline](../../equations/trusted-expression-pipeline.md)
- [Finite Difference](../../core/finite-difference.md)
- Contrasts with: [CarpetX GR BSSN RHS, Ricci, Constraints, And Validation](../carpetx/gr-bssn-rhs-ricci-constraints-and-validation.md)
