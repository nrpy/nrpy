# CarpetX GR BSSN RHS, Ricci, Constraints, And Validation

> CarpetX registration path for generated BSSN Ricci, RHS, constraints, and RHS trusted-expression evidence. · Status: confirmed · Last reconciled: 07-20-2026
> Up: [CarpetX](index.md)

## Summary

This page owns the CarpetX CFunction registration layer for the generated BSSN
Ricci, RHS, constraints, and RHS trusted-expression checks. Symbolic BSSN
equation construction stays with
[BSSN Family](../../equations/general-relativity/bssn-family.md), and common
trusted-value mechanics stay with
[Trusted Expression Pipeline](../../equations/trusted-expression-pipeline.md).

`register_CFunction_Ricci_eval()` emits the Ricci auxiliary kernel and schedules
it before the BSSN RHS kernel in `ODESolvers_RHS`. `register_CFunction_rhs_eval()`
emits the evolved RHS kernel with gauge RHSs and optional matter, dissipation,
and improvement terms. `register_CFunction_BSSN_constraints()` emits
Hamiltonian, momentum, and `MSQUARED` diagnostics in `ODESolvers_PostStep`.

## Detail

All three CarpetX GR registration functions use the same parallel-codegen
entry pattern: during `pcg.pcg_registration_phase()` they record the call with
`parallel_codegen.register_func_call()` and return `None`; during codegen they
temporarily replace the global `fd_order` parameter with the requested local
order, register the CFunction, restore the old `fd_order`, and return
`pcg.NRPyEnv()`.

`register_CFunction_Ricci_eval()` selects `BSSN_quantities` with the optional
`_rfm_precompute` key suffix and writes the six `RbarDD*` CarpetX
gridfunctions from `Bq.Ricci_exprs`. Its generated body declares inverse grid
spacings from `CCTK_DELTA_SPACE` and runs finite-difference codegen with
finite-difference helper functions and Golden Kernels inside an interior CarpetX
`simple_loop`. Although the source starts adding SIMD-specific declarations when
requested, the current CarpetX `simple_loop(enable_simd=True)` path raises
`ValueError`, so this registration path does not currently emit a usable SIMD
CarpetX loop kernel. The registered
prefunc rewrites finite-difference helper text from `NO_INLINE` to
`CCTK_ATTRIBUTE_NOINLINE`. Its schedule is guarded by `if(fd_order == <order>)`
and places the function in `ODESolvers_RHS as <thorn>_Ricci before
<thorn>_RHS`, reading `hDD*` and `lambdaU*` and writing `RbarDD*`.

`register_CFunction_rhs_eval()` selects `BSSN_RHSs` with key suffixes for
reference-metric precompute, `RbarDD` gridfunction input, and optional
`T4munu`, then calls `BSSN_gauge_RHSs()` with the requested lapse and shift
evolution options. It appends `alpha_rhs`, `vet_rhsU*`, and `bet_rhsU*` to a
local RHS dictionary, sorts that dictionary, and maps each RHS expression name
to the corresponding CarpetX RHS gridfunction access before codegen.

The RHS options that change generated behavior are `enable_T4munu`,
`enable_KreissOliger_dissipation`, `enable_CAKO`, `enable_CAHD`, `enable_SSL`,
`enable_rfm_precompute`, `enable_simd`, `LapseEvolutionOption`, and
`ShiftEvolutionOption`. `enable_T4munu` changes the BSSN RHS and gauge RHS
lookup keys. `enable_KreissOliger_dissipation` adds `_dKOD` derivative terms to
gauge and non-gauge RHSs; when `enable_CAKO` is true, separate gauge and
non-gauge dissipation strengths are registered and multiplied by the evolved
conformal-factor-derived `W`. `enable_CAHD` builds the matching
`BSSN_constraints` object and adds a Hamiltonian-damping term to `cf_rhs`.
`enable_SSL` registers slow-start lapse parameters and subtracts the
time-dependent lapse term from `alpha_rhs`. `enable_rfm_precompute` changes
reference-metric and BSSN lookup keys. `enable_simd` starts adding the SIMD
header and `CCTK_ParameterGet()` plus `ConstSIMD()` parameter reads, but this
path then passes `enable_simd=True` into CarpetX `simple_loop()`, which raises
`ValueError`; the current simple-loop-backed RHS kernel path is therefore not
usable with `enable_simd=True`.

After option handling, the RHS codegen constructs an upwind control vector
`betaU[i] = vetU[i] * rfm.ReU[i]` and calls `c_codegen()` with finite-difference
codegen, finite-difference helper functions, Golden Kernels, and
`upwind_control_vec=betaU`. The registered prefunc uses the same
`NO_INLINE` to `CCTK_ATTRIBUTE_NOINLINE` rewrite as Ricci. The schedule is
guarded by the finite-difference order and places the function in
`ODESolvers_RHS as <thorn>_RHS after <thorn>_Ricci`, reading
`evol_variables(everywhere)` and `auxevol_variables(interior)` and writing
`evol_variables_rhs(interior)`.

RHS parameter metadata is passed through `ET_current_thorn_CodeParams_used`.
The list always includes `eta` and `fd_order`; it includes either
`diss_strength_gauge` plus `diss_strength_nongauge` when `enable_CAKO` is true,
or `diss_strength` otherwise. It adds `SSL_h` and `SSL_sigma` only when
`enable_SSL` is true, adds `C_CAHD` and
`CFL_FACTOR__ignore_repeats_Carpet_timeref_factors` only when `enable_CAHD` is
true, and adds `PI` for the `Baikal` thorn path.

Two emitted-code caveats are source-observed and not interpreted here as fixes.
The CAHD body templates currently leave an empty operand in the SIMD
`MulSIMD(C_CAHD, MulSIMD(,))` path and in the scalar `C_CAHD * ;` path. The SSL
SIMD body constructs `noSIMD_SSL_Gaussian_prefactor` as a scalar value and then
passes `*noSIMD_SSL_Gaussian_prefactor` to `ConstSIMD()`.

`register_CFunction_BSSN_constraints()` selects `BSSN_constraints` with optional
reference-metric precompute and `T4munu` suffixes, computes `H`, `MU0`, `MU1`,
`MU2`, and `MSQUARED`, and writes them through CarpetX auxiliary gridfunction
accesses. Its generated body uses the same finite-difference helper functions,
Golden Kernels, and `CCTK_ATTRIBUTE_NOINLINE` helper rewrite. As with Ricci and
RHS, SIMD-specific declarations begin when requested, but the downstream
CarpetX `simple_loop(enable_simd=True)` call raises `ValueError`, so this
registration path does not currently emit a usable SIMD CarpetX loop kernel.
Its schedule is guarded by
`fd_order` and places `<thorn>_BSSN_constraints` in `ODESolvers_PostStep`,
reading BSSN state and optional stress-energy fields, writing `aux_variables`,
and syncing `aux_variables`.

RHS trusted-expression validation belongs here because `rhs_eval.py` validates
the CarpetX-specific assembled RHS dictionary after CarpetX option handling and
before final C code emission. The normal in-function path calls
`process_dictionary_of_expressions(..., fixed_mpfs_for_free_symbols=True)` and
then `compare_or_generate_trusted_results()` with a basename containing the
lapse option, shift option, coordinate system, `T4munu`, `KO`, and
`improvements{enable_SSL}`. In this normal path, `improvements` reflects
`enable_SSL`, not a combined CAKO/CAHD/SSL state. With
`validate_expressions=True`, the function returns the processed dictionary
before that in-function comparison. The module `__main__` path drives Cartesian
OnePlusLog/Gamma-driver trusted checks for both `T4munu` states and both
all-improvements states; because it performs the trusted comparison outside
`register_CFunction_rhs_eval()`, those basenames omit the `KO...` segment and
use the local `enable_improvements` loop variable.

CarpetX has six backend-local RHS trusted baselines: four covariant cases span
both `T4munu` states and both improvements states, while two noncovariant
KO-enabled cases span both `T4munu` states with improvements disabled.

## Sources

- [rhs_eval.py](../../../nrpy/infrastructures/CarpetX/general_relativity/rhs_eval.py) - `register_CFunction_rhs_eval`, `validate_expressions`, `__main__`
- [Ricci_eval.py](../../../nrpy/infrastructures/CarpetX/general_relativity/Ricci_eval.py) - `register_CFunction_Ricci_eval`
- [BSSN_constraints.py](../../../nrpy/infrastructures/CarpetX/general_relativity/BSSN_constraints.py) - `register_CFunction_BSSN_constraints`
- [rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuFalse_improvementsFalse.py](../../../nrpy/infrastructures/CarpetX/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuFalse_improvementsFalse.py) - `trusted_dict`
- [rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuFalse_improvementsTrue.py](../../../nrpy/infrastructures/CarpetX/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuFalse_improvementsTrue.py) - `trusted_dict`
- [rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuTrue_improvementsFalse.py](../../../nrpy/infrastructures/CarpetX/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuTrue_improvementsFalse.py) - `trusted_dict`
- [rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuTrue_improvementsTrue.py](../../../nrpy/infrastructures/CarpetX/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuTrue_improvementsTrue.py) - `trusted_dict`
- [rhs_eval_OnePlusLog_GammaDriving2ndOrder_NoCovariant_Cartesian_T4munuFalse_KOTrue_improvementsFalse.py](../../../nrpy/infrastructures/CarpetX/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_NoCovariant_Cartesian_T4munuFalse_KOTrue_improvementsFalse.py) - `trusted_dict`
- [rhs_eval_OnePlusLog_GammaDriving2ndOrder_NoCovariant_Cartesian_T4munuTrue_KOTrue_improvementsFalse.py](../../../nrpy/infrastructures/CarpetX/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_NoCovariant_Cartesian_T4munuTrue_KOTrue_improvementsFalse.py) - `trusted_dict`

## See Also

- [CarpetX](index.md)
- [Code Parameters, Includes, And Loops](code-parameters-includes-and-loops.md)
- [GR ADM/BSSN And Matter Coupling](gr-adm-bssn-and-matter-coupling.md)
- Depends on: [BSSN Family](../../equations/general-relativity/bssn-family.md)
- [Trusted Expression Pipeline](../../equations/trusted-expression-pipeline.md)
- [C Codegen](../../core/c-codegen.md)
- [Finite Difference](../../core/finite-difference.md)
- [Reference Metrics](../../core/reference-metrics.md)
- Contrasts with: [ETLegacy GR BSSN RHS, Ricci, Constraints, And Validation](../etlegacy/gr-bssn-rhs-ricci-constraints-and-validation.md)
