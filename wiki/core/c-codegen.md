# C Codegen

> Core route for turning SymPy expressions into generated C text. · Status: contested
> Up: [Core APIs](index.md)

## Summary

`nrpy.c_codegen` is the main bridge from symbolic NRPy expressions to C code strings. The public `c_codegen()` function normalizes expression and target-name inputs, builds a `CCodeGen` option object, emits C assignments, and delegates to the finite-difference/gridfunction path when requested. That path is deliberately tied to `nrpy.finite_difference`, so generated kernels can read registered gridfunctions, compute derivative stencils, handle upwinded derivatives, and then write expression results. The scalar CSE/substitution combination has a known failure described below.

## Detail

`CCodeGen` stores the option matrix used by `c_codegen()`: formatting (`prestring`, `poststring`, `include_braces`, `enable_clang_format`, and the `clang_format_options` parameter), floating-point type and alias, verbose comments, CSE switches and sorting, optional CSE preprocessing, SIMD rewrite flags, GoldenKernels mode, scalar temporaries, postprocessing substitution suffixes, finite-difference flags, memory layout style, upwind control vector, rational-constant handling, and the accepted `enforce_c_parameters_must_be_defined` flag. It also reads `finite_difference::fd_order` and rejects non-positive or odd finite-difference orders.

Some accepted options have narrower current behavior than their names suggest. `enable_clang_format` is consulted only at the end of `gridfunction_management_and_FD_codegen()`. Public `c_codegen()` does not forward this option to that function, and ordinary non-gridfunction emission does not apply it either. Only a direct call to the finite-difference helper can activate its final formatting branch. On the `enable_fd_codegen` and `automatically_read_gf_data_from_memory` paths `c_codegen()` always passes the assembled finite-difference body through `order_statements_for_register_pressure()` (see [CSE And Printer Support](helpers/cse-and-printer-support.md)) before returning it, so whenever the body is scheduled the generated kernel no longer carries the read/derivative/expression phase comments; the helper returns the body unchanged, comments included, when it contains anything other than single-assignment `const` definitions and stores, or a load that follows a store. `enforce_c_parameters_must_be_defined` is stored but is not read elsewhere in `nrpy/c_codegen.py`; it currently enforces nothing. `SCALAR_TMP_varnames` and `SCALAR_TMP_sympyexprs` are caller-coordinated positional lists, and the implementation does not validate their paired lengths before indexing. Treat these as current implementation limits, not guarantees implied by option names.

Claim evidence:
- Claim: Public `c_codegen()` does not apply `enable_clang_format`, including its finite-difference route; the formatting branch exists only in `gridfunction_management_and_FD_codegen()` when called directly with that option.
- Role: descriptive behavior
- Deciding authority: [c_codegen.py](../../nrpy/c_codegen.py), `c_codegen` and `gridfunction_management_and_FD_codegen`
- Corroboration: none available; the public call's keyword list and the helper's final branch define this behavior.

Floating-point aliases are infrastructure-aware. Without SIMD, raw `NRPy` code uses the selected C type such as `double`, BHaH uses `REAL`, and ETLegacy or CarpetX use `CCTK_REAL`, unless the caller explicitly supplies `fp_type_alias`. With `enable_simd`, only double precision is supported; BHaH, CarpetX, ETLegacy, and Dendro use `REAL_SIMD_ARRAY`. Dendro scalar callers must explicitly pass `fp_type_alias=gri.DENDRO_SCALAR_TYPE` to obtain `DendroScalar`; the constructor otherwise leaves the scalar alias empty for that infrastructure.

Claim evidence:
- Claim: Dendro scalar code generation requires an explicit `fp_type_alias` to obtain `DendroScalar`; `CCodeGen` does not choose that alias automatically.
- Role: descriptive behavior
- Deciding authority: [c_codegen.py](../../nrpy/c_codegen.py), `CCodeGen.__init__`
- Corroboration: [Dendro Ricci_eval.py](../../nrpy/infrastructures/Dendro/general_relativity/Ricci_eval.py), `register_CFunction_Ricci_eval`, explicitly supplies the scalar type.

GoldenKernels is a convenience mode layered on the normal options. `enable_GoldenKernels` forces CSE preprocessing on and enables extra SIMD substitution/FMA discovery by setting `enable_cse_preprocess`, `simd_find_more_subs`, and `simd_find_more_FMAsFMSs` true.

`c_codegen()` accepts one SymPy expression or a list of expressions, paired with one output variable name or a list. It rejects tuple inputs, checks expression/name length agreement, and can emit verbose comments containing the original symbolic expressions. The non-FD path either prints direct assignments when `enable_cse` is false or runs SymPy CSE when it is true. CSE may first call `cse_preprocess()` to factor expressions and convert rationals into symbols; declared rational constants use `rational_const_alias` and the active floating-point alias, while SIMD constants also get `double` backing values and `ConstSIMD()` arrays. When `cse_sorting="none"`, the CSE output passes through `sort_cse_output_deterministically()` afterward; tied structural keys can retain input order, and generated names must avoid input-symbol collisions. Scalar temporaries use `SCALAR_TMP_varnames` and `SCALAR_TMP_sympyexprs`: matching outputs are inserted into the CSE input as `sp.Eq()` objects so postprocessing can preserve their ordering constraints.

Claim status: contested; contradiction: CONTR-0014.
See [CONTR-0014](../contradictions.md#contr-0014) and
[CSE And Printer Support](helpers/cse-and-printer-support.md) for ordering and naming limits.

Claim evidence:
- Claim: The post-sort for `cse_sorting="none"` does not guarantee permutation-invariant output for tied structural keys or collision-free generated names.
- Role: descriptive behavior
- Deciding authority: [cse_preprocess_postprocess.py](../../nrpy/helpers/cse_preprocess_postprocess.py), `cse_preprocess` and `sort_cse_output_deterministically`
- Corroboration: [c_codegen.py](../../nrpy/c_codegen.py), `c_codegen`, selects these helpers during CSE emission.

`apply_substitution_dict()` appends configured suffixes to matching free-symbol names and returns the substituted SymPy expression. It is a shared postprocessing helper used by C codegen and Python-side codegen routes when symbols must be renamed consistently after expression construction.

Ordinary non-FD, non-SIMD CSE emission has a substitution limitation shared with
Python codegen. When a postprocessed CSE temporary survives and
`postproc_substitution_dict` is nonempty, the scalar C loop assigns to a tuple
entry and raises `TypeError`, even if no substitution key matches. That option
combination is not guaranteed to return C code.

Claim status: contested; contradiction: CONTR-0012.
See [CONTR-0012](../contradictions.md#contr-0012) for the deciding code and
resolution check.

`nrpyAbs` is an unevaluated SymPy function used to avoid hangs from SymPy `Abs()` evaluation on complicated expressions while still printing the intended C absolute-value function through the custom C printer mappings.

When `automatically_read_gf_data_from_memory` or `enable_fd_codegen` is enabled, `c_codegen()` extracts derivative symbols from expression free symbols, maps derivative symbols to base gridfunctions and derivative operators, computes finite-difference coefficients and stencils, and calls `gridfunction_management_and_FD_codegen()`. `enable_fd_codegen` is stronger than `automatically_read_gf_data_from_memory`: the constructor turns on automatic gridfunction reads, and the public `c_codegen()` entry clears `FDFunctions_dict` at the start of that FD codegen call. `automatically_read_gf_data_from_memory` alone uses the same read/planning route but does not reset the FD helper registry.

`stored_first_derivatives` is an optional sequence of registered gridfunction names that the caller makes available to this kernel. Before constructing FD prototypes and planning reads, the FD path calls `select_stored_first_derivatives()` to select point reads for centered first derivatives and first-derivative stencils for canonical mixed second derivatives. The original derivative temporary names and equation expressions remain unchanged. Merely registering storage does not activate it: the consumer must explicitly select the fields, supply their array, and guarantee current values and stencil halos. An omitted or empty sequence retains ordinary lowering. Diagonal, upwind, KO, and unselected derivatives keep their operators. See [Finite Difference](finite-difference.md) for the naming contract.

Only a list-valued `upwind_control_vec` activates direction selection and control-vector expressions. Current constructor logic does not comprehensively validate list length or every element, so callers must supply an entry for each derivative direction used. Non-list values do not establish an upwind control vector even if accepted by construction. With the default `sp.Symbol("unset")` and a `dupD` request, generation warns and the zero control selector chooses the downwind input; it does not simply return the positive-shift stencil.

Claim evidence:
- Claim: With the default symbolic `upwind_control_vec` and a requested `dupD`, generation warns and uses a zero selector, which selects the downwind input. A list-valued control vector is needed for direction-dependent selection.
- Role: descriptive behavior
- Deciding authority: [c_codegen.py](../../nrpy/c_codegen.py), `gridfunction_management_and_FD_codegen`, warning, `upwindU`, and selection expression
- Corroboration: [finite_difference.py](../../nrpy/finite_difference.py), `extract_list_of_deriv_var_strings_from_sympyexpr_list`, supplies the implicit downwind derivative for a non-string control value.

`gridfunction_management_and_FD_codegen()` plans memory reads after eliminating reused implicit downwind derivatives. Pointer/stride mixed helpers and CUDA KO helpers own their stencil loads; direct field references and other derivatives still contribute caller reads. It emits up to three phase comment sections, which `c_codegen()` removes when it re-emits the scheduled body. Step 1 reads gridfunctions from memory and computes finite-difference stencils; Step 2 implements the upwind algorithm when an upwind control vector and `dupD` operators require it; Step 3 evaluates the remaining SymPy expressions and writes results. The upwind path extracts the used directions, emits `UpwindControlVectorU*`, and computes `Upwind* = UPWIND_ALG(UpwindControlVectorU*)`. Normally, selection combines `UpwindAlgInput*_ddnD*` and `UpwindAlgInput*_dupD*` into the final `dupD` derivative. For BHaH with `fp_type="double"` and a list-valued control vector, a requested `dKOD` derivative for the same field and direction as `dupD0`, `dupD1`, or `dupD2` replaces the implicit downwind computation, provided `ko_fd_order` is omitted or equals `fd_order` and the original expressions do not explicitly request that `ddnD` derivative. Selection then uses `dup + (Upwind - 1) * A(n) * raw_KO`, where `dup` is the upwind input, `n = fd_order`, and `A(n) = 2**(n+3) / (n * binomial(n, n/2))`. Raw KO is the derivative before any dissipation strength or equation prefactor and remains available to the RHS. This applies to scalar/SIMD and inline/helper emission; explicit downwind requests, unmatched field/direction pairs, other infrastructures, and other precisions retain the existing construction. The stencil identity is exact, but reassociated floating-point evaluation need not be bitwise identical.

Claim evidence:
- Claim: Only BHaH double-precision generation with list-valued upwind control and a same-field/same-direction KO request at the regular derivative order replaces an implicit `ddnD` computation with raw KO; explicitly requested downwind derivatives are preserved.
- Role: descriptive behavior
- Deciding authority: [nrpy/c_codegen.py](../../nrpy/c_codegen.py) - `gridfunction_management_and_FD_codegen`, `upwind_ko_symbols` eligibility/filtering and final upwind selection; `c_codegen` scalar/SIMD and inline/helper doctest examples.
- Corroboration: [nrpy/finite_difference.py](../../nrpy/finite_difference.py) - `compute_fdcoeffs_fdstencl` defines the shifted upwind/downwind and raw KO stencils used by the identity.

SIMD Step 3 assigns each RHS to a temporary `REAL_SIMD_ARRAY __RHS_exp_*` and emits `WriteSIMD(&output, __RHS_exp_*)`.

Finite-difference codegen recursively calls `c_codegen()` for each generated part with local options: `FDPart1` for stencil arithmetic, `FDPart2` for upwind arithmetic, and `FDPart3` for final expression evaluation. Those recursive calls explicitly disable automatic memory reads and FD recursion for the stencil part and pass rational-symbol dictionaries from finite-difference preprocessing. With `enable_fd_functions`, Step 1 emits helper calls instead of inline stencil formulas, prepares upwind-control expressions separately, and builds each `FDFunction.CFunction` so `construct_FD_functions_prefunc()` can later collect the static helper functions for a C function `prefunc`.

An explicit `ko_fd_order` must be a positive, even Python integer, with Boolean
values rejected. `CCodeGen` rejects any explicit `ko_fd_order` when
`enable_fd_functions=True`, including a value equal to `fd_order`; leave it
unset for helper generation.

Claim evidence:
- Claim: `CCodeGen` validates explicit KO order as a positive, even, non-Boolean Python integer and rejects combining any explicit KO order with FD helper generation.
- Role: descriptive behavior
- Deciding authority: [c_codegen.py](../../nrpy/c_codegen.py), `CCodeGen.__init__`
- Corroboration: none available; the constructor's explicit validation branches define these option restrictions.

Helper pages own the lower-level transformation mechanics used by this interface. CSE preprocessing, postprocessing, temporary sorting, custom C printer mappings, and JAX printer behavior live under [CSE And Printer Support](helpers/cse-and-printer-support.md); symbolic SIMD rewrites and handwritten intrinsic headers live under [SIMD And Intrinsic Support](helpers/simd-and-intrinsic-support.md). Use this page for the codegen interface itself. Use the finite-difference page when the question is about derivative naming, stencil construction, prototype operators, or memory-read ordering.

Import `CCodeGen` and `c_codegen` from `nrpy.c_codegen`. The empty package initializer `nrpy/__init__.py` does not re-export them.

## Sources

- [Dendro Ricci_eval.py](../../nrpy/infrastructures/Dendro/general_relativity/Ricci_eval.py) - `register_CFunction_Ricci_eval`

- [nrpy/c_codegen.py](../../nrpy/c_codegen.py) - `CCodeGen`, `c_codegen`, `apply_substitution_dict`, `nrpyAbs`, `gridfunction_management_and_FD_codegen`
- [nrpy/py_codegen.py](../../nrpy/py_codegen.py) - `py_codegen`
- [nrpy/finite_difference.py](../../nrpy/finite_difference.py) - `compute_fdcoeffs_fdstencl`, `read_gfs_from_memory`, `FDFunction`
- [nrpy/__init__.py](../../nrpy/__init__.py) - empty package initializer; no core-codegen re-exports

## See Also

- Parent: [Core APIs](index.md)
- Depends on: [Finite Difference](finite-difference.md)
- Depends on: [CSE And Printer Support](helpers/cse-and-printer-support.md)
- Depends on: [SIMD And Intrinsic Support](helpers/simd-and-intrinsic-support.md)
- Contrasts with: [Python Codegen](python-codegen.md)
- See also: [C Function Registry](c-function-registry.md)
