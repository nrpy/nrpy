# C Codegen

> Core route for turning SymPy expressions into generated C text. · Status: contested
> Up: [Core APIs](index.md)

## Summary

`nrpy.c_codegen` lowers symbolic NRPy expressions into C code strings. The public `c_codegen()` function normalizes expression and target-name inputs, builds a `CCodeGen` option object, emits C assignments, and delegates to the finite-difference/gridfunction path when requested. That path is deliberately tied to `nrpy.finite_difference`, so generated kernels can read registered gridfunctions, compute derivative stencils, handle upwinded derivatives, and then write expression results. The scalar CSE/substitution combination has a known failure described below.

## Detail

### `CCodeGen` Options And Scalar CSE Emission

`CCodeGen` stores the option matrix used by `c_codegen()`: formatting (`prestring`, `poststring`, `include_braces`, `enable_clang_format`, and the `clang_format_options` parameter), floating-point type and alias, verbose comments, CSE switches and sorting, optional CSE preprocessing, SIMD rewrite flags, GoldenKernels mode, scalar temporaries, postprocessing substitution suffixes, finite-difference flags, memory layout style, upwind control vector, rational-constant handling, and the accepted `enforce_c_parameters_must_be_defined` flag. It also reads `finite_difference::fd_order` and rejects non-positive or odd finite-difference orders.

Some accepted options have narrower current behavior than their names suggest. `enable_clang_format` is consulted only at the end of `gridfunction_management_and_FD_codegen()`. Public `c_codegen()` does not forward this option to that function, and ordinary non-gridfunction emission does not apply it either; only a direct call to the finite-difference helper can activate its final formatting branch.
On the `enable_fd_codegen` and `automatically_read_gf_data_from_memory` paths `c_codegen()` always passes the assembled finite-difference body through `order_statements_for_register_pressure()` (see [CSE And Printer Support](helpers/cse-and-printer-support.md)) before returning it, so whenever the body is scheduled the generated kernel no longer carries the read/derivative/expression phase comments; the helper returns the body unchanged, comments included, when it contains anything other than single-assignment `const` definitions and stores, or a `const` definition that follows a store.
`enforce_c_parameters_must_be_defined` is stored but is not read elsewhere in `nrpy/c_codegen.py`; it currently enforces nothing. `SCALAR_TMP_varnames` and `SCALAR_TMP_sympyexprs` are caller-coordinated positional lists, and the implementation does not validate their paired lengths before indexing. Treat these as current implementation limits, not guarantees implied by option names.

Evidence (descriptive behavior): `c_codegen.py::c_codegen`; `c_codegen.py::gridfunction_management_and_FD_codegen`.

Floating-point aliases are infrastructure-aware. Without SIMD, raw `NRPy` code uses the selected C type such as `double`, BHaH uses `REAL`, and ETLegacy or CarpetX use `CCTK_REAL`, unless the caller explicitly supplies `fp_type_alias`. With `enable_simd`, only double precision is supported; BHaH, CarpetX, ETLegacy, and Dendro use `REAL_SIMD_ARRAY`. Dendro scalar callers must explicitly pass `fp_type_alias=gri.DENDRO_SCALAR_TYPE` to obtain `DendroScalar`; the constructor otherwise leaves the scalar alias empty for that infrastructure.

Evidence (descriptive behavior): `c_codegen.py::CCodeGen.__init__`; corroboration: `Dendro/general_relativity/Ricci_eval.py::register_CFunction_Ricci_eval`, explicitly supplies the scalar type.

GoldenKernels is a convenience mode layered on the normal options. `enable_GoldenKernels` sets `enable_cse_preprocess` (default false), `simd_find_more_subs` (default false), and `simd_find_more_FMAsFMSs` (default true) to true, even when the caller passed false.

`c_codegen()` accepts one SymPy expression or a list of expressions, paired with one output variable name or a list. It rejects tuple inputs, checks expression/name length agreement, and can emit verbose comments containing the original symbolic expressions. The non-FD path either prints direct assignments when `enable_cse` is false or runs SymPy CSE when it is true.
CSE may first call `cse_preprocess()` to factor expressions and convert rationals into symbols; declared rational constants use `rational_const_alias` and the active floating-point alias, while SIMD constants also get `double` backing values and `ConstSIMD()` arrays. Temporaries from `sp.cse` are numbered `<cse_varprefix>tmpN` and skip every identifier that occurs in the output strings. When `cse_sorting="none"`, the CSE output passes through `sort_cse_output_deterministically()` afterward, which renames temporaries without taking a name from the output strings or the input symbols; tied structural keys can retain input order. Scalar temporaries use `SCALAR_TMP_varnames` and `SCALAR_TMP_sympyexprs`: matching outputs are inserted into the CSE input as `sp.Eq()` objects so postprocessing can preserve their ordering constraints.

Claim status: contested; contradiction: CONTR-0014.
See [CONTR-0014](../contradictions.md#contr-0014) and
[CSE And Printer Support](helpers/cse-and-printer-support.md) for the ordering limit.

Evidence (descriptive behavior): `cse_preprocess_postprocess.py::sort_cse_output_deterministically`; corroboration: `c_codegen.py::c_codegen`.

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

### Finite-Difference, Stored-Derivative, And SIMD Lowering

When `automatically_read_gf_data_from_memory` or `enable_fd_codegen` is enabled, `c_codegen()` extracts derivative symbols from expression free symbols, maps derivative symbols to base gridfunctions and derivative operators, computes finite-difference coefficients and stencils, and calls `gridfunction_management_and_FD_codegen()`. `enable_fd_codegen` is stronger than `automatically_read_gf_data_from_memory`: the constructor turns on automatic gridfunction reads, and the public `c_codegen()` entry clears `FDFunctions_dict` at the start of that FD codegen call. `automatically_read_gf_data_from_memory` alone uses the same read/planning route but does not reset the FD helper registry.

`stored_first_derivatives` is an optional sequence of registered gridfunction names that the caller makes available to this kernel. Before constructing FD prototypes and planning reads, the FD path calls `select_stored_first_derivatives()` to select point reads for centered first derivatives and first-derivative stencils for canonical mixed second derivatives. The original derivative temporary names and equation expressions remain unchanged. Merely registering storage does not activate it: the calling module must explicitly select the fields, supply their array, and guarantee current values and stencil halos. An omitted or empty sequence retains ordinary lowering. Diagonal, upwind, KO, and unselected derivatives keep their operators. See [Finite Difference](finite-difference.md) for the derivative naming rules.

Only a list-valued `upwind_control_vec` activates direction selection and control-vector expressions. Current constructor logic does not comprehensively validate list length or every element, so callers must supply an entry for each derivative direction used. Non-list values do not establish an upwind control vector even if accepted by construction. With the default `sp.Symbol("unset")` and a `dupD` request, generation warns and the zero control selector chooses the downwind input; it does not simply return the positive-shift stencil.

Evidence (descriptive behavior): `c_codegen.py::gridfunction_management_and_FD_codegen` (`upwindU` warning and selection); corroboration: `finite_difference.py::extract_list_of_deriv_var_strings_from_sympyexpr_list` (implicit downwind derivative).

`gridfunction_management_and_FD_codegen()` plans memory reads after eliminating reused implicit downwind derivatives. Pointer/stride mixed helpers and CUDA KO helpers load their own stencils; direct field references and other derivatives still contribute caller reads. It emits up to three phase comment sections, which `c_codegen()` removes when it re-emits the scheduled body. Each section is labeled `Step N of M`, where `N` counts the sections in order and `M` is the number emitted: the read section reads gridfunctions from memory and computes finite-difference stencils when reads or helper calls are present, the upwind section implements the upwind algorithm when an upwind control vector and `dupD` operators require it, and the evaluation section evaluates the remaining SymPy expressions and writes results, always appearing last (and numbered Step 3 only when the other two sections appear).
The upwind path extracts the used directions, emits `UpwindControlVectorU*`, and computes `Upwind* = UPWIND_ALG(UpwindControlVectorU*)`. Normally, selection combines `UpwindAlgInput*_ddnD*` and `UpwindAlgInput*_dupD*` into the final `dupD` derivative. For BHaH with `fp_type="double"` and a list-valued control vector, a requested `dKOD` derivative for the same field and direction as `dupD0`, `dupD1`, or `dupD2` replaces the implicit downwind computation, provided `ko_fd_order` is omitted or equals `fd_order` and the original expressions do not explicitly request that `ddnD` derivative.
Selection then uses `dup + (Upwind - 1) * A(n) * raw_KO`, where `dup` is the upwind input, `n = fd_order`, and `A(n) = 2**(n+3) / (n * binomial(n, n/2))`; raw KO is the derivative before any dissipation strength or equation prefactor and remains available to the RHS. This applies to scalar/SIMD and inline/helper emission, whereas explicit downwind requests, unmatched field/direction pairs, other infrastructures, and other precisions retain the existing construction. The stencil identity is exact, though reassociated floating-point evaluation need not be bitwise identical.

Evidence (descriptive behavior): `c_codegen.py::gridfunction_management_and_FD_codegen`; `c_codegen.py::c_codegen`.

With SIMD enabled, the evaluation section of `gridfunction_management_and_FD_codegen()` (Step 3 when all three phase comment sections appear) assigns each RHS to a temporary `REAL_SIMD_ARRAY __RHS_exp_*` and emits `WriteSIMD(&output, __RHS_exp_*)`.

### Recursive FD Emission And `ko_fd_order`

Finite-difference codegen recursively calls `c_codegen()` for each generated part with local options: `FDPart1` for stencil arithmetic, `FDPart2` for upwind arithmetic, and `FDPart3` for final expression evaluation. Those recursive calls explicitly disable automatic memory reads and FD recursion for the stencil part and pass rational-symbol dictionaries from finite-difference preprocessing. With `enable_fd_functions` and at least one derivative operator, Step 1 of `gridfunction_management_and_FD_codegen()` emits helper calls instead of inline stencil formulas, prepares upwind-control expressions separately, and builds each `FDFunction.CFunction` so `construct_FD_functions_prefunc()` can later collect the static helper functions for a C function `prefunc`.

An explicit `ko_fd_order` must be a positive, even Python integer, with Boolean
values rejected. `CCodeGen` rejects any explicit `ko_fd_order` when
`enable_fd_functions=True`, including a value equal to `fd_order`; leave it
unset for helper generation.

Evidence (descriptive behavior): `c_codegen.py::CCodeGen.__init__`.

Import `CCodeGen` and `c_codegen` from `nrpy.c_codegen`. The empty package initializer `nrpy/__init__.py` does not re-export them. For CSE preprocessing/sorting and printer mappings see [CSE And Printer Support](helpers/cse-and-printer-support.md); for SIMD intrinsics see [SIMD And Intrinsic Support](helpers/simd-and-intrinsic-support.md).

## Sources

- [Dendro Ricci_eval.py](../../nrpy/infrastructures/Dendro/general_relativity/Ricci_eval.py) - `register_CFunction_Ricci_eval`
- [cse_preprocess_postprocess.py](../../nrpy/helpers/cse_preprocess_postprocess.py) - `cse_preprocess`, `sort_cse_output_deterministically`
- [nrpy/c_codegen.py](../../nrpy/c_codegen.py) - `CCodeGen`, `c_codegen`, `apply_substitution_dict`, `nrpyAbs`, `gridfunction_management_and_FD_codegen`
- [nrpy/finite_difference.py](../../nrpy/finite_difference.py) - `compute_fdcoeffs_fdstencl`, `read_gfs_from_memory`, `FDFunction`
- [nrpy/__init__.py](../../nrpy/__init__.py) - empty package initializer; no core-codegen re-exports

## See Also

- Depends on: [Finite Difference](finite-difference.md)
- Depends on: [CSE And Printer Support](helpers/cse-and-printer-support.md)
- Depends on: [SIMD And Intrinsic Support](helpers/simd-and-intrinsic-support.md)
