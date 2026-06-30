# C Codegen

> Core route for turning SymPy expressions into generated C text. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Core APIs](index.md)

## Summary

`nrpy.c_codegen` is the main bridge from symbolic NRPy expressions to C code strings. The public `c_codegen()` function normalizes expression and target-name inputs, builds a `CCodeGen` option object, emits C assignments, and delegates to the finite-difference/gridfunction path when requested. That path is deliberately tied to `nrpy.finite_difference`, so generated kernels can read registered gridfunctions, compute derivative stencils, handle upwinded derivatives, and then write expression results.

## Detail

`CCodeGen` stores the option matrix used by `c_codegen()`: formatting (`prestring`, `poststring`, `include_braces`, `enable_clang_format`, and the `clang_format_options` parameter), floating-point type and alias, verbose comments, CSE switches and sorting, optional CSE preprocessing, SIMD rewrite flags, GoldenKernels mode, scalar temporaries, postprocessing substitution suffixes, finite-difference flags, memory layout style, upwind control vector, rational-constant handling, and runtime parameter enforcement. It also reads `finite_difference::fd_order` and rejects non-positive or odd finite-difference orders.

Floating-point aliases are infrastructure-aware. Without SIMD, raw `NRPy` code uses the selected C type such as `double`, BHaH uses `REAL`, and ETLegacy or CarpetX use `CCTK_REAL`, unless the caller explicitly supplies `fp_type_alias`. With `enable_simd`, only double precision is supported; BHaH, CarpetX, and ETLegacy use `REAL_SIMD_ARRAY`, while other infrastructures error because no SIMD floating-point alias is defined.

GoldenKernels is a convenience mode layered on the normal options. `enable_GoldenKernels` forces CSE preprocessing on and enables extra SIMD substitution/FMA discovery by setting `enable_cse_preprocess`, `simd_find_more_subs`, and `simd_find_more_FMAsFMSs` true.

`c_codegen()` accepts one SymPy expression or a list of expressions, paired with one output variable name or a list. It rejects tuple inputs, checks expression/name length agreement, and can emit verbose comments containing the original symbolic expressions. The non-FD path either prints direct assignments when `enable_cse` is false or runs SymPy CSE when it is true. CSE may first call `cse_preprocess()` to factor expressions and convert rationals into symbols; declared rational constants use `rational_const_alias` and the active floating-point alias, while SIMD constants also get `double` backing values and `ConstSIMD()` arrays. When `cse_sorting="none"`, the CSE output is sorted deterministically afterward. Scalar temporaries use `SCALAR_TMP_varnames` and `SCALAR_TMP_sympyexprs`: matching outputs are inserted into the CSE input as `sp.Eq()` objects so postprocessing can preserve their ordering constraints.

`apply_substitution_dict()` appends configured suffixes to matching free-symbol names and returns the substituted SymPy expression. It is a shared postprocessing helper used by C codegen and Python-side codegen routes when symbols must be renamed consistently after expression construction.

`nrpyAbs` is an unevaluated SymPy function used to avoid hangs from SymPy `Abs()` evaluation on complicated expressions while still printing the intended C absolute-value function through the custom C printer mappings.

When `automatically_read_gf_data_from_memory` or `enable_fd_codegen` is enabled, `c_codegen()` extracts derivative symbols from expression free symbols, maps derivative symbols to base gridfunctions and derivative operators, computes finite-difference coefficients and stencils, asks `read_gfs_from_memory()` to generate the needed memory reads, and calls `gridfunction_management_and_FD_codegen()`. `enable_fd_codegen` is stronger than `automatically_read_gf_data_from_memory`: the constructor turns on automatic gridfunction reads, and the public `c_codegen()` entry clears `FDFunctions_dict` at the start of that FD codegen call. `automatically_read_gf_data_from_memory` alone uses the same read/planning route but does not reset the FD helper registry.

`gridfunction_management_and_FD_codegen()` emits up to three generated comment sections. Step 1 reads gridfunctions from memory and computes finite-difference stencils; Step 2 implements the upwind algorithm when an upwind control vector and `dupD` operators require it; Step 3 evaluates the remaining SymPy expressions and writes results. The upwind path extracts the used directions, emits `UpwindControlVectorU*`, computes `UPWIND_ALG(UpwindControlVectorU*)`, and combines `UpwindAlgInput*_ddnD*` and `UpwindAlgInput*_dupD*` into the final `dupD` derivative. SIMD Step 3 assigns each RHS to a temporary `REAL_SIMD_ARRAY __RHS_exp_*` and emits `WriteSIMD(&output, __RHS_exp_*)`.

Finite-difference codegen recursively calls `c_codegen()` for each generated part with local options: `FDPart1` for stencil arithmetic, `FDPart2` for upwind arithmetic, and `FDPart3` for final expression evaluation. Those recursive calls explicitly disable automatic memory reads and FD recursion for the stencil part and pass rational-symbol dictionaries from finite-difference preprocessing. With `enable_fd_functions`, Step 1 emits helper calls instead of inline stencil formulas, prepares upwind-control expressions separately, and builds each `FDFunction.CFunction` so `construct_FD_functions_prefunc()` can later collect the static helper functions for a C function `prefunc`.

Helper pages own the lower-level transformation mechanics used by this interface. CSE preprocessing, postprocessing, deterministic temporary sorting, custom C printer mappings, and JAX printer behavior live under [CSE And Printer Support](helpers/cse-and-printer-support.md); symbolic SIMD rewrites and handwritten intrinsic headers live under [SIMD And Intrinsic Support](helpers/simd-and-intrinsic-support.md). Use this page for the codegen interface itself. Use the finite-difference page when the question is about derivative naming, stencil construction, prototype operators, or memory-read ordering.

## Sources

- [nrpy/c_codegen.py](../../nrpy/c_codegen.py) - `CCodeGen`, `c_codegen`, `apply_substitution_dict`, `nrpyAbs`, `gridfunction_management_and_FD_codegen`
- [nrpy/py_codegen.py](../../nrpy/py_codegen.py) - `py_codegen`
- [nrpy/finite_difference.py](../../nrpy/finite_difference.py) - `compute_fdcoeffs_fdstencl`, `read_gfs_from_memory`, `FDFunction`

## See Also

- Parent: [Core APIs](index.md)
- Depends on: [Finite Difference](finite-difference.md)
- Depends on: [CSE And Printer Support](helpers/cse-and-printer-support.md)
- Depends on: [SIMD And Intrinsic Support](helpers/simd-and-intrinsic-support.md)
- Contrasts with: [Python Codegen](python-codegen.md)
- See also: [C Function Registry](c-function-registry.md)
