# SIMD And Intrinsic Support

> Helper route for symbolic SIMD rewrites and handwritten intrinsic macro headers. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [Helper APIs](index.md)

## Summary

`nrpy.helpers.simd` rewrites SymPy expressions into symbolic function calls with names such as `AddSIMD`, `MulSIMD`, `DivSIMD`, and `FusedMulAddSIMD`. Those symbolic calls are still generator-side expressions; the handwritten headers `simd_intrinsics.h` and `cuda_intrinsics.h` provide the C/CUDA macro contracts that generated code can include when those names need concrete scalar, SIMD, or CUDA behavior.

## Detail

`expr_convert_to_simd_intrins()` is the transformation entry point. If the caller does not supply a `symbol_to_Rational_dict`, it first calls `cse_preprocess()` so rational constants can be represented by generated symbols before the SIMD pass. The function then builds an `ExprTree` over the expression, creates SymPy function objects for the SIMD spelling, and rewrites the tree in stages: absolute value and transcendental functions, powers and selected roots, subtraction patterns, binary addition and multiplication trees, division patterns, fused multiply-add/subtract patterns, negative fused patterns, and optional cleanup of an unused prefixed `_NegativeOne_` constant.

The SIMD pass is algebraic and symbolic. It does not choose an instruction set, allocate vector registers, or include a header. It produces expression text containing names such as `SqrtSIMD`, `CosSIMD`, `PowSIMD`, `FusedMulSubSIMD`, `NegFusedMulAddSIMD`, and `NegFusedMulSubSIMD`; later C code generation prints those names into generated code. `nrpy.c_codegen.CCodeGen` is the integration point that enables this path through `enable_simd`, restricts SIMD output to double precision, uses `REAL_SIMD_ARRAY` for supported infrastructures, calls `expr_convert_to_simd_intrins()` on common subexpressions and final expressions, and emits SIMD constant declarations through `ConstSIMD`.

Debug mode is an internal validation path for the symbolic rewrite. `_debug_eval_namespace()` builds an explicit namespace containing SymPy objects, free symbols, and `*_SIMD_check` functions such as `AbsSIMD_check`, `ConstSIMD_check`, `FusedMulAddSIMD_check`, and `NegFusedMulSubSIMD_check`. With `debug=True`, `expr_convert_to_simd_intrins()` evaluates the transformed spelling through that restricted namespace, back-substitutes rational symbols, compares against the original expression, and raises `ValueError` if simplification cannot prove the difference is zero. This documents the helper's own validation behavior, not a general recommendation to evaluate arbitrary expression strings.

`simd_intrinsics.h` is the C-side source contract for SIMD names. It defines `REAL_SIMD_ARRAY`, `SIMD_WIDTH`, load/store macros, constant construction, arithmetic, selected math functions, fused multiply-add/subtract spellings, absolute value, zero initialization, horizontal addition, and `UPWIND_ALG`. The header selects among AVX512F, AVX, SSE2/SSE3, and scalar fallback branches using compiler feature macros such as `__AVX512F__`, `__AVX__`, `__SSE2__`, `__SSE3__`, and `__FMA__`. The scalar branch keeps the same macro names with width one, so generated code has a fallback contract even when the vector feature macros are absent.

`cuda_intrinsics.h` is the CUDA-oriented scalar/intrinsic contract. When `__CUDACC__` is defined, macros such as `AddCUDA`, `SubCUDA`, `MulCUDA`, `DivCUDA`, `FusedMulAddCUDA`, `SqrtCUDA`, and `ReadCUDA` use CUDA device intrinsics or device read helpers where the header defines them. The non-CUDA branch keeps the same `REAL_CUDA_ARRAY`, `CUDA_WIDTH`, arithmetic, load/store, horizontal-add, and `UPWIND_ALG` macro names as scalar C expressions. This header is separate from the SIMD symbolic transformer: it defines CUDA macro spellings and fallback behavior, but `expr_convert_to_simd_intrins()` itself emits `*SIMD` names, not `*CUDA` names.

Generated projects may contain copied versions of these headers, but those copies are generated output. For KB purposes, cite the handwritten sources under `nrpy/helpers/` and the codegen integration in `nrpy/c_codegen.py`, not copies under `project/`.

## Sources

- [nrpy/helpers/simd.py](../../../nrpy/helpers/simd.py) - `expr_convert_to_simd_intrins`, `_debug_eval_namespace`, `AbsSIMD_check`, `ConstSIMD_check`, `FusedMulAddSIMD_check`, `FusedMulSubSIMD_check`, `NegFusedMulAddSIMD_check`, `NegFusedMulSubSIMD_check`
- [nrpy/helpers/simd_intrinsics.h](../../../nrpy/helpers/simd_intrinsics.h) - `REAL_SIMD_ARRAY`, `SIMD_WIDTH`, `UPWIND_ALG`, `FusedMulAddSIMD`, `HorizAddSIMD`
- [nrpy/helpers/cuda_intrinsics.h](../../../nrpy/helpers/cuda_intrinsics.h) - `REAL_CUDA_ARRAY`, `CUDA_WIDTH`, `UPWIND_ALG`, `FusedMulAddCUDA`, `ReadCUDA`
- [nrpy/c_codegen.py](../../../nrpy/c_codegen.py) - `CCodeGen`, `c_codegen`, `gridfunction_management_and_FD_codegen`
- [Intel Intrinsics Guide](https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html) - SSE, AVX, AVX512, and intrinsic-family terminology; accessed 07-12-2026
- [NVIDIA CUDA Programming Guide](https://docs.nvidia.com/cuda/cuda-programming-guide/index.html) - CUDA programming-model terminology; version 13.3 page accessed 07-12-2026

## See Also

- [Helper APIs](index.md)
- [CSE And Printer Support](cse-and-printer-support.md)
- [Loop Kernel And Device Helpers](loop-kernel-and-device-helpers.md)
- [C Codegen](../c-codegen.md)
- [Finite Difference](../finite-difference.md)
- [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
