# SIMD And Intrinsic Support

> Helper route for symbolic SIMD rewrites and handwritten intrinsic macro headers. · Status: contested
> Up: [Helper APIs](index.md)

## Summary

`nrpy.helpers.simd` rewrites SymPy expressions into symbolic calls such as `AddSIMD`, `MulSIMD`, `DivSIMD`, and `FusedMulAddSIMD`. Handwritten headers supply C/CUDA implementations for supported names. The transformation and header must be checked together: some emitted names and constants are missing definitions, and some header branches differ from the symbolic arithmetic.

## Detail

`expr_convert_to_simd_intrins()` is the transformation entry point. If the caller does not supply a `symbol_to_Rational_dict`, it first calls `cse_preprocess()` so rational constants can be represented by generated symbols before the SIMD pass. The function then builds an `ExprTree` over the expression, creates SymPy function objects for the SIMD spelling, and rewrites the tree in stages: absolute value and transcendental functions, powers and selected roots, subtraction patterns, binary addition and multiplication trees, division patterns, fused multiply-add/subtract patterns, negative fused patterns, and optional cleanup of an unused prefixed `_NegativeOne_` constant.

The SIMD pass is algebraic and symbolic. It does not choose an instruction set, allocate vector registers, or include a header. It produces expression text containing names such as `SqrtSIMD`, `CosSIMD`, `PowSIMD`, `FusedMulSubSIMD`, `NegFusedMulAddSIMD`, and `NegFusedMulSubSIMD`; later C code generation prints those names into generated code. `nrpy.c_codegen.CCodeGen` is the integration point that enables this path through `enable_simd`, restricts SIMD output to double precision, uses `REAL_SIMD_ARRAY` for supported infrastructures, calls `expr_convert_to_simd_intrins()` on common subexpressions and final expressions, and emits SIMD constant declarations through `ConstSIMD`.

Reciprocal half-power rewrites for exponents `-1/2`, `-3/2`, and `-5/2` introduce the prefixed `_Integer_1` numerator without registering it in `symbol_to_Rational_dict`. If it was not already present, declarations derived from that dictionary omit it; the expression requires a matching numerator constant. Separately, `PowSIMD`, `CbrtSIMD`, `LogSIMD`, and `SignSIMD` can be emitted but have no definitions in `simd_intrinsics.h`. Including that header alone is insufficient for those expressions; supply suitable definitions or use a supported expression form.

Claim status: contested; contradiction: CONTR-0015. See [CONTR-0015](../../contradictions.md#contr-0015).

Claim evidence:
- Claim: Reciprocal half-powers can omit numerator registration, and emitted Pow/Cbrt/Log/Sign SIMD names lack definitions in the SIMD header.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/simd.py](../../../nrpy/helpers/simd.py), `expr_convert_to_simd_intrins`
- Corroboration: `none available`; the named implementation decides this behavior.
- Additional deciding authority: [nrpy/helpers/simd_intrinsics.h](../../../nrpy/helpers/simd_intrinsics.h), all instruction-set branches.

Debug mode is an internal validation path for the symbolic rewrite. `_debug_eval_namespace()` builds an explicit namespace containing SymPy objects, free symbols, and `*_SIMD_check` functions such as `AbsSIMD_check`, `ConstSIMD_check`, `FusedMulAddSIMD_check`, and `NegFusedMulSubSIMD_check`. With `debug=True`, `expr_convert_to_simd_intrins()` evaluates the transformed spelling through that restricted namespace, back-substitutes rational symbols, compares against the original expression, and raises `ValueError` if simplification cannot prove the difference is zero. This documents the helper's own validation behavior, not a general recommendation to evaluate arbitrary expression strings.

`simd_intrinsics.h` defines vector types, widths, loads/stores, broadcasts, arithmetic, selected math functions, horizontal addition, and `UPWIND_ALG`. Compiler macros select AVX512F, AVX, SSE2/SSE3, or scalar fallback. The supported scalar macros use width one. Header arithmetic still needs separate validation: without FMA, the AVX and SSE `NegFusedMulSubSIMD(a,b,c)` expansion computes `c-(c+(a*b+c))`, which can overflow when `-(a*b+c)` remains finite. The SSE branch with `__FMA__` defined instead computes `a*b-c`, contrary to the symbolic checker's `-(a*b+c)`; this branch is not ordinarily selected when AVX is enabled. Symbolic debug checks do not execute these headers or establish floating-point equivalence.

Claim evidence:
- Claim: SIMD header coverage and arithmetic differ from the symbolic checker; non-FMA negated subtraction can overflow differently and the SSE/FMA form has a different sign expression.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/simd_intrinsics.h](../../../nrpy/helpers/simd_intrinsics.h), `NegFusedMulSubSIMD`
- Corroboration: `none available`; the named implementation decides this behavior.
- Corroborating comparison: [nrpy/helpers/simd.py](../../../nrpy/helpers/simd.py), `NegFusedMulSubSIMD_check`, defines the symbolic intended expression.

`cuda_intrinsics.h` is the CUDA-oriented scalar/intrinsic contract. When `__CUDACC__` is defined, `AddCUDA`, `SubCUDA`, and `MulCUDA` are plain `+`, `-`, and `*`, so that nvcc's default `-fmad=true` contraction can fuse multiply-add pairs the symbolic lowering did not fuse itself; the opaque `__dadd_rn`/`__dsub_rn`/`__dmul_rn` intrinsics they replaced are never contracted, and each plain operator is still IEEE round-to-nearest. `DivCUDA`, `FusedMulAddCUDA`, `SqrtCUDA`, and `ReadCUDA` use CUDA device intrinsics or device read helpers, and the fused subtract and negated forms fold their exact negations into a single `__fma_rn`. The non-CUDA branch keeps the same `REAL_CUDA_ARRAY`, `CUDA_WIDTH`, arithmetic, load/store, horizontal-add, and `UPWIND_ALG` macro names as scalar C expressions. This header is separate from the SIMD symbolic transformer: it defines CUDA macro spellings and fallback behavior, but `expr_convert_to_simd_intrins()` itself emits `*SIMD` names, not `*CUDA` names.

Generated projects may contain copied versions of these headers, but those copies are generated output. For KB purposes, cite the handwritten sources under `nrpy/helpers/` and the codegen integration in `nrpy/c_codegen.py`, not copies under `project/`.

## Sources

- [nrpy/helpers/simd.py](../../../nrpy/helpers/simd.py) - `expr_convert_to_simd_intrins`, `_debug_eval_namespace`, `AbsSIMD_check`, `ConstSIMD_check`, `FusedMulAddSIMD_check`, `FusedMulSubSIMD_check`, `NegFusedMulAddSIMD_check`, `NegFusedMulSubSIMD_check`
- [nrpy/helpers/simd_intrinsics.h](../../../nrpy/helpers/simd_intrinsics.h) - `REAL_SIMD_ARRAY`, `SIMD_WIDTH`, `UPWIND_ALG`, `FusedMulAddSIMD`, `HorizAddSIMD`
- [nrpy/helpers/cuda_intrinsics.h](../../../nrpy/helpers/cuda_intrinsics.h) - `REAL_CUDA_ARRAY`, `CUDA_WIDTH`, `UPWIND_ALG`, `FusedMulAddCUDA`, `ReadCUDA`
- [nrpy/c_codegen.py](../../../nrpy/c_codegen.py) - `CCodeGen`, `c_codegen`, `gridfunction_management_and_FD_codegen`
- [Intel Intrinsics Guide](https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html) - SSE, AVX, AVX512, and intrinsic-family terminology
- [NVIDIA CUDA Programming Guide](https://docs.nvidia.com/cuda/cuda-programming-guide/index.html) - CUDA programming-model terminology from version 13.3

## See Also

- Parent: [Helper APIs](index.md)
- See also: [CSE And Printer Support](cse-and-printer-support.md)
- See also: [Loop Kernel And Device Helpers](loop-kernel-and-device-helpers.md)
- See also: [C Codegen](../c-codegen.md)
- See also: [Finite Difference](../finite-difference.md)
- See also: [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
