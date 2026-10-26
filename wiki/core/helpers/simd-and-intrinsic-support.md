# SIMD And Intrinsic Support

> Helper route for symbolic SIMD rewrites and handwritten intrinsic macro headers. · Status: contested
> Up: [Helper APIs](index.md)

## Summary

`nrpy.helpers.simd` rewrites SymPy expressions into symbolic calls such as `AddSIMD`, `MulSIMD`, `DivSIMD`, and `FusedMulAddSIMD`. Handwritten headers supply C/CUDA implementations for supported names. The transformation and header must be checked together: some emitted names and constants are missing definitions. Fused and ordinary arithmetic can round differently.

## Detail

`expr_convert_to_simd_intrins()` is the transformation entry point. If the caller does not supply a `symbol_to_Rational_dict`, it first calls `cse_preprocess()` so rational constants can be represented by generated symbols before the SIMD pass. The function then builds an `ExprTree` over the expression, creates SymPy function objects for the SIMD spelling, broadcasts each real number that preprocessing leaves as a scalar (a float, `pi`, `sqrt(2)`) as `ConstSIMD(<C double literal>)` in Step 0 of `simd.py`, except that Step 0 skips the exponent of each power that is neither skipped nor wrapped itself. The function next rewrites the tree in later stages: Step 1 replaces absolute value (`Abs` and `nrpyAbs` become `AbsSIMD`), the functions `exp`, `log`, `sin`, `cos`, and `sign`, and powers and selected roots, and the stages after Step 1 replace subtraction patterns, binary addition and multiplication trees, division patterns, fused multiply-add/subtract patterns, and negative fused patterns, followed by optional cleanup of an unused prefixed `_NegativeOne_` constant. `c_codegen()` skips that cleanup for the final expressions when any common subexpression reads `_NegativeOne_`. A subtraction rewrite removes the two matched summands by position, so a repeated summand is kept. The option `simd_find_more_FMAsFMSs`, true by default in both `expr_convert_to_simd_intrins()` and `CCodeGen`, enables extra searches for fused multiply-add/subtract patterns; for example, `a + b*c + d*e` gives `FusedMulAddSIMD(b, c, FusedMulAddSIMD(d, e, a))` with the option true, whereas with it false the product `d*e` stays unfused, as in `AddSIMD(a, FusedMulAddSIMD(b, c, MulSIMD(d, e)))`.

Claim evidence:
- Claim: `expr_convert_to_simd_intrins()` and `CCodeGen` take the option `simd_find_more_FMAsFMSs`, true by default; it enables extra fused multiply-add/subtract searches in `simd.py`, so `a + b*c + d*e` gives `FusedMulAddSIMD(b, c, FusedMulAddSIMD(d, e, a))` with the option true and `AddSIMD(a, FusedMulAddSIMD(b, c, MulSIMD(d, e)))` with it false; the option does not guarantee more fused operations for other expressions.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/simd.py](../../../nrpy/helpers/simd.py), `expr_convert_to_simd_intrins`; [nrpy/c_codegen.py](../../../nrpy/c_codegen.py), `CCodeGen`
- Corroboration: `none available`; the named implementation decides this behavior.

The SIMD pass is algebraic and symbolic. It does not choose an instruction set, allocate vector registers, or include a header. It produces expression text containing names such as `SqrtSIMD`, `CosSIMD`, `PowSIMD`, `FusedMulSubSIMD`, `NegFusedMulAddSIMD`, and `NegFusedMulSubSIMD`; later C code generation prints those names into generated code. `nrpy.c_codegen.CCodeGen` is the integration point that enables this path through `enable_simd`, restricts SIMD output to double precision, uses `REAL_SIMD_ARRAY` for supported infrastructures, calls `expr_convert_to_simd_intrins()` on common subexpressions and final expressions, and emits SIMD constant declarations through `ConstSIMD`.

Reciprocal half-power rewrites for exponents `-1/2`, `-3/2`, and `-5/2`, and negative integer powers, use the prefixed `_Integer_1` numerator and register it in `symbol_to_Rational_dict` when it is not already present, so declarations derived from that dictionary include it. The Step 0 pass skips the exponent of each power that is neither skipped nor wrapped itself, so that exponent is not wrapped as one constant: when Step 1 emits `PowSIMD` for such a power, it emits the exponent as a C double literal if the exponent is a Float, a Rational, a SymPy `NumberSymbol` such as `pi` or `E`, or a real root of a rational that SymPy stores as one power, such as `sqrt(2)`, as in `PowSIMD(x, 0.3)` for `x**0.3`; SymPy stores `sqrt(8)` as the product `2*sqrt(2)`, which is rewritten as in the next sentence. The nodes inside a skipped exponent are still processed, and a skipped exponent that combines numbers with other operations or functions, such as `-pi`, `2*pi`, `sqrt(pi)`, or `log(2)`, is not wrapped as one constant and its parts are rewritten as in any other expression, as in `PowSIMD(x, MulSIMD(_Integer_2, ConstSIMD(<pi as a C double literal>)))` for `x**(2*pi)`. A skipped power does not skip its own exponent, so `x**(y**pi)` gives `PowSIMD(x, PowSIMD(y, ConstSIMD(<pi as a C double literal>)))`. A rational or integer base under a symbolic exponent stays as written, as in `PowSIMD(1/2, x)`.

Claim evidence:
- Claim: `expr_convert_to_simd_intrins()` registers `_Integer_1` for the reciprocal half-power and negative-integer-power rewrites, and wraps real non-rational numbers as `ConstSIMD(<C double literal>)`, except that it skips the exponent of each power that is neither skipped nor wrapped itself.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/simd.py](../../../nrpy/helpers/simd.py), `expr_convert_to_simd_intrins`
- Corroboration: `none available`; the named implementation decides this behavior.

`PowSIMD`, `CbrtSIMD`, `LogSIMD`, and `SignSIMD` can be emitted but have no definitions in `simd_intrinsics.h`. Including that header alone is insufficient for those expressions; supply suitable definitions or use a supported expression form.

Claim status: contested; contradiction: CONTR-0015. See [CONTR-0015](../../contradictions.md#contr-0015).

Claim evidence:
- Claim: Emitted Pow/Cbrt/Log/Sign SIMD names lack definitions in the SIMD header.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/simd.py](../../../nrpy/helpers/simd.py), `expr_convert_to_simd_intrins`; [nrpy/helpers/simd_intrinsics.h](../../../nrpy/helpers/simd_intrinsics.h), all instruction-set branches.
- Corroboration: `none available`; the named implementation and header decide this behavior.

Debug mode is an internal validation path for the symbolic rewrite. `_debug_eval_namespace()` builds an explicit namespace containing SymPy objects, free symbols, and `*_SIMD_check` functions such as `AbsSIMD_check`, `ConstSIMD_check`, `FusedMulAddSIMD_check`, and `NegFusedMulSubSIMD_check`. With `debug=True`, `expr_convert_to_simd_intrins()` evaluates the transformed spelling through that restricted namespace, back-substitutes rational symbols, compares against the original expression, and raises `ValueError` if simplification cannot prove the difference is zero. Real constants such as `pi` or `sqrt(2)` enter the transformed expression as rounded double literals, so debug mode reports their difference from the exact value as `ValueError`. This documents the helper's own validation behavior, not a general recommendation to evaluate arbitrary expression strings.

`simd_intrinsics.h` defines vector types, widths, loads/stores, broadcasts, arithmetic, selected math functions, horizontal addition, and `UPWIND_ALG`. Compiler macros select AVX512F, AVX, SSE2/SSE3, or scalar fallback. The supported scalar macros use width one.

`NegFusedMulSubSIMD(a,b,c)` represents `-(a*b+c)`. The non-FMA SSE/AVX forms multiply, add `c`, and flip the sum's sign bit, preserving unary negation of signed zero without redundant arithmetic involving `c`. The SSE/FMA, AVX/FMA, and AVX512 forms use their corresponding `fnmsub` intrinsics, equivalent to `fma(-a,b,-c)`. With contraction disabled, ordinary arithmetic rounds the product and sum separately; fused arithmetic rounds once. These forms can differ in rounding, overflow, and the zero sign after exact cancellation. The SSE/FMA branch is normally bypassed because AVX takes precedence in header selection. Symbolic debug checks do not execute these headers or establish floating-point equivalence.

Claim evidence:
- Claim: NegFusedMulSubSIMD represents -(a*b+c): non-FMA SSE/AVX multiply, add c, and flip the sign bit; SSE/FMA, AVX/FMA, and AVX512 use fnmsub, equivalent to fma(-a,b,-c). Fused and ordinary forms can differ in rounding, overflow, and the zero sign after exact cancellation.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/simd_intrinsics.h](../../../nrpy/helpers/simd_intrinsics.h), `NegFusedMulSubSIMD`
- Corroboration: [nrpy/helpers/simd.py](../../../nrpy/helpers/simd.py), `NegFusedMulSubSIMD_check`, supplies the symbolic sign expression; it does not execute the header or validate floating-point results.

`cuda_intrinsics.h` is the header that defines the CUDA-oriented scalar and intrinsic macros. When `__CUDACC__` is defined, `AddCUDA`, `SubCUDA`, and `MulCUDA` are plain `+`, `-`, and `*`, so that nvcc's default `-fmad=true` contraction can fuse multiply-add pairs the symbolic lowering did not fuse itself; the opaque `__dadd_rn`/`__dsub_rn`/`__dmul_rn` intrinsics they replaced are never contracted, and each plain operator is still IEEE round-to-nearest. `DivCUDA`, `FusedMulAddCUDA`, `SqrtCUDA`, and `ReadCUDA` use CUDA device intrinsics or device read helpers, and the fused subtract and negated forms fold their exact negations into a single `__fma_rn`. The non-CUDA branch keeps the same `REAL_CUDA_ARRAY`, `CUDA_WIDTH`, arithmetic, load/store, horizontal-add, and `UPWIND_ALG` macro names as scalar C expressions. This header is separate from the SIMD symbolic transformer: it defines CUDA macro spellings and fallback behavior, but `expr_convert_to_simd_intrins()` itself emits `*SIMD` names, not `*CUDA` names.

Generated projects may contain copied versions of these headers, but those copies are generated output. For KB purposes, cite the handwritten sources under `nrpy/helpers/` and the codegen integration in `nrpy/c_codegen.py`, not copies under `project/`.

## Sources

- [nrpy/helpers/simd.py](../../../nrpy/helpers/simd.py) - `expr_convert_to_simd_intrins`, `_debug_eval_namespace`, `AbsSIMD_check`, `ConstSIMD_check`, `FusedMulAddSIMD_check`, `FusedMulSubSIMD_check`, `NegFusedMulAddSIMD_check`, `NegFusedMulSubSIMD_check`
- [nrpy/helpers/simd_intrinsics.h](../../../nrpy/helpers/simd_intrinsics.h) - `REAL_SIMD_ARRAY`, `SIMD_WIDTH`, `UPWIND_ALG`, `FusedMulAddSIMD`, `NegFusedMulSubSIMD`, `HorizAddSIMD`
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
