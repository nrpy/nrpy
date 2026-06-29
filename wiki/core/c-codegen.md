# C Codegen

> Core route for turning SymPy expressions into generated C text. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Core APIs](index.md)

## Summary

`nrpy.c_codegen` is the main bridge from symbolic NRPy expressions to C code strings. The public `c_codegen()` function normalizes expression and target-name inputs, builds a `CCodeGen` option object, emits C assignments, and delegates to the finite-difference/gridfunction path when requested. That path is deliberately tied to `nrpy.finite_difference`, so generated kernels can read registered gridfunctions, compute derivative stencils, handle upwinded derivatives, and then write expression results.

## Detail

`CCodeGen` stores the codegen configuration used by `c_codegen()`: floating-point type and alias, CSE behavior, SIMD flags, finite-difference flags, memory layout style, upwind control vector, rational-constant handling, and optional clang-format output. The constructor validates supported floating-point modes and switches aliases according to the active infrastructure, for example using `REAL` for BHaH and `REAL_SIMD_ARRAY` for SIMD-capable BHaH, CarpetX, or ETLegacy targets.

`c_codegen()` accepts one SymPy expression or a list of expressions, paired with one output variable name or a list. It rejects tuple inputs, checks expression/name length agreement, and can emit verbose comments containing the original symbolic expressions. With ordinary symbolic codegen it applies CSE and SymPy C printing according to `CCodeGen` options.

When `automatically_read_gf_data_from_memory` or `enable_fd_codegen` is enabled, `c_codegen()` first extracts derivative symbols from expression free symbols, maps derivative symbols to base gridfunctions and derivative operators, computes finite-difference coefficients and stencils, asks `read_gfs_from_memory()` to generate the needed memory reads, and calls `gridfunction_management_and_FD_codegen()`. That helper emits up to three generated sections: gridfunction reads plus stencil arithmetic, upwind-selection arithmetic, and final expression evaluation. Optional `enable_fd_functions` populates finite-difference helper functions through `FDFunction` objects instead of inlining every derivative expression.

Use this page for the codegen interface itself. Use the finite-difference page when the question is about derivative naming, stencil construction, or memory-read ordering.

## Sources

- [nrpy/c_codegen.py](../../nrpy/c_codegen.py) - `CCodeGen`, `c_codegen`, `gridfunction_management_and_FD_codegen`
- [nrpy/finite_difference.py](../../nrpy/finite_difference.py) - `compute_fdcoeffs_fdstencl`, `read_gfs_from_memory`, `FDFunction`
- [coding_style.md](../../coding_style.md) - SymPy to C and BHaH symbolic codegen rules

## See Also

- [Core APIs](index.md)
- [Finite Difference](finite-difference.md)
- [C Function Registry](c-function-registry.md)
