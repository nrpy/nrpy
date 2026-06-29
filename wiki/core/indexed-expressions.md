# Indexed Expressions

> Core route for tensor-like symbolic arrays and symmetry-aware declarations. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Core APIs](index.md)

## Summary

`nrpy.indexedexp` is NRPy's standard tensor-expression helper module. It creates nested SymPy-symbol arrays for rank-1 through rank-4 indexed expressions, applies supported symmetry strings, initializes zero tensors, provides symmetric matrix inversion helpers, and builds Levi-Civita symbols or tensors. Equation and gridfunction code use it instead of ad hoc list construction so names, rank, dimension, and symmetry stay consistent.

## Detail

The central implementation is `declare_indexedexp()`. It validates the requested rank and dimension, builds symbols with component suffixes, parses optional symmetry strings, applies symmetry or antisymmetry, and then applies symmetry-axis derivative zeroing when the `indexedexp::symmetry_axes` parameter is set. The public wrappers `declarerank1()`, `declarerank2()`, `declarerank3()`, and `declarerank4()` set the rank and delegate to `declare_indexedexp()`.

Supported symmetry syntax is token based. Rank-2 declarations can use strings such as `sym01`, `anti01`, or `nosym`. Higher ranks accept underscore-separated tokens such as `sym12`, `sym012`, or `sym01_sym23`; malformed concatenations are rejected with errors that point to the expected token form. The coding rules use these strings for common cases: `sym01` for metrics, `sym12` for derivatives of symmetric tensors, `sym01_sym23` for Riemann-like tensors, and `nosym` when no symmetry is intended.

Zero initializers `zerorank1()` through `zerorank4()` call `declare_indexedexp(None, ...)`, producing mutable nested arrays initialized to symbolic zero. These are the preferred accumulators for loop-built tensor expressions.

Matrix inversion helpers include symmetric inverters such as `symm_matrix_inverter3x3()`, which computes inverse and determinant formulas directly for symmetric matrices and raises `NonInvertibleMatrixError` on zero determinant. Levi-Civita helpers include `LeviCivitaSymbol_dim3_rank3()` and `LeviCivitaTensorUUU_dim3_rank3()`, the latter dividing the symbol by the supplied metric determinant square root.

Use these helpers with explicit nested loops. The project style explicitly avoids Einstein summation notation and matrix-multiply operators in equation construction.

## Sources

- [nrpy/indexedexp.py](../../nrpy/indexedexp.py) - `declare_indexedexp`, `declarerank1`, `declarerank2`, `declarerank3`, `declarerank4`
- [nrpy/indexedexp.py](../../nrpy/indexedexp.py) - `zerorank1`, `zerorank2`, `zerorank3`, `zerorank4`, `symm_matrix_inverter3x3`
- [nrpy/indexedexp.py](../../nrpy/indexedexp.py) - `LeviCivitaSymbol_dim3_rank3`, `LeviCivitaTensorUUU_dim3_rank3`
- [coding_style.md](../../coding_style.md) - indexed-expression conventions and expression construction pattern

## See Also

- [Core APIs](index.md)
- [Gridfunctions And Parameters](gridfunctions-and-parameters.md)
- [Reference Metrics](reference-metrics.md)
