# Indexed Expressions

> Core route for tensor-like symbolic arrays and symmetry-aware declarations. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [Core APIs](index.md)

## Summary

`nrpy.indexedexp` is NRPy's standard tensor-expression helper module. It creates nested SymPy-symbol arrays for rank-1 through rank-4 indexed expressions, applies supported symmetry strings, initializes zero tensors, infers nested-list rank, provides explicit 2x2 through 4x4 matrix inversion helpers, and builds Levi-Civita symbols or tensors. Equation and gridfunction code use it instead of ad hoc list construction so names, rank, dimension, and symmetry stay consistent.

## Detail

`create_tensor_symbolic()` is the low-level constructor. Its `shape` argument may be a single integer for rank-1 output or a list of per-axis extents for nested output. When `symbol` is nonempty, each component is named from the base plus component suffixes; `preindex` prepends numeric suffix pieces during recursive construction, and `character_zero_index` switches suffixes from digits to alphabetic labels starting at the requested character, for example `x`, `y`, `z` or `A`, `B`, `C`. When no symbol is supplied, components are initialized to symbolic zero.

The central public implementation is `declare_indexedexp()`. It validates rank and dimension, rejects deprecated `DIM=` in favor of `dimension=`, accepts only symbol basenames made from alphanumerics plus `_`, `-`, and `>`, builds symbols with `create_tensor_symbolic()`, parses optional symmetry strings, applies symmetry or antisymmetry, and then applies symmetry-axis derivative zeroing. The public wrappers `declarerank1()`, `declarerank2()`, `declarerank3()`, and `declarerank4()` set the rank and delegate to `declare_indexedexp()`.

Supported symmetry syntax is deliberately strict. The accepted tokens are `nosym`, `sym...`, and `anti...`; `nosym` must appear alone. Non-`nosym` tokens are joined by underscores, so `sym01_sym23` is valid and concatenated forms such as `sym01sym23` are rejected. Token indices must be digits, unique, ascending, and inside the rank range. Rank-2 tokens use exactly two indices; higher-rank tokens may use multi-index groups such as `sym012` or `anti0123`.

Symmetry expansion reduces multi-index groups to pairwise operations. For example, a group over indices `012` expands through adjacent pair operations so the rank-specific symmetrizers can reuse the same pairwise cases. Antisymmetric operations mirror components with a negative sign and set repeated-index components to zero.

The `indexedexp::symmetry_axes` parameter controls post-declaration derivative zeroing. `zero_out_derivatives_across_symmetry_axes()` scans component names for `_dD` or `_dDD` suffixes and zeros first- or second-derivative components whose derivative direction lies across a configured symmetry axis. This pass supports ranks 1 through 4. Names that indicate derivative order greater than 2 raise `ValueError`.

Zero initializers `zerorank1()` through `zerorank4()` call `declare_indexedexp(None, ...)`, producing mutable nested arrays initialized to symbolic zero. `get_rank()` infers rank by walking nested lists until it reaches a non-list element, so callers can classify tensor-like structures without separate metadata.

Matrix inversion helpers are split into symmetric and generic families: `symm_matrix_inverter2x2()`, `symm_matrix_inverter3x3()`, `symm_matrix_inverter4x4()`, and `generic_matrix_inverter2x2()`, `generic_matrix_inverter3x3()`, `generic_matrix_inverter4x4()`. They compute explicit determinant and inverse formulas and raise `NonInvertibleMatrixError` when the determinant is zero.

Levi-Civita helpers include `LeviCivitaSymbol_dim3_rank3()`, `LeviCivitaTensorUUU_dim3_rank3()`, and `LeviCivitaTensorDDD_dim3_rank3()`. The symbol helper builds the rank-3 3D alternating symbol. The UUU tensor divides by the supplied determinant square root, while the DDD tensor multiplies by it.

Use these helpers with explicit nested loops. The project style explicitly avoids Einstein summation notation and matrix-multiply operators in equation construction.

## Sources

- [nrpy/indexedexp.py](../../nrpy/indexedexp.py) - `create_tensor_symbolic`, `declare_indexedexp`, `_parse_symmetry_tokens`, `_expanded_symmetry_tokens`
- [nrpy/indexedexp.py](../../nrpy/indexedexp.py) - `declarerank1`, `declarerank2`, `declarerank3`, `declarerank4`, `zero_out_derivatives_across_symmetry_axes`, `get_rank`
- [nrpy/indexedexp.py](../../nrpy/indexedexp.py) - `zerorank1`, `zerorank2`, `zerorank3`, `zerorank4`
- [nrpy/indexedexp.py](../../nrpy/indexedexp.py) - `symm_matrix_inverter2x2`, `symm_matrix_inverter3x3`, `symm_matrix_inverter4x4`, `generic_matrix_inverter2x2`, `generic_matrix_inverter3x3`, `generic_matrix_inverter4x4`, `NonInvertibleMatrixError`
- [nrpy/indexedexp.py](../../nrpy/indexedexp.py) - `LeviCivitaSymbol_dim3_rank3`, `LeviCivitaTensorUUU_dim3_rank3`, `LeviCivitaTensorDDD_dim3_rank3`

## See Also

- [Core APIs](index.md)
- [Gridfunctions And Parameters](gridfunctions-and-parameters.md)
- [Finite Difference](finite-difference.md)
- [Reference Metrics](reference-metrics.md)
