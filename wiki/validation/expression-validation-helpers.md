# Expression Validation Helpers

> Numerical SymPy-expression validation helpers for equality, zero checks, and trusted dictionaries. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Validation](index.md)

## Summary

`nrpy/validate_expressions/validate_expressions.py` validates SymPy expressions
by evaluating them at high precision with `mpmath`. This is high-precision
numerical validation, not a formal symbolic proof. The same helpers support
equality assertions, zero checks, trusted dictionary generation, and trusted
dictionary comparison.

## Detail

The module uses `precision = 30` as the standard decimal precision. Relative
comparison tolerance is `10 ** (-4.0 / 5.0 * precision)` where the helper uses
the module constant, or `10 ** (-4.0 / 5.0 * mp.dps)` where it uses the current
`mpmath` decimal precision.

`assert_equal()` accepts either dictionaries or single expressions. Non-dict
inputs are `sympify`-wrapped into one-entry dictionaries keyed by `""`, both
sides are processed with fixed `mpf` substitutions, and each paired numerical
result is compared by relative error. A mismatch prints the failing key and
error bound, then raises `AssertionError`; a passing check prints
`Assertion Passed!` unless messages are suppressed.

`check_zero()` sends one expression through `convert_one_expression_to_mpfmpc()`
and returns whether the final numerical result is exactly `mp.mpf("0.0")`.
That exact-zero comparison happens after the conversion helper has applied its
near-zero retry behavior.

Free-symbol substitution has two modes. Fixed mode seeds Python `random` per
symbol with
`random.seed(int(hashlib.md5(str(var).encode()).hexdigest(), 16 + hex_offset))`.
Here `hex_offset` changes the integer base argument passed to `int(...)`; it
does not change the MD5 string or the symbol string being hashed. Non-fixed
mode uses ordinary random values. Symbols named `PI` or `M_PI` are assigned
`mp.pi`; symbols named `SQRT1_2` or `M_SQRT1_2` are assigned `1 / sqrt(2)`.

Single-expression conversion first sets `mp.dps = precision` and returns
`mpf(0.0)` immediately for exact SymPy zero. Otherwise it runs
`sp.cse(expr, order="none")`, collects free symbols from the reduced expression
and common subexpressions, injects `mpf` values with `xreplace`, and evaluates
common-subexpression replacements at the current `mp.dps`. The reduced
expression is converted to `mpf` first. If that raises `TypeError`, an
`sp.nan` result takes the diagnostic print path; other non-real results fall
back to `mpc(sp.N(..., mp.dps))`.

Near-zero handling is a retry, not symbolic simplification. If a nonzero result
has magnitude below `10 ** (-4.0 / 5.0 * mp.dps)`, the helper reruns the same
common-subexpression form at twice the standard precision. It collapses the
result to `mp.mpf("0.0")` only if the higher-precision result is still below
the tolerance at that higher `mp.dps`; otherwise it keeps the original nonzero
value and restores `mp.dps` to `precision`.

`process_dictionary_of_expressions()` sorts dictionary items, skips keys
containing `funcform`, converts scalar SymPy expressions directly, and flattens
tensor lists recursively. Flattened tensor entries are named
`<lhs>_<flat_index>` before conversion to `mpf` or `mpc`.

Trusted output is caller-relative. `output_trusted()` derives a sibling
`tests/<trusted_file_basename>.py` path from `Path(os_path_abspath)` relative to
`os_getcwd`, writes only the `mpmath` imports needed by the generated values,
defines `trusted_dict`, formats the file with Black, and creates the `tests/`
directory if needed. `compare_or_generate_trusted_results()` writes that file
when missing and compares against it when present.

Trusted comparison imports the caller-relative module path with `importlib`
and expects a `trusted_dict`. It first checks only the dictionary lengths; when
the lengths differ, it reports the missing key set and tells maintainers how to
regenerate a trusted file. When the lengths match, it iterates trusted keys and
compares each value by relative error against `results_dict[key]`; a changed
key with unchanged count will surface as `KeyError`.

## Sources

- [nrpy/validate_expressions/validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `precision`, `assert_equal`, `check_zero`
- [nrpy/validate_expressions/validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `convert_free_symbols_set_to_mpf_dict`, `convert_one_expression_to_mpfmpc`, `inject_mpfs_into_cse_expression`
- [nrpy/validate_expressions/validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `process_dictionary_of_expressions`, `output_trusted`, `compare_against_trusted`, `compare_or_generate_trusted_results`

## See Also

- [Validation](index.md)
- [Trusted Expression Pipeline](../equations/trusted-expression-pipeline.md)
- [Static Analysis](static-analysis.md)
- [Reference Metrics](../core/reference-metrics.md)
