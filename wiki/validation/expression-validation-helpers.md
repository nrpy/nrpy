# Expression Validation Helpers

> Numerical SymPy-expression validation helpers for equality, zero checks, and trusted dictionaries. · Status: confirmed · Last reconciled: 07-13-2026
> Up: [Validation](index.md)

## Summary

`nrpy/validate_expressions/validate_expressions.py` validates SymPy expressions
by evaluating them at high precision with `mpmath`. This is high-precision
numerical spot validation, not a formal symbolic proof, runtime test, generated
project test, or general scientific-accuracy guarantee. The same helpers support
equality assertions, zero checks, trusted dictionary generation, and trusted
dictionary comparison under the limits below.

## Detail

Prefer an exact analytic, symbolic, or semantic invariant when practical. When
sampled regression fits the contract, follow the caller safeguards in [Test
Oracles And Safe Updates](test-oracles-and-safe-updates.md); this page retains
the current helper mechanics and their limits.

The module uses `precision = 30` as the standard decimal precision. Relative
comparison tolerance is `10 ** (-4.0 / 5.0 * precision)` where the helper uses
the module constant, or `10 ** (-4.0 / 5.0 * mp.dps)` where it uses the current
`mpmath` decimal precision.

The same module contains doctests for `assert_equal()`, `check_zero()`, and
exact-zero conversion, and its `__main__` block runs them. No separate dedicated
test module is checked in under `nrpy/validate_expressions/`; broader equation
modules exercise these helpers through their own trusted-result paths.

`assert_equal()` accepts either dictionaries or single expressions. Non-dict
inputs are `sympify`-wrapped into one-entry dictionaries keyed by `""`, both
sides are processed with fixed `mpf` substitutions, and values at each common
processed key are compared by relative error. Dictionary values must be SymPy
expressions or recursively nested lists of them; unsupported leaves raise
`TypeError`. Dictionary inputs must have identical raw key sets and list
nesting; extra, missing, renamed, or differently shaped entries raise
`AssertionError` before value comparison. Among entries retained for numerical
processing, a scalar name that would collide with a flattened tensor-leaf name
also raises. Keys containing `funcform` participate in raw key, nesting, and
leaf validation but are intentionally omitted from numerical comparison and
processed-name collision checks. There is no positional dictionary mode.

For finite values, a mismatch prints the failing key and error bound, then
raises `AssertionError`; a passing check prints `Assertion Passed!` unless
messages are suppressed. Non-finite components are handled before subtraction:
matching NaNs in the same real or imaginary component and same-signed
infinities compare equal, while one-sided NaNs, finite-versus-non-finite values,
opposite infinities, or different finite companion components fail.

Claim evidence:
- Claim: `assert_equal()` requires SymPy-expression leaves, supports recursively nested lists, rejects non-identical raw key sets and different list nesting, rejects flattened-name collisions among numerically processed entries, omits structurally valid `funcform` keys from numerical comparison and collision checks, and explicitly distinguishes matching from mismatched non-finite components.
- Role: descriptive behavior
- Deciding authority: [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py), `assert_equal` and `_nonfinite_values_match`
- Corroboration: [test_parse_BSSN.py](../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py), `test_example_BSSN`, exercises direct dictionary comparison across scalar, vector, and matrix expression values; filter and error-path tests remain colocated with the deciding helper
- Validation: `inspected=pass; generated=not-run; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0, mpmath 1.3.0; backend=not-applicable; precision=30 decimal digits; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=pass; options=68 validator doctests, explicit funcform collision and leaf-type probes, and BSSN cross-representation run; date=07-13-2026`

`check_zero()` sends one expression through `convert_one_expression_to_mpfmpc()`
and returns whether the final numerical result is exactly `mp.mpf("0.0")`.
That exact-zero comparison happens after the conversion helper has applied its
near-zero retry behavior. Regression calls with free symbols must pass
`fixed_mpfs_for_free_symbols=True`; symbol-free results remain numerical
evaluations, not formal identity proofs.

Free-symbol substitution has two modes. Fixed mode seeds Python `random` per
symbol with
`random.seed(int(hashlib.md5(str(var).encode()).hexdigest(), 16 + hex_offset))`.
Here `hex_offset` changes the integer base argument passed to `int(...)`; it
does not change the MD5 string or the symbol string being hashed. Non-fixed
mode uses ordinary random values. Symbols named `PI` or `M_PI` are assigned
`mp.pi`; symbols named `SQRT1_2` or `M_SQRT1_2` are assigned `1 / sqrt(2)`.
Fixed mode is repeatable only under documented current runtime assumptions; it
does not guarantee bitwise stability across Python, SymPy, mpmath, platform, or
codegen versions.

Single-expression conversion first sets `mp.dps = precision` and returns
`mpf(0.0)` immediately for exact SymPy zero. Otherwise it runs
`sp.cse(expr, order="none")`, collects free symbols from the reduced expression
and common subexpressions, injects `mpf` values with `xreplace`, and evaluates
common-subexpression replacements at the current `mp.dps`. The reduced
expression is converted to `mpf` first. If that raises `TypeError`, an
`sp.nan` result takes the diagnostic print path, positive and negative SymPy
infinities convert explicitly to their `mpmath` counterparts, and other
non-real results are split into real and imaginary components before each
finite, NaN, or infinity component is converted and combined as `mpc`. On the
recorded SymPy 1.14.0/mpmath 1.3.0 stack, both signed SymPy infinities raise
`TypeError` in the initial `mpf` conversion and exercise the explicit fallback.

Claim evidence:
- Claim: `inject_mpfs_into_cse_expression()` first attempts `mpf` conversion; SymPy NaN uses the diagnostic fallback, on the recorded SymPy 1.14.0/mpmath 1.3.0 stack signed SymPy infinities raise `TypeError` and use the explicit signed-infinity fallback, and other non-real fallback values are converted component by component before constructing an `mpc` value.
- Role: descriptive behavior
- Deciding authority: [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py), `inject_mpfs_into_cse_expression`
- Corroboration: none available; direct signed-infinity and `assert_equal()` doctests are colocated with the deciding helper, and no separate source exercises every conversion branch
- Validation: `inspected=pass; generated=not-run; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0, mpmath 1.3.0; backend=not-applicable; precision=30 decimal digits; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=pass; options=direct signed-infinity fallback doctests plus aggregate NaN and complex-fallback cases in 68 total validator doctests; date=07-13-2026`

Near-zero handling is a retry, not symbolic simplification. If a nonzero result
has magnitude below `10 ** (-4.0 / 5.0 * mp.dps)`, the helper reruns the same
common-subexpression form at twice the standard precision. The retry does not
forward `hex_offset`; in non-fixed mode it also draws new random substitutions.
It therefore need not evaluate the same numerical point. The helper collapses
the result to `mp.mpf("0.0")` only if the higher-precision result is still below
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

Missing-file generation creates a baseline; it does not independently establish
that baseline as scientifically correct or count as a successful comparison.
Candidate review and acceptance follow the isolated two-process procedure in
[Test Oracles And Safe Updates](test-oracles-and-safe-updates.md). Existing-file
comparison checks regression against stored numbers only under the chosen
substitutions and tolerance.

Trusted comparison imports the caller-relative module path with `importlib`
and expects a `trusted_dict`. It first checks dictionary lengths; when the
lengths differ, it reports the missing key set and tells maintainers how to
regenerate a trusted file. When lengths match, it iterates trusted keys and
compares each value against `results_dict[key]`; a changed key with unchanged
count surfaces as `KeyError`. Finite values use relative error. Non-finite
values use the same component-aware rule as `assert_equal()`: intentional NaN
sentinels or same-signed infinities pass only when both sides match, and a
computed NaN against a finite trusted value, or a finite result against a
trusted NaN, raises `ValueError`.

Claim evidence:
- Claim: `compare_against_trusted()` uses relative error for finite values, accepts matching NaN components and same-signed infinities, and raises `ValueError` for finite/non-finite or differently structured non-finite mismatches.
- Role: descriptive behavior
- Deciding authority: [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py), `compare_against_trusted` and `_nonfinite_values_match`
- Corroboration: [reference_metric_GeneralRFM_fisheyeN2.py](../../nrpy/tests/reference_metric_GeneralRFM_fisheyeN2.py), `trusted_dict`, supplies three intentional trusted NaN sentinels exercised by the reference-metric owner run
- Validation: `inspected=pass; generated=not-run; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0, mpmath 1.3.0; backend=not-applicable; precision=30 decimal digits; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=pass; options=full reference_metric validation plus both trusted-NaN mismatch doctests; date=07-13-2026`

All expression checks sample numerical substitutions rather than proving an
identity over a domain. A sampled point can miss a discrepancy or encounter a
singularity/branch-sensitive case. The helpers do not exercise generated C/CUDA,
time evolution, boundary conditions, executable output, or convergence.

## Sources

- [nrpy/validate_expressions/validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `precision`, `assert_equal`, `check_zero`
- [nrpy/validate_expressions/validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `convert_free_symbols_set_to_mpf_dict`, `convert_one_expression_to_mpfmpc`, `inject_mpfs_into_cse_expression`
- [nrpy/validate_expressions/validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `process_dictionary_of_expressions`, `output_trusted`, `compare_against_trusted`, `compare_or_generate_trusted_results`

## See Also

- Parent: [Validation](index.md)
- Depends on: [Test Oracles And Safe Updates](test-oracles-and-safe-updates.md)
- See also: [Trusted Expression Pipeline](../equations/trusted-expression-pipeline.md)
- Validated by: [Static Analysis](static-analysis.md)
- See also: [Reference Metrics](../core/reference-metrics.md)
