# Expression Validation Helpers

> High-precision `mpmath` evaluation helpers for SymPy equality checks, zero tests, and `tests/*.py` trusted dictionaries. · Status: confirmed
> Up: [Validation](index.md)

## Summary

`nrpy/validate_expressions/validate_expressions.py` evaluates SymPy expressions in `mpmath` (`precision = 30`, relative tolerance `10 ** (-4.0 / 5.0 * precision)`) for equality assertions (`assert_equal`), zero tests (`check_zero`), and `tests/<basename>.py` `trusted_dict` regression across NRPy's equation and reference-metric modules.

## Detail

### Equation Caller Procedure

Prefer exact analytic or symbolic identities when practical ([Test Oracles And Safe Updates](test-oracles-and-safe-updates.md#expression-oracle-selection)). When sampled numerical regression is appropriate (used across GR, GRHD, GRMHD, wave, elliptic, TOV, SEOBNR/BOB, and geometry modules; `quaternion_rotations/tensor_rotation.py` is a doctest-only exception):

1. Build a dictionary with stable descriptive keys from module attributes or expression names.
2. Call `process_dictionary_of_expressions(..., fixed_mpfs_for_free_symbols=True)` explicitly (`fixed_mpfs_for_free_symbols` defaults to `False`, which draws non-repeatable random samples). The processor sorts items, flattens nested tensor lists to `<name>_<flat_index>`, skips keys containing `funcform`, and **silently overwrites** colliding processed names (e.g., `{"v": [sp.Integer(1)], "v_0": sp.Integer(2)}` keeps only `v_0 = 2`), so callers must choose collision-free names.
3. Call `compare_or_generate_trusted_results(os.path.abspath(__file__), os.getcwd(), basename, results_dict)` to compare against or initialize `tests/<basename>.py`.

Evidence (descriptive behavior): `validate_expressions.py::process_dictionary_of_expressions` and `convert_free_symbols_set_to_mpf_dict`; corroboration: `WaveEquation_RHSs.py::__main__` (`fixed_mpfs_for_free_symbols=True`).

### Equality Assertions And Zero Checks

`assert_equal()` accepts two dictionaries or two SymPy expressions (wrapping non-dict inputs into `{ "": sp.sympify(x) }`), evaluates both with fixed `mpf` substitutions, and checks relative error at every processed key. Type annotations cover scalar SymPy expressions and tensor ranks `1`–`4`; non-SymPy leaves raise `TypeError`. Dictionary inputs require identical raw key sets and list nesting (`AssertionError` on mismatch or flattened-name collisions among numerical entries); keys containing `funcform` still undergo key, nesting, and SymPy-leaf checks while skipping numerical evaluation and collision checks. Matching `NaN` real/imaginary parts and same-signed infinities compare equal (`_nonfinite_values_match`), whereas one-sided `NaN`, finite-vs-infinite, or opposite-sign infinities fail.

Evidence (descriptive behavior): `validate_expressions.py::assert_equal` and `_nonfinite_values_match`; corroboration: `nrpylatex/test_parse_BSSN.py::test_example_BSSN`.

`check_zero()` evaluates one expression through `convert_one_expression_to_mpfmpc()` (including its doubled-precision retry) and checks exact equality with `mp.mpf("0.0")`; callers with free symbols must pass `fixed_mpfs_for_free_symbols=True`.

### Substitution Seeding, CSE Evaluation, And Doubled-Precision Retry

In fixed mode, `convert_free_symbols_set_to_mpf_dict` seeds `random.seed(int(hashlib.md5(str(var).encode()).hexdigest(), 16 + hex_offset))` per symbol (`hex_offset` shifts the integer conversion base, not the MD5 digest), assigning `mp.pi` to `PI`/`M_PI` and `1/sqrt(2)` to `SQRT1_2`/`M_SQRT1_2` (repeatable only within matching Python/SymPy/`mpmath` versions).

`convert_one_expression_to_mpfmpc` returns `mpf(0.0)` immediately for exact SymPy `0`; otherwise `inject_mpfs_into_cse_expression` runs `sp.cse(expr, order="none")`, substitutes `mpf` symbols via `xreplace`, and attempts `mpf` conversion—falling back on `TypeError` to diagnostic `sp.nan` handling, explicit `mpf("+inf")`/`mpf("-inf")`, or component-wise `mpc` construction. If `0 < |result| < 10 ** (-4.0 / 5.0 * mp.dps)`, `convert_one_expression_to_mpfmpc` clears SymPy's process-wide cache, re-evaluates at `2 * precision` (without forwarding `hex_offset`), clears the cache again so later standard-precision calls do not reuse `2 * precision` values, collapses results below the doubled-precision tolerance to `mp.mpf("0.0")`, and restores `mp.dps = precision`.

Evidence (descriptive behavior): `validate_expressions.py::inject_mpfs_into_cse_expression` and `convert_one_expression_to_mpfmpc`; corroboration: `nrpy/tests/reference_metric_GeneralRFM_fisheyeN2.py::trusted_dict` (`Jac_dUrfm_dDSphUD_8 == 0.0`).

### Trusted-File Generation And Comparison

`output_trusted()` writes `tests/<trusted_file_basename>.py` relative to caller path and `os_getcwd`, emitting required `mpmath` imports and `trusted_dict` formatted with Black; missing-file creation initializes a candidate file rather than validating correctness (use the two-process update workflow in [Test Oracles And Safe Updates](test-oracles-and-safe-updates.md)). `compare_against_trusted()` imports `trusted_dict`, checks length equality first (printing missing keys on mismatch), looks up each key in `results_dict` (`KeyError` if renamed at equal length), checks relative error on finite entries, and raises `ValueError` on finite/non-finite or mismatched non-finite components.

Evidence (descriptive behavior): `validate_expressions.py::compare_against_trusted` and `_nonfinite_values_match`; corroboration: `nrpy/tests/reference_metric_GeneralRFM_fisheyeN2.py::trusted_dict` (`NaN` sentinels).

## Sources

- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `precision`, `assert_equal`, `check_zero`, `convert_free_symbols_set_to_mpf_dict`, `convert_one_expression_to_mpfmpc`, `inject_mpfs_into_cse_expression`, `process_dictionary_of_expressions`, `output_trusted`, `compare_against_trusted`, `compare_or_generate_trusted_results`, `_nonfinite_values_match`
- [reference_metric_GeneralRFM_fisheyeN2.py](../../nrpy/tests/reference_metric_GeneralRFM_fisheyeN2.py) - `trusted_dict`
- [test_parse_BSSN.py](../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py) - `test_example_BSSN`
- [WaveEquation_RHSs.py](../../nrpy/equations/wave_equation/WaveEquation_RHSs.py) - `__main__`
- [tensor_rotation.py](../../nrpy/equations/quaternion_rotations/tensor_rotation.py) - `rotate` doctests

## See Also

- Depends on: [Test Oracles And Safe Updates](test-oracles-and-safe-updates.md)
- Depends on: [Code Test Policy](code-test-policy.md)
- Validated by: [Static Analysis](static-analysis.md)
