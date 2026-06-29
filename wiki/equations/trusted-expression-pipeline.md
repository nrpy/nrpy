# Trusted Expression Pipeline

> Explain how symbolic equation outputs become trusted numerical validation files. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Equations](index.md)

## Summary

Equation modules validate symbolic expressions by building a results dictionary,
converting it to deterministic high-precision numerical values, and comparing
or generating a sibling trusted file under `tests/`. Trusted files are generated
evidence, not hand-authored code, and ordinary per-file static analysis is
skipped for them.

## Detail

The common equation-module flow is:

1. Build a dictionary from object state or explicit expression names.
2. Call `process_dictionary_of_expressions(...)`, usually with
   `fixed_mpfs_for_free_symbols=True`.
3. Call `compare_or_generate_trusted_results(...)` with the owning module path,
   working directory, trusted-file basename, and processed results.

`process_dictionary_of_expressions` sorts dictionary items, ignores keys
containing `funcform`, flattens tensor lists, and converts every SymPy expression
to an `mpf` or `mpc` result. `compare_or_generate_trusted_results` derives the
owning module's `tests/<trusted_file_basename>.py` path. If the file exists, it
loads `trusted_dict` and compares expression count and values; if it does not
exist, it writes a new trusted dictionary.

`output_trusted` writes only the needed `mpmath` imports plus `trusted_dict`, and
formats the file with Black. `compare_against_trusted` raises on missing keys or
value mismatch and tells maintainers to delete the stale trusted file and rerun
the owning module only when the new result is trusted.

Trusted-value files under `*/tests/*.py` are treated specially. They should
contain only generated `mpf` or `mpc` data, no module docstrings, no functions,
and no classes. The preserved agent rules and coding style both say not to
hand-edit trusted values; regenerate them from the owning module and explain the
reason in the commit message.

## Sources

- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `process_dictionary_of_expressions`, `compare_or_generate_trusted_results`
- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `output_trusted`, `compare_against_trusted`
- [coding_style.md](../../coding_style.md) - `Expression Validation via Trusted Dictionaries`, `Trusted Vector File Contract`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `Required Checks`, `Expression Validation`
- [BSSN_RHSs.py](../../nrpy/equations/general_relativity/BSSN_RHSs.py) - `BSSNRHSs`, `BSSN_RHSs`
- [BSSN_RHSs_Cartesian.py](../../nrpy/equations/general_relativity/tests/BSSN_RHSs_Cartesian.py) - `trusted_dict`

## See Also

- [Equations](index.md)
- [BSSN Family](bssn-family.md)
