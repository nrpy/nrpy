# Test Oracles And Safe Updates

> Prospective rules for selecting, storing, reviewing, and updating NRPy `trusted_dict` and `validate_strings` reference files safely. · Status: provisional
> Up: [Validation](index.md)

## Summary

Stored reference files (`trusted_dict` modules and `validate_strings` C/C++/CUDA baselines) become test evidence only when an executable defining-module `__main__` path generates a candidate and compares it against an independently specified mathematical or structural requirement. Prefer exact symbolic or analytic invariants first; create or update reference files only in an isolated worktree copy using two fresh Python processes and full diff review.

## Detail

### `tests/` Store And Reference-File Rules

Outside [Narrow Examples Product Exception](examples-product-exception.md), sibling `nrpy/**/tests/` directories are non-executable reference stores holding only: (1) generated `trusted_dict` modules (`mpf`/`mpc` imports + `trusted_dict` only—never hand-edited, no docstrings, functions, classes, or assertions), (2) generated C/C++/CUDA baselines, (3) empty `__init__.py` markers, and (4) maintainer-approved plain-text scientific provenance cited directly by a production module to document shipped coefficients without commands or assertions (`nrpy/equations/seobnr/tests/BOB_v2_fit_sim_list.md` is nonprecedential legacy text).

### Generation Check Requirements

Generation checks must assert file structure, metadata, or bug-specific output invariants specified independently (via API, docstring, or acceptance criteria) before inspecting candidate output—never completion, exit status, or file existence alone ([meaningful assertion requirements](code-test-policy.md#meaningful-assertion-requirements)).

### Expression Oracle Selection

Prefer exact analytic or symbolic identities (`check_zero(..., fixed_mpfs_for_free_symbols=True)` when free symbols exist; symbol-free expressions need no substitution). When sampled numerical regression is appropriate, the defining module's `__main__` must build a stable-keyed dictionary, evaluate `process_dictionary_of_expressions(..., fixed_mpfs_for_free_symbols=True)`, and compare via `compare_or_generate_trusted_results()` ([Expression Validation Helpers](expression-validation-helpers.md)). Fixed `mpf`/`mpc` sampling gives repeatable high-precision checks under current SymPy/mpmath versions rather than formal proof; meaningful doctest-only equation modules need no artificial `trusted_dict`.

### Emitted Source Baselines, Focused Assertions, And Variant Selection

Use `validate_strings()` (`file_ext` in `{"c", "cpp", "cu"}`, normalized with `clang_format()` from clean registry/parameter state via public `CFunction.full_function` or header/file emitters; untouched legacy baselines stay unformatted) only for small (`~1–80` KB), predominantly handwritten infrastructure control code: initial-data fills, algebraic projection, pointer/loop scaffolding, headers, and parameter files. Never golden-test large SymPy/FD kernels (RHS, Ricci, or constraint evaluations)—validate those symbolically in `nrpy/equations/`.

Once any variant of a `CFunction` family has a `validate_strings()` baseline, add no substring, line, count, or ordering fragment checks on any family variant; cover extra variants with full baselines or compiled integration checks, and reserve focused assertions for metadata, API errors, or properties stripped by formatting. Validate every modified variant branch plus one documented representative per shared-logic equivalence class (`reference_metric.py::unittest_CoordSystems` is a module fact, not a universal template).

Evidence (normative rule): `coding_style.md` (`#### validate_strings pattern`, `### Trusted Vector File Contract`); `reference_metric.py::unittest_CoordSystems`.

### Two-Process Update Procedure, State Cleanup, And Subprocesses

Reference creation or regeneration requires eight steps:
1. Inspect the defining module's `__main__`, `cwd`, cache, registry, and file writes.
2. Work from repository root in an isolated copy containing only intended task edits.
3. Run a first fresh process to emit the candidate (removing only the exact old reference file when absence is required—never directories or unrelated baselines).
4. Review the full diff, `git status`, and underlying equations or templates independently.
5. Run a second fresh process with the reviewed file in place so the final run **compares** rather than generates.
6. Retain new files only in explicit update mode; ordinary compare/doctest runs leave zero new files, restore `cwd`/environment (and registries if same-process callers follow), stay offline, use disposable caches, and never register `CodeParameter` objects at import time (`coding_style.md` `### CodeParameter Registration Scope`).
7. Never hand-edit or bulk-regenerate unrelated reference files.
8. Explain every intentional reference change in the commit message.

Ordinary module checks may call `clang-format`; compilers and runtimes belong in scoped CI. Every subprocess helper requires `shell=False` argument vectors, explicit timeouts, checked exit codes, small fixed job counts, bounded log capture, and `timeout-minutes` in CI YAML, with pinned dependency commits or logged versions.

### Known Gaps And Non-Precedents

`nrpy/helpers/generic.py::clang_format` calls `Popen.communicate()` without a timeout (bound enclosing commands instead), and `BHaH/Makefile_helpers.py::compile_Makefile` uses `shell=True` without timeouts or doctests (do not copy it). Substantive edits to legacy compile doctests must migrate them to CI—or, if out of scope, bound them in place with `TemporaryDirectory`, `shell=False`, fixed `-j`, explicit timeout, checked exit, and `try/finally` `cwd` cleanup without broadening coverage. In the shared `/work` tree, run only side-effect-free checks; run generators, builds, reference updates, and `black .` only in isolated copies.

Evidence (descriptive behavior): `nrpy/helpers/generic.py::clang_format`; `BHaH/Makefile_helpers.py::compile_Makefile`; `coding_style.md` (`### CodeParameter Registration Scope`).

### Contributor Checklist

- [ ] Keep `nrpy/**/tests/` non-executable; never hand-edit `trusted_dict` or source baselines, and distinguish candidate creation from second-process comparison.
- [ ] Prefer exact invariants; use `fixed_mpfs_for_free_symbols=True` for free-symbol `check_zero` / `process_dictionary_of_expressions` checks.
- [ ] Apply `clang_format()` + `validate_strings()` only to small structural C/C++/CUDA emitters—never RHS/Ricci/constraint kernels or substring fragments alongside a family baseline.
- [ ] Cover every modified branch plus one documented representative per equivalence class.
- [ ] Run updates, builds, and `black .` in an isolated copy with two fresh processes, `shell=False` bounded subprocesses, offline deterministic state, and side-effect-free shared-tree checks; follow [Code Test Policy](code-test-policy.md#contributor-checklist).

## Sources

- [coding_style.md](../../coding_style.md) - `### Expression Validation via Trusted Dictionaries`, `### Trusted Vector File Contract`, `### CodeParameter Registration Scope`, `#### validate_strings pattern`
- [generic.py](../../nrpy/helpers/generic.py) - `clang_format`, `validate_strings`
- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `assert_equal`, `check_zero`, `process_dictionary_of_expressions`, `compare_or_generate_trusted_results`
- [WaveEquation_RHSs.py](../../nrpy/equations/wave_equation/WaveEquation_RHSs.py) - `__main__` expression check
- [reference_metric.py](../../nrpy/reference_metric.py) - `unittest_CoordSystems`
- [Makefile_helpers.py](../../nrpy/infrastructures/BHaH/Makefile_helpers.py) - `compile_Makefile`

## See Also

- Depends on: [Code Test Policy](code-test-policy.md)
- Depends on: [Expression Validation Helpers](expression-validation-helpers.md)
- Validated by: [Static Analysis](static-analysis.md)
