# Code Test Policy

> Prospective placement, regression-value criteria, integration boundaries, and proof-scope rules for tests of `nrpy/**` code. · Status: provisional
> Up: [Validation](index.md)

## Summary

New or modified code under `nrpy/**` uses the lowest-cost validation layer that fails on a realistic regression and verifies an explicit mathematical, API, or output requirement. In-process Python, symbolic-expression, and structural C/C++/CUDA checks live in the defining module's docstrings and `__main__` block, whereas compiler, build, and runtime checks belong in CI. When no nontrivial behavioral assertion exists, add no behavioral test while still running static analysis.

## Detail

### Module-Local Placement And Integration Boundary

In `nrpy/**`, place Python API, symbolic-expression, and `validate_strings()` checks in the defining module's docstrings executed via `if __name__ == "__main__":` (`doctest.testmod()` failure runner; infrastructure modules place `Doctests:` immediately before `>>>` prompts and keep test-only imports local). `static-analysis` in `.github/workflows/main.yml` skips `*/tests/*` and executes non-example modules directly. Never add `pytest`/`unittest` dependencies, standalone test runners, executable sibling test modules, `# FIXME` doctest placeholders, or wrapper tests duplicating defining-module checks.

External compilers, build systems, and runtimes belong in CI workflows rather than module doctests (`clang-format` remains permitted for `validate_strings()`). The sole standalone C/C++ exception is `coding_style.md` `### External-Host Test Harnesses` (`tests_infra/`, such as `Dendro/tests_infra/dendrolib_capability_test.cpp`), built and run manually to verify external-host layout assumptions absent from NRPy CI without recording volatile run snapshots; legacy compile doctests follow [Known Gaps And Non-Precedents](test-oracles-and-safe-updates.md#known-gaps-and-non-precedents).

Evidence (normative rule): `coding_style.md` (`### if __name__ == "__main__": Block`, `### Doctest Conventions`, `### External-Host Test Harnesses`); `.github/workflows/main.yml`.

### Meaningful Assertion Requirements

Every behavioral or regression test must fail on a plausible bug and verify a concrete requirement: a return value or state transition, a documented exception, a mathematical/symbolic identity, exact registration metadata (required names, fields, or relations—never bare registration success), small handwritten emitted control code, or an explicit CI numerical tolerance ([Generation Check Requirements](test-oracles-and-safe-updates.md#generation-check-requirements)). Symbolic equations are validated under `nrpy/equations/`, never via kernel-scale infrastructure expression files.

Prohibited low-signal tests: import-only, call-only, bare registry-presence, generation-completed, file-count/existence, generic exit-status, incidental codegen assignment substrings, and generated `__version__` import checks.

### Status-Only Build Or Crash Conditions

A compile, link, or crash test asserting only exit `0` is permitted only when all six conditions hold: (1) a concrete compile/link compatibility requirement or crash bug is documented, (2) target/configuration/input are minimal, (3) the pre-fix revision fails to build or crashes while the fix succeeds identically, (4) compiler/runtime versions are named, (5) no stronger numerical output exists, and (6) claims are restricted to that exact build or crash check. New standalone CLI helpers under `nrpy/examples/**` additionally require all twelve conditions in [Narrow Examples Product Exception](examples-product-exception.md).

### Validation Layers, Legacy Outliers, And Decision Sequence

| Layer | What it establishes (and never establishes later layers) |
| --- | --- |
| Doctest | In-process Python return value, state, exception, or small algebraic identity. |
| Sampled expression | High-precision numerical agreement at fixed `mpf`/`mpc` points (not formal proof). |
| Emitted source (`validate_strings`) | Exact `clang-format` text match on small structural code (neither compiles nor runs). |
| Generation | Independently specified file layout, metadata, or output invariant. |
| Build / Runtime crash | Toolchain compile/link or non-crashing execution under the 6-condition check above. |
| Regression / Convergence | Output matches reference within tolerance, or multi-resolution error scales at expected order. |

Decision order: (1) Python/API -> module doctest; (2) symbolic equation -> exact identity or `trusted_dict` in `nrpy/equations/`; (3) small structural emitter -> `validate_strings()`; (4) large FD kernel -> upstream symbolic check; (5) compiler/toolchain -> CI; (6) standalone example numerical test -> 12-condition [examples exception](examples-product-exception.md); (7) no nontrivial invariant -> no behavioral test; (8) edited handwritten Python -> [Static Analysis](static-analysis.md).

Retained nonprecedential outliers: `nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py` (standalone sampled cross-representation script) and `JAX/jax_project_generator.py`'s emitted downstream `pytest` `__version__` file (unrun in CI).

Evidence (descriptive behavior): `nrpylatex/test_parse_BSSN.py::test_example_BSSN`; `JAX/jax_project_generator.py::output_PyFunction_files_and_construct_project`; `.github/workflows/main.yml`.

### Contributor Checklist

- [ ] Name the defining module, mathematical/API requirement, validation layer, and pre-fix failure mode before writing a test; add no test when no nontrivial invariant exists.
- [ ] Keep ordinary checks inside defining-module docstrings/`__main__`; add no `pytest`/`unittest`, `# FIXME` placeholders, compile doctests, or duplicate wrapper tests.
- [ ] Satisfy all 6 conditions above for any status-only build/crash check and all 12 conditions in [Narrow Examples Product Exception](examples-product-exception.md) for any `nrpy/examples/**` test script.
- [ ] Run [Static Analysis](static-analysis.md) on every edited handwritten Python file and follow the [reference-update checklist](test-oracles-and-safe-updates.md#contributor-checklist).

## Sources

- [coding_style.md](../../coding_style.md) - `### if __name__ == "__main__": Block`, `### Doctest Conventions`, `### External-Host Test Harnesses`, `## Static Analysis Configuration`
- [README.md](../../nrpy/infrastructures/Dendro/tests_infra/README.md) - external-host `dendrolib_capability_test.cpp` scope
- [main.yml](../../.github/workflows/main.yml) - `static-analysis`, `codegen-ubuntu`, `codegen-mac`
- [generic.py](../../nrpy/helpers/generic.py) - `validate_strings`
- [test_parse_BSSN.py](../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py) - `test_example_BSSN`
- [jax_project_generator.py](../../nrpy/infrastructures/JAX/jax_project_generator.py) - `_generate_project_metadata`, `output_PyFunction_files_and_construct_project`

## See Also

- Depends on: [Test Oracles And Safe Updates](test-oracles-and-safe-updates.md)
- Validated by: [Static Analysis](static-analysis.md)
- See also: [Narrow Examples Product Exception](examples-product-exception.md)
