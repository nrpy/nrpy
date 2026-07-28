# Test Oracles And Safe Updates

> Prospective rules for selecting, storing, reviewing, and safely updating NRPy test oracles. · Status: provisional · Last reconciled: 07-20-2026
> Up: [Validation](index.md)

## Summary

A test oracle is reviewed expected data, not executable coverage. The owning
production validation path must generate or observe a candidate and compare it
against an independently defined contract. Choose exact analytic or semantic
checks first; use sampled expression or normalized emitted-source oracles only
when their proof limits fit the contract.

Oracle creation and regeneration happen in an isolated user-owned intended-
change tree, use two fresh processes, and end with a comparison against the
reviewed candidate. This page sets prospective policy. Its provisional status
does not claim that legacy stores, baselines, helpers, or subprocesses already
conform.

## Detail

### Scope And Authority

This page owns prospective policy for oracle stores, oracle selection, variant
coverage, update review, mutable state, dependencies, subprocesses, and
workspace safety. [Code Test Policy](code-test-policy.md) owns placement,
meaningful-contract gates, validation-layer names, integration boundaries, and
the examples exception. Repository behavior corroborates mechanics and known
gaps; retained exceptions do not authorize new ones.

### Oracle Versus Executable Test

An oracle is stored data or an expected artifact. It becomes evidence only when
an executable owner path generates or observes a candidate and performs a
meaningful comparison or assertion. Oracle presence alone is not a test.

Missing-expected-file creation is candidate generation, not a passing
comparison, and cannot establish that the candidate is true. [Expression
Validation Helpers](expression-validation-helpers.md) owns trusted-expression
helper behavior; [Maintenance And Validation
Helpers](../core/helpers/maintenance-and-validation-helpers.md) owns emitted-
source helper behavior.

### Ordinary `tests/` Store Contract

Outside the narrow generated-product examples exception, sibling
`nrpy/**/tests/` directories are non-executable oracle and support-evidence
stores, not runner roots. They may contain only:

- generated trusted numerical dictionaries;
- generated exact C or CUDA source baselines;
- required empty package markers; and
- narrowly admitted scientific provenance or support text.

Every new support category or file must satisfy all six admissions:

1. It is plain, non-executable text.
2. Its owning production module directly cites or consumes it.
3. Its source or provenance is necessary to interpret a shipped coefficient or
   oracle.
4. Its colocation reason is documented.
5. It contains no commands, harness logic, executable snippets, or test
   assertions.
6. A maintainer explicitly decides to admit that category and file.

New handwritten executable logic, runners, assertions, fixtures, and framework
suites are prohibited in ordinary stores. A trusted Python file contains only
needed `mpf` or `mpc` imports plus `trusted_dict`: no module docstring,
function, class, runner, fixture, or assertion. Never hand-edit trusted values.
Directory placement never exempts legacy or exceptional handwritten Python
from static analysis, but that rule does not admit new executable Python.

Retained `nrpy/equations/seobnr/tests/BOB_v2_fit_sim_list.md` scientific
provenance is legacy support, not generated or validated test output and not
precedent for another support file. The retained classification is prospective
governance; that unregistered support file is not cited as evidence here.

### Generation Contract

Generation validation must assert product structure, schema, content, metadata,
meaningful artifact contents, or a bug-specific invariant defined before the
candidate is inspected. The contract may pre-exist or be introduced
independently in the same reviewed feature change as an API, docstring,
specification, or explicit acceptance criteria. Candidate and golden output
cannot author their own contract. Completion, status, existence, and file count
alone do not pass the [meaningful-contract gate](code-test-policy.md#meaningful-contract-gate).

### Expression Oracle Selection

Prefer an exact analytic, symbolic, or semantic invariant when practical. When
sampled trusted-result regression is appropriate, the owner path must:

1. build a dictionary with stable descriptive keys;
2. process trusted results with `fixed_mpfs_for_free_symbols=True`; and
3. compare through `compare_or_generate_trusted_results()` in the established
   owner path.

Direct dictionary `assert_equal()` calls must use stable keys. [Expression
Validation Helpers](expression-validation-helpers.md) owns helper mechanics,
structure checks, comparison behavior, and limitations.

Regression uses of `check_zero()` with free symbols must pass
`fixed_mpfs_for_free_symbols=True`; truly symbol-free expressions do not need a
substitution. Fixed sampling means repeatable sampled high-precision evidence
under documented current runtime assumptions when free symbols exist. It is
not formal identity or guaranteed bitwise stability across Python, SymPy,
mpmath, platform, or codegen versions. A symbol-free result is a numerical
evaluation, still not a formal proof. Add analytic, branch, singularity, zero,
and boundary cases when the contract needs them. Do not create golden data to
satisfy a test count; meaningful doctest-only equation modules are allowed.

### Emitted Source Oracle Selection

Use an exact emitted-source oracle only when full text is a useful stable
contract, normally predominantly handwritten infrastructure control logic.
Supported candidates are complete `CFunction.full_function` text retrieved
after public registration, or complete output from a public header or file
generator.

The owner must establish clean explicit registry and parameter state, generate
through the public path, capture complete output, normalize every new or
regenerated full text with `clang_format()`, and call `validate_strings()` with
a stable whitespace-free description and correct `file_ext="c"` or
`file_ext="cu"`. Untouched raw baselines are not bulk-normalized.

Do not golden-test a large SymPy/codegen-dominated kernel. Validate upstream
expressions, dimensions, key sets, semantic invariants, or a small stable
wrapper instead. Incidental generated assignments are not contracts, and exact
source comparison neither compiles nor executes output.

### Focused Assertions

A focused assertion may replace an unnecessary full golden with a narrower
durable semantic contract, or test metadata, API behavior, runtime behavior, or
information that full-source equality omits or formatting normalizes away.
Once an owner has a complete emitted-source oracle for any member of a
`CFunction` family, that owner must not also inspect generated-text fragments
for any family variant through substring membership, counts, line matching, or
position/order checks. Such fragment checks duplicate the oracle mechanism and
must not substitute for variant coverage. Cover a variant gap with a complete
variant oracle or a semantic compile/runtime test. Incidental assignment checks
remain prohibited.

### Variant Selection

Validate every modified variant-specific branch plus one representative from
each documented equivalence class that exercises shared logic. If no defensible
partition exists, validate all variants. Record every omitted variant, its
class, the chosen representative, and the shared-code-path rationale.

`reference_metric.py` contains an existing `unittest_CoordSystems` subset.

Claim evidence:
- Claim: `reference_metric.py` contains an existing `unittest_CoordSystems` subset.
- Role: descriptive behavior
- Deciding authority: [reference_metric.py](../../nrpy/reference_metric.py), `unittest_CoordSystems`
- Corroboration: none available; the subset is local to its owning module
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-run; tool_version=not-run; backend=not-applicable; precision=not-run; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=not-run; date=07-13-2026`

That subset is a repository fact, not universal precedent for selecting
variants.

### Two-Process Oracle Update

Oracle creation or update must follow all eight steps:

1. Inspect the owner `__main__` path, generator, cwd, cache, registry, and
   output effects before execution.
2. Work from repository root in an isolated user-owned worktree or copy that
   contains only intended task changes and no unrelated modifications. Preserve
   the old tracked artifact through Git and, when useful, a read-only comparison
   copy.
3. Run the owner in a first fresh process to produce the candidate. Remove only
   the exact old artifact if the helper requires absence; never delete a
   directory, support data, shared output, or unrelated baseline.
4. Inspect the complete candidate, repository status, and full old/new diff.
   Review the mathematics, algorithm, or template independently; a candidate
   cannot attest itself.
5. Accept the candidate explicitly, then run the owner in a second fresh
   process with the reviewed candidate present. The final observed path must
   compare, not generate.
6. Only explicit oracle-update mode may retain the exact candidate through
   review and the second run. Clean every other incidental filesystem effect;
   ordinary comparison or doctest mode retains nothing new.
7. Never hand-edit oracle data to force a match or bulk-regenerate unrelated
   output.
8. Explain every intentional oracle change in the commit message.

### Retention, Cleanup, And State

Ordinary compare or doctest mode preserves no new output and restores every
filesystem effect. Explicit oracle-update mode alone may retain the exact
reviewed candidate; all other effects are cleaned unconditionally.

Run ordinary owner checks from documented repository-root cwd and prefer a
fresh process matching CI direct execution. Initialize relevant registries and
parameters explicitly. Restore cwd and environment unconditionally with a
context manager or `finally`, and restore mutable registries when later same-
process checks can observe mutation; a fresh process with explicit
initialization is sufficient otherwise. Use an owned disposable cache when
needed, make results independent of pre-existing cache contents, and never
clear or overwrite ambient or shared cache. Seed or fix unavoidable samples,
sort unordered emission inputs when order affects the contract, remove owned
temporary outputs, and keep ordinary checks offline.

New or modified code must not register `CodeParameter` objects at import time,
because import-time registration mutates shared generated-parameter state.

Claim evidence:
- Claim: New or modified code must not register `CodeParameter` objects at import time, because import-time registration mutates shared generated-parameter state.
- Role: normative rule
- Deciding authority: [coding_style.md](../../coding_style.md), `### CodeParameter Registration Scope`
- Corroboration: none available; current coding-style governance is the owning authority

### Dependencies And Subprocesses

Ordinary owner checks may depend on repository-configured `clang-format` for
source normalization; it is not compiler/build/runtime integration. Compiler,
build, and runtime toolchains belong only in scoped integration with authority
to install or access them and CI-owned output. Controllable dependencies use an
exact pin, tag, or commit. For an uncontrollable dependency, record the resolved
version and reason and label evidence non-reproducible beyond the named run.

Every new or modified subprocess wrapper and every examples or integration
helper subprocess needs an explicit timeout, checked failure, small fixed
resource or job bounds, and bounded captured output or logging. Its CI job needs
`timeout-minutes` or equivalent. Python always uses `shell=False` argument
vectors; unavoidable fixed shell orchestration stays in audited CI YAML and
never interpolates helper or user data. Document cleanup and every network,
installation, or external-access need.

### Known Gaps And Non-Precedents

Current `clang_format()` calls `Popen.communicate()` without a timeout.

Claim evidence:
- Claim: Current `clang_format()` calls `Popen.communicate()` without a timeout.
- Role: descriptive behavior
- Deciding authority: [generic.py](../../nrpy/helpers/generic.py), `clang_format`
- Corroboration: none available; the subprocess implementation is local to the helper
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-run; tool_version=not-run; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=not-run; date=07-13-2026`

Ordinary callers remain allowed, but their enclosing validation command or job
must be bounded until the helper is hardened.

Existing raw full-text baselines are historical facts. Apply normalization only
when a baseline is regenerated through a scoped reviewed change; bulk migration
is separate work. Retained BOB support is a non-precedent under the store
admission rule.

BHaH `compile_Makefile()` currently runs an external compilation doctest through
`shell=True`, broad host-derived parallelism, and no explicit timeout.

Claim evidence:
- Claim: BHaH `compile_Makefile()` currently runs an external compilation doctest through `shell=True`, broad host-derived parallelism, and no explicit timeout.
- Role: descriptive behavior
- Deciding authority: [Makefile_helpers.py](../../nrpy/infrastructures/BHaH/Makefile_helpers.py), `compile_Makefile`
- Corroboration: none available; the compile doctest and subprocess implementation share one owner
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-run; tool_version=not-run; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=not-run; date=07-13-2026`

It is an unsafe legacy outlier, not a bounded template.

A substantive change to that compile doctest must migrate it to scoped CI. Only
when migration is outside the authorized task scope may the existing surface be
hardened in place with `TemporaryDirectory`, `shell=False` argument vectors,
small fixed parallelism, explicit timeout, checked failure, cwd restoration,
and `try/finally` cleanup. The fallback must not broaden the tested matrix,
generated project, toolchain work, or subprocess surface.

### Workspace Safety

In the shared repository working tree, run only demonstrably side-effect-free
scoped checks. Run generators, builds, oracle creation, and `black .` only in an
isolated user-owned intended-change worktree or copy containing no unrelated
modifications. No coordination exception permits mutation in the shared working
tree; never reset or delete shared generated output.

### Contributor Checklist

- [ ] Keep ordinary `tests/` stores non-executable and apply all six admissions
      to new support data.
- [ ] Distinguish stored oracle from the executable owner path; missing-file
      creation is candidate generation, not comparison.
- [ ] Prefer exact expression invariants; otherwise check dictionary lengths
      and keys and use fixed substitutions where free symbols exist.
- [ ] Normalize complete stable public-generator output; do not golden-test a
      large kernel or incidental substring.
- [ ] For a `CFunction` family with any complete emitted-source oracle, add no
      generated-fragment membership, count, line, or position/order assertions;
      cover variant gaps with a complete variant oracle or semantic
      compile/runtime test.
- [ ] Use focused assertions only to replace an unnecessary golden, test
      metadata, API, or runtime behavior, or detect information equality or
      formatting omits.
- [ ] Cover modified variant branches and each documented class, or all
      variants; record every omission and rationale.
- [ ] Use an isolated owned intended-change tree and two fresh processes;
      inspect the full candidate, status, and diff before the final comparison.
- [ ] Delete only the exact old artifact when required; never delete a
      directory, support/shared data, or unrelated baseline, and never hand-
      edit or bulk-regenerate oracles.
- [ ] Restore cwd, environment, and filesystem effects unconditionally; restore
      registries when same-process observers remain; use owned cache and temp
      output and keep ordinary checks offline.
- [ ] Bound subprocess time, resources, failures, and output; use `shell=False`
      argument vectors, exact pins or qualified resolved versions, and explicit
      access authority.
- [ ] In the shared repository working tree, run only side-effect-free scoped
      checks.
- [ ] Follow [Code Test Policy](code-test-policy.md#contributor-checklist) for
      placement, meaningfulness, integration, examples, claims, and static
      gates.

## Sources

- [coding_style.md](../../coding_style.md) - `### Expression Validation via Trusted Dictionaries`, `### Trusted Vector File Contract`, `### CodeParameter Registration Scope`, `#### validate_strings pattern`
- [generic.py](../../nrpy/helpers/generic.py) - `clang_format`, `validate_strings`
- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `assert_equal`, `check_zero`, `process_dictionary_of_expressions`, `compare_or_generate_trusted_results`
- [WaveEquation_RHSs.py](../../nrpy/equations/wave_equation/WaveEquation_RHSs.py) - `WaveEquation_RHSs`, module `__main__` path
- [reference_metric.py](../../nrpy/reference_metric.py) - `unittest_CoordSystems`
- [Makefile_helpers.py](../../nrpy/infrastructures/BHaH/Makefile_helpers.py) - `compile_Makefile`

## See Also

- Parent: [Validation](index.md)
- Depends on: [Code Test Policy](code-test-policy.md)
- Depends on: [Expression Validation Helpers](expression-validation-helpers.md)
- See also: [Trusted Expression Pipeline](../equations/trusted-expression-pipeline.md)
- See also: [Infrastructure Code Style](../infrastructures/infrastructure-code-style.md)
- See also: [Maintenance And Validation Helpers](../core/helpers/maintenance-and-validation-helpers.md)
- Validated by: [Static Analysis](static-analysis.md)
- See also: [Workflows](../workflows.md)
