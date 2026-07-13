# Code Test Policy

> Prospective placement, meaningfulness, and proof rules for tests of NRPy code. · Status: provisional · Last reconciled: 07-13-2026
> Up: [Validation](index.md)

## Summary

New or modified code under `nrpy/**` uses the cheapest validation layer that
can fail on a plausible production regression and prove a durable contract.
Ordinary Python, symbolic-expression, and emitted-source checks stay with their
owning production module; compiler, build, runtime, and numerical-product
integration belongs in scoped CI. No meaningful contract means no behavioral
test, while required static quality gates still apply.

This page sets prospective policy. Its provisional status reflects incomplete
durable support for the commissioned norms; it does not claim that every legacy
test shape or current automation path conforms.

## Detail

### Scope And Authority

This policy governs new and modified code under `nrpy/**`. **Must**, **must
not**, and **only** state requirements. **Should** states the default and needs
a concrete documented reason to override. Explicit requirements for a task
control that task, followed by current coding-style and routed KB governance,
executable local or CI configuration, and repeated source practice. A legacy
outlier is evidence about current repository contents, not precedent for new
work.

### Owner-Local Architecture

For ordinary in-process Python API behavior, symbolic-expression checks, and
emitted-source checks, put the smallest meaningful assertion in the production
module that owns the contract. Use its function, class, or module docstring and
run it through the owner's `__main__` validation path. Runnable equation
modules use the canonical `doctest.testmod()` failure runner; runnable
infrastructure modules should use that shape and place the final `Doctests:`
label immediately before prompts. Keep test-only imports inside the doctest or
`__main__` block.

The configured `static-analysis` job excludes every `*/tests/*` path and
directly executes selected non-example source modules; this describes job
configuration, not a latest successful run.

Claim evidence:
- Claim: The configured `static-analysis` job excludes every `*/tests/*` path and directly executes selected non-example source modules; this describes job configuration, not a latest successful run.
- Role: CI behavior
- Deciding authority: [main.yml](../../.github/workflows/main.yml), `static-analysis`
- Corroboration: [single_file_static_analysis.sh](../../.github/single_file_static_analysis.sh), `run_test_step`, corroborates direct per-file execution but not CI exclusions

New core code must not add pytest or unittest dependencies, configuration,
collectors, fixtures, test classes, executable sibling test modules, standalone
runners, test-only abstractions, or wrapper tests that duplicate the owner
doctest. A runner is useful only when the module contains prompts or performs
meaningful subsequent owner validation, such as an established trusted-
expression comparison. A new runner with neither is prohibited. Zero attempted
doctests is not coverage.

### Meaningful Contract Gate

Every behavioral or regression test must fail for a plausible production
regression and assert a durable contract. Valid contracts include:

- returned value or documented state transition;
- documented error type or condition;
- nontrivial mathematical or semantic property;
- exact documented registration metadata, such as required names, fields,
  values, or relationships, whose mismatch would violate a durable contract;
  bare registration success or presence of some registry entry is not a
  meaningful contract;
- stable, predominantly handwritten emitted control logic; or
- explicit build, runtime, or numerical criterion in scoped integration CI.

A generation check must assert independently documented product structure,
schema, content, metadata, meaningful artifact content, or a bug-specific
invariant. That contract may predate the change or be introduced independently
in the same reviewed API, docstring, specification, or acceptance criteria.
Candidate or golden output cannot define its own contract; see [Test Oracles And
Safe Updates](test-oracles-and-safe-updates.md#generation-contract).

Prohibited coverage includes import-only, call-only, bare registration-success
or unqualified registry-presence-only, generation-completed-only, generic exit-
status, file-existence, and file-count assertions; tests added only because a
function exists; empty runners and placeholder doctests; duplicate owner-plus-
wrapper coverage; lint, discovery, fixture, or baseline plumbing tested without
a production contract; incidental SymPy/codegen assignment substrings; and
downstream smoke scaffolds that only import a package or inspect `__version__`.
A narrow status-only compile, link, or crash regression is allowed only under
the next gate. If no cheap meaningful assertion exists, add no behavioral test.
Static quality gates remain mandatory and separate.

### Integration Boundary

External compiler, build, runtime, and numerical-product validation belongs in
scoped CI, not ordinary owner doctests. Repository-configured `clang-format`
remains an allowed dependency for owner emitted-source normalization. New core
standalone harnesses and compile doctests are prohibited. A substantive change
to an existing compile doctest must follow the migration or tightly bounded
fallback in [Test Oracles And Safe Updates](test-oracles-and-safe-updates.md#known-gaps-and-non-precedents).

### Status-Only Build Or Crash Gate

Generic status-only coverage is prohibited. A compile, link, or crash
regression whose only assertion is success is allowed only when every condition
holds:

1. A production compile/link compatibility contract or crash bug is documented.
2. Target, configuration, and input are minimal and exact.
3. The old revision demonstrably fails to compile/link or crashes, while the
   new revision succeeds under the same setup.
4. Compiler, build, or runtime toolchain and resolved version are named.
5. No stronger semantic result is available from the run.
6. The claim is limited to the named compile/link compatibility or crash
   regression.

### Validation Layers And Claims

Use the cheapest layer that proves the contract and name it accurately:

| Layer | What it may establish |
| --- | --- |
| Doctest | Python/API value, state, error, or small semantic invariant. |
| Sampled expression | Numerical agreement at named substitutions; not formal identity. |
| Trusted source | Exact normalized emitted-text agreement with a reviewed oracle; not compilation or execution. |
| Generation | Independently documented output structure, schema, content, metadata, or bug-specific artifact invariant. Completion alone does not pass. |
| Build | Compile/link compatibility for a named generated product and toolchain. |
| Runtime compatibility/crash | Start or completion for a named input only when the six-condition status gate applies; otherwise runtime needs concrete output, state, or semantic assertion. |
| Regression/numerical | Parsed result satisfies an explicit reference and tolerance. |
| Scientific/convergence | Dedicated evidence designed for the named scientific claim. |

An earlier layer proves nothing about a later one. In particular,
`validate_strings()` neither compiles nor executes emitted code. A generated-
project matrix can establish only named generation and named toolchain
compile/link compatibility unless it runs and checks a stronger result. CI
configuration proves configured commands and assertions, not a latest pass. A
named successful run proves only the commands, assertions, inputs, platforms,
versions, backends, and options it actually exercised; it does not establish
omitted backends, general convergence, or broad scientific correctness.

### Narrow Examples Product Exception

Outside this exception, examples are not ordinary policy evidence. A new
command-line helper under `nrpy/examples/**` is allowed only when all twelve
conditions hold:

1. The target is an NRPy-generated standalone product. Invoking an example
   generator may be setup only; generator argument behavior belongs in an owner
   API/doctest unless inseparable from the product-process contract.
2. An external compiler, build/runtime toolchain, or launcher is used only when
   the generated-product contract requires it.
3. The contract is observable only at the product process boundary: product
   arguments, parameter file, runtime state, or parsed output.
4. No cheaper owner API, symbolic, source-oracle, or semantic check proves the
   same regression.
5. The same change wires an explicit CI invocation and target build/run.
6. Python helper subprocesses use `shell=False` argument vectors. Any
   unavoidable fixed shell orchestration stays in audited CI YAML and never
   interpolates helper or user data.
7. Every helper subprocess has an explicit timeout, small fixed resource or job
   limits, checked failure, and bounded or capped captured output. The CI job
   has `timeout-minutes` or equivalent.
8. Cwd and environment restoration and output cleanup are unconditional.
9. Inputs and randomness are deterministic; generated output is CI-owned or
   temporary.
10. Network, installation, and external access are explicitly authorized.
    Controllable dependencies use an exact pin, tag, or commit; otherwise
    record resolved version and reason and label evidence non-reproducible
    beyond the named run.
11. The helper asserts concrete output, state, semantic property, or numerical
    criterion beyond status, file existence, file count, or import success. A
    status-only product build/crash assertion must also satisfy all six status
    gates and keep the narrow compatibility/crash claim.
12. The assertion names the exact validation layer proved.

If any condition fails, the helper is prohibited. This exception never
authorizes general pytest or unittest adoption in examples.

The checked-in SEOB helpers are finite-input numerical-and-argument precedents:
they use deterministic samples, argument-vector subprocesses with checked
failure, parsed stdout, and a numerical failure criterion, and the workflow
invokes them explicitly. They lack explicit subprocess timeouts, bounded or
capped output capture, unconditional cwd restoration and cleanup, and parts of
the current resource-limit policy; their CI jobs also lack `timeout-minutes`.
None of those gaps may be copied.

Claim evidence:
- Claim: The checked-in SEOB helpers use deterministic samples, argument-vector
  subprocesses with checked failure, parsed stdout, and a numerical failure
  criterion, and the workflow invokes them explicitly; they lack explicit
  subprocess timeouts, bounded or capped output capture, unconditional cwd
  restoration and cleanup, and parts of the current resource-limit policy, and
  their CI jobs lack `timeout-minutes`.
- Role: descriptive behavior
- Deciding authority: [sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py),
  `run_sebob`, `process_input_set`;
  [sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py),
  `run_sebobv2`, `process_input_set`; [main.yml](../../.github/workflows/main.yml),
  `sebob-consistency-test` and `sebobv2-consistency-test`
- Corroboration: none available; helper and workflow configuration jointly decide the distributed claim
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-run; tool_version=not-run; backend=not-run;
  precision=not-run; GPU=not-applicable; restart=not-applicable;
  distributed=not-applicable; error_path=not-run; options=not-run;
  date=07-13-2026`

### Static Summary

Run exact governing command `black .` from repository root in an isolated,
user-owned worktree or copy containing only intended task changes and no
unrelated modifications, then inspect its diff. Run
`./.github/single_file_static_analysis.sh <path.py>` for each modified
handwritten file. Inspect the wrapper's command construction, the target's
direct-execution effects, and every invoked tool's cache and filesystem output
effects. The current wrapper dispatches interpolated command strings through
`eval`. Before invocation, require a repository-relative argument whose resolved
target stays in the repository, is a regular non-symlink file, contains no `..`
component, does not begin with `-`, and matches `^[A-Za-z0-9_./-]+$`. A failing
path blocks invocation until separately authorized wrapper hardening. Run an
accepted path only in the isolated intended-change tree or copy; disable or
redirect every writable tool cache and filesystem output to an owned disposable
location, then inspect repository status. Each file must pass the full gate and
report Pylint **10.00/10.00**. Generated trusted `*/tests/*.py` data is
exempt from per-file analysis; its owning handwritten module is not. Current
automation must not be assumed to guarantee this policy. See [Static
Analysis](static-analysis.md) for command mechanics and enforcement details.

### Decision Tree And Legacy Shapes

1. Ordinary Python/API behavior or error: use an owner doctest.
2. Symbolic property: prefer an exact invariant; otherwise use a fixed sampled
   owner comparison with explicit dictionary-structure checks.
3. Stable handwritten emitted function, header, or file: normalize complete
   public-generator output and compare it with `validate_strings()`.
4. Large generated kernel: validate expressions or semantics upstream; do not
   use a full-text golden or incidental substrings.
5. Required build/toolchain: use scoped CI; add no compile doctest.
6. Runtime or numerical standalone generated example product: use the examples
   exception only when all twelve gates apply; status-only cases also need all
   six status conditions.
7. No meaningful contract: add no behavioral test.
8. Handwritten Python changed: use the isolated exact formatter and full per-
   file static gate, including Pylint 10.00.

Three checked-in shapes are retained facts, not templates. This policy treats
the NRPyLaTeX BSSN file as a unique retained direct-execution sampled cross-
representation harness. The JAX
project generator emits a downstream pytest import/`__version__` scaffold that
configured generated-project CI does not execute; it is low-signal and its
removal or replacement needs a separate JAX decision. BHaH `compile_Makefile()`
contains an external-compilation doctest that is unsafe under this policy.

Claim evidence:
- Claim: The NRPyLaTeX BSSN file directly executes a sampled cross-representation comparison; the JAX project generator emits a downstream pytest import/`__version__` scaffold that configured generated-project CI does not execute; and BHaH `compile_Makefile()` contains an external-compilation doctest.
- Role: descriptive behavior
- Deciding authority: [test_parse_BSSN.py](../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py), `test_example_BSSN`; [jax_project_generator.py](../../nrpy/infrastructures/JAX/jax_project_generator.py), `_generate_project_metadata` and `output_PyFunction_files_and_construct_project`; [Makefile_helpers.py](../../nrpy/infrastructures/BHaH/Makefile_helpers.py), `compile_Makefile`
- Corroboration: [main.yml](../../.github/workflows/main.yml), `codegen-ubuntu` and `codegen-mac`, corroborates that JAX output is generated without executing its scaffold; no independent corroboration for the other two outliers
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-run; tool_version=not-run; backend=not-run; precision=not-run; GPU=not-run; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=not-run; date=07-13-2026`

Do not expand any of those shapes. Whether an unrelated touch must remove a
legacy empty runner remains maintainer judgment. Meaningful doctest-only
equation modules remain valid exceptions to broad historical wording about
trusted dictionaries; meaningful validation wins over adding a golden merely
for test count.

### Contributor Checklist

- [ ] Identify owner, durable contract, validation layer, and old-behavior
      failure.
- [ ] Use owner doctest for ordinary in-process behavior; count meaningful
      subsequent owner validation separately from attempted prompts.
- [ ] Add no core standalone runner, pytest/unittest surface, test-only
      abstraction, or duplicate wrapper.
- [ ] Add no test when no meaningful behavioral contract can be asserted.
- [ ] Keep compiler/build/runtime integration in scoped CI; add no compile
      doctest.
- [ ] Apply all six status-only conditions before making a narrow compile/link
      or crash claim.
- [ ] Apply all twelve examples-product gates, including same-change CI,
      concrete result, subprocess bounds, access authority, and exact claim.
- [ ] State only what configuration or a named successful run establishes.
- [ ] For every handwritten Python change, use isolated `black .`, inspect its
      diff, pass the full per-file gate, and report Pylint 10.00.
- [ ] Follow the separate [oracle checklist](test-oracles-and-safe-updates.md#contributor-checklist)
      for stored evidence, variants, updates, and mutable state.

## Sources

- [coding_style.md](../../coding_style.md) - `### if __name__ == "__main__": Block`, `### Doctest Conventions`, `#### Doctest placeholders`, `## Static Analysis Configuration`
- [main.yml](../../.github/workflows/main.yml) - `static-analysis`, `codegen-ubuntu`, `codegen-mac`, `sebob-consistency-test`, `sebobv2-consistency-test`
- [single_file_static_analysis.sh](../../.github/single_file_static_analysis.sh) - `run_test_step`
- [.pylintrc](../../.pylintrc) - `[MASTER]`; [.pylintrc_python36](../../.pylintrc_python36) - `[MASTER]`
- [generic.py](../../nrpy/helpers/generic.py) - `validate_strings`
- [test_parse_BSSN.py](../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py) - `test_example_BSSN`
- [jax_project_generator.py](../../nrpy/infrastructures/JAX/jax_project_generator.py) - `_generate_project_metadata`, `output_PyFunction_files_and_construct_project`
- [Makefile_helpers.py](../../nrpy/infrastructures/BHaH/Makefile_helpers.py) - `compile_Makefile`
- [sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py) - `run_sebob`, `process_input_set`
- [sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py) - `run_sebobv2`, `process_input_set`

## See Also

- Parent: [Validation](index.md)
- Depends on: [Test Oracles And Safe Updates](test-oracles-and-safe-updates.md)
- Validated by: [Static Analysis](static-analysis.md)
- See also: [Expression Validation Helpers](expression-validation-helpers.md)
- See also: [Generated Project CI](generated-project-ci.md)
