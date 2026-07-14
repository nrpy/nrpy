# Workflows

> Procedures for ingesting, querying, and maintaining the KB. · Status: confirmed · Last reconciled: 07-13-2026

## Summary

Start from `AGENTS.md`, follow catalog and router paths to compiled leaves, and
answer from sourced pages. Durable answers are filed back into the wiki when
they should compound. Maintenance work registers sources first, updates the
owning leaf, then fixes nearby links, catalog entries, source-map rows, and
glossary terms.
NRPy code changes keep their normal project workflow: direct example runs need
`PYTHONPATH=.` when there is no editable install. Modified handwritten Python
requires `black .` in an isolated user-owned intended-change worktree or copy,
then individual single-file static analysis. New handwritten files require
Pylint **10.00/10.00**; existing tracked handwritten files must not regress
from their grandfathered pre-change scores.

## Search Order

Use this order for agent navigation:

1. Read [AGENTS.md](../AGENTS.md) for root routes and governance pointers.
2. Read [catalog.md](catalog.md) for ambiguous or cross-branch queries.
3. Follow branch `index.md` routers from root to owner.
4. Read the owning leaves and their typed `See Also` neighbors.
5. Use `rg` over `wiki/` only when catalog and router paths fail.
6. Use `rg` over source code only to verify or update sourced facts.

## Ingest

1. Register the source in [raw/SOURCES.md](../raw/SOURCES.md) with provenance,
   `frozen` or `living` status, and ingest state. Exact cited-file rows may
   abbreviate to source and status per [SCHEMA.md](SCHEMA.md); external-source
   rows may also carry an accessed date and notes. Do not record or compute
   source-tracking digests or timestamps.
2. Decide which branch owns the compiled facts.
3. Update the owning leaf in synthesized prose and cite exact files plus stable
   symbols or headings.
4. From glossary entities touched by the source, follow existing links one hop
   and update affected neighbors.
5. Update a router only when a new leaf or route exists.
6. Run the checks in [lint/CHECKS.md](lint/CHECKS.md) on the touched subgraph.

## Query

Start with the search order above. Answer from compiled pages whenever they are
current and sufficiently sourced. If a normal question requires broad tree
search, repair the relevant router, catalog row, or focused sourced leaf instead
of making broad search the permanent workflow.

## Query Filing

Do not file one-off answers.

File durable answers when they span three or more leaves, resolve recurring
ambiguity, compare subsystems, record a gotcha, or change how future queries
should route.

Use this filing target:

1. Durable single-topic synthesis updates the existing owning leaf.
2. Durable narrow new topic creates a leaf under the owning branch.
3. Durable cross-branch synthesis goes under an existing branch when one branch
   owns it cleanly.
4. Durable cross-branch synthesis with no clean owner creates
   `wiki/syntheses/index.md` plus the first `wiki/syntheses/*.md` leaf in the
   same change.

When filing, update [catalog.md](catalog.md), source dependencies in
[source-map.md](source-map.md) when they change, and typed neighbor links. Use
plain markdown; do not require Obsidian metadata or frontmatter.

## Living-Source Change

Review the change dependency-aware: inspect the changed paths, confirm the
source's [raw/SOURCES.md](../raw/SOURCES.md) row (provenance, status, ingest
state) still describes it, follow [source-map.md](source-map.md) ownership rows
and maintainer signals to the affected compiled pages, then re-ingest affected
leaves and update one-hop neighbors that share changed glossary entities. Mark
a page `stale` only when it cannot be reconciled immediately. Record durable
source drift and reconciliation decisions in the affected pages,
[source-map.md](source-map.md), [contradictions.md](contradictions.md), commit
messages, or PR descriptions as appropriate.

## External Source Trust Change

Classify external sources as `background` unless a durable KB claim depends on
them. Promote an external source to `external-spec` only when the claim needs
the source as authority and one of these is true:

1. A frozen markdown excerpt or version note exists under `raw/source-docs/`.
2. The upstream source is stable by version, commit, or immutable URL.
3. A documented exception explains why the live external source is acceptable.

Any external source trust change requires updates to
[source-map.md](source-map.md). If the change affects a contested or stale
claim, also update [contradictions.md](contradictions.md). Repository facts
still follow the authority order in [SCHEMA.md](SCHEMA.md): code, tests,
generated evidence, docs, external specs, then background sources.

## Contradiction And Stale Claims

When sources disagree, create or update a structured
`CONTR-0001`-form row in [contradictions.md](contradictions.md) before changing
the affected claim's status to `contested`. Record exact claim/status,
competing sources, authority decision, complete affected-page links,
page-status rationale, owner/trigger, resolution test, opened/resolved dates,
and notes. Put `Claim status: contested; contradiction: CONTR-0001.` on every
active affected page; use `stale` in the same form when applicable.

When a living source has moved and reconciliation cannot finish in the same
change, mark affected pages `stale`, add or update the contradiction/staleness
row, and record the blocked reconciliation reason there. When a claim is
reconciled, update affected pages and close or revise the row in
[contradictions.md](contradictions.md).

Before closure, review source-map reverse dependents, all pages citing either
source, catalog aliases/key symbols, current affected pages, typed neighbors,
and targeted exact/key-phrase wiki hits. Remove active markers only after every
affected page is reconciled and the resolution test passes.

## Claim Adjudication

The [Claim And Evidence Contract](SCHEMA.md#claim-and-evidence-contract) is
prospective after its 07-13-2026 adoption change. High-risk claims predating it
and claims changed in that same adoption change remain baseline-uncovered unless
an exact block is present; that adoption change asserts no completed block
coverage. When a baseline claim or its deciding source is next materially
changed after adoption, add the exact block immediately after the claim in its
owning `Detail` section. For an active contradiction, add it to the matching
`### CONTR-*` subsection, never the fixed register row. Add validation and
dimensions only for behavioral claims. Use code for descriptive behavior;
owning
governance/configuration for normative rules; stable specification plus targeted
tests for intended public/scientific contracts; workflow/configuration for CI
job shape; and frozen generated evidence only for its pinned context. Synthesis
agreement is never authority. Navigation, structure, provenance,
status, symbolic definition, and normative rules do not receive behavioral
validation lines. In every behavioral dimension, use an exact value when
exercised, `not-run` when applicable but unexercised, and `not-applicable` only
when the dimension does not apply; a date value uses `MM-DD-YYYY`.

## Safe Reproduction

In the shared repository working tree, run only demonstrably side-effect-free
scoped checks. Run generators, builds, oracle creation or update, and blanket
formatters only in an isolated user-owned intended-change worktree or copy with
no unrelated changes.
Inspect side effects first, use owned disposable cache and output, and impose
time/resource limits. Retain no incidental output and never clear or overwrite
ambient or shared cache. No coordination exception permits mutation in the
shared working tree. Record command,
working directory, observed assertion, result, limits, cleanup, and behavioral
tuple when relevant. Network, installs, remote CI, and external toolchains need
user authority. Never reset or clean shared `project/` output.

## Deterministic Checks

Run:

```bash
python tools/kb_lint.py
git diff --check
```

Expected success is exit 0; linter prints `KB lint passed.` `--all` is an
identical compatibility alias, not stronger coverage. Inventory commands are
diagnostics, never semantic completeness proof.

## Page Move Or Delete

Move a leaf only when its owning route changes or it is promoted into a router.
Update parent routers, `Up:` headers, `See Also` links, and glossary references
in the same change. Delete a page only after no router or leaf links to it.

## Evidence-First Debugging

Do not propose a root cause or patch for runtime failures until the failure is
reproduced or explicitly marked unreproduced. If a fix depends on cached,
optimized, approximate, or generated helpers matching runtime behavior,
instrument the reproduced case and compare runtime-visible quantities first.
Durable evidence belongs under `raw/source-docs/` as a frozen source; the wiki
keeps only the durable lesson.

## NRPy Workflow Notes

For direct example runs without an editable install, append `.` to
`PYTHONPATH`. Run `black .` only in an isolated, user-owned intended-change
worktree or copy with no unrelated modifications, inspect its diff, and run
`./.github/single_file_static_analysis.sh <path.py>` for each modified
handwritten Python file. A new handwritten file must report Pylint
**10.00/10.00**. An existing tracked handwritten file must not regress from its
pre-change score, including a grandfathered score at or below `9.5`. Inspect the
wrapper's command construction, the target's direct-execution effects, and every
invoked tool's cache and filesystem output effects. The current wrapper
dispatches interpolated command strings through `eval`. Before invocation,
require a repository-relative argument whose resolved target stays in the
repository, is a regular non-symlink file, contains no `..` component, does not
begin with `-`, and matches `^[A-Za-z0-9_./-]+$`. A failing path blocks invocation
until separately authorized wrapper hardening. Run an accepted path only in the
isolated intended-change tree or copy; disable or redirect every writable tool
cache and filesystem output to an owned disposable location, then inspect
repository status. Generated trusted `*/tests/*.py` data is exempt from per-file
analysis, but its handwritten owner is not. [Code Test
Policy](validation/code-test-policy.md) owns test selection; [Test Oracles And
Safe Updates](validation/test-oracles-and-safe-updates.md) owns oracle format
and the two-process update procedure.
Generated C/CUDA projects, generated thorns, generated Charm++
projects, and generated JAX projects are normally products of Python
generators; cite the generator unless a generated artifact is deliberately
registered as frozen evidence.

## Sources

- [kb-instructions.md](../raw/source-docs/kb-instructions.md) - section 7.6 for `wiki/workflows.md`
- [coding_style.md](../coding_style.md) - `## Python Coding Style`, `### Formatting`
- [original-agents.md](../raw/source-docs/original-agents.md) - historical `## Required Checks`; current `coding_style.md` decides conflicts
- [README.md](../README.md) - `## Contributor Setup`
- [single_file_static_analysis.sh](../.github/single_file_static_analysis.sh) - `run_test_step`
- [main.yml](../.github/workflows/main.yml) - `static-analysis`

## See Also

- Depends on: [Schema](SCHEMA.md)
- See also: [Lint Checks](lint/CHECKS.md)
- See also: [Contribution Style And Static Analysis](architecture/contribution-style-and-static-analysis.md)
- See also: [Static Analysis](validation/static-analysis.md)
- See also: [Code Test Policy](validation/code-test-policy.md)
- See also: [Test Oracles And Safe Updates](validation/test-oracles-and-safe-updates.md)
