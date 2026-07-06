# Workflows

> Procedures for ingesting, querying, and maintaining the KB. · Status: confirmed · Last reconciled: 07-05-2026

## Summary

Start from `AGENTS.md`, follow catalog and router paths to compiled leaves, and
answer from sourced pages. Durable answers are filed back into the wiki when
they should compound. Maintenance work registers sources first, updates the
owning leaf, then fixes nearby links, catalog/log entries, and glossary terms.
NRPy code changes keep their normal project workflow: direct example runs need
`PYTHONPATH=.` when there is no editable install, Python edits require Black
before commit, and every modified Python file gets the single-file
static-analysis script unless it is a generated trusted-value file under
`*/tests/*.py`.

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

When filing, update [catalog.md](catalog.md), [log.md](log.md), source
dependencies in [source-map.md](source-map.md) when they change, and typed
neighbor links. Use plain markdown; do not require Obsidian metadata or
frontmatter.

## Living-Source Change

Review the change dependency-aware: inspect the changed paths, confirm the
source's [raw/SOURCES.md](../raw/SOURCES.md) row (provenance, status, ingest
state) still describes it, follow [source-map.md](source-map.md) ownership rows
and maintainer signals to the affected compiled pages, then re-ingest affected
leaves and update one-hop neighbors that share changed glossary entities. Mark
a page `stale` only when it cannot be reconciled immediately. Record durable
source drift and reconciliation decisions in [log.md](log.md).

## External Source Trust Change

Classify external sources as `background` unless a durable KB claim depends on
them. Promote an external source to `external-spec` only when the claim needs
the source as authority and one of these is true:

1. A frozen markdown excerpt or version note exists under `raw/source-docs/`.
2. The upstream source is stable by version, commit, or immutable URL.
3. A logged exception explains why the live external source is acceptable.

Any external source trust change requires updates to [source-map.md](source-map.md)
and [log.md](log.md). If the change affects a contested or stale claim, also
update [contradictions.md](contradictions.md). Repository facts still follow the
authority order in [SCHEMA.md](SCHEMA.md): code, tests, generated evidence,
docs, external specs, then background sources.

## Contradiction And Stale Claims

When sources disagree, create or update a row in
[contradictions.md](contradictions.md) before changing the affected claim's
status to `contested`. Record the competing sources, authority decision,
affected pages, opened date, and next action. Add a matching `reconcile` or
`source-drift` entry to [log.md](log.md).

When a living source has moved and reconciliation cannot finish in the same
change, mark affected pages `stale`, add or update the contradiction/staleness
row, and log the blocked reconciliation reason. When a claim is reconciled,
update affected pages, close or revise the row in [contradictions.md](contradictions.md),
and append the resolution to [log.md](log.md).

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
`PYTHONPATH`. Python code changes require `black .` before commit, and each
modified Python file requires `./.github/single_file_static_analysis.sh
<path.py>` before commit. Trusted values under `*/tests/*.py` are regenerated
from owning modules, not hand-edited, and are exempt from single-file static
analysis. Generated C/CUDA projects, generated thorns, generated Charm++
projects, and generated JAX projects are normally products of Python
generators; cite the generator unless a generated artifact is deliberately
registered as frozen evidence.

## Sources

- [kb-instructions.md](../raw/source-docs/kb-instructions.md) - section 7.6 for `wiki/workflows.md`
- [original-agents.md](../raw/source-docs/original-agents.md) - `## Required Checks`
- [README.md](../README.md) - `## Contributor Setup`
- [single_file_static_analysis.sh](../.github/single_file_static_analysis.sh) - `run_test_step`
- [main.yml](../.github/workflows/main.yml) - `static-analysis`

## See Also

- [Schema](SCHEMA.md)
- [Lint Checks](lint/CHECKS.md)
- [Contribution Style And Static Analysis](architecture/contribution-style-and-static-analysis.md)
