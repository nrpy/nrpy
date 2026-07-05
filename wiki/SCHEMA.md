# Knowledge Base Schema

> Rules for the self-maintaining markdown KB. · Status: confirmed · Last reconciled: 07-05-2026

## Summary

This KB has three ownership layers: immutable or living sources, compiled wiki
pages, and schema/governance pages. Those ownership layers are not navigation
depth layers. Navigation is an unbounded router tree that starts at
`AGENTS.md`, adds `index.md` routers only where real branching exists, and puts
durable implementation facts in leaves.

## Commitments

1. Sources are authoritative for wiki purposes; wrong or stale claims are fixed
   by updating or adding sources and reconciling pages.
2. Sources, compiled pages, and schema are separate layers with separate owners.
3. The model-owned compiled wiki summarizes, cross-references, files, updates
   neighbors, and lints sourced facts.
4. Compile once into structured pages; do not re-answer every normal query by
   searching the whole repository.
5. Ingest one source at a time and trace its implications across linked pages.
6. Use ordinary relative Markdown links, not Obsidian-style links.
7. Navigate from `AGENTS.md` through as many routers as needed; if a normal
   query needs broad repo search, the router tree has drifted.
8. Lint for broken links, orphan pages, stale claims, low-confidence claims,
   dead code references, and entity drift.
9. Start with a small, high-value source set and grow conventions deliberately.

## Ownership

Sources include code, config, docs, preserved root instructions, selected
evidence, and build inputs. Compiled wiki leaves synthesize those sources and
cite them. Governance pages define the page contract, source vocabulary, router
rules, maintenance procedures, and lint checks.

Source status is either `frozen` for preserved snapshots and selected evidence
or `living` for files that may change. Ingest state is `registered`, `partial`,
or `ingested`.

## Source-Tracking Metadata And Dates

A repository source-manifest row records source identity, provenance, status,
and ingest state as its only source-tracking fields. Exact cited-file rows may
abbreviate to source and status: the repository path is the provenance, and
ingest state is tracked through covering aggregate rows and
[source-map.md](source-map.md). External-source rows may also carry an
accessed date (`MM-DD-YYYY`) and notes; those are ordinary KB metadata, not
source-tracking metadata. Manifests and KB docs carry no source-tracking
hash columns or values of any kind (`sha256` or any other digest algorithm)
and no `mtime` columns or values; sources are never hashed for tracking.
Technical, non-source-tracking hash facts (cache hashing, coordinate hash
macros) remain allowed as reviewed domain facts, but never as stored digest
values.

Source drift and staleness are resolved by dependency-aware reconciliation, not
stored fingerprints: inspect the changed paths, the source's manifest row and
status, [source-map.md](source-map.md) ownership rows, maintainer signals, and
the affected compiled pages, then re-ingest and reconcile.

All retained KB dates use the `MM-DD-YYYY` format. Approved placeholders such
as `n/a` or `-` are allowed where no date applies.

## Support Pages

Support pages keep the wiki compounding without changing router semantics:

- `wiki/catalog.md` is the global content catalog. It lists pages with type,
  one-line answer, route, query terms, status, reconciliation date, source
  count, and concept-hub candidacy. It does not replace reading leaves.
- `wiki/log.md` is append-only chronological memory for durable KB operations:
  ingest, query filing, lint, source drift, reconciliation, page movement, and
  other maintenance decisions. It must stay free of chat transcripts, scratch
  logs, token reports, and planning dumps.
- `wiki/source-map.md` records source dependencies: source or aggregate,
  source authority tier, ingest status, dependent pages, covered subpaths,
  known gaps, last check, and next action.
- `wiki/contradictions.md` records contested or stale claims, competing
  sources, authority decisions, affected pages, and resolution state.

`index.md` files remain router-only. Do not turn router indexes into global
content catalogs.

Every durable KB operation updates `wiki/log.md`. Every page add, move, or
delete updates `wiki/catalog.md`. Every durable source dependency change
updates `wiki/source-map.md`.

## Page Contract

Leaves use this exact section order: title, two-line blockquote header,
`Summary`, `Detail`, `Sources`, and `See Also`. Page status is one of
`confirmed`, `provisional`, `contested`, or `stale`. Use `confirmed` only when
registered sources directly support the page. Use `provisional` when structure
is useful but not fully reconciled, `contested` when sources disagree, and
`stale` when a living source has moved since the page was checked.

Each leaf must cite at least one source. `See Also` must include the parent
router and at least one neighboring wiki page; the header `Up:` link does not
count.

Use typed `See Also` labels when the relation is known:

- `Parent:` owning router or parent concept.
- `Depends on:` prerequisite concept, API, or workflow.
- `Implements:` page, source, or workflow this page realizes.
- `Validated by:` tests, examples, CI, or evidence pages.
- `Example:` worked example or generated application that demonstrates the
  page.
- `Contrasts with:` related page that differs in an important way.
- `See also:` useful neighbor when a stronger typed relation does not apply.

New or substantially edited leaves should use typed labels for non-parent
links when doing so improves traversal.

## Citation Rules

Code citations use a relative path plus a stable symbol, function, class,
macro, or header name. Documentation citations use a relative path plus a
stable heading or document role. Generated line numbers are not stable
authorities. Every cited source belongs in `raw/SOURCES.md` directly or through
a clearly named aggregate source row.

Source authority tiers are:

- `primary-code`: repository source code, configuration, and checked-in build
  inputs.
- `primary-test`: repository tests and validation harnesses.
- `generated-evidence`: selected frozen generated files, run outputs, or build
  evidence deliberately registered as evidence.
- `primary-doc`: repository documentation and preserved root instructions.
- `external-spec`: stable external specification or upstream documentation
  used for a durable KB claim.
- `background`: orientation material that informs structure or interpretation
  but does not decide repository facts.

Default contradiction precedence is `primary-code` > `primary-test` >
`generated-evidence` > `primary-doc` > `external-spec` > `background`.
Exceptions are allowed only with source-specific reasoning in `wiki/log.md`
and, when a claim is contested, `wiki/contradictions.md`.

Avoid per-sentence citations. Add inline source tags only for high-risk claims:
public APIs, generated-output boundaries, CI guarantees, source authority
decisions, contradiction decisions, and claims where a nearby source paragraph
could plausibly be misread.

## Query Filing

One-off answers do not need filing. File durable answers when they span three
or more leaves, resolve recurring ambiguity, compare subsystems, record a
gotcha, or change how future queries should route.

File durable single-topic synthesis into the owning leaf. Create a new leaf
only for a durable narrow topic with a clear owning branch, or for a
cross-branch synthesis that cannot cleanly live under an existing branch. In
all cases, update neighboring links, `wiki/catalog.md`, `wiki/log.md`, and
`wiki/source-map.md` when source dependencies change.

## Router Rules

Routers are pointers only. They contain a title, a one-line blockquote, a
`Page | Go here when...` table, and root or parent navigation. Routers do not
have `Detail` sections, ordinary `Sources` sections, source dumps, or unique
implementation facts. A router links only to immediate children plus parent,
root, or global navigation.

Create navigation depth only from real fan-out. Do not create empty leaves,
padding routers, or single-child padding directories. Promote a broad leaf
into a router only when there is enough sourced content to split into children.

## Canonical Terms

Canonical spellings and recurring entity names live in
[glossary.md](glossary.md). When a term recurs across branches or has multiple
spellings, update the glossary and one-hop neighbors that rely on it.

## Not In The KB

Generated project outputs under `project/` are excluded unless a selected
generated file is deliberately registered as frozen evidence. The KB also
excludes binaries, images, archives, compiled artifacts, rendered outputs,
scratch logs, token reports, latest-snapshot reports, and planning or ranking
artifacts.

## Sources

- [kb-instructions.md](../raw/source-docs/kb-instructions.md) - `# Instructions: Build a Self-Maintaining Markdown Knowledge Base From Scratch`

## See Also

- [Workflows](workflows.md)
- [Glossary](glossary.md)
- [Lint Checks](lint/CHECKS.md)
