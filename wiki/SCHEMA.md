# Knowledge Base Schema

> Rules for the self-maintaining markdown KB. · Status: confirmed · Last reconciled: 2026-06-29

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

## Citation Rules

Code citations use a relative path plus a stable symbol, function, class,
macro, or header name. Documentation citations use a relative path plus a
stable heading or document role. Generated line numbers are not stable
authorities. Every cited source belongs in `raw/SOURCES.md` directly or through
a clearly named aggregate source row.

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
