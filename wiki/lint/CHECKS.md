# Lint Checks

> Mechanical and review checks for the Markdown KB. · Status: confirmed · Last reconciled: 07-05-2026

## Summary

Run these checks before considering KB work complete. The goal is a reachable,
source-backed tree: root `AGENTS.md` routes, routers stay fact-light, leaves
carry facts and sources, and excluded artifacts stay out of the KB.

## Mechanical Checks

- `git status --short` shows only intended documentation/source-doc changes and
  no generated artifacts.
- `git diff --check` reports no whitespace errors.
- `find wiki raw -type f ! -name '*.md' -print` returns no files.
- `find . -maxdepth 1 -type f -name '*kb*instructions*.md' -print` returns no
  root KB instruction file.
- `rg -n '\[\[' AGENTS.md wiki raw` returns no Obsidian-style links.
- A marker scan returns no unfilled dates, to-do markers, or template text.
- No removed source-tracking metadata appears in governed KB files, including
  `raw/source-docs/**/*.md`: no hash digest values of any algorithm (`sha256`
  or otherwise), no `Mtime`/`Hash` manifest table columns, and no `mtime`
  values. Mentions of the removed metadata are allowed only in prohibition or
  supersession statements.
- No `YYYY-MM-DD` date literals remain; retained KB dates use `MM-DD-YYYY` or
  approved placeholders such as `n/a`/`-`.
- Agents are never told to compute, recompute, compare, or lint source-tracking
  digests or timestamps; source drift is dependency-aware review per
  [Workflows](../workflows.md).
- Relative Markdown links from `AGENTS.md`, `wiki/**/*.md`, and `raw/**/*.md`
  resolve within the repository.

## Router Checks

- Root `AGENTS.md` is the only root-level KB document.
- Every content-bearing `wiki/` directory has an `index.md` router, except the
  allowed `wiki/lint/` utility directory or another explicit schema exemption.
- Every new `wiki/` directory with content pages has an `index.md` router in
  the same change, unless explicitly exempted by [Schema](../SCHEMA.md).
- No router has a `Detail` section.
- Routers generally have no `Sources` section and carry no unique durable
  implementation facts.
- Routers link only to immediate children plus parent, root, or global hubs.

## Compounding Memory Checks

- Every `wiki/**/*.md` page appears exactly once in [catalog.md](../catalog.md),
  including routers, leaves, and governance/support pages.
- Catalog rows are one-line navigation aids only; they do not replace leaf
  reading, leaf `Sources`, or source-backed detail.
- Every durable ingest, query filing, lint pass, reconciliation, source-drift
  decision, page move, page add, or page delete has a parseable
  [log.md](../log.md) entry.
- `wiki/log.md` contains decisions and checks only; it has no chat transcripts,
  scratch output, token reports, planning dumps, or routine lint reports.
- Every glossary term has an owner page link, an explicit external/background
  marker, or a matching concept-hub candidate row in `wiki/catalog.md`.
- [source-map.md](../source-map.md) covers seeded high-value dependencies,
  especially exact cited files and aggregate `partial` rows from
  [raw/SOURCES.md](../../raw/SOURCES.md).

## Leaf Checks

- Every non-router leaf has the required title, two-line header, and section
  order: `Summary`, `Detail`, `Sources`, `See Also`.
- Every leaf has non-empty sources and cites files by stable symbol or heading,
  not generated line numbers.
- Every leaf's `See Also` includes its parent router and at least one
  neighboring wiki page.
- Every cited file or symbol exists, and every cited authority is registered in
  [raw/SOURCES.md](../../raw/SOURCES.md) directly or through a clearly named
  aggregate row.

## NRPy Baseline

- No generated project output, binaries, images, archives, rendered artifacts,
  scratch reports, token/latest reports, or planning/ranking artifacts are
  added.
- Trusted-value files are source evidence, not hand-authored prose pages.
- Python edits, if any accidentally occur, require `black .` and the
  single-file static-analysis script for each modified Python file, except
  trusted-value files under `*/tests/*.py`.

## Sources

- [kb-instructions.md](../../raw/source-docs/kb-instructions.md) - section 7.8 for `wiki/lint/CHECKS.md`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `## Required Checks`

## See Also

- [Schema](../SCHEMA.md)
- [Workflows](../workflows.md)
