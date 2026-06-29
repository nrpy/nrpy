# Lint Checks

> Mechanical and review checks for the Markdown KB. · Status: confirmed · Last reconciled: 2026-06-29

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
- A marker scan returns no unfilled dates, hashes, to-do markers, or template
  text.
- Relative Markdown links from `AGENTS.md`, `wiki/**/*.md`, and `raw/**/*.md`
  resolve within the repository.

## Router Checks

- Root `AGENTS.md` is the only root-level KB document.
- Every content-bearing `wiki/` directory has an `index.md` router, except the
  allowed `wiki/lint/` utility directory.
- No router has a `Detail` section.
- Routers generally have no `Sources` section and carry no unique durable
  implementation facts.
- Routers link only to immediate children plus parent, root, or global hubs.

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
- [coding_style.md](../../coding_style.md) - `## Static Analysis Configuration`

## See Also

- [Schema](../SCHEMA.md)
- [Workflows](../workflows.md)
