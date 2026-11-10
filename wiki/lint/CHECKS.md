# Lint Checks

> Mechanical and review checks for the Markdown KB. · Status: confirmed

## Summary

Run these checks before considering KB work complete. The goal is a reachable,
source-backed tree: root `AGENTS.md` routes, routers stay fact-light, leaves
carry facts and sources, and excluded files stay out of the KB.

## Deterministic Checks

- `python tools/kb_lint.py` prints `KB lint passed.` It checks Markdown links
  and anchors, Obsidian-style links, router sections and edges, leaf section
  order, `Up:` headers, `See Also` parents and neighbors, source registration,
  reachability from `AGENTS.md`, one catalog row per normal leaf and catalog
  status, a defining-page link or external/background marker for each glossary
  term, status vocabularies, source-map columns and targets, contradiction rows
  and markers, an `index.md` router in each content directory, non-Markdown or
  maintenance files in `wiki/` and `raw/`, root-level maintenance files, and
  volatile metadata.
- `git diff --check` reports no whitespace errors.

## Review Checks

- `git status --short` shows only intended documentation/source-doc changes and
  no generated files.
- A marker scan returns no to-do markers or template text.
- Every added or substantially revised passage has completed the review in
  [Language For Computational Physicists](../computational-physics-language.md).
- Routers carry no `Sources` section and no unique durable implementation facts.
- Catalog rows are one-line navigation aids only.
- [source-map.md](../source-map.md) covers seeded high-value dependencies,
  especially exact cited files and aggregate `partial` rows from
  [raw/SOURCES.md](../../raw/SOURCES.md).
- Cited files are cited by stable symbol or heading, not generated line
  numbers; cited symbols exist; semantic truth, dynamic or generated symbols,
  and C/CUDA macros are checked by reading the source.
- No generated project output, binaries, images, archives, rendered files,
  scratch reports, or planning/ranking files are added. Root planning and task
  files that follow the
  [Schema naming grammar](../SCHEMA.md#coordination-artifacts) stay outside KB
  scope.
- Trusted-value files are source evidence, not hand-authored prose pages.
- Handwritten Python edits follow
  [Static Analysis](../validation/static-analysis.md).

## Sources

- [kb-instructions.md](../../raw/source-docs/kb-instructions.md) - section 7.8 for `wiki/lint/CHECKS.md`

## See Also

- [KB Page Format](../SCHEMA.md)
- [Workflows](../workflows.md)
- [Language For Computational Physicists](../computational-physics-language.md)
