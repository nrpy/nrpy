# Knowledge Base Layout And Rules

> Page layout, source provenance, router bounds, and checker rules for the Markdown KB. · Status: confirmed

## Summary

The NRPy knowledge base separates authoritative sources (`nrpy/`, `coding_style.md`, `.github/`, `raw/`), compiled `wiki/` pages, and governance rules (`AGENTS.md`, `SCHEMA.md`, `workflows.md`). Navigation starts at `AGENTS.md`, branches through `index.md` routers where fan-out warrants a directory (`2` to `14` children), and stores each durable technical fact once in its defining leaf.

## Volatile Information Policy

Source-manifest rows in [raw/SOURCES.md](../raw/SOURCES.md) record source identity, provenance, status (`living`, `frozen`, or `frozen-historical`), ingest state (`registered`, `partial`, or `ingested`), and—for external URLs—authority tier (`external-spec` or `background`) plus concise notes (`<= 80` B). Add an exact local file row only when no aggregate row in that file covers the path.

Authored KB files (`AGENTS.md`, `wiki/**/*.md`, and `raw/SOURCES.md`) forbid maintenance or runtime snapshots: no dates, times, timestamps, `mtime` values, commit or file hashes/digests, inventory/file/page/source/job counts, audit/reconciliation/resolution fields, environment tuples, recorded run results, or separate maintenance logs (`wiki/log.md` is retired; Git history records changes). Stable scientific constants, stencil orders, tensor component counts, versioned specification URLs, arXiv identifiers, and code identifiers such as `CoordSystem_hash` remain valid. Preserved files under `raw/source-docs/` remain untouched verbatim snapshots on disk, while their `raw/SOURCES.md` rows obey this policy.

Universal symbolic-versus-runtime scope rule: symbolic expressions, doctests, and code-generation tests establish symbolic identities and emitted source structure only; compilation, runtime behavior, and numerical convergence require executed build or solver checks. Do not repeat generic symbolic-only disclaimers across individual leaves.

## Support Pages

Every table row across `catalog.md`, `examples/example-generator-catalog.md`, `glossary.md`, `source-map.md`, `contradictions.md`, and `raw/SOURCES.md` must stay `<= 450` B:

- `wiki/catalog.md` and `wiki/examples/example-generator-catalog.md`: `<= 240` B/row; `catalog.md` uses the 3-column leaf index (`Page | Scope | Query terms`) listing each normal leaf once with non-duplicate symbol and topic search terms (no `Type` or `Status` columns; update it only when adding, moving, renaming, or deleting a leaf).
- `wiki/source-map.md`: aggregate ingest table (`Source or aggregate | Glob pattern | Authority tier | Status | Ingest | Uncovered scope`, `<= 450` B/row) covering only the multi-file aggregates from `raw/SOURCES.md`. Single-file rows and manual `Dependent pages` lists are prohibited; compute reverse citations dynamically with `rg -l -F '<path>' wiki/`.
- `wiki/contradictions.md`: active conflict register (`ID | Claim status | Affected pages | Summary`, `<= 220` B/row) paired with self-contained `### CONTR-XXXX` subsections.
- `wiki/glossary.md`: canonical `Term | Meaning` table (`<= 450` B/row, `Meaning <= 180` B) where every `Meaning` cell links its defining leaf or marks `external`/`background`, carries no `## Sources` section, and uses typed `See Also` bullets.
- `wiki/workflows.md` and `wiki/computational-physics-language.md`: procedures, deterministic/review checks, and computational-physics prose rules.

## Page Layout

Normal leaves stay `<= 20,000` B and use this exact structure: `# Title`, two-line blockquote (`> One-line scope. · Status: <confirmed|provisional|contested|stale>` and `> Up: [Parent](...)`), `## Summary`, `## Detail`, `## Sources`, and `## See Also`.

- **Subsection partitioning and line bounds**: Whenever `## Detail` exceeds `4,000` B, partition it with one or more `###` headings naming concrete functions, structs, parameters, or equations so no contiguous `## Detail` or `###` chunk exceeds `6,000` B (measured on raw UTF-8 text excluding HTML comments). Keep every `wiki/` line `<= 1,000` B, avoid double blank lines (`\n\n\n`), and never restate table cells in adjacent prose paragraphs.
- **Parent and neighbor links**: Header `> Up:` is the sole parent-router link (`../AGENTS.md` for leaves directly in `wiki/`, otherwise `index.md`). `## See Also` must hold `1` to `6` non-parent neighbor `- ` bullets (targeting `<= 3` high-signal links, `1` unique link per bullet, no self-links), never duplicating `> Up:`, and every bullet must start with `Depends on:`, `Implements:`, `Validated by:`, `Example:`, `Contrasts with:`, or `See also:`.
- **Page status**: Use `contested` when an unresolved source conflict alters the `Summary`, public interface, or multiple core `Detail` claims; `stale` when living-source drift leaves principal claims unreconciled; `provisional` when synthesis is incomplete without known conflict; and `confirmed` when registered sources directly support the page (allowing at most one bounded, noncentral `Detail` defect marked inline with its `CONTR-XXXX` link).

## Citation Rules

Source authority precedence is `primary-code` > `primary-test` > `generated-evidence` > `primary-doc` > `external-spec` > `background`; exceptions require explicit justification in the defining leaf or [contradictions.md](contradictions.md).

- **Unique `## Sources` entries**: List each cited repository file or external URL at most once as a clickable Markdown link in `## Sources`, grouping multiple functions, classes, macros, or headings on that bullet. `confirmed` pages must not cite `frozen-historical` sources unless the bullet explicitly notes historical provenance.
- **Single Defining-Page Rule**: State and evidence each repository fact in one defining leaf. Neighboring leaves link to that defining leaf in one sentence rather than copying claims or `Evidence` lines.

### Claim And Evidence Rules

Legacy multi-line `Claim evidence:` blocks are prohibited. Add a single-line evidence tag (`<= 300` B) immediately after new or materially revised high-risk claims (public interfaces, mathematical/numerical guarantees, CLI options, generated-file boundaries, CI job behavior, authority/contradiction decisions, and easily misread source passages):

```text
Evidence (<role>): `<subpath.py::symbol>`; corroboration: `<source>`
```

Allowed roles are `descriptive behavior`, `normative rule`, `scientific/interface requirement`, `CI behavior`, and `generated evidence`. Cite sources in `## Detail` via shortest-unambiguous `` `subpath.py::symbol` `` tags matching `## Sources` (adding parent directories when basenames collide in that leaf), omitting `corroboration:` when none exists. Inside tables, attach compact ``(Ev: `module.py::symbol`)`` tags directly in cells instead of repeating table rows in prose.

## Router Rules

Routers (`AGENTS.md` and `wiki/**/index.md`) contain only a title, one-line blockquote, `Page | Go here when...` table, and parent/root links—never `## Detail`, `## Sources`, or standalone implementation facts. `AGENTS.md` stays `<= 4,800` B. Every `wiki/**/index.md` router must link `2` to `14` immediate children (leaves or child `index.md` routers); single-child padding directories are prohibited. `AGENTS.md` may link directly to `wiki/*.md` leaves, branch `wiki/*/index.md` routers, interior `wiki/**/index.md` sub-routers, and support pages, while every normal leaf and `index.md` sub-router must also be linked from its immediate parent router.

## Contradiction Rules

Active rows in [contradictions.md](contradictions.md) `## Register` use the 4-column header:

`ID | Claim status | Affected pages | Summary`

Rows stay `<= 220` B with `ID` in `CONTR-0001` form and `Claim status` in `{"contested", "stale"}`. Each active ID has a matching `### CONTR-XXXX` subsection recording competing sources, authority decision, `Evidence (<role>): ...` (`<= 300` B), recheck trigger, and `Resolution check:`. Every active affected leaf carries `Claim status: contested; contradiction: CONTR-0001.` (or `stale`) outside code fences with a link to `contradictions.md#contr-0001`. When the resolution check passes and affected leaves are reconciled, remove all `#contr-*` links, delete the register row and `### CONTR-*` subsection in the same commit, and append the ID to `Retired IDs:` in `contradictions.md` (retired IDs are never reused).

## Checker And Frozen Scope

`PYTHONDONTWRITEBYTECODE=1 python3 -B tools/kb_lint.py` checks only that local links and heading anchors resolve in `AGENTS.md`, `wiki/**/*.md`, and `raw/SOURCES.md` (no wikilinks), that every wiki page is reachable from `AGENTS.md`, and that each backticked source file in an `Evidence (...)` line or `(Ev: ...)` span names exactly one git-tracked file present on disk (the page's `## Sources` links break ties) whose cited `::` symbols appear in it; stage new source files with `git add` before citing them. Every other rule on this page and in [Language For Computational Physicists](computational-physics-language.md) is checked in review, not by the lint. Files under `raw/source-docs/**/*.md` remain immutable verbatim snapshots.

### Root Planning And Task Files

At repository root only, commissioned planning or task Markdown files whose lowercase basename starts with `plan`, `task`, or `tasks` (optionally followed by digits and hyphen/underscore lowercase alphanumeric tags, with no segment matching `audit`, `scan`, `report`, `log`, `latest`, `token`, `rank`, or `ranking` in singular, plural, or numbered form) are exempt from KB routing checks and must never be filed, staged, moved, or deleted as KB content unless the user directs.

## Not In The KB

Generated directories (`project/`, `build/`, `dist/`), caches, binaries, images, archives, scratch logs, token/audit/scan/latest reports, and uncommissioned planning or ranking files remain excluded from `wiki/` and `raw/`.

## Sources

- [AGENTS.md](../AGENTS.md) - root router, computational-physics prose mandate, and protected-file scope
- [tools/kb_lint.py](../tools/kb_lint.py) - link, reachability, and Evidence-citation checks

## See Also

- Depends on: [Language For Computational Physicists](computational-physics-language.md)
- See also: [Workflows](workflows.md)
- See also: [Glossary](glossary.md)
