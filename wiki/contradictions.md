# Contradictions

> Register for contested and stale KB claims. Plain Markdown only.

Known contested/stale claims as of 06-30-2026 are tracked below.

## Register

| Claim | Status | Source A | Source B | Authority decision | Affected pages | Opened | Resolved | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `manga_bhah_lib` final print suggests running `./bhah_lib`, but the generated target is a library build. | stale | `nrpy/examples/manga_bhah_lib.py` final print | [Matter TOV Workflows](examples/matter-tov-workflows.md) source/Makefile interpretation | Page treats print as stale generic messaging until generator text is reconciled. | [Matter TOV Workflows](examples/matter-tov-workflows.md) | 06-30-2026 | - | Recheck when example generator output text changes. |

## Rules

- Open a row when registered sources disagree, when a living source supersedes a KB claim, or when a page is marked `contested` or `stale`.
- Update a row when new evidence, source authority, affected pages, or reconciliation status changes.
- Resolve a row only after affected pages have been reconciled or the claim has been removed; set `Resolved` to the resolution date and summarize the decision in `Notes`.
- Link affected pages to the active row while the claim remains `contested` or `stale`; remove or update those links when the row is resolved.
