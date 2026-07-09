# Contradictions

> Register for contested and stale KB claims. Plain Markdown only. · Status: confirmed · Last reconciled: 07-06-2026

Known contested/stale claims as of 07-06-2026 are tracked below.

## Register

| Claim | Status | Source A | Source B | Authority decision | Affected pages | Opened | Resolved | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `manga_bhah_lib` final print suggests running `./bhah_lib`, but the generated target is a library build. | stale | `nrpy/examples/manga_bhah_lib.py` final print | [Matter TOV Workflows](examples/matter-tov-workflows.md) source/Makefile interpretation | Page treats print as stale generic messaging until generator text is reconciled. | [Matter TOV Workflows](examples/matter-tov-workflows.md) | 06-30-2026 | - | Recheck when example generator output text changes. |
| `sebobv1_jax` emits a `Commondata(...)` return that includes `a_f`, but its batch commondata registration supplies fourteen names and only thirteen dtypes/defaults, so `zip()` truncation drops `a_f` from the generated dataclass. | contested | `nrpy/infrastructures/JAX/sebob/SEOBNRv5_aligned_spin_coefficients.py` emitted return | `nrpy/examples/sebobv1_jax.py` `register_commondata` batch call via `nrpy/infrastructures/JAX/commondata.py` `zip()` | Pages document the mismatch as a known upstream inconsistency until `a_f` is registered or the return signature changes. | [SEBOBv1 JAX Workflow](infrastructures/jax/sebobv1-jax-workflow.md); [Commondata And PyFunction Registry](infrastructures/jax/commondata-and-pyfunction-registry.md) | 07-06-2026 | - | `register_commondata()` does not length-check its four input lists; truncation count defines registered fields. |

## Rules

- Open a row when registered sources disagree, when a living source supersedes a KB claim, or when a page is marked `contested` or `stale`.
- Update a row when new evidence, source authority, affected pages, or reconciliation status changes.
- Resolve a row only after affected pages have been reconciled or the claim has been removed; set `Resolved` to the resolution date and summarize the decision in `Notes`.
- Link affected pages to the active row while the claim remains `contested` or `stale`; remove or update those links when the row is resolved.
