# Workflows

> Procedures for ingesting, querying, and maintaining the KB. · Status: confirmed · Last reconciled: 2026-06-29

## Summary

Start from `AGENTS.md`, follow routers to one or two compiled leaves, and answer
from sourced pages. Maintenance work registers sources first, updates the owning
leaf, then fixes nearby links and glossary terms. NRPy code changes keep their
normal project workflow: direct example runs need `PYTHONPATH=.` when there is
no editable install, Python edits require Black before commit, and every
modified Python file gets the single-file static-analysis script unless it is a
generated trusted-value file under `*/tests/*.py`.

## Ingest

1. Register the source in [raw/SOURCES.md](../raw/SOURCES.md) with provenance,
   `frozen` or `living` status, mtime, hash, and ingest state.
2. Decide which branch owns the compiled facts.
3. Update the owning leaf in synthesized prose and cite exact files plus stable
   symbols or headings.
4. From glossary entities touched by the source, follow existing links one hop
   and update affected neighbors.
5. Update a router only when a new leaf or route exists.
6. Run the checks in [lint/CHECKS.md](lint/CHECKS.md) on the touched subgraph.

## Query

Start at [AGENTS.md](../AGENTS.md). Follow the narrow route through branch
routers to the owning leaf, read one or two leaves, then answer from compiled
pages. If a normal question requires broad tree search, repair the relevant
router or add a focused sourced leaf.

## Living-Source Change

Recompute the changed source's mtime and hash, re-ingest affected leaves, and
update one-hop neighbors that share changed glossary entities. Mark a page
`stale` only when it cannot be reconciled immediately.

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
- [coding_style.md](../coding_style.md) - `## Static Analysis Configuration`
- [coding_style.md](../coding_style.md) - `### Trusted Vector File Contract`
- [README.md](../README.md) - `## Contributor Setup`
- [single_file_static_analysis.sh](../.github/single_file_static_analysis.sh) - `run_test_step`
- [main.yml](../.github/workflows/main.yml) - `static-analysis`

## See Also

- [Schema](SCHEMA.md)
- [Lint Checks](lint/CHECKS.md)
- [Contribution Style And Static Analysis](architecture/contribution-style-and-static-analysis.md)
