# Workflows

> Procedures for querying, ingesting, reconciling, and checking the NRPy KB. · Status: confirmed

## Summary

Start at [AGENTS.md](../AGENTS.md), jump directly to the relevant branch or sub-branch `index.md` router, and read only the needed leaf or `###` subsection. State each repository fact in one defining leaf, compute reverse dependencies dynamically with `rg -l -F`, and verify changes with `tools/kb_lint.py` and [Language For Computational Physicists](computational-physics-language.md).

## Search Order

1. Read [AGENTS.md](../AGENTS.md) and follow the direct branch or sub-router link to the target `index.md` (or top-level leaf).
2. For ambiguous or cross-branch topics, run `rg -n -i -F '<term>' wiki/catalog.md` and read only matching rows.
3. Read compact leaves (`< 8,000` B) directly in one `read_file` call. On leaves `>= 8,000` B, run `rg -n '^### ' <page>` first and read only the needed `###` line range. Follow a `## See Also` link only when its relation prefix (`Depends on:`, `Implements:`, `Validated by:`, `Example:`, `Contrasts with:`, `See also:`) identifies a missing prerequisite or check.
4. Search support tables selectively without reading whole files:
   - Reverse file citations: `rg -l -F '<path>' wiki/`
   - Exact source or aggregate keys: ``rg -n -F '| `<key>` |' raw/SOURCES.md`` or ``rg -n -F '| `<aggregate>` |' wiki/source-map.md``
   - Contradiction entries: `rg -n -F '| <CONTR-ID> |' wiki/contradictions.md` plus the matching `### CONTR-XXXX` subsection
   - Canonical term spellings: `rg -n -i -F '<term>' wiki/glossary.md`
5. Search `nrpy/` source files directly when verifying or updating code-backed claims.

## Ingest

1. Confirm the source is covered in [raw/SOURCES.md](../raw/SOURCES.md) via an aggregate row, an exact local row (`Source | Status`), or an external row (`Source | Tier | Provenance | Status | Ingest | Notes`, `Notes <= 80` B). Never store hashes, dates, timestamps, `mtime` values, file/page/job counts, or run results.
2. Select the single defining leaf for the topic under the appropriate branch router.
3. Synthesize durable facts in direct computational-physics prose, add `<= 300` B `Evidence (<role>):` lines using shortest-unambiguous `` `subpath.py::symbol` `` tags for high-risk claims, and list each cited file or URL once in `## Sources`.
4. Partition `## Detail` with descriptive `###` headings whenever `## Detail` exceeds `4,000` B (`<= 6,000` B per chunk, `<= 20,000` B per leaf) and keep `1` to `6` prefixed non-parent `## See Also` bullets.
5. Update the immediate parent `index.md` and [catalog.md](catalog.md) only when adding, moving, renaming, or deleting a leaf; update [source-map.md](source-map.md) only when an aggregate's scope changes.

## Query Filing

Do not file one-off answers. File durable syntheses only when they resolve recurring cross-module questions, record non-obvious numerical or build constraints, or fix missing router paths:

1. Single-topic updates go into the existing defining leaf (splitting into two focused leaves if the leaf would exceed `20,000` B).
2. Narrow new topics add a leaf under the branch router (`2` to `14` children per `index.md`).
3. Cross-branch comparisons with a single leaf sit directly under `wiki/` (`> Up: [NRPy Knowledge Base](../AGENTS.md)`, such as [Generated Backend Comparison](generated-backend-comparison.md)); create `wiki/syntheses/index.md` only when `>= 2` cross-branch synthesis leaves exist.

## Living-Source Change

1. Inspect the changed code or configuration diff and verify its [raw/SOURCES.md](../raw/SOURCES.md) coverage.
2. Locate citing wiki pages dynamically with `rg -l -F '<changed_path>' wiki/`.
3. Run `rg -l -F '<changed_symbol>' wiki/` across other wiki pages **only** when a public function signature, `CodeParameter`, gridfunction group, or mathematical convention changed—never expand 1-hop across generic [glossary.md](glossary.md) terms.
4. Reconcile the defining leaf; if a neighboring page duplicates the changed claim instead of linking to the defining leaf, replace the duplicate prose with a 1-sentence link. Mark a leaf `stale` (with a `CONTR-XXXX` entry) only when immediate reconciliation is blocked.

## External Source And Contradiction Updates

Classify external URLs in [raw/SOURCES.md](../raw/SOURCES.md) as `background` unless a durable KB claim depends on a version-pinned paper, tag, or frozen excerpt (`external-spec`). Repository code, tests, generated evidence, and repository docs still take precedence per [SCHEMA.md](SCHEMA.md#citation-rules).

When sources conflict or living-source drift cannot be reconciled in the same change, add a 4-column row (`<= 220` B) and `### CONTR-XXXX` subsection in [contradictions.md](contradictions.md) and place `Claim status: contested; contradiction: CONTR-0001.` (or `stale`) with a `#contr-0001` link on every affected leaf. When resolved and verified by its resolution check, reconcile affected leaves, remove all `#contr-*` links, delete the row and subsection, and append the ID to `Retired IDs:`.

## Safe Reproduction

In the shared `/work/` checkout, run only side-effect-free scoped checks (`PYTHONDONTWRITEBYTECODE=1 python3 -B tools/kb_lint.py`, `git diff --check`). Run code generators, C/CUDA/Fortran/Charm++ builds, test-oracle updates, and static-analysis wrappers (`Black`, `isort`, `mypy`, `Pylint`, `pydocstyle`, `darglint`) inside an isolated copy with redirected cache and output directories (`HOME`, `XDG_CACHE_HOME`, `PYLINTHOME`, `MYPY_CACHE_DIR`). Inspect scripts for hardcoded paths, deletions, network calls, or package installs before execution; never clean or overwrite shared `project/` outputs or modify `.github/workflows/main.yml` without explicit user authorization.

## Evidence-First Debugging

Follow the runtime failure inspection rule in [AGENTS.md](../AGENTS.md#agent-working-rules), then record only the durable mathematical or implementation lesson in `wiki/`.

## Deterministic And Review Checks

Run before completing KB or Python edits:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -B tools/kb_lint.py
git diff --check
```

`tools/kb_lint.py` must exit `0` and print `KB lint passed.` For modified handwritten Python files (including `tools/kb_lint.py`; generated `*/tests/*.py` trusted dictionaries are exempt, though their handwritten callers are not), run `.github/single_file_static_analysis.sh <file.py>` in an isolated copy per [Static Analysis](validation/static-analysis.md), [Code Test Policy](validation/code-test-policy.md), and [Test Oracles And Safe Updates](validation/test-oracles-and-safe-updates.md).

Semantic review checklist:
- `git status --short` shows only intended changes, with `raw/source-docs/*.md` untouched and no generated `project/`, cache, or scratch report files added.
- Every added or revised passage satisfies [KB Rules](SCHEMA.md) and [Language For Computational Physicists](computational-physics-language.md) with zero vague metaphors, no prose restating adjacent table cells, and no template/to-do markers.
- Cited symbols, headings, and equations exist in the cited files, and each fact lives on a single defining leaf.

## Sources

- [AGENTS.md](../AGENTS.md) - root navigation, command conventions, runtime debugging rule, and protected workflow scope
- [tools/kb_lint.py](../tools/kb_lint.py) - KB link, reachability, and Evidence-citation checks
- [coding_style.md](../coding_style.md) - `## Python Coding Style`, `## Static Analysis Configuration`
- [README.md](../README.md) - `## Contributor Setup`

## See Also

- Depends on: [KB Rules](SCHEMA.md)
- Depends on: [Language For Computational Physicists](computational-physics-language.md)
- Validated by: [Static Analysis](validation/static-analysis.md)
- See also: [Code Test Policy](validation/code-test-policy.md)
- See also: [Test Oracles And Safe Updates](validation/test-oracles-and-safe-updates.md)
