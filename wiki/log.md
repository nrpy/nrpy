# KB Log

> Append-only chronological memory for durable KB operations.

This file records operations that should compound into future KB state. It is
not a transcript archive or a scratchpad. Keep entries concise, factual, and
parseable with simple Unix tools.

## Entry Headings

Use one of these exact heading forms:

- `## [MM-DD-YYYY] ingest | <title>`
- `## [MM-DD-YYYY] query-filed | <title>`
- `## [MM-DD-YYYY] lint | <scope>`
- `## [MM-DD-YYYY] source-drift | <source>`
- `## [MM-DD-YYYY] reconcile | <page-or-source>`
- `## [MM-DD-YYYY] page-add | <page>`
- `## [MM-DD-YYYY] page-move | <page>`
- `## [MM-DD-YYYY] page-delete | <page>`

## Entry Fields

Each entry uses these fields, in this order:

- `Sources:` source files, source manifests, or external ground truths used.
- `Pages touched:` wiki pages changed or reviewed for the operation.
- `Decision:` durable change, conclusion, or reason no page changed.
- `Checks:` validation commands, lint checks, review checks, or `not run`.
- `Follow-up:` next durable action, or `none`.

## Exclusions

Do not store chat transcripts, scratch logs, token reports, routine lint dumps,
temporary planning notes, command output dumps, or agent-status chatter here.
Only record durable operations whose result should help future agents understand
what changed, why it changed, and what remains to reconcile.

## [06-30-2026] reconcile | compounding-memory support

- `Sources:` [AGENTS.md](../AGENTS.md); [SCHEMA.md](SCHEMA.md); [Workflows](workflows.md); [Sources](../raw/SOURCES.md); Karpathy LLM Wiki approach (`https://gist.githubusercontent.com/karpathy/442a6bf555914893e9891c11519de94f/raw/ac46de1ad27f92b28ac95459c782c07f6b8c964a/llm-wiki.md`).
- `Pages touched:` [KB Log](log.md).
- `Decision:` created append-only chronological memory file for durable KB operations, with parseable headings and fixed entry fields.
- `Checks:` `git diff --check -- wiki/log.md`.
- `Follow-up:` update schema, workflows, catalog, and lint rules in their scoped tasks so future durable operations consistently update this log.

## [06-30-2026] reconcile | governance support pages

- `Sources:` [AGENTS.md](../AGENTS.md); [SCHEMA.md](SCHEMA.md); [Workflows](workflows.md); [Contradictions](contradictions.md); [Sources](../raw/SOURCES.md); Karpathy LLM Wiki approach (`https://gist.githubusercontent.com/karpathy/442a6bf555914893e9891c11519de94f/raw/ac46de1ad27f92b28ac95459c782c07f6b8c964a/llm-wiki.md`).
- `Pages touched:` [AGENTS.md](../AGENTS.md); [SCHEMA.md](SCHEMA.md); [Workflows](workflows.md); [Contradictions](contradictions.md); [Sources](../raw/SOURCES.md).
- `Decision:` added governance routes and support-page rules for log, contradictions, catalog, and source map; defined authority tiers, query filing, external-source trust handling, and an initial contradiction register.
- `Checks:` scoped `git diff --check` commands reported clean by the implementing tasks.
- `Follow-up:` closed by later root routing and source-map reconciliation.

## [06-30-2026] reconcile | global catalog

- `Sources:` [Global Catalog](catalog.md); current `wiki/**/*.md` inventory.
- `Pages touched:` [Global Catalog](catalog.md).
- `Decision:` created global page inventory with page type, route, query aliases, status, last reconciled date, source count, and concept-hub candidate signal.
- `Checks:` catalog task reported `git diff --check -- wiki/catalog.md` clean and row count matched current markdown page count at creation time.
- `Follow-up:` refresh catalog after synthesis pages or new support pages are added.

## [06-30-2026] reconcile | source map

- `Sources:` [Sources](../raw/SOURCES.md); [Source Map](source-map.md); Karpathy LLM Wiki approach as background only.
- `Pages touched:` [Source Map](source-map.md).
- `Decision:` created seed source-dependency map with aggregate rows, exact seed rows, authority tiers, partial-ingest gaps, and next actions for staleness control.
- `Checks:` source-map task reported `git diff --check -- wiki/source-map.md` clean, with 8 aggregate rows and 20 exact-source rows at creation time. Later reconciliation left the source map with 8 aggregate rows and 21 exact-source rows.
- `Follow-up:` add exact dependency rows as high-risk pages change.

## [06-30-2026] lint | compounding-memory checks

- `Sources:` [Lint Checks](lint/CHECKS.md); [kb_lint.py](../tools/kb_lint.py); Karpathy LLM Wiki approach as lint/health-check background.
- `Pages touched:` [Lint Checks](lint/CHECKS.md); [kb_lint.py](../tools/kb_lint.py).
- `Decision:` added mechanical checks for catalog membership, log presence, support-page hygiene, router/leaf contracts, links, local source paths, glossary ownership signals, synthesis routers, and banned wikilinks.
- `Checks:` implementing task reported `git diff --check -- wiki/lint/CHECKS.md`, `git diff --check -- tools/kb_lint.py`, and `python -m py_compile tools/kb_lint.py` clean. Later cleanup reconciled remaining lint failures.
- `Follow-up:` none.

## [06-30-2026] reconcile | typed See Also migration

- `Sources:` [SCHEMA.md](SCHEMA.md); target leaf pages.
- `Pages touched:` [C Codegen](core/c-codegen.md); [Python Codegen](core/python-codegen.md); [C Function Registry](core/c-function-registry.md); [Gridfunctions And Parameters](core/gridfunctions-and-parameters.md); [Reference Metrics](core/reference-metrics.md); [Finite Difference](core/finite-difference.md); [Generated Output Boundaries](architecture/generated-output-boundaries.md); [Generated Project CI](validation/generated-project-ci.md); [Example Generator Catalog](examples/example-generator-catalog.md); [Trusted Expression Pipeline](equations/trusted-expression-pipeline.md); [BHaH Lifecycle](infrastructures/bhah-lifecycle.md); [GR Application Wiring](infrastructures/bhah/gr-application-wiring.md).
- `Decision:` migrated high-traffic `See Also` sections to typed labels such as `Parent`, `Depends on`, `Implements`, `Validated by`, `Example`, `Contrasts with`, and `See also`.
- `Checks:` scoped `git diff --check` commands reported clean by the implementing tasks.
- `Follow-up:` extend typed labels opportunistically when other leaves receive substantial edits.

## [06-30-2026] query-filed | generated backend comparison

- `Sources:` [Generated Backend Comparison](syntheses/generated-backend-comparison.md); [Generated Output Boundaries](architecture/generated-output-boundaries.md); [Generated Project CI](validation/generated-project-ci.md); infrastructure branch pages cited by the synthesis; Karpathy LLM Wiki approach as query-filing background.
- `Pages touched:` [Syntheses](syntheses/index.md); [Generated Backend Comparison](syntheses/generated-backend-comparison.md).
- `Decision:` filed the cross-branch generated-backend comparison as a durable synthesis page and routed it through `wiki/syntheses/index.md`.
- `Checks:` synthesis task reported `git diff --check -- wiki/syntheses/index.md wiki/syntheses/generated-backend-comparison.md` clean and no Obsidian-style wikilinks, `TODO`, or `TBD` matches under `wiki/syntheses`.
- `Follow-up:` keep [Global Catalog](catalog.md) and [Source Map](source-map.md) synchronized when the synthesis page's cited sources change.

## [06-30-2026] lint | final compounding-memory audit

- `Sources:` [Lint Checks](lint/CHECKS.md); [kb_lint.py](../tools/kb_lint.py).
- `Pages touched:` [KB Log](log.md); [kb_lint.py](../tools/kb_lint.py).
- `Decision:` final audit result is recorded in the append-only log, not as a root-level KB markdown artifact.
- `Checks:` `git diff --check`; `rg -n '\[\[' AGENTS.md wiki raw`; `python tools/kb_lint.py`.
- `Follow-up:` remove generated bytecode such as `tools/__pycache__/kb_lint.cpython-312.pyc` before commit.

## [06-30-2026] reconcile | final compounding-memory cleanup

- `Sources:` [AGENTS.md](../AGENTS.md); [Source Map](source-map.md); [kb_lint.py](../tools/kb_lint.py); [Lint Checks](lint/CHECKS.md).
- `Pages touched:` [AGENTS.md](../AGENTS.md); [Source Map](source-map.md); [kb_lint.py](../tools/kb_lint.py); [KB Log](log.md).
- `Decision:` closed earlier follow-ups: root routes now include catalog/source-map, source-map reconciliation left 8 aggregate rows and 21 exact-source rows, and final KB lint passes.
- `Checks:` `git diff --check`; `rg -n '\[\[' AGENTS.md wiki raw`; `python tools/kb_lint.py`.
- `Follow-up:` none from the final mechanical pass.

## [07-02-2026] reconcile | BHaH SO(3) rotation provenance

- `Sources:` [Sources](../raw/SOURCES.md); [Source Map](source-map.md); `nrpy/infrastructures/BHaH/rotation/register_all.py`; `nrpy/infrastructures/BHaH/rotation/*.py`; `nrpy/infrastructures/BHaH/rotation/tests/so3_*__openmp.c`.
- `Pages touched:` [Sources](../raw/SOURCES.md); [Source Map](source-map.md); [KB Log](log.md).
- `Decision:` registered exact raw/source rows for the BHaH rotation package init, bulk registrar, seven SO(3) helper Python modules, and the seven trusted `so3_*__openmp.c` validation files currently read by `register_all.py`; deferred BHaH router, owner leaf, equation handoff, and catalog edits to later Issue 1 phases.
- `Checks:` `git diff --check -- raw/SOURCES.md wiki/source-map.md wiki/log.md`; `python - <<'PY' ...` source-row verification for registered rotation files.
- `Follow-up:` create and route the BHaH SO(3) Rotation Helpers leaf, then update source-map dependent page status from registered to exact seed.

## [07-02-2026] page-add | BHaH SO(3) Rotation Helpers

- `Sources:` [SO(3) Rotation Helpers](infrastructures/bhah/so3-rotation-helpers.md); [Sources](../raw/SOURCES.md); `nrpy/equations/rotation/SO3_rotations.py`; `nrpy/infrastructures/BHaH/rotation/register_all.py`; `nrpy/infrastructures/BHaH/rotation/*.py`; `nrpy/infrastructures/BHaH/rotation/tests/so3_*__openmp.c`.
- `Pages touched:` [SO(3) Rotation Helpers](infrastructures/bhah/so3-rotation-helpers.md); [BHaH](infrastructures/bhah/index.md); [Geometry And Special-Function Support](equations/geometry-and-special-function-support.md); [Global Catalog](catalog.md); [Source Map](source-map.md); [KB Log](log.md).
- `Decision:` added and routed the BHaH SO(3) Rotation Helpers leaf, cited the equation-side `SO3_rotations.py` source, linked the equation-side `SO3Expressions` page to the generated-C helper owner, added the global catalog row, and reconciled the source-map row from registered provenance to exact-seed owner coverage.
- `Checks:` `git diff --check -- wiki/infrastructures/bhah/so3-rotation-helpers.md wiki/infrastructures/bhah/index.md wiki/equations/geometry-and-special-function-support.md wiki/catalog.md raw/SOURCES.md wiki/source-map.md wiki/log.md`; refreshed exact SO(3) and aggregate source rows; `python tools/kb_lint.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache python nrpy/infrastructures/BHaH/rotation/register_all.py`.
- `Follow-up:` none for Issue 1; the two `register_all_unrotate_*__openmp.c` files remain outside this source-map seed row because `register_all.py` does not validate them.

## [07-02-2026] reconcile | BHaH SO(3) tensor helper documentation

- `Sources:` [SO(3) Rotation Helpers](infrastructures/bhah/so3-rotation-helpers.md); [Sources](../raw/SOURCES.md); [Source Map](source-map.md); [Global Catalog](catalog.md); `nrpy/infrastructures/BHaH/rotation/register_all.py`; `nrpy/infrastructures/BHaH/rotation/so3_apply_R_to_tensorDD.py`; `nrpy/infrastructures/BHaH/rotation/tests/so3_apply_R_to_tensorDD_so3_apply_R_to_tensorDD__openmp.c`.
- `Pages touched:` [SO(3) Rotation Helpers](infrastructures/bhah/so3-rotation-helpers.md); [Sources](../raw/SOURCES.md); [Source Map](source-map.md); [Global Catalog](catalog.md); [KB Log](log.md).
- `Decision:` reconciled the BHaH SO(3) helper owner page and exact source rows with the tensor helper now registered by `register_all.py`, raising the documented helper count and trusted generated-C source count from seven to eight.
- `Checks:` `git diff --check`; `python tools/kb_lint.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. python nrpy/infrastructures/BHaH/rotation/so3_apply_R_to_tensorDD.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. python nrpy/infrastructures/BHaH/rotation/register_all.py`.
- `Follow-up:` none.

## [07-02-2026] reconcile | spinning black-hole vector-spin docs and sources

- `Sources:` [Standalone GR/BHaH](examples/standalone-gr-bhah.md); [Sources](../raw/SOURCES.md); [Source Map](source-map.md); [Global Catalog](catalog.md); `nrpy/examples/spinning_blackhole.py`; `nrpy/equations/general_relativity/InitialData_Spherical.py`.
- `Pages touched:` [Standalone GR/BHaH](examples/standalone-gr-bhah.md); [Sources](../raw/SOURCES.md); [Source Map](source-map.md); [Global Catalog](catalog.md); [KB Log](log.md).
- `Decision:` reconciled the standalone GR/BHaH page with public UIUC vector-spin defaults and refreshed stale source-manifest rows for changed equation, infrastructure, and example sources.
- `Checks:` `black nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py nrpy/infrastructures/BHaH/general_relativity/initial_data.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. ./.github/single_file_static_analysis.sh nrpy/infrastructures/BHaH/general_relativity/initial_data.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. ./.github/single_file_static_analysis.sh nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py`; `python tools/kb_lint.py`; `git diff --check`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. python -m nrpy.examples.spinning_blackhole`; `make -C project/spinning_blackhole -j2`.
- `Follow-up:` none.

## [07-02-2026] reconcile | spinning black-hole vector-spin cleanup

- `Sources:` [Standalone GR/BHaH](examples/standalone-gr-bhah.md); [Sources](../raw/SOURCES.md); `nrpy/examples/spinning_blackhole.py`; `nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py`; `nrpy/infrastructures/BHaH/rotation/so3_apply_R_to_tensorDD.py`.
- `Pages touched:` [Standalone GR/BHaH](examples/standalone-gr-bhah.md); [Sources](../raw/SOURCES.md); [KB Log](log.md).
- `Decision:` aligned the standalone page reconciliation date, kept SO(3) helper imports scoped to runtime needs, inlined the single-use Python wrapper for the emitted ADM rotation helper, and constrained the inactive OffsetKerrSchild example branch to signed z-aligned spin.
- `Checks:` `black nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py nrpy/infrastructures/BHaH/rotation/so3_apply_R_to_tensorDD.py nrpy/examples/spinning_blackhole.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. ./.github/single_file_static_analysis.sh nrpy/examples/spinning_blackhole.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. ./.github/single_file_static_analysis.sh nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. ./.github/single_file_static_analysis.sh nrpy/infrastructures/BHaH/rotation/so3_apply_R_to_tensorDD.py`; `python tools/kb_lint.py`; `git diff --check`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. python nrpy/infrastructures/BHaH/rotation/so3_apply_R_to_tensorDD.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. python -m nrpy.examples.spinning_blackhole`; `make -C project/spinning_blackhole -j2`.
- `Follow-up:` none.

## [07-02-2026] reconcile | spin-vector style and generated comments

- `Sources:` [Sources](../raw/SOURCES.md); `nrpy/infrastructures/BHaH/general_relativity/initial_data.py`; `nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py`; `nrpy/infrastructures/BHaH/general_relativity/spin_vector.py`.
- `Pages touched:` [Sources](../raw/SOURCES.md); [KB Log](log.md).
- `Decision:` clarified spin-aligned sampling comments and ValueError documentation, normalized touched docstring return style, matched C-function registration variable structure, removed empty checkpoint-placeholder output for non-spin restarts, and recorded the aligned-frame polar-axis caveat without adding coordinate clamping.
- `Checks:` `black nrpy/infrastructures/BHaH/general_relativity/initial_data.py nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py nrpy/infrastructures/BHaH/general_relativity/spin_vector.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. ./.github/single_file_static_analysis.sh nrpy/infrastructures/BHaH/general_relativity/initial_data.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. ./.github/single_file_static_analysis.sh nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. ./.github/single_file_static_analysis.sh nrpy/infrastructures/BHaH/general_relativity/spin_vector.py`; `python tools/kb_lint.py`; `git diff --check`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. python -m nrpy.examples.spinning_blackhole`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. python - <<'PY' ...` checkpoint placeholder check; `make -C project/spinning_blackhole -j2`.
- `Follow-up:` none.

## [07-02-2026] reconcile | spin-vector registration comment cleanup

- `Sources:` [Sources](../raw/SOURCES.md); `nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py`; `nrpy/infrastructures/BHaH/general_relativity/spin_vector.py`.
- `Pages touched:` [Sources](../raw/SOURCES.md); [KB Log](log.md).
- `Decision:` removed trailing whitespace from the generated validation helper description, ordered C-function registration variables before the emitted body, and clarified the Spherical tilted-spin coordinate-singularity caveat.
- `Checks:` `black nrpy/infrastructures/BHaH/general_relativity/spin_vector.py nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. ./.github/single_file_static_analysis.sh nrpy/infrastructures/BHaH/general_relativity/spin_vector.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. ./.github/single_file_static_analysis.sh nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py`; `python tools/kb_lint.py`; `git diff --check`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. python -m nrpy.examples.spinning_blackhole`; `python - <<'PY' ...` generated validation-helper header check; `make -C project/spinning_blackhole -j2`.
- `Follow-up:` none.

## [07-03-2026] reconcile | inline UIUC spin-vector registration into initial_data

- `Sources:` [Sources](../raw/SOURCES.md); `nrpy/infrastructures/BHaH/general_relativity/initial_data.py`; `nrpy/infrastructures/BHaH/general_relativity/__init__.py`.
- `Pages touched:` [Sources](../raw/SOURCES.md); [KB Log](log.md).
- `Decision:` deleted the single-consumer `spin_vector.py` module and inlined its spin-vector CodeParameter registration and UIUC validation-helper C-function registration into the spin branch of `register_CFunction_initial_data`, so the entire UIUC spin flow reads top-to-bottom in one file, the CodeParameters are registered by the routine that owns the generated code consuming them, and the package no longer exports a two-step registration API that is only safe in one call order; generated project output is unchanged except CodeParameter module-attribution comments now cite `initial_data`.
- `Checks:` `black nrpy/infrastructures/BHaH/general_relativity/initial_data.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. ./.github/single_file_static_analysis.sh nrpy/infrastructures/BHaH/general_relativity/initial_data.py`; `python tools/kb_lint.py`; `git diff --check`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. python -m nrpy.examples.spinning_blackhole`; before/after diff of all generated `project/spinning_blackhole` `.c`/`.h`/`.par` files (attribution comments only); `make -C project/spinning_blackhole -j4`.
- `Follow-up:` none.

## [07-05-2026] reconcile | BHaH MoL diagonal-RK3 hooks

- `Sources:` [MoL Time Integration](infrastructures/bhah/mol-time-integration.md); [Sources](../raw/SOURCES.md); [Source Map](source-map.md); `nrpy/infrastructures/BHaH/MoLtimestepping/MoL_step_forward_in_time.py`.
- `Pages touched:` [MoL Time Integration](infrastructures/bhah/mol-time-integration.md); [Sources](../raw/SOURCES.md); [Source Map](source-map.md); [KB Log](log.md).
- `Decision:` reconciled BHaH MoL documentation and infrastructure aggregate source metadata with the diagonal-RK3 k2 post-hook fix; kept compact-binary example documentation on RK4 because `two_blackholes_collide.py` defaults to RK4.
- `Checks:` `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. python nrpy/infrastructures/BHaH/MoLtimestepping/MoL_step_forward_in_time.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. ./.github/single_file_static_analysis.sh nrpy/infrastructures/BHaH/MoLtimestepping/MoL_step_forward_in_time.py`; `XDG_CACHE_HOME=/tmp/nrpy-cache PYTHONPATH=. python nrpy/infrastructures/BHaH/MoLtimestepping/register_all.py`; `python tools/kb_lint.py`; `git diff --check`.
- `Follow-up:` none.

## [07-05-2026] reconcile | remove source-tracking hash and timestamp metadata

- `Sources:` [AGENTS.md](../AGENTS.md); [Schema](SCHEMA.md); [Workflows](workflows.md); [Sources](../raw/SOURCES.md); [kb-instructions.md](../raw/source-docs/kb-instructions.md); [kb_lint.py](../tools/kb_lint.py).
- `Pages touched:` [AGENTS.md](../AGENTS.md); [Schema](SCHEMA.md); [Workflows](workflows.md); [Lint Checks](lint/CHECKS.md); [Source Map](source-map.md); [Sources](../raw/SOURCES.md); [kb-instructions.md](../raw/source-docs/kb-instructions.md); [KB Log](log.md); every `wiki/**/*.md` page for date normalization, including [Global Catalog](catalog.md) and [Contradictions](contradictions.md).
- `Decision:` source tracking no longer stores or computes hash or `mtime` metadata: sources are never hashed for tracking with any algorithm, and manifest rows are source identity, provenance, status, and ingest state only, with no digest values, no `Mtime`/`Hash` columns, and no `mtime` values; a root source-tracking and date policy section was added to `AGENTS.md`; all retained KB dates now use `MM-DD-YYYY`; source drift moved from stored-fingerprint comparison to dependency-aware review; and the frozen `kb-instructions.md` scaffold was narrowly edited, as an explicit logged exception to frozen-source preservation, so its superseded source-tracking instructions cannot keep the removed behavior alive.
- `Checks:` `git status --short`; `git diff --check`; `python -m py_compile tools/kb_lint.py`; `python tools/kb_lint.py`; `find wiki raw -type f ! -name '*.md' -print`; `rg -n '\[\[' AGENTS.md wiki raw`; targeted `rg` scans over `AGENTS.md`, `wiki`, `raw`, and `tools` for removed source-tracking metadata (digest values, manifest timestamp/digest columns), `YYYY-MM-DD` date literals, and remaining technical `hash` facts.
- `Follow-up:` none.

## [07-05-2026] reconcile | metadata lint allowlist and manifest contract

- `Sources:` [kb_lint.py](../tools/kb_lint.py); [Schema](SCHEMA.md); [Workflows](workflows.md); [Sources](../raw/SOURCES.md); code-review findings on the metadata migration.
- `Pages touched:` [Schema](SCHEMA.md); [Workflows](workflows.md); [Sources](../raw/SOURCES.md); [kb_lint.py](../tools/kb_lint.py); [KB Log](log.md).
- `Decision:` tightened the removed-metadata lint exemption so a mention passes only when a prohibition or supersession keyword precedes it with no sentence boundary in between, making positive restore-the-metadata instructions fail; and clarified the manifest contract: repository rows carry identity, provenance, status, and ingest state as their only source-tracking fields, exact cited-file rows abbreviate to source and status with provenance implied by the repository path and ingest state tracked through covering aggregate rows and the source map, and external-source accessed dates and notes are ordinary KB metadata.
- `Checks:` `python -m py_compile tools/kb_lint.py`; `python tools/kb_lint.py`; unit exercise of the exemption helper against positive-instruction and prohibition sample lines; `git diff --check`.
- `Follow-up:` none.

## [07-06-2026] page-add | coding style leaves

- `Sources:` [Schema](SCHEMA.md); [Workflows](workflows.md); [Source Map](source-map.md); [Global Catalog](catalog.md); project style-guide source used during initial filing.
- `Pages touched:` [Python Coding Style](architecture/python-coding-style.md); [C And Embedded C Style](architecture/c-and-embedded-c-style.md); [Equation Setup Style](equations/equation-setup-style.md); [Infrastructure Code Style](infrastructures/infrastructure-code-style.md); [Contribution Style And Static Analysis](architecture/contribution-style-and-static-analysis.md); [Architecture](architecture/index.md); [Equations](equations/index.md); [Infrastructures](infrastructures/index.md); [Global Catalog](catalog.md); [Source Map](source-map.md); [KB Log](log.md).
- `Decision:` filed the project style-guide contents into focused KB leaves for Python, C/embedded C, equation setup, and infrastructure code style; kept the contribution page as the quick route; routed new leaves through owning branch indexes; and expanded source dependency routing for affected style pages.
- `Checks:` `git diff --check`; `python tools/kb_lint.py`.
- `Follow-up:` none.

## [07-06-2026] reconcile | remove root style-guide references

- `Sources:` [Schema](SCHEMA.md); [Workflows](workflows.md); [Sources](../raw/SOURCES.md); [Source Map](source-map.md); [Global Catalog](catalog.md); [original-agents.md](../raw/source-docs/original-agents.md).
- `Pages touched:` [Sources](../raw/SOURCES.md); [Source Map](source-map.md); [Global Catalog](catalog.md); [Contribution Style And Static Analysis](architecture/contribution-style-and-static-analysis.md); [Generated Output Boundaries](architecture/generated-output-boundaries.md); [Python Coding Style](architecture/python-coding-style.md); [C And Embedded C Style](architecture/c-and-embedded-c-style.md); [Equation Setup Style](equations/equation-setup-style.md); [BSSN Family](equations/general-relativity/bssn-family.md); [Trusted Expression Pipeline](equations/trusted-expression-pipeline.md); [Infrastructure Code Style](infrastructures/infrastructure-code-style.md); [Static Analysis](validation/static-analysis.md); [Workflows](workflows.md); [Lint Checks](lint/CHECKS.md); [original-agents.md](../raw/source-docs/original-agents.md); [KB Log](log.md).
- `Decision:` removed all AGENTS/wiki/raw references to the deprecated root style-guide filename, replaced style citations with the preserved agent-rule source and direct tool/config sources, removed the exact source-manifest row for the root style-guide file, and refreshed catalog/source-map metadata for the affected pages.
- `Checks:` targeted filename scan over `AGENTS.md`, `wiki`, and `raw` returned no matches; `python tools/kb_lint.py`; `git diff --check`.
- `Follow-up:` none.
