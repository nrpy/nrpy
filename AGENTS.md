# NRPy Knowledge Base

This repository carries a plain-markdown knowledge base for NRPy. Start here,
follow a few links, and synthesize from the compiled pages instead of grepping
the whole tree first.

## Router

| Go to | Use it for |
| --- | --- |
| [Architecture](wiki/architecture/index.md) | Project purpose, build/run paths, generated-output boundaries, and contribution rules. |
| [Core APIs](wiki/core/index.md) | Core codegen APIs, parameters, gridfunctions, indexed expressions, reference metrics, and finite difference support. |
| [Equations](wiki/equations/index.md) | Symbolic equation families, GR and SEOBNR routes, support helpers, and trusted expression validation. |
| [Infrastructures](wiki/infrastructures/index.md) | Generated-backend lifecycle and infrastructure routes for BHaH, ETLegacy, CarpetX, superB, and JAX. |
| [Examples](wiki/examples/index.md) | First wave-equation run and black-hole evolution examples. |
| [Validation](wiki/validation/index.md) | Static analysis and generated-project CI. |
| [Glossary](wiki/glossary.md) | Canonical terms. |
| [Catalog](wiki/catalog.md) | Global page inventory and query-routing terms. |
| [Workflows](wiki/workflows.md) | KB ingest, query, and maintenance procedures. |
| [Lint Checks](wiki/lint/CHECKS.md) | Mechanical and review checks for the KB. |
| [Sources](raw/SOURCES.md) | Source manifest with provenance, status, and ingest state. |
| [Source Map](wiki/source-map.md) | Source-to-page dependency seed map and drift follow-up. |
| [Contradictions](wiki/contradictions.md) | Contested, stale, and reconciled claims. |
| [Syntheses](wiki/syntheses/index.md) | Cross-branch filed syntheses. |
| [Schema](wiki/SCHEMA.md) | Page contracts and governance. |

## Where Do I Start?

| Task | Read first |
| --- | --- |
| Build, run, inspect code generation, generated outputs, or contribution rules | [Architecture](wiki/architecture/index.md) |
| Find `CFunction`, `c_codegen`, parameters, or gridfunctions | [Core APIs](wiki/core/index.md) |
| Change equation modules or expression validation | [Equations](wiki/equations/index.md) |
| Work on BHaH, ETLegacy, CarpetX, superB, or JAX generation | [Infrastructures](wiki/infrastructures/index.md) |
| Run or compare example generators | [Examples](wiki/examples/index.md) |
| Inspect static analysis, expression validation, or generated-project CI | [Validation](wiki/validation/index.md) |
| Update KB pages | [Workflows](wiki/workflows.md) |

## Source-Tracking And Date Policy

These rules bind every KB manifest and doc under `AGENTS.md`, `wiki/`, and
`raw/`:

- No source-tracking hash columns or values of any kind — `sha256` or any
  other digest — and no hashing of sources at all.
- No `mtime` columns or values.
- Retained KB dates use `MM-DD-YYYY`.
- Do not output KB maintenance notes to a separate log file. This KB already
  lives in a git repo: commit history records durable operations, so separate
  logs are redundant and wasteful.

Source drift is handled by dependency-aware review of changed paths, source
status, [Source Map](wiki/source-map.md) rows, and affected compiled pages —
not by stored fingerprints.

Rules for maintaining this KB live in [wiki/SCHEMA.md](wiki/SCHEMA.md).

## KB Checker And Coordination Scope

`python tools/kb_lint.py` is the canonical deterministic KB structure check.
`--all` is a compatibility alias with identical coverage, not a stronger mode.

Commissioned `plan1.md` through `plan5.md`, `plan_synth.md`, and `tasks1.md`
through `tasks6.md` are coordination artifacts. They may remain untracked and
are exempt from KB routing/catalog checks. Never file, stage, move, or delete
them as KB content unless the user directs that action.
