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
| [Sources](raw/SOURCES.md) | Source manifest with mtimes and hashes. |
| [Source Map](wiki/source-map.md) | Source-to-page dependency seed map and drift follow-up. |
| [Log](wiki/log.md) | Chronological KB maintenance history. |
| [Contradictions](wiki/contradictions.md) | Contested, stale, and reconciled claims. |
| [Syntheses](wiki/syntheses/index.md) | Cross-branch filed syntheses. |
| [Schema](wiki/SCHEMA.md) | Page contracts and governance. |

## Where Do I Start?

| Task | Read first |
| --- | --- |
| Build or run NRPy | [Build And Run](wiki/architecture/build-and-run.md) |
| Understand symbolic-to-C generation | [Symbolic Codegen Lifecycle](wiki/architecture/symbolic-codegen-lifecycle.md) |
| Understand generated outputs | [Generated Output Boundaries](wiki/architecture/generated-output-boundaries.md) |
| Check style and static-analysis rules | [Contribution Style And Static Analysis](wiki/architecture/contribution-style-and-static-analysis.md) |
| Find `CFunction`, `c_codegen`, parameters, or gridfunctions | [Core APIs](wiki/core/index.md) |
| Change equation modules or expression validation | [Equations](wiki/equations/index.md) |
| Work on BHaH generation | [BHaH](wiki/infrastructures/bhah/index.md) |
| Run the first example | [First Wave Equation Run](wiki/examples/first-wave-equation-run.md) |
| Inspect CI coverage | [Generated Project CI](wiki/validation/generated-project-ci.md) |
| Update KB pages | [Workflows](wiki/workflows.md) |

Rules for maintaining this KB live in [wiki/SCHEMA.md](wiki/SCHEMA.md).
