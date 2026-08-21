---
name: dialectic
description: Use two fresh-context independent agents to create or improve non-KB files, cross-critique complete drafts, synthesize and validate one candidate, and install it only after unanimous approval. Invoke for “engage dialectic,” “use dialectic,” “run dialectic,” two-agent consensus editing, or paired independent review. In NRPy, never use it for `AGENTS.md`, `wiki/**`, or `raw/**`.
---

# Dialectic

Produce two independent complete drafts, cross-critique them, synthesize the simplest sufficient candidate, validate it, and change live targets only after both agents accept the same candidate and evidence.

This `SKILL.md` is the entire Dialectic implementation. Never require or add Dialectic-owned companion scripts, tests, metadata, templates, state files, or other skill files. Use existing project tools and task-relevant tests when warranted.

## Freeze scope

Before drafting, give both agents one brief containing the exact request, acceptance criteria, target paths and create/revise/delete intent, initial target contents or absence, exact authorized implementation inputs, applicable ambient repository rules, and proportionate validation commands. Authorize individual paths: no globs, recursive discovery, neighboring-file reads, or writes outside targets. Preserve unrelated and concurrent work.

Never target, copy, inventory, or hash `AGENTS.md`, `wiki/**`, or `raw/**`; use ordinary NRPy KB workflow for those paths. Treat harness-supplied repository instructions as ambient policy, not Dialectic inputs, and obey all protected-file rules.

Use isolated temporary work areas outside live target tree. If request, scope, inputs, assumptions, target contents, or validation changes, discard current work and restart two fresh drafts from one revised brief.

## Choose roles

Treat work as software engineering when it creates, changes, reviews, tests, builds, configures, or documents executable software, code-facing interfaces, or developer tooling. For such work, both agents must receive the same exact authorized inputs and applicable standards. Both must also be NRPy coding-style experts: enforce architecture and module conventions, public interfaces, symbolic and code-generation idioms, generated-output ownership, infrastructure boundaries, validation and oracle policy, static-analysis expectations, compatibility, and the simplest-sufficient implementation rule in addition to their primary characteristics.

1. **Software engineer:** Check correctness, architecture, interfaces, failure cases, maintainability, tests, integration, and executable behavior.
2. **Code-simplification agent:** Aggressively reject over-engineering, speculative generality, duplication, needless indirection, and unjustified layers, wrappers, helpers, dependencies, configuration, or abstractions. Inspect every relevant exact authorized NRPy implementation input for established functions, classes, APIs, patterns, and components. Name what was considered; reuse a sufficient established component instead of reinventing it. Allow new machinery only when one concrete unmet requirement proves reuse insufficient. If evidence needs another path, identify that exact path and restart with authorization; never browse beyond scope or use missing access to justify reinvention.

For other work, use one first-principles analyst and one domain and implementation expert. Both retain NRPy coding-style expertise and apply it whenever NRPy conventions govern the artifact.

## Draft and cross-critique

Spawn exactly two fresh-context agents concurrently with separate work areas. Give them identical briefs plus their distinct roles and shared NRPy coding-style responsibilities. During independence, permit only exact authorized reads and each agent's own temporary writes; forbid peer reads and live writes. Require each complete draft to include all requested target operations and a short report with rationale, assumptions, inventory, role evidence, NRPy standards evidence, and final exact line `DRAFT COMPLETE`. Draft authors never approve their own draft.

After both drafts complete, freeze them and verify completeness and scope before peer access. Do not run expensive full validation on every draft unless needed to make it reviewable. Then have each agent critique the other agent's draft. Every material finding must name the exact location, violated requirement or standard, evidence, consequence, and smallest correction. Reviewers may expose defects but may not enlarge acceptance criteria.

## Synthesize and validate

Create one candidate outside live tree by resolving every material critique. Prefer the smallest candidate that meets all requirements and reuses established NRPy components. Represent requested deletion by absence. Confirm deleted resources have no stale references.

Validate exact candidate with existing checks appropriate to its format and risk. Inspect the diff for target confinement and unrelated changes; run applicable syntax, format, lint, build, focused behavior, documentation, or structural checks. Record commands, outcomes, manual limits, established components reused, and why any new abstraction is necessary. Any failed required check or unresolved material critique blocks approval.

Treat review and validation evidence as provisional until root mechanically verifies every declared candidate/input/validation binding against the exact brief, candidate bytes, and recorded check output. Before this preflight succeeds and evidence freezes, root may correct only mechanical binding or transcription defects; review reasoning, decisions, and conclusions cannot change under this exception.

After preflight succeeds, freeze candidate bytes and all review and validation evidence. Any later change to candidate bytes or evidence creates a new candidate round and requires fresh validation and two fresh decisions.

## Require unanimity and install

Give both agents the same frozen candidate, critique resolution, and validation evidence. Each independently returns either acceptance explicitly bound to those bytes and evidence or one concrete blocking defect. Silence, majority, stale review, or conditional approval does not pass.

Immediately before installation, confirm live targets and authorized inputs still match the shared brief. Drift restarts the workflow. After both acceptances, apply only approved writes and deletions, verify installed state equals the accepted candidate, and rerun applicable installed-tree checks. Report changed and deleted targets, validation results and limits, both decisions, and temporary run location.
