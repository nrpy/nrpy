---
name: trialectic
description: Use when explicitly asked for NRPy trialectic, "engage tri", "run tri", "use tri", or three-agent independent review, or required by repository policy. Do not activate merely to mention, audit, or edit this skill. Three independent seats; delta-and-impact follow-ups; verified /work/ delivery.
---

# NRPy Trialectic

Read [the required shared protocol](../../review-protocol.md) for modes, independence,
coverage, budgets, validation, and delivery. Load this skill and the protocol only;
distribute them together.

Use exactly three independent seats. Root owns live writes; seats never delegate
or edit live targets. Default to `review`, including documents and all KB work.
Use `design` or `tri-build` only under the shared criteria and limits.

## Roles

**Scientific-software verifier.** Check mathematics, executable behavior,
assumptions, indices, invariants, interfaces, failure cases, numerical/oracle
adequacy, and affected generated/backend contracts.

**Integration and simplification reviewer.** Check NRPy conventions, reusable
components, dependency direction, public interfaces, source/generated ownership,
docs/tests, and stale references. Prefer the simplest sufficient implementation;
reject unsupported abstractions and compatibility machinery.

**Standards and release auditor.** Trace criteria to authority/evidence. Check
omissions in validation, packaging, API/source/generated boundaries, docs, affected
backends, and delivery; audit delta effects, not unrelated repository areas.

For non-software work, use a first-principles analyst, a domain/implementation expert, and an adversarial editor/verifier.
Roles add complementary scrutiny, not duplicate test suites.

## Completion

Follow-ups cover the delta, open findings, and code/docs impacts only. Root must
[deliver and verify](../../review-protocol.md#delivery) the selected result at ALL
intended authorized `/work/` paths; a verdict or `DRAFT COMPLETE` is not delivery.
Preserve blocked drafts separately, never as accepted changes.

If the protocol is missing, do not invent approval. Preserve produced drafts in a
fresh `/work/review-results/<task>/unapproved/` area, verify copies, and report the
missing file. Explicit user write/output restrictions control.
