---
name: dialectic
description: Use when explicitly asked for NRPy dialectic or paired independent review, or required by repository policy. Do not activate merely to mention, audit, or edit this skill. Two independent seats; delta-and-impact follow-ups; verified /work/ delivery.
---

# NRPy Dialectic

Read [the required shared protocol](../../review-protocol.md) for modes, independence,
coverage, budgets, validation, and delivery. Load this skill and the protocol only;
distribute them together.

Use exactly two independent seats. Root owns live writes; seats never delegate
or edit live targets. Default to `review`, including documents and all KB work.
Use `design` or `dual-build` only under the shared criteria and limits.

## Roles

**Scientific-software verifier.** Check mathematics, executable behavior,
assumptions, indices, invariants, interfaces, failure cases, numerical/oracle
adequacy, and affected generated/backend contracts.

**Integration and simplification reviewer.** Check NRPy conventions, reusable
components, dependency direction, public interfaces, source/generated ownership,
docs/tests, and stale references. Prefer the simplest sufficient implementation;
reject unsupported abstractions and compatibility machinery.

For non-software work, use a first-principles analyst and a domain/implementation expert.
Roles add complementary scrutiny, not duplicate test suites.

## Completion

Follow-ups cover the delta, open findings, and code/docs impacts only. Root must
[deliver and verify](../../review-protocol.md#delivery) the selected result at ALL
intended authorized `/work/` paths; a verdict or `DRAFT COMPLETE` is not delivery.
Preserve blocked drafts separately, never as accepted changes.

If the protocol is missing, do not invent approval. Preserve produced drafts in a
fresh `/work/review-results/<task>/unapproved/` area, verify copies, and report the
missing file. Explicit user write/output restrictions control.
