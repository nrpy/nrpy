# Contribution Style And Static Analysis

> Preserve the practical style, artifact, trusted-value, and static-analysis rules for code changes. · Status: confirmed · Last reconciled: 07-20-2026
> Up: [Architecture](index.md)

## Summary

The KB files contributor rules for new or modified code into focused leaves for Python style, C/embedded-C style, equation setup style, infrastructure code style, defensive-guard evidence, and static analysis. Coordinate-admission checks are prohibited unless the current task's user-authored request expressly requests coordinate bounds checking and names the exact target; existing checks remain unchanged and are nonprecedential. API and ABI compatibility shims must not be added or expanded, and interface changes must keep the current commit self-consistent by atomically migrating all affected surfaces and removing the superseded surface. These rules are not claims that every legacy file already conforms. Handwritten Python changes follow [Static Analysis](../validation/static-analysis.md); documentation-only Markdown changes do not invoke Python checks unless Python files are edited.

## Detail

For Python source changes, follow [Static
Analysis](../validation/static-analysis.md), which owns command mechanics,
applicability, classification, and current enforcement gaps.

Trusted numerical value files are a special source class. They live under
`*/tests/`, contain generated reference values, and are regenerated from their
owning module rather than hand-edited. Generated trusted data is exempt from the
single-file script, but its handwritten owner is not. [Test Oracles And Safe
Updates](../validation/test-oracles-and-safe-updates.md) owns store admission
and update safety.

Python conventions are owned by [Python Coding
Style](python-coding-style.md), generated and handwritten C/CUDA/H conventions
by [C And Embedded C Style](c-and-embedded-c-style.md), equation-module rules
by [Equation Setup Style](../equations/equation-setup-style.md), and
backend-specific deltas by [Infrastructure Code
Style](../infrastructures/infrastructure-code-style.md).

Coordinate-admission checks are request-gated. Do not add, expand, strengthen,
duplicate, relocate, or factor logic into a helper when its purpose or effect is
to validate, bound, clamp, reject, abort, skip, return status for, alter control
flow for, select interpolation for, or classify coordinates, logical
coordinates, or coordinate-derived indices, candidates, or stencils by
finiteness, domain, patch or grid-interior membership, representability, or
range. This applies in all handwritten, emitted, and generated languages and to
functions, methods, callers, helpers, wrappers, lambdas, inline functions,
templates, and function-like macros. Authorization requires the current task's
user-authored request to expressly request coordinate bounds checking and name
the exact target. Prior reviews or plans, generic fix, hardening, review, or
undefined-behavior language, theoretical undefined behavior, and inferred need
are insufficient. Existing checks are grandfathered unchanged and are not
precedent.

The prohibition does not cover genuinely non-coordinate storage safety over
already-produced discrete values, such as serialized indices, array or buffer
extents, counts, allocation overflow, storage capacity, or distributed
ownership. Existing interpolation-stencil storage, distributed-ownership, and
boundary-classification checks may remain unchanged, but new or expanded
coordinate-derived stencil-membership checks still require explicit
authorization. This rule overrides ordinary defensive-guard
evidence: realistic caller evidence or a reproducer supports an authorized
check but does not itself authorize one.

Before adding a defensive guard, trace all in-repository callers and identify a
concrete path by which the guarded state can occur. A public untrusted-input
boundary, caller contract, reproduced failure, or realistic runtime path must
establish strong likelihood; a theoretically invalid value or possible
language-level undefined behavior alone does not justify guarding an internal
trusted path. Preserve `*_assume_valid` contracts when callers establish their
preconditions, and fix cross-stage failure propagation at its owning boundary
instead of duplicating conversion or validation logic downstream. A guard test
must exercise that realistic caller path; a synthetic direct call with
otherwise unreachable values is not evidence that a production guard is
needed.

Claim evidence:
- Claim: Coordinate-admission checks require the current task's user-authored request to expressly request bounds checking and name the exact target, regardless of ordinary defensive-guard evidence; existing checks are grandfathered unchanged and nonprecedential, while other defensive guards require a caller-traced realistic failure path and a test of that path.
- Role: normative rule
- Deciding authority: `coding_style.md` - `## Coordinate Bounds-Check Prohibition`, `### Defensive Guard Evidence`
- Corroboration: none available; the contributor guide is the owning authority.

Never add or expand an API or ABI compatibility shim. Prohibited shims include
aliases or re-exports, deprecated forwarders, obsolete generated symbols,
prototypes, or header declarations, dual registrations, old-signature or
old-option translation, and CLI, configuration, or schema bridges whose purpose
or effect is to keep a superseded interface usable.

An interface change must atomically update every affected current-commit definition,
caller, registry entry, generated header and prototype, document, test, and
golden, then remove the superseded surface in the same change. NRPy's only
interface-stability requirement is self-consistency within the current commit.
Existing untouched shims may remain as factual legacy only: they are
nonprecedential, provide no implicit authorization, and must not be copied or
expanded. This policy does not authorize unrelated cleanup, mass API deletion,
or other changes outside the current task.

A minimal internal adapter remains allowed when it connects distinct active
interfaces genuinely required by producers and consumers in the current
commit. It must not expose or translate an obsolete name, signature, option,
symbol, or version.

Claim evidence:
- Claim: API and ABI compatibility shims must never be added or expanded; interface changes must atomically migrate all affected current-commit surfaces and remove the superseded surface, while untouched legacy shims are factual and nonprecedential and distinct active current-commit interfaces may use minimal non-legacy adapters.
- Role: normative rule
- Deciding authority: `coding_style.md` - `## API And ABI Compatibility-Shim Prohibition`
- Corroboration: none available; the contributor guide is the owning authority.

Artifact handling follows [Generated Output
Boundaries](generated-output-boundaries.md): ordinary work does not add binary
or generated products without maintainer approval, and generated trusted values
remain evidence rather than prose.

## Sources

- [coding_style.md](../../coding_style.md) - `#### Embedded C Code String Conventions`, `## C/H Coding Style`, `## Coordinate Bounds-Check Prohibition`, `## API And ABI Compatibility-Shim Prohibition`, `### Defensive Guard Evidence`, `## Python Coding Style`, `### Formatting`
- [raw/source-docs/original-agents.md](../../raw/source-docs/original-agents.md) - historical `## Required Checks`, plus `## Equation Setup Rules` and `## Quick Reference`; current `coding_style.md` decides conflicts

## See Also

- Parent: [Architecture](index.md)
- See also: [Build And Run](build-and-run.md)
- See also: [Python Coding Style](python-coding-style.md)
- See also: [C And Embedded C Style](c-and-embedded-c-style.md)
- See also: [Equation Setup Style](../equations/equation-setup-style.md)
- See also: [Infrastructure Code Style](../infrastructures/infrastructure-code-style.md)
- Depends on: [Code Test Policy](../validation/code-test-policy.md)
- Depends on: [Test Oracles And Safe Updates](../validation/test-oracles-and-safe-updates.md)
- Validated by: [Static Analysis](../validation/static-analysis.md)
- Depends on: [Generated Output Boundaries](generated-output-boundaries.md)
- See also: [Symbolic Codegen Lifecycle](symbolic-codegen-lifecycle.md)
