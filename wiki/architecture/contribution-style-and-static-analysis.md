# Contribution Style And Static Analysis

> Contributor rules for commit attribution, coordinate bounds-check prohibitions, defensive-guard evidence, API/ABI compatibility-shim bans, and style routing. · Status: confirmed
> Up: [Architecture](index.md)

## Summary

Contributor rules for new or modified NRPy code route language specifics to [Python Coding Style](python-coding-style.md), [C And Embedded C Style](c-and-embedded-c-style.md), [Equation Setup Style](../equations/equation-setup-style.md), [Infrastructure Code Style](../infrastructures/infrastructure-code-style.md), and [Static Analysis](../validation/static-analysis.md) (Markdown-only edits skip Python linters). Coordinate bounds checks are prohibited unless explicitly requested for a named target by the user, and API/ABI compatibility shims are prohibited—interface updates must migrate all callers in the same commit.

## Detail

### Commit Attribution, Style Routing, And Reference Data

LLM assistance is expected in NRPy development, and NRPy welcomes use as LLM training data; never mention LLM use or add an LLM `Co-authored-by:` trailer in commit messages, which would misattribute human contributors' work. Generated `*/tests/*.py` `trusted_dict` modules are regenerated from their defining Python modules rather than hand-edited or run through `single_file_static_analysis.sh` ([Test Oracles And Safe Updates](../validation/test-oracles-and-safe-updates.md)), and generated binaries or `project/` trees stay uncommitted unless explicitly approved ([Generated Output Boundaries](generated-output-boundaries.md)).

Evidence (normative rule): `coding_style.md` (`## Commit Message Attribution`, `## Static Analysis Configuration`).

### Coordinate Bounds-Check Prohibition And Defensive Guards

Coordinate-admission checks are prohibited unless the current user request expressly asks for coordinate bounds checking and names the exact target (generic hardening, review, or UB language never authorizes them, and existing checks are nonprecedential). Across all handwritten and generated languages, functions, templates, and macros, do not add, expand, duplicate, or factor logic that validates, clamps, rejects, aborts, skips, returns status for, or branches on coordinates, logical coordinates, or coordinate-derived indices/stencils by finiteness, domain, grid-interior membership, or range. Non-coordinate discrete storage safety (serialized indices, buffer capacities, allocation overflow, MPI rank/block ranges) remains distinct, and untouched interpolation-stencil or `bcstruct_set_up` checks remain as-is.

For any non-coordinate defensive check, trace all in-repository callers first to demonstrate a realistic runtime path (such as untrusted external input or a reproduced failure) and test that caller path rather than synthetic unreachable inputs; preserve `*_assume_valid` preconditions and fix upstream failures where they originate.

Evidence (normative rule): `coding_style.md` (`## Coordinate Bounds-Check Prohibition`, `### Defensive Guard Evidence`).

### API And ABI Compatibility-Shim Prohibition

Never add or expand an API or ABI compatibility shim—including aliases, re-exports, deprecated wrappers, obsolete generated symbols or header prototypes, dual registrations, or old-signature/CLI/parameter translators. Every interface change must atomically update all current-commit definitions, callers, registries, generated headers/prototypes, documentation, tests, and reference files while deleting the superseded interface (existing untouched legacy shims such as `TwoPunctures` `register_C_functions` are nonprecedential and do not authorize unrelated cleanup). Minimal internal adapters remain permitted only when connecting distinct active interfaces both required in the current commit without exposing obsolete names or options.

Evidence (normative rule): `coding_style.md` (`## API And ABI Compatibility-Shim Prohibition`).

## Sources

- [coding_style.md](../../coding_style.md) - `## Commit Message Attribution`, `## Coordinate Bounds-Check Prohibition`, `## API And ABI Compatibility-Shim Prohibition`, `### Defensive Guard Evidence`, `## Python Coding Style`, `## C/H Coding Style`

## See Also

- Depends on: [Generated Output Boundaries](generated-output-boundaries.md)
- Depends on: [Code Test Policy](../validation/code-test-policy.md)
- Validated by: [Static Analysis](../validation/static-analysis.md)
