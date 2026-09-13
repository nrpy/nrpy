# Response to the PR #211 code review

Thank you for the detailed review. We re-evaluated the published report and its
supporting logs against the current tree using two independent dialectic seats.
All seven published issue classes were valid and have been addressed. The final
post-edit review concluded `ACCEPT` from both seats.

## Resolved findings

1. **Coordinate-admission control flow** — Removed floating coordinate
   finiteness, range, tolerance, and domain-membership predicates from mesh
   admission, block geometry, exterior writes, and the two host harnesses.
   Exterior ownership now uses Dendrolib's discrete block-boundary flags and
   padded storage indices. Coordinate arithmetic remains only where it computes
   numerical reference values, not where it admits, rejects, skips, or writes
   points. The affected validation page now describes this boundary mechanism
   accurately while retaining general application boundary conditions as an
   open gate.

2. **Unsafe `char[N]` defaults** — Dendro's CodeParameter emitter now requires
   a string without embedded nulls and emits its UTF-8 bytes using safe literal
   characters or fixed-width three-digit octal escapes. A compiled C++17 probe
   verified byte-for-byte preservation for quotes, backslashes, control
   characters, and non-ASCII input.

3. **Missing semantic `END` markers** — Added markers for the enclosing
   standalone `-r` branch, the real-host `try`, and the capability harness's
   padded-point branch. The previously present inner-block markers did not mark
   the two enclosing constructs, so those two template changes were necessary.

4. **Multiprecision `Any` annotations** — Replaced the explicit `Any` input,
   result, and accumulator annotations with the concrete accepted union and
   `mpmath.mpf`. Exact float/rational conversion and both precision passes are
   unchanged. Current mypy accepts the result.

5. **Incorrect `KO_ENABLED` documentation** — The docstring now states that
   real-host profile validation compares `ko_enabled` with `KO_ENABLED`, while
   emitted numerical code does not use it to derive stencil reach.

6. **BSSN parameter classification** — Both documentation occurrences now call
   `register_MU_gridfunctions` and
   `register_M_and_LAMBDA_CONSTRAINT_gridfunctions` generation-time NRPy
   parameters. Their names, defaults, memoization role, and scoped
   save/set/restore behavior remain unchanged.

7. **Transient qualification snapshots** — Removed dates, platform/toolchain
   tuples, run results, counts, measurements, and instructions to rewrite those
   receipts from `tests_infra/README.md`. Stable repository pins, setup and run
   commands, tested axes, assertions, tolerances, and expected markers remain.

## Verification

- Affected Python owners pass the repository single-file static-analysis gate
  under Python 3.12.3 and mypy 2.3.1.
- Repository-wide `black --check .` passes: 875 files unchanged.
- `python tools/kb_lint.py` and `git diff --check` pass.
- Generated default BSSN and fCCZ4 projects build and pass 28/28 CTests,
  including both non-flat GR reference tests.
- The pinned Dendrolib capability harness passes on one rank, two ranks, and
  order 10 with padding 5.
- The pinned real Dendro host builds and passes two-rank transport, parameter,
  and 100-step Minkowski checks.
- Deliberate offset, halo, and nonfinite-state corruptions each produce a
  nonzero failure.
- No new test files were added. Existing harnesses were adjusted where required.
- `raw/source-docs/` and `.github/workflows/main.yml` are unchanged.

## Rejected ancillary candidates

The supporting logs contained additional candidates that do not describe
current-tree defects:

- `_codeparameter_tail` and `_term` have multiple production calls;
  `gf_array_name`, ADM conversion, role metadata, and the cited public helpers
  have active consumers or supported direct use.
- Generic and GR Dendro emitters have separate backend and formulation
  responsibilities. Dendro consumes the canonical CodeParameter and finite-
  difference registries rather than replacing them.
- No equation, constraint, connection, projection, initial-data, cache-key,
  BHaH-factory, component-offset, x-fastest-layout, halo, padding-reach, or
  derivative-token defect was reproduced. Generated symbolic and runtime
  oracles passed.
- Default and FD6+KO generation, real-host context initialization, runtime
  parameter parsing, and the pinned-host integration are not broadly broken.
- FD8, output production, physical boundary conditions, remeshing, local time
  stepping, restart, GPU, threading, and optional upstream configurations are
  not promises of the reviewed interface.
- The profile label, projection diagnostic, standalone maximum extent, upstream
  revision handling, and direct-header preferences have no demonstrated harmful
  behavior under the governing contract.
- Blanket marker claims are overbroad: trivial braces, helper-emitted loop
  footers, the short refinement lambda, and the compact one-line `StencilTerm`
  aggregate are not established violations. Only the three published closers
  required changes; no separate marker-test requirement exists.
- Project policy does not require eliminating every third-party-boundary `Any`,
  accepting default `clang-format` output as an oracle, adding a local pytest,
  or treating code not reached by CI as unused.
- The unchanged BHaH character-default emitter and permissive malformed-array
  parsing are pre-existing and outside this Dendro-scoped remediation.
- Historical claims concerning `kernel_lowering`, superseded helpers and future
  slots, uninitialized pointers, empty runners, point-loop defaults, authorship,
  fail-open setup, obsolete type names, `FCCZ4RHSBuild`, tracked work records,
  workflow permissions, and earlier marker locations are absent or already
  corrected in the current tree.

Coverage ideas such as same-process registry reuse, Python 3.7 execution,
real-host FD2/FD6, optional host modes, the complete capability fault matrix,
multiprocessing, and nonlinear evolution remain unproven extensions rather than
review findings. Under the minimal-sufficient-work contract, they were not
pursued as fixes.
