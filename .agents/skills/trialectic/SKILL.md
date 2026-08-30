---
name: trialectic
description: Use when the user explicitly asks to engage tri, use tri, run tri, or invoke the NRPy trialectic, requests three-agent independent NRPy review, or repository policy requires it. Produce or review an NRPy change through three fresh-context agent seats, repository evidence, proportionate validation, at most three candidate cycles, and either unanimous acceptance or one tightly bounded lead-finalization pass for minor residual defects. Do not invoke for routine NRPy edits without an explicit trigger or policy requirement. NRPy KB work is review-only.
---

# NRPy Trialectic

Use three independent agent seats to create or review one NRPy candidate. Preserve independent evidence long enough to expose defects, then resolve findings against the frozen brief, authoritative repository evidence, and executable checks.

Root is the lead agent. Root owns discovery, candidate preparation, synthesis, mutation, validation, installation, and reporting. Delegated agents never write to live targets.

Unanimous evidence-bound acceptance is the normal release gate. The only exception is the terminal lead-finalization rule after the third candidate cycle. That exception may close only minor, deterministic residual defects and must never be reported as unanimous acceptance.

This workflow has no required runtime companion files. Use existing repository and environment tools, and keep run evidence in temporary work areas or the completion report.

## Invariants

- Preserve the user's request, authorization boundaries, unrelated work, and concurrent changes.
- Treat the user request, applicable repository instructions, authoritative project sources, executable contracts, and validation results as stronger evidence than agent agreement.
- Use exactly three independent seats in each delegated phase and no more than three delegated agents concurrently.
- Neither root nor a delegated agent may recursively invoke Dialectic or Trialectic. Delegated agents may not create subagents. Root launches every delegated wave directly, and every such wave counts against this invocation's budget.
- Start each independent triad in new agent threads with the same frozen brief and evidence but distinct roles. Do not give any occupant prior agent reasoning, peer output, or peer conclusions before all initial results are frozen.
- If the environment cannot provide three complete independent results, stop rather than silently degrade to fewer seats.
- Permit `NO ACCEPTABLE CANDIDATE`. Never force synthesis, weaken requirements, manufacture consensus, or convert silence into approval.
- Never let agreement or lead finalization override a failed required check, unsupported scientific claim, repository contradiction, protected-file rule, source-authority rule, or unresolved evidence gap.
- Treat applicable `AGENTS.md` files as harness-supplied ambient instructions. Do not target, modify, copy into evidence packets, inventory for this workflow, quote, summarize, or hash them.

## Select the least expensive useful mode

Record the mode and rationale before delegation.

- `review` - Default. Root prepares one complete candidate; three independent agents review and decide on it. Use for bounded implementations, bug fixes, documentation, configuration, deletion, established patterns, and all NRPy KB work.
- `design` - Three independent agents propose designs; root prepares one implementation candidate; a new independent triad reviews and decides on it. Use only when materially different valid architectures remain after repository discovery.
- `tri-build` - Three independent agents build isolated candidates; root selects or synthesizes one; a new independent triad reviews and decides on it. Use only when independent implementations materially improve confidence, provide an oracle, or reduce high technical risk.

Do not choose `design` or `tri-build` merely because a task is important. Prefer `review` when repository patterns and objective validators substantially determine the implementation.

## Hard effort budget and stop conditions

One invocation has one frozen scope round. Root may refine a draft brief during root-only discovery, but must freeze it before candidate-mutating preparation or the first delegated phase. After that freeze, a material scope change, required scope expansion, contract defect, or material pre-install drift ends the invocation. Report `RESTART REQUIRED` or the more specific blocker; do not automatically begin another scope round. A later explicit invocation starts with fresh counters and agents.

Use these global limits across the entire invocation:

- at most one design-proposal or `tri-build` phase;
- at most three candidate-review cycles total;
- at most one combined review-and-decision wave per candidate cycle;
- at most one root-owned correction batch between consecutive candidate cycles;
- no peer cross-critique, free-form debate, or rebuttal waves;
- at most one agent-recovery wave for the entire invocation; and
- at most one validation retry for the entire invocation, limited to a diagnosed transient infrastructure or service failure.

The pathological maximum is five delegated waves: one proposal or build wave, three review-and-decision waves, and one recovery wave. Ordinary `review` mode normally uses one delegated wave. Skip every wave that is not necessary. Root-only correction batches do not add delegated waves, but each is limited to one batch between cycles and must be based on already-known findings.

Root records counters in the run evidence and increments a counter before launching a wave or binding a cycle. If the increment would exceed a limit, stop without launching that work. Counters never reset within an invocation.

Two candidate cycles are the normal ceiling. Open a third and final cycle only when cycle two made concrete progress, the scope and validation plan remain fixed, and one bounded correction or evidence acquisition has a credible path to success. Record the remaining issue, planned change, and deciding check before opening cycle three. Never open a cycle merely to seek different votes.

Stop early when any of these occurs:

- the same blocker or evidence gap recurs in two consecutive cycles without materially new evidence or a materially different correction;
- a proposed correction would violate scope, source authority, preservation rules, or another acceptance criterion;
- required evidence, fresh-context delegation, candidate identity, or a required validator is unavailable;
- a deterministic required check fails and no bounded in-scope correction is identified;
- material drift is detected after scope freeze; or
- the next required wave would exceed a global limit.

Run each planned validation command once per candidate identity. The single validation retry may rerun only a command with a diagnosed transient failure. Never rerun an unchanged deterministic failure merely to seek a pass.

On a corrected candidate, rerun every check whose read set may intersect changed bytes. Carry forward a prior passing result only when root can mechanically establish that its relevant inputs and environment are unchanged; record that justification.

Keep reports compact. Agents cite exact locations and material evidence instead of restating the brief, reproducing candidate files, or repeating equivalent findings. Root deduplicates materially identical findings.

## Freeze the scope round

Before delegation, root performs enough read-only discovery to prepare one frozen brief. Do not modify live targets during discovery.

The brief states:

- exact user request, acceptance criteria, and explicit non-goals;
- exact target paths and create, revise, or delete intent;
- material target baselines, including object type, contents or required absence, permissions, and symlink target when relevant;
- bounded implementation inputs needed to understand the change;
- applicable repository rules and authoritative sources;
- material assumptions, compatibility requirements, and preservation rules;
- generated-output ownership and source-of-truth rules;
- selected mode and rationale;
- proportionate validation obligations and known unavailable checks; and
- a separate temporary area for every agent plus one root-owned candidate area outside the live target tree.

Bind writes to exact target paths. Bound reads to decision-relevant paths or repository areas; do not inspect unrelated regions. Validation may read transitive dependencies required by imports, compilers, generators, test collection, builds, or link checking.

An agent that needs material evidence outside the brief returns `SCOPE EXPANSION REQUEST` with the exact path or bounded area, why it is needed, and how it may affect the result. Root denies an unnecessary request. If the evidence is necessary, stop with `SCOPE INCOMPLETE` and report the exact expansion needed; do not authorize it mid-run.

After scope freeze, a material change to the request, acceptance criteria, target set or intent, target baseline, authoritative inputs, assumptions, repository rules, selected mode, or validation obligations ends the invocation with `RESTART REQUIRED`.

## Assign complementary roles

For software, symbolic, numerical, generated-code, build, configuration, interface, or developer-tooling work, use these roles.

### Seat 1: NRPy scientific-software verifier

Check mathematical and executable correctness, assumptions, invariants, interfaces, failure cases, generated-output ownership, compatibility, numerical behavior, backend implications, and validation oracles. Select applicable checks such as symbolic equivalence, deterministic regeneration, generated-code compilation and execution, reference-oracle comparison, justified numerical tolerances, convergence or order verification, and affected-backend coverage.

### Seat 2: NRPy integration and simplification reviewer

Check architecture, module conventions, public interfaces, symbolic and code-generation idioms, infrastructure boundaries, generated-output ownership, compatibility, maintainability, and the simplest-sufficient implementation rule. Inspect authorized implementation inputs for established NRPy functions, classes, APIs, patterns, and components. Name what was considered, reuse sufficient machinery, and reject speculative generality, duplication, needless indirection, unjustified layers, wrappers, helpers, dependencies, configuration, and abstractions.

### Seat 3: NRPy standards and release auditor

Apply NRPy coding and scientific-software standards with the deepest cross-cutting scrutiny. Trace acceptance criteria to candidate evidence and implementation choices to repository authority. Audit file placement, dependency direction, public and private API boundaries, symbolic and code-generation idioms, source-versus-generated ownership, deterministic regeneration, validation-oracle adequacy, static-analysis expectations, compatibility and transition handling, documentation and test integration, build integration, backend coverage, packaging, and stale references. Expose omissions that local correctness or simplification review can miss.

For non-software work, use one first-principles analyst, one domain-and-implementation expert, and one adversarial editor and verifier. All three still apply applicable NRPy conventions.

Role labels do not establish correctness. Every material conclusion requires repository evidence, derivation, executable validation, an authoritative source, or an explicit limitation.

## Produce independent work

Give all occupants the same frozen brief, authorized evidence, acceptance criteria, and validation plan. Vary only role instructions. Do not expose peer output or conclusions.

### Review mode

Root prepares one complete candidate in the root-owned area. Complete candidate-mutating preparation, including applicable formatting, regeneration, and fixers, before binding a candidate cycle.

### Design mode

Each proposal agent independently returns one complete design covering all target operations, architecture and interfaces, assumptions, established components to reuse, scientific and executable risks, compatibility, validation strategy, and unresolved questions. End with exact line `PROPOSAL COMPLETE`.

Freeze all proposals before comparison. Root selects one coherent proposal, combines only compatible evidence-supported elements, or returns `NO ACCEPTABLE CANDIDATE`. Root prepares one complete candidate. Retire proposal occupants and use a new independent triad for candidate review.

### Tri-build mode

Each builder independently creates one complete candidate in an isolated area and reports rationale, assumptions, authorized inputs inspected, established components considered and reused, validation performed, limitations, and final exact line `DRAFT COMPLETE`.

Freeze all candidates before comparison. Root verifies scope and completeness, then selects one coherent candidate or synthesizes only compatible evidence-supported elements. Do not average incompatible implementations or include an element merely to represent every draft. If none is sufficient, return `NO ACCEPTABLE CANDIDATE`. Retire builder occupants and use a new independent triad for candidate review.

There is no peer cross-critique phase. Root resolves differences directly against the brief, repository evidence, and validators.

## Bind and validate a candidate cycle

A candidate-review cycle binds one canonical candidate within the frozen scope round.

For each cycle, root:

1. finishes all intended candidate-mutating preparation;
2. assigns the cycle identifier;
3. freezes one canonical candidate snapshot;
4. records every target operation and required presence or absence;
5. establishes candidate identity mechanically;
6. runs all required validation that should inform review; and
7. launches one combined review-and-decision wave only if required validation passes.

For ordinary files, use an immutable snapshot, direct byte comparison, content digest, or an equivalent existing mechanism. For NRPy KB sources and candidates, never hash; use frozen copies, non-hash snapshots, or direct byte equality permitted by the normal KB workflow.

If exact candidate identity cannot bind review, validation, decisions, and installation to one state, stop.

Validators must not mutate the canonical snapshot. Run an inherently mutating validator on an expendable byte-identical copy. Adopting any intended validator delta creates a new candidate and consumes the next cycle.

As applicable, validate target confinement, unrelated changes, stale references, generated-output ownership, syntax, formatting, lint, static analysis, imports, documentation structure, links, focused behavior, failure cases, symbolic identities, assumptions, indices, reference quantities, deterministic regeneration, generated-target builds and execution, numerical order or convergence, justified reference oracles and tolerances, and affected backends.

Record each command's working directory, invocation, exit status, relevant output or explicit truncation, candidate identity before and after, unavailable dependencies, untested backends, manual limits, and tolerance rationale. Never record secrets.

Any candidate mutation ends the current cycle. Any substantive evidence change after review invalidates that cycle's decisions. Corrections normally consume the next cycle and require a new independent triad. Count cycles globally; there is never a fourth review cycle. Between consecutive cycles, root performs at most one bounded correction batch based on the frozen findings. If that batch cannot produce a validatable candidate, stop rather than iterating through uncounted root-only edit/check loops.

## Review and decide in one wave

After required validation passes, each reviewer independently inspects the same frozen candidate and evidence and returns one compact report containing:

- rationale and material assumptions;
- exact inventory reviewed;
- role-specific and applicable NRPy evidence;
- every material finding, or an explicit statement that none was found; and
- final exact line `DECISION: ACCEPT` or `DECISION: BLOCK`.

Every material finding states:

1. exact location;
2. violated acceptance criterion, repository rule, invariant, or authoritative standard;
3. concrete evidence;
4. consequence; and
5. smallest sufficient correction or missing evidence.

The only valid decisions are:

- `ACCEPT` - No evidence-backed material defect or required evidence gap remains for the frozen brief, candidate, and supplied validation evidence.
- `BLOCK` - One or more concrete material defects or required evidence gaps remain, each using the five-part format.

Preferences, unsupported suspicions, enlarged scope, conditional approval, silence, timeout, stale review, or approval of another candidate do not pass. Do not ask reviewers to agree with one another.

Root classifies each finding as exactly one of:

- `CONTRACT DEFECT` - The frozen brief is materially ambiguous, contradictory, incomplete, or wrong. Stop with `RESTART REQUIRED`.
- `BLOCKER` - Evidence demonstrates a candidate defect. Correct it before normal acceptance.
- `EVIDENCE GAP` - A required property remains unresolved. Obtain already-authorized evidence or stop.
- `MINOR FINALIZATION DEFECT` - Available only after cycle three. The issue meets every terminal criterion below and is eligible for the one lead-owned correction batch.
- `NONBLOCKING TRADEOFF` - The observation is real but compatible with the brief. Record it.
- `NOISE` - The finding is an unsupported preference, invented requirement, or claim contradicted by authoritative evidence.

Do not reject a valid finding because another reviewer disagrees. Do not accept a finding merely because every reviewer repeats it. Root may not overrule an evidence-backed substantive `BLOCK`. `MINOR FINALIZATION DEFECT` is unavailable in cycles one and two and may not be used when materiality is uncertain.

If all three reviewers return valid `ACCEPT` decisions for the same cycle, the candidate may proceed to installation. Otherwise root resolves blockers and evidence gaps with the smallest evidence-supported correction and opens the next cycle if one remains. Never replace, retry, or resample a substantive `BLOCK` in search of approval.

## Terminal lead finalization after cycle three

There is no fourth candidate-review cycle and no delegated call after this terminal step.

If cycle three does not obtain unanimous acceptance, root may perform exactly one lead-owned correction batch only when every remaining valid issue is a `MINOR FINALIZATION DEFECT`. An issue qualifies only when all of these are true:

- its exact location and correction are already determined by the frozen brief, authoritative repository evidence, or a required validator;
- the correction is local, bounded, and deterministic rather than a design choice;
- it does not change scope, target set, dependencies, architecture, public interfaces, mathematical or scientific behavior, algorithms, source authority, generated-output ownership, compatibility promises, backend behavior, or validation obligations;
- it requires no new implementation evidence and leaves no unresolved evidence gap; and
- existing required checks can conclusively verify it.

Root must not relabel a substantive blocker as minor to obtain completion. Typical qualifying defects are mechanical consistency errors, behavior-neutral formatting or spelling defects, mechanically identified stale-reference cleanup that cannot alter behavior or meaning, and other repository-prescribed edits with no semantic discretion.

For a qualifying terminal correction, root:

1. records each residual finding and why it is minor;
2. applies all final corrections once, in one batch, to a mutable copy of the cycle-three candidate;
3. freezes and identifies one terminal candidate;
4. reruns every affected required check exactly once, using the unused transient retry only when applicable;
5. confirms no protected rule, source-authority rule, scope binding, or live baseline changed; and
6. either installs that exact terminal candidate or stops.

If every required check passes and installation preflight succeeds, root may install and report `LEAD-FINALIZED AFTER CYCLE CAP`. The report must state that the terminal bytes were lead-finalized and were not unanimously re-reviewed. This is the sole exception to unanimous candidate approval.

If any residual issue is not minor, any required evidence remains missing, any required check fails, or preflight detects drift, do not install. Report `NO ACCEPTABLE CANDIDATE` and stop. Never begin another cycle, seek more votes, or perform a second lead-correction batch.

## NRPy KB review-only mode

Any candidate involving `wiki/**` or KB-maintenance files under `raw/**` uses `review` mode for the complete atomic KB-workflow change.

Before candidate preparation, root identifies the repository-defined KB schema, source-authority rules, and normal KB workflow. If they cannot be established from authorized repository evidence, stop rather than infer them.

Root creates the sole candidate through that workflow and completes every applicable source registration, dependency update, source-map update, catalog update, claim-evidence operation, link check, KB lint check, and other required step before review.

Additional invariants:

- Repository evidence and the KB schema determine source authority; agent agreement never does.
- All three reviewers inspect the same root-prepared candidate and evidence. They do not create alternatives, edit candidate files, choose source authority, or install changes.
- Preserved sources, including `raw/source-docs/**`, are immutable.
- Exclude `AGENTS.md` from candidate data and review packets.
- Never hash KB sources or candidates.
- Root applies corrections only through the normal KB workflow.
- Any KB candidate mutation before cycle three normally consumes the next cycle and uses a new independent triad.
- Installation occurs only through the normal KB workflow.

For KB work, terminal lead finalization is allowed only through the normal KB workflow and only for mechanical formatting, an exact repository-authoritative reference correction, or mechanically prescribed metadata. It may not alter claims, evidence relationships, source authority, schema, preserved sources, dependency meaning, catalog meaning, or workflow obligations. Rerun every affected KB workflow check before installation.

## Install and verify

Immediately before the first live write, root confirms:

1. the mutation remains authorized;
2. live targets and decision-relevant inputs have not materially drifted from the frozen scope;
3. the release-candidate identity and its validation evidence are unchanged since that release candidate was frozen; and
4. either all three valid `ACCEPT` decisions bind to that candidate, or the terminal candidate satisfies every lead-finalization requirement above.

Material drift ends the invocation with `DRIFT DETECTED`; do not rebase or restart automatically. Never overwrite concurrent work because an earlier candidate was accepted.

Apply only the release candidate's authorized writes and deletions. Prefer repository-native atomic or recoverable operations. For KB work, use only the normal KB workflow.

Mechanically verify installed state equals the release candidate, including exact file bytes, required absence, object type, permissions, and symlink target when relevant. Rerun applicable installed-tree checks and inspect for unauthorized outputs or validator mutations.

If installation is partial, validation fails, installed state differs, or unexpected drift appears, stop and report the exact live state. Use only an already-authorized recovery mechanism; do not improvise destructive rollback.

## Agent failure

The entire invocation permits one agent-recovery wave. It may ask mechanically failed or malformed seats from one phase to correct their output or replace those seats with fresh-context occupants using unchanged inputs. Retire failed occupants before replacement. If recovery is unavailable or fails, stop. Recovery is unavailable after terminal lead finalization begins.

A substantive `BLOCK` is not an agent failure. Never resample it.

## Report

Report:

- selected mode, scope-round identifier, and final cycle identifier;
- candidate-cycle, root-correction-batch, and delegated-wave counts;
- whether agent recovery or the transient validation retry was used;
- created, changed, and deleted targets;
- validation and installed-tree verification outcomes;
- unavailable checks, untested backends, redactions, and other limitations;
- established NRPy components reused and justification for new machinery;
- material nonblocking tradeoffs;
- the final review decisions and their candidate binding; for lead finalization, preserve all three cycle-three decisions that preceded the terminal correction;
- whether the outcome was `UNANIMOUSLY ACCEPTED`, `LEAD-FINALIZED AFTER CYCLE CAP`, or `NO ACCEPTABLE CANDIDATE`;
- for lead finalization, the exact minor corrections and the fact that terminal bytes were not unanimously re-reviewed;
- whether isolation was mechanical or instruction-enforced;
- temporary-area cleanup status; and
- the exact blocker when installation did not occur.
