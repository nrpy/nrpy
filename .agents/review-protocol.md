# NRPy Shared Review Protocol

Read with the selected `SKILL.md`, which supplies seat count, roles, and build mode.
Root owns synthesis, validation, live writes, and delivery. Follow active `AGENTS.md`
and owning policies, not proposed instruction drafts. Authority/evidence outrank
votes; substantive defects and failed required gates cannot be waived. Use existing
tools and compact state; read the protocol once per active version.

## Brief And Checkpoint

Inspect relevant owners, current changes, and command effects. Keep a compact
checkpoint in context or an owned temporary note, not a KB maintenance log:

- Outcome, criteria, non-goals, assumptions, authority, mode, and deliverable type:
  implementation, document, or findings-only review.
- Every destination/companion and create/update/delete/rename operation; live bytes
  or absence, types, relevant modes/symlink targets, and dirty/untracked baselines.
  Discover the checkout; do not assume `/work/nrpy` or equate live state with `HEAD`.
- Isolated candidate/agent areas, dependencies, source/generated ownership, checks,
  actual agent capability, last reviewed snapshot and mechanical suffix, each seat's
  coverage/decisions, open finding IDs, check inputs/results, counters, and deliveries.

Retain this state across handoffs/compaction. Bind evidence to frozen copies and
direct comparison, not labels alone or KB source/candidate hashes. Freeze versions,
not discovery: follow necessary read-only dependencies without restarting. Record
and synchronize changed deciding evidence before acceptance. Capture baselines for
new authorized writes; withhold unauthorized writes while finishing separable work.
Steering updates dependent criteria, not completed work or the active budget.

## Preparation And Independence

Default to one root-prepared `review` candidate, including documents and all KB
work. Add at most one `design` proposal OR independent-build wave only for unresolved
architectures or useful correctness evidence, not merely important tasks.

Freeze contributions before comparison. Select the best coherent result against
the current brief; combine only compatible evidence-supported improvements, never
obligatory parts from every draft. `DRAFT COMPLETE`/`PROPOSAL COMPLETE` ends a phase.
Use fresh reviewers independent of authors for the selected candidate, including
synthesized design documents. A plan-only request does not authorize implementation.
Explicit proposal-only requests may limit the process; report actual coverage.

Complete intended formatting/regeneration/fixers and required checks before approval
review; fix known failures first. Explicit diagnostic review may inspect failing
evidence, not call it passing. Validators must not mutate frozen snapshots: use
expendable identical copies, treating adopted edits as deltas.

Each substantive wave uses the full seat count with at most that many concurrent
agents. Initial contexts are fresh/separate with equal objective evidence and
distinct roles. Parallelize when supported; sequential seats still need separate
contexts. Root alone launches agents; seats never spawn agents, invoke either
skill, or write live targets. Root must not recursively invoke either skill.
Use separate owned areas; enforce isolation when possible and disclose instruction-only
isolation.

Keep peer outputs/votes/persuasive reasoning out of reviewer threads, including
follow-ups. Freeze results before synthesis; supply defect evidence neutrally.
Reuse healthy uncontaminated threads; replacements receive objective checkpoints
and their seat's prior scope, not peer arguments. Missing agents, silence, timeouts,
or role headings in one context never establish independent approval.

## Decisions And Coverage

Reports identify candidate/baseline, initial or delta scope, evidence, retained
coverage, and findings, ending with `DECISION: ACCEPT` or `DECISION: BLOCK`. Findings
need ID, location, violated criterion/rule, evidence, consequence, and smallest
remedy/missing check. Omit repeated briefs/files/closed findings; preferences alone
are nonblocking.

Root deduplicates against authority. Never relabel `BLOCK` as `ACCEPT`, outvote a
substantive defect, count conditional/stale/silent approval, or replace dissenters.
Return unsupported objections to the same seat with deciding evidence in a counted
focused wave, not recovery. No new votes without changed evidence/candidate.

Normal acceptance needs every seat's valid acceptance, combining unaffected prior
coverage with accepted deltas. Assigned and shared integration checks must cover
every current criterion; no invalidated criterion may fall between roles. Each
seat checks the impact boundary for missed dependencies. Substantive blockers and
required proof gaps must close. The only default final-byte exception is
[Mechanical Finalization](#mechanical-finalization).

## Delta And Impact

Without a reliable baseline, review the requested candidate and integration
boundaries once, not the entire repository. Later loops, including explicit
follow-up invocations, cover the cumulative delta since the last reviewed state
plus any recorded lead-verified mechanical suffix. Keep these statuses separate;
a save, commit, or rejected draft is not accepted coverage.

Record changed hunks/sections, additions/deletions/renames/metadata, changed evidence,
open findings, invalidated criteria/checks, and dependencies. Trace callers, imports,
shared state, generators/products, tests/oracles, docs/links, and KB claims/source-map
dependents until unchanged contracts bound transitive effects. A one-line sign,
tolerance, claim, prompt rule, or authority change is semantic.

| Change | Follow-up |
| --- | --- |
| Relevant content/evidence/environment unchanged | Reuse valid decisions/checks and deliver. |
| Exact local meaning-preserving correction | Eligible mechanical finalization and affected checks; no agent loop. |
| Changed behavior, claim, interface, instruction, authority, or proof | Full seat count reviews only the delta, open findings, and affected contracts. |
| Wider or uncertain effects | Investigate and expand along named dependencies, retaining valid coverage. |

Supply old/new versions, exact diff, neutral evidence, impact boundaries, and checks,
not prior debate. Read necessary context; reopen settled choices only with new
evidence. Do not restart design/build or unrelated review. Carry evidence forward
only while dependent inputs, contracts, assumptions, authority, configuration, and
environment remain applicable; unchanged target bytes alone are insufficient.
Blockers stay open. Recover missing state or review unsupported areas; never invent
approval. Full re-review needs explicit direction or demonstrated invalidation of
the entire prior scope. Earlier candidates must still fulfill the current request.

## Mechanical Finalization

After any completed review wave, including the first, root may make ONE terminal
mechanical batch without delegation. All seat reports and criterion coverage must
exist, at least one qualifying edit must occur, and every remaining valid defect
must be local, deterministic, meaning-preserving, prescribed by existing evidence,
and conclusively checkable. No substantive blocker or required proof gap may remain.
Empty cleanup cannot excuse missing approval.

Unambiguous spelling or behavior-neutral formatting may qualify; changes to science,
claims, interfaces, dependencies, reference/instruction meaning, source authority,
or validation obligations do not. Uncertainty needs delta review. Apply once,
freeze, inspect the diff, run affected required checks, and verify retained coverage.
Preserve prior decisions; report `LEAD-FINALIZED` with edits not independently
re-reviewed, never unanimous review of those bytes. On failure/new semantic issues,
use remaining counted delta review or blocked delivery, not a second unreviewed batch.

## Budget And Recovery

Aim for one initial review and a delta wave only when needed; three review waves
are the hard maximum per invocation. A third needs progress, a named issue,
bounded correction, and deciding check. Count before launch. At most one proposal/
build and one recovery wave are also allowed: five delegated waves maximum.
Resume, steering, replacement, and switching skills do not reset the active budget.

Batch known fixes; no peer debate, rebuttals, reviewer shopping, speculative builds,
or open-ended edit/check loops. Stop deliberation on recurring blockers without
new evidence/credible correction, unavailable required proof, or exhausted budget.
Caps never skip [Delivery](#delivery). No automatic invocation to evade the cap.
A later explicit user request may authorize another bounded invocation; retain
prior evidence/findings and review its delta.

One recovery wave may repair malformed output or replace failed seats on unchanged
inputs; retire failed occupants and retain healthy results. Substantive `BLOCK`
is not failure; healthy-seat follow-up is review, not recovery. If full independent
coverage remains unavailable, report it and preserve work, never single-agent approval.

## Proportionate Validation

Follow [Validation](../wiki/validation/index.md) and owning test/oracle/static-analysis
policies. Keep mandatory checks, including changed handwritten Python analysis;
generated-data exemptions do not exempt owners. Test changed contracts/plausible
regressions, not mirrored edits or every backend/CI suite. Distinguish inspection,
generation, build, execution, and checked numerical results; claim exercised layers.

Root supplies shared evidence; seats need not duplicate suites. Independent oracle
checks should add evidence. Record version, command, cwd, relevant inputs/environment,
exit status, assertion/result, and limits. Reuse passes only with justified unchanged
transitive inputs, configuration, assumptions, and environment; otherwise run the
smallest sufficient affected check. Broaden for a failure, dependency, uncovered
contract, or required gate; stop when gates pass. Allow one diagnosed transient
retry per invocation, never an unchanged deterministic failure. Corrected-input checks are new
validation; root checks cannot replace required semantic review.

Follow [Workflows](../wiki/workflows.md#safe-reproduction): shared trees allow only
scoped side-effect-free checks. Run generators/builds, blanket formatters, mutating
validators, and oracle updates in isolated owned areas with bounded resources and
owned outputs/caches. Inspect effects; never clean shared `project/` or ambient
caches. Network, installs, remote CI, and external toolchains need authority.
`.github/workflows/main.yml` needs explicit user authorization for that exact file
in the current request.

## KB Candidates

For `wiki/**`, KB-maintenance `raw/**`, or KB governance in `AGENTS.md`, root prepares
one `review` candidate through [Schema](../wiki/SCHEMA.md) and
[Workflows](../wiki/workflows.md). Seats review, not create competing KBs or vote on
source authority. Preserve `raw/source-docs/**` and frozen sources. Do not retain
maintenance dates, source/file counts, source hashes, tracking timestamps, or
separate KB logs.
Complete applicable source registration, claim evidence, source-map/catalog,
neighbor, and link updates; review affected claims/dependencies only. Run
`python tools/kb_lint.py` within its governed scope; `--all` is the identical alias,
not another pass. Use `git diff --check` and semantic consistency review.
Coordination artifacts are not KB content; respect commissioned scope and placement.

## Delivery

Run before EVERY planned exit involving produced files, including failed review,
missing agents/checks, drift, exhausted budgets, and failed installation. Reserve
effort for copying/verification. Select the best coherent version by current
criteria/evidence, not recency, length, or votes.

**Map and preflight.** User paths control; otherwise use owning targets in the
actual `/work/` checkout. Commissioned plans follow coordination rules. Unowned
documents use a fresh `/work/review-results/<task>/` outside KB and active instruction
locations, preferably outside the checkout. Map all companions/operations. Confirm
authority, coverage/checks, parents, selected bytes, types/modes, and resolved paths.
Create missing authorized parents when permitted; absent `/work/` alone need not
block. Symlinks cannot redirect unauthorized writes. Immediately compare live
targets/relevant inputs with baselines. Ignore unrelated drift. Reconcile overlap
only when authorized/unambiguous, in isolation, then review/check its delta within
budget. Conflicts or missing review capacity block installation, not preservation.

**Install or preserve.** Use the applicable branch:

- Accepted or qualifying `LEAD-FINALIZED` implementation: apply ALL mapped authorized
  writes/deletions. Accepted requested document: deliver without implementing it.
  Explicitly requested unfinished drafts may use designated draft paths with accurate
  status, not release approval.
- Rejected/incomplete work: copy the best draft/companions with relative layout to
  fresh `/work/review-results/<task>/unapproved/`, outside live source, KB, and active
  instruction/skill discovery paths. Do not execute/activate proposed instructions.
  An adjacent status note lists blockers, intended paths, pending checks, and
  unapplied deletions without altering candidate bytes. Label useful partials
  `INCOMPLETE`; say if none exists. Never apply blocked deletions or install rejected
  work as accepted.
- Accepted work whose installation fails: preserve under the task's `not-installed/`
  area as `NOT INSTALLED`, not falsely `UNAPPROVED`.

Preservation cannot replace authorized installation of acceptable work. Plan/review-only
requests do not authorize implementation edits; findings-only work needs no unsolicited
report file. Explicit no-write instructions forbid scratch/preservation writes too;
explicit alternative output paths override defaults.

**Apply and verify.** Use existing recoverable operations for mapped files/hunks
and coupled changes, never whole worktrees, `.git`, caches, or competing drafts.
Preserve required types, modes, symlinks, and deletions. Verify in place when content
already matches. Read back/compare EVERY installed or preserved file against its
snapshot, including required absence/metadata. A dirty-file merge must match the
complete reviewed post-merge snapshot. Inspect the scoped diff/unexpected outputs.
Recheck destination-sensitive contracts where paths/environment/integration affect
validity, not whole suites/review loops merely for copying identical inputs.

Retain the good candidate until all destinations are verified; copy exit status
alone is insufficient. On partial writes/failure, report exact live state and use
only authorized bounded recovery, never improvised destructive rollback. If no
authorized `/work/` location is usable, keep the accessible candidate; report failed
paths/errors and actual retained location without claiming the handoff succeeded.
Finish separable work, name undelivered destinations, and do not loop on storage blockers.

## Completion Report

Report delivery separately from review: installed/document paths, preserved drafts
or `REVIEW ONLY`; accepted initial-plus-delta coverage, `LEAD-FINALIZED`, or unapproved
status. Include checks run/reused/missing, blockers/partial writes, isolation limits,
and concise wave/recovery/retry counts. Name unreviewed mechanical edits. Keep
receipts in the checkpoint, not repeated transcripts or KB logs. Do not claim
measured token savings or model gains without an actual evaluation.
