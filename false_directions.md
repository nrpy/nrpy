# False positives and false directions

## Current runtime implementation — 09-07-2026

CONTR-0011 is resolved by one builder and one registrar per RHS/constraint
module. Full registered CFunction text, evolution/diagnostic ordering, and
Dendro registration metadata match the prior implementation for both shipped
profiles; existing trusted-expression checks pass.

The generated fCCZ4 solver now builds against Dendro-GR
`b3261e2a0d3457781b11d63ac5ab38375ffab93b` and Dendrolib
`246043709e806021fcfc011fe657b8bf964cae4c`. Real block offsets and padded
origins, distributed unzip–RHS–zip, ETS callbacks, TOML parameters, and
collective failure handling are implemented. Two active ranks pass the
independent affine transport and analytic eta-response checks, rank-local
negative cases, and a checked 100-step fixed-mesh Minkowski run. Reproduction
and numerical limits are in `nrpy/infrastructures/Dendro/tests_infra/README.md`.

This closes the requested bounded PR 7 host portion. General boundary
conditions, remeshing/LTS, checkpoint/restart, output, GPU/thread execution,
BSSN real-host qualification, and real-host CI remain unqualified. The existing
workflow is unchanged. Earlier checkpoints below retain their historical scope;
they do not describe the current runtime implementation or review status.


## Previous cleanup checkpoint — 09-07-2026

The earlier sections are historical records, not the current acceptance or
working-tree status. The reviewed implementation is now committed as
`4d97621e` (39 files, +696/-454); earlier staged/uncommitted descriptions and
wave counts refer to their original snapshots. Historical votes and build
claims have not been independently replayed by this cleanup.

The current register leaves CONTR-0011 open for the private builder and
registrar per formulation in Dendro RHS/constraints. Only the BSSN leaf is
`contested`; the fCCZ4 leaf is `provisional`. The pure build/register pairing
objection was withdrawn because established BHaH counterparts exist. This
supersedes the earlier two-provisional-leaves/two-divergences account in
progress §29.5 and false_directions §5b.2–§5b.3. The Schema's page-status matrix
still governs status; historical accounts cannot override it.

The cleanup corrects the stale quotation, stencil comparison, source coverage,
namespace count, and review rationales; adds a full-downwind reach doctest;
and makes small prose/comment clarifications. See `issues_todo.md` for each
disposition. CONTR-0011 consolidation and real-host integration remain open;
this cleanup does not implement either.

Date: 09-07-2026. Branch: `dendro-infra`, candidate staged on baseline
`93290e10`.

## Purpose and status

This is an untracked working document, in the same class as
`inconsistencies.md` and `progress.md`. It is never committed. It records every
claim, measurement, remedy or process step in the Dendro review effort that
turned out to be wrong, together with what disproved it, so the same wrong turn
is not taken again and so a reader of the review record knows which parts of it
not to trust.

It covers the effort that applied eleven entries of `inconsistencies.md` (I1,
I12, K1, K2, S15, S11, S2, I14, S7, S1, D2), across four review waves and
twelve independent seat reports, plus the earlier effort that landed the
infrastructure and the D5 trusted-baseline work.

Three kinds of entry appear here, and the distinction matters:

- **False positive** - something was reported as a defect and was not one.
- **False direction** - a claim, measurement or remedy that was itself wrong,
  usually mine, and that other work was then built on.
- **Declined, not false** - a real observation deliberately not acted on. These
  are recorded at the end so they are not re-litigated as though they were
  oversights.

---

## 1. Errors in `inconsistencies.md` itself

That document was produced by an independent panel and is careful, but it is
not authoritative and it does not authorize its own remedies. Two of its
entries are wrong, and one of its suggested remedies would have produced a
tautology. Everything else in the eleven verified cleanly.

### 1.1 K2's premise is false - two leaves never contradicted each other

`inconsistencies.md` K2 reads:

> `fccz4-application-wiring.md:18-19` states the opposite correctly, so two
> Dendro leaves contradict each other - a `wiki/contradictions.md` matter under
> `wiki/SCHEMA.md`.

They do not contradict each other. `fccz4-application-wiring.md:18` contains
the page's only sentence with the word "split", and it is about a *different*
split - the pairing of a pure `build_*` function with a `register_CFunctions_*`
function:

> Each family pairs a pure `build_*` function with a `register_CFunctions_*`
> function so one profile can assemble a different subset without the builders
> knowing about each other. That split is Dendro's own: the established
> infrastructures build and register inside a single `register_CFunction_*`.

The per-formulation private-builder split K2 is about is recorded on neither
leaf. Nor is the duplication it leaves. Found by the first wave's verifier seat
and confirmed independently by two later seats.

The half of K2 that *is* true, and that the remedy correctly fixed, is the
narrow one: the BSSN leaf called the private-builder-per-formulation layout
"BHaH's own arrangement", and `BHaH/general_relativity/rhs_eval.py` holds one
public function and no private builders.

**Consequence.** I took the premise on trust and wrote a remedy asserting that
the fCCZ4 leaf "describes the split the same way" and that the duplication "is
recorded there". Both false. One false claim was replaced by two, and
`CONTR-0010` inherited the error in five of its thirteen columns before it was
corrected. `CONTR-0010`'s Notes now record that the survey was wrong here.

**Lesson.** A findings document's *premise* needs verifying as much as its
conclusion. Where an entry says "page X says Y", open page X.

### 1.2 I14's alternative remedy would have been a tautology

I14 offers, as its fallback shape:

> If the emitted C++ case is to stay independent of `REQUIRED_PADDING`, emit a
> second generated constant carrying the exact reach and compare two
> registry-supplied numbers - which is what the comment claims it does.

Both numbers would come from the same source - `stencil_reach_per_axis` through
`padding_from_derivative_operators` through `roles.set_required_padding` - so
the comparison could not fail through any path the generator can produce. It
would have replaced a wrong check with an unfalsifiable one. The route taken
instead keeps a single authority for the reach and states the emitted case's
limit honestly. Three seats across two waves examined the substitution and
accepted it; the third wave's verifier confirmed the surviving bound is
one-sided but sound rather than vacuous.

### 1.3 Measurements the document itself withheld or corrected

To its credit, the panel recorded these rather than asserting them, and each
is a false positive avoided:

- The peer `__main__` `else:` ratio in S15: two review seats measured 98/111
  and 97/113, so no peer count is asserted and the finding rests on the rule
  text plus the Dendro-side 23-of-23 count.
- The BHaH `register_CFunction` field ratio in S1: disputed, so omitted.
- I10's remedy space: the document closes one of its own options mid-flight -
  core cannot absorb the generation-parameter check, because `NRPyParameter`
  has no allowed-values concept and the Dendro table is *narrower* than the
  parameter's real domain, so adding one would make `phi` illegal for BHaH.
- The `standalone_host/` question, which stood as a dispute and was resolved
  against the `BHaH/interpolation/` precedent: the directory is conformant and
  only the `__init__.py` import remained.
- One survey's clean list credited Dendro's `generate_default_parfile` as
  mirroring BHaH's function of the same name, which conflicts with I11's
  module-name finding. The document discloses the partial clearing rather than
  hiding it.

---

## 2. Mechanical tests that report false positives

`wiki/infrastructures/new-infrastructure-conformance.md` carries mechanical
greps as conformance aids. Two of them fire on the Dendro tree and both are
false positives on inspection. They are recorded in `inconsistencies.md` §4 and
re-confirmed here because a future audit will hit them again.

- **The formulation-name test.** "Expect zero" case-insensitive `fccz4` hits in
  the generic `Dendro/*.py` layer. The tree has six, and all six are doctest
  fixtures. The generic layer carries no formulation name in any load-bearing
  position.
- **The `snapshot|Frozen[A-Z]` test.** Reports 13 hits over `*.py`, 14 counting
  a prose line in `dendro_standalone_host.h`. Every one is emitted C++ for the
  lifecycle drift gate (`Ctx::snapshot_state`, `max_drift_from_snapshot`). No
  Python duplicate registry survives.

A third mechanical check is worth naming as a *structural* false negative
rather than a false positive: `python tools/kb_lint.py` passed at every single
stage of this effort, including while `wiki/contradictions.md` contained a dead
symbol citation, a retracted premise and an absolute its own resolution test
refuted. The linter checks structure and link resolution. It does not check
whether a cited symbol exists or whether a claim is true. Three separate
blocking findings landed in that one file with the linter green throughout.

---

## 3. My own false directions

These are the substantive ones. They are grouped by kind because the kinds
recur.

### 3.1 Process failures

**A patch script aborted mid-way and silently skipped every later edit.**
The wave-two batch applied nine substitutions in sequence. The seventh raised
an `AssertionError` because I had guessed the target text wrong, and Python
exited - so the eighth and ninth, including the `c-and-embedded-c-style.md`
"both forms" correction, never ran. I reported the batch as applied.

This is the second occurrence of exactly this failure mode in this effort. The
first, during the D5 work, left `infrastructure-code-style.md` asserting that
the RHS sweep "captures nothing" after a later assertion aborted the batch, and
the lesson I drew then was to apply substitutions one file at a time. I did not
carry it forward. The final batch was rewritten so each edit is independent,
reports `applied` or `SKIPPED` per edit, and exits non-zero if any edit did not
land.

**Editing `/work` while seats were still reading.** During the earlier D5
effort my brief told three seats the candidate was frozen, and I then edited the
tree while two of them were mid-review; one intermediate state raised
`NameError`. Flagged by a seat as review drift. Fixed by holding all edits until
every seat reports, which is what this effort did across all four waves.

**`pkill -f "<logname>"` matched my own command line**, killing my own shell and
a static-analysis run (exit 144). Replaced with pid loops that exclude `$$` and
`$PPID`.

### 3.2 Verification that could not detect the failure it was checking for

**A `grep -c` that could never match.** To confirm the "both forms" fix had
landed I ran `grep -c "both the guard and .pragma once"` and got `0`, which I
read as "closed". The phrase spans a line break in the file, so a line-oriented
grep could not match it whether the fix had applied or not. The check was
incapable of failing informatively, and it certified an edit that had never
run. Replaced with Python containment checks over the whole file text, which
are line-break-safe, and with paired positive/negative assertions so a missing
edit shows up as a failed *positive* check rather than a passed negative one.

**Reporting a filtered result as if the recorded command produced it.**
`CONTR-0009`'s resolution test names the command
`grep -rn "glb_gridfcs_dict" nrpy/infrastructures/Dendro`. In my own checking I
ran that command *plus* `grep -v ">>>"` to drop doctest lines, saw one hit, and
wrote "returns exactly one line" into both the register row and my report to
the user. Run as recorded, it returns 14 `clear()` hits, thirteen of them
doctest fixtures. Two separate seats in two separate waves blocked on this, the
second time after I had already "fixed" it once - my first correction scoped
`del` and `pop` but left a new absolute about `clear`.

**A header sweep with the wrong scope.** Checking the "both forms" claim I swept
`find Dendro-GR -name '*.h' -o -name '*.hpp'` and found 336 headers, 40 of them
carrying both a guard and `#pragma once` - which would have made the disputed
sentence defensible. Two seats reported 81 headers and zero carrying both. The
seats were right: my sweep included `build*/` and `_deps/` trees, so 255 of my
336 were build artifacts and vendored third-party copies rather than
Dendro-GR's own source. Restricted to non-build trees: 81 headers, 49
guard-only, 10 `#pragma once`-only, 22 neither, **0** carrying both.

### 3.3 False claims I wrote into the tree

Each of these was a statement in code, a docstring, an emitted comment or a KB
page that a seat then disproved. All are fixed; they are recorded because the
pattern - asserting a stronger version of something true - recurs.

| Claim I wrote | Why it is false | Caught by |
| --- | --- | --- |
| The fCCZ4 leaf "describes the split the same way", and the duplication "is recorded there" | Neither is true; see 1.1 | Wave 1, verifier |
| BHaH "branches inline on `enable_fCCZ4` at four sites" | Four *occurrences*; two are branches, the others are the parameter declaration and its `:param` line. `inconsistencies.md` correctly said "four mentions"; I converted mentions to branches | Waves 1 and 2 |
| "No `glb_gridfcs_dict` deletion exists anywhere in the Dendro infrastructure" | `trusted_capture.py:70` calls `gri.glb_gridfcs_dict.clear()` | Wave 2, standards |
| "the subpackage's only registry clear is ... `trusted_capture`" (the first correction) | Thirteen further `clear()` calls exist as doctest fixtures across five modules | Wave 3, all three seats |
| "A padding wired from somewhere other than the builder's record fails here" (emitted C++) | `required_padding = 2` at `fd_order 4` satisfies `2 >= 2` and passes. The deleted floor was also 2 with KO off, so no coverage was lost - only the comment overclaimed | Waves 1 and 3 |
| "every caller in this repository reads it back with `CFunction_roles.required_padding()`" | The module's own doctest, sixteen lines below, passes a literal `3` | Wave 2, integration |
| "correct only for the families that have a C-code path today" | Such a floor is wrong *within* the C-code-path families: at `fd_order` 4 with KO off it gives 2 against a true reach of 3, which is the configuration both applications ship | Wave 3, verifier |
| `dfullupD`/`dfulldnD` "have stencils but no such path, and would break it" (emitted C++) | Those families cannot reach emitted C at all - `c_codegen` defaults to `C_CODEGEN_DERIVATIVE_FAMILIES`, and forcing one through raises in `proto_FD_operators_to_sympy_expressions`. The sentence restored a live-hazard reading that the parallel docstring sentence denies | Wave 3, integration and standards |
| The host's headers "use both the guard and `#pragma once`" | True of zero headers. The deleted design record said "use both **forms**", meaning the corpus contains both styles; the paraphrase dropped the word and turned a true statement about a corpus into a false one about each header | Wave 3, all three seats |
| "each pair of private builders sharing a long verbatim tail" | Well supported for the constraint pair (33 + 20 + 12-line runs, 8-line identical tail); the RHS pair's duplication is scattered, longest contiguous run 11 lines | Wave 3, integration and standards |

Two further mechanical misses in the same class:

- **A dangling symbol my own rename created.** S2 deleted
  `BSSN_CONSTRAINTS_EVAL_ALL_BLOCKS_CFUNCTION`, and `CONTR-0008`'s
  claim-evidence block still cited it. I renamed constants without sweeping the
  KB for citations of the old names. All three wave-one seats found it
  independently, which is how obvious it was.
- **A catalog source count left stale.** I12 removed the ADR citation from the
  naming leaf's `Sources` list, taking it from seven entries to six. I updated
  that catalog row's date and its query keyword and left the count at `7`.

### 3.4 Incomplete work reported as complete

- **S1 at three sites.** I told the seats that `desc`, `cfunc_type`, `name`,
  `subdirectory` and `includes` were named locals "at all 21 production sites".
  At `CodeParameters.py`'s three sites only `desc`, `subdirectory` and
  `includes` were hoisted; `cfunc_type`, `name`, `params` and `body` were still
  inline. Two seats caught it, one calling it blocking and one nonblocking.
- **Reconciliation dates.** I bumped the dates on six touched pages and left
  `wiki/catalog.md`'s and `wiki/source-map.md`'s own header dates and their own
  catalog rows at `09-06-2026`, although the candidate edited both.
- **A propagated defect.** While correcting the "both forms" wording elsewhere I
  wrote the identical defective construction into `raw/SOURCES.md`'s coverage
  note, so a batch intended to fix the claim in one place doubled it.

### 3.5 A disproved rationale from the earlier effort

`trusted_capture.reset_generation_state`'s docstring originally justified
clearing the equations factories' memos with a cache-key argument. That was
disproved: the factories *do* rebuild on a changed conformal factor, because
`EvolvedConformalFactor_cf` is part of the memo key. The accurate reason, which
the docstring now records, is that the factory constructors register the evolved
state, so clearing the gridfunction registry without clearing the memos leaves
the next build with nothing registered.

---

## 4. Reviewer misses and disagreements

Seats are not infallible either, and the record should say so.

- **A miss, caught by a peer.** The second wave's verifier seat examined the
  `:param required_padding:` provenance sentence and explicitly accepted it
  ("the `:param`'s 'every caller in this repository reads it back' holds"). The
  integration seat in the same wave found the counterexample in the module's own
  doctest. Two seats, same sentence, opposite conclusions; the one that opened
  the file's doctests was right.
- **A severity disagreement, not an error.** S1's incompleteness at three sites
  was blocking to the second wave's integration seat and nonblocking to its
  standards seat. Resolved by taking the stricter reading.
- **A scope disagreement, resolved in the seats' favour.** The header-count
  discrepancy in 3.2: 336/40 (mine) against 81/0 (theirs). Theirs was the
  meaningful scope.
- **Correct arithmetic, wrong attribution, mine.** I reported the changed
  modules' doctest counts as 9/10/3 when they are 10/9/3. Caught by a seat in
  passing. No consequence beyond the report.

---

## 5. Declined, not false

These are real observations that were deliberately not acted on. They are
recorded so nobody mistakes them for oversights, and so they are not
re-litigated without new evidence. A standards seat assessed each decline
against the deciding authorities and found each defensible.

- **Removing the keyword-only `*` from 17 remaining Dendro functions** whose
  keyword-only parameters are optional profile knobs. No rule in `CLAUDE.md`,
  `coding_style.md` or the routed KB mandates positional parameters; S7 named
  only the name-last problem; and removing them would move argument-order risk
  into parameters where a wrong-order call is much harder to spot than the
  `solver_stem` one was.
- **`cmake_helpers.py`'s stem-versus-prefix argument order** across three
  neighbouring emitters. Unchanged from baseline and outside the eleven.
- **Hoisting `params=build.block_params` and `body=build.block_body`** into
  locals at the 21 registration sites. They are attribute reads off frozen
  build records. `coding_style.md:745` does list `params` and `body`, so this
  decline rests on scope rather than on the rule's text, and it leaves two
  patterns in the package: three sites name all six fields, eighteen name four.
- **`cfunc_type` in the two `register_CFunction` calls inside
  `CFunction_roles.py` doctests.** Outside the 21 production sites; the
  parameter defaults to `"void"`; the rule governs registration functions
  rather than doctest fixtures.
- **A `CONTR-0011` row** filing the private-builder-per-formulation split as a
  divergence `new-infrastructure-conformance.md` does not permit. The BSSN leaf
  names that rule inline and records the cost instead. Three seats accepted the
  inline treatment as satisfying the design-record obligation; the seat that
  proposed the row did not block when it was declined.
- **`KO_ENABLED` has no emitted consumer** after I14 removed the formula that
  read it. Kept as the only place the generated project records which
  dissipation the kernel was lowered with. A seat confirmed the precedent -
  `NUM_CODE_PARAMETERS` and `NUM_AUX_GFS` are already emitted without
  consumers - and found no rule requiring one.
- **The exact stencil reach is no longer asserted in emitted C++.** I14's
  emitted case is a one-sided bound. The exact reach is pinned instead by owner
  doctests (`rhs_eval.py` asserts `padding == 3` for both formulations) and by
  the emitted parameter file, and the host's real block padding is still gated
  against `REQUIRED_PADDING` by `Ctx::startup_checks`. A seat judged the
  residual loss to be defence-in-depth against a mis-wired padding, with nothing
  left unpinned.

---

## 5a. Added by the fresh review invocation

A fourth invocation of three fresh seats reviewed the candidate after the
wave budget was exhausted. Two blocked, one accepted. It produced six more
entries for this document.

### 5a.1 A prescribed remedy that was itself wrong in context

Three wave-three seats independently prescribed "use both forms" for the
upstream header sentence, and I applied it verbatim. A fresh seat then showed it
is **false under its nearest antecedent**: two sentences earlier the same
paragraph reads "Both simple names such as `BHAHAHA_HEADER_H` and legacy
double-underscore forms such as `__SIMD_INTRINSICS_H__` appear in the codebase",
so "both forms" naturally means the two *guard-macro naming* forms - and the
upstream checkout contains no double-underscore guard at all. The deleted design
record's original sentence was unambiguous because its paragraph never mentioned
macro-naming forms; I12's relocation created the ambiguity, and the prescribing
seats did not read the surrounding paragraph.

**Lesson.** A remedy prescribed by a reviewer is still my claim once I write it.
Three seats agreeing on a wording is not evidence the wording is true in its
destination. The fix now names both mechanisms explicitly and cites one header
of each kind.

### 5a.2 Three wrong counts in my own review brief

- I wrote that "three `CodeParameters.py` sites name all six fields, the other
  eighteen name four". An AST audit shows the keyword set passed to
  `register_CFunction` is **uniform** across all 21 sites -
  `(subdirectory, includes, desc, cfunc_type, name, params, body)`. The real
  distinction is narrower and different: three sites pass named locals for all
  five of `desc`/`cfunc_type`/`name`/`params`/`body`, the other eighteen pass
  locals for three and attribute reads for `params` and `body`. I conflated
  keyword sets with named-local counts and inflated both numbers.
- I wrote "62 files per project". It is 31 per project, 62 across both.
- I wrote "all 23 Dendro top-level modules and all 5 `general_relativity`
  modules", implying 28. There are 18 top-level plus 5, and 23 runnable in
  total.

None changed a conclusion, but a brief is evidence handed to reviewers, and two
seats spent effort correcting it.

### 5a.3 Doctest values I guessed instead of measured

Closing a real coverage gap, I added a doctest pinning the upwinded stencil
reach at orders 4 and 6 and asserted `3` and `4` for
`padding_from_derivative_operators([cf_dupD[1]], "unset", N)`. Both failed.
`padding_from_derivative_operators` raises unless the expressions reach every
axis, so a single-component upwind expression raises rather than returning a
reach. Measured with a three-axis expression the values are 3 and 4 as expected,
but the call shape I wrote could never have produced them.

**Lesson.** In the middle of an exercise about unverified claims I wrote two
unverified numbers. Measure, then assert - including in a doctest whose whole
purpose is to pin a measurement.

### 5a.4 A mechanism claim that named the wrong mechanism

I wrote that `dfullupD`/`dfulldnD` "have stencils and no C-code path and so
cannot appear in an emitted kernel today". That conclusion is too strong:
full-upwind lowering is unsupported, but the mechanism is not a rejection: `extract_list_of_deriv_var_strings_from_sympyexpr_list`
defaults to `C_CODEGEN_DERIVATIVE_FAMILIES` and drops the other two families
with an explicit `pass`, so `c_codegen` classifies such a symbol as an ordinary
one and emits it. A seat demonstrated it: `c_codegen([uu_dfullupD0], ...)`
returns `const double x = uu_dfullupD0;` with no error. The `ValueError` a
neighbouring module comment advertises never fires on that path. The docstring
now says the kernel emits an undeclared identifier and fails to compile.

### 5a.5 An emitted-constant claim with no consumer

The `KO_ENABLED` docstring said the constant is recorded "so a reader of the
generated header, **and the host**, can tell which dissipation the kernel was
lowered with". After I14 removed its last consumer, nothing reads it - no
emitted code, no self-test, and no host, since a real-host build is a named
deferral gate on the validation leaf. The docstring asserted a consumer the
KB's own gates say does not exist yet.

### 5a.6 A KB marker that failed lint, and a correction to section 2

Opening `CONTR-0011` I wrote the inline `Claim status: contested;
contradiction: CONTR-0011.` markers without the backlink the register contract
requires. `python tools/kb_lint.py` **caught it**:
"marker lacks backlink for CONTR-0011", on both affected leaves.

This corrects section 2 of this document, which said the linter "checks
structure and link resolution" and does not check accuracy. That is true of
symbol existence and of factual claims, but the linter does enforce parts of
the contradiction contract, including marker backlinks. It is a narrower gap
than I described.

### 5a.7 A disagreement between seats, resolved on text

The fresh verifier seat wrote "I have no clause in `wiki/SCHEMA.md` that
*requires* a `CONTR-0011` row for a code-versus-governance divergence"; the
fresh standards seat quoted `wiki/SCHEMA.md:213`, "normative KB or contributor
rule: current owning governance or configuration decides; **implementation
divergence opens a contradiction**". I read the clause: it says what the
standards seat says it says. The row was opened. Recorded here because the
earlier decision to decline that row - accepted by a whole wave - rested on the
verifier's reading, and the clause was there the whole time.

---

## 5b. Added while closing the round

### 5b.1 The reconciliation-date omission, third recurrence

Adding the `CONTR-0011` marker to `fccz4-application-wiring.md` edited a page
whose header date and `wiki/catalog.md` row I then left at `09-06-2026`. This is
the third instance of the same failure in this effort: section 3.4 records
bumping six pages and leaving `catalog.md`'s and `source-map.md`'s own dates,
and a wave before that found the naming leaf's source count left at 7 after its
`Sources` list dropped to 6.

The pattern is specific and mechanical: I treat "edit the page" and "update the
page's metadata" as separate steps and complete only the first. `kb_lint` does
not check date freshness, so nothing catches it. The fix is to treat any page
edit as a three-part unit - content, header date, catalog row (date, status,
source count) - and to check the row every time, not only when the edit felt
metadata-shaped. I caught this one myself before the wave, which is the only
reason it is not a fourth blocking finding.

### 5b.2 A resolution test that could not detect half of what it covered

`CONTR-0011` covers two divergences. Its Resolution test named
`grep -n "^def " nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py`,
and its second clause - "each public `register_CFunctions_*` builds and
registers in one routine" - is about a pairing that spans every kernel family in
the subpackage, not one file. The inspection as written could not observe the
resolution of the thing it claimed to cover.

This is the same class as the `CONTR-0009` failure in section 3.2, where the
recorded command did not produce the result the row stated. There the recorded
command was too narrow because I had silently filtered its output; here it was
too narrow because I scoped it to the file I happened to be looking at. Both
come from writing the test after the claim instead of running it first. The
inspection now spans `general_relativity/*.py` and names both structures.

### 5b.3 Escalating a question the repository had already answered

I put two questions to a review wave as things I "could not settle alone":
whether an active `contested` row obliges a page-status change on the pages it
names, and whether one row or two should cover two divergences of the same
class. Both were answered by `wiki/contradictions.md` itself, in rows that had
been sitting in the file the entire time:

- CONTR-0001 and CONTR-0003 leave every affected page `confirmed`; CONTR-0002
  downgrades only the page whose sole generated interface the conflict
  invalidates. So an active row does not by itself lower a page's status.
- CONTR-0004, 0005 and 0006 are three separate rows for three discrepancies in
  one source paper, each with a different deciding passage and resolution test -
  which is the test for when to split, and both of mine share theirs.

Finding this took one command each once I looked. The standing rule for the user
is not to ask them design questions the tree already answers; the same rule
applies to reviewers, whose time is the same resource. Read the neighbours
before escalating a governance question - the KB is precedent, not just rules.

### 5b.4 An honest note on the shape of this round

Across four completed waves and twelve seat reports, **not one blocking finding
landed on the code or on the generated output**. All twelve were documentation
and governance-record defects, and all but two were in text this effort itself
wrote while fixing earlier defects of the same kind. The eleven substantive
remedies were verified sound at every wave and are provably output-neutral
apart from I14's intended change.

That asymmetry is the round's real result. The code changes were small,
mechanical and well-covered by existing oracles; the risk was concentrated
entirely in the prose written *about* them, where no automated gate exists. A
future round of this shape should budget its review effort accordingly - and
should expect the correction batches, not the original change, to be where new
defects enter.

---

## 6. Patterns worth carrying forward

1. **Verify the premise, not just the conclusion.** The single most expensive
   error here was trusting `inconsistencies.md`'s claim about what another page
   said. One sentence of checking would have prevented a retracted remedy and
   two rounds of register corrections.
2. **Scope governance claims to their evidence.** "No X exists anywhere"
   cost three blocking findings across two waves. Absolute claims require
   exhaustive evidence within an explicit scope; record the actual check.
3. **A recorded test must be the command that was run.** If a filter is needed
   to get the stated result, the filter belongs in the recorded command.
4. **A negative check that cannot distinguish "fixed" from "never applied" is
   not an adequate regression check.** Show that the check distinguishes the
   relevant good and bad states. Add positive assertions where needed, and
   make text checks robust to line breaks.
5. **Batch edits must report per-edit status and fail loudly.** An aborting
   script that has already written some files is worse than one that writes
   none, because the partial state looks like success.
6. **Renaming a symbol obliges a KB citation sweep**, not just a code sweep.
7. **`kb_lint` passing means the structure is intact, nothing more.** Factual
   accuracy in `wiki/` has no automated gate; it needs a reader.
8. **Scope every measurement before quoting it.** Build directories and
   vendored dependency trees are not the project.
