# Review findings and dispositions

Date: 09-07-2026. Branch: `dendro-infra`. The earlier candidate is committed as
`4d97621e` (39 files, +696/-454). This document records the subsequent cleanup.

Untracked working document, like `progress.md`, `inconsistencies.md` and
`false_directions.md`. Never committed.

## Status and scope

The previous version reported three `ACCEPT`s in an eighth review wave. That
is a historical report, not independently reproduced acceptance of this
cleanup. The dispositions below replace the former open list and distinguish
factual corrections from optional clarity edits. They describe the cleanup's
changes, not a claim that all Dendro work is finished.

The subsequent requested implementation resolves CONTR-0011: each RHS and
constraint module has one builder and one registrar, with unchanged registered
results for both shipped formulations. The real fCCZ4 host now builds at the
pinned Dendro-GR/Dendrolib revisions and passes two-rank transport, parameter,
collective-failure, and 100-step fixed-mesh Minkowski checks. See
`nrpy/infrastructures/Dendro/tests_infra/README.md` for reproduction and limits.

Remaining work: general physical boundaries, remeshing/LTS, checkpoint/restart,
output integration, GPU/thread qualification, BSSN real-host execution, and
real-host CI. The dispositions below describe the earlier cleanup; their
historical source inspections alone did not establish runtime integration.

## 1. Resolved: mark the deleted quotation as historical

In `wiki/contradictions.md`, CONTR-0010's Notes now identify the then-current
"That split is Dendro's own" sentence as superseded wording and use past tense.
The current fCCZ4 leaf no longer contains that sentence.

The former list overstated this defect's priority: the quotation was not the
only explanation for the BSSN/fCCZ4 status split. CONTR-0011's Notes and the
fCCZ4 Summary already explain the withdrawn pure build/register objection.

## 2. Resolved: remove the false full-upwind comparison

`nrpy/infrastructures/Dendro/constants_h.py` now states that `dfullupD` and
`dfulldnD` reach `fd_order`, without comparing that reach to single-point
upwinding. At orders 2/4/6, centered reach is 1/2/3, single-point upwinded reach
is 2/3/4, and both full-upwind reaches are 2/4/6.

The original proposed remedy was insufficient: "rather than one point past
the centered radius" is also false at order 2. Both comparative clauses are
removed.

## 3. Clarified: identify the BHaH RHS module precisely

CONTR-0010's Authority decision now uses "BHaH's only module holding both
formulations' RHS systems", matching its Claim. The old "two-formulation
module" was understandable in context; this is a precision edit, not another
proven contradiction.

## 4. Resolved: complete the cited upstream coverage

`raw/SOURCES.md` now includes root Dendro-GR `CMakeLists.txt`'s Dendrolib tag
ordering and `BSSN_GR` namespace declarations, alongside the four previously
listed files. `wiki/source-map.md` records the corresponding host-source
scope and dependent pages. These source observations do not prove real-host
integration.

## 5. Resolved: pin full-downwind reach

`nrpy/finite_difference.py` now declares `uu_dfulldnD` beside `uu_dfullupD`
and checks `stencil_reach_per_axis([uu_dfulldnD[1]], "unset", 4) == (0, 4, 0)`.
This closes the missing full-downwind reach assertion. It does not test the
coefficient signs or stencil orientation.

The old assertion that no emitted kernel can contain a full-upwind operator
was too strong. Supported Dendro lowering does not implement these families;
`c_codegen` can treat such a symbol as ordinary and emit an undeclared
identifier. This is not a guaranteed generator rejection.

## 6. Clarified: name the reach owner and maximum

`block_kernel_helpers.py` now repeats the explicit
`nrpy.finite_difference.stencil_reach_per_axis` reference and names the
maximum across axes. The reference was already present earlier in the
docstring; it had not disappeared entirely. The companion reach paragraph in
`finite_difference.py` is reflowed. These are readability edits.

## 7. Simplified: group RHS suffix comments

`general_relativity/rhs_eval.py` now has one shared comment, says "both
formulations", and retains the LTS flat-layout explanation without redundant
per-constant labels. This follows the three Dendro siblings as a local style
choice. The infrastructure-conformance three-example rule counts established
infrastructures, not three modules inside Dendro; the former rationale
misapplied that rule.

## 8. Simplified: remove the rejected-draft detour

The BSSN leaf's Summary now states the adopted module layout directly and then
describes the outstanding intra-module divergence. Its opening already gave
the disposition; removing the earlier-draft rationale makes that clearer.
No new claim about a repository-wide absence of subpackages is needed.

## 9. Resolved: remove the unreproducible namespace count

`wiki/infrastructures/new-infrastructure-conformance.md` retains the host's
lowercase `bssn` namespace observation and drops "14 occurrences". The count
depends on matching and file scope and is unnecessary for the naming rule.
The page's evidence block and catalog reconciliation date are updated.
Earlier counts in `progress.md` remain historical measurements, not current
claims.

## 10. Clarified and normalized: prose versus emitted spelling

The four prose occurrences of `behaviour` are normalized to `behavior`: the
`nrpy/examples/dendro_bssn.py` module docstring and three Dendro leaves.
The nine `centre` occurrences in `solver_context.py` remain unchanged. They
are emitted text, including comments as well as identifiers; describing all
nine as identifiers was inaccurate. The earlier broad reason for declining
all spelling edits did not apply to the four prose occurrences.

## Related corrections in the working records

`progress.md` and `false_directions.md` now begin with a current checkpoint
superseding stale working-tree and contradiction-status descriptions while
preserving the historical sections. `false_directions.md` also corrects the
full-upwind emission mechanism and scopes its general lessons: absolute
claims need adequate scoped evidence, and regression checks must distinguish
the relevant good and bad states. Positive assertions are useful where
needed, not mandatory partners for every negative assertion.

The existing statements that single-point upwinded and Kreiss-Oliger families
reach one point beyond centered stencils remain correctly scoped. The
full-upwind pair has a different reach contract and no supported C lowering.
