# ADR-001: Dendro fCCZ4 state names, order, and native semantics

- Date: 09-03-2026
- Amended: 09-06-2026 - decision point 3 corrected to describe the shipped
  whole-registry emission.  An earlier draft of this branch did compute a
  used-parameter closure; that design was discarded when the infrastructure was
  rebuilt on direct registry reads, and the ADR is only now catching up.
- Status: approved
- Scope: `nrpy/infrastructures/Dendro/`

## Decision

1. The 25 fCCZ4 EVOL gridfunctions keep their exact NRPy names and the
   `GridFunction.gridfunction_lists()` order. No BSSN aliases
   (`cf` stays `cf`, never `chi`; `lambdaU` never `Gt`;
   `Theta_fCCZ4` never `Theta`).
2. The state stores native `hDD` (reference-metric conformal-metric
   perturbation, Cartesian: $\bar{\gamma}_{ij} = \delta_{ij} + h_{ij}$)
   and native evolved `lambdaU` (fCCZ4 conformal connection quantity, not
   a BSSN contracted connection). No full-metric field replaces `hDD`.
3. Runtime physics parameters are the registered `CodeParameter` objects,
   emitted whole rather than filtered to a use closure; no Dendro physics
   table duplicates them.
4. A kernel is a registered `CFunction` plus non-authoritative Dendro
   role metadata; no second body registry exists.

5. Emitted headers carry traditional `#ifndef`/`#define` header guards, as
   `coding_style.md` section 4 requires and as NRPy's other infrastructures
   emit. The macro is derived from the file name of each generated header,
   which already carries the solver stem, so two generated solvers in one
   Dendro-GR checkout cannot collide. Dendro-GR's own headers use both forms, so nothing about the host
   argued for `#pragma once`.

## Consequences

- Generated machine identities preserve registered names byte-for-byte
  (`in_<name>`, `rhs_<name>` are reversible decorations), so a generated
  identifier always traces back to exactly one registered gridfunction.
- An intentional rename or reorder therefore changes the emitted `EvolVar`
  enum and the component order the generated state header declares. The
  checkpoint ABI is a separate qualified profile and is deferred, so nothing
  currently rejects an incompatible checkpoint.

## Deferred (not dropped)

- The `BlockGeometry` adapter proof: the single auditable host function
  normalizing `pmin_padded`/`component_offset`, plus two-block and offset
  sentinel tests. The Dendrolib commit and the capability proof it waited on have
  landed: the generated project's `CMakeLists.txt` records that commit and
  compares it against the tag a real-host build configures, and the proof is
  recorded on the Dendro validation page. What remains deferred is a build
  against a real Dendro-GR checkout. Only the
  `standalone_host/dendro_standalone_host.h` struct exists so far, and the
  adapter signatures stay frozen until that build runs.
- The checkpoint ABI and physical boundaries are separate qualified profiles,
  gated on that same real-host build.
