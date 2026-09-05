# ADR-001: Dendro fCCZ4 state names, order, and native semantics

- Date: 09-03-2026
- Status: approved
- Scope: `nrpy/infrastructures/Dendro/` (whitepaper section 17, Phase 0)

## Decision

1. The 25 fCCZ4 EVOL gridfunctions keep their exact NRPy names and the
   `GridFunction.gridfunction_lists()` order. No BSSN aliases
   (`cf` stays `cf`, never `chi`; `lambdaU` never `Gt`;
   `Theta_fCCZ4` never `Theta`).
2. The state stores native `hDD` (reference-metric conformal-metric
   perturbation, Cartesian: $\bar{\gamma}_{ij} = \delta_{ij} + h_{ij}$)
   and native evolved `lambdaU` (fCCZ4 conformal connection quantity, not
   a BSSN contracted connection). No full-metric field replaces `hDD`.
3. Runtime physics parameters are the used registered `CodeParameter`
   objects only; no Dendro physics table duplicates them.
4. A kernel is a registered `CFunction` plus non-authoritative Dendro
   role metadata; no second body registry exists.

## Consequences

- Generated machine identities preserve registered names byte-for-byte
  (`in_<name>`, `rhs_<name>` are reversible decorations).
- `state_schema_hash` changes on any intentional rename/reorder, and
  incompatible checkpoints are rejected.

## Deferred (not dropped)

- I3-2 `BlockGeometry` adapter proof (whitepaper section 7.4): the single
  auditable host function normalizing `pmin_padded`/`component_offset`
  plus two-block/offset sentinel tests await the I0-1 Dendrolib pin and
  capability proof. Only the `dendro_mock.hpp` struct exists. This ADR
  records the deferral; adapter signatures stay frozen until the pin lands.
- PR5 `used_codeparameters` scan still uses substring matching
  (`Dendro/general_relativity/rhs_eval.py`); the freeze-side closure is
  token-aware. Unifying PR5 on the shared helper is PR5-scoped follow-up
  (frozen R2 scope excludes the PR5 builder).
