# ADR-001: Dendro fCCZ4 state names, order, and native semantics

- Date: 09-03-2026
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
3. Runtime physics parameters are the used registered `CodeParameter`
   objects only; no Dendro physics table duplicates them.
4. A kernel is a registered `CFunction` plus non-authoritative Dendro
   role metadata; no second body registry exists.

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
  sentinel tests, awaits the Dendrolib pin and capability proof recorded in
  `dendrolib_pin.json` and `dendrolib_capabilities.json`. Only the
  `dendro_mock.h` struct exists. This ADR records the deferral; adapter
  signatures stay frozen until the pin lands.
- The checkpoint ABI and physical boundaries are separate qualified
  profiles, gated on the same pin.
