# BSSN Application Wiring

> Explain the Dendro BSSN builders, what they reuse from the fCCZ4 profile, and what adding a second formulation proved about the generic layer. · Status: provisional · Last reconciled: 09-05-2026
> Up: [Dendro](index.md)

## Summary

BSSN is the second formulation lowered through the Dendro infrastructure. It
lives under `general_relativity/BSSN/`, matching how BHaH puts variant families
in a subdirectory (`Kasner/`, `TOVola/`, `TwoPunctures/`) beside the default at
the package top level. Two modules are new — the right-hand side and the
constraint diagnostics — and everything else is shared with fCCZ4 unchanged:
the initial data, the algebraic projection, the whole generic layer, and the
formulation-agnostic lowering in `kernel_lowering`.

## Detail

### What the builders do

`BSSN/rhs_eval.py` assembles the evolution system the way ETLegacy and BHaH do:
the non-gauge equations from the cached `BSSN_RHSs` object, the lapse and shift
from `BSSN_gauge_RHSs` added to a *copy* of its dictionary, Kreiss-Oliger terms
through the shared `add_KreissOliger_dissipation_terms` helper with
`include_Theta_fCCZ4=False`, and the upwind control vector as the rescaled
shift `betaU[i] = vetU[i] * ReU[i]`. It then lowers all of it through
`kernel_lowering`, asserts the 24-field bijection against the registry, and
records the padding and the upwind control set.

`BSSN/diagnostics.py` emits the Hamiltonian constraint and the three momentum
constraint components from the established `BSSN_constraints` projector.

### The DIAG-before-projector ordering, and why

`BSSN_constraints` registers `H`, `M`, `LAMBDA_CONSTRAINT` and `MU` into the
**AUX** group when it is constructed. Dendro's diagnostics contract is the
**DIAG** group, which is the settled infrastructure convention: BHaH registers
31 DIAG gridfunctions across its wave-equation, elliptic and GR diagnostics,
while AUX appears once. Every registration in the projector is guarded by
`if <name> not in gri.glb_gridfcs_dict`, so the Dendro builder registers the
names it writes as DIAG *before* constructing the projector, which makes the
projector's own registration a no-op. No change to the shared equations module
and no change to the Dendro generic layer.

Claim evidence:
- Claim: `BSSN_constraints` registers its diagnostic gridfunctions into the AUX group, each guarded by an existence check, so a caller that registers those exact names first determines their group; the Dendro BSSN builder uses this to keep them in DIAG.
- Role: descriptive behavior
- Deciding authority: [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), the `group="AUX"` registrations and their `not in gri.glb_gridfcs_dict` guards
- Corroboration: [diagnostics.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN/diagnostics.py), `build_diagnostics` step 1
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3, GCC 13.3.0, CMake 3.28.3, OpenMPI 4.1.6; backend=Dendro; precision=double; GPU=not-applicable; restart=not-applicable; distributed=1 and 2 MPI ranks; error_path=not-run; options=--fd-order 4 --no-ko; date=09-05-2026`

### Dendro's own vocabulary

The emitted names follow Dendro-GR's BSSN solver rather than NRPy's: the solver
directory is `BSSN_GR`, CMake variables carry the `BSSN_` prefix, the object
library is `bssn_common`, the executable is `bssnSolver`, the namespace is
`bssn`, and the context source is `bssnCtx.cpp`. The emitted CFunctions are
`bssn_rhs*` and `bssn_constraints*`, matching `bssneqs.cpp` and
`bssn_constraints.cpp` upstream.

### What the port found

Adding the second formulation was the cheap test of whether the abstraction
exists, and it found three real defects that a single-formulation tree could
not expose:

- The "shared" ADM-to-evolved conversion asserted that exactly one evolved
  field was left undefined — true only for fCCZ4's Z4 scalar. The rule is now
  that *at most* one constraint scalar is zeroed.
- The diagnostics accessed-set intersected raw free symbols with the registry,
  so a field that a kernel only *differentiates* was never bound. The BSSN
  momentum constraint takes derivatives of `lambdaU`, and the emitted kernel
  read a pointer nothing declared. `kernel_lowering.base_gridfunction_of` now
  resolves derivative symbols back to their field.
- `naming.aux_pointer` had been deleted as dead code in an earlier review
  round. It was dead only because the tree had one formulation.

Nothing outside `general_relativity/` needed a change, which is the result the
port was run to obtain.

## Sources

- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN/rhs_eval.py) - `bssn_rhs_expressions`, `build_bssn_rhs`, `register_CFunctions_rhs_eval`
- [diagnostics.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN/diagnostics.py) - `build_diagnostics`, `register_CFunctions_diagnostics`
- [kernel_lowering.py](../../../nrpy/infrastructures/Dendro/kernel_lowering.py) - `base_gridfunction_of`, `accessed_gridfunctions`, `padding_from_operators`
- [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py) - `BSSNRHSs`
- [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py) - `BSSN_gauge_RHSs`
- [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py) - `BSSNconstraints`
- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - generation entry point and command-line profile

## See Also

- Parent: [Dendro](index.md)
- Depends on: [BSSN Family](../../equations/general-relativity/bssn-family.md)
- Implements: [Gridfunctions, Naming, And Loops](gridfunctions-naming-and-loops.md)
- Contrasts with: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
- Validated by: [Validation, Host Mock, And Deferral Gates](validation-host-mock-and-deferral-gates.md)
