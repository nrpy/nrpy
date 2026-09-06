# BSSN Application Wiring

> Explain the Dendro BSSN builders, what they reuse from the fCCZ4 profile, and what adding a second formulation proved about the generic layer. · Status: provisional · Last reconciled: 09-06-2026
> Up: [Dendro](index.md)

## Summary

BSSN is the second formulation lowered through the Dendro infrastructure, and
it lives under `general_relativity/BSSN/`. Two modules are new — the right-hand
side and the constraint diagnostics — and the initial data, the det(gammabar)/tr(Abar)
enforcement and the generic layer are shared with fCCZ4.

The layout is a deliberate divergence, not an imitation. NRPy's one existing
two-formulation infrastructure is BHaH, which lowers BSSN and fCCZ4 through a
single `general_relativity/rhs_eval.py` on an `enable_fCCZ4` boolean; BHaH's
subdirectories (`Kasner/`, `TOVola/`, `TwoPunctures/`, `psi4/`) are diagnostics
and initial-data providers, not second implementations of a top-level module's
role. Dendro diverges because its two formulations differ in evolved-field
count, in constraint set and in which shared equations module supplies the
expressions, so a single module would branch on formulation in every one of
those places. That is one instance against one instance, which the conformance
page rates as weak evidence either way; collapsing onto BHaH's shape remains a
live option.

## Detail

### What the builders do

`BSSN/rhs_eval.py` assembles the evolution system the way ETLegacy and BHaH do:
the non-gauge equations from the cached `BSSN_RHSs` object, the lapse and shift
from `BSSN_gauge_RHSs` added to a *copy* of its dictionary, Kreiss-Oliger terms
through the shared `add_KreissOliger_dissipation_terms` helper with
`include_Theta_fCCZ4=False`, and the upwind control vector as the rescaled
shift `betaU[i] = vetU[i] * ReU[i]`. It then lowers all of it through
`block_kernel_helpers`, asserts the 24-field bijection against the registry, and
records the padding and the upwind control set.

`BSSN/constraints_eval.py` emits the Hamiltonian constraint and the three momentum
constraint components from the established `BSSN_constraints` factory.

### The DIAG-before-factory ordering, and why

`BSSN_constraints` registers `H`, `M` and `LAMBDA_CONSTRAINT` into the **AUX**
group when it is constructed, each guarded by an existence check. `MU` is
different: the factory registers it only under the `register_MU_gridfunctions`
CodeParameter, which defaults to `False` and which no Dendro module sets, so
nothing competes for those names. Dendro's diagnostics contract is the
**DIAG** group, which is the settled infrastructure convention: BHaH registers
31 DIAG gridfunctions across its wave-equation, elliptic and GR diagnostics,
while AUX appears once. Every registration in the factory is guarded by
`if <name> not in gri.glb_gridfcs_dict`, so the Dendro builder registers the
names it writes as DIAG *before* constructing the factory. That is
load-bearing for `H`; `MU0`-`MU2` are DIAG simply because this builder is their
sole registrant. The factory still registers `M` and `LAMBDA_CONSTRAINT`,
which this kernel does not compute, so the builder removes those two
afterwards — restricted to newly added AUX names, because the factory's
construction also pulls in the evolved state.
Without that the generated state header would advertise two variables no kernel
writes and no vector backs. No change to the shared equations module.

Claim evidence:
- Claim: `BSSN_constraints` registers `H`, `M` and `LAMBDA_CONSTRAINT` into the AUX group, each guarded by an existence check, so a caller that registers `H` first determines its group; it registers `MU` only under the `register_MU_gridfunctions` CodeParameter, which defaults to `False` and which no Dendro module sets. The Dendro BSSN builder pre-registers the names it writes and afterwards deletes the newly added AUX names it does not write, leaving `NUM_AUX_GFS = 0` in the emitted state header.
- Role: descriptive behavior
- Deciding authority: [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), the `group="AUX"` registrations and their `not in gri.glb_gridfcs_dict` guards
- Corroboration: [constraints_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN/constraints_eval.py), `build_constraints_eval` step 1
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3, GCC 13.3.0, CMake 3.28.3, OpenMPI 4.1.6; backend=Dendro; precision=double; GPU=not-applicable; restart=not-applicable; distributed=1 and 2 MPI ranks; error_path=not-run; options=--fd-order 4 --no-ko; date=09-05-2026`

### Dendro's own vocabulary

The emitted names follow Dendro-GR's BSSN solver rather than NRPy's: the solver
directory is `BSSN_GR`, CMake variables carry the `BSSN_` prefix, the object
library is `bssn_common`, the executable is `bssnSolver`, the namespace is
`bssn`, and the context source is `bssnCtx.cpp`. The emitted kernel names are
NRPy's, not the host's: `bssn_rhs_eval*`, `bssn_constraints_eval*` and
`bssn_enforce_detgbar_equals_detghat_trAzero*` carry the stem the way ETLegacy
carries its thorn name, and take the operation names BHaH and ETLegacy already
use. Dendro-GR's own `bssn_constraints.cpp` is the det/trace *enforcement*, and
its H/M diagnostics live in `physcon.cpp`, so upstream vocabulary would have
pointed the constraint-diagnostics kernel at the wrong operation.

### What the port found

Adding the second formulation required generic-layer work: the
formulation-agnostic lowering was extracted into `block_kernel_helpers` and
`tensor_family_of` moved into `gridfunction_name_decorations`. No formulation-specific content entered
the generic layer and no existing emitter changed behaviour.

The port was the cheap test of whether the abstraction exists, and it found
defects that a single-formulation tree could not expose:

- The "shared" ADM-to-evolved conversion asserted that exactly one evolved
  field was left undefined — true only for fCCZ4's Z4 scalar. The rule is now
  that *at most* one constraint scalar is zeroed.
- The diagnostics accessed-set intersected raw free symbols with the registry,
  so a field that a kernel only *differentiates* was never bound. The BSSN
  momentum constraint takes derivatives of `lambdaU`, and the emitted kernel
  read a pointer nothing declared. `block_kernel_helpers.accessed_gridfunctions` now
  composes the canonical NRPy derivative extraction to resolve each derivative
  symbol back to the field it differentiates.
- The diagnostic role pointer (`aux_pointer` then, `diag_pointer` now) had
  been deleted as dead code in an earlier review round. It was dead only
  because the tree had one formulation.
- The shared initial-data and enforcement builders hardcoded `fccz4_` into the
  CFunction names they registered, so the first BSSN solver shipped eight
  `fccz4_`-named sources. The stem is now threaded from the caller.
- The perturbation added an identical profile to every field, and 22 of 24
  fields have an asymptotic value of zero, so the probe state carried two
  distinct component values and the `FLATADAPTER` gate could not see a
  component bound to the wrong flat-layout slab. Each component is now scaled
  by one plus its registry position.

## Sources

- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN/rhs_eval.py) - `BSSN_rhs_expressions`, `build_rhs_eval`, `register_CFunctions_rhs_eval`
- [constraints_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN/constraints_eval.py) - `build_constraints_eval`, `register_CFunctions_constraints_eval`
- [block_kernel_helpers.py](../../../nrpy/infrastructures/Dendro/block_kernel_helpers.py) - `accessed_gridfunctions`, `padding_from_operators`, `emitted_operators`
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
