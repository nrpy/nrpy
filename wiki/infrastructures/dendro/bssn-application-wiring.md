# BSSN Application Wiring

> Explain the Dendro BSSN builders, what they reuse from the fCCZ4 profile, and what adding a second formulation proved about the generic layer. · Status: contested · Last reconciled: 09-07-2026
> Up: [Dendro](index.md)

## Summary

BSSN is the second formulation lowered through the Dendro infrastructure, and
it lives beside the fCCZ4 branch in `general_relativity/`. No module is new: the
shared right-hand-side and constraint-diagnostic modules each gained a BSSN
builder beside the fCCZ4 one, and the initial data, the det(gammabar)/tr(Abar)
enforcement and the generic layer are shared with fCCZ4.

The layout follows BHaH. NRPy's established two-formulation infrastructure
emits BSSN and fCCZ4 from a single `general_relativity/rhs_eval.py` on an
`enable_fCCZ4` boolean, and Dendro follows it as far as the module layout:
one module per artifact, taking the same boolean. An earlier draft of this
branch put the BSSN builders in their own `general_relativity/BSSN/`
subpackage, on the argument that a single module would branch on formulation in
every place the two differ: evolved-field count, constraint set, and which
shared equations module supplies the expressions.

Where Dendro departs from BHaH is inside the module. BHaH's
`general_relativity/rhs_eval.py` is a single public
`register_CFunction_rhs_eval` that branches inline on `enable_fCCZ4` in the few
places the formulations differ. Dendro instead holds a private builder and a
private registrar per formulation, and the public registrar dispatches once on
the boolean. That arrangement is Dendro's own, and no host requirement drove
it, so it is a divergence [New Infrastructure
Conformance](../new-infrastructure-conformance.md) does not permit. It leaves
the constraint builders sharing a long verbatim tail and the right-hand-side
builders substantially duplicated. Consolidating them is open.

Claim status: contested; contradiction: CONTR-0011.
See [CONTR-0011](../../contradictions.md#contr-0011) for the deciding authority
and the inspection that would resolve it.

## Detail

### What the builders do

`general_relativity/rhs_eval.py`'s BSSN builder assembles the evolution system
the way ETLegacy and BHaH do: the non-gauge equations from the cached
`BSSN_RHSs` object, the lapse and shift from `BSSN_gauge_RHSs` added to a *copy*
of its dictionary, Kreiss-Oliger terms through the shared
`add_KreissOliger_dissipation_terms` helper with `include_Theta_fCCZ4=False`,
and the upwind control vector as the rescaled shift
`betaU[i] = vetU[i] * ReU[i]`. It then emits all of it through
`block_kernel_helpers`, asserts the 24-field bijection against the registry, and
records the padding and the upwind control set.

`general_relativity/constraints_eval.py`'s BSSN builder emits the Hamiltonian
constraint and the three momentum constraint components from the established
`BSSN_constraints` factory.

### The DIAG-before-factory ordering, and why

`BSSN_constraints` registers `H`, `M` and `LAMBDA_CONSTRAINT` into the **AUX**
group when it is constructed, each guarded by an existence check. Two
CodeParameters gate those registrations. `M` and `LAMBDA_CONSTRAINT` are gated
by `register_M_and_LAMBDA_CONSTRAINT_gridfunctions`, which defaults to `True`,
so a caller that wants neither must ask; `H` is unconditional. `MU` is
different again: the factory registers it only under
`register_MU_gridfunctions`, which defaults to `False` and which no Dendro
module sets, so nothing competes for those names. Both gates are keys of the
`BSSNconstraints_dict` memo, so flipping either forces a rebuild rather than
returning an object constructed under the previous setting.

Dendro's diagnostics contract is the
**DIAG** group, which is the settled infrastructure convention: BHaH registers
31 DIAG gridfunctions across its wave-equation, elliptic and GR diagnostics,
while AUX appears once. Every registration in the factory is guarded by
`if <name> not in gri.glb_gridfcs_dict`, so the Dendro builder registers the
names it writes as DIAG *before* constructing the factory. That is
load-bearing for `H`; `MU0`-`MU2` are DIAG simply because this builder is their
sole registrant. The factory would otherwise also register `M` and
`LAMBDA_CONSTRAINT`, which this kernel does not compute, so the builder
suppresses those two at the source: it reads
`register_M_and_LAMBDA_CONSTRAINT_gridfunctions`, sets it to `False` for the
duration of the factory construction, and restores the previous value in a
`finally` block so no other caller inherits this builder's choice. Nothing is
deleted from `glb_gridfcs_dict` afterwards -- suppressing the registration is
what makes the deletion unnecessary, and a deletion pass would have to tell the
two names apart from the evolved state the same construction pulls in.
Without the suppression the generated state header would advertise two
variables no kernel writes and no vector backs. No change to the shared
equations module.

Claim evidence:
- Claim: `BSSN_constraints` registers `H`, `M` and `LAMBDA_CONSTRAINT` into the AUX group, each guarded by an existence check, so a caller that registers `H` first determines its group; `M` and `LAMBDA_CONSTRAINT` are additionally gated by the `register_M_and_LAMBDA_CONSTRAINT_gridfunctions` CodeParameter, which defaults to `True`, and `MU` by `register_MU_gridfunctions`, which defaults to `False` and which no Dendro module sets. Both gates are keys of the `BSSNconstraints_dict` construction-parameter memo, so flipping either forces a rebuild. The Dendro BSSN builder pre-registers the names it writes as DIAG and suppresses the two AUX names it does not write by setting `register_M_and_LAMBDA_CONSTRAINT_gridfunctions` to `False` across the factory construction and restoring it afterwards, rather than deleting anything from `gri.glb_gridfcs_dict`, leaving `NUM_AUX_GFS = 0` in the emitted state header.
- Role: descriptive behavior
- Deciding authority: [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), the `group="AUX"` registrations and their `not in gri.glb_gridfcs_dict` guards
- Corroboration: [constraints_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/constraints_eval.py), the `register_M_and_LAMBDA_CONSTRAINT_gridfunctions` save/set/restore around the `BSSN_constraints` construction in the BSSN builder
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3, GCC 13.3.0, CMake 3.28.3, OpenMPI 4.1.6; backend=Dendro; precision=double; GPU=not-applicable; restart=not-applicable; distributed=1 and 2 MPI ranks; error_path=not-run; options=--fd-order 4 --no-ko; date=09-07-2026`

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
formulation-agnostic kernel emission was extracted into `block_kernel_helpers`
and `tensor_family_of` moved into `gridfunction_name_decorations`. No
formulation-specific content entered the generic layer and no existing emitter
changed behaviour.

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

- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - `BSSN_rhs_expressions`, `build_rhs_eval`, `register_CFunctions_rhs_eval`
- [constraints_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/constraints_eval.py) - `build_constraints_eval`, `register_CFunctions_constraints_eval`
- [block_kernel_helpers.py](../../../nrpy/infrastructures/Dendro/block_kernel_helpers.py) - `accessed_gridfunctions`, `padding_from_derivative_operators`, `emitted_derivative_operators`
- [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py) - `BSSNRHSs`
- [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py) - `BSSN_gauge_RHSs`
- [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py) - `BSSNconstraints`
- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - generation entry point and command-line profile

## See Also

- Parent: [Dendro](index.md)
- Depends on: [BSSN Family](../../equations/general-relativity/bssn-family.md)
- Implements: [Gridfunctions, Naming, And Loops](gridfunctions-naming-and-loops.md)
- Contrasts with: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
- Validated by: [Validation, Standalone Host, And Deferral Gates](validation-standalone-host-and-deferral-gates.md)
