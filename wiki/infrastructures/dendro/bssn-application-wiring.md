# BSSN Application Wiring

> Explain the Dendro BSSN builders, what they reuse from the fCCZ4 profile, and what adding a second formulation proved about the generic layer. · Status: provisional
> Up: [Dendro](index.md)

## Summary

BSSN is the second formulation lowered through the Dendro infrastructure, and
it lives beside the fCCZ4 branch in `general_relativity/`. The shared
right-hand-side and constraint-diagnostic modules each use one builder with
formulation-specific expression assembly. Initial data and the
det(gammabar)/tr(Abar) enforcement are shared with fCCZ4, while the GR-owned
assembly modules supply both formulations' scientific content to the generic
Dendro emitters.

The layout follows BHaH. NRPy's established two-formulation infrastructure
emits BSSN and fCCZ4 from a single `general_relativity/rhs_eval.py` on an
`enable_fCCZ4` boolean, and Dendro follows it as far as the module layout:
one module per generated file, taking the same boolean.

The per-formulation builder and registrar duplication is consolidated. Each
module now has one public builder and one public registrar, with the
formulation branch confined to expression assembly and its field-count check.

## Detail

The example now assembles GR-owned context policy, lifecycle CTest statements,
projection status, and scientific self-tests explicitly into generic Dendro
emitters. BSSN/W is checked at regular finite-difference orders 4, 6, and 8 on
a nonflat fixed block, component by component, for the RHS block kernel, flat
adapter, and constraint diagnostics.
The expected values use the same exactly emitted binary64 samples as the kernels
and an independent 80/100-digit stencil-and-CSE evaluation. Order 6 is the
configured default, and KO-off is the BSSN default. [Finite-Difference Profiles
And Dendro Conformance](finite-difference-profiles-and-dendro-conformance.md)
defines the supported order pairs, stencil reach, Dendro-GR comparison, and
qualification limits.

### Shared assembly and registration

`build_rhs_eval` and `build_constraints_eval` share lowering, wrapper emission,
and metadata handling across BSSN and fCCZ4. Their public registrars share the
registration path. The registered interfaces and generated kernels are
unchanged. This closes [CONTR-0011](../../contradictions.md#contr-0011).

Claim evidence:
- Claim: each Dendro RHS and constraint module has one builder and one registrar shared across formulations; consolidation preserves the registered interfaces and generated kernels for the shipped BSSN and fCCZ4 profiles.
- Role: descriptive behavior
- Deciding authority: [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py), `build_rhs_eval` and `register_CFunctions_rhs_eval`; [constraints_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/constraints_eval.py), `build_constraints_eval` and `register_CFunctions_constraints_eval`
- Corroboration: complete generated BSSN and fCCZ4 C++ project tests, including independent nonflat RHS and diagnostic references, exercise both shared registration paths

### What the builders do

`general_relativity/rhs_eval.py`'s BSSN builder assembles the evolution system
the way ETLegacy and BHaH do: the non-gauge equations from the cached
`BSSN_RHSs` object, the lapse and shift from `BSSN_gauge_RHSs` added to a *copy*
of its dictionary, Kreiss-Oliger terms through the shared
`add_KreissOliger_dissipation_terms` helper with `include_Theta_fCCZ4=False`,
and centered advection after normalizing the equation factory's directional
derivative symbols. It then emits all of it through `block_kernel_helpers`,
asserts the 24-field bijection against the registry, and records the exact
padding implied by the selected regular and KO operators. NRPy-generated
kernels do not implement Dendro-GR's `bflag`-dependent physical-boundary
derivative closures; the conformance page owns the exact comparison.

Claim evidence:
- Claim: NRPy-generated Dendro BSSN kernels use padded centered stencils and do not implement Dendro-GR's `bflag`-dependent physical-boundary derivative closures.
- Role: descriptive behavior
- Deciding authority: [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py), `build_rhs_eval`; [simple_loop.py](../../../nrpy/infrastructures/Dendro/simple_loop.py), `loop`
- Corroboration: [Dendro-GR derivs.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/derivs.cpp), derivative routines that consume `bflag`

`general_relativity/constraints_eval.py`'s BSSN builder emits the Hamiltonian
constraint and the three momentum constraint components from the established
`BSSN_constraints` factory.

### Production Dendro-GR application adapter

`nrpy.examples.dendro_bssn --dendro-gr-host` selects `chi` as the evolved
conformal factor and emits an opt-in `nrpy_bssnSolver` target for the public
[Dendro-GR repository](https://github.com/paralab/Dendro-GR). The target
compiles Dendro-GR's BSSN application and context, but renames that context's
`bssnRHS` reference to `nrpy_bssnRHS`. The adapter therefore replaces only the
right-hand-side evaluation. Dendro-GR still owns parameter input,
TwoPunctures or approximate initial data, octree construction, mesh changes,
Runge--Kutta evolution, physical-boundary treatment, algebraic projection,
diagnostics, checkpoints, and output.

Claim evidence:
- Claim: `--dendro-gr-host` emits an opt-in `nrpy_bssnSolver` that retains the Dendro-GR BSSN application and replaces its `bssnRHS` call with the generated block kernel.
- Role: public generated application interface
- Deciding authority: [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py), `parse_args` and `main`; [bssn_host_adapter.py](../../../nrpy/infrastructures/Dendro/general_relativity/bssn_host_adapter.py), `output_bssn_host_files`
- Corroboration: [host test README](../../../nrpy/infrastructures/Dendro/tests_infra/README.md#full-dendro-gr-bssn-application), generation, build, and short-run procedure against the public host

Dendro-GR stores the full conformal metric in its legacy 24-component order.
The generated BSSN kernel stores `hDD = gammabarDD - deltaDD` in NRPy's
registry order. The adapter remaps component pointers and subtracts one from
only the three diagonal metric inputs. Their right-hand sides need no value
conversion because the time derivative of the Cartesian reference metric is
zero. It maps `ETA_CONST` to generated `eta` and maps `KO_DISS_SIGMA` to both
generated KO strengths.

The adapter checks the host element order, block padding, scalar type, field
count, block geometry, and generated parameters before evaluation. It applies
the generated centered block kernel, then Dendro-GR's radiative physical
boundary routine. A KO-enabled profile adds Dendro-GR's boundary KO values
only at physical-boundary points because generated KO already covers the
remaining block interior. Conformal-factor-scaled KO and CUDA host execution
are rejected.

This path changes the equations evaluated by the Dendro-GR application. It
uses NRPy's centered-advection BSSN equations and lower-base-order KO profile.
Dendro-GR compile definitions for legacy SSL or CAHD right-hand-side terms do
not add those terms to the generated equations. The host's initial-lapse
choice and all retained mesh and output operations remain active. Exact build,
input, and run commands live in the [host test README](../../../nrpy/infrastructures/Dendro/tests_infra/README.md#full-dendro-gr-bssn-application).

### The DIAG-before-factory ordering, and why

`BSSN_constraints` registers `H`, `M` and `LAMBDA_CONSTRAINT` into the **AUX**
group when it is constructed, each guarded by an existence check. Two
Generation-time NRPy parameters gate those registrations. `M` and
`LAMBDA_CONSTRAINT` are gated
by `register_M_and_LAMBDA_CONSTRAINT_gridfunctions`, which defaults to `True`,
so a caller that wants neither must ask; `H` is unconditional. `MU` is
different again: the factory registers it only under
`register_MU_gridfunctions`, which defaults to `False` and which no Dendro
module sets, so nothing competes for those names. Both gates are keys of the
`BSSNconstraints_dict` memo, so flipping either forces a rebuild rather than
returning an object constructed under the previous setting.

Dendro's diagnostics contract is the
**DIAG** group, which is the settled infrastructure convention: BHaH registers
DIAG gridfunctions across its wave-equation, elliptic and GR diagnostics,
while AUX appears in a separate helper role. Every registration in the factory is guarded by
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
- Claim: `BSSN_constraints` registers `H`, `M` and `LAMBDA_CONSTRAINT` into the AUX group, each guarded by an existence check, so a caller that registers `H` first determines its group; `M` and `LAMBDA_CONSTRAINT` are additionally gated by the generation-time NRPy parameter `register_M_and_LAMBDA_CONSTRAINT_gridfunctions`, which defaults to `True`, and `MU` by the generation-time NRPy parameter `register_MU_gridfunctions`, which defaults to `False` and which no Dendro module sets. Both gates are keys of the `BSSNconstraints_dict` construction-parameter memo, so flipping either forces a rebuild. The Dendro BSSN builder pre-registers the names it writes as DIAG and suppresses the two AUX names it does not write by setting `register_M_and_LAMBDA_CONSTRAINT_gridfunctions` to `False` across the factory construction and restoring it afterwards, rather than deleting anything from `gri.glb_gridfcs_dict`, leaving `NUM_AUX_GFS = 0` in the emitted state header.
- Role: descriptive behavior
- Deciding authority: [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), the `group="AUX"` registrations and their `not in gri.glb_gridfcs_dict` guards
- Corroboration: [constraints_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/constraints_eval.py), the `register_M_and_LAMBDA_CONSTRAINT_gridfunctions` save/set/restore around the `BSSN_constraints` construction in the BSSN builder

### NRPy module identity and BSSN names

NRPy authors the generated module, so its directory and CMake project are
`nrpy_bssn`, and its C++ namespace is `nrpy::bssn`. Names inside the module
retain the conventional BSSN stem: CMake variables carry the `BSSN_` prefix,
the production library is `nrpy_bssn_dendro`, the qualification executable is
`nrpy_bssn_dendro_qualify`, and the context source is `bssnCtx.cpp`. The kernel
names `bssn_rhs_eval*`, `bssn_constraints_eval*` and
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
changed behavior.

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
  CFunction names they registered, so the BSSN solver shipped
  `fccz4_`-named sources. The stem is now threaded from the caller.
- The perturbation added an identical profile to every field, and most fields
  have an asymptotic value of zero, so the probe state carried too few
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
- [bssn_host_adapter.py](../../../nrpy/infrastructures/Dendro/general_relativity/bssn_host_adapter.py) - Dendro-GR application target, component mapping, physical boundaries, and host parameter mapping
- [tests_infra/README.md](../../../nrpy/infrastructures/Dendro/tests_infra/README.md) - Dendro-GR application generation, build, and run procedure

## See Also

- Parent: [Dendro](index.md)
- Depends on: [BSSN Family](../../equations/general-relativity/bssn-family.md)
- Depends on: [Finite-Difference Profiles And Dendro Conformance](finite-difference-profiles-and-dendro-conformance.md)
- Implements: [Gridfunctions, Naming, And Loops](gridfunctions-naming-and-loops.md)
- Contrasts with: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
- Validated by: [Validation, Standalone Host, And Deferred Tests](validation-standalone-host-and-deferral-gates.md)
