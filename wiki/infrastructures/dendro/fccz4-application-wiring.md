# fCCZ4 Application Wiring

> Explain the Dendro fCCZ4 builders: the direct finite-difference right-hand side, the algebraic projection, initial-data conversion, and constraint diagnostics. · Status: provisional · Last reconciled: 09-04-2026
> Up: [Dendro](index.md)

## Summary

Four builders lower the shared fCCZ4 expression factory into Dendro kernels.
The scientific formulation comes from the shared factory and from established
NRPy projectors and conversions, and every field name is read back from the
registry rather than written down. What the builders do author is bounded and
deliberate: the analytic test perturbation profile, the two algebraic residuals
the projection reports, and the rescaled connection relation — all assembled
from registered quantities, and all registered as CFunctions rather than written
into a fixed template. Each builder
returns CFunction bodies and metadata; the caller registers them, freezes the
environment, and exports. Keeping registration out of the builders is what lets
one profile assemble a different subset without the builders knowing about each
other.

## Detail

### Right-hand side

The RHS builder maps each RHS symbol to its registered EVOL gridfunction name
algorithmically and asserts the bijection against the registry, derives the
output lvalues and the `in_` and `rhs_` pointer bindings from the gridfunction
registry, and runs `c_codegen` with direct finite differences inside a scoped
access-capture context using the shared factory's upwind control vector and the
`DendroScalar` alias. Point and block loops are emitted through the Dendro loop
helpers. It provides three registered CFunction bodies — per-block, all-block,
and a local-time-stepping flat-block adapter that reuses the same numerical
body — and exposes the operator and access manifests together with the
access-derived padding.

### Algebraic projection

The projection restores the two algebraic constraints of the conformal
decomposition at every point of a block: the conformal metric determinant ratio
returns to one, and the trace of the conformal traceless extrinsic curvature
returns to zero. The projected values come from the established NRPy projector
`BSSN_algebraic_constraints`, so the module contributes no new formulation
content; it lowers those expressions into a Dendro point loop and adds the
structured status record the generated host lifecycle consumes. The kernel never
calls `exit()`.

Because the projection is pointwise, it needs no padding, and computing every
projected value into a local before storing is what makes an in-place call
safe: the input and output pointers may alias the same block arrays. That
argument is specific to a pointwise kernel and does not extend to the stencil
kernels on this page.

Every written field name is read back from the registered BSSN quantities, so
no field name is hardcoded. Only `hDD` and `aDD` are written: the projection is
purely algebraic in those two families, and `lambdaU` and `Theta_fCCZ4` are
never touched by it.

### Initial data

The Minkowski fill writes every EVOL field to its asymptotic value through
exact-name `out_` bindings derived from the frozen EVOL registry. The role is
`out_` rather than `rhs_` because this writer produces state, not a right-hand
side; reusing `rhs_` would make the generated signature claim it fills the
right-hand-side vector, which is how a caller silently zeroes the state. No field name,
count, or asymptotic value is hardcoded: all three come from the registered
records, which makes this the first end-to-end proof that state allocation,
field order, loops, and CFunctions agree. The module's own docstring still says
`rhs_<name>` here; that wording is stale, and the pointer-binding helper is the
deciding source.

Claim evidence:
- Claim: the generated Minkowski initial-data writer binds its outputs through `out_<name>` pointers, not `rhs_<name>`, because it produces state rather than a right-hand side.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/general_relativity/initial_data.py`, `_block_pointer_bindings`
- Corroboration: `nrpy/infrastructures/Dendro/naming.py`, `out_pointer` docstring recording why `rhs_` was rejected
- Validation: `inspected=pass; generated=pass; built=pass; run=not-run; result_checked=not-run`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3, GCC 13.3.0, CMake 3.28.3; backend=Dendro; precision=double; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=--fd-order 4 --no-ko; date=09-04-2026`

A separate builder emits the smooth analytic perturbation that makes the
stencils observable. Its profile is authored in NRPy and lowered by `c_codegen`,
so it is formulation content and therefore belongs in a registered CFunction
rather than in a fixed template. It is the only physics *profile* a builder
authors, as distinct from the two algebraic residuals and the connection
relation named above. The generated Minkowski lifecycle uses it for
the perturbed-RHS and convergence-order gates, which a flat state cannot
exercise.

The smooth ADM-to-fCCZ4 conversion reuses the established `ADM_to_BSSN` map, so
the registered `EvolvedConformalFactor_cf` choice governs the conformal
convention in exactly one place. A separate pass initializes the connection as
`lambdaU^i = DeltaGamma^i / ReU^i`, which is the statement that the connection
constraint holds at the initial slice. Splitting the connection pass out is
deliberate: it depends on derivatives of quantities the conversion pass has just
written, so it cannot share the conversion's point loop.

### Constraint diagnostics

The first diagnostic set is the Hamiltonian constraint of the fCCZ4 system and
the spatial Z4 connection constraint. Both come from the shared expression
factory with diagnostics enabled, so a Dendro kernel and any other
infrastructure lower the same expressions. The diagnostic gridfunctions are
registered in the DIAG group: they are recomputed from the evolved state and are
never authoritative checkpoint state. The kernel uses the same
finite-difference order and the same memory-access mechanism as the RHS, so no
separate diagnostic profile exists.

### Shared factory boundary

The formulation itself — evolution equations, gauge, algebraic constraints,
diagnostics, and the Kreiss-Oliger terms — belongs to the equations layer and is
validated there against trusted expression dictionaries. This page covers only
the lowering. For the equations themselves see
[Fully Covariant Conformal Z4](../../equations/general-relativity/fccz4.md).

## Sources

- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - `build_fccz4_rhs`, `FCCZ4RHSBuild`
- [projection.py](../../../nrpy/infrastructures/Dendro/general_relativity/projection.py) - `build_projection`, `register_projection_CFunctions`
- [initial_data.py](../../../nrpy/infrastructures/Dendro/general_relativity/initial_data.py) - `build_minkowski_initial_data`, `build_adm_to_evolved`, `build_lambda_initialization`
- [diagnostics.py](../../../nrpy/infrastructures/Dendro/general_relativity/diagnostics.py) - `build_diagnostics`, `tensor_family_of`, `register_diagnostics_CFunctions`
- [fCCZ4_system.py](../../../nrpy/equations/general_relativity/fCCZ4_system.py) - `build_fccz4_expression_bundle`
- [kreiss_oliger_terms.py](../../../nrpy/equations/general_relativity/kreiss_oliger_terms.py) - Kreiss-Oliger dissipation terms
- [BSSN_algebraic_constraints.py](../../../nrpy/equations/general_relativity/BSSN_algebraic_constraints.py) - `BSSN_algebraic_constraints`
- [ADM_to_BSSN.py](../../../nrpy/equations/general_relativity/ADM_to_BSSN.py) - `ADM_to_BSSN`

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Fully Covariant Conformal Z4](../../equations/general-relativity/fccz4.md)
- Implements: [Gridfunctions, Naming, Access Capture, And Loops](gridfunctions-naming-access-capture-and-loops.md)
- Contrasts with: [GR Application Wiring](../bhah/gr-application-wiring.md)
- Validated by: [Validation, Host Mock, And Deferral Gates](validation-host-mock-and-deferral-gates.md)
- See also: [C Codegen](../../core/c-codegen.md)
