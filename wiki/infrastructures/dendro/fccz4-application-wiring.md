# fCCZ4 Application Wiring

> Describe fCCZ4 equations, initial data, conversions, and runtime services in `NRPy_fCCZ4_GR`. · Status: provisional
> Up: [Dendro](index.md)

## Summary

`nrpy.examples.dendro_fccz4` emits a complete fCCZ4 application with W or chi
selected at code-generation time (W by default). Its
24-component BSSN-compatible state is followed by `Theta_fCCZ4`. It uses SSL,
constant default `eta=1`, centered derivatives, its fCCZ4 CAHD contribution,
and separate conformal Ricci and RHS kernels.

## Detail

The fCCZ4 CAHD contribution is

```text
W_rhs   += 2 C_CAHD W H dt,    for W evolution;
chi_rhs += 4 C_CAHD chi H dt,  for chi evolution.
Default C_CAHD = 0.15.
```

`C_CAHD` is a runtime parameter. The chi term follows from `chi=W^2` and
`chi_rhs=2 W W_rhs`. Here `H` is the BSSN-shaped Hamiltonian expression used
by the RHS adjustment and now also reported separately; it is not the distinct
`H_Z4` constraint. This equation is not the BSSN CAHD equation. fCCZ4 retains
its formulation's Lambda RHS rather than adding Brown's BSSN adjustment. Order-specific
constraint kernels are named `fCCZ4_constraints_order_N`.

Claim evidence:
- Claim: fCCZ4's W and chi CAHD terms use the BSSN-shaped Hamiltonian expression with a runtime coefficient, not `H_Z4`; changing the generated conformal factor changes the coefficient from 2 to 4.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py`, `register_CFunction_rhs_eval`.
- Corroboration: `nrpy/infrastructures/Dendro/general_relativity/fCCZ4_constraints.py`, `register_CFunction_fCCZ4_constraints`, emits both `H` and formulation-specific `H_Z4` diagnostics.

The initial octree uses the same analytic Dendro-GR puncture seed as BSSN.
Evolved initial data uses the same TwoPunctures ADM data, full
`psi=psi_background+u`, `alpha=W=psi^(-2)`, ADM conversion, halo exchange,
physical-boundary fill, separate `initial_data_lambdaU`, and algebraic
projection sequence as BSSN. `Theta_fCCZ4` receives its formulation-defined
initial value. W and chi builds have distinct checkpoint formulation IDs and
reject cross-formulation restores.

After initialization, every RK stage, and AMR transfer, `alpha` is floored at
`CHI_FLOOR` and the evolved conformal factor at `CHI_FLOOR` for chi or
`sqrt(CHI_FLOOR)` for W before algebraic projection.

Solver context evaluates Ricci and RHS consecutively for each block after one
halo exchange and physical exterior-ghost fill. Six `RbarDD` values remain scratch storage,
not evolved state. Apparent-horizon searches interpolate the evolved
fCCZ4 fields and convert each search point to ADM variables. ADM surface
quantities use grid fields from `BSSN_to_ADM`; waveform extraction evaluates
Psi4 directly from the evolved fields.

The diagnostic kernel retains `H_Z4` and the three `Z4constraintU` components
as its first four fields. It then reports `H`, physical lower-index momentum
components `MU0..MU2`, `M_CONSTRAINT=sqrt(gamma_ij M^i M^j)`, and
`LAMBDA_CONSTRAINT=sqrt(gammabar_ij Z4constraintU^i Z4constraintU^j)`.
The connection residual is reported, not enforced. Both excised physical-volume
and unique-node RMS files include all ten fields and are emitted after remeshing.

The shared runtime parses puncture-centered AMR controls, tracks centers with
`vetU`, retains puncture history in checkpoints, and computes diagnostics
after any scheduled remesh. This shared path does not make the fCCZ4
constraint or CAHD equations identical to BSSN's.

## Sources

- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - fCCZ4 generation profile.
- [param_toml.py](../../../nrpy/infrastructures/Dendro/param_toml.py) - sample parameter file and TwoPunctures lapse default.
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - fCCZ4 RHS and CAHD registration.
- [fCCZ4_constraints.py](../../../nrpy/infrastructures/Dendro/general_relativity/fCCZ4_constraints.py) - fCCZ4 constraint registration.
- [Ricci_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/Ricci_eval.py) - separate conformal Ricci kernel.
- [ADM_to_BSSN.py](../../../nrpy/infrastructures/Dendro/general_relativity/ADM_to_BSSN.py) - ADM conversion.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - traversal and runtime scheduling.
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - analytic initial-grid seed and startup.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [fCCZ4 System](../../equations/general-relativity/fccz4.md)
- Contrasts with: [BSSN Application Wiring](bssn-application-wiring.md)
- See also: [Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md)
