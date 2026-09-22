# fCCZ4 Application Wiring

> Describe fCCZ4 equations, initial data, conversions, and runtime services in `NRPy_fCCZ4_GR`. · Status: provisional
> Up: [Dendro](index.md)

## Summary

`nrpy.examples.dendro_fccz4` emits a complete W-form fCCZ4 application. Its
24-component BSSN-compatible state is followed by `Theta_fCCZ4`. It uses SSL,
constant default `eta=1`, centered derivatives, its fCCZ4 CAHD contribution,
and separate conformal Ricci and RHS kernels.

## Detail

The fCCZ4 CAHD contribution is

```text
W_rhs += 2 C_CAHD W H dt,  C_CAHD = 0.15.
```

This equation is not the BSSN CAHD equation. fCCZ4 retains its formulation's
Lambda RHS rather than adding Brown's BSSN adjustment. Order-specific
constraint kernels are named `fCCZ4_constraints_order_N`.

Initial data uses the same TwoPunctures ADM data, full
`psi=psi_background+u`, `alpha=W=psi^(-2)`, ADM conversion, halo exchange,
physical-boundary fill, separate `initial_data_lambdaU`, and algebraic
projection sequence as BSSN. `Theta_fCCZ4` receives its formulation-defined
initial value.

After initialization, every RK stage, and AMR transfer, `alpha` is floored at
`CHI_FLOOR` and W at `sqrt(CHI_FLOOR)` before algebraic projection.

Solver context evaluates Ricci and RHS consecutively for each block after one
halo exchange and boundary fill. Six `RbarDD` values remain scratch storage,
not evolved state. Apparent-horizon searches interpolate selected evolved
W-BSSN fields and convert each search point to ADM variables. ADM surface
quantities use grid fields from `BSSN_to_ADM`; waveform extraction evaluates
Psi4 directly from the evolved fields.

## Sources

- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - fCCZ4 generation profile.
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - fCCZ4 RHS and CAHD registration.
- [fCCZ4_constraints.py](../../../nrpy/infrastructures/Dendro/general_relativity/fCCZ4_constraints.py) - fCCZ4 constraint registration.
- [Ricci_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/Ricci_eval.py) - separate conformal Ricci kernel.
- [ADM_to_BSSN.py](../../../nrpy/infrastructures/Dendro/general_relativity/ADM_to_BSSN.py) - ADM conversion.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - traversal and runtime scheduling.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [fCCZ4 System](../../equations/general-relativity/fccz4.md)
- Contrasts with: [BSSN Application Wiring](bssn-application-wiring.md)
