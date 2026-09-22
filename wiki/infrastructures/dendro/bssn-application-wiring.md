# BSSN Application Wiring

> Describe BSSN equations, initial data, conversions, and runtime services in `NRPy_BSSN_GR`. · Status: provisional
> Up: [Dendro](index.md)

## Summary

`nrpy.examples.dendro_bssn` emits a complete W-form BSSN application. It uses
centered derivatives, SSL, constant default `eta=1`, Dendro-style CAHD,
Brown's covariant Lambda adjustment, and separate conformal Ricci and RHS
kernels. Runtime parameter input can override the generated eta default.

## Detail

The BSSN CAHD contribution is

```text
W_rhs += (C_CAHD/2) W H dx^2 / (dt (1 + 10 dx^2)),  C_CAHD = 0.06.
```

`dx` and `dt` are the current block spacing and timestep. Brown's adjustment
appears exactly once in the covariant `lambdaU` RHS. The generated constraint
kernels are named `BSSN_constraints_order_N`; there is no generic constraint
catch-all source.

TwoPunctures produces ADM data and full `psi=psi_background+u`. Initialization
sets `alpha=W=psi^(-2)` and follows this sequence:

```text
TwoPunctures -> ADM_to_BSSN -> zip/exchange/unzip -> physical boundary
             -> initial_data_lambdaU -> algebraic projection
```

Apparent-horizon searches interpolate selected evolved W-BSSN fields and then
convert each search point to ADM variables. `BSSN_to_ADM` reconstructs grid
fields used by ADM surface quantities; waveform extraction evaluates Psi4
directly from W-BSSN fields. Generated runtime services also provide
physical boundaries, volume-normalized constraint diagnostics, checkpoint and
restart, and scheduled output.

After one evolved-state halo exchange and physical-boundary fill, solver
context traverses blocks once. For each block it calls `Ricci_eval_order_N`
immediately followed by `rhs_eval_order_N`. No synchronization or second block
traversal occurs between them. Algebraic projection runs after initialization,
each RK stage before the next exchange, and AMR transfer.
At the same three points, `alpha` is floored at `CHI_FLOOR` and W at
`sqrt(CHI_FLOOR)` before projection. This preserves Dendro-BSSN's chi-floor
scale when the evolved conformal factor is W.

## Sources

- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - BSSN generation profile.
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - BSSN RHS, SSL, CAHD, and Brown adjustment registration.
- [BSSN_constraints.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py) - BSSN constraint kernel registration.
- [ADM_to_BSSN.py](../../../nrpy/infrastructures/Dendro/general_relativity/ADM_to_BSSN.py) - ADM conversion.
- [initial_data_lambdaU.py](../../../nrpy/infrastructures/Dendro/general_relativity/initial_data_lambdaU.py) - connection initialization.
- [twopunctures.py](../../../nrpy/infrastructures/Dendro/general_relativity/twopunctures.py) - TwoPunctures data and interpolation.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - traversal and runtime scheduling.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [BSSN Family](../../equations/general-relativity/bssn-family.md)
- Contrasts with: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
