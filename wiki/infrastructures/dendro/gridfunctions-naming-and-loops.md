# Gridfunctions, Naming, And Loops

> Define Dendro storage names, fixed component order, and block-local loop generation. · Status: provisional
> Up: [Dendro](index.md)

## Summary

`DendroGridFunction` supplies the Dendro C++ storage expression directly. The
Dendro connector maps NRPy's canonical conformal-factor gridfunction `cf` to
`cf_W_or_chi`; equation modules retain `cf`. Numerical registrars bind pointers
from fixed component lists and emit direct `ot::Block` loops.

## Detail

BSSN evolved storage has this fixed order:

```text
alpha, cf_W_or_chi, trK,
lambdaU0, lambdaU1, lambdaU2,
vetU0, vetU1, vetU2,
betU0, betU1, betU2,
hDD00, hDD01, hDD02, hDD11, hDD12, hDD22,
aDD00, aDD01, aDD02, aDD11, aDD12, aDD22
```

fCCZ4 appends `Theta_fCCZ4`. Conformal Ricci storage is a separate six-component
scratch list: `RbarDD00`, `RbarDD01`, `RbarDD02`, `RbarDD11`, `RbarDD12`, and
`RbarDD22`. Enum values, names, pointer indices, transfer, output, and checkpoint
metadata use these same lists. Alphabetical registry sorting never defines
runtime storage.

Each registrar emits its pointer declarations once in registry order. It does
not infer inputs by scanning expressions or rewrite identifiers with strings or
regular expressions. `simple_loop.py` emits x-fastest padded-block loops without
a nested OpenMP region. Numerical kernels receive `ot::Block` directly and
derive offsets, dimensions, spacing, padding, and interior bounds from it.

## Sources

- [grid.py](../../../nrpy/grid.py) - `DendroGridFunction` and `dendro_name`.
- [state_h.py](../../../nrpy/infrastructures/Dendro/state_h.py) - state lists and enumeration emission.
- [simple_loop.py](../../../nrpy/infrastructures/Dendro/simple_loop.py) - padded block point loops.
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - direct input and output pointer binding.
- [Ricci_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/Ricci_eval.py) - separate Ricci scratch binding.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)
- Used by: [BSSN Application Wiring](bssn-application-wiring.md)
