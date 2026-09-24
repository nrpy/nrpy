# Gridfunctions, Naming, And Loops

> Define Dendro storage names, fixed component order, and block-local loop generation. · Status: provisional
> Up: [Dendro](index.md)

## Summary

`DendroGridFunction` supplies the Dendro C++ storage expression directly. The
Dendro connector maps NRPy's canonical conformal-factor gridfunction `cf` to
`cf_W_or_chi`; equation modules retain `cf`. Stencil registrars bind pointers
from fixed component lists and emit direct `ot::Block` loops. Algebraic floors
and projection instead operate on owned nodes of zipped evolved vectors.

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
a nested OpenMP region. Stencil kernels receive `ot::Block` directly and
derive offsets, dimensions, spacing, padding, and interior bounds from it.
The floor and determinant/trace projection kernels take zipped field pointers
and the mesh-owned node range instead of block geometry.

Claim evidence:
- Claim: Dendro's algebraic floor and determinant/trace projection kernels take zipped field pointers and an owned-node range, while stencil kernels retain padded-block geometry.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py` and `enforce_detgbar_equals_detghat_trAzero.py`, CFunction registrations.
- Corroboration: `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::post_timestep` within `output_solver_context_cpp`, passes the zipped stage and owned-node range.

## Sources

- [grid.py](../../../nrpy/grid.py) - `DendroGridFunction` and `dendro_name`.
- [state_h.py](../../../nrpy/infrastructures/Dendro/state_h.py) - state lists and enumeration emission.
- [simple_loop.py](../../../nrpy/infrastructures/Dendro/simple_loop.py) - padded block point loops.
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - direct input and output pointer binding.
- [Ricci_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/Ricci_eval.py) - separate Ricci scratch binding.
- [floor_the_lapse_and_conformal_factor.py](../../../nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py) - owned-node floor.
- [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/Dendro/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - owned-node projection.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)
- Used by: [BSSN Application Wiring](bssn-application-wiring.md)
