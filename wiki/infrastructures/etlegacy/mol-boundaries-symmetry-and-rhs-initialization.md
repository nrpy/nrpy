# MoL, Boundaries, Symmetry, And RHS Initialization

> ETLegacy leaf for MoL registration, old CartGrid3D symmetry parity, boundary setup, NewRad RHS boundaries, and RHS zeroing. Status: confirmed. Last reconciled: 2026-06-30
> Up: [ETLegacy](index.md)

## Summary

ETLegacy registers runtime services around generated equation kernels instead of
owning the Einstein Toolkit runtime itself. MoL registration binds the generated
`evol_variables` group to the generated `evol_variables_rhs` group. Symmetry
registration emits old CartGrid3D parity declarations from registered
gridfunction names. RHS initialization zeros every evolved RHS gridfunction at
`BASEGRID`. Boundary setup is split between Driver, Boundary, auxiliary,
evolved, and NewRad registration paths.

## Detail

`register_CFunction_MoL_registration()` emits one thorn-local C function named
`<thorn>_MoL_registration`. The body looks up
`<thorn>::evol_variables` and `<thorn>::evol_variables_rhs` with
`CCTK_GroupIndex()` and passes the pair to `MoLRegisterEvolvedGroup()`, so MoL
knows which RHS group advances which evolved-variable group. The generated
function is scheduled in `MoL_Register` with `OPTIONS: META`.

`register_CFunction_Symmetry_registration_oldCartGrid3D()` emits
`<thorn>_Symmetry_registration_oldCartGrid3D` at `BASEGRID` as
`Symmetry_registration`. It iterates over registered `EVOL`, `AUXEVOL`, and
`AUX` gridfunctions, defaults parity to scalar `+1` across the three coordinate
planes, then flips parity entries from name suffix digits. One trailing digit is
treated as rank 1, two trailing digits as rank 2, and a name with a digit three
places from the end is rejected as unsupported rank 3 or higher. For
four-dimensional rank-2 naming, the emitted parity code subtracts one from each
component index before mapping to the three spatial parity slots. Each accepted
gridfunction ends with a `SetCartSymVN(cctkGH, sym, "<thorn>::<gfname>GF")`
call.

`register_CFunction_zero_rhss()` emits `<thorn>_zero_rhss` at `BASEGRID` after
`Symmetry_registration`. For each registered `EVOL` gridfunction, it writes the
matching `<gfname>_rhs` ETLegacy gridfunction access to `0.0`. The body is
wrapped in ETLegacy `simple_loop()` with `loop_region="all points"`,
OpenMP enabled, and SIMD disabled, and the schedule declares
`WRITES: evol_variables_rhs(everywhere)`.

`boundary_conditions.py` has a dispatcher, `register_CFunctions()`, that
registers four separate paths. `register_CFunction_specify_Driver_BoundaryConditions()`
schedules a meta function in `Driver_BoundarySelect`. It selects every `EVOL`
gridfunction for Driver boundary type `"none"` with width `1`, and every `AUX`
gridfunction for type `"flat"` with `bndsize = fd_order / 2 + 1`. This is the
only boundary path in this module that combines evolved and auxiliary variables
in one generated function.

`register_CFunction_specify_evol_BoundaryConditions()` schedules evolved
Boundary registration in `MoL_PostStep` with `OPTIONS: level` and
`SYNC: evol_variables`, then schedules Cactus `ApplyBCs` after that function as
`<thorn>_evol_ApplyBCs`. It selects each `EVOL` gridfunction for Boundary type
`"none"` with width `1`; the local description explains that evolved outer
boundary handling is delegated to NewRad because NewRad modifies RHSs.

`register_CFunction_specify_aux_BoundaryConditions()` schedules auxiliary
Boundary registration in `MoL_PseudoEvolution` after
`<thorn>_BSSN_constraints`, with `OPTIONS: level` and `SYNC: aux_variables`.
It then schedules `ApplyBCs` as `<thorn>_aux_ApplyBCs` after the registration
function. It selects each `AUX` gridfunction for Boundary type `"flat"` with
`bndsize = fd_order / 2 + 1`. `AUXEVOL` gridfunctions are used by the symmetry
registration path above, but this boundary module does not select `AUXEVOL`
gridfunctions in its Driver, Boundary, or NewRad loops.

`register_CFunction_specify_NewRad_BoundaryConditions_parameters()` schedules
NewRad setup in `MoL_CalcRHS` after `<thorn>_RHS`. For each registered `EVOL`
gridfunction that is an `ETLegacyGridFunction`, it emits
`NewRad_Apply(cctkGH, <gfname>GF, <gfname>_rhsGF, f_infinity, wavespeed, 1.0)`.
The infinity value and characteristic speed come from the gridfunction's
`f_infinity` and `wavespeed` metadata, and the radial power is fixed to `1.0`
in this generator. The schedule declares reads from `evol_variables(everywhere)`
and writes to `evol_variables_rhs(boundary)`.

## Sources

- [MoL_registration.py](../../../nrpy/infrastructures/ETLegacy/MoL_registration.py) - `register_CFunction_MoL_registration`
- [Symmetry_registration.py](../../../nrpy/infrastructures/ETLegacy/Symmetry_registration.py) - `register_CFunction_Symmetry_registration_oldCartGrid3D`
- [zero_rhss.py](../../../nrpy/infrastructures/ETLegacy/zero_rhss.py) - `register_CFunction_zero_rhss`
- [boundary_conditions.py](../../../nrpy/infrastructures/ETLegacy/boundary_conditions.py) - `register_CFunction_specify_Driver_BoundaryConditions`, `register_CFunction_specify_evol_BoundaryConditions`, `register_CFunction_specify_aux_BoundaryConditions`, `register_CFunction_specify_NewRad_BoundaryConditions_parameters`, `register_CFunctions`
- [grid.py](../../../nrpy/grid.py) - `GridFunction`, `ETLegacyGridFunction`

## See Also

- [ETLegacy](index.md)
- [Code Parameters, Includes, And Loops](code-parameters-includes-and-loops.md)
- [GR BSSN RHS, Ricci, Constraints, And Validation](gr-bssn-rhs-ricci-constraints-and-validation.md)
- [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)
- [Generated Project CI](../../validation/generated-project-ci.md)
