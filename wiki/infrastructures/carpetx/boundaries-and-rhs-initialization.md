# CarpetX Boundaries And RHS Initialization

> CarpetX-local NewRadX boundary registration and RHS zeroing behavior. · Status: confirmed · Last reconciled: 2026-06-30
> Up: [CarpetX](index.md)

## Summary

CarpetX boundary support in the local infrastructure is source-observed as a
NewRadX registration path for evolved variables plus a separate RHS zeroing
function. NewRadX setup includes `newradx.hxx`, applies
`NewRadX_Apply(cctkGH, var, rhs, f_infinity, wavespeed, 1.0)` to CarpetX
`EVOL` gridfunctions, and schedules immediately after the thorn RHS in the
`ODESolvers_RHS` bin. RHS zeroing registers `<thorn>_zero_rhss` at `BASEGRID`
after `Symmetry_registration` and writes the evolved RHS group everywhere.

## Detail

`register_CFunction_specify_NewRad_BoundaryConditions_parameters()` builds a
Cactus C function named
`<thorn>_specify_NewRad_BoundaryConditions_parameters`. Its include list adds
`newradx.hxx` alongside `math.h`, `cctk.h`, `cctk_Arguments.h`, and
`cctk_Parameters.h`. The generated body declares Cactus arguments and
parameters, enters `using namespace NewRadX`, and then iterates
`gri.glb_gridfcs_dict`.

The NewRadX loop selects gridfunctions whose group is `EVOL` and whose object
is a `CarpetXGridFunction`. For each such gridfunction, it emits
`NewRadX_Apply(cctkGH, <gfname>GF, <gfname>_rhsGF, <f_infinity>, <wavespeed>, 1.0)`.
The infinity value and wave speed come from the CarpetX gridfunction metadata;
the radial power argument is fixed to `1.0` in this generator.

The NewRadX schedule entry is registered in `ODESolvers_RHS after <thorn>_RHS`.
It declares `READS: evol_variables(everywhere)` and
`WRITES: evol_variables_rhs(boundary)`. The schedule bin ordering context is
owned by `construct_schedule_ccl()`, which emits `ODESolvers_RHS` after
`STARTUP`, `BASEGRID`, and `CCTK_INITIAL`, and before
`ODESolvers_PostStep` and remaining bins.

`register_CFunction_zero_rhss()` registers `<thorn>_zero_rhss`. It iterates
registered `EVOL` gridfunctions and emits assignments setting each matching
`<gfname>_rhs` access to `0.0`. The schedule entry places the function at
`BASEGRID after Symmetry_registration` and declares
`WRITES: evol_variables_rhs(everywhere)`.

Observed weak behavior: `zero_rhss.py` is inside the CarpetX infrastructure but
imports `nrpy.infrastructures.ETLegacy.simple_loop as lp` and uses
`gri.ETLegacyGridFunction.access_gf()` for RHS accesses. This page documents
that source behavior only; it is not a code-fix task.

A scoped inventory of `nrpy/infrastructures/CarpetX/` on 2026-06-30 shows no
local CarpetX files named for MoL registration, Symmetry registration, Driver
boundary selection, or separate Boundary registration beyond
`boundary_conditions.py`. The local CarpetX boundary-specific source in this
scope is the NewRadX-oriented `boundary_conditions.py` file described above.

## Sources

- [boundary_conditions.py](../../../nrpy/infrastructures/CarpetX/boundary_conditions.py) - `register_CFunction_specify_NewRad_BoundaryConditions_parameters`, `register_CFunctions`
- [zero_rhss.py](../../../nrpy/infrastructures/CarpetX/zero_rhss.py) - `register_CFunction_zero_rhss`
- [schedule_ccl.py](../../../nrpy/infrastructures/CarpetX/schedule_ccl.py) - `construct_schedule_ccl`
- [Source Manifest](../../../raw/SOURCES.md) - `carpetx-package-inventory` scoped inventory for absence claims

## See Also

- [CarpetX](index.md)
- [Thorn Assembly, Configuration, And CCL Files](thorn-assembly-configuration-and-ccl-files.md)
- [Code Parameters, Includes, And Loops](code-parameters-includes-and-loops.md)
- [ETLegacy MoL, Boundaries, Symmetry, And RHS Initialization](../etlegacy/mol-boundaries-symmetry-and-rhs-initialization.md)
