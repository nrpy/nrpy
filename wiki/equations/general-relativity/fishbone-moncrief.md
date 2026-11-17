# Fishbone-Moncrief

> Map Fishbone-Moncrief torus initial data for GRHD and GRMHD calculations. · Status: confirmed
> Up: [General Relativity](index.md)

## Summary

`FishboneMoncriefID` constructs spherical-coordinate Fishbone-Moncrief torus
initial data around a spinning black hole. It stores ADM fields, fluid fields,
Valencia velocity, and magnetic-field quantities that are natural inputs for
GRHD or GRMHD initial data, while leaving coordinate-system and Jacobian
mapping outside this module.

## Detail

The class registers runtime code parameters `r_in`, `r_at_max_density`, `a`,
`M`, `kappa`, `Gamma`, and `A_b`. It defines spherical symbols `r`, `th`, and
`ph`, then initializes output containers for ADM data (`gammaDD`, `KDD`,
`alpha`, `betaU`, and `BU`), hydrodynamic data (`rho_initial`,
`Pressure_initial`, `LorentzFactor`, and `Valencia3velocityU`), and magnetic
data (`BmagU`, `BtildeU`, and `smallb2`).

`_calculate_l_at_r` computes the specific angular momentum at a supplied
Boyer-Lindquist radius. `_compute_initial_data` then builds the torus enthalpy,
density, pressure, Boyer-Lindquist four-velocity, transformed Kerr-Schild
four-velocity, Kerr-Schild ADM metric fields, extrinsic curvature, vector
potential, magnetic field, comoving magnetic invariant, and Valencia velocity.
It stores these outputs as class attributes and sets `BU` to zero. The density
expression evaluates
`exp(log(hm1*(Gamma-1)/(kappa*Gamma))/(Gamma-1))` directly. It supplies no
torus-support mask, atmosphere, or positivity enforcement. Callers must select
admissible parameters and torus coordinates and supply exterior treatment;
using the unmasked expressions outside their physical domain can yield
negative or complex fluid quantities.

Evidence (descriptive behavior): `general_relativity/fishbone_moncrief/fishbone_moncrief.py::FishboneMoncriefID._compute_initial_data`.

`FishboneMoncriefID` constructs equilibrium torus initial data in spherical `(r, th, ph)` coordinates; external coordinate transformations map those fields to target evolution grids.

`FishboneMoncriefID.BtildeU` is the coordinate curl of the vector potential
`A_i`, the densitized field `sqrt(gamma) B^i`, and
`FishboneMoncriefID.BmagU = BtildeU / sqrt(gamma)`. Both hold spherical
Kerr-Schild components without division by `sqrt(4*pi)`. Use `BmagU`, not
`BtildeU`, as GRMHD input. Inside `_compute_initial_data`, the local comoving
four-vector `smallbU` includes that division; the stored `smallb2` follows the
same convention. To use `FishboneMoncriefID.BmagU` with the
[GRMHD](../grmhd.md) speed or HLL functions, transform its contravariant
components into the same basis as the supplied face metric and
four-velocities, then divide by `sqrt(4*pi)`. For `GRMHD_Equations`, transform
to the evolved `CoordSystem` basis. If `F^i` denotes the transformed
Fishbone-Moncrief field, substitute `F^i / (sqrt(4*pi) * ReU[i])` for each
`rescaledBmagU` symbol. To assign expressions after constructing the class,
instead set `BmagU[i] = F^i / sqrt(4*pi)` before computing `T4UU`; changing
`rescaledBmagU` after construction does not recompute `BmagU`. For the same
metric, shift, lapse, four-velocity, and field at the same point, all in one
basis, GRMHD `smallb2` should equal the stored Fishbone-Moncrief `smallb2`; no
trusted test performs this comparison.

Evidence (descriptive behavior): `general_relativity/fishbone_moncrief/fishbone_moncrief.py::FishboneMoncriefID._compute_initial_data`; `grmhd/GRMHD_equations.py::GRMHD_Equations.__init__` and `compute_smallb4U`; corroboration: `grmhd/HLL_fluxes.py::calculate_HLL_fluxes` magnetic input scaling.

The representative trusted dictionary records outputs such as `rho_initial`,
`Pressure_initial`, `alpha`, `betaU`, `gammaDD`, `KDD`,
`Valencia3velocityU`, `BtildeU`, `BmagU`, `LorentzFactor`, and `smallb2` from a
`FishboneMoncriefID` instance.

## Sources

- [fishbone_moncrief.py](../../../nrpy/equations/general_relativity/fishbone_moncrief/fishbone_moncrief.py) - `FishboneMoncriefID`, `_calculate_l_at_r`, `_compute_initial_data`
- [fishbone_moncrief.py](../../../nrpy/equations/general_relativity/fishbone_moncrief/tests/fishbone_moncrief.py) - `trusted_dict`
- [GRMHD_equations.py](../../../nrpy/equations/grmhd/GRMHD_equations.py) - `GRMHD_Equations.__init__`, `compute_smallb4U`
- [HLL_fluxes.py](../../../nrpy/equations/grmhd/HLL_fluxes.py) - `calculate_HLL_fluxes`

## See Also

- See also: [Initial Data](initial-data.md)
- Depends on: [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)
- Validated by: [Expression Validation Helpers](../../validation/expression-validation-helpers.md)
