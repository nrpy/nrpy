# Fishbone-Moncrief

> Map Fishbone-Moncrief torus initial data for GRMHD-style consumers. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [General Relativity](index.md)

## Summary

`FishboneMoncriefID` constructs spherical-coordinate Fishbone-Moncrief torus
initial data around a spinning black hole. It stores ADM fields, fluid fields,
Valencia velocity, and magnetic-field quantities that are natural inputs for
GRHD or GRMHD initial-data consumers, while leaving coordinate-system and
Jacobian mapping outside this module.

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
At the end of the pipeline it stores all public outputs as class attributes and
sets `BU` to zero for this initial-data family.

The module-level note is part of the implementation boundary: data are
constructed in spherical `(r, th, ph)` coordinates, and coordinate mapping is
performed externally. This makes the page adjacent to GR initial data and GRHD
coverage, but the file itself is a torus initial-data provider rather than a
general GRHD evolution module.

The representative trusted dictionary records outputs such as `rho_initial`,
`Pressure_initial`, `alpha`, `betaU`, `gammaDD`, `KDD`,
`Valencia3velocityU`, `BtildeU`, `BmagU`, `LorentzFactor`, and `smallb2` from a
`FishboneMoncriefID` instance.

## Sources

- [fishbone_moncrief.py](../../../nrpy/equations/general_relativity/fishbone_moncrief/fishbone_moncrief.py) - `FishboneMoncriefID`, `_calculate_l_at_r`, `_compute_initial_data`
- [fishbone_moncrief.py](../../../nrpy/equations/general_relativity/fishbone_moncrief/tests/fishbone_moncrief.py) - `trusted_dict`

## See Also

- [General Relativity](index.md)
- [Initial Data](initial-data.md)
- [Metric Conversions And Matter](metric-conversions-and-matter.md)
- [GRHD](../grhd.md)
- [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
