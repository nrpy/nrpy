# Initial Data

> Map Cartesian and spherical ADM initial-data providers used by GR equation consumers. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [General Relativity](index.md)

## Summary

NRPy's general-relativity initial-data modules construct ADM quantities, not
evolved BSSN RHSs. Both providers expose `gammaDD`, `KDD`, `alpha`, `betaU`,
and `BU`, plus the coordinate symbols used by the data. Downstream consumers can
feed these ADM fields into `ADM_to_BSSN` or generated application setup code.

The Cartesian provider supports `BrillLindquist` and `Kasner`. The spherical
provider supports `UIUCBlackHole`, `StaticTrumpet`, and `OffsetKerrSchild`.

## Detail

`InitialData_Cartesian` initializes `gammaDD`, `KDD`, `alpha`, `betaU`, `BU`,
and the Cartesian symbols `x`, `y`, and `z`. `BrillLindquist` registers two
black-hole positions and masses, builds a conformally flat spatial metric, and
sets zero extrinsic curvature. `Kasner` delegates to
`kasner_adm_quantities`, which constructs diagonal exact-Kasner ADM fields from
`KASNER_t0`, `KASNER_p1`, `KASNER_p2`, and `KASNER_p3`.

`InitialData_Spherical` initializes the same ADM outputs with spherical symbols
`r`, `th`, and `ph`. `UIUCBlackHole` registers `M` and `chi` and returns the
spinning-black-hole spatial metric and extrinsic curvature in the UIUC slicing.
`StaticTrumpet` returns spherical trumpet data with its own lapse and radial
shift. `OffsetKerrSchild` registers `M`, `a`, and `r0`, uses the offset
Kerr-Schild radius `r + r0`, and returns ADM metric, extrinsic curvature,
lapse, radial shift, and zero `BU`.

Gauge handling is explicit. If an initial-data family does not define gauge
quantities, or if `override_gauge_with_standard=True`, the provider sets
`betaU` and `BU` to zero and computes the standard lapse from
`ADM_to_BSSN(..., compute_cf_only=True)`. The lapse branch follows
`EvolvedConformalFactor_cf`: `phi` uses `exp(-2*cf)`, `chi` uses `sqrt(cf)`,
and `W` uses `cf`.

The trusted files validate the public object state for each supported family:
Cartesian Brill-Lindquist and Kasner, plus spherical UIUCBlackHole,
StaticTrumpet, and OffsetKerrSchild. The separate `ADM_to_BSSN_StaticTrumpet`
trusted file covers the handoff from spherical ADM initial data into BSSN
variables.

## Sources

- [InitialData_Cartesian.py](../../../nrpy/equations/general_relativity/InitialData_Cartesian.py) - `InitialData_Cartesian`, `BrillLindquist`, `Kasner`, `kasner_adm_quantities`
- [InitialData_Spherical.py](../../../nrpy/equations/general_relativity/InitialData_Spherical.py) - `InitialData_Spherical`, `UIUCBlackHole`, `StaticTrumpet`, `OffsetKerrSchild`
- [ADM_to_BSSN.py](../../../nrpy/equations/general_relativity/ADM_to_BSSN.py) - `ADM_to_BSSN`, `compute_cf_only`
- [InitialData_Cartesian_BrillLindquist.py](../../../nrpy/equations/general_relativity/tests/InitialData_Cartesian_BrillLindquist.py) - `trusted_dict`
- [InitialData_Cartesian_Kasner.py](../../../nrpy/equations/general_relativity/tests/InitialData_Cartesian_Kasner.py) - `trusted_dict`
- [InitialData_Spherical_UIUCBlackHole.py](../../../nrpy/equations/general_relativity/tests/InitialData_Spherical_UIUCBlackHole.py) - `trusted_dict`
- [InitialData_Spherical_StaticTrumpet.py](../../../nrpy/equations/general_relativity/tests/InitialData_Spherical_StaticTrumpet.py) - `trusted_dict`
- [InitialData_Spherical_OffsetKerrSchild.py](../../../nrpy/equations/general_relativity/tests/InitialData_Spherical_OffsetKerrSchild.py) - `trusted_dict`
- [ADM_to_BSSN_StaticTrumpet.py](../../../nrpy/equations/general_relativity/tests/ADM_to_BSSN_StaticTrumpet.py) - `trusted_dict`

## See Also

- [General Relativity](index.md)
- [Metric Conversions And Matter](metric-conversions-and-matter.md)
- [BSSN Family](bssn-family.md)
- [Fishbone-Moncrief](fishbone-moncrief.md)
- [Conformally Flat Elliptic](../conformally-flat-elliptic.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
