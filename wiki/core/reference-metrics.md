# Reference Metrics

> Core route for coordinate-system reference metrics and BHaH precompute support. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Core APIs](index.md)

## Summary

`ReferenceMetric` builds and stores coordinate-system-specific reference-metric data used by NRPy equations and infrastructures. It sets coordinate maps, scale factors, reference metric tensors, inverse metrics, determinant derivatives, rescaling factors, and reference Christoffel symbols. The same source supports optional BHaH reference-metric precomputation, and representative trusted dictionaries under `nrpy/tests` validate generated reference-metric quantities for supported coordinate systems.

## Detail

`ReferenceMetric` is initialized with a coordinate-system name and an `enable_rfm_precompute` flag. It registers shared coordinate-system parameters, creates symbolic grid coordinates `xx`, declares Cartesian symbols, and initializes mappings such as `xx_to_Cart`, `Cart_to_xx`, `xxSph`, orthogonal scale factors, unit vectors, and rescaling arrays. Coordinate-system dispatch selects Cartesian-like, spherical-like, wedge, prolate-spheroidal, cylindrical, or `GeneralRFM` setup paths, raising `ValueError` for unknown coordinate systems.

After coordinate setup, the class computes Jacobians to Cartesian and spherical coordinates, the reference metric `ghatDD`, inverse metric `ghatUU`, determinant `detgammahat`, determinant derivatives, rescaling quantities and their derivatives, reference-metric derivatives, `GammahatUDD`, and `GammahatUDDdD`. `GeneralRFM` differs from ordinary diagonal reference metrics: it registers `ghatDD`, `ghatDDdD`, and `ghatDDdDD` as `AUXEVOL` gridfunctions and then computes derived reference-metric quantities algebraically.

When precompute is enabled for non-`GeneralRFM` systems, expensive hatted quantities are rewritten from generic SymPy function forms into NRPy variable names and collected for generated storage. `ReferenceMetricPrecompute` consumes `refmetric.reference_metric[CoordSystem + "_rfm_precompute"]`, discovers precomputed expressions, emits struct members, allocation/free specifications, reader strings, and host/CUDA populate code. `register_CFunctions_rfm_precompute()` registers coordinate-specific `rfm_precompute_malloc`, `rfm_precompute_defines`, and `rfm_precompute_free` functions and adds the `rfm_struct` typedef to BHaH defines.

The files under `nrpy/tests/reference_metric_*.py` are representative trusted-value dictionaries. They contain only an `mpf` import and `trusted_dict`, matching the project rule that trusted test-vector files are generated evidence rather than prose or hand-authored logic.

## Sources

- [nrpy/reference_metric.py](../../nrpy/reference_metric.py) - `ReferenceMetric`, `reference_metric`, `supported_CoordSystems`
- [nrpy/infrastructures/BHaH/rfm_precompute.py](../../nrpy/infrastructures/BHaH/rfm_precompute.py) - `ReferenceMetricPrecompute`, `register_CFunctions_rfm_precompute`
- [nrpy/tests/reference_metric_Cartesian.py](../../nrpy/tests/reference_metric_Cartesian.py) - `trusted_dict`
- [nrpy/tests/reference_metric_Spherical.py](../../nrpy/tests/reference_metric_Spherical.py) - `trusted_dict`
- [nrpy/tests/reference_metric_GeneralRFM_fisheyeN2.py](../../nrpy/tests/reference_metric_GeneralRFM_fisheyeN2.py) - `trusted_dict`
- [coding_style.md](../../coding_style.md) - trusted vector file contract

## See Also

- [Core APIs](index.md)
- [Indexed Expressions](indexed-expressions.md)
- [Gridfunctions And Parameters](gridfunctions-and-parameters.md)
