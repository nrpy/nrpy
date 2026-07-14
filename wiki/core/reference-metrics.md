# Reference Metrics

> Core route for coordinate-system reference metrics and precompute support. · Status: confirmed · Last reconciled: 07-13-2026
> Up: [Core APIs](index.md)

## Summary

`ReferenceMetric` builds and stores coordinate-system-specific reference-metric data used by NRPy equations and infrastructures. It sets coordinate maps, scale factors, reference metric tensors, inverse metrics, determinant derivatives, rescaling factors, and reference Christoffel symbols. `reference_metric` is a lazy `rfm_dict` cache that builds base and precompute variants on demand. The module-level validation path processes every supported coordinate system into trusted dictionaries and includes a separate inverse-map check for `SinhSymTP`.

## Detail

`ReferenceMetric` is initialized with a coordinate-system name and an `enable_rfm_precompute` flag. The module parameter `CoordSystem_to_register_CodeParameters` controls constructor CodeParameter registration: when this parameter equals one coordinate-system name, constructor-registered coordinate-system CodeParameters are added to the global CodeParameter dictionary only for that system, and when it equals `All`, all systems register globally. This guard prevents multipatch or multi-system setup from polluting the global CodeParameter namespace unintentionally.

Every constructor registers or uses the shared reference-metric parameters `enable_grid_physical_size`, `grid_physical_size`, `grid_hole_radius`, and `CoordSystemName`. It creates symbolic grid coordinates `xx`, declares Cartesian symbols, and initializes mappings such as `xx_to_Cart`, `Cart_to_xx`, `xxSph`, orthogonal scale factors, unit vectors, grid bounds, hole-radius dictionaries, physical-size dictionaries, and rescaling arrays.

Coordinate-system dispatch is by coordinate family. Names beginning with `GeneralRFM` use `general_rfm_like()`. Names containing `Cartesian`, `Spherical`, `Wedge`, `SymTP`, or `Cylindrical` route to Cartesian-like, spherical-like, wedge-like, SymTP/prolate-spheroidal, or cylindrical-like setup paths. Unknown names raise `ValueError`.

After coordinate setup, the class computes Jacobians to and from Cartesian and spherical coordinates; rescaling arrays `ReU`, `ReD`, and `ReDD`; the reference metric `ghatDD`; inverse metric `ghatUU`; determinant `detgammahat`; determinant derivatives `detgammahatdD` and `detgammahatdDD`; derivatives of the rescaling arrays; reference-metric derivatives `ghatDDdD` and `ghatDDdDD`; and reference Christoffel quantities `GammahatUDD` and `GammahatUDDdD`.

Supported coordinate systems are `Spherical`, `SinhSpherical`, `SinhSphericalv2n2`, `Cartesian`, `SinhCartesian`, `Cylindrical`, `SinhCylindrical`, `SinhCylindricalv2n2`, `SymTP`, `SinhSymTP`, `LWedgeHSinhSph`, `UWedgeHSinhSph`, `RingHoleySinhSpherical`, `HoleySinhSpherical`, and `GeneralRFM_fisheyeN2`. `unittest_CoordSystems` is a smaller regression subset: `SinhSymTP`, `HoleySinhSpherical`, `Cartesian`, and `SinhCylindricalv2n2`.

`rfm_dict.__getitem__()` lazily builds reference metrics. If a requested key is missing, it strips any `_rfm_precompute` suffix, builds the base `ReferenceMetric(CoordSystem, enable_rfm_precompute=False)`, caches it under the base key, and also caches a precompute key. Non-`GeneralRFM` systems get a separate `ReferenceMetric(CoordSystem, enable_rfm_precompute=True)` under `CoordSystem + "_rfm_precompute"`; `GeneralRFM` precompute keys map to the same base object because GeneralRFM data is already gridfunction-provided.

`GeneralRFM` behavior is distinct from diagonal coordinate families. `general_rfm_like()` sets neutral bounds and identity Cartesian labels by default, marks the analytic inverse map as unavailable with `nan` sentinels, stores provider metadata, and recognizes the `GeneralRFM_fisheyeN*` pattern by building a fisheye provider. For GeneralRFM systems, `ghatDD`, `ghatDDdD`, and `ghatDDdDD` are registered as `AUXEVOL` gridfunctions, rescaling factors are set to one, and derived quantities such as `ghatUU`, `detgammahat`, determinant derivatives, and `Gammahat*` are computed algebraically from those gridfunctions.

The `GeneralRFM_fisheyeN2` trusted dictionary contains three intentional NaN
sentinels, at `Cart_to_xx_0` through `Cart_to_xx_2`, for that unavailable
analytic inverse map. Trusted comparison accepts them only when the computed
values contain matching NaN components. A finite value at one of those keys, or
a computed NaN at a finite trusted key, is a validation failure; NaN is not a
blanket wildcard.

Claim evidence:
- Claim: The `GeneralRFM_fisheyeN2` trusted dictionary contains exactly three intentional NaN inverse-map sentinels, and validation accepts them only against matching NaN components.
- Role: descriptive behavior
- Deciding authority: [reference_metric_GeneralRFM_fisheyeN2.py](../../nrpy/tests/reference_metric_GeneralRFM_fisheyeN2.py), `trusted_dict`, and [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py), `compare_against_trusted` and `_nonfinite_values_match`
- Corroboration: [reference_metric.py](../../nrpy/reference_metric.py), `general_rfm_like`, assigns unavailable analytic inverse-map entries to `nan`
- Validation: `inspected=pass; generated=not-run; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0, mpmath 1.3.0; backend=not-applicable; precision=30 decimal digits; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=pass; options=full reference_metric validation including both trusted-NaN mismatch doctests; date=07-13-2026`

When precompute is enabled for non-`GeneralRFM` systems, hatted quantities are first expressed through generic function forms such as `f0_of_xx0_funcform`. The replacement pass rewrites functions and derivatives into rigid NRPy variable names such as base function names plus `__D...` derivative suffixes, folds derivatives that evaluate to coordinate-independent values, and skips creating a separate precompute object for `GeneralRFM`.

`register_pi()` and `register_sqrt1_2()` provide special constants as CodeParameters named `PI` and `SQRT1_2`, with values emitted into generated CodeParameters headers and excluded from parfiles.

In the module `__main__` validation path, doctests run first. Then every entry in `supported_CoordSystems` is processed through `validate_expressions.process_dictionary_of_expressions()` and compared or regenerated through `compare_or_generate_trusted_results()`. A separate `SinhSymTP` validation substitutes `Cart_to_xx` into `xx_to_Cart` and checks each Cartesian residual with `check_zero()`.

## Sources

- [nrpy/reference_metric.py](../../nrpy/reference_metric.py) - `ReferenceMetric`, `CoordSystem_to_register_CodeParameters`, `reference_metric`, `rfm_dict`
- [nrpy/reference_metric.py](../../nrpy/reference_metric.py) - `supported_CoordSystems`, `unittest_CoordSystems`
- [nrpy/reference_metric.py](../../nrpy/reference_metric.py) - `Sinhv1`, `Sinhv2`, `cartesian_like`, `spherical_like`, `spherical_wedge_like`, `prolate_spheroidal_like`, `cylindrical_like`, `general_rfm_like`
- [nrpy/reference_metric.py](../../nrpy/reference_metric.py) - `register_pi`, `register_sqrt1_2`
- [nrpy/reference_metric.py](../../nrpy/reference_metric.py) - `__main__` validation loop and `SinhSymTP` `check_zero` inverse check
- [nrpy/tests/reference_metric_GeneralRFM_fisheyeN2.py](../../nrpy/tests/reference_metric_GeneralRFM_fisheyeN2.py) - `trusted_dict`

## See Also

- Parent: [Core APIs](index.md)
- Depends on: [Indexed Expressions](indexed-expressions.md)
- Depends on: [Gridfunctions And Parameters](gridfunctions-and-parameters.md)
- Validated by: [Expression Validation Helpers](../validation/expression-validation-helpers.md)
