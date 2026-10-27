# Source Map

> Aggregate ingest coverage and remaining scope for NRPy source trees. · Status: provisional
> Up: [NRPy Knowledge Base](../AGENTS.md)

## Scope

This table tracks ingest coverage across the multi-file source trees registered in [raw/SOURCES.md](../raw/SOURCES.md), following the authority tiers in [KB Rules](SCHEMA.md#citation-rules) and [Claim And Evidence Rules](SCHEMA.md#claim-and-evidence-rules). Reverse file-to-page dependencies are computed dynamically with `rg -l -F '<path>' wiki/`; external source tiers live in [raw/SOURCES.md](../raw/SOURCES.md).

## Aggregate Rows

| Source or aggregate | Glob pattern | Authority tier | Status | Ingest | Uncovered scope |
| --- | --- | --- | --- | --- | --- |
| `core-top-level-package-modules` | `nrpy/*.py`, `nrpy/*.txt`, `nrpy/py.typed` | primary-code | living | partial | Minor internal utility branches in top-level `nrpy/*.py` modules outside public codegen, parameter, grid, reference-metric, and FD entrypoints. |
| `helpers-validation-and-reference-metric-tests` | `nrpy/helpers/**`, `nrpy/validate_expressions/**`, `nrpy/tests/**` | primary-code; primary-test | living | partial | Individual `nrpy/tests/` reference-metric regression scripts not cited directly by [Reference Metrics](core/reference-metrics.md). |
| `equation-modules-and-trusted-values` | `nrpy/equations/**` | primary-code; primary-test | living | partial | Generated `*/tests/*.py` `trusted_dict` reference modules (cited via their calling equation modules) and lower-traffic Post-Newtonian helper subfunctions. |
| `infrastructure-modules-and-embedded-headers` | `nrpy/infrastructures/**` | primary-code; primary-test | living | partial | Auxiliary backend formatting helpers and standalone C/C++ micro-test files not individually cited by backend leaves. |
| `example-generators-and-companion-scripts` | `nrpy/examples/**` | primary-code; primary-test | living | partial | Non-CI plotting or post-processing helper scripts under `nrpy/examples/` beyond the top-level generators in [Example Generator Catalog](examples/example-generator-catalog.md). |
| `ci-and-local-automation` | `.github/**` | primary-test | living | partial | Auxiliary workflow steps outside `.github/workflows/main.yml`, `.github/workflows/dendrolib-canary.yml`, `.github/workflows/sebobv2-projected-attachment.yml`, `.github/workflows/spheroidal_fisheye_roundtrip.yml`, `.github/full_nrpy_local_ci.sh`, and `.github/single_file_static_analysis.sh`. |
| `nrpy/equations/generalrfm/inverse_spheroidal_fisheye.py` | primary-code | living | exact seed | [Geometry Support](equations/geometry-and-special-function-support.md); [Reference Metrics](core/reference-metrics.md); [BHaH Grids](infrastructures/bhah/grids-coordinates-reference-metrics-and-boundaries.md). | Inverse spheroidal map, physical normalization, BHaH geometry, and scoped numerical validation. | No CUDA support or evolution convergence claim. | 10-16-2026 | Recheck these owners when the map or integration contracts change. |
| `nrpy/infrastructures/BHaH/inverse_spheroidal_fisheye.py` | primary-code | living | exact seed | [Geometry Support](equations/geometry-and-special-function-support.md); [Reference Metrics](core/reference-metrics.md); [BHaH Grids](infrastructures/bhah/grids-coordinates-reference-metrics-and-boundaries.md). | Inverse spheroidal map, physical normalization, BHaH geometry, and scoped numerical validation. | No CUDA support or evolution convergence claim. | 10-16-2026 | Recheck these owners when the map or integration contracts change. |
| `.github/scripts/check_inverse_spheroidal_fisheye.py` | primary-test | living | exact seed | [Geometry Support](equations/geometry-and-special-function-support.md); [Reference Metrics](core/reference-metrics.md); [BHaH Grids](infrastructures/bhah/grids-coordinates-reference-metrics-and-boundaries.md). | Inverse spheroidal map, physical normalization, BHaH geometry, and scoped numerical validation. | No CUDA support or evolution convergence claim. | 10-16-2026 | Recheck these owners when the map or integration contracts change. |
