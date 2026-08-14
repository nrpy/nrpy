# Apparent Horizon Library

> Route BHaHAHA library generation and keep horizon-enabled examples linked to, not owners of, the apparent-horizon internals. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [Examples](index.md)

## Summary

`python -m nrpy.examples.bhahaha` is the owning example leaf for generating the
BHaHAHA apparent-horizon library. It builds the BHaHAHA C function registry,
writes a library project under the requested output root, emits the public
`BHaHAHA.h` header, and creates a static archive target named from the
finite-difference order.

The generator deletes and recreates its selected output directory. Preserve
prior library output before rerunning. Configured Ubuntu/macOS codegen jobs
generate and build the default static library; they do not run a solver or
establish horizon-finding accuracy, and workflow YAML is not latest-pass proof.

Horizon-enabled evolution examples call this generator as a dependency, then
register BHaH-facing horizon orchestration in their own generated applications.
They do not own BHaHAHA's public ABI, interpolation stages, relaxation loop, or
diagnostic vocabulary.

## Detail

The command shape is:

```bash
python -m nrpy.examples.bhahaha --fdorder 6 --outrootdir project
```

`--fdorder` sets the finite-difference order, with interpolation order one less
than that value. The default is `6`; non-default orders change the generated
project name from `BHaHAHA` to `BHaHAHA-<fdorder>o`, so `--fdorder 8` writes
`BHaHAHA-8o`. `--outrootdir` sets the output root directory, defaulting to
`project`. `--cpp` wraps the emitted public header in `extern "C"`, defines
`restrict` as `__restrict__`, and changes the fixed radii array parameter in
the header to an unsized array for C++ compatibility. `--no-openmp` is the
source-backed option that disables OpenMP flags in the generated Makefile.

The generator fixes `CoordSystem = "Spherical"`, uses `SSPRK33`, inner
boundaries only, reference-metric precomputation, no SIMD, and
`enable_fd_functions = False`. It registers BHaHAHA setup, poisoning,
relaxation, radial-grid, interpolation, diagnostics, error-message, local
boundary, RHS, KO, and BHaH define/header functions, then removes
`__rfm__Spherical` wrapper suffixes so the library exports plain `bah_`
functions.

The Makefile helper is called with `exec_or_library_name = "lib<project_name>"`,
`lib_function_prefix = "bah_"`, `create_lib = True`, and `static_lib = True`.
The resulting generated project is therefore a static library workflow, not a
standalone evolution executable. The final generator message tells users to
build in the generated BHaHAHA directory and start linking guidance from
`BHaHAHA.h`.

`BHaHAHA.h` is assembled from the checked-in header template plus generator
additions. It defines `REAL` when absent, declares the 12 Cartesian ADM input
gridfunction slots, provides `IDX2(itheta, iphi)` for flattened horizon-surface
storage, and appends `BHAHAHA_NGHOSTS` from half the finite-difference order.
It also appends `bhahaha_error_codes` from `error_code_msg_tuples_list` and the
public `bah_error_message()` prototype. The public header documents the caller
input struct, diagnostics struct, poisoning helpers, radial-grid setup,
center/radius helper, `bah_find_horizon()`, and diagnostic-file output entry
point; full field-level ownership stays in
[BHaHAHA Public API And Input Contract](../infrastructures/bhah/bhahaha-public-api-and-input-contract.md).

The horizon-enabled BHaH examples relate to this library by composition.
`two_blackholes_collide.py`, `blackhole_spectroscopy.py`,
`spinning_blackhole.py`, `superB_two_blackholes_collide.py`, and
`superB_blackhole_spectroscopy.py` are cross-link contexts for generated
applications that use or mirror apparent-horizon support. The standalone BHaH
black-hole examples run `bhahaha.py` into their project directory, link the
matching `BHaHAHA` or `BHaHAHA-<fdorder>o` archive, include `BHaHAHA.h`, and
register `bhahaha_find_horizons` with a source-specific horizon count. Their
initial-data, gauge, diagnostics, checkpointing, Psi4, and runtime ownership
belongs in [Standalone GR/BHaH](standalone-gr-bhah.md), while BHaHAHA solver
internals belong under [BHaHAHA Horizon Runtime](../infrastructures/bhah/bhahaha-horizon-runtime.md).

## Sources

- [bhahaha.py](../../nrpy/examples/bhahaha.py) - `--fdorder`, `--outrootdir`, `--cpp`, `--no-openmp`, `project_name`, `output_CFunctions_function_prototypes_and_construct_Makefile`
- [BHaHAHA_header.h](../../nrpy/infrastructures/BHaH/BHaHAHA/BHaHAHA_header.h) - `bhahaha_params_and_data_struct`, `bhahaha_diagnostics_struct`, `NUM_EXT_INPUT_CARTESIAN_GFS`, `IDX2`, public `bah_*` prototypes
- [error_message.py](../../nrpy/infrastructures/BHaH/BHaHAHA/error_message.py) - `error_code_msg_tuples_list`, `register_CFunction_error_message`
- [two_blackholes_collide.py](../../nrpy/examples/two_blackholes_collide.py) - `BHaH.BHaHAHA.BHaH_implementation.register_CFunction_bhahaha_find_horizons`
- [spinning_blackhole.py](../../nrpy/examples/spinning_blackhole.py) - `BHaH.BHaHAHA.BHaH_implementation.register_CFunction_bhahaha_find_horizons`
- [blackhole_spectroscopy.py](../../nrpy/examples/blackhole_spectroscopy.py) - `BHaH.BHaHAHA.BHaH_implementation.register_CFunction_bhahaha_find_horizons`
- [superB_two_blackholes_collide.py](../../nrpy/examples/superB_two_blackholes_collide.py) - apparent-horizon cross-link context
- [superB_blackhole_spectroscopy.py](../../nrpy/examples/superB_blackhole_spectroscopy.py) - apparent-horizon cross-link context

## See Also

- [Examples](index.md)
- [Standalone GR/BHaH](standalone-gr-bhah.md)
- [BHaHAHA Horizon Runtime](../infrastructures/bhah/bhahaha-horizon-runtime.md)
- [BHaHAHA Public API And Input Contract](../infrastructures/bhah/bhahaha-public-api-and-input-contract.md)
- [BHaHAHA Grid, Interpolation, And Boundaries](../infrastructures/bhah/bhahaha-grid-interpolation-and-boundaries.md)
