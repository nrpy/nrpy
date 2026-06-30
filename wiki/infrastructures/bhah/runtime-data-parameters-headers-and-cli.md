# Runtime Data, Parameters, Headers, And CLI

> BHaH route for generated runtime structs, parameter defaults, headers, parfiles, and command-line overrides. Status: confirmed. Last reconciled: 2026-06-29
> Up: [BHaH](index.md)

## Summary

BHaH runtime state is assembled from NRPy `CodeParameter` registrations, extra
struct contributions, generated headers, default-setter C functions, and a
parameter-file parser. `commondata_struct` stores simulation-wide runtime data,
`params_struct` stores per-grid parameters, and `griddata_struct` stores per-grid
pointers and module-owned substructures such as MoL, reference-metric
precompute, and curvilinear boundary data.

## Detail

`output_BHaH_defines_h()` builds `BHaH_defines.h`. It registers general includes
and macros, emits `REAL` and `DOUBLE`, converts `#define`-typed CodeParameters
into C preprocessor definitions, registers `params_struct` and
`commondata_struct`, adds finite-difference macros such as `NGHOSTS` when finite
differencing is present, and appends grid macros plus `griddata_struct`. The
struct split is controlled by each CodeParameter's `commondata` flag:
non-`#define` parameters with `commondata=True` become members of
`commondata_struct`; the rest become members of `params_struct`.

`griddata_struct` always contains `REAL *xx[3]` and then collects registered
per-grid declarations from `griddata_commondata`. The BHaH core registers
`params_struct params`; optional modules register their own entries, such as
`rfm_struct* rfmstruct`, `bc_struct bcstruct`, and
`MoL_gridfunctions_struct gridfuncs`. `griddata_commondata` also owns
`griddata_free()` and, for non-OpenMP builds, `griddata_free_device()`, freeing
reference-metric tables, curvilinear boundary arrays, MoL gridfunction storage,
auxevol storage, coordinate arrays, and the `griddata` array according to enabled
features.

`CodeParameters.py` registers two generated default setters:
`commondata_struct_set_to_default()` zeroes `commondata` and fills defaults for
simulation-wide parameters, while `params_struct_set_to_default()` loops over
`MAXNUMGRIDS` and fills each grid's `params`. Scalars, strings, and `REAL` or
`int` arrays have separate assignment paths. The same module writes generated
parameter-access fragments: `set_CodeParameters.h`,
`set_CodeParameters-nopointer.h`, and `set_CodeParameters-simd.h`. These are the
BHaH source-backed parameter include files used by generated C functions; the
local source does not show a separate BHaH writer for a literal
`CodeParameters.h` file.

`cmdline_input_and_parfiles.py` generates the runtime parser
`cmdline_input_and_parfile_parser(commondata, argc, argv)`. The parser considers
only CodeParameters that are both `commondata=True` and `add_to_parfile=True`.
It builds a parameter descriptor table, strips comments and whitespace, accepts
scalar and array values, rejects duplicate definitions, warns on unrecognized
parameters, and writes parsed values into `commondata`. The default parfile
writer emits a project-named parfile, grouped by CodeParameter module, with only
commondata parfile parameters. Command-line steerable parameters are explicitly
validated: duplicate names, unsupported names, and array parameters are rejected.
Runtime invocation may read the default parfile, read a named parfile, or read a
parfile and then override approved steerable scalar/string parameters from the
last command-line arguments.

CUDA builds add device-facing runtime state. `output_device_headers()` writes
`BHaH_device_defines.h`, `BHaH_global_device_defines.h`, and
`BHaH_CUDA_global_init.h` only when `parallelization == "cuda"`. The device
header declares `__constant__ params_struct d_params[NUM_STREAMS]`,
`__constant__ commondata_struct d_commondata`, CUDA streams, parity arrays, and
optional wavespeed or `f_infinity` arrays. CUDA runtime helpers copy
`commondata` and per-grid `params` to constant memory, copy gridfunction buffers
host to device or device to host asynchronously, copy coordinate arrays back to a
host mirror grid, and copy `bc_struct` metadata and inner/outer boundary arrays
to device storage.

## Sources

- [BHaH_defines_h.py](../../../nrpy/infrastructures/BHaH/BHaH_defines_h.py) - `output_BHaH_defines_h`, `register_BHaH_defines`, `parse_cparam_type`, `_register_param_structs`, `register_griddata_struct_and_return_griddata_struct_str`
- [BHaH_device_defines_h.py](../../../nrpy/infrastructures/BHaH/BHaH_device_defines_h.py) - `CUDA_BHaH_device_defines_h`, `BHaH_CUDA_global_init_h`, `BHaH_CUDA_global_defines_h`, `output_device_headers`
- [CodeParameters.py](../../../nrpy/infrastructures/BHaH/CodeParameters.py) - `register_CFunctions_params_commondata_struct_set_to_default`, `write_CodeParameters_h_files`
- [cmdline_input_and_parfiles.py](../../../nrpy/infrastructures/BHaH/cmdline_input_and_parfiles.py) - `register_CFunction_cmdline_input_and_parfile_parser`, `generate_default_parfile`
- [griddata_commondata.py](../../../nrpy/infrastructures/BHaH/griddata_commondata.py) - `GridCommonData`, `register_griddata_commondata`, `register_CFunction_griddata_free`, `register_CFunction_griddata_free__device`
- [cuda_utilities.py](../../../nrpy/infrastructures/BHaH/parallelization/cuda_utilities.py) - `register_CFunctions_HostDevice__operations`, `register_CFunction_cpyHosttoDevice_commondata__constant`, `register_CFunction_cpyHosttoDevice_params__constant`, `register_CFunction_cpyHosttoDevice_bc_struct`, `register_CFunction_cpyDevicetoHost__grid`, `register_CFunction_cpyDevicetoHost__gf`, `register_CFunction_cpyHosttoDevice__gf`
- [c_function.py](../../../nrpy/c_function.py) - `CFunction`, `generate_full_function`

## See Also

- [BHaH](index.md)
- [Lifecycle And Project Assembly](lifecycle-and-project-assembly.md)
- [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)
