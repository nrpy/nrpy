# CarpetX Code Parameters, Includes, And Loops

> CarpetX-local helpers for Cactus parameter reads, standard includes, and generated Loop kernels. · Status: confirmed · Last reconciled: 2026-06-30
> Up: [CarpetX](index.md)

## Summary

CarpetX generated C functions use infrastructure-local helpers for Cactus
runtime parameters, standard include lists, and CarpetX Loop wrapper text. The
parameter helper has two distinct paths: scalar Cactus parameters use
`DECLARE_CCTK_PARAMETERS`, while SIMD-compatible real parameters are read
through `CCTK_ParameterGet()` and converted with `ConstSIMD()`. Loop wrapper
generation is CarpetX-specific and emits `grid.loop_*` calls with
`Loop::PointDesc` lambdas, not ETLegacy-style nested C loops.

## Detail

`read_CodeParameters()` receives `(thorn_name, code_parameter_name)` string
tuples, sorts them, and checks each code-parameter name against
`par.glb_code_params_dict`. In the scalar
non-SIMD code-parameter path it emits `DECLARE_CCTK_PARAMETERS` before the
tuple-driven declarations. In the SIMD path for `CCTK_REAL` and `REAL`
parameters, it emits a scalar pointer from
`CCTK_ParameterGet("<name>", "<thorn>", NULL)`, then wraps the dereferenced
value as a `REAL_SIMD_ARRAY` with `ConstSIMD()`. Non-real, non-char registered
types are emitted as `const <type> <name> = CCTK_ParameterGet(...)` using the
registered `CodeParameter.module`.

The inverse grid spacings `invdxx0`, `invdxx1`, and `invdxx2` are special
cases. They are generated from `1.0/CCTK_DELTA_SPACE(<dir>)` rather than from
`CCTK_ParameterGet()`. With SIMD enabled, each inverse spacing also gets a
`NOSIMD*` scalar and a `ConstSIMD()` vector declaration. Any registered
CodeParameter whose type string contains `char` raises `ValueError` before
either scalar or SIMD declarations are emitted.

`define_standard_includes()` centralizes the default CarpetX include list:
`loop_device.hxx`, `math.h`, `cctk.h`, `cctk_Arguments.h`, and
`cctk_Parameters.h`. Registration sites can request this list instead of
duplicating the common Cactus and CarpetX Loop headers.

`simple_loop()` emits CarpetX Loop calls. For `loop_region="all points"` the
base call is `grid.loop_all`; with the default `run_on_device=True`, the emitted
call is `grid.loop_all_device`. For `loop_region="interior"` the base call is
`grid.loop_int`; with the default device path it emits `grid.loop_int_device`,
while `run_on_device=False` emits `grid.loop_int`. Unsupported regions raise
`ValueError`. The emitted lambda captures by value and receives
`const Loop::PointDesc &p`, with `CCTK_DEVICE` for device loops and `CCTK_HOST`
for host loops.

Loop centering is a three-integer template argument list. `None` defaults to
`[0, 0, 0]`; otherwise callers must provide exactly three integers. Each entry
must be `0` for vertex centering or `1` for cell centering. CarpetX SIMD loop
generation is not implemented in this helper: `simple_loop(enable_simd=True)`
raises `ValueError` and asks callers to generate with `enable_simd=False`.

## Sources

- [CodeParameters.py](../../../nrpy/infrastructures/CarpetX/CodeParameters.py) - `read_CodeParameters`
- [CarpetX_include_header.py](../../../nrpy/infrastructures/CarpetX/CarpetX_include_header.py) - `define_standard_includes`
- [simple_loop.py](../../../nrpy/infrastructures/CarpetX/simple_loop.py) - `simple_loop`

## See Also

- [CarpetX](index.md)
- [Thorn Assembly, Configuration, And CCL Files](thorn-assembly-configuration-and-ccl-files.md)
- [Boundaries And RHS Initialization](boundaries-and-rhs-initialization.md)
- [C Codegen](../../core/c-codegen.md)
- [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)
- [Loop Kernel And Device Helpers](../../core/helpers/loop-kernel-and-device-helpers.md)
- [SIMD And Intrinsic Support](../../core/helpers/simd-and-intrinsic-support.md)
- [ETLegacy Code Parameters, Includes, And Loops](../etlegacy/code-parameters-includes-and-loops.md)
