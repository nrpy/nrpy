# SO(3) Rotation Helpers

> Map the BHaH generated-C SO(3) helper layer for rotation matrices, hat vectors, and axis-angle recovery. Status: confirmed. Last reconciled: 07-02-2026
> Up: [BHaH](index.md)

## Summary

The BHaH rotation package wraps equation-side `SO3Expressions` from
`nrpy/equations/rotation/SO3_rotations.py` in eight public generated
CFunctions. `register_all.register_CFunctions()` registers helpers for building
rotation matrices from cumulative hats, converting axis-angle data to a matrix,
applying `R` or `R^T` to vectors, applying `R` to rank-2 covariant tensors,
computing a relative rotation, left-multiplying cumulative hats, and recovering
an axis-angle pair that maps one unit vector to another.

These helpers are generated under the BHaH CFunction `rotation/` subdirectory,
include `BHaH_defines.h`, and deliberately omit `set_CodeParameters.h` because
their APIs are small pure helper routines using `REAL`, fixed-size arrays, and
local scalar math.

## Detail

`register_CFunctions()` is the router for the generated helper set. It calls
the eight helper-specific registration functions:
`register_CFunction_so3_build_R_from_hats`,
`register_CFunction_so3_axis_angle_to_R`,
`register_CFunction_so3_apply_R_to_vector`,
`register_CFunction_so3_apply_RT_to_vector`,
`register_CFunction_so3_apply_R_to_tensorDD`,
`register_CFunction_so3_relative_R_dst_from_src`,
`register_CFunction_so3_left_multiply_hats_with_R`, and
`register_CFunction_so3_find_nU_and_dphi_from_unit_vectors`.

The common convention is that `R` maps rotating-frame components to fixed-frame
components. Thus `v_fixed = R v_rot`, while `R^T` maps fixed-frame components
back to the rotating frame. C matrix references are written as `R[i][j]`, where
`i` is the row and `j` is the column. Hat-vector routines treat `xhatU`,
`yhatU`, and `zhatU` as cumulative rotating-basis vectors expressed in fixed
components, and `so3_build_R_from_hats` packs those hats as the three columns of
`R`.

`so3_axis_angle_to_R` normalizes the input axis when possible and builds `R`
with the Rodrigues formula from `SO3Expressions`. Degenerate input is
deterministic: near-zero angle or near-zero axis norm returns the identity
matrix. `so3_find_nU_and_dphi_from_unit_vectors` performs the inverse
unit-vector task. It clamps the dot product, returns `(1,0,0), 0` for parallel
input, and chooses a deterministic orthogonal axis with angle `pi` for
antiparallel input.

The in-place vector, rank-2 tensor, and hat updates are alias-safe.
`so3_apply_R_to_vector`, `so3_apply_RT_to_vector`, and
`so3_apply_R_to_tensorDD` snapshot their input arrays before overwriting them.
The tensor helper applies `T_out = R T_in R^T` to all nine components. The
cumulative-hat update snapshots all three hats, interprets them as columns of
`R_old`, then applies `R_new = DeltaR * R_old`. Relative rotations use
`DeltaR_dst_from_src = R_dst^T R_src`, which maps source rotating-basis
components into destination rotating-basis components.

Each helper registration passes `subdirectory="rotation"` to
`register_CFunction`, so generated project output places these C functions
under the generated rotation helper directory. The individual modules validate
their generated OpenMP strings through `validate_strings` against trusted
`so3_*__openmp.c` files. The aggregate `register_all.py` script also clears
`CFunction_dict`, sets `parallelization` to `openmp`, registers all eight
helpers, and directly compares each generated `full_function` with its trusted
file under `nrpy/infrastructures/BHaH/rotation/tests/`.

## Sources

- [register_all.py](../../../nrpy/infrastructures/BHaH/rotation/register_all.py) - `register_CFunctions`
- [SO3_rotations.py](../../../nrpy/equations/rotation/SO3_rotations.py) - `SO3Expressions`, `rodrigues_matrix_from_axis_angle`
- [so3_build_R_from_hats.py](../../../nrpy/infrastructures/BHaH/rotation/so3_build_R_from_hats.py) - `register_CFunction_so3_build_R_from_hats`
- [so3_axis_angle_to_R.py](../../../nrpy/infrastructures/BHaH/rotation/so3_axis_angle_to_R.py) - `register_CFunction_so3_axis_angle_to_R`
- [so3_apply_R_to_vector.py](../../../nrpy/infrastructures/BHaH/rotation/so3_apply_R_to_vector.py) - `register_CFunction_so3_apply_R_to_vector`
- [so3_apply_RT_to_vector.py](../../../nrpy/infrastructures/BHaH/rotation/so3_apply_RT_to_vector.py) - `register_CFunction_so3_apply_RT_to_vector`
- [so3_apply_R_to_tensorDD.py](../../../nrpy/infrastructures/BHaH/rotation/so3_apply_R_to_tensorDD.py) - `register_CFunction_so3_apply_R_to_tensorDD`
- [so3_relative_R_dst_from_src.py](../../../nrpy/infrastructures/BHaH/rotation/so3_relative_R_dst_from_src.py) - `register_CFunction_so3_relative_R_dst_from_src`
- [so3_left_multiply_hats_with_R.py](../../../nrpy/infrastructures/BHaH/rotation/so3_left_multiply_hats_with_R.py) - `register_CFunction_so3_left_multiply_hats_with_R`
- [so3_find_nU_and_dphi_from_unit_vectors.py](../../../nrpy/infrastructures/BHaH/rotation/so3_find_nU_and_dphi_from_unit_vectors.py) - `register_CFunction_so3_find_nU_and_dphi_from_unit_vectors`
- [so3_build_R_from_hats_so3_build_R_from_hats__openmp.c](../../../nrpy/infrastructures/BHaH/rotation/tests/so3_build_R_from_hats_so3_build_R_from_hats__openmp.c) - `so3_build_R_from_hats`
- [so3_axis_angle_to_R_so3_axis_angle_to_R__openmp.c](../../../nrpy/infrastructures/BHaH/rotation/tests/so3_axis_angle_to_R_so3_axis_angle_to_R__openmp.c) - `so3_axis_angle_to_R`
- [so3_apply_R_to_vector_so3_apply_R_to_vector__openmp.c](../../../nrpy/infrastructures/BHaH/rotation/tests/so3_apply_R_to_vector_so3_apply_R_to_vector__openmp.c) - `so3_apply_R_to_vector`
- [so3_apply_RT_to_vector_so3_apply_RT_to_vector__openmp.c](../../../nrpy/infrastructures/BHaH/rotation/tests/so3_apply_RT_to_vector_so3_apply_RT_to_vector__openmp.c) - `so3_apply_RT_to_vector`
- [so3_apply_R_to_tensorDD_so3_apply_R_to_tensorDD__openmp.c](../../../nrpy/infrastructures/BHaH/rotation/tests/so3_apply_R_to_tensorDD_so3_apply_R_to_tensorDD__openmp.c) - `so3_apply_R_to_tensorDD`
- [so3_relative_R_dst_from_src_so3_relative_R_dst_from_src__openmp.c](../../../nrpy/infrastructures/BHaH/rotation/tests/so3_relative_R_dst_from_src_so3_relative_R_dst_from_src__openmp.c) - `so3_relative_R_dst_from_src`
- [so3_left_multiply_hats_with_R_so3_left_multiply_hats_with_R__openmp.c](../../../nrpy/infrastructures/BHaH/rotation/tests/so3_left_multiply_hats_with_R_so3_left_multiply_hats_with_R__openmp.c) - `so3_left_multiply_hats_with_R`
- [so3_find_nU_and_dphi_from_unit_vectors_so3_find_nU_and_dphi_from_unit_vectors__openmp.c](../../../nrpy/infrastructures/BHaH/rotation/tests/so3_find_nU_and_dphi_from_unit_vectors_so3_find_nU_and_dphi_from_unit_vectors__openmp.c) - `so3_find_nU_and_dphi_from_unit_vectors`

## See Also

- Parent: [BHaH](index.md)
- Depends on: [Geometry And Special-Function Support](../../equations/geometry-and-special-function-support.md)
- Depends on: [C Function Registry](../../core/c-function-registry.md)
- See also: [GR Application Wiring](gr-application-wiring.md)
- See also: [Lifecycle And Project Assembly](lifecycle-and-project-assembly.md)
