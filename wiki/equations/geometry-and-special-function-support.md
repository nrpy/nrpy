# Geometry And Special-Function Support

> Map equation-tree helpers for basis transforms, GeneralRFM fisheye maps, rotations, and spin-weighted spherical harmonics. · Status: confirmed · Last reconciled: 07-02-2026
> Up: [Equations](index.md)

## Summary

This page owns support modules under `nrpy/equations` that are not themselves
evolution systems: reference-metric basis transforms, GeneralRFM fisheye
coordinate maps, SO(3) rotation expressions, quaternion tensor rotation, and
spin-weighted spherical harmonics. The generic tensor containers, parameters,
and reference-metric registry remain owned by the Core API pages.

## Detail

`BasisTransforms` wraps an existing coordinate-system `ReferenceMetric` entry
and applies its stored Jacobians. It transforms contravariant vectors,
covariant vectors, rank-2 covariant tensors, and time-independent four-tensors
between the reference-metric basis and Cartesian or spherical bases. The lazy
`basis_transforms` dictionary constructs one `BasisTransforms` object per
coordinate system, and the module validation emits trusted dictionaries for
each supported coordinate system.

`GeneralRFMFisheye` builds an N-transition radial fisheye map from raw
Cartesian coordinates `xx[i]` to physical Cartesian coordinates. It registers
the plateau, transition-center, width, and global-scale code parameters,
stores `xx_to_CartU` and `dCart_dxxUD`, and constructs the induced flat
reference metric `ghatDD` plus first and second derivatives in the raw
coordinates. `build_fisheye` is the public constructor wrapper, and trusted
files cover the N=1 and N=2 variants.

`SO3Expressions` stores symbolic matrix-rotation expressions used by equation
validation and by the generated BHaH SO(3) helper layer. Its helpers build a
rotation matrix from orthonormal frame hats, apply `R` or `R^T` to vectors and
rank-2 covariant tensors, compute relative rotations, dot products, norms,
cross products, right-handed-frame invariants, Rodrigues axis-angle matrices,
and axis recovery branches from a rotation matrix. Its script entry point sends
the object dictionary through the trusted-expression pipeline. Generated BHaH C
helpers that wrap this symbolic layer are routed from
[BHaH SO(3) Rotation Helpers](../infrastructures/bhah/so3-rotation-helpers.md).

`tensor_rotation.rotate` is a lighter quaternion path for symbolic vectors and
3x3 tensors. It constructs a SymPy quaternion from an axis and angle, applies
`q v q*` for vectors, and applies the corresponding row and column rotations
for matrices. The module is doctest-only: its `__main__` block runs doctests
and does not generate a sibling trusted dictionary.

`spin_weighted_spherical_harmonics.Y` implements the Goldberg-formula path for
spin-weighted spherical harmonics. It registers `M_PI` for generated C contexts,
builds the finite sum with a Mathematica-compatible cotangent option or a
C-code-friendly tangent reciprocal, and returns the cached-simplified SymPy
expression. Its validation runs doctests and then generates trusted values for
representative `s=-2` harmonics.

## Sources

- [jacobians.py](../../nrpy/equations/basis_transforms/jacobians.py) - `BasisTransforms`, `basis_transforms`, `basis_transform_vectorU_from_rfmbasis_to_Cartesian`, `basis_transform_4tensorUU_from_Cartesian_to_time_indep_rfmbasis`
- [jacobians_Cartesian.py](../../nrpy/equations/basis_transforms/tests/jacobians_Cartesian.py) - `trusted_dict`
- [jacobians_Spherical.py](../../nrpy/equations/basis_transforms/tests/jacobians_Spherical.py) - `trusted_dict`
- [jacobians_GeneralRFM_fisheyeN2.py](../../nrpy/equations/basis_transforms/tests/jacobians_GeneralRFM_fisheyeN2.py) - `trusted_dict`
- [fisheye.py](../../nrpy/equations/generalrfm/fisheye.py) - `GeneralRFMFisheye`, `build_fisheye`, `xx_to_CartU`, `ghatDD`, `ghatDDdD`, `ghatDDdDD`
- [fisheye_N1.py](../../nrpy/equations/generalrfm/tests/fisheye_N1.py) - `trusted_dict`
- [fisheye_N2.py](../../nrpy/equations/generalrfm/tests/fisheye_N2.py) - `trusted_dict`
- [SO3_rotations.py](../../nrpy/equations/rotation/SO3_rotations.py) - `SO3Expressions`, `apply_R_to_vector`, `apply_R_to_tensorDD`, `rodrigues_matrix_from_axis_angle`
- [SO3_rotations.py](../../nrpy/equations/rotation/tests/SO3_rotations.py) - `trusted_dict`
- [tensor_rotation.py](../../nrpy/equations/quaternion_rotations/tensor_rotation.py) - `rotate`
- [spin_weighted_spherical_harmonics.py](../../nrpy/equations/special_functions/spin_weighted_spherical_harmonics.py) - `Y`
- [spin_weighted_spherical_harmonics.py](../../nrpy/equations/special_functions/tests/spin_weighted_spherical_harmonics.py) - `trusted_dict`
- [Goldberg formula reference](https://web2.ph.utexas.edu/~gsudama/pub/1967_008.pdf) - mathematical background for spin-weighted spherical harmonics
- [Spin-weighted functions with quaternions](https://pubs.aip.org/aip/jmp/article/57/9/092504/648118/How-should-spin-weighted-spherical-functions-be) - mathematical background
- [SO(3) rotation background](https://rotations.berkeley.edu/geodesics-of-the-rotation-group-so3/) - mathematical background

## See Also

- [Equations](index.md)
- [Wave Equation](wave-equation.md)
- [Conformally Flat Elliptic](conformally-flat-elliptic.md)
- [SEOBNR And BOB](seobnr/index.md)
- [Reference Metrics](../core/reference-metrics.md)
- [Indexed Expressions](../core/indexed-expressions.md)
- [Gridfunctions And Parameters](../core/gridfunctions-and-parameters.md)
- [Trusted Expression Pipeline](trusted-expression-pipeline.md)
