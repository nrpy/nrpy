"""
Register C helper for matrix-first unrotation of Cartesian vectors.

This hot path uses cumulative hats -> R directly and avoids
axis-angle two-stage solve machinery.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_function as cfc
from nrpy.infrastructures.BHaH.rotation.so3_matrix_ops import (
    register_CFunction_build_R_from_cumulative_hats,
    register_CFunction_so3_apply_R_to_vector,
    register_CFunction_so3_validate_and_optionally_fix_hats,
    register_rotation_commondata_CodeParameters,
)


def verify_unrotate_hot_path_antiparallel_regression() -> None:
    r"""
    Check the direct-matrix hot path on a deterministic anti-parallel basis case.

    This verifies numerical behavior for a 180-degree basis flip without any
    axis-angle decomposition in the hot path.

    Doctests:
    >>> from nrpy.equations.quaternion_rotations.tensor_rotation import rotate
    >>> xhat = [-1.0, 0.0, 0.0]
    >>> yhat = [0.0, 1.0, 0.0]
    >>> zhat = [0.0, 0.0, -1.0]
    >>> R = [[xhat[0], yhat[0], zhat[0]], [xhat[1], yhat[1], zhat[1]], [xhat[2], yhat[2], zhat[2]]]
    >>> v_rot = [1.25, -0.5, 2.0]
    >>> v_fixed = [sum(R[i][j] * v_rot[j] for j in range(3)) for i in range(3)]
    >>> v_fixed
    [-1.25, -0.5, -2.0]
    >>> v_expected = [float(v) for v in rotate(v_rot, [0.0, 1.0, 0.0], 3.141592653589793)]
    >>> max(abs(a - b) for a, b in zip(v_fixed, v_expected)) < 1e-14
    True
    """


def register_CFunction_unrotate_xCart_to_fixed_frame() -> None:
    r"""
    Register C function ``unrotate_xCart_to_fixed_frame``.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> from nrpy.infrastructures.BHaH.rotation.so3_matrix_ops import assert_so3_convention_in_text
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_unrotate_xCart_to_fixed_frame()
    >>> generated_str = cfc.CFunction_dict["unrotate_xCart_to_fixed_frame"].full_function
    >>> assert_so3_convention_in_text(generated_str, "unrotate_xCart_to_fixed_frame")
    >>> "unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame(" in generated_str
    False
    >>> "build_R_from_cumulative_hats(" in generated_str and "so3_apply_R_to_vector(" in generated_str
    True
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "openmp", file_ext="c")
    """
    register_rotation_commondata_CodeParameters()
    if "build_R_from_cumulative_hats" not in cfc.CFunction_dict:
        register_CFunction_build_R_from_cumulative_hats()
    if "so3_validate_and_optionally_fix_hats" not in cfc.CFunction_dict:
        register_CFunction_so3_validate_and_optionally_fix_hats()
    if "so3_apply_R_to_vector" not in cfc.CFunction_dict:
        register_CFunction_so3_apply_R_to_vector()

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r"""
@brief Map a Cartesian vector from rotating-frame coordinates to fixed-frame coordinates.

Convention:
- R maps rotating-frame components -> fixed-frame components.
- R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
- v_fixed = R v_rot, v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
- DeltaR_dst_from_src = R_dst^T R_src.
- C layout statement for helper calls: R[i][j] is row i, column j.
- This routine applies x_fixed = R * x_rot directly.

@param[in] commondata Commondata structure with cumulative regrid basis vectors.
@param[in,out] xCart Cartesian point/vector to be unrotated in place.
"""
    cfunc_type = "void"
    name = "unrotate_xCart_to_fixed_frame"
    params = "const commondata_struct commondata, REAL xCart[3]"
    body = r"""
  REAL xhatU[3] = {commondata.cumulative_regrid_xhatU[0], commondata.cumulative_regrid_xhatU[1], commondata.cumulative_regrid_xhatU[2]};
  REAL yhatU[3] = {commondata.cumulative_regrid_yhatU[0], commondata.cumulative_regrid_yhatU[1], commondata.cumulative_regrid_yhatU[2]};
  REAL zhatU[3] = {commondata.cumulative_regrid_zhatU[0], commondata.cumulative_regrid_zhatU[1], commondata.cumulative_regrid_zhatU[2]};
  if (so3_validate_and_optionally_fix_hats(xhatU, yhatU, zhatU, 0) != 0) {
    fprintf(stderr,
            "ERROR in %s: cumulative hats are invalid at the consuming boundary. "
            "Repair must be applied at the producing/update boundary.\n",
            __func__);
    exit(1);
  }

  REAL R[3][3];
  build_R_from_cumulative_hats(xhatU, yhatU, zhatU, R);
  so3_apply_R_to_vector(R, xCart, xCart);
"""

    cfc.register_CFunction(
        subdirectory="rotation",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=False,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
