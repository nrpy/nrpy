"""
Register a helper that returns two axis-angle output slots.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_function as cfc
from nrpy.infrastructures.BHaH.rotation.so3_matrix_ops import (
    register_CFunction_build_R_from_cumulative_hats,
    register_CFunction_so3_matrix_to_axis_angle,
    register_CFunction_so3_validate_and_optionally_fix_hats,
    register_rotation_commondata_CodeParameters,
)


def register_CFunction_unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame() -> (
    None
):
    r"""
    Register helper for two axis-angle output slots.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> from nrpy.infrastructures.BHaH.rotation.so3_matrix_ops import assert_so3_convention_in_text
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame()
    >>> generated_str = cfc.CFunction_dict["unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame"].full_function
    >>> assert_so3_convention_in_text(generated_str, "unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame")
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "openmp", file_ext="c")
    """
    register_rotation_commondata_CodeParameters()
    if "so3_validate_and_optionally_fix_hats" not in cfc.CFunction_dict:
        register_CFunction_so3_validate_and_optionally_fix_hats()
    if "build_R_from_cumulative_hats" not in cfc.CFunction_dict:
        register_CFunction_build_R_from_cumulative_hats()
    if "so3_matrix_to_axis_angle" not in cfc.CFunction_dict:
        register_CFunction_so3_matrix_to_axis_angle()

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r"""
@brief Compute axis-angle outputs from cumulative hats.

This routine computes a single unrotation from the cumulative basis matrix and converts it to
one axis-angle pair in (@p nU_part1, @p dphi_part1). The second pair is set to
identity: (@p nU_part2=(1,0,0), @p dphi_part2=0).

Convention:
- R maps rotating-frame components -> fixed-frame components.
- R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
- v_fixed = R v_rot, v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
- DeltaR_dst_from_src = R_dst^T R_src.
- C layout statement for helper calls: R[i][j] is row i, column j.
- Axis-angle bridge is produced from this R using so3_matrix_to_axis_angle().

@param[in] commondata Commondata structure containing cumulative regrid basis vectors.
@param[out] nU_part1 Rotation axis for first rotation.
@param[out] dphi_part1 Rotation angle for first rotation.
@param[out] nU_part2 Rotation axis for second rotation.
@param[out] dphi_part2 Rotation angle for second rotation.
"""
    cfunc_type = "void"
    name = "unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame"
    params = (
        "const commondata_struct commondata, REAL nU_part1[3], "
        "REAL *restrict dphi_part1, REAL nU_part2[3], REAL *restrict dphi_part2"
    )
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

  so3_matrix_to_axis_angle(R, nU_part1, dphi_part1);

  // Second output slot is deterministic identity.
  nU_part2[0] = 1.0;
  nU_part2[1] = 0.0;
  nU_part2[2] = 0.0;
  *dphi_part2 = 0.0;
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
