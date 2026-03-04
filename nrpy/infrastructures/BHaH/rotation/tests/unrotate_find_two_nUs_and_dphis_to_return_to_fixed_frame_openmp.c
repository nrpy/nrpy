#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * @brief Compute axis-angle outputs from cumulative hats.
 *
 * This routine computes a single unrotation from the cumulative basis matrix and converts it to
 * one axis-angle pair in (@p nU_part1, @p dphi_part1). The second pair is set to
 * identity: (@p nU_part2=(1,0,0), @p dphi_part2=0).
 *
 * Convention:
 * - R maps rotating-frame components -> fixed-frame components.
 * - R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
 * - v_fixed = R v_rot, v_rot = R^T v_fixed.
 * - T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
 * - DeltaR_dst_from_src = R_dst^T R_src.
 * - C layout statement for helper calls: R[i][j] is row i, column j.
 * - Axis-angle bridge is produced from this R using SO3_matrix_to_axis_angle().
 * - Einstein notation for the matrix map represented by part 1:
 *   v^i_fixed = R^i{}_j v^j_rot.
 *
 * @param[in] commondata Commondata structure containing cumulative regrid basis vectors.
 * @param[out] nU_part1 Rotation axis for first rotation.
 * @param[out] dphi_part1 Rotation angle for first rotation.
 * @param[out] nU_part2 Rotation axis for second rotation.
 * @param[out] dphi_part2 Rotation angle for second rotation.
 */
void unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame(const commondata_struct commondata, REAL nU_part1[3], REAL *restrict dphi_part1,
                                                              REAL nU_part2[3], REAL *restrict dphi_part2) {
  REAL xhatU[3] = {commondata.cumulative_regrid_xhatU[0], commondata.cumulative_regrid_xhatU[1], commondata.cumulative_regrid_xhatU[2]};
  REAL yhatU[3] = {commondata.cumulative_regrid_yhatU[0], commondata.cumulative_regrid_yhatU[1], commondata.cumulative_regrid_yhatU[2]};
  REAL zhatU[3] = {commondata.cumulative_regrid_zhatU[0], commondata.cumulative_regrid_zhatU[1], commondata.cumulative_regrid_zhatU[2]};
  if (SO3_validate_and_optionally_fix_hats(xhatU, yhatU, zhatU, 0) != 0) {
    fprintf(stderr,
            "ERROR in %s: cumulative hats are invalid at the consuming boundary. "
            "Repair must be applied at the producing/update boundary.\n",
            __func__);
    exit(1);
  }

  REAL R[3][3];
  build_R_from_cumulative_hats(xhatU, yhatU, zhatU, R);

  SO3_matrix_to_axis_angle(R, nU_part1, dphi_part1);

  // Second output slot is deterministic identity.
  nU_part2[0] = 1.0;
  nU_part2[1] = 0.0;
  nU_part2[2] = 0.0;
  *dphi_part2 = 0.0;
} // END FUNCTION unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame
