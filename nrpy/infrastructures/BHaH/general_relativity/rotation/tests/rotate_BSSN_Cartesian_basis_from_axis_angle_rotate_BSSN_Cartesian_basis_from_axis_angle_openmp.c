#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * @brief Rotate BSSN Cartesian-basis fields from axis-angle input.
 *
 * This routine converts (nU,dphi) to DeltaR_dst_from_src and then calls
 * rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src().
 *
 * Convention:
 * - R maps rotating-frame components -> fixed-frame components.
 * - R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
 * - v_fixed = R v_rot, v_rot = R^T v_fixed.
 * - T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
 * - DeltaR_dst_from_src = R_dst^T R_src.
 * - C layout statement for helper calls: R[i][j] is row i, column j.
 * - Einstein notation for the callee update:
 *   v^i_dst = (DeltaR)^i{}_j v^j_src,
 *   T^dst_{ij} = (DeltaR)_i{}^k (DeltaR)_j{}^l T^src_{kl}.
 * - Symmetric tensor storage contract inherited from
 *   rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src(): only upper-triangular components
 *   (i <= j) are updated.
 *
 * @param[in,out] vetU BSSN rescaled shift vector, rotated in place.
 * @param[in,out] betU BSSN rescaled driver vector, rotated in place.
 * @param[in,out] lambdaU BSSN conformal connection vector, rotated in place.
 * @param[in,out] hDD BSSN conformal metric perturbation tensor, rotated in place.
 * @param[in,out] aDD BSSN trace-free extrinsic curvature tensor, rotated in place.
 * @param[in] nU Rotation axis (need not be normalized).
 * @param[in] dphi Rotation angle in radians.
 */
void rotate_BSSN_Cartesian_basis_from_axis_angle(REAL vetU[3], REAL betU[3], REAL lambdaU[3], REAL hDD[3][3], REAL aDD[3][3], const REAL nU[3],
                                                 const REAL dphi) {
  const REAL n0 = nU[0];
  const REAL n1 = nU[1];
  const REAL n2 = nU[2];
  const REAL nnorm = sqrt(n0 * n0 + n1 * n1 + n2 * n2);

  REAL R[3][3];
  if (nnorm < 1e-300) {
    // Deterministic fallback for degenerate axis input.
    R[0][0] = 1.0;
    R[0][1] = 0.0;
    R[0][2] = 0.0;
    R[1][0] = 0.0;
    R[1][1] = 1.0;
    R[1][2] = 0.0;
    R[2][0] = 0.0;
    R[2][1] = 0.0;
    R[2][2] = 1.0;
  } else {
    const REAL nx = n0 / nnorm;
    const REAL ny = n1 / nnorm;
    const REAL nz = n2 / nnorm;
    const REAL tmp0 = cos(dphi);
    const REAL tmp2 = sin(dphi);
    const REAL tmp1 = 1 - tmp0;
    const REAL tmp5 = nx * ny * tmp1;
    const REAL tmp7 = nx * nz * tmp1;
    R[0][0] = ((nx) * (nx)) * tmp1 + tmp0;
    R[0][1] = -nz * tmp2 + tmp5;
    R[0][2] = ny * tmp2 + tmp7;
    R[1][0] = nz * tmp2 + tmp5;
    R[1][1] = ((ny) * (ny)) * tmp1 + tmp0;
    R[1][2] = -nx * tmp2 + ny * nz * tmp1;
    R[2][0] = -ny * tmp2 + tmp7;
    R[2][1] = nx * tmp2 + ny * nz * tmp1;
    R[2][2] = ((nz) * (nz)) * tmp1 + tmp0;
  }

  rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src(vetU, betU, lambdaU, hDD, aDD, R);
} // END FUNCTION rotate_BSSN_Cartesian_basis_from_axis_angle
