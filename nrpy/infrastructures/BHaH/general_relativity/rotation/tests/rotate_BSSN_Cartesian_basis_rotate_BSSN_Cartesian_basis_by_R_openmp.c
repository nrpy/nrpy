#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * @brief Rotate Cartesian-basis BSSN vectors and symmetric tensors using a matrix.
 *
 * Convention:
 * - R maps rotating-frame components -> fixed-frame components.
 * - R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
 * - v_fixed = R v_rot, v_rot = R^T v_fixed.
 * - T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
 * - DeltaR_dst_from_src = R_dst^T R_src.
 * - C layout statement for helper calls: R[i][j] is row i, column j.
 * - DeltaR_dst_from_src maps source rotating-basis components to destination
 *   rotating-basis components.
 * - Vector update in Einstein notation:
 *   v^i_dst = (DeltaR)^i{}_j v^j_src.
 * - Covariant rank-2 update in Einstein notation:
 *   T^dst_{ij} = (DeltaR)_i{}^k (DeltaR)_j{}^l T^src_{kl}.
 * - Symmetric tensor storage contract:
 *   only upper-triangular components (i <= j) are authoritative on input
 *   and are overwritten on output; lower-triangular entries are untouched.
 * - Internally, full symmetric local tensors are reconstructed from upper
 *   components before calling SO3_apply_R_to_tensorDD().
 *
 * @param[in,out] vetU BSSN rescaled shift vector, rotated in place.
 * @param[in,out] betU BSSN rescaled driver vector, rotated in place.
 * @param[in,out] lambdaU BSSN conformal connection vector, rotated in place.
 * @param[in,out] hDD BSSN conformal metric perturbation tensor, rotated in place.
 * @param[in,out] aDD BSSN trace-free extrinsic curvature tensor, rotated in place.
 * @param[in] DeltaR_dst_from_src Relative basis rotation matrix.
 */
void rotate_BSSN_Cartesian_basis_by_R(REAL vetU[3], REAL betU[3], REAL lambdaU[3], REAL hDD[3][3], REAL aDD[3][3],
                                      const REAL DeltaR_dst_from_src[3][3]) {
  REAL vetU_out[3], betU_out[3], lambdaU_out[3];
  REAL hDD_sym[3][3], aDD_sym[3][3];
  REAL hDD_out[3][3], aDD_out[3][3];

  // Reconstruct full symmetric tensors from authoritative upper-triangular input.
  hDD_sym[0][0] = hDD[0][0];
  hDD_sym[0][1] = hDD[0][1];
  hDD_sym[0][2] = hDD[0][2];
  hDD_sym[1][0] = hDD[0][1];
  hDD_sym[1][1] = hDD[1][1];
  hDD_sym[1][2] = hDD[1][2];
  hDD_sym[2][0] = hDD[0][2];
  hDD_sym[2][1] = hDD[1][2];
  hDD_sym[2][2] = hDD[2][2];

  aDD_sym[0][0] = aDD[0][0];
  aDD_sym[0][1] = aDD[0][1];
  aDD_sym[0][2] = aDD[0][2];
  aDD_sym[1][0] = aDD[0][1];
  aDD_sym[1][1] = aDD[1][1];
  aDD_sym[1][2] = aDD[1][2];
  aDD_sym[2][0] = aDD[0][2];
  aDD_sym[2][1] = aDD[1][2];
  aDD_sym[2][2] = aDD[2][2];

  SO3_apply_R_to_vector(DeltaR_dst_from_src, vetU, vetU_out);
  SO3_apply_R_to_vector(DeltaR_dst_from_src, betU, betU_out);
  SO3_apply_R_to_vector(DeltaR_dst_from_src, lambdaU, lambdaU_out);
  SO3_apply_R_to_tensorDD(DeltaR_dst_from_src, hDD_sym, hDD_out);
  SO3_apply_R_to_tensorDD(DeltaR_dst_from_src, aDD_sym, aDD_out);

  // Write back vectors completely.
  vetU[0] = vetU_out[0];
  vetU[1] = vetU_out[1];
  vetU[2] = vetU_out[2];
  betU[0] = betU_out[0];
  betU[1] = betU_out[1];
  betU[2] = betU_out[2];
  lambdaU[0] = lambdaU_out[0];
  lambdaU[1] = lambdaU_out[1];
  lambdaU[2] = lambdaU_out[2];

  // Write back only upper-triangular tensor components.
  hDD[0][0] = hDD_out[0][0];
  hDD[0][1] = hDD_out[0][1];
  hDD[0][2] = hDD_out[0][2];
  hDD[1][1] = hDD_out[1][1];
  hDD[1][2] = hDD_out[1][2];
  hDD[2][2] = hDD_out[2][2];

  aDD[0][0] = aDD_out[0][0];
  aDD[0][1] = aDD_out[0][1];
  aDD[0][2] = aDD_out[0][2];
  aDD[1][1] = aDD_out[1][1];
  aDD[1][2] = aDD_out[1][2];
  aDD[2][2] = aDD_out[2][2];
} // END FUNCTION rotate_BSSN_Cartesian_basis_by_R
