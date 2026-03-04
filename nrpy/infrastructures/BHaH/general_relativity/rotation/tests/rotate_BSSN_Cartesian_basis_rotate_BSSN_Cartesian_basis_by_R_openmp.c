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
 *   components before applying equation-derived SO(3) tensor-rotation expressions.
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

  // Alias-safe vector rotations: v_dst = DeltaR * v_src.
  const REAL vetU_in[3] = {vetU[0], vetU[1], vetU[2]};
  const REAL betU_in[3] = {betU[0], betU[1], betU[2]};
  const REAL lambdaU_in[3] = {lambdaU[0], lambdaU[1], lambdaU[2]};
  const REAL tmp0 = DeltaR_dst_from_src[0][0] * DeltaR_dst_from_src[0][1];
  const REAL tmp1 = DeltaR_dst_from_src[0][0] * DeltaR_dst_from_src[0][2];
  const REAL tmp2 = DeltaR_dst_from_src[0][1] * DeltaR_dst_from_src[0][2];
  const REAL tmp17 = DeltaR_dst_from_src[1][1] * hDD_sym[1][0];
  const REAL tmp18 = DeltaR_dst_from_src[1][2] * hDD_sym[2][0];
  const REAL tmp19 = DeltaR_dst_from_src[1][0] * hDD_sym[0][1];
  const REAL tmp20 = DeltaR_dst_from_src[1][2] * hDD_sym[2][1];
  const REAL tmp21 = DeltaR_dst_from_src[1][0] * hDD_sym[0][2];
  const REAL tmp22 = DeltaR_dst_from_src[1][1] * hDD_sym[1][2];
  const REAL tmp30 = DeltaR_dst_from_src[2][1] * hDD_sym[1][0];
  const REAL tmp31 = DeltaR_dst_from_src[2][2] * hDD_sym[2][0];
  const REAL tmp32 = DeltaR_dst_from_src[2][0] * hDD_sym[0][1];
  const REAL tmp33 = DeltaR_dst_from_src[2][2] * hDD_sym[2][1];
  const REAL tmp34 = DeltaR_dst_from_src[2][0] * hDD_sym[0][2];
  const REAL tmp35 = DeltaR_dst_from_src[2][1] * hDD_sym[1][2];
  const REAL tmp50 = DeltaR_dst_from_src[1][1] * aDD_sym[1][0];
  const REAL tmp51 = DeltaR_dst_from_src[1][2] * aDD_sym[2][0];
  const REAL tmp52 = DeltaR_dst_from_src[1][0] * aDD_sym[0][1];
  const REAL tmp53 = DeltaR_dst_from_src[1][2] * aDD_sym[2][1];
  const REAL tmp54 = DeltaR_dst_from_src[1][0] * aDD_sym[0][2];
  const REAL tmp55 = DeltaR_dst_from_src[1][1] * aDD_sym[1][2];
  const REAL tmp57 = DeltaR_dst_from_src[2][1] * aDD_sym[1][0];
  const REAL tmp58 = DeltaR_dst_from_src[2][2] * aDD_sym[2][0];
  const REAL tmp59 = DeltaR_dst_from_src[2][0] * aDD_sym[0][1];
  const REAL tmp60 = DeltaR_dst_from_src[2][2] * aDD_sym[2][1];
  const REAL tmp61 = DeltaR_dst_from_src[2][0] * aDD_sym[0][2];
  const REAL tmp62 = DeltaR_dst_from_src[2][1] * aDD_sym[1][2];
  const REAL tmp15 = DeltaR_dst_from_src[0][0] * DeltaR_dst_from_src[1][0] * hDD_sym[0][0] +
                     DeltaR_dst_from_src[0][1] * DeltaR_dst_from_src[1][1] * hDD_sym[1][1] +
                     DeltaR_dst_from_src[0][2] * DeltaR_dst_from_src[1][2] * hDD_sym[2][2];
  const REAL tmp16 = DeltaR_dst_from_src[0][0] * DeltaR_dst_from_src[2][0] * hDD_sym[0][0] +
                     DeltaR_dst_from_src[0][1] * DeltaR_dst_from_src[2][1] * hDD_sym[1][1] +
                     DeltaR_dst_from_src[0][2] * DeltaR_dst_from_src[2][2] * hDD_sym[2][2];
  const REAL tmp29 = DeltaR_dst_from_src[1][0] * DeltaR_dst_from_src[2][0] * hDD_sym[0][0] +
                     DeltaR_dst_from_src[1][1] * DeltaR_dst_from_src[2][1] * hDD_sym[1][1] +
                     DeltaR_dst_from_src[1][2] * DeltaR_dst_from_src[2][2] * hDD_sym[2][2];
  const REAL tmp48 = DeltaR_dst_from_src[0][0] * DeltaR_dst_from_src[1][0] * aDD_sym[0][0] +
                     DeltaR_dst_from_src[0][1] * DeltaR_dst_from_src[1][1] * aDD_sym[1][1] +
                     DeltaR_dst_from_src[0][2] * DeltaR_dst_from_src[1][2] * aDD_sym[2][2];
  const REAL tmp49 = DeltaR_dst_from_src[0][0] * DeltaR_dst_from_src[2][0] * aDD_sym[0][0] +
                     DeltaR_dst_from_src[0][1] * DeltaR_dst_from_src[2][1] * aDD_sym[1][1] +
                     DeltaR_dst_from_src[0][2] * DeltaR_dst_from_src[2][2] * aDD_sym[2][2];
  const REAL tmp56 = DeltaR_dst_from_src[1][0] * DeltaR_dst_from_src[2][0] * aDD_sym[0][0] +
                     DeltaR_dst_from_src[1][1] * DeltaR_dst_from_src[2][1] * aDD_sym[1][1] +
                     DeltaR_dst_from_src[1][2] * DeltaR_dst_from_src[2][2] * aDD_sym[2][2];
  vetU_out[0] = DeltaR_dst_from_src[0][0] * vetU_in[0] + DeltaR_dst_from_src[0][1] * vetU_in[1] + DeltaR_dst_from_src[0][2] * vetU_in[2];
  vetU_out[1] = DeltaR_dst_from_src[1][0] * vetU_in[0] + DeltaR_dst_from_src[1][1] * vetU_in[1] + DeltaR_dst_from_src[1][2] * vetU_in[2];
  vetU_out[2] = DeltaR_dst_from_src[2][0] * vetU_in[0] + DeltaR_dst_from_src[2][1] * vetU_in[1] + DeltaR_dst_from_src[2][2] * vetU_in[2];
  betU_out[0] = DeltaR_dst_from_src[0][0] * betU_in[0] + DeltaR_dst_from_src[0][1] * betU_in[1] + DeltaR_dst_from_src[0][2] * betU_in[2];
  betU_out[1] = DeltaR_dst_from_src[1][0] * betU_in[0] + DeltaR_dst_from_src[1][1] * betU_in[1] + DeltaR_dst_from_src[1][2] * betU_in[2];
  betU_out[2] = DeltaR_dst_from_src[2][0] * betU_in[0] + DeltaR_dst_from_src[2][1] * betU_in[1] + DeltaR_dst_from_src[2][2] * betU_in[2];
  lambdaU_out[0] = DeltaR_dst_from_src[0][0] * lambdaU_in[0] + DeltaR_dst_from_src[0][1] * lambdaU_in[1] + DeltaR_dst_from_src[0][2] * lambdaU_in[2];
  lambdaU_out[1] = DeltaR_dst_from_src[1][0] * lambdaU_in[0] + DeltaR_dst_from_src[1][1] * lambdaU_in[1] + DeltaR_dst_from_src[1][2] * lambdaU_in[2];
  lambdaU_out[2] = DeltaR_dst_from_src[2][0] * lambdaU_in[0] + DeltaR_dst_from_src[2][1] * lambdaU_in[1] + DeltaR_dst_from_src[2][2] * lambdaU_in[2];
  hDD_out[0][0] = ((DeltaR_dst_from_src[0][0]) * (DeltaR_dst_from_src[0][0])) * hDD_sym[0][0] +
                  ((DeltaR_dst_from_src[0][1]) * (DeltaR_dst_from_src[0][1])) * hDD_sym[1][1] +
                  ((DeltaR_dst_from_src[0][2]) * (DeltaR_dst_from_src[0][2])) * hDD_sym[2][2] + hDD_sym[0][1] * tmp0 + hDD_sym[0][2] * tmp1 +
                  hDD_sym[1][0] * tmp0 + hDD_sym[1][2] * tmp2 + hDD_sym[2][0] * tmp1 + hDD_sym[2][1] * tmp2;
  hDD_out[0][1] =
      DeltaR_dst_from_src[0][0] * DeltaR_dst_from_src[1][1] * hDD_sym[0][1] + DeltaR_dst_from_src[0][0] * DeltaR_dst_from_src[1][2] * hDD_sym[0][2] +
      DeltaR_dst_from_src[0][1] * DeltaR_dst_from_src[1][0] * hDD_sym[1][0] + DeltaR_dst_from_src[0][1] * DeltaR_dst_from_src[1][2] * hDD_sym[1][2] +
      DeltaR_dst_from_src[0][2] * DeltaR_dst_from_src[1][0] * hDD_sym[2][0] + DeltaR_dst_from_src[0][2] * DeltaR_dst_from_src[1][1] * hDD_sym[2][1] +
      tmp15;
  hDD_out[0][2] =
      DeltaR_dst_from_src[0][0] * DeltaR_dst_from_src[2][1] * hDD_sym[0][1] + DeltaR_dst_from_src[0][0] * DeltaR_dst_from_src[2][2] * hDD_sym[0][2] +
      DeltaR_dst_from_src[0][1] * DeltaR_dst_from_src[2][0] * hDD_sym[1][0] + DeltaR_dst_from_src[0][1] * DeltaR_dst_from_src[2][2] * hDD_sym[1][2] +
      DeltaR_dst_from_src[0][2] * DeltaR_dst_from_src[2][0] * hDD_sym[2][0] + DeltaR_dst_from_src[0][2] * DeltaR_dst_from_src[2][1] * hDD_sym[2][1] +
      tmp16;
  hDD_out[1][0] = DeltaR_dst_from_src[0][0] * tmp17 + DeltaR_dst_from_src[0][0] * tmp18 + DeltaR_dst_from_src[0][1] * tmp19 +
                  DeltaR_dst_from_src[0][1] * tmp20 + DeltaR_dst_from_src[0][2] * tmp21 + DeltaR_dst_from_src[0][2] * tmp22 + tmp15;
  hDD_out[1][1] = ((DeltaR_dst_from_src[1][0]) * (DeltaR_dst_from_src[1][0])) * hDD_sym[0][0] + DeltaR_dst_from_src[1][0] * tmp17 +
                  DeltaR_dst_from_src[1][0] * tmp18 + ((DeltaR_dst_from_src[1][1]) * (DeltaR_dst_from_src[1][1])) * hDD_sym[1][1] +
                  DeltaR_dst_from_src[1][1] * tmp19 + DeltaR_dst_from_src[1][1] * tmp20 +
                  ((DeltaR_dst_from_src[1][2]) * (DeltaR_dst_from_src[1][2])) * hDD_sym[2][2] + DeltaR_dst_from_src[1][2] * tmp21 +
                  DeltaR_dst_from_src[1][2] * tmp22;
  hDD_out[1][2] = DeltaR_dst_from_src[2][0] * tmp17 + DeltaR_dst_from_src[2][0] * tmp18 + DeltaR_dst_from_src[2][1] * tmp19 +
                  DeltaR_dst_from_src[2][1] * tmp20 + DeltaR_dst_from_src[2][2] * tmp21 + DeltaR_dst_from_src[2][2] * tmp22 + tmp29;
  hDD_out[2][0] = DeltaR_dst_from_src[0][0] * tmp30 + DeltaR_dst_from_src[0][0] * tmp31 + DeltaR_dst_from_src[0][1] * tmp32 +
                  DeltaR_dst_from_src[0][1] * tmp33 + DeltaR_dst_from_src[0][2] * tmp34 + DeltaR_dst_from_src[0][2] * tmp35 + tmp16;
  hDD_out[2][1] = DeltaR_dst_from_src[1][0] * tmp30 + DeltaR_dst_from_src[1][0] * tmp31 + DeltaR_dst_from_src[1][1] * tmp32 +
                  DeltaR_dst_from_src[1][1] * tmp33 + DeltaR_dst_from_src[1][2] * tmp34 + DeltaR_dst_from_src[1][2] * tmp35 + tmp29;
  hDD_out[2][2] = ((DeltaR_dst_from_src[2][0]) * (DeltaR_dst_from_src[2][0])) * hDD_sym[0][0] + DeltaR_dst_from_src[2][0] * tmp30 +
                  DeltaR_dst_from_src[2][0] * tmp31 + ((DeltaR_dst_from_src[2][1]) * (DeltaR_dst_from_src[2][1])) * hDD_sym[1][1] +
                  DeltaR_dst_from_src[2][1] * tmp32 + DeltaR_dst_from_src[2][1] * tmp33 +
                  ((DeltaR_dst_from_src[2][2]) * (DeltaR_dst_from_src[2][2])) * hDD_sym[2][2] + DeltaR_dst_from_src[2][2] * tmp34 +
                  DeltaR_dst_from_src[2][2] * tmp35;
  aDD_out[0][0] = ((DeltaR_dst_from_src[0][0]) * (DeltaR_dst_from_src[0][0])) * aDD_sym[0][0] +
                  ((DeltaR_dst_from_src[0][1]) * (DeltaR_dst_from_src[0][1])) * aDD_sym[1][1] +
                  ((DeltaR_dst_from_src[0][2]) * (DeltaR_dst_from_src[0][2])) * aDD_sym[2][2] + aDD_sym[0][1] * tmp0 + aDD_sym[0][2] * tmp1 +
                  aDD_sym[1][0] * tmp0 + aDD_sym[1][2] * tmp2 + aDD_sym[2][0] * tmp1 + aDD_sym[2][1] * tmp2;
  aDD_out[0][1] =
      DeltaR_dst_from_src[0][0] * DeltaR_dst_from_src[1][1] * aDD_sym[0][1] + DeltaR_dst_from_src[0][0] * DeltaR_dst_from_src[1][2] * aDD_sym[0][2] +
      DeltaR_dst_from_src[0][1] * DeltaR_dst_from_src[1][0] * aDD_sym[1][0] + DeltaR_dst_from_src[0][1] * DeltaR_dst_from_src[1][2] * aDD_sym[1][2] +
      DeltaR_dst_from_src[0][2] * DeltaR_dst_from_src[1][0] * aDD_sym[2][0] + DeltaR_dst_from_src[0][2] * DeltaR_dst_from_src[1][1] * aDD_sym[2][1] +
      tmp48;
  aDD_out[0][2] =
      DeltaR_dst_from_src[0][0] * DeltaR_dst_from_src[2][1] * aDD_sym[0][1] + DeltaR_dst_from_src[0][0] * DeltaR_dst_from_src[2][2] * aDD_sym[0][2] +
      DeltaR_dst_from_src[0][1] * DeltaR_dst_from_src[2][0] * aDD_sym[1][0] + DeltaR_dst_from_src[0][1] * DeltaR_dst_from_src[2][2] * aDD_sym[1][2] +
      DeltaR_dst_from_src[0][2] * DeltaR_dst_from_src[2][0] * aDD_sym[2][0] + DeltaR_dst_from_src[0][2] * DeltaR_dst_from_src[2][1] * aDD_sym[2][1] +
      tmp49;
  aDD_out[1][0] = DeltaR_dst_from_src[0][0] * tmp50 + DeltaR_dst_from_src[0][0] * tmp51 + DeltaR_dst_from_src[0][1] * tmp52 +
                  DeltaR_dst_from_src[0][1] * tmp53 + DeltaR_dst_from_src[0][2] * tmp54 + DeltaR_dst_from_src[0][2] * tmp55 + tmp48;
  aDD_out[1][1] = ((DeltaR_dst_from_src[1][0]) * (DeltaR_dst_from_src[1][0])) * aDD_sym[0][0] + DeltaR_dst_from_src[1][0] * tmp50 +
                  DeltaR_dst_from_src[1][0] * tmp51 + ((DeltaR_dst_from_src[1][1]) * (DeltaR_dst_from_src[1][1])) * aDD_sym[1][1] +
                  DeltaR_dst_from_src[1][1] * tmp52 + DeltaR_dst_from_src[1][1] * tmp53 +
                  ((DeltaR_dst_from_src[1][2]) * (DeltaR_dst_from_src[1][2])) * aDD_sym[2][2] + DeltaR_dst_from_src[1][2] * tmp54 +
                  DeltaR_dst_from_src[1][2] * tmp55;
  aDD_out[1][2] = DeltaR_dst_from_src[2][0] * tmp50 + DeltaR_dst_from_src[2][0] * tmp51 + DeltaR_dst_from_src[2][1] * tmp52 +
                  DeltaR_dst_from_src[2][1] * tmp53 + DeltaR_dst_from_src[2][2] * tmp54 + DeltaR_dst_from_src[2][2] * tmp55 + tmp56;
  aDD_out[2][0] = DeltaR_dst_from_src[0][0] * tmp57 + DeltaR_dst_from_src[0][0] * tmp58 + DeltaR_dst_from_src[0][1] * tmp59 +
                  DeltaR_dst_from_src[0][1] * tmp60 + DeltaR_dst_from_src[0][2] * tmp61 + DeltaR_dst_from_src[0][2] * tmp62 + tmp49;
  aDD_out[2][1] = DeltaR_dst_from_src[1][0] * tmp57 + DeltaR_dst_from_src[1][0] * tmp58 + DeltaR_dst_from_src[1][1] * tmp59 +
                  DeltaR_dst_from_src[1][1] * tmp60 + DeltaR_dst_from_src[1][2] * tmp61 + DeltaR_dst_from_src[1][2] * tmp62 + tmp56;
  aDD_out[2][2] = ((DeltaR_dst_from_src[2][0]) * (DeltaR_dst_from_src[2][0])) * aDD_sym[0][0] + DeltaR_dst_from_src[2][0] * tmp57 +
                  DeltaR_dst_from_src[2][0] * tmp58 + ((DeltaR_dst_from_src[2][1]) * (DeltaR_dst_from_src[2][1])) * aDD_sym[1][1] +
                  DeltaR_dst_from_src[2][1] * tmp59 + DeltaR_dst_from_src[2][1] * tmp60 +
                  ((DeltaR_dst_from_src[2][2]) * (DeltaR_dst_from_src[2][2])) * aDD_sym[2][2] + DeltaR_dst_from_src[2][2] * tmp61 +
                  DeltaR_dst_from_src[2][2] * tmp62;

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
