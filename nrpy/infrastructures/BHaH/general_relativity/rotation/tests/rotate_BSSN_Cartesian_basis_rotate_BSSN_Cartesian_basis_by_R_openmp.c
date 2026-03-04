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
 * - Vectors: v_dst = DeltaR_dst_from_src * v_src.
 * - Rank-2 tensors: T_dst = DeltaR_dst_from_src * T_src *
 *   DeltaR_dst_from_src^T.
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
  // Enforce symmetry by mirroring upper-triangular entries to lower-triangular.
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < i; j++) {
      hDD[i][j] = hDD[j][i];
      aDD[i][j] = aDD[j][i];
    } // END LOOP over strictly-lower-triangular entries.
  } // END LOOP over tensor rows.

  REAL vetU_out[3], betU_out[3], lambdaU_out[3];
  REAL hDD_out[3][3], aDD_out[3][3];
  so3_apply_R_to_vector(DeltaR_dst_from_src, vetU, vetU_out);
  so3_apply_R_to_vector(DeltaR_dst_from_src, betU, betU_out);
  so3_apply_R_to_vector(DeltaR_dst_from_src, lambdaU, lambdaU_out);
  so3_apply_R_to_tensorDD(DeltaR_dst_from_src, hDD, hDD_out);
  so3_apply_R_to_tensorDD(DeltaR_dst_from_src, aDD, aDD_out);

  for (int i = 0; i < 3; i++) {
    vetU[i] = vetU_out[i];
    betU[i] = betU_out[i];
    lambdaU[i] = lambdaU_out[i];
    for (int j = 0; j < 3; j++) {
      hDD[i][j] = hDD_out[i][j];
      aDD[i][j] = aDD_out[i][j];
    }
  }

  // Re-enforce symmetry by averaging mirrored tensor entries.
  for (int i = 0; i < 3; i++) {
    for (int j = i + 1; j < 3; j++) {
      const REAL hsym = 0.5 * (hDD[i][j] + hDD[j][i]);
      const REAL asym = 0.5 * (aDD[i][j] + aDD[j][i]);
      hDD[i][j] = hsym;
      hDD[j][i] = hsym;
      aDD[i][j] = asym;
      aDD[j][i] = asym;
    } // END LOOP over strictly-upper-triangular entries.
  } // END LOOP over tensor rows.
} // END FUNCTION rotate_BSSN_Cartesian_basis_by_R
