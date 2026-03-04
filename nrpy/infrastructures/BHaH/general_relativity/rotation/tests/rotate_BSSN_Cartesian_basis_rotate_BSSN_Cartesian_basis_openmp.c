#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * @brief Rotate BSSN Cartesian-basis fields from axis-angle input.
 *
 * This routine converts (nU,dphi) to DeltaR_dst_from_src and then calls
 * rotate_BSSN_Cartesian_basis_by_R().
 *
 * Convention:
 * - R maps rotating-frame components -> fixed-frame components.
 * - R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
 * - v_fixed = R v_rot, v_rot = R^T v_fixed.
 * - T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
 * - DeltaR_dst_from_src = R_dst^T R_src.
 * - C layout statement for helper calls: R[i][j] is row i, column j.
 *
 * @param[in,out] vetU BSSN rescaled shift vector, rotated in place.
 * @param[in,out] betU BSSN rescaled driver vector, rotated in place.
 * @param[in,out] lambdaU BSSN conformal connection vector, rotated in place.
 * @param[in,out] hDD BSSN conformal metric perturbation tensor, rotated in place.
 * @param[in,out] aDD BSSN trace-free extrinsic curvature tensor, rotated in place.
 * @param[in] nU Rotation axis (need not be normalized).
 * @param[in] dphi Rotation angle in radians.
 */
void rotate_BSSN_Cartesian_basis(REAL vetU[3], REAL betU[3], REAL lambdaU[3], REAL hDD[3][3], REAL aDD[3][3], const REAL nU[3], const REAL dphi) {
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
    const REAL c = cos(dphi);
    const REAL s = sin(dphi);
    const REAL one_minus_c = 1.0 - c;

    // Rodrigues formula: R = c I + (1-c) n n^T + s [n]_x.
    R[0][0] = c + one_minus_c * nx * nx;
    R[0][1] = one_minus_c * nx * ny - s * nz;
    R[0][2] = one_minus_c * nx * nz + s * ny;
    R[1][0] = one_minus_c * ny * nx + s * nz;
    R[1][1] = c + one_minus_c * ny * ny;
    R[1][2] = one_minus_c * ny * nz - s * nx;
    R[2][0] = one_minus_c * nz * nx - s * ny;
    R[2][1] = one_minus_c * nz * ny + s * nx;
    R[2][2] = c + one_minus_c * nz * nz;
  }

  rotate_BSSN_Cartesian_basis_by_R(vetU, betU, lambdaU, hDD, aDD, R);
} // END FUNCTION rotate_BSSN_Cartesian_basis
