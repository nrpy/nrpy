#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * @brief Transform one BSSN point from Cartesian basis to a destination rfm basis.
 *
 * This private kernel assumes one destination coordinate system. The public
 * unsuffixed runtime-dispatch wrapper is emitted separately by
 * ``rfm_wrapper_functions``.
 *
 * @param[in] params Destination-grid parameters used for indexing and coordinate-dependent expressions.
 * @param[in] idx3 Destination-grid flattened point index.
 * @param[in] xx_dst Destination-grid local coordinates for this point.
 * @param[in,out] y_n_gfs Gridfunction storage updated in place.
 */
void basis_transform_BSSN_Cartesian_to_rfm_single_point(const params_struct *restrict params, const int idx3, const REAL xx_dst[3],
                                                        REAL *restrict y_n_gfs) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    basis_transform_BSSN_Cartesian_to_rfm_single_point__rfm__SinhSpherical(params, idx3, xx_dst, y_n_gfs);
    break;
  case SINHSYMTP:
    basis_transform_BSSN_Cartesian_to_rfm_single_point__rfm__SinhSymTP(params, idx3, xx_dst, y_n_gfs);
    break;
  default:
    fprintf(stderr, "ERROR in basis_transform_BSSN_Cartesian_to_rfm_single_point(): CoordSystem hash = %d not #define'd!\n",
            params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION basis_transform_BSSN_Cartesian_to_rfm_single_point
