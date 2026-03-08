#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * @brief Transform one BSSN point from a source rfm basis to Cartesian basis.
 *
 * This private kernel assumes source-grid coordinate-dependent expressions for a
 * single coordinate system. The public unsuffixed runtime-dispatch wrapper is
 * emitted separately by ``rfm_wrapper_functions``.
 *
 * @param[in] params Source-grid parameters used for coordinate-dependent expressions.
 * @param[in] paramsdst Destination-grid parameters used only for indexing at ``idx3``.
 * @param[in] idx3 Destination-grid flattened point index.
 * @param[in] xx_src Source-grid local coordinates for this point.
 * @param[in,out] y_n_gfs Gridfunction storage updated in place.
 */
void basis_transform_BSSN_rfm_to_Cartesian_single_point(const params_struct *restrict params, const params_struct *restrict paramsdst, const int idx3,
                                                        const REAL xx_src[3], REAL *restrict y_n_gfs) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    basis_transform_BSSN_rfm_to_Cartesian_single_point__rfm__SinhSpherical(params, paramsdst, idx3, xx_src, y_n_gfs);
    break;
  case SINHSYMTP:
    basis_transform_BSSN_rfm_to_Cartesian_single_point__rfm__SinhSymTP(params, paramsdst, idx3, xx_src, y_n_gfs);
    break;
  default:
    fprintf(stderr, "ERROR in basis_transform_BSSN_rfm_to_Cartesian_single_point(): CoordSystem hash = %d not #define'd!\n",
            params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION basis_transform_BSSN_rfm_to_Cartesian_single_point
