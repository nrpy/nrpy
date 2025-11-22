#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * @file sqrt_detgammahat_d3xx_volume_element.c
 * @brief Compute the local 3D volume element sqrt(detgammahat) * d3xx at a point.
 *
 * This routine evaluates the positive volume element used when integrating scalar fields on a
 * 3D grid. The evaluation uses reference-metric data and grid spacings stored in the params
 * structure. The computed value is returned by reference via the dV pointer.
 *
 * The expression implemented is:
 *   dV = sqrt(detgammahat(xx0, xx1, xx2)) * abs(dxx0 * dxx1 * dxx2)
 * The absolute value ensures a positive volume element regardless of coordinate orientation.
 *
 * @param[in]  params  Pointer to parameter struct (reference-metric data, grid spacings, and sizes).
 * @param[in]  xx0     Local coordinate 0 at which to evaluate the volume element.
 * @param[in]  xx1     Local coordinate 1 at which to evaluate the volume element.
 * @param[in]  xx2     Local coordinate 2 at which to evaluate the volume element.
 * @param[out] dV      Pointer to the location where the volume element will be stored.
 *
 * @return     void. The result is written to *dV.
 *
 * Note: If a user-editable block is provided in the implementation, users may add custom logic,
 * such as scaling or additional diagnostics, prior to writing the result.
 *
 */
void sqrt_detgammahat_d3xx_volume_element(const params_struct *restrict params, const REAL xx0, const REAL xx1, const REAL xx2, REAL *restrict dV) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    sqrt_detgammahat_d3xx_volume_element__rfm__SinhSpherical(params, xx0, xx1, xx2, dV);
    break;
  default:
    fprintf(stderr, "ERROR in sqrt_detgammahat_d3xx_volume_element(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION sqrt_detgammahat_d3xx_volume_element
