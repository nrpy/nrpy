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
void sqrt_detgammahat_d3xx_volume_element__rfm__SinhSpherical(const params_struct *restrict params, const REAL xx0, const REAL xx1, const REAL xx2,
                                                              REAL *restrict dV) {
  const REAL AMPL = params->AMPL;
  const REAL SINHW = params->SINHW;
  const REAL dxx0 = params->dxx0;
  const REAL dxx1 = params->dxx1;
  const REAL dxx2 = params->dxx2;
  /*
   *  Original SymPy expression:
   *  "*dV = AMPL**2*sqrt((1 - exp(2*xx0/SINHW))**4*(exp(2*xx0/SINHW) + 1)**2*exp(6*(1 - xx0)/SINHW)*sin(xx1)**2/(SINHW**2*(1 -
   * exp(2/SINHW))**6))*Abs(AMPL)*Abs(dxx0*dxx1*dxx2)"
   */
  const REAL tmp0 = (1.0 / (SINHW));
  const REAL tmp2 = exp(2 * tmp0 * xx0);
  *dV = ((AMPL) * (AMPL)) *
        sqrt(((1 - tmp2) * (1 - tmp2) * (1 - tmp2) * (1 - tmp2)) * ((tmp2 + 1) * (tmp2 + 1)) * exp(6 * tmp0 * (1 - xx0)) * ((sin(xx1)) * (sin(xx1))) /
             (((SINHW) * (SINHW)) * pow(1 - exp(2 * tmp0), 6))) *
        fabs(AMPL) * fabs(dxx0 * dxx1 * dxx2);
} // END FUNCTION sqrt_detgammahat_d3xx_volume_element__rfm__SinhSpherical
