#include "BHaH_defines.h"

/**
 * @brief Computes conserved quantities for KerrSchild_Cartesian (massive).
 *
 *         Computed: E, Lx, Ly, Lz, Q.
 *
 *         Expects state vector f[8] containing position and 4-velocity u^mu.
 *
 *         Note: Input vector components are mapped to local variables p0..p3 to match
 *         the symbolic variable names used in the equation generation module.
 */
void conserved_quantities_KerrSchild_Cartesian_massive(const commondata_struct *restrict commondata, const double f[8], double *restrict E,
                                                       double *restrict Lx, double *restrict Ly, double *restrict Lz, double *restrict Q) {
#include "set_CodeParameters.h"
  // Unpack position coordinates from f[0]..f[3]
  const double t = f[0];
  const double x = f[1];
  const double y = f[2];
  const double z = f[3];

  // Unpack 4-momentum components from f[4]..f[7]
  // If particle is massive the mass is assumed to be m=1
  const double p0 = f[4];
  const double p1 = f[5];
  const double p2 = f[6];
  const double p3 = f[7];

  const REAL tmp0 = ((a_spin) * (a_spin));
  const REAL tmp4 = ((z) * (z));
  const REAL tmp6 = ((x) * (x)) + ((y) * (y));
  const REAL tmp7 = (1.0 / 2.0) * tmp4 + (1.0 / 2.0) * ((x) * (x)) + (1.0 / 2.0) * ((y) * (y)) +
                    (1.0 / 2.0) * sqrt(4 * tmp0 * tmp4 + ((-tmp0 + tmp4 + tmp6) * (-tmp0 + tmp4 + tmp6)));
  const REAL tmp30 = sqrt(tmp6);
  const REAL tmp8 = -1.0 / 2.0 * tmp0 + tmp7;
  const REAL tmp13 = (1.0 / 2.0) * tmp0 + tmp7;
  const REAL tmp9 = 2 * M_scale / (tmp0 * tmp4 + ((tmp8) * (tmp8)));
  const REAL tmp14 = (1.0 / (tmp13));
  const REAL tmp16 = sqrt(tmp8);
  const REAL tmp31 = sqrt(tmp13);
  const REAL tmp10 = tmp8 * tmp9 * z;
  const REAL tmp12 = pow(tmp8, 3.0 / 2.0) * tmp9;
  const REAL tmp17 = a_spin * y + tmp16 * x;
  const REAL tmp19 = -a_spin * x + tmp16 * y;
  const REAL tmp15 = tmp12 * tmp14;
  const REAL tmp18 = p1 * tmp17;
  const REAL tmp26 = tmp12 / ((tmp13) * (tmp13));
  const REAL tmp24 = p0 * tmp10 + p2 * tmp10 * tmp14 * tmp19 + p3 * (tmp16 * tmp4 * tmp9 + 1) + tmp10 * tmp14 * tmp18;
  const REAL tmp25 = p3 * tmp10 * tmp14;
  const REAL tmp21 = -p0 * (tmp12 - 1) - p2 * tmp15 * tmp19 - p3 * tmp10 - tmp15 * tmp18;
  const REAL tmp27 = p0 * tmp15 * tmp19 + p2 * (((tmp19) * (tmp19)) * tmp26 + 1) + tmp18 * tmp19 * tmp26 + tmp19 * tmp25;
  const REAL tmp28 = p0 * tmp15 * tmp17 + p1 * (((tmp17) * (tmp17)) * tmp26 + 1) + p2 * tmp17 * tmp19 * tmp26 + tmp17 * tmp25;
  const REAL tmp29 = tmp27 * x - tmp28 * y;
  *E = tmp21;
  *Lx = tmp24 * y - tmp27 * z;
  *Ly = -tmp24 * x + tmp28 * z;
  *Lz = tmp29;
  *Q = tmp4 * (tmp0 * (1 - ((tmp21) * (tmp21))) + tmp13 * ((tmp29) * (tmp29)) / tmp6) / tmp8 +
       ((-tmp16 * tmp24 * tmp30 / tmp31 + tmp31 * z * (tmp27 * y + tmp28 * x) / (tmp16 * tmp30)) *
        (-tmp16 * tmp24 * tmp30 / tmp31 + tmp31 * z * (tmp27 * y + tmp28 * x) / (tmp16 * tmp30)));
} // END FUNCTION conserved_quantities_KerrSchild_Cartesian_massive
