#include "bssn_defines.h"

/**
 * Per-block Minkowski initial data fill (all EVOL fields to their asymptotic values).
 */
void bssn_minkowski_initial_data_block(const block_geometry_struct &geom, DendroScalar *const *out_gfs) {
  DendroScalar *const out_aDD00 = out_gfs[0] + geom.component_offset;
  DendroScalar *const out_aDD01 = out_gfs[1] + geom.component_offset;
  DendroScalar *const out_aDD02 = out_gfs[2] + geom.component_offset;
  DendroScalar *const out_aDD11 = out_gfs[3] + geom.component_offset;
  DendroScalar *const out_aDD12 = out_gfs[4] + geom.component_offset;
  DendroScalar *const out_aDD22 = out_gfs[5] + geom.component_offset;
  DendroScalar *const out_alpha = out_gfs[6] + geom.component_offset;
  DendroScalar *const out_betU0 = out_gfs[7] + geom.component_offset;
  DendroScalar *const out_betU1 = out_gfs[8] + geom.component_offset;
  DendroScalar *const out_betU2 = out_gfs[9] + geom.component_offset;
  DendroScalar *const out_cf = out_gfs[10] + geom.component_offset;
  DendroScalar *const out_hDD00 = out_gfs[11] + geom.component_offset;
  DendroScalar *const out_hDD01 = out_gfs[12] + geom.component_offset;
  DendroScalar *const out_hDD02 = out_gfs[13] + geom.component_offset;
  DendroScalar *const out_hDD11 = out_gfs[14] + geom.component_offset;
  DendroScalar *const out_hDD12 = out_gfs[15] + geom.component_offset;
  DendroScalar *const out_hDD22 = out_gfs[16] + geom.component_offset;
  DendroScalar *const out_lambdaU0 = out_gfs[17] + geom.component_offset;
  DendroScalar *const out_lambdaU1 = out_gfs[18] + geom.component_offset;
  DendroScalar *const out_lambdaU2 = out_gfs[19] + geom.component_offset;
  DendroScalar *const out_trK = out_gfs[20] + geom.component_offset;
  DendroScalar *const out_vetU0 = out_gfs[21] + geom.component_offset;
  DendroScalar *const out_vetU1 = out_gfs[22] + geom.component_offset;
  DendroScalar *const out_vetU2 = out_gfs[23] + geom.component_offset;
  const std::ptrdiff_t nx = static_cast<std::ptrdiff_t>(geom.nx);
  const std::ptrdiff_t ny = static_cast<std::ptrdiff_t>(geom.ny);
  const std::ptrdiff_t nz = static_cast<std::ptrdiff_t>(geom.nz);
  [[maybe_unused]] const std::ptrdiff_t nxy = nx * ny;
  const std::ptrdiff_t padding = static_cast<std::ptrdiff_t>(0);
  [[maybe_unused]] const DendroScalar invdxx0 = static_cast<DendroScalar>(1) / geom.dx[0];
  [[maybe_unused]] const DendroScalar invdxx1 = static_cast<DendroScalar>(1) / geom.dx[1];
  [[maybe_unused]] const DendroScalar invdxx2 = static_cast<DendroScalar>(1) / geom.dx[2];
  for (int i2 = static_cast<int>(padding); i2 < static_cast<int>(nz - padding); i2++) {
    for (int i1 = static_cast<int>(padding); i1 < static_cast<int>(ny - padding); i1++) {
      for (int i0 = static_cast<int>(padding); i0 < static_cast<int>(nx - padding); i0++) {
        const std::ptrdiff_t pp = i0 + nx * (i1 + ny * i2);
        [[maybe_unused]] const DendroScalar xx0 = geom.pmin_padded[0] + static_cast<DendroScalar>(i0) * geom.dx[0];
        [[maybe_unused]] const DendroScalar xx1 = geom.pmin_padded[1] + static_cast<DendroScalar>(i1) * geom.dx[1];
        [[maybe_unused]] const DendroScalar xx2 = geom.pmin_padded[2] + static_cast<DendroScalar>(i2) * geom.dx[2];
        out_aDD00[pp] = DendroScalar{0.0};
        out_aDD01[pp] = DendroScalar{0.0};
        out_aDD02[pp] = DendroScalar{0.0};
        out_aDD11[pp] = DendroScalar{0.0};
        out_aDD12[pp] = DendroScalar{0.0};
        out_aDD22[pp] = DendroScalar{0.0};
        out_alpha[pp] = DendroScalar{1.0};
        out_betU0[pp] = DendroScalar{0.0};
        out_betU1[pp] = DendroScalar{0.0};
        out_betU2[pp] = DendroScalar{0.0};
        out_cf[pp] = DendroScalar{1.0};
        out_hDD00[pp] = DendroScalar{0.0};
        out_hDD01[pp] = DendroScalar{0.0};
        out_hDD02[pp] = DendroScalar{0.0};
        out_hDD11[pp] = DendroScalar{0.0};
        out_hDD12[pp] = DendroScalar{0.0};
        out_hDD22[pp] = DendroScalar{0.0};
        out_lambdaU0[pp] = DendroScalar{0.0};
        out_lambdaU1[pp] = DendroScalar{0.0};
        out_lambdaU2[pp] = DendroScalar{0.0};
        out_trK[pp] = DendroScalar{0.0};
        out_vetU0[pp] = DendroScalar{0.0};
        out_vetU1[pp] = DendroScalar{0.0};
        out_vetU2[pp] = DendroScalar{0.0};
      } // END LOOP: for i0 over [static_cast<int>(padding), static_cast<int>(nx - padding))
    } // END LOOP: for i1 over [static_cast<int>(padding), static_cast<int>(ny - padding))
  } // END LOOP: for i2 over [static_cast<int>(padding), static_cast<int>(nz - padding))
} // END FUNCTION: bssn_minkowski_initial_data_block
