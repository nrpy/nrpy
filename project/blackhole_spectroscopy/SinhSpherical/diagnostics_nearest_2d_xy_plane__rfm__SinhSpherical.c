#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/*
 * Output diagnostic quantities at gridpoints closest to xy plane.
 */
void diagnostics_nearest_2d_xy_plane__rfm__SinhSpherical(commondata_struct *restrict commondata, const params_struct *restrict params,
                                                         REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs) {
#include "../set_CodeParameters.h"

  // Unpack gridfuncs struct:
  const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
  const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
  const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;

  // 1D output
  char filename[256];
  sprintf(filename, "out2d-xy-conv_factor%.2f-t%08.2f.txt", convergence_factor, time);
  FILE *outfile = (nn == 0) ? fopen(filename, "w") : fopen(filename, "a");
  if (!outfile) {
    fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
    exit(1);
  }

  // Output data in xy-plane in SinhSpherical coordinates.
  const int numpts_i0 = Nxx0, numpts_i1 = 1, numpts_i2 = Nxx2;
  int i0_pts[numpts_i0], i1_pts[numpts_i1], i2_pts[numpts_i2];
#pragma omp parallel for
  for (int i0 = NGHOSTS; i0 < Nxx0 + NGHOSTS; i0++)
    i0_pts[i0 - NGHOSTS] = i0;
  i1_pts[0] = (int)((1.0 / 2.0) * Nxx_plus_2NGHOSTS1);
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx2 + NGHOSTS; i2++)
    i2_pts[i2 - NGHOSTS] = i2;
  // Main loop:
  LOOP_NOOMP(i0_pt, 0, numpts_i0, i1_pt, 0, numpts_i1, i2_pt, 0, numpts_i2) {
    const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
    const int idx3 = IDX3(i0, i1, i2);
    REAL xCart[3];
    xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);
    {
      const REAL log10HL = log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16));
      const REAL log10sqrtM2L = log10(sqrt(diagnostic_output_gfs[IDX4pt(MSQUAREDGF, idx3)]) + 1e-16);
      const REAL cfL = y_n_gfs[IDX4pt(CFGF, idx3)];
      const REAL alphaL = y_n_gfs[IDX4pt(ALPHAGF, idx3)];
      const REAL trKL = y_n_gfs[IDX4pt(TRKGF, idx3)];
      fprintf(outfile, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e\n", xCart[0], xCart[1], log10HL, log10sqrtM2L, cfL, alphaL, trKL);
    }
  }

  fclose(outfile);
}
