#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"

typedef struct {
  REAL xCart_axis;
  REAL log10HL;
  REAL log10sqrtM2L;
  REAL cfL;
  REAL alphaL;
  REAL trKL;
} data_point_1d_struct;

// qsort() comparison function for 1D output.
static int compare(const void *a, const void *b) {
  REAL l = ((data_point_1d_struct *)a)->xCart_axis;
  REAL r = ((data_point_1d_struct *)b)->xCart_axis;
  return (l > r) - (l < r);
}

/*
 * Output diagnostic quantities at gridpoints closest to y axis.
 */
void diagnostics_nearest_1d_y_axis__rfm__SinhSpherical(commondata_struct *restrict commondata, const params_struct *restrict params,
                                                       REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs) {
#include "../set_CodeParameters.h"

  // Unpack gridfuncs struct:
  const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
  const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
  const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;

  // 1D output
  char filename[256];
  sprintf(filename, "out1d-y-conv_factor%.2f-t%08.2f.txt", convergence_factor, time);
  FILE *outfile = (nn == 0) ? fopen(filename, "w") : fopen(filename, "a");
  if (!outfile) {
    fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
    exit(1);
  }

  // Output along y-axis in SinhSpherical coordinates.
  const int numpts_i0 = Nxx0, numpts_i1 = 1, numpts_i2 = 2;
  int i0_pts[numpts_i0], i1_pts[numpts_i1], i2_pts[numpts_i2];

  data_point_1d_struct data_points[numpts_i0 * numpts_i1 * numpts_i2];
  int data_index = 0;
#pragma omp parallel for
  for (int i0 = NGHOSTS; i0 < Nxx0 + NGHOSTS; i0++)
    i0_pts[i0 - NGHOSTS] = i0;
  i1_pts[0] = (int)((1.0 / 2.0) * Nxx_plus_2NGHOSTS1);
  i2_pts[0] = (int)(NGHOSTS + (1.0 / 4.0) * Nxx2 - 1.0 / 2.0);
  i2_pts[1] = (int)(NGHOSTS + (3.0 / 4.0) * Nxx2 - 1.0 / 2.0);
  // Main loop:
  LOOP_NOOMP(i0_pt, 0, numpts_i0, i1_pt, 0, numpts_i1, i2_pt, 0, numpts_i2) {
    const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
    const int idx3 = IDX3(i0, i1, i2);
    REAL xCart[3];
    xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);

    {
      data_point_1d_struct dp1d;
      dp1d.xCart_axis = xCart[1];
      dp1d.log10HL = log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16));
      dp1d.log10sqrtM2L = log10(sqrt(diagnostic_output_gfs[IDX4pt(MSQUAREDGF, idx3)]) + 1e-16);
      dp1d.cfL = y_n_gfs[IDX4pt(CFGF, idx3)];
      dp1d.alphaL = y_n_gfs[IDX4pt(ALPHAGF, idx3)];
      dp1d.trKL = y_n_gfs[IDX4pt(TRKGF, idx3)];
      data_points[data_index] = dp1d;
      data_index++;
    }
  }

  qsort(data_points, data_index, sizeof(data_point_1d_struct), compare);

  for (int i = 0; i < data_index; i++) {
    fprintf(outfile, "%.15e %.15e %.15e %.15e %.15e %.15e\n", data_points[i].xCart_axis, data_points[i].log10HL, data_points[i].log10sqrtM2L,
            data_points[i].cfL, data_points[i].alphaL, data_points[i].trKL);
  }

  fclose(outfile);
}
