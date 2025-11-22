#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "diagnostics/diagnostics_nearest_common.h"

// Data point for sorting by physical axis coordinate
typedef struct {
  REAL coord; // physical y or z
  int idx3;   // 3D index
  int i0;
  int i1;
  int i2;
} data_point_1d_struct;

// qsort comparator (file scope to satisfy -std=c11)
static int compare_by_coord(const void *a, const void *b) {
  const REAL lv = ((const data_point_1d_struct *)a)->coord;
  const REAL rv = ((const data_point_1d_struct *)b)->coord;
  return (lv > rv) - (lv < rv);
} // END FUNCTION compare_by_coord

void diagnostics_nearest_1d_y_and_z_axes__rfm__SinhSpherical(
    commondata_struct *restrict commondata, const int grid,
    const params_struct *restrict params,
    const params_struct *restrict params_chare,
    const REAL *restrict xx[3], const REAL *restrict xx_chare[3],
    const int NUM_GFS_NEAREST, const int which_gfs[],
    const char **diagnostic_gf_names,
    const REAL *restrict gridfuncs_diags[],
    const charecomm_struct *restrict charecommstruct,
    diagnostic_struct *restrict diagnosticstruct,
    const int chare_index[3], Ck::IO::Session token,
    const int which_diagnostics_part) {

#include "set_CodeParameters.h"

  switch (which_diagnostics_part) {

    case DIAGNOSTICS_SETUP_1D: {

      const int Nxx0chare = params_chare->Nxx0;
      const int Nxx1chare = params_chare->Nxx1;
      const int Nxx2chare = params_chare->Nxx2;

      // Build filename component with runtime coordinate system name and grid number
      char coordsys_with_grid[128];
      snprintf(coordsys_with_grid, sizeof(coordsys_with_grid),
               "%s-grid%02d", params->CoordSystemName, grid);
      strcpy(diagnosticstruct->filename_2d_xy, coordsys_with_grid);
      strcpy(diagnosticstruct->filename_2d_yz, coordsys_with_grid);

      diagnosticstruct->num_output_quantities = NUM_GFS_NEAREST;

      // Interior counts
      const int N0int = params->Nxx_plus_2NGHOSTS0 - 2 * NGHOSTS;
      const int N1int = params->Nxx_plus_2NGHOSTS1 - 2 * NGHOSTS;
      const int N2int = params->Nxx_plus_2NGHOSTS2 - 2 * NGHOSTS;

      // Common fixed-point helpers
      MAYBE_UNUSED const int i0_mid = params->Nxx_plus_2NGHOSTS0 / 2;
      MAYBE_UNUSED const int i1_mid = params->Nxx_plus_2NGHOSTS1 / 2;
      MAYBE_UNUSED const int i2_mid = params->Nxx_plus_2NGHOSTS2 / 2;

      MAYBE_UNUSED const int i1_min = NGHOSTS;
      MAYBE_UNUSED const int i1_max =
          params->Nxx_plus_2NGHOSTS1 - NGHOSTS - 1;
      MAYBE_UNUSED const int i2_min = NGHOSTS;

      MAYBE_UNUSED const int i0_rmin = NGHOSTS; // rho = 0 (Cylindrical)
      MAYBE_UNUSED const int i1_pmin = NGHOSTS; // phi = -pi

      // Quarter-plane indices for cell-centered grids
      MAYBE_UNUSED const int i2_q1 =
          (int)(NGHOSTS + 0.25 * (REAL)N2int - 0.5);
      MAYBE_UNUSED const int i2_q3 =
          (int)(NGHOSTS + 0.75 * (REAL)N2int - 0.5);
      MAYBE_UNUSED const int i1_q1 =
          (int)(NGHOSTS + 0.25 * (REAL)N1int - 0.5);
      MAYBE_UNUSED const int i1_q3 =
          (int)(NGHOSTS + 0.75 * (REAL)N1int - 0.5);

      // Allocate buffers (tight upper bounds)
      const int max_y = N0int + N0int;
      const int max_z = N0int + N0int;
      data_point_1d_struct *data_points_y =
          max_y > 0 ? (data_point_1d_struct *)malloc(
                          sizeof(data_point_1d_struct) * (size_t)max_y)
                    : NULL;
      data_point_1d_struct *data_points_z =
          max_z > 0 ? (data_point_1d_struct *)malloc(
                          sizeof(data_point_1d_struct) * (size_t)max_z)
                    : NULL;

      // ----------------------
      // Build y-axis samples
      // ----------------------
      int count_y = 0;
      for (int i0 = NGHOSTS;
           i0 < params->Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        const int i1 = i1_mid;
        const int i2 = i2_q1;
        const int idx3 = IDX3P(params, i0, i1, i2);
        REAL xCart[3], xOrig[3] = {xx[0][i0], xx[1][i1], xx[2][i2]};
        xx_to_Cart(params, xOrig, xCart);
        data_points_y[count_y].coord = xCart[1];
        data_points_y[count_y].idx3 = idx3;
        data_points_y[count_y].i0 = i0;
        data_points_y[count_y].i1 = i1;
        data_points_y[count_y].i2 = i2;
        count_y++;
      }
      for (int i0 = NGHOSTS;
           i0 < params->Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        const int i1 = i1_mid;
        const int i2 = i2_q3;
        const int idx3 = IDX3P(params, i0, i1, i2);
        REAL xCart[3], xOrig[3] = {xx[0][i0], xx[1][i1], xx[2][i2]};
        xx_to_Cart(params, xOrig, xCart);
        data_points_y[count_y].coord = xCart[1];
        data_points_y[count_y].idx3 = idx3;
        data_points_y[count_y].i0 = i0;
        data_points_y[count_y].i1 = i1;
        data_points_y[count_y].i2 = i2;
        count_y++;
      }

      if (count_y > 1)
        qsort(data_points_y, (size_t)count_y,
              sizeof(data_point_1d_struct), compare_by_coord);

      diagnosticstruct->tot_num_diagnostic_1d_y_pts = count_y;

      int num_diagnostics_chare = 0;
      for (int i = 0; i < count_y; i++) {
        const int i0 = data_points_y[i].i0;
        const int i1 = data_points_y[i].i1;
        const int i2 = data_points_y[i].i2;
        const int idx3 = IDX3(i0, i1, i2);
        if (charecommstruct->globalidx3pt_to_chareidx3[idx3] ==
            IDX3_OF_CHARE(chare_index[0], chare_index[1],
                          chare_index[2])) {
          num_diagnostics_chare++;
        }
      }
      diagnosticstruct->num_diagnostic_1d_y_pts = num_diagnostics_chare;
      diagnosticstruct->localidx3_diagnostic_1d_y_pt =
          (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->locali0_diagnostic_1d_y_pt =
          (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->locali1_diagnostic_1d_y_pt =
          (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->locali2_diagnostic_1d_y_pt =
          (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->offset_diagnostic_1d_y_pt =
          (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);

      int sizeinbytes = 24 * (diagnosticstruct->num_output_quantities + 1);
      int which_diagnostics_chare = 0;
      int which_diagnostic_global = 0;
      for (int i = 0; i < count_y; i++) {
        const int i0 = data_points_y[i].i0;
        const int i1 = data_points_y[i].i1;
        const int i2 = data_points_y[i].i2;
        const int idx3 = IDX3(i0, i1, i2);
        if (charecommstruct->globalidx3pt_to_chareidx3[idx3] ==
            IDX3_OF_CHARE(chare_index[0], chare_index[1],
                          chare_index[2])) {
          int localidx3 =
              charecommstruct->globalidx3pt_to_localidx3pt[idx3];
          diagnosticstruct
              ->localidx3_diagnostic_1d_y_pt[which_diagnostics_chare] =
              localidx3;
          diagnosticstruct
              ->locali0_diagnostic_1d_y_pt[which_diagnostics_chare] =
              MAP_GLOBAL_TO_LOCAL_IDX0(chare_index[0], i0, Nxx0chare);
          diagnosticstruct
              ->locali1_diagnostic_1d_y_pt[which_diagnostics_chare] =
              MAP_GLOBAL_TO_LOCAL_IDX1(chare_index[1], i1, Nxx1chare);
          diagnosticstruct
              ->locali2_diagnostic_1d_y_pt[which_diagnostics_chare] =
              MAP_GLOBAL_TO_LOCAL_IDX2(chare_index[2], i2, Nxx2chare);
          diagnosticstruct
              ->offset_diagnostic_1d_y_pt[which_diagnostics_chare] =
              which_diagnostic_global * sizeinbytes;
          which_diagnostics_chare++;
        }
        which_diagnostic_global++;
      }

      // ----------------------
      // Build z-axis samples
      // ----------------------
      int count_z = 0;
      for (int i0 = NGHOSTS;
           i0 < params->Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        const int i1 = i1_min;
        const int i2 = i2_min;
        const int idx3 = IDX3P(params, i0, i1, i2);
        REAL xCart[3], xOrig[3] = {xx[0][i0], xx[1][i1], xx[2][i2]};
        xx_to_Cart(params, xOrig, xCart);
        data_points_z[count_z].coord = xCart[2];
        data_points_z[count_z].idx3 = idx3;
        data_points_z[count_z].i0 = i0;
        data_points_z[count_z].i1 = i1;
        data_points_z[count_z].i2 = i2;
        count_z++;
      }
      for (int i0 = NGHOSTS;
           i0 < params->Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        const int i1 = i1_max;
        const int i2 = i2_min;
        const int idx3 = IDX3P(params, i0, i1, i2);
        REAL xCart[3], xOrig[3] = {xx[0][i0], xx[1][i1], xx[2][i2]};
        xx_to_Cart(params, xOrig, xCart);
        data_points_z[count_z].coord = xCart[2];
        data_points_z[count_z].idx3 = idx3;
        data_points_z[count_z].i0 = i0;
        data_points_z[count_z].i1 = i1;
        data_points_z[count_z].i2 = i2;
        count_z++;
      }

      if (count_z > 1)
        qsort(data_points_z, (size_t)count_z,
              sizeof(data_point_1d_struct), compare_by_coord);

      num_diagnostics_chare = 0;
      for (int i = 0; i < count_z; i++) {
        const int i0 = data_points_z[i].i0;
        const int i1 = data_points_z[i].i1;
        const int i2 = data_points_z[i].i2;
        const int idx3 = IDX3(i0, i1, i2);
        if (charecommstruct->globalidx3pt_to_chareidx3[idx3] ==
            IDX3_OF_CHARE(chare_index[0], chare_index[1],
                          chare_index[2])) {
          num_diagnostics_chare++;
        }
      }
      diagnosticstruct->num_diagnostic_1d_z_pts = num_diagnostics_chare;
      diagnosticstruct->localidx3_diagnostic_1d_z_pt =
          (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->locali0_diagnostic_1d_z_pt =
          (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->locali1_diagnostic_1d_z_pt =
          (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->locali2_diagnostic_1d_z_pt =
          (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->offset_diagnostic_1d_z_pt =
          (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);

      int sizeinbztes =
          24 * (diagnosticstruct->num_output_quantities + 1);
      which_diagnostics_chare = 0;
      which_diagnostic_global = 0;
      for (int i = 0; i < count_z; i++) {
        const int i0 = data_points_z[i].i0;
        const int i1 = data_points_z[i].i1;
        const int i2 = data_points_z[i].i2;
        const int idx3 = IDX3(i0, i1, i2);
        if (charecommstruct->globalidx3pt_to_chareidx3[idx3] ==
            IDX3_OF_CHARE(chare_index[0], chare_index[1],
                          chare_index[2])) {
          int localidx3 =
              charecommstruct->globalidx3pt_to_localidx3pt[idx3];
          diagnosticstruct
              ->localidx3_diagnostic_1d_z_pt[which_diagnostics_chare] =
              localidx3;
          diagnosticstruct
              ->locali0_diagnostic_1d_z_pt[which_diagnostics_chare] =
              MAP_GLOBAL_TO_LOCAL_IDX0(chare_index[0], i0, Nxx0chare);
          diagnosticstruct
              ->locali1_diagnostic_1d_z_pt[which_diagnostics_chare] =
              MAP_GLOBAL_TO_LOCAL_IDX1(chare_index[1], i1, Nxx1chare);
          diagnosticstruct
              ->locali2_diagnostic_1d_z_pt[which_diagnostics_chare] =
              MAP_GLOBAL_TO_LOCAL_IDX2(chare_index[2], i2, Nxx2chare);
          diagnosticstruct
              ->offset_diagnostic_1d_z_pt[which_diagnostics_chare] =
              which_diagnostic_global * sizeinbztes;
          which_diagnostics_chare++;
        }
        which_diagnostic_global++;
      }

      free(data_points_y);
      free(data_points_z);

      break;
    }

    case DIAGNOSTICS_WRITE_Y: {

      const REAL *restrict src = gridfuncs_diags[grid];
      const int NUM_COLS = 1 + NUM_GFS_NEAREST;
      REAL *row =
          (REAL *)malloc(sizeof(REAL) * (size_t)(NUM_COLS));
      if (!row) {
        fprintf(stderr,
                "Error: Failed to allocate memory for row buffer.\n");
        exit(1);
      }

      const int num_diagnostic_pts =
          diagnosticstruct->num_diagnostic_1d_y_pts;
      const int *restrict idx3_diagnostic_pt =
          diagnosticstruct->localidx3_diagnostic_1d_y_pt;
      const int *restrict i0_diagnostic_pt =
          diagnosticstruct->locali0_diagnostic_1d_y_pt;
      const int *restrict i1_diagnostic_pt =
          diagnosticstruct->locali1_diagnostic_1d_y_pt;
      const int *restrict i2_diagnostic_pt =
          diagnosticstruct->locali2_diagnostic_1d_y_pt;
      const int *restrict offsetpt_firstfield =
          diagnosticstruct->offset_diagnostic_1d_y_pt;

      int sizeinbytes =
          24 * (diagnosticstruct->num_output_quantities + 1);

      //~ fprintf(stderr,
              //~ "[WRITE_Y] grid=%d NUM_GFS_NEAREST=%d "
              //~ "num_output_quantities=%d num_pts=%d sizeinbytes=%d\n",
              //~ grid, NUM_GFS_NEAREST,
              //~ diagnosticstruct->num_output_quantities,
              //~ num_diagnostic_pts, sizeinbytes);
      //~ fflush(stderr);

      for (int which_pt = 0; which_pt < num_diagnostic_pts;
           which_pt++) {
        const int idx3 = idx3_diagnostic_pt[which_pt];
        const int i0 = i0_diagnostic_pt[which_pt];
        const int i1 = i1_diagnostic_pt[which_pt];
        const int i2 = i2_diagnostic_pt[which_pt];
        REAL xCart[3];
        REAL xOrig[3] = {xx_chare[0][i0], xx_chare[1][i1],
                         xx_chare[2][i2]};

        // minimal debug: before xx_to_Cart
        fprintf(stderr,
                "[WRITE_Y] pt=%d idx3=%d i=(%d,%d,%d) "
                "xOrig=(%.6e,%.6e,%.6e) params=%p\n",
                which_pt, idx3, i0, i1, i2,
                xOrig[0], xOrig[1], xOrig[2],
                (const void *)params);
        fflush(stderr);

        xx_to_Cart(params, xOrig, xCart);

        char out[sizeinbytes + 1];
        row[0] = xCart[1];
        for (int gf_i = 0; gf_i < NUM_GFS_NEAREST; ++gf_i) {
          const int gf = which_gfs[gf_i];
          row[1 + gf_i] =
              src[IDX4Ppt(params_chare, gf, idx3)];
        }

        int n = 0;
        n += snprintf(out + n,
                      (size_t)(sizeinbytes + 1 - n),
                      "% .15e", row[0]);
        for (int col = 1; col < NUM_COLS; col++) {
          n += snprintf(out + n,
                        (size_t)(sizeinbytes + 1 - n),
                        " % .15e", row[col]);
        }
        out[n++] = '\n';
        Ck::IO::write(token, out, sizeinbytes,
                      offsetpt_firstfield[which_pt]);
      }

      free(row);
      break;
    }

    case DIAGNOSTICS_WRITE_Z: {

      const REAL *restrict src = gridfuncs_diags[grid];
      const int NUM_COLS = 1 + NUM_GFS_NEAREST;
      REAL *row =
          (REAL *)malloc(sizeof(REAL) * (size_t)(NUM_COLS));
      if (!row) {
        fprintf(stderr,
                "Error: Failed to allocate memory for row buffer.\n");
        exit(1);
      }

      const int num_diagnostic_pts =
          diagnosticstruct->num_diagnostic_1d_z_pts;
      const int *restrict idx3_diagnostic_pt =
          diagnosticstruct->localidx3_diagnostic_1d_z_pt;
      const int *restrict i0_diagnostic_pt =
          diagnosticstruct->locali0_diagnostic_1d_z_pt;
      const int *restrict i1_diagnostic_pt =
          diagnosticstruct->locali1_diagnostic_1d_z_pt;
      const int *restrict i2_diagnostic_pt =
          diagnosticstruct->locali2_diagnostic_1d_z_pt;
      const int *restrict offsetpt_firstfield =
          diagnosticstruct->offset_diagnostic_1d_z_pt;

      int sizeinbytes =
          24 * (diagnosticstruct->num_output_quantities + 1);

      fprintf(stderr,
              "[WRITE_Z] grid=%d NUM_GFS_NEAREST=%d "
              "num_output_quantities=%d num_pts=%d sizeinbytes=%d\n",
              grid, NUM_GFS_NEAREST,
              diagnosticstruct->num_output_quantities,
              num_diagnostic_pts, sizeinbytes);
      fflush(stderr);

      for (int which_pt = 0; which_pt < num_diagnostic_pts;
           which_pt++) {
        const int idx3 = idx3_diagnostic_pt[which_pt];
        const int i0 = i0_diagnostic_pt[which_pt];
        const int i1 = i1_diagnostic_pt[which_pt];
        const int i2 = i2_diagnostic_pt[which_pt];
        REAL xCart[3];
        REAL xOrig[3] = {xx_chare[0][i0], xx_chare[1][i1],
                         xx_chare[2][i2]};

        fprintf(stderr,
                "[WRITE_Z] pt=%d idx3=%d i=(%d,%d,%d) "
                "xOrig=(%.6e,%.6e,%.6e) params=%p\n",
                which_pt, idx3, i0, i1, i2,
                xOrig[0], xOrig[1], xOrig[2],
                (const void *)params);
        fflush(stderr);

        xx_to_Cart(params, xOrig, xCart);

        char out[sizeinbytes + 1];
        row[0] = xCart[2];
        for (int gf_i = 0; gf_i < NUM_GFS_NEAREST; ++gf_i) {
          const int gf = which_gfs[gf_i];
          row[1 + gf_i] =
              src[IDX4Ppt(params_chare, gf, idx3)];
        }

        int n = 0;
        n += snprintf(out + n,
                      (size_t)(sizeinbytes + 1 - n),
                      "% .15e", row[0]);
        for (int col = 1; col < NUM_COLS; col++) {
          n += snprintf(out + n,
                        (size_t)(sizeinbytes + 1 - n),
                        " % .15e", row[col]);
        }
        out[n++] = '\n';
        Ck::IO::write(token, out, sizeinbytes,
                      offsetpt_firstfield[which_pt]);
      }

      free(row);
      break;
    }
  }

} // END FUNCTION diagnostics_nearest_1d_y_and_z_axes__rfm__SinhSpherical
