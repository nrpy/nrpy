#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "diagnostics/diagnostics_nearest_common.h"
#include "stdlib.h"

/**
 * @file diagnostics_nearest_2d_xy_and_yz_planes.c
 * @brief Sample and write 2D diagnostics on the nearest xy and yz planes for a single SinhSpherical grid.
 *
 * Overview:
 * For the specified grid at the current time, this routine:
 *  - Locates interior slices that are nearest to the xy and yz planes, with selection rules
 *    specialized by the runtime coordinate system.
 *  - Converts native coordinates xx to Cartesian coordinates (x, y, z) using xx_to_Cart so that
 *    the first columns of each row contain mapped coordinates.
 *  - Writes two files per call, one for the xy plane and one for the yz plane. Each file begins
 *    with a time comment and a header, followed by one row per interior point that contains the
 *    mapped coordinates and sampled diagnostic values.
 *  - Performs sampling without interpolation; values are read directly from gridfuncs_diags[grid].
 *
 * Plane selection notes (examples, not exhaustive):
 *  - Cartesian: xy fixes i2 at mid; yz fixes i0 at mid.
 *  - Cylindrical: xy fixes z at mid; yz emits two phi slices near +/- pi/2 to realize x=0.
 *  - Spherical and SymTP: xy fixes the polar-like angle at mid; yz emits two phi-like slices near +/- pi/2.
 *  - Wedge: xy may emit two z-like slices at quarter indices; yz fixes the across-wedge index at mid.
 *
 * If a user-editable block is provided in the implementation, users may insert custom logic such as
 * adding extra columns or filtering before rows are written.
 *
 * @param[in,out] commondata            Pointer to common runtime data used for time and I/O.
 * @param[in]     grid                  Grid index to process.
 * @param[in]     params                Pointer to simulation and grid parameters (sizes, names, strides).
 * @param[in]     xx                    Native grid coordinates; xx[d][i_d] gives the coordinate along dimension d.
 * @param[in]     NUM_GFS_NEAREST       Number of diagnostic gridfunctions to sample at each interior point.
 * @param[in]     which_gfs             Array of length NUM_GFS_NEAREST specifying which gridfunctions to sample.
 * @param[in]     diagnostic_gf_names   Array of length NUM_GFS_NEAREST with human-readable names for headers.
 * @param[in]     gridfuncs_diags       Array of pointers; gridfuncs_diags[grid] points to this grid's diagnostic data.
 *
 * @return        void                  No return value. On success two text files are written and closed. Fatal I/O
 *                                      or allocation failures result in program termination.
 */
void diagnostics_nearest_2d_xy_and_yz_planes__rfm__SinhSpherical(commondata_struct *restrict commondata, const int grid,
                                                                 const params_struct *restrict params, const params_struct *restrict params_chare, const REAL *restrict xx[3], const REAL *restrict xx_chare[3],
                                                                 const int NUM_GFS_NEAREST, const int which_gfs[], const char **diagnostic_gf_names,
                                                                 const REAL *restrict gridfuncs_diags[],
                                                                 const charecomm_struct *restrict charecommstruct, diagnostic_struct *restrict diagnosticstruct,
                                                                 const int chare_index[3], Ck::IO::Session token, const int which_diagnostics_part) {

#include "set_CodeParameters.h"


  switch (which_diagnostics_part) {

    case DIAGNOSTICS_SETUP_2D: {

      const int Nxx0chare = params_chare->Nxx0;
      const int Nxx1chare = params_chare->Nxx1;
      const int Nxx2chare = params_chare->Nxx2;

      // Build filename component with runtime coordinate system name and grid number
      char coordsys_with_grid[128];
      snprintf(coordsys_with_grid, sizeof(coordsys_with_grid), "%s-grid%02d", params->CoordSystemName, grid);

      char filename[256];

      build_outfile_name(filename, sizeof filename, "out2d-xy", coordsys_with_grid, commondata, /*include_time=*/0);
      strcpy(diagnosticstruct->filename_2d_xy, filename);

      build_outfile_name(filename, sizeof filename, "out2d-yz", coordsys_with_grid, commondata, /*include_time=*/0);
      strcpy(diagnosticstruct->filename_2d_yz, filename);

      diagnosticstruct->num_output_quantities = NUM_GFS_NEAREST;

      // compute offset in bytes for first field for each diagnostic pt
      int sizeinbytes = 23 * (diagnosticstruct->num_output_quantities + 2);

      // Interior grid counts and loop bounds
      MAYBE_UNUSED const int N0int = params->Nxx_plus_2NGHOSTS0 - 2 * NGHOSTS;
      MAYBE_UNUSED const int N1int = params->Nxx_plus_2NGHOSTS1 - 2 * NGHOSTS;
      MAYBE_UNUSED const int N2int = params->Nxx_plus_2NGHOSTS2 - 2 * NGHOSTS;
      const int i0_end = params->Nxx_plus_2NGHOSTS0 - NGHOSTS;
      const int i1_end = params->Nxx_plus_2NGHOSTS1 - NGHOSTS;
      const int i2_end = params->Nxx_plus_2NGHOSTS2 - NGHOSTS;

      // Fixed-point index helpers
      MAYBE_UNUSED const int i0_mid = params->Nxx_plus_2NGHOSTS0 / 2;
      MAYBE_UNUSED const int i1_mid = params->Nxx_plus_2NGHOSTS1 / 2;
      MAYBE_UNUSED const int i2_mid = params->Nxx_plus_2NGHOSTS2 / 2;
      MAYBE_UNUSED const int i1_q1 = (int)(NGHOSTS + 0.25 * (REAL)N1int - 0.5);
      MAYBE_UNUSED const int i1_q3 = (int)(NGHOSTS + 0.75 * (REAL)N1int - 0.5);
      MAYBE_UNUSED const int i2_q1 = (int)(NGHOSTS + 0.25 * (REAL)N2int - 0.5);
      MAYBE_UNUSED const int i2_q3 = (int)(NGHOSTS + 0.75 * (REAL)N2int - 0.5);

      // --- Sample and setup data for the xy-plane---
      {
        diagnosticstruct->tot_num_diagnostic_2d_xy_pts = N0int * N2int;

        int num_diagnostics_chare = 0;
        int i1 = i1_mid;
        for (int i2 = NGHOSTS; i2 < i2_end; i2++) {
          for (int i0 = NGHOSTS; i0 < i0_end; i0++) {
            const int idx3 = IDX3P(params, i0, i1, i2);
            if (charecommstruct->globalidx3pt_to_chareidx3[idx3] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])) {
              num_diagnostics_chare++;
            }
          } // END LOOP over i0
        } // END LOOP over i2

        diagnosticstruct->num_diagnostic_2d_xy_pts = num_diagnostics_chare;
        diagnosticstruct->locali0_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->locali1_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->locali2_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->localidx3_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->offset_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);




        int which_diagnostics_chare = 0;
        int which_diagnostic_global = 0;
        i1 = i1_mid;
        for (int i2 = NGHOSTS; i2 < i2_end; i2++) {
          for (int i0 = NGHOSTS; i0 < i0_end; i0++) {
            const int idx3 = IDX3P(params, i0, i1, i2);
            if (charecommstruct->globalidx3pt_to_chareidx3[idx3] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])) {
              // store the local idx3 of diagnostic point
              int localidx3 = charecommstruct->globalidx3pt_to_localidx3pt[idx3];
              diagnosticstruct->localidx3_diagnostic_2d_xy_pt[which_diagnostics_chare] = localidx3;
              diagnosticstruct->locali0_diagnostic_2d_xy_pt[which_diagnostics_chare] = MAP_GLOBAL_TO_LOCAL_IDX0(chare_index[0], i0, Nxx0chare);
              diagnosticstruct->locali1_diagnostic_2d_xy_pt[which_diagnostics_chare] = MAP_GLOBAL_TO_LOCAL_IDX1(chare_index[1], i1, Nxx1chare);
              diagnosticstruct->locali2_diagnostic_2d_xy_pt[which_diagnostics_chare] = MAP_GLOBAL_TO_LOCAL_IDX2(chare_index[2], i2, Nxx2chare);
              diagnosticstruct->offset_diagnostic_2d_xy_pt[which_diagnostics_chare] = which_diagnostic_global * sizeinbytes;
              which_diagnostics_chare++;
            }
            which_diagnostic_global++;
          }
        }
      }

      // --- Sample and setup data for the yz-plane---
      {
        diagnosticstruct->tot_num_diagnostic_2d_yz_pts = 2* N0int * N1int; //2 slices??

        int num_diagnostics_chare = 0;
        int i2_slices[2] = {i2_q1, i2_q3};
        for (int slice = 0; slice < 2; slice++) {
          const int i2 = i2_slices[slice];
          for (int i1 = NGHOSTS; i1 < i1_end; i1++) {
            for (int i0 = NGHOSTS; i0 < i0_end; i0++) {
              const int idx3 = IDX3P(params, i0, i1, i2);
              if (charecommstruct->globalidx3pt_to_chareidx3[idx3] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])) {
                num_diagnostics_chare++;
              }
            } // END LOOP over i0
          } // END LOOP over i1
        } // END LOOP over slices

        diagnosticstruct->num_diagnostic_2d_yz_pts = num_diagnostics_chare;
        diagnosticstruct->locali0_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->locali1_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->locali2_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->localidx3_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->offset_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);


        int which_diagnostics_chare = 0;
        int which_diagnostic_global = 0;
        for (int slice = 0; slice < 2; slice++) {
          const int i2 = i2_slices[slice];
          for (int i1 = NGHOSTS; i1 < i1_end; i1++) {
            for (int i0 = NGHOSTS; i0 < i0_end; i0++) {
              const int idx3 = IDX3P(params, i0, i1, i2);
              if (charecommstruct->globalidx3pt_to_chareidx3[idx3] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])) {
                // store the local idx3 of diagnostic point
                int localidx3 = charecommstruct->globalidx3pt_to_localidx3pt[idx3];
                diagnosticstruct->localidx3_diagnostic_2d_yz_pt[which_diagnostics_chare] = localidx3;
                diagnosticstruct->locali0_diagnostic_2d_yz_pt[which_diagnostics_chare] = MAP_GLOBAL_TO_LOCAL_IDX0(chare_index[0], i0, Nxx0chare);
                diagnosticstruct->locali1_diagnostic_2d_yz_pt[which_diagnostics_chare] = MAP_GLOBAL_TO_LOCAL_IDX1(chare_index[1], i1, Nxx1chare);
                diagnosticstruct->locali2_diagnostic_2d_yz_pt[which_diagnostics_chare] = MAP_GLOBAL_TO_LOCAL_IDX2(chare_index[2], i2, Nxx2chare);
                diagnosticstruct->offset_diagnostic_2d_yz_pt[which_diagnostics_chare] = which_diagnostic_global * sizeinbytes;
                which_diagnostics_chare++;
              }
              which_diagnostic_global++;
            }// END LOOP over i0
          } // END LOOP over i1
        } // END LOOP over slices
      }

      break;
    }

    case DIAGNOSTICS_WRITE_XY: {

      // Active grid data pointer and reusable row buffer
      const REAL *restrict src = gridfuncs_diags[grid];
      const int NUM_COLS = 2 + NUM_GFS_NEAREST;
      REAL *row = (REAL *)malloc(sizeof(REAL) * (size_t)NUM_COLS);
      if (!row) {
        fprintf(stderr, "Error: Failed to allocate memory for row buffer.\n");
        exit(1);
      } // END IF row allocation failure

      // Unpack diagnosticptoffset struct:
      const int num_diagnostic_pts = diagnosticstruct->num_diagnostic_2d_xy_pts;
      const int *restrict idx3_diagnostic_pt = diagnosticstruct->localidx3_diagnostic_2d_xy_pt;
      const int *restrict i0_diagnostic_pt = diagnosticstruct->locali0_diagnostic_2d_xy_pt;
      const int *restrict i1_diagnostic_pt = diagnosticstruct->locali1_diagnostic_2d_xy_pt;
      const int *restrict i2_diagnostic_pt = diagnosticstruct->locali2_diagnostic_2d_xy_pt;
      const int *restrict offsetpt_firstfield = diagnosticstruct->offset_diagnostic_2d_xy_pt;

      for (int which_pt = 0; which_pt < num_diagnostic_pts; which_pt++) {
        const int idx3 = idx3_diagnostic_pt[which_pt];
        const int i0 = i0_diagnostic_pt[which_pt];
        const int i1 = i1_diagnostic_pt[which_pt];
        const int i2 = i2_diagnostic_pt[which_pt];
        REAL xCart[3];

        //~ REAL xOrig[3] = {xx[0][i0], xx[1][i1], xx[2][i2]};
        //~ xx_to_Cart(params, xOrig, xCart);
        REAL xOrig[3] = {xx_chare[0][i0], xx_chare[1][i1], xx_chare[2][i2]};
        xx_to_Cart(params_chare, xOrig, xCart);

        int sizeinbytes = 23 * (diagnosticstruct->num_output_quantities + 2);
        char out[sizeinbytes + 1];
        row[0] = xCart[0];
        row[1] = xCart[1];
        for (int gf_idx = 0; gf_idx < NUM_GFS_NEAREST; gf_idx++) {
          const int gf = which_gfs[gf_idx];
          row[2 + gf_idx] = src[IDX4Ppt(params_chare, gf, idx3)];
        } // END LOOP over gridfunctions

        int n = 0;
        n += sprintf(out + n, "% .15e", row[0]);
        for (int col = 1; col < NUM_COLS; col++) {
          n += sprintf(out + n, " % .15e", row[col]);
        }
        out[sizeinbytes - 1] = '\n';
        Ck::IO::write(token, out, sizeinbytes, offsetpt_firstfield[which_pt]);
      }
      // Finalize
      free(row);
      break;
    }

    case DIAGNOSTICS_WRITE_YZ: {

      // Active grid data pointer and reusable row buffer
      const REAL *restrict src = gridfuncs_diags[grid];
      const int NUM_COLS = 2 + NUM_GFS_NEAREST;
      REAL *row = (REAL *)malloc(sizeof(REAL) * (size_t)NUM_COLS);
      if (!row) {
        fprintf(stderr, "Error: Failed to allocate memory for row buffer.\n");
        exit(1);
      } // END IF row allocation failure

      // Unpack diagnosticptoffset struct:
      const int num_diagnostic_pts = diagnosticstruct->num_diagnostic_2d_yz_pts;
      const int *restrict idx3_diagnostic_pt = diagnosticstruct->localidx3_diagnostic_2d_yz_pt;
      const int *restrict i0_diagnostic_pt = diagnosticstruct->locali0_diagnostic_2d_yz_pt;
      const int *restrict i1_diagnostic_pt = diagnosticstruct->locali1_diagnostic_2d_yz_pt;
      const int *restrict i2_diagnostic_pt = diagnosticstruct->locali2_diagnostic_2d_yz_pt;
      const int *restrict offsetpt_firstfield = diagnosticstruct->offset_diagnostic_2d_yz_pt;

      for (int which_pt = 0; which_pt < num_diagnostic_pts; which_pt++) {
        const int idx3 = idx3_diagnostic_pt[which_pt];
        const int i0 = i0_diagnostic_pt[which_pt];
        const int i1 = i1_diagnostic_pt[which_pt];
        const int i2 = i2_diagnostic_pt[which_pt];
        REAL xCart[3];
        //~ REAL xOrig[3] = {xx[0][i0], xx[1][i1], xx[2][i2]};
        //~ xx_to_Cart(params, xOrig, xCart);
        REAL xOrig[3] = {xx_chare[0][i0], xx_chare[1][i1], xx_chare[2][i2]};
        xx_to_Cart(params_chare, xOrig, xCart);


        int sizeinbytes = 23 * (diagnosticstruct->num_output_quantities + 2);
        char out[sizeinbytes + 1];
        row[0] = xCart[1];
        row[1] = xCart[2];
        for (int gf_idx = 0; gf_idx < NUM_GFS_NEAREST; gf_idx++) {
          const int gf = which_gfs[gf_idx];
          row[2 + gf_idx] = src[IDX4Ppt(params_chare, gf, idx3)];
        } // END LOOP over gridfunctions

        int n = 0;
        n += sprintf(out + n, "% .15e", row[0]);
        for (int col = 1; col < NUM_COLS; col++) {
          n += sprintf(out + n, " % .15e", row[col]);
        }
        out[sizeinbytes - 1] = '\n';
        Ck::IO::write(token, out, sizeinbytes, offsetpt_firstfield[which_pt]);
      }
      // Finalize
      free(row);
      break;
    }
  }




  // Open output files (one file per timestep, per plane, for this grid)
  //~ FILE *out_xy = open_outfile("out2d-xy", coordsys_with_grid, commondata, /*include_time=*/1);
  //~ FILE *out_yz = open_outfile("out2d-yz", coordsys_with_grid, commondata, /*include_time=*/1);

  //~ if (!out_xy || !out_yz) {
    //~ if (out_xy)
      //~ fclose(out_xy);
    //~ if (out_yz)
      //~ fclose(out_yz);
    //~ fprintf(stderr, "Error: Cannot open output files for grid %d.\n", grid);
    //~ exit(1);
  //~ } // END IF cannot open output files

  //~ // Write time comment and headers
  //~ diag_write_time_comment(out_xy, commondata->time);
  //~ diag_write_header(out_xy, "x y", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);
  //~ diag_write_time_comment(out_yz, commondata->time);
  //~ diag_write_header(out_yz, "y z", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);


} // END FUNCTION
