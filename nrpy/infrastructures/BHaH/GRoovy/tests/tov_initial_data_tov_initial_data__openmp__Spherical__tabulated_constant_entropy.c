#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Read a fixed-entropy, three-column EOS slice from disk.
 *
 * @param[in] filename Path to the slice file.
 * @param[in] num_rows Number of rows expected in the file.
 * @return Dynamically allocated array with shape [num_rows][3], or NULL on failure.
 */
static double (*eos_slice_file(const char *filename, int num_rows))[3] {
  FILE *file_ptr = NULL;
  double (*data_array)[3] = malloc(num_rows * sizeof(double[3]));
  char line_buffer[256];

  if (num_rows <= 0) {
    fprintf(stderr, "Error: Number of rows must be positive.\n");
    return NULL;
  } // END IF: invalid row count
  if (data_array == NULL) {
    fprintf(stderr, "Error: Could not allocate memory for %d rows.\n", num_rows);
    return NULL;
  } // END IF: allocation failure

  file_ptr = fopen(filename, "r");
  if (file_ptr == NULL) {
    fprintf(stderr, "Error opening file '%s': %s\n", filename, strerror(errno));
    free(data_array);
    return NULL;
  } // END IF: failed file open

  for (int i = 0; i < num_rows; ++i) {
    if (fgets(line_buffer, sizeof(line_buffer), file_ptr) == NULL) {
      fprintf(stderr, "Error: Unexpected EOF or read error at row %d.\n", i + 1);
      free(data_array);
      fclose(file_ptr);
      return NULL;
    } // END IF: failed line read

    if (sscanf(line_buffer, "%lf %lf %lf", &data_array[i][0], &data_array[i][1], &data_array[i][2]) != 3) {
      fprintf(stderr, "Error: Could not parse 3 doubles on line %d.\n", i + 1);
      free(data_array);
      fclose(file_ptr);
      return NULL;
    } // END IF: failed line parse
  } // END LOOP: for i over EOS slice rows

  fclose(file_ptr);
  return data_array;
} // END FUNCTION: eos_slice_file

/**
 * Spherically symmetric TOV initial data with tabulated EOS data at fixed entropy
 */
void initial_data(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  const char *eos_slice_filename = "constant_ent_Beq_slice.txt";
  double (*constant_ent_Beq_slice)[3] = eos_slice_file(eos_slice_filename, commondata->eos.N_rho);

  commondata->eos.lp_of_lr = (double *)malloc(sizeof(double) * commondata->eos.N_rho);
  commondata->eos.le_of_lr = (double *)malloc(sizeof(double) * commondata->eos.N_rho);
  commondata->eos.Ye_of_lr = (double *)malloc(sizeof(double) * commondata->eos.N_rho);

  for (int i = 0; i < commondata->eos.N_rho; ++i) {
    commondata->eos.Ye_of_lr[i] = constant_ent_Beq_slice[i][0];
    commondata->eos.lp_of_lr[i] = constant_ent_Beq_slice[i][1];
    commondata->eos.le_of_lr[i] = constant_ent_Beq_slice[i][2];
  } // END LOOP: for i over EOS entropy-slice rows

  free(constant_ent_Beq_slice);
  const double beq_entropy = 0.3;

  ID_persist_struct ID_persist;
  TOV_read_data_file_set_ID_persist("TOVdata.txt", &ID_persist);

  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    // Unpack the current grid and its gridfunctions.
    params_struct *restrict params = &griddata[grid].params;
    bc_struct *restrict bcstruct = &griddata[grid].bcstruct;
#include "set_CodeParameters.h"
    REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++) {
      xx[ww] = griddata[grid].xx[ww];
    } // END LOOP: for ww over coordinate directions
    REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
    REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;

    // Import the TOV spacetime data into the BSSN gridfunctions.
    initial_data_reader__convert_ADM_Spherical_to_BSSN(commondata, params, xx, bcstruct, &griddata[grid].gridfuncs, &ID_persist, TOV_ID_function);
#pragma omp parallel for
    for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
      MAYBE_UNUSED const REAL xx2 = xx[2][i2];
      for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
        MAYBE_UNUSED const REAL xx1 = xx[1][i1];
        for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
          MAYBE_UNUSED const REAL xx0 = xx[0][i0];

          // Interpolate the TOV matter profile at this grid point.
          REAL rbar, rho, rho_baryon, press, M, expnu, exp4phi;
          /*
           *  Original SymPy expression:
           *  "rbar = xx0"
           */
          rbar = xx0;

          TOV_interpolate_1D(rbar, &ID_persist, &rho, &rho_baryon, &press, &M, &expnu, &exp4phi);

          // Seed the tabulated primitive variables at fixed entropy.
          ghl_primitive_quantities prims;
          ghl_metric_quantities metric;
          ghl_conservative_quantities cons;
          double dummy0, dummy1, dummy2, dummy3, dummy4;

          prims.rho = rho_baryon;
          prims.press = press;
          prims.vU[0] = 0.0;
          prims.vU[1] = 0.0;
          prims.vU[2] = 0.0;

          const double rhoL = prims.rho;
          if (rhoL <= 1.01 * commondata->eos.rho_atm) {
            prims.rho = commondata->eos.rho_atm;
            prims.press = commondata->eos.press_atm;
            prims.eps = commondata->eos.eps_atm;
            prims.Y_e = commondata->eos.Y_e_atm;
            prims.temperature = commondata->eos.T_atm;
          } else {
            prims.Y_e = ghl_tabulated_compute_Ye_from_rho(&commondata->eos, rhoL);

            // Supply a temperature guess before inverting for the fixed-entropy state.
            prims.temperature = 1.0;
            ghl_tabulated_compute_P_T_from_S(&commondata->eos, prims.rho, prims.Y_e, beq_entropy, &prims.press, &prims.temperature);
          } // END IF: atmosphere reset versus tabulated-entropy seed

          ghl_return_primitives(&prims, &auxevol_gfs[IDX4(RHOBGF, i0, i1, i2)], &auxevol_gfs[IDX4(PGF, i0, i1, i2)], &dummy0,
                                &auxevol_gfs[IDX4(RESCALEDVU0GF, i0, i1, i2)], &auxevol_gfs[IDX4(RESCALEDVU1GF, i0, i1, i2)],
                                &auxevol_gfs[IDX4(RESCALEDVU2GF, i0, i1, i2)], &dummy1, &dummy2, &dummy3, &dummy4,
                                &auxevol_gfs[IDX4(YEGF, i0, i1, i2)], &auxevol_gfs[IDX4(TEMPERATUREGF, i0, i1, i2)]);

          // Rebuild GRHayL structs from the TOV-seeded gridfunction data.
          basis_transform_rfm_basis_to_Cartesian(commondata, params, &prims, &cons, &metric, i0, i1, i2, xx, auxevol_gfs, in_gfs);

          bool speed_limited;
          ghl_enforce_primitive_limits_and_compute_u0(&commondata->ghl_params, &commondata->eos, &metric, &prims, &speed_limited);

          // Store the limited primitive variables and reconstructed conservatives.
          basis_transform_Cartesian_to_rfm_basis(commondata, params, &prims, &cons, i0, i1, i2, xx, auxevol_gfs, in_gfs);

          auxevol_gfs[IDX4(U4UTGF, i0, i1, i2)] = prims.u0;

        } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
      } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
    } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)

  } // END LOOP: for grid over grids

  // Release the interpolated TOV profile storage after initial-data setup.
  free(ID_persist.r_Schw_arr);
  free(ID_persist.rho_arr);
  free(ID_persist.rho_baryon_arr);
  free(ID_persist.P_arr);
  free(ID_persist.M_arr);
  free(ID_persist.expnu_arr);
  free(ID_persist.exp4phi_arr);
  free(ID_persist.rbar_arr);

  free(commondata->eos.Ye_of_lr);
  commondata->eos.Ye_of_lr = NULL;
  free(commondata->eos.lp_of_lr);
  commondata->eos.lp_of_lr = NULL;
  free(commondata->eos.le_of_lr);
  commondata->eos.le_of_lr = NULL;
} // END FUNCTION initial_data
