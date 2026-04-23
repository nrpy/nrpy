#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Spherically symmetric TOV initial data with tabulated EOS data at fixed temperature
 */
void initial_data(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  const double beq_temperature = 0.3;
  ghl_tabulated_compute_Ye_of_rho_beq_constant_T(beq_temperature, &commondata->eos);

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

          // Seed the tabulated primitive variables at fixed temperature.
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
            const double tempL = beq_temperature;
            const double YeL = ghl_tabulated_compute_Ye_from_rho(&commondata->eos, rhoL);

            double pressL, epsL;
            ghl_tabulated_compute_P_eps_from_T(&commondata->eos, rhoL, YeL, tempL, &pressL, &epsL);

            prims.press = pressL;
            prims.eps = epsL;
            prims.Y_e = YeL;
            prims.temperature = tempL;
          } // END IF: atmosphere reset versus tabulated-temperature seed

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
} // END FUNCTION initial_data
