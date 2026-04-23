#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Initial data for GRoovy shock tests
 */
void initial_data(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    // Unpack the current grid and its gridfunctions.
    params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
    REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++) {
      xx[ww] = griddata[grid].xx[ww];
    } // END LOOP: for ww over coordinate directions
    REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
    REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
#pragma omp parallel for
    for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
      MAYBE_UNUSED const REAL xx2 = xx[2][i2];
      for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
        MAYBE_UNUSED const REAL xx1 = xx[1][i1];
        for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
          MAYBE_UNUSED const REAL xx0 = xx[0][i0];

          /*
           *  Original SymPy expressions:
           *  "[auxevol_gfs[IDX4(RHOBGF, i0, i1, i2)] = (TINYDOUBLE/2 + xx0*cos(xx1)/2 + nrpyAbs(TINYDOUBLE + xx0*cos(xx1) - 1)/2 -
           * 1/2)/(8*(TINYDOUBLE + xx0*cos(xx1) - 1)) + (xx0*cos(xx1)/2 - nrpyAbs(xx0*cos(xx1) - 1)/2 - 1/2)/(-TINYDOUBLE + xx0*cos(xx1) - 1)]"
           *  "[auxevol_gfs[IDX4(PGF, i0, i1, i2)] = (TINYDOUBLE/2 + xx0*cos(xx1)/2 + nrpyAbs(TINYDOUBLE + xx0*cos(xx1) - 1)/2 - 1/2)/(10*(TINYDOUBLE
           * + xx0*cos(xx1) - 1)) + (xx0*cos(xx1)/2 - nrpyAbs(xx0*cos(xx1) - 1)/2 - 1/2)/(-TINYDOUBLE + xx0*cos(xx1) - 1)]"
           *  "[auxevol_gfs[IDX4(RESCALEDVU0GF, i0, i1, i2)] = 0]"
           *  "[auxevol_gfs[IDX4(RESCALEDVU1GF, i0, i1, i2)] = 0]"
           *  "[auxevol_gfs[IDX4(RESCALEDVU2GF, i0, i1, i2)] = 0]"
           *  "[auxevol_gfs[IDX4(HDD00GF, i0, i1, i2)] = (xx0**4*sin(xx1)**2/((xx0**2*sin(xx1)**2*sin(xx2)**2 +
           * xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0**2*sin(xx1)**2 + xx0**2*sin(xx2)**2*cos(xx1)**2 +
           * xx0**2*cos(xx1)**2*cos(xx2)**2)*(sin(xx1)**2*sin(xx2)**2 + sin(xx1)**2*cos(xx2)**2 + cos(xx1)**2) - (xx0**2*sin(xx1)**2*sin(xx2)**2 +
           * xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0*sin(xx1)*cos(xx1)*cos(xx2)**2 -
           * xx0*sin(xx1)*cos(xx1))**2))**(1/3)*(sin(xx1)**2*sin(xx2)**2 + sin(xx1)**2*cos(xx2)**2 + cos(xx1)**2) - 1]"
           *  "[auxevol_gfs[IDX4(HDD01GF, i0, i1, i2)] = (xx0**4*sin(xx1)**2/((xx0**2*sin(xx1)**2*sin(xx2)**2 +
           * xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0**2*sin(xx1)**2 + xx0**2*sin(xx2)**2*cos(xx1)**2 +
           * xx0**2*cos(xx1)**2*cos(xx2)**2)*(sin(xx1)**2*sin(xx2)**2 + sin(xx1)**2*cos(xx2)**2 + cos(xx1)**2) - (xx0**2*sin(xx1)**2*sin(xx2)**2 +
           * xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0*sin(xx1)*cos(xx1)*cos(xx2)**2 -
           * xx0*sin(xx1)*cos(xx1))**2))**(1/3)*(xx0*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0*sin(xx1)*cos(xx1)*cos(xx2)**2 - xx0*sin(xx1)*cos(xx1))/xx0]"
           *  "[auxevol_gfs[IDX4(HDD02GF, i0, i1, i2)] = 0]"
           *  "[auxevol_gfs[IDX4(HDD11GF, i0, i1, i2)] = (-xx0**2 + (xx0**4*sin(xx1)**2/((xx0**2*sin(xx1)**2*sin(xx2)**2 +
           * xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0**2*sin(xx1)**2 + xx0**2*sin(xx2)**2*cos(xx1)**2 +
           * xx0**2*cos(xx1)**2*cos(xx2)**2)*(sin(xx1)**2*sin(xx2)**2 + sin(xx1)**2*cos(xx2)**2 + cos(xx1)**2) - (xx0**2*sin(xx1)**2*sin(xx2)**2 +
           * xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0*sin(xx1)*cos(xx1)*cos(xx2)**2 -
           * xx0*sin(xx1)*cos(xx1))**2))**(1/3)*(xx0**2*sin(xx1)**2 + xx0**2*sin(xx2)**2*cos(xx1)**2 + xx0**2*cos(xx1)**2*cos(xx2)**2))/xx0**2]"
           *  "[auxevol_gfs[IDX4(HDD12GF, i0, i1, i2)] = 0]"
           *  "[auxevol_gfs[IDX4(HDD22GF, i0, i1, i2)] = (-xx0**2*sin(xx1)**2 + (xx0**4*sin(xx1)**2/((xx0**2*sin(xx1)**2*sin(xx2)**2 +
           * xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0**2*sin(xx1)**2 + xx0**2*sin(xx2)**2*cos(xx1)**2 +
           * xx0**2*cos(xx1)**2*cos(xx2)**2)*(sin(xx1)**2*sin(xx2)**2 + sin(xx1)**2*cos(xx2)**2 + cos(xx1)**2) - (xx0**2*sin(xx1)**2*sin(xx2)**2 +
           * xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0*sin(xx1)*cos(xx1)*cos(xx2)**2 -
           * xx0*sin(xx1)*cos(xx1))**2))**(1/3)*(xx0**2*sin(xx1)**2*sin(xx2)**2 + xx0**2*sin(xx1)**2*cos(xx2)**2))/(xx0**2*sin(xx1)**2)]"
           *  "[auxevol_gfs[IDX4(ADD00GF, i0, i1, i2)] = 0]"
           *  "[auxevol_gfs[IDX4(ADD01GF, i0, i1, i2)] = 0]"
           *  "[auxevol_gfs[IDX4(ADD02GF, i0, i1, i2)] = 0]"
           *  "[auxevol_gfs[IDX4(ADD11GF, i0, i1, i2)] = 0]"
           *  "[auxevol_gfs[IDX4(ADD12GF, i0, i1, i2)] = 0]"
           *  "[auxevol_gfs[IDX4(ADD22GF, i0, i1, i2)] = 0]"
           *  "[auxevol_gfs[IDX4(VETU0GF, i0, i1, i2)] = 0]"
           *  "[auxevol_gfs[IDX4(VETU1GF, i0, i1, i2)] = 0]"
           *  "[auxevol_gfs[IDX4(VETU2GF, i0, i1, i2)] = 0]"
           *  "[auxevol_gfs[IDX4(ALPHAGF, i0, i1, i2)] = 1]"
           *  "[auxevol_gfs[IDX4(CFGF, i0, i1, i2)] = (((xx0**2*sin(xx1)**2*sin(xx2)**2 + xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0**2*sin(xx1)**2 +
           * xx0**2*sin(xx2)**2*cos(xx1)**2 + xx0**2*cos(xx1)**2*cos(xx2)**2)*(sin(xx1)**2*sin(xx2)**2 + sin(xx1)**2*cos(xx2)**2 + cos(xx1)**2) -
           * (xx0**2*sin(xx1)**2*sin(xx2)**2 + xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0*sin(xx1)*cos(xx1)*cos(xx2)**2
           * - xx0*sin(xx1)*cos(xx1))**2)/((xx0**4*sin(xx1)**2/((xx0**2*sin(xx1)**2*sin(xx2)**2 + xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0**2*sin(xx1)**2
           * + xx0**2*sin(xx2)**2*cos(xx1)**2 + xx0**2*cos(xx1)**2*cos(xx2)**2)*(sin(xx1)**2*sin(xx2)**2 + sin(xx1)**2*cos(xx2)**2 + cos(xx1)**2) -
           * (xx0**2*sin(xx1)**2*sin(xx2)**2 + xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0*sin(xx1)*cos(xx1)*cos(xx2)**2
           * - xx0*sin(xx1)*cos(xx1))**2))*(xx0**2*sin(xx1)**2*sin(xx2)**2 + xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0**2*sin(xx1)**2 +
           * xx0**2*sin(xx2)**2*cos(xx1)**2 + xx0**2*cos(xx1)**2*cos(xx2)**2)*(sin(xx1)**2*sin(xx2)**2 + sin(xx1)**2*cos(xx2)**2 + cos(xx1)**2) -
           * xx0**4*sin(xx1)**2/((xx0**2*sin(xx1)**2*sin(xx2)**2 + xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0**2*sin(xx1)**2 +
           * xx0**2*sin(xx2)**2*cos(xx1)**2 + xx0**2*cos(xx1)**2*cos(xx2)**2)*(sin(xx1)**2*sin(xx2)**2 + sin(xx1)**2*cos(xx2)**2 + cos(xx1)**2) -
           * (xx0**2*sin(xx1)**2*sin(xx2)**2 + xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0*sin(xx1)*cos(xx1)*cos(xx2)**2
           * - xx0*sin(xx1)*cos(xx1))**2)*(xx0**2*sin(xx1)**2*sin(xx2)**2 + xx0**2*sin(xx1)**2*cos(xx2)**2)*(xx0*sin(xx1)*sin(xx2)**2*cos(xx1) +
           * xx0*sin(xx1)*cos(xx1)*cos(xx2)**2 - xx0*sin(xx1)*cos(xx1))**2))**(-1/6)]"
           *  "[auxevol_gfs[IDX4(TRKGF, i0, i1, i2)] = 0]"
           */
          {
            const REAL tmp0 = cos(xx1);
            const REAL tmp7 = sin(xx1);
            const REAL tmp9 = ((sin(xx2)) * (sin(xx2)));
            const REAL tmp11 = ((cos(xx2)) * (cos(xx2)));
            const REAL tmp14 = ((xx0) * (xx0));
            const REAL tmp1 = tmp0 * xx0;
            const REAL tmp8 = ((tmp7) * (tmp7));
            const REAL tmp25 = (1.0 / (tmp14));
            const REAL tmp16 = tmp1 * tmp7;
            const REAL tmp20 = ((tmp0) * (tmp0)) * tmp14;
            const REAL tmp3 = ((1.0 / 2.0) * tmp0 * xx0 - 1.0 / 2.0 * fabs(tmp1 - 1) - 1.0 / 2.0) / (-TINYDOUBLE + tmp0 * xx0 - 1);
            const REAL tmp4 = TINYDOUBLE + tmp1 - 1;
            const REAL tmp13 = ((tmp0) * (tmp0)) + tmp11 * tmp8 + tmp8 * tmp9;
            const REAL tmp15 = tmp11 * tmp14 * tmp8 + tmp14 * tmp8 * tmp9;
            const REAL tmp17 = tmp11 * tmp16 + tmp16 * tmp9 - tmp16;
            const REAL tmp21 = tmp11 * tmp20 + tmp14 * tmp8 + tmp20 * tmp9;
            const REAL tmp5 = ((1.0 / 2.0) * TINYDOUBLE + (1.0 / 2.0) * tmp1 + (1.0 / 2.0) * fabs(tmp4) - 1.0 / 2.0) / tmp4;
            const REAL tmp22 = tmp13 * tmp15 * tmp21 - tmp15 * ((tmp17) * (tmp17));
            const REAL tmp23 = tmp8 * ((xx0) * (xx0) * (xx0) * (xx0)) / tmp22;
            const REAL tmp24 = cbrt(tmp23);
            auxevol_gfs[IDX4(RHOBGF, i0, i1, i2)] = tmp3 + (1.0 / 8.0) * tmp5;
            auxevol_gfs[IDX4(PGF, i0, i1, i2)] = tmp3 + (1.0 / 10.0) * tmp5;
            auxevol_gfs[IDX4(RESCALEDVU0GF, i0, i1, i2)] = 0;
            auxevol_gfs[IDX4(RESCALEDVU1GF, i0, i1, i2)] = 0;
            auxevol_gfs[IDX4(RESCALEDVU2GF, i0, i1, i2)] = 0;
            auxevol_gfs[IDX4(HDD00GF, i0, i1, i2)] = tmp13 * tmp24 - 1;
            auxevol_gfs[IDX4(HDD01GF, i0, i1, i2)] = tmp17 * tmp24 / xx0;
            auxevol_gfs[IDX4(HDD02GF, i0, i1, i2)] = 0;
            auxevol_gfs[IDX4(HDD11GF, i0, i1, i2)] = tmp25 * (-tmp14 + tmp21 * tmp24);
            auxevol_gfs[IDX4(HDD12GF, i0, i1, i2)] = 0;
            auxevol_gfs[IDX4(HDD22GF, i0, i1, i2)] = tmp25 * (-tmp14 * tmp8 + tmp15 * tmp24) / tmp8;
            auxevol_gfs[IDX4(ADD00GF, i0, i1, i2)] = 0;
            auxevol_gfs[IDX4(ADD01GF, i0, i1, i2)] = 0;
            auxevol_gfs[IDX4(ADD02GF, i0, i1, i2)] = 0;
            auxevol_gfs[IDX4(ADD11GF, i0, i1, i2)] = 0;
            auxevol_gfs[IDX4(ADD12GF, i0, i1, i2)] = 0;
            auxevol_gfs[IDX4(ADD22GF, i0, i1, i2)] = 0;
            auxevol_gfs[IDX4(VETU0GF, i0, i1, i2)] = 0;
            auxevol_gfs[IDX4(VETU1GF, i0, i1, i2)] = 0;
            auxevol_gfs[IDX4(VETU2GF, i0, i1, i2)] = 0;
            auxevol_gfs[IDX4(ALPHAGF, i0, i1, i2)] = 1;
            auxevol_gfs[IDX4(CFGF, i0, i1, i2)] = (1.0 / (sqrt(cbrt(tmp22 / (tmp13 * tmp15 * tmp21 * tmp23 - tmp15 * ((tmp17) * (tmp17)) * tmp23)))));
            auxevol_gfs[IDX4(TRKGF, i0, i1, i2)] = 0;
          }

          ghl_primitive_quantities prims;
          ghl_metric_quantities metric;
          ghl_conservative_quantities cons;

          // Rebuild GRHayL structs from the freshly seeded gridfunction data.
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
} // END FUNCTION initial_data
