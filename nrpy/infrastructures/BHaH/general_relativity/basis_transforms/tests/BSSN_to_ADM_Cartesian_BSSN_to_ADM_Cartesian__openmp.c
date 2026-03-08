#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

#define ADMALPHA 0
#define ADMBETAX 1
#define ADMBETAY 2
#define ADMBETAZ 3
#define ADMGXX 4
#define ADMGXY 5
#define ADMGXZ 6
#define ADMGYY 7
#define ADMGYZ 8
#define ADMGZZ 9
#define ADMKXX 10
#define ADMKXY 11
#define ADMKXZ 12
#define ADMKYY 13
#define ADMKYZ 14
#define ADMKZZ 15

/**
 * @brief Convert Cartesian-basis BSSN gridfunctions to Cartesian-basis ADM gridfunctions.
 *
 * The input is read from ``in_gfs`` and the ADM output is written to ``out_gfs``
 * at matching grid points using the requested ADM slot names.
 *
 * @param[in] params Grid parameters used for loop extents.
 * @param[in] in_gfs Cartesian-basis BSSN input gridfunctions.
 * @param[out] out_gfs Cartesian-basis ADM output gridfunctions.
 */
void BSSN_to_ADM_Cartesian(const params_struct *restrict params, const REAL *restrict in_gfs, REAL *restrict out_gfs) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

#pragma omp parallel for
  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
      for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {

        const REAL alpha = in_gfs[IDX4(ALPHAGF, i0, i1, i2)];
        const REAL cf = in_gfs[IDX4(CFGF, i0, i1, i2)];
        const REAL trK = in_gfs[IDX4(TRKGF, i0, i1, i2)];
        const REAL vetU0 = in_gfs[IDX4(VETU0GF, i0, i1, i2)];
        const REAL vetU1 = in_gfs[IDX4(VETU1GF, i0, i1, i2)];
        const REAL vetU2 = in_gfs[IDX4(VETU2GF, i0, i1, i2)];
        const REAL hDD00 = in_gfs[IDX4(HDD00GF, i0, i1, i2)];
        const REAL hDD01 = in_gfs[IDX4(HDD01GF, i0, i1, i2)];
        const REAL hDD02 = in_gfs[IDX4(HDD02GF, i0, i1, i2)];
        const REAL hDD11 = in_gfs[IDX4(HDD11GF, i0, i1, i2)];
        const REAL hDD12 = in_gfs[IDX4(HDD12GF, i0, i1, i2)];
        const REAL hDD22 = in_gfs[IDX4(HDD22GF, i0, i1, i2)];
        const REAL aDD00 = in_gfs[IDX4(ADD00GF, i0, i1, i2)];
        const REAL aDD01 = in_gfs[IDX4(ADD01GF, i0, i1, i2)];
        const REAL aDD02 = in_gfs[IDX4(ADD02GF, i0, i1, i2)];
        const REAL aDD11 = in_gfs[IDX4(ADD11GF, i0, i1, i2)];
        const REAL aDD12 = in_gfs[IDX4(ADD12GF, i0, i1, i2)];
        const REAL aDD22 = in_gfs[IDX4(ADD22GF, i0, i1, i2)];
        const REAL tmp0 = (1.0 / ((cf) * (cf)));
        const REAL tmp7 = (1.0 / 3.0) * trK;
        const REAL tmp1 = tmp0 * (hDD00 + 1);
        const REAL tmp4 = tmp0 * (hDD11 + 1);
        const REAL tmp6 = tmp0 * (hDD22 + 1);
        out_gfs[IDX4(ADMALPHA, i0, i1, i2)] = alpha;
        out_gfs[IDX4(ADMBETAX, i0, i1, i2)] = vetU0;
        out_gfs[IDX4(ADMBETAY, i0, i1, i2)] = vetU1;
        out_gfs[IDX4(ADMBETAZ, i0, i1, i2)] = vetU2;
        out_gfs[IDX4(ADMGXX, i0, i1, i2)] = tmp1;
        out_gfs[IDX4(ADMGXY, i0, i1, i2)] = hDD01 * tmp0;
        out_gfs[IDX4(ADMGXZ, i0, i1, i2)] = hDD02 * tmp0;
        out_gfs[IDX4(ADMGYY, i0, i1, i2)] = tmp4;
        out_gfs[IDX4(ADMGYZ, i0, i1, i2)] = hDD12 * tmp0;
        out_gfs[IDX4(ADMGZZ, i0, i1, i2)] = tmp6;
        out_gfs[IDX4(ADMKXX, i0, i1, i2)] = aDD00 * tmp0 + tmp1 * tmp7;
        out_gfs[IDX4(ADMKXY, i0, i1, i2)] = aDD01 * tmp0 + hDD01 * tmp0 * tmp7;
        out_gfs[IDX4(ADMKXZ, i0, i1, i2)] = aDD02 * tmp0 + hDD02 * tmp0 * tmp7;
        out_gfs[IDX4(ADMKYY, i0, i1, i2)] = aDD11 * tmp0 + tmp4 * tmp7;
        out_gfs[IDX4(ADMKYZ, i0, i1, i2)] = aDD12 * tmp0 + hDD12 * tmp0 * tmp7;
        out_gfs[IDX4(ADMKZZ, i0, i1, i2)] = aDD22 * tmp0 + tmp6 * tmp7;

      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
} // END FUNCTION BSSN_to_ADM_Cartesian
