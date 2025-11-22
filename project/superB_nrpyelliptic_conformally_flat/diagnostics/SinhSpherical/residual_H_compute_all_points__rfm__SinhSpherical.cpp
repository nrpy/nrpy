#include "BHaH_defines.h"

/**
 * @file residual_H_compute_all_points.c
 * @brief Compute the Hamiltonian-constraint residual on all interior grid points and store the result.
 *
 * The generated function "residual_H_compute_all_points" iterates over the interior region of the grid
 * for the selected coordinate system and evaluates the Hamiltonian-constraint residual used by the
 * hyperbolic relaxation scheme. Results are written into the destination gridfunction buffer.
 *
 * Reference-metric precomputation is currently disabled at code-generation time, so this routine
 * accepts xx coordinate arrays and does not take an rfm_struct parameter.
 *
 * If a user-editable block is present in the implementation, users may insert custom logic such as
 * additional diagnostics or instrumentation without changing the function interface.
 *
 * @param[in]  commondata        Pointer to read-only global simulation metadata (e.g., time, step counters); may be unused.
 * @param[in]  params            Pointer to read-only per-grid parameters (sizes, ghost zones, strides, names).
 * @param[in]  xx                Array of three coordinate arrays used for coordinate-dependent operations.
 * @param[in]  auxevol_gfs       Pointer to read-only auxiliary evolution gridfunctions required by the residual.
 * @param[in]  in_gfs            Pointer to read-only input gridfunctions (e.g., current solution fields).
 * @param[out] dest_gf_address   Pointer to the destination gridfunction buffer where the residual is stored.
 *
 * @return     void.
 */
void residual_H_compute_all_points__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                       REAL *restrict xx[3], const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs,
                                                       REAL *restrict dest_gf_address) {
#include "set_CodeParameters.h"
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    MAYBE_UNUSED const REAL xx2 = xx[2][i2];
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      MAYBE_UNUSED const REAL xx1 = xx[1][i1];
      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        MAYBE_UNUSED const REAL xx0 = xx[0][i0];

        /*
         * NRPy-Generated GF Access/FD Code, Step 1 of 2:
         * Read gridfunction(s) from main memory and compute FD stencils as needed.
         */
        const REAL ADD_times_AUU = auxevol_gfs[IDX4(ADD_TIMES_AUUGF, i0, i1, i2)];
        const REAL psi_background = auxevol_gfs[IDX4(PSI_BACKGROUNDGF, i0, i1, i2)];
        const REAL uu_i2m5 = in_gfs[IDX4(UUGF, i0, i1, i2 - 5)];
        const REAL uu_i2m4 = in_gfs[IDX4(UUGF, i0, i1, i2 - 4)];
        const REAL uu_i2m3 = in_gfs[IDX4(UUGF, i0, i1, i2 - 3)];
        const REAL uu_i2m2 = in_gfs[IDX4(UUGF, i0, i1, i2 - 2)];
        const REAL uu_i2m1 = in_gfs[IDX4(UUGF, i0, i1, i2 - 1)];
        const REAL uu_i1m5 = in_gfs[IDX4(UUGF, i0, i1 - 5, i2)];
        const REAL uu_i1m4 = in_gfs[IDX4(UUGF, i0, i1 - 4, i2)];
        const REAL uu_i1m3 = in_gfs[IDX4(UUGF, i0, i1 - 3, i2)];
        const REAL uu_i1m2 = in_gfs[IDX4(UUGF, i0, i1 - 2, i2)];
        const REAL uu_i1m1 = in_gfs[IDX4(UUGF, i0, i1 - 1, i2)];
        const REAL uu_i0m5 = in_gfs[IDX4(UUGF, i0 - 5, i1, i2)];
        const REAL uu_i0m4 = in_gfs[IDX4(UUGF, i0 - 4, i1, i2)];
        const REAL uu_i0m3 = in_gfs[IDX4(UUGF, i0 - 3, i1, i2)];
        const REAL uu_i0m2 = in_gfs[IDX4(UUGF, i0 - 2, i1, i2)];
        const REAL uu_i0m1 = in_gfs[IDX4(UUGF, i0 - 1, i1, i2)];
        const REAL uu = in_gfs[IDX4(UUGF, i0, i1, i2)];
        const REAL uu_i0p1 = in_gfs[IDX4(UUGF, i0 + 1, i1, i2)];
        const REAL uu_i0p2 = in_gfs[IDX4(UUGF, i0 + 2, i1, i2)];
        const REAL uu_i0p3 = in_gfs[IDX4(UUGF, i0 + 3, i1, i2)];
        const REAL uu_i0p4 = in_gfs[IDX4(UUGF, i0 + 4, i1, i2)];
        const REAL uu_i0p5 = in_gfs[IDX4(UUGF, i0 + 5, i1, i2)];
        const REAL uu_i1p1 = in_gfs[IDX4(UUGF, i0, i1 + 1, i2)];
        const REAL uu_i1p2 = in_gfs[IDX4(UUGF, i0, i1 + 2, i2)];
        const REAL uu_i1p3 = in_gfs[IDX4(UUGF, i0, i1 + 3, i2)];
        const REAL uu_i1p4 = in_gfs[IDX4(UUGF, i0, i1 + 4, i2)];
        const REAL uu_i1p5 = in_gfs[IDX4(UUGF, i0, i1 + 5, i2)];
        const REAL uu_i2p1 = in_gfs[IDX4(UUGF, i0, i1, i2 + 1)];
        const REAL uu_i2p2 = in_gfs[IDX4(UUGF, i0, i1, i2 + 2)];
        const REAL uu_i2p3 = in_gfs[IDX4(UUGF, i0, i1, i2 + 3)];
        const REAL uu_i2p4 = in_gfs[IDX4(UUGF, i0, i1, i2 + 4)];
        const REAL uu_i2p5 = in_gfs[IDX4(UUGF, i0, i1, i2 + 5)];
        static const REAL FDPart1_Rational_5_6 = 5.0 / 6.0;
        static const REAL FDPart1_Rational_5_21 = 5.0 / 21.0;
        static const REAL FDPart1_Rational_5_84 = 5.0 / 84.0;
        static const REAL FDPart1_Rational_5_504 = 5.0 / 504.0;
        static const REAL FDPart1_Rational_1_1260 = 1.0 / 1260.0;
        static const REAL FDPart1_Rational_5269_1800 = 5269.0 / 1800.0;
        static const REAL FDPart1_Rational_5_1008 = 5.0 / 1008.0;
        static const REAL FDPart1_Rational_1_3150 = 1.0 / 3150.0;
        static const REAL FDPart1_Rational_5_3 = 5.0 / 3.0;
        static const REAL FDPart1_Rational_5_126 = 5.0 / 126.0;
        const REAL FDPart1tmp0 = -FDPart1_Rational_5269_1800 * uu;
        const REAL uu_dD0 = invdxx0 * (FDPart1_Rational_1_1260 * (-uu_i0m5 + uu_i0p5) + FDPart1_Rational_5_21 * (uu_i0m2 - uu_i0p2) +
                                       FDPart1_Rational_5_504 * (uu_i0m4 - uu_i0p4) + FDPart1_Rational_5_6 * (-uu_i0m1 + uu_i0p1) +
                                       FDPart1_Rational_5_84 * (-uu_i0m3 + uu_i0p3));
        const REAL uu_dD1 = invdxx1 * (FDPart1_Rational_1_1260 * (-uu_i1m5 + uu_i1p5) + FDPart1_Rational_5_21 * (uu_i1m2 - uu_i1p2) +
                                       FDPart1_Rational_5_504 * (uu_i1m4 - uu_i1p4) + FDPart1_Rational_5_6 * (-uu_i1m1 + uu_i1p1) +
                                       FDPart1_Rational_5_84 * (-uu_i1m3 + uu_i1p3));
        const REAL uu_dDD00 =
            ((invdxx0) * (invdxx0)) * (FDPart1_Rational_1_3150 * (uu_i0m5 + uu_i0p5) + FDPart1_Rational_5_1008 * (-uu_i0m4 - uu_i0p4) +
                                       FDPart1_Rational_5_126 * (uu_i0m3 + uu_i0p3) + FDPart1_Rational_5_21 * (-uu_i0m2 - uu_i0p2) +
                                       FDPart1_Rational_5_3 * (uu_i0m1 + uu_i0p1) + FDPart1tmp0);
        const REAL uu_dDD11 =
            ((invdxx1) * (invdxx1)) * (FDPart1_Rational_1_3150 * (uu_i1m5 + uu_i1p5) + FDPart1_Rational_5_1008 * (-uu_i1m4 - uu_i1p4) +
                                       FDPart1_Rational_5_126 * (uu_i1m3 + uu_i1p3) + FDPart1_Rational_5_21 * (-uu_i1m2 - uu_i1p2) +
                                       FDPart1_Rational_5_3 * (uu_i1m1 + uu_i1p1) + FDPart1tmp0);
        const REAL uu_dDD22 =
            ((invdxx2) * (invdxx2)) * (FDPart1_Rational_1_3150 * (uu_i2m5 + uu_i2p5) + FDPart1_Rational_5_1008 * (-uu_i2m4 - uu_i2p4) +
                                       FDPart1_Rational_5_126 * (uu_i2m3 + uu_i2p3) + FDPart1_Rational_5_21 * (-uu_i2m2 - uu_i2p2) +
                                       FDPart1_Rational_5_3 * (uu_i2m1 + uu_i2p1) + FDPart1tmp0);

        /*
         * NRPy-Generated GF Access/FD Code, Step 2 of 2:
         * Evaluate SymPy expressions and write to main memory.
         */
        const REAL FDPart3tmp0 = (1.0 / (SINHW));
        const REAL FDPart3tmp7 = sin(xx1);
        const REAL FDPart3tmp12 = 2 / ((SINHW) * (SINHW));
        const REAL FDPart3tmp5 = ((exp(FDPart3tmp0) - exp(-FDPart3tmp0)) * (exp(FDPart3tmp0) - exp(-FDPart3tmp0))) / ((AMPL) * (AMPL));
        const REAL FDPart3tmp2 = exp(FDPart3tmp0 * xx0);
        const REAL FDPart3tmp3 = exp(-FDPart3tmp0 * xx0);
        const REAL FDPart3tmp4 = FDPart3tmp2 - FDPart3tmp3;
        const REAL FDPart3tmp6 = FDPart3tmp5 / ((FDPart3tmp4) * (FDPart3tmp4));
        const REAL FDPart3tmp10 = FDPart3tmp0 * FDPart3tmp2 + FDPart3tmp0 * FDPart3tmp3;
        const REAL FDPart3tmp11 = FDPart3tmp5 / ((FDPart3tmp10) * (FDPart3tmp10));
        dest_gf_address[IDX3(i0, i1, i2)] = (1.0 / 8.0) * ADD_times_AUU / pow(psi_background + uu, 7) + FDPart3tmp11 * uu_dDD00 +
                                            FDPart3tmp6 * uu_dDD11 + FDPart3tmp6 * uu_dD1 * cos(xx1) / FDPart3tmp7 +
                                            FDPart3tmp6 * uu_dDD22 / ((FDPart3tmp7) * (FDPart3tmp7)) -
                                            uu_dD0 * (-FDPart3tmp11 * (2 * FDPart3tmp0 * FDPart3tmp2 + 2 * FDPart3tmp0 * FDPart3tmp3) / FDPart3tmp4 +
                                                      (1.0 / 2.0) * FDPart3tmp5 * (FDPart3tmp12 * FDPart3tmp2 - FDPart3tmp12 * FDPart3tmp3) /
                                                          ((FDPart3tmp10) * (FDPart3tmp10) * (FDPart3tmp10)));

      } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
} // END FUNCTION residual_H_compute_all_points__rfm__SinhSpherical
