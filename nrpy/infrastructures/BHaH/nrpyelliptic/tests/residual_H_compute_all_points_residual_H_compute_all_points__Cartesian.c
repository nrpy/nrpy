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
void residual_H_compute_all_points__rfm__Cartesian(const commondata_struct *restrict commondata, const params_struct *restrict params,
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
         * NRPy+-Generated GF Access/FD Code, Step 1 of 2:
         * Read gridfunction(s) from main memory and compute FD stencils as needed.
         */
        const REAL ADD_times_AUU = auxevol_gfs[IDX4(ADD_TIMES_AUUGF, i0, i1, i2)];
        const REAL psi_background = auxevol_gfs[IDX4(PSI_BACKGROUNDGF, i0, i1, i2)];
        const REAL uu_i2m2 = in_gfs[IDX4(UUGF, i0, i1, i2 - 2)];
        const REAL uu_i2m1 = in_gfs[IDX4(UUGF, i0, i1, i2 - 1)];
        const REAL uu_i1m2 = in_gfs[IDX4(UUGF, i0, i1 - 2, i2)];
        const REAL uu_i1m1 = in_gfs[IDX4(UUGF, i0, i1 - 1, i2)];
        const REAL uu_i0m2 = in_gfs[IDX4(UUGF, i0 - 2, i1, i2)];
        const REAL uu_i0m1 = in_gfs[IDX4(UUGF, i0 - 1, i1, i2)];
        const REAL uu = in_gfs[IDX4(UUGF, i0, i1, i2)];
        const REAL uu_i0p1 = in_gfs[IDX4(UUGF, i0 + 1, i1, i2)];
        const REAL uu_i0p2 = in_gfs[IDX4(UUGF, i0 + 2, i1, i2)];
        const REAL uu_i1p1 = in_gfs[IDX4(UUGF, i0, i1 + 1, i2)];
        const REAL uu_i1p2 = in_gfs[IDX4(UUGF, i0, i1 + 2, i2)];
        const REAL uu_i2p1 = in_gfs[IDX4(UUGF, i0, i1, i2 + 1)];
        const REAL uu_i2p2 = in_gfs[IDX4(UUGF, i0, i1, i2 + 2)];
        static const REAL FDPart1_Rational_5_2 = 5.0 / 2.0;
        static const REAL FDPart1_Rational_1_12 = 1.0 / 12.0;
        static const REAL FDPart1_Rational_4_3 = 4.0 / 3.0;
        const REAL FDPart1tmp0 = -FDPart1_Rational_5_2 * uu;
        const REAL uu_dDD00 =
            ((invdxx0) * (invdxx0)) * (FDPart1_Rational_1_12 * (-uu_i0m2 - uu_i0p2) + FDPart1_Rational_4_3 * (uu_i0m1 + uu_i0p1) + FDPart1tmp0);
        const REAL uu_dDD11 =
            ((invdxx1) * (invdxx1)) * (FDPart1_Rational_1_12 * (-uu_i1m2 - uu_i1p2) + FDPart1_Rational_4_3 * (uu_i1m1 + uu_i1p1) + FDPart1tmp0);
        const REAL uu_dDD22 =
            ((invdxx2) * (invdxx2)) * (FDPart1_Rational_1_12 * (-uu_i2m2 - uu_i2p2) + FDPart1_Rational_4_3 * (uu_i2m1 + uu_i2p1) + FDPart1tmp0);

        /*
         * NRPy+-Generated GF Access/FD Code, Step 2 of 2:
         * Evaluate SymPy expressions and write to main memory.
         */
        dest_gf_address[IDX3(i0, i1, i2)] = (1.0 / 8.0) * ADD_times_AUU / pow(psi_background + uu, 7) + uu_dDD00 + uu_dDD11 + uu_dDD22;

      } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
} // END FUNCTION residual_H_compute_all_points__rfm__Cartesian
