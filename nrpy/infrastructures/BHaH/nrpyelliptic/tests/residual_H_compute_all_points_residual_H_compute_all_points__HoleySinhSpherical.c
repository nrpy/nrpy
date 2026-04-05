#include "BHaH_defines.h"
#include "intrinsics/simd_intrinsics.h"

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
void residual_H_compute_all_points__rfm__HoleySinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                            const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs,
                                                            const REAL *restrict in_gfs, REAL *restrict dest_gf_address) {
#include "set_CodeParameters-simd.h"
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      const double NOSIMDf1_of_xx1 = rfmstruct->f1_of_xx1[i1];
      MAYBE_UNUSED const REAL_SIMD_ARRAY f1_of_xx1 = ConstSIMD(NOSIMDf1_of_xx1);
      const double NOSIMDf1_of_xx1__D1 = rfmstruct->f1_of_xx1__D1[i1];
      MAYBE_UNUSED const REAL_SIMD_ARRAY f1_of_xx1__D1 = ConstSIMD(NOSIMDf1_of_xx1__D1);
      const double NOSIMDf1_of_xx1__DD11 = rfmstruct->f1_of_xx1__DD11[i1];
      MAYBE_UNUSED const REAL_SIMD_ARRAY f1_of_xx1__DD11 = ConstSIMD(NOSIMDf1_of_xx1__DD11);

      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += simd_width) {
        MAYBE_UNUSED const REAL_SIMD_ARRAY f0_of_xx0 = ReadSIMD(&rfmstruct->f0_of_xx0[i0]);
        MAYBE_UNUSED const REAL_SIMD_ARRAY f0_of_xx0__D0 = ReadSIMD(&rfmstruct->f0_of_xx0__D0[i0]);
        MAYBE_UNUSED const REAL_SIMD_ARRAY f0_of_xx0__DD00 = ReadSIMD(&rfmstruct->f0_of_xx0__DD00[i0]);
        MAYBE_UNUSED const REAL_SIMD_ARRAY f0_of_xx0__DDD000 = ReadSIMD(&rfmstruct->f0_of_xx0__DDD000[i0]);

        /*
         * NRPy-Generated GF Access/FD Code, Step 1 of 2:
         * Read gridfunction(s) from main memory and compute FD stencils as needed.
         */
        const REAL_SIMD_ARRAY ADD_times_AUU = ReadSIMD(&auxevol_gfs[IDX4(ADD_TIMES_AUUGF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY psi_background = ReadSIMD(&auxevol_gfs[IDX4(PSI_BACKGROUNDGF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY uu_i2m2 = ReadSIMD(&in_gfs[IDX4(UUGF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY uu_i2m1 = ReadSIMD(&in_gfs[IDX4(UUGF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY uu_i1m2 = ReadSIMD(&in_gfs[IDX4(UUGF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY uu_i1m1 = ReadSIMD(&in_gfs[IDX4(UUGF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY uu_i0m2 = ReadSIMD(&in_gfs[IDX4(UUGF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY uu_i0m1 = ReadSIMD(&in_gfs[IDX4(UUGF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY uu = ReadSIMD(&in_gfs[IDX4(UUGF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY uu_i0p1 = ReadSIMD(&in_gfs[IDX4(UUGF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY uu_i0p2 = ReadSIMD(&in_gfs[IDX4(UUGF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY uu_i1p1 = ReadSIMD(&in_gfs[IDX4(UUGF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY uu_i1p2 = ReadSIMD(&in_gfs[IDX4(UUGF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY uu_i2p1 = ReadSIMD(&in_gfs[IDX4(UUGF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY uu_i2p2 = ReadSIMD(&in_gfs[IDX4(UUGF, i0, i1, i2 + 2)]);
        static const double dblFDPart1_NegativeOne_ = -1.0;
        MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart1_NegativeOne_ = ConstSIMD(dblFDPart1_NegativeOne_);

        static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
        const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

        static const double dblFDPart1_Rational_2_3 = 2.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart1_Rational_2_3 = ConstSIMD(dblFDPart1_Rational_2_3);

        static const double dblFDPart1_Rational_4_3 = 4.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart1_Rational_4_3 = ConstSIMD(dblFDPart1_Rational_4_3);

        static const double dblFDPart1_Rational_5_2 = 5.0 / 2.0;
        const REAL_SIMD_ARRAY FDPart1_Rational_5_2 = ConstSIMD(dblFDPart1_Rational_5_2);

        const REAL_SIMD_ARRAY FDPart1tmp0 = MulSIMD(FDPart1_Rational_5_2, uu);
        const REAL_SIMD_ARRAY uu_dD0 = MulSIMD(
            invdxx0, FusedMulAddSIMD(FDPart1_Rational_1_12, SubSIMD(uu_i0m2, uu_i0p2), MulSIMD(FDPart1_Rational_2_3, SubSIMD(uu_i0p1, uu_i0m1))));
        const REAL_SIMD_ARRAY uu_dD1 = MulSIMD(
            invdxx1, FusedMulAddSIMD(FDPart1_Rational_1_12, SubSIMD(uu_i1m2, uu_i1p2), MulSIMD(FDPart1_Rational_2_3, SubSIMD(uu_i1p1, uu_i1m1))));
        const REAL_SIMD_ARRAY uu_dDD00 =
            MulSIMD(MulSIMD(invdxx0, invdxx0), FusedMulSubSIMD(FDPart1_Rational_4_3, AddSIMD(uu_i0m1, uu_i0p1),
                                                               FusedMulAddSIMD(FDPart1_Rational_1_12, AddSIMD(uu_i0m2, uu_i0p2), FDPart1tmp0)));
        const REAL_SIMD_ARRAY uu_dDD11 =
            MulSIMD(MulSIMD(invdxx1, invdxx1), FusedMulSubSIMD(FDPart1_Rational_4_3, AddSIMD(uu_i1m1, uu_i1p1),
                                                               FusedMulAddSIMD(FDPart1_Rational_1_12, AddSIMD(uu_i1m2, uu_i1p2), FDPart1tmp0)));
        const REAL_SIMD_ARRAY uu_dDD22 =
            MulSIMD(MulSIMD(invdxx2, invdxx2), FusedMulSubSIMD(FDPart1_Rational_4_3, AddSIMD(uu_i2m1, uu_i2p1),
                                                               FusedMulAddSIMD(FDPart1_Rational_1_12, AddSIMD(uu_i2m2, uu_i2p2), FDPart1tmp0)));

        /*
         * NRPy-Generated GF Access/FD Code, Step 2 of 2:
         * Evaluate SymPy expressions and write to main memory.
         */
        static const double dblFDPart3_Integer_1 = 1.0;
        MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart3_Integer_1 = ConstSIMD(dblFDPart3_Integer_1);

        static const double dblFDPart3_Integer_2 = 2.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_2 = ConstSIMD(dblFDPart3_Integer_2);

        static const double dblFDPart3_NegativeOne_ = -1.0;
        MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart3_NegativeOne_ = ConstSIMD(dblFDPart3_NegativeOne_);

        static const double dblFDPart3_Rational_1_8 = 1.0 / 8.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_8 = ConstSIMD(dblFDPart3_Rational_1_8);

        const REAL_SIMD_ARRAY FDPart3tmp0 = DivSIMD(FDPart3_Integer_1, MulSIMD(f0_of_xx0, f0_of_xx0));
        const REAL_SIMD_ARRAY __RHS_exp_0 = FusedMulAddSIMD(
            FDPart3tmp0, DivSIMD(uu_dDD22, MulSIMD(f1_of_xx1, f1_of_xx1)),
            AddSIMD(FusedMulAddSIMD(FDPart3_Rational_1_8,
                                    DivSIMD(ADD_times_AUU,
                                            MulSIMD(MulSIMD(MulSIMD(MulSIMD(MulSIMD(MulSIMD(AddSIMD(psi_background, uu), AddSIMD(psi_background, uu)),
                                                                                    AddSIMD(psi_background, uu)),
                                                                            AddSIMD(psi_background, uu)),
                                                                    AddSIMD(psi_background, uu)),
                                                            AddSIMD(psi_background, uu)),
                                                    AddSIMD(psi_background, uu))),
                                    DivSIMD(uu_dDD00, MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp0, f1_of_xx1__D1), DivSIMD(uu_dD1, f1_of_xx1),
                        FusedMulSubSIMD(
                            FDPart3tmp0, uu_dDD11,
                            MulSIMD(uu_dD0,
                                    NegFusedMulAddSIMD(FDPart3_Integer_2, DivSIMD(DivSIMD(FDPart3_Integer_1, f0_of_xx0__D0), f0_of_xx0),
                                                       DivSIMD(f0_of_xx0__DD00, MulSIMD(MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0), f0_of_xx0__D0))))))));

        WriteSIMD(&dest_gf_address[IDX3(i0, i1, i2)], __RHS_exp_0);

      } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += simd_width)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
} // END FUNCTION residual_H_compute_all_points__rfm__HoleySinhSpherical
