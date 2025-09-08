#include "BHaH_defines.h"
#include "intrinsics/simd_intrinsics.h"
/**
 * Kernel: compute_residual_all_points_host.
 * Kernel to compute the residual throughout the grid.
 */
static void compute_residual_all_points_host(const params_struct *restrict params, const rfm_struct *restrict rfmstruct,
                                             const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict aux_gfs) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  const REAL NOSIMDinvdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL_SIMD_ARRAY invdxx0 = ConstSIMD(NOSIMDinvdxx0);
  const REAL NOSIMDinvdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL_SIMD_ARRAY invdxx1 = ConstSIMD(NOSIMDinvdxx1);
  const REAL NOSIMDinvdxx2 = params->invdxx2;
  MAYBE_UNUSED const REAL_SIMD_ARRAY invdxx2 = ConstSIMD(NOSIMDinvdxx2);

#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += simd_width) {

        /*
         * NRPy+-Generated GF Access/FD Code, Step 1 of 2:
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

        static const double dblFDPart1_Rational_4_3 = 4.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart1_Rational_4_3 = ConstSIMD(dblFDPart1_Rational_4_3);

        static const double dblFDPart1_Rational_5_2 = 5.0 / 2.0;
        const REAL_SIMD_ARRAY FDPart1_Rational_5_2 = ConstSIMD(dblFDPart1_Rational_5_2);

        const REAL_SIMD_ARRAY FDPart1tmp0 = MulSIMD(FDPart1_Rational_5_2, uu);
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
         * NRPy+-Generated GF Access/FD Code, Step 2 of 2:
         * Evaluate SymPy expressions and write to main memory.
         */
        static const double dblFDPart3_Integer_1 = 1.0;
        MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart3_Integer_1 = ConstSIMD(dblFDPart3_Integer_1);

        static const double dblFDPart3_Rational_1_8 = 1.0 / 8.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_8 = ConstSIMD(dblFDPart3_Rational_1_8);

        const REAL_SIMD_ARRAY __RHS_exp_0 =
            AddSIMD(AddSIMD(uu_dDD00, uu_dDD11),
                    FusedMulAddSIMD(FDPart3_Rational_1_8,
                                    DivSIMD(ADD_times_AUU,
                                            MulSIMD(MulSIMD(MulSIMD(MulSIMD(MulSIMD(MulSIMD(AddSIMD(psi_background, uu), AddSIMD(psi_background, uu)),
                                                                                    AddSIMD(psi_background, uu)),
                                                                            AddSIMD(psi_background, uu)),
                                                                    AddSIMD(psi_background, uu)),
                                                            AddSIMD(psi_background, uu)),
                                                    AddSIMD(psi_background, uu))),
                                    uu_dDD22));

        WriteSIMD(&aux_gfs[IDX4(RESIDUAL_HGF, i0, i1, i2)], __RHS_exp_0);

      } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += simd_width)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
} // END FUNCTION compute_residual_all_points_host

/**
 * Compute residual of the Hamiltonian constraint for the hyperbolic relaxation equation.
 */
void compute_residual_all_points__rfm__Cartesian(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                 const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs,
                                                 REAL *restrict aux_gfs) {
  compute_residual_all_points_host(params, rfmstruct, auxevol_gfs, in_gfs, aux_gfs);

} // END FUNCTION compute_residual_all_points__rfm__Cartesian
