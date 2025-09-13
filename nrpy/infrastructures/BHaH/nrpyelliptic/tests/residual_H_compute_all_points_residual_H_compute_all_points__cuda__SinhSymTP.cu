#include "BHaH_defines.h"
#include "intrinsics/cuda_intrinsics.h"
/**
 * Kernel: residual_H_compute_all_points_gpu.
 * Kernel to compute the residual throughout the grid.
 */
__global__ static void residual_H_compute_all_points_gpu(const size_t streamid, const rfm_struct *restrict rfmstruct,
                                                         const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict aux_gfs) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = d_params[streamid].Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = d_params[streamid].Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = d_params[streamid].Nxx_plus_2NGHOSTS2;

  const REAL NOCUDAinvdxx0 = d_params[streamid].invdxx0;
  MAYBE_UNUSED const REAL_CUDA_ARRAY invdxx0 = ConstCUDA(NOCUDAinvdxx0);
  const REAL NOCUDAinvdxx1 = d_params[streamid].invdxx1;
  MAYBE_UNUSED const REAL_CUDA_ARRAY invdxx1 = ConstCUDA(NOCUDAinvdxx1);
  const REAL NOCUDAinvdxx2 = d_params[streamid].invdxx2;
  MAYBE_UNUSED const REAL_CUDA_ARRAY invdxx2 = ConstCUDA(NOCUDAinvdxx2);

  MAYBE_UNUSED const int tid0 = blockIdx.x * blockDim.x + threadIdx.x;
  MAYBE_UNUSED const int tid1 = blockIdx.y * blockDim.y + threadIdx.y;
  MAYBE_UNUSED const int tid2 = blockIdx.z * blockDim.z + threadIdx.z;

  MAYBE_UNUSED const int stride0 = blockDim.x * gridDim.x;
  MAYBE_UNUSED const int stride1 = blockDim.y * gridDim.y;
  MAYBE_UNUSED const int stride2 = blockDim.z * gridDim.z;

  for (int i2 = tid2 + NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2 += stride2) {
    for (int i1 = tid1 + NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1 += stride1) {
      const double NOCUDAf1_of_xx1 = rfmstruct->f1_of_xx1[i1];
      MAYBE_UNUSED const REAL_CUDA_ARRAY f1_of_xx1 = ConstCUDA(NOCUDAf1_of_xx1);
      const double NOCUDAf1_of_xx1__D1 = rfmstruct->f1_of_xx1__D1[i1];
      MAYBE_UNUSED const REAL_CUDA_ARRAY f1_of_xx1__D1 = ConstCUDA(NOCUDAf1_of_xx1__D1);
      const double NOCUDAf1_of_xx1__DD11 = rfmstruct->f1_of_xx1__DD11[i1];
      MAYBE_UNUSED const REAL_CUDA_ARRAY f1_of_xx1__DD11 = ConstCUDA(NOCUDAf1_of_xx1__DD11);
      const double NOCUDAf4_of_xx1 = rfmstruct->f4_of_xx1[i1];
      MAYBE_UNUSED const REAL_CUDA_ARRAY f4_of_xx1 = ConstCUDA(NOCUDAf4_of_xx1);
      const double NOCUDAf4_of_xx1__D1 = rfmstruct->f4_of_xx1__D1[i1];
      MAYBE_UNUSED const REAL_CUDA_ARRAY f4_of_xx1__D1 = ConstCUDA(NOCUDAf4_of_xx1__D1);
      const double NOCUDAf4_of_xx1__DD11 = rfmstruct->f4_of_xx1__DD11[i1];
      MAYBE_UNUSED const REAL_CUDA_ARRAY f4_of_xx1__DD11 = ConstCUDA(NOCUDAf4_of_xx1__DD11);

      for (int i0 = tid0 + NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += stride0) {
        MAYBE_UNUSED const REAL_CUDA_ARRAY f0_of_xx0 = ReadCUDA(&rfmstruct->f0_of_xx0[i0]);
        MAYBE_UNUSED const REAL_CUDA_ARRAY f0_of_xx0__D0 = ReadCUDA(&rfmstruct->f0_of_xx0__D0[i0]);
        MAYBE_UNUSED const REAL_CUDA_ARRAY f0_of_xx0__DD00 = ReadCUDA(&rfmstruct->f0_of_xx0__DD00[i0]);
        MAYBE_UNUSED const REAL_CUDA_ARRAY f0_of_xx0__DDD000 = ReadCUDA(&rfmstruct->f0_of_xx0__DDD000[i0]);
        MAYBE_UNUSED const REAL_CUDA_ARRAY f2_of_xx0 = ReadCUDA(&rfmstruct->f2_of_xx0[i0]);
        MAYBE_UNUSED const REAL_CUDA_ARRAY f2_of_xx0__D0 = ReadCUDA(&rfmstruct->f2_of_xx0__D0[i0]);
        MAYBE_UNUSED const REAL_CUDA_ARRAY f2_of_xx0__DD00 = ReadCUDA(&rfmstruct->f2_of_xx0__DD00[i0]);

        /*
         * NRPy+-Generated GF Access/FD Code, Step 1 of 2:
         * Read gridfunction(s) from main memory and compute FD stencils as needed.
         */
        const REAL_CUDA_ARRAY ADD_times_AUU = ReadCUDA(&auxevol_gfs[IDX4(ADD_TIMES_AUUGF, i0, i1, i2)]);
        const REAL_CUDA_ARRAY psi_background = ReadCUDA(&auxevol_gfs[IDX4(PSI_BACKGROUNDGF, i0, i1, i2)]);
        const REAL_CUDA_ARRAY uu_i2m2 = ReadCUDA(&in_gfs[IDX4(UUGF, i0, i1, i2 - 2)]);
        const REAL_CUDA_ARRAY uu_i2m1 = ReadCUDA(&in_gfs[IDX4(UUGF, i0, i1, i2 - 1)]);
        const REAL_CUDA_ARRAY uu_i1m2 = ReadCUDA(&in_gfs[IDX4(UUGF, i0, i1 - 2, i2)]);
        const REAL_CUDA_ARRAY uu_i1m1 = ReadCUDA(&in_gfs[IDX4(UUGF, i0, i1 - 1, i2)]);
        const REAL_CUDA_ARRAY uu_i0m2 = ReadCUDA(&in_gfs[IDX4(UUGF, i0 - 2, i1, i2)]);
        const REAL_CUDA_ARRAY uu_i0m1 = ReadCUDA(&in_gfs[IDX4(UUGF, i0 - 1, i1, i2)]);
        const REAL_CUDA_ARRAY uu = ReadCUDA(&in_gfs[IDX4(UUGF, i0, i1, i2)]);
        const REAL_CUDA_ARRAY uu_i0p1 = ReadCUDA(&in_gfs[IDX4(UUGF, i0 + 1, i1, i2)]);
        const REAL_CUDA_ARRAY uu_i0p2 = ReadCUDA(&in_gfs[IDX4(UUGF, i0 + 2, i1, i2)]);
        const REAL_CUDA_ARRAY uu_i1p1 = ReadCUDA(&in_gfs[IDX4(UUGF, i0, i1 + 1, i2)]);
        const REAL_CUDA_ARRAY uu_i1p2 = ReadCUDA(&in_gfs[IDX4(UUGF, i0, i1 + 2, i2)]);
        const REAL_CUDA_ARRAY uu_i2p1 = ReadCUDA(&in_gfs[IDX4(UUGF, i0, i1, i2 + 1)]);
        const REAL_CUDA_ARRAY uu_i2p2 = ReadCUDA(&in_gfs[IDX4(UUGF, i0, i1, i2 + 2)]);
        static const double dblFDPart1_NegativeOne_ = -1.0;
        MAYBE_UNUSED const REAL_CUDA_ARRAY FDPart1_NegativeOne_ = ConstCUDA(dblFDPart1_NegativeOne_);

        static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
        const REAL_CUDA_ARRAY FDPart1_Rational_1_12 = ConstCUDA(dblFDPart1_Rational_1_12);

        static const double dblFDPart1_Rational_2_3 = 2.0 / 3.0;
        const REAL_CUDA_ARRAY FDPart1_Rational_2_3 = ConstCUDA(dblFDPart1_Rational_2_3);

        static const double dblFDPart1_Rational_4_3 = 4.0 / 3.0;
        const REAL_CUDA_ARRAY FDPart1_Rational_4_3 = ConstCUDA(dblFDPart1_Rational_4_3);

        static const double dblFDPart1_Rational_5_2 = 5.0 / 2.0;
        const REAL_CUDA_ARRAY FDPart1_Rational_5_2 = ConstCUDA(dblFDPart1_Rational_5_2);

        const REAL_CUDA_ARRAY FDPart1tmp0 = MulCUDA(FDPart1_Rational_5_2, uu);
        const REAL_CUDA_ARRAY uu_dD0 = MulCUDA(
            invdxx0, FusedMulAddCUDA(FDPart1_Rational_1_12, SubCUDA(uu_i0m2, uu_i0p2), MulCUDA(FDPart1_Rational_2_3, SubCUDA(uu_i0p1, uu_i0m1))));
        const REAL_CUDA_ARRAY uu_dD1 = MulCUDA(
            invdxx1, FusedMulAddCUDA(FDPart1_Rational_1_12, SubCUDA(uu_i1m2, uu_i1p2), MulCUDA(FDPart1_Rational_2_3, SubCUDA(uu_i1p1, uu_i1m1))));
        const REAL_CUDA_ARRAY uu_dDD00 =
            MulCUDA(MulCUDA(invdxx0, invdxx0), FusedMulSubCUDA(FDPart1_Rational_4_3, AddCUDA(uu_i0m1, uu_i0p1),
                                                               FusedMulAddCUDA(FDPart1_Rational_1_12, AddCUDA(uu_i0m2, uu_i0p2), FDPart1tmp0)));
        const REAL_CUDA_ARRAY uu_dDD11 =
            MulCUDA(MulCUDA(invdxx1, invdxx1), FusedMulSubCUDA(FDPart1_Rational_4_3, AddCUDA(uu_i1m1, uu_i1p1),
                                                               FusedMulAddCUDA(FDPart1_Rational_1_12, AddCUDA(uu_i1m2, uu_i1p2), FDPart1tmp0)));
        const REAL_CUDA_ARRAY uu_dDD22 =
            MulCUDA(MulCUDA(invdxx2, invdxx2), FusedMulSubCUDA(FDPart1_Rational_4_3, AddCUDA(uu_i2m1, uu_i2p1),
                                                               FusedMulAddCUDA(FDPart1_Rational_1_12, AddCUDA(uu_i2m2, uu_i2p2), FDPart1tmp0)));

        /*
         * NRPy+-Generated GF Access/FD Code, Step 2 of 2:
         * Evaluate SymPy expressions and write to main memory.
         */
        static const double dblFDPart3_Integer_1 = 1.0;
        MAYBE_UNUSED const REAL_CUDA_ARRAY FDPart3_Integer_1 = ConstCUDA(dblFDPart3_Integer_1);

        static const double dblFDPart3_Integer_2 = 2.0;
        const REAL_CUDA_ARRAY FDPart3_Integer_2 = ConstCUDA(dblFDPart3_Integer_2);

        static const double dblFDPart3_NegativeOne_ = -1.0;
        MAYBE_UNUSED const REAL_CUDA_ARRAY FDPart3_NegativeOne_ = ConstCUDA(dblFDPart3_NegativeOne_);

        static const double dblFDPart3_Rational_1_2 = 1.0 / 2.0;
        const REAL_CUDA_ARRAY FDPart3_Rational_1_2 = ConstCUDA(dblFDPart3_Rational_1_2);

        static const double dblFDPart3_Rational_1_8 = 1.0 / 8.0;
        const REAL_CUDA_ARRAY FDPart3_Rational_1_8 = ConstCUDA(dblFDPart3_Rational_1_8);

        const REAL_CUDA_ARRAY FDPart3tmp4 = MulCUDA(f2_of_xx0, f2_of_xx0);
        const REAL_CUDA_ARRAY FDPart3tmp1 = FusedMulAddCUDA(f0_of_xx0, f0_of_xx0, MulCUDA(f4_of_xx1, f4_of_xx1));
        const REAL_CUDA_ARRAY FDPart3tmp8 = DivCUDA(FDPart3_Integer_2, FDPart3tmp4);
        const REAL_CUDA_ARRAY FDPart3tmp2 = DivCUDA(FDPart3_Integer_1, FDPart3tmp1);
        const REAL_CUDA_ARRAY FDPart3tmp6 = DivCUDA(FDPart3_Integer_1, MulCUDA(FDPart3tmp1, FDPart3tmp1));
        const REAL_CUDA_ARRAY __RHS_exp_0 = FusedMulAddCUDA(
            MulCUDA(FDPart3tmp2, FDPart3tmp4), DivCUDA(uu_dDD00, MulCUDA(f0_of_xx0__D0, f0_of_xx0__D0)),
            FusedMulAddCUDA(
                FDPart3_Rational_1_8,
                DivCUDA(ADD_times_AUU, MulCUDA(MulCUDA(MulCUDA(MulCUDA(MulCUDA(MulCUDA(AddCUDA(psi_background, uu), AddCUDA(psi_background, uu)),
                                                                               AddCUDA(psi_background, uu)),
                                                                       AddCUDA(psi_background, uu)),
                                                               AddCUDA(psi_background, uu)),
                                                       AddCUDA(psi_background, uu)),
                                               AddCUDA(psi_background, uu))),
                FusedMulAddCUDA(
                    uu_dDD22, DivCUDA(DivCUDA(FDPart3_Integer_1, MulCUDA(f1_of_xx1, f1_of_xx1)), MulCUDA(f0_of_xx0, f0_of_xx0)),
                    FusedMulAddCUDA(
                        MulCUDA(FDPart3tmp2, f1_of_xx1__D1), DivCUDA(uu_dD1, f1_of_xx1),
                        FusedMulSubCUDA(
                            FDPart3tmp2, uu_dDD11,
                            MulCUDA(
                                uu_dD0,
                                FusedMulSubCUDA(
                                    FDPart3tmp6, MulCUDA(MulCUDA(FDPart3_NegativeOne_, FDPart3tmp4), DivCUDA(f0_of_xx0, f0_of_xx0__D0)),
                                    FusedMulSubCUDA(
                                        MulCUDA(FDPart3tmp2, FDPart3tmp4), DivCUDA(DivCUDA(FDPart3_Integer_1, f0_of_xx0__D0), f0_of_xx0),
                                        DivCUDA(
                                            MulCUDA(MulCUDA(FDPart3_Rational_1_2, FDPart3tmp6),
                                                    MulCUDA(MulCUDA(MulCUDA(MulCUDA(f2_of_xx0, f2_of_xx0), f2_of_xx0), f2_of_xx0),
                                                            FusedMulAddCUDA(
                                                                MulCUDA(FDPart3tmp1, FDPart3tmp8), MulCUDA(f0_of_xx0__D0, f0_of_xx0__DD00),
                                                                FusedMulSubCUDA(
                                                                    FDPart3tmp8,
                                                                    MulCUDA(f0_of_xx0, MulCUDA(MulCUDA(f0_of_xx0__D0, f0_of_xx0__D0), f0_of_xx0__D0)),
                                                                    MulCUDA(f2_of_xx0__D0,
                                                                            MulCUDA(MulCUDA(FDPart3_Integer_2, FDPart3tmp1),
                                                                                    DivCUDA(MulCUDA(f0_of_xx0__D0, f0_of_xx0__D0),
                                                                                            MulCUDA(MulCUDA(f2_of_xx0, f2_of_xx0), f2_of_xx0)))))))),
                                            MulCUDA(MulCUDA(MulCUDA(f0_of_xx0__D0, f0_of_xx0__D0), f0_of_xx0__D0), f0_of_xx0__D0))))))))));

        WriteCUDA(&aux_gfs[IDX4(RESIDUAL_HGF, i0, i1, i2)], __RHS_exp_0);

      } // END LOOP: for (int i0 = tid0+NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += stride0)
    } // END LOOP: for (int i1 = tid1+NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1 += stride1)
  } // END LOOP: for (int i2 = tid2+NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2 += stride2)
} // END FUNCTION residual_H_compute_all_points_gpu

/**
 * Compute residual of the Hamiltonian constraint for the hyperbolic relaxation equation.
 */
void residual_H_compute_all_points__rfm__SinhSymTP(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                   const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs,
                                                   const REAL *restrict in_gfs, REAL *restrict aux_gfs) {

  const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_NELL_H;
  const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_NELL_H;
  const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_NELL_H;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((params->Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
                       (params->Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
                       (params->Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  residual_H_compute_all_points_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, auxevol_gfs, in_gfs, aux_gfs);
  cudaCheckErrors(cudaKernel, "residual_H_compute_all_points_gpu failure");

} // END FUNCTION residual_H_compute_all_points__rfm__SinhSymTP
