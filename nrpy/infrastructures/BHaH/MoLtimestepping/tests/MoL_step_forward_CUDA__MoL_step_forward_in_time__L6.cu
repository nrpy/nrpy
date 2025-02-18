#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "intrinsics/cuda_intrinsics.h"

#define LOOP_ALL_GFS_GPS(ii)                                                                                                                         \
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;                                                                                            \
  const int stride0 = blockDim.x * gridDim.x;                                                                                                        \
  for (int(ii) = (tid0);                                                                                                                             \
       (ii) < d_params[streamid].Nxx_plus_2NGHOSTS0 * d_params[streamid].Nxx_plus_2NGHOSTS1 * d_params[streamid].Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;  \
       (ii) += (stride0))
/**
 * Kernel: rk_substep_1_gpu.
 * Compute RK substep 1.
 */
__global__ static void rk_substep_1_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs,
                                        const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    const REAL_CUDA_ARRAY __rk_exp_0 = FusedMulAddCUDA(k1_gfsL, dt, y_n_gfsL);
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_1_gpu

/**
 * Runge-Kutta function for substep 1.
 */
static void rk_substep_1__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs,
                                   const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  {

    const size_t threads_in_x_dir = 32;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rk_substep_1_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, y_n_gfs, next_y_input_gfs, dt);
    cudaCheckErrors(cudaKernel, "rk_substep_1_gpu failure");
  }
} // END FUNCTION rk_substep_1__launcher

/**
 * Kernel: rk_substep_2_gpu.
 * Compute RK substep 2.
 */
__global__ static void rk_substep_2_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k2_gfs, REAL *restrict y_n_gfs,
                                        REAL *restrict next_y_input_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k2_gfsL = k2_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dblRK_Rational_1_8 = 1.0 / 8.0;
    const REAL_CUDA_ARRAY RK_Rational_1_8 = ConstCUDA(dblRK_Rational_1_8);

    static const double dblRK_Rational_3_8 = 3.0 / 8.0;
    const REAL_CUDA_ARRAY RK_Rational_3_8 = ConstCUDA(dblRK_Rational_3_8);

    const REAL_CUDA_ARRAY __rk_exp_0 =
        FusedMulAddCUDA(RK_Rational_1_8, MulCUDA(k2_gfsL, dt), FusedMulAddCUDA(RK_Rational_3_8, MulCUDA(k1_gfsL, dt), y_n_gfsL));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_2_gpu

/**
 * Runge-Kutta function for substep 2.
 */
static void rk_substep_2__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict k2_gfs, REAL *restrict y_n_gfs,
                                   REAL *restrict next_y_input_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  {

    const size_t threads_in_x_dir = 32;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rk_substep_2_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k2_gfs, y_n_gfs, next_y_input_gfs, dt);
    cudaCheckErrors(cudaKernel, "rk_substep_2_gpu failure");
  }
} // END FUNCTION rk_substep_2__launcher

/**
 * Kernel: rk_substep_3_gpu.
 * Compute RK substep 3.
 */
__global__ static void rk_substep_3_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k2_gfs, REAL *restrict k3_gfs,
                                        REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k2_gfsL = k2_gfs[i];
    const REAL k3_gfsL = k3_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dblRK_Rational_2_27 = 2.0 / 27.0;
    const REAL_CUDA_ARRAY RK_Rational_2_27 = ConstCUDA(dblRK_Rational_2_27);

    static const double dblRK_Rational_8_27 = 8.0 / 27.0;
    const REAL_CUDA_ARRAY RK_Rational_8_27 = ConstCUDA(dblRK_Rational_8_27);

    const REAL_CUDA_ARRAY __rk_exp_0 = FusedMulAddCUDA(RK_Rational_8_27, FusedMulAddCUDA(k1_gfsL, dt, MulCUDA(k3_gfsL, dt)),
                                                       FusedMulAddCUDA(RK_Rational_2_27, MulCUDA(k2_gfsL, dt), y_n_gfsL));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_3_gpu

/**
 * Runge-Kutta function for substep 3.
 */
static void rk_substep_3__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict k2_gfs, REAL *restrict k3_gfs,
                                   REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  {

    const size_t threads_in_x_dir = 32;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rk_substep_3_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k2_gfs, k3_gfs, y_n_gfs, next_y_input_gfs, dt);
    cudaCheckErrors(cudaKernel, "rk_substep_3_gpu failure");
  }
} // END FUNCTION rk_substep_3__launcher

/**
 * Kernel: rk_substep_4_gpu.
 * Compute RK substep 4.
 */
__global__ static void rk_substep_4_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k2_gfs, REAL *restrict k3_gfs,
                                        REAL *restrict k4_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k2_gfsL = k2_gfs[i];
    const REAL k3_gfsL = k3_gfs[i];
    const REAL k4_gfsL = k4_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dbl_NegativeOne_ = -1.0;
    MAYBE_UNUSED const REAL_CUDA_ARRAY _NegativeOne_ = ConstCUDA(dbl_NegativeOne_);

    static const double dblRK_Rational_1_49 = 1.0 / 49.0;
    const REAL_CUDA_ARRAY RK_Rational_1_49 = ConstCUDA(dblRK_Rational_1_49);

    static const double dblRK_Rational_1_7 = 1.0 / 7.0;
    const REAL_CUDA_ARRAY RK_Rational_1_7 = ConstCUDA(dblRK_Rational_1_7);

    static const double dblRK_Rational_3_392 = 3.0 / 392.0;
    const REAL_CUDA_ARRAY RK_Rational_3_392 = ConstCUDA(dblRK_Rational_3_392);

    static const double dblRK_Rational_3_56 = 3.0 / 56.0;
    const REAL_CUDA_ARRAY RK_Rational_3_56 = ConstCUDA(dblRK_Rational_3_56);

    static const double dblRK_Rational_6_49 = 6.0 / 49.0;
    const REAL_CUDA_ARRAY RK_Rational_6_49 = ConstCUDA(dblRK_Rational_6_49);

    static const double dblRK_Rational_6_7 = 6.0 / 7.0;
    const REAL_CUDA_ARRAY RK_Rational_6_7 = ConstCUDA(dblRK_Rational_6_7);

    static const double dblRK_Rational_9_392 = 9.0 / 392.0;
    const REAL_CUDA_ARRAY RK_Rational_9_392 = ConstCUDA(dblRK_Rational_9_392);

    static const double dblRK_Rational_9_56 = 9.0 / 56.0;
    const REAL_CUDA_ARRAY RK_Rational_9_56 = ConstCUDA(dblRK_Rational_9_56);

    const REAL_CUDA_ARRAY __rk_exp_0 = FusedMulAddCUDA(
        k2_gfsL, MulCUDA(dt, FusedMulSubCUDA(RK_Rational_1_49, SqrtCUDA(21), RK_Rational_1_7)),
        FusedMulAddCUDA(
            k3_gfsL, MulCUDA(dt, NegFusedMulAddCUDA(RK_Rational_6_49, SqrtCUDA(21), RK_Rational_6_7)),
            FusedMulAddCUDA(k4_gfsL, MulCUDA(dt, FusedMulSubCUDA(RK_Rational_3_392, SqrtCUDA(21), RK_Rational_9_56)),
                            FusedMulAddCUDA(k1_gfsL, MulCUDA(dt, FusedMulSubCUDA(RK_Rational_9_392, SqrtCUDA(21), RK_Rational_3_56)), y_n_gfsL))));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_4_gpu

/**
 * Runge-Kutta function for substep 4.
 */
static void rk_substep_4__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict k2_gfs, REAL *restrict k3_gfs,
                                   REAL *restrict k4_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  {

    const size_t threads_in_x_dir = 32;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rk_substep_4_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k2_gfs, k3_gfs, k4_gfs, y_n_gfs,
                                                                                    next_y_input_gfs, dt);
    cudaCheckErrors(cudaKernel, "rk_substep_4_gpu failure");
  }
} // END FUNCTION rk_substep_4__launcher

/**
 * Kernel: rk_substep_5_gpu.
 * Compute RK substep 5.
 */
__global__ static void rk_substep_5_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k2_gfs, REAL *restrict k3_gfs,
                                        REAL *restrict k4_gfs, REAL *restrict k5_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs,
                                        const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k2_gfsL = k2_gfs[i];
    const REAL k3_gfsL = k3_gfs[i];
    const REAL k4_gfsL = k4_gfs[i];
    const REAL k5_gfsL = k5_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dbl_NegativeOne_ = -1.0;
    MAYBE_UNUSED const REAL_CUDA_ARRAY _NegativeOne_ = ConstCUDA(dbl_NegativeOne_);

    static const double dblRK_Rational_1_49 = 1.0 / 49.0;
    const REAL_CUDA_ARRAY RK_Rational_1_49 = ConstCUDA(dblRK_Rational_1_49);

    static const double dblRK_Rational_1_5 = 1.0 / 5.0;
    const REAL_CUDA_ARRAY RK_Rational_1_5 = ConstCUDA(dblRK_Rational_1_5);

    static const double dblRK_Rational_1_7 = 1.0 / 7.0;
    const REAL_CUDA_ARRAY RK_Rational_1_7 = ConstCUDA(dblRK_Rational_1_7);

    static const double dblRK_Rational_33_56 = 33.0 / 56.0;
    const REAL_CUDA_ARRAY RK_Rational_33_56 = ConstCUDA(dblRK_Rational_33_56);

    static const double dblRK_Rational_363_1960 = 363.0 / 1960.0;
    const REAL_CUDA_ARRAY RK_Rational_363_1960 = ConstCUDA(dblRK_Rational_363_1960);

    static const double dblRK_Rational_51_392 = 51.0 / 392.0;
    const REAL_CUDA_ARRAY RK_Rational_51_392 = ConstCUDA(dblRK_Rational_51_392);

    static const double dblRK_Rational_6_5 = 6.0 / 5.0;
    const REAL_CUDA_ARRAY RK_Rational_6_5 = ConstCUDA(dblRK_Rational_6_5);

    static const double dblRK_Rational_8_49 = 8.0 / 49.0;
    const REAL_CUDA_ARRAY RK_Rational_8_49 = ConstCUDA(dblRK_Rational_8_49);

    static const double dblRK_Rational_9_280 = 9.0 / 280.0;
    const REAL_CUDA_ARRAY RK_Rational_9_280 = ConstCUDA(dblRK_Rational_9_280);

    const REAL_CUDA_ARRAY tmp1 = MulCUDA(_NegativeOne_, SqrtCUDA(21));
    const REAL_CUDA_ARRAY __rk_exp_0 = FusedMulAddCUDA(
        k2_gfsL, MulCUDA(dt, FusedMulSubCUDA(RK_Rational_1_49, tmp1, RK_Rational_1_7)),
        FusedMulAddCUDA(
            k4_gfsL, MulCUDA(dt, FusedMulAddCUDA(RK_Rational_363_1960, SqrtCUDA(21), RK_Rational_9_280)),
            FusedMulAddCUDA(
                k5_gfsL, MulCUDA(dt, FusedMulAddCUDA(RK_Rational_1_5, SqrtCUDA(21), RK_Rational_6_5)),
                FusedMulAddCUDA(MulCUDA(RK_Rational_8_49, k3_gfsL), MulCUDA(tmp1, dt),
                                FusedMulAddCUDA(k1_gfsL, MulCUDA(dt, FusedMulSubCUDA(RK_Rational_51_392, tmp1, RK_Rational_33_56)), y_n_gfsL)))));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_5_gpu

/**
 * Runge-Kutta function for substep 5.
 */
static void rk_substep_5__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict k2_gfs, REAL *restrict k3_gfs,
                                   REAL *restrict k4_gfs, REAL *restrict k5_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs,
                                   const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  {

    const size_t threads_in_x_dir = 32;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rk_substep_5_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k2_gfs, k3_gfs, k4_gfs, k5_gfs, y_n_gfs,
                                                                                    next_y_input_gfs, dt);
    cudaCheckErrors(cudaKernel, "rk_substep_5_gpu failure");
  }
} // END FUNCTION rk_substep_5__launcher

/**
 * Kernel: rk_substep_6_gpu.
 * Compute RK substep 6.
 */
__global__ static void rk_substep_6_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k2_gfs, REAL *restrict k3_gfs,
                                        REAL *restrict k4_gfs, REAL *restrict k5_gfs, REAL *restrict k6_gfs, REAL *restrict y_n_gfs,
                                        REAL *restrict next_y_input_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k2_gfsL = k2_gfs[i];
    const REAL k3_gfsL = k3_gfs[i];
    const REAL k4_gfsL = k4_gfs[i];
    const REAL k5_gfsL = k5_gfs[i];
    const REAL k6_gfsL = k6_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dbl_NegativeOne_ = -1.0;
    MAYBE_UNUSED const REAL_CUDA_ARRAY _NegativeOne_ = ConstCUDA(dbl_NegativeOne_);

    static const double dblRK_Rational_10_9 = 10.0 / 9.0;
    const REAL_CUDA_ARRAY RK_Rational_10_9 = ConstCUDA(dblRK_Rational_10_9);

    static const double dblRK_Rational_11_6 = 11.0 / 6.0;
    const REAL_CUDA_ARRAY RK_Rational_11_6 = ConstCUDA(dblRK_Rational_11_6);

    static const double dblRK_Rational_14_9 = 14.0 / 9.0;
    const REAL_CUDA_ARRAY RK_Rational_14_9 = ConstCUDA(dblRK_Rational_14_9);

    static const double dblRK_Rational_21_20 = 21.0 / 20.0;
    const REAL_CUDA_ARRAY RK_Rational_21_20 = ConstCUDA(dblRK_Rational_21_20);

    static const double dblRK_Rational_2_3 = 2.0 / 3.0;
    const REAL_CUDA_ARRAY RK_Rational_2_3 = ConstCUDA(dblRK_Rational_2_3);

    static const double dblRK_Rational_343_90 = 343.0 / 90.0;
    const REAL_CUDA_ARRAY RK_Rational_343_90 = ConstCUDA(dblRK_Rational_343_90);

    static const double dblRK_Rational_49_18 = 49.0 / 18.0;
    const REAL_CUDA_ARRAY RK_Rational_49_18 = ConstCUDA(dblRK_Rational_49_18);

    static const double dblRK_Rational_7_10 = 7.0 / 10.0;
    const REAL_CUDA_ARRAY RK_Rational_7_10 = ConstCUDA(dblRK_Rational_7_10);

    static const double dblRK_Rational_7_12 = 7.0 / 12.0;
    const REAL_CUDA_ARRAY RK_Rational_7_12 = ConstCUDA(dblRK_Rational_7_12);

    static const double dblRK_Rational_7_18 = 7.0 / 18.0;
    const REAL_CUDA_ARRAY RK_Rational_7_18 = ConstCUDA(dblRK_Rational_7_18);

    const REAL_CUDA_ARRAY tmp1 = MulCUDA(_NegativeOne_, SqrtCUDA(21));
    const REAL_CUDA_ARRAY __rk_exp_0 = FusedMulAddCUDA(
        k4_gfsL, MulCUDA(dt, FusedMulAddCUDA(RK_Rational_21_20, tmp1, RK_Rational_7_10)),
        FusedMulAddCUDA(
            k1_gfsL, MulCUDA(dt, FusedMulAddCUDA(RK_Rational_7_12, SqrtCUDA(21), RK_Rational_11_6)),
            FusedMulAddCUDA(k3_gfsL, MulCUDA(dt, FusedMulSubCUDA(RK_Rational_14_9, SqrtCUDA(21), RK_Rational_10_9)),
                            FusedMulAddCUDA(k5_gfsL, MulCUDA(dt, FusedMulSubCUDA(RK_Rational_7_10, tmp1, RK_Rational_343_90)),
                                            FusedMulAddCUDA(k6_gfsL, MulCUDA(dt, FusedMulAddCUDA(RK_Rational_7_18, tmp1, RK_Rational_49_18)),
                                                            FusedMulAddCUDA(RK_Rational_2_3, MulCUDA(k2_gfsL, dt), y_n_gfsL))))));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_6_gpu

/**
 * Runge-Kutta function for substep 6.
 */
static void rk_substep_6__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict k2_gfs, REAL *restrict k3_gfs,
                                   REAL *restrict k4_gfs, REAL *restrict k5_gfs, REAL *restrict k6_gfs, REAL *restrict y_n_gfs,
                                   REAL *restrict next_y_input_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  {

    const size_t threads_in_x_dir = 32;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rk_substep_6_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k2_gfs, k3_gfs, k4_gfs, k5_gfs, k6_gfs, y_n_gfs,
                                                                                    next_y_input_gfs, dt);
    cudaCheckErrors(cudaKernel, "rk_substep_6_gpu failure");
  }
} // END FUNCTION rk_substep_6__launcher

/**
 * Kernel: rk_substep_7_gpu.
 * Compute RK substep 7.
 */
__global__ static void rk_substep_7_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k3_gfs, REAL *restrict k5_gfs,
                                        REAL *restrict k6_gfs, REAL *restrict k7_gfs, REAL *restrict y_n_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k3_gfsL = k3_gfs[i];
    const REAL k5_gfsL = k5_gfs[i];
    const REAL k6_gfsL = k6_gfs[i];
    const REAL k7_gfsL = k7_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dblRK_Rational_16_45 = 16.0 / 45.0;
    const REAL_CUDA_ARRAY RK_Rational_16_45 = ConstCUDA(dblRK_Rational_16_45);

    static const double dblRK_Rational_1_20 = 1.0 / 20.0;
    const REAL_CUDA_ARRAY RK_Rational_1_20 = ConstCUDA(dblRK_Rational_1_20);

    static const double dblRK_Rational_49_180 = 49.0 / 180.0;
    const REAL_CUDA_ARRAY RK_Rational_49_180 = ConstCUDA(dblRK_Rational_49_180);

    const REAL_CUDA_ARRAY __rk_exp_0 =
        FusedMulAddCUDA(RK_Rational_49_180, FusedMulAddCUDA(k5_gfsL, dt, MulCUDA(k6_gfsL, dt)),
                        FusedMulAddCUDA(RK_Rational_16_45, MulCUDA(k3_gfsL, dt),
                                        FusedMulAddCUDA(RK_Rational_1_20, FusedMulAddCUDA(k1_gfsL, dt, MulCUDA(k7_gfsL, dt)), y_n_gfsL)));
    WriteCUDA(&y_n_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_7_gpu

/**
 * Runge-Kutta function for substep 7.
 */
static void rk_substep_7__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict k3_gfs, REAL *restrict k5_gfs,
                                   REAL *restrict k6_gfs, REAL *restrict k7_gfs, REAL *restrict y_n_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  {

    const size_t threads_in_x_dir = 32;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rk_substep_7_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k3_gfs, k5_gfs, k6_gfs, k7_gfs, y_n_gfs, dt);
    cudaCheckErrors(cudaKernel, "rk_substep_7_gpu failure");
  }
} // END FUNCTION rk_substep_7__launcher

/**
 * Method of Lines (MoL) for "L6" method: Step forward one full timestep.
 *
 */
void MoL_step_forward_in_time(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {

  // C code implementation of -={ L6 }=- Method of Lines timestepping.

  // First set the initial time:
  const REAL time_start = commondata->time;
  // -={ START k1 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 0.00000000000000000e+00 * commondata->dt;
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict next_y_input_gfs = griddata[grid].gridfuncs.next_y_input_gfs;
    MAYBE_UNUSED REAL *restrict k1_gfs = griddata[grid].gridfuncs.k1_gfs;
    MAYBE_UNUSED REAL *restrict k2_gfs = griddata[grid].gridfuncs.k2_gfs;
    MAYBE_UNUSED REAL *restrict k3_gfs = griddata[grid].gridfuncs.k3_gfs;
    MAYBE_UNUSED REAL *restrict k4_gfs = griddata[grid].gridfuncs.k4_gfs;
    MAYBE_UNUSED REAL *restrict k5_gfs = griddata[grid].gridfuncs.k5_gfs;
    MAYBE_UNUSED REAL *restrict k6_gfs = griddata[grid].gridfuncs.k6_gfs;
    MAYBE_UNUSED REAL *restrict k7_gfs = griddata[grid].gridfuncs.k7_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, y_n_gfs, k1_gfs);
    rk_substep_1__launcher(params, k1_gfs, y_n_gfs, next_y_input_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k1 substep }=-

  // -={ START k2 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 1.00000000000000000e+00 * commondata->dt;
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict next_y_input_gfs = griddata[grid].gridfuncs.next_y_input_gfs;
    MAYBE_UNUSED REAL *restrict k1_gfs = griddata[grid].gridfuncs.k1_gfs;
    MAYBE_UNUSED REAL *restrict k2_gfs = griddata[grid].gridfuncs.k2_gfs;
    MAYBE_UNUSED REAL *restrict k3_gfs = griddata[grid].gridfuncs.k3_gfs;
    MAYBE_UNUSED REAL *restrict k4_gfs = griddata[grid].gridfuncs.k4_gfs;
    MAYBE_UNUSED REAL *restrict k5_gfs = griddata[grid].gridfuncs.k5_gfs;
    MAYBE_UNUSED REAL *restrict k6_gfs = griddata[grid].gridfuncs.k6_gfs;
    MAYBE_UNUSED REAL *restrict k7_gfs = griddata[grid].gridfuncs.k7_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k2_gfs);
    rk_substep_2__launcher(params, k1_gfs, k2_gfs, y_n_gfs, next_y_input_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k2 substep }=-

  // -={ START k3 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 5.00000000000000000e-01 * commondata->dt;
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict next_y_input_gfs = griddata[grid].gridfuncs.next_y_input_gfs;
    MAYBE_UNUSED REAL *restrict k1_gfs = griddata[grid].gridfuncs.k1_gfs;
    MAYBE_UNUSED REAL *restrict k2_gfs = griddata[grid].gridfuncs.k2_gfs;
    MAYBE_UNUSED REAL *restrict k3_gfs = griddata[grid].gridfuncs.k3_gfs;
    MAYBE_UNUSED REAL *restrict k4_gfs = griddata[grid].gridfuncs.k4_gfs;
    MAYBE_UNUSED REAL *restrict k5_gfs = griddata[grid].gridfuncs.k5_gfs;
    MAYBE_UNUSED REAL *restrict k6_gfs = griddata[grid].gridfuncs.k6_gfs;
    MAYBE_UNUSED REAL *restrict k7_gfs = griddata[grid].gridfuncs.k7_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k3_gfs);
    rk_substep_3__launcher(params, k1_gfs, k2_gfs, k3_gfs, y_n_gfs, next_y_input_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k3 substep }=-

  // -={ START k4 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 6.66666666666666630e-01 * commondata->dt;
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict next_y_input_gfs = griddata[grid].gridfuncs.next_y_input_gfs;
    MAYBE_UNUSED REAL *restrict k1_gfs = griddata[grid].gridfuncs.k1_gfs;
    MAYBE_UNUSED REAL *restrict k2_gfs = griddata[grid].gridfuncs.k2_gfs;
    MAYBE_UNUSED REAL *restrict k3_gfs = griddata[grid].gridfuncs.k3_gfs;
    MAYBE_UNUSED REAL *restrict k4_gfs = griddata[grid].gridfuncs.k4_gfs;
    MAYBE_UNUSED REAL *restrict k5_gfs = griddata[grid].gridfuncs.k5_gfs;
    MAYBE_UNUSED REAL *restrict k6_gfs = griddata[grid].gridfuncs.k6_gfs;
    MAYBE_UNUSED REAL *restrict k7_gfs = griddata[grid].gridfuncs.k7_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k4_gfs);
    rk_substep_4__launcher(params, k1_gfs, k2_gfs, k3_gfs, k4_gfs, y_n_gfs, next_y_input_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k4 substep }=-

  // -={ START k5 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 1.72673164646011429e-01 * commondata->dt;
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict next_y_input_gfs = griddata[grid].gridfuncs.next_y_input_gfs;
    MAYBE_UNUSED REAL *restrict k1_gfs = griddata[grid].gridfuncs.k1_gfs;
    MAYBE_UNUSED REAL *restrict k2_gfs = griddata[grid].gridfuncs.k2_gfs;
    MAYBE_UNUSED REAL *restrict k3_gfs = griddata[grid].gridfuncs.k3_gfs;
    MAYBE_UNUSED REAL *restrict k4_gfs = griddata[grid].gridfuncs.k4_gfs;
    MAYBE_UNUSED REAL *restrict k5_gfs = griddata[grid].gridfuncs.k5_gfs;
    MAYBE_UNUSED REAL *restrict k6_gfs = griddata[grid].gridfuncs.k6_gfs;
    MAYBE_UNUSED REAL *restrict k7_gfs = griddata[grid].gridfuncs.k7_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k5_gfs);
    rk_substep_5__launcher(params, k1_gfs, k2_gfs, k3_gfs, k4_gfs, k5_gfs, y_n_gfs, next_y_input_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k5 substep }=-

  // -={ START k6 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 8.27326835353988543e-01 * commondata->dt;
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict next_y_input_gfs = griddata[grid].gridfuncs.next_y_input_gfs;
    MAYBE_UNUSED REAL *restrict k1_gfs = griddata[grid].gridfuncs.k1_gfs;
    MAYBE_UNUSED REAL *restrict k2_gfs = griddata[grid].gridfuncs.k2_gfs;
    MAYBE_UNUSED REAL *restrict k3_gfs = griddata[grid].gridfuncs.k3_gfs;
    MAYBE_UNUSED REAL *restrict k4_gfs = griddata[grid].gridfuncs.k4_gfs;
    MAYBE_UNUSED REAL *restrict k5_gfs = griddata[grid].gridfuncs.k5_gfs;
    MAYBE_UNUSED REAL *restrict k6_gfs = griddata[grid].gridfuncs.k6_gfs;
    MAYBE_UNUSED REAL *restrict k7_gfs = griddata[grid].gridfuncs.k7_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k6_gfs);
    rk_substep_6__launcher(params, k1_gfs, k2_gfs, k3_gfs, k4_gfs, k5_gfs, k6_gfs, y_n_gfs, next_y_input_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k6 substep }=-

  // -={ START k7 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 1.00000000000000000e+00 * commondata->dt;
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict next_y_input_gfs = griddata[grid].gridfuncs.next_y_input_gfs;
    MAYBE_UNUSED REAL *restrict k1_gfs = griddata[grid].gridfuncs.k1_gfs;
    MAYBE_UNUSED REAL *restrict k2_gfs = griddata[grid].gridfuncs.k2_gfs;
    MAYBE_UNUSED REAL *restrict k3_gfs = griddata[grid].gridfuncs.k3_gfs;
    MAYBE_UNUSED REAL *restrict k4_gfs = griddata[grid].gridfuncs.k4_gfs;
    MAYBE_UNUSED REAL *restrict k5_gfs = griddata[grid].gridfuncs.k5_gfs;
    MAYBE_UNUSED REAL *restrict k6_gfs = griddata[grid].gridfuncs.k6_gfs;
    MAYBE_UNUSED REAL *restrict k7_gfs = griddata[grid].gridfuncs.k7_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k7_gfs);
    rk_substep_7__launcher(params, k1_gfs, k3_gfs, k5_gfs, k6_gfs, k7_gfs, y_n_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, y_n_gfs);
  }
  // -={ END k7 substep }=-

  // Adding dt to commondata->time many times will induce roundoff error,
  // so here we set time based on the iteration number:
  commondata->time = (REAL)(commondata->nn + 1) * commondata->dt;

  // Increment the timestep n:
  commondata->nn++;
} // END FUNCTION MoL_step_forward_in_time
