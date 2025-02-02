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
 * GPU Kernel: rk_substep_1_gpu.
 * GPU Kernel to compute RK substep 1.
 */
__global__ static void rk_substep_1_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs,
                                        const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dblRK_Rational_1_18 = 1.0 / 18.0;
    const REAL_CUDA_ARRAY RK_Rational_1_18 = ConstCUDA(dblRK_Rational_1_18);

    const REAL_CUDA_ARRAY __rk_exp_0 = FusedMulAddCUDA(RK_Rational_1_18, MulCUDA(k1_gfsL, dt), y_n_gfsL);
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

  const size_t threads_in_x_dir = 32;
  const size_t threads_in_y_dir = 1;
  const size_t threads_in_z_dir = 1;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_1_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, y_n_gfs, next_y_input_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_1_gpu failure");
} // END FUNCTION rk_substep_1__launcher

/**
 * GPU Kernel: rk_substep_2_gpu.
 * GPU Kernel to compute RK substep 2.
 */
__global__ static void rk_substep_2_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k2_gfs, REAL *restrict y_n_gfs,
                                        REAL *restrict next_y_input_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k2_gfsL = k2_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dblRK_Rational_1_16 = 1.0 / 16.0;
    const REAL_CUDA_ARRAY RK_Rational_1_16 = ConstCUDA(dblRK_Rational_1_16);

    static const double dblRK_Rational_1_48 = 1.0 / 48.0;
    const REAL_CUDA_ARRAY RK_Rational_1_48 = ConstCUDA(dblRK_Rational_1_48);

    const REAL_CUDA_ARRAY __rk_exp_0 =
        FusedMulAddCUDA(RK_Rational_1_16, MulCUDA(k2_gfsL, dt), FusedMulAddCUDA(RK_Rational_1_48, MulCUDA(k1_gfsL, dt), y_n_gfsL));
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

  const size_t threads_in_x_dir = 32;
  const size_t threads_in_y_dir = 1;
  const size_t threads_in_z_dir = 1;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_2_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k2_gfs, y_n_gfs, next_y_input_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_2_gpu failure");
} // END FUNCTION rk_substep_2__launcher

/**
 * GPU Kernel: rk_substep_3_gpu.
 * GPU Kernel to compute RK substep 3.
 */
__global__ static void rk_substep_3_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k3_gfs, REAL *restrict y_n_gfs,
                                        REAL *restrict next_y_input_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k3_gfsL = k3_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dblRK_Rational_1_32 = 1.0 / 32.0;
    const REAL_CUDA_ARRAY RK_Rational_1_32 = ConstCUDA(dblRK_Rational_1_32);

    static const double dblRK_Rational_3_32 = 3.0 / 32.0;
    const REAL_CUDA_ARRAY RK_Rational_3_32 = ConstCUDA(dblRK_Rational_3_32);

    const REAL_CUDA_ARRAY __rk_exp_0 =
        FusedMulAddCUDA(RK_Rational_1_32, MulCUDA(k1_gfsL, dt), FusedMulAddCUDA(RK_Rational_3_32, MulCUDA(k3_gfsL, dt), y_n_gfsL));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_3_gpu

/**
 * Runge-Kutta function for substep 3.
 */
static void rk_substep_3__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict k3_gfs, REAL *restrict y_n_gfs,
                                   REAL *restrict next_y_input_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  const size_t threads_in_x_dir = 32;
  const size_t threads_in_y_dir = 1;
  const size_t threads_in_z_dir = 1;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_3_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k3_gfs, y_n_gfs, next_y_input_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_3_gpu failure");
} // END FUNCTION rk_substep_3__launcher

/**
 * GPU Kernel: rk_substep_4_gpu.
 * GPU Kernel to compute RK substep 4.
 */
__global__ static void rk_substep_4_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k3_gfs, REAL *restrict k4_gfs,
                                        REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k3_gfsL = k3_gfs[i];
    const REAL k4_gfsL = k4_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dbl_NegativeOne_ = -1.0;
    MAYBE_UNUSED const REAL_CUDA_ARRAY _NegativeOne_ = ConstCUDA(dbl_NegativeOne_);

    static const double dblRK_Rational_5_16 = 5.0 / 16.0;
    const REAL_CUDA_ARRAY RK_Rational_5_16 = ConstCUDA(dblRK_Rational_5_16);

    static const double dblRK_Rational_75_64 = 75.0 / 64.0;
    const REAL_CUDA_ARRAY RK_Rational_75_64 = ConstCUDA(dblRK_Rational_75_64);

    const REAL_CUDA_ARRAY __rk_exp_0 = FusedMulAddCUDA(RK_Rational_75_64, FusedMulSubCUDA(k4_gfsL, dt, MulCUDA(k3_gfsL, dt)),
                                                       FusedMulAddCUDA(RK_Rational_5_16, MulCUDA(k1_gfsL, dt), y_n_gfsL));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_4_gpu

/**
 * Runge-Kutta function for substep 4.
 */
static void rk_substep_4__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict k3_gfs, REAL *restrict k4_gfs,
                                   REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  const size_t threads_in_x_dir = 32;
  const size_t threads_in_y_dir = 1;
  const size_t threads_in_z_dir = 1;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_4_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k3_gfs, k4_gfs, y_n_gfs, next_y_input_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_4_gpu failure");
} // END FUNCTION rk_substep_4__launcher

/**
 * GPU Kernel: rk_substep_5_gpu.
 * GPU Kernel to compute RK substep 5.
 */
__global__ static void rk_substep_5_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k4_gfs, REAL *restrict k5_gfs,
                                        REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k4_gfsL = k4_gfs[i];
    const REAL k5_gfsL = k5_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dblRK_Rational_3_16 = 3.0 / 16.0;
    const REAL_CUDA_ARRAY RK_Rational_3_16 = ConstCUDA(dblRK_Rational_3_16);

    static const double dblRK_Rational_3_20 = 3.0 / 20.0;
    const REAL_CUDA_ARRAY RK_Rational_3_20 = ConstCUDA(dblRK_Rational_3_20);

    static const double dblRK_Rational_3_80 = 3.0 / 80.0;
    const REAL_CUDA_ARRAY RK_Rational_3_80 = ConstCUDA(dblRK_Rational_3_80);

    const REAL_CUDA_ARRAY __rk_exp_0 =
        FusedMulAddCUDA(RK_Rational_3_20, MulCUDA(k5_gfsL, dt),
                        FusedMulAddCUDA(RK_Rational_3_80, MulCUDA(k1_gfsL, dt), FusedMulAddCUDA(RK_Rational_3_16, MulCUDA(k4_gfsL, dt), y_n_gfsL)));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_5_gpu

/**
 * Runge-Kutta function for substep 5.
 */
static void rk_substep_5__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict k4_gfs, REAL *restrict k5_gfs,
                                   REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  const size_t threads_in_x_dir = 32;
  const size_t threads_in_y_dir = 1;
  const size_t threads_in_z_dir = 1;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_5_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k4_gfs, k5_gfs, y_n_gfs, next_y_input_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_5_gpu failure");
} // END FUNCTION rk_substep_5__launcher

/**
 * GPU Kernel: rk_substep_6_gpu.
 * GPU Kernel to compute RK substep 6.
 */
__global__ static void rk_substep_6_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k4_gfs, REAL *restrict k5_gfs,
                                        REAL *restrict k6_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k4_gfsL = k4_gfs[i];
    const REAL k5_gfsL = k5_gfs[i];
    const REAL k6_gfsL = k6_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dbl_NegativeOne_ = -1.0;
    MAYBE_UNUSED const REAL_CUDA_ARRAY _NegativeOne_ = ConstCUDA(dbl_NegativeOne_);

    static const double dblRK_Rational_23124283_1800000000 = 23124283.0 / 1800000000.0;
    const REAL_CUDA_ARRAY RK_Rational_23124283_1800000000 = ConstCUDA(dblRK_Rational_23124283_1800000000);

    static const double dblRK_Rational_28693883_1125000000 = 28693883.0 / 1125000000.0;
    const REAL_CUDA_ARRAY RK_Rational_28693883_1125000000 = ConstCUDA(dblRK_Rational_28693883_1125000000);

    static const double dblRK_Rational_29443841_614563906 = 29443841.0 / 614563906.0;
    const REAL_CUDA_ARRAY RK_Rational_29443841_614563906 = ConstCUDA(dblRK_Rational_29443841_614563906);

    static const double dblRK_Rational_77736538_692538347 = 77736538.0 / 692538347.0;
    const REAL_CUDA_ARRAY RK_Rational_77736538_692538347 = ConstCUDA(dblRK_Rational_77736538_692538347);

    const REAL_CUDA_ARRAY __rk_exp_0 =
        FusedMulAddCUDA(RK_Rational_23124283_1800000000, MulCUDA(k6_gfsL, dt),
                        FusedMulAddCUDA(RK_Rational_29443841_614563906, MulCUDA(k1_gfsL, dt),
                                        FusedMulAddCUDA(RK_Rational_77736538_692538347, MulCUDA(k4_gfsL, dt),
                                                        NegFusedMulAddCUDA(RK_Rational_28693883_1125000000, MulCUDA(k5_gfsL, dt), y_n_gfsL))));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_6_gpu

/**
 * Runge-Kutta function for substep 6.
 */
static void rk_substep_6__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict k4_gfs, REAL *restrict k5_gfs,
                                   REAL *restrict k6_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  const size_t threads_in_x_dir = 32;
  const size_t threads_in_y_dir = 1;
  const size_t threads_in_z_dir = 1;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_6_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k4_gfs, k5_gfs, k6_gfs, y_n_gfs, next_y_input_gfs,
                                                                                  dt);
  cudaCheckErrors(cudaKernel, "rk_substep_6_gpu failure");
} // END FUNCTION rk_substep_6__launcher

/**
 * GPU Kernel: rk_substep_7_gpu.
 * GPU Kernel to compute RK substep 7.
 */
__global__ static void rk_substep_7_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k4_gfs, REAL *restrict k5_gfs,
                                        REAL *restrict k6_gfs, REAL *restrict k7_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs,
                                        const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k4_gfsL = k4_gfs[i];
    const REAL k5_gfsL = k5_gfs[i];
    const REAL k6_gfsL = k6_gfs[i];
    const REAL k7_gfsL = k7_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dbl_NegativeOne_ = -1.0;
    MAYBE_UNUSED const REAL_CUDA_ARRAY _NegativeOne_ = ConstCUDA(dbl_NegativeOne_);

    static const double dblRK_Rational_16016141_946692911 = 16016141.0 / 946692911.0;
    const REAL_CUDA_ARRAY RK_Rational_16016141_946692911 = ConstCUDA(dblRK_Rational_16016141_946692911);

    static const double dblRK_Rational_180193667_1043307555 = 180193667.0 / 1043307555.0;
    const REAL_CUDA_ARRAY RK_Rational_180193667_1043307555 = ConstCUDA(dblRK_Rational_180193667_1043307555);

    static const double dblRK_Rational_22789713_633445777 = 22789713.0 / 633445777.0;
    const REAL_CUDA_ARRAY RK_Rational_22789713_633445777 = ConstCUDA(dblRK_Rational_22789713_633445777);

    static const double dblRK_Rational_545815736_2771057229 = 545815736.0 / 2771057229.0;
    const REAL_CUDA_ARRAY RK_Rational_545815736_2771057229 = ConstCUDA(dblRK_Rational_545815736_2771057229);

    static const double dblRK_Rational_61564180_158732637 = 61564180.0 / 158732637.0;
    const REAL_CUDA_ARRAY RK_Rational_61564180_158732637 = ConstCUDA(dblRK_Rational_61564180_158732637);

    const REAL_CUDA_ARRAY __rk_exp_0 = FusedMulAddCUDA(
        RK_Rational_545815736_2771057229, MulCUDA(k6_gfsL, dt),
        FusedMulAddCUDA(RK_Rational_16016141_946692911, MulCUDA(k1_gfsL, dt),
                        FusedMulAddCUDA(RK_Rational_22789713_633445777, MulCUDA(k5_gfsL, dt),
                                        FusedMulAddCUDA(RK_Rational_61564180_158732637, MulCUDA(k4_gfsL, dt),
                                                        NegFusedMulAddCUDA(RK_Rational_180193667_1043307555, MulCUDA(k7_gfsL, dt), y_n_gfsL)))));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_7_gpu

/**
 * Runge-Kutta function for substep 7.
 */
static void rk_substep_7__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict k4_gfs, REAL *restrict k5_gfs,
                                   REAL *restrict k6_gfs, REAL *restrict k7_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs,
                                   const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  const size_t threads_in_x_dir = 32;
  const size_t threads_in_y_dir = 1;
  const size_t threads_in_z_dir = 1;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_7_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k4_gfs, k5_gfs, k6_gfs, k7_gfs, y_n_gfs,
                                                                                  next_y_input_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_7_gpu failure");
} // END FUNCTION rk_substep_7__launcher

/**
 * GPU Kernel: rk_substep_8_gpu.
 * GPU Kernel to compute RK substep 8.
 */
__global__ static void rk_substep_8_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k4_gfs, REAL *restrict k5_gfs,
                                        REAL *restrict k6_gfs, REAL *restrict k7_gfs, REAL *restrict k8_gfs, REAL *restrict y_n_gfs,
                                        REAL *restrict next_y_input_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k4_gfsL = k4_gfs[i];
    const REAL k5_gfsL = k5_gfs[i];
    const REAL k6_gfsL = k6_gfs[i];
    const REAL k7_gfsL = k7_gfs[i];
    const REAL k8_gfsL = k8_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dbl_NegativeOne_ = -1.0;
    MAYBE_UNUSED const REAL_CUDA_ARRAY _NegativeOne_ = ConstCUDA(dbl_NegativeOne_);

    static const double dblRK_Rational_100302831_723423059 = 100302831.0 / 723423059.0;
    const REAL_CUDA_ARRAY RK_Rational_100302831_723423059 = ConstCUDA(dblRK_Rational_100302831_723423059);

    static const double dblRK_Rational_39632708_573591083 = 39632708.0 / 573591083.0;
    const REAL_CUDA_ARRAY RK_Rational_39632708_573591083 = ConstCUDA(dblRK_Rational_39632708_573591083);

    static const double dblRK_Rational_421739975_2616292301 = 421739975.0 / 2616292301.0;
    const REAL_CUDA_ARRAY RK_Rational_421739975_2616292301 = ConstCUDA(dblRK_Rational_421739975_2616292301);

    static const double dblRK_Rational_433636366_683701615 = 433636366.0 / 683701615.0;
    const REAL_CUDA_ARRAY RK_Rational_433636366_683701615 = ConstCUDA(dblRK_Rational_433636366_683701615);

    static const double dblRK_Rational_790204164_839813087 = 790204164.0 / 839813087.0;
    const REAL_CUDA_ARRAY RK_Rational_790204164_839813087 = ConstCUDA(dblRK_Rational_790204164_839813087);

    static const double dblRK_Rational_800635310_3783071287 = 800635310.0 / 3783071287.0;
    const REAL_CUDA_ARRAY RK_Rational_800635310_3783071287 = ConstCUDA(dblRK_Rational_800635310_3783071287);

    const REAL_CUDA_ARRAY __rk_exp_0 =
        FusedMulAddCUDA(RK_Rational_790204164_839813087, MulCUDA(k7_gfsL, dt),
                        FusedMulAddCUDA(RK_Rational_800635310_3783071287, MulCUDA(k8_gfsL, dt),
                                        FusedMulAddCUDA(RK_Rational_100302831_723423059, MulCUDA(k6_gfsL, dt),
                                                        FusedMulAddCUDA(RK_Rational_39632708_573591083, MulCUDA(k1_gfsL, dt),
                                                                        NegFusedMulAddCUDA(RK_Rational_433636366_683701615, MulCUDA(k4_gfsL, dt),
                                                                                           NegFusedMulAddCUDA(RK_Rational_421739975_2616292301,
                                                                                                              MulCUDA(k5_gfsL, dt), y_n_gfsL))))));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_8_gpu

/**
 * Runge-Kutta function for substep 8.
 */
static void rk_substep_8__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict k4_gfs, REAL *restrict k5_gfs,
                                   REAL *restrict k6_gfs, REAL *restrict k7_gfs, REAL *restrict k8_gfs, REAL *restrict y_n_gfs,
                                   REAL *restrict next_y_input_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  const size_t threads_in_x_dir = 32;
  const size_t threads_in_y_dir = 1;
  const size_t threads_in_z_dir = 1;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_8_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k4_gfs, k5_gfs, k6_gfs, k7_gfs, k8_gfs, y_n_gfs,
                                                                                  next_y_input_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_8_gpu failure");
} // END FUNCTION rk_substep_8__launcher

/**
 * GPU Kernel: rk_substep_9_gpu.
 * GPU Kernel to compute RK substep 9.
 */
__global__ static void rk_substep_9_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k4_gfs, REAL *restrict k5_gfs,
                                        REAL *restrict k6_gfs, REAL *restrict k7_gfs, REAL *restrict k8_gfs, REAL *restrict k9_gfs,
                                        REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k4_gfsL = k4_gfs[i];
    const REAL k5_gfsL = k5_gfs[i];
    const REAL k6_gfsL = k6_gfs[i];
    const REAL k7_gfsL = k7_gfs[i];
    const REAL k8_gfsL = k8_gfs[i];
    const REAL k9_gfsL = k9_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dbl_NegativeOne_ = -1.0;
    MAYBE_UNUSED const REAL_CUDA_ARRAY _NegativeOne_ = ConstCUDA(dbl_NegativeOne_);

    static const double dblRK_Rational_123872331_1001029789 = 123872331.0 / 1001029789.0;
    const REAL_CUDA_ARRAY RK_Rational_123872331_1001029789 = ConstCUDA(dblRK_Rational_123872331_1001029789);

    static const double dblRK_Rational_12992083_490766935 = 12992083.0 / 490766935.0;
    const REAL_CUDA_ARRAY RK_Rational_12992083_490766935 = ConstCUDA(dblRK_Rational_12992083_490766935);

    static const double dblRK_Rational_246121993_1340847787 = 246121993.0 / 1340847787.0;
    const REAL_CUDA_ARRAY RK_Rational_246121993_1340847787 = ConstCUDA(dblRK_Rational_246121993_1340847787);

    static const double dblRK_Rational_309121744_1061227803 = 309121744.0 / 1061227803.0;
    const REAL_CUDA_ARRAY RK_Rational_309121744_1061227803 = ConstCUDA(dblRK_Rational_309121744_1061227803);

    static const double dblRK_Rational_37695042795_15268766246 = 37695042795.0 / 15268766246.0;
    const REAL_CUDA_ARRAY RK_Rational_37695042795_15268766246 = ConstCUDA(dblRK_Rational_37695042795_15268766246);

    static const double dblRK_Rational_393006217_1396673457 = 393006217.0 / 1396673457.0;
    const REAL_CUDA_ARRAY RK_Rational_393006217_1396673457 = ConstCUDA(dblRK_Rational_393006217_1396673457);

    static const double dblRK_Rational_6005943493_2108947869 = 6005943493.0 / 2108947869.0;
    const REAL_CUDA_ARRAY RK_Rational_6005943493_2108947869 = ConstCUDA(dblRK_Rational_6005943493_2108947869);

    const REAL_CUDA_ARRAY __rk_exp_0 = NegFusedMulAddCUDA(
        RK_Rational_309121744_1061227803, MulCUDA(k5_gfsL, dt),
        FusedMulAddCUDA(RK_Rational_393006217_1396673457, MulCUDA(k8_gfsL, dt),
                        FusedMulAddCUDA(RK_Rational_6005943493_2108947869, MulCUDA(k7_gfsL, dt),
                                        FusedMulAddCUDA(RK_Rational_123872331_1001029789, MulCUDA(k9_gfsL, dt),
                                                        FusedMulAddCUDA(RK_Rational_246121993_1340847787, MulCUDA(k1_gfsL, dt),
                                                                        NegFusedMulAddCUDA(RK_Rational_37695042795_15268766246, MulCUDA(k4_gfsL, dt),
                                                                                           NegFusedMulAddCUDA(RK_Rational_12992083_490766935,
                                                                                                              MulCUDA(k6_gfsL, dt), y_n_gfsL)))))));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_9_gpu

/**
 * Runge-Kutta function for substep 9.
 */
static void rk_substep_9__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict k4_gfs, REAL *restrict k5_gfs,
                                   REAL *restrict k6_gfs, REAL *restrict k7_gfs, REAL *restrict k8_gfs, REAL *restrict k9_gfs, REAL *restrict y_n_gfs,
                                   REAL *restrict next_y_input_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  const size_t threads_in_x_dir = 32;
  const size_t threads_in_y_dir = 1;
  const size_t threads_in_z_dir = 1;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_9_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k4_gfs, k5_gfs, k6_gfs, k7_gfs, k8_gfs, k9_gfs,
                                                                                  y_n_gfs, next_y_input_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_9_gpu failure");
} // END FUNCTION rk_substep_9__launcher

/**
 * GPU Kernel: rk_substep_10_gpu.
 * GPU Kernel to compute RK substep 10.
 */
__global__ static void rk_substep_10_gpu(const size_t streamid, REAL *restrict k10_gfs, REAL *restrict k1_gfs, REAL *restrict k4_gfs,
                                         REAL *restrict k5_gfs, REAL *restrict k6_gfs, REAL *restrict k7_gfs, REAL *restrict k8_gfs,
                                         REAL *restrict k9_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k10_gfsL = k10_gfs[i];
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k4_gfsL = k4_gfs[i];
    const REAL k5_gfsL = k5_gfs[i];
    const REAL k6_gfsL = k6_gfs[i];
    const REAL k7_gfsL = k7_gfs[i];
    const REAL k8_gfsL = k8_gfs[i];
    const REAL k9_gfsL = k9_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dbl_NegativeOne_ = -1.0;
    MAYBE_UNUSED const REAL_CUDA_ARRAY _NegativeOne_ = ConstCUDA(dbl_NegativeOne_);

    static const double dblRK_Rational_1028468189_846180014 = 1028468189.0 / 846180014.0;
    const REAL_CUDA_ARRAY RK_Rational_1028468189_846180014 = ConstCUDA(dblRK_Rational_1028468189_846180014);

    static const double dblRK_Rational_10304129995_1701304382 = 10304129995.0 / 1701304382.0;
    const REAL_CUDA_ARRAY RK_Rational_10304129995_1701304382 = ConstCUDA(dblRK_Rational_10304129995_1701304382);

    static const double dblRK_Rational_1311729495_1432422823 = 1311729495.0 / 1432422823.0;
    const REAL_CUDA_ARRAY RK_Rational_1311729495_1432422823 = ConstCUDA(dblRK_Rational_1311729495_1432422823);

    static const double dblRK_Rational_15336726248_1032824649 = 15336726248.0 / 1032824649.0;
    const REAL_CUDA_ARRAY RK_Rational_15336726248_1032824649 = ConstCUDA(dblRK_Rational_15336726248_1032824649);

    static const double dblRK_Rational_3065993473_597172653 = 3065993473.0 / 597172653.0;
    const REAL_CUDA_ARRAY RK_Rational_3065993473_597172653 = ConstCUDA(dblRK_Rational_3065993473_597172653);

    static const double dblRK_Rational_45442868181_3398467696 = 45442868181.0 / 3398467696.0;
    const REAL_CUDA_ARRAY RK_Rational_45442868181_3398467696 = ConstCUDA(dblRK_Rational_45442868181_3398467696);

    static const double dblRK_Rational_48777925059_3047939560 = 48777925059.0 / 3047939560.0;
    const REAL_CUDA_ARRAY RK_Rational_48777925059_3047939560 = ConstCUDA(dblRK_Rational_48777925059_3047939560);

    static const double dblRK_Rational_8478235783_508512852 = 8478235783.0 / 508512852.0;
    const REAL_CUDA_ARRAY RK_Rational_8478235783_508512852 = ConstCUDA(dblRK_Rational_8478235783_508512852);

    const REAL_CUDA_ARRAY __rk_exp_0 = NegFusedMulAddCUDA(
        RK_Rational_10304129995_1701304382, MulCUDA(k6_gfsL, dt),
        NegFusedMulAddCUDA(
            RK_Rational_45442868181_3398467696, MulCUDA(k9_gfsL, dt),
            FusedMulAddCUDA(
                RK_Rational_3065993473_597172653, MulCUDA(k10_gfsL, dt),
                FusedMulAddCUDA(RK_Rational_8478235783_508512852, MulCUDA(k4_gfsL, dt),
                                FusedMulAddCUDA(RK_Rational_1311729495_1432422823, MulCUDA(k5_gfsL, dt),
                                                FusedMulAddCUDA(RK_Rational_15336726248_1032824649, MulCUDA(k8_gfsL, dt),
                                                                NegFusedMulAddCUDA(RK_Rational_48777925059_3047939560, MulCUDA(k7_gfsL, dt),
                                                                                   NegFusedMulAddCUDA(RK_Rational_1028468189_846180014,
                                                                                                      MulCUDA(k1_gfsL, dt), y_n_gfsL))))))));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_10_gpu

/**
 * Runge-Kutta function for substep 10.
 */
static void rk_substep_10__launcher(params_struct *restrict params, REAL *restrict k10_gfs, REAL *restrict k1_gfs, REAL *restrict k4_gfs,
                                    REAL *restrict k5_gfs, REAL *restrict k6_gfs, REAL *restrict k7_gfs, REAL *restrict k8_gfs, REAL *restrict k9_gfs,
                                    REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  const size_t threads_in_x_dir = 32;
  const size_t threads_in_y_dir = 1;
  const size_t threads_in_z_dir = 1;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_10_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k10_gfs, k1_gfs, k4_gfs, k5_gfs, k6_gfs, k7_gfs, k8_gfs,
                                                                                   k9_gfs, y_n_gfs, next_y_input_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_10_gpu failure");
} // END FUNCTION rk_substep_10__launcher

/**
 * GPU Kernel: rk_substep_11_gpu.
 * GPU Kernel to compute RK substep 11.
 */
__global__ static void rk_substep_11_gpu(const size_t streamid, REAL *restrict k10_gfs, REAL *restrict k11_gfs, REAL *restrict k1_gfs,
                                         REAL *restrict k4_gfs, REAL *restrict k5_gfs, REAL *restrict k6_gfs, REAL *restrict k7_gfs,
                                         REAL *restrict k8_gfs, REAL *restrict k9_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs,
                                         const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k10_gfsL = k10_gfs[i];
    const REAL k11_gfsL = k11_gfs[i];
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k4_gfsL = k4_gfs[i];
    const REAL k5_gfsL = k5_gfs[i];
    const REAL k6_gfsL = k6_gfs[i];
    const REAL k7_gfsL = k7_gfs[i];
    const REAL k8_gfsL = k8_gfs[i];
    const REAL k9_gfsL = k9_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dbl_NegativeOne_ = -1.0;
    MAYBE_UNUSED const REAL_CUDA_ARRAY _NegativeOne_ = ConstCUDA(dbl_NegativeOne_);

    static const double dblRK_Rational_185892177_718116043 = 185892177.0 / 718116043.0;
    const REAL_CUDA_ARRAY RK_Rational_185892177_718116043 = ConstCUDA(dblRK_Rational_185892177_718116043);

    static const double dblRK_Rational_3185094517_667107341 = 3185094517.0 / 667107341.0;
    const REAL_CUDA_ARRAY RK_Rational_3185094517_667107341 = ConstCUDA(dblRK_Rational_3185094517_667107341);

    static const double dblRK_Rational_3962137247_1805957418 = 3962137247.0 / 1805957418.0;
    const REAL_CUDA_ARRAY RK_Rational_3962137247_1805957418 = ConstCUDA(dblRK_Rational_3962137247_1805957418);

    static const double dblRK_Rational_4093664535_808688257 = 4093664535.0 / 808688257.0;
    const REAL_CUDA_ARRAY RK_Rational_4093664535_808688257 = ConstCUDA(dblRK_Rational_4093664535_808688257);

    static const double dblRK_Rational_477755414_1098053517 = 477755414.0 / 1098053517.0;
    const REAL_CUDA_ARRAY RK_Rational_477755414_1098053517 = ConstCUDA(dblRK_Rational_477755414_1098053517);

    static const double dblRK_Rational_5232866602_850066563 = 5232866602.0 / 850066563.0;
    const REAL_CUDA_ARRAY RK_Rational_5232866602_850066563 = ConstCUDA(dblRK_Rational_5232866602_850066563);

    static const double dblRK_Rational_5731566787_1027545527 = 5731566787.0 / 1027545527.0;
    const REAL_CUDA_ARRAY RK_Rational_5731566787_1027545527 = ConstCUDA(dblRK_Rational_5731566787_1027545527);

    static const double dblRK_Rational_65686358_487910083 = 65686358.0 / 487910083.0;
    const REAL_CUDA_ARRAY RK_Rational_65686358_487910083 = ConstCUDA(dblRK_Rational_65686358_487910083);

    static const double dblRK_Rational_703635378_230739211 = 703635378.0 / 230739211.0;
    const REAL_CUDA_ARRAY RK_Rational_703635378_230739211 = ConstCUDA(dblRK_Rational_703635378_230739211);

    const REAL_CUDA_ARRAY __rk_exp_0 = NegFusedMulAddCUDA(
        RK_Rational_477755414_1098053517, MulCUDA(k5_gfsL, dt),
        FusedMulAddCUDA(
            RK_Rational_65686358_487910083, MulCUDA(k11_gfsL, dt),
            NegFusedMulAddCUDA(
                RK_Rational_4093664535_808688257, MulCUDA(k9_gfsL, dt),
                FusedMulAddCUDA(
                    RK_Rational_5232866602_850066563, MulCUDA(k8_gfsL, dt),
                    FusedMulAddCUDA(RK_Rational_5731566787_1027545527, MulCUDA(k7_gfsL, dt),
                                    FusedMulAddCUDA(RK_Rational_185892177_718116043, MulCUDA(k1_gfsL, dt),
                                                    FusedMulAddCUDA(RK_Rational_3962137247_1805957418, MulCUDA(k10_gfsL, dt),
                                                                    NegFusedMulAddCUDA(RK_Rational_703635378_230739211, MulCUDA(k6_gfsL, dt),
                                                                                       NegFusedMulAddCUDA(RK_Rational_3185094517_667107341,
                                                                                                          MulCUDA(k4_gfsL, dt), y_n_gfsL)))))))));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_11_gpu

/**
 * Runge-Kutta function for substep 11.
 */
static void rk_substep_11__launcher(params_struct *restrict params, REAL *restrict k10_gfs, REAL *restrict k11_gfs, REAL *restrict k1_gfs,
                                    REAL *restrict k4_gfs, REAL *restrict k5_gfs, REAL *restrict k6_gfs, REAL *restrict k7_gfs, REAL *restrict k8_gfs,
                                    REAL *restrict k9_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  const size_t threads_in_x_dir = 32;
  const size_t threads_in_y_dir = 1;
  const size_t threads_in_z_dir = 1;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_11_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k10_gfs, k11_gfs, k1_gfs, k4_gfs, k5_gfs, k6_gfs, k7_gfs,
                                                                                   k8_gfs, k9_gfs, y_n_gfs, next_y_input_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_11_gpu failure");
} // END FUNCTION rk_substep_11__launcher

/**
 * GPU Kernel: rk_substep_12_gpu.
 * GPU Kernel to compute RK substep 12.
 */
__global__ static void rk_substep_12_gpu(const size_t streamid, REAL *restrict k10_gfs, REAL *restrict k11_gfs, REAL *restrict k1_gfs,
                                         REAL *restrict k4_gfs, REAL *restrict k5_gfs, REAL *restrict k6_gfs, REAL *restrict k7_gfs,
                                         REAL *restrict k8_gfs, REAL *restrict k9_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs,
                                         const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k10_gfsL = k10_gfs[i];
    const REAL k11_gfsL = k11_gfs[i];
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k4_gfsL = k4_gfs[i];
    const REAL k5_gfsL = k5_gfs[i];
    const REAL k6_gfsL = k6_gfs[i];
    const REAL k7_gfsL = k7_gfs[i];
    const REAL k8_gfsL = k8_gfs[i];
    const REAL k9_gfsL = k9_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dbl_NegativeOne_ = -1.0;
    MAYBE_UNUSED const REAL_CUDA_ARRAY _NegativeOne_ = ConstCUDA(dbl_NegativeOne_);

    static const double dblRK_Rational_11173962825_925320556 = 11173962825.0 / 925320556.0;
    const REAL_CUDA_ARRAY RK_Rational_11173962825_925320556 = ConstCUDA(dblRK_Rational_11173962825_925320556);

    static const double dblRK_Rational_13158990841_6184727034 = 13158990841.0 / 6184727034.0;
    const REAL_CUDA_ARRAY RK_Rational_13158990841_6184727034 = ConstCUDA(dblRK_Rational_13158990841_6184727034);

    static const double dblRK_Rational_160528059_685178525 = 160528059.0 / 685178525.0;
    const REAL_CUDA_ARRAY RK_Rational_160528059_685178525 = ConstCUDA(dblRK_Rational_160528059_685178525);

    static const double dblRK_Rational_248638103_1413531060 = 248638103.0 / 1413531060.0;
    const REAL_CUDA_ARRAY RK_Rational_248638103_1413531060 = ConstCUDA(dblRK_Rational_248638103_1413531060);

    static const double dblRK_Rational_3936647629_1978049680 = 3936647629.0 / 1978049680.0;
    const REAL_CUDA_ARRAY RK_Rational_3936647629_1978049680 = ConstCUDA(dblRK_Rational_3936647629_1978049680);

    static const double dblRK_Rational_403863854_491063109 = 403863854.0 / 491063109.0;
    const REAL_CUDA_ARRAY RK_Rational_403863854_491063109 = ConstCUDA(dblRK_Rational_403863854_491063109);

    static const double dblRK_Rational_411421997_543043805 = 411421997.0 / 543043805.0;
    const REAL_CUDA_ARRAY RK_Rational_411421997_543043805 = ConstCUDA(dblRK_Rational_411421997_543043805);

    static const double dblRK_Rational_5068492393_434740067 = 5068492393.0 / 434740067.0;
    const REAL_CUDA_ARRAY RK_Rational_5068492393_434740067 = ConstCUDA(dblRK_Rational_5068492393_434740067);

    static const double dblRK_Rational_652783627_914296604 = 652783627.0 / 914296604.0;
    const REAL_CUDA_ARRAY RK_Rational_652783627_914296604 = ConstCUDA(dblRK_Rational_652783627_914296604);

    const REAL_CUDA_ARRAY __rk_exp_0 = NegFusedMulAddCUDA(
        RK_Rational_411421997_543043805, MulCUDA(k5_gfsL, dt),
        FusedMulAddCUDA(
            RK_Rational_652783627_914296604, MulCUDA(k6_gfsL, dt),
            NegFusedMulAddCUDA(
                RK_Rational_160528059_685178525, MulCUDA(k10_gfsL, dt),
                FusedMulAddCUDA(
                    RK_Rational_3936647629_1978049680, MulCUDA(k9_gfsL, dt),
                    FusedMulAddCUDA(RK_Rational_403863854_491063109, MulCUDA(k1_gfsL, dt),
                                    FusedMulAddCUDA(RK_Rational_11173962825_925320556, MulCUDA(k7_gfsL, dt),
                                                    FusedMulAddCUDA(RK_Rational_248638103_1413531060, MulCUDA(k11_gfsL, dt),
                                                                    NegFusedMulAddCUDA(RK_Rational_5068492393_434740067, MulCUDA(k4_gfsL, dt),
                                                                                       NegFusedMulAddCUDA(RK_Rational_13158990841_6184727034,
                                                                                                          MulCUDA(k8_gfsL, dt), y_n_gfsL)))))))));
    WriteCUDA(&next_y_input_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_12_gpu

/**
 * Runge-Kutta function for substep 12.
 */
static void rk_substep_12__launcher(params_struct *restrict params, REAL *restrict k10_gfs, REAL *restrict k11_gfs, REAL *restrict k1_gfs,
                                    REAL *restrict k4_gfs, REAL *restrict k5_gfs, REAL *restrict k6_gfs, REAL *restrict k7_gfs, REAL *restrict k8_gfs,
                                    REAL *restrict k9_gfs, REAL *restrict y_n_gfs, REAL *restrict next_y_input_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  const size_t threads_in_x_dir = 32;
  const size_t threads_in_y_dir = 1;
  const size_t threads_in_z_dir = 1;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_12_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k10_gfs, k11_gfs, k1_gfs, k4_gfs, k5_gfs, k6_gfs, k7_gfs,
                                                                                   k8_gfs, k9_gfs, y_n_gfs, next_y_input_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_12_gpu failure");
} // END FUNCTION rk_substep_12__launcher

/**
 * GPU Kernel: rk_substep_13_gpu.
 * GPU Kernel to compute RK substep 13.
 */
__global__ static void rk_substep_13_gpu(const size_t streamid, REAL *restrict k10_gfs, REAL *restrict k11_gfs, REAL *restrict k12_gfs,
                                         REAL *restrict k13_gfs, REAL *restrict k1_gfs, REAL *restrict k6_gfs, REAL *restrict k7_gfs,
                                         REAL *restrict k8_gfs, REAL *restrict k9_gfs, REAL *restrict y_n_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k10_gfsL = k10_gfs[i];
    const REAL k11_gfsL = k11_gfs[i];
    const REAL k12_gfsL = k12_gfs[i];
    const REAL k13_gfsL = k13_gfs[i];
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k6_gfsL = k6_gfs[i];
    const REAL k7_gfsL = k7_gfs[i];
    const REAL k8_gfsL = k8_gfs[i];
    const REAL k9_gfsL = k9_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static const double dbl_NegativeOne_ = -1.0;
    MAYBE_UNUSED const REAL_CUDA_ARRAY _NegativeOne_ = ConstCUDA(dbl_NegativeOne_);

    static const double dblRK_Rational_1041891430_1371343529 = 1041891430.0 / 1371343529.0;
    const REAL_CUDA_ARRAY RK_Rational_1041891430_1371343529 = ConstCUDA(dblRK_Rational_1041891430_1371343529);

    static const double dblRK_Rational_118820643_751138087 = 118820643.0 / 751138087.0;
    const REAL_CUDA_ARRAY RK_Rational_118820643_751138087 = ConstCUDA(dblRK_Rational_118820643_751138087);

    static const double dblRK_Rational_14005451_335480064 = 14005451.0 / 335480064.0;
    const REAL_CUDA_ARRAY RK_Rational_14005451_335480064 = ConstCUDA(dblRK_Rational_14005451_335480064);

    static const double dblRK_Rational_181606767_758867731 = 181606767.0 / 758867731.0;
    const REAL_CUDA_ARRAY RK_Rational_181606767_758867731 = ConstCUDA(dblRK_Rational_181606767_758867731);

    static const double dblRK_Rational_1_4 = 1.0 / 4.0;
    const REAL_CUDA_ARRAY RK_Rational_1_4 = ConstCUDA(dblRK_Rational_1_4);

    static const double dblRK_Rational_528747749_2220607170 = 528747749.0 / 2220607170.0;
    const REAL_CUDA_ARRAY RK_Rational_528747749_2220607170 = ConstCUDA(dblRK_Rational_528747749_2220607170);

    static const double dblRK_Rational_561292985_797845732 = 561292985.0 / 797845732.0;
    const REAL_CUDA_ARRAY RK_Rational_561292985_797845732 = ConstCUDA(dblRK_Rational_561292985_797845732);

    static const double dblRK_Rational_59238493_1068277825 = 59238493.0 / 1068277825.0;
    const REAL_CUDA_ARRAY RK_Rational_59238493_1068277825 = ConstCUDA(dblRK_Rational_59238493_1068277825);

    static const double dblRK_Rational_760417239_1151165299 = 760417239.0 / 1151165299.0;
    const REAL_CUDA_ARRAY RK_Rational_760417239_1151165299 = ConstCUDA(dblRK_Rational_760417239_1151165299);

    const REAL_CUDA_ARRAY __rk_exp_0 = NegFusedMulAddCUDA(
        RK_Rational_528747749_2220607170, MulCUDA(k12_gfsL, dt),
        FusedMulAddCUDA(
            RK_Rational_561292985_797845732, MulCUDA(k8_gfsL, dt),
            FusedMulAddCUDA(
                RK_Rational_760417239_1151165299, MulCUDA(k10_gfsL, dt),
                FusedMulAddCUDA(
                    RK_Rational_181606767_758867731, MulCUDA(k7_gfsL, dt),
                    FusedMulAddCUDA(RK_Rational_1_4, MulCUDA(k13_gfsL, dt),
                                    FusedMulAddCUDA(RK_Rational_118820643_751138087, MulCUDA(k11_gfsL, dt),
                                                    FusedMulAddCUDA(RK_Rational_14005451_335480064, MulCUDA(k1_gfsL, dt),
                                                                    NegFusedMulAddCUDA(RK_Rational_59238493_1068277825, MulCUDA(k6_gfsL, dt),
                                                                                       NegFusedMulAddCUDA(RK_Rational_1041891430_1371343529,
                                                                                                          MulCUDA(k9_gfsL, dt), y_n_gfsL)))))))));
    WriteCUDA(&y_n_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_13_gpu

/**
 * Runge-Kutta function for substep 13.
 */
static void rk_substep_13__launcher(params_struct *restrict params, REAL *restrict k10_gfs, REAL *restrict k11_gfs, REAL *restrict k12_gfs,
                                    REAL *restrict k13_gfs, REAL *restrict k1_gfs, REAL *restrict k6_gfs, REAL *restrict k7_gfs,
                                    REAL *restrict k8_gfs, REAL *restrict k9_gfs, REAL *restrict y_n_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  const size_t threads_in_x_dir = 32;
  const size_t threads_in_y_dir = 1;
  const size_t threads_in_z_dir = 1;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_13_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k10_gfs, k11_gfs, k12_gfs, k13_gfs, k1_gfs, k6_gfs,
                                                                                   k7_gfs, k8_gfs, k9_gfs, y_n_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_13_gpu failure");
} // END FUNCTION rk_substep_13__launcher

/**
 * Method of Lines (MoL) for "DP8" method: Step forward one full timestep.
 *
 */
void MoL_step_forward_in_time(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {

  // C code implementation of -={ DP8 }=- Method of Lines timestepping.

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
    MAYBE_UNUSED REAL *restrict k8_gfs = griddata[grid].gridfuncs.k8_gfs;
    MAYBE_UNUSED REAL *restrict k9_gfs = griddata[grid].gridfuncs.k9_gfs;
    MAYBE_UNUSED REAL *restrict k10_gfs = griddata[grid].gridfuncs.k10_gfs;
    MAYBE_UNUSED REAL *restrict k11_gfs = griddata[grid].gridfuncs.k11_gfs;
    MAYBE_UNUSED REAL *restrict k12_gfs = griddata[grid].gridfuncs.k12_gfs;
    MAYBE_UNUSED REAL *restrict k13_gfs = griddata[grid].gridfuncs.k13_gfs;
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
    commondata->time = time_start + 5.55555555555555525e-02 * commondata->dt;
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
    MAYBE_UNUSED REAL *restrict k8_gfs = griddata[grid].gridfuncs.k8_gfs;
    MAYBE_UNUSED REAL *restrict k9_gfs = griddata[grid].gridfuncs.k9_gfs;
    MAYBE_UNUSED REAL *restrict k10_gfs = griddata[grid].gridfuncs.k10_gfs;
    MAYBE_UNUSED REAL *restrict k11_gfs = griddata[grid].gridfuncs.k11_gfs;
    MAYBE_UNUSED REAL *restrict k12_gfs = griddata[grid].gridfuncs.k12_gfs;
    MAYBE_UNUSED REAL *restrict k13_gfs = griddata[grid].gridfuncs.k13_gfs;
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
    commondata->time = time_start + 8.33333333333333287e-02 * commondata->dt;
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
    MAYBE_UNUSED REAL *restrict k8_gfs = griddata[grid].gridfuncs.k8_gfs;
    MAYBE_UNUSED REAL *restrict k9_gfs = griddata[grid].gridfuncs.k9_gfs;
    MAYBE_UNUSED REAL *restrict k10_gfs = griddata[grid].gridfuncs.k10_gfs;
    MAYBE_UNUSED REAL *restrict k11_gfs = griddata[grid].gridfuncs.k11_gfs;
    MAYBE_UNUSED REAL *restrict k12_gfs = griddata[grid].gridfuncs.k12_gfs;
    MAYBE_UNUSED REAL *restrict k13_gfs = griddata[grid].gridfuncs.k13_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k3_gfs);
    rk_substep_3__launcher(params, k1_gfs, k3_gfs, y_n_gfs, next_y_input_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k3 substep }=-

  // -={ START k4 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 1.25000000000000000e-01 * commondata->dt;
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
    MAYBE_UNUSED REAL *restrict k8_gfs = griddata[grid].gridfuncs.k8_gfs;
    MAYBE_UNUSED REAL *restrict k9_gfs = griddata[grid].gridfuncs.k9_gfs;
    MAYBE_UNUSED REAL *restrict k10_gfs = griddata[grid].gridfuncs.k10_gfs;
    MAYBE_UNUSED REAL *restrict k11_gfs = griddata[grid].gridfuncs.k11_gfs;
    MAYBE_UNUSED REAL *restrict k12_gfs = griddata[grid].gridfuncs.k12_gfs;
    MAYBE_UNUSED REAL *restrict k13_gfs = griddata[grid].gridfuncs.k13_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k4_gfs);
    rk_substep_4__launcher(params, k1_gfs, k3_gfs, k4_gfs, y_n_gfs, next_y_input_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k4 substep }=-

  // -={ START k5 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 3.12500000000000000e-01 * commondata->dt;
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
    MAYBE_UNUSED REAL *restrict k8_gfs = griddata[grid].gridfuncs.k8_gfs;
    MAYBE_UNUSED REAL *restrict k9_gfs = griddata[grid].gridfuncs.k9_gfs;
    MAYBE_UNUSED REAL *restrict k10_gfs = griddata[grid].gridfuncs.k10_gfs;
    MAYBE_UNUSED REAL *restrict k11_gfs = griddata[grid].gridfuncs.k11_gfs;
    MAYBE_UNUSED REAL *restrict k12_gfs = griddata[grid].gridfuncs.k12_gfs;
    MAYBE_UNUSED REAL *restrict k13_gfs = griddata[grid].gridfuncs.k13_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k5_gfs);
    rk_substep_5__launcher(params, k1_gfs, k4_gfs, k5_gfs, y_n_gfs, next_y_input_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k5 substep }=-

  // -={ START k6 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 3.75000000000000000e-01 * commondata->dt;
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
    MAYBE_UNUSED REAL *restrict k8_gfs = griddata[grid].gridfuncs.k8_gfs;
    MAYBE_UNUSED REAL *restrict k9_gfs = griddata[grid].gridfuncs.k9_gfs;
    MAYBE_UNUSED REAL *restrict k10_gfs = griddata[grid].gridfuncs.k10_gfs;
    MAYBE_UNUSED REAL *restrict k11_gfs = griddata[grid].gridfuncs.k11_gfs;
    MAYBE_UNUSED REAL *restrict k12_gfs = griddata[grid].gridfuncs.k12_gfs;
    MAYBE_UNUSED REAL *restrict k13_gfs = griddata[grid].gridfuncs.k13_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k6_gfs);
    rk_substep_6__launcher(params, k1_gfs, k4_gfs, k5_gfs, k6_gfs, y_n_gfs, next_y_input_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k6 substep }=-

  // -={ START k7 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 1.47499999999999992e-01 * commondata->dt;
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
    MAYBE_UNUSED REAL *restrict k8_gfs = griddata[grid].gridfuncs.k8_gfs;
    MAYBE_UNUSED REAL *restrict k9_gfs = griddata[grid].gridfuncs.k9_gfs;
    MAYBE_UNUSED REAL *restrict k10_gfs = griddata[grid].gridfuncs.k10_gfs;
    MAYBE_UNUSED REAL *restrict k11_gfs = griddata[grid].gridfuncs.k11_gfs;
    MAYBE_UNUSED REAL *restrict k12_gfs = griddata[grid].gridfuncs.k12_gfs;
    MAYBE_UNUSED REAL *restrict k13_gfs = griddata[grid].gridfuncs.k13_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k7_gfs);
    rk_substep_7__launcher(params, k1_gfs, k4_gfs, k5_gfs, k6_gfs, k7_gfs, y_n_gfs, next_y_input_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k7 substep }=-

  // -={ START k8 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 4.65000000000000024e-01 * commondata->dt;
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
    MAYBE_UNUSED REAL *restrict k8_gfs = griddata[grid].gridfuncs.k8_gfs;
    MAYBE_UNUSED REAL *restrict k9_gfs = griddata[grid].gridfuncs.k9_gfs;
    MAYBE_UNUSED REAL *restrict k10_gfs = griddata[grid].gridfuncs.k10_gfs;
    MAYBE_UNUSED REAL *restrict k11_gfs = griddata[grid].gridfuncs.k11_gfs;
    MAYBE_UNUSED REAL *restrict k12_gfs = griddata[grid].gridfuncs.k12_gfs;
    MAYBE_UNUSED REAL *restrict k13_gfs = griddata[grid].gridfuncs.k13_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k8_gfs);
    rk_substep_8__launcher(params, k1_gfs, k4_gfs, k5_gfs, k6_gfs, k7_gfs, k8_gfs, y_n_gfs, next_y_input_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k8 substep }=-

  // -={ START k9 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 5.64865451382259520e-01 * commondata->dt;
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
    MAYBE_UNUSED REAL *restrict k8_gfs = griddata[grid].gridfuncs.k8_gfs;
    MAYBE_UNUSED REAL *restrict k9_gfs = griddata[grid].gridfuncs.k9_gfs;
    MAYBE_UNUSED REAL *restrict k10_gfs = griddata[grid].gridfuncs.k10_gfs;
    MAYBE_UNUSED REAL *restrict k11_gfs = griddata[grid].gridfuncs.k11_gfs;
    MAYBE_UNUSED REAL *restrict k12_gfs = griddata[grid].gridfuncs.k12_gfs;
    MAYBE_UNUSED REAL *restrict k13_gfs = griddata[grid].gridfuncs.k13_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k9_gfs);
    rk_substep_9__launcher(params, k1_gfs, k4_gfs, k5_gfs, k6_gfs, k7_gfs, k8_gfs, k9_gfs, y_n_gfs, next_y_input_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k9 substep }=-

  // -={ START k10 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 6.50000000000000022e-01 * commondata->dt;
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
    MAYBE_UNUSED REAL *restrict k8_gfs = griddata[grid].gridfuncs.k8_gfs;
    MAYBE_UNUSED REAL *restrict k9_gfs = griddata[grid].gridfuncs.k9_gfs;
    MAYBE_UNUSED REAL *restrict k10_gfs = griddata[grid].gridfuncs.k10_gfs;
    MAYBE_UNUSED REAL *restrict k11_gfs = griddata[grid].gridfuncs.k11_gfs;
    MAYBE_UNUSED REAL *restrict k12_gfs = griddata[grid].gridfuncs.k12_gfs;
    MAYBE_UNUSED REAL *restrict k13_gfs = griddata[grid].gridfuncs.k13_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k10_gfs);
    rk_substep_10__launcher(params, k10_gfs, k1_gfs, k4_gfs, k5_gfs, k6_gfs, k7_gfs, k8_gfs, k9_gfs, y_n_gfs, next_y_input_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k10 substep }=-

  // -={ START k11 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 9.24656277640504398e-01 * commondata->dt;
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
    MAYBE_UNUSED REAL *restrict k8_gfs = griddata[grid].gridfuncs.k8_gfs;
    MAYBE_UNUSED REAL *restrict k9_gfs = griddata[grid].gridfuncs.k9_gfs;
    MAYBE_UNUSED REAL *restrict k10_gfs = griddata[grid].gridfuncs.k10_gfs;
    MAYBE_UNUSED REAL *restrict k11_gfs = griddata[grid].gridfuncs.k11_gfs;
    MAYBE_UNUSED REAL *restrict k12_gfs = griddata[grid].gridfuncs.k12_gfs;
    MAYBE_UNUSED REAL *restrict k13_gfs = griddata[grid].gridfuncs.k13_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k11_gfs);
    rk_substep_11__launcher(params, k10_gfs, k11_gfs, k1_gfs, k4_gfs, k5_gfs, k6_gfs, k7_gfs, k8_gfs, k9_gfs, y_n_gfs, next_y_input_gfs,
                            commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k11 substep }=-

  // -={ START k12 substep }=-
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
    MAYBE_UNUSED REAL *restrict k8_gfs = griddata[grid].gridfuncs.k8_gfs;
    MAYBE_UNUSED REAL *restrict k9_gfs = griddata[grid].gridfuncs.k9_gfs;
    MAYBE_UNUSED REAL *restrict k10_gfs = griddata[grid].gridfuncs.k10_gfs;
    MAYBE_UNUSED REAL *restrict k11_gfs = griddata[grid].gridfuncs.k11_gfs;
    MAYBE_UNUSED REAL *restrict k12_gfs = griddata[grid].gridfuncs.k12_gfs;
    MAYBE_UNUSED REAL *restrict k13_gfs = griddata[grid].gridfuncs.k13_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k12_gfs);
    rk_substep_12__launcher(params, k10_gfs, k11_gfs, k1_gfs, k4_gfs, k5_gfs, k6_gfs, k7_gfs, k8_gfs, k9_gfs, y_n_gfs, next_y_input_gfs,
                            commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, next_y_input_gfs);
  }
  // -={ END k12 substep }=-

  // -={ START k13 substep }=-
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
    MAYBE_UNUSED REAL *restrict k8_gfs = griddata[grid].gridfuncs.k8_gfs;
    MAYBE_UNUSED REAL *restrict k9_gfs = griddata[grid].gridfuncs.k9_gfs;
    MAYBE_UNUSED REAL *restrict k10_gfs = griddata[grid].gridfuncs.k10_gfs;
    MAYBE_UNUSED REAL *restrict k11_gfs = griddata[grid].gridfuncs.k11_gfs;
    MAYBE_UNUSED REAL *restrict k12_gfs = griddata[grid].gridfuncs.k12_gfs;
    MAYBE_UNUSED REAL *restrict k13_gfs = griddata[grid].gridfuncs.k13_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k13_gfs);
    rk_substep_13__launcher(params, k10_gfs, k11_gfs, k12_gfs, k13_gfs, k1_gfs, k6_gfs, k7_gfs, k8_gfs, k9_gfs, y_n_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, y_n_gfs);
  }
  // -={ END k13 substep }=-

  // Adding dt to commondata->time many times will induce roundoff error,
  // so here we set time based on the iteration number:
  commondata->time = (REAL)(commondata->nn + 1) * commondata->dt;

  // Increment the timestep n:
  commondata->nn++;
} // END FUNCTION MoL_step_forward_in_time
