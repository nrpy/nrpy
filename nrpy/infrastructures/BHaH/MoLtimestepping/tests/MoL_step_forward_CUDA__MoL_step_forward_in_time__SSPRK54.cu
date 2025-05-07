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
    static constexpr double dblRK_Rational_27567112243069_70368744177664 = 27567112243069.0 / 70368744177664.0;
    const REAL_CUDA_ARRAY RK_Rational_27567112243069_70368744177664 = ConstCUDA(dblRK_Rational_27567112243069_70368744177664);

    const REAL_CUDA_ARRAY __rk_exp_0 = FusedMulAddCUDA(RK_Rational_27567112243069_70368744177664, MulCUDA(k1_gfsL, dt), y_n_gfsL);
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

  const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_MOL_SUBSTEP;
  const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_MOL_SUBSTEP;
  const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_MOL_SUBSTEP;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_1_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, y_n_gfs, next_y_input_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_1_gpu failure");
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
    static constexpr double dblRK_Rational_1659173807685965_4503599627370496 = 1659173807685965.0 / 4503599627370496.0;
    const REAL_CUDA_ARRAY RK_Rational_1659173807685965_4503599627370496 = ConstCUDA(dblRK_Rational_1659173807685965_4503599627370496);

    static constexpr double dblRK_Rational_7842355689270681_36028797018963968 = 7842355689270681.0 / 3.602879701896397e+16;
    const REAL_CUDA_ARRAY RK_Rational_7842355689270681_36028797018963968 = ConstCUDA(dblRK_Rational_7842355689270681_36028797018963968);

    const REAL_CUDA_ARRAY __rk_exp_0 =
        FusedMulAddCUDA(RK_Rational_1659173807685965_4503599627370496, MulCUDA(k2_gfsL, dt),
                        FusedMulAddCUDA(RK_Rational_7842355689270681_36028797018963968, MulCUDA(k1_gfsL, dt), y_n_gfsL));
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

  const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_MOL_SUBSTEP;
  const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_MOL_SUBSTEP;
  const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_MOL_SUBSTEP;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_2_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k2_gfs, y_n_gfs, next_y_input_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_2_gpu failure");
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
    static constexpr double dblRK_Rational_2521268231078959_18014398509481984 = 2521268231078959.0 / 1.8014398509481984e+16;
    const REAL_CUDA_ARRAY RK_Rational_2521268231078959_18014398509481984 = ConstCUDA(dblRK_Rational_2521268231078959_18014398509481984);

    static constexpr double dblRK_Rational_4537678802552775_18014398509481984 = 4537678802552775.0 / 1.8014398509481984e+16;
    const REAL_CUDA_ARRAY RK_Rational_4537678802552775_18014398509481984 = ConstCUDA(dblRK_Rational_4537678802552775_18014398509481984);

    static constexpr double dblRK_Rational_5958592814262287_72057594037927936 = 5958592814262287.0 / 7.205759403792794e+16;
    const REAL_CUDA_ARRAY RK_Rational_5958592814262287_72057594037927936 = ConstCUDA(dblRK_Rational_5958592814262287_72057594037927936);

    const REAL_CUDA_ARRAY __rk_exp_0 =
        FusedMulAddCUDA(RK_Rational_4537678802552775_18014398509481984, MulCUDA(k3_gfsL, dt),
                        FusedMulAddCUDA(RK_Rational_5958592814262287_72057594037927936, MulCUDA(k1_gfsL, dt),
                                        FusedMulAddCUDA(RK_Rational_2521268231078959_18014398509481984, MulCUDA(k2_gfsL, dt), y_n_gfsL)));
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

  const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_MOL_SUBSTEP;
  const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_MOL_SUBSTEP;
  const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_MOL_SUBSTEP;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_3_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k2_gfs, k3_gfs, y_n_gfs, next_y_input_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_3_gpu failure");
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
    static constexpr double dblRK_Rational_4144561800390297_36028797018963968 = 4144561800390297.0 / 3.602879701896397e+16;
    const REAL_CUDA_ARRAY RK_Rational_4144561800390297_36028797018963968 = ConstCUDA(dblRK_Rational_4144561800390297_36028797018963968);

    static constexpr double dblRK_Rational_4897486879351823_72057594037927936 = 4897486879351823.0 / 7.205759403792794e+16;
    const REAL_CUDA_ARRAY RK_Rational_4897486879351823_72057594037927936 = ConstCUDA(dblRK_Rational_4897486879351823_72057594037927936);

    static constexpr double dblRK_Rational_4908696163965517_9007199254740992 = 4908696163965517.0 / 9007199254740992.0;
    const REAL_CUDA_ARRAY RK_Rational_4908696163965517_9007199254740992 = ConstCUDA(dblRK_Rational_4908696163965517_9007199254740992);

    static constexpr double dblRK_Rational_7459218339277047_36028797018963968 = 7459218339277047.0 / 3.602879701896397e+16;
    const REAL_CUDA_ARRAY RK_Rational_7459218339277047_36028797018963968 = ConstCUDA(dblRK_Rational_7459218339277047_36028797018963968);

    const REAL_CUDA_ARRAY __rk_exp_0 = FusedMulAddCUDA(
        RK_Rational_4897486879351823_72057594037927936, MulCUDA(k1_gfsL, dt),
        FusedMulAddCUDA(RK_Rational_4908696163965517_9007199254740992, MulCUDA(k4_gfsL, dt),
                        FusedMulAddCUDA(RK_Rational_7459218339277047_36028797018963968, MulCUDA(k3_gfsL, dt),
                                        FusedMulAddCUDA(RK_Rational_4144561800390297_36028797018963968, MulCUDA(k2_gfsL, dt), y_n_gfsL))));
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

  const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_MOL_SUBSTEP;
  const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_MOL_SUBSTEP;
  const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_MOL_SUBSTEP;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_4_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k2_gfs, k3_gfs, k4_gfs, y_n_gfs, next_y_input_gfs,
                                                                                  dt);
  cudaCheckErrors(cudaKernel, "rk_substep_4_gpu failure");
} // END FUNCTION rk_substep_4__launcher

/**
 * Kernel: rk_substep_5_gpu.
 * Compute RK substep 5.
 */
__global__ static void rk_substep_5_gpu(const size_t streamid, REAL *restrict k1_gfs, REAL *restrict k2_gfs, REAL *restrict k3_gfs,
                                        REAL *restrict k4_gfs, REAL *restrict k5_gfs, REAL *restrict y_n_gfs, const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k1_gfsL = k1_gfs[i];
    const REAL k2_gfsL = k2_gfs[i];
    const REAL k3_gfsL = k3_gfs[i];
    const REAL k4_gfsL = k4_gfs[i];
    const REAL k5_gfsL = k5_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    static constexpr double dblRK_Rational_1235962931917479_4503599627370496 = 1235962931917479.0 / 4503599627370496.0;
    const REAL_CUDA_ARRAY RK_Rational_1235962931917479_4503599627370496 = ConstCUDA(dblRK_Rational_1235962931917479_4503599627370496);

    static constexpr double dblRK_Rational_2644727643550321_18014398509481984 = 2644727643550321.0 / 1.8014398509481984e+16;
    const REAL_CUDA_ARRAY RK_Rational_2644727643550321_18014398509481984 = ConstCUDA(dblRK_Rational_2644727643550321_18014398509481984);

    static constexpr double dblRK_Rational_3756320236709225_36028797018963968 = 3756320236709225.0 / 3.602879701896397e+16;
    const REAL_CUDA_ARRAY RK_Rational_3756320236709225_36028797018963968 = ConstCUDA(dblRK_Rational_3756320236709225_36028797018963968);

    static constexpr double dblRK_Rational_4476270149944963_18014398509481984 = 4476270149944963.0 / 1.8014398509481984e+16;
    const REAL_CUDA_ARRAY RK_Rational_4476270149944963_18014398509481984 = ConstCUDA(dblRK_Rational_4476270149944963_18014398509481984);

    static constexpr double dblRK_Rational_8142777736761735_36028797018963968 = 8142777736761735.0 / 3.602879701896397e+16;
    const REAL_CUDA_ARRAY RK_Rational_8142777736761735_36028797018963968 = ConstCUDA(dblRK_Rational_8142777736761735_36028797018963968);

    const REAL_CUDA_ARRAY __rk_exp_0 = FusedMulAddCUDA(
        RK_Rational_2644727643550321_18014398509481984, MulCUDA(k1_gfsL, dt),
        FusedMulAddCUDA(
            RK_Rational_3756320236709225_36028797018963968, MulCUDA(k3_gfsL, dt),
            FusedMulAddCUDA(RK_Rational_4476270149944963_18014398509481984, MulCUDA(k2_gfsL, dt),
                            FusedMulAddCUDA(RK_Rational_8142777736761735_36028797018963968, MulCUDA(k5_gfsL, dt),
                                            FusedMulAddCUDA(RK_Rational_1235962931917479_4503599627370496, MulCUDA(k4_gfsL, dt), y_n_gfsL)))));
    WriteCUDA(&y_n_gfs[i], __rk_exp_0);
  }
} // END FUNCTION rk_substep_5_gpu

/**
 * Runge-Kutta function for substep 5.
 */
static void rk_substep_5__launcher(params_struct *restrict params, REAL *restrict k1_gfs, REAL *restrict k2_gfs, REAL *restrict k3_gfs,
                                   REAL *restrict k4_gfs, REAL *restrict k5_gfs, REAL *restrict y_n_gfs, const REAL dt) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;

  const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_MOL_SUBSTEP;
  const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_MOL_SUBSTEP;
  const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_MOL_SUBSTEP;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((Ntot + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  rk_substep_5_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, k1_gfs, k2_gfs, k3_gfs, k4_gfs, k5_gfs, y_n_gfs, dt);
  cudaCheckErrors(cudaKernel, "rk_substep_5_gpu failure");
} // END FUNCTION rk_substep_5__launcher

/**
 * Method of Lines (MoL) for "SSPRK54" method: Step forward one full timestep.
 *
 */
void MoL_step_forward_in_time(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {

  // C code implementation of -={ SSPRK54 }=- Method of Lines timestepping.

  // First set the initial time:
  const REAL time_start = commondata->time;
  // -={ START k1 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 0.00000000000000000e+00 * commondata->dt;
    cpyHosttoDevice_params__constant(&griddata[grid].params, griddata[grid].params.grid_idx % NUM_STREAMS);
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict next_y_input_gfs = griddata[grid].gridfuncs.next_y_input_gfs;
    MAYBE_UNUSED REAL *restrict k1_gfs = griddata[grid].gridfuncs.k1_gfs;
    MAYBE_UNUSED REAL *restrict k2_gfs = griddata[grid].gridfuncs.k2_gfs;
    MAYBE_UNUSED REAL *restrict k3_gfs = griddata[grid].gridfuncs.k3_gfs;
    MAYBE_UNUSED REAL *restrict k4_gfs = griddata[grid].gridfuncs.k4_gfs;
    MAYBE_UNUSED REAL *restrict k5_gfs = griddata[grid].gridfuncs.k5_gfs;
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
    commondata->time = time_start + 3.91752227003919984e-01 * commondata->dt;
    cpyHosttoDevice_params__constant(&griddata[grid].params, griddata[grid].params.grid_idx % NUM_STREAMS);
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict next_y_input_gfs = griddata[grid].gridfuncs.next_y_input_gfs;
    MAYBE_UNUSED REAL *restrict k1_gfs = griddata[grid].gridfuncs.k1_gfs;
    MAYBE_UNUSED REAL *restrict k2_gfs = griddata[grid].gridfuncs.k2_gfs;
    MAYBE_UNUSED REAL *restrict k3_gfs = griddata[grid].gridfuncs.k3_gfs;
    MAYBE_UNUSED REAL *restrict k4_gfs = griddata[grid].gridfuncs.k4_gfs;
    MAYBE_UNUSED REAL *restrict k5_gfs = griddata[grid].gridfuncs.k5_gfs;
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
    commondata->time = time_start + 5.86079688967789947e-01 * commondata->dt;
    cpyHosttoDevice_params__constant(&griddata[grid].params, griddata[grid].params.grid_idx % NUM_STREAMS);
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict next_y_input_gfs = griddata[grid].gridfuncs.next_y_input_gfs;
    MAYBE_UNUSED REAL *restrict k1_gfs = griddata[grid].gridfuncs.k1_gfs;
    MAYBE_UNUSED REAL *restrict k2_gfs = griddata[grid].gridfuncs.k2_gfs;
    MAYBE_UNUSED REAL *restrict k3_gfs = griddata[grid].gridfuncs.k3_gfs;
    MAYBE_UNUSED REAL *restrict k4_gfs = griddata[grid].gridfuncs.k4_gfs;
    MAYBE_UNUSED REAL *restrict k5_gfs = griddata[grid].gridfuncs.k5_gfs;
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
    commondata->time = time_start + 4.74542363026869973e-01 * commondata->dt;
    cpyHosttoDevice_params__constant(&griddata[grid].params, griddata[grid].params.grid_idx % NUM_STREAMS);
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict next_y_input_gfs = griddata[grid].gridfuncs.next_y_input_gfs;
    MAYBE_UNUSED REAL *restrict k1_gfs = griddata[grid].gridfuncs.k1_gfs;
    MAYBE_UNUSED REAL *restrict k2_gfs = griddata[grid].gridfuncs.k2_gfs;
    MAYBE_UNUSED REAL *restrict k3_gfs = griddata[grid].gridfuncs.k3_gfs;
    MAYBE_UNUSED REAL *restrict k4_gfs = griddata[grid].gridfuncs.k4_gfs;
    MAYBE_UNUSED REAL *restrict k5_gfs = griddata[grid].gridfuncs.k5_gfs;
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
    commondata->time = time_start + 9.35010631009239979e-01 * commondata->dt;
    cpyHosttoDevice_params__constant(&griddata[grid].params, griddata[grid].params.grid_idx % NUM_STREAMS);
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict next_y_input_gfs = griddata[grid].gridfuncs.next_y_input_gfs;
    MAYBE_UNUSED REAL *restrict k1_gfs = griddata[grid].gridfuncs.k1_gfs;
    MAYBE_UNUSED REAL *restrict k2_gfs = griddata[grid].gridfuncs.k2_gfs;
    MAYBE_UNUSED REAL *restrict k3_gfs = griddata[grid].gridfuncs.k3_gfs;
    MAYBE_UNUSED REAL *restrict k4_gfs = griddata[grid].gridfuncs.k4_gfs;
    MAYBE_UNUSED REAL *restrict k5_gfs = griddata[grid].gridfuncs.k5_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(commondata, params, rfmstruct, auxevol_gfs, next_y_input_gfs, k5_gfs);
    rk_substep_5__launcher(params, k1_gfs, k2_gfs, k3_gfs, k4_gfs, k5_gfs, y_n_gfs, commondata->dt);
    if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
      apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, y_n_gfs);
  }
  // -={ END k5 substep }=-

  // Adding dt to commondata->time many times will induce roundoff error,
  // so here we set time based on the iteration number:
  commondata->time = (REAL)(commondata->nn + 1) * commondata->dt;

  // Increment the timestep n:
  commondata->nn++;
} // END FUNCTION MoL_step_forward_in_time
