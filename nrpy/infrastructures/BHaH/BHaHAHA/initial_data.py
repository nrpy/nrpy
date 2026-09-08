"""
Register C functions for setting up initial data for the 4-metric g_{mu nu}.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_function as cfc


def register_CFunction_initial_data() -> None:
    """
    Register the C function responsible for initializing data related to the 4-metric g_{mu nu}.

    This function sets up the necessary parameters, includes, and body of the C function that reads
    3D metric data, performs basis transformations, applies boundary conditions, computes derivatives,
    and initializes the field h(theta, phi) on all grids.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = r"""
/**
 * Reads 3D metric data (in Cartesian basis) from file, basis transform, apply BCs, compute h_{ij,k}, then set up initial guess h(r,theta)
 */
#ifdef __CUDACC__
__global__
#endif
void set_initial_data(commondata_struct *restrict commondata, griddata_struct *restrict griddata, REAL *coarse_to_fine, REAL dst_pts[][2])
{
#ifdef __CUDACC__
  // Set up cooperative group
  namespace cg = cooperative_groups;
  cg::grid_group gpu_grid = cg::this_grid();
#endif

  int grid = 0;
  params_struct *params = &griddata[0].params;
#include "set_CodeParameters.h"
  const int NUM_THETA = Nxx1;
  if (commondata->use_coarse_horizon) {
    const int num_dst_pts = Nxx1 * Nxx2;
    PARALLEL_2D_LOOP(i1, NGHOSTS, Nxx1 + NGHOSTS, i2, NGHOSTS, Nxx2 + NGHOSTS) {
          dst_pts[IDX2(i1 - NGHOSTS, i2 - NGHOSTS)][0] = griddata->xx[1][i1];
          dst_pts[IDX2(i1 - NGHOSTS, i2 - NGHOSTS)][1] = griddata->xx[2][i2];
    } END_PARALLEL_2D_LOOP
      // Choosing an NGHOSTS stencil half-width significantly speeds up BHaHAHA finds.

    bah_interpolation_2d_general__uniform_src_grid(NGHOSTS, commondata->coarse_horizon_dxx1, commondata->coarse_horizon_dxx2,
                                                   commondata->coarse_horizon_Nxx_plus_2NGHOSTS1, commondata->coarse_horizon_Nxx_plus_2NGHOSTS2,
                                                   commondata->coarse_horizon_r_theta_phi, commondata->coarse_horizon, num_dst_pts, dst_pts,
                                                   coarse_to_fine, &commondata->error_flag);
#ifdef __CUDACC__
    gpu_grid.sync();
#endif
    if (commondata->error_flag != BHAHAHA_SUCCESS)
      return;
  }

  // Step 2.a: Use CUDA or OpenMP to parallelize the loop over the entire grid, initializing the h(theta, phi) scalar.
  const REAL times[3] = {commondata->bhahaha_params_and_data->t_m1, //
                         commondata->bhahaha_params_and_data->t_m2, //
                         commondata->bhahaha_params_and_data->t_m3};
  PARALLEL_LOOP(
           i0, NGHOSTS, Nxx0 + NGHOSTS,   // Loop over radial grid points
           i1, NGHOSTS, Nxx1 + NGHOSTS,   // Loop over polar grid points
           i2, NGHOSTS, Nxx2 + NGHOSTS) { // Loop over azimuthal grid points
    // Step 2.b: Set an initial guess value for the h(theta, phi) field at each grid point.
    //      NOTE: coarse_to_fine (interpolated above) & horizon_guess contain Nxx1 x Nxx2 points.
    if (commondata->use_coarse_horizon) {
      griddata[grid].gridfuncs.y_n_gfs[IDX4(HHGF, i0, i1, i2)] = coarse_to_fine[IDX2(i1 - NGHOSTS, i2 - NGHOSTS)];
    } else if (commondata->bhahaha_params_and_data->use_fixed_radius_guess_on_full_sphere) {
      // r_max_interior = r_min_external_input + ((Nr_external_input-BHAHAHA_NGHOSTS) + 0.5) * dr
      const REAL r_max_interior =
          commondata->bhahaha_params_and_data->r_min_external_input +
          ((commondata->bhahaha_params_and_data->Nr_external_input - BHAHAHA_NGHOSTS) + 0.5) * commondata->bhahaha_params_and_data->dr_external_input;
      griddata[grid].gridfuncs.y_n_gfs[IDX4(HHGF, i0, i1, i2)] = 0.8 * r_max_interior;
    } else {
      griddata[grid].gridfuncs.y_n_gfs[IDX4(HHGF, i0, i1, i2)] =
          bah_quadratic_extrapolation(times, //
                                      commondata->bhahaha_params_and_data->prev_horizon_m1[IDX2(i1 - NGHOSTS, i2 - NGHOSTS)],
                                      commondata->bhahaha_params_and_data->prev_horizon_m2[IDX2(i1 - NGHOSTS, i2 - NGHOSTS)],
                                      commondata->bhahaha_params_and_data->prev_horizon_m3[IDX2(i1 - NGHOSTS, i2 - NGHOSTS)],
                                      commondata->bhahaha_params_and_data->time_external_input);
    }
    // set VVGF = eta * HHGF,
    //  so that partial_t h = VVGF - eta * HHGF = 0 at t=0. Otherwise we get really ugly dynamics.
    griddata[grid].gridfuncs.y_n_gfs[IDX4(VVGF, i0, i1, i2)] = eta_damping * griddata[grid].gridfuncs.y_n_gfs[IDX4(HHGF, i0, i1, i2)];
  } END_PARALLEL_LOOP // END LOOP over all gridpoints

#ifdef __CUDACC__
  gpu_grid.sync();
#endif

  bah_apply_bcs_inner_only(commondata, &griddata[grid].params, &griddata[grid].bcstruct, griddata[grid].gridfuncs.y_n_gfs);

  return;
}
"""
    desc = "Read 3D metric data (in Cartesian basis) from file, basis transform, apply BCs, compute h_{ij,k}, then set up initial guess h(r,theta)"
    cfunc_type = "void"
    name = "initial_data"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    body = r"""
  const int grid = 0;
  const params_struct *restrict params = &griddata[grid].params;

  // Allocate memory for set_initial_data kernel
#include "set_CodeParameters.h"
  REAL *restrict coarse_to_fine = NULL;
  REAL(*dst_pts)[2] = NULL;
  if (commondata->use_coarse_horizon) {
    const int num_dst_pts = Nxx1 * Nxx2;
    #ifdef __CUDACC__
    cudaMalloc((void**)&coarse_to_fine, sizeof(REAL) * Nxx1 * Nxx2);
    cudaMalloc((void**)&dst_pts, sizeof(REAL)*num_dst_pts*2);
    #else
    coarse_to_fine = (double *)malloc(sizeof(REAL) * Nxx1 * Nxx2);
    dst_pts = (double (*)[2])malloc(num_dst_pts * sizeof(*dst_pts));
    #endif
    if (dst_pts == NULL || coarse_to_fine == NULL) {
      commondata->error_flag = INITIAL_DATA_MALLOC_ERROR;
      return;
    }
  }
#ifdef __CUDACC__
  // Create device side copies of commondata and griddata
  commondata_struct *d_commondata = NULL;
  cudaMalloc((void**)&d_commondata, sizeof(commondata_struct));
  cudaMemcpy(d_commondata, commondata, sizeof(commondata_struct), cudaMemcpyHostToDevice);

  griddata_struct *d_griddata = NULL;
  cudaMalloc((void**)&d_griddata, sizeof(griddata_struct));
  cudaMemcpy(d_griddata, griddata, sizeof(griddata_struct), cudaMemcpyHostToDevice);
#endif
    
#ifdef __CUDACC__
  void *Args[] = {&d_commondata, &d_griddata, (void *)&coarse_to_fine, &dst_pts};
  int INTERP_ORDER = (2 * NGHOSTS + 1);
  int sharesize = sizeof(REAL)*(INTERP_ORDER + (THREADSPERBLOCK*4)*INTERP_ORDER + THREADSPERBLOCK*INTERP_ORDER*INTERP_ORDER);
  COOPERATIVE_KERNEL(set_initial_data, Args, sharesize)
#else
  set_initial_data(commondata, griddata, coarse_to_fine, dst_pts); 
#endif

#ifdef __CUDACC__
  //Retrieve commondata
  gpuErrchk( cudaMemcpy(commondata, d_commondata, sizeof(commondata_struct), cudaMemcpyDeviceToHost) );
  //Free device side copies of commondata and griddata
  cudaFree(d_commondata);
  cudaFree(d_griddata);
#endif

  // Free memory
  if (commondata->use_coarse_horizon) {
    FREE(coarse_to_fine);
    FREE(dst_pts);
    FREE(commondata->coarse_horizon);
    for (int ii = 0; ii < 3; ii++)
      FREE(commondata->coarse_horizon_r_theta_phi[ii]);
  }

  commondata->use_coarse_horizon = 1; // for next time initial_data() is called

  return;
"""
    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,  # params not passed to function
        body=body,
    )
