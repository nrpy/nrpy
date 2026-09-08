"""
Module for registering the CFL-limited timestep C function.

Registers a C function that computes the timestep based on the minimum grid spacing
in a 2D spherical numerical grid.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_function as cfc


def register_CFunction_cfl_limited_timestep_based_on_h_equals_r() -> None:
    """
    Register the C function for computing the CFL-limited timestep.

    The registered C function calculates the timestep using:
        dt = CFL_FACTOR * ds_min
    where ds_min is the smallest grid spacing in the numerical grid.

    This ensures stability by adhering to the CFL condition.
    """
    includes = ["BHaH_defines.h"]
    description = "Compute minimum timestep dt = CFL_FACTOR * ds_min on a 2D spherical numerical grid."
    cfunc_type = "void"
    name = "cfl_limited_timestep_based_on_h_equals_r"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    cfunc_decorators = r"""
#ifdef __CUDACC__
__device__
#endif
"""
    body = r"""
#ifdef __CUDACC__
  //Set up cooperative group
  namespace cg = cooperative_groups;
  cg::grid_group gpu_grid = cg::this_grid();

  //Setup shared memory
  extern __shared__ REAL s[];
#endif
#ifdef __CUDACC__
  CUDA_ONE_THREAD(gpu_grid) {
#endif
  commondata->dt = 1e30;
#ifdef __CUDACC__
  } END_CUDA_ONE_THREAD //End global variable modification
#endif
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    const params_struct *restrict params = &griddata[grid].params;
    const REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
    REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++) {
      xx[ww] = griddata[grid].xx[ww];
    }

#include "set_CodeParameters.h"

#ifdef __CUDACC__
    REAL *ds_min = s;
    int tid = threadIdx.x;
    ds_min[tid] = 1e38;
#else
    REAL ds_min = 1e38;
#endif
#ifdef __CUDACC__
    PARALLEL_LOOP(i2, NGHOSTS, Nxx2 + NGHOSTS, i1, NGHOSTS, Nxx1 + NGHOSTS, i0, NGHOSTS, Nxx0 + NGHOSTS) {
#else
#pragma omp parallel for reduction(min : ds_min)
    LOOP_NOOMP(i0, NGHOSTS, Nxx0 + NGHOSTS, i1, NGHOSTS, Nxx1 + NGHOSTS, i2, NGHOSTS, Nxx2 + NGHOSTS) {
#endif
      const REAL hh = in_gfs[IDX4(HHGF, i0, i1, i2)];
      const REAL xx1 = xx[1][i1];
      REAL dsmin1, dsmin2;

      dsmin1 = fabs(hh * dxx1);
      dsmin2 = fabs(hh * dxx2 * sin(xx1));
#ifdef __CUDACC__
      ds_min[tid] = NRPYMIN(ds_min[tid], NRPYMIN(dsmin1, dsmin2));
#else
      ds_min = NRPYMIN(ds_min, NRPYMIN(dsmin1, dsmin2));
#endif
#ifndef __CUDACC__
    }
#else
    } END_PARALLEL_LOOP
    gpu_grid.sync();
#endif
#ifdef __CUDACC__
    //Reduction
    unsigned int halfsize = blockDim.x>>1;
    while (halfsize > 0) {
      if (tid < halfsize) {
        if ((blockDim.x*blockIdx.x + tid) < Nxx0*Nxx1*Nxx2 && (blockDim.x*blockIdx.x + (tid + halfsize)) < Nxx0*Nxx1*Nxx2 ) {
          if (ds_min[tid] > ds_min[tid + halfsize])
            ds_min[tid] = ds_min[tid + halfsize];
        }
      }
      __syncthreads();
      halfsize = halfsize>>1;
    }
    if (blockDim.x*blockIdx.x < Nxx0*Nxx1*Nxx2 && threadIdx.x == 0 ) {
      atomicMin_double(&commondata->dt, ds_min[0] * commondata->CFL_FACTOR);
    }
    gpu_grid.sync();
#else 
    commondata->dt = NRPYMIN(commondata->dt, ds_min * commondata->CFL_FACTOR);
#endif
  }
"""

    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        desc=description,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
