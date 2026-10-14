"""
Register the C function for computing minimum, maximum, and mean radius *from the coordinate centroid* of the horizon.
Needs: coordinate centroid of the horizon.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.general_relativity.bhahaha.area as bhahaha_area
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_diagnostics_min_max_mean_radii_wrt_centroid(
    enable_fd_functions: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for computing minimum, maximum, and mean radius *from the coordinate centroid* of the horizon.
    Needs: coordinate centroid of the horizon.

    :param enable_fd_functions: Whether to enable finite difference functions, defaults to True.
    :return: An NRPyEnv_type object if registration is successful, otherwise None.

    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = r"""
#ifdef __CUDACC__
/**
 * Function that leverages CUDA's atomicCAS to compute the minimum between 
 * a double stored at the initial address and a new double. Akin to CUDA's 
 * atomicMin function but for doubles, this function returns the previous 
 * value stored at the original address.
 */
__device__ double atomicMin_double(double* address, double val)
{
  unsigned long long int* address_as_ull = (unsigned long long int*) address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed, 
        __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
  } while (assumed != old);
  return __longlong_as_double(old);
}

/**
 * Function that leverages CUDA's atomicCAS to compute the maximium between 
 * a double stored at the initial address and a new double. Akin to CUDA's 
 * atomicMax function but for doubles, this function returns the previous 
 * value stored at the original address.
 */
__device__ double atomicMax_double(double* address, double val)
{
  unsigned long long int* address_as_ull = (unsigned long long int*) address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed, 
        __double_as_longlong(fmax(val, __longlong_as_double(assumed))));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif
"""
    desc = "BHaHAHA apparent horizon diagnostics: Compute Theta L2 and Linfinity norms."
    cfunc_type = "void"
    name = "diagnostics_min_max_mean_radii_wrt_centroid"
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
  // Set up cooperative group
  namespace cg = cooperative_groups;
  cg::grid_group gpu_grid = cg::this_grid();
  
  // Set up shared memory
  extern __shared__ REAL s[];
  
  // Set global space for intermediate norm calculation and storage
  diag_norms_struct *norms = commondata->norms;
  #endif
  const int grid = 0;
  bhahaha_diagnostics_struct *restrict bhahaha_diags = commondata->bhahaha_diagnostics;
  const params_struct *restrict params = &griddata[grid].params;
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
  const REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];
#include "set_CodeParameters.h"

  // Set integration weights.
  const REAL *restrict weights;
  int weight_stencil_size;
  bah_diagnostics_integration_weights(Nxx1, Nxx2, &weights, &weight_stencil_size);

  // Compute radii quantities. Mean radius is area-weighted, since physical gridspacing is quite uneven.
  #ifdef __CUDACC__
  CUDA_ONE_THREAD(gpu_grid) {
    norms->sum_curr_area = 0.0;
    norms->sum_mean_radius = 0.0;
    norms->min_radius_squared = +1e30;
    norms->max_radius_squared = -1e30;
  } END_CUDA_ONE_THREAD;
  REAL *s_sum_curr_area = &s[0];
  REAL *s_sum_mean_radius = &s[1*blockDim.x];
  REAL *s_min_radius_squared = &s[2*blockDim.x];
  REAL *s_max_radius_squared = &s[3*blockDim.x];

  int tid = threadIdx.x;
  s_sum_curr_area[threadIdx.x] = 0.0;
  s_sum_mean_radius[threadIdx.x] = 0.0;
  s_min_radius_squared[threadIdx.x] = +1e30;
  s_max_radius_squared[threadIdx.x] = -1e30;
  #endif
  REAL sum_curr_area = 0.0;
  REAL sum_mean_radius = 0.0;
  REAL min_radius_squared = +1e30;
  REAL max_radius_squared = -1e30;


#ifdef __CUDACC__
  CUDA_3D_LOOP(i0, NGHOSTS, NGHOSTS + Nxx0, i1, NGHOSTS, NGHOSTS + Nxx1, i2, NGHOSTS, NGHOSTS + Nxx2) {
      const REAL weight2 = weights[(i2 - NGHOSTS) % weight_stencil_size];
      MAYBE_UNUSED const REAL xx2 = xx[2][i2];
      const REAL weight1 = weights[(i1 - NGHOSTS) % weight_stencil_size];
      MAYBE_UNUSED const REAL xx1 = xx[1][i1];
#else
#pragma omp parallel
  {
#pragma omp for
    for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
      const REAL weight2 = weights[(i2 - NGHOSTS) % weight_stencil_size];
      MAYBE_UNUSED const REAL xx2 = xx[2][i2];
      for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
        const REAL weight1 = weights[(i1 - NGHOSTS) % weight_stencil_size];
        MAYBE_UNUSED const REAL xx1 = xx[1][i1];
        for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
#endif          
"""
    body += (
        ccg.c_codegen(
            bhahaha_area.area3(),
            "const REAL area_element",
            enable_fd_codegen=True,
            enable_fd_functions=enable_fd_functions,
        )
        + """
#ifndef __CUDACC__
#pragma omp critical
          {
            sum_curr_area += area_element * weight1 * weight2;
            const REAL r = in_gfs[IDX4(HHGF, NGHOSTS, i1, i2)];
            const REAL theta = commondata->interp_src_r_theta_phi[1][i1];
            const REAL phi = commondata->interp_src_r_theta_phi[2][i2];
            const REAL xx = r * sin(theta) * cos(phi);
            const REAL yy = r * sin(theta) * sin(phi);
            const REAL zz = r * cos(theta);
            // Radius as measured from AH centroid:
            REAL radius_squared = ((xx - bhahaha_diags->x_centroid_wrt_coord_origin) * (xx - bhahaha_diags->x_centroid_wrt_coord_origin) + //
                                   (yy - bhahaha_diags->y_centroid_wrt_coord_origin) * (yy - bhahaha_diags->y_centroid_wrt_coord_origin) + //
                                   (zz - bhahaha_diags->z_centroid_wrt_coord_origin) * (zz - bhahaha_diags->z_centroid_wrt_coord_origin));
            sum_mean_radius += sqrt(radius_squared) * area_element * weight1 * weight2;
            if (radius_squared < min_radius_squared)
              min_radius_squared = radius_squared;
            if (radius_squared > max_radius_squared)
              max_radius_squared = radius_squared;
          } // END OMP CRITICAL
        } // END LOOP over i0
      } // END LOOP over i1
    } // END LOOP over i2
  } // END OMP PARALLEL
#else
//Prep for Reduction
          s_sum_curr_area[tid] += area_element *weight1 * weight2;
          const REAL r = in_gfs[IDX4(HHGF, NGHOSTS, i1, i2)];
          const REAL theta = commondata->interp_src_r_theta_phi[1][i1];
          const REAL phi = commondata->interp_src_r_theta_phi[2][i2];
          const REAL xx = r * sin(theta) * cos(phi);
          const REAL yy = r * sin(theta) * sin(phi);
          const REAL zz = r * cos(theta);
          REAL radius_squared = ((xx - bhahaha_diags->x_centroid_wrt_coord_origin) * (xx - bhahaha_diags->x_centroid_wrt_coord_origin) + //
                                 (yy - bhahaha_diags->y_centroid_wrt_coord_origin) * (yy - bhahaha_diags->y_centroid_wrt_coord_origin) + //
                                 (zz - bhahaha_diags->z_centroid_wrt_coord_origin) * (zz - bhahaha_diags->z_centroid_wrt_coord_origin));
          s_sum_mean_radius[tid] += sqrt(radius_squared) * area_element * weight1 * weight2;
          if (radius_squared < s_min_radius_squared[tid])
            s_min_radius_squared[tid] = radius_squared;
          if (radius_squared > s_max_radius_squared[tid])
            s_max_radius_squared[tid]  = radius_squared;
} END_CUDA_3D_LOOP          
#endif

#ifdef __CUDACC__
  gpu_grid.sync();
  // Block level Reduction to thread 0 in shared memory
  unsigned int halfsize = blockDim.x>>1;
  while (halfsize > 0) {
    if (tid < halfsize) {
      if ((blockDim.x*blockIdx.x + tid) < Nxx0*Nxx1*Nxx2 && (blockDim.x*blockIdx.x + (tid + halfsize)) < Nxx0*Nxx1*Nxx2 ) {
        s_sum_curr_area[tid] += s_sum_curr_area[tid + halfsize];
        s_sum_mean_radius[tid] += s_sum_mean_radius[tid + halfsize];
        if (s_min_radius_squared[tid] > s_min_radius_squared[tid + halfsize])
          s_min_radius_squared[tid] = s_min_radius_squared[tid + halfsize];
        if (s_max_radius_squared[tid] < s_max_radius_squared[tid + halfsize])
          s_max_radius_squared[tid] = s_max_radius_squared[tid + halfsize];
      }
    }
    __syncthreads();
    halfsize = halfsize>>1;
  }
  // Grid level reduction to norms struct
  if (blockDim.x*blockIdx.x < Nxx0*Nxx1*Nxx2 && threadIdx.x == 0) {
    atomicAdd(&norms->sum_curr_area, s_sum_curr_area[0]);
    atomicAdd(&norms->sum_mean_radius, s_sum_mean_radius[0]);
    atomicMin_double(&norms->min_radius_squared, s_min_radius_squared[0]);
    atomicMax_double(&norms->max_radius_squared, s_max_radius_squared[0]);
  }
  gpu_grid.sync();
#endif

#ifdef __CUDACC__
  CUDA_ONE_THREAD(gpu_grid)
#endif
  {
#ifdef __CUDACC__
    sum_curr_area = norms->sum_curr_area;
    sum_mean_radius = norms->sum_mean_radius;
    min_radius_squared = norms->min_radius_squared;
    max_radius_squared = norms->max_radius_squared;
#endif

    // Store commondata diagnostic parameters.
    const REAL curr_area = sum_curr_area * params->dxx1 * params->dxx2;
    bhahaha_diags->min_coord_radius_wrt_centroid = sqrt(min_radius_squared);
    bhahaha_diags->max_coord_radius_wrt_centroid = sqrt(max_radius_squared);
    bhahaha_diags->mean_coord_radius_wrt_centroid = sum_mean_radius * params->dxx1 * params->dxx2 / curr_area;
  }
#ifdef __CUDACC__
  END_CUDA_ONE_THREAD;
#endif
"""
    )
    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        cfunc_decorators=cfunc_decorators,
        body=body,
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
