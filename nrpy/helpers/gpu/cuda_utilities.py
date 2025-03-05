"""
Module that provides CUDA specific utilities.

Authors: Samuel D. Tootle; sdtootle **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.params as par  # NRPy+: Parameter interface
from nrpy.helpers.gpu.gpu_kernel import GPU_Kernel

"""
    Define the default launch dictionary for CUDA kernels.

    threads_per_block: Number of threads per block. Default is 32 in the x direction, 1 in the y and z directions.
    blocks_per_grid: Number of blocks needed in the [x,y,z] direction.  Default is computed based on Nxx and threads_per_block.
    stream: Stream ID for the kernel.  Empty string defaults to param_streamid % NUM_STREAMS. Exclude from dictionary to use null stream.
"""
default_launch_dictionary = {
    "blocks_per_grid": [],
    "threads_per_block": ["32"],
    "stream": "",
}


# Define functions to copy params to device
def register_CFunction_cpyHosttoDevice_params__constant() -> None:
    """
    Register C function for copying params to __constant__ space on device.

    DOCTEST:
    >>> import nrpy.c_function as cfc
    >>> register_CFunction_cpyHosttoDevice_params__constant()
    >>> print(cfc.CFunction_dict['cpyHosttoDevice_params__constant'].full_function)
    #include "../BHaH_defines.h"
    /**
     * Copy parameters to GPU __constant__.
     */
    __host__ void cpyHosttoDevice_params__constant(const params_struct *restrict params, const int streamid) {
      cudaMemcpyToSymbol(d_params, params, sizeof(params_struct), streamid * sizeof(params_struct), cudaMemcpyHostToDevice);
    } // END FUNCTION cpyHosttoDevice_params__constant
    <BLANKLINE>
    """
    includes = ["BHaH_defines.h"]

    desc = r"""Copy parameters to GPU __constant__."""
    cfunc_type = "__host__ void"
    name = "cpyHosttoDevice_params__constant"
    params = r"""const params_struct *restrict params, const int streamid"""
    body = "cudaMemcpyToSymbol(d_params, params, sizeof(params_struct), streamid * sizeof(params_struct), cudaMemcpyHostToDevice);"
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        subdirectory="CUDA_utils",
    )


# Define functions to copy params to device
def register_CFunction_cpyHosttoDevice_bc_struct() -> None:
    """
    Register C function for copying relevant aspects of bc_struct to the device.

    Currently bc_struct is allocated on the host with inner/outer_bc arrays being
    allocated on the device.

    DOCTEST:
    >>> import nrpy.c_function as cfc
    >>> register_CFunction_cpyHosttoDevice_bc_struct()
    >>> print(cfc.CFunction_dict['cpyHosttoDevice_bc_struct'].full_function)
    #include "../BHaH_defines.h"
    /**
     * Copy parameters to GPU __constant__.
     */
    __host__ void cpyHosttoDevice_bc_struct(const bc_struct *restrict host_src, bc_struct *restrict device_dst, const int streamid) {
    <BLANKLINE>
      // Copy the bc_info structure (basic metadata)
      // This is always stored on the host currently.
      memcpy(&device_dst->bc_info, &host_src->bc_info, sizeof(bc_info_struct));
    <BLANKLINE>
      // Copy inner boundary array if it exists
      int num_inner = host_src->bc_info.num_inner_boundary_points;
      if (num_inner > 0) {
        cudaMalloc(&device_dst->inner_bc_array, sizeof(innerpt_bc_struct) * num_inner);
        cudaCheckErrors(cudaMalloc, "inner_bc_array malloc failed");
    <BLANKLINE>
        cudaMemcpy(device_dst->inner_bc_array, host_src->inner_bc_array, sizeof(innerpt_bc_struct) * num_inner, cudaMemcpyHostToDevice);
        cudaCheckErrors(cudaMemcpy, "inner_bc_array Memcpy failed.");
      } else {
        device_dst->inner_bc_array = NULL;
      }
    <BLANKLINE>
      // Copy outer boundary arrays
      for (int gz = 0; gz < NGHOSTS; gz++) {
        for (int dirn = 0; dirn < 3; dirn++) {
          int num_outer = host_src->bc_info.num_pure_outer_boundary_points[gz][dirn];
          if (num_outer > 0) {
            int idx = gz * 3 + dirn; // Calculate the index in the array
    <BLANKLINE>
            cudaMalloc(&device_dst->pure_outer_bc_array[idx], sizeof(outerpt_bc_struct) * num_outer);
            cudaCheckErrors(cudaMalloc, "inner_bc_array malloc failed");
    <BLANKLINE>
            cudaMemcpy(device_dst->pure_outer_bc_array[idx], host_src->pure_outer_bc_array[idx], sizeof(outerpt_bc_struct) * num_outer,
                       cudaMemcpyHostToDevice);
            cudaCheckErrors(cudaMemcpy, "inner_bc_array Memcpy failed.");
          } else {
            device_dst->pure_outer_bc_array[gz * 3 + dirn] = NULL;
          }
        }
      }
    } // END FUNCTION cpyHosttoDevice_bc_struct
    <BLANKLINE>
    """
    includes = ["BHaH_defines.h"]

    desc = r"""Copy parameters to GPU __constant__."""
    cfunc_type = "__host__ void"
    name = "cpyHosttoDevice_bc_struct"
    params = r"""const bc_struct *restrict host_src, bc_struct *restrict device_dst, const int streamid"""
    body = """
  // Copy the bc_info structure (basic metadata)
  // This is always stored on the host currently.
  memcpy(&device_dst->bc_info, &host_src->bc_info, sizeof(bc_info_struct));

  // Copy inner boundary array if it exists
  int num_inner = host_src->bc_info.num_inner_boundary_points;
  if (num_inner > 0) {
      cudaMalloc(&device_dst->inner_bc_array, sizeof(innerpt_bc_struct) * num_inner);
      cudaCheckErrors(cudaMalloc, "inner_bc_array malloc failed");

      cudaMemcpy(device_dst->inner_bc_array, host_src->inner_bc_array, sizeof(innerpt_bc_struct) * num_inner, cudaMemcpyHostToDevice);
      cudaCheckErrors(cudaMemcpy, "inner_bc_array Memcpy failed.");
  } else {
      device_dst->inner_bc_array = NULL;
  }

  // Copy outer boundary arrays
  for (int gz = 0; gz < NGHOSTS; gz++) {
      for (int dirn = 0; dirn < 3; dirn++) {
          int num_outer = host_src->bc_info.num_pure_outer_boundary_points[gz][dirn];
          if (num_outer > 0) {
              int idx = gz * 3 + dirn; // Calculate the index in the array

              cudaMalloc(&device_dst->pure_outer_bc_array[idx], sizeof(outerpt_bc_struct) * num_outer);
              cudaCheckErrors(cudaMalloc, "inner_bc_array malloc failed");

              cudaMemcpy(device_dst->pure_outer_bc_array[idx], host_src->pure_outer_bc_array[idx],sizeof(outerpt_bc_struct) * num_outer, cudaMemcpyHostToDevice);
              cudaCheckErrors(cudaMemcpy, "inner_bc_array Memcpy failed.");
          } else {
              device_dst->pure_outer_bc_array[gz * 3 + dirn] = NULL;
          }
      }
  }
    """
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        subdirectory="CUDA_utils",
    )


# Define functions to copy params to device
def register_CFunction_cpyHosttoDevice_commondata__constant() -> None:
    """
    Register C function for copying commondata to __constant__ space on device.

    DOCTEST:
    >>> import nrpy.c_function as cfc
    >>> register_CFunction_cpyHosttoDevice_commondata__constant()
    >>> print(cfc.CFunction_dict['cpyHosttoDevice_commondata__constant'].full_function)
    #include "../BHaH_defines.h"
    /**
     * Copy parameters to GPU __constant__.
     */
    __host__ void cpyHosttoDevice_commondata__constant(const commondata_struct *restrict commondata) {
      cudaMemcpyToSymbol(d_commondata, commondata, sizeof(commondata_struct));
    } // END FUNCTION cpyHosttoDevice_commondata__constant
    <BLANKLINE>
    """
    includes = ["BHaH_defines.h"]

    desc = r"""Copy parameters to GPU __constant__."""
    cfunc_type = "__host__ void"
    name = "cpyHosttoDevice_commondata__constant"
    params = r"""const commondata_struct *restrict commondata"""
    body = "cudaMemcpyToSymbol(d_commondata, commondata, sizeof(commondata_struct));"
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        subdirectory="CUDA_utils",
    )


def generate_CFunction_mallocHostgrid() -> str:
    """
    Generate the kernel that allocates Host side storage.
    :returns: Full kernel function string

    DOCTEST:
    >>> kernel = generate_CFunction_mallocHostgrid()
    >>> print(kernel)
    /**
     * Kernel: mallocHostgrid.
     * Allocate griddata_struct[grid].xx for host.
     */
    __host__ static void mallocHostgrid(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                        griddata_struct *restrict gd_host, const griddata_struct *restrict gd_gpu) {
    <BLANKLINE>
      int const &Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
      int const &Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
      int const &Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
    <BLANKLINE>
      // Set up cell-centered Cartesian coordinate grid, centered at the origin.
      gd_host->xx[0] = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS0);
      gd_host->xx[1] = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS1);
      gd_host->xx[2] = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS2);
    } // END FUNCTION mallocHostgrid
    <BLANKLINE>
    """
    desc = r"""Allocate griddata_struct[grid].xx for host."""
    name = "mallocHostgrid"
    params_dict = {
        "commondata": "const commondata_struct *restrict",
        "params": "const params_struct *restrict",
        "gd_host": "griddata_struct *restrict",
        "gd_gpu": "const griddata_struct *restrict",
    }
    body = """
  int const& Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  int const& Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  int const& Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  // Set up cell-centered Cartesian coordinate grid, centered at the origin.
  gd_host->xx[0] = (REAL*) malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS0);
  gd_host->xx[1] = (REAL*) malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS1);
  gd_host->xx[2] = (REAL*) malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS2);
"""
    kernel = GPU_Kernel(body, params_dict, name, decorators="__host__", comments=desc)
    return kernel.CFunction.full_function


def register_CFunction_cpyDevicetoHost__grid() -> None:
    """
    Register C function for copying grid from device to host.

    DOCTEST:
    >>> import nrpy.c_function as cfc
    >>> register_CFunction_cpyDevicetoHost__grid()
    >>> print(cfc.CFunction_dict['cpyDevicetoHost__grid'].full_function)
    #include "../BHaH_defines.h"
    /**
     * Kernel: mallocHostgrid.
     * Allocate griddata_struct[grid].xx for host.
     */
    __host__ static void mallocHostgrid(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                        griddata_struct *restrict gd_host, const griddata_struct *restrict gd_gpu) {
    <BLANKLINE>
      int const &Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
      int const &Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
      int const &Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
    <BLANKLINE>
      // Set up cell-centered Cartesian coordinate grid, centered at the origin.
      gd_host->xx[0] = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS0);
      gd_host->xx[1] = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS1);
      gd_host->xx[2] = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS2);
    } // END FUNCTION mallocHostgrid
    <BLANKLINE>
    /**
     * Copy griddata_struct[grid].xx from GPU to host.
     */
    __host__ void cpyDevicetoHost__grid(const commondata_struct *restrict commondata, griddata_struct *restrict gd_host,
                                        const griddata_struct *restrict gd_gpu) {
      for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
        const params_struct *restrict params = &gd_gpu[grid].params;
        int const &Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
        int const &Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
        int const &Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
    <BLANKLINE>
        mallocHostgrid(commondata, params, &gd_host[grid], gd_gpu);
        cudaMemcpy(gd_host[grid].xx[0], gd_gpu[grid].xx[0], sizeof(REAL) * Nxx_plus_2NGHOSTS0, cudaMemcpyDeviceToHost);
        cudaMemcpy(gd_host[grid].xx[1], gd_gpu[grid].xx[1], sizeof(REAL) * Nxx_plus_2NGHOSTS1, cudaMemcpyDeviceToHost);
        cudaMemcpy(gd_host[grid].xx[2], gd_gpu[grid].xx[2], sizeof(REAL) * Nxx_plus_2NGHOSTS2, cudaMemcpyDeviceToHost);
      } // END LOOP over grids
    } // END FUNCTION cpyDevicetoHost__grid
    <BLANKLINE>
    """
    includes = ["BHaH_defines.h"]
    prefunc = generate_CFunction_mallocHostgrid()
    desc = r"""Copy griddata_struct[grid].xx from GPU to host."""
    cfunc_type = "__host__ void"
    name = "cpyDevicetoHost__grid"
    params = "const commondata_struct *restrict commondata, griddata_struct *restrict gd_host, "
    params += "const griddata_struct *restrict gd_gpu"
    body = """for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
  const params_struct *restrict params = &gd_gpu[grid].params;
  int const& Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  int const& Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  int const& Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  mallocHostgrid(commondata, params, &gd_host[grid], gd_gpu);
  cudaMemcpy(gd_host[grid].xx[0], gd_gpu[grid].xx[0], sizeof(REAL) * Nxx_plus_2NGHOSTS0, cudaMemcpyDeviceToHost);
  cudaMemcpy(gd_host[grid].xx[1], gd_gpu[grid].xx[1], sizeof(REAL) * Nxx_plus_2NGHOSTS1, cudaMemcpyDeviceToHost);
  cudaMemcpy(gd_host[grid].xx[2], gd_gpu[grid].xx[2], sizeof(REAL) * Nxx_plus_2NGHOSTS2, cudaMemcpyDeviceToHost);
} // END LOOP over grids
"""
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        subdirectory="CUDA_utils",
    )


def register_CFunction_CUDA__malloc_host_gfs() -> None:
    """
    Register C function for allocating sufficient Host storage for diagnostics GFs.

    DOCTEST:
    >>> import nrpy.c_function as cfc
    >>> register_CFunction_CUDA__malloc_host_gfs()
    >>> print(cfc.CFunction_dict['CUDA__malloc_host_gfs'].full_function)
    #include "../BHaH_defines.h"
    /**
     * Allocate Host storage for diagnostics GFs.
     */
    __host__ void CUDA__malloc_host_gfs(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                        MoL_gridfunctions_struct *restrict gridfuncs) {
    <BLANKLINE>
      int const &Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
      int const &Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
      int const &Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
      const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
    <BLANKLINE>
      cudaMallocHost((void **)&gridfuncs->y_n_gfs, sizeof(REAL) * Nxx_plus_2NGHOSTS_tot * NUM_EVOL_GFS);
      cudaCheckErrors(cudaMallocHost, "Malloc y_n diagnostic GFs failed.");
    <BLANKLINE>
      cudaMallocHost((void **)&gridfuncs->diagnostic_output_gfs, sizeof(REAL) * Nxx_plus_2NGHOSTS_tot * NUM_AUX_GFS);
      cudaCheckErrors(cudaMallocHost, "Malloc diagnostic GFs failed.")
    } // END FUNCTION CUDA__malloc_host_gfs
    <BLANKLINE>
    """
    includes = ["BHaH_defines.h"]
    desc = r"""Allocate Host storage for diagnostics GFs."""
    cfunc_type = "__host__ void"
    name = "CUDA__malloc_host_gfs"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, "
    params += "MoL_gridfunctions_struct *restrict gridfuncs"
    body = """
  int const& Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  int const& Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  int const& Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

  cudaMallocHost((void**)&gridfuncs->y_n_gfs, sizeof(REAL) * Nxx_plus_2NGHOSTS_tot * NUM_EVOL_GFS);
  cudaCheckErrors(cudaMallocHost, "Malloc y_n diagnostic GFs failed.");

  cudaMallocHost((void**)&gridfuncs->diagnostic_output_gfs, sizeof(REAL) * Nxx_plus_2NGHOSTS_tot * NUM_AUX_GFS);
  cudaCheckErrors(cudaMallocHost, "Malloc diagnostic GFs failed.")
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        subdirectory="CUDA_utils",
    )


def register_CFunction_CUDA__free_host_gfs() -> None:
    """
    Register C function for allocating sufficient Host storage for diagnostics GFs.

    DOCTEST:
    >>> import nrpy.c_function as cfc
    >>> register_CFunction_CUDA__free_host_gfs()
    >>> print(cfc.CFunction_dict['CUDA__free_host_gfs'].full_function)
    #include "../BHaH_defines.h"
    /**
     * Free Host storage for diagnostics GFs.
     */
    __host__ void CUDA__free_host_gfs(MoL_gridfunctions_struct *gridfuncs) {
    <BLANKLINE>
      cudaFreeHost(gridfuncs->y_n_gfs);
      cudaCheckErrors(free, "Host-ynFree failed");
      cudaFreeHost(gridfuncs->diagnostic_output_gfs);
      cudaCheckErrors(free, "Host-non-ynFree failed");
    } // END FUNCTION CUDA__free_host_gfs
    <BLANKLINE>
    """
    includes = ["BHaH_defines.h"]
    desc = r"""Free Host storage for diagnostics GFs."""
    cfunc_type = "__host__ void"
    name = "CUDA__free_host_gfs"
    params = "MoL_gridfunctions_struct * gridfuncs"
    body = """
    cudaFreeHost(gridfuncs->y_n_gfs);
    cudaCheckErrors(free, "Host-ynFree failed");
    cudaFreeHost(gridfuncs->diagnostic_output_gfs);
    cudaCheckErrors(free, "Host-non-ynFree failed");
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        subdirectory="CUDA_utils",
    )


def register_CFunction_cpyDevicetoHost__gf() -> None:
    """
    Register C function for asynchronously copying data from device to host.

    DOCTEST:
    >>> import nrpy.c_function as cfc
    >>> register_CFunction_cpyDevicetoHost__gf()
    >>> print(cfc.CFunction_dict['cpyDevicetoHost__gf'].full_function)
    #include "../BHaH_defines.h"
    /**
     * Asynchronously copying a grid function from device to host.
     */
    __host__ size_t cpyDevicetoHost__gf(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *gf_host,
                                        const REAL *gf_gpu, const int host_GF_IDX, const int gpu_GF_IDX) {
    <BLANKLINE>
      int const Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
      int const Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
      int const Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
      const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
    <BLANKLINE>
      size_t streamid = (params->grid_idx + gpu_GF_IDX) % NUM_STREAMS;
      int offset_gpu = Nxx_plus_2NGHOSTS_tot * gpu_GF_IDX;
      int offset_host = Nxx_plus_2NGHOSTS_tot * host_GF_IDX;
      cudaMemcpyAsync(&gf_host[offset_host], &gf_gpu[offset_gpu], sizeof(REAL) * Nxx_plus_2NGHOSTS_tot, cudaMemcpyDeviceToHost, streams[streamid]);
      cudaCheckErrors(cudaMemcpyAsync, "Copy of gf data failed");
      return streamid;
    } // END FUNCTION cpyDevicetoHost__gf
    <BLANKLINE>
    """
    includes = ["BHaH_defines.h"]
    desc = r"""Asynchronously copying a grid function from device to host."""
    cfunc_type = "__host__ size_t"
    name = "cpyDevicetoHost__gf"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, "
    params += "REAL * gf_host, const REAL * gf_gpu, const int host_GF_IDX, const int gpu_GF_IDX"
    body = """
  int const Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  int const Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  int const Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

  size_t streamid = (params->grid_idx + gpu_GF_IDX) % NUM_STREAMS;
  int offset_gpu  = Nxx_plus_2NGHOSTS_tot * gpu_GF_IDX;
  int offset_host = Nxx_plus_2NGHOSTS_tot * host_GF_IDX;
  cudaMemcpyAsync(&gf_host[offset_host],
                  &gf_gpu[offset_gpu],
                  sizeof(REAL) * Nxx_plus_2NGHOSTS_tot,
                  cudaMemcpyDeviceToHost, streams[streamid]);
  cudaCheckErrors(cudaMemcpyAsync, "Copy of gf data failed");
  return streamid;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        subdirectory="CUDA_utils",
    )


def register_CFunction_cpyHosttoDevice__gf() -> None:
    """
    Register C function for asynchronously copying data from host to device.

    DOCTEST:
    >>> import nrpy.c_function as cfc
    >>> register_CFunction_cpyHosttoDevice__gf()
    >>> print(cfc.CFunction_dict['cpyDevicetoHost__gf'].full_function)
    #include "../BHaH_defines.h"
    /**
     * Asynchronously copying a grid function from device to host.
     */
    __host__ size_t cpyDevicetoHost__gf(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *gf_host,
                                        const REAL *gf_gpu, const int host_GF_IDX, const int gpu_GF_IDX) {
    <BLANKLINE>
      int const Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
      int const Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
      int const Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
      const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
    <BLANKLINE>
      size_t streamid = (params->grid_idx + gpu_GF_IDX) % NUM_STREAMS;
      int offset_gpu = Nxx_plus_2NGHOSTS_tot * gpu_GF_IDX;
      int offset_host = Nxx_plus_2NGHOSTS_tot * host_GF_IDX;
      cudaMemcpyAsync(&gf_host[offset_host], &gf_gpu[offset_gpu], sizeof(REAL) * Nxx_plus_2NGHOSTS_tot, cudaMemcpyDeviceToHost, streams[streamid]);
      cudaCheckErrors(cudaMemcpyAsync, "Copy of gf data failed");
      return streamid;
    } // END FUNCTION cpyDevicetoHost__gf
    <BLANKLINE>
    """
    includes = ["BHaH_defines.h"]
    desc = r"""Asynchronously copying a grid function from host to device."""
    cfunc_type = "__host__ size_t"
    name = "cpyHosttoDevice__gf"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, "
    params += "const REAL * gf_host, REAL * gf_gpu, const int host_GF_IDX, const int gpu_GF_IDX"
    body = """
  int const Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  int const Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  int const Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

  size_t streamid = (params->grid_idx + gpu_GF_IDX) % NUM_STREAMS;
  int offset_gpu  = Nxx_plus_2NGHOSTS_tot * gpu_GF_IDX;
  int offset_host = Nxx_plus_2NGHOSTS_tot * host_GF_IDX;
  cudaMemcpyAsync(&gf_gpu[offset_host],
                  &gf_host[offset_gpu],
                  sizeof(REAL) * Nxx_plus_2NGHOSTS_tot,
                  cudaMemcpyHostToDevice, streams[streamid]);
  cudaCheckErrors(cudaMemcpyAsync, "Copy of gf data failed");
  return streamid;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        subdirectory="CUDA_utils",
    )


def register_CFunctions_HostDevice__operations() -> None:
    """Generate all of the Host to/from Device specific functions."""
    register_CFunction_cpyHosttoDevice_commondata__constant()
    register_CFunction_cpyHosttoDevice_params__constant()
    register_CFunction_cpyHosttoDevice__gf()

    register_CFunction_cpyDevicetoHost__grid()
    register_CFunction_cpyDevicetoHost__gf()
    register_CFunction_CUDA__malloc_host_gfs()
    register_CFunction_CUDA__free_host_gfs()


# Define implemented reductions
def minimum_reduction(cmp_var: str, reduction_var: str = "local_reduced") -> str:
    """
    Operator for minimum reduction.

    :param cmp_var: String of variable to compare against
    :param reduction_var: String of reduction variable
    :returns: Reduction operation string
    """
    return f"""
if({reduction_var} > {cmp_var}) {{
    {reduction_var} = {cmp_var};
}}
"""


def maximum_reduction(cmp_var: str, reduction_var: str = "local_reduced") -> str:
    """
    Operator for maximum reduction.

    :param cmp_var: String of variable to compare against
    :param reduction_var: String of reduction variable
    :returns: Reduction operation string
    """
    return f"""
if({reduction_var} < {cmp_var}) {{
    {reduction_var} = {cmp_var};
}}
"""


def sum_reduction(cmp_var: str, reduction_var: str = "local_reduced") -> str:
    """
    Operator for sum reduction.

    :param cmp_var: String of variable to compare against
    :param reduction_var: String of reduction variable
    :returns: Reduction operation string
    """
    return f"{reduction_var} += {cmp_var};"


# End reduction definitions

# Store local reductions for later use in substitution
implemented_reduction_dict = {
    "sum": sum_reduction,
    "minimum": minimum_reduction,
    "maximum": maximum_reduction,
}


class CUDA_reductions:
    r"""
    Provides a template for generating CUDA compatible reductions.

    :param reduction_type: Specify the reduction type to generate
    :param cfunc_decorators: Set decorators for the reduction (e.g. template, __host__, inline)
    :param cfunc_type: Return type of the function

    DOCTEST:
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.params as par
    >>> for fp_type in ['float', 'double']:
    ...     par.set_parval_from_str("fp_type", fp_type)
    ...     for reduction_type, _reduction_func in implemented_reduction_dict.items():
    ...         reduction = CUDA_reductions(reduction_type=reduction_type)
    ...         reduction.generate_CFunction()
    ...         generated_str = reduction.CFunction.full_function
    ...         validate_strings(generated_str, f"CUDA_reductions_{fp_type}__{reduction_type}")
    """

    def __init__(
        self,
        reduction_type: str = "minimum",
        cfunc_decorators: str = "__host__",
        cfunc_type: str = "REAL",
    ) -> None:
        self.reduction_type = reduction_type
        self.cfunc_decorators = cfunc_decorators
        self.cfunc_type = cfunc_type
        self.fp_type = par.parval_from_str("fp_type")

        self.type_dict = {
            "double": "unsigned long long",
            "float": "unsigned",
        }

        self.initial_value_dict = {
            "sum": "0U",
            "maximum": "0U",
            "minimum": "0xFFFFFFFFU",
        }

        self.reduction_dict = implemented_reduction_dict

        self.atomic_operations = {
            "sum": "atomicAdd",
            "minimum": "atomicMin",
            "maximum": "atomicMax",
        }

        if self.reduction_type not in self.initial_value_dict:
            raise ValueError(
                f"{self.reduction_type} is not a defined reduction. Choose from {self.initial_value_dict.keys()}"
            )
        self.includes = ["BHaH_defines.h"]
        self.desc = f"""Find array global {self.reduction_type}."""
        self.name = f"find_global__{self.reduction_type}"
        self.params = "REAL * data, uint const data_length"

        self.kernel_name = f"{self.name}__cuda"
        self.recast_type = (
            f"({self.type_dict[self.fp_type]} int *)"
            if self.reduction_type != "sum"
            else ""
        )
        self.prefunc = f"""
__global__
static void {self.kernel_name}(REAL * data, REAL * min, uint const data_length) {{
    // shared data between all warps
    // Assumes one block = 32 warps = 32 * 32 threads
    // As of today, the standard maximum threads per
    // block is 1024 = 32 * 32
    __shared__ REAL shared_data[32];

    // Some initial value
    REAL REDUCTION_LIMIT = (REAL) {self.initial_value_dict[self.reduction_type]};

    // Global data index - expecting a 1D dataset
    uint idx = threadIdx.x + blockDim.x * blockIdx.x;

    // thread index
    uint tid = threadIdx.x;

    // local thread reduction
    REAL local_reduced = REDUCTION_LIMIT;

    // warp mask - says all threads are involved in shuffle
    // 0xFFFFFFFFU in binary is 32 1's.
    unsigned mask = 0xFFFFFFFFU;

    // lane = which thread am I in the warp
    uint lane = threadIdx.x % warpSize;
    // warpID = which warp am I in the block
    uint warpID = threadIdx.x / warpSize;

    // Stride through data for each thread
    while(idx < data_length) {{
        {self.reduction_dict[self.reduction_type]("data[idx]")}
        // idx stride
        idx += gridDim.x * blockDim.x;
    }}
    // Shuffle down kernel
    for(int offset = warpSize / 2; offset > 0; offset >>= 1) {{
        REAL shfl = __shfl_down_sync(mask, local_reduced, offset);
        {self.reduction_dict[self.reduction_type]("shfl")}
    }}
    // Shuffle results in lane 0 have the shuffle result
    if(lane == 0) {{
        shared_data[warpID] = local_reduced;
    }}

    // Make sure all warps in the block are synchronized
    __syncthreads();

    // Since there is only 32 partial reductions, we only
    // have one warp worth of work
    if(warpID == 0) {{
        // Check to make sure we had 32 blocks of data
        if(tid < blockDim.x / warpSize) {{
            local_reduced = shared_data[lane];
        }} else {{
            local_reduced = REDUCTION_LIMIT;
        }}

        // Shuffle down kernel
        for(int offset = warpSize / 2; offset > 0; offset >>= 1) {{
            REAL shfl = __shfl_down_sync(mask, local_reduced, offset);
            {self.reduction_dict[self.reduction_type]("shfl")}
        }}
        if(tid == 0) {{
            {self.atomic_operations[self.reduction_type]}({self.recast_type}min, *({self.recast_type}&local_reduced));
        }}
    }}
}}
"""

        self.body = f"""
    // This can be tested up to 1024
    uint threadCount = 32;

    // Number of blocks
    uint blockCount = (data_length + threadCount - 1) / threadCount;

    // CUDA atomics other than cas are only
    // compatible with (u)int.  To be generic
    // we use {self.type_dict[self.fp_type]} to be able to handle
    // {self.fp_type} precision variables
    using ull = REAL;
    ull * h_reduced = (ull*)malloc(sizeof(ull));
    ull * d_reduced;
    *h_reduced = (ull){self.initial_value_dict[self.reduction_type]};

    cudaMalloc(&d_reduced, sizeof(ull));
    cudaCheckErrors(cudaMalloc, "cudaMalloc failure"); // error checking

    cudaMemcpy(d_reduced, h_reduced, sizeof(ull), cudaMemcpyHostToDevice);
    cudaCheckErrors(cudaMemcpy, "cudaCopyTo failure"); // error checking

    {self.name}__cuda<<<blockCount, threadCount>>>(data, d_reduced, data_length);
    cudaCheckErrors(find_min_cu, "cudaKernel - find_min_cu failed"); // error checking

    cudaMemcpy(h_reduced, d_reduced, sizeof(REAL), cudaMemcpyDeviceToHost);
    cudaCheckErrors(cudaMemcpy, "cudaCopyFrom failure"); // error checking

    cudaFree(d_reduced);
    cudaCheckErrors(cudaFree, "cudaFree failure"); // error checking

    // Recast back to result pointer type
    REAL * res = (REAL *) h_reduced;
    return *res;
"""
        self.CFunction: cfc.CFunction

    def generate_CFunction(self) -> None:
        """
        Generate CFunction from initialized parameters.
        For testing only.
        """
        self.CFunction = cfc.CFunction(
            prefunc=self.prefunc,
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            CoordSystem_for_wrapper_func="",
            name=self.name,
            params=self.params,
            include_CodeParameters_h=False,
            body=self.body,
            subdirectory="CUDA_utils",
        )


def register_CFunction_find_global_minimum() -> None:
    """Register C function for finding the global minimum of an array."""
    reduction = CUDA_reductions(
        reduction_type="minimum",
        cfunc_decorators="__host__",
        cfunc_type="REAL",
    )

    cfc.register_CFunction(
        prefunc=reduction.prefunc,
        includes=reduction.includes,
        desc=reduction.desc,
        cfunc_type=reduction.cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=reduction.name,
        params=reduction.params,
        include_CodeParameters_h=False,
        body=reduction.body,
        subdirectory="CUDA_utils",
    )


def register_CFunction_find_global_sum() -> None:
    """Register C function for finding the global sum of an array."""
    reduction = CUDA_reductions(
        reduction_type="sum",
        cfunc_decorators="__host__",
        cfunc_type="REAL",
    )

    cfc.register_CFunction(
        prefunc=reduction.prefunc,
        includes=reduction.includes,
        desc=reduction.desc,
        cfunc_type=reduction.cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=reduction.name,
        params=reduction.params,
        include_CodeParameters_h=False,
        body=reduction.body,
        subdirectory="CUDA_utils",
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
