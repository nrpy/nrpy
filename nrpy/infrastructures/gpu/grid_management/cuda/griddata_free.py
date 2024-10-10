"""
Class to register griddata_free method for CUDA based applications.

Author: Samuel D. Tootle
        sdtootle **at** gmail **dot** com
"""

import nrpy.infrastructures.gpu.grid_management.base_griddata_free as base_free


class register_CFunction_griddata_free(base_free.base_register_CFunction_griddata_free):
    """
    Register the C function griddata_free() to free all memory within the griddata struct.
    Overload based on CUDA parallelization

    :param enable_rfm_precompute: A flag to enable/disable rfm_precompute_free within the C function body.
    :param enable_CurviBCs: A flag to enable/disable freeing CurviBCs within the C function body.
    """

    def __init__(
        self,
        enable_rfm_precompute: bool,
        enable_CurviBCs: bool,
    ) -> None:
        super().__init__()
        self.params = "const commondata_struct *restrict commondata, griddata_struct *restrict griddata, griddata_struct *restrict griddata_host, const bool free_non_y_n_gfs_and_core_griddata_pointers"
        self.body = r"""for(int grid=0;grid<commondata->NUMGRIDS;grid++) {
"""
        if enable_rfm_precompute:
            self.body += "  rfm_precompute_free(commondata, &griddata[grid].params, &griddata[grid].rfmstruct);\n"
            self.body += '  cudaCheckErrors(free, "rfmFree failed")\n;'

        if enable_CurviBCs:
            self.body += r"""
  cudaFree(griddata[grid].bcstruct.inner_bc_array);
  cudaCheckErrors(free, "bcstruct inner_bc_arrayFree failed");
  for(int ng=0;ng<NGHOSTS*3;ng++) {
    cudaFree(griddata[grid].bcstruct.pure_outer_bc_array[ng]);
    cudaCheckErrors(free, "bcstruct pure_outer_bc_arrayFree failed");
  }
"""
        self.body += r"""

  MoL_free_memory_y_n_gfs(&griddata[grid].gridfuncs);
  cudaCheckErrors(free, "MoLFree failed");
  CUDA__free_host_gfs(&griddata_host[grid].gridfuncs);
  if(free_non_y_n_gfs_and_core_griddata_pointers) {
    MoL_free_memory_non_y_n_gfs(&griddata[grid].gridfuncs);
  }
  for(int i=0;i<3;i++) {
    cudaFree(griddata[grid].xx[i]);
    cudaCheckErrors(free, "griddata-xxFree failed");
    free(griddata_host[grid].xx[i]);
  }

} // END for(int grid=0;grid<commondata->NUMGRIDS;grid++)
"""
        self.body += """
  if(free_non_y_n_gfs_and_core_griddata_pointers) {
    free(griddata);
    free(griddata_host);
  }
"""
        self.register()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
