"""
Class to register griddata_free method for OpenMP based applications.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot** com
        Samuel D. Tootle
        sdtootle **at** gmail **dot** com
"""

import nrpy.infrastructures.BHaH.grid_management.base_griddata_free as base_free
class register_CFunction_griddata_free(base_free.base_register_CFunction_griddata_free):
    """
    Register the C function griddata_free() to free all memory within the griddata struct.
    Overload based on OpenMP parallelization

    :param enable_rfm_precompute: A flag to enable/disable rfm_precompute_free within the C function body.
    :param enable_CurviBCs: A flag to enable/disable freeing CurviBCs within the C function body.
    """    
    def __init__(
        self,
        enable_rfm_precompute: bool,
        enable_CurviBCs: bool,
    ) -> None:
        super().__init__()
        self.body = r"""for(int grid=0;grid<commondata->NUMGRIDS;grid++) {
"""
        if enable_rfm_precompute:
            self.body += "  rfm_precompute_free(commondata, &griddata[grid].params, &griddata[grid].rfmstruct);\n"
        
        if enable_CurviBCs:
            self.body += r"""
  free(griddata[grid].bcstruct.inner_bc_array);
  for(int ng=0;ng<NGHOSTS*3;ng++) free(griddata[grid].bcstruct.pure_outer_bc_array[ng]);
"""
        self.body += r"""

  MoL_free_memory_y_n_gfs(&griddata[grid].gridfuncs);
  if(enable_free_non_y_n_gfs)
  for(int i=0;i<3;i++) free(griddata[grid].xx[i]);
} // END for(int grid=0;grid<commondata->NUMGRIDS;grid++)
"""
        self.body += "free(griddata);\n"
        
        self.register_CFunction()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
