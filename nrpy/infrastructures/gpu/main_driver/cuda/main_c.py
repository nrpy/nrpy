"""
Generate the C main() function for all codes in the BHaH infrastructure using CUDA.

Author: Samuel D. Tootle
        sdtootle **at** gmail **dot** com
        Zachariah B. Etienne
        zachetie **at** gmail **dot** com
"""

import nrpy.infrastructures.gpu.main_driver.base_main_c as base_main


class register_CFunction_main_c(base_main.base_register_CFunction_main_c):
    """
    Generate the "generic" C main() function for all simulation codes in the BHaH infrastructure.

    :param MoL_method: Method of Lines algorithm used to step forward in time.
    :param initial_data_desc: Description for initial data, default is an empty string.
    :param set_initial_data_after_auxevol_malloc: Set to True if initial data for y_n_gfs requires auxevol variables be set, e.g., in case of GRHD, initial data must set primitives.
    :param boundary_conditions_desc: Description of the boundary conditions, default is an empty string.
    :param prefunc: String that appears before main(). DO NOT populate this, except when debugging, default is an empty string.
    :param post_non_y_n_auxevol_mallocs: Function calls after memory is allocated for non y_n and auxevol gridfunctions, default is an empty string.
    :param pre_diagnostics: Function calls prior to diagnostics; e.g., regridding. Default is an empty string.
    :param pre_MoL_step_forward_in_time: Code for handling pre-right-hand-side operations, default is an empty string.
    :param post_MoL_step_forward_in_time: Code for handling post-right-hand-side operations, default is an empty string.
    :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".
    :raises ValueError: Raised if any required function for BHaH main() is not registered.
    """

    def __init__(
        self,
        MoL_method: str,
        initial_data_desc: str = "",
        set_initial_data_after_auxevol_malloc: bool = False,
        boundary_conditions_desc: str = "",
        prefunc: str = "",
        post_non_y_n_auxevol_mallocs: str = "",
        pre_diagnostics: str = "",
        pre_MoL_step_forward_in_time: str = "",
        post_MoL_step_forward_in_time: str = "",
        clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
    ) -> None:
        super().__init__(
            MoL_method,
            initial_data_desc=initial_data_desc,
            set_initial_data_after_auxevol_malloc=set_initial_data_after_auxevol_malloc,
            boundary_conditions_desc=boundary_conditions_desc,
            prefunc=prefunc,
            post_non_y_n_auxevol_mallocs=post_non_y_n_auxevol_mallocs,
            pre_diagnostics=pre_diagnostics,
            pre_MoL_step_forward_in_time=pre_MoL_step_forward_in_time,
            post_MoL_step_forward_in_time=post_MoL_step_forward_in_time,
            clang_format_options=clang_format_options,
        )
        self.includes += ["BHaH_gpu_global_defines.h"]
        self.body = r"""
#include "BHaH_gpu_global_init.h"
commondata_struct commondata; // commondata contains parameters common to all grids.
griddata_struct *restrict griddata; // griddata contains data specific to an individual grid.
griddata_struct *restrict griddata_host; // stores only the host data needed for diagnostics

// Step 1.a: Initialize each CodeParameter in the commondata struc to its default value.
commondata_struct_set_to_default(&commondata);

// Step 1.b: Overwrite the default values with those from the parameter file.
//           Then overwrite the parameter file values with those provided via command line arguments.
cmdline_input_and_parfile_parser(&commondata, argc, argv);

// Step 1.c: Allocate memory for MAXNUMGRIDS griddata structs,
//           where each structure contains data specific to an individual grid.
griddata = (griddata_struct *)malloc(sizeof(griddata_struct) * commondata.NUMGRIDS);
griddata_host = (griddata_struct *)malloc(sizeof(griddata_struct) * commondata.NUMGRIDS);

// Step 1.d: Set each CodeParameter in griddata.params to default.
params_struct_set_to_default(&commondata, griddata);

// Step 1.e: Set up numerical grids, including parameters such as NUMGRIDS, xx[3], masks, Nxx, dxx, invdxx,
//           bcstruct, rfm_precompute, timestep, and others.
{
  // If this function is being called for the first time, initialize commondata time, nn, t_0, and nn_0 to 0.
  const bool calling_for_first_time = true;
  numerical_grids_and_timestep(&commondata, griddata, griddata_host, calling_for_first_time);
}

for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
  // Step 2.a: Allocate storage for the initial data (y_n_gfs gridfunctions) on each grid.
  MoL_malloc_y_n_gfs(&commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
  // Step 2.b: Allocate host storage for diagnostics
  CUDA__malloc_host_gfs(&commondata, &griddata[grid].params, &griddata_host[grid].gridfuncs);
}
"""
        setup_initial_data_code = """Set up initial data.
if (!read_checkpoint(&commondata, griddata_host, griddata))
  initial_data(&commondata, griddata);
"""
        allocate_storage_code = """Allocate storage for non-y_n gridfunctions, needed for the Runge-Kutta-like timestepping.
for(int grid=0; grid<commondata.NUMGRIDS; grid++)
  MoL_malloc_non_y_n_gfs(&commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
"""
        step3code = setup_initial_data_code
        step4code = allocate_storage_code
        if self.set_initial_data_after_auxevol_malloc:
            step3code = allocate_storage_code
            step4code = setup_initial_data_code
        self.body += f"""
// Step 3: {step3code}

// Step 4: {step4code}
"""
        if self.post_non_y_n_auxevol_mallocs:
            self.body += f"""// Step 4.a: Functions called after memory for non-y_n and auxevol gridfunctions is allocated.
{self.post_non_y_n_auxevol_mallocs}
"""
        self.body += """
// Step 5: MAIN SIMULATION LOOP
while(commondata.time < commondata.t_final) { // Main loop to progress forward in time.
"""
        if self.pre_diagnostics:
            self.body += self.pre_diagnostics
        else:
            self.body += "// (nothing here; specify by setting pre_diagnostics string in register_CFunction_main_c().)\n"
        self.body += """
  // Step 5.b: Main loop, part 2: Output diagnostics
  diagnostics(&commondata, griddata, griddata_host);

  // Step 5.c: Main loop, part 3 (pre_MoL_step_forward_in_time): Prepare to step forward in time
"""
        if self.pre_MoL_step_forward_in_time != "":
            self.body += self.pre_MoL_step_forward_in_time
        else:
            self.body += "  // (nothing here; specify by setting pre_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"
        self.body += f"""
  // Step 5.d: Main loop, part 4: Step forward in time using Method of Lines with {MoL_method} algorithm,
  //           applying {self.boundary_conditions_desc} boundary conditions.
  MoL_step_forward_in_time(&commondata, griddata);

  // Step 5.e: Main loop, part 5 (post_MoL_step_forward_in_time): Finish up step in time
"""
        if self.post_MoL_step_forward_in_time != "":
            self.body += self.post_MoL_step_forward_in_time
        else:
            self.body += "  // (nothing here; specify by setting post_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"
        self.body += r"""
} // End main loop to progress forward in time.
// Make sure all workers are done
cudaDeviceSynchronize();
for(int i = 0; i < nstreams; ++i) {
    cudaStreamDestroy(streams[i]);
}
// Step 6: Free all allocated memory
{
  const bool enable_free_non_y_n_gfs=true;
  griddata_free(&commondata, griddata, griddata_host, enable_free_non_y_n_gfs);
}
return 0;
"""

        self.register()
