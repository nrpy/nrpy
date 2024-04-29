"""
Base class for generating the C main() function for all codes in the BHaH infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot** com
        Samuel D. Tootle
        sdtootle **at** gmail **dot** com
"""

from typing import List, Tuple
import nrpy.c_function as cfc

class base_register_CFunction_main_c:
    """
    Base class for generating the "generic" C main() function for all simulation codes in the BHaH infrastructure.

    :param MoL_method: Method of Lines algorithm used to step forward in time.
    :param initial_data_desc: Description for initial data, default is an empty string.
    :param boundary_conditions_desc: Description of the boundary conditions, default is an empty string.
    :param prefunc: String that appears before main(). DO NOT populate this, except when debugging, default is an empty string.
    :param initialize_constant_auxevol: If set to True, `initialize_constant_auxevol` function will be called during the simulation initialization phase to set these constants. Default is False.
    :param pre_MoL_step_forward_in_time: Code for handling pre-right-hand-side operations, default is an empty string.
    :param post_MoL_step_forward_in_time: Code for handling post-right-hand-side operations, default is an empty string.
    :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".
    :raises ValueError: Raised if any required function for BHaH main() is not registered.
    """    
    def __init__(
        self,
        MoL_method: str,
        initial_data_desc: str = "",
        boundary_conditions_desc: str = "",
        prefunc: str = "",
        initialize_constant_auxevol: bool = False,
        pre_MoL_step_forward_in_time: str = "",
        post_MoL_step_forward_in_time: str = "",
        clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
    ) -> None:
        self.MoL_method=MoL_method
        self.initial_data_desc=initial_data_desc
        self.boundary_conditions_desc=boundary_conditions_desc
        self.prefunc=prefunc
        self.initialize_constant_auxevol=initialize_constant_auxevol
        self.pre_MoL_step_forward_in_time=pre_MoL_step_forward_in_time
        self.post_MoL_step_forward_in_time=post_MoL_step_forward_in_time
        self.clang_format_options=clang_format_options
        
        self.initial_data_desc += " "
        # Make sure all required C functions are registered
        missing_functions: List[Tuple[str, str]] = []
        for func_tuple in [
            ("params_struct_set_to_default", "CodeParameters.py"),
            ("commondata_struct_set_to_default", "CodeParameters.py"),
            ("cmdline_input_and_parfile_parser", "cmdline_input_and_parfiles.py"),
            (
                "numerical_grids_and_timestep",
                "e.g., numerical_grids_and_timestep.py or user defined",
            ),
            ("MoL_malloc_y_n_gfs", "MoL.py"),
            ("MoL_malloc_non_y_n_gfs", "MoL.py"),
            ("initial_data", "initial_data.py"),
            ("MoL_step_forward_in_time", "MoL.py"),
            ("diagnostics", "diagnostics.py"),
            ("MoL_free_memory_y_n_gfs", "MoL.py"),
            ("MoL_free_memory_non_y_n_gfs", "MoL.py"),
        ]:
            if func_tuple[0] not in cfc.CFunction_dict:
                missing_functions += [func_tuple]
        if missing_functions:
            error_msg = "Error: These functions are required for all BHaH main() functions, and are not registered.\n"
            for func_tuple in missing_functions:
                error_msg += (
                    f'  {func_tuple[0]}, registered by function within "{func_tuple[1]}"\n'
                )
            raise ValueError(error_msg)

        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
        self.desc = """-={ main() function }=-
Step 1.a: Set each commondata CodeParameter to default.
Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
Step 1.c: Allocate NUMGRIDS griddata structs, each containing data specific to an individual grid.
Step 1.d: Set each CodeParameter in griddata.params to default.
Step 1.e: Set up numerical grids: xx[3], masks, Nxx, dxx, invdxx, bcstruct, rfm_precompute, timestep, etc.
Step 2: Initial data are set on y_n_gfs gridfunctions. Allocate storage for them first.
Step 3: Finalize initialization: set up {self.initial_data_desc}initial data, etc.
Step 4: Allocate storage for non-y_n gridfunctions, needed for the Runge-Kutta-like timestepping.
"""
        if self.initialize_constant_auxevol:
            self.desc += "Step 4.a: Set AUXEVOL gridfunctions that will never change in time."
        self.desc += f"""Step 5: MAIN SIMULATION LOOP
- Step 5.a: Output diagnostics.
- Step 5.b: Prepare to step forward in time.
- Step 5.c: Step forward in time using Method of Lines with {MoL_method} algorithm, applying {boundary_conditions_desc} boundary conditions.
- Step 5.d: Finish up step in time.
Step 6: Free all allocated memory."""
        self.cfunc_type = "int"
        self.name = "main"
        self.params = "int argc, const char *argv[]"
        self.body = r"""  commondata_struct commondata; // commondata contains parameters common to all grids.
  griddata_struct *restrict griddata; // griddata contains data specific to an individual grid.

// Step 1.a: Set each commondata CodeParameter to default.
commondata_struct_set_to_default(&commondata);

// Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
cmdline_input_and_parfile_parser(&commondata, argc, argv);

// Step 1.c: Allocate NUMGRIDS griddata arrays, each containing data specific to an individual grid.
griddata = (griddata_struct *restrict)malloc(sizeof(griddata_struct)*commondata.NUMGRIDS);

// Step 1.d: Set each CodeParameter in griddata.params to default.
params_struct_set_to_default(&commondata, griddata);

// Step 1.e: Set up numerical grids: xx[3], masks, Nxx, dxx, invdxx, bcstruct, rfm_precompute, timestep, etc.
{
  // if calling_for_first_time, then initialize commondata time=nn=t_0=nn_0 = 0
  const bool calling_for_first_time = true;
  numerical_grids_and_timestep(&commondata, griddata, calling_for_first_time);
}

for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
  // Step 2: Initial data are set on y_n_gfs gridfunctions. Allocate storage for them first.
  MoL_malloc_y_n_gfs(&commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
}

// Step 3: Finalize initialization: set up initial data, etc.
initial_data(&commondata, griddata);

// Step 4: Allocate storage for non-y_n gridfunctions, needed for the Runge-Kutta-like timestepping
for(int grid=0; grid<commondata.NUMGRIDS; grid++)
  MoL_malloc_non_y_n_gfs(&commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
"""
        if self.initialize_constant_auxevol:
            self.body += """// Step 4.a: Set AUXEVOL gridfunctions that will never change in time.
initialize_constant_auxevol(&commondata, griddata);
"""
        self.body += """
// Step 5: MAIN SIMULATION LOOP
while(commondata.time < commondata.t_final) { // Main loop to progress forward in time.
  // Step 5.a: Main loop, part 1: Output diagnostics
  diagnostics(&commondata, griddata);

  // Step 5.b: Main loop, part 2 (pre_MoL_step_forward_in_time): Prepare to step forward in time
"""
        if self.pre_MoL_step_forward_in_time != "":
            self.body += self.pre_MoL_step_forward_in_time
        else:
            self.body += "  // (nothing here; specify by setting pre_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"
        self.body += f"""
  // Step 5.c: Main loop, part 3: Step forward in time using Method of Lines with {MoL_method} algorithm,
  //           applying {self.boundary_conditions_desc} boundary conditions.
  MoL_step_forward_in_time(&commondata, griddata);

  // Step 5.d: Main loop, part 4 (post_MoL_step_forward_in_time): Finish up step in time
"""
        if self.post_MoL_step_forward_in_time != "":
            self.body += self.post_MoL_step_forward_in_time
        else:
            self.body += "  // (nothing here; specify by setting post_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"
        self.body += r"""
} // End main loop to progress forward in time.

// Step 6: Free all allocated memory
{
  const bool enable_free_non_y_n_gfs=true;
  griddata_free(&commondata, griddata, enable_free_non_y_n_gfs);
}
return 0;
"""