"""
Generate the C main() function for all codes in the BHaH infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List, Tuple

import nrpy.c_function as cfc


def register_CFunction_main_c(
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
    """
    Generate the "generic" C main() function for all simulation codes in the BHaH infrastructure.

    :param MoL_method: Method of Lines algorithm used to step forward in time.
    :param initial_data_desc: Description for initial data, default is an empty string.
    :param set_initial_data_after_auxevol_malloc: Set to True if initial data for y_n_gfs requires auxevol variables be set, e.g., in case of GRHD, initial data must set primitives.
    :param boundary_conditions_desc: Description of the boundary conditions, default is an empty string.
    :param prefunc: String that appears before main(). DO NOT populate this, EXCEPT when debugging, default is an empty string.
    :param post_non_y_n_auxevol_mallocs: Function calls after memory is allocated for non y_n and auxevol gridfunctions, default is an empty string.
    :param pre_diagnostics: Function calls prior to diagnostics; e.g., regridding. Default is an empty string.
    :param pre_MoL_step_forward_in_time: Function calls AFTER diagnostics and prior to each right-hand-side update, default is an empty string.
    :param post_MoL_step_forward_in_time: Function calls after each right-hand-side update, default is an empty string.
    :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".
    :raises ValueError: Raised if any required function for BHaH main() is not registered.
    """
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

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    allocate_auxevol_desc = "Allocate storage for non-y_n gridfunctions needed for the Runge-Kutta-like time-stepping."
    set_initial_data_desc = f"Set {initial_data_desc} initial data."
    step3desc = set_initial_data_desc
    step4desc = allocate_auxevol_desc
    if set_initial_data_after_auxevol_malloc:
        step3desc = allocate_auxevol_desc
        step4desc = set_initial_data_desc
    desc = f"""-={{ main() function }}=-
Step 1.a: Initialize each CodeParameter in the commondata struct to its default value.
Step 1.b: Overwrite the default values with those from the parameter file.
          Then overwrite the parameter file values with those provided via command line arguments.
Step 1.c: Allocate memory for MAXNUMGRIDS griddata structs,
          where each structure contains data specific to an individual grid.
Step 1.d: Initialize each CodeParameter in griddata.params to its default value.
Step 1.e: Set up numerical grids, including parameters such as NUMGRIDS, xx[3], masks, Nxx, dxx, invdxx,
          bcstruct, rfm_precompute, timestep, and others.
Step 2: Allocate storage for the initial data (y_n_gfs gridfunctions) on each grid.
Step 3: {step3desc}
Step 4: {step4desc}
"""
    if post_non_y_n_auxevol_mallocs:
        desc += "Step 4.a: Post-initial-data functions called after memory for non-y_n and auxevol gridfunctions is allocated."
    desc += f"""Step 5: MAIN SIMULATION LOOP
- Step 5.a: Functions to run prior to diagnostics. E.g., regridding.
- Step 5.b: Output diagnostics.
- Step 5.c: Prepare to step forward in time.
- Step 5.d: Step forward in time using Method of Lines with {MoL_method} algorithm, applying {boundary_conditions_desc} boundary conditions.
- Step 5.e: Finish up step in time.
Step 6: Free all allocated memory."""
    cfunc_type = "int"
    name = "main"
    params = "int argc, const char *argv[]"
    body = r"""  commondata_struct commondata; // commondata contains parameters common to all grids.
  griddata_struct *restrict griddata; // griddata contains data specific to an individual grid.

// Step 1.a: Initialize each CodeParameter in the commondata struct to its default value.
commondata_struct_set_to_default(&commondata);

// Step 1.b: Overwrite the default values with those from the parameter file.
//           Then overwrite the parameter file values with those provided via command line arguments.
cmdline_input_and_parfile_parser(&commondata, argc, argv);

// Step 1.c: Allocate memory for MAXNUMGRIDS griddata structs,
//           where each structure contains data specific to an individual grid.
griddata = (griddata_struct *restrict)malloc(sizeof(griddata_struct)*MAXNUMGRIDS);

// Step 1.d: Initialize each CodeParameter in griddata.params to its default value.
params_struct_set_to_default(&commondata, griddata);

// Step 1.e: Set up numerical grids, including parameters such as NUMGRIDS, xx[3], masks, Nxx, dxx, invdxx,
//           bcstruct, rfm_precompute, timestep, and others.
{
  // If this function is being called for the first time, initialize commondata.time, nn, t_0, and nn_0 to 0.
  const bool calling_for_first_time = true;
  numerical_grids_and_timestep(&commondata, griddata, calling_for_first_time);
}

for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
  // Step 2: Allocate storage for the initial data (y_n_gfs gridfunctions) on each grid.
  MoL_malloc_y_n_gfs(&commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
}
"""
    setup_initial_data_code = """Set up initial data.
initial_data(&commondata, griddata);
"""
    allocate_storage_code = """Allocate storage for non-y_n gridfunctions, needed for the Runge-Kutta-like timestepping.
for(int grid=0; grid<commondata.NUMGRIDS; grid++)
  MoL_malloc_non_y_n_gfs(&commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
"""
    step3code = setup_initial_data_code
    step4code = allocate_storage_code
    if set_initial_data_after_auxevol_malloc:
        step3code = allocate_storage_code
        step4code = setup_initial_data_code
    body += f"""
// Step 3: {step3code}

// Step 4: {step4code}
"""
    if post_non_y_n_auxevol_mallocs:
        body += "// Step 4.a: Functions called after memory for non-y_n and auxevol gridfunctions is allocated.\n"
        body += post_non_y_n_auxevol_mallocs
    body += """
// Step 5: MAIN SIMULATION LOOP
while(commondata.time < commondata.t_final) { // Main loop to progress forward in time.
  // Step 5.a: Main loop, part 1 (pre_diagnostics): Functions to run prior to diagnostics. E.g., regridding.
"""
    if pre_diagnostics:
        body += pre_diagnostics
    else:
        body += "// (nothing here; specify by setting pre_diagnostics string in register_CFunction_main_c().)\n"
    body += """
  // Step 5.b: Main loop, part 2: Output diagnostics
  diagnostics(&commondata, griddata);

  // Step 5.c: Main loop, part 3 (pre_MoL_step_forward_in_time): Prepare to step forward in time
"""
    if pre_MoL_step_forward_in_time:
        body += pre_MoL_step_forward_in_time
    else:
        body += "// (nothing here; specify by setting pre_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"
    body += f"""
  // Step 5.d: Main loop, part 4: Step forward in time using Method of Lines with {MoL_method} algorithm,
  //           applying {boundary_conditions_desc} boundary conditions.
  MoL_step_forward_in_time(&commondata, griddata);

  // Step 5.e: Main loop, part 5 (post_MoL_step_forward_in_time): Finish up step in time
"""
    if post_MoL_step_forward_in_time != "":
        body += post_MoL_step_forward_in_time
    else:
        body += "  // (nothing here; specify by setting post_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"
    body += r"""
} // End main loop to progress forward in time.

// Step 6: Free all allocated memory
{
  const bool free_non_y_n_gfs_and_core_griddata_pointers=true;
  griddata_free(&commondata, griddata, free_non_y_n_gfs_and_core_griddata_pointers);
}
return 0;
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
        clang_format_options=clang_format_options,
    )
