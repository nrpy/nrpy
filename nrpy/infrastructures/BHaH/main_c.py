"""
Generates the C main() function for all codes in the BHaH infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List, Tuple
import nrpy.c_function as cfc


def register_CFunction_main_c(
    MoL_method: str,
    initial_data_desc: str = "",
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
    boundary_conditions_desc: str = "",
    prefunc: str = "",
    pre_MoL_step_forward_in_time: str = "",
    post_MoL_step_forward_in_time: str = "",
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 0}",
) -> None:
    """
    Generates the C main() function for all codes in the BHaH infrastructure.

    :param MoL_method: Method of Lines algorithm used to step forward in time.
    :param initial_data_desc: Description for initial data, default is an empty string.
    :param enable_rfm_precompute: Enable rfm precomputation, default is False.
    :param enable_CurviBCs: Enable CurviBCs, default is False.
    :param boundary_conditions_desc: Description of the boundary conditions, default is an empty string.
    :param prefunc: String that appears before main(). DO NOT populate this, except when debugging, default is an empty string.
    :param pre_MoL_step_forward_in_time: Code for handling pre-right-hand-side operations, default is an empty string.
    :param post_MoL_step_forward_in_time: Code for handling post-right-hand-side operations, default is an empty string.
    :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 0}".
    """
    initial_data_desc += " "
    # Make sure all required C functions are registered
    missing_functions: List[Tuple[str, str]] = []
    for func_tuple in [
        ("params_struct_set_to_default", "CodeParameters.py"),
        ("commondata_struct_set_to_default", "CodeParameters.py"),
        ("cmdline_input_and_parfile_parser", "cmdline_input_and_parfiles.py"),
        (
            "numerical_grids_and_timestep_setup",
            "e.g., reference_metric.py or user defined",
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
    desc = """-={ main() function }=-
Step 1.a: Set each commondata CodeParameter to default.
Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
Step 1.c: Allocate NUMGRIDS griddata arrays, each containing data specific to an individual grid.
Step 1.d: Set each CodeParameter in griddata.params to default.
"""
    if enable_CurviBCs:
        desc += "Step 1.e: Set non-parfile parameters related to numerical grid, then set up numerical grids and CFL-limited timestep.\n"
    if enable_rfm_precompute:
        desc += "Step 1.f: Set up boundary condition struct (bcstruct)\n"
    desc += f"""Step 2: Declare and allocate memory for gridfunctions
Step 3: Finalize initialization: set up {initial_data_desc}initial data, etc.
Step 4: MAIN SIMULATION LOOP
- Step 4.a: Output diagnostics
- Step 4.b: Prepare to step forward in time
- Step 4.c: Step forward in time using Method of Lines with {MoL_method} algorithm, applying {boundary_conditions_desc} boundary conditions.
- Step 4.d: Finish up step in time
Step 5: Free all allocated memory"""
    c_type = "int"
    name = "main"
    params = "int argc, const char *argv[]"
    body = r"""  commondata_struct commondata; // commondata contains parameters common to all grids.
  griddata_struct *restrict griddata; // griddata contains data specific to an individual grid.

// Step 1.a: Set each commondata CodeParameter to default.
commondata_struct_set_to_default(&commondata);
// Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
cmdline_input_and_parfile_parser(&commondata, argc, argv);

// Step 1.c: Allocate NUMGRIDS griddata arrays, each containing data specific to an individual grid.
griddata = (griddata_struct *restrict)malloc(sizeof(griddata_struct)*commondata.NUMGRIDS);
// Step 1.d: Set each CodeParameter in griddata.params to default.
params_struct_set_to_default(&commondata, griddata);
for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
  // Step 1.e: Set non-parfile parameters related to numerical grid, then set up numerical grids and CFL-limited timestep.
  numerical_grids_and_timestep_setup(&commondata, &griddata[grid].params, &griddata[grid]);
"""
    if enable_CurviBCs:
        body += r"""
  // Step 1.f: Set up boundary condition struct (bcstruct)
  bcstruct_set_up(&commondata, &griddata[grid].params, griddata[grid].xx, &griddata[grid].bcstruct);
"""
    if enable_rfm_precompute:
        body += r"""
  // Step 1.g: Allocate memory for and define reference-metric precomputation lookup arrays
  rfm_precompute_malloc(&commondata, &griddata[grid].params, &griddata[grid].rfmstruct);
  rfm_precompute_defines(&commondata, &griddata[grid].params, &griddata[grid].rfmstruct, griddata[grid].xx);
"""
    body += r"""
  // Step 2: Declare and allocate memory for gridfunctions
  MoL_malloc_y_n_gfs(&commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
  MoL_malloc_non_y_n_gfs(&commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
}
// Step 3: Finalize initialization: set up initial data, etc.
initial_data(&commondata, griddata);

// Step 4: MAIN SIMULATION LOOP
while(commondata.time < commondata.t_final) { // Main loop to progress forward in time.
  // Step 4.a: Main loop, part 1: Output diagnostics
  diagnostics(&commondata, griddata);

  // Step 4.b: Main loop, part 2 (pre_MoL_step_forward_in_time): Prepare to step forward in time
"""
    if pre_MoL_step_forward_in_time != "":
        body += pre_MoL_step_forward_in_time
    else:
        body += "  // (nothing here; specify by setting pre_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"
    body += f"""
  // Step 4.c: Main loop, part 3: Step forward in time using Method of Lines with {MoL_method} algorithm,
  //           applying {boundary_conditions_desc} boundary conditions.
  MoL_step_forward_in_time(&commondata, griddata);

  // Step 4.d: Main loop, part 4 (post_MoL_step_forward_in_time): Finish up step in time
"""
    if post_MoL_step_forward_in_time != "":
        body += post_MoL_step_forward_in_time
    else:
        body += "  // (nothing here; specify by setting post_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"
    body += r"""
} // End main loop to progress forward in time.

// Step 5: Free all allocated memory
for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
  MoL_free_memory_y_n_gfs(&griddata[grid].gridfuncs);
  MoL_free_memory_non_y_n_gfs(&griddata[grid].gridfuncs);"""
    if enable_rfm_precompute:
        body += r"""
  rfm_precompute_free(&griddata[grid].rfmstruct);"""
    if enable_CurviBCs:
        body += r"""
  free(griddata[grid].bcstruct.inner_bc_array);
  for(int ng=0;ng<NGHOSTS*3;ng++) free(griddata[grid].bcstruct.pure_outer_bc_array[ng]);
"""
    body += r"""
  for(int i=0;i<3;i++) free(griddata[grid].xx[i]);
}
free(griddata);
return 0;
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
        clang_format_options=clang_format_options,
    )
