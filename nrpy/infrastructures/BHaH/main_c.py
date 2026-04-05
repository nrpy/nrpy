# nrpy/infrastructures/BHaH/main_c.py
"""
Generates the main C `main()` function for simulations within the BHaH infrastructure.

This module acts as the final assembler for a BHaH-based executable, creating
the top-level `main()` function that orchestrates the entire simulation.

The `register_CFunction_main_c` function constructs a standard simulation
lifecycle, which includes:
- Parsing command-line arguments and parameter files.
- Setting up numerical grids and the time step.
- Allocating memory for all required gridfunctions.
- Setting the initial data.
- Executing the main time-evolution loop.
- Freeing memory upon completion.

A key feature is its configurability. The main loop is designed with specific
hooks (e.g., `pre_diagnostics`, `post_MoL_step_forward_in_time`) that allow
custom C code to be injected at various stages, enabling flexible simulation
logic without altering the core structure. The module also validates that all
necessary C functions from other components have been registered before
generating the `main` function that calls them.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List, Tuple

import nrpy.c_function as cfc
import nrpy.params as par


def _check_required_functions() -> None:
    """
    Check for all required C functions for BHaH main().

    :raises ValueError: If any required C functions are not registered.
    """
    required_functions: List[Tuple[str, str]] = [
        ("params_struct_set_to_default", "CodeParameters.py"),
        ("commondata_struct_set_to_default", "CodeParameters.py"),
        ("cmdline_input_and_parfile_parser", "cmdline_input_and_parfiles.py"),
        (
            "numerical_grids_and_timestep",
            "e.g., numerical_grids_and_timestep.py or user defined",
        ),
        ("MoL_free_intermediate_stage_gfs", "MoLtimestepping/"),
        ("MoL_malloc_intermediate_stage_gfs", "MoLtimestepping/"),
        ("MoL_step_forward_in_time", "MoLtimestepping/"),
        ("diagnostics", "log10_L2norm_gf.py"),
        ("initial_data", "initial_data.py"),
    ]

    missing_functions = [
        func_tuple
        for func_tuple in required_functions
        if func_tuple[0] not in cfc.CFunction_dict
    ]

    if missing_functions:
        error_msg = "Error: These functions are required for all BHaH main() functions, and are not registered.\n"
        for name, origin in missing_functions:
            error_msg += f'  {name}, registered by function within "{origin}"\n'
        raise ValueError(error_msg)


def _generate_main_desc(
    MoL_method: str,
    initial_data_desc: str,
    set_initial_data_after_auxevol_malloc: bool,
    boundary_conditions_desc: str,
    post_non_y_n_auxevol_mallocs: str,
) -> str:
    """
    Generate the descriptive comment block for the main() C function.

    :param MoL_method: Method of Lines algorithm.
    :param initial_data_desc: Description for initial data.
    :param set_initial_data_after_auxevol_malloc: Flag to set initial data after auxevol malloc.
    :param boundary_conditions_desc: Description of the boundary conditions.
    :param post_non_y_n_auxevol_mallocs: String of post-malloc function calls.
    :return: The formatted description string for the C function.
    """
    # Define descriptions for steps that can be reordered.
    allocate_auxevol_desc = "Allocate storage for non-y_n gridfunctions needed for the Runge-Kutta-like time-stepping."
    set_initial_data_desc = f"Set {initial_data_desc} initial data."

    # Determine the order of Step 3 and 4 descriptions based on the input flag.
    if set_initial_data_after_auxevol_malloc:
        step3_desc = allocate_auxevol_desc
        step4_desc = set_initial_data_desc
    else:
        step3_desc = set_initial_data_desc
        step4_desc = allocate_auxevol_desc

    # Build the main description string.
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
Step 3: {step3_desc}
Step 4: {step4_desc}
"""
    if post_non_y_n_auxevol_mallocs:
        # This string is appended without a preceding newline to match original output.
        desc += "Step 4.a: Post-initial-data functions called after memory for non-y_n and auxevol gridfunctions is allocated."

    desc += f"""Step 5: MAIN SIMULATION LOOP
- Step 5.a: Functions to run prior to diagnostics. E.g., regridding.
- Step 5.b: Output diagnostics.
- Step 5.c: Prepare to step forward in time.
- Step 5.d: Step forward in time using Method of Lines with {MoL_method} algorithm, applying {boundary_conditions_desc} boundary conditions.
- Step 5.e: Finish up step in time.
Step 6: Free all allocated memory."""
    return desc


def _generate_main_body(
    MoL_method: str,
    set_initial_data_after_auxevol_malloc: bool,
    boundary_conditions_desc: str,
    post_non_y_n_auxevol_mallocs: str,
    pre_diagnostics: str,
    pre_MoL_step_forward_in_time: str,
    post_MoL_step_forward_in_time: str,
) -> str:
    """
    Generate the C code for the body of the main() function.

    This function constructs the main C function body by appending logical C code
    blocks to a list, exactly replicating the original's proven-correct assembly
    logic. This approach guarantees byte-for-byte compatibility with the trusted
    output while improving readability over the original.

    :param MoL_method: Method of Lines algorithm.
    :param set_initial_data_after_auxevol_malloc: Flag to set initial data after auxevol malloc.
    :param boundary_conditions_desc: Description of the boundary conditions.
    :param post_non_y_n_auxevol_mallocs: String of post-malloc function calls.
    :param pre_diagnostics: String of pre-diagnostics function calls.
    :param pre_MoL_step_forward_in_time: String of pre-MoL step function calls.
    :param post_MoL_step_forward_in_time: String of post-MoL step function calls.
    :return: The formatted C code for the function body.
    """
    # Prepare frequently used strings and logic based on parallelization settings.
    parallelization = par.parval_from_str("parallelization")
    is_cuda = parallelization == "cuda"
    compute_griddata = "griddata_device" if is_cuda else "griddata"

    # The body of the main() function is built as a list of C code strings.
    body_parts = []

    # Step 0: Variable declarations
    declarations_c_code = "commondata_struct commondata; // commondata contains parameters common to all grids.\n"
    declarations_c_code += f"griddata_struct *{compute_griddata}; // griddata contains data specific to an individual grid.\n"
    if is_cuda:
        declarations_c_code += r"""griddata_struct *griddata_host; // stores only the host data needed for diagnostics
#include "BHaH_CUDA_global_init.h"
"""
    body_parts.append(declarations_c_code)

    # Step 1: Initialization
    griddata_malloc_code = f"{compute_griddata} = (griddata_struct *)malloc(sizeof(griddata_struct) * MAXNUMGRIDS);"
    if is_cuda:
        griddata_malloc_code += "\ngriddata_host = (griddata_struct *)malloc(sizeof(griddata_struct) * MAXNUMGRIDS);"
    numerical_grids_args = (
        f"&commondata, {compute_griddata}"
        + (", griddata_host" if is_cuda else "")
        + ", calling_for_first_time"
    )

    step1_c_code = f"""
// Step 1.a: Initialize each CodeParameter in the commondata struct to its default value.
commondata_struct_set_to_default(&commondata);

// Step 1.b: Overwrite the default values with those from the parameter file.
//           Then overwrite the parameter file values with those provided via command line arguments.
cmdline_input_and_parfile_parser(&commondata, argc, argv);

// Step 1.c: Allocate memory for MAXNUMGRIDS griddata structs,
//           where each structure contains data specific to an individual grid.
{griddata_malloc_code}

// Step 1.d: Initialize each CodeParameter in {compute_griddata}.params to its default value.
params_struct_set_to_default(&commondata, {compute_griddata});

// Step 1.e: Set up numerical grids, including parameters such as NUMGRIDS, xx[3], masks, Nxx, dxx, invdxx,
//           bcstruct, rfm_precompute, timestep, and others.
{{
  IFCUDARUN(for (int grid = 0; grid < MAXNUMGRIDS; grid++) griddata_device[grid].params.is_host = false;);
  // If this function is being called for the first time, initialize commondata.time, nn, t_0, and nn_0 to 0.
  const bool calling_for_first_time = true;
  numerical_grids_and_timestep({numerical_grids_args});
}} // END setup of numerical & temporal grids.
"""
    body_parts.append(step1_c_code)

    allocator_macro = (
        "BHAH_MALLOC_DEVICE" if parallelization == "cuda" else "BHAH_MALLOC"
    )
    # Step 2: Allocate memory for evolved gridfunctions (y_n_gfs).
    step2_c_code = f"""
// Step 2: Allocate storage for the initial data (y_n_gfs gridfunctions) on each grid.
for(int grid=0; grid<commondata.NUMGRIDS; grid++) {{
  const int Nxx_plus_2NGHOSTS_tot = ({compute_griddata}[grid].params.Nxx_plus_2NGHOSTS0 * //
                                     {compute_griddata}[grid].params.Nxx_plus_2NGHOSTS1 * //
                                     {compute_griddata}[grid].params.Nxx_plus_2NGHOSTS2);
  {allocator_macro}({compute_griddata}[grid].gridfuncs.y_n_gfs, sizeof(REAL) * Nxx_plus_2NGHOSTS_tot * NUM_EVOL_GFS);
  if (NUM_AUXEVOL_GFS > 0) {{
    {allocator_macro}({compute_griddata}[grid].gridfuncs.auxevol_gfs, sizeof(REAL) * Nxx_plus_2NGHOSTS_tot * NUM_AUXEVOL_GFS);
    IFCUDARUN(BHAH_MALLOC(griddata_host[grid].gridfuncs.auxevol_gfs, sizeof(REAL) * Nxx_plus_2NGHOSTS_tot * NUM_AUXEVOL_GFS););
  }} // END IF NUM_AUXEVOL_GFS > 0

  // On GPU, separately allocate y_n_gfs on the host, for diagnostics purposes.
  IFCUDARUN(BHAH_MALLOC_PINNED(griddata_host[grid].gridfuncs.y_n_gfs, sizeof(REAL) * Nxx_plus_2NGHOSTS_tot * NUM_EVOL_GFS););
}} // END LOOP over grids
"""
    body_parts.append(step2_c_code)

    # Steps 3 & 4: Set initial data and allocate memory for auxiliary gridfunctions.
    initial_data_call = f"initial_data(&commondata, {f'griddata_host, {compute_griddata}' if is_cuda else compute_griddata});"
    setup_initial_data_code = f"Set up initial data.\n{initial_data_call}\n"
    allocate_storage_code = f"""Allocate storage for non-y_n gridfunctions, needed for the Runge-Kutta-like timestepping.
for(int grid=0; grid<commondata.NUMGRIDS; grid++)
  MoL_malloc_intermediate_stage_gfs(&commondata, &{compute_griddata}[grid].params, &{compute_griddata}[grid].gridfuncs);\n"""

    if set_initial_data_after_auxevol_malloc:
        step3_code, step4_code = (allocate_storage_code, setup_initial_data_code)
    else:
        step3_code, step4_code = (setup_initial_data_code, allocate_storage_code)
    body_parts.append(f"\n// Step 3: {step3_code}\n// Step 4: {step4_code}")

    if post_non_y_n_auxevol_mallocs:
        body_parts.append(
            "// Step 4.a: Functions called after memory for non-y_n and auxevol gridfunctions is allocated.\n"
        )
        body_parts.append(post_non_y_n_auxevol_mallocs)

    # Step 5: Main simulation loop.
    diagnostics_call_args = f"&commondata, {f'{compute_griddata}, griddata_host' if is_cuda else compute_griddata}"
    body_parts.append("""
// Step 5: MAIN SIMULATION LOOP
while(commondata.time < commondata.t_final) { // Main loop to progress forward in time.
  // Step 5.a: Main loop, part 1 (pre_diagnostics): Functions to run prior to diagnostics. E.g., regridding.
""")
    body_parts.append(
        pre_diagnostics
        or "// (nothing here; specify by setting pre_diagnostics string in register_CFunction_main_c().)\n"
    )

    body_parts.append(f"""
  // Step 5.b: Main loop, part 2: Output diagnostics
  diagnostics({diagnostics_call_args});

  // Step 5.c: Main loop, part 3 (pre_MoL_step_forward_in_time): Prepare to step forward in time
""")
    body_parts.append(
        pre_MoL_step_forward_in_time
        or "// (nothing here; specify by setting pre_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"
    )

    body_parts.append(f"""
  // Step 5.d: Main loop, part 4: Step forward in time using Method of Lines with {MoL_method} algorithm,
  //           applying {boundary_conditions_desc} boundary conditions.
  MoL_step_forward_in_time(&commondata, {compute_griddata});

  // Step 5.e: Main loop, part 5 (post_MoL_step_forward_in_time): Finish up step in time
""")
    body_parts.append(
        post_MoL_step_forward_in_time
        or "  // (nothing here; specify by setting post_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"
    )

    # Step 6: Free all allocated memory
    device_sync = "BHAH_DEVICE_SYNC();" if is_cuda else ""
    if not is_cuda:
        free_memory_code = rf"""
  const bool free_non_y_n_gfs_and_core_griddata_pointers=true;
  griddata_free(&commondata, {compute_griddata}, free_non_y_n_gfs_and_core_griddata_pointers);
}}
        """
    else:
        free_memory_code = rf"""
  const bool free_non_y_n_gfs_and_core_griddata_pointers=true;
  griddata_free_device(&commondata, {compute_griddata}, free_non_y_n_gfs_and_core_griddata_pointers);
  griddata_free(&commondata, griddata_host, free_non_y_n_gfs_and_core_griddata_pointers);
}}
for (int i = 0; i < NUM_STREAMS; ++i) {{
  cudaStreamDestroy(streams[i]);
}}
BHAH_DEVICE_SYNC();
cudaDeviceReset();
"""
    body_parts.append(f"""
}} // End main loop to progress forward in time.
{device_sync}
// Step 6: Free all allocated memory
{{{free_memory_code}""")
    body_parts.append(r"""return 0;
""")

    # Construct the final body string and perform necessary replacements.
    body = "".join(body_parts)
    return body


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

    Doctest:
    >>> import nrpy.c_function as cfc
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.params as par
    >>> from nrpy.infrastructures import BHaH
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> set_of_coordsys = {"Cartesian"}
    >>> project_name = "main_test"
    >>> Nxx_dict = {"Cartesian" : [1, 1, 1]}
    >>> for parallelization in supported_Parallelizations:
    ...    par.glb_extras_dict.clear()
    ...    cfc.CFunction_dict.clear()
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    BHaH.CodeParameters.register_CFunctions_params_commondata_struct_set_to_default()
    ...    BHaH.cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(project_name)
    ...    _ = BHaH.diagnostics.diagnostics.register_all_diagnostics(project_dir="/tmp", set_of_CoordSystems=["Spherical"], default_diagnostics_out_every=100, enable_nearest_diagnostics=True, enable_interp_diagnostics=False, enable_volume_integration_diagnostics=False, enable_free_auxevol=False)
    ...    _ = BHaH.numerical_grids_and_timestep.register_CFunctions(set_of_coordsys, [5], Nxx_dict)
    ...    _ = BHaH.MoLtimestepping.register_all.register_CFunctions()
    ...    _ = BHaH.wave_equation.initial_data_exact_soln.register_CFunction_initial_data()
    ...    register_CFunction_main_c("RK4")
    ...    generated_str = cfc.CFunction_dict["main"].full_function
    ...    validation_desc = f"_{parallelization}".replace(" ", "_")
    ...    validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    Setting up reference_metric[Cartesian]...
    """
    # Step 1: Check for all required functions.
    _check_required_functions()

    # Step 2: Set up parameters and include files based on parallelization.
    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if parallelization == "cuda":
        includes.append("BHaH_global_device_defines.h")

    # Step 3: Generate the description and body of the main function using helpers.
    desc = _generate_main_desc(
        MoL_method,
        initial_data_desc,
        set_initial_data_after_auxevol_malloc,
        boundary_conditions_desc,
        post_non_y_n_auxevol_mallocs,
    )
    body = _generate_main_body(
        MoL_method,
        set_initial_data_after_auxevol_malloc,
        boundary_conditions_desc,
        post_non_y_n_auxevol_mallocs,
        pre_diagnostics,
        pre_MoL_step_forward_in_time,
        post_MoL_step_forward_in_time,
    )

    # Step 4: Register the CFunction.
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type="int",
        name="main",
        params="int argc, const char *argv[]",
        body=body,
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
