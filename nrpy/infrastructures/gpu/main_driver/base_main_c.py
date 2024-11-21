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
        self.MoL_method = MoL_method
        self.initial_data_desc = initial_data_desc
        self.set_initial_data_after_auxevol_malloc = (
            set_initial_data_after_auxevol_malloc
        )
        self.boundary_conditions_desc = boundary_conditions_desc
        self.prefunc = prefunc
        self.post_non_y_n_auxevol_mallocs = post_non_y_n_auxevol_mallocs
        self.pre_diagnostics = pre_diagnostics
        self.pre_MoL_step_forward_in_time = pre_MoL_step_forward_in_time
        self.post_MoL_step_forward_in_time = post_MoL_step_forward_in_time
        self.clang_format_options = clang_format_options

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
                error_msg += f'  {func_tuple[0]}, registered by function within "{func_tuple[1]}"\n'
            raise ValueError(error_msg)

        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

        set_initial_data_desc = f"Set {initial_data_desc} initial data."
        allocate_auxevol_desc = "Allocate storage for non-y_n gridfunctions needed for the Runge-Kutta-like time-stepping."

        step3desc = set_initial_data_desc
        step4desc = allocate_auxevol_desc
        if set_initial_data_after_auxevol_malloc:
            step3desc = allocate_auxevol_desc
            step4desc = set_initial_data_desc

        self.desc = f"""-={{ main() function }}=-
Step 1.a: Initialize each CodeParameter in the commondata struc to its default value.
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
        if self.post_non_y_n_auxevol_mallocs:
            self.desc += "Step 4.a: Post-initial-data functions called after memory for non-y_n and auxevol gridfunctions is allocated."
        self.desc += f"""Step 5: MAIN SIMULATION LOOP
- Step 5.a: Functions to run prior to diagnostics. E.g., regridding.
- Step 5.b: Output diagnostics.
- Step 5.c: Prepare to step forward in time.
- Step 5.d: Step forward in time using Method of Lines with {self.MoL_method} algorithm, applying {self.boundary_conditions_desc} boundary conditions.
- Step 5.e: Finish up step in time.
Step 6: Free all allocated memory."""
        self.cfunc_type = "int"
        self.name = "main"
        self.params = "int argc, const char *argv[]"
        self.body = ""

    def register(self) -> None:
        """Register CFunction."""
        cfc.register_CFunction(
            prefunc=self.prefunc,
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            body=self.body,
            clang_format_options=self.clang_format_options,
        )
