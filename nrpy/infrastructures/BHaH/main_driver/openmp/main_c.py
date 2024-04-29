"""
Generate the C main() function for all codes in the BHaH infrastructure using OpenMP.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot** com
        Samuel D. Tootle
        sdtootle **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.main_driver.base_main_c as base_main

class register_CFunction_main_c(base_main.base_register_CFunction_main_c):
    """
    Generate the "generic" C main() function for all simulation codes in the BHaH infrastructure.

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
    def __init__ (
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
        super().__init__(
            MoL_method,
            initial_data_desc=initial_data_desc,
            boundary_conditions_desc=boundary_conditions_desc,
            prefunc=prefunc,
            initialize_constant_auxevol=initialize_constant_auxevol,
            pre_MoL_step_forward_in_time=pre_MoL_step_forward_in_time,
            post_MoL_step_forward_in_time=post_MoL_step_forward_in_time,
            clang_format_options=clang_format_options
        )
    
        cfc.register_CFunction(
            includes=self.includes,
            prefunc=self.prefunc,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            body=self.body,
            clang_format_options=clang_format_options,
        )
