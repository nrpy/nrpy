# nrpy/infrastructures/BHaH/MoLtimestepping/MoL_rk_substep.py
"""
Functions and classes for constructing individual Runge-Kutta substeps and the dictionary that stores them.

Authors: Zachariah B. Etienne (lead maintainer)
         zachetie **at** gmail **dot* com
         Samuel Tootle (GPU support in NRPy2)
         Brandon Clark (original, NRPy1 version)
"""

from typing import Dict, List, Union

import sympy as sp

import nrpy.c_function as cfc
from nrpy.c_codegen import c_codegen
from nrpy.helpers.generic import superfast_uniq
from nrpy.helpers.gpu.utilities import generate_kernel_and_launch_code

##############################################################################
# Parallelization check and set:

supported_parallelization = {"cuda", "openmp"}


def check_supported_parallelization(
    function_name: str,
    parallelization: str,
) -> None:
    """
    Check if the given parallelization is supported.

    :param function_name: Name of the function where the check is performed.
    :param parallelization: Parameter to specify parallelization (openmp or cuda).
    :raises ValueError: If the parallelization is not supported.
    """
    if parallelization not in supported_parallelization:
        raise ValueError(
            f"ERROR ({function_name}): {parallelization} is not supported."
        )


##############################################################################
# The RKFunction class and a dictionary storing substep logic:


class RKFunction:
    """
    A class to represent Runge-Kutte (RK) substep functions in C/C++.

    :param RK_lhs_list: List of LHS expressions for RK substep.
    :param RK_rhs_list: List of RHS expressions for RK substep.
    :param enable_intrinsics: A flag to specify if hardware instructions should be used.
    :param cfunc_type: decorators and return type for the RK substep function
    :param rk_step: current step (> 0).  Default (None) assumes Euler step
    :param rational_const_alias: Overload const specifier for Rational definitions
    :param parallelization: Parameter to specify parallelization (openmp or cuda).
    """

    def __init__(
        self,
        RK_lhs_list: List[sp.Basic],
        RK_rhs_list: List[sp.Basic],
        enable_intrinsics: bool = False,
        cfunc_type: str = "static void",
        rk_step: Union[int, None] = None,
        rational_const_alias: str = "const",
        parallelization: str = "openmp",
    ) -> None:
        check_supported_parallelization("RKFunction", parallelization)
        self.rk_step = rk_step
        self.enable_intrinsics = enable_intrinsics
        self.intrinsics_str = "CUDA" if parallelization == "cuda" else "SIMD"
        self.RK_rhs_list = RK_rhs_list
        self.RK_lhs_list = RK_lhs_list
        self.name: str = ""
        self.params: str = "params_struct *restrict params, "
        self.body: str = ""

        self.cfunc_type = cfunc_type
        self.rational_const_alias = rational_const_alias

        self.CFunction: cfc.CFunction
        self.desc = f"Runge-Kutta function for substep {self.rk_step}."

        # Save variables that appear in function arguments to be used
        # to generate a function call
        self.param_vars: List[str] = []
        self.includes: List[str] = []

        self.RK_lhs_str_list = [
            (
                f"const REAL_{self.intrinsics_str}_ARRAY __rk_exp_{i}"
                if self.enable_intrinsics
                else f"{str(el).replace('gfsL', 'gfs[i]')}"
            )
            for i, el in enumerate(self.RK_lhs_list)
        ]

        self.loop_body = c_codegen(
            self.RK_rhs_list,
            self.RK_lhs_str_list,
            include_braces=False,
            verbose=False,
            enable_simd=self.enable_intrinsics,
            enable_cse_preprocess=True,
            rational_const_alias=self.rational_const_alias,
        ).replace("SIMD", self.intrinsics_str)
        # Give rationals a better name
        self.loop_body = self.loop_body.replace("_Rational", "RK_Rational")

        if enable_intrinsics:
            for i, el in enumerate(self.RK_lhs_list):
                self.loop_body += f"Write{self.intrinsics_str}(&{str(el).replace('gfsL', 'gfs[i]')}, __rk_exp_{i});\n"

        self.name = f"rk_substep_{self.rk_step}"

        # Populate build and populate self.CFunction
        kernel_body: str = ""
        self.kernel_params = {}
        self.kernel_params["params"] = "params_struct *restrict"

        if parallelization == "cuda":
            self.body += (
                "\n".join(
                    f"const int Nxx_plus_2NGHOSTS{X} = params->Nxx_plus_2NGHOSTS{X};"
                    for X in ["0", "1", "2"]
                )
                + "\n"
            )
            self.body += "MAYBE_UNUSED const int Ntot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;\n\n"

        kernel_body += "LOOP_ALL_GFS_GPS(i) {\n"

        read_list = [
            read
            for el in self.RK_rhs_list
            for read in list(sp.ordered(el.free_symbols))
        ]
        read_list_unique = superfast_uniq(read_list)

        for el in read_list_unique:
            if str(el) != "commondata->dt":
                gfs_el = str(el).replace("gfsL", "gfs[i]")
                key = gfs_el[:-3]
                if self.enable_intrinsics and parallelization == "openmp":
                    self.kernel_params[key] = "REAL *restrict"
                    self.params += f"REAL *restrict {key},"
                    kernel_body += f"const REAL_SIMD_ARRAY {el} = Read{self.intrinsics_str}(&{gfs_el});\n"
                else:
                    self.kernel_params[key] = "REAL *restrict"
                    self.params += f"REAL *restrict {key},"
                    kernel_body += f"const REAL {el} = {gfs_el};\n"
        for el in self.RK_lhs_list:
            lhs_var = str(el).replace("_gfsL", "_gfs")
            if lhs_var not in self.params:
                self.kernel_params[lhs_var] = "REAL *restrict"
                self.params += f"REAL *restrict {lhs_var},"
        self.params += f"const REAL dt"
        self.kernel_params["dt"] = "const REAL"

        kernel_body += self.loop_body.replace("commondata->dt", "dt") + "\n}\n"
        self.kernel_params_lst = [f"{v} {k}" for k, v in self.kernel_params.items()]

        comments = f"Compute RK substep {self.rk_step}."
        # Prepare the argument dicts
        arg_dict_cuda = self.kernel_params.copy()
        arg_dict_cuda.pop("params")
        if self.enable_intrinsics and parallelization == "openmp":
            kernel_body = kernel_body.replace("dt", "DT")
            kernel_body = "const REAL_SIMD_ARRAY DT = ConstSIMD(dt);\n" + kernel_body
        prefunc, new_body = generate_kernel_and_launch_code(
            self.name,
            kernel_body,
            arg_dict_cuda,
            self.kernel_params,
            parallelization=parallelization,
            comments=comments,
            launch_dict={
                "blocks_per_grid": ["(Ntot + threads_in_x_dir - 1) / threads_in_x_dir"],
                "threads_per_block": ["32"],
                "stream": "params->grid_idx % NUM_STREAMS",
            },
        )
        self.body += new_body

        self.name += "__launcher"
        self.CFunction = cfc.CFunction(
            prefunc=prefunc,
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            body=self.body,
        )

    def c_function_call(self) -> str:
        """
        Generate the C function call for a given RK substep.

        :return: The C function call as a string.
        """
        c_function_call: str = self.name + "("
        for p in self.kernel_params:
            if p != "dt":
                c_function_call += f"{p}, "

        c_function_call += "commondata->dt);\n"
        return c_function_call


MoL_Functions_dict: Dict[str, RKFunction] = {}


def construct_RK_functions_prefunc() -> str:
    """
    Construct the prefunc (CFunction) strings for all RK functions stored in MoL_Functions_dict.

    :return: The concatenated prefunc (CFunction) strings as a single string.
    :raises ValueError: If the MoL_Functions_dict is empty
    """
    if len(MoL_Functions_dict.values()) == 0:
        raise ValueError("ERROR: MoL_Functions_dict is empty")

    prefunc = ""
    for fd_func in MoL_Functions_dict.values():
        prefunc += fd_func.CFunction.full_function + "\n\n"
    return prefunc


def single_RK_substep_input_symbolic(
    substep_time_offset_dt: Union[sp.Basic, int, str],
    rhs_str: str,
    rhs_input_expr: sp.Basic,
    rhs_output_expr: sp.Basic,
    RK_lhs_list: Union[sp.Basic, List[sp.Basic]],
    RK_rhs_list: Union[sp.Basic, List[sp.Basic]],
    post_rhs_list: Union[str, List[str]],
    post_rhs_output_list: Union[sp.Basic, List[sp.Basic]],
    additional_comments: str = "",
    enable_intrinsics: bool = False,
    gf_aliases: str = "",
    post_post_rhs_string: str = "",
    rational_const_alias: str = "const",
    rk_step: Union[int, None] = None,
    parallelization: str = "openmp",
) -> str:
    """
    Generate C code for a given Runge-Kutta substep.

    :param substep_time_offset_dt: Time offset for the RK substep.
    :param rhs_str: Right-hand side string of the C code.
    :param rhs_input_expr: Input expression for the RHS.
    :param rhs_output_expr: Output expression for the RHS.
    :param RK_lhs_list: List of LHS expressions for RK.
    :param RK_rhs_list: List of RHS expressions for RK.
    :param post_rhs_list: List of post-RHS expressions.
    :param post_rhs_output_list: List of outputs for post-RHS expressions.
    :param post_post_rhs_string: String to be used after the post-RHS phase.
    :param additional_comments: additional comments to append to auto-generated comment block.
    :param enable_intrinsics: A flag to specify if hardware instructions should be used.
    :param gf_aliases: Additional aliases for grid functions.
    :param rational_const_alias: Provide additional/alternative alias to const for rational definitions
    :param rk_step: Optional integer representing the current RK step.
    :param parallelization: Parameter to specify parallelization (openmp or cuda).

    :return: A string containing the generated C code.
    :raises ValueError: If substep_time_offset_dt cannot be extracted from the Butcher table.
    """
    check_supported_parallelization("single_RK_substep_input_symbolic", parallelization)
    # Ensure lists are lists
    if not isinstance(RK_lhs_list, list):
        RK_lhs_list = [RK_lhs_list]
    if not isinstance(RK_rhs_list, list):
        RK_rhs_list = [RK_rhs_list]
    if not isinstance(post_rhs_list, list):
        post_rhs_list = [post_rhs_list]
    if not isinstance(post_rhs_output_list, list):
        post_rhs_output_list = [post_rhs_output_list]

    comment_block = (
        "// ***Euler timestepping only requires one RHS evaluation***"
        if rk_step is None
        else f"// -={{ START k{rk_step} substep }}=-"
    )
    comment_block += additional_comments
    body = f"{comment_block}\n"

    if isinstance(substep_time_offset_dt, (int, sp.Rational, sp.Mul, sp.Float)):
        substep_time_offset_str = f"{float(substep_time_offset_dt):.17e}"
    else:
        raise ValueError(
            f"Could not extract substep_time_offset_dt={substep_time_offset_dt} from Butcher table"
        )

    body += "for(int grid=0; grid<commondata->NUMGRIDS; grid++) {\n"
    body += (
        f"commondata->time = time_start + {substep_time_offset_str} * commondata->dt;\n"
    )
    body += gf_aliases

    # (1) RHS
    updated_rhs_str = (
        str(rhs_str)
        .replace("RK_INPUT_GFS", str(rhs_input_expr).replace("gfsL", "gfs"))
        .replace("RK_OUTPUT_GFS", str(rhs_output_expr).replace("gfsL", "gfs"))
    )
    body += updated_rhs_str + "\n"

    # (2) RK update
    RK_key = f"RK_STEP{rk_step}"
    MoL_Functions_dict[RK_key] = RKFunction(
        RK_lhs_list,
        RK_rhs_list,
        rk_step=rk_step,
        enable_intrinsics=enable_intrinsics,
        rational_const_alias=rational_const_alias,
        parallelization=parallelization,
    )
    body += MoL_Functions_dict[RK_key].c_function_call()

    # (3) Post-RHS
    for post_rhs, post_rhs_output in zip(post_rhs_list, post_rhs_output_list):
        body += post_rhs.replace(
            "RK_OUTPUT_GFS", str(post_rhs_output).replace("gfsL", "gfs")
        )

    body += "}\n"

    for post_rhs, post_rhs_output in zip(post_rhs_list, post_rhs_output_list):
        body += post_post_rhs_string.replace(
            "RK_OUTPUT_GFS", str(post_rhs_output).replace("gfsL", "gfs")
        )

    return body


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
