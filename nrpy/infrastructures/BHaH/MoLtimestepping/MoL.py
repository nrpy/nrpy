"""
Module for producing C codes related to MoL timestepping within the BHaH infrastructure.
This includes implementation details and functions for allocating and deallocating the necessary memory.

Authors: Brandon Clark
         Zachariah B. Etienne (maintainer)
         zachetie **at** gmail **dot* com
"""

import os  # Standard Python module for multiplatform OS-level functions
import warnings
from typing import Dict, List, Tuple, Union

import sympy as sp  # Import SymPy, a computer algebra system written entirely in Python

import nrpy.c_function as cfc
import nrpy.helpers.gpu.gpu_kernel as gputils
import nrpy.params as par  # NRPy+: Parameter interface
from nrpy.c_codegen import c_codegen
from nrpy.helpers.generic import superfast_uniq
from nrpy.infrastructures.BHaH import BHaH_defines_h, griddata_commondata
from nrpy.infrastructures.BHaH.MoLtimestepping.RK_Butcher_Table_Dictionary import (
    generate_Butcher_tables,
)

# fmt: off
_ = par.CodeParameter("int", __name__, "nn_0", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("int", __name__, "nn", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("REAL", __name__, "CFL_FACTOR", 0.5, commondata=True)
_ = par.CodeParameter("REAL", __name__, "dt", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("REAL", __name__, "t_0", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("REAL", __name__, "time", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("REAL", __name__, "t_final", 10.0, commondata=True)
# fmt: on

supported_parallelization = {"cuda", "openmp"}
def check_supported_parallelization(function_name: str, parallelization: str,
) -> None:
    """
    Check if the given parallelization is supported.

    :param parallelization: Parameter to specify parallelization (openmp or cuda).
    :raises ValueError: If the parallelization is not supported.
    """
    if parallelization not in supported_parallelization:
        raise ValueError(
            f"ERROR ({function_name}): {parallelization} is not supported."
        )

class RKFunction:
    """
    A class to represent Runge-Kutte (RK) substep functions in C/C++.

    :param operator: The operator with respect to which the derivative is taken.
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
            # Add points in launcher body to compute Block/Grid kernel launch parameters
            self.body += "\n".join(
                f"const int Nxx_plus_2NGHOSTS{X} = params->Nxx_plus_2NGHOSTS{X};"
                for X in ["0", "1", "2"]
            ) + "\n"
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
            if not lhs_var in self.params:
                self.kernel_params[lhs_var] = "REAL *restrict"
                self.params += f"REAL *restrict {lhs_var},"
        self.params += f"const REAL dt"
        self.kernel_params["dt"] = "const REAL"

        kernel_body += self.loop_body.replace("commondata->dt", "dt") + "\n}\n"
        self.kernel_params_lst = [f"{v} {k}" for k, v in self.kernel_params.items()]

        if parallelization == "cuda":
            gpu_kernel_params = self.kernel_params.copy()
            gpu_kernel_params.pop("params")
            device_kernel = gputils.GPU_Kernel(
                kernel_body,
                gpu_kernel_params,
                f"{self.name}_gpu",
                launch_dict={
                    "blocks_per_grid": [
                        "(Ntot + threads_in_x_dir - 1) / threads_in_x_dir"
                    ],
                    "threads_per_block": ["32"],
                    "stream": "params->grid_idx % NUM_STREAMS",
                },
                comments=f"GPU Kernel to compute RK substep {self.rk_step}.",
            )
            self.body += device_kernel.launch_block
            self.body += device_kernel.c_function_call()
            prefunc = device_kernel.CFunction.full_function
        else:
            if self.enable_intrinsics:
                kernel_body = kernel_body.replace("dt", "DT")
                kernel_body = "const REAL_SIMD_ARRAY DT = ConstSIMD(dt);\n" + kernel_body
            rk_substep_prefunc = cfc.CFunction(
                desc=self.desc,
                cfunc_type=self.cfunc_type,
                name=self.name,
                params=",".join(self.kernel_params_lst),
                body=kernel_body,
            )
            c_function_call: str = f"{self.name}("
            for p in self.kernel_params:
                c_function_call += f"{p}, "
            c_function_call = c_function_call[:-2] + ")"
            prefunc = rk_substep_prefunc.full_function
            self.body += f"{c_function_call};\n"

        # Store CFunction
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


# Store RK substeps as RKFunctions
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


def is_diagonal_Butcher(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    key: str,
) -> bool:
    """
    Check if a given Butcher table (from Butcher_dict) is diagonal.

    :param Butcher_dict: A dictionary containing Butcher tables. Each key maps to a table (list of lists).
    :param key: The key to retrieve the Butcher table from the Butcher_dict.

    :return: True if the table is diagonal, False otherwise.

    Note:
    A Butcher table is diagonal if all its non-diagonal elements are zero.
    """
    # Get the Butcher table corresponding to the provided key
    Butcher = Butcher_dict[key][0]

    # Establish the number of rows to check for diagonal trait, all but the last row
    L = len(Butcher) - 1

    # Check all the desired rows
    for i in range(L):
        # Check each element before the diagonal element in a row
        for j in range(1, i):
            # If any non-diagonal coefficient is non-zero, then the table is not diagonal
            if Butcher[i][j] != sp.sympify(0):
                return False

    # If all non-diagonal elements are zero, the table is diagonal
    return True


def generate_gridfunction_names(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str = "RK4",
) -> Tuple[str, List[str], str, Union[str, None]]:
    """
    Generate gridfunction names for the specified Method of Lines (MoL) method.
    Used for setting up MoL_malloc, MoL_step_forward_in_time, MoL_free, and BHaH_defines.h

    :param Butcher_dict: Dictionary of Butcher tables for the MoL method.
    :param MoL_method: The MoL method to generate gridfunction names for.
    :return: A tuple containing y_n_gridfunctions, non_y_n_gridfunctions_list,
             diagnostic_gridfunctions_point_to, and diagnostic_gridfunctions2_point_to.

    Doctests:
    >>> Butcher_dict = generate_Butcher_tables()
    >>> generate_gridfunction_names(Butcher_dict, "RK2 Heun")
    ('y_n_gfs', ['y_nplus1_running_total_gfs', 'k_odd_gfs', 'k_even_gfs', 'auxevol_gfs'], 'y_nplus1_running_total_gfs', 'k_odd_gfs')
    >>> generate_gridfunction_names(Butcher_dict, "RK3")
    ('y_n_gfs', ['next_y_input_gfs', 'k1_gfs', 'k2_gfs', 'k3_gfs', 'auxevol_gfs'], 'k1_gfs', 'k2_gfs')
    >>> generate_gridfunction_names(Butcher_dict, "RK3 Ralston")
    ('y_n_gfs', ['k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs', 'k2_or_y_nplus_a32_k2_gfs', 'auxevol_gfs'], 'k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs', 'k2_or_y_nplus_a32_k2_gfs')
    >>> generate_gridfunction_names(Butcher_dict, "RK4")
    ('y_n_gfs', ['y_nplus1_running_total_gfs', 'k_odd_gfs', 'k_even_gfs', 'auxevol_gfs'], 'y_nplus1_running_total_gfs', 'k_odd_gfs')
    """
    # y_n_gridfunctions store data for the vector of gridfunctions y_i at t_n,
    # the start of each MoL timestep.
    y_n_gridfunctions = "y_n_gfs"

    # non_y_n_gridfunctions are needed to compute the data at t_{n+1}.
    # Often labeled with k_i in the name, these gridfunctions are *not*
    # needed at the start of each timestep.
    non_y_n_gridfunctions_list = []

    # Diagnostic output gridfunctions diagnostic_output_gfs & diagnostic_output_gfs2.
    if is_diagonal_Butcher(Butcher_dict, MoL_method) and "RK3" in MoL_method:
        non_y_n_gridfunctions_list.extend(
            [
                "k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs",
                "k2_or_y_nplus_a32_k2_gfs",
            ]
        )
        diagnostic_gridfunctions_point_to = (
            "k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs"
        )
        diagnostic_gridfunctions2_point_to = "k2_or_y_nplus_a32_k2_gfs"
    else:
        if not is_diagonal_Butcher(Butcher_dict, MoL_method):
            # Allocate memory for non-diagonal Butcher tables
            # Determine the number of k_i steps based on length of Butcher Table
            num_k = len(Butcher_dict[MoL_method][0]) - 1
            # For non-diagonal tables, an intermediate gridfunction "next_y_input" is used for rhs evaluations
            non_y_n_gridfunctions_list.append("next_y_input_gfs")
            for i in range(num_k):  # Need to allocate all k_i steps for a given method
                non_y_n_gridfunctions_list.append(f"k{i + 1}_gfs")
            diagnostic_gridfunctions_point_to = "k1_gfs"
            if "k2_gfs" in non_y_n_gridfunctions_list:
                diagnostic_gridfunctions2_point_to = "k2_gfs"
            else:
                print(
                    "MoL WARNING: No gridfunction group available for diagnostic_output_gfs2"
                )
                diagnostic_gridfunctions2_point_to = None
        else:
            # Allocate memory for diagonal Butcher tables, which use a "y_nplus1_running_total gridfunction"
            non_y_n_gridfunctions_list.append("y_nplus1_running_total_gfs")
            if MoL_method != "Euler":
                # Allocate memory for diagonal Butcher tables that aren't Euler
                # Need k_odd for k_1,3,5... and k_even for k_2,4,6...
                non_y_n_gridfunctions_list.extend(["k_odd_gfs", "k_even_gfs"])
            diagnostic_gridfunctions_point_to = "y_nplus1_running_total_gfs"
            diagnostic_gridfunctions2_point_to = "k_odd_gfs"

    non_y_n_gridfunctions_list.append("auxevol_gfs")

    return (
        y_n_gridfunctions,
        non_y_n_gridfunctions_list,
        diagnostic_gridfunctions_point_to,
        diagnostic_gridfunctions2_point_to,
    )


def register_CFunction_MoL_malloc(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
    which_gfs: str,
    parallelization: str = "openmp",
) -> None:
    """
    Register MoL_malloc_y_n_gfs() and MoL_malloc_non_y_n_gfs(), allocating memory for the gridfunctions indicated.

    :param Butcher_dict: Dictionary of Butcher tables for the MoL method.
    :param MoL_method: Method for the Method of Lines.
    :param which_gfs: Specifies which gridfunctions to consider ("y_n_gfs" or "non_y_n_gfs").
    :param parallelization: Parameter to specify parallelization (openmp or cuda).

    :raises ValueError: If the which_gfs parameter is neither "y_n_gfs" nor "non_y_n_gfs".

    Doctest: FIXME
    # >>> register_CFunction_MoL_malloc("Euler", "y_n_gfs")
    """
    check_supported_parallelization("register_CFunction_MoL_malloc", parallelization)
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    # Create a description for the function
    desc = f"""Method of Lines (MoL) for "{MoL_method}" method: Allocate memory for "{which_gfs}" gridfunctions
   - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep
   - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method"""

    cfunc_type = "void"

    (
        y_n_gridfunctions,
        non_y_n_gridfunctions_list,
        diagnostic_gridfunctions_point_to,
        diagnostic_gridfunctions2_point_to,
    ) = generate_gridfunction_names(Butcher_dict, MoL_method=MoL_method)

    # Determine which gridfunctions to allocate memory for
    if which_gfs == "y_n_gfs":
        gridfunctions_list = [y_n_gridfunctions]
    elif which_gfs == "non_y_n_gfs":
        gridfunctions_list = non_y_n_gridfunctions_list
    else:
        raise ValueError(f'ERROR: which_gfs = "{which_gfs}" unrecognized.')

    name = f"MoL_malloc_{which_gfs}"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs"

    # Generate the body of the function
    body = "const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;\n"
    for gridfunctions in gridfunctions_list:
        num_gfs = (
            "NUM_EVOL_GFS" if gridfunctions != "auxevol_gfs" else "NUM_AUXEVOL_GFS"
        )
        # Don't malloc a zero-sized array.
        if num_gfs == "NUM_AUXEVOL_GFS":
            body += "  if(NUM_AUXEVOL_GFS > 0) "
        if parallelization == "cuda":
            body += (
                f"cudaMalloc(&gridfuncs->{gridfunctions}, sizeof(REAL) * {num_gfs} * "
                "Nxx_plus_2NGHOSTS_tot);\n"
                'cudaCheckErrors(malloc, "Malloc failed");\n'
            )
        else:
            body += (
                f"gridfuncs->{gridfunctions} = (REAL *restrict)malloc(sizeof(REAL) * {num_gfs} * "
                "Nxx_plus_2NGHOSTS_tot);\n"
            )

    body += f"\ngridfuncs->diagnostic_output_gfs  = gridfuncs->{diagnostic_gridfunctions_point_to};\n"
    body += f"gridfuncs->diagnostic_output_gfs2 = gridfuncs->{diagnostic_gridfunctions2_point_to};\n"

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )


# single_RK_substep_input_symbolic() performs necessary replacements to
#   define C code for a single RK substep
#   (e.g., computing k_1 and then updating the outer boundaries)
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
    # Ensure all input lists are lists
    RK_lhs_list = [RK_lhs_list] if not isinstance(RK_lhs_list, list) else RK_lhs_list
    RK_rhs_list = [RK_rhs_list] if not isinstance(RK_rhs_list, list) else RK_rhs_list
    post_rhs_list = (
        [post_rhs_list] if not isinstance(post_rhs_list, list) else post_rhs_list
    )
    post_rhs_output_list = (
        [post_rhs_output_list]
        if not isinstance(post_rhs_output_list, list)
        else post_rhs_output_list
    )
    comment_block = (
        "// ***Euler timestepping only requires one RHS evaluation***"
        if rk_step is None
        else f"// -={{ START k{rk_step} substep }}=-"
    )
    comment_block += additional_comments
    body = f"{comment_block}\n"

    if isinstance(substep_time_offset_dt, (int, sp.Rational, sp.Mul)):
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

    # Part 1: RHS evaluation
    updated_rhs_str = (
        str(rhs_str)
        .replace("RK_INPUT_GFS", str(rhs_input_expr).replace("gfsL", "gfs"))
        .replace("RK_OUTPUT_GFS", str(rhs_output_expr).replace("gfsL", "gfs"))
    )
    body += updated_rhs_str + "\n"

    # Part 2: RK update
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

    # Part 3: Call post-RHS functions
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


########################################################################################################################
# EXAMPLE
# ODE: y' = f(t,y), y(t_0) = y_0
# Starting at time t_n with solution having value y_n and trying to update to y_nplus1 with timestep dt

# Example of scheme for RK4 with k_1, k_2, k_3, k_4 (Using non-diagonal algorithm) Notice this requires storage of
# y_n, y_nplus1, k_1 through k_4

# k_1      = dt*f(t_n, y_n)
# k_2      = dt*f(t_n + 1/2*dt, y_n + 1/2*k_1)
# k_3      = dt*f(t_n + 1/2*dt, y_n + 1/2*k_2)
# k_4      = dt*f(t_n + dt, y_n + k_3)
# y_nplus1 = y_n + 1/3k_1 + 1/6k_2 + 1/6k_3 + 1/3k_4

# Example of scheme RK4 using only k_odd and k_even (Diagonal algroithm) Notice that this only requires storage


# k_odd     = dt*f(t_n, y_n)
# y_nplus1  = 1/3*k_odd
# k_even    = dt*f(t_n + 1/2*dt, y_n + 1/2*k_odd)
# y_nplus1 += 1/6*k_even
# k_odd     = dt*f(t_n + 1/2*dt, y_n + 1/2*k_even)
# y_nplus1 += 1/6*k_odd
# k_even    = dt*f(t_n + dt, y_n + k_odd)
# y_nplus1 += 1/3*k_even
########################################################################################################################
def register_CFunction_MoL_step_forward_in_time(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
    rhs_string: str = "",
    post_rhs_string: str = "",
    post_post_rhs_string: str = "",
    enable_rfm_precompute: bool = False,
    enable_curviBCs: bool = False,
    enable_intrinsics: bool = False,
    parallelization: str = "openmp",
    rational_const_alias: str = "const",
) -> None:
    """
    Register MoL_step_forward_in_time() C function, which is the core driver for time evolution in BHaH codes.

    :param Butcher_dict: A dictionary containing the Butcher tables for various RK-like methods.
    :param MoL_method: The method of lines (MoL) used for time-stepping.
    :param rhs_string: Right-hand side string of the C code.
    :param post_rhs_string: Input string for post-RHS phase in the C code.
    :param post_post_rhs_string: String to be used after the post-RHS phase.
    :param enable_rfm_precompute: Flag to enable reference metric functionality.
    :param enable_curviBCs: Flag to enable curvilinear boundary conditions.
    :param enable_intrinsics: A flag to specify if hardware instructions should be used.
    :param parallelization: Parameter to specify parallelization (openmp or cuda).
    :param rational_const_alias: Overload const specifier for Rational definitions
    :raises ValueError: If unsupported Butcher table specified since adaptive RK steps are not implemented in MoL.

    Doctest:
    >>> import nrpy.c_function as cfc
    >>> from nrpy.infrastructures.BHaH.MoLtimestepping.MoL import register_CFunction_MoL_step_forward_in_time
    >>> from nrpy.infrastructures.BHaH.MoLtimestepping.RK_Butcher_Table_Dictionary import (
    ...     generate_Butcher_tables,
    ... )
    >>> Butcher_dict = generate_Butcher_tables()
    >>> expected_str_dict=dict()
    >>> try:
    ...     register_CFunction_MoL_step_forward_in_time(Butcher_dict, "AHE")
    ... except ValueError as e:
    ...     print(f"ValueError: {e.args[0]}")
    ValueError: Adaptive order Butcher tables are currently not supported in MoL.
    """
    check_supported_parallelization("register_CFunction_MoL_step_forward_in_time", parallelization)
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if enable_intrinsics and parallelization == "cuda":
        includes += [os.path.join("intrinsics", "cuda_intrinsics.h")]
    elif enable_intrinsics:
        includes += [os.path.join("intrinsics", "simd_intrinsics.h")]

    desc = f'Method of Lines (MoL) for "{MoL_method}" method: Step forward one full timestep.\n'
    cfunc_type = "void"
    name = "MoL_step_forward_in_time"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    # Code body
    body = f"""
// C code implementation of -={{ {MoL_method} }}=- Method of Lines timestepping.

// First set the initial time:
const REAL time_start = commondata->time;
"""
    (
        y_n_gridfunctions,
        non_y_n_gridfunctions_list,
        _throwaway,
        _throwaway2,
    ) = generate_gridfunction_names(Butcher_dict, MoL_method)

    gf_prefix = "griddata[grid].gridfuncs."

    gf_aliases = f"""// Set gridfunction aliases from gridfuncs struct
// y_n gridfunctions
MAYBE_UNUSED REAL *restrict {y_n_gridfunctions} = {gf_prefix}{y_n_gridfunctions};
// Temporary timelevel & AUXEVOL gridfunctions:\n"""
    for gf in non_y_n_gridfunctions_list:
        gf_aliases += f"MAYBE_UNUSED REAL *restrict {gf} = {gf_prefix}{gf};\n"

    gf_aliases += (
        "MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;\n"
    )
    if enable_rfm_precompute:
        gf_aliases += "MAYBE_UNUSED const rfm_struct *restrict rfmstruct = &griddata[grid].rfmstruct;\n"
    else:
        gf_aliases += "MAYBE_UNUSED REAL *restrict xx[3]; for(int ww=0;ww<3;ww++) xx[ww] = griddata[grid].xx[ww];\n"

    if enable_curviBCs:
        gf_aliases += "MAYBE_UNUSED const bc_struct *restrict bcstruct = &griddata[grid].bcstruct;\n"
    for i in ["0", "1", "2"]:
        gf_aliases += f"MAYBE_UNUSED const int Nxx_plus_2NGHOSTS{i} = griddata[grid].params.Nxx_plus_2NGHOSTS{i};\n"

    # Implement Method of Lines (MoL) Timestepping
    rk_step_body_dict: Dict[str, str] = {}
    Butcher = Butcher_dict[MoL_method][
        0
    ]  # Get the desired Butcher table from the dictionary

    if not Butcher[-1][0] == "":
        raise ValueError(
            "Adaptive order Butcher tables are currently not supported in MoL."
        )

    # Note: this calculation does not work for adaptive steps
    num_steps = (
        len(Butcher) - 1
    )  # Specify the number of required steps to update solution

    dt = sp.Symbol("commondata->dt", real=True)

    if is_diagonal_Butcher(Butcher_dict, MoL_method) and "RK3" in MoL_method:
        # Diagonal RK3 only!!!
        #  In a diagonal RK3 method, only 3 gridfunctions need be defined. Below implements this approach.
        y_n_gfs = sp.Symbol("y_n_gfsL", real=True)
        k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs = sp.Symbol(
            "k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfsL", real=True
        )
        k2_or_y_nplus_a32_k2_gfs = sp.Symbol("k2_or_y_nplus_a32_k2_gfsL", real=True)
        # k_1
        body += """
// In a diagonal RK3 method like this one, only 3 gridfunctions need be defined. Below implements this approach.
// Using y_n_gfs as input, k1 and apply boundary conditions\n"""
        rk_step_body_dict["RK_SUBSTEP_K1"] = (
            single_RK_substep_input_symbolic(
                additional_comments="""
// RHS evaluation:
//  1. We will store k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs now as
//     ...  the update for the next rhs evaluation y_n + a21*k1*dt
// Post-RHS evaluation:
//  1. Apply post-RHS to y_n + a21*k1*dt""",
                substep_time_offset_dt=Butcher[0][0],
                rhs_str=rhs_string,
                rhs_input_expr=y_n_gfs,
                rhs_output_expr=k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs,
                RK_lhs_list=[k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs],
                RK_rhs_list=[
                    Butcher[1][1]
                    * k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs
                    * dt
                    + y_n_gfs
                ],
                rk_step=1,
                post_rhs_list=[post_rhs_string],
                post_rhs_output_list=[
                    k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs
                ],
                enable_intrinsics=enable_intrinsics,
                gf_aliases=gf_aliases,
                post_post_rhs_string=post_post_rhs_string,
                rational_const_alias=rational_const_alias,
                parallelization=parallelization,
            )
            + "// -={ END k1 substep }=-\n\n"
        )

        # k_2
        rk_step_body_dict["RK_SUBSTEP_K2"] = (
            single_RK_substep_input_symbolic(
                additional_comments="""
// RHS evaluation:
//    1. Reassign k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs to be the running total y_{n+1}; a32*k2*dt to the running total
//    2. Store k2_or_y_nplus_a32_k2_gfs now as y_n + a32*k2*dt
// Post-RHS evaluation:
//    1. Apply post-RHS to both y_n + a32*k2 (stored in k2_or_y_nplus_a32_k2_gfs)
//       ... and the y_{n+1} running total, as they have not been applied yet to k2-related gridfunctions""",
                rhs_str=rhs_string,
                substep_time_offset_dt=Butcher[1][0],
                rhs_input_expr=k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs,
                rhs_output_expr=k2_or_y_nplus_a32_k2_gfs,
                RK_lhs_list=[
                    k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs,
                    k2_or_y_nplus_a32_k2_gfs,
                ],
                RK_rhs_list=[
                    Butcher[3][1]
                    * (k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs - y_n_gfs)
                    / Butcher[1][1]
                    + y_n_gfs
                    + Butcher[3][2] * k2_or_y_nplus_a32_k2_gfs * dt,
                    Butcher[2][2] * k2_or_y_nplus_a32_k2_gfs * dt + y_n_gfs,
                ],
                rk_step=2,
                post_rhs_list=[post_rhs_string, post_rhs_string],
                post_rhs_output_list=[
                    k2_or_y_nplus_a32_k2_gfs,
                    k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs,
                ],
                enable_intrinsics=enable_intrinsics,
                gf_aliases=gf_aliases,
                post_post_rhs_string=post_post_rhs_string,
                rational_const_alias=rational_const_alias,
                parallelization=parallelization,
            )
            + "// -={ END k2 substep }=-\n\n"
        )

        # k_3
        rk_step_body_dict["RK_SUBSTEP_K3"] = (
            single_RK_substep_input_symbolic(
                additional_comments="""
// RHS evaluation:
//    1. Add k3 to the running total and save to y_n
// Post-RHS evaluation:
//    1. Apply post-RHS to y_n""",
                substep_time_offset_dt=Butcher[2][0],
                rhs_str=rhs_string,
                rhs_input_expr=k2_or_y_nplus_a32_k2_gfs,
                rhs_output_expr=y_n_gfs,
                RK_lhs_list=[y_n_gfs],
                RK_rhs_list=[
                    k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs
                    + Butcher[3][3] * y_n_gfs * dt
                ],
                rk_step=3,
                post_rhs_list=[post_rhs_string],
                post_rhs_output_list=[y_n_gfs],
                enable_intrinsics=enable_intrinsics,
                gf_aliases=gf_aliases,
                post_post_rhs_string=post_post_rhs_string,
                rational_const_alias=rational_const_alias,
                parallelization=parallelization,
            )
            + "// -={ END k3 substep }=-\n\n"
        )
    else:
        y_n = sp.Symbol("y_n_gfsL", real=True)
        if not is_diagonal_Butcher(Butcher_dict, MoL_method):
            for s in range(num_steps):
                next_y_input = sp.Symbol("next_y_input_gfsL", real=True)

                # If we're on the first step (s=0), we use y_n gridfunction as input.
                #      Otherwise next_y_input is input. Output is just the reverse.
                if s == 0:  # If on first step:
                    rhs_input = y_n
                else:  # If on second step or later:
                    rhs_input = next_y_input
                rhs_output = sp.Symbol("k" + str(s + 1) + "_gfs", real=True)
                if s == num_steps - 1:  # If on final step:
                    RK_lhs = y_n
                else:  # If on anything but the final step:
                    RK_lhs = next_y_input
                RK_rhs = y_n
                for m in range(s + 1):
                    k_mp1_gfs = sp.Symbol("k" + str(m + 1) + "_gfsL")
                    if Butcher[s + 1][m + 1] != 0:
                        if Butcher[s + 1][m + 1] != 1:
                            RK_rhs += dt * k_mp1_gfs * Butcher[s + 1][m + 1]
                        else:
                            RK_rhs += dt * k_mp1_gfs

                post_rhs = post_rhs_string
                if s == num_steps - 1:  # If on final step:
                    post_rhs_output = y_n
                else:  # If on anything but the final step:
                    post_rhs_output = next_y_input

                rk_step_body_dict[
                    f"RK_SUBSTEP_K{s+1}"
                ] = f"""{single_RK_substep_input_symbolic(
                    substep_time_offset_dt=Butcher[s][0],
                    rhs_str=rhs_string,
                    rhs_input_expr=rhs_input,
                    rhs_output_expr=rhs_output,
                    RK_lhs_list=[RK_lhs],
                    RK_rhs_list=[RK_rhs],
                    rk_step=s+1,
                    post_rhs_list=[post_rhs],
                    post_rhs_output_list=[post_rhs_output],
                    enable_intrinsics=enable_intrinsics,
                    gf_aliases=gf_aliases,
                    post_post_rhs_string=post_post_rhs_string,
                    rational_const_alias=rational_const_alias,
                    parallelization=parallelization
                )}// -={{ END k{str(s + 1)} substep }}=-\n\n"""

        else:
            y_n = sp.Symbol("y_n_gfsL", real=True)
            y_nplus1_running_total = sp.Symbol("y_nplus1_running_total_gfsL", real=True)
            if (
                MoL_method == "Euler"
            ):  # Euler's method doesn't require any k_i, and gets its own unique algorithm
                rk_step_body_dict[f"{MoL_method}"] = single_RK_substep_input_symbolic(
                    additional_comments="// ***Euler timestepping only requires one RHS evaluation***",
                    substep_time_offset_dt=Butcher[0][0],
                    rhs_str=rhs_string,
                    rhs_input_expr=y_n,
                    rhs_output_expr=y_nplus1_running_total,
                    RK_lhs_list=[y_n],
                    RK_rhs_list=[y_n + y_nplus1_running_total * dt],
                    post_rhs_list=[post_rhs_string],
                    post_rhs_output_list=[y_n],
                    enable_intrinsics=enable_intrinsics,
                    gf_aliases=gf_aliases,
                    post_post_rhs_string=post_post_rhs_string,
                    rational_const_alias=rational_const_alias,
                    parallelization=parallelization,
                )
            else:
                for s in range(num_steps):
                    # If we're on the first step (s=0), we use y_n gridfunction as input.
                    # and k_odd as output.
                    if s == 0:
                        rhs_input = sp.Symbol("y_n_gfsL", real=True)
                        rhs_output = sp.Symbol("k_odd_gfsL", real=True)
                    # For the remaining steps the inputs and ouputs alternate between k_odd and k_even
                    elif s % 2 == 0:
                        rhs_input = sp.Symbol("k_even_gfsL", real=True)
                        rhs_output = sp.Symbol("k_odd_gfsL", real=True)
                    else:
                        rhs_input = sp.Symbol("k_odd_gfsL", real=True)
                        rhs_output = sp.Symbol("k_even_gfsL", real=True)
                    RK_lhs_list: List[sp.Basic] = []
                    RK_rhs_list = []
                    if s != num_steps - 1:  # For anything besides the final step
                        if s == 0:  # The first RK step
                            RK_lhs_list.append(y_nplus1_running_total)
                            RK_rhs_list.append(
                                rhs_output * dt * Butcher[num_steps][s + 1]
                            )

                            RK_lhs_list.append(rhs_output)
                            RK_rhs_list.append(
                                y_n + rhs_output * dt * Butcher[s + 1][s + 1]
                            )
                        else:
                            if Butcher[num_steps][s + 1] != 0:
                                RK_lhs_list.append(y_nplus1_running_total)
                                if Butcher[num_steps][s + 1] != 1:
                                    RK_rhs_list.append(
                                        y_nplus1_running_total
                                        + rhs_output * dt * Butcher[num_steps][s + 1]
                                    )
                                else:
                                    RK_rhs_list.append(
                                        y_nplus1_running_total + rhs_output * dt
                                    )
                            if Butcher[s + 1][s + 1] != 0:
                                RK_lhs_list.append(rhs_output)
                                if Butcher[s + 1][s + 1] != 1:
                                    RK_rhs_list.append(
                                        y_n + rhs_output * dt * Butcher[s + 1][s + 1]
                                    )
                                else:
                                    RK_rhs_list.append(y_n + rhs_output * dt)
                        post_rhs_output = rhs_output
                    if s == num_steps - 1:  # If on the final step
                        if Butcher[num_steps][s + 1] != 0:
                            RK_lhs_list.append(y_n)
                            if Butcher[num_steps][s + 1] != 1:
                                RK_rhs_list.append(
                                    y_n
                                    + y_nplus1_running_total
                                    + rhs_output * dt * Butcher[num_steps][s + 1]
                                )
                            else:
                                RK_rhs_list.append(
                                    y_n + y_nplus1_running_total + rhs_output * dt
                                )
                        post_rhs_output = y_n
                    rk_step_body_dict[f"RK_SUBSTEP_K{s+1}"] = (
                        single_RK_substep_input_symbolic(
                            substep_time_offset_dt=Butcher[s][0],
                            rhs_str=rhs_string,
                            rhs_input_expr=rhs_input,
                            rhs_output_expr=rhs_output,
                            RK_lhs_list=RK_lhs_list,
                            RK_rhs_list=RK_rhs_list,
                            post_rhs_list=[post_rhs_string],
                            post_rhs_output_list=[post_rhs_output],
                            rk_step=s + 1,
                            enable_intrinsics=enable_intrinsics,
                            gf_aliases=gf_aliases,
                            post_post_rhs_string=post_post_rhs_string,
                            rational_const_alias=rational_const_alias,
                            parallelization=parallelization,
                        )
                        + f"// -={{ END k{s + 1} substep }}=-\n\n"
                    )

    if parallelization == "cuda":
        prefunc = """
#define LOOP_ALL_GFS_GPS(ii) \
const int tid0 = threadIdx.x + blockIdx.x*blockDim.x; \
const int stride0 = blockDim.x * gridDim.x; \
  for(int (ii)=(tid0);(ii)<d_params[streamid].Nxx_plus_2NGHOSTS0*d_params[streamid].Nxx_plus_2NGHOSTS1*d_params[streamid].Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;(ii)+=(stride0))
"""
    else:
        prefunc = """
#define LOOP_ALL_GFS_GPS(ii) \
_Pragma("omp parallel for") \
  for(int (ii)=0;(ii)<params->Nxx_plus_2NGHOSTS0*params->Nxx_plus_2NGHOSTS1*params->Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;(ii)++)
"""
        if enable_intrinsics:
            warnings.warn(
                "SIMD intrinsics in MoL is not properly supported -- MoL update loops are not properly bounds checked."
            )
            prefunc = prefunc.replace("(ii)++", "(ii) += (simd_width)")
    prefunc += construct_RK_functions_prefunc()

    for _, v in rk_step_body_dict.items():
        body += v
    body += """
// Adding dt to commondata->time many times will induce roundoff error,
//   so here we set time based on the iteration number.
commondata->time = (REAL)(commondata->nn + 1) * commondata->dt;

// Finally, increment the timestep n:
commondata->nn++;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        prefunc=prefunc,
    )


# register_CFunction_MoL_free_memory() registers
#           MoL_free_memory_y_n_gfs() and
#           MoL_free_memory_non_y_n_gfs(), which free memory for
#           the indicated sets of gridfunctions
def register_CFunction_MoL_free_memory(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
    which_gfs: str,
    parallelization: str = "openmp",
) -> None:
    """
    Free memory for the specified Method of Lines (MoL) gridfunctions, given an MoL_method.

    :param Butcher_dict: Dictionary containing Butcher tableau for MoL methods.
    :param MoL_method: The Method of Lines method.
    :param which_gfs: The gridfunctions to be freed, either 'y_n_gfs' or 'non_y_n_gfs'.
    :param parallelization: Parameter to specify parallelization (openmp or cuda).

    :raises ValueError: If the 'which_gfs' argument is unrecognized.
    """
    check_supported_parallelization("register_CFunction_MoL_free_memory", parallelization)
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f'Method of Lines (MoL) for "{MoL_method}" method: Free memory for "{which_gfs}" gridfunctions\n'
    desc += "   - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep\n"
    desc += "   - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method\n"
    cfunc_type = "void"
    mem_free_func = "free" if parallelization != "cuda" else "cudaFree"

    (
        y_n_gridfunctions,
        non_y_n_gridfunctions_list,
        _diagnostic_gridfunctions_point_to,
        _diagnostic_gridfunctions2_point_to,
    ) = generate_gridfunction_names(Butcher_dict, MoL_method)

    if which_gfs == "y_n_gfs":
        gridfunctions_list = [y_n_gridfunctions]
    elif which_gfs == "non_y_n_gfs":
        gridfunctions_list = non_y_n_gridfunctions_list
    else:
        raise ValueError(f'ERROR: which_gfs = "{which_gfs}" unrecognized.')

    name: str = f"MoL_free_memory_{which_gfs}"
    params = "MoL_gridfunctions_struct *restrict gridfuncs"
    body = ""
    for gridfunction in gridfunctions_list:
        # Don't free a zero-sized array.
        if gridfunction == "auxevol_gfs":
            body += (
                f"  if(NUM_AUXEVOL_GFS > 0) {mem_free_func}(gridfuncs->{gridfunction});"
            )
        else:
            body += f"  {mem_free_func}(gridfuncs->{gridfunction});"
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )


# Register all the CFunctions and NRPy basic defines
def register_CFunctions(
    MoL_method: str = "RK4",
    rhs_string: str = "rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, RK_INPUT_GFS, RK_OUTPUT_GFS);",
    post_rhs_string: str = "apply_bcs(Nxx, Nxx_plus_2NGHOSTS, RK_OUTPUT_GFS);",
    post_post_rhs_string: str = "",
    enable_rfm_precompute: bool = False,
    enable_curviBCs: bool = False,
    enable_intrinsics: bool = False,
    register_MoL_step_forward_in_time: bool = True,
    parallelization: str = "openmp",
    rational_const_alias: str = "const",
) -> None:
    r"""
    Register all MoL C functions and NRPy basic defines.

    :param MoL_method: The method to be used for MoL. Default is 'RK4'.
    :param rhs_string: RHS function call as string. Default is "rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, RK_INPUT_GFS, RK_OUTPUT_GFS);"
    :param post_rhs_string: Post-RHS function call as string. Default is "apply_bcs(Nxx, Nxx_plus_2NGHOSTS, RK_OUTPUT_GFS);"
    :param post_post_rhs_string: Post-post-RHS function call as string. Default is an empty string.
    :param enable_rfm_precompute: Enable reference metric support. Default is False.
    :param enable_curviBCs: Enable curvilinear boundary conditions. Default is False.
    :param enable_intrinsics: A flag to specify if hardware instructions should be used.
    :param register_MoL_step_forward_in_time: Whether to register the MoL step forward function. Default is True.
    :param parallelization: Parameter to specify parallelization (openmp or cuda).
    :param rational_const_alias: Overload const specifier for Rational definitions

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunctions()
    >>> validate_strings(cfc.CFunction_dict["MoL_step_forward_in_time"].full_function, "MoL_step_forward_in_time")
    >>> sorted(cfc.CFunction_dict.keys())
    ['MoL_free_memory_non_y_n_gfs', 'MoL_free_memory_y_n_gfs', 'MoL_malloc_non_y_n_gfs', 'MoL_malloc_y_n_gfs', 'MoL_step_forward_in_time']
    >>> print(cfc.CFunction_dict["MoL_free_memory_non_y_n_gfs"].full_function)
    #include "BHaH_defines.h"
    #include "BHaH_function_prototypes.h"
    /**
     * Method of Lines (MoL) for "RK4" method: Free memory for "non_y_n_gfs" gridfunctions
     * - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep
     * - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method
     *
     */
    void MoL_free_memory_non_y_n_gfs(MoL_gridfunctions_struct *restrict gridfuncs) {
      free(gridfuncs->y_nplus1_running_total_gfs);
      free(gridfuncs->k_odd_gfs);
      free(gridfuncs->k_even_gfs);
      if (NUM_AUXEVOL_GFS > 0)
        free(gridfuncs->auxevol_gfs);
    } // END FUNCTION MoL_free_memory_non_y_n_gfs
    <BLANKLINE>
    >>> print(cfc.CFunction_dict["MoL_malloc_non_y_n_gfs"].full_function)
    #include "BHaH_defines.h"
    #include "BHaH_function_prototypes.h"
    /**
     * Method of Lines (MoL) for "RK4" method: Allocate memory for "non_y_n_gfs" gridfunctions
     * - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep
     * - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method
     */
    void MoL_malloc_non_y_n_gfs(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                MoL_gridfunctions_struct *restrict gridfuncs) {
    #include "set_CodeParameters.h"
      const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
      gridfuncs->y_nplus1_running_total_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
      gridfuncs->k_odd_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
      gridfuncs->k_even_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
      if (NUM_AUXEVOL_GFS > 0)
        gridfuncs->auxevol_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_AUXEVOL_GFS * Nxx_plus_2NGHOSTS_tot);
    <BLANKLINE>
      gridfuncs->diagnostic_output_gfs = gridfuncs->y_nplus1_running_total_gfs;
      gridfuncs->diagnostic_output_gfs2 = gridfuncs->k_odd_gfs;
    } // END FUNCTION MoL_malloc_non_y_n_gfs
    <BLANKLINE>
    """
    check_supported_parallelization("register_CFunctions", parallelization)
    Butcher_dict = generate_Butcher_tables()

    for which_gfs in ["y_n_gfs", "non_y_n_gfs"]:
        register_CFunction_MoL_malloc(
            Butcher_dict, MoL_method, which_gfs, parallelization=parallelization
        )
        register_CFunction_MoL_free_memory(
            Butcher_dict, MoL_method, which_gfs, parallelization=parallelization
        )
    if register_MoL_step_forward_in_time:
        register_CFunction_MoL_step_forward_in_time(
            Butcher_dict,
            MoL_method,
            rhs_string,
            post_rhs_string,
            post_post_rhs_string,
            enable_rfm_precompute=enable_rfm_precompute,
            enable_curviBCs=enable_curviBCs,
            enable_intrinsics=enable_intrinsics,
            parallelization=parallelization,
            rational_const_alias=rational_const_alias,
        )

    griddata_commondata.register_griddata_commondata(
        __name__, "MoL_gridfunctions_struct gridfuncs", "MoL gridfunctions"
    )

    # Generating gridfunction names based on the given MoL method
    (
        y_n_gridfunctions,
        non_y_n_gridfunctions_list,
        _diagnostic_gridfunctions_point_to,
        _diagnostic_gridfunctions2_point_to,
    ) = generate_gridfunction_names(Butcher_dict, MoL_method=MoL_method)

    # Step 3.b: Create MoL_timestepping struct:
    BHaH_MoL_body: str = (
        "typedef struct __MoL_gridfunctions_struct__ {\n"
        + f"REAL *restrict {y_n_gridfunctions};\n"
        + "".join(f"REAL *restrict {gfs};\n" for gfs in non_y_n_gridfunctions_list)
        + r"""REAL *restrict diagnostic_output_gfs;
REAL *restrict diagnostic_output_gfs2;
} MoL_gridfunctions_struct;
"""
    )

    if parallelization != "openmp":
        BHaH_MoL_body = BHaH_MoL_body.replace("REAL *restrict ", "REAL * ")
    BHaH_defines_h.register_BHaH_defines(__name__, BHaH_MoL_body)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
