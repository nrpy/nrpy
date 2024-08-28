"""
Module for producing C codes related to MoL timestepping within the superB infrastructure.
This includes implementation details and functions for allocating and deallocating the necessary memory.

Authors: Brandon Clark
         Zachariah B. Etienne (maintainer)
         zachetie **at** gmail **dot* com
         Nishita Jadoo
         njadoo **at** uidaho **dot* edu

superB changes/additions to nrpy.infrastructures.BHaH.MoLtimestepping.MoL.py:
-added time_start as parameter
-added which RK stage as parameter
-added switch case depending on RK stage
-allocate memory to diagnostic output gfs
"""

import os  # Standard Python module for multiplatform OS-level functions
import warnings
from typing import Dict, List, Tuple, Union

import sympy as sp  # Import SymPy, a computer algebra system written entirely in Python

import nrpy.c_function as cfc
import nrpy.params as par  # NRPy+: Parameter interface
from nrpy.c_codegen import c_codegen
from nrpy.helpers.generic import superfast_uniq
from nrpy.infrastructures.BHaH import BHaH_defines_h, griddata_commondata
from nrpy.infrastructures.BHaH.MoLtimestepping.MoL import (
    generate_gridfunction_names,
    is_diagonal_Butcher,
    register_CFunction_MoL_free_memory,
    register_CFunction_MoL_malloc,
)
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


# single_RK_substep_input_symbolic() performs necessary replacements to
#   define C code for a single RK substep
#   (e.g., computing k_1 and then updating the outer boundaries)
def single_RK_substep_input_symbolic(
    comment_block: str,
    substep_time_offset_dt: Union[sp.Basic, int, str],
    rhs_str: str,
    rhs_input_expr: sp.Basic,
    rhs_output_expr: sp.Basic,
    RK_lhs_list: Union[sp.Basic, List[sp.Basic]],
    RK_rhs_list: Union[sp.Basic, List[sp.Basic]],
    post_rhs_list: Union[str, List[str]],
    post_rhs_output_list: Union[sp.Basic, List[sp.Basic]],
    enable_simd: bool = False,
    gf_aliases: str = "",
    post_post_rhs_string: str = "",
    fp_type: str = "double",
) -> str:
    """
    Generate C code for a given Runge-Kutta substep.

    :param comment_block: Block of comments for the generated code.
    :param substep_time_offset_dt: Time offset for the RK substep.
    :param rhs_str: Right-hand side string of the C code.
    :param rhs_input_expr: Input expression for the RHS.
    :param rhs_output_expr: Output expression for the RHS.
    :param RK_lhs_list: List of LHS expressions for RK.
    :param RK_rhs_list: List of RHS expressions for RK.
    :param post_rhs_list: List of post-RHS expressions.
    :param post_rhs_output_list: List of outputs for post-RHS expressions.
    :param enable_simd: Whether SIMD optimization is enabled.
    :param gf_aliases: Additional aliases for grid functions.
    :param post_post_rhs_string: String to be used after the post-RHS phase.
    :param fp_type: Floating point type, e.g., "double".

    :return: A string containing the generated C code.

    :raises ValueError: If substep_time_offset_dt cannot be extracted from the Butcher table.
    """
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

    return_str = f"{comment_block}\n"
    if isinstance(substep_time_offset_dt, (int, sp.Rational, sp.Mul)):
        substep_time_offset_str = f"{float(substep_time_offset_dt):.17e}"
    else:
        raise ValueError(
            f"Could not extract substep_time_offset_dt={substep_time_offset_dt} from Butcher table"
        )
    return_str += "for(int grid=0; grid<commondata->NUMGRIDS; grid++) {\n"
    return_str += (
        f"commondata->time = time_start + {substep_time_offset_str} * commondata->dt;\n"
    )
    return_str += gf_aliases

    return_str += """
switch (which_MOL_part) {
  case MOL_PART_1: {"""

    # Part 1: RHS evaluation
    updated_rhs_str = (
        str(rhs_str)
        .replace("RK_INPUT_GFS", str(rhs_input_expr).replace("gfsL", "gfs"))
        .replace("RK_OUTPUT_GFS", str(rhs_output_expr).replace("gfsL", "gfs"))
    )
    return_str += updated_rhs_str + "\n"

    return_str += """
     break;
  }
  case MOL_PART_2: {"""

    # Part 2: RK update
    if enable_simd:
        warnings.warn(
            "enable_simd in MoL is not properly supported -- MoL update loops are not properly bounds checked."
        )
        return_str += "#pragma omp parallel for\n"
        return_str += "for(int i=0;i<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;i+=simd_width) {{\n"
    else:
        return_str += "LOOP_ALL_GFS_GPS(i) {\n"

    var_type = "REAL_SIMD_ARRAY" if enable_simd else "REAL"

    RK_lhs_str_list = [
        (
            f"const REAL_SIMD_ARRAY __rhs_exp_{i}"
            if enable_simd
            else f"{str(el).replace('gfsL', 'gfs[i]')}"
        )
        for i, el in enumerate(RK_lhs_list)
    ]

    read_list = [
        read for el in RK_rhs_list for read in list(sp.ordered(el.free_symbols))
    ]
    read_list_unique = superfast_uniq(read_list)

    for el in read_list_unique:
        if str(el) != "commondata->dt":
            if enable_simd:
                simd_el = str(el).replace("gfsL", "gfs[i]")
                return_str += f"const {var_type} {el} = ReadSIMD(&{simd_el});\n"
            else:
                return_str += (
                    f"const {var_type} {el} = {str(el).replace('gfsL', 'gfs[i]')};\n"
                )

    if enable_simd:
        return_str += "const REAL_SIMD_ARRAY DT = ConstSIMD(commondata->dt);\n"

    kernel = c_codegen(
        RK_rhs_list,
        RK_lhs_str_list,
        include_braces=False,
        verbose=False,
        enable_simd=enable_simd,
        fp_type=fp_type,
    )

    if enable_simd:
        return_str += kernel.replace("commondata->dt", "DT")
        for i, el in enumerate(RK_lhs_list):
            return_str += (
                f"  WriteSIMD(&{str(el).replace('gfsL', 'gfs[i]')}, __rhs_exp_{i});\n"
            )

    else:
        return_str += kernel

    return_str += "}\n"
    return_str += """
    break;
  }
"""
    return_str += """
  case MOL_PART_3_APPLY_BCS: {
"""
    # Part 3: Call post-RHS functions
    for post_rhs, post_rhs_output in zip(post_rhs_list, post_rhs_output_list):
        parts = [part.strip() for part in post_rhs.split(";") if part.strip()]
        part1 = parts[0] + ";"
        return_str += part1.replace(
            "RK_OUTPUT_GFS", str(post_rhs_output).replace("gfsL", "gfs")
        )
        return_str += "\n"

    return_str += """
    break;
  }
"""
    return_str += """
  case MOL_PART_3_AFTER_APPLY_BCS: {
"""
    for post_rhs, post_rhs_output in zip(post_rhs_list, post_rhs_output_list):
        parts = [part.strip() for part in post_rhs.split(";") if part.strip()]
        if len(parts) > 1:
            part2 = parts[1] + ";"
            return_str += part2.replace(
                "RK_OUTPUT_GFS", str(post_rhs_output).replace("gfsL", "gfs")
            )

    for post_rhs, post_rhs_output in zip(post_rhs_list, post_rhs_output_list):
        return_str += post_post_rhs_string.replace(
            "RK_OUTPUT_GFS", str(post_rhs_output).replace("gfsL", "gfs")
        )

    return_str += """
    break;
  }
}
}
"""

    return return_str


def generate_post_rhs_output_list(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
    rk_substep: int,
) -> List[str]:
    """
    Generate post_rhs_output_list based on the Method of Lines (MoL) method and the current RK substep.

    :param Butcher_dict: A dictionary containing the Butcher tableau and its order.
    :param MoL_method: The method of lines method name.
    :param rk_substep: The current Runge-Kutta substep.
    :return: A list of strings representing the post RHS output.
    """
    num_steps = len(Butcher_dict[MoL_method][0]) - 1
    s = rk_substep - 1  # Convert to 0-indexed
    post_rhs = []

    if is_diagonal_Butcher(Butcher_dict, MoL_method) and "RK3" in MoL_method:
        y_n_gfs = "Y_N_GFS"
        k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs = (
            "K1_OR_Y_NPLUS_A21_K1_OR_Y_NPLUS1_RUNNING_TOTAL_GFS"
        )
        k2_or_y_nplus_a32_k2_gfs = "K2_OR_Y_NPLUS_A32_K2_GFS"

        if s == 0:
            post_rhs = [k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs]
        elif s == 1:
            post_rhs = [
                k2_or_y_nplus_a32_k2_gfs,
                k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs,
            ]
        elif s == 2:
            post_rhs = [y_n_gfs]
    else:
        y_n = "Y_N_GFS"
        if not is_diagonal_Butcher(Butcher_dict, MoL_method):
            next_y_input = "NEXT_Y_INPUT_GFS"
            if s == num_steps - 1:  # If on final step
                post_rhs = [y_n]
            else:  # If on anything but the final step
                post_rhs = [next_y_input]
        else:
            if MoL_method == "Euler":
                post_rhs = [y_n]
            else:
                if s % 2 == 0:
                    rhs_output = "K_ODD_GFS"
                else:
                    rhs_output = "K_EVEN_GFS"

                if s == num_steps - 1:  # If on the final step
                    post_rhs = [y_n]
                else:  # For anything besides the final step
                    post_rhs = [rhs_output]

    return post_rhs  # Return the list of strings representing the post RHS output


def generate_rhs_output_exprs(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
    rk_substep: int,
) -> List[str]:
    """
    Generate rhs_output_exprs based on the Method of Lines (MoL) method and the current RK substep.

    :param Butcher_dict: A dictionary containing the Butcher tableau and its order.
    :param MoL_method: The method of lines method name.
    :param rk_substep: The current Runge-Kutta substep (1-indexed).
    :return: A list of sympy expressions representing the RHS output expressions.
    """
    s = rk_substep - 1  # Convert to 0-indexed
    rhs_output_expr = []

    if is_diagonal_Butcher(Butcher_dict, MoL_method) and "RK3" in MoL_method:
        y_n_gfs = "Y_N_GFS"
        k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs = (
            "K1_OR_Y_NPLUS_A21_K1_OR_Y_NPLUS1_RUNNING_TOTAL_GFS"
        )
        k2_or_y_nplus_a32_k2_gfs = "K2_OR_Y_NPLUS_A32_K2_GFS"

        if s == 0:
            rhs_output_expr = [k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs]
        elif s == 1:
            rhs_output_expr = [k2_or_y_nplus_a32_k2_gfs]
        elif s == 2:
            rhs_output_expr = [y_n_gfs]
    else:
        if not is_diagonal_Butcher(Butcher_dict, MoL_method):
            rhs_output_expr = [f"K{rk_substep}_GFS"]
        else:
            if MoL_method == "Euler":
                rhs_output_expr = ["Y_NPLUS1_RUNNING_TOTAL_GFS"]
            else:
                if s == 0:
                    rhs_output_expr = ["K_ODD_GFS"]
                # For the remaining steps the inputs and ouputs alternate between k_odd and k_even
                elif s % 2 == 0:
                    rhs_output_expr = ["K_ODD_GFS"]
                else:
                    rhs_output_expr = ["K_EVEN_GFS"]

    return (
        rhs_output_expr  # Return the list of strings representing the RHS output expr
    )


def register_CFunction_MoL_malloc_diagnostic_gfs() -> None:
    """
    Register the CFunction 'MoL_malloc_diagnostic_gfs'.
    This function allocates memory for diagnostic grid functions.
    :raises RuntimeError: If an error occurs while registering the CFunction
    :return None
    """
    includes: List[str] = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc: str = "Allocate memory for diagnostic gfs"
    cfunc_type: str = "void"
    name: str = "MoL_malloc_diagnostic_gfs"
    params: str = (
        "const commondata_struct *restrict commondata, const params_struct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs"
    )
    body: str = """
const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
gridfuncs->diagnostic_output_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
"""
    try:
        cfc.register_CFunction(
            includes=includes,
            desc=desc,
            cfunc_type=cfunc_type,
            name=name,
            params=params,
            include_CodeParameters_h=True,
            body=body,
        )
    except Exception as e:
        raise RuntimeError(
            f"Error registering CFunction 'MoL_malloc_diagnostic_gfs': {str(e)}"
        ) from e


def register_CFunction_MoL_free_memory_diagnostic_gfs() -> None:
    """
    Register the CFunction 'MoL_free_memory_diagnostic_gfs'.

    This function frees the memory allocated for diagnostic grid functions.

    :raises RuntimeError: If an error occurs while registering the CFunction

    :return None
    """
    includes: List[str] = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc: str = "Free memory for diagnostic gfs"
    cfunc_type: str = "void"
    name: str = "MoL_free_memory_diagnostic_gfs"
    params: str = "MoL_gridfunctions_struct *restrict gridfuncs"
    body: str = """
  free(gridfuncs->diagnostic_output_gfs);
"""
    try:
        cfc.register_CFunction(
            includes=includes,
            desc=desc,
            cfunc_type=cfunc_type,
            name=name,
            params=params,
            body=body,
        )
    except Exception as e:
        raise RuntimeError(
            f"Error registering CFunction 'MoL_free_memory_diagnostic_gfs': {str(e)}"
        ) from e


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
    enable_simd: bool = False,
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
    :param enable_simd: Flag to enable SIMD functionality.

    :return None

    Doctest:
    # FIXME
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if enable_simd:
        includes += [os.path.join("SIMD", "SIMD_intrinsics.h")]

    desc = f'Method of Lines (MoL) for "{MoL_method}" method: Step forward one full timestep.\n'
    cfunc_type = "void"
    name = "MoL_step_forward_in_time"
    params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata, const REAL time_start, const int which_RK_substep, const int which_MOL_part"

    # Code body
    body = f"""
// C code implementation of -={{ {MoL_method} }}=- Method of Lines timestepping.


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
REAL *restrict {y_n_gridfunctions} = {gf_prefix}{y_n_gridfunctions};
// Temporary timelevel & AUXEVOL gridfunctions:\n"""
    for gf in non_y_n_gridfunctions_list:
        gf_aliases += f"REAL *restrict {gf} = {gf_prefix}{gf};\n"

    gf_aliases += "params_struct *restrict params = &griddata[grid].params;\n"
    if enable_rfm_precompute:
        gf_aliases += (
            "const rfm_struct *restrict rfmstruct = &griddata[grid].rfmstruct;\n"
        )
    else:
        gf_aliases += "REAL *restrict xx[3]; for(int ww=0;ww<3;ww++) xx[ww] = griddata[grid].xx[ww];\n"

    if enable_curviBCs:
        gf_aliases += "const bc_struct *restrict bcstruct = &griddata[grid].bcstruct;\n"
    for i in ["0", "1", "2"]:
        gf_aliases += f"const int Nxx_plus_2NGHOSTS{i} = griddata[grid].params.Nxx_plus_2NGHOSTS{i};\n"

    # Implement Method of Lines (MoL) Timestepping
    Butcher = Butcher_dict[MoL_method][
        0
    ]  # Get the desired Butcher table from the dictionary
    num_steps = (
        len(Butcher) - 1
    )  # Specify the number of required steps to update solution

    dt = sp.Symbol("commondata->dt", real=True)

    body += """
  switch (which_RK_substep) {
"""

    if is_diagonal_Butcher(Butcher_dict, MoL_method) and "RK3" in MoL_method:
        # Diagonal RK3 only!!!
        #  In a diagonal RK3 method, only 3 gridfunctions need be defined. Below implements this approach.
        y_n_gfs = sp.Symbol("y_n_gfsL", real=True)
        k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs = sp.Symbol(
            "k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfsL", real=True
        )
        k2_or_y_nplus_a32_k2_gfs = sp.Symbol("k2_or_y_nplus_a32_k2_gfsL", real=True)

        # k_1
        body += r"""
          case RK_SUBSTEP_K1:
          {
"""

        body += """
// In a diagonal RK3 method like this one, only 3 gridfunctions need be defined. Below implements this approach.
// Using y_n_gfs as input, k1 and apply boundary conditions\n"""
        body += (
            single_RK_substep_input_symbolic(
                comment_block="""// -={ START k1 substep }=-
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
                post_rhs_list=[post_rhs_string],
                post_rhs_output_list=[
                    k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs
                ],
                enable_simd=enable_simd,
                gf_aliases=gf_aliases,
                post_post_rhs_string=post_post_rhs_string,
            )
            + "// -={ END k1 substep }=-\n\n"
        )
        body += """
        break;
        }
"""

        body += r"""
          case RK_SUBSTEP_K2:
          {
"""
        body += (
            single_RK_substep_input_symbolic(
                comment_block="""// -={ START k2 substep }=-
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
                post_rhs_list=[post_rhs_string, post_rhs_string],
                post_rhs_output_list=[
                    k2_or_y_nplus_a32_k2_gfs,
                    k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs,
                ],
                enable_simd=enable_simd,
                gf_aliases=gf_aliases,
                post_post_rhs_string=post_post_rhs_string,
            )
            + "// -={ END k2 substep }=-\n\n"
        )
        body += """
        break;
        }
"""

        body += r"""
          case RK_SUBSTEP_K3:
          {
"""

        body += (
            single_RK_substep_input_symbolic(
                comment_block="""// -={ START k3 substep }=-
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
                post_rhs_list=[post_rhs_string],
                post_rhs_output_list=[y_n_gfs],
                enable_simd=enable_simd,
                gf_aliases=gf_aliases,
                post_post_rhs_string=post_post_rhs_string,
            )
            + "// -={ END k3 substep }=-\n\n"
        )
        body += """
        break;
        }
"""

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

                # ~ body += f"""
                # ~ case RK_SUBSTEP_K{str(s + 1)}:
                # ~ """
                body += f"""
          case RK_SUBSTEP_K{str(s + 1)}:
          {{
"""

                body += f"""{single_RK_substep_input_symbolic(
                    comment_block=f"// -={{ START k{str(s + 1)} substep }}=-",
                    substep_time_offset_dt=Butcher[s][0],
                    rhs_str=rhs_string,
                    rhs_input_expr=rhs_input,
                    rhs_output_expr=rhs_output,
                    RK_lhs_list=[RK_lhs],
                    RK_rhs_list=[RK_rhs],
                    post_rhs_list=[post_rhs],
                    post_rhs_output_list=[post_rhs_output],
                    enable_simd=enable_simd,
                    gf_aliases=gf_aliases,
                    post_post_rhs_string=post_post_rhs_string,
                )}// -={{ END k{str(s + 1)} substep }}=-\n\n"""
                body += """
          break;
          }"""

        else:
            y_n = sp.Symbol("y_n_gfsL", real=True)
            y_nplus1_running_total = sp.Symbol("y_nplus1_running_total_gfsL", real=True)
            if (
                MoL_method == "Euler"
            ):  # Euler's method doesn't require any k_i, and gets its own unique algorithm
                body += single_RK_substep_input_symbolic(
                    comment_block="// ***Euler timestepping only requires one RHS evaluation***",
                    substep_time_offset_dt=Butcher[0][0],
                    rhs_str=rhs_string,
                    rhs_input_expr=y_n,
                    rhs_output_expr=y_nplus1_running_total,
                    RK_lhs_list=[y_n],
                    RK_rhs_list=[y_n + y_nplus1_running_total * dt],
                    post_rhs_list=[post_rhs_string],
                    post_rhs_output_list=[y_n],
                    enable_simd=enable_simd,
                    gf_aliases=gf_aliases,
                    post_post_rhs_string=post_post_rhs_string,
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
                    body += f"""
          case RK_SUBSTEP_K{str(s + 1)}:
          {{
"""
                    body += (
                        single_RK_substep_input_symbolic(
                            comment_block=f"// -={{ START k{s + 1} substep }}=-",
                            substep_time_offset_dt=Butcher[s][0],
                            rhs_str=rhs_string,
                            rhs_input_expr=rhs_input,
                            rhs_output_expr=rhs_output,
                            RK_lhs_list=RK_lhs_list,
                            RK_rhs_list=RK_rhs_list,
                            post_rhs_list=[post_rhs_string],
                            post_rhs_output_list=[post_rhs_output],
                            enable_simd=enable_simd,
                            gf_aliases=gf_aliases,
                            post_post_rhs_string=post_post_rhs_string,
                        )
                        + f"// -={{ END k{s + 1} substep }}=-\n\n"
                    )
                    body += """
          break;
          }"""

    body += """
  }
"""

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


def create_gf_constants(gf_list: List[str]) -> str:
    """
    Create C preprocessor constant definitions from a list of grid functions.

    :param gf_list: A list of grid function names.
    :return: A string with C preprocessor #define statements for each grid function.
    """
    return "\n".join(f"#define {gf.upper()} {i}" for i, gf in enumerate(gf_list))


def create_rk_substep_constants(num_steps: int) -> str:
    """
    Create C preprocessor constant definitions for Runge-Kutta substeps.

    :param num_steps: The number of Runge-Kutta substeps.
    :return: A string with C preprocessor #define statements for each substep.
    """
    return "\n".join(f"#define RK_SUBSTEP_K{s+1} {s+1}" for s in range(num_steps))


# Register all the CFunctions and NRPy basic defines
def register_CFunctions(
    MoL_method: str = "RK4",
    rhs_string: str = "rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, RK_INPUT_GFS, RK_OUTPUT_GFS);",
    post_rhs_string: str = "apply_bcs(Nxx, Nxx_plus_2NGHOSTS, RK_OUTPUT_GFS);",
    post_post_rhs_string: str = "",
    enable_rfm_precompute: bool = False,
    enable_curviBCs: bool = False,
    enable_simd: bool = False,
    register_MoL_step_forward_in_time: bool = True,
) -> None:
    r"""
    Register all MoL C functions and NRPy basic defines.

    :param MoL_method: The method to be used for MoL. Default is 'RK4'.
    :param rhs_string: RHS function call as string. Default is "rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, RK_INPUT_GFS, RK_OUTPUT_GFS);"
    :param post_rhs_string: Post-RHS function call as string. Default is "apply_bcs(Nxx, Nxx_plus_2NGHOSTS, RK_OUTPUT_GFS);"
    :param post_post_rhs_string: Post-post-RHS function call as string. Default is an empty string.
    :param enable_rfm_precompute: Enable reference metric support. Default is False.
    :param enable_curviBCs: Enable curvilinear boundary conditions. Default is False.
    :param enable_simd: Enable Single Instruction, Multiple Data (SIMD). Default is False.
    :param register_MoL_step_forward_in_time: Whether to register the MoL step forward function. Default is True.

    :return None

    Doctests:
    >>> from nrpy.helpers.generic import compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunctions()
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4B/RA7ddABGaScZHDxOiAHcc/CXr2duHb+UiUyv83OVALtvJ+o7uK/PoSVGe7rPuvil8asOnzsX/43MrS1REEi/tau4rRkS3klwMCWne6D351BIv83jxwuBwBgfb9aLOiuMaxdzlpat7M5Zzy6cqD3qxMNABQOc2xVV5NC/sFWryHJK7NLtTQZSJAkfrM9dF6qg6pG5p6oN+o9MOcVuOHCVrZ0lCxYx6wuKz2IJ/mMdvxb9kpOoc+n71ZJsMV7tA14+9i8TawSx62Kef1R0clKDrO9YH+vibd56jTMWlJJj7qwA4ejblO0802o+9UA00dZsIZRIq9k/LjwpM473QxNNtZt03nkrs/BShBdq2ZyQpIOcw4mZH+AJkh40wl5nKuTRIlLFG4Q6NBh276EAjeZ2AuRBWcD9KjPiSnm6E53GBzACBbVbUxEVCiFzWJsjBwNoZj+o2133GqXbzjLmvBUd7vlc8BVGk7x7yQAmajdQHpOKiPAu1g1Ch753e50ASUNrK27jpZmaGBGoIomk+YaSm9q5qS+kEx597LllKKlVnGRwinxSiP2AGfto3ST6yNHdH81ISIrL/6a3I7xZbc5mqdhbzNABxa3Z5BMbWcIHHuS3FVQ06uehvm7rSoMNgfYFRVPZSR6oO1KdJkTKVK1HF1Di7riialwoMpPoP9hkaRNko3+VQsXDfD8zZD6wI3dyRjWxc7mvzHBwx0cY3LMQdzaJUI5wPMv+zt67Ss3boY3mmwDorWTZPgtpO1oIklYhzAPfebO0fc+Za54sKFZl2yiHAsfuQAcLorBWaWP0LvYxFo6Xd0JOxuCjnuDMO+UIBcHHjtaxxErecmOOKXpIHLh7gQVgj/S7iSe1t8omN0Hf5S6JxQQuBAX9xr5LyU6E5INDGSciZhysRx+RWycRTsZ+QYOWVDM2oYdjcKxz8DXe7E24kj2bL7y17VdoQa76t8IALsAt8p2krjxFZsRwiCk9o/xmuCoNLDNsMWXktDEOEv5Rvra2g6gxG7enKHHQkXggrnH8gRIPuk6miLu3g8v9Maua5um7sJGBXaHr9QjfEN/pToqqo+vLQnBqP/vukzhmNo7GkL9x5cPiZ5pbngo9cm1+bad3DG0EwUmGCmeFRIbPZpOcRoM9+RDz+uKcGc5ICRRe0H7yvOxhpwK7YU1WLWlfUA7NO56FenA6xc54N6mW+Ot/tRgF6w4kzANqpGGwkRoP0g2h+ufyLNzmzE50KK0kib3HI0fqOBnNUFp0dF6Df57M9OaQAAAAAA6HVu59NvyUAAdMH0j8AAPO1tUaxxGf7AgAAAAAEWVo=")
    >>> returned_string = cfc.CFunction_dict["MoL_step_forward_in_time"].full_function
    >>> if returned_string != expected_string:
    ...    compressed_str = compress_string_to_base64(returned_string)
    ...    error_message = "Trusted MoL_step_forward_in_time.full_function string changed!\n Here's the diff:\n"
    ...    error_message += "Here's the diff:\n" + diff_strings(expected_string, returned_string) + "\n"
    ...    raise ValueError(error_message + f"base64-encoded output: {compressed_str}")
    """
    Butcher_dict = generate_Butcher_tables()

    for which_gfs in ["y_n_gfs", "non_y_n_gfs"]:
        register_CFunction_MoL_malloc(Butcher_dict, MoL_method, which_gfs)
        register_CFunction_MoL_free_memory(Butcher_dict, MoL_method, which_gfs)

    register_CFunction_MoL_malloc_diagnostic_gfs()
    register_CFunction_MoL_free_memory_diagnostic_gfs()

    if register_MoL_step_forward_in_time:
        register_CFunction_MoL_step_forward_in_time(
            Butcher_dict,
            MoL_method,
            rhs_string,
            post_rhs_string,
            post_post_rhs_string,
            enable_rfm_precompute=enable_rfm_precompute,
            enable_curviBCs=enable_curviBCs,
            enable_simd=enable_simd,
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

    # Convert y_n_gridfunctions to a list if it's a string
    gf_list = (
        [y_n_gridfunctions] if isinstance(y_n_gridfunctions, str) else y_n_gridfunctions
    )
    gf_list.extend(non_y_n_gridfunctions_list)

    gf_constants = create_gf_constants(gf_list)

    rk_substep_constants = create_rk_substep_constants(
        len(Butcher_dict[MoL_method][0]) - 1
    )

    # Step 3.b: Create MoL_timestepping struct:
    BHaH_defines_h.register_BHaH_defines(
        "nrpy.infrastructures.BHaH.MoLtimestepping.MoL",
        f"typedef struct __MoL_gridfunctions_struct__ {{\n"
        f"REAL *restrict {y_n_gridfunctions};\n"
        + "".join(f"REAL *restrict {gfs};\n" for gfs in non_y_n_gridfunctions_list)
        + r"""REAL *restrict diagnostic_output_gfs;
REAL *restrict diagnostic_output_gfs2;
} MoL_gridfunctions_struct;
"""
        + rf"""// Define constants for accessing gridfunction types
{gf_constants}"""
        + r"""
#define LOOP_ALL_GFS_GPS(ii) \
_Pragma("omp parallel for") \
  for(int (ii)=0;(ii)<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;(ii)++)
"""
        + rf"""// Define constants for rk substeps
{rk_substep_constants}""",
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
