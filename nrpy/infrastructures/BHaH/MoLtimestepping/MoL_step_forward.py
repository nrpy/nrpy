# nrpy/infrastructures/BHaH/MoLtimestepping/MoL_step_forward.py
"""
Register the core function that performs one full MoL time step, assembling substeps as needed.

Authors: Zachariah B. Etienne (lead maintainer)
         zachetie **at** gmail **dot* com
         Samuel Tootle (GPU support in NRPy2)
         Brandon Clark (original, NRPy1 version)
"""

import os
import warnings
from typing import Dict, List, Tuple, Union

import sympy as sp

import nrpy.c_function as cfc
from nrpy.infrastructures.BHaH.MoLtimestepping.MoL_gridfunction_names import (
    generate_gridfunction_names,
    is_diagonal_Butcher,
)
from nrpy.infrastructures.BHaH.MoLtimestepping.MoL_rk_substep import (
    check_supported_parallelization,
    construct_RK_functions_prefunc,
    single_RK_substep_input_symbolic,
)


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
    r"""
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
    >>> from nrpy.infrastructures.BHaH.MoLtimestepping.MoL_step_forward import register_CFunction_MoL_step_forward_in_time
    >>> from nrpy.infrastructures.BHaH.MoLtimestepping.MoL_rk_substep import MoL_Functions_dict
    >>> from nrpy.helpers.generic import validate_strings
    >>> from nrpy.infrastructures.BHaH.MoLtimestepping.RK_Butcher_Table_Dictionary import generate_Butcher_tables
    >>> Butcher_dict = generate_Butcher_tables()
    >>> rhs_string = "rhs_eval(commondata, params, rfmstruct,  auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);"
    >>> post_rhs_string = (
    ... "if (strncmp(commondata->outer_bc_type, \"extrapolation\", 50) == 0)\n"
    ... "    apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);"
    ... )
    >>> # We'll run a quick check for each known Butcher table that is non-adaptive:
    >>> for k, v in Butcher_dict.items():
    ...     Butcher = Butcher_dict[k][0]
    ...     cfc.CFunction_dict.clear()
    ...     MoL_Functions_dict.clear()
    ...     if Butcher[-1][0] != "":
    ...         continue  # skip adaptive methods
    ...     register_CFunction_MoL_step_forward_in_time(
    ...         Butcher_dict,
    ...         k,
    ...         rhs_string=rhs_string,
    ...         post_rhs_string=post_rhs_string,
    ...         enable_intrinsics=True,
    ...         parallelization="cuda",
    ...         rational_const_alias="static constexpr"
    ...     )
    ...     generated_str = cfc.CFunction_dict["MoL_step_forward_in_time"].full_function
    ...     validation_desc = f"CUDA__MoL_step_forward_in_time__{k}".replace(" ", "_")
    ...     validate_strings(generated_str, validation_desc, file_ext="cu")
    >>> cfc.CFunction_dict.clear()
    >>> MoL_Functions_dict.clear()
    >>> try:
    ...     register_CFunction_MoL_step_forward_in_time(Butcher_dict, "AHE")
    ... except ValueError as e:
    ...     print(f"ValueError: {e.args[0]}")
    ValueError: Adaptive order Butcher tables are currently not supported in MoL.
    """
    check_supported_parallelization(
        "register_CFunction_MoL_step_forward_in_time", parallelization
    )
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

    (
        y_n_gridfunctions,
        non_y_n_gridfunctions_list,
        _discard_pt,
        _discard_pt2,
    ) = generate_gridfunction_names(Butcher_dict, MoL_method)

    gf_prefix = "griddata[grid].gridfuncs."

    gf_aliases = f"""// Set gridfunction aliases, from griddata[].gridfuncs.
MAYBE_UNUSED REAL *restrict {y_n_gridfunctions} = {gf_prefix}{y_n_gridfunctions};
"""
    for gf in non_y_n_gridfunctions_list:
        gf_aliases += f"MAYBE_UNUSED REAL *restrict {gf} = {gf_prefix}{gf};\n"
    gf_aliases += """// Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
"""
    if enable_rfm_precompute:
        gf_aliases += "MAYBE_UNUSED const rfm_struct *restrict rfmstruct = griddata[grid].rfmstruct;\n"
    else:
        gf_aliases += "MAYBE_UNUSED REAL *restrict xx[3]; for(int ww=0;ww<3;ww++) xx[ww] = griddata[grid].xx[ww];\n"
    if enable_curviBCs:
        gf_aliases += "MAYBE_UNUSED const bc_struct *restrict bcstruct = &griddata[grid].bcstruct;\n"

    body = f"""
// C code implementation of -={{ {MoL_method} }}=- Method of Lines timestepping.

// First set the initial time:
const REAL time_start = commondata->time;
"""

    Butcher = Butcher_dict[MoL_method][0]
    if Butcher[-1][0] != "":
        raise ValueError(
            "Adaptive order Butcher tables are currently not supported in MoL."
        )
    num_steps = len(Butcher) - 1

    rk_step_body_dict: Dict[str, str] = {}

    # Diagonal RK3 check
    if is_diagonal_Butcher(Butcher_dict, MoL_method) and "RK3" in MoL_method:
        body += """
// In a diagonal RK3 method, we only need 3 gridfunctions (y_n, k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs, k2_or_y_nplus_a32_k2_gfs).
"""

        y_n_gfs = sp.Symbol("y_n_gfsL", real=True)
        k1_or_yrun = sp.Symbol(
            "k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfsL", real=True
        )
        k2_or_yrun = sp.Symbol("k2_or_y_nplus_a32_k2_gfsL", real=True)

        # k1
        rk_step_body_dict["RK_SUBSTEP_K1"] = (
            single_RK_substep_input_symbolic(
                additional_comments="""
// RHS evaluation:
//  1. We will store k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs
// Post-RHS evaluation: apply boundaries""",
                substep_time_offset_dt=Butcher[0][0],
                rhs_str=rhs_string,
                rhs_input_expr=y_n_gfs,
                rhs_output_expr=k1_or_yrun,
                RK_lhs_list=[k1_or_yrun],
                RK_rhs_list=[
                    Butcher[1][1] * k1_or_yrun * sp.Symbol("commondata->dt", real=True)
                    + y_n_gfs
                ],
                rk_step=1,
                post_rhs_list=[post_rhs_string],
                post_rhs_output_list=[k1_or_yrun],
                enable_intrinsics=enable_intrinsics,
                gf_aliases=gf_aliases,
                post_post_rhs_string=post_post_rhs_string,
                rational_const_alias=rational_const_alias,
                parallelization=parallelization,
            )
            + "// -={ END k1 substep }=-\n\n"
        )

        # k2
        rk_step_body_dict["RK_SUBSTEP_K2"] = (
            single_RK_substep_input_symbolic(
                additional_comments="""
// RHS evaluation:
//    1. Reassign k1_or_yrun to the running total
// Post-RHS evaluation:
//    1. Apply post-RHS""",
                substep_time_offset_dt=Butcher[1][0],
                rhs_str=rhs_string,
                rhs_input_expr=k1_or_yrun,
                rhs_output_expr=k2_or_yrun,
                RK_lhs_list=[k1_or_yrun, k2_or_yrun],
                RK_rhs_list=[
                    Butcher[3][1] * (k1_or_yrun - y_n_gfs) / Butcher[1][1]
                    + y_n_gfs
                    + Butcher[3][2]
                    * k2_or_yrun
                    * sp.Symbol("commondata->dt", real=True),
                    Butcher[2][2] * k2_or_yrun * sp.Symbol("commondata->dt", real=True)
                    + y_n_gfs,
                ],
                rk_step=2,
                post_rhs_list=[post_rhs_string, post_rhs_string],
                post_rhs_output_list=[k2_or_yrun, k1_or_yrun],
                enable_intrinsics=enable_intrinsics,
                gf_aliases=gf_aliases,
                post_post_rhs_string=post_post_rhs_string,
                rational_const_alias=rational_const_alias,
                parallelization=parallelization,
            )
            + "// -={ END k2 substep }=-\n\n"
        )

        # k3
        rk_step_body_dict["RK_SUBSTEP_K3"] = (
            single_RK_substep_input_symbolic(
                additional_comments="""
// RHS evaluation:
//    1. Add k3 to the running total and save to y_n
// Post-RHS evaluation:
//    1. Apply post-RHS to y_n""",
                substep_time_offset_dt=Butcher[2][0],
                rhs_str=rhs_string,
                rhs_input_expr=k2_or_yrun,
                rhs_output_expr=y_n_gfs,
                RK_lhs_list=[y_n_gfs],
                RK_rhs_list=[
                    k1_or_yrun
                    + Butcher[3][3] * y_n_gfs * sp.Symbol("commondata->dt", real=True)
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
        # Non-diagonal or other diagonal
        y_n = sp.Symbol("y_n_gfsL", real=True)
        if not is_diagonal_Butcher(Butcher_dict, MoL_method):
            # Non-diagonal
            for s in range(num_steps):
                next_y_input = sp.Symbol("next_y_input_gfsL", real=True)
                if s == 0:
                    rhs_input = y_n
                else:
                    rhs_input = next_y_input

                rhs_output = sp.Symbol(f"k{s + 1}_gfsL", real=True)

                if s == num_steps - 1:
                    RK_lhs = y_n
                else:
                    RK_lhs = next_y_input

                RK_rhs = y_n
                for m in range(s + 1):
                    k_mp1_gfs = sp.Symbol("k" + str(m + 1) + "_gfsL")
                    if Butcher[s + 1][m + 1] != 0:
                        if Butcher[s + 1][m + 1] != 1:
                            RK_rhs += (
                                sp.Symbol("commondata->dt", real=True)
                                * k_mp1_gfs
                                * Butcher[s + 1][m + 1]
                            )
                        else:
                            RK_rhs += sp.Symbol("commondata->dt", real=True) * k_mp1_gfs

                if s == num_steps - 1:
                    post_rhs_output = y_n
                else:
                    post_rhs_output = next_y_input

                rk_step_body_dict[f"RK_SUBSTEP_K{s+1}"] = (
                    single_RK_substep_input_symbolic(
                        substep_time_offset_dt=Butcher[s][0],
                        rhs_str=rhs_string,
                        rhs_input_expr=rhs_input,
                        rhs_output_expr=rhs_output,
                        RK_lhs_list=[RK_lhs],
                        RK_rhs_list=[RK_rhs],
                        rk_step=s + 1,
                        post_rhs_list=[post_rhs_string],
                        post_rhs_output_list=[post_rhs_output],
                        enable_intrinsics=enable_intrinsics,
                        gf_aliases=gf_aliases,
                        post_post_rhs_string=post_post_rhs_string,
                        rational_const_alias=rational_const_alias,
                        parallelization=parallelization,
                    )
                    + f"// -={{ END k{str(s + 1)} substep }}=-\n\n"
                )
        else:
            # Diagonal Butcher, e.g., Euler or standard diagonal RK4
            y_nplus1_running_total = sp.Symbol("y_nplus1_running_total_gfsL", real=True)
            if MoL_method == "Euler":
                rk_step_body_dict[f"{MoL_method}"] = single_RK_substep_input_symbolic(
                    additional_comments="// ***Euler timestepping only requires one RHS evaluation***",
                    substep_time_offset_dt=Butcher[0][0],
                    rhs_str=rhs_string,
                    rhs_input_expr=y_n,
                    rhs_output_expr=y_nplus1_running_total,
                    RK_lhs_list=[y_n],
                    RK_rhs_list=[
                        y_n
                        + y_nplus1_running_total
                        * sp.Symbol("commondata->dt", real=True)
                    ],
                    post_rhs_list=[post_rhs_string],
                    post_rhs_output_list=[y_n],
                    enable_intrinsics=enable_intrinsics,
                    gf_aliases=gf_aliases,
                    post_post_rhs_string=post_post_rhs_string,
                    rational_const_alias=rational_const_alias,
                    parallelization=parallelization,
                    rk_step=None,
                )
            else:
                k_odd = sp.Symbol("k_odd_gfsL", real=True)
                k_even = sp.Symbol("k_even_gfsL", real=True)
                for s in range(num_steps):
                    if s == 0:
                        rhs_input = y_n
                        rhs_output = k_odd
                    elif s % 2 == 0:
                        rhs_input = k_even
                        rhs_output = k_odd
                    else:
                        rhs_input = k_odd
                        rhs_output = k_even

                    RK_lhs_list = []
                    RK_rhs_list = []
                    if s != num_steps - 1:
                        if s == 0:
                            RK_lhs_list.append(y_nplus1_running_total)
                            RK_rhs_list.append(
                                rhs_output
                                * sp.Symbol("commondata->dt", real=True)
                                * Butcher[num_steps][s + 1]
                            )
                            RK_lhs_list.append(rhs_output)
                            RK_rhs_list.append(
                                y_n
                                + rhs_output
                                * sp.Symbol("commondata->dt", real=True)
                                * Butcher[s + 1][s + 1]
                            )
                        else:
                            if Butcher[num_steps][s + 1] != 0:
                                RK_lhs_list.append(y_nplus1_running_total)
                                if Butcher[num_steps][s + 1] != 1:
                                    RK_rhs_list.append(
                                        y_nplus1_running_total
                                        + rhs_output
                                        * sp.Symbol("commondata->dt", real=True)
                                        * Butcher[num_steps][s + 1]
                                    )
                                else:
                                    RK_rhs_list.append(
                                        y_nplus1_running_total
                                        + rhs_output
                                        * sp.Symbol("commondata->dt", real=True)
                                    )
                            if Butcher[s + 1][s + 1] != 0:
                                RK_lhs_list.append(rhs_output)
                                if Butcher[s + 1][s + 1] != 1:
                                    RK_rhs_list.append(
                                        y_n
                                        + rhs_output
                                        * sp.Symbol("commondata->dt", real=True)
                                        * Butcher[s + 1][s + 1]
                                    )
                                else:
                                    RK_rhs_list.append(
                                        y_n
                                        + rhs_output
                                        * sp.Symbol("commondata->dt", real=True)
                                    )
                        post_rhs_output = rhs_output
                    else:
                        # final step
                        if Butcher[num_steps][s + 1] != 0:
                            RK_lhs_list.append(y_n)
                            if Butcher[num_steps][s + 1] != 1:
                                RK_rhs_list.append(
                                    y_n
                                    + y_nplus1_running_total
                                    + rhs_output
                                    * sp.Symbol("commondata->dt", real=True)
                                    * Butcher[num_steps][s + 1]
                                )
                            else:
                                RK_rhs_list.append(
                                    y_n
                                    + y_nplus1_running_total
                                    + rhs_output
                                    * sp.Symbol("commondata->dt", real=True)
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
        prefunc = r"""
#define LOOP_ALL_GFS_GPS(ii) \
const int tid0 = threadIdx.x + blockIdx.x*blockDim.x; \
const int stride0 = blockDim.x * gridDim.x; \
  for(int (ii)=(tid0);(ii)<d_params[streamid].Nxx_plus_2NGHOSTS0*d_params[streamid].Nxx_plus_2NGHOSTS1*d_params[streamid].Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;(ii)+=(stride0))
"""
    else:
        prefunc = r"""
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
// so here we set time based on the iteration number:
commondata->time = (REAL)(commondata->nn + 1) * commondata->dt;

// Increment the timestep n:
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


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
