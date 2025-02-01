# nrpy/infrastructures/BHaH/MoLtimestepping/MoL_allocators.py
"""
Memory allocators and deallocators for y_n_gfs and non_y_n_gfs in MoL-based codes.

Authors: Zachariah B. Etienne (lead maintainer)
         zachetie **at** gmail **dot* com
         Samuel Tootle (GPU support in NRPy2)
         Brandon Clark (original, NRPy1 version)
"""

from typing import Dict, List, Tuple, Union

import nrpy.c_function as cfc
from nrpy.infrastructures.BHaH.MoLtimestepping.MoL_gridfunction_names import (
    generate_gridfunction_names,
)
from nrpy.infrastructures.BHaH.MoLtimestepping.MoL_rk_substep import (
    check_supported_parallelization,
)


def register_CFunction_MoL_malloc(
    Butcher_dict: Dict[str, Tuple[List[List[Union[int, str]]], int]],
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
    """
    check_supported_parallelization("register_CFunction_MoL_malloc", parallelization)
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    (
        y_n_gridfunctions,
        non_y_n_gridfunctions_list,
        diagnostic_gridfunctions_point_to,
        diagnostic_gridfunctions2_point_to,
    ) = generate_gridfunction_names(Butcher_dict, MoL_method=MoL_method)

    if which_gfs == "y_n_gfs":
        gridfunctions_list = [y_n_gridfunctions]
    elif which_gfs == "non_y_n_gfs":
        gridfunctions_list = non_y_n_gridfunctions_list
    else:
        raise ValueError(f'ERROR: which_gfs = "{which_gfs}" unrecognized.')

    desc = f"""Method of Lines (MoL) for "{MoL_method}" method: Allocate memory for "{which_gfs}" gridfunctions
   - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep
   - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method"""
    cfunc_type = "void"
    name = f"MoL_malloc_{which_gfs}"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs"

    body = "const int Nxx_plus_2NGHOSTS_tot = params->Nxx_plus_2NGHOSTS0 * params->Nxx_plus_2NGHOSTS1 * params->Nxx_plus_2NGHOSTS2;\n"
    for gridfunctions in gridfunctions_list:
        num_gfs = (
            "NUM_EVOL_GFS" if gridfunctions != "auxevol_gfs" else "NUM_AUXEVOL_GFS"
        )
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
        include_CodeParameters_h=False,
        body=body,
    )


def register_CFunction_MoL_free_memory(
    Butcher_dict: Dict[str, Tuple[List[List[Union[int, str]]], int]],
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
    check_supported_parallelization(
        "register_CFunction_MoL_free_memory", parallelization
    )
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f'Method of Lines (MoL) for "{MoL_method}" method: Free memory for "{which_gfs}" gridfunctions\n'
    desc += "   - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep\n"
    desc += "   - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method\n"
    cfunc_type = "void"

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
        if gridfunction == "auxevol_gfs":
            body += "  if(NUM_AUXEVOL_GFS > 0)"
        if parallelization == "cuda":
            body += f" cudaFree(gridfuncs->{gridfunction});\n"
        else:
            body += f" free(gridfuncs->{gridfunction});\n"

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
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
