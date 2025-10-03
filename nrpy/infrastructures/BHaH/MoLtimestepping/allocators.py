# nrpy/infrastructures/BHaH/MoLtimestepping/allocators.py
"""
Generates C functions for memory management within the MoL time-stepping framework.

This module provides functions to create C code for allocating and deallocating
memory for the gridfunctions required by the Method of Lines (MoL) integrators.
It distinguishes between two categories of gridfunctions:

- `y_n_gfs`: The set of gridfunctions that store the state of the system at the
  beginning of a time step (t_n).
- `non_y_n_gfs`: Auxiliary gridfunctions needed for intermediate calculations
  within a time step, such as the stages (k_i) in a Runge-Kutta method.

The main functions are:
- `register_CFunction_MoL_malloc()`: Registers a C function to allocate memory.
- `register_CFunction_MoL_free_memory()`: Registers a C function to free memory.

These are specialized based on the chosen MoL algorithm (e.g., "RK4") to
ensure the correct number of auxiliary arrays are managed.

Authors: Zachariah B. Etienne (lead maintainer)
         zachetie **at** gmail **dot* com
         Samuel Tootle (GPU support in NRPy2)
         Brandon Clark (original, NRPy1 version)
"""

from typing import Dict, List, Tuple, Union

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures import BHaH


def register_CFunction_MoL_malloc(
    Butcher_dict: Dict[str, Tuple[List[List[Union[int, str]]], int]],
    MoL_method: str,
    which_gfs: str,
) -> None:
    """
    Register MoL_malloc_y_n_gfs() and MoL_malloc_non_y_n_gfs(), allocating memory for the gridfunctions indicated.

    :param Butcher_dict: Dictionary of Butcher tables for the MoL method.
    :param MoL_method: Method for the Method of Lines.
    :param which_gfs: Specifies which gridfunctions to consider ("y_n_gfs" or "non_y_n_gfs").

    :raises ValueError: If the which_gfs parameter is neither "y_n_gfs" nor "non_y_n_gfs".

    Doctest: FIXME
    """
    parallelization = par.parval_from_str("parallelization")
    BHaH.MoLtimestepping.rk_substep.check_supported_parallelization(
        "register_CFunction_MoL_malloc"
    )
    includes = ["BHaH_defines.h"]

    (
        y_n_gridfunctions,
        non_y_n_gridfunctions_list,
        diagnostic_gridfunctions_point_to,
        diagnostic_gridfunctions2_point_to,
    ) = BHaH.MoLtimestepping.gridfunction_names.generate_gridfunction_names(
        Butcher_dict, MoL_method=MoL_method
    )

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
        body += f"BHAH_MALLOC(gridfuncs->{gridfunctions}, sizeof(REAL) * {num_gfs} * Nxx_plus_2NGHOSTS_tot);\n".replace(
            "BHAH_MALLOC",
            "BHAH_MALLOC_DEVICE" if parallelization in ["cuda"] else "BHAH_MALLOC",
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
) -> None:
    """
    Free memory for the specified Method of Lines (MoL) gridfunctions, given an MoL_method.

    :param Butcher_dict: Dictionary containing Butcher tableau for MoL methods.
    :param MoL_method: The Method of Lines method.
    :param which_gfs: The gridfunctions to be freed, either 'y_n_gfs' or 'non_y_n_gfs'.

    :raises ValueError: If the 'which_gfs' argument is unrecognized.
    """
    parallelization = par.parval_from_str("parallelization")
    BHaH.MoLtimestepping.rk_substep.check_supported_parallelization(
        "register_CFunction_MoL_free_memory"
    )
    includes = ["BHaH_defines.h"]
    desc = f'Method of Lines (MoL) for "{MoL_method}" method: Free memory for "{which_gfs}" gridfunctions\n'
    desc += "   - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep\n"
    desc += "   - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method\n"
    cfunc_type = "void"

    (
        y_n_gridfunctions,
        non_y_n_gridfunctions_list,
        _diagnostic_gridfunctions_point_to,
        _diagnostic_gridfunctions2_point_to,
    ) = BHaH.MoLtimestepping.gridfunction_names.generate_gridfunction_names(
        Butcher_dict, MoL_method
    )

    if which_gfs == "y_n_gfs":
        gridfunctions_list = [y_n_gridfunctions]
    elif which_gfs == "non_y_n_gfs":
        gridfunctions_list = non_y_n_gridfunctions_list
    else:
        raise ValueError(f'ERROR: which_gfs = "{which_gfs}" unrecognized.')

    name: str = f"MoL_free_memory_{which_gfs}"
    params = "MoL_gridfunctions_struct *restrict gridfuncs"
    if which_gfs == "non_y_n_gfs":
        params += ", const bool free_auxevol_gfs_if_exist"
    body = ""
    for gridfunction in gridfunctions_list:
        if gridfunction == "auxevol_gfs":
            body += "  if(NUM_AUXEVOL_GFS > 0 && free_auxevol_gfs_if_exist)"
        if parallelization == "cuda":
            body += f" BHAH_FREE_DEVICE(gridfuncs->{gridfunction});\n"
        else:
            body += f" BHAH_FREE(gridfuncs->{gridfunction});\n"

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
