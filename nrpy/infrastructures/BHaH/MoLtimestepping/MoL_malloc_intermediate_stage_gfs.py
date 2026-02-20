# nrpy/infrastructures/BHaH/MoLtimestepping/MoL_malloc_intermediate_stage_gfs.py
"""
C function registration for allocating Method of Lines (MoL) intermediate-level storage.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Samuel Tootle (GPU support in NRPy2)
        Brandon Clark (original, NRPy1 version)
"""

from typing import Dict, List, Tuple, Union

import sympy as sp

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures import BHaH


def register_CFunction_MoL_malloc_intermediate_stage_gfs(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
) -> None:
    """
    Construct and register a C function that allocates intermediate-level (k_i) storage for the chosen Method of Lines scheme.

    This function inspects the provided Butcher table and MoL method to determine which
    intermediate-level gridfunction arrays must be allocated (with "auxevol_gfs" explicitly
    excluded). It then registers a C function named "MoL_malloc_intermediate_stage_gfs" that,
    at runtime, computes:
        Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2
    and allocates, for each required intermediate-level array:
        sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot
    The allocation macro is chosen based on the NRPy "parallelization" parameter:
    "BHAH_MALLOC_DEVICE" if "cuda", else "BHAH_MALLOC". The registered function includes
    "BHaH_defines.h" and uses the common NRPy/BHaH data structures.

    :param Butcher_dict: Dictionary mapping MoL scheme names to their Butcher tables and
                         number of intermediate stages.
    :param MoL_method: Name of the Method of Lines scheme whose intermediate arrays should
                       be allocated (for example, "RK4", "RK3", or "RK2").

    Doctests:
    TBD
    """
    parallelization = par.parval_from_str("parallelization")
    BHaH.MoLtimestepping.rk_substep.check_supported_parallelization(
        "register_CFunction_MoL_malloc_intermediate_stage_gfs"
    )

    includes = ["BHaH_defines.h"]

    desc = f"""
 * Allocate intermediate-level (k_i) storage for the "{MoL_method}" Method of Lines (MoL) scheme.
 *
 * This routine is registered as "MoL_malloc_intermediate_stage_gfs". At runtime it computes
 * the total number of grid points including ghost zones as:
 *     Nxx_plus_2NGHOSTS_tot = params->Nxx_plus_2NGHOSTS0
 *                            * params->Nxx_plus_2NGHOSTS1
 *                            * params->Nxx_plus_2NGHOSTS2;
 * For each required intermediate-level gridfunction array (with "auxevol_gfs" excluded),
 * the routine allocates:
 *     sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot
 * The allocation macro is chosen by the build/runtime configuration:
 *     BHAH_MALLOC_DEVICE for CUDA builds, otherwise BHAH_MALLOC.
 *
 * @param commondata Pointer to common, read-only runtime data. Included for a uniform
 *                   function signature across BHaH routines; it is not modified here.
 * @param params     Pointer to grid parameter struct providing Nxx_plus_2NGHOSTS0/1/2 and
 *                   related metadata needed to compute allocation sizes.
 * @param gridfuncs  Pointer to the MoL gridfunctions struct whose intermediate-level
 *                   arrays will be allocated by this routine.
 *
 * @return void
"""
    cfunc_type = "void"
    name = "MoL_malloc_intermediate_stage_gfs"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs"

    body = "const int Nxx_plus_2NGHOSTS_tot = params->Nxx_plus_2NGHOSTS0 * params->Nxx_plus_2NGHOSTS1 * params->Nxx_plus_2NGHOSTS2;\n"
    intermediate_stage_gfs_gridfunctions_list = BHaH.MoLtimestepping.rk_butcher_table_dictionary.intermediate_stage_gf_names_list(
        Butcher_dict, MoL_method=MoL_method
    )
    allocator_macro = (
        "BHAH_MALLOC_DEVICE" if parallelization == "cuda" else "BHAH_MALLOC"
    )
    for gridfunctions in intermediate_stage_gfs_gridfunctions_list:
        body += f"{allocator_macro}(gridfuncs->{gridfunctions}, sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);\n"

    cfc.register_CFunction(
        subdirectory="MoL",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
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
