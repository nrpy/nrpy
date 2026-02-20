# nrpy/infrastructures/BHaH/MoLtimestepping/MoL_free_intermediate_stage_gfs.py
"""
C function registration for freeing Method of Lines (MoL) intermediate-level storage.

Authors: Zachariah B. Etienne (lead maintainer)
         zachetie **at** gmail **dot* com
         Samuel Tootle (GPU support in NRPy2)
         Brandon Clark (original, NRPy1 version)
"""

from typing import Dict, List, Tuple, Union

import sympy as sp

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures import BHaH


def register_CFunction_MoL_free_intermediate_levels(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
) -> None:
    """
    Construct and register a C function that frees intermediate-level (k_i) storage for the chosen Method of Lines scheme.

    This function inspects the provided Butcher table and MoL method to determine which
    intermediate-level gridfunction arrays must be freed (with "auxevol_gfs" explicitly
    excluded). It then registers a C function named "MoL_malloc_intermediate_stage_gfs" that,
    frees intermediate-level gridfunction storage.

    :param Butcher_dict: Dictionary mapping MoL scheme names to their Butcher tables and
                         number of intermediate stages.
    :param MoL_method: Name of the Method of Lines scheme whose intermediate arrays should
                       be freed (for example, "RK4", "RK3", or "RK2").

    Doctests:
    TBD
    """
    parallelization = par.parval_from_str("parallelization")
    BHaH.MoLtimestepping.rk_substep.check_supported_parallelization(
        "register_CFunction_MoL_free_intermediate_stage_gfs"
    )
    includes = ["BHaH_defines.h"]
    desc = f"""
 * Free intermediate-level (k_i) storage for the "{MoL_method}" Method of Lines (MoL) scheme.
 *
 * This routine is registered as "MoL_free_intermediate_stage_gfs".
 * It frees intermediate-level gridfunction storage.
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
    name = "MoL_free_intermediate_stage_gfs"
    params = "MoL_gridfunctions_struct *restrict gridfuncs"
    intermediate_stage_gfs_gridfunctions_list = BHaH.MoLtimestepping.rk_butcher_table_dictionary.intermediate_stage_gf_names_list(
        Butcher_dict, MoL_method=MoL_method
    )
    body = ""
    free_macro = "BHAH_FREE_DEVICE" if parallelization == "cuda" else "BHAH_FREE"
    for gridfunctions in intermediate_stage_gfs_gridfunctions_list:
        body += f"{free_macro}(gridfuncs->{gridfunctions});\n"

    cfc.register_CFunction(
        subdirectory="MoL",
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
