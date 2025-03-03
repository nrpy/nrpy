"""
C function for setting constant auxiliary terms for hyperbolic relaxation.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.simple_loop as lp
from nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms import (
    compute_psi_background_and_ADD_times_AUU,
)


# Define functions to set AUXEVOL gridfunctions
def register_CFunction_auxevol_gfs_single_point(
    CoordSystem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for the AUXEVOL grid functions at a single point.

    :param CoordSystem: The coordinate system to use in setting up the AUXEVOL gridfunctions.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Compute psi_background and ADD_times_AUU
    psi_background, ADD_times_AUU = compute_psi_background_and_ADD_times_AUU(
        CoordSystem
    )

    includes = ["BHaH_defines.h"]

    desc = r"""Compute AUXEVOL grid functions at a single point."""
    cfunc_type = "void"
    name = "auxevol_gfs_single_point"
    params = r"""const commondata_struct *restrict commondata, const params_struct *restrict params,
    const REAL xx0, const REAL xx1, const REAL xx2,  REAL *restrict psi_background, REAL *restrict ADD_times_AUU
"""
    body = ccg.c_codegen(
        [psi_background, ADD_times_AUU],
        ["*psi_background", "*ADD_times_AUU"],
        verbose=False,
        include_braces=False,
    )
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_auxevol_gfs_all_points(
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for the AUXEVOL grid functions at all points.

    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""Set AUXEVOL gridfunctions at all points."""
    cfunc_type = "void"
    name = "auxevol_gfs_all_points"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    psi_background_memaccess = gri.BHaHGridFunction.access_gf("psi_background")
    ADD_times_AUU_memaccess = gri.BHaHGridFunction.access_gf("ADD_times_AUU")

    body = r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];
  REAL *restrict in_gfs = griddata[grid].gridfuncs.auxevol_gfs;
"""
    body += lp.simple_loop(
        loop_body="auxevol_gfs_single_point(commondata, params, xx0,xx1,xx2,"
        f"&{psi_background_memaccess},"
        f"&{ADD_times_AUU_memaccess});",
        read_xxs=True,
        loop_region="all points",
        OMP_collapse=OMP_collapse,
    )
    body += "}\n"
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_initialize_constant_auxevol() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function to call all functions that set up AUXEVOL gridfunctions.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""Call functions that set up all AUXEVOL gridfunctions."""
    cfunc_type = "void"
    name = "initialize_constant_auxevol"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    body = r"""
    // Set up variable wavespeed
    variable_wavespeed_gfs_all_points(commondata, griddata);

    // Set up all other AUXEVOL gridfunctions
    auxevol_gfs_all_points(commondata, griddata);
    """

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
