"""
C functions for setting the initial data of the wave equation to an exact solution at a particular time.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.simple_loop as lp
from nrpy.equations.wave_equation.WaveEquation_Solutions_InitialData import (
    WaveEquation_solution_Cartesian,
)


def register_CFunction_exact_solution_single_Cartesian_point(
    WaveType: str = "SphericalGaussian",
    default_sigma: float = 3.0,
    default_k0: float = 1.0,
    default_k1: float = 1.0,
    default_k2: float = 1.0,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for the exact solution at a single point.

    :param WaveType: The type of wave: SphericalGaussian or PlaneWave
    :param default_sigma: The default value for the Gaussian width (sigma).
    :param default_k0: The default value for the plane wave wavenumber k in the x-direction.
    :param default_k1: The default value for the plane wave wavenumber k in the y-direction.
    :param default_k2: The default value for the plane wave wavenumber k in the z-direction.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Populate uu_ID, vv_ID
    exactsoln = WaveEquation_solution_Cartesian(
        WaveType=WaveType,
        default_sigma=default_sigma,
        default_k0=default_k0,
        default_k1=default_k1,
        default_k2=default_k2,
    )

    includes = ["BHaH_defines.h"]

    desc = r"""Exact solution at a single Cartesian point (x, y, z) = (xCart0, xCart1, xCart2)."""
    cfunc_type = "void"
    name = "exact_solution_single_Cartesian_point"
    params = r"""const commondata_struct *restrict commondata, const params_struct *restrict params,
    const REAL xCart0, const REAL xCart1, const REAL xCart2,  REAL *restrict exact_soln_UUGF, REAL *restrict exact_soln_VVGF
"""
    body = ccg.c_codegen(
        [exactsoln.uu_exactsoln, exactsoln.vv_exactsoln],
        ["*exact_soln_UUGF", "*exact_soln_VVGF"],
        verbose=False,
        include_braces=False,
    )
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_initial_data(
    OMP_collapse: int,
    enable_checkpointing: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the initial data function for the wave equation with specific parameters.

    :param OMP_collapse: Degree of OpenMP loop collapsing.
    :param enable_checkpointing: Attempt to read from a checkpoint file before generating initial data.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""Set initial data to params.time==0 corresponds to the initial data."""
    cfunc_type = "void"
    name = "initial_data"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    uu_gf_memaccess = gri.BHaHGridFunction.access_gf("uu")
    vv_gf_memaccess = gri.BHaHGridFunction.access_gf("vv")
    body = ""
    if enable_checkpointing:
        body += """// Attempt to read checkpoint file. If it doesn't exist, then continue. Otherwise return.
if( read_checkpoint(commondata, griddata) ) return;
"""
    body += r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];
  REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
"""
    body += lp.simple_loop(
        loop_body="REAL xCart[3]; xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);\n"
        "exact_solution_single_Cartesian_point(commondata, params, xCart[0], xCart[1], xCart[2],"
        f"&{uu_gf_memaccess},"
        f"&{vv_gf_memaccess});",
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
