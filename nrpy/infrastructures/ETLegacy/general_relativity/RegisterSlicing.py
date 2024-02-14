"""
Generates function to register the space-time 3+1 slicing condition with Cactus.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Samuel Cupp
"""

from typing import Union, cast
from inspect import currentframe as cfr
from types import FrameType as FT

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_RegisterSlicing(
    thorn_name: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the slicing registration function.

    :param thorn_name: The Einstein Toolkit thorn name.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    desc = f"""Register slicing condition for NRPy+-generated thorn {thorn_name}."""
    name = f"{thorn_name}_RegisterSlicing"
    body = f"""
    Einstein_RegisterSlicing ("{thorn_name}");
  return 0;"""

    schedule = """
schedule FUNC_NAME at STARTUP
{
  LANG: C
  OPTIONS: meta
} "Register 3+1 slicing condition"
"""

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=["Slicing.h", "cctk.h"],
        desc=desc,
        c_type="int",
        name=name,
        params="",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[("STARTUP", schedule)],
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
