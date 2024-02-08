"""
Generates a function to enforce a minimum on the lapse grid function.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Samuel Cupp
"""

from typing import Union, cast
from inspect import currentframe as cfr
from types import FrameType as FT

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.ETLegacy.simple_loop as lp
from nrpy.infrastructures.ETLegacy.ETLegacy_include_header import (
    define_standard_includes,
)


def register_CFunction_floor_the_lapse(
    thorn_name: str,
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the function that enforces the det(gammabar) = det(gammahat) constraint.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    desc = """Apply floor to the lapse."""
    name = f"{thorn_name}_floor_the_lapse"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
  DECLARE_CCTK_PARAMETERS;

#ifndef MAX
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#endif

"""
    lapse_access_gfs = [gri.ETLegacyGridFunction.access_gf(gf_name="alpha")]
    loop_body = (
        lapse_access_gfs[0] + " = MAX(" + lapse_access_gfs[0] + ", lapse_floor);"
    )

    body += lp.simple_loop(
        loop_body=loop_body,
        loop_region="all points",
        enable_simd=False,
        OMP_collapse=OMP_collapse,
    )

    schedule = f"""
schedule FUNC_NAME in MoL_PostStep before {thorn_name}_enforce_detgammahat_constraint
{{
  LANG: C
  READS:  alphaGF(everywhere)
  WRITES: alphaGF(everywhere)
}} "Set lapse = max(lapse_floor, lapse)"
"""

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=define_standard_includes(),
        desc=desc,
        c_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[("MoL_PostStep", schedule)],
        ET_current_thorn_CodeParams_used=["lapse_floor"],
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
