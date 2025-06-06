"""
Generates a function to enforce a minimum on the lapse grid function.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Samuel Cupp
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.CarpetX.simple_loop as lp
from nrpy.infrastructures.CarpetX.CarpetX_include_header import define_standard_includes


def register_CFunction_floor_the_lapse(
    thorn_name: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the function that enforces the det(gammabar) = det(gammahat) constraint.

    :param thorn_name: The Einstein Toolkit thorn name.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    desc = """Apply floor to the lapse."""
    name = f"{thorn_name}_floor_the_lapse"
    body = f"""  DECLARE_CCTK_ARGUMENTSX_{name};
  DECLARE_CCTK_PARAMETERS;

#ifndef MAX
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#endif

"""
    lapse_access_gfs = [gri.CarpetXGridFunction.access_gf(gf_name="alpha")]
    loop_body = (
        lapse_access_gfs[0] + " = MAX(" + lapse_access_gfs[0] + ", lapse_floor);"
    )

    body += lp.simple_loop(
        loop_body=loop_body,
        loop_region="all points",
        enable_simd=False,
    )

    schedule = f"""
schedule FUNC_NAME in ODESolvers_PostStep before {thorn_name}_enforce_detgammahat_constraint
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
        cfunc_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[("ODESolvers_PostStep", schedule)],
        ET_current_thorn_CodeParams_used=["lapse_floor"],
    )
    return pcg.NRPyEnv()
