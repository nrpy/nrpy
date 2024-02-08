"""
Generates function to compute the ADM variables from the BSSN variables.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Samuel Cupp
"""

from typing import Union, cast
from inspect import currentframe as cfr
from types import FrameType as FT

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg

import nrpy.infrastructures.CarpetX.simple_loop as lp
from nrpy.infrastructures.CarpetX.CarpetX_include_header import define_standard_includes
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM


def register_CFunction_BSSN_to_ADM(
    thorn_name: str,
    CoordSystem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Convert BSSN variables in the Cartesian basis to ADM variables in the Cartesian basis.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = define_standard_includes()
    desc = r"""Perform BSSN-to-ADM conversion. Useful for diagnostics."""
    name = f"{thorn_name}_BSSN_to_ADM"
    body = f"""  DECLARE_CCTK_ARGUMENTSX_{name};
  DECLARE_CCTK_PARAMETERS;

"""

    bssn2adm = BSSN_to_ADM(CoordSystem)

    # lapse, shift, and dtshift are just straight copies
    lapse_gf_access = gri.CarpetXGridFunction.access_gf(
        gf_name="alp", use_GF_suffix=False
    )
    lapse2_gf_access = gri.CarpetXGridFunction.access_gf(gf_name="alpha")

    loop_body = lapse_gf_access + " = " + lapse2_gf_access + ";\n"

    coord_name = ["x", "y", "z"]
    for i in range(3):
        bssn_shift_gf_access = gri.CarpetXGridFunction.access_gf(
            gf_name="vetU" + str(i)
        )
        shift_gf_access = gri.CarpetXGridFunction.access_gf(
            gf_name="beta" + coord_name[i], use_GF_suffix=False
        )
        loop_body += f"{shift_gf_access} = {bssn_shift_gf_access};\n"

    for i in range(3):
        bssn_dtshift_gf_access = gri.CarpetXGridFunction.access_gf(
            gf_name="betU" + str(i)
        )
        dtshift_gf_access = gri.CarpetXGridFunction.access_gf(
            gf_name="dtbeta" + coord_name[i], use_GF_suffix=False
        )
        loop_body += f"{dtshift_gf_access} = {bssn_dtshift_gf_access};\n"

    list_of_output_exprs = []
    list_of_output_varnames = []

    cf_gf_access = gri.CarpetXGridFunction.access_gf(gf_name="cf")
    trK_gf_access = gri.CarpetXGridFunction.access_gf(gf_name="trK")
    loop_body += f"const CCTK_REAL cf = {cf_gf_access};\n"
    loop_body += f"const CCTK_REAL trK = {trK_gf_access};\n"

    for i in range(3):
        for j in range(i, 3):
            hDD_gf_access = gri.CarpetXGridFunction.access_gf(
                gf_name="hDD" + str(i) + str(j)
            )
            loop_body += f"const CCTK_REAL hDD{i}{j} = {hDD_gf_access};\n"

            gamma_gf_access = gri.CarpetXGridFunction.access_gf(
                gf_name="g" + coord_name[i] + coord_name[j], use_GF_suffix=False
            )
            list_of_output_exprs += [bssn2adm.gammaDD[i][j]]
            list_of_output_varnames += [gamma_gf_access]

    for i in range(3):
        for j in range(i, 3):
            aDD_gf_access = gri.CarpetXGridFunction.access_gf(
                gf_name="aDD" + str(i) + str(j)
            )
            loop_body += f"const CCTK_REAL aDD{i}{j} = {aDD_gf_access};\n"

            curv_gf_access = gri.CarpetXGridFunction.access_gf(
                gf_name="k" + coord_name[i] + coord_name[j], use_GF_suffix=False
            )
            list_of_output_exprs += [bssn2adm.KDD[i][j]]
            list_of_output_varnames += [curv_gf_access]

    loop_body += ccg.c_codegen(
        list_of_output_exprs,
        list_of_output_varnames,
        verbose=False,
        include_braces=False,
        enable_simd=False,
    )
    loop_body = loop_body.rstrip()

    body += lp.simple_loop(
        loop_body=loop_body,
        loop_region="all points",
        enable_simd=False,
    )

    schedule = (
        "ODESolvers_PostStep",
        f"""
schedule FUNC_NAME in ODESolvers_PostStep after {thorn_name}_enforce_detgammahat_constraint before ADMBaseX_SetADMVars
{{
  LANG: C
  READS:  aDD00GF, aDD01GF, aDD02GF, aDD11GF, aDD12GF, aDD22GF,
          hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF,
          vetU0GF, vetU1GF, vetU2GF, betU0GF, betU1GF, betU2GF,
          cfGF, trKGF, alphaGF
  WRITES: ADMBaseX::metric(everywhere),
          ADMBaseX::shift(everywhere),
          ADMBaseX::curv(everywhere),
          ADMBaseX::dtshift(everywhere),
          ADMBaseX::lapse(everywhere)
}} "Perform BSSN-to-ADM conversion. Useful for diagnostics."
""",
    )

    ET_schedule_bins_entries = [schedule]

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=includes,
        desc=desc,
        c_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=ET_schedule_bins_entries,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
