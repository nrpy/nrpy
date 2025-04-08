"""
Generate C function to set cahdprefactor gridfunction when solving the BSSN equations in curvilinear coordinates.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Set, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.reference_metric as refmetric


def register_CFunction_cahdprefactor_auxevol_gridfunction(
    set_of_CoordSystems: Set[str],
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Add function that sets cahdprefactor gridfunction = C_CAHD * CFL_FACTOR * dsmin to Cfunction dictionary.

    :param set_of_CoordSystems: Coordinate systems used.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    for CoordSystem in set_of_CoordSystems:
        desc = "cahdprefactor_auxevol_gridfunction(): Initialize CAHD prefactor (auxevol) gridfunction."
        name = "cahdprefactor_auxevol_gridfunction"
        params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], REAL *restrict auxevol_gfs"

        rfm = refmetric.reference_metric[CoordSystem]
        dxx0, dxx1, dxx2 = sp.symbols("dxx0 dxx1 dxx2", real=True)
        loop_body = r"""  // Compute cahdprefactor gridfunction = C_CAHD * CFL_FACTOR * dsmin.
REAL dsmin0, dsmin1, dsmin2;
"""
        loop_body += ccg.c_codegen(
            [
                sp.Abs(rfm.scalefactor_orthog[0] * dxx0),
                sp.Abs(rfm.scalefactor_orthog[1] * dxx1),
                sp.Abs(rfm.scalefactor_orthog[2] * dxx2),
            ],
            ["dsmin0", "dsmin1", "dsmin2"],
            include_braces=False,
        )
        loop_body += """auxevol_gfs[IDX4(CAHDPREFACTORGF, i0, i1, i2)] = C_CAHD * CFL_FACTOR * MIN(dsmin0, MIN(dsmin1, dsmin2));"""

        cfc.register_CFunction(
            includes=["BHaH_defines.h"],
            desc=desc,
            name=name,
            params=params,
            include_CodeParameters_h=True,
            body=lp.simple_loop(
                loop_body=loop_body,
                loop_region="all points",
                read_xxs=True,
            ),
            CoordSystem_for_wrapper_func=CoordSystem,
        )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
