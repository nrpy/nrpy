"""
Library of C functions for solving the BSSN equations in curvilinear coordinates, using a reference-metric formalism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
from nrpy.equations.general_relativity import psi4_tetrads

def register_CFunction_psi4_tetrad(
    CoordSystem: str,
    tetrad: str = "quasiKinnersley",
    use_metric_to_construct_unit_normal: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C function for psi4 tetrad computations.

    :param CoordSystem: The coordinate system to be used.
    :param tetrad: The type of tetrad. Defaults to "quasiKinnersley".
    :param use_metric_to_construct_unit_normal: Whether to use the metric to construct the unit normal. Defaults to False.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"Compute {tetrad} tetrad for psi4, with use_metric_to_construct_unit_normal={use_metric_to_construct_unit_normal}"
    name = "psi4_tetrad"

    # Initialize psi4 tetrad
    psi4tet_class = psi4_tetrads.Psi4Tetrads(
        CoordSystem,
        enable_rfm_precompute=False,
        tetrad=tetrad,
        use_metric_to_construct_unit_normal=use_metric_to_construct_unit_normal,
    )
    list_of_vrnms = []
    list_of_exprs = []

    for i in range(4):
        list_of_vrnms.append("*mre4U" + str(i))
        list_of_exprs.append(psi4tet_class.mre4U[i])
    for i in range(4):
        list_of_vrnms.append("*mim4U" + str(i))
        list_of_exprs.append(psi4tet_class.mim4U[i])
    for i in range(4):
        list_of_vrnms.append("*n4U" + str(i))
        list_of_exprs.append(psi4tet_class.n4U[i])

    params = "const commondata_struct *restrict commondata, const params_struct *restrict params,"
    list_of_metricvarnames = ["cf"] + [
        f"hDD{i}{j}" for i in range(3) for j in range(i, 3)
    ]
    params += ", ".join([f"const REAL {var}" for var in list_of_metricvarnames])
    params += ", " + ", ".join([f"REAL {var}" for var in list_of_vrnms])
    params += ", REAL *restrict xx[3], const int i0, const int i1, const int i2"

    body = ""

    for i in range(3):
        body += f"  const REAL xx{i} = xx[{i}][i{i}];\n"
    body += "  // Compute tetrads:\n"

    # Sort the lhss list alphabetically, and rhss to match:
    lhss, rhss = [
        list(x)
        for x in zip(
            *sorted(zip(list_of_vrnms, list_of_exprs), key=lambda pair: pair[0])
        )
    ]
    body += ccg.c_codegen(
        rhss,
        lhss,
        verbose=False,
        enable_cse=True,
        include_braces=False,
    )

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
