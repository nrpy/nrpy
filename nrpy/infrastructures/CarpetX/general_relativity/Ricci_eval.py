"""
Generates function to compute the Ricci tensor for use in the BSSN evolution equations.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Samuel Cupp
"""

from typing import Union, cast, List
from inspect import currentframe as cfr
from types import FrameType as FT

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par
import nrpy.helpers.parallel_codegen as pcg
import nrpy.finite_difference as fin

import nrpy.infrastructures.CarpetX.simple_loop as lp
from nrpy.infrastructures.CarpetX.CarpetX_include_header import define_standard_includes
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities


def register_CFunction_Ricci_eval(
    thorn_name: str,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    fd_order: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the right-hand side evaluation function for the BSSN equations.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param enable_simd: Whether or not to enable SIMD (Single Instruction, Multiple Data).
    :param fd_order: Order of finite difference method

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    old_fd_order = par.parval_from_str("fd_order")
    # Set this because parallel codegen needs the correct local values
    par.set_parval_from_str("fd_order", fd_order)
    par.set_parval_from_str("enable_RbarDD_gridfunctions", False)

    includes = define_standard_includes()
    if enable_simd:
        includes += [("./simd/simd_intrinsics.h")]
    desc = r"""Compute Ricci tensor for the BSSN evolution equations."""
    name = f"{thorn_name}_Ricci_eval_order_{fd_order}"
    body = f"""  DECLARE_CCTK_ARGUMENTSX_{name};
"""
    if enable_simd:
        body += """
  const REAL_SIMD_ARRAY invdxx0 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(0));
  const REAL_SIMD_ARRAY invdxx1 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(1));
  const REAL_SIMD_ARRAY invdxx2 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(2));

"""
    else:
        body += """  const CCTK_REAL invdxx0 CCTK_ATTRIBUTE_UNUSED = 1.0/CCTK_DELTA_SPACE(0);
  const CCTK_REAL invdxx1 CCTK_ATTRIBUTE_UNUSED = 1.0/CCTK_DELTA_SPACE(1);
  const CCTK_REAL invdxx2 CCTK_ATTRIBUTE_UNUSED = 1.0/CCTK_DELTA_SPACE(2);
  DECLARE_CCTK_PARAMETERS;

"""

    Bq = BSSN_quantities[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    # Populate Ricci tensor
    Ricci_access_gfs: List[str] = []
    for var in Bq.Ricci_varnames:
        Ricci_access_gfs += [gri.CarpetXGridFunction.access_gf(gf_name=var)]
    body += lp.simple_loop(
        loop_body=ccg.c_codegen(
            Bq.Ricci_exprs,
            Ricci_access_gfs,
            enable_fd_codegen=True,
            enable_simd=enable_simd,
            enable_fd_functions=True,
            enable_GoldenKernels=True,
        ),
        loop_region="interior",
        enable_simd=enable_simd,
    )

    schedule = f"""
if(FD_order == {fd_order}) {{
  schedule FUNC_NAME in ODESolvers_RHS as {thorn_name}_Ricci before {thorn_name}_RHS
  {{
    LANG: C
    READS:  hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF,
            lambdaU0GF, lambdaU1GF, lambdaU2GF
    WRITES: RbarDD00GF, RbarDD01GF, RbarDD02GF, RbarDD11GF, RbarDD12GF, RbarDD22GF
  }} "Compute Ricci tensor, needed for BSSN RHSs, at finite-differencing order {fd_order}"
}}
"""

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=includes,
        desc=desc,
        c_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        prefunc=fin.construct_FD_functions_prefunc(),
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[("ODESolvers_RHS", schedule)],
    )
    # Reset to the initial values
    par.set_parval_from_str("fd_order", old_fd_order)
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
