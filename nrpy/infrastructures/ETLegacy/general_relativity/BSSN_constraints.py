"""
Generates function to compute the constraints H, MU, and MSQUARED.

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

import nrpy.infrastructures.ETLegacy.simple_loop as lp
from nrpy.infrastructures.ETLegacy.ETLegacy_include_header import (
    define_standard_includes,
)
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints


def register_CFunction_BSSN_constraints(
    thorn_name: str,
    enable_T4munu: bool,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    fd_order: int,
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the BSSN constraints evaluation function.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_T4munu: Whether to include the stress-energy tensor.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param enable_simd: Whether or not to enable SIMD instructions.
    :param fd_order: Order of finite difference method
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    old_fd_order = par.parval_from_str("fd_order")
    old_enable_T4munu = par.parval_from_str("enable_T4munu")
    # Set this because parallel codegen needs the correct local values
    par.set_parval_from_str("fd_order", fd_order)
    par.set_parval_from_str("enable_T4munu", enable_T4munu)

    includes = define_standard_includes()
    if enable_simd:
        includes += [("./simd/simd_intrinsics.h")]
    desc = r"""Evaluate BSSN constraints."""
    name = f"{thorn_name}_BSSN_constraints_order_{fd_order}"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
"""
    if enable_simd:
        body += f"""
  const REAL_SIMD_ARRAY invdxx0 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(0));
  const REAL_SIMD_ARRAY invdxx1 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(1));
  const REAL_SIMD_ARRAY invdxx2 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(2));
  const CCTK_REAL *param_PI CCTK_ATTRIBUTE_UNUSED = CCTK_ParameterGet("PI", "{thorn_name}", NULL);
  const REAL_SIMD_ARRAY PI CCTK_ATTRIBUTE_UNUSED = ConstSIMD(*param_PI);

"""
    else:
        body += """  DECLARE_CCTK_PARAMETERS;

"""

    Bcon = BSSN_constraints[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    list_of_output_exprs = [Bcon.H]
    Constraints_access_gfs: List[str] = [
        gri.ETLegacyGridFunction.access_gf(gf_name="H")
    ]

    for index in range(3):
        list_of_output_exprs += [Bcon.MU[index]]
        Constraints_access_gfs += [
            gri.ETLegacyGridFunction.access_gf(gf_name="MU" + str(index))
        ]
    body += lp.simple_loop(
        loop_body=ccg.c_codegen(
            list_of_output_exprs,
            Constraints_access_gfs,
            enable_fd_codegen=True,
            enable_simd=enable_simd,
            enable_fd_functions=True,
            enable_GoldenKernels=True,
        ),
        loop_region="interior",
        enable_simd=enable_simd,
        OMP_collapse=OMP_collapse,
    )

    schedule = f"""
if(FD_order == {fd_order}) {{
  schedule FUNC_NAME in MoL_PseudoEvolution as {thorn_name}_BSSN_constraints
  {{
    LANG: C
    READS:  aDD00GF, aDD01GF, aDD02GF, aDD11GF, aDD12GF, aDD22GF,
            hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF,
            trKGF, cfGF, lambdaU0GF, lambdaU1GF, lambdaU2GF"""
    if enable_T4munu:
        schedule += """,
            alphaGF, vetU0GF, vetU1GF, vetU2GF,
            T4UU00GF, T4UU01GF, T4UU02GF, T4UU03GF"""
    schedule += f"""
    WRITES: aux_variables
  }} "Compute BSSN (Hamiltonian and momentum) constraints, at finite-differencing order {fd_order}"
}}
"""

    params = None
    if enable_T4munu:
        params = ["PI"]

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
        ET_schedule_bins_entries=[("MoL_PseudoEvolution", schedule)],
        ET_current_thorn_CodeParams_used=params,
    )

    # Reset to the initial values
    par.set_parval_from_str("fd_order", old_fd_order)
    par.set_parval_from_str("enable_T4munu", old_enable_T4munu)

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
