"""
Generates function to raise the stress-energy tensor using the BSSN variables.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Samuel Cupp
"""

from typing import Union, cast, List
from inspect import currentframe as cfr
from types import FrameType as FT
import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.indexedexp as ixp
import nrpy.helpers.parallel_codegen as pcg

import nrpy.infrastructures.ETLegacy.simple_loop as lp
from nrpy.infrastructures.ETLegacy.ETLegacy_include_header import (
    define_standard_includes,
)
import nrpy.equations.general_relativity.g4munu_conversions as g4conv


def register_CFunction_T4DD_to_T4UU(
    thorn_name: str,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the function that enforces the det(gammabar) = det(gammahat) constraint.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    desc = r"""Compute T4UU from T4DD (provided by TmunuBase),
using BSSN quantities as inputs for the 4D raising operation

WARNING: Do not enable SIMD here, as it is not guaranteed that
         cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2] is a multiple of
         SIMD_width!"""
    name = f"{thorn_name}_T4DD_to_T4UU"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
  DECLARE_CCTK_PARAMETERS;

"""

    g4UU = g4conv.BSSN_to_g4UU(
        CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
    )

    T4DD = ixp.declarerank2("T4DD", dimension=4, symmetry="sym01")
    T4UU = ixp.zerorank2(dimension=4)
    for mu in range(4):
        for nu in range(4):
            for alpha in range(4):
                for beta in range(4):
                    T4UU[mu][nu] += T4DD[alpha][beta] * g4UU[mu][alpha] * g4UU[nu][beta]

    T4DD_access_gfs: List[str] = []
    T4UU_expr_list: List[sp.Expr] = []
    for i in range(4):
        for j in range(i, 4):
            T4DD_access_gfs += [
                gri.ETLegacyGridFunction.access_gf(gf_name=f"T4UU{i}{j}")
            ]
            T4UU_expr_list += [T4UU[i][j]]

    loop_body = ""
    lapse_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name="alpha")
    loop_body += f"const CCTK_REAL alpha = {lapse_gf_access};\n"
    cf_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name="cf")
    loop_body += f"const CCTK_REAL cf = {cf_gf_access};\n"

    for i in range(3):
        vet_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name=f"vetU{i}")
        loop_body += f"const CCTK_REAL vetU{i} = {vet_gf_access};\n"

    for i in range(3):
        for j in range(i, 3):
            hDD_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name=f"hDD{i}{j}")
            loop_body += f"const CCTK_REAL hDD{i}{j} = {hDD_gf_access};\n"

    coord_name = ["t", "x", "y", "z"]
    for i in range(4):
        for j in range(i, 4):
            Tmunu_gf_access = gri.ETLegacyGridFunction.access_gf(
                gf_name="eT" + coord_name[i] + coord_name[j], use_GF_suffix=False
            )
            loop_body += f"const CCTK_REAL T4DD{i}{j} = {Tmunu_gf_access};\n"

    loop_body += ccg.c_codegen(
        T4UU_expr_list,
        T4DD_access_gfs,
        enable_simd=False,
    )

    body += lp.simple_loop(
        loop_body=loop_body,
        loop_region="all points",
        enable_simd=False,
        OMP_collapse=OMP_collapse,
    )

    schedule1 = f"""
schedule FUNC_NAME in MoL_CalcRHS before {thorn_name}_RHS
{{
  LANG: C
  READS:  TmunuBase::stress_energy_scalar,
          TmunuBase::stress_energy_vector,
          TmunuBase::stress_energy_tensor,
          hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF,
          alphaGF, cfGF, vetU0GF, vetU1GF, vetU2GF
  WRITES: T4UU00GF(everywhere), T4UU01GF(everywhere), T4UU02GF(everywhere), T4UU03GF(everywhere),
          T4UU11GF(everywhere), T4UU12GF(everywhere), T4UU13GF(everywhere),
          T4UU22GF(everywhere), T4UU23GF(everywhere), T4UU33GF(everywhere)
}} "Compute T4UU from T4DD (provided in eT?? from TmunuBase), needed for BSSN RHSs"
"""
    schedule2 = f"""
schedule FUNC_NAME in MoL_PseudoEvolution before {thorn_name}_BSSN_constraints
{{
  LANG: C
  READS:  TmunuBase::stress_energy_scalar,
          TmunuBase::stress_energy_vector,
          TmunuBase::stress_energy_tensor,
          hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF,
          alphaGF, cfGF, vetU0GF, vetU1GF, vetU2GF
  WRITES: T4UU00GF(everywhere), T4UU01GF(everywhere), T4UU02GF(everywhere), T4UU03GF(everywhere),
          T4UU11GF(everywhere), T4UU12GF(everywhere), T4UU13GF(everywhere),
          T4UU22GF(everywhere), T4UU23GF(everywhere), T4UU33GF(everywhere)
}} "Compute T4UU from T4DD (provided in eT?? from TmunuBase), needed for BSSN constraints"
"""
    schedule_RHS = ("MoL_CalcRHS", schedule1)
    schedule_constraints = ("MoL_PseudoEvolution", schedule2)

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=define_standard_includes(),
        desc=desc,
        c_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[schedule_RHS, schedule_constraints],
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
