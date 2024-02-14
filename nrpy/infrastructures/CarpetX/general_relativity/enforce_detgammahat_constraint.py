"""
Generates function to enforce the det(gammabar) = det(gammahat) constraint..

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

import nrpy.infrastructures.CarpetX.simple_loop as lp
from nrpy.infrastructures.CarpetX.CarpetX_include_header import define_standard_includes
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support


def register_CFunction_enforce_detgammahat_constraint(
    thorn_name: str,
    CoordSystem: str,
    enable_rfm_precompute: bool,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the function that enforces the det(gammabar) = det(gammahat) constraint.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    Bq = BSSN_quantities[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    rfm = refmetric.reference_metric[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    desc = r"""Enforce det(gammabar) = det(gammahat) constraint. Required for strong hyperbolicity."""
    name = f"{thorn_name}_enforce_detgammahat_constraint"
    body = f"""  DECLARE_CCTK_ARGUMENTSX_{name};
  DECLARE_CCTK_PARAMETERS;

"""

    # First define the Kronecker delta:
    KroneckerDeltaDD = ixp.zerorank2()
    for i in range(3):
        KroneckerDeltaDD[i][i] = sp.sympify(1)

    # The detgammabar in BSSN_RHSs is set to detgammahat when BSSN_RHSs::detgbarOverdetghat_equals_one=True (default),
    #    so we manually compute it here:
    dummygammabarUU, detgammabar = ixp.symm_matrix_inverter3x3(Bq.gammabarDD)

    # Next apply the constraint enforcement equation above.
    nrpyAbs = sp.Function("nrpyAbs")
    hprimeDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            hprimeDD[i][j] = (nrpyAbs(rfm.detgammahat) / detgammabar) ** (
                sp.Rational(1, 3)
            ) * (KroneckerDeltaDD[i][j] + Bq.hDD[i][j]) - KroneckerDeltaDD[i][j]

    hDD_access_gfs: List[str] = []
    hprimeDD_expr_list: List[sp.Expr] = []
    for i in range(3):
        for j in range(i, 3):
            hDD_access_gfs += [gri.CarpetXGridFunction.access_gf(gf_name=f"hDD{i}{j}")]
            hprimeDD_expr_list += [hprimeDD[i][j]]

    # To evaluate the cube root, SIMD support requires e.g., SLEEF.
    #   Also need to be careful to not access memory out of bounds!
    #   After all this is a loop over ALL POINTS.
    #   Exercise to the reader: prove that for any reasonable grid,
    #   SIMD loops over grid interiors never write data out of bounds
    #   and are threadsafe for any reasonable number of threads.
    body += lp.simple_loop(
        loop_body=ccg.c_codegen(
            hprimeDD_expr_list,
            hDD_access_gfs,
            enable_simd=False,
            automatically_read_gf_data_from_memory=True,
            enable_fd_codegen=True,
            enable_fd_functions=True,
        ),
        loop_region="all points",
        enable_simd=False,
    )

    schedule = """
schedule FUNC_NAME in ODESolvers_PostStep
{
  LANG: C
  READS:  hDD00GF(everywhere), hDD01GF(everywhere), hDD02GF(everywhere),
          hDD11GF(everywhere), hDD12GF(everywhere), hDD22GF(everywhere)
  WRITES: hDD00GF(everywhere), hDD01GF(everywhere), hDD02GF(everywhere),
          hDD11GF(everywhere), hDD12GF(everywhere), hDD22GF(everywhere)
} "Enforce detgammabar = detgammahat (= 1 in Cartesian)"
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
        ET_schedule_bins_entries=[("ODESolvers_PostStep", schedule)],
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
