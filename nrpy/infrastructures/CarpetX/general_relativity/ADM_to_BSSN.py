"""
Generates function to compute the BSSN variables from the ADM variables.

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
import nrpy.params as par
import nrpy.indexedexp as ixp
import nrpy.helpers.parallel_codegen as pcg

from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
import nrpy.infrastructures.CarpetX.simple_loop as lp
from nrpy.infrastructures.CarpetX.CarpetX_include_header import define_standard_includes
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support


def register_CFunction_ADM_to_BSSN(
    thorn_name: str,
    CoordSystem: str,
    fd_order: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Convert ADM variables in the Cartesian basis to BSSN variables in the Cartesian basis.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param fd_order: Order of finite difference method

    :return: A string representing the full C function.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    old_fd_order = par.parval_from_str("fd_order")
    # Set this because parallel codegen needs the correct local values
    par.set_parval_from_str("fd_order", fd_order)

    desc = f"""Converting from ADM to BSSN quantities is required in the Einstein Toolkit,
as initial data are given in terms of ADM quantities, and {thorn_name} evolves the BSSN quantities."""
    name = f"{thorn_name}_ADM_to_BSSN_order_{fd_order}"
    body = f"""  DECLARE_CCTK_ARGUMENTSX_{name};
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL invdxx0 = 1.0 / CCTK_DELTA_SPACE(0);
  const CCTK_REAL invdxx1 = 1.0 / CCTK_DELTA_SPACE(1);
  const CCTK_REAL invdxx2 = 1.0 / CCTK_DELTA_SPACE(2);

"""
    lapse_gf_access = gri.CarpetXGridFunction.access_gf(
        gf_name="alp", use_GF_suffix=False
    )
    lapse2_gf_access = gri.CarpetXGridFunction.access_gf(gf_name="alpha")
    loop_body = lapse2_gf_access + " = " + lapse_gf_access + ";\n"

    coord_name = ["x", "y", "z"]
    for i in range(3):
        shift_gf_access = gri.CarpetXGridFunction.access_gf(
            gf_name="beta" + coord_name[i], use_GF_suffix=False
        )
        loop_body += f"const CCTK_REAL local_betaU{i} = {shift_gf_access};\n"

    for i in range(3):
        dtshift_gf_access = gri.CarpetXGridFunction.access_gf(
            gf_name="dtbeta" + coord_name[i], use_GF_suffix=False
        )
        loop_body += f"const CCTK_REAL local_BU{i} = {dtshift_gf_access};\n"

    for i in range(3):
        for j in range(i, 3):
            gamma_gf_access = gri.CarpetXGridFunction.access_gf(
                gf_name="g" + coord_name[i] + coord_name[j], use_GF_suffix=False
            )
            loop_body += f"const CCTK_REAL local_gDD{i}{j} = {gamma_gf_access};\n"

    for i in range(3):
        for j in range(i, 3):
            curv_gf_access = gri.CarpetXGridFunction.access_gf(
                gf_name="k" + coord_name[i] + coord_name[j], use_GF_suffix=False
            )
            loop_body += f"const CCTK_REAL local_kDD{i}{j} = {curv_gf_access};\n"

    gammaCartDD = ixp.declarerank2("local_gDD", symmetry="sym01")
    KCartDD = ixp.declarerank2("local_kDD", symmetry="sym01")

    betaU = ixp.declarerank1("local_betaU")
    BU = ixp.declarerank1("local_BU")
    adm2bssn = ADM_to_BSSN(gammaCartDD, KCartDD, betaU, BU, "Cartesian")

    cf_gf_access = gri.CarpetXGridFunction.access_gf(gf_name="cf")
    trK_gf_access = gri.CarpetXGridFunction.access_gf(gf_name="trK")
    list_of_output_exprs = [adm2bssn.cf, adm2bssn.trK]
    list_of_output_varnames = [cf_gf_access, trK_gf_access]
    for i in range(3):
        vetU_gf_access = gri.CarpetXGridFunction.access_gf(gf_name=f"vetU{i}")
        list_of_output_exprs += [adm2bssn.vetU[i]]
        list_of_output_varnames += [vetU_gf_access]
    for i in range(3):
        betU_gf_access = gri.CarpetXGridFunction.access_gf(gf_name=f"betU{i}")
        list_of_output_exprs += [adm2bssn.betU[i]]
        list_of_output_varnames += [betU_gf_access]
    for i in range(3):
        for j in range(i, 3):
            hDD_gf_access = gri.CarpetXGridFunction.access_gf(gf_name=f"hDD{i}{j}")
            list_of_output_exprs += [adm2bssn.hDD[i][j]]
            list_of_output_varnames += [hDD_gf_access]
    for i in range(3):
        for j in range(i, 3):
            aDD_gf_access = gri.CarpetXGridFunction.access_gf(gf_name=f"aDD{i}{j}")
            list_of_output_exprs += [adm2bssn.aDD[i][j]]
            list_of_output_varnames += [aDD_gf_access]

    loop_body += ccg.c_codegen(
        list_of_output_exprs,
        list_of_output_varnames,
        verbose=False,
        include_braces=False,
    )
    loop_body = loop_body.rstrip()

    body += lp.simple_loop(
        loop_body=loop_body,
        enable_simd=False,
        loop_region="all points",
    )

    body += "\n"

    Bq = BSSN_quantities[CoordSystem]
    gammabarUU = Bq.gammabarUU
    GammabarUDD = Bq.GammabarUDD

    # Next evaluate \bar{\Lambda}^i, based on GammabarUDD above and GammahatUDD
    # (from the reference metric):
    rfm = refmetric.reference_metric[CoordSystem]
    LambdabarU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                LambdabarU[i] += gammabarUU[j][k] * (
                    GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k]
                )

    # Finally apply rescaling:
    # lambda^i = Lambdabar^i/\text{ReU[i]}
    lambdaU = ixp.zerorank1()
    for i in range(3):
        lambdaU[i] = LambdabarU[i] / rfm.ReU[i]

    body += lp.simple_loop(
        ccg.c_codegen(
            lambdaU,
            [
                gri.CarpetXGridFunction.access_gf("lambdaU0"),
                gri.CarpetXGridFunction.access_gf("lambdaU1"),
                gri.CarpetXGridFunction.access_gf("lambdaU2"),
            ],
            verbose=False,
            include_braces=False,
            enable_fd_codegen=True,
        ),
        loop_region="interior",
    )

    body += """
  //ExtrapolateGammas(cctkGH, lambdaU0GF);
  //ExtrapolateGammas(cctkGH, lambdaU1GF);
  //ExtrapolateGammas(cctkGH, lambdaU2GF);
"""

    schedule = f"""
if(FD_order == {fd_order}) {{
  schedule FUNC_NAME at CCTK_INITIAL after ADMBaseX_PostInitial
  {{
    LANG: C
    READS:  ADMBaseX::metric,
            ADMBaseX::shift,
            ADMBaseX::curv,
            ADMBaseX::dtshift,
            ADMBaseX::lapse
    WRITES: evol_variables
    SYNC: evol_variables
  }} "Convert initial data into BSSN variables"
}}
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
        ET_schedule_bins_entries=[("CCTK_INITIAL", schedule)],
    )

    # Reset to the initial values
    par.set_parval_from_str("fd_order", old_fd_order)

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
