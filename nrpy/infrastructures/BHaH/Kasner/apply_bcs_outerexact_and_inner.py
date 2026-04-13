"""
Register a Kasner-specific exact outer-BC C function for BHaH GR evolutions.

This module generates a C function that sets outer boundary points to the exact
Kasner BSSN solution (for selected evolved fields), then applies standard inner
boundary conditions.

Author: Nishita Jadoo
        njadoo **at** uidaho **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
from nrpy.equations.general_relativity.kasner_exact import kasner_exact_bssn_exprs


def register_CFunction_apply_bcs_outerexact_kasner_and_inner(
    CoordSystem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C function to set exact Kasner outer BCs and then apply inner BCs.

    :param CoordSystem: Coordinate system used to build exact Kasner expressions.
    :return: None during registration phase, else the NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    exact_bssn = kasner_exact_bssn_exprs(CoordSystem)
    exact_AbarDD = exact_bssn["AbarDD"]
    exact_hDD = exact_bssn["hDD"]
    exact_lambdaU = exact_bssn["lambdaU"]

    ghat_reads = ""
    if CoordSystem.startswith("GeneralRFM"):
        ghat_reads = """
          const REAL ghatDD00 = auxevol_gfs[IDX4pt(GHATDD00GF, idx3)];
          const REAL ghatDD01 = auxevol_gfs[IDX4pt(GHATDD01GF, idx3)];
          const REAL ghatDD02 = auxevol_gfs[IDX4pt(GHATDD02GF, idx3)];
          const REAL ghatDD11 = auxevol_gfs[IDX4pt(GHATDD11GF, idx3)];
          const REAL ghatDD12 = auxevol_gfs[IDX4pt(GHATDD12GF, idx3)];
          const REAL ghatDD22 = auxevol_gfs[IDX4pt(GHATDD22GF, idx3)];
          const REAL ghatDDdD000 = auxevol_gfs[IDX4pt(GHATDDDD000GF, idx3)];
          const REAL ghatDDdD001 = auxevol_gfs[IDX4pt(GHATDDDD001GF, idx3)];
          const REAL ghatDDdD002 = auxevol_gfs[IDX4pt(GHATDDDD002GF, idx3)];
          const REAL ghatDDdD010 = auxevol_gfs[IDX4pt(GHATDDDD010GF, idx3)];
          const REAL ghatDDdD011 = auxevol_gfs[IDX4pt(GHATDDDD011GF, idx3)];
          const REAL ghatDDdD012 = auxevol_gfs[IDX4pt(GHATDDDD012GF, idx3)];
          const REAL ghatDDdD020 = auxevol_gfs[IDX4pt(GHATDDDD020GF, idx3)];
          const REAL ghatDDdD021 = auxevol_gfs[IDX4pt(GHATDDDD021GF, idx3)];
          const REAL ghatDDdD022 = auxevol_gfs[IDX4pt(GHATDDDD022GF, idx3)];
          const REAL ghatDDdD110 = auxevol_gfs[IDX4pt(GHATDDDD110GF, idx3)];
          const REAL ghatDDdD111 = auxevol_gfs[IDX4pt(GHATDDDD111GF, idx3)];
          const REAL ghatDDdD112 = auxevol_gfs[IDX4pt(GHATDDDD112GF, idx3)];
          const REAL ghatDDdD120 = auxevol_gfs[IDX4pt(GHATDDDD120GF, idx3)];
          const REAL ghatDDdD121 = auxevol_gfs[IDX4pt(GHATDDDD121GF, idx3)];
          const REAL ghatDDdD122 = auxevol_gfs[IDX4pt(GHATDDDD122GF, idx3)];
          const REAL ghatDDdD220 = auxevol_gfs[IDX4pt(GHATDDDD220GF, idx3)];
          const REAL ghatDDdD221 = auxevol_gfs[IDX4pt(GHATDDDD221GF, idx3)];
          const REAL ghatDDdD222 = auxevol_gfs[IDX4pt(GHATDDDD222GF, idx3)];
"""

    assign_code = ccg.c_codegen(
        [
            exact_AbarDD[0][0],
            exact_AbarDD[0][1],
            exact_AbarDD[0][2],
            exact_AbarDD[1][1],
            exact_AbarDD[1][2],
            exact_AbarDD[2][2],
            exact_bssn["cf"],
            exact_bssn["trK"],
            exact_hDD[0][0],
            exact_hDD[0][1],
            exact_hDD[0][2],
            exact_hDD[1][1],
            exact_hDD[1][2],
            exact_hDD[2][2],
            exact_lambdaU[0],
            exact_lambdaU[1],
            exact_lambdaU[2],
        ],
        [
            "gfs[IDX4pt(ADD00GF, idx3)]",
            "gfs[IDX4pt(ADD01GF, idx3)]",
            "gfs[IDX4pt(ADD02GF, idx3)]",
            "gfs[IDX4pt(ADD11GF, idx3)]",
            "gfs[IDX4pt(ADD12GF, idx3)]",
            "gfs[IDX4pt(ADD22GF, idx3)]",
            "gfs[IDX4pt(CFGF, idx3)]",
            "gfs[IDX4pt(TRKGF, idx3)]",
            "gfs[IDX4pt(HDD00GF, idx3)]",
            "gfs[IDX4pt(HDD01GF, idx3)]",
            "gfs[IDX4pt(HDD02GF, idx3)]",
            "gfs[IDX4pt(HDD11GF, idx3)]",
            "gfs[IDX4pt(HDD12GF, idx3)]",
            "gfs[IDX4pt(HDD22GF, idx3)]",
            "gfs[IDX4pt(LAMBDAU0GF, idx3)]",
            "gfs[IDX4pt(LAMBDAU1GF, idx3)]",
            "gfs[IDX4pt(LAMBDAU2GF, idx3)]",
        ],
        include_braces=False,
        verbose=False,
    ).replace("\n", "\n          ")

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Apply exact Kasner outer BCs, then apply inner BCs."
    cfunc_type = "void"
    name = "apply_bcs_outerexact_kasner_and_inner"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL *restrict xx[3], const bc_struct *restrict bcstruct, REAL *restrict gfs, const REAL *restrict auxevol_gfs"
    body = f"""
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  const REAL KASNER_t_phys = commondata->KASNER_t0 + commondata->time;
  const REAL KASNER_p1 = commondata->KASNER_p1;
  const REAL KASNER_p2 = commondata->KASNER_p2;
  const REAL KASNER_p3 = commondata->KASNER_p3;

  const bc_info_struct *bc_info = &bcstruct->bc_info;
  for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {{
    for (int dirn = 0; dirn < 3; dirn++) {{
      const int num_pure_outer_boundary_points = bc_info->num_pure_outer_boundary_points[which_gz][dirn];
      if (num_pure_outer_boundary_points <= 0) continue;

      const size_t gz_idx = (size_t)dirn + (size_t)3 * (size_t)which_gz;
      const outerpt_bc_struct *restrict pure_outer_bc_array = bcstruct->pure_outer_bc_array[gz_idx];

#pragma omp parallel for
      for (int idx2d = 0; idx2d < num_pure_outer_boundary_points; idx2d++) {{
        const short i0 = pure_outer_bc_array[idx2d].i0;
        const short i1 = pure_outer_bc_array[idx2d].i1;
        const short i2 = pure_outer_bc_array[idx2d].i2;
        const int idx3 = IDX3(i0, i1, i2);
        const REAL xx0 = xx[0][i0];
        const REAL xx1 = xx[1][i1];
        const REAL xx2 = xx[2][i2];
#include "set_CodeParameters.h"
{ghat_reads}
          {assign_code}
      }}
    }}
  }}

  apply_bcs_inner_only(commondata, params, bcstruct, gfs);
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
