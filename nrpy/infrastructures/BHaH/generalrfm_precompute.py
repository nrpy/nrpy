# nrpy/infrastructures/BHaH/generalrfm_precompute.py
"""
Generate C functions that precompute GeneralRFM quantities into AUXEVOL gridfunctions.

Currently supports GeneralRFM_fisheyeN* via equations/generalrfm/fisheye.py.

Design:
  * For GeneralRFM, ghat and derivatives are treated as fixed AUXEVOL gridfunctions.
  * Derived quantities (e.g., Gammahat*) are computed algebraically from ghat
    expressions and are not stored as AUXEVOL gridfunctions.
  * This module generates a C function that evaluates the analytic expressions
    (from the GeneralRFM provider) once per gridpoint and stores the results.

Note: CUDA is not yet supported for this precompute path.

Author: Nishita Jadoo
        njadoo@uidaho.edu
"""
from __future__ import annotations

import re
from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, List, Tuple, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.equations.generalrfm import fisheye as generalrfm_fisheye
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)
import nrpy.helpers.parallel_codegen as pcg


def _parse_fisheye_num_transitions(CoordSystem: str) -> int:
    match = re.match(r"GeneralRFM_fisheyeN(\d+)$", CoordSystem)
    if not match:
        raise ValueError(
            f"GeneralRFM CoordSystem {CoordSystem} not supported (expected GeneralRFM_fisheyeN*)."
        )
    return int(match.group(1))


def _build_fisheye_exprs(CoordSystem: str):
    num_transitions = _parse_fisheye_num_transitions(CoordSystem)
    fisheye = generalrfm_fisheye.build_fisheye(num_transitions)
    return fisheye


def _access_auxevol_gf(name: str, i0: str = "i0", i1: str = "i1", i2: str = "i2") -> str:
    return gri.BHaHGridFunction.access_gf(name, 0, 0, 0, gf_array_name="auxevol_gfs").replace(
        "i0", i0
    ).replace("i1", i1).replace("i2", i2)


def _unique_assignments(assignments: List[Tuple[str, sp.Expr]]) -> List[Tuple[str, sp.Expr]]:
    seen: Dict[str, bool] = {}
    uniq: List[Tuple[str, sp.Expr]] = []
    for varname, expr in assignments:
        if varname in seen:
            continue
        seen[varname] = True
        uniq.append((varname, expr))
    return uniq


def register_CFunction_generalrfm_precompute(
    CoordSystem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register a C function that precomputes GeneralRFM quantities into AUXEVOL gridfunctions.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    parallelization = par.parval_from_str("parallelization")
    if parallelization == "cuda":
        raise ValueError(
            "GeneralRFM precompute to AUXEVOL gridfunctions does not yet support CUDA."
        )

    fisheye = _build_fisheye_exprs(CoordSystem)

    ghatDD = fisheye.ghatDD
    ghatDDdD = fisheye.ghatDDdD
    ghatDDdDD = fisheye.ghatDDdDD

    ghatUU, detgammahat = ixp.symm_matrix_inverter3x3(ghatDD)
    # Register gridfunctions (AUXEVOL) only if they are not already present.
    # ReferenceMetric (GeneralRFM branch) should have already registered these.
    if "ghatDD00" not in gri.glb_gridfcs_dict:
        _ = gri.register_gridfunctions_for_single_rankN(
            "ghatDD", rank=2, symmetry="sym01", group="AUXEVOL"
        )
    if "ghatUU00" not in gri.glb_gridfcs_dict:
        _ = gri.register_gridfunctions_for_single_rankN(
            "ghatUU", rank=2, symmetry="sym01", group="AUXEVOL"
        )
    if "detgammahat" not in gri.glb_gridfcs_dict:
        _ = gri.register_gridfunctions(["detgammahat"], group="AUXEVOL")
    if "ghatDDdD000" not in gri.glb_gridfcs_dict:
        _ = gri.register_gridfunctions_for_single_rankN(
            "ghatDDdD", rank=3, symmetry="sym01", group="AUXEVOL"
        )
    if "ghatDDdDD0000" not in gri.glb_gridfcs_dict:
        _ = gri.register_gridfunctions_for_single_rankN(
            "ghatDDdDD", rank=4, symmetry="sym01sym23", group="AUXEVOL"
        )
    assignments: List[Tuple[str, sp.Expr]] = []

    # ghatDD, ghatUU (sym01)
    for i in range(3):
        for j in range(i, 3):
            assignments.append((_access_auxevol_gf(f"ghatDD{i}{j}"), ghatDD[i][j]))
            assignments.append((_access_auxevol_gf(f"ghatUU{i}{j}"), ghatUU[i][j]))

    # detgammahat
    assignments.append((_access_auxevol_gf("detgammahat"), detgammahat))

    # ghatDD derivatives
    for i in range(3):
        for j in range(i, 3):
            for k in range(3):
                assignments.append(
                    (_access_auxevol_gf(f"ghatDDdD{i}{j}{k}"), ghatDDdD[i][j][k])
                )
                for l in range(k, 3):
                    assignments.append(
                        (
                            _access_auxevol_gf(f"ghatDDdDD{i}{j}{k}{l}"),
                            ghatDDdDD[i][j][k][l],
                        )
                    )

    assignments = _unique_assignments(assignments)
    expr_list = [expr for _, expr in assignments]
    var_list = [var for var, _ in assignments]

    # Parameter symbols
    # CodeParameters.h provides local definitions for these parameters.
    param_defs = ""

    fp_type = par.parval_from_str("fp_type")
    fp_type_alias = "DOUBLE" if fp_type == "float" else "REAL"

    ccode = ccg.c_codegen(
        expr_list,
        var_list,
        fp_type_alias=fp_type_alias,
        verbose=False,
        include_braces=False,
    )

    body = f"""
{param_defs}
for (int i2 = 0; i2 < params->Nxx_plus_2NGHOSTS2; i2++) {{
  for (int i1 = 0; i1 < params->Nxx_plus_2NGHOSTS1; i1++) {{
    for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {{
      const REAL xx0 = xx[0][i0];
      const REAL xx1 = xx[1][i1];
      const REAL xx2 = xx[2][i2];
{ccode}
    }}
  }}
}}
"""

    name = f"generalrfm_precompute__{CoordSystem}"
    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        desc=f"Precompute GeneralRFM quantities for {CoordSystem}.",
        cfunc_type="void",
        name=name,
        params="const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL *restrict xx[3], REAL *restrict auxevol_gfs",
        include_CodeParameters_h=True,
        body=body,
    )
    return pcg.NRPyEnv()


def register_CFunctions_generalrfm_support(CoordSystem: str) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register all GeneralRFM support C functions (precompute and Cart_to_xx).
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Precompute.
    register_CFunction_generalrfm_precompute(CoordSystem)

    # Cart_to_xx numeric inverse.
    from nrpy.infrastructures.BHaH import generalrfm_cart_to_xx

    generalrfm_cart_to_xx.register_CFunction_generalrfm_Cart_to_xx(CoordSystem)

    return pcg.NRPyEnv()
