"""
Register a C helper that numerically inverts GeneralRFM xx->Cart maps.

Currently supports GeneralRFM_fisheyeN* via equations/generalrfm/fisheye.py.

Author: Nishita Jadoo
        njadoo@uidaho.edu
"""

from __future__ import annotations

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.infrastructures.BHaH.xx_tofrom_Cart import (
    _generate_bracketed_radial_inverse_body,
)


def register_CFunction_generalrfm_Cart_to_xx(
    CoordSystem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register a C function to numerically invert the GeneralRFM map Cart -> xx.

    :param CoordSystem: GeneralRFM coordinate system name.
    :raises ValueError: If CoordSystem is not GeneralRFM or its provider is unsupported/missing.
    :return: None during parallel-codegen registration phase, otherwise the NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    if par.parval_from_str("parallelization") == "cuda":
        raise ValueError("GeneralRFM Cart_to_xx does not support CUDA parallelization.")

    rfm = refmetric.reference_metric[CoordSystem]
    if not CoordSystem.startswith("GeneralRFM"):
        raise ValueError(f"{CoordSystem} is not a GeneralRFM coordinate system.")
    provider_name = getattr(rfm, "general_rfm_provider_name", "")
    if provider_name != "fisheye":
        raise ValueError(
            f"GeneralRFM provider '{provider_name}' for {CoordSystem} is not yet supported in generalrfm_Cart_to_xx."
        )
    fisheye = getattr(rfm, "general_rfm_provider", None)
    if fisheye is None:
        raise ValueError(f"GeneralRFM provider object missing for {CoordSystem}.")

    r_local = sp.Symbol("r", real=True, nonnegative=True)
    rbar_expr, drbar_dr_expr, _, _ = fisheye.radius_map_and_derivs_for_inverse(r_local)
    asymptotic_scale_expr = fisheye.c * fisheye.a_list[-1]

    origin_body = """    xx[0] = (REAL)0.0;
    xx[1] = (REAL)0.0;
    xx[2] = (REAL)0.0;
    return 0;"""
    success_body = """    const REAL inv_rCart = (REAL)1.0 / rCart;
    xx[0] = Cart[0] * radial_seed * inv_rCart;
    xx[1] = Cart[1] * radial_seed * inv_rCart;
    xx[2] = Cart[2] * radial_seed * inv_rCart;
    return 0;"""
    failure_body = """      return 1;"""
    body = _generate_bracketed_radial_inverse_body(
        r_local,
        rbar_expr,
        drbar_dr_expr,
        asymptotic_scale_expr,
        ("Cart[0]", "Cart[1]", "Cart[2]"),
        origin_body,
        success_body,
        failure_body,
    )

    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        desc=f"Numerically invert Cart->xx for {CoordSystem}.",
        cfunc_type="int",
        name=f"generalrfm_Cart_to_xx__{CoordSystem}",
        params="const params_struct *restrict params, const REAL Cart[3], REAL xx[3]",
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
