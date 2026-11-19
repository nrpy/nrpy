"""
Register a C helper for GeneralRFM physical-to-computational conversion.

Supports forward fisheye and inverse spheroidal fisheye GeneralRFM providers.

Author: Nishita Jadoo
        njadoo@uidaho.edu
"""

from __future__ import annotations

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.infrastructures.BHaH import inverse_fisheye
from nrpy.infrastructures.BHaH.xx_tofrom_Cart import (
    _generate_bracketed_radial_inverse_body,
)


def register_CFunction_generalrfm_Cart_to_xx(
    CoordSystem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register a C function for the GeneralRFM map Cart -> xx.

    :param CoordSystem: GeneralRFM coordinate system name.
    :raises ValueError: If CoordSystem is not GeneralRFM, its provider is
                        unsupported or missing, or CUDA parallelization is active.
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
    if provider_name == "inverse_fisheye":
        n = rfm.general_rfm_provider_meta["num_transitions"]
        inverse_fisheye.register_CFunctions(n)
        body = f"""
if (!isfinite(Cart[0]) || !isfinite(Cart[1]) || !isfinite(Cart[2])) return 1;
REAL lambda;
inverse_fisheye_scale__N{n}(params, Cart, 0, &lambda, NULL, NULL, NULL);
for (int i=0; i<3; i++) xx[i]=lambda*Cart[i];
return !(isfinite(xx[0]) && isfinite(xx[1]) && isfinite(xx[2]));
"""
        cfc.register_CFunction(
            includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
            desc=f"Evaluate the explicit physical-to-computational map for {CoordSystem}.",
            cfunc_type="int",
            name=f"generalrfm_Cart_to_xx__{CoordSystem}",
            params="const params_struct *restrict params, const REAL Cart[3], REAL xx[3]",
            include_CodeParameters_h=False,
            body=body,
        )
        return pcg.NRPyEnv()
    if provider_name not in {"fisheye", "spheroidal_fisheye"}:
        raise ValueError(
            f"GeneralRFM provider '{provider_name}' for {CoordSystem} is not yet supported in generalrfm_Cart_to_xx."
        )
    fisheye = getattr(rfm, "general_rfm_provider", None)
    if fisheye is None:
        raise ValueError(f"GeneralRFM provider object missing for {CoordSystem}.")

    origin_body = r"""    xx[0] = (REAL)0.0;
    xx[1] = (REAL)0.0;
    xx[2] = (REAL)0.0;
    return 0;"""
    failure_body = r"""      return 1;"""
    if provider_name == "fisheye":
        success_body = r"""    const REAL inv_rCart = (REAL)1.0 / rCart;
    xx[0] = Cart[0] * radial_seed * inv_rCart;
    xx[1] = Cart[1] * radial_seed * inv_rCart;
    xx[2] = Cart[2] * radial_seed * inv_rCart;
    return 0;"""
    else:
        success_body = r"""    xx[0] = Cart[0] / lam_xy;
    xx[1] = Cart[1] / lam_xy;
    xx[2] = Cart[2] / lam_z;
    return 0;"""

    body = _generate_bracketed_radial_inverse_body(
        fisheye,
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
