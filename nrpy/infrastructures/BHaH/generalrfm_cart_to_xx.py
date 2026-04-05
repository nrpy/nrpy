# nrpy/infrastructures/BHaH/generalrfm_cart_to_xx.py
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

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
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

    xx_to_Cart = fisheye.xx_to_CartU
    dCart_dxx = fisheye.dCart_dxxUD

    cart_exprs = [xx_to_Cart[0], xx_to_Cart[1], xx_to_Cart[2]]
    jac_exprs = [dCart_dxx[i][j] for i in range(3) for j in range(3)]

    param_symbols, _ = get_params_commondata_symbols_from_expr_list(
        cart_exprs + jac_exprs, exclude=[f"xx{j}" for j in range(3)]
    )
    param_defs = generate_definition_header(param_symbols, var_access="params->")

    fp_type = par.parval_from_str("fp_type")
    fp_type_alias = "DOUBLE" if fp_type == "float" else "REAL"

    cart_assigns = ccg.c_codegen(
        cart_exprs,
        ["Cart_guess[0]", "Cart_guess[1]", "Cart_guess[2]"],
        fp_type_alias=fp_type_alias,
        verbose=False,
        include_braces=False,
        cse_varprefix="cart_",
    )
    jac_assigns = ccg.c_codegen(
        jac_exprs,
        [f"J[{i}][{j}]" for i in range(3) for j in range(3)],
        fp_type_alias=fp_type_alias,
        verbose=False,
        include_braces=False,
        cse_varprefix="jac_",
    )

    body = f"""
{param_defs}
  // Newton solve for xx given Cart. Start with xx = Cart.
  REAL xx0 = Cart[0];
  REAL xx1 = Cart[1];
  REAL xx2 = Cart[2];
  const int max_iters = 50;
  const REAL tol = 1.0e-12;
  bool converged = false;

  for (int iter = 0; iter < max_iters; iter++) {{
    REAL Cart_guess[3];
    {cart_assigns}
    REAL F[3] = {{Cart_guess[0] - Cart[0], Cart_guess[1] - Cart[1], Cart_guess[2] - Cart[2]}};
    const REAL err = fabs(F[0]) + fabs(F[1]) + fabs(F[2]);
    if (err < tol) {{
      converged = true;
      break;
    }}

    REAL J[3][3];
    {jac_assigns}

    const REAL a = J[0][0], b = J[0][1], c = J[0][2];
    const REAL d = J[1][0], e = J[1][1], f = J[1][2];
    const REAL g = J[2][0], h = J[2][1], i = J[2][2];
    const REAL det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
    if (fabs(det) < 1.0e-14) {{
      return 1;
    }}
    const REAL invdet = 1.0 / det;
    REAL invJ[3][3];
    invJ[0][0] =  (e*i - f*h) * invdet;
    invJ[0][1] = -(b*i - c*h) * invdet;
    invJ[0][2] =  (b*f - c*e) * invdet;
    invJ[1][0] = -(d*i - f*g) * invdet;
    invJ[1][1] =  (a*i - c*g) * invdet;
    invJ[1][2] = -(a*f - c*d) * invdet;
    invJ[2][0] =  (d*h - e*g) * invdet;
    invJ[2][1] = -(a*h - b*g) * invdet;
    invJ[2][2] =  (a*e - b*d) * invdet;

    const REAL dxx0 = invJ[0][0]*F[0] + invJ[0][1]*F[1] + invJ[0][2]*F[2];
    const REAL dxx1 = invJ[1][0]*F[0] + invJ[1][1]*F[1] + invJ[1][2]*F[2];
    const REAL dxx2 = invJ[2][0]*F[0] + invJ[2][1]*F[1] + invJ[2][2]*F[2];

    xx0 -= dxx0;
    xx1 -= dxx1;
    xx2 -= dxx2;
  }}

  if (!converged) {{
    return 1;
  }}

  xx[0] = xx0;
  xx[1] = xx1;
  xx[2] = xx2;
  return 0;
"""

    name = f"generalrfm_Cart_to_xx__{CoordSystem}"
    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        desc=f"Numerically invert Cart->xx for {CoordSystem}.",
        cfunc_type="int",
        name=name,
        params="const params_struct *restrict params, const REAL Cart[3], REAL xx[3]",
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
