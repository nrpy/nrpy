"""
C function registration for converting between grid coordinate (xx0,xx1,xx2) (uniform grid spacing) to Cartesian coordinate (x,y,z).

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.reference_metric as refmetric
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri


# Construct Cart_to_xx_and_nearest_i0i1i2() C function for
# mapping from Cartesian->xx for the chosen CoordSystem.
def register_CFunction__Cart_to_xx_and_nearest_i0i1i2(
    CoordSystem: str, relative_to: str = "local_grid_center"
) -> None:
    """
    Construct a C function that maps from Cartesian coordinates to xx for the chosen coordinate system.
    Registers the C function for later use.

    :param CoordSystem: The coordinate system to use.
    :param relative_to: Whether the computation is relative to the "local_grid_center" or "global_grid_center".

    :raises ValueError: When the value of `relative_to` is not "local_grid_center" or "global_grid_center".
    """
    rfm = refmetric.reference_metric[CoordSystem]
    prefunc = ""
    desc = """Given Cartesian point (x,y,z), this function outputs the corresponding
  (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid"""

    namesuffix = f"_{relative_to}" if relative_to == "global_grid_center" else ""
    name = f"Cart_to_xx_and_nearest_i0i1i2{namesuffix}"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3], REAL xx[3], int Cart_to_i0i1i2[3]"

    body = ""
    if relative_to == "local_grid_center":
        body += """
  // See comments for description on how coordinates are computed relative to the local grid center.
  const REAL Cartx = xCart[0] - Cart_originx;
  const REAL Carty = xCart[1] - Cart_originy;
  const REAL Cartz = xCart[2] - Cart_originz;
"""
    elif relative_to == "global_grid_center":
        body += """
  // Coordinates are global and no transformation is needed.
  const REAL Cartx = xCart[0];
  const REAL Carty = xCart[1];
  const REAL Cartz = xCart[2];
"""
    else:
        raise ValueError(
            f"Error: relative_to must be set to either local_grid_center or global_grid_center.\n"
            f"{relative_to} was chosen."
        )

    if "theta_adj" in CoordSystem:
        body += ccg.c_codegen(
            [rfm.Cart_to_xx[0], rfm.Cart_to_xx[1], rfm.Cart_to_xx[2]],
            ["xx[0]", "const REAL target_th", "xx[2]"],
            include_braces=False,
        )
        body += "xx[1] = NewtonRaphson_get_xx1_from_th(params, target_th);\n"
    else:
        body += ccg.c_codegen(
            [rfm.Cart_to_xx[0], rfm.Cart_to_xx[1], rfm.Cart_to_xx[2]],
            ["xx[0]", "xx[1]", "xx[2]"],
            include_braces=False,
        )

    body += f"""
  // Then find the nearest index (i0,i1,i2) on underlying grid to (x,y,z)
  Cart_to_i0i1i2[0] = (int)( ( xx[0] - ({rfm.xxmin[0]})) / params->dxx0 + 0.5 + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
  Cart_to_i0i1i2[1] = (int)( ( xx[1] - ({rfm.xxmin[1]})) / params->dxx1 + 0.5 + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
  Cart_to_i0i1i2[2] = (int)( ( xx[2] - ({rfm.xxmin[2]})) / params->dxx2 + 0.5 + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
"""

    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        prefunc=prefunc,
        desc=desc,
        c_type="void",
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )


def register_CFunction_xx_to_Cart(CoordSystem: str) -> None:
    """
    Register a C function to convert arbitrary NRPy+ coordinates to Cartesian coordinates.

    :param CoordSystem: The coordinate system name as a string.
    """
    # Arbitrary-coordinate NRPy+ file output, Part 1: output the conversion from (x0,x1,x2) to Cartesian (x,y,z)
    # Suppose grid origin is at 1,1,1. Then the Cartesian gridpoint at 1,2,3 will be 2,3,4; hence
    # the xx_to_Cart[i]+gri.Cart_origin[i] below:
    rfm = refmetric.reference_metric[CoordSystem]

    body = """
REAL xx0 = xx[0][i0];
REAL xx1 = xx[1][i1];
REAL xx2 = xx[2][i2];
""" + ccg.c_codegen(
        [
            rfm.xx_to_Cart[0] + gri.Cart_origin[0],
            rfm.xx_to_Cart[1] + gri.Cart_origin[1],
            rfm.xx_to_Cart[2] + gri.Cart_origin[2],
        ],
        ["xCart[0]", "xCart[1]", "xCart[2]"],
    )

    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        desc="Compute Cartesian coordinates given local grid coordinate (xx0,xx1,xx2), "
        "  accounting for the origin of this grid being possibly off-center.",
        c_type="void",
        CoordSystem_for_wrapper_func=CoordSystem,
        name="xx_to_Cart",
        params="const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3]",
        include_CodeParameters_h=True,
        body=body,
    )
