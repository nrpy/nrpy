"""
C function registration for converting between grid coordinate (xx0,xx1,xx2) (uniform grid spacing) to Cartesian coordinate (x,y,z).

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.reference_metric as refmetric


# Construct Cart_to_xx_and_nearest_i0i1i2() C function for
# mapping from Cartesian->xx for the chosen CoordSystem.
def register_CFunction__Cart_to_xx_and_nearest_i0i1i2(
    CoordSystem: str,
    relative_to: str = "local_grid_center",
    gridding_approach: str = "independent grid(s)",
) -> None:
    """
    Construct and register a C function that maps Cartesian coordinates to xx and finds the nearest grid indices.

    This function generates a C function which, given Cartesian coordinates (x, y, z),
    computes the corresponding (xx0, xx1, xx2) coordinates and determines the "closest"
    grid indices (i0, i1, i2) for the specified coordinate system. The C function is
    then registered for later use.

    :param CoordSystem: The coordinate system for the local grid patch.
    :param relative_to: Whether the computation is relative to the "local_grid_center"
                        (default) or "global_grid_center".
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :raises ValueError: When the value of `gridding_approach` is not "independent grid(s)"
                        or "multipatch".
    """
    if gridding_approach not in {"independent grid(s)", "multipatch"}:
        raise ValueError(
            "Invalid value for 'gridding_approach'. Must be 'independent grid(s)' or 'multipatch'."
        )

    rfm = refmetric.reference_metric[CoordSystem]
    desc = "Given Cartesian point (x,y,z), this function "
    if gridding_approach == "multipatch":
        desc += "does stuff needed for multipatch, and then "
    desc += """unshifts the grid back to the origin to output the corresponding
            (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid"""

    namesuffix = f"_{relative_to}" if relative_to == "global_grid_center" else ""
    name = f"Cart_to_xx_and_nearest_i0i1i2{namesuffix}"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3], REAL xx[3], int Cart_to_i0i1i2[3]"

    body = """
  // Set (Cartx, Carty, Cartz) relative to the global (as opposed to local) grid.
  //   This local grid may be offset from the origin by adjusting
  //   (Cart_originx, Cart_originy, Cart_originz) to nonzero values.
  REAL Cartx = xCart[0];
  REAL Carty = xCart[1];
  REAL Cartz = xCart[2];
"""
    if relative_to == "local_grid_center":
        body += """
  // Set the origin, (Cartx, Carty, Cartz) = (0, 0, 0), to the center of the local grid patch.
  Cartx -= Cart_originx;
  Carty -= Cart_originy;
  Cartz -= Cart_originz;
  {
"""
    if rfm.requires_NewtonRaphson_for_Cart_to_xx:
        body += "  // First compute analytical coordinate inversions:\n"
        Cart_to_xx_exprs: List[sp.Expr] = []
        Cart_to_xx_names: List[str] = []
        for i in range(3):
            if rfm.NewtonRaphson_f_of_xx[i] == sp.sympify(0):
                Cart_to_xx_exprs += [rfm.Cart_to_xx[i]]
                Cart_to_xx_names += [f"xx[{i}]"]
        body += ccg.c_codegen(Cart_to_xx_exprs, Cart_to_xx_names, include_braces=False)
        body += """
  // Next perform Newton-Raphson iterations as needed:
  const REAL XX_TOLERANCE = 1e-12;  // that's 1 part in 1e12 dxxi.
  const REAL F_OF_XX_TOLERANCE = 1e-12;  // tolerance of function for which we're finding the root.
  const int ITER_MAX = 100;
  int iter, tolerance_has_been_met=0;
"""
        for i in range(3):
            if rfm.NewtonRaphson_f_of_xx[i] != sp.sympify(0):
                body += f"""
  iter=0;
  REAL xx{i}  = 0.5 * (params->xxmin{i} + params->xxmax{i});
  while(iter < ITER_MAX && !tolerance_has_been_met) {{
    REAL f_of_xx{i}, fprime_of_xx{i};

{ccg.c_codegen([rfm.NewtonRaphson_f_of_xx[i], sp.diff(rfm.NewtonRaphson_f_of_xx[i], rfm.xx[i])],
[f'f_of_xx{i}', f'fprime_of_xx{i}'], include_braces=True, verbose=False)}
    const REAL xx{i}_np1 = xx{i} - f_of_xx{i} / fprime_of_xx{i};

    if( fabs(xx{i} - xx{i}_np1) <= XX_TOLERANCE * params->dxx{i} && fabs(f_of_xx{i}) <= F_OF_XX_TOLERANCE ) {{
      tolerance_has_been_met = 1;
    }}
    xx{i} = xx{i}_np1;
    iter++;
  }} // END Newton-Raphson iterations to compute xx{i}
  if(iter >= ITER_MAX) {{
    fprintf(stderr, "ERROR: Newton-Raphson failed for {CoordSystem}: xx{i}, x,y,z = %.15e %.15e %.15e\\n", Cartx,Carty,Cartz);
    exit(1);
  }}
  xx[{i}] = xx{i};
"""
    else:
        body += ccg.c_codegen(
            [rfm.Cart_to_xx[0], rfm.Cart_to_xx[1], rfm.Cart_to_xx[2]],
            ["xx[0]", "xx[1]", "xx[2]"],
            include_braces=False,
        )

    body += """
      // Find the nearest grid indices (i0, i1, i2) for the given Cartesian coordinates (x, y, z).
      // Assuming a cell-centered grid, which follows the pattern:
      //   xx0[i0] = params->xxmin0 + ((REAL)(i0 - NGHOSTS) + 0.5) * params->dxx0
      // The index i0 can be derived as:
      //   i0 = (xx0[i0] - params->xxmin0) / params->dxx0 - 0.5 + NGHOSTS
      // Now, including typecasts:
      //   i0 = (int)((xx[0] - params->xxmin0) / params->dxx0 - 0.5 + (REAL)NGHOSTS)
      // The (int) typecast always rounds down, so we add 0.5 inside the outer parenthesis:
      //   i0 = (int)((xx[0] - params->xxmin0) / params->dxx0 - 0.5 + (REAL)NGHOSTS + 0.5)
      // The 0.5 values cancel out:
      //   i0 =           (int)( ( xx[0] - params->xxmin0 ) / params->dxx0 + (REAL)NGHOSTS )
      Cart_to_i0i1i2[0] = (int)( ( xx[0] - params->xxmin0 ) / params->dxx0 + (REAL)NGHOSTS );
      Cart_to_i0i1i2[1] = (int)( ( xx[1] - params->xxmin1 ) / params->dxx1 + (REAL)NGHOSTS );
      Cart_to_i0i1i2[2] = (int)( ( xx[2] - params->xxmin2 ) / params->dxx2 + (REAL)NGHOSTS );
  }
"""
    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        desc=desc,
        cfunc_type="void",
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )


def register_CFunction_xx_to_Cart(
    CoordSystem: str,
    gridding_approach: str = "independent grid(s)",
) -> None:
    """
    Convert uniform-grid coordinate (xx[0], xx[1], xx[2]) to the corresponding Cartesian coordinate.

    :param CoordSystem: The coordinate system name as a string.
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".

    :raises ValueError: If an invalid gridding_approach is provided.
    """
    if gridding_approach not in {"independent grid(s)", "multipatch"}:
        raise ValueError(
            "Invalid value for 'gridding_approach'. Must be 'independent grid(s)' or 'multipatch'."
        )

    # Description of the conversion process
    desc = """Compute Cartesian coordinates {x, y, z} = {xCart[0], xCart[1], xCart[2]} in terms of
              local grid coordinates {xx[0][i0], xx[1][i1], xx[2][i2]} = {xx0, xx1, xx2},
              taking into account the possibility that the origin of this grid is off-center."""

    cfunc_type = "void"
    name = "xx_to_Cart"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3]"""

    rfm = refmetric.reference_metric[CoordSystem]

    # ** Code body for the conversion process **
    # Suppose grid origin is at (1,1,1). Then the Cartesian gridpoint at (1,2,3) will be (2,3,4);
    # hence the xx_to_Cart[i] + gri.Cart_origin[i] below:
    body = """
const REAL xx0 = xx[0][i0];
const REAL xx1 = xx[1][i1];
const REAL xx2 = xx[2][i2];
""" + ccg.c_codegen(
        [
            rfm.xx_to_Cart[0] + gri.Cart_origin[0],
            rfm.xx_to_Cart[1] + gri.Cart_origin[1],
            rfm.xx_to_Cart[2] + gri.Cart_origin[2],
        ],
        ["xCart[0]", "xCart[1]", "xCart[2]"],
    )

    # Register the C function with the provided details
    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
