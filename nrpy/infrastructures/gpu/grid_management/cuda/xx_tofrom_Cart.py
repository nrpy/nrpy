"""
C function registration for converting between grid coordinate (xx0,xx1,xx2) (uniform grid spacing) to Cartesian coordinate (x,y,z).

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_codegen as ccg
import nrpy.infrastructures.gpu.grid_management.base_xx_tofrom_Cart as base_xx_classes


# Construct Cart_to_xx_and_nearest_i0i1i2() C function for
# mapping from Cartesian->xx for the chosen CoordSystem.
class register_CFunction__Cart_to_xx_and_nearest_i0i1i2(
    base_xx_classes.base_register_CFunction__Cart_to_xx_and_nearest_i0i1i2
):
    """
    Construct a C function that maps from Cartesian coordinates to xx for the chosen coordinate system.
    Registers the C function for later use.

    :param CoordSystem: The coordinate system to use.
    :param relative_to: Whether the computation is relative to the "local_grid_center" or "global_grid_center".
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :param fp_type: Floating point type, e.g., "double".
    :raises ValueError: When the value of `gridding_approach` is not "independent grid(s)"
                        or "multipatch".
    """

    def __init__(
        self,
        CoordSystem: str,
        relative_to: str = "local_grid_center",
        fp_type: str = "double",
    ) -> None:
        super().__init__(
            CoordSystem,
            relative_to=relative_to,
            fp_type=fp_type,
        )

        self.body += """
  // Set (Cartx, Carty, Cartz) relative to the global (as opposed to local) grid.
  //   This local grid may be offset from the origin by adjusting
  //   (Cart_originx, Cart_originy, Cart_originz) to nonzero values.
  REAL Cartx = xCart[0];
  REAL Carty = xCart[1];
  REAL Cartz = xCart[2];
"""
        if self.relative_to == "local_grid_center":
            self.body += """
  // Set the origin, (Cartx, Carty, Cartz) = (0, 0, 0), to the center of the local grid patch.
  Cartx -= Cart_originx;
  Carty -= Cart_originy;
  Cartz -= Cart_originz;
  {
"""

        if "theta_adj" in self.CoordSystem:
            self.body += ccg.c_codegen(
                [
                    self.rfm.Cart_to_xx[0],
                    self.rfm.Cart_to_xx[1],
                    self.rfm.Cart_to_xx[2],
                ],
                ["xx[0]", "const REAL target_th", "xx[2]"],
                include_braces=False,
                fp_type=self.fp_type,
            )
            self.body += "xx[1] = NewtonRaphson_get_xx1_from_th(params, target_th);\n"
        else:
            self.body += ccg.c_codegen(
                [
                    self.rfm.Cart_to_xx[0],
                    self.rfm.Cart_to_xx[1],
                    self.rfm.Cart_to_xx[2],
                ],
                ["xx[0]", "xx[1]", "xx[2]"],
                include_braces=False,
                fp_type=self.fp_type,
            )

        self.body += f"""
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
  Cart_to_i0i1i2[0] = (int)( ( xx[0] - ({self.rfm.xxmin[0]})) / params->dxx0 + 0.5 + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
  Cart_to_i0i1i2[1] = (int)( ( xx[1] - ({self.rfm.xxmin[1]})) / params->dxx1 + 0.5 + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
  Cart_to_i0i1i2[2] = (int)( ( xx[2] - ({self.rfm.xxmin[2]})) / params->dxx2 + 0.5 + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
"""

        self.register()


class register_CFunction_xx_to_Cart(base_xx_classes.base_register_CFunction_xx_to_Cart):
    """
    Convert uniform-grid coordinate (xx[0], xx[1], xx[2]) to the corresponding Cartesian coordinate.

    :param CoordSystem: The coordinate system name as a string.
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :param fp_type: Floating point type, e.g., "double".

    :raises ValueError: If an invalid gridding_approach is provided.
    """

    def __init__(
        self,
        CoordSystem: str,
        gridding_approach: str = "independent grid(s)",
        fp_type: str = "double",
    ) -> None:
        super().__init__(
            CoordSystem,
            gridding_approach=gridding_approach,
            fp_type=fp_type,
        )
        self.cfunc_type = "void"
        self.cfunc_decorators = "__host__"
        kernel_body = ""
        for sym in self.unique_symbols:
            kernel_body += f"const REAL {sym} = params->{sym};\n"
        self.body = kernel_body + self.body

        self.register()
