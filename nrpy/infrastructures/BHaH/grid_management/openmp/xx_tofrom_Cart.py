"""
C function registration for converting between grid coordinate (xx0,xx1,xx2) (uniform grid spacing) to Cartesian coordinate (x,y,z).

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.grid_management.base_xx_tofrom_Cart as base_xx_classes


# Construct Cart_to_xx_and_nearest_i0i1i2() C function for
# mapping from Cartesian->xx for the chosen CoordSystem.
class register_CFunction__Cart_to_xx_and_nearest_i0i1i2(
    base_xx_classes.base_register_CFunction__Cart_to_xx_and_nearest_i0i1i2
):
    def __init__(
        self,
        CoordSystem: str,
        relative_to: str = "local_grid_center",
        fp_type: str = "double",
    ) -> None:
        """
        Construct a C function that maps from Cartesian coordinates to xx for the chosen coordinate system.
        Registers the C function for later use.

        :param CoordSystem: The coordinate system to use.
        :param relative_to: Whether the computation is relative to the "local_grid_center" or "global_grid_center".
        :param fp_type: Floating point type, e.g., "double".

        :raises ValueError: When the value of `relative_to` is not "local_grid_center" or "global_grid_center".
        """
        super().__init__(
            CoordSystem,
            relative_to=relative_to,
            fp_type=fp_type,
        )

        if self.relative_to == "local_grid_center":
            self.body += """
  // See comments for description on how coordinates are computed relative to the local grid center.
  const REAL Cartx = xCart[0] - Cart_originx;
  const REAL Carty = xCart[1] - Cart_originy;
  const REAL Cartz = xCart[2] - Cart_originz;
"""
        elif self.relative_to == "global_grid_center":
            self.body += """
  // Coordinates are global and no transformation is needed.
  const REAL Cartx = xCart[0];
  const REAL Carty = xCart[1];
  const REAL Cartz = xCart[2];
"""
        else:
            raise ValueError(
                f"Error: relative_to must be set to either local_grid_center or global_grid_center.\n"
                f"{self.relative_to} was chosen."
            )

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
  // Then find the nearest index (i0,i1,i2) on underlying grid to (x,y,z)
  Cart_to_i0i1i2[0] = (int)( ( xx[0] - ({self.rfm.xxmin[0]})) / params->dxx0 + 0.5 + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
  Cart_to_i0i1i2[1] = (int)( ( xx[1] - ({self.rfm.xxmin[1]})) / params->dxx1 + 0.5 + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
  Cart_to_i0i1i2[2] = (int)( ( xx[2] - ({self.rfm.xxmin[2]})) / params->dxx2 + 0.5 + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
"""

        cfc.register_CFunction(
            includes=self.includes,
            prefunc=self.prefunc,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            CoordSystem_for_wrapper_func=self.CoordSystem,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=True,
            body=self.body,
        )


class register_CFunction_xx_to_Cart(base_xx_classes.base_register_CFunction_xx_to_Cart):

    def __init__(self, CoordSystem: str, fp_type: str = "double") -> None:
        """
        Register a C function to convert arbitrary NRPy+ coordinates to Cartesian coordinates.

        :param CoordSystem: The coordinate system name as a string.
        :param fp_type: Floating point type, e.g., "double".
        """
        super().__init__(
            CoordSystem,
            fp_type=fp_type,
        )

        cfc.register_CFunction(
            includes=self.includes,
            desc=self.descr,
            cfunc_type=self.cfunc_type,
            CoordSystem_for_wrapper_func=self.CoordSystem,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=True,
            body=self.body,
        )
