"""
Base classes for Coordinate conversions.

The parallelization modules will generate functions
to convert between grid coordinate (xx0,xx1,xx2) (uniform grid spacing)
to Cartesian coordinate (x,y,z) within the BHaH infrastructure.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot** com
        Samuel D. Tootle
        sdtootle **at** gmail **dot** com
"""

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.reference_metric as refmetric
from nrpy.helpers.expr_tree import get_unique_expression_symbols


# Construct Cart_to_xx_and_nearest_i0i1i2() C function for
# mapping from Cartesian->xx for the chosen CoordSystem.
class base_register_CFunction__Cart_to_xx_and_nearest_i0i1i2:
    """
    Base to generate the C function that maps Cartesian coordinates to xx and finds the nearest grid indices.
    This function generates a C function which, given Cartesian coordinates (x, y, z),
    computes the corresponding (xx0, xx1, xx2) coordinates and determines the "closest"
    grid indices (i0, i1, i2) for the specified coordinate system. The C function is
    then registered for later use.

    :param CoordSystem: The coordinate system for the local grid patch.
    :param relative_to: Whether the computation is relative to the "local_grid_center"
                        (default) or "global_grid_center".
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :param fp_type: Floating point type, e.g., "double".
    :raises ValueError: When the value of `gridding_approach` is not "independent grid(s)"
                        or "multipatch".
    """

    def __init__(
        self,
        CoordSystem: str,
        relative_to: str = "local_grid_center",
        gridding_approach: str = "independent grid(s)",
        fp_type: str = "double",
    ) -> None:
        self.CoordSystem = CoordSystem
        self.relative_to = relative_to
        self.gridding_approach = gridding_approach
        self.fp_type = fp_type
        self.includes = ["BHaH_defines.h"]
        self.cfunc_type = "void"
        self.prefunc = ""
        self.desc = """Given Cartesian point (x,y,z), this function outputs the corresponding
    (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid"""

        self.namesuffix = (
            f"_{self.relative_to}" if self.relative_to == "global_grid_center" else ""
        )
        self.name = f"Cart_to_xx_and_nearest_i0i1i2{self.namesuffix}"
        self.params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3], REAL xx[3], int Cart_to_i0i1i2[3]"
        self.rfm = refmetric.reference_metric[self.CoordSystem]

        self.body = ""

    def register(self) -> None:
        """Register CFunction."""
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


class base_register_CFunction_xx_to_Cart:
    """
    Base for generating conversion from uniform-grid coordinate (xx[0], xx[1], xx[2]) to the corresponding Cartesian coordinate.

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
        self.gridding_approach = gridding_approach
        self.CoordSystem = CoordSystem
        self.fp_type = fp_type
        self.includes = ["BHaH_defines.h"]
        self.cfunc_type = "void"
        self.name = "xx_to_Cart"
        self.params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3]"
        self.desc = """Compute Cartesian coordinates {x, y, z} = {xCart[0], xCart[1], xCart[2]} in terms of
              local grid coordinates {xx[0][i0], xx[1][i1], xx[2][i2]} = {xx0, xx1, xx2},
              taking into account the possibility that the origin of this grid is off-center."""
        self.cfunc_decorators = ""
        if gridding_approach not in {"independent grid(s)", "multipatch"}:
            raise ValueError(
                "Invalid value for 'gridding_approach'. Must be 'independent grid(s)' or 'multipatch'."
            )

        # Arbitrary-coordinate NRPy+ file output, Part 1: output the conversion from (x0,x1,x2) to Cartesian (x,y,z)
        # Suppose grid origin is at 1,1,1. Then the Cartesian gridpoint at 1,2,3 will be 2,3,4; hence
        # the xx_to_Cart[i]+gri.Cart_origin[i] below:
        self.rfm = refmetric.reference_metric[self.CoordSystem]
        expr_list = [
            self.rfm.xx_to_Cart[0] + gri.Cart_origin[0],
            self.rfm.xx_to_Cart[1] + gri.Cart_origin[1],
            self.rfm.xx_to_Cart[2] + gri.Cart_origin[2],
        ]
        self.unique_symbols = []
        for expr in expr_list:
            self.unique_symbols += get_unique_expression_symbols(
                expr, exclude=[f"xx{i}" for i in range(3)]
            )
        self.unique_symbols = sorted(list(set(self.unique_symbols)))
        self.body = """
REAL xx0 = xx[0][i0];
REAL xx1 = xx[1][i1];
REAL xx2 = xx[2][i2];
""" + ccg.c_codegen(
            expr_list,
            ["xCart[0]", "xCart[1]", "xCart[2]"],
            fp_type=self.fp_type,
        )

    def register(self) -> None:
        """Register CFunction."""
        cfc.register_CFunction(
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            CoordSystem_for_wrapper_func=self.CoordSystem,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=False,
            body=self.body,
            cfunc_decorators=self.cfunc_decorators,
        )
