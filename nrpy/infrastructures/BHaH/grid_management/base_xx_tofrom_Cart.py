"""
Base classes for Coordinate conversions

The parallelization modules will generate functions 
to convert between grid coordinate (xx0,xx1,xx2) (uniform grid spacing) 
to Cartesian coordinate (x,y,z) within the BHaH infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Samuel D. Tootle
        sdtootle **at** gmail **dot* com        
"""

import nrpy.reference_metric as refmetric
import nrpy.c_codegen as ccg
import nrpy.grid as gri


# Construct Cart_to_xx_and_nearest_i0i1i2() C function for
# mapping from Cartesian->xx for the chosen CoordSystem.
class base_register_CFunction__Cart_to_xx_and_nearest_i0i1i2:
    def __init__(
        self,
        CoordSystem: str,
        relative_to: str = "local_grid_center",
        fp_type: str = "double",
    ) -> None:
        """
        Base to generate the function that maps from Cartesian coordinates to xx for the chosen coordinate system.
        Registers the C function for later use.

        :param CoordSystem: The coordinate system to use.
        :param relative_to: Whether the computation is relative to the "local_grid_center" or "global_grid_center".
        :param fp_type: Floating point type, e.g., "double".

        :raises ValueError: When the value of `relative_to` is not "local_grid_center" or "global_grid_center".
        """
        self.CoordSystem=CoordSystem
        self.relative_to=relative_to
        self.fp_type=fp_type
        self.includes=["BHaH_defines.h"]
        self.cfunc_type="void"
        self.prefunc = ""
        self.desc = """Given Cartesian point (x,y,z), this function outputs the corresponding
    (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid"""
        
        self.namesuffix = f"_{self.relative_to}" if self.relative_to == "global_grid_center" else ""
        self.name = f"Cart_to_xx_and_nearest_i0i1i2{self.namesuffix}"
        self.params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3], REAL xx[3], int Cart_to_i0i1i2[3]"
        self.rfm = refmetric.reference_metric[self.CoordSystem]        

        self.body = ""


class base_register_CFunction_xx_to_Cart:
    
    def __init__ (
        self,
        CoordSystem: str, 
        fp_type: str = "double"
    ) -> None:
        """
        Base for generating the function to convert arbitrary NRPy+ coordinates to Cartesian coordinates.

        :param CoordSystem: The coordinate system name as a string.
        :param fp_type: Floating point type, e.g., "double".
        """
        self.CoordSystem=CoordSystem
        self.fp_type=fp_type
        self.includes=["BHaH_defines.h"]
        self.cfunc_type="void"
        self.name="xx_to_Cart"
        self.params="const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3]"
        self.descr="Compute Cartesian coordinates given local grid coordinate (xx0,xx1,xx2), " \
            "  accounting for the origin of this grid being possibly off-center."
        
        # Arbitrary-coordinate NRPy+ file output, Part 1: output the conversion from (x0,x1,x2) to Cartesian (x,y,z)
        # Suppose grid origin is at 1,1,1. Then the Cartesian gridpoint at 1,2,3 will be 2,3,4; hence
        # the xx_to_Cart[i]+gri.Cart_origin[i] below:
        self.rfm = refmetric.reference_metric[self.CoordSystem]

        self.body = """
REAL xx0 = xx[0][i0];
REAL xx1 = xx[1][i1];
REAL xx2 = xx[2][i2];
""" + ccg.c_codegen(
        [
            self.rfm.xx_to_Cart[0] + gri.Cart_origin[0],
            self.rfm.xx_to_Cart[1] + gri.Cart_origin[1],
            self.rfm.xx_to_Cart[2] + gri.Cart_origin[2],
        ],
        ["xCart[0]", "xCart[1]", "xCart[2]"],
        fp_type=self.fp_type,
        )
