"""
Library of C functions for solving the hyperbolic relaxation equation in curvilinear coordinates, using a reference-metric formalism.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
         Nishita Jadoo; njadoo **at** uidaho **dot* edu

"""

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.reference_metric as refmetric


# Define function to compute the l^2 of a gridfunction between 2 radii
def register_CFunction_L2norm_gf_between_r1_and_r2(
    CoordSystem: str,
) -> None:
    """
    Register function to compute l2-norm of a gridfunction between 2 radii assuming a single grid.

    Note that parallel codegen is disabled for this function, as it sometimes causes a
    multiprocess race condition on Python 3.6.7

    :param CoordSystem: the rfm coordinate system.
    """
    includes = ["BHaH_defines.h"]
    desc = "Compute l2-norm of a gridfunction between 2 radii  assuming a single grid."
    cfunc_type = "void"
    name = "L2norm_gf_between_r1_and_r2"
    params = """commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                const REAL integration_radius1, const REAL integration_radius2, const int gf_index, const REAL *restrict in_gf, REAL localsums[2]"""

    rfm = refmetric.reference_metric[CoordSystem]

    loop_body = ccg.c_codegen(
        [
            rfm.xxSph[0],
            rfm.detgammahat,
        ],
        [
            "const REAL r",
            "const REAL sqrtdetgamma",
        ],
        include_braces=False,
    )

    loop_body += r"""
if(r > integration_radius1 && r < integration_radius2) {
  const REAL gf_of_x = in_gf[IDX4(gf_index, i0, i1, i2)];
  const REAL dV = sqrtdetgamma * dxx0 * dxx1 * dxx2;
  squared_sum += gf_of_x * gf_of_x * dV;
  volume_sum  += dV;
} // END if(r > integration_radius1 && r < integration_radius2)
"""
    body = r"""
  // Unpack grid parameters assuming a single grid
  const int grid = 0;
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"

  // Define reference metric grid
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];

  // Set summation variables to compute l2-norm
  REAL squared_sum = 0.0;
  REAL volume_sum  = 0.0;

"""

    body += lp.simple_loop(
        loop_body="\n" + loop_body,
        read_xxs=True,
        loop_region="interior",
        OMP_custom_pragma=r"#pragma omp parallel for reduction(+:squared_sum,volume_sum)",
    )

    body += r"""
  localsums[0] = squared_sum;
	localsums[1] = volume_sum;
"""

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,  # set_CodeParameters.h is manually included after the declaration of params_struct *restrict params
        body=body,
    )
