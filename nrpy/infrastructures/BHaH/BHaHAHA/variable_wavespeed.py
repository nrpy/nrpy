"""
Set up variable wavespeed for the computational grid, largely written by Thiago Assumpção.

Handles the configuration of a variable wavespeed for use in numerical relativity simulations,
including registering the appropriate functions and code generation tasks.

Authors: Zachariah B. Etienne; zachetie **at** gmail **dot* com
         Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par


def register_CFunction_variable_wavespeed(
    CoordSystem: str = "Spherical",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register a C function to set the variable wavespeed on the computational grid.

    Registers a function for calculating the variable wavespeed. The function is also registered
    with the parallel codegen system if the registration phase is active.

    :param CoordSystem: Coordinate system in use (e.g., "Spherical").
    :return: None if in registration phase, otherwise an instance of NRPyEnv_type.

    >>> register_CFunction_variable_wavespeed("Cartesian", "float")
    """
    # Suppress unused argument warnings by documenting that they are placeholders for future extensions.
    _ = CoordSystem

    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    _MINIMUM_GLOBAL_WAVESPEED = par.register_CodeParameter(
        "REAL",
        __name__,
        name="MINIMUM_GLOBAL_WAVESPEED",
        defaultvalue=0.2,
        commondata=True,
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Set variable wavespeed"
    cfunc_type = "void"
    name = "variable_wavespeed"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    body = r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  const params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];
  REAL *restrict evol_gfs = griddata[grid].gridfuncs.y_n_gfs;
  REAL *restrict in_gfs = griddata[grid].gridfuncs.auxevol_gfs;
#pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
      const REAL xx2 = xx[2][i2];
      for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
        const REAL xx1 = xx[1][i1];
        for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
          const REAL hh = evol_gfs[IDX4(HHGF, i0, i1, i2)];
          const REAL dsmin1 = dxx1 * hh;
          const REAL dsmin2 = dxx2 * hh * sin(xx1);

          // Set local wavespeed
          in_gfs[IDX4(ZZVARIABLE_WAVESPEEDGF, i0, i1, i2)] = MINIMUM_GLOBAL_WAVESPEED * MIN(dsmin1, dsmin2) / dt;

        } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
      } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
    } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
  } // END LOOP for(int grid=0; grid<commondata->NUMGRIDS; grid++)
"""

    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
