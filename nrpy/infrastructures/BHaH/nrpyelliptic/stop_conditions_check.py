"""
C function for hyperbolic relaxation diagnostics.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par


# Define function to evaluate stop conditions
def register_CFunction_stop_conditions_check() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function to evaluate stop conditions.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h"]
    desc = "Evaluate stop conditions."
    cfunc_type = "void"
    name = "stop_conditions_check"
    params = (
        """commondata_struct *restrict commondata, griddata_struct *restrict griddata"""
    )

    # Register a parameter to stop hyperbolic relaxation
    _stop_relaxation = par.register_CodeParameter(
        "bool", __name__, "stop_relaxation", False, commondata=True
    )

    # Register parameter that sets the total number of relaxation steps
    _nn_max = par.register_CodeParameter("int", __name__, "nn_max", 0, commondata=True)

    # Register parameter that sets the tolerance for log of residual
    _log10_residual_tolerance = par.register_CodeParameter(
        "REAL", __name__, "log10_residual_tolerance", -15.8, commondata=True
    )

    # Register parameter that sets log of residual to be updated at every time step
    _log10_current_residual = par.register_CodeParameter(
        "REAL", __name__, "log10_current_residual", 1.0, commondata=True
    )

    body = r"""  // Since this version of NRPyElliptic is unigrid, we simply set the grid index to 0
  const int grid = 0;

  // Set params
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"

  // Check if total number of iteration steps has been reached
  if ((nn >= nn_max) || (log10_current_residual < log10_residual_tolerance)){
    printf("\nExiting main loop after %8d iterations\n", nn);
    printf("The tolerance for the logarithmic residual is %.8e\n", log10_residual_tolerance);
    printf("Exiting relaxation with logarithmic residual of %.8e\n", log10_current_residual);
    commondata->stop_relaxation = true;
  }
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
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
