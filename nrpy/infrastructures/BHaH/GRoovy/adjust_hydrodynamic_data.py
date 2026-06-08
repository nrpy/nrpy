"""
C function to reapply GRHayL primitive-variable limits after spacetime updates.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.simple_loop as lp


def register_CFunction_adjust_hydrodynamic_data(
    CoordSystem: str,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register a routine that rebuilds and re-limits hydrodynamic primitives after initial data import.

    This routine is intended for use after the spacetime data have been
    defined, so that the GRHayL primitive variables and stored `u^0` remain
    consistent with the updated metric.

    :param CoordSystem: The coordinate system.
    :param OMP_collapse: Level of OpenMP loop collapsing.
    :return: None if in registration phase, else the updated NRPy environment.
    """
    # Step 1: Register the function with NRPy's parallel codegen infrastructure.
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 2: Set up C function headers, signature, and parameters.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Adjust hydrodynamic initial data after spacetime import"
    cfunc_type = "void"
    name = "adjust_hydrodynamic_data"
    params = (
        "const commondata_struct *restrict commondata, "
        "const params_struct *restrict params, "
        "REAL *restrict xx[3], "
        "REAL *restrict evol_gfs, "
        "REAL *restrict auxevol_gfs"
    )

    # Step 3: Define the pointwise hydrodynamic adjustment.
    loop_body = r"""
      ghl_primitive_quantities prims;
      ghl_metric_quantities metric;
      ghl_conservative_quantities cons;

      // Rebuild GRHayL structs from the current gridfunction data.
      basis_transform_rfm_basis_to_Cartesian(
          commondata, params, &prims, &cons, &metric, i0, i1, i2, xx,
          auxevol_gfs, evol_gfs);

      bool speed_limited;
      ghl_enforce_primitive_limits_and_compute_u0(
          &commondata->ghl_params, &commondata->eos, &metric, &prims,
          &speed_limited);

      // Store the limited variables back in the native gridfunction basis.
      basis_transform_Cartesian_to_rfm_basis(
          commondata, params, &prims, &cons, i0, i1, i2, xx,
          auxevol_gfs, evol_gfs);

      auxevol_gfs[IDX4(U4UTGF, i0, i1, i2)] = prims.u0;
"""

    # Step 4: Wrap the pointwise adjustment in a grid loop.
    body = lp.simple_loop(
        loop_body=loop_body,
        loop_region="all points",
        enable_intrinsics=False,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=False,
        read_xxs=False,
        OMP_collapse=OMP_collapse,
    )

    # Step 5: Register the final C function.
    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        CoordSystem_for_wrapper_func=CoordSystem,
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
