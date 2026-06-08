"""
C function to recompute entropy from the evolved GRHD primitive state.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.simple_loop as lp


def register_CFunction_compute_entropy(
    OMP_collapse: int,
    evolving_temperature: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register a routine that recomputes entropy from other GRHD variables.

    :param OMP_collapse: Level of OpenMP loop collapsing.
    :param evolving_temperature: Whether temperature is an evolved primitive.
    :return: None if in registration phase, else the updated NRPy environment.
    """
    # Step 1: Register the function with NRPy's parallel codegen infrastructure.
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 2: Set up C function headers, signature, and parameters.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Compute entropy using the current GRHD primitive data"
    cfunc_type = "void"
    name = "compute_entropy"
    params = (
        "const commondata_struct *restrict commondata, "
        "const params_struct *restrict params, "
        "const ghl_eos_parameters *restrict eos, "
        "REAL *restrict auxevol_gfs"
    )

    # Step 3: Define the pointwise entropy update.
    if evolving_temperature:
        loop_body = r"""
    REAL dummy;

    // Enforce EOS bounds before recomputing pressure and entropy from T.
    ghl_tabulated_enforce_bounds_rho_Ye_T(
        eos,
        &auxevol_gfs[IDX4(RHOBGF, i0, i1, i2)],
        &auxevol_gfs[IDX4(YEGF, i0, i1, i2)],
        &auxevol_gfs[IDX4(TEMPERATUREGF, i0, i1, i2)]);

    ghl_tabulated_compute_P_eps_S_from_T(
        eos,
        auxevol_gfs[IDX4(RHOBGF, i0, i1, i2)],
        auxevol_gfs[IDX4(YEGF, i0, i1, i2)],
        auxevol_gfs[IDX4(TEMPERATUREGF, i0, i1, i2)],
        &auxevol_gfs[IDX4(PGF, i0, i1, i2)],
        &dummy,
        &auxevol_gfs[IDX4(SGF, i0, i1, i2)]);
"""
    else:
        loop_body = r"""
    // For a hybrid EOS, entropy depends only on density and pressure.
    auxevol_gfs[IDX4(SGF, i0, i1, i2)] = ghl_hybrid_compute_entropy_function(
        eos,
        auxevol_gfs[IDX4(RHOBGF, i0, i1, i2)],
        auxevol_gfs[IDX4(PGF, i0, i1, i2)]);
"""

    # Step 4: Wrap the pointwise update in a grid loop.
    body = lp.simple_loop(
        loop_body=loop_body,
        loop_region="all points",
        enable_intrinsics=False,
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
