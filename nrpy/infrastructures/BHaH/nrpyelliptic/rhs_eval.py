"""
C function for evaluating the RHS of the hyperbolic relaxation equation in curvilinear coordinates, using a reference-metric formalism.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from pathlib import Path
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.simple_loop as lp
from nrpy.equations.nrpyelliptic.ConformallyFlat_RHSs import (
    HyperbolicRelaxationCurvilinearRHSs,
)


# Define function to evaluate RHSs
def register_CFunction_rhs_eval(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the right-hand side (RHS) evaluation function for the hyperbolic relaxation equation.

    This function sets the right-hand side of the hyperbolic relaxation equation according to the
    selected coordinate system and specified parameters.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_simd: Whether to enable SIMD.
    :param OMP_collapse: Level of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    if enable_simd:
        includes += [str(Path("intrinsics") / "simd_intrinsics.h")]
    desc = r"""Set RHSs for hyperbolic relaxation equation."""
    cfunc_type = "void"
    name = "rhs_eval"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict rhs_gfs"
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )
    # Populate uu_rhs, vv_rhs
    rhs = HyperbolicRelaxationCurvilinearRHSs(CoordSystem, enable_rfm_precompute)
    body = lp.simple_loop(
        loop_body=ccg.c_codegen(
            [rhs.uu_rhs, rhs.vv_rhs],
            [
                gri.BHaHGridFunction.access_gf("uu", gf_array_name="rhs_gfs"),
                gri.BHaHGridFunction.access_gf("vv", gf_array_name="rhs_gfs"),
            ],
            enable_fd_codegen=True,
            enable_simd=enable_simd,
        ),
        loop_region="interior",
        enable_intrinsics=enable_simd,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )

    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        body=body,
        enable_simd=enable_simd,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
