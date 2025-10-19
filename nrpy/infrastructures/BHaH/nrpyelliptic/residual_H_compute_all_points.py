"""
C function for hyperbolic relaxation diagnostics.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.equations.nrpyelliptic.ConformallyFlat_RHSs import (
    HyperbolicRelaxationCurvilinearRHSs,
)
from nrpy.infrastructures import BHaH


# Define function to compute residual the solution
def register_CFunction_residual_H_compute_all_points(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_intrinsics: bool,
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the residual evaluation function.

    This function sets the residual of the Hamiltonian constraint in the hyperbolic
    relaxation equation according to the selected coordinate system and specified
    parameters.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_intrinsics: Whether to enable hardware intrinsics (e.g. simd).
    :param OMP_collapse: Level of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.

    Doctest:
    >>> import nrpy.c_function as cfc
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import unittest_CoordSystems
    >>> name="residual_H_compute_all_points"
    ... for CoordSystem in unittest_CoordSystems:
    ...   cfc.CFunction_dict.clear()
    ...   _ = register_CFunction_residual_H_compute_all_points(CoordSystem, True, True)
    ...   generated_str = cfc.CFunction_dict[f'{name}__rfm__{CoordSystem}'].full_function
    ...   validation_desc = f"{name}__{CoordSystem}"
    ...   validate_strings(generated_str, validation_desc)
    Setting up reference_metric[SinhSymTP_rfm_precompute]...
    Setting up reference_metric[HoleySinhSpherical_rfm_precompute]...
    Setting up reference_metric[Cartesian_rfm_precompute]...
    Setting up reference_metric[SinhCylindricalv2n2_rfm_precompute]...
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    enable_rfm_precompute = False
    includes = ["BHaH_defines.h"]
    desc = r"""Compute residual of the Hamiltonian constraint for the hyperbolic relaxation equation."""
    cfunc_type = "void"
    name = "residual_H_compute_all_points"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
                REAL *restrict xx[3], const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs,
                REAL *restrict dest_gf_address"""
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )
    # Populate residual_H
    rhs = HyperbolicRelaxationCurvilinearRHSs(CoordSystem, enable_rfm_precompute)
    orig_par = par.parval_from_str("parallelization")
    par.set_parval_from_str("parallelization", "openmp")
    body = BHaH.simple_loop.simple_loop(
        loop_body=ccg.c_codegen(
            [rhs.residual],
            ["dest_gf_address[IDX3(i0,i1,i2)]"],
            enable_fd_codegen=True,
            enable_simd=False,
        ),
        loop_region="interior",
        enable_intrinsics=False,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )
    par.set_parval_from_str("parallelization", orig_par)

    cfc.register_CFunction(
        subdirectory="diagnostics",
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=body,
        enable_simd=False,
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
