"""
C function for hyperbolic relaxation diagnostics.

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
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.params as par
from nrpy.equations.nrpyelliptic.ConformallyFlat_RHSs import (
    HyperbolicRelaxationCurvilinearRHSs,
)
from nrpy.helpers.expression_utils import get_params_commondata_symbols_from_expr_list


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
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name="residual_H_compute_all_points"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in unittest_CoordSystems:
    ...       cfc.CFunction_dict.clear()
    ...       _ = register_CFunction_residual_H_compute_all_points(CoordSystem, True, True)
    ...       generated_str = cfc.CFunction_dict[f'{name}__rfm__{CoordSystem}'].full_function
    ...       validation_desc = f"{name}__{parallelization}__{CoordSystem}".replace(" ", "_")
    ...       validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    Setting up reference_metric[SinhSymTP_rfm_precompute]...
    Setting up reference_metric[HoleySinhSpherical_rfm_precompute]...
    Setting up reference_metric[Cartesian_rfm_precompute]...
    Setting up reference_metric[SinhCylindricalv2n2_rfm_precompute]...
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h"]
    if enable_intrinsics:
        includes += (
            [str(Path("intrinsics") / "cuda_intrinsics.h")]
            if parallelization in ["cuda"]
            else [str(Path("intrinsics") / "simd_intrinsics.h")]
        )
    desc = r"""Compute residual of the Hamiltonian constraint for the hyperbolic relaxation equation."""
    cfunc_type = "void"
    name = "residual_H_compute_all_points"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
                REAL *restrict xx[3], const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs,
                REAL *restrict aux_gfs"""
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )
    # Populate residual_H
    rhs = HyperbolicRelaxationCurvilinearRHSs(CoordSystem, enable_rfm_precompute)
    loop_body = lp.simple_loop(
        loop_body=ccg.c_codegen(
            [rhs.residual],
            [
                gri.BHaHGridFunction.access_gf("residual_H", gf_array_name="aux_gfs"),
            ],
            enable_fd_codegen=True,
            enable_simd=enable_intrinsics,
        ).replace("SIMD", "CUDA" if parallelization == "cuda" else "SIMD"),
        loop_region="interior",
        enable_intrinsics=enable_intrinsics,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )

    loop_params = parallel_utils.get_loop_parameters(
        parallelization, enable_intrinsics=enable_intrinsics
    )

    params_symbols, _ = get_params_commondata_symbols_from_expr_list(
        [rhs.residual], exclude=[f"xx{i}" for i in range(3)]
    )
    if len(params_symbols) > 0:
        loop_params += "// Load necessary parameters from params_struct\n"
        for param in params_symbols:
            loop_params += f"const REAL {param} = {parallel_utils.get_params_access(parallelization)}{param};\n"
        loop_params += "\n"

    comments = "Kernel to compute the residual throughout the grid."

    # Prepare the argument dicts
    arg_dict_cuda = (
        {"rfmstruct": "const rfm_struct *restrict"}
        if enable_rfm_precompute
        else {f"x{i}": "const REAL *restrict" for i in range(3)}
    )
    arg_dict_cuda.update(
        {
            "auxevol_gfs": "const REAL *restrict",
            "in_gfs": "const REAL *restrict",
            "aux_gfs": "REAL *restrict",
        }
    )
    arg_dict_host = {
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }

    kernel_body = f"{loop_params}\n{loop_body}"

    for i in range(3):
        kernel_body = kernel_body.replace(f"xx[{i}]", f"x{i}")
    prefunc, new_body = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body,
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=comments,
        launch_dict={
            "blocks_per_grid": [],
            "threads_per_block": ["32", "NGHOSTS"],
            "stream": "default",
        },
        thread_tiling_macro_suffix="NELL_H",
    )
    for i in range(3):
        new_body = new_body.replace(f"x{i}", f"xx[{i}]")
    body = f"{new_body}\n"

    cfc.register_CFunction(
        subdirectory="diagnostics",
        prefunc=prefunc,
        include_CodeParameters_h=False,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=body,
        enable_simd=enable_intrinsics,
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
