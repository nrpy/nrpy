"""
C function for setting constant auxiliary terms for hyperbolic relaxation.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Tuple, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms import (
    compute_psi_background_and_ADD_times_AUU,
)
from nrpy.helpers.expression_utils import get_params_commondata_symbols_from_expr_list


def generate_prefunc_variable_wavespeed_gfs_all_points(
    CoordSystem: str,
    enable_intrinsics: bool = False,
) -> Tuple[str, str]:
    """
    Generate function to compute variable wavespeed based on local grid spacing for a single coordinate system.

    :param CoordSystem: The coordinate system to use in the hyperbolic relaxation.
    :param enable_intrinsics: Whether to enable hardware intrinsics (default: False).

    :return: None if in registration phase, else the updated NRPy environment.
    """
    comments = r"""Kernel to compute variable wavespeed for all grids based on local grid spacing."""
    name = "variable_wavespeed_gfs_all_points"
    parallelization = par.parval_from_str("parallelization")
    enable_intrinsics = enable_intrinsics if parallelization == "cuda" else False
    # Prepare the argument dicts
    arg_dict_cuda = {
        "x0": "const REAL *restrict",
        "x1": "const REAL *restrict",
        "x2": "const REAL *restrict",
        "in_gfs": "REAL *restrict",
        "dt": "const REAL",
        "MINIMUM_GLOBAL_WAVESPEED": "const REAL",
    }
    arg_dict_host = {
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }

    rfm = refmetric.reference_metric[CoordSystem]
    dxx0, dxx1, dxx2 = sp.symbols("dxx0 dxx1 dxx2", real=True)
    expr_list = [
        rfm.scalefactor_orthog[0] * dxx0,
        rfm.scalefactor_orthog[1] * dxx1,
        rfm.scalefactor_orthog[2] * dxx2,
    ]
    dsmin_computation_str = ccg.c_codegen(
        expr_list,
        ["const REAL dsmin0", "const REAL dsmin1", "const REAL dsmin2"],
        include_braces=False,
        enable_simd=enable_intrinsics,
    ).replace("SIMD", "CUDA" if parallelization == "cuda" else "SIMD")

    variable_wavespeed_memaccess = gri.BHaHGridFunction.access_gf("variable_wavespeed")

    dsmin_computation_str += f"""\n// Set local wavespeed
        {variable_wavespeed_memaccess} = MINIMUM_GLOBAL_WAVESPEED * MIN(dsmin0, MIN(dsmin1, dsmin2)) / dt;\n"""

    loop_body = lp.simple_loop(
        loop_body="\n" + dsmin_computation_str,
        read_xxs=True,
        loop_region="interior",
        CoordSystem=CoordSystem,
        enable_intrinsics=enable_intrinsics,
    )
    for i in range(3):
        loop_body = loop_body.replace(f"xx[{i}]", f"x{i}")

    loop_params = parallel_utils.get_loop_parameters(
        parallelization,
        enable_intrinsics=enable_intrinsics,
    )

    params_symbols, _ = get_params_commondata_symbols_from_expr_list(
        expr_list, exclude=[f"xx{i}" for i in range(3)]
    )
    loop_params += "// Load necessary parameters from params_struct\n"
    for param in params_symbols:
        loop_params += f"const REAL {param} = {parallel_utils.get_params_access(parallelization)}{param};\n"
    loop_params += "\n"

    kernel_body = f"{loop_params}\n{loop_body}"
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
        thread_tiling_macro_suffix="NELL_WAVESPEED",
    )

    return prefunc, new_body


# Define functions to set AUXEVOL gridfunctions
def generate_prefunc_auxevol_gfs_single_point(
    CoordSystem: str,
    enable_intrinsics: bool = False,
) -> Tuple[str, str]:
    """
    Generate function for the AUXEVOL grid functions at a single point.

    :param CoordSystem: The coordinate system to use in setting up the AUXEVOL gridfunctions.
    :param enable_intrinsics: Whether to enable hardware intrinsics (default: False).

    :return: None if in registration phase, else the updated NRPy environment.
    """
    # Compute psi_background and ADD_times_AUU
    psi_background, ADD_times_AUU = compute_psi_background_and_ADD_times_AUU(
        CoordSystem
    )

    parallelization = par.parval_from_str("parallelization")

    comments = r"""Kernel to compute AUXEVOL grid functions at a single point."""
    name = "auxevol_gfs_single_point"
    # Prepare the argument dicts
    arg_dict_cuda = {
        "xx0": "const REAL",
        "xx1": "const REAL",
        "xx2": "const REAL",
        "psi_background": "REAL *restrict",
        "ADD_times_AUU": "REAL *restrict",
    }
    arg_dict_host = {
        "commondata": "const commondata_struct *restrict",
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }

    # Generate kernel body
    kernel_body = ""
    expr_list = [psi_background, ADD_times_AUU]
    params_symbols, commondata_symbols = get_params_commondata_symbols_from_expr_list(
        expr_list, exclude=[f"xx{i}" for i in range(3)]
    )
    kernel_body += "// Load necessary parameters from params_struct\n"
    for param in params_symbols:
        kernel_body += f"const REAL {param} = {parallel_utils.get_params_access(parallelization)}{param};\n"
    kernel_body += "\n// Load necessary parameters from commondata_struct\n"
    for param in commondata_symbols:
        kernel_body += f"const REAL {param} = {parallel_utils.get_commondata_access(parallelization)}{param};\n"
    kernel_body += "\n"
    kernel_body += ccg.c_codegen(
        expr_list,
        ["*psi_background", "*ADD_times_AUU"],
        verbose=False,
        include_braces=False,
        enable_simd=enable_intrinsics,
    ).replace("SIMD", "CUDA" if parallelization == "cuda" else "SIMD")

    prefunc, new_body = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body,
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=comments,
        launch_dict=None,
        cfunc_decorators="__device__",
    )
    return prefunc, new_body


def generate_prefunc_auxevol_gfs_all_points(
    CoordSystem: str,
    OMP_collapse: int = 1,
    enable_intrinsics: bool = False,
) -> Tuple[str, str]:
    """
    Register the C function for the AUXEVOL grid functions at all points.

    :param OMP_collapse: Degree of OpenMP loop collapsing.
    :param CoordSystem: The coordinate system to use in setting up the AUXEVOL gridfunctions.
    :param enable_intrinsics: Whether to enable hardware intrinsics (default: False).

    :return: None if in registration phase, else the updated NRPy environment.
    """
    parallelization = par.parval_from_str("parallelization")

    name = "auxevol_gfs_all_points"
    psi_background_memaccess = gri.BHaHGridFunction.access_gf("psi_background")
    ADD_times_AUU_memaccess = gri.BHaHGridFunction.access_gf("ADD_times_AUU")
    enable_intrinsics = enable_intrinsics if parallelization == "cuda" else False

    kernel_body = f"{parallel_utils.get_loop_parameters(parallelization, enable_intrinsics=enable_intrinsics)}\n"

    single_point_prefunc, single_point_launch = (
        generate_prefunc_auxevol_gfs_single_point(CoordSystem)
    )
    single_point_launch = single_point_launch.replace(
        "psi_background", f"&{psi_background_memaccess}"
    ).replace("ADD_times_AUU", f"&{ADD_times_AUU_memaccess}")
    kernel_body += lp.simple_loop(
        f"{single_point_launch}",
        read_xxs=True,
        loop_region="all points",
        OMP_collapse=OMP_collapse,
        enable_intrinsics=enable_intrinsics,
    )
    for i in range(3):
        kernel_body = kernel_body.replace(f"xx[{i}]", f"x{i}")

    comments = "Kernel to initialize auxillary grid functions at all grid points."
    # Prepare the argument dicts
    arg_dict_cuda = {
        "x0": "const REAL *restrict",
        "x1": "const REAL *restrict",
        "x2": "const REAL *restrict",
        "in_gfs": "REAL *restrict",
    }
    arg_dict_host = {
        "commondata": "const commondata_struct *restrict",
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }
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
        thread_tiling_macro_suffix="NELL_AUX",
    )
    prefunc = f"{single_point_prefunc}\n{prefunc}"
    return prefunc, new_body


def register_CFunction_auxevol_gfs_set_to_constant(
    CoordSystem: str,
    OMP_collapse: int = 1,
    enable_intrinsics: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function to call all functions that set up AUXEVOL gridfunctions.

    :param CoordSystem: The coordinate system to use in setting up the AUXEVOL gridfunctions.
    :param OMP_collapse: Degree of OpenMP loop collapsing.
    :param enable_intrinsics: Whether to enable hardware intrinsics (default: False).
    :return: None if in registration phase, else the updated NRPy environment.

    Doctest:
    >>> import nrpy.c_function as cfc
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import unittest_CoordSystems
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name="auxevol_gfs_set_to_constant"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in unittest_CoordSystems:
    ...       cfc.CFunction_dict.clear()
    ...       _ = register_CFunction_auxevol_gfs_set_to_constant(CoordSystem)
    ...       generated_str = cfc.CFunction_dict[f'{name}__rfm__{CoordSystem}'].full_function
    ...       validation_desc = f"{name}__{parallelization}__{CoordSystem}".replace(" ", "_")
    ...       validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    Setting up reference_metric[SinhSymTP]...
    Setting up reference_metric[HoleySinhSpherical]...
    Setting up reference_metric[Cartesian]...
    Setting up reference_metric[SinhCylindricalv2n2]...
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""Call functions that set up all AUXEVOL gridfunctions."""
    cfunc_type = "void"
    name = "auxevol_gfs_set_to_constant"
    params = "commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs"
    body = (
        "cpyHosttoDevice_commondata__constant(commondata);\n"
        if parallelization in ["cuda"]
        else ""
    )

    auxevol_gfs_all_points_prefunc, auxevol_gfs_all_points_launch = (
        generate_prefunc_auxevol_gfs_all_points(
            CoordSystem=CoordSystem,
            OMP_collapse=OMP_collapse,
            enable_intrinsics=enable_intrinsics,
        )
    )
    wavespeed_all_points_prefunc, wavespeed_all_points_launch = (
        generate_prefunc_variable_wavespeed_gfs_all_points(
            CoordSystem=CoordSystem, enable_intrinsics=enable_intrinsics
        )
    )
    auxevol_gfs_all_points_launch = (
        f"{{{auxevol_gfs_all_points_launch}}}"
        if parallelization in ["cuda"]
        else auxevol_gfs_all_points_launch
    )
    wavespeed_all_points_launch = (
        f"{{{wavespeed_all_points_launch}}}"
        if parallelization in ["cuda"]
        else wavespeed_all_points_launch
    )
    body += rf"""
    REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
    REAL *restrict x0 = xx[0];
    REAL *restrict x1 = xx[1];
    REAL *restrict x2 = xx[2];

    // Set up variable wavespeed
    {wavespeed_all_points_launch.replace("in_gfs", "auxevol_gfs")}

    // Set up all other AUXEVOL gridfunctions
    {auxevol_gfs_all_points_launch.replace("in_gfs", "auxevol_gfs")}
    """

    prefunc = rf"""
    {auxevol_gfs_all_points_prefunc}
    {wavespeed_all_points_prefunc}
    """

    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
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
