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
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.params as par
from nrpy.equations.nrpyelliptic.ConformallyFlat_RHSs import (
    HyperbolicRelaxationCurvilinearRHSs,
)
from nrpy.helpers.expression_utils import get_params_commondata_symbols_from_expr_list


# Define function to evaluate RHSs
def register_CFunction_rhs_eval(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_intrinsics: bool,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the right-hand side (RHS) evaluation function for the hyperbolic relaxation equation.

    This function sets the right-hand side of the hyperbolic relaxation equation according to the
    selected coordinate system and specified parameters.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_intrinsics: Whether to enable hardware intrinsics.
    :param OMP_collapse: Level of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
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
    desc = r"""Set RHSs for hyperbolic relaxation equation."""
    cfunc_type = "void"
    name = "rhs_eval"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict rhs_gfs"
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )
    body = ""
    if not enable_rfm_precompute:
        for i in range(3):
            body += f"const REAL *restrict x{i} = xx[{i}];\n"
    # Populate uu_rhs, vv_rhs
    rhs = HyperbolicRelaxationCurvilinearRHSs(CoordSystem, enable_rfm_precompute)
    expr_list = [rhs.uu_rhs, rhs.vv_rhs]
    loop_body = lp.simple_loop(
        loop_body=ccg.c_codegen(
            expr_list,
            [
                gri.BHaHGridFunction.access_gf("uu", gf_array_name="rhs_gfs"),
                gri.BHaHGridFunction.access_gf("vv", gf_array_name="rhs_gfs"),
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
    # Since simd intrinsics only supports double precision, CCG has constants hardcoded as dbl,
    # this is not the case for CUDA intrinsics, for example, so we can use REAL
    # to explore other precisions, e.g. single precision.
    # Note: dbl prefix is not changed to avoid potential replace issues in the calculations themselves
    loop_body = loop_body.replace(
        "static const double dbl",
        (
            "static constexpr REAL dbl"
            if parallelization in ["cuda"]
            else "static const double dbl"
        ),
    )
    loop_params = parallel_utils.get_loop_parameters(
        parallelization, enable_intrinsics=enable_intrinsics
    )
    params_symbols, _ = get_params_commondata_symbols_from_expr_list(
        [rhs.residual], exclude=[f"xx{i}" for i in range(3)]
    )
    loop_params += "// Load necessary parameters from params_struct\n"
    for param in params_symbols:
        loop_params += f"const REAL {param} = {parallel_utils.get_params_access(parallelization)}{param};\n"

    loop_params += "\n// Setup parameters from function arguments\n"
    loop_params += (
        "const REAL NOSIMDeta_dampening = eta_damping_in;\n"
        "MAYBE_UNUSED const REAL_SIMD_ARRAY eta_damping = ConstSIMD(NOSIMDeta_dampening);\n"
        if enable_intrinsics
        else "const REAL eta_damping = eta_damping_in;\n"
    )
    loop_params += "\n"

    if parallelization == "cuda":
        loop_params = loop_params.replace("SIMD", "CUDA")

    comments = "Kernel to evaluate RHS on the interior."

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
            "rhs_gfs": "REAL *restrict",
            "eta_damping_in": "const REAL",
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
            "threads_per_block": ["64", "NGHOSTS"],
            "stream": "default",
        },
        thread_tiling_macro_suffix="NELL_RHS",
    )

    new_body = new_body.replace("eta_damping_in", "eta_damping")
    for i in range(3):
        new_body = new_body.replace(f"x{i}", f"xx[{i}]")
    body = f"{new_body}\n"

    cfc.register_CFunction(
        prefunc=prefunc,
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
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
