"""
C functions for setting the initial data of the wave equation to an exact solution at a particular time.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Tuple, Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.params as par
from nrpy.equations.wave_equation.WaveEquation_Solutions_InitialData import (
    WaveEquation_solution_Cartesian,
)
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)


def generate_CFunction_exact_solution_single_Cartesian_point(
    WaveType: str = "SphericalGaussian",
    default_sigma: float = 3.0,
    default_k0: float = 1.0,
    default_k1: float = 1.0,
    default_k2: float = 1.0,
) -> Tuple[str, str]:
    """
    Register the C function for the exact solution at a single point.

    :param WaveType: The type of wave: SphericalGaussian or PlaneWave
    :param default_sigma: The default value for the Gaussian width (sigma).
    :param default_k0: The default value for the plane wave wavenumber k in the x-direction.
    :param default_k1: The default value for the plane wave wavenumber k in the y-direction.
    :param default_k2: The default value for the plane wave wavenumber k in the z-direction.

    :return: A tuple containing the kernel and launch body for the exact solution.
    """
    # Populate uu_ID, vv_ID
    exactsoln = WaveEquation_solution_Cartesian(
        WaveType=WaveType,
        default_sigma=default_sigma,
        default_k0=default_k0,
        default_k1=default_k1,
        default_k2=default_k2,
    )

    parallelization = par.parval_from_str("parallelization")

    cfunc_decorators = "__device__" if parallelization in ["cuda"] else ""
    desc = r"""Exact solution at a single Cartesian point (x, y, z) = (xCart0, xCart1, xCart2)."""
    cfunc_type = "static void"
    name = "exact_solution_single_Cartesian_point"

    arg_dict_cuda = {
        "params": "const params_struct *restrict",
        "xCart0": "const REAL",
        "xCart1": "const REAL",
        "xCart2": "const REAL",
        "exact_soln_UUGF": "REAL *restrict",
        "exact_soln_VVGF": "REAL *restrict",
    }

    arg_dict_host = {
        "commondata": "const commondata_struct *restrict",
        **arg_dict_cuda,
    }

    expr_list = [exactsoln.uu_exactsoln, exactsoln.vv_exactsoln]

    # Find symbols stored in params
    param_symbols, commondata_symbols = get_params_commondata_symbols_from_expr_list(
        expr_list, exclude=[f"xx{j}" for j in range(3)]
    )

    # Since xx_to_cart requires params as a host & device kernel we pass
    # params struct as a pointer rather than use, e.g. constant CUDA
    # memory
    params_definitions = generate_definition_header(
        param_symbols,
        var_access=parallel_utils.get_params_access("openmp"),
    )
    commondata_definitions = generate_definition_header(
        commondata_symbols,
        var_access=parallel_utils.get_commondata_access(parallelization),
    )

    kernel_body = ccg.c_codegen(
        expr_list,
        ["*exact_soln_UUGF", "*exact_soln_VVGF"],
        verbose=False,
        include_braces=False,
    )

    kernel_body = f"{params_definitions}\n{commondata_definitions}\n\n{kernel_body}"

    kernel, launch_body = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body.replace("SIMD", "CUDA" if parallelization == "cuda" else "SIMD"),
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=desc,
        cfunc_type=cfunc_type,
        launchblock_with_braces=False,
        cfunc_decorators=cfunc_decorators,
    )
    return kernel, launch_body


def generate_CFunction_initial_data_compute(
    OMP_collapse: int,
    WaveType: str = "SphericalGaussian",
    default_sigma: float = 3.0,
    default_k0: float = 1.0,
    default_k1: float = 1.0,
    default_k2: float = 1.0,
) -> Tuple[str, str]:
    """
    Generate the initial data compute kernel for the wave equation with specific parameters.

    :param OMP_collapse: Degree of OpenMP loop collapsing.
    :param WaveType: The type of wave: SphericalGaussian or PlaneWave
    :param default_sigma: The default value for the Gaussian width (sigma).
    :param default_k0: The default value for the plane wave wavenumber k in the x-direction.
    :param default_k1: The default value for the plane wave wavenumber k in the y-direction.
    :param default_k2: The default value for the plane wave wavenumber k in the z-direction.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    parallelization = par.parval_from_str("parallelization")

    desc = r"""Set initial data to params.time==0 corresponds to the initial data."""
    cfunc_type = "void"
    name = "initial_data"
    uu_gf_memaccess = gri.BHaHGridFunction.access_gf("uu")
    vv_gf_memaccess = gri.BHaHGridFunction.access_gf("vv")
    kernel_body = ""

    kernel_body += parallel_utils.get_loop_parameters(parallelization)

    arg_dict_cuda = {
        "params": "const params_struct *restrict",
        "x0": "const REAL *restrict",
        "x1": "const REAL *restrict",
        "x2": "const REAL *restrict",
        "in_gfs": "REAL *restrict",
    }

    arg_dict_host = {
        "commondata": "const commondata_struct *restrict",
        **arg_dict_cuda,
    }

    exact_singlept_kernel, exact_singlept_launch = (
        generate_CFunction_exact_solution_single_Cartesian_point(
            WaveType=WaveType,
            default_sigma=default_sigma,
            default_k0=default_k0,
            default_k1=default_k1,
            default_k2=default_k2,
        )
    )

    loop_body = (
        "REAL xCart[3]; REAL xOrig[3] = {xx0, xx1, xx2}; xx_to_Cart(params, xOrig, xCart);\n"
        f"{exact_singlept_launch}"
    )
    for i in range(3):
        loop_body = loop_body.replace(f"xCart{i}", f"xCart[{i}]")
    loop_body = loop_body.replace("exact_soln_UUGF", f"&{uu_gf_memaccess}")
    loop_body = loop_body.replace("exact_soln_VVGF", f"&{vv_gf_memaccess}")

    kernel_body += lp.simple_loop(
        loop_body=loop_body,
        read_xxs=True,
        loop_region="all points",
        OMP_collapse=OMP_collapse,
    )

    kernel, launch_body = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body.replace("SIMD", "CUDA" if parallelization == "cuda" else "SIMD"),
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=desc,
        cfunc_type=cfunc_type,
        launchblock_with_braces=False,
        thread_tiling_macro_suffix="WAVE_ID_EXACT",
    )

    for i in range(3):
        kernel = kernel.replace(f"xx[{i}]", f"x{i}")

    prefunc = exact_singlept_kernel + kernel

    return prefunc, launch_body


def register_CFunction_initial_data_exact(
    OMP_collapse: int,
    WaveType: str = "SphericalGaussian",
    default_sigma: float = 3.0,
    default_k0: float = 1.0,
    default_k1: float = 1.0,
    default_k2: float = 1.0,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the initial data function for the wave equation with specific parameters.

    :param OMP_collapse: Degree of OpenMP loop collapsing.
    :param WaveType: The type of wave: SphericalGaussian or PlaneWave
    :param default_sigma: The default value for the Gaussian width (sigma).
    :param default_k0: The default value for the plane wave wavenumber k in the x-direction.
    :param default_k1: The default value for the plane wave wavenumber k in the y-direction.
    :param default_k2: The default value for the plane wave wavenumber k in the z-direction.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""Set initial data to params.time==0 corresponds to the initial data."""
    cfunc_type = "void"
    name = "initial_data_exact"
    arg_dict_host = {
        "commondata": "const commondata_struct *restrict",
        "params_in": "const params_struct *restrict",
        "x0": "const REAL *restrict",
        "x1": "const REAL *restrict",
        "x2": "const REAL *restrict",
        "in_gfs": "REAL *restrict",
    }
    params = ", ".join([f"{ptype} {pname}" for pname, ptype in arg_dict_host.items()])
    body = ""

    if parallelization in ["cuda"]:
        body += "cpyHosttoDevice_params__constant(params_in, params_in->grid_idx % NUM_STREAMS);\n"
        body += "params_struct * params;\n"
        body += "BHAH_MALLOC_DEVICE(params, sizeof(params_struct));\n"
        body += (
            "BHAH_MEMCPY_HOST_TO_DEVICE(params, params_in, sizeof(params_struct));\n"
        )
    else:
        body += "const params_struct *restrict params = params_in;\n"

    prefunc, id_compute_launch = generate_CFunction_initial_data_compute(
        OMP_collapse=OMP_collapse,
        WaveType=WaveType,
        default_sigma=default_sigma,
        default_k0=default_k0,
        default_k1=default_k1,
        default_k2=default_k2,
    )

    body += id_compute_launch.replace("params->", "params_in->")
    if parallelization in ["cuda"]:
        body += "BHAH_FREE_DEVICE(params);\n"
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_initial_data(
    enable_checkpointing: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the initial data function for the wave equation with specific parameters.

    :param enable_checkpointing: Attempt to read from a checkpoint file before generating initial data.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""Set initial data to params.time==0 corresponds to the initial data."""
    cfunc_type = "void"
    name = "initial_data"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata_host, griddata_struct *restrict griddata"
        if parallelization in ["cuda"]
        else "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    body = ""
    host_griddata = "griddata_host" if parallelization in ["cuda"] else "griddata"
    if enable_checkpointing:
        body += """// Attempt to read checkpoint file. If it doesn't exist, then continue. Otherwise return.
if( read_checkpoint(commondata, griddata) ) return;
""".replace(
            "griddata",
            f"{host_griddata}, griddata" if parallelization in ["cuda"] else "griddata",
        )
    if parallelization in ["cuda"]:
        body += "cpyHosttoDevice_commondata__constant(commondata);\n"

    body += r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
"""
    for i in range(3):
        body += f"REAL *restrict x{i} = griddata[grid].xx[{i}];\n"
    body += "REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;\n"
    body += "params_struct *restrict params = &griddata[grid].params;\n"
    body += "initial_data_exact(commondata, params, x0, x1, x2, in_gfs);\n"
    body += "}\n"
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
