"""
C function for setting the initial data to be used in hyperbolic relaxation.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par
from nrpy.infrastructures import BHaH


# Define functions to set up initial guess
def register_CFunction_initial_guess_single_point() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for initial guess of solution at a single point.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h"]

    desc = r"""Compute initial guess at a single point."""
    cfunc_type = "void"
    name = "initial_guess_single_point"
    cfunc_decorators = "__device__ __host__" if parallelization in ["cuda"] else ""
    params = "const REAL xx0, const REAL xx1, const REAL xx2,  REAL *restrict uu_ID, REAL *restrict vv_ID"
    body = ccg.c_codegen(
        [sp.sympify(0), sp.sympify(0)],
        ["*uu_ID", "*vv_ID"],
        verbose=False,
        include_braces=False,
    )
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        cfunc_decorators=cfunc_decorators,
    )
    return pcg.NRPyEnv()


def register_CFunction_initial_guess_all_points(
    OMP_collapse: int,
    enable_checkpointing: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the initial guess function for the hyperbolic relaxation equation.

    :param enable_checkpointing: Attempt to read from a checkpoint file before generating initial guess.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""Set initial guess to solutions of hyperbolic relaxation equation at all points."""
    cfunc_type = "void"
    name = "initial_data"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata_host, griddata_struct *restrict griddata"
        if parallelization in ["cuda"]
        else "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    uu_gf_memaccess = gri.BHaHGridFunction.access_gf("uu")
    vv_gf_memaccess = gri.BHaHGridFunction.access_gf("vv")
    body = ""
    if enable_checkpointing:
        body += """// Attempt to read checkpoint file. If it doesn't exist, then continue. Otherwise return.
if( read_checkpoint(commondata, griddata) ) return;
""".replace(
            "griddata",
            "griddata_host, griddata" if parallelization in ["cuda"] else "griddata",
        )
    body += r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
  REAL *restrict x0 = griddata[grid].xx[0];
  REAL *restrict x1 = griddata[grid].xx[1];
  REAL *restrict x2 = griddata[grid].xx[2];
  REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
"""
    kernel_body = f"{parallel_utils.get_loop_parameters(parallelization)}\n"
    kernel_body += BHaH.simple_loop.simple_loop(
        loop_body="initial_guess_single_point(xx0,xx1,xx2,"
        f"&{uu_gf_memaccess},"
        f"&{vv_gf_memaccess});",
        read_xxs=True,
        loop_region="all points",
        OMP_collapse=OMP_collapse,
    )
    for i in range(3):
        kernel_body = kernel_body.replace(f"xx[{i}]", f"x{i}")

    comments = "Kernel to initialize all grid points."
    # Prepare the argument dicts
    arg_dict_cuda = {
        "x0": "const REAL *restrict",
        "x1": "const REAL *restrict",
        "x2": "const REAL *restrict",
        "in_gfs": "REAL *restrict",
    }
    arg_dict_host = {
        # "commondata": "const commondata_struct *restrict",
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
        thread_tiling_macro_suffix="NELL_ID",
    )

    body += f"{new_body}\n"
    body += "}\n"
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
    return pcg.NRPyEnv()
