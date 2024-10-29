"""
Library of C functions for solving the hyperbolic relaxation equation in curvilinear coordinates, using a reference-metric formalism.
Here we specify overloaded methods and classes utilizing GPU parallelization

Authors: Samuel D. Tootle; sdtootle **at** gmail **dot** com
         Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cf
from pathlib import Path
from types import FrameType as FT
from typing import Any, Dict, Tuple, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.gpu_kernels.kernel_base as gputils
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.diagnostics.output_0d_1d_2d_nearest_gridpoint_slices as out012d
import nrpy.infrastructures.gpu.loop_utilities.cuda.simple_loop as lp
import nrpy.infrastructures.gpu.nrpyelliptic.base_conformally_flat_C_codegen_library as base_npe_classes


# Define functions to set up initial guess
class gpu_register_CFunction_initial_guess_single_point(
    base_npe_classes.base_register_CFunction_initial_guess_single_point
):
    """
    GPU overload to generate a function that computes the initial guess at a single point.
    :param fp_type: floating point precision of sympy expressions translated to C code.
    :return None.
    """

    def __init__(self, fp_type: str = "double") -> None:
        super().__init__(fp_type=fp_type)
        self.cfunc_type = """__device__ __host__ void"""
        self.params = r"""const REAL xx0, const REAL xx1, const REAL xx2,  REAL *restrict uu_ID, REAL *restrict vv_ID
"""

        self.body = ccg.c_codegen(
            [sp.sympify(0), sp.sympify(0)],
            ["*uu_ID", "*vv_ID"],
            verbose=False,
            include_braces=False,
            fp_type=self.fp_type,
        )

        self.register()


def register_CFunction_initial_guess_single_point(
    fp_type: str = "double",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Support function for gpu_register_CFunction_initial_guess_single_point.
    Facilitates generating a function to compute the initial guess at a
    single point and retain compatibility with parallel_codegen.

    :param fp_type: floating point precision of sympy expressions translated to C code.
    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None

    gpu_register_CFunction_initial_guess_single_point(fp_type=fp_type)

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


class gpu_register_CFunction_initial_guess_all_points(
    base_npe_classes.base_register_CFunction_initial_guess_all_points
):
    """
    GPU overload to register the initial guess function for the hyperbolic relaxation equation.

    :param enable_checkpointing: Attempt to read from a checkpoint file before generating initial guess.
    :param fp_type: Floating point type, e.g., "double".

    :return: None.
    """

    def __init__(
        self,
        enable_checkpointing: bool = False,
        fp_type: str = "double",
        **_: Any,
    ) -> None:

        super().__init__(enable_checkpointing, fp_type=fp_type)
        self.loop_body = lp.simple_loop(
            loop_body="initial_guess_single_point(xx0,xx1,xx2,"
            f"&{self.uu_gf_memaccess},"
            f"&{self.vv_gf_memaccess});",
            read_xxs=True,
            loop_region="all points",
            fp_type=self.fp_type,
        ).full_loop_body

        # Put loop_body into a device kernel
        self.device_kernel = gputils.GPU_Kernel(
            self.loop_body,
            {
                "x0": "const REAL *restrict",
                "x1": "const REAL *restrict",
                "x2": "const REAL *restrict",
                "in_gfs": "REAL *restrict",
            },
            "initial_guess_all_points_gpu",
            launch_dict={
                "blocks_per_grid": [],
                "threads_per_block": ["32", "NGHOSTS"],
                "stream": "default",
            },
            fp_type=self.fp_type,
            comments="GPU Kernel to initialize all grid points.",
        )

        self.body = r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
  REAL *restrict x0 = griddata[grid].xx[0];
  REAL *restrict x1 = griddata[grid].xx[1];
  REAL *restrict x2 = griddata[grid].xx[2];
  REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
#include "set_CodeParameters.h"
"""
        self.body += f"{self.device_kernel.launch_block}"
        self.body += f"{self.device_kernel.c_function_call()}"
        self.body += "}\n"

        self.prefunc = self.device_kernel.CFunction.full_function
        self.include_CodeParameters_h = False
        self.register()


def register_CFunction_initial_guess_all_points(
    enable_checkpointing: bool = False,
    fp_type: str = "double",
    **_kwargs: Any,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Support function for gpu_register_CFunction_initial_guess_all_points.
    Facilitates generating a function to compute the initial guess at all
    points and retain compatibility with parallel_codegen.

    :param enable_checkpointing: Attempt to read from a checkpoint file before generating initial guess.
    :param fp_type: Floating point type, e.g., "double".
    :param _kwargs: capture unused arguments from openmp-like calls

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None
    gpu_register_CFunction_initial_guess_all_points(
        enable_checkpointing, fp_type=fp_type
    )

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


# Define functions to set AUXEVOL gridfunctions
class gpu_register_CFunction_auxevol_gfs_single_point(
    base_npe_classes.base_register_CFunction_auxevol_gfs_single_point
):
    """
    Register the C function for the AUXEVOL grid functions at a single point.

    :param CoordSystem: The coordinate system to use in setting up the AUXEVOL gridfunctions.
    :param fp_type: Floating point type, e.g., "double".

    :return: None.
    """

    def __init__(
        self,
        CoordSystem: str,
        fp_type: str = "double",
    ) -> None:
        super().__init__(CoordSystem, fp_type=fp_type)

        self.cfunc_type = """__device__ void"""
        self.params = r"""const size_t streamid, const REAL xx0, const REAL xx1, const REAL xx2, REAL *restrict psi_background, REAL *restrict ADD_times_AUU
"""
        # Is there a better way to do this?
        commondata_refs = [
            "zPunc",
            "P0_x",
            "P0_y",
            "P0_z",
            "P1_x",
            "P1_y",
            "P1_z",
            "S0_x",
            "S0_y",
            "S0_z",
            "S1_x",
            "S1_y",
            "S1_z",
            "bare_mass_0",
            "bare_mass_1",
        ]

        self.body = "// Temporary variables to needed parameters and commondata.\n"
        for p in self.unique_symbols:
            if p not in commondata_refs:
                self.body += f"const REAL {p} = d_params[streamid].{p};\n"
        for cd in commondata_refs:
            self.body += f"const REAL {cd} = d_commondata.{cd};\n"
        self.body += "\n\n"
        self.body += ccg.c_codegen(
            [self.psi_background, self.ADD_times_AUU],
            ["*psi_background", "*ADD_times_AUU"],
            verbose=False,
            include_braces=False,
            fp_type=fp_type,
        )
        self.include_CodeParameters_h = False

        self.register()


def register_CFunction_auxevol_gfs_single_point(
    CoordSystem: str,
    fp_type: str = "double",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Support function for gpu_register_CFunction_auxevol_gfs_single_point.
    Facilitates generating a function for the AUXEVOL grid functions at a single point
    and retain compatibility with parallel_codegen.

    :param CoordSystem: The coordinate system to use in setting up the AUXEVOL gridfunctions.
    :param fp_type: Floating point type, e.g., "double".

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None
    gpu_register_CFunction_auxevol_gfs_single_point(CoordSystem, fp_type=fp_type)

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


class gpu_register_CFunction_auxevol_gfs_all_points(
    base_npe_classes.base_register_CFunction_auxevol_gfs_all_points
):
    """
    GPU overload to generate the C function for the AUXEVOL grid functions at all points.

    :param OMP_collapse: Degree of GPU loop collapsing.
    :param fp_type: Floating point type, e.g., "double".

    :return: None.
    """

    def __init__(
        self,
        fp_type: str = "double",
        **_: Any,
    ) -> None:
        super().__init__(fp_type=fp_type)
        self.loop_body = lp.simple_loop(
            loop_body="auxevol_gfs_single_point(streamid, xx0,xx1,xx2,"
            f"&{self.psi_background_memaccess},"
            f"&{self.ADD_times_AUU_memaccess});",
            read_xxs=True,
            loop_region="all points",
            fp_type=self.fp_type,
        ).full_loop_body

        # Put loop_body into a device kernel
        self.device_kernel = gputils.GPU_Kernel(
            self.loop_body,
            {
                "x0": "const REAL *restrict",
                "x1": "const REAL *restrict",
                "x2": "const REAL *restrict",
                "in_gfs": "REAL *restrict",
            },
            f"{self.name}_gpu",
            launch_dict={
                "blocks_per_grid": [],
                "threads_per_block": ["32", "NGHOSTS"],
                "stream": "default",
            },
            fp_type=self.fp_type,
            comments="GPU Kernel to initialize auxillary grid functions at all grid points.",
        )

        self.body = r"""cpyHosttoDevice_commondata__constant(commondata);
        for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
  REAL *restrict x0 = griddata[grid].xx[0];
  REAL *restrict x1 = griddata[grid].xx[1];
  REAL *restrict x2 = griddata[grid].xx[2];
  REAL *restrict in_gfs = griddata[grid].gridfuncs.auxevol_gfs;
"""
        self.body += f"{self.device_kernel.launch_block}"
        self.body += f"{self.device_kernel.c_function_call()}"
        self.body += "}\n"
        self.prefunc = self.device_kernel.CFunction.full_function
        self.include_CodeParameters_h = False

        self.register()


def register_CFunction_auxevol_gfs_all_points(
    fp_type: str = "double",
    **_kwargs: Any,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Support function for gpu_register_CFunction_auxevol_gfs_single_point.
    Facilitates generating a function for the AUXEVOL grid functions at all points
    and retain compatibility with parallel_codegen.

    :param fp_type: Floating point type, e.g., "double".
    :param _kwargs: capture unused arguments from openmp-like calls

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None
    gpu_register_CFunction_auxevol_gfs_all_points(fp_type=fp_type)

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


class gpu_register_CFunction_variable_wavespeed_gfs_all_points(
    base_npe_classes.base_register_CFunction_variable_wavespeed_gfs_all_points
):
    """
    GPU overload to generate function to compute variable wavespeed based on local grid spacing for a single coordinate system.

    :param CoordSystem: The coordinate system to use in the hyperbolic relaxation.
    :param fp_type: Floating point type, e.g., "double".

    :return: None.
    """

    def __init__(
        self,
        CoordSystem: str,
        fp_type: str = "double",
    ) -> None:
        super().__init__(CoordSystem, fp_type=fp_type)

        self.body = r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
  REAL *restrict x0 = griddata[grid].xx[0];
  REAL *restrict x1 = griddata[grid].xx[1];
  REAL *restrict x2 = griddata[grid].xx[2];
  REAL *restrict in_gfs = griddata[grid].gridfuncs.auxevol_gfs;
"""
        self.loop_body = lp.simple_loop(
            loop_body="\n" + self.dsmin_computation_str,
            read_xxs=True,
            loop_region="interior",
            CoordSystem=self.CoordSystem,
            fp_type=self.fp_type,
        ).full_loop_body
        kernel_body = "// Temporary parameters\n"
        for sym in self.unique_symbols:
            kernel_body += f"const REAL {sym} = d_params[streamid].{sym};\n"
        # Put loop_body into a device kernel
        self.device_kernel = gputils.GPU_Kernel(
            kernel_body + self.loop_body,
            {
                "x0": "const REAL *restrict",
                "x1": "const REAL *restrict",
                "x2": "const REAL *restrict",
                "in_gfs": "REAL *restrict",
                "dt": "const REAL",
                "MINIMUM_GLOBAL_WAVESPEED": "const REAL",
            },
            f"{self.name}_gpu",
            launch_dict={
                "blocks_per_grid": [],
                "threads_per_block": ["32", "NGHOSTS"],
                "stream": "default",
            },
            fp_type=self.fp_type,
            comments="GPU Kernel to initialize auxillary grid functions at all grid points.",
        )
        self.body += f"{self.device_kernel.launch_block}"
        self.body += f"{self.device_kernel.c_function_call()}"
        # We must close the loop that was opened in the line 'for(int grid=0; grid<commondata->NUMGRIDS; grid++) {'
        self.body += r"""} // END LOOP for(int grid=0; grid<commondata->NUMGRIDS; grid++)
                """
        self.prefunc = self.device_kernel.CFunction.full_function
        self.include_CodeParameters_h = False
        self.register()


def register_CFunction_variable_wavespeed_gfs_all_points(
    CoordSystem: str,
    fp_type: str = "double",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Support function for gpu_register_CFunction_variable_wavespeed_gfs_all_points.
    Facilitates generating a function to compute variable wavespeed based on local
    grid spacing for a single coordinate system and retain compatibility with parallel_codegen.

    :param CoordSystem: The coordinate system to use in the hyperbolic relaxation.
    :param fp_type: Floating point type, e.g., "double".

    :return: None if in registration phase, else the updated NRPy environment.
    """
    gpu_register_CFunction_variable_wavespeed_gfs_all_points(
        CoordSystem, fp_type=fp_type
    )

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


class gpu_register_CFunction_initialize_constant_auxevol(
    base_npe_classes.base_register_CFunction_initialize_constant_auxevol
):
    """
    Register function to call all functions that set up AUXEVOL gridfunctions.

    :return: None.
    """

    def __init__(self) -> None:
        super().__init__()
        self.register()


def register_CFunction_initialize_constant_auxevol() -> Union[None, pcg.NRPyEnv_type]:
    """
    Support function for gpu_register_CFunction_initialize_constant_auxevol.
    Facilitates generating a function to call all functions that set up AUXEVOL gridfunctions
    and retain compatibility with parallel_codegen.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None
    gpu_register_CFunction_initialize_constant_auxevol()

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


# Define function to compute the l^2 of a gridfunction
class register_CFunction_compute_L2_norm_of_gridfunction(
    base_npe_classes.base_register_CFunction_compute_L2_norm_of_gridfunction
):
    """
    Register function to compute l2-norm of a gridfunction assuming a single grid.

    Note that parallel codegen is disabled for this function, as it sometimes causes a
    multiprocess race condition on Python 3.6.7

    :param CoordSystem: the rfm coordinate system.
    :param fp_type: Floating point type, e.g., "double".
    :return: None.
    """

    def __init__(
        self,
        CoordSystem: str,
        fp_type: str = "double",
    ) -> None:

        super().__init__(CoordSystem, fp_type=fp_type)
        self.fp_type = "double"
        self.includes += ["BHaH_function_prototypes.h"]
        reduction_loop_body = ccg.c_codegen(
            self.expr_list,
            [
                "const REAL r",
                "const REAL sqrtdetgamma",
            ],
            include_braces=False,
            fp_type=self.fp_type,
        )

        l2_squared_dv_memaccess = "aux_gfs[IDX4(L2_SQUARED_DVGF, i0, i1, i2)]"
        l2_dv_memaccess = "aux_gfs[IDX4(L2_DVGF, i0, i1, i2)]"

        reduction_loop_body += rf"""
if(r < integration_radius) {{
  const REAL gf_of_x = in_gfs[IDX4(gf_index, i0, i1, i2)];
  const REAL dV = sqrtdetgamma * dxx0 * dxx1 * dxx2;
  {l2_squared_dv_memaccess} = gf_of_x * gf_of_x * dV;
  {l2_dv_memaccess} = dV;
}} // END if(r < integration_radius)
"""
        if fp_type == "float":
            reduction_loop_body = reduction_loop_body.replace("REAL", "double")
        self.body += r"""
  params_struct *restrict params = &griddata->params;
#include "set_CodeParameters.h"
  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
  REAL *restrict x0 = griddata->xx[0];
  REAL *restrict x1 = griddata->xx[1];
  REAL *restrict x2 = griddata->xx[2];
  REAL *restrict in_gfs = griddata->gridfuncs.diagnostic_output_gfs;
  REAL *restrict aux_gfs = griddata->gridfuncs.diagnostic_output_gfs2;

  // Since we're performing sums, make sure arrays are zero'd
  cudaMemset(aux_gfs, 0, sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);

  // Set summation variables to compute l2-norm


"""

        self.loop_body = lp.simple_loop(
            loop_body="\n" + reduction_loop_body,
            read_xxs=True,
            loop_region="interior",
            fp_type=self.fp_type,
        ).full_loop_body

        kernel_body = "// Temporary parameters\n"

        # Add symbols that are hard coded into reduction_loop_body
        self.unique_symbols += [f"dxx{i}" for i in range(3)]
        for sym in self.unique_symbols:
            kernel_body += f"const REAL {sym} = d_params[streamid].{sym};\n"

        # Put loop_body into a device kernel
        self.device_kernel = gputils.GPU_Kernel(
            kernel_body + self.loop_body,
            {
                "x0": "const REAL *restrict",
                "x1": "const REAL *restrict",
                "x2": "const REAL *restrict",
                "in_gfs": "const REAL *restrict",
                "aux_gfs": "REAL *restrict",
                "integration_radius": "const REAL",
                "gf_index": "const int",
            },
            f"{self.name}_gpu",
            launch_dict={
                "blocks_per_grid": [],
                "threads_per_block": ["32", "NGHOSTS"],
                "stream": "default",
            },
            fp_type=self.fp_type,
            comments="GPU Kernel to compute L2 quantities pointwise (not summed).",
        )
        self.body += f"{self.device_kernel.launch_block}"
        self.body += f"{self.device_kernel.c_function_call()}"

        self.body += r"""
  REAL squared_sum = find_global__sum(&aux_gfs[IDX4(L2_SQUARED_DVGF, 0, 0, 0)], Nxx_plus_2NGHOSTS_tot);
  REAL volume_sum = find_global__sum(&aux_gfs[IDX4(L2_DVGF, 0, 0, 0)], Nxx_plus_2NGHOSTS_tot);
  // Compute and output the log of the l2-norm.
  REAL local_norm = log10(1e-16 + sqrt(squared_sum / volume_sum));  // 1e-16 + ... avoids log10(0)
  return local_norm;
"""

        cfc.register_CFunction(
            prefunc=self.device_kernel.CFunction.full_function,
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            CoordSystem_for_wrapper_func="",
            name=self.name,
            params=self.params,
            include_CodeParameters_h=False,  # set_CodeParameters.h is manually included after the declaration of params_struct *restrict params
            body=self.body,
        )


# Define diagnostics function
class gpu_register_CFunction_diagnostics(
    base_npe_classes.base_register_CFunction_diagnostics
):
    """
    GPU overload to facilitate generating the C function for simulation diagnostics.

    :param CoordSystem: Coordinate system used.
    :param default_diagnostics_out_every: Specifies the default diagnostics output frequency.
    :param enable_progress_indicator: Whether to enable the progress indicator.
    :param axis_filename_tuple: Tuple containing filename and variables for axis output.
    :param plane_filename_tuple: Tuple containing filename and variables for plane output.
    :param out_quantities_dict: Dictionary or string specifying output quantities.
    :return: None.
    """

    def __init__(
        self,
        CoordSystem: str,
        default_diagnostics_out_every: int,
        enable_progress_indicator: bool = False,
        axis_filename_tuple: Tuple[str, str] = (
            "out1d-AXIS-n-%08d.txt",
            "nn",
        ),
        plane_filename_tuple: Tuple[str, str] = (
            "out2d-PLANE-n-%08d.txt",
            "nn",
        ),
        out_quantities_dict: Union[str, Dict[Tuple[str, str], str]] = "default",
    ) -> None:
        super().__init__(
            default_diagnostics_out_every,
            out_quantities_dict=out_quantities_dict,
        )
        self.params += ", griddata_struct *restrict griddata_host"

        # This has to be here to avoid type issues with mypy
        # An error will throw in super().__init__() if out_quantities_dict != dict
        self.out_quantities_dict: Dict[Tuple[str, str], str] = self.out_quantities_dict

        for axis in ["y", "z"]:
            out012d.register_CFunction_diagnostics_nearest_1d_axis(
                CoordSystem=CoordSystem,
                out_quantities_dict=self.out_quantities_dict,
                filename_tuple=axis_filename_tuple,
                axis=axis,
                pointer_decorator="[[maybe_unused]] ",
            )
        for plane in ["xy", "yz"]:
            out012d.register_CFunction_diagnostics_nearest_2d_plane(
                CoordSystem=CoordSystem,
                out_quantities_dict=self.out_quantities_dict,
                filename_tuple=plane_filename_tuple,
                plane=plane,
                pointer_decorator="[[maybe_unused]] ",
            )

        self.body = r"""  // Output progress to stderr
  progress_indicator(commondata, griddata);

  // Grid data output
  const int n_step = commondata->nn, outevery = commondata->diagnostics_output_every;

  REAL global_norm = 0.0;
  for(int grid = 0; grid < commondata->NUMGRIDS; ++grid) {

  // Set gridfunctions aliases
  REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
  REAL *restrict diagnostic_output_gfs = griddata[grid].gridfuncs.diagnostic_output_gfs;

  // Set params and rfm_struct
  params_struct *restrict params = &griddata[grid].params;
  const rfm_struct *restrict rfmstruct = &griddata[grid].rfmstruct;
#include "set_CodeParameters.h"
  REAL *restrict host_y_n_gfs = griddata_host[grid].gridfuncs.y_n_gfs;
  REAL *restrict host_diag_gfs = griddata_host[grid].gridfuncs.diagnostic_output_gfs;
  if (n_step % outevery == 0) {
    size_t streamid = cpyDevicetoHost__gf(commondata, params, host_y_n_gfs, y_n_gfs, UUGF, UUGF);
  }

  // Compute Hamiltonian constraint violation and store it at diagnostic_output_gfs
  compute_residual_all_points(commondata, params, rfmstruct, auxevol_gfs, y_n_gfs, diagnostic_output_gfs);
  cudaDeviceSynchronize();
  if (n_step % outevery == 0) {
    size_t streamid = cpyDevicetoHost__gf(commondata, params, host_diag_gfs, diagnostic_output_gfs, RESIDUAL_HGF, RESIDUAL_HGF);
  }

  // Set integration radius for l2-norm computation
  const REAL integration_radius = 1000;

  // Compute l2-norm of Hamiltonian constraint violation
  const REAL residual_H = compute_L2_norm_of_gridfunction(commondata, &griddata[grid], integration_radius, RESIDUAL_HGF, diagnostic_output_gfs);
  global_norm = MAX(global_norm, residual_H);
  }
  // Update residual to be used in stop condition
  commondata->log10_current_residual = global_norm;

  // Output l2-norm of Hamiltonian constraint violation to file
  {
    // Only consider a single grid for now.
    const int grid = 0;
    params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
    char filename[256];
    sprintf(filename, "residual_l2_norm.txt");
    FILE *outfile = (nn == 0) ? fopen(filename, "w") : fopen(filename, "a");
    if (!outfile) {
      fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
      exit(1);
    }
    fprintf(outfile, "%6d %10.4e %.17e\n", nn, time, global_norm);
    fclose(outfile);
  }


  if (n_step % outevery == 0) {
    // Only consider a single grid for now.
    const int grid = 0;
    params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
    // Set reference metric grid xx
    REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
        xx[ww] = griddata_host[grid].xx[ww];

    // Ensure all device workers are done
    cudaDeviceSynchronize();

    // 1D output
    diagnostics_nearest_1d_y_axis(commondata, params, xx, &griddata_host[grid].gridfuncs);
    diagnostics_nearest_1d_z_axis(commondata, params, xx, &griddata_host[grid].gridfuncs);

    // 2D output
    diagnostics_nearest_2d_xy_plane(commondata, params, xx, &griddata_host[grid].gridfuncs);
    diagnostics_nearest_2d_yz_plane(commondata, params, xx, &griddata_host[grid].gridfuncs);
  }
"""
        if enable_progress_indicator:
            self.body += "progress_indicator(commondata, griddata);"
        self.body += r"""
  if (commondata->time + commondata->dt > commondata->t_final)
    printf("\n");
"""
        cfc.register_CFunction(
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=False,
            body=self.body,
        )


def register_CFunction_diagnostics(
    CoordSystem: str,
    default_diagnostics_out_every: int,
    enable_progress_indicator: bool = False,
    axis_filename_tuple: Tuple[str, str] = (
        "out1d-AXIS-n-%08d.txt",
        "nn",
    ),
    plane_filename_tuple: Tuple[str, str] = (
        "out2d-PLANE-n-%08d.txt",
        "nn",
    ),
    out_quantities_dict: Union[str, Dict[Tuple[str, str], str]] = "default",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Support function for gpu_register_CFunction_diagnostics.
    Facilitates generating a function to compute simulation diagnostics
    and retain compatibility with parallel_codegen.

    :param CoordSystem: Coordinate system used.
    :param default_diagnostics_out_every: Specifies the default diagnostics output frequency.
    :param enable_progress_indicator: Whether to enable the progress indicator.
    :param axis_filename_tuple: Tuple containing filename and variables for axis output.
    :param plane_filename_tuple: Tuple containing filename and variables for plane output.
    :param out_quantities_dict: Dictionary or string specifying output quantities.
    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None
    gpu_register_CFunction_diagnostics(
        CoordSystem,
        default_diagnostics_out_every,
        enable_progress_indicator=enable_progress_indicator,
        axis_filename_tuple=axis_filename_tuple,
        plane_filename_tuple=plane_filename_tuple,
        out_quantities_dict=out_quantities_dict,
    )

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


# Define function to evaluate stop conditions
def register_CFunction_check_stop_conditions() -> Union[None, pcg.NRPyEnv_type]:
    """
    Support function for base_register_CFunction_check_stop_conditions.
    Facilitates generating a function to evaluate stop conditions
    and retain compatibility with parallel_codegen.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None
    base_npe_classes.base_register_CFunction_check_stop_conditions()

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


# Define function to evaluate RHSs
class gpu_register_CFunction_rhs_eval(
    base_npe_classes.base_register_CFunction_rhs_eval
):
    """
    GPU overload to facilitate generating the right-hand side (RHS) evaluation function.

    This function sets the right-hand side of the hyperbolic relaxation equation according to the
    selected coordinate system and specified parameters.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_simd: Whether to enable SIMD.
    :param OMP_collapse: Level of GPU loop collapsing.
    :param fp_type: Floating point type, e.g., "double".
    :param enable_intrinsics: Toggle using CUDA intrinsics for calculations.

    :return: None.
    """

    def __init__(
        self,
        CoordSystem: str,
        enable_rfm_precompute: bool,
        fp_type: str = "double",
        enable_intrinsics: bool = False,
    ) -> None:

        super().__init__(CoordSystem, enable_rfm_precompute)

        if enable_intrinsics:
            self.includes += [str(Path("simd") / "simd_intrinsics.h")]

        self.simple_loop = lp.simple_loop(
            loop_body=ccg.c_codegen(
                [self.rhs.uu_rhs, self.rhs.vv_rhs],
                [
                    gri.BHaHGridFunction.access_gf("uu", gf_array_name="rhs_gfs"),
                    gri.BHaHGridFunction.access_gf("vv", gf_array_name="rhs_gfs"),
                ],
                enable_fd_codegen=True,
                fp_type=fp_type,
                rational_const_alias="static constexpr",
                enable_simd=enable_intrinsics,
            ),
            loop_region="interior",
            CoordSystem=CoordSystem,
            enable_rfm_precompute=enable_rfm_precompute,
            read_xxs=not enable_rfm_precompute,
            fp_type=fp_type,
            enable_intrinsics=enable_intrinsics,
        )
        self.loop_body = self.simple_loop.full_loop_body.replace(
            "const REAL f", "[[maybe_unused]] const REAL f"
        )
        self.loop_body = self.loop_body.replace(
            "const REAL_SIMD_ARRAY f", "[[maybe_unused]] const REAL_SIMD_ARRAY f"
        )
        self.loop_body = self.loop_body.replace(
            "const double dbl", "static constexpr double dbl"
        )
        self.kernel_comments = "GPU Kernel to evaluate RHS on the interior."
        self.params_dict_coord = {f"x{i}": "const REAL *restrict" for i in range(3)}
        self.body = ""

        if enable_rfm_precompute:

            self.params_dict_coord = {
                f'{rfm_f.replace("REAL *restrict ","")[:-1]}': "const REAL *restrict"
                for rfm_f in self.simple_loop.rfmp.BHaH_defines_list
            }
            params_dict = {f"rfm_{k}": v for k, v in self.params_dict_coord.items()}
            params_dict["auxevol_gfs"] = "const REAL *restrict"
            params_dict["in_gfs"] = "const REAL *restrict"
            params_dict["rhs_gfs"] = "REAL *restrict"
            params_dict["eta_damping"] = "const REAL"

            # Put loop_body into a device kernel
            self.device_kernel = gputils.GPU_Kernel(
                self.loop_body,
                params_dict,
                "rhs_eval_gpu",
                launch_dict={
                    "blocks_per_grid": [],
                    "threads_per_block": ["32", "NGHOSTS"],
                    "stream": "default",
                },
                fp_type=fp_type,
                comments=self.kernel_comments,
            )
            for k, v in self.params_dict_coord.items():
                self.body += f"{v} rfm_{k} = rfmstruct->{k};\n"
        else:
            raise ValueError(
                "rhs_eval without rfm_precompute has not been implemented."
            )

        self.body += f"{self.device_kernel.launch_block}"
        self.body += f"{self.device_kernel.c_function_call()}"
        cfc.register_CFunction(
            prefunc=self.device_kernel.CFunction.full_function,
            include_CodeParameters_h=True,
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            CoordSystem_for_wrapper_func="",
            name=self.name,
            params=self.params,
            body=self.body,
        )


def register_CFunction_rhs_eval(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    fp_type: str = "double",
    enable_intrinsics: bool = False,
    **_kwargs: Any,
) -> Union[None, pcg.NRPyEnv_type]:
    """

    Support function for gpu_register_CFunction_rhs_eval.
    Facilitates generating the right-hand side (RHS) evaluation function
    for the hyperbolic relaxation equation and retain compatibility with parallel_codegen.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param fp_type: Floating point type, e.g., "double".
    :param enable_intrinsics: Toogle using CUDA intrinsics for calculations.
    :param _kwargs: capture unused arguments from openmp-like calls

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None
    gpu_register_CFunction_rhs_eval(
        CoordSystem,
        enable_rfm_precompute,
        fp_type=fp_type,
        enable_intrinsics=enable_intrinsics,
    )

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


# Define function to compute residual the solution
class gpu_register_CFunction_compute_residual_all_points(
    base_npe_classes.base_register_CFunction_compute_residual_all_points
):
    """
    GPU overload to facilitate generating the residual evaluation function.

    This function sets the residual of the Hamiltonian constraint in the hyperbolic
    relaxation equation according to the selected coordinate system and specified
    parameters.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_intrinsicsv: Whether to enable cuda intrinsics.
    :param OMP_collapse: Level of GPU loop collapsing.
    :param fp_type: Floating point type, e.g., "double".

    :return: None.
    """

    def __init__(
        self,
        CoordSystem: str,
        enable_rfm_precompute: bool,
        enable_intrinsics: bool,
        fp_type: str = "double",
    ) -> None:
        super().__init__(
            CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
        )
        self.body = ""
        if enable_intrinsics:
            self.includes += [str(Path("simd") / "simd_intrinsics.h")]

        self.simple_loop = lp.simple_loop(
            loop_body=ccg.c_codegen(
                [self.rhs.residual],
                [
                    gri.BHaHGridFunction.access_gf(
                        "residual_H", gf_array_name="aux_gfs"
                    ),
                ],
                enable_fd_codegen=True,
                enable_simd=enable_intrinsics,
                fp_type=fp_type,
            ),
            loop_region="interior",
            enable_intrinsics=enable_intrinsics,
            CoordSystem=CoordSystem,
            enable_rfm_precompute=enable_rfm_precompute,
            read_xxs=not enable_rfm_precompute,
            fp_type=fp_type,
        )
        self.kernel_body = self.simple_loop.full_loop_body
        self.params_dict_coord = {
            f'{rfm_f.replace("REAL *restrict ","")[:-1]}': "const REAL *restrict"
            for rfm_f in self.simple_loop.rfmp.BHaH_defines_list
        }
        params_dict = {f"rfm_{k}": v for k, v in self.params_dict_coord.items()}
        params_dict["auxevol_gfs"] = "const REAL *restrict"
        params_dict["in_gfs"] = "const REAL *restrict"
        params_dict["aux_gfs"] = "REAL *restrict"
        self.kernel_body = self.kernel_body.replace(
            "const REAL f", "[[maybe_unused]] const REAL f"
        )
        self.kernel_body = self.kernel_body.replace(
            "const REAL_SIMD_ARRAY f", "[[maybe_unused]] const REAL_SIMD_ARRAY f"
        )
        self.kernel_body = self.kernel_body.replace(
            "const double dbl", "static constexpr double dbl"
        )
        self.device_kernel = gputils.GPU_Kernel(
            self.kernel_body,
            params_dict,
            f"{self.name}_gpu",
            launch_dict={
                "blocks_per_grid": [],
                "threads_per_block": ["32", "NGHOSTS"],
                "stream": "default",
            },
            fp_type=fp_type,
            comments="GPU Kernel to compute the residual throughout the grid.",
        )
        self.prefunc = self.device_kernel.CFunction.full_function
        for k, v in self.params_dict_coord.items():
            self.body += f"{v} rfm_{k} = rfmstruct->{k};\n"
        self.body += self.device_kernel.launch_block + "\n\n"
        self.body += self.device_kernel.c_function_call()
        cfc.register_CFunction(
            prefunc=self.prefunc,
            include_CodeParameters_h=True,
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            CoordSystem_for_wrapper_func="",
            name=self.name,
            params=self.params,
            body=self.body,
        )


def register_CFunction_compute_residual_all_points(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    fp_type: str = "double",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Support function for gpu_register_CFunction_compute_residual_all_points.
    Facilitates generating the function to evaluate the residual at all points
    and retain compatibility with parallel_codegen.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_simd: Whether to enable SIMD.
    :param fp_type: Floating point type, e.g., "double".

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None
    gpu_register_CFunction_compute_residual_all_points(
        CoordSystem, enable_rfm_precompute, enable_simd, fp_type=fp_type
    )

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
