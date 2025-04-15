"""
C function for hyperbolic relaxation diagnostics.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from pathlib import Path
from types import FrameType as FT
from typing import Dict, Tuple, Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.diagnostics.output_0d_1d_2d_nearest_gridpoint_slices as out012d
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.equations.nrpyelliptic.ConformallyFlat_RHSs import (
    HyperbolicRelaxationCurvilinearRHSs,
)
from nrpy.helpers.expression_utils import get_params_commondata_symbols_from_expr_list


# Define function to compute the l^2 of a gridfunction
def register_CFunction_compute_L2_norm_of_gridfunction(
    CoordSystem: str,
) -> None:
    """
    Register function to compute l2-norm of a gridfunction assuming a single grid.

    Note that parallel codegen is disabled for this function, as it sometimes causes a
    multiprocess race condition on Python 3.6.7

    :param CoordSystem: the rfm coordinate system.

    Doctest:
    >>> import nrpy.c_function as cfc
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import unittest_CoordSystems
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name="compute_L2_norm_of_gridfunction"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in unittest_CoordSystems:
    ...       cfc.CFunction_dict.clear()
    ...       _ = register_CFunction_compute_L2_norm_of_gridfunction(CoordSystem)
    ...       generated_str = cfc.CFunction_dict[f'{name}__rfm__{CoordSystem}'].full_function
    ...       validation_desc = f"{name}__{parallelization}__{CoordSystem}".replace(" ", "_")
    ...       validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    Setting up reference_metric[SinhSymTP]...
    Setting up reference_metric[HoleySinhSpherical]...
    Setting up reference_metric[Cartesian]...
    Setting up reference_metric[SinhCylindricalv2n2]...
    """
    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Compute l2-norm of a gridfunction assuming a single grid."
    cfunc_type = "void"
    name = "compute_L2_norm_of_gridfunction"
    params = """commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3],
                const REAL integration_radius, const int gf_index, REAL * l2norm, const REAL *restrict in_gfs"""
    params += ", REAL *restrict aux_gfs" if parallelization in ["cuda"] else ""

    rfm = refmetric.reference_metric[CoordSystem]

    fp_type = par.parval_from_str("fp_type")

    # Enforce at least double precision since calcuations
    # can go beyond single precision numerical limits for some
    # coordinate systems
    fp_type_alias = "DOUBLE" if fp_type == "float" else "REAL"
    ccg_fp_type = "double" if fp_type == "float" else fp_type
    expr_list = [
        rfm.xxSph[0],
        rfm.detgammahat,
    ]

    # Define the norm calculations within the loop
    reduction_loop_body = ccg.c_codegen(
        expr_list,
        [
            "const DOUBLE r",
            "const DOUBLE sqrtdetgamma",
        ],
        include_braces=False,
        fp_type=ccg_fp_type,
        fp_type_alias=fp_type_alias,
    )

    reduction_loop_body += r"""
if(r < integration_radius) {
  const DOUBLE gf_of_x = in_gfs[IDX4(gf_index, i0, i1, i2)];
  const DOUBLE dV = sqrtdetgamma * dxx0 * dxx1 * dxx2;
"""
    reduction_loop_body += (
        r"""
  aux_gfs[IDX4(L2_SQUARED_DVGF, i0, i1, i2)] = gf_of_x * gf_of_x * dV;
  aux_gfs[IDX4(L2_DVGF, i0, i1, i2)] = dV;
} // END if(r < integration_radius)
"""
        if parallelization in ["cuda"]
        else r"""
  squared_sum += gf_of_x * gf_of_x * dV;
  volume_sum  += dV;
} // END if(r < integration_radius)
"""
    )
    reduction_loop_body = reduction_loop_body.replace("REAL", fp_type_alias)

    OMP_custom_pragma = (
        r"#pragma omp parallel for reduction(+:squared_sum,volume_sum)"
        if parallelization not in ["cuda"]
        else ""
    )

    # Generate the loop for the reduction_loop_body
    loop_body = lp.simple_loop(
        loop_body="\n" + reduction_loop_body,
        read_xxs=True,
        loop_region="interior",
        OMP_custom_pragma=OMP_custom_pragma,
    )

    # Device code computes the local L2 quantities which are reduced
    # in a separate algorithm, find_global__sum.
    # Host code uses OpenMP reduction #pragma to compute the
    # L2 quantities and the global norms in a single kernel
    comments = (
        "Kernel to compute L2 quantities pointwise (not summed)."
        if parallelization in ["cuda"]
        else "Kernel to compute L2 quantities pointwise (summed)."
    )
    # Prepare the argument dictionaries
    arg_dict_cuda = {
        "x0": "const REAL *restrict",
        "x1": "const REAL *restrict",
        "x2": "const REAL *restrict",
        "in_gfs": "const REAL *restrict",
        "aux_gfs": "REAL *restrict",
        "integration_radius": "const REAL",
        "gf_index": "const int",
    }
    arg_dict_host = {
        "params": "const params_struct *restrict",
        "x0": "const REAL *restrict",
        "x1": "const REAL *restrict",
        "x2": "const REAL *restrict",
        "in_gfs": "const REAL *restrict",
        "squared_sum_final": "REAL *restrict",
        "volume_sum_final": "REAL *restrict",
        "integration_radius": "const REAL",
        "gf_index": "const int",
    }
    loop_params = parallel_utils.get_loop_parameters(parallelization)
    params_symbols, _ = get_params_commondata_symbols_from_expr_list(
        expr_list, exclude=[f"xx{i}" for i in range(3)]
    )
    loop_params += "// Load necessary parameters from params_struct\n"
    # We have to manually add dxx{i} here since they are not in the SymPy
    # expression list, but, instead, are manually used in calculations
    # in reduction_loop_body above
    for param in params_symbols + [f"dxx{i}" for i in range(3)]:
        loop_params += f"const REAL {param} = {parallel_utils.get_params_access(parallelization)}{param};\n"
    loop_params += (
        "\n"
        if parallelization in ["cuda"]
        else r"""
  DOUBLE squared_sum = 0.0;
  DOUBLE volume_sum  = 0.0;
        """
    )

    kernel_body = f"{loop_params}\n{loop_body}"
    for i in range(3):
        kernel_body = kernel_body.replace(f"xx[{i}]", f"x{i}")

    if parallelization not in ["cuda"]:
        kernel_body += r"""*squared_sum_final = squared_sum;
*volume_sum_final = volume_sum;
        """

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
        thread_tiling_macro_suffix="NELL_GRIDL2",
    )

    # Define launch kernel body
    body = r"""
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
  REAL *restrict x0 = xx[0];
  REAL *restrict x1 = xx[1];
  REAL *restrict x2 = xx[2];
"""

    body += (
        r"""
  // Since we're performing sums, make sure arrays are zero'd
  cudaMemset(aux_gfs, 0, sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
"""
        if parallelization in ["cuda"]
        else r"""
  // Set summation variables to compute l2-norm
  DOUBLE squared_sum = 0.0;
  DOUBLE volume_sum  = 0.0;
"""
    )

    body += f"{new_body}\n"
    body += (
        r"""
  // Set summation variables to compute l2-norm
  REAL squared_sum = find_global__sum(&aux_gfs[IDX4(L2_SQUARED_DVGF, 0, 0, 0)], Nxx_plus_2NGHOSTS_tot);
  REAL volume_sum = find_global__sum(&aux_gfs[IDX4(L2_DVGF, 0, 0, 0)], Nxx_plus_2NGHOSTS_tot);
  // Compute and output the log of the l2-norm.
  REAL local_norm = log10(1e-16 + sqrt(squared_sum / volume_sum));  // 1e-16 + ... avoids log10(0)
"""
        if parallelization in ["cuda"]
        else ""
    )

    # For host reduction code.
    body = body.replace("squared_sum_final", "&squared_sum").replace(
        "volume_sum_final", "&volume_sum"
    )

    body += r"""
  // Compute and output the log of the l2-norm.
  *l2norm = log10(1e-16 + sqrt(squared_sum / volume_sum));  // 1e-16 + ... avoids log10(0)
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


# Define function to compute residual the solution
def register_CFunction_compute_residual_all_points(
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
    >>> name="compute_residual_all_points"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in unittest_CoordSystems:
    ...       cfc.CFunction_dict.clear()
    ...       _ = register_CFunction_compute_residual_all_points(CoordSystem, True, True)
    ...       generated_str = cfc.CFunction_dict[f'{name}__rfm__{CoordSystem}'].full_function
    ...       validation_desc = f"{name}__{parallelization}__{CoordSystem}".replace(" ", "_")
    ...       validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
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
    name = "compute_residual_all_points"
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
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


# Define diagnostics function
def register_CFunction_diagnostics(
    CoordSystem: str,
    enable_rfm_precompute: bool,
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
    Register C function for simulation diagnostics.

    :param CoordSystem: Coordinate system used.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param default_diagnostics_out_every: Specifies the default diagnostics output frequency.
    :param enable_progress_indicator: Whether to enable the progress indicator.
    :param axis_filename_tuple: Tuple containing filename and variables for axis output.
    :param plane_filename_tuple: Tuple containing filename and variables for plane output.
    :param out_quantities_dict: Dictionary or string specifying output quantities.
    :raises TypeError: If `out_quantities_dict` is not a dictionary and not set to "default".
    :return: None if in registration phase, else the updated NRPy environment.

    Doctest:
    >>> import nrpy.c_function as cfc
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.params as par
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name="diagnostics"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    cfc.CFunction_dict.clear()
    ...    _ = register_CFunction_diagnostics("Cartesian", True, 100)
    ...    generated_str = cfc.CFunction_dict[f'{name}'].full_function
    ...    validation_desc = f"{name}__{parallelization}".replace(" ", "_")
    ...    validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    _ = par.CodeParameter(
        "int",
        __name__,
        "diagnostics_output_every",
        default_diagnostics_out_every,
        commondata=True,
    )

    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Diagnostics."
    cfunc_type = "void"
    name = "diagnostics"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    params += (
        ", griddata_struct *restrict griddata_host"
        if parallelization in ["cuda"]
        else ""
    )

    # fmt: off
    if out_quantities_dict == "default":
        out_quantities_dict = {
            ("REAL", "numUU"): "y_n_gfs[IDX4pt(UUGF, idx3)]",
            ("REAL", "log10ResidualH"): "log10(fabs(diagnostic_output_gfs[IDX4pt(RESIDUAL_HGF, idx3)] + 1e-16))",
        }
    if not isinstance(out_quantities_dict, dict):
        raise TypeError(f"out_quantities_dict was initialized to {out_quantities_dict}, which is not a dictionary!")
    # fmt: on
    for axis in ["y", "z"]:
        out012d.register_CFunction_diagnostics_nearest_1d_axis(
            CoordSystem=CoordSystem,
            out_quantities_dict=out_quantities_dict,
            axis=axis,
            filename_tuple=axis_filename_tuple,
        )
    for plane in ["xy", "yz"]:
        out012d.register_CFunction_diagnostics_nearest_2d_plane(
            CoordSystem=CoordSystem,
            out_quantities_dict=out_quantities_dict,
            plane=plane,
            filename_tuple=plane_filename_tuple,
        )

    body = r"""  // Output progress to stderr
  progress_indicator(commondata, griddata);

  // Grid data output
  const int n_step = commondata->nn, outevery = commondata->diagnostics_output_every;
  const REAL time = commondata->time;

  REAL global_norm = -1e9;
  for(int grid = 0; grid < commondata->NUMGRIDS; ++grid) {

  // Set gridfunctions aliases
  REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
  REAL *restrict diagnostic_output_gfs = griddata[grid].gridfuncs.diagnostic_output_gfs;
  MAYBE_UNUSED REAL *restrict diagnostic_output_gfs2 = griddata[grid].gridfuncs.diagnostic_output_gfs2;
  // Set params
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
"""
    if parallelization in ["cuda"]:
        body += r"""
  REAL *restrict host_y_n_gfs = griddata_host[grid].gridfuncs.y_n_gfs;
  REAL *restrict host_diag_gfs = griddata_host[grid].gridfuncs.diagnostic_output_gfs;
  if (n_step % outevery == 0) {
    size_t streamid = params->grid_idx % NUM_STREAMS;
    cpyDevicetoHost__gf(commondata, params, host_y_n_gfs, y_n_gfs, UUGF, UUGF, streamid);
  }
"""
    if enable_rfm_precompute:
        body += r"""  // Set rfm_struct
  const rfm_struct *restrict rfmstruct = griddata[grid].rfmstruct;

  // Compute Hamiltonian constraint violation and store it at diagnostic_output_gfs
  compute_residual_all_points(commondata, params, rfmstruct, auxevol_gfs, y_n_gfs, diagnostic_output_gfs);
"""
    else:
        body += r"""
  // Set xx
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];

  // Compute Hamiltonian constraint violation and store it at diagnostic_output_gfs
  compute_residual_all_points(commondata, params, xx, auxevol_gfs, y_n_gfs, diagnostic_output_gfs);
"""

    if parallelization in ["cuda"]:
        body += r"""if (n_step % outevery == 0) {
    BHAH_DEVICE_SYNC();
    size_t streamid = params->grid_idx % NUM_STREAMS;
    cpyDevicetoHost__gf(commondata, params, host_diag_gfs, diagnostic_output_gfs, RESIDUAL_HGF, RESIDUAL_HGF, streamid);
  }"""

    body += r"""
  // Set integration radius for l2-norm computation
  const REAL integration_radius = 1000;

  // Compute l2-norm of Hamiltonian constraint violation
  REAL residual_H;
  compute_L2_norm_of_gridfunction(commondata, &griddata[grid].params, griddata[grid].xx, integration_radius, RESIDUAL_HGF, &residual_H, diagnostic_output_gfs);
  global_norm = MAX(global_norm, residual_H);
  } // END for(grid=0; grid<commondata->NUMGRIDS; ++grid)

  // Update residual to be used in stop condition
  commondata->log10_current_residual = global_norm;

  // Output l2-norm of Hamiltonian constraint violation to file
  {
    char filename[256];
    sprintf(filename, "residual_l2_norm.txt");
    FILE *outfile = (n_step == 0) ? fopen(filename, "w") : fopen(filename, "a");
    if (!outfile) {
      fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
      exit(1);
    }
    fprintf(outfile, "%6d %10.4e %.17e\n", n_step, time, global_norm);
    fclose(outfile);
  }

  if (n_step % outevery == 0) {
    // Only consider a single grid for now.
    const int grid = 0;
    params_struct *restrict params = &griddata[grid].params;
""".replace(
        "diagnostic_output_gfs)",
        (
            "diagnostic_output_gfs, diagnostic_output_gfs2)"
            if parallelization in ["cuda"]
            else "diagnostic_output_gfs)"
        ),
    )

    host_griddata = "griddata_host" if parallelization in ["cuda"] else "griddata"
    body += rf"""// Set reference metric grid xx
    REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
        xx[ww] = {host_griddata}[grid].xx[ww];
"""
    if parallelization in ["cuda"]:
        body += r"""
    // Sync with device when relevant
    BHAH_DEVICE_SYNC();
"""
    body += rf"""
    // 1D output
    diagnostics_nearest_1d_y_axis(commondata, params, xx, &{host_griddata}[grid].gridfuncs);
    diagnostics_nearest_1d_z_axis(commondata, params, xx, &{host_griddata}[grid].gridfuncs);

    // 2D output
    diagnostics_nearest_2d_xy_plane(commondata, params, xx, &{host_griddata}[grid].gridfuncs);
    diagnostics_nearest_2d_yz_plane(commondata, params, xx, &{host_griddata}[grid].gridfuncs);
  }}
"""
    if enable_progress_indicator:
        body += "progress_indicator(commondata, griddata);"
    body += r"""
  if (commondata->time + commondata->dt > commondata->t_final)
    printf("\n");
"""

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


# Define function to evaluate stop conditions
def register_CFunction_check_stop_conditions() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function to evaluate stop conditions.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h"]
    desc = "Evaluate stop conditions."
    cfunc_type = "void"
    name = "check_stop_conditions"
    params = (
        """commondata_struct *restrict commondata, griddata_struct *restrict griddata"""
    )

    # Register a parameter to stop hyperbolic relaxation
    _stop_relaxation = par.register_CodeParameter(
        "bool", __name__, "stop_relaxation", False, commondata=True
    )

    # Register parameter that sets the total number of relaxation steps
    _nn_max = par.register_CodeParameter("int", __name__, "nn_max", 0, commondata=True)

    # Register parameter that sets the tolerance for log of residual
    _log10_residual_tolerance = par.register_CodeParameter(
        "REAL", __name__, "log10_residual_tolerance", -15.8, commondata=True
    )

    # Register parameter that sets log of residual to be updated at every time step
    _log10_current_residual = par.register_CodeParameter(
        "REAL", __name__, "log10_current_residual", 1.0, commondata=True
    )

    body = r"""  // Since this version of NRPyElliptic is unigrid, we simply set the grid index to 0
  const int grid = 0;

  // Set params
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"

  // Check if total number of iteration steps has been reached
  if ((nn >= nn_max) || (log10_current_residual < log10_residual_tolerance)){
    printf("\nExiting main loop after %8d iterations\n", nn);
    printf("The tolerance for the logarithmic residual is %.8e\n", log10_residual_tolerance);
    printf("Exiting relaxation with logarithmic residual of %.8e\n", log10_current_residual);
    commondata->stop_relaxation = true;
  }
"""

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,  # set_CodeParameters.h is manually included after the declaration of params_struct *restrict params
        body=body,
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
