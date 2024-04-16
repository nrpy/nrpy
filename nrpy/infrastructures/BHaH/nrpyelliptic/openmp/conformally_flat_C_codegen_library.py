"""
Library of C functions for solving the hyperbolic relaxation equation in curvilinear coordinates, using a reference-metric formalism.
using openmp parallelization

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
         Samuel D. Tootle; sdtootle **at** gmail **dot** com
"""

from typing import Union, cast, Tuple, Dict
from inspect import currentframe as cfr
from types import FrameType as FT
from pathlib import Path
from inspect import currentframe as cf

import sympy as sp
import nrpy.grid as gri
import nrpy.params as par
import nrpy.reference_metric as refmetric
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc

import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.nrpyelliptic.base_conformally_flat_C_codegen_library as base_npe_classes

import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.infrastructures.BHaH.diagnostics.output_0d_1d_2d_nearest_gridpoint_slices as out012d

# Define functions to set up initial guess


class openmp_register_CFunction_initial_guess_single_point(
    base_npe_classes.base_register_CFunction_initial_guess_single_point
):

    def __init__(self, fp_type: str = "double") -> None:
        super().__init__(fp_type=fp_type)

        self.body = ccg.c_codegen(
            [sp.sympify(0), sp.sympify(0)],
            ["*uu_ID", "*vv_ID"],
            verbose=False,
            include_braces=False,
            fp_type=self.fp_type,
        )


def register_CFunction_initial_guess_single_point(fp_type: str = "double") -> Union[None, pcg.NRPyEnv_type]:
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None

    self_class = openmp_register_CFunction_initial_guess_single_point(fp_type=fp_type)
    cfc.register_CFunction(
        includes=self_class.includes,
        desc=self_class.desc,
        cfunc_type=self_class.cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=self_class.name,
        params=self_class.params,
        include_CodeParameters_h=True,
        body=self_class.body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

class openmp_register_CFunction_initial_guess_all_points(
    base_npe_classes.base_register_CFunction_initial_guess_all_points
):
    """
    Register the initial guess function for the hyperbolic relaxation equation.

    :param enable_checkpointing: Attempt to read from a checkpoint file before generating initial guess.
    :param OMP_collapse: Degree of OpenMP loop collapsing.
    :param fp_type: Floating point type, e.g., "double".

    :return: None if in registration phase, else the updated NRPy environment.
    """

    def __init__(
        self,
        OMP_collapse: int,
        enable_checkpointing: bool = False,
        fp_type: str = "double",
    ) -> None:
        super().__init__(enable_checkpointing, fp_type=fp_type)

        self.body += r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];
  REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
"""
        self.body += lp.simple_loop(
            loop_body="initial_guess_single_point(commondata, params, xx0,xx1,xx2,"
            f"&{self.uu_gf_memaccess},"
            f"&{self.vv_gf_memaccess});",
            read_xxs=True,
            loop_region="all points",
            OMP_collapse=OMP_collapse,
            fp_type=self.fp_type,
        )
        self.body += "}\n"


def register_CFunction_initial_guess_all_points(
    OMP_collapse: int,
    enable_checkpointing: bool = False,
    fp_type: str = "double",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the initial guess function for the hyperbolic relaxation equation.

    :param enable_checkpointing: Attempt to read from a checkpoint file before generating initial guess.
    :param OMP_collapse: Degree of OpenMP loop collapsing.
    :param fp_type: Floating point type, e.g., "double".

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None    
    self_class = openmp_register_CFunction_initial_guess_all_points(OMP_collapse, enable_checkpointing, fp_type=fp_type)
    cfc.register_CFunction(
        includes=self_class.includes,
        desc=self_class.desc,
        cfunc_type=self_class.cfunc_type,
        name=self_class.name,
        params=self_class.params,
        include_CodeParameters_h=False,
        body=self_class.body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
# Define functions to set AUXEVOL gridfunctions
class openmp_register_CFunction_auxevol_gfs_single_point(
    base_npe_classes.base_register_CFunction_auxevol_gfs_single_point
):
    def __init__(
        self,
        CoordSystem: str,
        fp_type: str = "double",
    ) -> None:
        """
        Register the C function for the AUXEVOL grid functions at a single point.

        :param CoordSystem: The coordinate system to use in setting up the AUXEVOL gridfunctions.
        :param fp_type: Floating point type, e.g., "double".

        :return: None if in registration phase, else the updated NRPy environment.
        """
        super().__init__(CoordSystem, fp_type=fp_type)

        self.body = ccg.c_codegen(
            [self.psi_background, self.ADD_times_AUU],
            ["*psi_background", "*ADD_times_AUU"],
            verbose=False,
            include_braces=False,
            fp_type=fp_type,
        )

def register_CFunction_auxevol_gfs_single_point(
    CoordSystem: str,
    fp_type: str = "double",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for the AUXEVOL grid functions at a single point.

    :param CoordSystem: The coordinate system to use in setting up the AUXEVOL gridfunctions.
    :param fp_type: Floating point type, e.g., "double".

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None
    self_class = openmp_register_CFunction_auxevol_gfs_single_point(CoordSystem, fp_type=fp_type)
    cfc.register_CFunction(
        includes=self_class.includes,
        desc=self_class.desc,
        cfunc_type=self_class.cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=self_class.name,
        params=self_class.params,
        include_CodeParameters_h=True,
        body=self_class.body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
class openmp_register_CFunction_auxevol_gfs_all_points(
    base_npe_classes.base_register_CFunction_auxevol_gfs_all_points
):
    def __init__(
        self,
        OMP_collapse: int,
        fp_type: str = "double",
    ) -> None:
        """
        OpenMP overload to generate the C function for the AUXEVOL grid functions at all points.

        :param OMP_collapse: Degree of OpenMP loop collapsing.
        :param fp_type: Floating point type, e.g., "double".

        :return: None if in registration phase, else the updated NRPy environment.
        """
        super().__init__(fp_type=fp_type)

        self.body = r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];
  REAL *restrict in_gfs = griddata[grid].gridfuncs.auxevol_gfs;
"""
        self.body += lp.simple_loop(
            loop_body="auxevol_gfs_single_point(commondata, params, xx0,xx1,xx2,"
            f"&{self.psi_background_memaccess},"
            f"&{self.ADD_times_AUU_memaccess});",
            read_xxs=True,
            loop_region="all points",
            OMP_collapse=OMP_collapse,
            fp_type=self.fp_type,
        )
        self.body += "}\n"

def register_CFunction_auxevol_gfs_all_points(
    OMP_collapse: int,
    fp_type: str = "double",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for the AUXEVOL grid functions at all points.

    :param OMP_collapse: Degree of OpenMP loop collapsing.
    :param fp_type: Floating point type, e.g., "double".

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None
    self_class = openmp_register_CFunction_auxevol_gfs_all_points(OMP_collapse, fp_type=fp_type)
    cfc.register_CFunction(
        includes=self_class.includes,
        desc=self_class.desc,
        cfunc_type=self_class.cfunc_type,
        name=self_class.name,
        params=self_class.params,
        include_CodeParameters_h=False,
        body=self_class.body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
    
class openmp_register_CFunction_variable_wavespeed_gfs_all_points(
    base_npe_classes.base_register_CFunction_variable_wavespeed_gfs_all_points
):
    def __init__(
        self,
        CoordSystem: str,
        fp_type: str = "double",
    ) -> None:
        """
        OpenMP overload to generate function to compute variable wavespeed based on local grid spacing for a single coordinate system.

        :param CoordSystem: The coordinate system to use in the hyperbolic relaxation.
        :param fp_type: Floating point type, e.g., "double".

        :return: None if in registration phase, else the updated NRPy environment.
        """
        super().__init__(CoordSystem, fp_type=fp_type)

        self.body = r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];
  REAL *restrict in_gfs = griddata[grid].gridfuncs.auxevol_gfs;
"""

        self.body += lp.simple_loop(
            loop_body="\n" + self.dsmin_computation_str,
            read_xxs=True,
            loop_region="interior",
            fp_type=self.fp_type,
        )

        # We must close the loop that was opened in the line 'for(int grid=0; grid<commondata->NUMGRIDS; grid++) {'
        self.body += r"""} // END LOOP for(int grid=0; grid<commondata->NUMGRIDS; grid++)
                """

def register_CFunction_variable_wavespeed_gfs_all_points(
    CoordSystem: str,
    fp_type: str = "double",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function to compute variable wavespeed based on local grid spacing for a single coordinate system.

    :param CoordSystem: The coordinate system to use in the hyperbolic relaxation.
    :param fp_type: Floating point type, e.g., "double".

    :return: None if in registration phase, else the updated NRPy environment.
    """
    self_class=openmp_register_CFunction_variable_wavespeed_gfs_all_points(CoordSystem, fp_type=fp_type)

    cfc.register_CFunction(
        includes=self_class.includes,
        desc=self_class.desc,
        cfunc_type=self_class.cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=self_class.name,
        params=self_class.params,
        include_CodeParameters_h=False,  # set_CodeParameters.h is manually included after the declaration of params_struct *restrict params
        body=self_class.body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

class openmp_register_CFunction_initialize_constant_auxevol(
    base_npe_classes.base_register_CFunction_initialize_constant_auxevol
):

    def __init__(self) -> None:
        """
        Register function to call all functions that set up AUXEVOL gridfunctions.

        :return: None if in registration phase, else the updated NRPy environment.
        """
        super().__init__()

        self.body = r"""
        // Set up variable wavespeed
        variable_wavespeed_gfs_all_points(commondata, griddata);

        // Set up all other AUXEVOL gridfunctions
        auxevol_gfs_all_points(commondata, griddata);
        """


def register_CFunction_initialize_constant_auxevol() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function to call all functions that set up AUXEVOL gridfunctions.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None
    self_class=openmp_register_CFunction_initialize_constant_auxevol()
    
    cfc.register_CFunction(
        includes=self_class.includes,
        desc=self_class.desc,
        cfunc_type=self_class.cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=self_class.name,
        params=self_class.params,
        include_CodeParameters_h=False,
        body=self_class.body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
    
# Define function to compute the l^2 of a gridfunction
class register_CFunction_compute_L2_norm_of_gridfunction(
    base_npe_classes.base_register_CFunction_compute_L2_norm_of_gridfunction
):
    def __init__(
        self,
        CoordSystem: str,
        fp_type: str = "double",
    ) -> None:
        """
        Register function to compute l2-norm of a gridfunction assuming a single grid.

        Note that parallel codegen is disabled for this function, as it sometimes causes a
        multiprocess race condition on Python 3.6.7

        :param CoordSystem: the rfm coordinate system.
        :param fp_type: Floating point type, e.g., "double".
        """
        super().__init__(CoordSystem, fp_type=fp_type)

        loop_body = ccg.c_codegen(
            [
                self.rfm.xxSph[0],
                self.rfm.detgammahat,
            ],
            [
                "const REAL r",
                "const REAL sqrtdetgamma",
            ],
            include_braces=False,
            fp_type=self.fp_type,
        )

        loop_body += r"""
if(r < integration_radius) {
  const REAL gf_of_x = in_gf[IDX4(gf_index, i0, i1, i2)];
  const REAL dV = sqrtdetgamma * dxx0 * dxx1 * dxx2;
  squared_sum += gf_of_x * gf_of_x * dV;
  volume_sum  += dV;
} // END if(r < integration_radius)
"""
        self.body += r"""
  // Unpack grid parameters assuming a single grid
  const int grid = 0;
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"

  // Define reference metric grid
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];

  // Set summation variables to compute l2-norm
  REAL squared_sum = 0.0;
  REAL volume_sum  = 0.0;

"""

        self.body += lp.simple_loop(
            loop_body="\n" + loop_body,
            read_xxs=True,
            loop_region="interior",
            OMP_custom_pragma=r"#pragma omp parallel for reduction(+:squared_sum,volume_sum)",
            fp_type=self.fp_type,
        )

        self.body += r"""
  // Compute and output the log of the l2-norm.
  return log10(1e-16 + sqrt(squared_sum / volume_sum));  // 1e-16 + ... avoids log10(0)
"""

        cfc.register_CFunction(
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
class openmp_register_CFunction_diagnostics(
    base_npe_classes.base_register_CFunction_diagnostics
):

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
        """
        OpenMP overload to facilitate generating the C function for simulation diagnostics.

        :param CoordSystem: Coordinate system used.
        :param default_diagnostics_out_every: Specifies the default diagnostics output frequency.
        :param enable_progress_indicator: Whether to enable the progress indicator.
        :param axis_filename_tuple: Tuple containing filename and variables for axis output.
        :param plane_filename_tuple: Tuple containing filename and variables for plane output.
        :param out_quantities_dict: Dictionary or string specifying output quantities.
        :return: None if in registration phase, else the updated NRPy environment.
        """
        super().__init__(
            default_diagnostics_out_every,
            out_quantities_dict=out_quantities_dict,
        )
        
        # This has to be here to avoid type issues with mypy
        # An error will throw in super().__init__() if out_quantities_dict != dict
        self.out_quantities_dict: Dict[Tuple[str, str], str] = self.out_quantities_dict
        
        for axis in ["y", "z"]:
            out012d.register_CFunction_diagnostics_nearest_1d_axis(
                CoordSystem=CoordSystem,
                out_quantities_dict=self.out_quantities_dict,
                filename_tuple=axis_filename_tuple,
                axis=axis,
            )
        for plane in ["xy", "yz"]:
            out012d.register_CFunction_diagnostics_nearest_2d_plane(
                CoordSystem=CoordSystem,
                out_quantities_dict=self.out_quantities_dict,
                filename_tuple=plane_filename_tuple,
                plane=plane,
            )

        self.body = r"""  // Output progress to stderr
  progress_indicator(commondata, griddata);

  // Since this version of NRPyElliptic is unigrid, we simply set the grid index to 0
  const int grid = 0;

  // Set gridfunctions aliases
  REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
  REAL *restrict diagnostic_output_gfs = griddata[grid].gridfuncs.diagnostic_output_gfs;

  // Set params and rfm_struct
  params_struct *restrict params = &griddata[grid].params;
  const rfm_struct *restrict rfmstruct = &griddata[grid].rfmstruct;
#include "set_CodeParameters.h"

  // Compute Hamiltonian constraint violation and store it at diagnostic_output_gfs
  compute_residual_all_points(commondata, params, rfmstruct, auxevol_gfs, y_n_gfs, diagnostic_output_gfs);

  // Set integration radius for l2-norm computation
  const REAL integration_radius = 1000;

  // Compute l2-norm of Hamiltonian constraint violation
  const REAL residual_H = compute_L2_norm_of_gridfunction(commondata, griddata, integration_radius, RESIDUAL_HGF, diagnostic_output_gfs);

  // Update residual to be used in stop condition
  commondata->log10_current_residual = residual_H;

  // Output l2-norm of Hamiltonian constraint violation to file
  {
    char filename[256];
    sprintf(filename, "residual_l2_norm.txt");
    FILE *outfile = (nn == 0) ? fopen(filename, "w") : fopen(filename, "a");
    if (!outfile) {
      fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
      exit(1);
    }
    fprintf(outfile, "%6d %10.4e %.17e\n", nn, time, residual_H);
    fclose(outfile);
  }

  // Grid data output
  const int n_step = commondata->nn, outevery = commondata->diagnostics_output_every;
  if (n_step % outevery == 0) {
    // Set reference metric grid xx
    REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
        xx[ww] = griddata[grid].xx[ww];

    // 1D output
    diagnostics_nearest_1d_y_axis(commondata, params, xx, &griddata[grid].gridfuncs);
    diagnostics_nearest_1d_z_axis(commondata, params, xx, &griddata[grid].gridfuncs);

    // 2D output
    diagnostics_nearest_2d_xy_plane(commondata, params, xx, &griddata[grid].gridfuncs);
    diagnostics_nearest_2d_yz_plane(commondata, params, xx, &griddata[grid].gridfuncs);
  }
"""
        if enable_progress_indicator:
            self.body += "progress_indicator(commondata, griddata);"
        self.body += r"""
  if (commondata->time + commondata->dt > commondata->t_final)
    printf("\n");
"""

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
    Register C function for simulation diagnostics.

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
    self_class = openmp_register_CFunction_diagnostics(
        CoordSystem,
        default_diagnostics_out_every,
        enable_progress_indicator=enable_progress_indicator,
        axis_filename_tuple=axis_filename_tuple,
        plane_filename_tuple=plane_filename_tuple,
        out_quantities_dict=out_quantities_dict
    )
    cfc.register_CFunction(
        includes=self_class.includes,
        desc=self_class.desc,
        cfunc_type=self_class.cfunc_type,
        name=self_class.name,
        params=self_class.params,
        include_CodeParameters_h=False,
        body=self_class.body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

# Define function to evaluate stop conditions
def register_CFunction_check_stop_conditions() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function to evaluate stop conditions.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None
    self_class=base_npe_classes.base_register_CFunction_check_stop_conditions()
    cfc.register_CFunction(
        includes=self_class.includes,
        desc=self_class.desc,
        cfunc_type=self_class.cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=self_class.name,
        params=self_class.params,
        include_CodeParameters_h=False,  # set_CodeParameters.h is manually included after the declaration of params_struct *restrict params
        body=self_class.body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

# Define function to evaluate RHSs
class openmp_register_CFunction_rhs_eval(base_npe_classes.base_register_CFunction_rhs_eval):

    def __init__(
        self,
        CoordSystem: str,
        enable_rfm_precompute: bool,
        enable_simd: bool,
        OMP_collapse: int,
        fp_type: str = "double",
    ) -> None:
        """
        OpenMP overload to facilitate generating the right-hand side (RHS) evaluation function for the hyperbolic relaxation equation.

        This function sets the right-hand side of the hyperbolic relaxation equation according to the
        selected coordinate system and specified parameters.

        :param CoordSystem: The coordinate system.
        :param enable_rfm_precompute: Whether to enable reference metric precomputation.
        :param enable_simd: Whether to enable SIMD.
        :param OMP_collapse: Level of OpenMP loop collapsing.
        :param fp_type: Floating point type, e.g., "double".

        :return: None if in registration phase, else the updated NRPy environment.
        """
        super().__init__(
            CoordSystem, enable_rfm_precompute
        )

        if enable_simd:
            self.includes += [str(Path("simd") / "simd_intrinsics.h")]

        self.body = lp.simple_loop(
            loop_body=ccg.c_codegen(
                [self.rhs.uu_rhs, self.rhs.vv_rhs],
                [
                    gri.BHaHGridFunction.access_gf("uu", gf_array_name="rhs_gfs"),
                    gri.BHaHGridFunction.access_gf("vv", gf_array_name="rhs_gfs"),
                ],
                enable_fd_codegen=True,
                enable_simd=enable_simd,
                fp_type=fp_type,
            ),
            loop_region="interior",
            enable_simd=enable_simd,
            CoordSystem=CoordSystem,
            enable_rfm_precompute=enable_rfm_precompute,
            read_xxs=not enable_rfm_precompute,
            OMP_collapse=OMP_collapse,
            fp_type=fp_type,
        )


def register_CFunction_rhs_eval(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    OMP_collapse: int,
    fp_type: str = "double",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the right-hand side (RHS) evaluation function for the hyperbolic relaxation equation.

    This function sets the right-hand side of the hyperbolic relaxation equation according to the
    selected coordinate system and specified parameters.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_simd: Whether to enable SIMD.
    :param OMP_collapse: Level of OpenMP loop collapsing.
    :param fp_type: Floating point type, e.g., "double".

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None
    self_class=openmp_register_CFunction_rhs_eval(
        CoordSystem, enable_rfm_precompute, enable_simd, OMP_collapse, fp_type=fp_type
    )

    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=self_class.includes,
        desc=self_class.desc,
        cfunc_type=self_class.cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=self_class.name,
        params=self_class.params,
        body=self_class.body,
        enable_simd=enable_simd,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

# Define function to compute residual the solution
class openmp_register_CFunction_compute_residual_all_points(
    base_npe_classes.base_register_CFunction_compute_residual_all_points
):
    def __init__(
        self,
        CoordSystem: str,
        enable_rfm_precompute: bool,
        enable_simd: bool,
        OMP_collapse: int,
        fp_type: str = "double",
    ) -> None:
        """
        OpenMP overload to facilitate generating the residual evaluation function.

        This function sets the residual of the Hamiltonian constraint in the hyperbolic
        relaxation equation according to the selected coordinate system and specified
        parameters.

        :param CoordSystem: The coordinate system.
        :param enable_rfm_precompute: Whether to enable reference metric precomputation.
        :param enable_simd: Whether to enable SIMD.
        :param OMP_collapse: Level of OpenMP loop collapsing.
        :param fp_type: Floating point type, e.g., "double".

        :return: None if in registration phase, else the updated NRPy environment.
        """
        super().__init__(
            CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
        )
        if enable_simd:
            self.includes += [str(Path("simd") / "simd_intrinsics.h")]

        self.body = lp.simple_loop(
            loop_body=ccg.c_codegen(
                [self.rhs.residual],
                [
                    gri.BHaHGridFunction.access_gf(
                        "residual_H", gf_array_name="aux_gfs"
                    ),
                ],
                enable_fd_codegen=True,
                enable_simd=enable_simd,
                fp_type=fp_type,
            ),
            loop_region="interior",
            enable_simd=enable_simd,
            CoordSystem=CoordSystem,
            enable_rfm_precompute=enable_rfm_precompute,
            read_xxs=not enable_rfm_precompute,
            OMP_collapse=OMP_collapse,
            fp_type=fp_type,
        )

def register_CFunction_compute_residual_all_points(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    OMP_collapse: int,
    fp_type: str = "double",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the residual evaluation function.

    This function sets the residual of the Hamiltonian constraint in the hyperbolic
    relaxation equation according to the selected coordinate system and specified
    parameters.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_simd: Whether to enable SIMD.
    :param OMP_collapse: Level of OpenMP loop collapsing.
    :param fp_type: Floating point type, e.g., "double".

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None
    self_class = openmp_register_CFunction_compute_residual_all_points(CoordSystem,
                                                                       enable_rfm_precompute, 
                                                                       enable_simd, OMP_collapse,
                                                                       fp_type=fp_type)
    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=self_class.includes,
        desc=self_class.desc,
        cfunc_type=self_class.cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=self_class.name,
        params=self_class.params,
        body=self_class.body,
        enable_simd=enable_simd,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
