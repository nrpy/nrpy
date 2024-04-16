"""
Library of parallelization paradigm agnostic base classes.
Base classes are used to support generating C functions for solving the 
hyperbolic relaxation equation in curvilinear coordinates, using a reference-metric formalism.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
         Samuel D. Tootle; sdtootle **at** gmail **dot** com
"""

from typing import Union, Tuple, Dict

import sympy as sp
import nrpy.grid as gri
import nrpy.params as par
import nrpy.reference_metric as refmetric
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc

from nrpy.equations.nrpyelliptic.ConformallyFlat_RHSs import (
    HyperbolicRelaxationCurvilinearRHSs,
)
from nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms import (
    compute_psi_background_and_ADD_times_AUU,
)


# Define functions to set up initial guess
class base_register_CFunction_initial_guess_single_point:
    """
    Base class for generating the initial guess of solution at a single point.

    :param fp_type: Floating point type, e.g., "double".
    :return: None.
    """

    def __init__(self, fp_type: str = "double") -> None:
        self.fp_type = fp_type
        self.includes = ["BHaH_defines.h"]

        self.desc = r"""Compute initial guess at a single point."""
        self.cfunc_type = "void"
        self.name = "initial_guess_single_point"
        self.params = r"""const commondata_struct *restrict commondata, const params_struct *restrict params,
        const REAL xx0, const REAL xx1, const REAL xx2,  REAL *restrict uu_ID, REAL *restrict vv_ID
"""


class base_register_CFunction_initial_guess_all_points:
    """
    Base class for generating the initial guess function for the hyperbolic relaxation equation.

    :param enable_checkpointing: Attempt to read from a checkpoint file before generating initial guess.
    :param fp_type: Floating point type, e.g., "double".
    :return: None.
    """

    def __init__(
        self, enable_checkpointing: bool = False, fp_type: str = "double"
    ) -> None:
        self.fp_type = fp_type
        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

        self.desc = r"""Set initial guess to solutions of hyperbolic relaxation equation at all points."""
        self.cfunc_type = "void"
        self.name = "initial_data"
        self.params = (
            "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
        )
        self.uu_gf_memaccess = gri.BHaHGridFunction.access_gf("uu")
        self.vv_gf_memaccess = gri.BHaHGridFunction.access_gf("vv")
        self.body = ""
        if enable_checkpointing:
            self.body += """// Attempt to read checkpoint file. If it doesn't exist, then continue. Otherwise return.
if( read_checkpoint(commondata, griddata) ) return;
"""


# Define functions to set AUXEVOL gridfunctions


class base_register_CFunction_auxevol_gfs_single_point:
    """
    Base class for generating the function to compute AUXEVOL grid functions at a single point.

    :param CoordSystem: The coordinate system to use in setting up the AUXEVOL gridfunctions.
    :param fp_type: Floating point type, e.g., "double".
    :return: None.
    """

    def __init__(
        self,
        CoordSystem: str,
        fp_type: str = "double",
    ) -> None:
        self.fp_type = fp_type
        # Compute psi_background and ADD_times_AUU
        self.psi_background, self.ADD_times_AUU = (
            compute_psi_background_and_ADD_times_AUU(CoordSystem)
        )

        self.includes = ["BHaH_defines.h"]

        self.desc = r"""Compute AUXEVOL grid functions at a single point."""
        self.cfunc_type = "void"
        self.name = "auxevol_gfs_single_point"
        self.params = r"""const commondata_struct *restrict commondata, const params_struct *restrict params,
    const REAL xx0, const REAL xx1, const REAL xx2,  REAL *restrict psi_background, REAL *restrict ADD_times_AUU
"""


class base_register_CFunction_auxevol_gfs_all_points:
    """
    Base class for generating the function for the AUXEVOL grid functions at all points.

    :param fp_type: Floating point type, e.g., "double".
    :return: None.
    """

    def __init__(
        self,
        fp_type: str = "double",
    ) -> None:

        self.fp_type = fp_type

        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

        self.desc = r"""Set AUXEVOL gridfunctions at all points."""
        self.cfunc_type = "void"
        self.name = "auxevol_gfs_all_points"
        self.params = (
            "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
        )
        self.psi_background_memaccess = gri.BHaHGridFunction.access_gf("psi_background")
        self.ADD_times_AUU_memaccess = gri.BHaHGridFunction.access_gf("ADD_times_AUU")


class base_register_CFunction_variable_wavespeed_gfs_all_points:
    """
    Base class for generating the function to compute variable wavespeed.
    Calculation is based on the local grid spacing for a single coordinate system.

    :param CoordSystem: The coordinate system to use in the hyperbolic relaxation.
    :param fp_type: Floating point type, e.g., "double".

    :return: None.
    """

    def __init__(
        self,
        CoordSystem: str,
        fp_type: str = "double",
    ) -> None:
        self.fp_type = fp_type
        self.CoordSystem = CoordSystem

        self.includes = ["BHaH_defines.h"]
        self.desc = (
            "Compute variable wavespeed for all grids based on local grid spacing."
        )
        self.cfunc_type = "void"
        self.name = "variable_wavespeed_gfs_all_points"
        self.params = (
            "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
        )

        rfm = refmetric.reference_metric[self.CoordSystem]
        dxx0, dxx1, dxx2 = sp.symbols("dxx0 dxx1 dxx2", real=True)
        self.dsmin_computation_str = ccg.c_codegen(
            [
                rfm.scalefactor_orthog[0] * dxx0,
                rfm.scalefactor_orthog[1] * dxx1,
                rfm.scalefactor_orthog[2] * dxx2,
            ],
            ["const REAL dsmin0", "const REAL dsmin1", "const REAL dsmin2"],
            include_braces=False,
            fp_type=self.fp_type,
        )

        variable_wavespeed_memaccess = gri.BHaHGridFunction.access_gf(
            "variable_wavespeed"
        )

        self.dsmin_computation_str += f"""\n// Set local wavespeed
            {variable_wavespeed_memaccess} = MINIMUM_GLOBAL_WAVESPEED * MIN(dsmin0, MIN(dsmin1, dsmin2)) / dt;\n"""


class base_register_CFunction_initialize_constant_auxevol:
    """
    Base class for generating the function to call all functions that set up AUXEVOL gridfunctions.

    :return: None
    """

    def __init__(self) -> None:
        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

        self.desc = r"""Call functions that set up all AUXEVOL gridfunctions."""
        self.cfunc_type = "void"
        self.name = "initialize_constant_auxevol"
        self.params = (
            "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
        )


# Define function to compute the l^2 of a gridfunction
class base_register_CFunction_compute_L2_norm_of_gridfunction:
    """
    Base class for generating the function to compute l2-norm of a gridfunction assuming a single grid.

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
        self.fp_type = fp_type
        self.includes = ["BHaH_defines.h"]
        self.desc = "Compute l2-norm of a gridfunction assuming a single grid."
        self.cfunc_type = "REAL"
        self.name = "compute_L2_norm_of_gridfunction"
        self.params = """commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                    const REAL integration_radius, const int gf_index, const REAL *restrict in_gf"""

        self.rfm = refmetric.reference_metric[CoordSystem]
        self.body = ""


# Define diagnostics function
class base_register_CFunction_diagnostics:
    """
    Base class for generating the function for simulation diagnostics.

    :param default_diagnostics_out_every: Specifies the default diagnostics output frequency.
    :param out_quantities_dict: Dictionary or string specifying output quantities.
    :raises TypeError: If `out_quantities_dict` is not a dictionary and not set to "default".
    """

    def __init__(
        self,
        default_diagnostics_out_every: int,
        out_quantities_dict: Union[str, Dict[Tuple[str, str], str]] = "default",
    ) -> None:
        _ = par.CodeParameter(
            "int",
            __name__,
            "diagnostics_output_every",
            default_diagnostics_out_every,
            commondata=True,
        )

        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
        self.desc = "Diagnostics."
        self.cfunc_type = "void"
        self.name = "diagnostics"
        self.params = (
            "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
        )
        self.out_quantities_dict = out_quantities_dict

        # fmt: off
        if self.out_quantities_dict == "default":
            self.out_quantities_dict = {
                ("REAL", "numUU"): "y_n_gfs[IDX4pt(UUGF, idx3)]",
                ("REAL", "log10ResidualH"): "log10(fabs(diagnostic_output_gfs[IDX4pt(RESIDUAL_HGF, idx3)] + 1e-16))",
            }
        if not isinstance(self.out_quantities_dict, dict):
            raise TypeError(f"out_quantities_dict was initialized to {self.out_quantities_dict}, which is not a dictionary!")
        # fmt: on


# Define function to evaluate stop conditions
class base_register_CFunction_check_stop_conditions:
    """
    Base class for generating the function to evaluate stop conditions.

    :return: None.
    """

    def __init__(self) -> None:
        self.includes = ["BHaH_defines.h"]
        self.desc = "Evaluate stop conditions."
        self.cfunc_type = "void"
        self.name = "check_stop_conditions"
        self.params = """commondata_struct *restrict commondata, griddata_struct *restrict griddata"""

        # Register a parameter to stop hyperbolic relaxation
        _stop_relaxation = par.register_CodeParameter(
            "bool", __name__, "stop_relaxation", False, commondata=True
        )

        # Register parameter that sets the total number of relaxation steps
        _nn_max = par.register_CodeParameter(
            "int", __name__, "nn_max", 0, commondata=True
        )

        # Register parameter that sets the tolerance for log of residual
        _log10_residual_tolerance = par.register_CodeParameter(
            "REAL", __name__, "log10_residual_tolerance", -15.8, commondata=True
        )

        # Register parameter that sets log of residual to be updated at every time step
        _log10_current_residual = par.register_CodeParameter(
            "REAL", __name__, "log10_current_residual", 1.0, commondata=True
        )

        self.body = r"""  // Since this version of NRPyElliptic is unigrid, we simply set the grid index to 0
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
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            CoordSystem_for_wrapper_func="",
            name=self.name,
            params=self.params,
            include_CodeParameters_h=False,  # set_CodeParameters.h is manually included after the declaration of params_struct *restrict params
            body=self.body,
        )


# Define function to evaluate RHSs
class base_register_CFunction_rhs_eval:
    """
    Base class for generating the right-hand side (RHS) evaluation function for the hyperbolic relaxation equation.

    This function sets the right-hand side of the hyperbolic relaxation equation according to the
    selected coordinate system and specified parameters.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :return: None.
    """

    def __init__(
        self,
        CoordSystem: str,
        enable_rfm_precompute: bool,
    ) -> None:
        self.includes = ["BHaH_defines.h"]
        self.desc = r"""Set RHSs for hyperbolic relaxation equation."""
        self.cfunc_type = "void"
        self.name = "rhs_eval"
        self.params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict rhs_gfs"
        if enable_rfm_precompute:
            self.params = self.params.replace(
                "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
            )
        # Populate uu_rhs, vv_rhs
        self.rhs = HyperbolicRelaxationCurvilinearRHSs(
            CoordSystem, enable_rfm_precompute
        )


# Define function to compute residual the solution
class base_register_CFunction_compute_residual_all_points:
    """
    Base class for generating the residual evaluation function.

    This function sets the residual of the Hamiltonian constraint in the hyperbolic
    relaxation equation according to the selected coordinate system and specified
    parameters.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :return: None.
    """

    def __init__(
        self,
        CoordSystem: str,
        enable_rfm_precompute: bool,
    ) -> None:
        self.includes = ["BHaH_defines.h"]

        self.desc = r"""Compute residual of the Hamiltonian constraint for the hyperbolic relaxation equation."""
        self.cfunc_type = "void"
        self.name = "compute_residual_all_points"
        self.params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
                    REAL *restrict xx[3], const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs,
                    REAL *restrict aux_gfs"""
        if enable_rfm_precompute:
            self.params = self.params.replace(
                "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
            )
        # Populate residual_H
        self.rhs = HyperbolicRelaxationCurvilinearRHSs(
            CoordSystem, enable_rfm_precompute
        )
