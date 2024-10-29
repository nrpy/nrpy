"""
Base classes for numerical_grids_and_timestep() C function.

The parallelization modules will generate functions
to set up numerical grids for use within the BHaH infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Samuel D. Tootle
        sdtootle **at** gmail **dot* com
"""

import os
from typing import Dict, List

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.helpers.expr_tree import get_unique_expression_symbols


class base_register_CFunction_numerical_grid_params_Nxx_dxx_xx:
    """
    Base class for generating the function to Set up a cell-centered grid of size grid_physical_size.
    Set params: Nxx, Nxx_plus_2NGHOSTS, dxx, invdxx, and xx.

    :param CoordSystem: The coordinate system used for the simulation.
    :param Nxx_dict: A dictionary that maps coordinate systems to lists containing the number of grid points along each direction.

    :raises ValueError: If CoordSystem is not in Nxx_dict.
    """

    def __init__(
        self,
        CoordSystem: str,
        Nxx_dict: Dict[str, List[int]],
    ) -> None:
        self.prefunc = ""
        self.CoordSystem = CoordSystem
        self.Nxx_dict = Nxx_dict
        if self.CoordSystem not in self.Nxx_dict:
            raise ValueError(
                f"{CoordSystem} is not in Nxx_dict = {self.Nxx_dict}. Please add it."
            )
        for dirn in range(3):
            par.adjust_CodeParam_default(f"Nxx{dirn}", Nxx_dict[CoordSystem][dirn])
        par.adjust_CodeParam_default("CoordSystemName", CoordSystem)

        for dirn in range(3):
            par.adjust_CodeParam_default(
                f"Nxx{dirn}", self.Nxx_dict[self.CoordSystem][dirn]
            )
        par.adjust_CodeParam_default("CoordSystemName", self.CoordSystem)
        self.rfm = refmetric.reference_metric[self.CoordSystem]
        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
        self.desc = f"""Initializes a cell-centered grid in {CoordSystem} coordinates based on physical dimensions (grid_physical_size).

Inputs:
- Nx[] inputs: Specifies new grid dimensions, if needed.
- params.convergence_factor (set to 1.0 by default): Factor by which grid resolution is increased; set to 1.0 by default.
- set_xxmin_xxmax_to_defaults: Whether to set xxmin[3], xxmax[3] to default values set in reference_metric.py.

Parameter outputs:
- Nxx: Number of grid points in each direction.
- Nxx_plus_2NGHOSTS: Total grid points including ghost zones.
- dxx: Grid spacing.
- invdxx: Inverse of grid spacing.

Grid setup output:
- xx: Coordinate values for each (cell-centered) grid point.
"""
        self.cfunc_type = "void"
        self.name = "numerical_grid_params_Nxx_dxx_xx"
        self.params = "const commondata_struct *restrict commondata, params_struct *restrict params, REAL * xx[3], const int Nx[3], const bool set_xxmin_xxmax_to_defaults"
        self.body = "// Set default values for the grid resolution in each dimension.\n"
        for dirn in range(3):
            self.body += (
                f"params->Nxx{dirn} = {self.Nxx_dict[self.CoordSystem][dirn]};\n"
            )
        self.body += """
// If all components of Nx[] are set to valid values (i.e., not -1), override the default values with Nx[].
if( Nx[0]!=-1 && Nx[1]!=-1 && Nx[2]!=-1 ) {
"""

        for dirn in range(3):
            self.body += f"params->Nxx{dirn} = Nx[{dirn}];\n"
        self.body += f"""}}
snprintf(params->CoordSystemName, 50, "{CoordSystem}");

// Resize grid by convergence_factor; used for convergence testing.
{{
  // convergence_factor does not increase resolution across an axis of symmetry (Nxx == 2):
  if(params->Nxx0 != 2) params->Nxx0 *= commondata->convergence_factor;
  if(params->Nxx1 != 2) params->Nxx1 *= commondata->convergence_factor;
  if(params->Nxx2 != 2) params->Nxx2 *= commondata->convergence_factor;
}}

// Set the full grid size; including the ghostzones (of width NGHOSTS) on the boundaries.
params->Nxx_plus_2NGHOSTS0 = params->Nxx0 + 2*NGHOSTS;
params->Nxx_plus_2NGHOSTS1 = params->Nxx1 + 2*NGHOSTS;
params->Nxx_plus_2NGHOSTS2 = params->Nxx2 + 2*NGHOSTS;
"""

        # Set grid_physical_size & grid_hole_radius
        self.body += """{
#include "../set_CodeParameters.h"
// Set grid size to a function of grid_physical_size (grid_physical_size set in set_CodeParameters.h above):
"""
        for key, value in self.rfm.grid_physical_size_dict.items():
            self.body += f"params->{key} = {value};\n"
        self.body += "}\n"

        # Set grid_hole_radius
        if "Holey" in self.CoordSystem or "Wedge" in self.CoordSystem:
            self.body += """{
#include "../set_CodeParameters.h"
// Set grid hole radius to a function of grid_hole_radius (grid_hole_radius set in set_CodeParameters.h above):
"""
            for key, value in self.rfm.grid_hole_radius_dict.items():
                self.body += f"params->{key} = {value};\n"
            self.body += "}\n"

        # Set minimum and maximum values of xx[][] for each grid.
        self.body += """if (set_xxmin_xxmax_to_defaults) {
#include "../set_CodeParameters.h"
// Set {xxmin[], xxmax[]} to default values, which could be functions of other rfm params (set in set_CodeParameters.h above):
"""
        for minmax in ["min", "max"]:
            for dirn in range(3):
                rfm_value = (
                    self.rfm.xxmin[dirn] if minmax == "min" else self.rfm.xxmax[dirn]
                )
                self.body += f"params->xx{minmax}{dirn} = {rfm_value};\n"
        self.body += "}\n"

        # Set quantities that depend on Nxx and {xxmin, xxmax}, then set up coordinate arrays xx[3][Nxxi].
        self.body += """
// Set quantities that depend on Nxx and {xxmin, xxmax}: dxx, invdxx.
params->dxx0 = (params->xxmax0 - params->xxmin0) / ((REAL)params->Nxx0);
params->dxx1 = (params->xxmax1 - params->xxmin1) / ((REAL)params->Nxx1);
params->dxx2 = (params->xxmax2 - params->xxmin2) / ((REAL)params->Nxx2);

params->invdxx0 = ((REAL)params->Nxx0) / (params->xxmax0 - params->xxmin0);
params->invdxx1 = ((REAL)params->Nxx1) / (params->xxmax1 - params->xxmin1);
params->invdxx2 = ((REAL)params->Nxx2) / (params->xxmax2 - params->xxmin2);
#include "../set_CodeParameters.h"
"""

    def register(self) -> None:
        """Register CFunction."""
        _, actual_name = cfc.function_name_and_subdir_with_CoordSystem(
            os.path.join("."), self.name, self.CoordSystem
        )
        if not actual_name in cfc.CFunction_dict:
            cfc.register_CFunction(
                prefunc=self.prefunc,
                includes=self.includes,
                desc=self.desc,
                cfunc_type=self.cfunc_type,
                CoordSystem_for_wrapper_func=self.CoordSystem,
                name=self.name,
                params=self.params,
                include_CodeParameters_h=False,  # keep this False or regret having to debug the mess.
                body=self.body,
            )


class base_register_CFunction_cfl_limited_timestep:
    """
    Base class for generating the function to find the CFL-limited timestep dt on a numerical grid.

    The timestep is determined by the relation dt = CFL_FACTOR * ds_min, where ds_min
    is the minimum spacing between neighboring gridpoints on a numerical grid.

    :param CoordSystem: The coordinate system used for the simulation.
    :param fp_type: Floating point type, e.g., "double".
    """

    def __init__(self, CoordSystem: str, fp_type: str = "double") -> None:

        self.CoordSystem = CoordSystem
        self.fp_type = fp_type

        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
        self.desc = (
            f"Output minimum gridspacing ds_min on a {CoordSystem} numerical grid."
        )
        self.cfunc_type = "void"
        self.name = "cfl_limited_timestep"
        self.params = "commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3]"
        self.prefunc = ""

        self.rfm = refmetric.reference_metric[CoordSystem]
        dxx0, dxx1, dxx2 = sp.symbols("dxx0 dxx1 dxx2", real=True)
        self.body = ""
        self.min_expressions = [
            sp.Abs(self.rfm.scalefactor_orthog[0] * dxx0),
            sp.Abs(self.rfm.scalefactor_orthog[1] * dxx1),
            sp.Abs(self.rfm.scalefactor_orthog[2] * dxx2),
        ]

        # Save unique symbols not related to coordinates in case
        # we need to include them in the function body
        self.unique_symbols = []
        for expr in self.min_expressions:
            sub_list = get_unique_expression_symbols(
                expr, exclude=[f"xx{i}" for i in range(3)]
            )
            self.unique_symbols += sub_list
        self.unique_symbols = sorted(list(set(self.unique_symbols)))

        self.min_body_compute = ccg.c_codegen(
            self.min_expressions,
            ["dsmin0", "dsmin1", "dsmin2"],
            include_braces=False,
            fp_type=fp_type,
        )
        self.min_body = self.min_body_compute
        self.min_body += """
  ds_min = MIN(ds_min, MIN(dsmin0, MIN(dsmin1, dsmin2)));
}
commondata->dt = MIN(commondata->dt, ds_min * commondata->CFL_FACTOR);
"""

    def register(self) -> None:
        """Register CFunction."""
        _, actual_name = cfc.function_name_and_subdir_with_CoordSystem(
            os.path.join("."), self.name, self.CoordSystem
        )
        if not actual_name in cfc.CFunction_dict:
            cfc.register_CFunction(
                prefunc=self.prefunc,
                includes=self.includes,
                desc=self.desc,
                cfunc_type=self.cfunc_type,
                CoordSystem_for_wrapper_func=self.CoordSystem,
                name=self.name,
                params=self.params,
                include_CodeParameters_h=True,
                body=self.body,
            )


class base_register_CFunction_numerical_grids_and_timestep:
    """
    Base class for generating the function to set up all numerical grids and timestep.

    The function configures the numerical grids based on given parameters, specifically
    focusing on the usage of reference metric precomputations and curvilinear boundary
    conditions.

    :param list_of_CoordSystems: List of CoordSystems
    :param list_of_grid_physical_sizes: List of grid_physical_size for each CoordSystem; needed for Independent grids.
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :param enable_rfm_precompute: Whether to enable reference metric precomputation (default: False).
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions (default: False).

    :raises ValueError: If invalid gridding_approach selected.
    """

    def __init__(
        self,
        list_of_CoordSystems: List[str],
        list_of_grid_physical_sizes: List[float],
        gridding_approach: str = "independent grid(s)",
        enable_rfm_precompute: bool = False,
        enable_CurviBCs: bool = False,
    ) -> None:

        self.list_of_CoordSystems = list_of_CoordSystems
        self.list_of_grid_physical_sizes = list_of_grid_physical_sizes
        self.enable_rfm_precompute = enable_rfm_precompute
        self.enable_CurviBCs = enable_CurviBCs
        self.gridding_approach = gridding_approach

        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
        self.desc = "Set up numerical grids and timestep."
        self.cfunc_type = "void"
        self.name = "numerical_grids_and_timestep"
        self.params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata, bool calling_for_first_time"
        self.body = r"""
  // Step 1.a: Set each CodeParameter in griddata.params to default, for MAXNUMGRIDS grids.
  params_struct_set_to_default(commondata, griddata);"""
        self.body += rf"""
  // Independent grids
  int Nx[3] = {{ -1, -1, -1 }};

  // Step 1.b: Set commondata->NUMGRIDS to number of CoordSystems we have
  commondata->NUMGRIDS = {len(list_of_CoordSystems)};
"""
        if self.gridding_approach == "independent grid(s)":
            self.body += """
  {
    // Step 1.c: For each grid, set Nxx & Nxx_plus_2NGHOSTS, as well as dxx, invdxx, & xx based on grid_physical_size
    const bool set_xxmin_xxmax_to_defaults = true;
    int grid=0;
"""
            for which_CoordSystem, CoordSystem in enumerate(self.list_of_CoordSystems):
                self.body += (
                    f"  griddata[grid].params.CoordSystem_hash = {CoordSystem.upper()};\n"
                    f"  griddata[grid].params.grid_idx = grid;\n"
                )
                self.body += f"griddata[grid].params.grid_physical_size = {list_of_grid_physical_sizes[which_CoordSystem]};\n"
                self.body += "numerical_grid_params_Nxx_dxx_xx(commondata, &griddata[grid].params, griddata[grid].xx, Nx, set_xxmin_xxmax_to_defaults);\n"
                self.body += "grid++;\n\n"
            self.body += "}\n"
        elif self.gridding_approach == "multipatch":
            # fmt: off
            _ = par.CodeParameter("char[200]", __name__, "multipatch_choice", "", commondata=True, add_to_parfile=True)
            # fmt: on
            for dirn in ["x", "y", "z"]:
                # Direction of unit vectors relative to original, accounting for accumulation of regrids."
                _ = par.register_CodeParameter(
                    "REAL",
                    __name__,
                    f"cumulative_regrid_{dirn}hatU[3]",
                    "unset",  # Set below in C code when calling_for_first_time.
                    commondata=True,
                    add_to_parfile=False,
                    add_to_set_CodeParameters_h=False,
                )
            self.body += """
  if(calling_for_first_time) {
    // Initialize rotation unit vectors
    // Set the x-hat unit vector (1, 0, 0)
    commondata->cumulative_regrid_xhatU[0] = 1;
    commondata->cumulative_regrid_xhatU[1] = 0;
    commondata->cumulative_regrid_xhatU[2] = 0;

    // Set the y-hat unit vector (0, 1, 0)
    commondata->cumulative_regrid_yhatU[0] = 0;
    commondata->cumulative_regrid_yhatU[1] = 1;
    commondata->cumulative_regrid_yhatU[2] = 0;

    // Set the z-hat unit vector (0, 0, 1)
    commondata->cumulative_regrid_zhatU[0] = 0;
    commondata->cumulative_regrid_zhatU[1] = 0;
    commondata->cumulative_regrid_zhatU[2] = 1;
  }
  // Step 1.c: Multipatch grid structures are set up algorithmically.
  multipatch_grids_set_up(commondata, griddata);
  """
        else:
            raise ValueError(
                f"""gridding_approach == "{gridding_approach}" not supported.
            Supported approaches include: "independent grid(s)" (default) and "multipatch"."""
            )
        self.body += "\n// Step 1.d: Allocate memory for and define reference-metric precomputation lookup tables\n"

    def register(self) -> None:
        """Register C function."""
        cfc.register_CFunction(
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=False,
            body=self.body,
        )
