"""
Base classes for numerical_grids_and_timestep() C function.

The parallelization modules will generate functions 
to set up numerical grids for use within the BHaH infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Samuel D. Tootle
        sdtootle **at** gmail **dot* com
"""

from typing import Dict, List
import sympy as sp
import nrpy.params as par
import nrpy.reference_metric as refmetric
import nrpy.c_codegen as ccg
from nrpy.helpers.expr_tree import get_unique_expression_symbols

class base_register_CFunction_numerical_grid_params_Nxx_dxx_xx:
    """
    Base class for generating the function to Set up a cell-centered grid of size grid_physical_size.
    Set params: Nxx, Nxx_plus_2NGHOSTS, dxx, invdxx, and xx.

    :param CoordSystem: The coordinate system used for the simulation.
    :param grid_physical_size: The physical size of the grid.
    :param Nxx_dict: A dictionary that maps coordinate systems to lists containing the number of grid points along each direction.

    :raises ValueError: If CoordSystem is not in Nxx_dict.
    """

    def __init__(
        self,
        CoordSystem: str,
        grid_physical_size: float,
        Nxx_dict: Dict[str, List[int]],
    ) -> None:

        self.CoordSystem = CoordSystem
        self.grid_physical_size = grid_physical_size
        self.Nxx_dict = Nxx_dict
        if self.CoordSystem not in self.Nxx_dict:
            raise ValueError(
                f"{CoordSystem} is not in Nxx_dict = {self.Nxx_dict}. Please add it."
            )
        for dirn in range(3):
            par.adjust_CodeParam_default(
                f"Nxx{dirn}", self.Nxx_dict[self.CoordSystem][dirn]
            )
        par.adjust_CodeParam_default("grid_physical_size", self.grid_physical_size)
        par.adjust_CodeParam_default("CoordSystemName", self.CoordSystem)
        self.rfm = refmetric.reference_metric[self.CoordSystem]
        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
        self.desc = f"""
Initializes a cell-centered grid in {self.CoordSystem} coordinates based on physical dimensions (grid_physical_size).

Inputs:
- Nx[] inputs: Specifies new grid dimensions, if needed.
- grid_is_resized: Indicates if the grid has been manually resized, triggering adjustments to grid parameters.
- convergence_factor: Multiplier for grid dimensions to refine resolution, applied only if grid hasn't been resized.

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
        self.params = "const commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3], const int Nx[3], const bool grid_is_resized"
        self.body = "// Start by setting default values for Nxx.\n"
        for dirn in range(3):
            self.body += (
                f"params->Nxx{dirn} = {self.Nxx_dict[self.CoordSystem][dirn]};\n"
            )
        self.body += """
// If all components of Nx[] are set to reasonable values (i.e., not -1), then set params->Nxx{} to Nx[].
if( !(Nx[0]==-1 || Nx[1]==-1 || Nx[2]==-1) ) {
"""

        for dirn in range(3):
            self.body += f"params->Nxx{dirn} = Nx[{dirn}];\n"
        self.body += f"""}}
snprintf(params->CoordSystemName, 50, "{CoordSystem}");

if( !grid_is_resized ) {{
  // convergence_factor does not increase resolution across an axis of symmetry:
  if(params->Nxx0 != 2) params->Nxx0 *= commondata->convergence_factor;
  if(params->Nxx1 != 2) params->Nxx1 *= commondata->convergence_factor;
  if(params->Nxx2 != 2) params->Nxx2 *= commondata->convergence_factor;
}}

params->Nxx_plus_2NGHOSTS0 = params->Nxx0 + 2*NGHOSTS;
params->Nxx_plus_2NGHOSTS1 = params->Nxx1 + 2*NGHOSTS;
params->Nxx_plus_2NGHOSTS2 = params->Nxx2 + 2*NGHOSTS;
#include "../set_CodeParameters.h"

// Set grid size to grid_physical_size (set above, based on params->grid_physical_size):
"""

        # Set grid_physical_size & grid_hole_radius
        self.body += """{
        // Set grid size to a function of grid_physical_size, just set in set_CodeParameters.h above:
        """
        for key, value in self.rfm.grid_physical_size_dict.items():
            self.body += f"params->{key} = {value};\n"
        self.body += "}\n"

        # Set xxmin and xxmax
        self.body += """if( !grid_is_resized ) {
    """
        self.body += "// Set xxmin, xxmax\n"
        if "Cartesian" in CoordSystem:
            self.body += """
    REAL factor_rescale_minmax_xx0 = 1.0;
    REAL factor_rescale_minmax_xx1 = 1.0;
    REAL factor_rescale_minmax_xx2 = 1.0;
    if(params->Nxx0 < params->Nxx1 || params->Nxx0 < params->Nxx2)
    factor_rescale_minmax_xx0 = ((REAL)params->Nxx0) / ((REAL)MAX(params->Nxx1, params->Nxx2));
    if(params->Nxx1 < params->Nxx2 || params->Nxx1 < params->Nxx0)
    factor_rescale_minmax_xx1 = ((REAL)params->Nxx1) / ((REAL)MAX(params->Nxx2, params->Nxx0));
    if(params->Nxx2 < params->Nxx0 || params->Nxx2 < params->Nxx1)
    factor_rescale_minmax_xx2 = ((REAL)params->Nxx2) / ((REAL)MAX(params->Nxx0, params->Nxx1));
    """
        for minmax in ["min", "max"]:
            for dirn in range(3):
                rfm_value = (
                    self.rfm.xxmin[dirn] if minmax == "min" else self.rfm.xxmax[dirn]
                )
                str_rfm_value = str(rfm_value)
                param_dir = f"params->xx{minmax}{dirn}"

                if str_rfm_value in par.glb_code_params_dict:
                    cparam_type = par.glb_code_params_dict[str_rfm_value].cparam_type
                    if cparam_type != "#define":
                        self.body += f"{param_dir} = params->{rfm_value};\n"
                        continue

                if "Cartesian" in self.CoordSystem:
                    self.body += (
                        f"{param_dir} = {rfm_value} * factor_rescale_minmax_xx{dirn};\n"
                    )
                else:
                    self.body += f"{param_dir} = {rfm_value};\n"

        self.body += """}

    params->dxx0 = (params->xxmax0 - params->xxmin0) / ((REAL)params->Nxx0);
    params->dxx1 = (params->xxmax1 - params->xxmin1) / ((REAL)params->Nxx1);
    params->dxx2 = (params->xxmax2 - params->xxmin2) / ((REAL)params->Nxx2);

    params->invdxx0 = ((REAL)params->Nxx0) / (params->xxmax0 - params->xxmin0);
    params->invdxx1 = ((REAL)params->Nxx1) / (params->xxmax1 - params->xxmin1);
    params->invdxx2 = ((REAL)params->Nxx2) / (params->xxmax2 - params->xxmin2);
"""


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
        self.params = "commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3], bc_struct *restrict bcstruct"

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
        self.unique_symbols = list()
        for expr in self.min_expressions:
            sub_list = get_unique_expression_symbols(expr, exclude=[f'xx{i}' for i in range(3)])
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


class base_register_CFunction_numerical_grids_and_timestep:
    """
    Base class for generating the function to set up all numerical grids and timestep.

    The function configures the numerical grids based on given parameters, specifically
    focusing on the usage of reference metric precomputations and curvilinear boundary
    conditions.

    :param list_of_CoordSystems: List of CoordSystems
    :param enable_rfm_precompute: Whether to enable reference metric precomputation (default: False).
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions (default: False).
    """

    def __init__(
        self,
        list_of_CoordSystems: List[str],
        enable_rfm_precompute: bool = False,
        enable_CurviBCs: bool = False,
    ) -> None:

        self.list_of_CoordSystems = list_of_CoordSystems
        self.enable_rfm_precompute = enable_rfm_precompute
        self.enable_CurviBCs = enable_CurviBCs

        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
        self.desc = "Set up a cell-centered grids of size grid_physical_size."
        self.cfunc_type = "void"
        self.name = "numerical_grids_and_timestep"
        self.params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata, bool calling_for_first_time"
        self.body = ""
