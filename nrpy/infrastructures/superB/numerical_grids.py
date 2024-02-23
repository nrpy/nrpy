"""
Register numerical_grids_chare() C function, as well as functions called by this one.

These functions set up numerical grids for use within the superB infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from typing import Dict, List
import sympy as sp
import nrpy.c_function as cfc
import nrpy.reference_metric as refmetric
import nrpy.c_codegen as ccg

def register_CFunction_numerical_grid_params_Nxx_dxx_xx(
    CoordSystem: str, grid_physical_size: float, Nxx_dict: Dict[str, List[int]], Nchare_dict: Dict[str, List[int]]
) -> None:
    """
    Register a C function to Set up a cell-centered grid of size grid_physical_size.
       Set params: Nxx, Nxx_plus_2NGHOSTS, dxx, invdxx, and xx.

    :param CoordSystem: The coordinate system used for the simulation.
    :param grid_physical_size: The physical size of the grid.
    :param Nxx_dict: A dictionary that maps coordinate systems to lists containing the number of grid points along each direction.
    """        
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"Set up a cell-centered {CoordSystem} grid of size grid_physical_size. Set params: Nxx, Nxx_plus_2NGHOSTS, dxx, invdxx, and xx."
    c_type = "void"
    name = "numerical_grid_params_Nxx_dxx_xx_chare"
    params = "commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3], const int chare_index[3]"
    body = ""
    for dirn in range(3):
        NxxoverNchare = Nxx_dict[CoordSystem][dirn] // Nchare_dict[CoordSystem][dirn] if Nxx_dict[CoordSystem][dirn] % Nchare_dict[CoordSystem][dirn] == 0 else ValueError("Result of dividing Nxx by Nchare is not an integer")
        body += f"params->Nxx{dirn} = {NxxoverNchare};\n"
    body += rf"""
const REAL grid_physical_size = params->grid_physical_size;
snprintf(params->CoordSystemName, 50, "{CoordSystem}");

// convergence_factor does not increase resolution across an axis of symmetry:
if(params->Nxx0 != 2) params->Nxx0 *= commondata->convergence_factor;
if(params->Nxx1 != 2) params->Nxx1 *= commondata->convergence_factor;
if(params->Nxx2 != 2) params->Nxx2 *= commondata->convergence_factor;

params->Nxx_plus_2NGHOSTS0 = params->Nxx0 + 2*NGHOSTS;
params->Nxx_plus_2NGHOSTS1 = params->Nxx1 + 2*NGHOSTS;
params->Nxx_plus_2NGHOSTS2 = params->Nxx2 + 2*NGHOSTS;

// Set grid size to grid_physical_size (set above, based on params->grid_physical_size):
"""
    rfm = refmetric.reference_metric[CoordSystem]
    for key, value in rfm.grid_physical_size_dict.items():
        body += f"params->{key} = {value};\n"
    body += "\n// Set xxmin, xxmax\n"
    for minmax in ["min", "max"]:
        for dirn in range(3):
            dx = (rfm.xxmax[dirn] - rfm.xxmin[dirn])/Nxx_dict[CoordSystem][dirn]            
            str_chare_value = (
              f"{rfm.xxmin[dirn]} + ({dx} * ((REAL)params->Nxx{dirn}) * ((REAL)chare_index[{dirn}]))" 
              if minmax == "min"
              else f"{rfm.xxmin[dirn]} + ({dx} * ((REAL)params->Nxx{dirn}) * ((REAL)(chare_index[{dirn}] + 1)))"
            )
            param_dir = f"params->xx{minmax}{dirn}"          
            body += f"{param_dir} = {str_chare_value};\n"

    body += """
params->dxx0 = (params->xxmax0 - params->xxmin0) / ((REAL)params->Nxx0);
params->dxx1 = (params->xxmax1 - params->xxmin1) / ((REAL)params->Nxx1);
params->dxx2 = (params->xxmax2 - params->xxmin2) / ((REAL)params->Nxx2);

params->invdxx0 = ((REAL)params->Nxx0) / (params->xxmax0 - params->xxmin0);
params->invdxx1 = ((REAL)params->Nxx1) / (params->xxmax1 - params->xxmin1);
params->invdxx2 = ((REAL)params->Nxx2) / (params->xxmax2 - params->xxmin2);

// Set up cell-centered Cartesian coordinate grid, centered at the origin.
xx[0] = (REAL *restrict)malloc(sizeof(REAL)*params->Nxx_plus_2NGHOSTS0);
xx[1] = (REAL *restrict)malloc(sizeof(REAL)*params->Nxx_plus_2NGHOSTS1);
xx[2] = (REAL *restrict)malloc(sizeof(REAL)*params->Nxx_plus_2NGHOSTS2);
for (int j = 0; j < params->Nxx_plus_2NGHOSTS0; j++) xx[0][j] = params->xxmin0 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * params->dxx0;
for (int j = 0; j < params->Nxx_plus_2NGHOSTS1; j++) xx[1][j] = params->xxmin1 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * params->dxx1;
for (int j = 0; j < params->Nxx_plus_2NGHOSTS2; j++) xx[2][j] = params->xxmin2 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * params->dxx2;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        c_type=c_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,  # keep this False or regret having to debug the mess.
        body=body,
    )


def register_CFunction_numerical_grids(
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
) -> None:
    """
    Register a C function to set up all numerical grids and timestep.

    The function configures the numerical grids based on given parameters, specifically
    focusing on the usage of reference metric precomputations and curvilinear boundary
    conditions.

    :param enable_rfm_precompute: Whether to enable reference metric precomputation (default: False).
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions (default: False).
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Set up a cell-centered grids of size grid_physical_size."
    c_type = "void"
    name = "numerical_grids_chare"
    params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    body = r"""// Step 1.b: Set Nxx & Nxx_plus_2NGHOSTS, as well as dxx, invdxx, & xx based on grid_physical_size
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  numerical_grid_params_Nxx_dxx_xx_chare(commondata, &griddata[grid].params, griddata[grid].xx);
}

// Step 1.c: Allocate memory for and define reference-metric precomputation lookup tables
"""
    if enable_rfm_precompute:
        body += r"""
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  rfm_precompute_malloc(commondata, &griddata[grid].params, &griddata[grid].rfmstruct);
  rfm_precompute_defines(commondata, &griddata[grid].params, &griddata[grid].rfmstruct, griddata[grid].xx);
}
"""
    else:
        body += "// (reference-metric precomputation disabled)\n"
    body += "\n// Step 1.d: Set up curvilinear boundary condition struct (bcstruct)\n"
    if enable_CurviBCs:
        body += r"""
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  bcstruct_set_up_chare(commondata, &griddata[grid].params, griddata[grid].xx, &griddata[grid].bcstruct);
}
"""
    else:
        body += "// (curvilinear boundary conditions bcstruct disabled)\n"
    
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


def register_CFunctions(
    list_of_CoordSystems: List[str],
    grid_physical_size: float,
    Nxx_dict: Dict[str, List[int]],
    Nchare_dict: Dict[str, List[int]],
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
) -> None:
    """
    Register C functions related to coordinate systems and grid parameters.

    :param list_of_CoordSystems: List of CoordSystems
    :param grid_physical_size: Physical size of the grid.
    :param Nxx_dict: Dictionary containing number of grid points.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions.
    """
    for CoordSystem in list_of_CoordSystems:
        register_CFunction_numerical_grid_params_Nxx_dxx_xx(
            CoordSystem=CoordSystem,
            grid_physical_size=grid_physical_size,
            Nxx_dict=Nxx_dict,
            Nchare_dict=Nchare_dict,
        )
    register_CFunction_numerical_grids(
        enable_rfm_precompute=enable_rfm_precompute,
        enable_CurviBCs=enable_CurviBCs,
    )
