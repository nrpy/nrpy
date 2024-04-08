"""
Register numerical_grids_and_timestep() C function, as well as functions called by this one.

These functions set up numerical grids for use within the BHaH infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Dict, List
import sympy as sp
import nrpy.c_function as cfc
import nrpy.params as par
import nrpy.reference_metric as refmetric
import nrpy.c_codegen as ccg

# fmt: off
for i in range(3):
    _ = par.CodeParameter("int", __name__, f"Nxx_plus_2NGHOSTS{i}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("int", __name__, f"Nxx{i}", 64)
    # reference_metric sets xxmin and xxmax below.
    _ = par.CodeParameter("REAL", __name__, f"xxmin{i}", -10.0, add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("REAL", __name__, f"xxmax{i}", 10.0, add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("REAL", __name__, f"invdxx{i}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("REAL", __name__, f"dxx{i}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
_ = par.CodeParameter("REAL", __name__, "convergence_factor", 1.0, commondata=True)
_ = par.CodeParameter("int", __name__, "CoordSystem_hash", commondata=False, add_to_parfile=False)
_ = par.CodeParameter("char[200]", __name__, "gridding_choice", "independent grid(s)", commondata=True, add_to_parfile=True)
# fmt: on


def register_CFunction_numerical_grid_params_Nxx_dxx_xx(
    CoordSystem: str, grid_physical_size: float, Nxx_dict: Dict[str, List[int]]
) -> None:
    """
    Register a C function to Set up a cell-centered grid of size grid_physical_size.
       Set params: Nxx, Nxx_plus_2NGHOSTS, dxx, invdxx, and xx.

    :param CoordSystem: The coordinate system used for the simulation.
    :param grid_physical_size: The physical size of the grid.
    :param Nxx_dict: A dictionary that maps coordinate systems to lists containing the number of grid points along each direction.

    :raises ValueError: If CoordSystem is not in Nxx_dict.
    """
    if CoordSystem not in Nxx_dict:
        raise ValueError(
            f"{CoordSystem} is not in Nxx_dict = {Nxx_dict}. Please add it."
        )
    for dirn in range(3):
        par.adjust_CodeParam_default(f"Nxx{dirn}", Nxx_dict[CoordSystem][dirn])
    par.adjust_CodeParam_default("grid_physical_size", grid_physical_size)
    par.adjust_CodeParam_default("CoordSystemName", CoordSystem)

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"""Initializes a cell-centered grid in {CoordSystem} coordinates based on physical dimensions (grid_physical_size).

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
    cfunc_type = "void"
    name = "numerical_grid_params_Nxx_dxx_xx"
    params = "const commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3], const int Nx[3], const bool grid_is_resized"
    body = "// Start by setting default values for Nxx.\n"
    for dirn in range(3):
        body += f"params->Nxx{dirn} = {Nxx_dict[CoordSystem][dirn]};\n"
    body += """
// If all components of Nx[] are set to reasonable values (i.e., not -1), then set params->Nxx{} to Nx[].
if( !(Nx[0]==-1 || Nx[1]==-1 || Nx[2]==-1) ) {
"""
    for dirn in range(3):
        body += f"params->Nxx{dirn} = Nx[{dirn}];\n"
    body += f"""}}
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

// Set grid size to grid_physical_size (set above, based on params->grid_physical_size):
"""
    rfm = refmetric.reference_metric[CoordSystem]

    # Set grid_physical_size & grid_hole_radius
    body += """{
    #include "../set_CodeParameters.h"
    // Set grid size to a function of grid_physical_size, just set in set_CodeParameters.h above:
    """
    for key, value in rfm.grid_physical_size_dict.items():
        body += f"params->{key} = {value};\n"
    body += "}\n"

    # Set xxmin and xxmax
    body += """if( !grid_is_resized ) {
#include "../set_CodeParameters.h"
"""
    body += "// Set xxmin, xxmax\n"
    if "Cartesian" in CoordSystem:
        body += """
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
            rfm_value = rfm.xxmin[dirn] if minmax == "min" else rfm.xxmax[dirn]
            str_rfm_value = str(rfm_value)
            param_dir = f"params->xx{minmax}{dirn}"

            if str_rfm_value in par.glb_code_params_dict:
                c_type_alias = par.glb_code_params_dict[str_rfm_value].cparam_type
                if c_type_alias != "#define":
                    body += f"{param_dir} = params->{rfm_value};\n"
                    continue

            if "Cartesian" in CoordSystem:
                body += f"{param_dir} = {rfm_value} * factor_rescale_minmax_xx{dirn};\n"
            else:
                body += f"{param_dir} = {rfm_value};\n"

    body += """}

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
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,  # keep this False or regret having to debug the mess.
        body=body,
    )


def register_CFunction_cfl_limited_timestep(
    CoordSystem: str, fp_type: str = "double"
) -> None:
    """
    Register a C function to find the CFL-limited timestep dt on a numerical grid.

    The timestep is determined by the relation dt = CFL_FACTOR * ds_min, where ds_min
    is the minimum spacing between neighboring gridpoints on a numerical grid.

    :param CoordSystem: The coordinate system used for the simulation.
    :param fp_type: Floating point type, e.g., "double".
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"Output minimum gridspacing ds_min on a {CoordSystem} numerical grid."
    cfunc_type = "void"
    name = "cfl_limited_timestep"
    params = "commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3], bc_struct *restrict bcstruct"
    body = r"""
REAL ds_min = 1e38;
#pragma omp parallel for reduction(min:ds_min)
LOOP_NOOMP(i0, 0, Nxx_plus_2NGHOSTS0,
           i1, 0, Nxx_plus_2NGHOSTS1,
           i2, 0, Nxx_plus_2NGHOSTS2) {
    const REAL xx0 = xx[0][i0];
    const REAL xx1 = xx[1][i1];
    const REAL xx2 = xx[2][i2];
    REAL dsmin0, dsmin1, dsmin2;
"""
    rfm = refmetric.reference_metric[CoordSystem]
    dxx0, dxx1, dxx2 = sp.symbols("dxx0 dxx1 dxx2", real=True)
    body += ccg.c_codegen(
        [
            sp.Abs(rfm.scalefactor_orthog[0] * dxx0),
            sp.Abs(rfm.scalefactor_orthog[1] * dxx1),
            sp.Abs(rfm.scalefactor_orthog[2] * dxx2),
        ],
        ["dsmin0", "dsmin1", "dsmin2"],
        include_braces=False,
        fp_type=fp_type,
    )
    body += """ds_min = MIN(ds_min, MIN(dsmin0, MIN(dsmin1, dsmin2)));
}
commondata->dt = MIN(commondata->dt, ds_min * commondata->CFL_FACTOR);
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )


def register_CFunction_numerical_grids_and_timestep(
    list_of_CoordSystems: List[str],
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
) -> None:
    """
    Register a C function to set up all numerical grids and timestep.

    The function configures the numerical grids based on given parameters, specifically
    focusing on the usage of reference metric precomputations and curvilinear boundary
    conditions.

    :param list_of_CoordSystems: List of CoordSystems
    :param enable_rfm_precompute: Whether to enable reference metric precomputation (default: False).
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions (default: False).
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Set up a cell-centered grids of size grid_physical_size."
    cfunc_type = "void"
    name = "numerical_grids_and_timestep"
    params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata, bool calling_for_first_time"
    body = r"""
    // Step 1.a: Set each CodeParameter in griddata.params to default, for MAXNUMGRIDS grids.
    params_struct_set_to_default(commondata, griddata);"""
    body += r"""
      if(strncmp(commondata->gridding_choice, "independent grid(s)", 200) == 0) {
        // Independent grids
        bool grid_is_resized=false;
        int Nx[3] = { -1, -1, -1 };


        // Step 1.b: For each grid, set Nxx & Nxx_plus_2NGHOSTS, as well as dxx, invdxx, & xx based on grid_physical_size
        int grid=0;
    """
    for CoordSystem in list_of_CoordSystems:
        body += f"griddata[grid].params.CoordSystem_hash = {CoordSystem.upper()};\n"
        body += "numerical_grid_params_Nxx_dxx_xx(commondata, &griddata[grid].params, griddata[grid].xx, Nx, grid_is_resized);\n"
        body += "grid++;\n\n"
    body += r"""}

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
  bcstruct_set_up(commondata, &griddata[grid].params, griddata[grid].xx, &griddata[grid].bcstruct);
}
"""
    else:
        body += "// (curvilinear boundary conditions bcstruct disabled)\n"
    body += r"""
// Step 1.e: Set timestep based on minimum spacing between neighboring gridpoints.
commondata->dt = 1e30;
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  cfl_limited_timestep(commondata, &griddata[grid].params, griddata[grid].xx, &griddata[grid].bcstruct);
}

// Step 1.f: Initialize timestepping parameters to zero if this is the first time this function is called.
if(calling_for_first_time) {
  commondata->nn = 0;
  commondata->nn_0 = 0;
  commondata->t_0 = 0.0;
  commondata->time = 0.0;
}
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


def register_CFunctions(
    list_of_CoordSystems: List[str],
    grid_physical_size: float,
    Nxx_dict: Dict[str, List[int]],
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
    fp_type: str = "double",
) -> None:
    """
    Register C functions related to coordinate systems and grid parameters.

    :param list_of_CoordSystems: List of CoordSystems
    :param grid_physical_size: Physical size of the grid.
    :param Nxx_dict: Dictionary containing number of grid points.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions.
    :param fp_type: Floating point type, e.g., "double".
    """
    for CoordSystem in list_of_CoordSystems:
        register_CFunction_numerical_grid_params_Nxx_dxx_xx(
            CoordSystem=CoordSystem,
            grid_physical_size=grid_physical_size,
            Nxx_dict=Nxx_dict,
        )
        register_CFunction_cfl_limited_timestep(
            CoordSystem=CoordSystem, fp_type=fp_type
        )
    register_CFunction_numerical_grids_and_timestep(
        list_of_CoordSystems=list_of_CoordSystems,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_CurviBCs=enable_CurviBCs,
    )
