"""
Register numerical_grids_and_timestep() C function, as well as functions called by this one.

These functions set up numerical grids for use within the BHaH infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Dict, List

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.params as par
import nrpy.reference_metric as refmetric

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
_ = par.CodeParameter("int", __name__, "grid_idx", commondata=False, add_to_parfile=False)
# fmt: on


def register_CFunction_numerical_grid_params_Nxx_dxx_xx(
    CoordSystem: str, Nxx_dict: Dict[str, List[int]]
) -> None:
    """
    Register a C function to set up a cell-centered grid of size grid_physical_size.
       Set params: Nxx, Nxx_plus_2NGHOSTS, dxx, invdxx, and xx.

    :param CoordSystem: The coordinate system used for the simulation.
    :param Nxx_dict: A dictionary that maps coordinate systems to lists containing the number of grid points along each direction.

    :raises ValueError: If CoordSystem is not in Nxx_dict.
    """
    if CoordSystem not in Nxx_dict:
        raise ValueError(
            f"{CoordSystem} is not in Nxx_dict = {Nxx_dict}. Please add it."
        )
    for dirn in range(3):
        par.adjust_CodeParam_default(f"Nxx{dirn}", Nxx_dict[CoordSystem][dirn])
    par.adjust_CodeParam_default("CoordSystemName", CoordSystem)

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"""Initializes a cell-centered grid in {CoordSystem} coordinates based on physical dimensions (grid_physical_size).

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
    cfunc_type = "void"
    name = "numerical_grid_params_Nxx_dxx_xx"
    params = "const commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3], const int Nx[3], const bool set_xxmin_xxmax_to_defaults"
    body = "// Set default values for the grid resolution in each dimension.\n"
    for dirn in range(3):
        body += f"params->Nxx{dirn} = {Nxx_dict[CoordSystem][dirn]};\n"
    body += """
// If all components of Nx[] are set to valid values (i.e., not -1), override the default values with Nx[].
if( Nx[0]!=-1 && Nx[1]!=-1 && Nx[2]!=-1 ) {
"""
    for dirn in range(3):
        body += f"params->Nxx{dirn} = Nx[{dirn}];\n"
    body += f"""}}
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
    rfm = refmetric.reference_metric[CoordSystem]

    # Set grid_physical_size & grid_hole_radius
    body += """{
#include "../set_CodeParameters.h"
// Set grid size to a function of grid_physical_size (grid_physical_size set in set_CodeParameters.h above):
"""
    for key, value in rfm.grid_physical_size_dict.items():
        body += f"params->{key} = {value};\n"
    body += "}\n"

    # Set grid_hole_radius
    if "Holey" in CoordSystem or "Wedge" in CoordSystem:
        body += """{
#include "../set_CodeParameters.h"
// Set grid hole radius to a function of grid_hole_radius (grid_hole_radius set in set_CodeParameters.h above):
"""
        for key, value in rfm.grid_hole_radius_dict.items():
            body += f"params->{key} = {value};\n"
        body += "}\n"

    # Set minimum and maximum values of xx[][] for each grid.
    body += """if (set_xxmin_xxmax_to_defaults) {
#include "../set_CodeParameters.h"
// Set {xxmin[], xxmax[]} to default values, which could be functions of other rfm params (set in set_CodeParameters.h above):
"""
    for minmax in ["min", "max"]:
        for dirn in range(3):
            rfm_value = rfm.xxmin[dirn] if minmax == "min" else rfm.xxmax[dirn]
            body += f"params->xx{minmax}{dirn} = {rfm_value};\n"
    body += "}\n"

    # Set quantities that depend on Nxx and {xxmin, xxmax}, then set up coordinate arrays xx[3][Nxxi].
    body += """
// Set quantities that depend on Nxx and {xxmin, xxmax}: dxx, invdxx.
params->dxx0 = (params->xxmax0 - params->xxmin0) / ((REAL)params->Nxx0);
params->dxx1 = (params->xxmax1 - params->xxmin1) / ((REAL)params->Nxx1);
params->dxx2 = (params->xxmax2 - params->xxmin2) / ((REAL)params->Nxx2);

params->invdxx0 = ((REAL)params->Nxx0) / (params->xxmax0 - params->xxmin0);
params->invdxx1 = ((REAL)params->Nxx1) / (params->xxmax1 - params->xxmin1);
params->invdxx2 = ((REAL)params->Nxx2) / (params->xxmax2 - params->xxmin2);

// Set up uniform, cell-centered, topologically Cartesian numerical grid,
//   centered at (xxmin[i] + xxmax[i])/2 in direction i, and store
//   {xx[0], xx[1], xx[2]} arrays.
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
    CoordSystem: str,
) -> None:
    """
    Register a C function to find the CFL-limited timestep dt on a numerical grid.

    The timestep is determined by the relation dt = CFL_FACTOR * ds_min, where ds_min
    is the minimum spacing between neighboring gridpoints on a numerical grid.

    :param CoordSystem: The coordinate system used for the simulation.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"Compute minimum timestep dt = CFL_FACTOR * ds_min on a {CoordSystem} numerical grid."
    cfunc_type = "void"
    name = "cfl_limited_timestep"
    params = "commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3]"
    body = r"""
REAL ds_min = 1e38;
#pragma omp parallel for reduction(min:ds_min)
LOOP_NOOMP(i0, 0, Nxx_plus_2NGHOSTS0,
           i1, 0, Nxx_plus_2NGHOSTS1,
           i2, 0, Nxx_plus_2NGHOSTS2) {
    MAYBE_UNUSED const REAL xx0 = xx[0][i0];
    MAYBE_UNUSED const REAL xx1 = xx[1][i1];
    MAYBE_UNUSED const REAL xx2 = xx[2][i2];
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
    list_of_grid_physical_sizes: List[float],
    gridding_approach: str = "independent grid(s)",
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
    enable_set_cfl_timestep: bool = True,
) -> None:
    """
    Register a C function to set up all numerical grids and timestep.

    The function configures the numerical grids based on given parameters, specifically
    focusing on the usage of reference metric precomputations and curvilinear boundary
    conditions.

    :param list_of_CoordSystems: List of CoordSystems
    :param list_of_grid_physical_sizes: List of grid_physical_size for each CoordSystem; needed for Independent grids.
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :param enable_rfm_precompute: Whether to enable reference metric precomputation (default: False).
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions (default: False).
    :param enable_set_cfl_timestep: Whether to enable computation of dt, the CFL timestep. A custom version can be implemented later.

    :raises ValueError: If invalid gridding_approach selected.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Set up numerical grids and timestep."
    cfunc_type = "void"
    name = "numerical_grids_and_timestep"
    params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata, bool calling_for_first_time"
    body = r"""
  // Step 1.a: Set each CodeParameter in griddata.params to default, for MAXNUMGRIDS grids.
  params_struct_set_to_default(commondata, griddata);"""
    body += rf"""
  // Independent grids
  int Nx[3] = {{ -1, -1, -1 }};

  // Step 1.b: Set commondata->NUMGRIDS to number of CoordSystems we have
  commondata->NUMGRIDS = {len(list_of_CoordSystems)};
"""
    if gridding_approach == "independent grid(s)":
        body += """
  {
    // Step 1.c: For each grid, set Nxx & Nxx_plus_2NGHOSTS, as well as dxx, invdxx, & xx based on grid_physical_size
    const bool set_xxmin_xxmax_to_defaults = true;
    int grid=0;
"""
        for which_CoordSystem, CoordSystem in enumerate(list_of_CoordSystems):
            body += (
                f"  griddata[grid].params.CoordSystem_hash = {CoordSystem.upper()};\n"
            )
            body += f"  griddata[grid].params.grid_physical_size = {list_of_grid_physical_sizes[which_CoordSystem]};\n"
            body += "  numerical_grid_params_Nxx_dxx_xx(commondata, &griddata[grid].params, griddata[grid].xx, Nx, set_xxmin_xxmax_to_defaults);\n"
            body += "  grid++;\n\n"
        body += "}\n"
    elif gridding_approach == "multipatch":
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
        body += """
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
    body += "\n// Step 1.d: Allocate memory for and define reference-metric precomputation lookup tables\n"
    if enable_rfm_precompute:
        body += r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  griddata[grid].rfmstruct = (rfm_struct *)malloc(sizeof(rfm_struct));
  rfm_precompute_malloc(commondata, &griddata[grid].params, griddata[grid].rfmstruct);
  rfm_precompute_defines(commondata, &griddata[grid].params, griddata[grid].rfmstruct, griddata[grid].xx);
}
"""
    else:
        body += "// (reference-metric precomputation disabled)\n"
    body += "\n// Step 1.e: Set up curvilinear boundary condition struct (bcstruct)\n"
    if enable_CurviBCs:
        body += r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  bcstruct_set_up(commondata, &griddata[grid].params, griddata[grid].xx, &griddata[grid].bcstruct);
}
"""
    else:
        body += "// (curvilinear boundary conditions bcstruct disabled)\n"
    if enable_set_cfl_timestep:
        body += r"""
// Step 1.f: Set timestep based on minimum spacing between neighboring gridpoints.
commondata->dt = 1e30;
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  cfl_limited_timestep(commondata, &griddata[grid].params, griddata[grid].xx);
}"""
    body += r"""
// Step 1.g: Initialize timestepping parameters to zero if this is the first time this function is called.
if(calling_for_first_time) {
  commondata->nn = 0;
  commondata->nn_0 = 0;
  commondata->t_0 = 0.0;
  commondata->time = 0.0;
}

// Step 1.h: Set grid_idx for each grid.
for(int grid=0;grid<commondata->NUMGRIDS;grid++) {
   griddata[grid].params.grid_idx = grid;
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
    list_of_grid_physical_sizes: List[float],
    Nxx_dict: Dict[str, List[int]],
    gridding_approach: str = "independent grid(s)",
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
    enable_set_cfl_timestep: bool = True,
) -> None:
    """
    Register C functions related to coordinate systems and grid parameters.

    :param list_of_CoordSystems: List of CoordSystems
    :param list_of_grid_physical_sizes: List of grid_physical_size for each CoordSystem; needed for Independent grids.
    :param Nxx_dict: Dictionary containing number of grid points.
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions.
    :param enable_set_cfl_timestep: Whether to enable computation of dt, the CFL timestep. A custom version can be implemented later.
    """
    for CoordSystem in list_of_CoordSystems:
        register_CFunction_numerical_grid_params_Nxx_dxx_xx(
            CoordSystem=CoordSystem,
            Nxx_dict=Nxx_dict,
        )
        if enable_set_cfl_timestep:
            register_CFunction_cfl_limited_timestep(
                CoordSystem=CoordSystem,
            )
    register_CFunction_numerical_grids_and_timestep(
        list_of_CoordSystems=list_of_CoordSystems,
        list_of_grid_physical_sizes=list_of_grid_physical_sizes,
        gridding_approach=gridding_approach,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_CurviBCs=enable_CurviBCs,
        enable_set_cfl_timestep=enable_set_cfl_timestep,
    )

    if gridding_approach == "multipatch":
        # Register regrid & masking functions
        pass
