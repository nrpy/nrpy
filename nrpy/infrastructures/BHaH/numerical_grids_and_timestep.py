"""
Register numerical_grids_and_timestep_setup() C function,
  which sets up a single numerical grid for use within
  the BHaH infrastructure.

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
# fmt: on


def register_CFunction_numerical_grids_and_timestep_setup(
    CoordSystem: str, grid_physical_size: float, Nxx_dict: Dict[str, List[int]]
) -> None:
    """
    Registers a C function to set up a numerical grid and timestep.

    :param CoordSystem: The coordinate system used for the simulation.
    :param grid_physical_size: The physical size of the grid.
    :param Nxx_dict: A dictionary that maps coordinate systems to lists containing the number of grid points along each direction.
    """
    if CoordSystem not in Nxx_dict:
        raise ValueError(
            f"{CoordSystem} is not in Nxx_dict = {str(Nxx_dict)}. Please add it."
        )
    for dirn in range(3):
        par.adjust_CodeParam_default(f"Nxx{dirn}", Nxx_dict[CoordSystem][dirn])
    par.adjust_CodeParam_default("grid_physical_size", grid_physical_size)

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"Set up a cell-centered {CoordSystem} grid of size grid_physical_size."
    c_type = "void"
    name = "numerical_grids_and_timestep_setup"
    params = "commondata_struct *restrict commondata, params_struct *restrict params, griddata_struct *restrict griddata"
    body = r"""const REAL grid_physical_size = params->grid_physical_size;
const REAL t_final = commondata->t_final;
// Don't increase resolution across an axis of symmetry:
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
            rfm_value = rfm.xxmin[dirn] if minmax == "min" else rfm.xxmax[dirn]
            str_rfm_value = str(rfm_value)
            param_dir = f"params->xx{minmax}{dirn}"

            if str_rfm_value in par.glb_code_params_dict:
                c_type_alias = par.glb_code_params_dict[str_rfm_value].c_type_alias
                if c_type_alias != "#define":
                    body += f"{param_dir} = params->{rfm_value};\n"
                    continue

            body += f"{param_dir} = {rfm_value};\n"

    body += """
params->dxx0 = (params->xxmax0 - params->xxmin0) / ((REAL)params->Nxx0);
params->dxx1 = (params->xxmax1 - params->xxmin1) / ((REAL)params->Nxx1);
params->dxx2 = (params->xxmax2 - params->xxmin2) / ((REAL)params->Nxx2);

params->invdxx0 = ((REAL)params->Nxx0) / (params->xxmax0 - params->xxmin0);
params->invdxx1 = ((REAL)params->Nxx1) / (params->xxmax1 - params->xxmin1);
params->invdxx2 = ((REAL)params->Nxx2) / (params->xxmax2 - params->xxmin2);

// Time parameters:
commondata->nn = 0;
commondata->nn_0 = 0;
commondata->t_0 = 0.0;
commondata->time = 0.0;

// Set up cell-centered Cartesian coordinate grid, centered at the origin.
griddata->xx[0] = (REAL *restrict)malloc(sizeof(REAL)*params->Nxx_plus_2NGHOSTS0);
griddata->xx[1] = (REAL *restrict)malloc(sizeof(REAL)*params->Nxx_plus_2NGHOSTS1);
griddata->xx[2] = (REAL *restrict)malloc(sizeof(REAL)*params->Nxx_plus_2NGHOSTS2);
for (int j = 0; j < params->Nxx_plus_2NGHOSTS0; j++) griddata->xx[0][j] = params->xxmin0 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * params->dxx0;
for (int j = 0; j < params->Nxx_plus_2NGHOSTS1; j++) griddata->xx[1][j] = params->xxmin1 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * params->dxx1;
for (int j = 0; j < params->Nxx_plus_2NGHOSTS2; j++) griddata->xx[2][j] = params->xxmin2 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * params->dxx2;

REAL ds_min = 1e38;
LOOP_NOOMP(i0, NGHOSTS, params->Nxx_plus_2NGHOSTS0-NGHOSTS,
           i1, NGHOSTS, params->Nxx_plus_2NGHOSTS1-NGHOSTS,
           i2, NGHOSTS, params->Nxx_plus_2NGHOSTS2-NGHOSTS) {
#include "set_CodeParameters.h"
    const REAL xx0 = griddata->xx[0][i0];
    const REAL xx1 = griddata->xx[1][i1];
    const REAL xx2 = griddata->xx[2][i2];
    REAL dsmin0, dsmin1, dsmin2;
"""
    dxx0, dxx1, dxx2 = sp.symbols("dxx0 dxx1 dxx2", real=True)
    body += ccg.c_codegen(
        [
            rfm.scalefactor_orthog[0] * dxx0,
            rfm.scalefactor_orthog[1] * dxx1,
            rfm.scalefactor_orthog[2] * dxx2,
        ],
        ["dsmin0", "dsmin1", "dsmin2"],
        include_braces=False,
    )
    body += """ds_min = MIN(ds_min, MIN(dsmin0, MIN(dsmin1, dsmin2)));
}
commondata->dt = commondata->CFL_FACTOR * ds_min;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,  # keep this False or regret having to debug the mess.
        body=body,
    )
