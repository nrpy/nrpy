"""
Register numerical_grids_chare() C function, as well as functions called by this one.

These functions set up local chare numerical grid using the global grid for use within the superB infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from typing import List

import nrpy.c_function as cfc
import nrpy.reference_metric as refmetric


def register_CFunction_numerical_grid_params_Nxx_dxx_xx_chare(
    CoordSystem: str,
) -> None:
    """
    Register a C function to Set up a cell-centered grid of size grid_physical_size.
       Set params: Nxx, Nxx_plus_2NGHOSTS, dxx, invdxx, and xx.`

    :param CoordSystem: The coordinate system used for the simulation.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"Set up a cell-centered {CoordSystem} grid of size grid_physical_size. Set params: Nxx, Nxx_plus_2NGHOSTS, dxx, invdxx, and xx."
    cfunc_type = "void"
    name = "numerical_grid_params_Nxx_dxx_xx_chare"
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, params_struct *restrict params_chare, REAL *restrict xx[3], const int chare_index[3]"
    body = r"""
const int Nchare0 = commondata->Nchare0;
const int Nchare1 = commondata->Nchare1;
const int Nchare2 = commondata->Nchare2;
if (params->Nxx0 % Nchare0 != 0) {
  fprintf(stderr, "Error: Division does not result in an integer value: Nxx0 %% Nchare0 = %d %% %d = %d\n", params->Nxx0, Nchare0, params->Nxx0 % Nchare0);
  exit(1);
}
if (params->Nxx1 % Nchare1 != 0) {
  fprintf(stderr, "Error: Division does not result in an integer value: Nxx1 %% Nchare1 = %d %% %d = %d\n", params->Nxx1, Nchare1, params->Nxx1 % Nchare1);
  exit(1);
}
if (params->Nxx2 % Nchare2 != 0) {
  fprintf(stderr, "Error: Division does not result in an integer value: Nxx2 %% Nchare2 = %d %% %d = %d\n", params->Nxx2, Nchare2, params->Nxx2 % Nchare2);
  exit(1);
}
if (Nchare0 > 1 && params->Nxx0 / Nchare0 < NGHOSTS) {
  fprintf(stderr, "Error: params->Nxx0 / Nchare0 is less than NGHOSTS: %d / %d = %d < %d\n", params->Nxx0, Nchare0, params->Nxx0 / Nchare0, NGHOSTS);
  exit(1);
}
if (Nchare1 > 1 && params->Nxx1 / Nchare1 < NGHOSTS) {
  fprintf(stderr, "Error: params->Nxx1 / Nchare1 is less than NGHOSTS: %d / %d = %d < %d\n", params->Nxx1, Nchare1, params->Nxx1 / Nchare1, NGHOSTS);
  exit(1);
}
if (Nchare2 > 1 && params->Nxx2 / Nchare2 < NGHOSTS) {
  fprintf(stderr, "Error: params->Nxx2 / Nchare2 is less than NGHOSTS: %d / %d = %d < %d\n", params->Nxx2, Nchare2, params->Nxx2 / Nchare2, NGHOSTS);
  exit(1);
}
"""
    for dirn in range(3):
        body += f"params_chare->Nxx{dirn} = params->Nxx{dirn}/Nchare{dirn};\n"
    body += rf"""
const REAL grid_physical_size = params_chare->grid_physical_size;
snprintf(params_chare->CoordSystemName, 50, "{CoordSystem}");

params_chare->Nxx_plus_2NGHOSTS0 = params_chare->Nxx0 + 2*NGHOSTS;
params_chare->Nxx_plus_2NGHOSTS1 = params_chare->Nxx1 + 2*NGHOSTS;
params_chare->Nxx_plus_2NGHOSTS2 = params_chare->Nxx2 + 2*NGHOSTS;

// Set grid size to grid_physical_size (set above, based on params->grid_physical_size):
"""
    rfm = refmetric.reference_metric[CoordSystem]
    for key, value in rfm.grid_physical_size_dict.items():
        body += f"params_chare->{key} = {value};\n"

    body += "\n// Set xxmin, xxmax\n"
    for minmax in ["min", "max"]:
        for dirn in range(3):
            param_dir = f"params_chare->xx{minmax}{dirn}"
            if minmax == "min":
                body += f"{param_dir} = params->xx{minmax}{dirn} + (params->dxx{dirn} * (REAL)(params_chare->Nxx{dirn} * chare_index[{dirn}]));\n"
            else:
                body += f"{param_dir} = params->xx{minmax}{dirn} - (params->dxx{dirn} * (REAL)(params_chare->Nxx{dirn} * (Nchare{dirn} - 1 - chare_index[{dirn}])));\n"

    body += r"""

params_chare->dxx0 = params->dxx0;
params_chare->dxx1 = params->dxx1;
params_chare->dxx2 = params->dxx2;

params_chare->invdxx0 = params->invdxx0;
params_chare->invdxx1 = params->invdxx1;
params_chare->invdxx2 = params->invdxx2;

// Set up cell-centered Cartesian coordinate grid, centered at the origin.
xx[0] = (REAL *restrict)malloc(sizeof(REAL)*params_chare->Nxx_plus_2NGHOSTS0);
xx[1] = (REAL *restrict)malloc(sizeof(REAL)*params_chare->Nxx_plus_2NGHOSTS1);
xx[2] = (REAL *restrict)malloc(sizeof(REAL)*params_chare->Nxx_plus_2NGHOSTS2);
for (int j = 0; j < params_chare->Nxx_plus_2NGHOSTS0; j++)
  xx[0][j] = params->xxmin0 + ((REAL)(j - NGHOSTS + (params_chare->Nxx0 * chare_index[0])) + (1.0 / 2.0)) * params_chare->dxx0;
for (int j = 0; j < params_chare->Nxx_plus_2NGHOSTS1; j++)
  xx[1][j] = params->xxmin1 + ((REAL)(j - NGHOSTS + (params_chare->Nxx1 * chare_index[1])) + (1.0 / 2.0)) * params_chare->dxx1;
for (int j = 0; j < params_chare->Nxx_plus_2NGHOSTS2; j++)
  xx[2][j] = params->xxmin2 + ((REAL)(j - NGHOSTS + (params_chare->Nxx2 * chare_index[2])) + (1.0 / 2.0)) * params_chare->dxx2;
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


def register_CFunction_numerical_grids_chare(
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
    enable_psi4_diagnostics: bool = False,
) -> None:
    """
    Register a C function to set up all numerical grids and timestep.

    The function configures the numerical grids based on given parameters, specifically
    focusing on the usage of reference metric precomputations and curvilinear boundary
    conditions.

    :param enable_rfm_precompute: Whether to enable reference metric precomputation (default: False).
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions (default: False).
    :param enable_psi4_diagnostics: Whether or not to enable psi4 diagnostics.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Set up a cell-centered grids of size grid_physical_size."
    cfunc_type = "void"
    name = "numerical_grids_chare"
    params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata, griddata_struct *restrict griddata_chare, const int chare_index[3]"
    body = r"""
        int grid=0;
        griddata_chare[grid].params.CoordSystem_hash = griddata[grid].params.CoordSystem_hash;
        griddata_chare[grid].params.grid_physical_size = griddata[grid].params.grid_physical_size;
    """
    body += r"""
    // Step 1.b: Set Nxx & Nxx_plus_2NGHOSTS, as well as dxx, invdxx, & xx based on grid_physical_size
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  numerical_grid_params_Nxx_dxx_xx_chare(commondata, &griddata[grid].params, &griddata_chare[grid].params, griddata_chare[grid].xx, chare_index);
}

// Step 1.c: Allocate memory for and define reference-metric precomputation lookup tables
"""
    if enable_rfm_precompute:
        body += r"""
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  griddata_chare[grid].rfmstruct = (rfm_struct *)malloc(sizeof(rfm_struct));
  rfm_precompute_malloc(commondata, &griddata_chare[grid].params, griddata_chare[grid].rfmstruct);
  rfm_precompute_defines(commondata, &griddata_chare[grid].params, griddata_chare[grid].rfmstruct, griddata_chare[grid].xx);
}
"""
    else:
        body += "// (reference-metric precomputation disabled)\n"
    body += r"""
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  charecommstruct_set_up(commondata, &griddata[grid].params, &griddata_chare[grid].params, &griddata_chare[grid].charecommstruct, chare_index);
"""
    body += "\n// Step 1.d: Set up curvilinear boundary condition struct (bcstruct)\n"
    if enable_CurviBCs:
        body += r"""
  bcstruct_chare_set_up(commondata, &griddata[grid].params, &griddata_chare[grid].params, &griddata_chare[grid].charecommstruct, griddata_chare[grid].xx, &griddata[grid].bcstruct, &griddata_chare[grid].bcstruct, &griddata_chare[grid].nonlocalinnerbcstruct, chare_index);
"""
    else:
        body += "// (curvilinear boundary conditions bcstruct disabled)\n"

    body += r"""
  // Initialize the diagnostics struct with zero
  griddata_chare[grid].diagnosticstruct = (diagnostic_struct){0};
  // 1D diagnostics set up
  diagnosticstruct_set_up_nearest_1d_y_axis(commondata, &griddata[grid].params, &griddata_chare[grid].params, &griddata_chare[grid].charecommstruct, griddata[grid].xx, chare_index, &griddata_chare[grid].diagnosticstruct);
  diagnosticstruct_set_up_nearest_1d_z_axis(commondata, &griddata[grid].params, &griddata_chare[grid].params, &griddata_chare[grid].charecommstruct, griddata[grid].xx, chare_index, &griddata_chare[grid].diagnosticstruct);
  // 2D diagnostics set up
  diagnosticstruct_set_up_nearest_2d_xy_plane(commondata, &griddata[grid].params, &griddata_chare[grid].params, &griddata_chare[grid].charecommstruct, griddata[grid].xx, chare_index, &griddata_chare[grid].diagnosticstruct);
  diagnosticstruct_set_up_nearest_2d_yz_plane(commondata, &griddata[grid].params, &griddata_chare[grid].params, &griddata_chare[grid].charecommstruct, griddata[grid].xx, chare_index, &griddata_chare[grid].diagnosticstruct);
    """
    if enable_psi4_diagnostics:
        body += r"""
  psi4_diagnostics_set_up(commondata, &griddata[grid].params, &griddata_chare[grid].params,
                                                &griddata_chare[grid].charecommstruct, griddata[grid].xx, chare_index,
                                                &griddata_chare[grid].diagnosticstruct);
"""
    body += r"""
}"""

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
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
    enable_psi4_diagnostics: bool = False,
) -> None:
    """
    Register C functions related to coordinate systems and grid parameters.

    :param list_of_CoordSystems: List of CoordSystems
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions.
    :param enable_psi4_diagnostics: Whether or not to enable psi4 diagnostics.
    """
    for CoordSystem in list_of_CoordSystems:
        register_CFunction_numerical_grid_params_Nxx_dxx_xx_chare(
            CoordSystem=CoordSystem,
        )
    register_CFunction_numerical_grids_chare(
        enable_rfm_precompute=enable_rfm_precompute,
        enable_CurviBCs=enable_CurviBCs,
        enable_psi4_diagnostics=enable_psi4_diagnostics,
    )
