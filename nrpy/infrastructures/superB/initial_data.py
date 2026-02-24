"""
Library of C functions for solving the BSSN equations in curvilinear coordinates, using a reference-metric formalism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Optional, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
from nrpy.equations.general_relativity.InitialData_Cartesian import (
    InitialData_Cartesian,
)
from nrpy.equations.general_relativity.InitialData_Spherical import (
    InitialData_Spherical,
)
from nrpy.infrastructures import BHaH


def register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN(
    CoordSystem: str,
    addl_includes: Optional[List[str]] = None,
    IDCoordSystem: str = "Spherical",
    enable_T4munu: bool = False,
    enable_fd_functions: bool = False,
    ID_persist_struct_str: str = "",
) -> None:
    """
    Register the CFunction for converting initial ADM data to BSSN variables.

    :param CoordSystem: Coordinate system for output BSSN variables.
    :param addl_includes: Additional header files to include.
    :param IDCoordSystem: Coordinate system for input ADM variables. Defaults to "Spherical".
    :param enable_T4munu: Whether to include stress-energy tensor components.
    :param enable_fd_functions: Whether to enable finite-difference functions.
    :param ID_persist_struct_str: String for persistent ID structure.
    """
    includes, prefunc, lambdaU_launch = (
        BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter.setup_ADM_initial_data_reader(
            ID_persist_struct_str=ID_persist_struct_str,
            enable_T4munu=enable_T4munu,
            enable_fd_functions=enable_fd_functions,
            addl_includes=addl_includes,
            CoordSystem=CoordSystem,
            IDCoordSystem=IDCoordSystem,
        )
    )

    desc = f"Read ADM data in the {IDCoordSystem} basis, and output rescaled BSSN data in the {CoordSystem} basis"
    cfunc_type = "void"
    name = f"initial_data_reader__convert_ADM_{IDCoordSystem}_to_BSSN"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    REAL *restrict xx[3], bc_struct *restrict bcstruct, MoL_gridfunctions_struct *restrict gridfuncs,
    ID_persist_struct *restrict ID_persist,
    void ID_function(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
                     const ID_persist_struct *restrict ID_persist,
                     initial_data_struct *restrict initial_data), const int initial_data_part"""

    body = r"""
  switch (initial_data_part) {
    case INITIALDATA_BIN_ONE: {"""

    body += BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter.build_initial_data_conversion_loop(
        enable_T4munu
    )

    body += (
        BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter.build_lambdaU_zeroing_block()
    )

    body += """
      break;
    }
"""
    body += f"""
    case INITIALDATA_BIN_TWO: {{
      {lambdaU_launch}
      break;
    }}
  }}
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


def register_CFunction_initial_data(
    CoordSystem: str,
    IDtype: str,
    IDCoordSystem: str,
    ID_persist_struct_str: str,
    enable_checkpointing: bool = False,
    populate_ID_persist_struct_str: str = "",
    free_ID_persist_struct_str: str = "",
    enable_T4munu: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C functions for converting ADM initial data to BSSN variables and applying boundary conditions.

    The function performs the following operations:
    1. Registers the exact ADM initial data function.
    2. Registers a function for converting ADM initial data to BSSN variables in the specified coordinate system.
    3. Generates C code for setting initial data and applying boundary conditions.

    :param CoordSystem: The coordinate system for the calculation.
    :param IDtype: The type of initial data.
    :param IDCoordSystem: The native coordinate system of the initial data.
    :param enable_checkpointing: Attempt to read from a checkpoint file before generating initial data.
    :param ID_persist_struct_str: A string representing the persistent structure for the initial data.
    :param populate_ID_persist_struct_str: Optional string to populate the persistent structure for initial data.
    :param free_ID_persist_struct_str: Optional string to free the persistent structure for initial data.
    :param enable_T4munu: Whether to include the stress-energy tensor. Defaults to False.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    try:
        ID: Union[InitialData_Cartesian, InitialData_Spherical]
        if IDCoordSystem == "Cartesian":
            ID = InitialData_Cartesian(IDtype=IDtype)
        else:
            ID = InitialData_Spherical(IDtype=IDtype)

        BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter.register_CFunction_exact_ADM_ID_function(
            IDCoordSystem,
            IDtype,
            ID.alpha,
            ID.betaU,
            ID.BU,
            ID.gammaDD,
            ID.KDD,
        )
    except (ValueError, RuntimeError):
        print(
            f"Warning: {IDtype} does not correspond to an implemented exact initial data type."
        )
        print("Assuming initial data functionality is implemented elsewhere.")

    register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN(
        CoordSystem,
        IDCoordSystem=IDCoordSystem,
        ID_persist_struct_str=ID_persist_struct_str,
        enable_T4munu=enable_T4munu,
    )
    if CoordSystem.startswith("GeneralRFM"):
        BHaH.generalrfm_precompute.register_CFunctions_generalrfm_support(CoordSystem)

    desc = "Set initial data."
    cfunc_type = "void"
    name = "initial_data"
    params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata, const int initial_data_part"

    body = ""
    if enable_checkpointing:
        body += """// Attempt to read checkpoint file. If it doesn't exist, then continue. Otherwise return.
if( read_checkpoint(commondata, griddata) ) return;
"""

    body += """
switch (initial_data_part) {
  case INITIALDATA_BIN_ONE: {"""
    body += "ID_persist_struct ID_persist;\n"

    if populate_ID_persist_struct_str:
        body += populate_ID_persist_struct_str
    body += """
    for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
      // Unpack griddata struct:
      params_struct *restrict params = &griddata[grid].params;
"""
    if CoordSystem.startswith("GeneralRFM"):
        body += """      generalrfm_precompute(commondata, params, (const REAL *restrict *)griddata[grid].xx, griddata[grid].gridfuncs.auxevol_gfs);
"""
    body += f"""initial_data_reader__convert_ADM_{IDCoordSystem}_to_BSSN(commondata, params,
griddata[grid].xx, &griddata[grid].bcstruct, &griddata[grid].gridfuncs, &ID_persist, {IDtype}, initial_data_part);
    }}"""
    if free_ID_persist_struct_str:
        body += free_ID_persist_struct_str
    body += """
    break;
  }"""

    apply_inner_bcs_block = (
        BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter.build_apply_inner_bcs_block()
    )
    body += f"""
  case INITIALDATA_APPLYBCS_INNERONLY: {{
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {{

      // Unpack griddata struct:
      params_struct *restrict params = &griddata[grid].params;
      bc_struct *restrict bcstruct = &griddata[grid].bcstruct;
      MoL_gridfunctions_struct *restrict gridfuncs = &griddata[grid].gridfuncs;

      {apply_inner_bcs_block}
    }}
    break;
  }}"""
    body += f"""
  case INITIALDATA_BIN_TWO: {{
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {{
      // Unpack griddata struct:
      params_struct *restrict params = &griddata[grid].params;
"""
    body += f"""      initial_data_reader__convert_ADM_{IDCoordSystem}_to_BSSN(commondata, params, griddata[grid].xx, &griddata[grid].bcstruct, &griddata[grid].gridfuncs,
                                                         NULL, {IDtype}, initial_data_part);
    }}
    break;
  }}"""
    body += """
  case INITIALDATA_APPLYBCS_OUTEREXTRAPANDINNER: {
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
      // Unpack griddata struct:
      params_struct *restrict params = &griddata[grid].params;
      apply_bcs_outerextrap_and_inner(commondata, params, &griddata[grid].bcstruct, griddata[grid].gridfuncs.y_n_gfs);
    }
    break;
  }
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
    return pcg.NRPyEnv()
