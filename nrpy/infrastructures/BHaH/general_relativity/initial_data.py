"""
Generate C functions for interfacing BSSN applications with initial data importers.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.equations.general_relativity.InitialData_Cartesian import (
    InitialData_Cartesian,
)
from nrpy.equations.general_relativity.InitialData_Spherical import (
    InitialData_Spherical,
)
from nrpy.infrastructures import BHaH


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
    parallelization = par.parval_from_str("parallelization")

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

    BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter.register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN(
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
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata_host, griddata_struct *restrict griddata"
        if parallelization in ["cuda"]
        else "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    body = ""
    host_griddata = "griddata_host" if parallelization in ["cuda"] else "griddata"
    if enable_checkpointing:
        body += """// Attempt to read checkpoint file. If it doesn't exist, then continue. Otherwise return.
if( read_checkpoint(commondata, griddata) ) return;
""".replace(
            "griddata",
            f"{host_griddata}, griddata" if parallelization in ["cuda"] else "griddata",
        )
    body += "ID_persist_struct ID_persist;\n"
    if populate_ID_persist_struct_str:
        body += populate_ID_persist_struct_str
    body += """
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
"""
    body += (
        f"initial_data_reader__convert_ADM_{IDCoordSystem}_to_BSSN(commondata, params,"
        f"(const REAL *restrict *) {host_griddata}[grid].xx, (const REAL *restrict *) griddata[grid].xx,"
        f"&griddata[grid].bcstruct, &{host_griddata}[grid].gridfuncs, &griddata[grid].gridfuncs, &ID_persist, {IDtype});"
        if parallelization in ["cuda"]
        else f"initial_data_reader__convert_ADM_{IDCoordSystem}_to_BSSN(commondata, params,"
        f"(const REAL *restrict *) {host_griddata}[grid].xx, &griddata[grid].bcstruct, &griddata[grid].gridfuncs, &ID_persist, {IDtype});"
    )
    if CoordSystem.startswith("GeneralRFM"):
        body += f"\n  generalrfm_precompute__{CoordSystem}(commondata, params, (const REAL *restrict *) {host_griddata}[grid].xx, {host_griddata}[grid].gridfuncs.auxevol_gfs);"
    body += """
  apply_bcs_outerextrap_and_inner(commondata, params, &griddata[grid].bcstruct, griddata[grid].gridfuncs.y_n_gfs);
}
"""
    if free_ID_persist_struct_str:
        body += free_ID_persist_struct_str

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
