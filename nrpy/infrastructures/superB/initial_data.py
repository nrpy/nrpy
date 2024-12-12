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
import nrpy.infrastructures.BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter as admid
from nrpy.equations.general_relativity.InitialData_Cartesian import (
    InitialData_Cartesian,
)
from nrpy.equations.general_relativity.InitialData_Spherical import (
    InitialData_Spherical,
)
from nrpy.infrastructures.BHaH import BHaH_defines_h


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

    :raises ValueError: If `addl_includes` is provided but is not a list, ensuring that additional includes are correctly formatted for inclusion.
    """
    # Step 1: construct this function's contribution to BHaH_defines.h:
    BHd = r"""typedef struct __initial_data_struct__ {
  REAL alpha;

  REAL betaSphorCartU0, betaSphorCartU1, betaSphorCartU2;
  REAL BSphorCartU0, BSphorCartU1, BSphorCartU2;

  REAL gammaSphorCartDD00, gammaSphorCartDD01, gammaSphorCartDD02;
  REAL gammaSphorCartDD11, gammaSphorCartDD12, gammaSphorCartDD22;

  REAL KSphorCartDD00, KSphorCartDD01, KSphorCartDD02;
  REAL KSphorCartDD11, KSphorCartDD12, KSphorCartDD22;
"""
    if enable_T4munu:
        BHd += """
  REAL T4SphorCartUU00,T4SphorCartUU01,T4SphorCartUU02,T4SphorCartUU03;
  REAL                 T4SphorCartUU11,T4SphorCartUU12,T4SphorCartUU13;
  REAL                                 T4SphorCartUU22,T4SphorCartUU23;
  REAL                                                 T4SphorCartUU33;
"""
    BHd += """
} initial_data_struct;
"""
    BHd += "typedef struct __ID_persist_struct__ {\n"
    BHd += ID_persist_struct_str + "\n"
    BHd += "} ID_persist_struct;\n"
    BHaH_defines_h.register_BHaH_defines(__name__, BHd)

    # Step 2: include BHaH_defines.h and register CFunction.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if enable_fd_functions:
        includes += ["finite_difference_functions.h"]

    if addl_includes is not None:
        if not isinstance(addl_includes, list):
            raise ValueError("Error: addl_includes must be a list.")
        includes += addl_includes

    def T4UU_prettyprint() -> str:
        """
        Return a pretty-printed string for T4UU variables in C code.

        :return: A string containing the C declarations for T4UU variables.
        """
        return r"""
  REAL T4UU00,T4UU01,T4UU02,T4UU03;
  REAL        T4UU11,T4UU12,T4UU13;
  REAL               T4UU22,T4UU23;
  REAL                      T4UU33;
"""

    prefunc = """
// ADM variables in the Cartesian basis:
typedef struct __ADM_Cart_basis_struct__ {
  REAL alpha, betaU0,betaU1,betaU2, BU0,BU1,BU2;
  REAL gammaDD00,gammaDD01,gammaDD02,gammaDD11,gammaDD12,gammaDD22;
  REAL KDD00,KDD01,KDD02,KDD11,KDD12,KDD22;
"""
    if enable_T4munu:
        prefunc += T4UU_prettyprint()
    prefunc += "} ADM_Cart_basis_struct;\n"
    ##############
    prefunc += """
// BSSN variables in the Cartesian basis:
typedef struct __BSSN_Cart_basis_struct__ {
  REAL alpha, betaU0,betaU1,betaU2, BU0,BU1,BU2;
  REAL cf, trK;
  REAL gammabarDD00,gammabarDD01,gammabarDD02,gammabarDD11,gammabarDD12,gammabarDD22;
  REAL AbarDD00,AbarDD01,AbarDD02,AbarDD11,AbarDD12,AbarDD22;
"""
    if enable_T4munu:
        prefunc += T4UU_prettyprint()
    prefunc += "} BSSN_Cart_basis_struct;\n"
    ##############
    prefunc += """
// Rescaled BSSN variables in the rfm basis:
typedef struct __rescaled_BSSN_rfm_basis_struct__ {
  REAL alpha, vetU0,vetU1,vetU2, betU0,betU1,betU2;
  REAL cf, trK;
  REAL hDD00,hDD01,hDD02,hDD11,hDD12,hDD22;
  REAL aDD00,aDD01,aDD02,aDD11,aDD12,aDD22;
"""
    if enable_T4munu:
        prefunc += T4UU_prettyprint()
    prefunc += "} rescaled_BSSN_rfm_basis_struct;\n"
    ##############
    ##############
    prefunc += admid.Cfunction_ADM_SphorCart_to_Cart(
        IDCoordSystem=IDCoordSystem,
        enable_T4munu=enable_T4munu,
    )
    prefunc += admid.Cfunction_ADM_Cart_to_BSSN_Cart(enable_T4munu=enable_T4munu)
    prefunc += admid.Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm(
        CoordSystem=CoordSystem,
        enable_T4munu=enable_T4munu,
    )
    prefunc += admid.Cfunction_initial_data_lambdaU_grid_interior(
        CoordSystem=CoordSystem
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
    case INITIALDATA_BIN_ONE: {
    const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

    LOOP_OMP("omp parallel for", i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
      // xxL are the local coordinates on the destination grid
      REAL xxL[3] = { xx[0][i0], xx[1][i1], xx[2][i2] };

      // xCart is the global Cartesian coordinate, which accounts for any grid offsets from the origin.
      REAL xCart[3];
      xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);

      // Read or compute initial data at destination point xCart
      initial_data_struct initial_data;
      ID_function(commondata, params, xCart, ID_persist, &initial_data);

      ADM_Cart_basis_struct ADM_Cart_basis;
      ADM_SphorCart_to_Cart(commondata, params, xCart, &initial_data, &ADM_Cart_basis);

      BSSN_Cart_basis_struct BSSN_Cart_basis;
      ADM_Cart_to_BSSN_Cart(commondata, params, xCart, &ADM_Cart_basis, &BSSN_Cart_basis);

      rescaled_BSSN_rfm_basis_struct rescaled_BSSN_rfm_basis;
      BSSN_Cart_to_rescaled_BSSN_rfm(commondata, params, xxL, &BSSN_Cart_basis, &rescaled_BSSN_rfm_basis);

      const int idx3 = IDX3(i0, i1, i2);
"""
    gf_list = ["alpha", "trK", "cf"]
    for i in range(3):
        gf_list += [f"vetU{i}", f"betU{i}"]
        for j in range(i, 3):
            gf_list += [f"hDD{i}{j}", f"aDD{i}{j}"]
    for gf in sorted(gf_list):
        body += f"gridfuncs->y_n_gfs[IDX4pt({gf.upper()}GF, idx3)] = rescaled_BSSN_rfm_basis.{gf};\n"
    if enable_T4munu:
        for mu in range(4):
            for nu in range(mu, 4):
                gf = f"T4UU{mu}{nu}"
                body += f"gridfuncs->auxevol_gfs[IDX4pt({gf.upper()}GF, idx3)] = rescaled_BSSN_rfm_basis.{gf};\n"
    body += """
        // Initialize lambdaU to zero
        gridfuncs->y_n_gfs[IDX4pt(LAMBDAU0GF, idx3)] = 0.0;
        gridfuncs->y_n_gfs[IDX4pt(LAMBDAU1GF, idx3)] = 0.0;
        gridfuncs->y_n_gfs[IDX4pt(LAMBDAU2GF, idx3)] = 0.0;
      } // END LOOP over all gridpoints on given grid
      break;
    }
"""
    body += """
    case INITIALDATA_BIN_TWO: {
      initial_data_lambdaU_grid_interior(commondata, params, xx, gridfuncs->y_n_gfs);
      break;
    }
  }
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

        admid.register_CFunction_exact_ADM_ID_function(
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
    body += f"""initial_data_reader__convert_ADM_{IDCoordSystem}_to_BSSN(commondata, params,
griddata[grid].xx, &griddata[grid].bcstruct, &griddata[grid].gridfuncs, &ID_persist, {IDtype}, initial_data_part);
    }}"""
    if free_ID_persist_struct_str:
        body += free_ID_persist_struct_str
    body += rf"""
    break;
  }}
  case INITIALDATA_APPLYBCS_INNERONLY: {{
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {{
      // Unpack griddata struct:
      params_struct *restrict params = &griddata[grid].params;
      apply_bcs_inner_only(commondata, params, &griddata[grid].bcstruct, griddata[grid].gridfuncs.y_n_gfs);
    }}
    break;
  }}
  case INITIALDATA_BIN_TWO: {{
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {{
      // Unpack griddata struct:
      params_struct *restrict params = &griddata[grid].params;
      initial_data_reader__convert_ADM_{IDCoordSystem}_to_BSSN(commondata, params, griddata[grid].xx, &griddata[grid].bcstruct, &griddata[grid].gridfuncs,
                                                         NULL, {IDtype}, initial_data_part);
    }}
    break;
  }}
  case INITIALDATA_APPLYBCS_OUTEREXTRAPANDINNER: {{
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {{
      // Unpack griddata struct:
      params_struct *restrict params = &griddata[grid].params;
      apply_bcs_outerextrap_and_inner(commondata, params, &griddata[grid].bcstruct, griddata[grid].gridfuncs.y_n_gfs);
    }}
    break;
  }}
}}
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
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
