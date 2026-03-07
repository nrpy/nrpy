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

_TWOPUNCTURES_SOLVE_ONCE_AND_BROADCAST_PREFUNC = r"""
static inline void tp_write_derivs(FILE *fp, const derivs *v, const int ntotal) {
#define TP_WRITE(ptr) fwrite((const void *)(ptr), sizeof(REAL), ntotal, fp)
  TP_WRITE(v->d0);
  TP_WRITE(v->d1);
  TP_WRITE(v->d2);
  TP_WRITE(v->d3);
  TP_WRITE(v->d11);
  TP_WRITE(v->d12);
  TP_WRITE(v->d13);
  TP_WRITE(v->d22);
  TP_WRITE(v->d23);
  TP_WRITE(v->d33);
#undef TP_WRITE
}

static inline void tp_read_derivs(FILE *fp, derivs *v, const int ntotal) {
#define TP_READ(ptr) fread((void *)(ptr), sizeof(REAL), ntotal, fp)
  TP_READ(v->d0);
  TP_READ(v->d1);
  TP_READ(v->d2);
  TP_READ(v->d3);
  TP_READ(v->d11);
  TP_READ(v->d12);
  TP_READ(v->d13);
  TP_READ(v->d22);
  TP_READ(v->d23);
  TP_READ(v->d33);
#undef TP_READ
}

static void tp_solve_once_and_broadcast(ID_persist_struct *restrict ID_persist) {
  if (CkNumPes() <= 1) {
    const double t0 = CkWallTimer();
    TP_solve(ID_persist);
    CkPrintf("[startup] TP_solve complete in %.3f s (single-PE)\n", CkWallTimer() - t0);
    return;
  }
  const int ntotal = ID_persist->npoints_A * ID_persist->npoints_B * ID_persist->npoints_phi;
  const char *tp_dump = "twopunctures_idpersist.bin";
  const char *tp_done = "twopunctures_idpersist.done";

  if (CkMyPe() == 0) {
    const double t0 = CkWallTimer();
    (void)remove(tp_done);
    (void)remove(tp_dump);
    ID_persist->verbose = false;
    TP_solve(ID_persist);
    CkPrintf("[startup] TP_solve complete in %.3f s on PE 0\n", CkWallTimer() - t0);
    FILE *fp = fopen(tp_dump, "wb");
    if (fp == NULL)
      CkAbort("ERROR: could not open twopunctures_idpersist.bin for writing.");
    if (fwrite((const void *)ID_persist, sizeof(ID_persist_struct), 1, fp) != 1)
      CkAbort("ERROR: writing ID_persist header failed.");
    tp_write_derivs(fp, &ID_persist->v, ntotal);
    tp_write_derivs(fp, &ID_persist->cf_v, ntotal);
    fclose(fp);
    FILE *done = fopen(tp_done, "wb");
    if (done == NULL)
      CkAbort("ERROR: could not create twopunctures_idpersist.done.");
    fclose(done);
    return;
  }

  extern void allocate_derivs(derivs * v, int n);
  allocate_derivs(&ID_persist->v, ntotal);
  allocate_derivs(&ID_persist->cf_v, ntotal);

  int tries = 0;
  const double t_wait_start = CkWallTimer();
  while (access(tp_done, F_OK) != 0) {
    usleep(100000);
    tries++;
    if (tries > 6000)
      CkAbort("ERROR: timed out waiting for twopunctures_idpersist.done");
  }

  FILE *fp = fopen(tp_dump, "rb");
  if (fp == NULL)
    CkAbort("ERROR: could not open twopunctures_idpersist.bin for reading.");
  const derivs v_local = ID_persist->v;
  const derivs cf_v_local = ID_persist->cf_v;
  ID_persist_struct tmp;
  if (fread((void *)&tmp, sizeof(ID_persist_struct), 1, fp) != 1)
    CkAbort("ERROR: reading ID_persist header failed.");
  *ID_persist = tmp;
  ID_persist->v = v_local;
  ID_persist->cf_v = cf_v_local;
  tp_read_derivs(fp, &ID_persist->v, ntotal);
  tp_read_derivs(fp, &ID_persist->cf_v, ntotal);
  fclose(fp);
  if (CkMyPe() == 1)
    CkPrintf("[startup] non-root PEs loaded TP data after %.3f s wait\n", CkWallTimer() - t_wait_start);
}
"""

_TWOPUNCTURES_SOLVE_ONCE_AND_BROADCAST_POPULATE = r"""
initialize_ID_persist_struct(commondata, &ID_persist);
tp_solve_once_and_broadcast(&ID_persist);
"""


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
    initial_data_addl_includes: Optional[List[str]] = None,
    initial_data_prefunc_str: str = "",
    enable_tp_solve_broadcast: bool = False,
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
    :param initial_data_addl_includes: Optional additional include headers for initial_data().
    :param initial_data_prefunc_str: Optional helper C code emitted before initial_data().
    :param enable_tp_solve_broadcast: If True, inject helper C code and default setup that solves TwoPunctures on PE 0 and broadcasts the result to other PEs.
    :param populate_ID_persist_struct_str: Optional string to populate the persistent structure for initial data.
    :param free_ID_persist_struct_str: Optional string to free the persistent structure for initial data.
    :param enable_T4munu: Whether to include the stress-energy tensor. Defaults to False.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if initial_data_addl_includes:
        includes.extend(initial_data_addl_includes)
    if enable_tp_solve_broadcast and "<unistd.h>" not in includes:
        includes.append("<unistd.h>")

    prefunc_blocks: List[str] = []
    if initial_data_prefunc_str:
        prefunc_blocks.append(initial_data_prefunc_str)
    if enable_tp_solve_broadcast:
        prefunc_blocks.append(_TWOPUNCTURES_SOLVE_ONCE_AND_BROADCAST_PREFUNC)
    effective_prefunc = "\n".join(prefunc_blocks) if prefunc_blocks else None

    effective_populate_ID_persist_struct_str = populate_ID_persist_struct_str
    if enable_tp_solve_broadcast and not effective_populate_ID_persist_struct_str:
        effective_populate_ID_persist_struct_str = (
            _TWOPUNCTURES_SOLVE_ONCE_AND_BROADCAST_POPULATE
        )

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
    body += """
    const bool is_root_pe = (CkMyPe() == 0);
    const double t_bin_one_start = CkWallTimer();
    if (is_root_pe) CkPrintf("[startup][initial_data] INITIALDATA_BIN_ONE begin (TwoPunctures setup + ADM->BSSN part 1)\\n");
"""
    body += "ID_persist_struct ID_persist;\n"

    if effective_populate_ID_persist_struct_str:
        body += """
    const double t_id_setup_start = CkWallTimer();
"""
        body += effective_populate_ID_persist_struct_str
        body += """
    if (is_root_pe) CkPrintf("[startup][initial_data] TwoPunctures ID setup/solve took %.3f s\\n", CkWallTimer() - t_id_setup_start);
"""
    body += """
    const double t_convert_start = CkWallTimer();
"""
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
    body += """
    if (is_root_pe) CkPrintf("[startup][initial_data] ADM->BSSN conversion (INITIALDATA_BIN_ONE) took %.3f s\\n", CkWallTimer() - t_convert_start);
"""
    if free_ID_persist_struct_str:
        body += """
    const double t_id_free_start = CkWallTimer();
"""
        body += free_ID_persist_struct_str
        body += """
    if (is_root_pe) CkPrintf("[startup][initial_data] Free TwoPunctures ID memory took %.3f s\\n", CkWallTimer() - t_id_free_start);
"""
    body += """
    if (is_root_pe) CkPrintf("[startup][initial_data] INITIALDATA_BIN_ONE total %.3f s\\n", CkWallTimer() - t_bin_one_start);
"""
    body += """
    break;
  }"""

    apply_inner_bcs_block = (
        BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter.build_apply_inner_bcs_block()
    )
    body += f"""
  case INITIALDATA_APPLYBCS_INNERONLY: {{
    const bool is_root_pe = (CkMyPe() == 0);
    const double t_inner_bcs_start = CkWallTimer();
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {{

      // Unpack griddata struct:
      params_struct *restrict params = &griddata[grid].params;
      bc_struct *restrict bcstruct = &griddata[grid].bcstruct;
      MoL_gridfunctions_struct *restrict gridfuncs = &griddata[grid].gridfuncs;

      {apply_inner_bcs_block}
    }}
    if (is_root_pe) CkPrintf("[startup][initial_data] INITIALDATA_APPLYBCS_INNERONLY took %.3f s\\n", CkWallTimer() - t_inner_bcs_start);
    break;
  }}"""
    body += """
  case INITIALDATA_BIN_TWO: {
    const bool is_root_pe = (CkMyPe() == 0);
    const double t_bin_two_start = CkWallTimer();
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
      // Unpack griddata struct:
      params_struct *restrict params = &griddata[grid].params;
"""
    body += f"""      initial_data_reader__convert_ADM_{IDCoordSystem}_to_BSSN(commondata, params, griddata[grid].xx, &griddata[grid].bcstruct, &griddata[grid].gridfuncs,
                                                         NULL, {IDtype}, initial_data_part);
    }}
    if (is_root_pe) CkPrintf("[startup][initial_data] ADM->BSSN conversion (INITIALDATA_BIN_TWO / lambdaU pass) took %.3f s\\n", CkWallTimer() - t_bin_two_start);
    break;
  }}"""
    body += """
  case INITIALDATA_APPLYBCS_OUTEREXTRAPANDINNER: {
    const bool is_root_pe = (CkMyPe() == 0);
    const double t_outer_inner_bcs_start = CkWallTimer();
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
      // Unpack griddata struct:
      params_struct *restrict params = &griddata[grid].params;
      apply_bcs_outerextrap_and_inner(commondata, params, &griddata[grid].bcstruct, griddata[grid].gridfuncs.y_n_gfs);
    }
    if (is_root_pe) CkPrintf("[startup][initial_data] INITIALDATA_APPLYBCS_OUTEREXTRAPANDINNER took %.3f s\\n", CkWallTimer() - t_outer_inner_bcs_start);
    break;
  }
}
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=effective_prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
