"""
Generate timestepping.cpp, timestepping.h and timestepping.ci for the superB infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from pathlib import Path
from typing import Dict, List, Tuple, Union

import sympy as sp  # Import SymPy, a computer algebra system written entirely in Python

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.helpers.generic import clang_format
from nrpy.infrastructures.BHaH import griddata_commondata
from nrpy.infrastructures.BHaH.MoLtimestepping.MoL import generate_gridfunction_names
from nrpy.infrastructures.BHaH.MoLtimestepping.RK_Butcher_Table_Dictionary import (
    generate_Butcher_tables,
)
from nrpy.infrastructures.superB.MoL import (
    generate_post_rhs_output_list,
    generate_rhs_output_exprs,
)


def generate_send_nonlocalinnerbc_data_code(which_gf: str) -> str:
    """
    Generate code for sending nonlocal inner bc data.

    :param which_gf: Which grid function
    :return: Code for sending nonlocal inner bc data
    """
    return f"""
    if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares > 0) {{
        serial {{
          send_nonlocalinnerbc_data({which_gf}, grid);
        }}
      }}
    """


def generate_process_nonlocalinnerbc_code() -> str:
    """
    Generate code for nonlocal inner bc processing.

    :return: Code for nonlocal inner bc processing.
    """
    return """
      if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares > 0) {
        for (iter = 0; iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares; iter++) {
          when receiv_nonlocalinnerbc_data(int src_chare_idx3, int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]) {
            serial {
              set_tmpBuffer_innerbc_receiv(src_chare_idx3, len_tmpBuffer, tmpBuffer, grid);
              type_gfs_nonlocal_innerbc = type_gfs;
            }
          }
        }
        serial {
          process_nonlocalinnerbc(type_gfs_nonlocal_innerbc, grid);
        }
      }
"""


def generate_mol_step_forward_code(
    rk_substep: str, rhs_output_exprs_list: List[str], post_rhs_output_list: List[str]
) -> str:
    """
    Generate code for MoL step forward in time.

    :param rk_substep: Runge-Kutta substep.
    :param rhs_output_exprs_list: List of output expression for the RHS.
    :param post_rhs_output_list: List of outputs for post-RHS expressions.
    :return: Code for MoL step forward in time.
    """
    return_str = f"""
    serial{{
        MoL_step_forward_in_time(&commondata, griddata_chare, time_start, {rk_substep},  MOL_PART_1);
    }}
"""
    return_str += """
    if (strncmp(commondata.outer_bc_type, "radiation", 50) == 0) {
"""
    for rhs_output_exprs in rhs_output_exprs_list:
        return_str += generate_send_nonlocalinnerbc_data_code(rhs_output_exprs)
        return_str += generate_process_nonlocalinnerbc_code()

    return_str += """
    }
"""
    return_str += f"""
    serial{{
        MoL_step_forward_in_time(&commondata, griddata_chare, time_start, {rk_substep}, MOL_PART_2);
        MoL_step_forward_in_time(&commondata, griddata_chare, time_start, {rk_substep}, MOL_PART_3_APPLY_BCS);
    }}
"""
    return_str += """
    if (strncmp(commondata.outer_bc_type, "extrapolation", 50) == 0) {
"""
    for post_rhs_output in post_rhs_output_list:
        return_str += generate_send_nonlocalinnerbc_data_code(post_rhs_output)
        return_str += generate_process_nonlocalinnerbc_code()

    return_str += """
    }
"""
    return_str += f"""
    serial{{
        MoL_step_forward_in_time(&commondata, griddata_chare, time_start, {rk_substep}, MOL_PART_3_AFTER_APPLY_BCS);
    }}
"""
    return return_str


def generate_send_neighbor_data_code(which_gf: str, neighbor_direction: str) -> str:
    """
    Generate code for sending neighbor data.

    :param which_gf: Which grid function
    :param neighbor_direction: Direction of neighboring grid
    :return: Code for sending neighbor data
    """
    return f"""
    serial {{
        send_neighbor_data({which_gf}, {neighbor_direction}, grid);
    }}
    """


def generate_process_ghost_code(
    axis: str, pos_ghost_type: str, neg_ghost_type: str, nchare_var: str
) -> str:
    """
    Generate code for ghost zone processing.

    :param axis: Axis index
    :param pos_ghost_type: Positive ghost zone type
    :param neg_ghost_type: Negative ghost zone type
    :param nchare_var: Number of charm grids
    :return: Code for ghost zone processing
    """
    this_index_var = f"thisIndex.{axis}"
    pos_ghost_func = f"{pos_ghost_type.lower()}"
    neg_ghost_func = f"{neg_ghost_type.lower()}"
    return f"""
    if({this_index_var} < {nchare_var} - 1) {{
        when {pos_ghost_func}(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]) {{
            serial {{
                process_ghost({pos_ghost_type}, type_gfs, len_tmpBuffer, tmpBuffer, grid);
            }}
        }}
    }}
    if({this_index_var} > 0) {{
        when {neg_ghost_func}(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]) {{
            serial {{
                process_ghost({neg_ghost_type}, type_gfs, len_tmpBuffer, tmpBuffer, grid);
            }}
        }}
    }}
    """


# Define ghost types for each axis
x_pos_ghost_type = "EAST_GHOST"
x_neg_ghost_type = "WEST_GHOST"
y_pos_ghost_type = "NORTH_GHOST"
y_neg_ghost_type = "SOUTH_GHOST"
z_pos_ghost_type = "TOP_GHOST"
z_neg_ghost_type = "BOTTOM_GHOST"


def generate_diagnostics_code(
    dimension: str, direction: str, num_fields: str, tot_num_diagnostic_pts: str
) -> str:
    """
    Generate code for diagnostics.

    :param dimension: Dimension index
    :param direction: Direction of diagnostics (e.g. "x", "y", "z")
    :param num_fields: Number of fields for diagnostics
    :param tot_num_diagnostic_pts: Total number of diagnostic points
    :return: Code for diagnostics
    """
    code = f"""
    // {dimension} {direction}
    when ready_{dimension}_{direction}(Ck::IO::FileReadyMsg *m_{dimension}_{direction}) {{
        serial {{
        f_{dimension}_{direction} = m_{dimension}_{direction}->file;
        CkCallback sessionStart_{dimension}_{direction}(CkIndex_Timestepping::start_write_{dimension}_{direction}(0), thisProxy);
        CkCallback sessionEnd_{dimension}_{direction}(CkIndex_Timestepping::test_written_{dimension}_{direction}(0), thisProxy);
        int num_fields = {num_fields};
        int tot_num_diagnostic_pts = {tot_num_diagnostic_pts};
        int totsizeinbytes = 23 * num_fields * tot_num_diagnostic_pts;
        Ck::IO::startSession(f_{dimension}_{direction}, totsizeinbytes, 0, sessionStart_{dimension}_{direction}, sessionEnd_{dimension}_{direction});
        delete m_{dimension}_{direction};
        }}
    }}
    when start_write_{dimension}_{direction}(Ck::IO::SessionReadyMsg *m_{dimension}_{direction}){{
        serial {{
        thisProxy.diagnostics_ckio(m_{dimension}_{direction}->session, OUTPUT_{dimension.upper()}_{direction.upper()});
        delete m_{dimension}_{direction};
        }}
    }}
    when test_written_{dimension}_{direction}(CkReductionMsg *m_{dimension}_{direction}) {{
        serial {{
        delete m_{dimension}_{direction};
        CkCallback cb_{dimension}_{direction}(CkIndex_Timestepping::closed_{dimension}_{direction}(0), thisProxy);
        Ck::IO::close(f_{dimension}_{direction}, cb_{dimension}_{direction});
        count_filewritten++;
        }}
    }}
    when closed_{dimension}_{direction}(CkReductionMsg *m_{dimension}_{direction}) {{
        serial {{
        delete m_{dimension}_{direction};
        }}
    }}
    """
    return code


def register_CFunction_timestepping_malloc() -> None:
    """
    Register a C function for timestepping malloc.

    :return None
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Allocate memory for temporary buffers used to communicate face data"
    cfunc_type = "void"
    name = "timestepping_malloc_tmpBuffer"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const nonlocalinnerbc_struct *restrict nonlocalinnerbcstruct, tmpBuffers_struct *restrict tmpBuffers"
    body = """
const int Nxx_plus_2NGHOSTS_face0 = Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
const int Nxx_plus_2NGHOSTS_face1 = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS2;
const int Nxx_plus_2NGHOSTS_face2 = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1;

tmpBuffers->tmpBuffer_EW = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * NGHOSTS * Nxx_plus_2NGHOSTS_face0);
tmpBuffers->tmpBuffer_NS = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * NGHOSTS * Nxx_plus_2NGHOSTS_face1);
tmpBuffers->tmpBuffer_TB = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * NGHOSTS * Nxx_plus_2NGHOSTS_face2);

// Unpack nonlocalinnerbcstruct
  const int tot_num_dst_chares = nonlocalinnerbcstruct->tot_num_dst_chares;
  const int *restrict num_srcpts_tosend_each_chare = nonlocalinnerbcstruct->num_srcpts_tosend_each_chare;
  const int tot_num_src_chares = nonlocalinnerbcstruct->tot_num_src_chares;
  const int *restrict num_srcpts_each_chare = nonlocalinnerbcstruct->num_srcpts_each_chare;

  tmpBuffers->tmpBuffer_innerbc_send = (REAL **)malloc(tot_num_dst_chares * sizeof(REAL *));
  for (int which_chare = 0; which_chare < tot_num_dst_chares; which_chare++) {
    tmpBuffers->tmpBuffer_innerbc_send[which_chare] = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * num_srcpts_tosend_each_chare[which_chare]);
  }
  tmpBuffers->tmpBuffer_innerbc_receiv = (REAL **)malloc(tot_num_src_chares * sizeof(REAL *));
  for (int which_chare = 0; which_chare < tot_num_src_chares; which_chare++) {
    tmpBuffers->tmpBuffer_innerbc_receiv[which_chare] = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * num_srcpts_each_chare[which_chare]);
  }

"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )


def register_CFunction_timestepping_free_memory() -> None:
    """
    Register a C function for timestepping free memory.

    :return None
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Free memory for temporary buffers used to communicate face data"
    cfunc_type = "void"
    name = "timestepping_free_memory_tmpBuffer"
    params = "const nonlocalinnerbc_struct *restrict nonlocalinnerbcstruct, tmpBuffers_struct *restrict tmpBuffers"
    body = """
if (tmpBuffers->tmpBuffer_EW != NULL) free(tmpBuffers->tmpBuffer_EW);
if (tmpBuffers->tmpBuffer_NS != NULL) free(tmpBuffers->tmpBuffer_NS);
if (tmpBuffers->tmpBuffer_TB != NULL) free(tmpBuffers->tmpBuffer_TB);

// Unpack nonlocalinnerbcstruct
const int tot_num_dst_chares = nonlocalinnerbcstruct->tot_num_dst_chares;
const int tot_num_src_chares = nonlocalinnerbcstruct->tot_num_src_chares;

// Free tmpBuffer_innerbc_send
if (tmpBuffers->tmpBuffer_innerbc_send != NULL) {
  for (int which_chare = 0; which_chare < tot_num_dst_chares; which_chare++) {
    if (tmpBuffers->tmpBuffer_innerbc_send[which_chare] != NULL) {
      free(tmpBuffers->tmpBuffer_innerbc_send[which_chare]);
    }
  }
  free(tmpBuffers->tmpBuffer_innerbc_send);
}
// Free tmpBuffer_innerbc_receiv
if (tmpBuffers->tmpBuffer_innerbc_receiv != NULL) {
  for (int which_chare = 0; which_chare < tot_num_src_chares; which_chare++) {
    if (tmpBuffers->tmpBuffer_innerbc_receiv[which_chare] != NULL) {
      free(tmpBuffers->tmpBuffer_innerbc_receiv[which_chare]);
    }
  }
  free(tmpBuffers->tmpBuffer_innerbc_receiv);
}
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )


def output_timestepping_h(
    project_dir: str,
    enable_residual_diagnostics: bool = False,
    enable_psi4_diagnostics: bool = False,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    """
    Generate timestepping.h.

    :param project_dir: Directory where the project C code is output.
    :param enable_residual_diagnostics: Flag to enable residual diagnostics, default is False.
    :param enable_psi4_diagnostics: Whether or not to enable psi4 diagnostics.
    :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = r"""#ifndef __TIMESTEPPING_H__
#define __TIMESTEPPING_H__
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "timestepping.decl.h"
"""
    if enable_psi4_diagnostics:
        file_output_str += r"""
struct sectionBcastMsg : public CkMcastBaseMsg, public CMessage_sectionBcastMsg {
  int k;
  sectionBcastMsg(int _k) : k(_k) {}
  void pup(PUP::er &p) {
    CMessage_sectionBcastMsg::pup(p);
    p|k;
  }
};
"""
    file_output_str += r"""
class Timestepping : public CBase_Timestepping {
  Timestepping_SDAG_CODE

  private:
    /// Member Variables (Object State) ///"""
    if enable_psi4_diagnostics:
        file_output_str += r"""
    CProxySection_Timestepping secProxy;
    CkSectionInfo cookie;"""
    file_output_str += r"""
    commondata_struct commondata;
    griddata_struct *griddata;
    griddata_struct *griddata_chare;
    bool is_boundarychare;
    REAL time_start;
    //bool contains_gridcenter;
    const int grid = 0;
    const int which_grid_diagnostics = 0;
    bool write_diagnostics_this_step;
    Ck::IO::File f_1d_y;
    Ck::IO::File f_1d_z;
    Ck::IO::File f_2d_xy;
    Ck::IO::File f_2d_yz;
    int count_filewritten = 0;
    const int expected_count_filewritten = 4;
    int iter = 0;
    int type_gfs_nonlocal_innerbc;
"""
    file_output_str += r"""
    /// Member Functions (private) ///
    void send_neighbor_data(const int type_gfs, const int dir, const int grid);
    void process_ghost(const int type_ghost, const int type_gfs, const int len_tmpBuffer, const REAL *restrict tmpBuffer, const int grid);
    void send_nonlocalinnerbc_idx3srcpts_toreceiv();
    void process_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, const int *restrict globalidx3_srcpts);
    void send_nonlocalinnerbc_data(const int type_gfs, const int grid);
    void set_tmpBuffer_innerbc_receiv(const int src_chare_idx3, const int len_tmpBuffer, const REAL *restrict vals, const int grid);
    void process_nonlocalinnerbc(const int type_gfs, const int grid);"""
    if enable_residual_diagnostics:
        file_output_str += r"""
    void contribute_localsums_for_residualH(REAL localsums_for_residualH[2]);
    void send_wavespeed_at_outer_boundary(const int grid);"""
    if enable_psi4_diagnostics:
        file_output_str += r"""
    void contribute_localsums_for_psi4_decomp(sectionBcastMsg *msg, const int grid);
"""
    file_output_str += r"""
  public:
    /// Constructors ///
    Timestepping(CommondataObject &&inData);
    Timestepping(CkMigrateMessage* msg);
    /// Destructor ///
    ~Timestepping();

    /// Entry Methods ///
};

#endif //__TIMESTEPPING_H__
"""
    timestepping_h_file = project_Path / "timestepping.h"
    with timestepping_h_file.open("w", encoding="utf-8") as file:
        file.write(
            clang_format(file_output_str, clang_format_options=clang_format_options)
        )


def generate_switch_statement_for_gf_types(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
) -> str:
    """
    Generate the switch statement for grid function types based on the given Method of Lines (MoL) method.

    :param Butcher_dict: Dictionary containing Butcher tableau data.
    :param MoL_method: Method of Lines (MoL) method name.
    :return: A string representing the switch statement for the grid function types.
    """
    # Generating gridfunction names based on the given MoL method
    (
        y_n_gridfunctions,
        non_y_n_gridfunctions_list,
        _diagnostic_gridfunctions_point_to,
        _diagnostic_gridfunctions2_point_to,
    ) = generate_gridfunction_names(Butcher_dict, MoL_method=MoL_method)

    # Convert y_n_gridfunctions to a list if it's a string
    gf_list = (
        [y_n_gridfunctions] if isinstance(y_n_gridfunctions, str) else y_n_gridfunctions
    )
    gf_list.extend(non_y_n_gridfunctions_list)

    switch_statement = """
switch (type_gfs) {
"""
    switch_cases = []
    for gf in gf_list:
        switch_cases.append(f"  case {gf.upper()}:")
        switch_cases.append(f"    gfs = griddata_chare[grid].gridfuncs.{gf.lower()};")
        switch_cases.append("    break;")
    switch_cases.append(
        """
  default:
    break;
}
"""
    )

    switch_body = "\n".join(switch_cases)
    return switch_statement + switch_body


def output_timestepping_cpp(
    project_dir: str,
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
    initial_data_desc: str = "",
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
    initialize_constant_auxevol: bool = False,
    enable_residual_diagnostics: bool = False,
    enable_psi4_diagnostics: bool = False,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    """
    Generate timestepping.cpp.

    :param project_dir: Directory where the project C code is output.
    :param Butcher_dict: Dictionary containing Butcher tableau for the MoL method.
    :param MoL_method: Method of Lines (MoL) method name.
    :param initial_data_desc: Description for initial data, default is an empty string.
    :param enable_rfm_precompute: Enable rfm precomputation, default is False.
    :param enable_CurviBCs: Enable CurviBCs, default is False.
    :param initialize_constant_auxevol: If set to True, `initialize_constant_auxevol` function will be called during the simulation initialization phase to set these constants. Default is False.
    :param enable_residual_diagnostics: Enable residual diagnostics, default is False.
    :param enable_psi4_diagnostics: Whether or not to enable psi4 diagnostics.
    :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".
    :raises ValueError: Raised if any required function is not registered.
    """
    initial_data_desc += " "
    # Make sure all required C functions are registered
    missing_functions: List[Tuple[str, str]] = []
    for func_tuple in [
        ("params_struct_set_to_default", "CodeParameters.py"),
        (
            "numerical_grids_and_timestep",
            "e.g., numerical_grids_and_timestep.py or user defined",
        ),
        ("MoL_malloc_y_n_gfs", "MoL.py"),
        ("MoL_malloc_non_y_n_gfs", "MoL.py"),
        ("initial_data", "initial_data.py"),
        ("MoL_step_forward_in_time", "MoL.py"),
        ("diagnostics", "diagnostics.py"),
        ("MoL_free_memory_y_n_gfs", "MoL.py"),
        ("MoL_free_memory_non_y_n_gfs", "MoL.py"),
    ]:
        if func_tuple[0] not in cfc.CFunction_dict:
            missing_functions += [func_tuple]
    if missing_functions:
        error_msg = "Error: These functions are required and are not registered.\n"
        for func_tuple in missing_functions:
            error_msg += (
                f'  {func_tuple[0]}, registered by function within "{func_tuple[1]}"\n'
            )
        raise ValueError(error_msg)

    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = r"""#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "timestepping.h"
#include "main.h"

extern /* readonly */ CProxy_Main mainProxy;

/*
*Step 1.c: Allocate NUMGRIDS griddata arrays, each containing data specific to an individual grid.
*Step 1.d: Set each CodeParameter in griddata.params to default.
"""
    if enable_CurviBCs:
        file_output_str += "*Step 1.e: Set non-parfile parameters related to numerical grid, then set up numerical grids and CFL-limited timestep.\n"
    if enable_rfm_precompute:
        file_output_str += "*Step 1.f: Set up boundary condition struct (bcstruct).\n"
    file_output_str += rf"""Step 2: Initial data are set on y_n_gfs gridfunctions. Allocate storage for them first.
*Step 3: Finalize initialization: set up {initial_data_desc}initial data, etc.
*Step 4: Allocate storage for non-y_n gridfunctions, needed for the Runge-Kutta-like timestepping.
"""
    if initialize_constant_auxevol:
        file_output_str += (
            "*Step 4.a: Set AUXEVOL gridfunctions that will never change in time."
        )
    file_output_str += r"""*/
"""
    file_output_str += r"""
Timestepping::Timestepping(CommondataObject &&inData) {

  commondata = inData.commondata;

  if (thisIndex.x == 0 || thisIndex.y == 0 || thisIndex.z == 0 || thisIndex.x == (commondata.Nchare0 - 1) || thisIndex.y == (commondata.Nchare1 - 1) ||
      thisIndex.z == (commondata.Nchare2 - 1)) {
    is_boundarychare = true;
  }

  // Step 1.c: Allocate NUMGRIDS griddata arrays, each containing data specific to an individual grid.
  griddata = (griddata_struct *restrict)malloc(sizeof(griddata_struct) * commondata.NUMGRIDS);
  griddata_chare = (griddata_struct *restrict)malloc(sizeof(griddata_struct) * commondata.NUMGRIDS);

  // Step 1.d: Set each CodeParameter in griddata.params to default.
  params_struct_set_to_default(&commondata, griddata);
  params_struct_set_to_default(&commondata, griddata_chare);

  // Step 1.e: Set up numerical grids: xx[3], masks, Nxx, dxx, invdxx, bcstruct, rfm_precompute, timestep, etc.
  {
    // if calling_for_first_time, then initialize commondata time=nn=t_0=nn_0 = 0
    const bool calling_for_first_time = true;
    numerical_grids_and_timestep(&commondata, griddata, calling_for_first_time);

    const int thisIndex_arr[3] = {thisIndex.x, thisIndex.y, thisIndex.z};
    numerical_grids_chare(&commondata, griddata, griddata_chare, thisIndex_arr);
  }

  // Step 2: Initial data are set on y_n_gfs gridfunctions. Allocate storage for them first.
  for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
    MoL_malloc_y_n_gfs(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].gridfuncs);
  }

  // Step 3: Allocate storage for non-y_n gridfunctions, needed for the Runge-Kutta-like timestepping
  for(int grid=0; grid<commondata.NUMGRIDS; grid++)
    MoL_malloc_non_y_n_gfs(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].gridfuncs);
"""

    file_output_str += """
  // Allocate storage for diagnostic gridfunctions
  for(int grid=0; grid<commondata.NUMGRIDS; grid++)
    MoL_malloc_diagnostic_gfs(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].gridfuncs);

  // Allocate storage for temporary buffers, needed for communicating face data
  for(int grid=0; grid<commondata.NUMGRIDS; grid++)
    timestepping_malloc_tmpBuffer(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].nonlocalinnerbcstruct, &griddata_chare[grid].tmpBuffers);

"""
    if initialize_constant_auxevol:
        file_output_str += r"""
  // Step 4.a: Set AUXEVOL gridfunctions that will never change in time.
  initialize_constant_auxevol(&commondata, griddata_chare);
"""
    file_output_str += """
}
"""
    file_output_str += r"""
// migration constructor
Timestepping::Timestepping(CkMigrateMessage *msg) { }

// destructor
Timestepping::~Timestepping() {
  // Step 5: Free all allocated memory
  for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
    MoL_free_memory_y_n_gfs(&griddata_chare[grid].gridfuncs);
    MoL_free_memory_non_y_n_gfs(&griddata_chare[grid].gridfuncs);
    MoL_free_memory_diagnostic_gfs(&griddata_chare[grid].gridfuncs);
    timestepping_free_memory_tmpBuffer(&griddata_chare[grid].nonlocalinnerbcstruct, &griddata_chare[grid].tmpBuffers);"""
    if enable_rfm_precompute:
        file_output_str += r"""
    rfm_precompute_free(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].rfmstruct);"""
    if enable_CurviBCs:
        file_output_str += r"""
    free(griddata[grid].bcstruct.inner_bc_array);
    free(griddata_chare[grid].bcstruct.inner_bc_array);
    free(griddata_chare[grid].bcstruct.inner_bc_array_nonlocal);
    for(int ng=0;ng<NGHOSTS*3;ng++) {
     free(griddata[grid].bcstruct.pure_outer_bc_array[ng]);
     free(griddata_chare[grid].bcstruct.pure_outer_bc_array[ng]);
    }
    for(int i=0;i<3;i++) {
      free(griddata[grid].xx[i]);
      free(griddata_chare[grid].xx[i]);
    }
    free(griddata_chare[grid].diagnosticstruct.localidx3_diagnostic_1d_y_pt);
    free(griddata_chare[grid].diagnosticstruct.locali0_diagnostic_1d_y_pt);
    free(griddata_chare[grid].diagnosticstruct.locali1_diagnostic_1d_y_pt);
    free(griddata_chare[grid].diagnosticstruct.locali2_diagnostic_1d_y_pt);
    free(griddata_chare[grid].diagnosticstruct.offset_diagnostic_1d_y_pt);
    free(griddata_chare[grid].diagnosticstruct.localidx3_diagnostic_1d_z_pt);
    free(griddata_chare[grid].diagnosticstruct.locali0_diagnostic_1d_z_pt);
    free(griddata_chare[grid].diagnosticstruct.locali1_diagnostic_1d_z_pt);
    free(griddata_chare[grid].diagnosticstruct.locali2_diagnostic_1d_z_pt);
    free(griddata_chare[grid].diagnosticstruct.offset_diagnostic_1d_z_pt);
    free(griddata_chare[grid].diagnosticstruct.localidx3_diagnostic_2d_xy_pt);
    free(griddata_chare[grid].diagnosticstruct.locali0_diagnostic_2d_xy_pt);
    free(griddata_chare[grid].diagnosticstruct.locali1_diagnostic_2d_xy_pt);
    free(griddata_chare[grid].diagnosticstruct.locali2_diagnostic_2d_xy_pt);
    free(griddata_chare[grid].diagnosticstruct.offset_diagnostic_2d_xy_pt);
    free(griddata_chare[grid].diagnosticstruct.localidx3_diagnostic_2d_yz_pt);
    free(griddata_chare[grid].diagnosticstruct.locali0_diagnostic_2d_yz_pt);
    free(griddata_chare[grid].diagnosticstruct.locali1_diagnostic_2d_yz_pt);
    free(griddata_chare[grid].diagnosticstruct.locali2_diagnostic_2d_yz_pt);
    free(griddata_chare[grid].diagnosticstruct.offset_diagnostic_2d_yz_pt);
    if (griddata_chare[grid].diagnosticstruct.list_of_R_exts_chare != NULL) {
      free(griddata_chare[grid].diagnosticstruct.list_of_R_exts_chare);
    }
    if (griddata_chare[grid].diagnosticstruct.localsums_for_psi4_decomp != NULL) {
      free(griddata_chare[grid].diagnosticstruct.localsums_for_psi4_decomp);
    }
    if (griddata_chare[grid].diagnosticstruct.globalsums_for_psi4_decomp != NULL) {
      free(griddata_chare[grid].diagnosticstruct.globalsums_for_psi4_decomp);
    }
    free(griddata_chare[grid].charecommstruct.globalidx3pt_to_chareidx3);
    free(griddata_chare[grid].charecommstruct.globalidx3pt_to_localidx3pt);
    free(griddata_chare[grid].charecommstruct.localidx3pt_to_globalidx3pt);
    free(griddata_chare[grid].nonlocalinnerbcstruct.idx3_of_src_chares);
    free(griddata_chare[grid].nonlocalinnerbcstruct.idx3chare_to_src_chare_id);
    free(griddata_chare[grid].nonlocalinnerbcstruct.num_srcpts_each_chare);
    for (int i = 0; i < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares; i++) {
      free(griddata_chare[grid].nonlocalinnerbcstruct.map_srcchare_and_srcpt_id_to_linear_id[i]);
    }
    free(griddata_chare[grid].nonlocalinnerbcstruct.map_srcchare_and_srcpt_id_to_linear_id);
    for (int i = 0; i < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares; i++) {
      free(griddata_chare[grid].nonlocalinnerbcstruct.globalidx3_srcpts[i]);
    }
    free(griddata_chare[grid].nonlocalinnerbcstruct.globalidx3_srcpts);
    free(griddata_chare[grid].nonlocalinnerbcstruct.idx3_of_dst_chares);
    free(griddata_chare[grid].nonlocalinnerbcstruct.idx3chare_to_dst_chare_id);
    free(griddata_chare[grid].nonlocalinnerbcstruct.num_srcpts_tosend_each_chare);
    for (int i = 0; i < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares; i++) {
      free(griddata_chare[grid].nonlocalinnerbcstruct.globalidx3_srcpts_tosend[i]);
    }
    free(griddata_chare[grid].nonlocalinnerbcstruct.globalidx3_srcpts_tosend);
  }
  free(griddata);
  free(griddata_chare);
}
"""

    switch_case_code = generate_switch_statement_for_gf_types(Butcher_dict, MoL_method)

    file_output_str += r"""
// send NGHOSTS number of interior faces with face extents that include ghosts
void Timestepping::send_neighbor_data(const int type_gfs, const int dir, const int grid) {
  const int Nchare0 = commondata.Nchare0;
  const int Nchare1 = commondata.Nchare1;
  const int Nchare2 = commondata.Nchare2;
  const int Nxx0 = griddata_chare[grid].params.Nxx0;
  const int Nxx1 = griddata_chare[grid].params.Nxx1;
  const int Nxx2 = griddata_chare[grid].params.Nxx2;
  const int Nxx_plus_2NGHOSTS0 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS2;
  REAL *restrict tmpBuffer_EW = griddata_chare[grid].tmpBuffers.tmpBuffer_EW;
  REAL *restrict tmpBuffer_NS = griddata_chare[grid].tmpBuffers.tmpBuffer_NS;
  REAL *restrict tmpBuffer_TB = griddata_chare[grid].tmpBuffers.tmpBuffer_TB;
  const REAL *restrict gfs = nullptr;
  """
    file_output_str += switch_case_code
    file_output_str += r"""
  switch (dir) {
    case EAST_WEST:
      //send to west
      if (thisIndex.x > 0) {
        for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
          int i0 = 2*NGHOSTS - 1;
          for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
            for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
              for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
                tmpBuffer_EW[IDXFACES0(which_gf, which_inner, i1, i2)] = gfs[IDX4(which_gf, i0, i1, i2)];
              }
            }
            i0--;
          }
        }
        thisProxy[CkArrayIndex3D(thisIndex.x - 1, thisIndex.y, thisIndex.z)].east_ghost(type_gfs, NUM_EVOL_GFS*NGHOSTS*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2, tmpBuffer_EW);
      }
      //send to east
      if (thisIndex.x < Nchare0 - 1) {
        for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
          int i0 = Nxx0;
          for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
            for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
              for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
                tmpBuffer_EW[IDXFACES0(which_gf, which_inner, i1, i2)] = gfs[IDX4(which_gf, i0, i1, i2)];
              }
            }
            i0++;
          }
        }
        thisProxy[CkArrayIndex3D(thisIndex.x + 1, thisIndex.y, thisIndex.z)].west_ghost(type_gfs, NUM_EVOL_GFS*NGHOSTS*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2, tmpBuffer_EW);
      }
      break;
    case NORTH_SOUTH:
      //send to south
      if (thisIndex.y > 0) {
        for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
          int i1 = 2*NGHOSTS - 1;
          for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
            for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
              for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
                tmpBuffer_NS[IDXFACES1(which_gf, which_inner, i0, i2)] = gfs[IDX4(which_gf, i0, i1, i2)];
              }
            }
            i1--;
          }
        }
        thisProxy[CkArrayIndex3D(thisIndex.x, thisIndex.y - 1, thisIndex.z)].north_ghost(type_gfs, NUM_EVOL_GFS*NGHOSTS*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS2, tmpBuffer_NS);
      }
      //send to north
      if (thisIndex.y < Nchare1 - 1) {
        for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
          int i1 = Nxx1;
          for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
            for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
              for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
                tmpBuffer_NS[IDXFACES1(which_gf, which_inner, i0, i2)] = gfs[IDX4(which_gf, i0, i1, i2)];
              }
            }
            i1++;
          }
        }
        thisProxy[CkArrayIndex3D(thisIndex.x, thisIndex.y + 1, thisIndex.z)].south_ghost(type_gfs, NUM_EVOL_GFS*NGHOSTS*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS2, tmpBuffer_NS);
      }
      break;
    case TOP_BOTTOM:
      //send to bottom
      if (thisIndex.z > 0) {
        for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
          int i2 = 2*NGHOSTS - 1;
          for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
            for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
              for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
                tmpBuffer_TB[IDXFACES2(which_gf, which_inner, i0, i1)] = gfs[IDX4(which_gf, i0, i1, i2)];
              }
            }
            i2--;
          }
        }
        thisProxy[CkArrayIndex3D(thisIndex.x, thisIndex.y, thisIndex.z - 1)].top_ghost(type_gfs, NUM_EVOL_GFS*NGHOSTS*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1, tmpBuffer_TB);
      }
      //send to top
      if (thisIndex.z < Nchare2 - 1) {
        for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
          int i2 = Nxx2;
          for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
            for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
              for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
                tmpBuffer_TB[IDXFACES2(which_gf, which_inner, i0, i1)] = gfs[IDX4(which_gf, i0, i1, i2)];
              }
            }
            i2++;
          }
        }
        thisProxy[CkArrayIndex3D(thisIndex.x, thisIndex.y, thisIndex.z + 1)].bottom_ghost(type_gfs, NUM_EVOL_GFS*NGHOSTS*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1, tmpBuffer_TB);
      }
      break;
    default:
      return;
    }
}
"""
    file_output_str += r"""
// process neighbor ghosts
void Timestepping::process_ghost(const int type_ghost, const int type_gfs, const int len_tmpBuffer, const REAL *restrict vals, const int grid) {
  const int Nxx0 = griddata_chare[grid].params.Nxx0;
  const int Nxx1 = griddata_chare[grid].params.Nxx1;
  const int Nxx2 = griddata_chare[grid].params.Nxx2;
  const int Nxx_plus_2NGHOSTS0 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS2;

  REAL *restrict gfs = nullptr;
"""
    switch_case_code = generate_switch_statement_for_gf_types(Butcher_dict, MoL_method)
    file_output_str += switch_case_code
    file_output_str += r"""
  switch (type_ghost) {
    case EAST_GHOST:
      for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
        int i0 = Nxx0 + (2 * NGHOSTS) - 1;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
            for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
              gfs[IDX4(which_gf, i0, i1, i2)] = vals[IDXFACES0(which_gf, which_inner, i1, i2)];
            }
          }
          i0--;
        }
      }
      break;
    case WEST_GHOST:
      for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
        int i0 = 0;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
            for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
              gfs[IDX4(which_gf, i0, i1, i2)] = vals[IDXFACES0(which_gf, which_inner, i1, i2)];
            }
          }
          i0++;
        }
      }
      break;
    case NORTH_GHOST:
      for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
        int i1 = Nxx1 + (2 * NGHOSTS) - 1;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
              gfs[IDX4(which_gf, i0, i1, i2)] = vals[IDXFACES1(which_gf, which_inner, i0, i2)];
            }
          }
          i1--;
        }
      }
      break;
    case SOUTH_GHOST:
      for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
        int i1 = 0;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
              gfs[IDX4(which_gf, i0, i1, i2)] = vals[IDXFACES1(which_gf, which_inner, i0, i2)];
            }
          }
          i1++;
        }
      }
      break;
    case TOP_GHOST:
      for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
        int i2 = Nxx2 + (2 * NGHOSTS) - 1;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
              gfs[IDX4(which_gf, i0, i1, i2)] = vals[IDXFACES2(which_gf, which_inner, i0, i1)];
            }
          }
          i2--;
        }
      }
      break;
    case BOTTOM_GHOST:
      for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
        int i2 = 0;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
              gfs[IDX4(which_gf, i0, i1, i2)] = vals[IDXFACES2(which_gf, which_inner, i0, i1)];
            }
          }
          i2++;
        }
      }
      break;
  }
}
"""

    if enable_residual_diagnostics:
        file_output_str += r"""
void Timestepping::contribute_localsums_for_residualH(REAL localsums_for_residualH[2]) {
  std::vector<double> outdoubles(2);
  outdoubles[0] = localsums_for_residualH[0];
  outdoubles[1] = localsums_for_residualH[1];
  CkCallback cb(CkIndex_Timestepping::report_sums_for_residualH(NULL), thisProxy);
  contribute(outdoubles, CkReduction::sum_double, cb);
}

void Timestepping::send_wavespeed_at_outer_boundary(const int grid) {
  const int Nchare0 = commondata.Nchare0;
  const int Nchare1 = commondata.Nchare1;
  const int Nchare2 = commondata.Nchare2;
  const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
  const int Nxx0chare = griddata_chare[grid].params.Nxx0;
  const int Nxx1chare = griddata_chare[grid].params.Nxx1;
  const int Nxx2chare = griddata_chare[grid].params.Nxx2;
  const int Nxx_plus_2NGHOSTS0chare = griddata_chare[grid].params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1chare = griddata_chare[grid].params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2chare = griddata_chare[grid].params.Nxx_plus_2NGHOSTS2;

  const int pt_at_outer_boundary_i0 = Nxx_plus_2NGHOSTS0 - NGHOSTS - 1;
  const int pt_at_outer_boundary_i1 = NGHOSTS;
  const int pt_at_outer_boundary_i2 = Nxx_plus_2NGHOSTS2 / 2 ;
  const int pt_at_outer_boundary_idx3 =  IDX3(pt_at_outer_boundary_i0, pt_at_outer_boundary_i1, pt_at_outer_boundary_i2);
  const int idx3_this_chare = IDX3_OF_CHARE(thisIndex.x, thisIndex.y, thisIndex.z);

  if (griddata_chare[grid].charecommstruct.globalidx3pt_to_chareidx3[pt_at_outer_boundary_idx3] == idx3_this_chare) {
    const int locali0 = MAP_GLOBAL_TO_LOCAL_IDX0(thisIndex.x, pt_at_outer_boundary_i0, Nxx0chare);
    const int locali1 = MAP_GLOBAL_TO_LOCAL_IDX1(thisIndex.y, pt_at_outer_boundary_i1, Nxx1chare);
    const int locali2 = MAP_GLOBAL_TO_LOCAL_IDX2(thisIndex.z, pt_at_outer_boundary_i2, Nxx2chare);
    const REAL wavespeed_at_outer_boundary = griddata_chare[grid].gridfuncs.auxevol_gfs[IDX4GENERAL(VARIABLE_WAVESPEEDGF, locali0, locali1, locali2, Nxx_plus_2NGHOSTS0chare, Nxx_plus_2NGHOSTS1chare, Nxx_plus_2NGHOSTS2chare)];
    thisProxy.receiv_wavespeed_at_outer_boundary(wavespeed_at_outer_boundary);
  }
}
"""
    if enable_psi4_diagnostics:
        file_output_str += r"""
void Timestepping::contribute_localsums_for_psi4_decomp(sectionBcastMsg *msg, const int grid) {
  // Unpack diagnosticptoffset struct:
  const int length_localsums_for_psi4_decomp = griddata_chare[grid].diagnosticstruct.length_localsums_for_psi4_decomp;
  const REAL *restrict localsums_for_psi4_decomp = griddata_chare[grid].diagnosticstruct.localsums_for_psi4_decomp;

  // Initialize outdoubles with the correct size
  std::vector<double> outdoubles(length_localsums_for_psi4_decomp);

  // Copy and convert data from localsums_for_psi4_decomp to outdoubles
  for (int i = 0; i < length_localsums_for_psi4_decomp; ++i) {
    outdoubles[i] = static_cast<double>(localsums_for_psi4_decomp[i]);
  }

  CkGetSectionInfo(cookie, msg);
  CkCallback cb(CkIndex_Timestepping::report_sums_for_psi4_diagnostics(NULL), thisProxy[CkArrayIndex3D(thisIndex.x, 0, 0)]);
  CProxySection_Timestepping::contribute(outdoubles, CkReduction::sum_double, cookie, cb);
  delete msg;
}
"""

    file_output_str += r"""
void Timestepping::send_nonlocalinnerbc_idx3srcpts_toreceiv() {
  // Unpack griddata_chare[grid].nonlocalinnerbcstruct
  int grid = 0;
  const int tot_num_src_chares = griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares;
  const int *restrict idx3_of_src_chares = griddata_chare[grid].nonlocalinnerbcstruct.idx3_of_src_chares;
  int **restrict globalidx3_srcpts = griddata_chare[grid].nonlocalinnerbcstruct.globalidx3_srcpts;
  const int *restrict num_srcpts_each_chare = griddata_chare[grid].nonlocalinnerbcstruct.num_srcpts_each_chare;

  const int Nchare0 = commondata.Nchare0;
  const int Nchare1 = commondata.Nchare1;
  const int Nchare2 = commondata.Nchare2;
  const int idx3_of_thischare = IDX3_OF_CHARE(thisIndex.x, thisIndex.y, thisIndex.z);

  for (int src_chare_id = 0; src_chare_id < tot_num_src_chares; src_chare_id++) {
    int idx3srcchare = idx3_of_src_chares[src_chare_id];
    int src_chare_index0;
    int src_chare_index1;
    int src_chare_index2;
    REVERSE_IDX3GENERAL(idx3srcchare, Nchare0, Nchare1, src_chare_index0, src_chare_index1, src_chare_index2);
    int num_srcpts  = num_srcpts_each_chare[src_chare_id];
    thisProxy[CkArrayIndex3D(src_chare_index0, src_chare_index1, src_chare_index2)].receiv_nonlocalinnerbc_idx3srcpt_tosend(idx3_of_thischare, num_srcpts, globalidx3_srcpts[src_chare_id]);
  }
}
"""

    file_output_str += r"""
void Timestepping::process_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, const int *restrict globalidx3_srcpts) {
  // Unpack griddata_chare[grid].nonlocalinnerbcstruct
  int grid = 0;
  const int *restrict idx3chare_to_dst_chare_id = griddata_chare[grid].nonlocalinnerbcstruct.idx3chare_to_dst_chare_id;
  int **restrict globalidx3_srcpts_tosend = griddata_chare[grid].nonlocalinnerbcstruct.globalidx3_srcpts_tosend;

  const int dst_chare_id_val = idx3chare_to_dst_chare_id[idx3_of_sendingchare];
  for (int count = 0; count < num_srcpts; count++) {
    globalidx3_srcpts_tosend[dst_chare_id_val][count] =  globalidx3_srcpts[count];
  }
}
"""

    file_output_str += r"""
void Timestepping::send_nonlocalinnerbc_data(const int type_gfs, const int grid) {
  const int Nchare0 = commondata.Nchare0;
  const int Nchare1 = commondata.Nchare1;
  const int Nchare2 = commondata.Nchare2;
  const int Nxx_plus_2NGHOSTS0 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS2;

  // Unpack nonlocalinnerbcstruct and charecommstruct
  const int tot_num_dst_chares = griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares;
  const int *restrict idx3_of_dst_chares = griddata_chare[grid].nonlocalinnerbcstruct.idx3_of_dst_chares;
  const int *restrict num_srcpts_tosend_each_chare = griddata_chare[grid].nonlocalinnerbcstruct.num_srcpts_tosend_each_chare;
  int **restrict globalidx3_srcpts_tosend = griddata_chare[grid].nonlocalinnerbcstruct.globalidx3_srcpts_tosend;
  const int *restrict globalidx3pt_to_localidx3pt = griddata_chare[grid].charecommstruct.globalidx3pt_to_localidx3pt;

  const REAL *restrict gfs = nullptr;
  """
    file_output_str += switch_case_code
    file_output_str += r"""
  const int idx3_this_chare = IDX3_OF_CHARE(thisIndex.x, thisIndex.y, thisIndex.z);

  for (int which_dst_chare = 0; which_dst_chare < tot_num_dst_chares; which_dst_chare++) {
    REAL *restrict tmpBuffer_innerbc_send = griddata_chare[grid].tmpBuffers.tmpBuffer_innerbc_send[which_dst_chare];
    for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
      for (int which_srcpt = 0; which_srcpt < num_srcpts_tosend_each_chare[which_dst_chare]; which_srcpt++) {
        const int globalidx3srcpt = globalidx3_srcpts_tosend[which_dst_chare][which_srcpt];
        const int localidx3srcpt = globalidx3pt_to_localidx3pt[globalidx3srcpt];
        const int idx2 = IDX2NONLOCALINNERBC(which_gf, which_srcpt, num_srcpts_tosend_each_chare[which_dst_chare]);
        tmpBuffer_innerbc_send[idx2] = gfs[IDX4pt(which_gf, localidx3srcpt)];
      }
    }
    int dst_chare_index0;
    int dst_chare_index1;
    int dst_chare_index2;
    REVERSE_IDX3GENERAL(idx3_of_dst_chares[which_dst_chare], Nchare0, Nchare1, dst_chare_index0, dst_chare_index1, dst_chare_index2);
    thisProxy[CkArrayIndex3D(dst_chare_index0, dst_chare_index1, dst_chare_index2)].receiv_nonlocalinnerbc_data(idx3_this_chare, type_gfs, NUM_EVOL_GFS * num_srcpts_tosend_each_chare[which_dst_chare], tmpBuffer_innerbc_send);
	}
}
"""

    file_output_str += r"""
void Timestepping::set_tmpBuffer_innerbc_receiv(const int src_chare_idx3, const int len_tmpBuffer, const REAL *restrict vals, const int grid) {
  const int src_chare_id = griddata_chare[grid].nonlocalinnerbcstruct.idx3chare_to_src_chare_id[src_chare_idx3];
  REAL *restrict tmpBuffer_innerbc_receiv = griddata_chare[grid].tmpBuffers.tmpBuffer_innerbc_receiv[src_chare_id];

  for (int i = 0; i < len_tmpBuffer; i++) {
    tmpBuffer_innerbc_receiv[i] = vals[i];
  }
}
"""

    file_output_str += r"""
void Timestepping::process_nonlocalinnerbc(const int type_gfs, const int grid) {
  REAL *restrict gfs = nullptr;
"""
    file_output_str += switch_case_code
    file_output_str += r"""
  apply_bcs_inner_only_nonlocal(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].bcstruct, &griddata_chare[grid].nonlocalinnerbcstruct, gfs, griddata_chare[grid].tmpBuffers.tmpBuffer_innerbc_receiv);
}
"""
    file_output_str += r"""
#include "timestepping.def.h"
"""

    timestepping_cpp_file = project_Path / "timestepping.cpp"
    with timestepping_cpp_file.open("w", encoding="utf-8") as file:
        file.write(
            clang_format(file_output_str, clang_format_options=clang_format_options)
        )


def output_timestepping_ci(
    project_dir: str,
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
    post_non_y_n_auxevol_mallocs: str,
    pre_MoL_step_forward_in_time: str = "",
    post_MoL_step_forward_in_time: str = "",
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
    enable_psi4_diagnostics: bool = False,
    enable_residual_diagnostics: bool = False,
) -> None:
    """
    Generate timestepping.ci.

    :param project_dir: Directory where the project C code is output.
    :param Butcher_dict: Dictionary containing Butcher tableau for the MoL method.
    :param MoL_method: Method of Lines (MoL) method name.
    :param post_non_y_n_auxevol_mallocs: Function calls after memory is allocated for non y_n and auxevol gridfunctions, default is an empty string.
    :param pre_MoL_step_forward_in_time: Code for handling pre-right-hand-side operations, default is an empty string.
    :param post_MoL_step_forward_in_time: Code for handling post-right-hand-side operations, default is an empty string.
    :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".
    :param enable_psi4_diagnostics: Whether or not to enable psi4 diagnostics.
    :param enable_residual_diagnostics: Whether or not to enable residual diagnostics.
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = r"""module timestepping {
  include "BHaH_defines.h";
  include "BHaH_function_prototypes.h";
  include "commondata_object.h";
  include "ckio.h";
  """
    if enable_psi4_diagnostics:
        file_output_str += r"""
  message sectionBcastMsg;
        """
    file_output_str += r"""
  array [3D] Timestepping {
    entry Timestepping(CommondataObject &inData);
    entry void ready_1d_y(Ck::IO::FileReadyMsg *m);
    entry void ready_1d_z(Ck::IO::FileReadyMsg *m);
    entry void ready_2d_xy(Ck::IO::FileReadyMsg *m);
    entry void ready_2d_yz(Ck::IO::FileReadyMsg *m);
    // Step 5: MAIN SIMULATION LOOP
    entry void start() {"""

    if enable_psi4_diagnostics:
        file_output_str += r"""
      serial {
        secProxy = CProxySection_Timestepping::ckNew(thisProxy.ckGetArrayID(), thisIndex.x, thisIndex.x, 1, 0, commondata.Nchare1 - 1, 1, 0, commondata.Nchare2 - 1, 1);
      }"""
    file_output_str += r"""
      if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares > 0) {
        serial {
          send_nonlocalinnerbc_idx3srcpts_toreceiv();
        }
      }
      if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares > 0) {
        for (iter = 0; iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares; iter++) {
          when receiv_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, int globalidx3_srcpts[num_srcpts]) {
            serial {
              process_nonlocalinnerbc_idx3srcpt_tosend(idx3_of_sendingchare, num_srcpts, globalidx3_srcpts);
            }
          }
        }
      }"""

    if not enable_residual_diagnostics:
        # If anything other than NRPy elliptic
        file_output_str += """
      serial {
        initial_data(&commondata, griddata_chare, INITIALDATA_LOOPOVERALLGRIDPTS);
        initial_data(&commondata, griddata_chare, INITIALDATA_APPLYBCSINNERONLY);
      }"""
        file_output_str += generate_send_nonlocalinnerbc_data_code("Y_N_GFS")
        file_output_str += generate_process_nonlocalinnerbc_code()
        file_output_str += """
      serial {
        initial_data(&commondata, griddata_chare, INITIALDATA_LAMBDAUGRIDINTERIOR);
        initial_data(&commondata, griddata_chare, INITIALDATA_APPLYBCSOUTEREXTRAPANDINNER);
      }"""
        file_output_str += generate_send_nonlocalinnerbc_data_code("Y_N_GFS")
        file_output_str += generate_process_nonlocalinnerbc_code()
    else:
        # If NRPy elliptic
        file_output_str += """
      serial {
        initial_data(&commondata, griddata_chare);
"""
        if post_non_y_n_auxevol_mallocs:
            file_output_str += (
                """   // Step 4.a: Functions called after memory for non-y_n and auxevol gridfunctions is allocated.
"""
                + post_non_y_n_auxevol_mallocs
            )
        file_output_str += """
        send_wavespeed_at_outer_boundary(grid);
      }
      when receiv_wavespeed_at_outer_boundary(REAL wavespeed_at_outer_boundary) {
        serial {
          griddata_chare[grid].params.wavespeed_at_outer_boundary = wavespeed_at_outer_boundary;
        }
      }
"""

    for loop_direction in ["x", "y", "z"]:
        # Determine ghost types and configuration based on the current axis
        if loop_direction == "x":
            pos_ghost_type = x_pos_ghost_type
            neg_ghost_type = x_neg_ghost_type
            nchare_var = "commondata.Nchare0"
            grid_split_direction = "EAST_WEST"
        elif loop_direction == "y":
            pos_ghost_type = y_pos_ghost_type
            neg_ghost_type = y_neg_ghost_type
            nchare_var = "commondata.Nchare1"
            grid_split_direction = "NORTH_SOUTH"
        else:  # loop_direction == "z"
            pos_ghost_type = z_pos_ghost_type
            neg_ghost_type = z_neg_ghost_type
            nchare_var = "commondata.Nchare2"
            grid_split_direction = "TOP_BOTTOM"

        file_output_str += generate_send_neighbor_data_code(
            "Y_N_GFS", grid_split_direction
        )
        file_output_str += generate_process_ghost_code(
            loop_direction, pos_ghost_type, neg_ghost_type, nchare_var
        )

    file_output_str += r"""
      while (commondata.time < commondata.t_final) { // Main loop to progress forward in time.
        serial {
          time_start = commondata.time;
        }
        """
    if enable_residual_diagnostics:
        file_output_str += r"""
        serial {
          Ck::IO::Session token;  //pass a null token
          const int thisIndex_arr[3] = {thisIndex.x, thisIndex.y, thisIndex.z};
          REAL localsums_for_residualH[2];
          diagnostics(&commondata, griddata_chare, griddata, token, OUTPUT_RESIDUAL, which_grid_diagnostics, thisIndex_arr, localsums_for_residualH);
          contribute_localsums_for_residualH(localsums_for_residualH);
        }
        when continue_after_residual_H_done() {
        """
    if enable_residual_diagnostics:
        file_output_str += r"""
            serial {
              const int n_step = commondata.nn;
              const int outevery = commondata.diagnostics_output_every;
              write_diagnostics_this_step = n_step % outevery == 0;
            }
            """
    else:
        file_output_str += r"""
            serial {
              write_diagnostics_this_step = fabs(round(commondata.time / commondata.diagnostics_output_every) * commondata.diagnostics_output_every - commondata.time) < 0.5 * commondata.dt;
            }
            """
    if enable_psi4_diagnostics:
        file_output_str += r"""
        if (write_diagnostics_this_step) {
          serial {
            if(thisIndex.y == 0 && thisIndex.z == 0) {
              sectionBcastMsg *msg = new sectionBcastMsg(1);
              secProxy.recvMsg_to_contribute_localsums_for_psi4_decomp(msg);
            }
          }
        }
        """
    if enable_residual_diagnostics:
        filename_format = "commondata.nn"
    else:
        filename_format = "commondata.convergence_factor, commondata.time"

    file_output_str += rf"""
        // Step 5.a: Main loop, part 1: Output diagnostics
        //serial {{
        //  if (write_diagnostics_this_step && contains_gridcenter) {{
        //    diagnostics(&commondata, griddata_chare, Ck::IO::Session(), OUTPUT_0D, which_grid_diagnostics);
        //  }}
        //}}
        // Create sessions for ckio file writing from first chare only
        if (thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {{
          serial {{
            progress_indicator(&commondata, griddata_chare);
            if (commondata.time + commondata.dt > commondata.t_final)
              printf("\n");
          }}
          if (write_diagnostics_this_step) {{
            serial {{
              count_filewritten = 0;
              {{
                char filename[256];
                sprintf(filename, griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_1d_y, {filename_format});
                Ck::IO::Options opts;
                CkCallback opened_1d_y(CkIndex_Timestepping::ready_1d_y(NULL), thisProxy);
                Ck::IO::open(filename, opened_1d_y, opts);
              }}
              {{
                char filename[256];
                sprintf(filename, griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_1d_z, {filename_format});
                Ck::IO::Options opts;
                CkCallback opened_1d_z(CkIndex_Timestepping::ready_1d_z(NULL), thisProxy);
                Ck::IO::open(filename, opened_1d_z, opts);
              }}
              {{
                char filename[256];
                sprintf(filename, griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_2d_xy, {filename_format});
                Ck::IO::Options opts;
                CkCallback opened_2d_xy(CkIndex_Timestepping::ready_2d_xy(NULL), thisProxy);
                Ck::IO::open(filename, opened_2d_xy, opts);
              }}
              {{
                char filename[256];
                sprintf(filename, griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_2d_yz, {filename_format});
                Ck::IO::Options opts;
                CkCallback opened_2d_yz(CkIndex_Timestepping::ready_2d_yz(NULL), thisProxy);
                Ck::IO::open(filename, opened_2d_yz, opts);
              }}
            }}
"""
    # Generate code for 1d y diagnostics
    file_output_str += generate_diagnostics_code(
        "1d",
        "y",
        "griddata_chare[which_grid_diagnostics].diagnosticstruct.num_output_quantities + 1",
        "griddata_chare[which_grid_diagnostics].diagnosticstruct.tot_num_diagnostic_1d_y_pts",
    )

    # Generate code for 1d z diagnostics
    file_output_str += generate_diagnostics_code(
        "1d",
        "z",
        "griddata_chare[which_grid_diagnostics].diagnosticstruct.num_output_quantities + 1",
        "griddata_chare[which_grid_diagnostics].diagnosticstruct.tot_num_diagnostic_1d_z_pts",
    )

    # Generate code for 2d xy diagnostics
    file_output_str += generate_diagnostics_code(
        "2d",
        "xy",
        "griddata_chare[which_grid_diagnostics].diagnosticstruct.num_output_quantities + 2",
        "griddata_chare[which_grid_diagnostics].diagnosticstruct.tot_num_diagnostic_2d_xy_pts",
    )

    # Generate code for 2d yz diagnostics
    file_output_str += generate_diagnostics_code(
        "2d",
        "yz",
        "griddata_chare[which_grid_diagnostics].diagnosticstruct.num_output_quantities + 2",
        "griddata_chare[which_grid_diagnostics].diagnosticstruct.tot_num_diagnostic_2d_yz_pts",
    )

    file_output_str += r"""
            if (count_filewritten == expected_count_filewritten) {
              serial {thisProxy.continue_timestepping(); }
            }
          } else {
            serial {thisProxy.continue_timestepping(); }
          }
        }
        when continue_timestepping() {
"""
    if pre_MoL_step_forward_in_time != "":
        file_output_str += pre_MoL_step_forward_in_time
    else:
        file_output_str += "  // (nothing here; specify by setting pre_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"

    file_output_str += r"""
    """

    num_steps = len(Butcher_dict[MoL_method][0]) - 1

    # Loop over RK substeps and loop directions.
    for s in range(num_steps):
        rhs_output_exprs_list = generate_rhs_output_exprs(
            Butcher_dict, MoL_method, s + 1
        )
        post_rhs_output_list = generate_post_rhs_output_list(
            Butcher_dict, MoL_method, s + 1
        )
        file_output_str += generate_mol_step_forward_code(
            f"RK_SUBSTEP_K{s+1}", rhs_output_exprs_list, post_rhs_output_list
        )
        for loop_direction in ["x", "y", "z"]:
            # Determine ghost types and configuration based on the current axis
            if loop_direction == "x":
                pos_ghost_type = x_pos_ghost_type
                neg_ghost_type = x_neg_ghost_type
                nchare_var = "commondata.Nchare0"
                grid_split_direction = "EAST_WEST"
            elif loop_direction == "y":
                pos_ghost_type = y_pos_ghost_type
                neg_ghost_type = y_neg_ghost_type
                nchare_var = "commondata.Nchare1"
                grid_split_direction = "NORTH_SOUTH"
            else:  # loop_direction == "z"
                pos_ghost_type = z_pos_ghost_type
                neg_ghost_type = z_neg_ghost_type
                nchare_var = "commondata.Nchare2"
                grid_split_direction = "TOP_BOTTOM"

            # Loop over each element in post_rhs_output_list
            for post_rhs_output in post_rhs_output_list:
                file_output_str += generate_send_neighbor_data_code(
                    post_rhs_output, grid_split_direction
                )
                file_output_str += generate_process_ghost_code(
                    loop_direction, pos_ghost_type, neg_ghost_type, nchare_var
                )

    file_output_str += r"""
        """

    if post_MoL_step_forward_in_time != "":
        file_output_str += post_MoL_step_forward_in_time
    else:
        file_output_str += "  // (nothing here; specify by setting post_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"

    if enable_residual_diagnostics:
        file_output_str += r"""
      }
        """

    file_output_str += r"""
        }
    """

    file_output_str += r"""
        serial {
          // Adding dt to commondata.time many times will induce roundoff error,
          //   so here we set time based on the iteration number.
          commondata.time = (REAL)(commondata.nn + 1) * commondata.dt;
          // Finally, increment the timestep n:
          commondata.nn++;
        }
      } // End main loop to progress forward in time.
"""
    file_output_str += r"""
      serial{
        mainProxy.done();
      }
    };
    entry void start_write_1d_y(Ck::IO::SessionReadyMsg *m);
    entry void start_write_1d_z(Ck::IO::SessionReadyMsg *m);
    entry void start_write_2d_xy(Ck::IO::SessionReadyMsg *m);
    entry void start_write_2d_yz(Ck::IO::SessionReadyMsg *m);
    entry void test_written_1d_y(CkReductionMsg *m);
    entry void test_written_1d_z(CkReductionMsg *m);
    entry void test_written_2d_xy(CkReductionMsg *m);
    entry void test_written_2d_yz(CkReductionMsg *m);
    entry void closed_1d_y(CkReductionMsg *m);
    entry void closed_1d_z(CkReductionMsg *m);
    entry void closed_2d_xy(CkReductionMsg *m);
    entry void closed_2d_yz(CkReductionMsg *m);
    entry void diagnostics_ckio(Ck::IO::Session token, int which_output) {
      serial {
        const int thisIndex_arr[3] = {thisIndex.x, thisIndex.y, thisIndex.z};"""
    if enable_residual_diagnostics:
        file_output_str += r"""
        REAL unused_var[2];
        diagnostics(&commondata, griddata_chare, griddata, token, which_output, which_grid_diagnostics, thisIndex_arr, unused_var);"""
    else:
        file_output_str += r"""
        diagnostics(&commondata, griddata_chare, griddata, token, which_output, which_grid_diagnostics, thisIndex_arr);"""
    file_output_str += r"""
      }
    }
    entry void east_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void west_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void north_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void south_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void top_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void bottom_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void continue_timestepping();
    entry void receiv_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, int globalidx3_srcpts[num_srcpts]);
    entry void receiv_nonlocalinnerbc_data(int src_chare_idx3, int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);"""
    if enable_residual_diagnostics:
        file_output_str += r"""
    entry void continue_after_residual_H_done();
    entry void report_sums_for_residualH(CkReductionMsg *msg) {
      serial {
        int reducedArrSize=msg->getSize()/sizeof(double);
        CkAssert(reducedArrSize == 2);
        double *output=(double *)msg->getData();
        // Update residual to be used in stop condition
        commondata.log10_current_residual = log10(1e-16 + sqrt(output[0] / output[1])); // 1e-16 + ... avoids log10(0)
        // Output l2-norm of Hamiltonian constraint violation to file
        if (thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
          char filename[256];
          sprintf(filename, "residual_l2_norm.txt");
          const int nn = commondata.nn;
          const REAL time = commondata.time;
          const REAL residual_H =  commondata.log10_current_residual ;
          FILE *outfile = (nn == 0) ? fopen(filename, "w") : fopen(filename, "a");
          if (!outfile) {
            fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
            exit(1);
          }
          fprintf(outfile, "%6d %10.4e %.17e\n", nn, time, residual_H);
          fclose(outfile);
        }
        delete msg;
        thisProxy[CkArrayIndex3D(thisIndex.x, thisIndex.y, thisIndex.z)].continue_after_residual_H_done();
      }
    }
    entry void receiv_wavespeed_at_outer_boundary(REAL wavespeed_at_outer_boundary);"""
    if enable_psi4_diagnostics:
        file_output_str += r"""
    entry void recvMsg_to_contribute_localsums_for_psi4_decomp(sectionBcastMsg *msg){
      serial {
        Ck::IO::Session token; // pass a null token
        const int thisIndex_arr[3] = {thisIndex.x, thisIndex.y, thisIndex.z};
        diagnostics(&commondata, griddata_chare, griddata, token, OUTPUT_PSI4, which_grid_diagnostics, thisIndex_arr);
        contribute_localsums_for_psi4_decomp(msg, which_grid_diagnostics);
      }
    }
    entry void report_sums_for_psi4_diagnostics(CkReductionMsg * msg) {
      serial {
        int reducedArrSize = msg->getSize() / sizeof(double);
        double *output = (double *)msg->getData();
        const int length_localsums_for_psi4_decomp = griddata_chare[which_grid_diagnostics].diagnosticstruct.length_localsums_for_psi4_decomp;
        for (int i = 0; i < length_localsums_for_psi4_decomp; i++) {
          griddata_chare[which_grid_diagnostics].diagnosticstruct.globalsums_for_psi4_decomp[i] = (REAL)output[i];
        }
        psi4_spinweightm2_decomposition_file_write(&commondata, &griddata_chare[which_grid_diagnostics].diagnosticstruct);
        delete msg;
      }
    }"""

    file_output_str += r"""
  };
};
"""

    timestepping_ci_file = project_Path / "timestepping.ci"
    with timestepping_ci_file.open("w", encoding="utf-8") as file:
        file.write(
            clang_format(file_output_str, clang_format_options=clang_format_options)
        )


def output_timestepping_h_cpp_ci_register_CFunctions(
    project_dir: str,
    MoL_method: str = "RK4",
    enable_rfm_precompute: bool = False,
    post_non_y_n_auxevol_mallocs: str = "",
    pre_MoL_step_forward_in_time: str = "",
    post_MoL_step_forward_in_time: str = "",
    enable_psi4_diagnostics: bool = False,
    enable_residual_diagnostics: bool = False,
) -> None:
    """
    Output timestepping h, cpp, and ci files and register C functions.

    :param project_dir: Directory where the project C code is output
    :param MoL_method: Method of Lines (MoL) method name, default is "RK4".
    :param enable_rfm_precompute: Enable RFM precompute, default is False.
    :param post_non_y_n_auxevol_mallocs: Function calls after memory is allocated for non y_n and auxevol gridfunctions, default is an empty string.
    :param pre_MoL_step_forward_in_time: Code for handling pre-right-hand-side operations, default is an empty string.
    :param post_MoL_step_forward_in_time: Code for handling post-right-hand-side operations, default is an empty string.
    :param enable_psi4_diagnostics: Whether or not to enable psi4 diagnostics.
    :param enable_residual_diagnostics: Whether or not to enable residual diagnostics.
    """
    # For NRPy elliptic: register parameter wavespeed at outer boundary
    if enable_residual_diagnostics:
        _wavespeed_at_outer_boundary = par.register_CodeParameter(
            "REAL", __name__, "wavespeed_at_outer_boundary", 0.0, commondata=False
        )

    output_timestepping_h(
        project_dir=project_dir,
        enable_residual_diagnostics=enable_residual_diagnostics,
        enable_psi4_diagnostics=enable_psi4_diagnostics,
    )

    Butcher_dict = generate_Butcher_tables()

    output_timestepping_cpp(
        project_dir=project_dir,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_CurviBCs=True,
        Butcher_dict=Butcher_dict,
        MoL_method=MoL_method,
        enable_residual_diagnostics=enable_residual_diagnostics,
        enable_psi4_diagnostics=enable_psi4_diagnostics,
    )

    output_timestepping_ci(
        project_dir=project_dir,
        MoL_method=MoL_method,
        post_non_y_n_auxevol_mallocs=post_non_y_n_auxevol_mallocs,
        pre_MoL_step_forward_in_time=pre_MoL_step_forward_in_time,
        post_MoL_step_forward_in_time=post_MoL_step_forward_in_time,
        enable_psi4_diagnostics=enable_psi4_diagnostics,
        enable_residual_diagnostics=enable_residual_diagnostics,
        Butcher_dict=Butcher_dict,
    )

    register_CFunction_timestepping_malloc()

    register_CFunction_timestepping_free_memory()

    # Register temporary buffers for face data communication to griddata_struct:
    griddata_commondata.register_griddata_commondata(
        __name__,
        "tmpBuffers_struct tmpBuffers",
        "temporary buffer for sending face data to neighbor chares",
    )
