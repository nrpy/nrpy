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
from nrpy.infrastructures import BHaH, superB


def generate_complete_gf_list(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
) -> List[str]:
    """
    Generate the complete list of grid function type names based on the given Method of Lines (MoL) method.

    :param Butcher_dict: Dictionary containing Butcher tableau data.
    :param MoL_method: Method of Lines (MoL) method name.
    :return: A list of grid function type names (intermediate-stage names plus
             'y_n_gfs', 'auxevol_gfs_gfs', and 'diagnostic_output_gfs').
    """
    gf_list = BHaH.MoLtimestepping.rk_butcher_table_dictionary.intermediate_stage_gf_names_list(
        Butcher_dict, MoL_method=MoL_method
    )
    gf_list.extend(["y_n_gfs", "auxevol_gfs", "diagnostic_output_gfs"])
    return gf_list


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


def generate_process_nonlocalinnerbc_code(which_gf: str) -> str:
    """
    Generate code for nonlocal inner bc processing.

    :param which_gf: Which grid function
    :return: Code for nonlocal inner bc processing.
    """
    return f"""
      if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares > 0) {{
        for (iter = 0; iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares; iter++) {{
          when receiv_nonlocalinnerbc_data_{which_gf.lower()}(int src_chare_idx3, int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]) {{
            serial {{
              set_tmpBuffer_innerbc_receiv(src_chare_idx3, len_tmpBuffer, tmpBuffer, grid);
              type_gfs_nonlocal_innerbc = type_gfs;
            }}
          }}
        }}
        serial {{
          process_nonlocalinnerbc(type_gfs_nonlocal_innerbc, grid);
        }}
      }}
"""


def generate_mol_step_forward_code(
    rk_substep: str,
    rhs_output_exprs_list: List[str],
    post_rhs_output_list: List[str],
    outer_bcs_type: str = "radiation",
) -> str:
    """
    Generate code for MoL step forward in time.

    :param rk_substep: Runge-Kutta substep.
    :param rhs_output_exprs_list: List of output expression for the RHS.
    :param post_rhs_output_list: List of outputs for post-RHS expressions.
    :param outer_bcs_type: type of outer boundary BCs to apply. Only options are radiation or extrapolation in superB.
    :return: Code for MoL step forward in time.
    """
    return_str = f"""
    serial{{
        MoL_step_forward_in_time(&commondata, griddata_chare, time_start, {rk_substep},  MOL_PRE_RK_UPDATE);
    }}
"""
    if outer_bcs_type == "radiation":
        for rhs_output_exprs in rhs_output_exprs_list:
            return_str += generate_send_nonlocalinnerbc_data_code(rhs_output_exprs)
            return_str += generate_process_nonlocalinnerbc_code(rhs_output_exprs)

    return_str += f"""
    serial{{
        MoL_step_forward_in_time(&commondata, griddata_chare, time_start, {rk_substep}, MOL_RK_UPDATE);
    }}
"""
    if outer_bcs_type == "extrapolation":
        return_str += f"""
    serial{{
        MoL_step_forward_in_time(&commondata, griddata_chare, time_start, {rk_substep}, MOL_POST_RK_UPDATE_APPLY_BCS);
    }}
"""
        for post_rhs_output in post_rhs_output_list:
            return_str += generate_send_nonlocalinnerbc_data_code(post_rhs_output)
            return_str += generate_process_nonlocalinnerbc_code(post_rhs_output)

    return_str += f"""
    serial{{
        MoL_step_forward_in_time(&commondata, griddata_chare, time_start, {rk_substep}, MOL_POST_RK_UPDATE);
    }}
"""
    return_str += """
    if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {"""
    return_str += generate_send_nonlocalinnerbc_data_code("AUXEVOL_GFS")
    return_str += generate_process_nonlocalinnerbc_code("AUXEVOL_GFS")
    return_str += """}
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
    dimension: str,
    direction: str,
    totsize_field: str,
    which_output: str,
) -> str:
    """
    Generate code for diagnostics.

    :param dimension: Dimension index
    :param direction: Direction of diagnostics (e.g. "x", "y", "z")
    :param totsize_field: Name of the diagnosticstruct field storing total output size in bytes.
    :param which_output: Constant selecting which diagnostics output to write.
    :return: Code for diagnostics
    """
    code = f"""
    // {dimension} {direction}
    when ready_{dimension}_{direction}(Ck::IO::FileReadyMsg *m_{dimension}_{direction}) {{
        serial {{
        f_{dimension}_{direction} = m_{dimension}_{direction}->file;
        CkCallback sessionStart_{dimension}_{direction}(CkIndex_Timestepping::start_write_{dimension}_{direction}(0), thisProxy);
        CkCallback sessionEnd_{dimension}_{direction}(CkIndex_Timestepping::test_written_{dimension}_{direction}(0), thisProxy);
        int totsizeinbytes = griddata_chare[which_grid_diagnostics].diagnosticstruct.{totsize_field};
        Ck::IO::startSession(f_{dimension}_{direction}, totsizeinbytes, 0, sessionStart_{dimension}_{direction}, sessionEnd_{dimension}_{direction});
        delete m_{dimension}_{direction};
        }}
    }}
    when start_write_{dimension}_{direction}(Ck::IO::SessionReadyMsg *m_{dimension}_{direction}){{
        serial {{
        thisProxy.diagnostics_ckio(m_{dimension}_{direction}->session, {which_output});
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


def generate_PUP_code() -> str:
    """
    Generate code for PUP routine for Timestepping class.

    :return: A string representing the PUP routine for Timestepping class.
    """
    code = r"""
// PUP routine for class Timestepping
// Note that the user can choose what to PUP, for example not all structs in griddata are PUP'ed
// gridddata has to be unpacked first since it is needed to unpack griddata_chare(for the size of certain arrays)
void Timestepping::pup(PUP::er &p) {
  CBase_Timestepping::pup(p);
  __sdag_pup(p);
  pup_commondata_struct(p, commondata);
  if (p.isUnpacking()) {
    griddata = (griddata_struct *restrict)malloc(sizeof(griddata_struct) * commondata.NUMGRIDS);
  }
  for (int i = 0; i < commondata.NUMGRIDS; i++) {
    pup_griddata(p, griddata[i]);
  }
  if (p.isUnpacking()) {
    griddata_chare = (griddata_struct *restrict)malloc(sizeof(griddata_struct) * commondata.NUMGRIDS);
  }
  for (int i = 0; i < commondata.NUMGRIDS; i++) {
    pup_griddata_chare(p, griddata_chare[i], griddata[i].params, commondata);
  }
  p | is_boundarychare;
  p | time_start;
  p | count_filewritten;
  p | iter;
  p | type_gfs_nonlocal_innerbc;
  p | write_diagnostics_this_step;
  p | f_1d_y;
  p | f_1d_z;
  p | f_2d_xy;
  p | f_2d_yz;
  p | const_cast<int&>(grid);
  p | const_cast<int&>(which_grid_diagnostics);
  p | const_cast<int&>(expected_count_filewritten);
}
"""
    return code


def register_CFunction_timestepping_malloc(
    enable_BHaHAHA: bool = False,
) -> None:
    """
    Register a C function for timestepping malloc.

    :param enable_BHaHAHA: If True, allocate memory of tmpBuffer_bhahaha_gfs..
    :return None
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Allocate memory for temporary buffers used to communicate face data"
    cfunc_type = "void"
    name = "timestepping_malloc_tmpBuffer"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const MoL_gridfunctions_struct *restrict gridfuncs, const nonlocalinnerbc_struct *restrict nonlocalinnerbcstruct, tmpBuffers_struct *restrict tmpBuffers"
    body = """
const int Nxx_plus_2NGHOSTS_face0 = Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
const int Nxx_plus_2NGHOSTS_face1 = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS2;
const int Nxx_plus_2NGHOSTS_face2 = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1;

const int max_sync_gfs = gridfuncs->max_sync_gfs;

tmpBuffers->tmpBuffer_EW = (REAL *restrict)malloc(sizeof(REAL) * max_sync_gfs * NGHOSTS * Nxx_plus_2NGHOSTS_face0);
tmpBuffers->tmpBuffer_NS = (REAL *restrict)malloc(sizeof(REAL) * max_sync_gfs * NGHOSTS * Nxx_plus_2NGHOSTS_face1);
tmpBuffers->tmpBuffer_TB = (REAL *restrict)malloc(sizeof(REAL) * max_sync_gfs * NGHOSTS * Nxx_plus_2NGHOSTS_face2);

// Unpack nonlocalinnerbcstruct
const int tot_num_dst_chares = nonlocalinnerbcstruct->tot_num_dst_chares;
const int *restrict num_srcpts_tosend_each_chare = nonlocalinnerbcstruct->num_srcpts_tosend_each_chare;
const int tot_num_src_chares = nonlocalinnerbcstruct->tot_num_src_chares;
const int *restrict num_srcpts_each_chare = nonlocalinnerbcstruct->num_srcpts_each_chare;

tmpBuffers->tmpBuffer_innerbc_send = (REAL **)malloc(tot_num_dst_chares * sizeof(REAL *));
for (int which_chare = 0; which_chare < tot_num_dst_chares; which_chare++) {
  tmpBuffers->tmpBuffer_innerbc_send[which_chare] = (REAL *restrict)malloc(sizeof(REAL) * max_sync_gfs * num_srcpts_tosend_each_chare[which_chare]);
}
tmpBuffers->tmpBuffer_innerbc_receiv = (REAL **)malloc(tot_num_src_chares * sizeof(REAL *));
for (int which_chare = 0; which_chare < tot_num_src_chares; which_chare++) {
  tmpBuffers->tmpBuffer_innerbc_receiv[which_chare] = (REAL *restrict)malloc(sizeof(REAL) * max_sync_gfs * num_srcpts_each_chare[which_chare]);
}
"""
    if enable_BHaHAHA:
        body += """
const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
tmpBuffers->tmpBuffer_bhahaha_gfs = (REAL *restrict)malloc(sizeof(REAL) * BHAHAHA_NUM_INTERP_GFS * Nxx_plus_2NGHOSTS_tot);
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


def register_CFunction_timestepping_free_memory(
    enable_BHaHAHA: bool = False,
) -> None:
    """
    Register a C function for timestepping free memory.

    :param enable_BHaHAHA: If True, free memory of tmpBuffer_bhahaha_gfs.
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
    if enable_BHaHAHA:
        body += """
if (tmpBuffers->tmpBuffer_bhahaha_gfs != NULL)
  free(tmpBuffers->tmpBuffer_bhahaha_gfs);
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
    nrpyelliptic_project: bool = False,
    enable_psi4: bool = False,
    enable_charm_checkpointing: bool = False,
    enable_BHaHAHA: bool = False,
) -> None:
    """
    Generate timestepping.h.

    :param project_dir: Directory where the project C code is output.
    :param nrpyelliptic_project: If True, enable NRPyElliptic project mode (enables residual
                    diagnostics and NRPyElliptic-specific behavior).
    :param enable_psi4: Whether or not to enable psi4 diagnostics.
    :param enable_charm_checkpointing: Enable checkpointing using Charm++.
    :param enable_BHaHAHA: If True, add creation of horizon_finder and interpolator3d chares and communication with them.
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = r"""#ifndef __TIMESTEPPING_H__
#define __TIMESTEPPING_H__
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "superB/superB_pup_function_prototypes.h"
#include "timestepping.decl.h"
#include "diagnostics/diagnostic_gfs.h"
#include "diagnostics/diagnostics_nearest_common.h"
#include "diagnostics/diagnostics_volume_integration_helpers.h"
"""
    file_output_str += r"""
/**
 * @brief Build the default set of recipes used by diagnostics_volume_integration and reductions.
 *
 * @param[out] recipes Array of length DIAGS_INTEGRATION_MAX_RECIPES to populate.
 * @return Number of recipes written into @p recipes.
 *
 * Note: Keep this function self-contained so callers in different translation units
 * can construct the same canonical recipe book without duplicating code.
 */
static inline int diags_integration_build_default_recipes(diags_integration_recipe_t recipes[DIAGS_INTEGRATION_MAX_RECIPES]) {
  diags_integration_initialize_recipes(recipes);
  int NUM_RECIPES = 0;

  if (NUM_RECIPES < DIAGS_INTEGRATION_MAX_RECIPES) {"""

    if nrpyelliptic_project:
        file_output_str += r"""
    recipes[NUM_RECIPES].name = "whole_domain";
    recipes[NUM_RECIPES].num_rules = 0; // No rules means the whole domain is used.

    // Important: is_squared=1 enables computation of L2 norm & RMS; RMS_f = sqrt(int f^2 dV / int dV)
    recipes[NUM_RECIPES].integrands[0] = (diags_integration_integrand_spec_t){.gf_index = DIAG_RESIDUALGF, .is_squared = 1};
    recipes[NUM_RECIPES].num_integrands = 1;"""

    else:
        file_output_str += r"""
    recipes[NUM_RECIPES].name = "whole_domain";
    recipes[NUM_RECIPES].num_rules = 0; // No rules means the whole domain is used.

    // Important: is_squared=1 enables computation of L2 norm & RMS; RMS_f = sqrt(int f^2 dV / int dV)
    recipes[NUM_RECIPES].integrands[0] = (diags_integration_integrand_spec_t){.gf_index = DIAG_HAMILTONIANGF, .is_squared = 1};
    recipes[NUM_RECIPES].integrands[1] = (diags_integration_integrand_spec_t){.gf_index = DIAG_MSQUAREDGF, .is_squared = 1};
    recipes[NUM_RECIPES].num_integrands = 2;"""

    file_output_str += r"""
    NUM_RECIPES++;
  }

  if (NUM_RECIPES < DIAGS_INTEGRATION_MAX_RECIPES) {"""

    if nrpyelliptic_project:
        file_output_str += r"""
    const REAL R_outer = 80;
    recipes[NUM_RECIPES].name = "sphere_R_80";
    recipes[NUM_RECIPES].num_rules = 1;

    // Important: exclude_inside=0 implies outer is excluded.
    recipes[NUM_RECIPES].rules[0] = (diags_integration_sphere_rule_t){.center_xyz = {0, 0, 0}, .radius = R_outer, .exclude_inside = 0};

    // Important: is_squared=1 enables computation of L2 norm & RMS; RMS_f = sqrt(int f^2 dV / int dV)
    recipes[NUM_RECIPES].integrands[0] = (diags_integration_integrand_spec_t){.gf_index = DIAG_RESIDUALGF, .is_squared = 1};
    recipes[NUM_RECIPES].num_integrands = 1;"""

    else:
        file_output_str += r"""
    const REAL R_outer = 8;
    recipes[NUM_RECIPES].name = "sphere_R_gt_8";
    recipes[NUM_RECIPES].num_rules = 1;

    // Important: exclude_inside=0 implies outer is excluded.
    recipes[NUM_RECIPES].rules[0] = (diags_integration_sphere_rule_t){.center_xyz = {0, 0, 0}, .radius = R_outer, .exclude_inside = 1};

    // Important: is_squared=1 enables computation of L2 norm & RMS; RMS_f = sqrt(int f^2 dV / int dV)
    recipes[NUM_RECIPES].integrands[0] = (diags_integration_integrand_spec_t){.gf_index = DIAG_HAMILTONIANGF, .is_squared = 1};
    recipes[NUM_RECIPES].integrands[1] = (diags_integration_integrand_spec_t){.gf_index = DIAG_MSQUAREDGF, .is_squared = 1};
    recipes[NUM_RECIPES].num_integrands = 2;"""

    file_output_str += r"""
    NUM_RECIPES++;
  }
  return NUM_RECIPES;
}
"""

    file_output_str += r"""
class Timestepping : public CBase_Timestepping {
  Timestepping_SDAG_CODE

  private:
    /// Member Variables (Object State) ///"""
    file_output_str += r"""
    commondata_struct commondata;
    griddata_struct *griddata;
    griddata_struct *griddata_chare;
    bool is_boundarychare;
    REAL time_start;
    //bool contains_gridcenter;
    const int grid = 0;
    REAL *diagnostic_gfs[MAXNUMGRIDS];
    const int which_grid_diagnostics = 0;
    const bool free_non_y_n_gfs_and_core_griddata_pointers=true;
    bool write_diagnostics_this_step;"""
    if enable_charm_checkpointing:
        file_output_str += r"""
    bool write_chckpt_this_step;"""
    if enable_BHaHAHA:
        file_output_str += r"""
    bool do_horizon_find;"""
    file_output_str += r"""
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
    if nrpyelliptic_project:
        file_output_str += r"""
    void send_wavespeed_at_outer_boundary(const int grid);"""
    file_output_str += r"""
    void contribute_localsums_for_diagnostic_volume_integ();"""

    if enable_BHaHAHA:
        file_output_str += r"""
    void send_interp_gfs_to_corresponding_interpolator_chare(const int grid, const REAL *restrict gfs, const int *restrict gf_indices,
                                                             const int num_gfs, const int request_type);
    void send_bhahaha_gfs_to_corresponding_interpolator_chare(const int grid);
"""
        if enable_psi4:
            file_output_str += r"""
    void send_psi4_gfs_to_corresponding_interpolator_chare(const int grid);
"""
    file_output_str += r"""
  public:
    /// Constructors ///
    Timestepping(CommondataObject &&inData);
    Timestepping(CkMigrateMessage* msg);
    /// Destructor ///
    ~Timestepping();"""
    if enable_charm_checkpointing:
        file_output_str += r"""
    void pup(PUP::er &p);"""
    file_output_str += r"""
};

#endif //__TIMESTEPPING_H__
"""
    timestepping_h_file = project_Path / "timestepping.h"
    with timestepping_h_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def generate_switch_statement_for_gf_types(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
    set_parity_types: bool = False,
) -> str:
    """
    Generate the switch statement for grid function types based on the given Method of Lines (MoL) method.

    :param Butcher_dict: Dictionary containing Butcher tableau data.
    :param MoL_method: Method of Lines (MoL) method name.
    :param set_parity_types: whether to set parity types in switch statements.
    :return: A string representing the switch statement for the grid function types.
    """
    gf_list = generate_complete_gf_list(Butcher_dict, MoL_method=MoL_method)

    switch_statement = """
switch (type_gfs) {
"""
    switch_cases = []
    for gf in gf_list:
        switch_cases.append(f"  case {gf.upper()}:")
        if gf == "y_n_gfs":
            switch_cases.append("  case Y_N_GFS_INITIALDATA_PART1:")
            switch_cases.append("  case Y_N_GFS_INITIALDATA_PART2:")
        if gf not in ("auxevol_gfs", "diagnostic_output_gfs"):
            switch_cases.append(
                f"    gfs = griddata_chare[grid].gridfuncs.{gf.lower()};"
            )
            switch_cases.append(
                "    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;"
            )
            switch_cases.append(
                "    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;"
            )
            if set_parity_types:
                switch_cases.append("    gf_parity_types = evol_gf_parity;")
        elif gf == "auxevol_gfs":
            switch_cases.append(
                f"    gfs = griddata_chare[grid].gridfuncs.{gf.lower()};"
            )
            switch_cases.append(
                "    NUM_GFS = griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync;"
            )
            switch_cases.append(
                "    gfs_to_sync = griddata_chare[grid].gridfuncs.auxevol_gfs_to_sync;"
            )
            if set_parity_types:
                switch_cases.append("    gf_parity_types = auxevol_gf_parity;")
        elif gf == "diagnostic_output_gfs":
            switch_cases.append(
                f"    gfs = griddata_chare[grid].gridfuncs.{gf.lower()};"
            )
            switch_cases.append(
                "    NUM_GFS = griddata_chare[grid].gridfuncs.num_aux_gfs_to_sync;"
            )
            switch_cases.append(
                "    gfs_to_sync = griddata_chare[grid].gridfuncs.aux_gfs_to_sync;"
            )
            if set_parity_types:
                switch_cases.append("    gf_parity_types = aux_gf_parity;")
        switch_cases.append("    break;")
    switch_cases.append("""
  default:
    break;
}
""")

    switch_body = "\n".join(switch_cases)

    return switch_statement + switch_body


def generate_switch_statement_for_gf_types_for_entry_method(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
) -> str:
    """
    Generate the switch statement for grid function types based on the given Method of Lines (MoL) method for calling entry method based on gf type.

    :param Butcher_dict: Dictionary containing Butcher tableau data.
    :param MoL_method: Method of Lines (MoL) method name.
    :return: A string representing the switch statement for the grid function types.
    """
    gf_list = generate_complete_gf_list(Butcher_dict, MoL_method=MoL_method)

    # Keep synching of y n gfs during initial data distinct to prevent mismatch of messages
    gf_list.append("y_n_gfs_initialdata_part1")
    gf_list.append("y_n_gfs_initialdata_part2")

    switch_statement = """
switch (type_gfs) {
"""
    switch_cases = []
    for gf in gf_list:
        switch_cases.append(f"  case {gf.upper()}:")
        switch_cases.append(
            f"    thisProxy[CkArrayIndex3D(dst_chare_index0, dst_chare_index1, dst_chare_index2)].receiv_nonlocalinnerbc_data_{gf.lower()}(idx3_this_chare, type_gfs, NUM_GFS * num_srcpts_tosend_each_chare[which_dst_chare], tmpBuffer_innerbc_send);"
        )
        switch_cases.append("    break;")
    switch_cases.append("""
  default:
    break;
}
""")
    switch_body = "\n".join(switch_cases)
    return switch_statement + switch_body


def generate_entry_methods_for_receiv_nonlocalinnerbc_for_gf_types(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
    outer_bcs_type: str = "radiation",
    enable_psi4: bool = False,
    nrpyelliptic_project: bool = False,
) -> str:
    """
    Generate entry method declarations based on grid function types.

    :param Butcher_dict: Dictionary containing Butcher tableau data.
    :param MoL_method: Method of Lines (MoL) method name.
    :param outer_bcs_type: type of outer boundary BCs to apply. Only options are radiation or extrapolation in superB.
    :param enable_psi4: Whether to enable psi4 diagnostics.
    :param nrpyelliptic_project: If True, enable NRPyElliptic project mode (enables residual
                    diagnostics and NRPyElliptic-specific behavior).
    :return: A string containing entry method declarations separated by newlines.
    :raises ValueError: If `outer_bcs_type` is not set to either 'radiation' or 'extrapolation'.
    """
    gf_list = generate_complete_gf_list(Butcher_dict, MoL_method=MoL_method)

    # need separate entry methods of y n gf during initial data set up to prevent mismatch of messages
    gf_list.append("Y_N_GFS_INITIALDATA_PART1")
    gf_list.append("Y_N_GFS_INITIALDATA_PART2")

    # Find the list of all gfs to which inner bc synching is performed in the .ci file
    # for gf not in the list above, finish entry method by "{}" instead of ";"
    rhs_output_exprs_list_all = []
    post_rhs_output_list_all = []
    num_steps = len(Butcher_dict[MoL_method][0]) - 1
    for s in range(num_steps):
        rhs_output_exprs_list = superB.MoL.generate_rhs_output_exprs(
            Butcher_dict, MoL_method, s + 1
        )
        post_rhs_output_list = superB.MoL.generate_post_rhs_output_list(
            Butcher_dict, MoL_method, s + 1
        )
        rhs_output_exprs_list_all.extend(rhs_output_exprs_list)
        post_rhs_output_list_all.extend(post_rhs_output_list)
    if outer_bcs_type == "radiation":
        inner_bc_synching_gfs = rhs_output_exprs_list_all
    elif outer_bcs_type == "extrapolation":
        inner_bc_synching_gfs = post_rhs_output_list_all
    else:
        raise ValueError(
            f"Invalid value for 'outer_bcs_type': '{outer_bcs_type}'. "
            "Expected 'radiation' or 'extrapolation'. Please provide a valid boundary condition type."
        )

    inner_bc_synching_gfs.append("AUXEVOL_GFS")

    if enable_psi4:
        inner_bc_synching_gfs.append("DIAGNOSTIC_OUTPUT_GFS")

    # If anything other than NRPy elliptic, in NRPy elliptic initial data is set up differently
    if not nrpyelliptic_project:
        inner_bc_synching_gfs.append("Y_N_GFS_INITIALDATA_PART1")
        inner_bc_synching_gfs.append("Y_N_GFS_INITIALDATA_PART2")

    entry_method_for_gf_types: List[str] = []
    for gf in gf_list:
        if gf.upper() in inner_bc_synching_gfs:
            entry_method = f"entry void receiv_nonlocalinnerbc_data_{gf.lower()}(int src_chare_idx3, int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);"
        else:
            entry_method = f"entry void receiv_nonlocalinnerbc_data_{gf.lower()}(int src_chare_idx3, int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]){{}}"
        entry_method_for_gf_types.append(entry_method)

    # Return the concatenated string
    return "\n".join(entry_method_for_gf_types)


def output_timestepping_cpp(
    project_dir: str,
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
    initial_data_desc: str = "",
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
    initialize_constant_auxevol: bool = False,
    nrpyelliptic_project: bool = False,
    enable_psi4: bool = False,
    enable_charm_checkpointing: bool = False,
    enable_BHaHAHA: bool = False,
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
    :param nrpyelliptic_project: If True, enable NRPyElliptic project mode (enables residual
                    diagnostics and NRPyElliptic-specific behavior).
    :param enable_psi4: Whether or not to enable psi4 diagnostics.
    :param enable_charm_checkpointing: Enable checkpointing using Charm++.
    :param enable_BHaHAHA: If True, add creation of horizon_finder and interpolator3d chares and communication with them.
    """
    initial_data_desc += " "

    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = r"""#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "timestepping.h"
#include "main.h"
#include "diagnostics/diagnostics_volume_integration_helpers.h"
"""

    if enable_BHaHAHA:
        file_output_str += r"""
#include "horizon_finder.h"
"""

    file_output_str += r"""
extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Timestepping timesteppingArray;"""

    if enable_BHaHAHA:
        file_output_str += r"""
extern /* readonly */ CProxy_Horizon_finder horizon_finderProxy;
extern /* readonly */ CProxy_Interpolator3d interpolator3dArray;"""

    file_output_str += r"""
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

  CkPrintf("Timestepping chare %d,%d,%d created on PE %d\n", thisIndex.x, thisIndex.y, thisIndex.z, CkMyPe());

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
    // Define data needed for syncing gfs across chares
    MoL_sync_data_defines(&griddata_chare[grid].gridfuncs);
  }
"""

    file_output_str += """
//Allocate storage for gridfunctions on each grid.
for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
  const int Nxx_plus_2NGHOSTS_tot = (griddata[grid].params.Nxx_plus_2NGHOSTS0 * //
                                     griddata[grid].params.Nxx_plus_2NGHOSTS1 * //
                                     griddata[grid].params.Nxx_plus_2NGHOSTS2);

  BHAH_MALLOC(griddata_chare[grid].gridfuncs.y_n_gfs, sizeof(REAL) * Nxx_plus_2NGHOSTS_tot * NUM_EVOL_GFS);

  MoL_malloc_intermediate_stage_gfs(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].gridfuncs);

  if (NUM_AUXEVOL_GFS > 0) {
    BHAH_MALLOC(griddata_chare[grid].gridfuncs.auxevol_gfs, sizeof(REAL) * Nxx_plus_2NGHOSTS_tot * NUM_AUXEVOL_GFS);
  } // END IF NUM_AUXEVOL_GFS > 0

} // END LOOP over grids
"""

    file_output_str += """
  // Initialize y n and non-y n gfs to NAN to avoid uninitialized memory errors
  for(int grid=0; grid<commondata.NUMGRIDS; grid++)
    initialize_yn_and_non_yn_gfs_to_nan(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].gridfuncs);
"""

    file_output_str += """
  // Allocate storage for temporary buffers, needed for communicating face data
  for(int grid=0; grid<commondata.NUMGRIDS; grid++)
    timestepping_malloc_tmpBuffer(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].gridfuncs, &griddata_chare[grid].nonlocalinnerbcstruct, &griddata_chare[grid].tmpBuffers);
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
Timestepping::Timestepping(CkMigrateMessage *msg): CBase_Timestepping(msg) { }


// destructor
Timestepping::~Timestepping() {

  griddata_free(&commondata, griddata, free_non_y_n_gfs_and_core_griddata_pointers);

  for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
    timestepping_free_memory_tmpBuffer(&griddata_chare[grid].nonlocalinnerbcstruct, &griddata_chare[grid].tmpBuffers);"""
    if enable_CurviBCs:
        file_output_str += r"""
    free(griddata_chare[grid].bcstruct.inner_bc_array);
    free(griddata_chare[grid].bcstruct.inner_bc_array_nonlocal);
    for(int ng=0;ng<NGHOSTS*3;ng++) {
     free(griddata_chare[grid].bcstruct.pure_outer_bc_array[ng]);
    }"""
    file_output_str += r"""
    for(int i=0;i<3;i++) {
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
    free(griddata_chare[grid].diagnosticstruct.offset_diagnostic_2d_yz_pt);"""

    file_output_str += r"""
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
  int NUM_GFS;
  const int* gfs_to_sync = nullptr;
  """
    file_output_str += switch_case_code
    file_output_str += r"""
  switch (dir) {
    case EAST_WEST:
      //send to west
      if (thisIndex.x > 0) {
        for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
          int i0 = 2*NGHOSTS - 1;
          for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
            for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
              for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
                tmpBuffer_EW[IDXFACES0(which_gf, which_inner, i1, i2)] = gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)];
              }
            }
            i0--;
          }
        }
        thisProxy[CkArrayIndex3D(thisIndex.x - 1, thisIndex.y, thisIndex.z)].east_ghost(type_gfs, NUM_GFS*NGHOSTS*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2, tmpBuffer_EW);
      }
      //send to east
      if (thisIndex.x < Nchare0 - 1) {
        for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
          int i0 = Nxx0;
          for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
            for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
              for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
                tmpBuffer_EW[IDXFACES0(which_gf, which_inner, i1, i2)] = gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)];
              }
            }
            i0++;
          }
        }
        thisProxy[CkArrayIndex3D(thisIndex.x + 1, thisIndex.y, thisIndex.z)].west_ghost(type_gfs, NUM_GFS*NGHOSTS*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2, tmpBuffer_EW);
      }
      break;
    case NORTH_SOUTH:
      //send to south
      if (thisIndex.y > 0) {
        for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
          int i1 = 2*NGHOSTS - 1;
          for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
            for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
              for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
                tmpBuffer_NS[IDXFACES1(which_gf, which_inner, i0, i2)] = gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)];
              }
            }
            i1--;
          }
        }
        thisProxy[CkArrayIndex3D(thisIndex.x, thisIndex.y - 1, thisIndex.z)].north_ghost(type_gfs, NUM_GFS*NGHOSTS*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS2, tmpBuffer_NS);
      }
      //send to north
      if (thisIndex.y < Nchare1 - 1) {
        for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
          int i1 = Nxx1;
          for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
            for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
              for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
                tmpBuffer_NS[IDXFACES1(which_gf, which_inner, i0, i2)] = gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)];
              }
            }
            i1++;
          }
        }
        thisProxy[CkArrayIndex3D(thisIndex.x, thisIndex.y + 1, thisIndex.z)].south_ghost(type_gfs, NUM_GFS*NGHOSTS*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS2, tmpBuffer_NS);
      }
      break;
    case TOP_BOTTOM:
      //send to bottom
      if (thisIndex.z > 0) {
        for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
          int i2 = 2*NGHOSTS - 1;
          for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
            for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
              for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
                tmpBuffer_TB[IDXFACES2(which_gf, which_inner, i0, i1)] = gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)];
              }
            }
            i2--;
          }
        }
        thisProxy[CkArrayIndex3D(thisIndex.x, thisIndex.y, thisIndex.z - 1)].top_ghost(type_gfs, NUM_GFS*NGHOSTS*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1, tmpBuffer_TB);
      }
      //send to top
      if (thisIndex.z < Nchare2 - 1) {
        for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
          int i2 = Nxx2;
          for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
            for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
              for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
                tmpBuffer_TB[IDXFACES2(which_gf, which_inner, i0, i1)] = gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)];
              }
            }
            i2++;
          }
        }
        thisProxy[CkArrayIndex3D(thisIndex.x, thisIndex.y, thisIndex.z + 1)].bottom_ghost(type_gfs, NUM_GFS*NGHOSTS*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1, tmpBuffer_TB);
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
  int NUM_GFS;
  const int* gfs_to_sync = nullptr;
"""
    switch_case_code = generate_switch_statement_for_gf_types(Butcher_dict, MoL_method)
    file_output_str += switch_case_code
    file_output_str += r"""
  switch (type_ghost) {
    case EAST_GHOST:
      for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
        int i0 = Nxx0 + (2 * NGHOSTS) - 1;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
            for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
              gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)] = vals[IDXFACES0(which_gf, which_inner, i1, i2)];
            }
          }
          i0--;
        }
      }
      break;
    case WEST_GHOST:
      for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
        int i0 = 0;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
            for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
              gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)] = vals[IDXFACES0(which_gf, which_inner, i1, i2)];
            }
          }
          i0++;
        }
      }
      break;
    case NORTH_GHOST:
      for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
        int i1 = Nxx1 + (2 * NGHOSTS) - 1;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
              gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)] = vals[IDXFACES1(which_gf, which_inner, i0, i2)];
            }
          }
          i1--;
        }
      }
      break;
    case SOUTH_GHOST:
      for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
        int i1 = 0;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
              gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)] = vals[IDXFACES1(which_gf, which_inner, i0, i2)];
            }
          }
          i1++;
        }
      }
      break;
    case TOP_GHOST:
      for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
        int i2 = Nxx2 + (2 * NGHOSTS) - 1;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
              gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)] = vals[IDXFACES2(which_gf, which_inner, i0, i1)];
            }
          }
          i2--;
        }
      }
      break;
    case BOTTOM_GHOST:
      for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
        int i2 = 0;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
              gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)] = vals[IDXFACES2(which_gf, which_inner, i0, i1)];
            }
          }
          i2++;
        }
      }
      break;
  }
}
"""

    if nrpyelliptic_project:
        file_output_str += r"""
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

    file_output_str += r"""
void Timestepping::contribute_localsums_for_diagnostic_volume_integ() {
  // Build the canonical recipe list using the shared builder.
  diags_integration_recipe_t recipes[DIAGS_INTEGRATION_MAX_RECIPES];
  const int NUM_RECIPES = diags_integration_build_default_recipes(recipes);

  // Compute local contributions for each recipe and flatten into double vector
  std::vector<double> outdoubles;
  for (int r = 0; r < NUM_RECIPES; r++) {
    const int prefix = recipes[r].num_rules;
    REAL selected_volume = 0.0;
    REAL integrals[DIAGS_INTEGRATION_MAX_INTEGRANDS] = {0.0};

    diagnostics_integration_apply_rules(&commondata, griddata_chare, (const REAL *restrict *)diagnostic_gfs, recipes[r].rules, prefix,
                      recipes[r].integrands, recipes[r].num_integrands, &selected_volume, integrals);

    // Push volume then integrand integrals in order
    outdoubles.push_back((double)selected_volume);
    for (int s = 0; s < recipes[r].num_integrands; s++)
      outdoubles.push_back((double)integrals[s]);
  }

  CkCallback cb(CkIndex_Timestepping::report_sums_for_volume(NULL), thisProxy[CkArrayIndex3D(0,0,0)]);
  contribute(outdoubles, CkReduction::sum_double, cb);
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
  int NUM_GFS;
  const int* gfs_to_sync = nullptr;
  """
    file_output_str += switch_case_code
    file_output_str += r"""
  const int idx3_this_chare = IDX3_OF_CHARE(thisIndex.x, thisIndex.y, thisIndex.z);

  for (int which_dst_chare = 0; which_dst_chare < tot_num_dst_chares; which_dst_chare++) {
    REAL *restrict tmpBuffer_innerbc_send = griddata_chare[grid].tmpBuffers.tmpBuffer_innerbc_send[which_dst_chare];
    for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
      for (int which_srcpt = 0; which_srcpt < num_srcpts_tosend_each_chare[which_dst_chare]; which_srcpt++) {
        const int globalidx3srcpt = globalidx3_srcpts_tosend[which_dst_chare][which_srcpt];
        const int localidx3srcpt = globalidx3pt_to_localidx3pt[globalidx3srcpt];
        const int idx2 = IDX2NONLOCALINNERBC(which_gf, which_srcpt, num_srcpts_tosend_each_chare[which_dst_chare]);
        tmpBuffer_innerbc_send[idx2] = gfs[IDX4pt(gfs_to_sync[which_gf], localidx3srcpt)];
      }
    }
    int dst_chare_index0;
    int dst_chare_index1;
    int dst_chare_index2;
    REVERSE_IDX3GENERAL(idx3_of_dst_chares[which_dst_chare], Nchare0, Nchare1, dst_chare_index0, dst_chare_index1, dst_chare_index2);"""

    switch_case_code_for_entry_method = (
        generate_switch_statement_for_gf_types_for_entry_method(
            Butcher_dict, MoL_method
        )
    )
    file_output_str += switch_case_code_for_entry_method
    file_output_str += """
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
  int NUM_GFS;
  const int* gfs_to_sync = nullptr;
  const int8_t* gf_parity_types = nullptr;
"""
    file_output_str += generate_switch_statement_for_gf_types(
        Butcher_dict, MoL_method, set_parity_types=True
    )
    file_output_str += r"""
  apply_bcs_inner_only_nonlocal(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].bcstruct, &griddata_chare[grid].nonlocalinnerbcstruct, NUM_GFS, gfs, gfs_to_sync, gf_parity_types, griddata_chare[grid].tmpBuffers.tmpBuffer_innerbc_receiv);
}
"""

    if enable_charm_checkpointing:
        file_output_str += generate_PUP_code()

    if enable_BHaHAHA:
        file_output_str += r"""
void Timestepping::send_bhahaha_gfs_to_corresponding_interpolator_chare(const int grid) {
  send_interp_gfs_to_corresponding_interpolator_chare(grid, griddata_chare[grid].gridfuncs.y_n_gfs, bhahaha_gf_interp_indices, BHAHAHA_NUM_INTERP_GFS,
                                                      INTERP_REQUEST_BHAHAHA);
}

void Timestepping::send_interp_gfs_to_corresponding_interpolator_chare(const int grid, const REAL *restrict gfs, const int *restrict gf_indices,
                                                                       const int num_gfs, const int request_type) {
  const int Nchare0 = commondata.Nchare0;
  const int Nchare1 = commondata.Nchare1;
  const int Nchare2 = commondata.Nchare2;
  const int Nxx_plus_2NGHOSTS0 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS2;
  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

  REAL *restrict tmpBuffer = griddata_chare[grid].tmpBuffers.tmpBuffer_bhahaha_gfs;


  int idx = 0;
  for (int which_gf = 0; which_gf < num_gfs; which_gf++) {
    for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
      for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
        for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
          tmpBuffer[idx] = gfs[IDX4(gf_indices[which_gf], i0, i1, i2)];
          idx++;
        }
      }
    }
  }
  interpolator3dArray[CkArrayIndex3D(thisIndex.x, thisIndex.y, thisIndex.z)].receiv_interp_gfs(request_type, num_gfs, commondata.nn,
                                                                                               num_gfs * Nxx_plus_2NGHOSTS_tot, tmpBuffer);
}

"""
        if enable_psi4:
            file_output_str += r"""
void Timestepping::send_psi4_gfs_to_corresponding_interpolator_chare(const int grid) {
  const int psi4_gf_indices[2] = {DIAG_PSI4_REGF, DIAG_PSI4_IMGF};
  send_interp_gfs_to_corresponding_interpolator_chare(grid, diagnostic_gfs[grid], psi4_gf_indices, 2, INTERP_REQUEST_PSI4);
}
"""

    file_output_str += r"""
#include "timestepping.def.h"
"""

    timestepping_cpp_file = project_Path / "timestepping.cpp"
    with timestepping_cpp_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_timestepping_ci(
    project_dir: str,
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str,
    post_non_y_n_auxevol_mallocs: str,
    pre_MoL_step_forward_in_time: str = "",
    post_MoL_step_forward_in_time: str = "",
    outer_bcs_type: str = "radiation",
    enable_psi4: bool = False,
    nrpyelliptic_project: bool = False,
    enable_charm_checkpointing: bool = False,
    enable_BHaHAHA: bool = False,
) -> None:
    """
    Generate timestepping.ci.

    :param project_dir: Directory where the project C code is output.
    :param Butcher_dict: Dictionary containing Butcher tableau for the MoL method.
    :param MoL_method: Method of Lines (MoL) method name.
    :param post_non_y_n_auxevol_mallocs: Function calls after memory is allocated for non y_n and auxevol gridfunctions, default is an empty string.
    :param pre_MoL_step_forward_in_time: Code for handling pre-right-hand-side operations, default is an empty string.
    :param post_MoL_step_forward_in_time: Code for handling post-right-hand-side operations, default is an empty string.
    :param outer_bcs_type: type of outer boundary BCs to apply. Only options are radiation or extrapolation in superB.
    :param enable_psi4: Whether or not to enable psi4 diagnostics.
    :param nrpyelliptic_project: If True, enable NRPyElliptic project mode (enables residual
                    diagnostics and NRPyElliptic-specific behavior).
    :param enable_charm_checkpointing: Enable checkpointing using Charm++.
    :param enable_BHaHAHA: If True, add creation of horizon_finder and interpolator3d chares and communication with them.
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = r"""module timestepping {
  include "BHaH_defines.h";
  include "BHaH_function_prototypes.h";
  include "commondata_object.h";
  include "ckio.h";
  include "pup_stl.h";
  """

    file_output_str += r"""
  array [3D] Timestepping {
    entry Timestepping(CommondataObject &inData);
    entry void ready_1d_y(Ck::IO::FileReadyMsg *m);
    entry void ready_1d_z(Ck::IO::FileReadyMsg *m);
    entry void ready_2d_xy(Ck::IO::FileReadyMsg *m);
    entry void ready_2d_yz(Ck::IO::FileReadyMsg *m);
    // Step 5: MAIN SIMULATION LOOP
    entry void start() {
"""

    if enable_BHaHAHA:
        file_output_str += r"""
      // Only Timestepping chare (0,0,0) launches the horizon finder chares.
      if (thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
        serial {
          GriddataObject griddataObj;
          griddataObj.size_griddata = commondata.NUMGRIDS;
          griddataObj.griddata = new griddata_struct[commondata.NUMGRIDS];
          for(int grid=0; grid<commondata.NUMGRIDS; grid++){
            griddataObj.griddata[grid] = griddata[grid];

          }
          CommondataObject commondataObj;
          commondataObj.commondata = commondata;
          horizon_finderProxy.start(commondataObj, griddataObj);
        }
      }

      // Each Timestepping chare starts its corresponding Interpolator3d chare.
      serial {
          GriddataObject griddataObj;
          griddataObj.size_griddata = commondata.NUMGRIDS;
          griddataObj.griddata = new griddata_struct[commondata.NUMGRIDS];
          for(int grid=0; grid<commondata.NUMGRIDS; grid++){
            griddataObj.griddata[grid] = griddata_chare[grid];
          }
          CommondataObject commondataObj;
          commondataObj.commondata = commondata;
          interpolator3dArray[CkArrayIndex3D(thisIndex.x, thisIndex.y, thisIndex.z)].start(commondataObj, griddataObj);
      }
"""

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

    # If anything other than NRPy elliptic
    if not nrpyelliptic_project:
        file_output_str += """
      serial {
        initial_data(&commondata, griddata_chare, INITIALDATA_BIN_ONE);
        initial_data(&commondata, griddata_chare, INITIALDATA_APPLYBCS_INNERONLY);
      }"""
        file_output_str += generate_send_nonlocalinnerbc_data_code(
            "Y_N_GFS_INITIALDATA_PART1"
        )
        file_output_str += generate_process_nonlocalinnerbc_code(
            "Y_N_GFS_INITIALDATA_PART1"
        )
        file_output_str += (
            """if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {"""
        )
        file_output_str += generate_send_nonlocalinnerbc_data_code("AUXEVOL_GFS")
        file_output_str += generate_process_nonlocalinnerbc_code("AUXEVOL_GFS")
        file_output_str += """}"""
        file_output_str += """
      serial {
        initial_data(&commondata, griddata_chare, INITIALDATA_BIN_TWO);
        initial_data(&commondata, griddata_chare, INITIALDATA_APPLYBCS_OUTEREXTRAPANDINNER);
      }"""
        if post_non_y_n_auxevol_mallocs:
            file_output_str += (
                """   // Step 4.a: Functions called after memory for non-y_n and auxevol gridfunctions is allocated.
      serial {
"""
                + post_non_y_n_auxevol_mallocs
                + """
      }
"""
            )

        file_output_str += generate_send_nonlocalinnerbc_data_code(
            "Y_N_GFS_INITIALDATA_PART2"
        )
        file_output_str += generate_process_nonlocalinnerbc_code(
            "Y_N_GFS_INITIALDATA_PART2"
        )
        file_output_str += (
            """if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {"""
        )
        file_output_str += generate_send_nonlocalinnerbc_data_code("AUXEVOL_GFS")
        file_output_str += generate_process_nonlocalinnerbc_code("AUXEVOL_GFS")
        file_output_str += """}"""
    # If NRPy elliptic
    else:
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
        file_output_str += (
            """if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {"""
        )
        file_output_str += generate_send_neighbor_data_code(
            "AUXEVOL_GFS", grid_split_direction
        )
        file_output_str += generate_process_ghost_code(
            loop_direction, pos_ghost_type, neg_ghost_type, nchare_var
        )
        file_output_str += """}"""

    file_output_str += r"""
      while (commondata.time < commondata.t_final) { // Main loop to progress forward in time.
        serial {
          time_start = commondata.time;
        }
      """

    file_output_str += r"""
        serial {
          const REAL currtime = commondata.time, currdt = commondata.dt, outevery = commondata.diagnostics_output_every;
          write_diagnostics_this_step = fabs(round(currtime / outevery) * outevery - currtime) < 0.5 * currdt;
        """
    if enable_charm_checkpointing:
        file_output_str += r"""
          write_chckpt_this_step = fabs(round(commondata.time / commondata.checkpoint_every) * commondata.checkpoint_every -
                                         commondata.time) < 0.5 * commondata.dt;"""
    file_output_str += r"""
        }"""

    file_output_str += r"""
            //START DIAGNOSTICS
            if (write_diagnostics_this_step) {
              serial {
                // Allocate temporary storage for diagnostic_gfs.
                for (int grid = 0; grid < commondata.NUMGRIDS; grid++) {
                  const int Nxx_plus_2NGHOSTS0 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS0;                                                                      \
                  const int Nxx_plus_2NGHOSTS1 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS1;                                                                      \
                  const int Nxx_plus_2NGHOSTS2 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS2;
                  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
                  BHAH_MALLOC(diagnostic_gfs[grid], TOTAL_NUM_DIAG_GFS * Nxx_plus_2NGHOSTS_tot * sizeof(REAL));

                  // Initialize to NaNs so any unfilled outer-boundary points are obvious in diagnostics.
                  for (int ii = 0; ii < TOTAL_NUM_DIAG_GFS * Nxx_plus_2NGHOSTS_tot; ii++) {
                    diagnostic_gfs[grid][ii] = NAN;
                  }
                } // END LOOP over grids

                // Set diagnostics_gfs -- see nrpy/infrastructures/BHaH/[project]/diagnostics/ for definition.
                diagnostic_gfs_set(&commondata, griddata_chare, diagnostic_gfs);

                // Alias MoL's diagnostic_output_gfs to Timestepping's per-grid diagnostic_gfs storage.
                for (int grid = 0; grid < commondata.NUMGRIDS; ++grid) {
                  griddata_chare[grid].gridfuncs.diagnostic_output_gfs = diagnostic_gfs[grid];
                }
              }
                """
    # Sync psi4
    if enable_psi4:
        file_output_str += generate_send_nonlocalinnerbc_data_code(
            "DIAGNOSTIC_OUTPUT_GFS"
        )
        file_output_str += generate_process_nonlocalinnerbc_code(
            "DIAGNOSTIC_OUTPUT_GFS"
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

            file_output_str += generate_send_neighbor_data_code(
                "DIAGNOSTIC_OUTPUT_GFS", grid_split_direction
            )
            file_output_str += generate_process_ghost_code(
                loop_direction, pos_ghost_type, neg_ghost_type, nchare_var
            )

    file_output_str += r"""
              serial {
"""
    if enable_psi4:
        file_output_str += r"""
                // Send psi4 diagnostic gfs to interpolator chares for psi4 decomposition.
                if (commondata.num_psi4_extraction_radii > 0) {
                  send_psi4_gfs_to_corresponding_interpolator_chare(grid);
                }
"""
    file_output_str += r"""
                // Diagnostics center
                diagnostics_ckio(Ck::IO::Session(), DIAGNOSTICS_WRITE_CENTER);

                // Execute volume-integration recipe for this chare and contribute results
                diagnostics_ckio(Ck::IO::Session(), DIAGNOSTICS_VOLUME);
              }"""

    file_output_str += r"""
            }
         """

    if enable_BHaHAHA:
        file_output_str += r"""
          serial {
            do_horizon_find = true;  // default: yes, find horizon
            // STEP 1: Check if horizon find is scheduled for the current iteration.
            if (commondata.diagnostics_output_every <= 0 ||
                (commondata.nn % (int)(commondata.diagnostics_output_every / commondata.dt + 0.5)) != 0) {
              int bah_find_every = 1;  // Placeholder: find every iteration. This should be a commondata param.
              if (commondata.diagnostics_output_every > commondata.dt) {
                // A basic way to get find_every from time interval
                bah_find_every = (int)(commondata.diagnostics_output_every / commondata.dt + 0.5);
                if (bah_find_every == 0)
                  bah_find_every = 1;
              }
              if (bah_find_every <= 0 || (commondata.nn % bah_find_every) != 0) {
                do_horizon_find = false;  // not scheduled this iteration
              }
            } // END IF: diagnostics_output_every > 0
          }
          // If do horison for this step, send the entire set of required BHaHAHA gridfunctions to the corresponding Interpolator3d chare
          if (do_horizon_find) {
              serial { send_bhahaha_gfs_to_corresponding_interpolator_chare(grid); }
          }
        """

    file_output_str += r"""
        if (thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
          serial {
            progress_indicator(&commondata, griddata_chare);
            if (commondata.time + commondata.dt > commondata.t_final)
              printf("\n");
          }
          if (write_diagnostics_this_step) {
            serial {
              count_filewritten = 0;
              {
                char filename[512];
                build_outfile_name(filename, sizeof filename, "out1d-y",
                                   griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_1d_y,
                                   &commondata, /*include_time=*/1);
                Ck::IO::Options opts;
                CkCallback opened_1d_y(CkIndex_Timestepping::ready_1d_y(NULL), thisProxy);
                Ck::IO::open(filename, opened_1d_y, opts);
              }
              {
                char filename[512];
                build_outfile_name(filename, sizeof filename, "out1d-z",
                                   griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_1d_z,
                                   &commondata, /*include_time=*/1);
                Ck::IO::Options opts;
                CkCallback opened_1d_z(CkIndex_Timestepping::ready_1d_z(NULL), thisProxy);
                Ck::IO::open(filename, opened_1d_z, opts);
              }
              {
                char filename[512];
                build_outfile_name(filename, sizeof filename, "out2d-xy",
                                   griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_2d_xy,
                                   &commondata, /*include_time=*/1);
                Ck::IO::Options opts;
                CkCallback opened_2d_xy(CkIndex_Timestepping::ready_2d_xy(NULL), thisProxy);
                Ck::IO::open(filename, opened_2d_xy, opts);
              }
              {
                char filename[512];
                build_outfile_name(filename, sizeof filename, "out2d-yz",
                                   griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_2d_yz,
                                   &commondata, /*include_time=*/1);
                Ck::IO::Options opts;
                CkCallback opened_2d_yz(CkIndex_Timestepping::ready_2d_yz(NULL), thisProxy);
                Ck::IO::open(filename, opened_2d_yz, opts);
              }
            }
"""

    file_output_str += generate_diagnostics_code(
        "1d", "y", "totsizeinbytes_1d_y", "DIAGNOSTICS_WRITE_Y"
    )
    file_output_str += generate_diagnostics_code(
        "1d", "z", "totsizeinbytes_1d_z", "DIAGNOSTICS_WRITE_Z"
    )
    file_output_str += generate_diagnostics_code(
        "2d", "xy", "totsizeinbytes_2d_xy", "DIAGNOSTICS_WRITE_XY"
    )
    file_output_str += generate_diagnostics_code(
        "2d", "yz", "totsizeinbytes_2d_yz", "DIAGNOSTICS_WRITE_YZ"
    )

    file_output_str += r"""
            if (count_filewritten == expected_count_filewritten) {
              serial {thisProxy.continue_timestepping(); }
            }
          } else {
            serial {thisProxy.continue_timestepping(); }
          }
        }
"""
    file_output_str += r"""
        when continue_timestepping() {
          //END DIAGNOSTICS
          if (write_diagnostics_this_step) {
            serial {
              // Free temporary storage allocated to diagnostic_gfs.
              for (int grid = 0; grid < commondata.NUMGRIDS; grid++)
                free(diagnostic_gfs[grid]);
            }
          }
        }
    """

    if enable_charm_checkpointing:
        file_output_str += r"""
        // periodically checkpointing
        if (write_chckpt_this_step) {
          serial {
            // coordinate to start checkpointing
            if (thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
              CkPrintf("[%d] CHECKPOINT at step %d\n", CkMyPe(), commondata.nn);
            }
            contribute(CkCallback(CkReductionTarget(Timestepping, startCheckpoint), thisProxy(0, 0, 0)));
          }
          if (thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
            when startCheckpoint() serial {
              CkCallback cb(CkIndex_Timestepping::recvCheckPointDone(), thisProxy);
              CkStartCheckpoint("log", cb);
            }
          }
          when recvCheckPointDone() { }
        }
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
        rhs_output_exprs_list = superB.MoL.generate_rhs_output_exprs(
            Butcher_dict, MoL_method, s + 1
        )
        post_rhs_output_list = superB.MoL.generate_post_rhs_output_list(
            Butcher_dict, MoL_method, s + 1
        )
        file_output_str += generate_mol_step_forward_code(
            f"RK_SUBSTEP_K{s+1}",
            rhs_output_exprs_list,
            post_rhs_output_list,
            outer_bcs_type=outer_bcs_type,
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
            file_output_str += (
                """if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {"""
            )
            file_output_str += generate_send_neighbor_data_code(
                "AUXEVOL_GFS", grid_split_direction
            )
            file_output_str += generate_process_ghost_code(
                loop_direction, pos_ghost_type, neg_ghost_type, nchare_var
            )
            file_output_str += """}"""

    file_output_str += r"""
        """

    if post_MoL_step_forward_in_time != "":
        file_output_str += post_MoL_step_forward_in_time
    else:
        file_output_str += "  // (nothing here; specify by setting post_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"

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
    if not enable_BHaHAHA:
        file_output_str += r"""
      serial{
        mainProxy.done();
      }"""
    file_output_str += r"""
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
    entry void diagnostics_ckio(Ck::IO::Session token, const int which_diagnostics_part) {
      serial {
        const int thisIndex_arr[3] = {thisIndex.x, thisIndex.y, thisIndex.z};
        if (which_diagnostics_part == DIAGNOSTICS_VOLUME) {
          // Compute local volume-integral contributions and participate in a global reduction
          contribute_localsums_for_diagnostic_volume_integ();
        } else {
          diagnostics(&commondata, griddata, griddata_chare, (const REAL *restrict *)diagnostic_gfs, thisIndex_arr, which_grid_diagnostics, token, which_diagnostics_part);
        }
      }
    }
    entry void east_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void west_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void north_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void south_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void top_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void bottom_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void continue_timestepping();
    entry void receiv_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, int globalidx3_srcpts[num_srcpts]);"""
    file_output_str += generate_entry_methods_for_receiv_nonlocalinnerbc_for_gf_types(
        Butcher_dict,
        MoL_method,
        outer_bcs_type,
        enable_psi4,
        nrpyelliptic_project,
    )

    if nrpyelliptic_project:
        file_output_str += r"""
    entry void receiv_wavespeed_at_outer_boundary(REAL wavespeed_at_outer_boundary);"""

    if enable_charm_checkpointing:
        file_output_str += r"""
    entry [reductiontarget] void startCheckpoint(); //reduction to start checkpointing
    entry void recvCheckPointDone();  //checkpointing done, resume application
    """

    file_output_str += r"""
    entry void report_sums_for_volume(CkReductionMsg * msg) {
      serial {
        // Build canonical recipes via the shared builder so parsing matches contribution order.
        diags_integration_recipe_t recipes[DIAGS_INTEGRATION_MAX_RECIPES];
        const int NUM_RECIPES = diags_integration_build_default_recipes(recipes);

        // Expected reduced array length: for each recipe, one volume + one entry per integrand
        int expectedLen = 0;
        for (int r = 0; r < NUM_RECIPES; r++)
          expectedLen += 1 + recipes[r].num_integrands;

        int reducedArrSize = msg->getSize() / sizeof(double);
        CkAssert(reducedArrSize == expectedLen);
        double *output = (double *)msg->getData();

        diags_integration_results_t integration_results;
        integration_results.num_recipe_results = 0;

        int idx = 0;
        for (int r = 0; r < NUM_RECIPES; r++) {
          diags_integration_recipe_result_t *result = &integration_results.recipe_results[integration_results.num_recipe_results];
          result->recipe_name = recipes[r].name;
          result->proper_volume = (REAL)output[idx++];
          result->num_integrand_results = recipes[r].num_integrands;

          // Fill integrand-level data
          for (int s = 0; s < recipes[r].num_integrands; s++) {
            diags_integration_integrand_result_t *int_result = &result->integrand_results[s];
            const diags_integration_integrand_spec_t *spec_ptr = &recipes[r].integrands[s];
            int_result->gf_index = spec_ptr->gf_index;
            int_result->is_squared = spec_ptr->is_squared;
            if (spec_ptr->gf_index == diags_integration_UNIT_INTEGRAND_GFINDEX)
              int_result->gf_name = "one";
            else
              int_result->gf_name = diagnostic_gf_names[spec_ptr->gf_index];

            int_result->integral = (REAL)output[idx++];
            if (spec_ptr->is_squared) {
              const double integral_value = (double)int_result->integral;
              int_result->L2 = sqrt(fmax(integral_value, 0.0));
              int_result->RMS = (result->proper_volume > 0.0) ? sqrt(fmax(integral_value, 0.0) / (double)result->proper_volume) : 0.0;
              int_result->L2_valid = 1;
              int_result->RMS_valid = 1;
            } else {
              int_result->L2 = 0.0;
              int_result->RMS = 0.0;
              int_result->L2_valid = 0;
              int_result->RMS_valid = 0;
            }
          }

          // Write per-recipe file on root chare only
          if (thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
            char filename[256];
            snprintf(filename, sizeof(filename), "out3d-integrals-conv_factor%.2f-%s.txt", commondata.convergence_factor, recipes[r].name);
            FILE *file_ptr = diags_integration_open_output_file(commondata.nn, filename);
            if (commondata.nn == 0)
              diags_integration_write_header(file_ptr, r, &recipes[r]);

            REAL integrals_arr[DIAGS_INTEGRATION_MAX_INTEGRANDS];
            for (int s = 0; s < recipes[r].num_integrands; s++)
              integrals_arr[s] = result->integrand_results[s].integral;

            diags_integration_write_row(file_ptr, (double)commondata.time, (double)result->proper_volume, recipes[r].integrands, integrals_arr,
                                        recipes[r].num_integrands);
            fclose(file_ptr);
          }

          integration_results.num_recipe_results++;
        }    """
    if nrpyelliptic_project:
        file_output_str += r"""
        // Update a convenience field (same as original user-edit capture)
        REAL rms_residual_r_80 = NAN;
        diags_integration_get_rms(&integration_results, "sphere_R_80", DIAG_RESIDUALGF, &rms_residual_r_80);
        commondata.log10_current_residual = log10(rms_residual_r_80);    """

    file_output_str += r"""

        delete msg;
      }
    }
"""
    file_output_str += r"""
  };
};
"""

    timestepping_ci_file = project_Path / "timestepping.ci"
    with timestepping_ci_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_timestepping_h_cpp_ci_register_CFunctions(
    project_dir: str,
    MoL_method: str = "RK4",
    enable_rfm_precompute: bool = False,
    post_non_y_n_auxevol_mallocs: str = "",
    pre_MoL_step_forward_in_time: str = "",
    post_MoL_step_forward_in_time: str = "",
    outer_bcs_type: str = "radiation",
    enable_psi4: bool = False,
    nrpyelliptic_project: bool = False,
    enable_charm_checkpointing: bool = False,
    enable_BHaHAHA: bool = False,
) -> None:
    """
    Output timestepping h, cpp, and ci files and register C functions.

    :param project_dir: Directory where the project C code is output
    :param MoL_method: Method of Lines (MoL) method name, default is "RK4".
    :param enable_rfm_precompute: Enable RFM precompute, default is False.
    :param post_non_y_n_auxevol_mallocs: Function calls after memory is allocated for non y_n and auxevol gridfunctions, default is an empty string.
    :param pre_MoL_step_forward_in_time: Code for handling pre-right-hand-side operations, default is an empty string.
    :param post_MoL_step_forward_in_time: Code for handling post-right-hand-side operations, default is an empty string.
    :param outer_bcs_type: type of outer boundary BCs to apply. Only options are radiation or extrapolation in superB.
    :param enable_psi4: Whether or not to enable psi4 diagnostics.
    :param nrpyelliptic_project: If True, enable NRPyElliptic project mode (enables residual
                    diagnostics and NRPyElliptic-specific behavior).
    :param enable_charm_checkpointing: Enable checkpointing using Charm++.
    :param enable_BHaHAHA: If True, add creation of horizon_finder and interpolator3d chares and communication with them.
    """
    # For NRPy elliptic: register parameter wavespeed at outer boundary
    if nrpyelliptic_project:
        _wavespeed_at_outer_boundary = par.register_CodeParameter(
            "REAL", __name__, "wavespeed_at_outer_boundary", 0.0, commondata=False
        )

    output_timestepping_h(
        project_dir=project_dir,
        nrpyelliptic_project=nrpyelliptic_project,
        enable_psi4=enable_psi4,
        enable_charm_checkpointing=enable_charm_checkpointing,
        enable_BHaHAHA=enable_BHaHAHA,
    )

    Butcher_dict = (
        BHaH.MoLtimestepping.rk_butcher_table_dictionary.generate_Butcher_tables()
    )

    output_timestepping_cpp(
        project_dir=project_dir,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_CurviBCs=True,
        Butcher_dict=Butcher_dict,
        MoL_method=MoL_method,
        nrpyelliptic_project=nrpyelliptic_project,
        enable_psi4=enable_psi4,
        enable_charm_checkpointing=enable_charm_checkpointing,
        enable_BHaHAHA=enable_BHaHAHA,
    )

    output_timestepping_ci(
        project_dir=project_dir,
        MoL_method=MoL_method,
        post_non_y_n_auxevol_mallocs=post_non_y_n_auxevol_mallocs,
        pre_MoL_step_forward_in_time=pre_MoL_step_forward_in_time,
        post_MoL_step_forward_in_time=post_MoL_step_forward_in_time,
        outer_bcs_type=outer_bcs_type,
        enable_psi4=enable_psi4,
        nrpyelliptic_project=nrpyelliptic_project,
        Butcher_dict=Butcher_dict,
        enable_charm_checkpointing=enable_charm_checkpointing,
        enable_BHaHAHA=enable_BHaHAHA,
    )

    register_CFunction_timestepping_malloc(
        enable_BHaHAHA=enable_BHaHAHA,
    )

    register_CFunction_timestepping_free_memory(
        enable_BHaHAHA=enable_BHaHAHA,
    )

    # Register temporary buffers for face data communication to griddata_struct:
    BHaH.griddata_commondata.register_griddata_commondata(
        __name__,
        "tmpBuffers_struct tmpBuffers",
        "temporary buffer for sending face data to neighbor chares",
    )
