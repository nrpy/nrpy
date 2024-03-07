"""

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from typing import List, Tuple, Optional, Dict
import itertools
import nrpy.c_function as cfc
import sys
from pathlib import Path
import nrpy.params as par
import nrpy.grid as gri
from nrpy.infrastructures.BHaH import griddata_commondata
from nrpy.helpers.generic import clang_format
from nrpy.infrastructures.BHaH import BHaH_defines_h

def generate_mol_step_forward_code(rk_substep):
    return f"""
    serial{{
        MoL_step_forward_in_time(&commondata, griddata_chare, time_start, {rk_substep});
    }}
    """

def generate_send_neighbor_data_code(which_gf, neighbor_direction):
    return f"""
    serial {{
        send_neighbor_data({which_gf}, {neighbor_direction}, grid);
    }}
    """

def generate_ghost_code(axis, axis_index, pos_ghost_type, neg_ghost_type, nchare_var):
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

def generate_diagnostics_code(dimension, direction, num_fields, tot_num_diagnostic_pts):
    code = f"""
    // {dimension} {direction}
    when ready_{dimension}_{direction}(Ck::IO::FileReadyMsg *m_{dimension}_{direction}) {{
        serial {{
        f_{dimension}_{direction} = m_{dimension}_{direction}->file;
        CkCallback sessionStart_{dimension}_{direction}(CkIndex_Timestepping::start_write_{dimension}_{direction}(0), thisProxy);
        CkCallback sessionEnd_{dimension}_{direction}(CkIndex_Timestepping::test_written_{dimension}_{direction}(0), thisProxy);
        int totsizeinbytes = 23 * {num_fields} * {tot_num_diagnostic_pts};
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
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Allocate memory for temporary buffers used to communicate face data"
    c_type = "void"
    name = "timestepping_malloc_tmpBuffer"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, tmpBuffers_struct *restrict tmpBuffers"
    body = """
const int Nxx_plus_2NGHOSTS_face0 = Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
const int Nxx_plus_2NGHOSTS_face1 = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS2;
const int Nxx_plus_2NGHOSTS_face2 = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1;

tmpBuffers->tmpBuffer_EW = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * NGHOSTS * Nxx_plus_2NGHOSTS_face0);
tmpBuffers->tmpBuffer_NS = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * NGHOSTS * Nxx_plus_2NGHOSTS_face1);
tmpBuffers->tmpBuffer_TB = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * NGHOSTS * Nxx_plus_2NGHOSTS_face2);
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )

def register_CFunction_timestepping_free_memory() -> None:
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Free memory for temporary buffers used to communicate face data"
    c_type = "void"
    name = "timestepping_free_memory_tmpBuffer"
    params = "tmpBuffers_struct *restrict tmpBuffers"
    body = """
if (tmpBuffers->tmpBuffer_EW != NULL) free(tmpBuffers->tmpBuffer_EW);
if (tmpBuffers->tmpBuffer_NS != NULL) free(tmpBuffers->tmpBuffer_NS);
if (tmpBuffers->tmpBuffer_TB != NULL) free(tmpBuffers->tmpBuffer_TB);
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
    )



def output_timestepping_h(
    project_dir: str,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    """
    Generate timestepping.h
    :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = r"""#ifndef __TIMESTEPPING_H__
#define __TIMESTEPPING_H__
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "timestepping.decl.h"

class Timestepping : public CBase_Timestepping {
  Timestepping_SDAG_CODE

  private:
    /// Member Variables (Object State) ///
    commondata_struct commondata;
    griddata_struct *griddata;
    griddata_struct *griddata_chare;
    bool is_boundarychare;
    REAL time_start;
    //bool contains_gridcenter;
    bool write_diagnostics_this_step;
    const int which_grid_diagnostics = 0;
    Ck::IO::File f_1d_y;
    Ck::IO::File f_1d_z;
    Ck::IO::File f_2d_xy;
    Ck::IO::File f_2d_yz;
    int grid; // grid in "for (grid = 0; grid < commondata.NUMGRIDS; grid++)" in .ci file must be declared here


    /// Member Functions (private) ///
    void send_neighbor_data(const int type_gfs, const int dir, const int grid);
    void process_ghost(const int type_ghost, const int type_gfs, const int len_tmpBuffer, const REAL *restrict tmpBuffer, const int grid);


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


def output_timestepping_cpp(
    project_dir: str,
    MoL_method: str,
    initial_data_desc: str = "",
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
    boundary_conditions_desc: str = "",
    prefunc: str = "",
    initialize_constant_auxevol: bool = False,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    """
    Generate timestepping.cpp

    :param MoL_method: Method of Lines algorithm used to step forward in time.
    :param initial_data_desc: Description for initial data, default is an empty string.
    :param enable_rfm_precompute: Enable rfm precomputation, default is False.
    :param enable_CurviBCs: Enable CurviBCs, default is False.
    :param boundary_conditions_desc: Description of the boundary conditions, default is an empty string.
    :param prefunc: String that appears before main(). DO NOT populate this, except when debugging, default is an empty string.
    :param initialize_constant_auxevol: If set to True, `initialize_constant_auxevol` function will be called during the simulation initialization phase to set these constants. Default is False.
    :param pre_MoL_step_forward_in_time: Code for handling pre-right-hand-side operations, default is an empty string.
    :param post_MoL_step_forward_in_time: Code for handling post-right-hand-side operations, default is an empty string.
    :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".
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
        file_output_str += "*Step 4.a: Set AUXEVOL gridfunctions that will never change in time."
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

  for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
    // Step 2: Initial data are set on y_n_gfs gridfunctions. Allocate storage for them first.
    MoL_malloc_y_n_gfs(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].gridfuncs);
  }

  // Step 3: Finalize initialization: set up initial data, etc.
  initial_data(&commondata, griddata_chare);

  // Step 4: Allocate storage for non-y_n gridfunctions, needed for the Runge-Kutta-like timestepping
  for(int grid=0; grid<commondata.NUMGRIDS; grid++)
    MoL_malloc_non_y_n_gfs(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].gridfuncs);

  // Allocate storage for temporary buffers, needed for communicating face data
  for(int grid=0; grid<commondata.NUMGRIDS; grid++)
    timestepping_malloc_tmpBuffer(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].tmpBuffers);

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
    timestepping_free_memory_tmpBuffer(&griddata_chare[grid].tmpBuffers);"""
    if enable_rfm_precompute:
        file_output_str += r"""
    rfm_precompute_free(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].rfmstruct);"""
    if enable_CurviBCs:
        file_output_str += r"""
    free(griddata[grid].bcstruct.inner_bc_array);
    free(griddata_chare[grid].bcstruct.inner_bc_array);
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
    free(griddata_chare[grid].charecommstruct.globalidx3pt_to_chareidx3);
    free(griddata_chare[grid].charecommstruct.globalidx3pt_to_localidx3pt);
    free(griddata_chare[grid].charecommstruct.localidx3pt_to_globalidx3pt);
  }
  free(griddata);
  free(griddata_chare);
}
"""
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
  switch (type_gfs) {
    case K_ODD:
      gfs = griddata_chare[grid].gridfuncs.k_odd_gfs;
      break;
    case K_EVEN:
      gfs = griddata_chare[grid].gridfuncs.k_even_gfs;
      break;
    case Y_N:
      gfs = griddata_chare[grid].gridfuncs.y_n_gfs;
      break;
    default:
      break;
  }
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
  switch (type_gfs) {
  case K_ODD:
    gfs = griddata_chare[grid].gridfuncs.k_odd_gfs;
    break;
  case K_EVEN:
    gfs = griddata_chare[grid].gridfuncs.k_even_gfs;
    break;
  case Y_N:
    gfs = griddata_chare[grid].gridfuncs.y_n_gfs;
    break;
  default:
    break;
  }
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
#include "timestepping.def.h"
"""


    timestepping_cpp_file = project_Path / "timestepping.cpp"
    with timestepping_cpp_file.open("w", encoding="utf-8") as file:
        file.write(
            clang_format(file_output_str, clang_format_options=clang_format_options)
        )


def output_timestepping_ci(
    project_dir: str,
    MoL_method: str,
    pre_MoL_step_forward_in_time: str = "",
    post_MoL_step_forward_in_time: str = "",
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    """
    Generate timestepping.ci

    :param MoL_method: Method of Lines algorithm used to step forward in time.
    :param pre_MoL_step_forward_in_time: Code for handling pre-right-hand-side operations, default is an empty string.
    :param post_MoL_step_forward_in_time: Code for handling post-right-hand-side operations, default is an empty string.
    :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".
    """

    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = r"""module timestepping {
  include "BHaH_defines.h";
  include "BHaH_function_prototypes.h";
  include "commondata_object.h";
  include "ckio.h";
  array [3D] Timestepping {
    entry Timestepping(CommondataObject &inData);
    entry void ready_1d_y(Ck::IO::FileReadyMsg *m);
    entry void ready_1d_z(Ck::IO::FileReadyMsg *m);
    entry void ready_2d_xy(Ck::IO::FileReadyMsg *m);
    entry void ready_2d_yz(Ck::IO::FileReadyMsg *m);
    // Step 5: MAIN SIMULATION LOOP
    entry void start() {
      while (commondata.time < commondata.t_final) { // Main loop to progress forward in time.
        serial {
          time_start = commondata.time;
        }
        serial {
          write_diagnostics_this_step = fabs(round(commondata.time / commondata.diagnostics_output_every) * commondata.diagnostics_output_every - commondata.time) < 0.5 * commondata.dt;
        }
        // Step 5.a: Main loop, part 1: Output diagnostics
        //serial {
        //  if (write_diagnostics_this_step && contains_gridcenter) {
        //    diagnostics(&commondata, griddata_chare, Ck::IO::Session(), OUTPUT_0D, which_grid_diagnostics);
        //  }
        //}
        // Create sessions for ckio file writing from first chare only
        if (write_diagnostics_this_step && thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
          serial {
            progress_indicator(&commondata, griddata_chare);
            if (commondata.time + commondata.dt > commondata.t_final)
              printf("\n");


            {
              char filename[256];
              sprintf(filename, griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_1d_y, commondata.convergence_factor, commondata.time);
              Ck::IO::Options opts;
              CkCallback opened_1d_y(CkIndex_Timestepping::ready_1d_y(NULL), thisProxy);
              Ck::IO::open(filename, opened_1d_y, opts);
            }
            {
              char filename[256];
              sprintf(filename, griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_1d_z, commondata.convergence_factor, commondata.time);
              Ck::IO::Options opts;
              CkCallback opened_1d_z(CkIndex_Timestepping::ready_1d_z(NULL), thisProxy);
              Ck::IO::open(filename, opened_1d_z, opts);
            }
            {
              char filename[256];
              sprintf(filename, griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_2d_xy, commondata.convergence_factor, commondata.time);
              Ck::IO::Options opts;
              CkCallback opened_2d_xy(CkIndex_Timestepping::ready_2d_xy(NULL), thisProxy);
              Ck::IO::open(filename, opened_2d_xy, opts);
            }
            {
              char filename[256];
              sprintf(filename, griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_2d_yz, commondata.convergence_factor, commondata.time);
              Ck::IO::Options opts;
              CkCallback opened_2d_yz(CkIndex_Timestepping::ready_2d_yz(NULL), thisProxy);
              Ck::IO::open(filename, opened_2d_yz, opts);
            }
          }
"""
    # Generate code for 1d y diagnostics
    file_output_str += generate_diagnostics_code("1d", "y", "griddata_chare[which_grid_diagnostics].diagnosticstruct.num_output_quantities + 1", "griddata_chare[which_grid_diagnostics].diagnosticstruct.tot_num_diagnostic_1d_y_pts")

    # Generate code for 1d z diagnostics
    file_output_str += generate_diagnostics_code("1d", "z", "griddata_chare[which_grid_diagnostics].diagnosticstruct.num_output_quantities + 1", "griddata_chare[which_grid_diagnostics].diagnosticstruct.tot_num_diagnostic_1d_z_pts")

    # Generate code for 2d xy diagnostics
    file_output_str += generate_diagnostics_code("2d", "xy", "griddata_chare[which_grid_diagnostics].diagnosticstruct.num_output_quantities + 2", "griddata_chare[which_grid_diagnostics].diagnosticstruct.tot_num_diagnostic_2d_xy_pts")

    # Generate code for 2d yz diagnostics
    file_output_str += generate_diagnostics_code("2d", "yz", "griddata_chare[which_grid_diagnostics].diagnosticstruct.num_output_quantities + 2", "griddata_chare[which_grid_diagnostics].diagnosticstruct.tot_num_diagnostic_2d_yz_pts")

    file_output_str += r"""
      }
"""
    if pre_MoL_step_forward_in_time != "":
        file_output_str += pre_MoL_step_forward_in_time
    else:
        file_output_str += "  // (nothing here; specify by setting pre_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"


    file_output_str += r"""
    """
    # Loop over RK substeps and axes
    for k in range(1, 5):
        rk_substep = f"RK_SUBSTEP_K{k}"
        file_output_str += generate_mol_step_forward_code(rk_substep)
        file_output_str += "for (grid = 0; grid < commondata.NUMGRIDS; grid++) {"
        for axis in ['x', 'y', 'z']:
            # do something with rk_substep and axis
            if axis == 'x':
                pos_ghost_type = x_pos_ghost_type
                neg_ghost_type = x_neg_ghost_type
                nchare_var = "commondata.Nchare0"
                direction = 'EAST_WEST'
            elif axis == 'y':
                pos_ghost_type = y_pos_ghost_type
                neg_ghost_type = y_neg_ghost_type
                nchare_var = "commondata.Nchare1"
                direction = 'NORTH_SOUTH'
            elif axis == 'z':
                pos_ghost_type = z_pos_ghost_type
                neg_ghost_type = z_neg_ghost_type
                nchare_var = "commondata.Nchare2"
                direction = 'TOP_BOTTOM'

            # Generate code for this RK substep and axis
            if rk_substep == 'RK_SUBSTEP_K1':
                which_gf = 'K_ODD'
            elif rk_substep == 'RK_SUBSTEP_K2':
                which_gf = 'K_EVEN'
            elif rk_substep == 'RK_SUBSTEP_K3':
                which_gf = 'K_ODD'
            elif rk_substep == 'RK_SUBSTEP_K4':
                which_gf = 'Y_N'
            else:
                raise ValueError(f"Unknown RK substep: {rk_substep}")

            file_output_str += generate_send_neighbor_data_code(which_gf, direction)
            file_output_str += generate_ghost_code(axis, 0, pos_ghost_type, neg_ghost_type, nchare_var)
        file_output_str += "}"

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
        diagnostics(&commondata, griddata_chare, token, which_output, which_grid_diagnostics);
      }
    }
    entry void east_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void west_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void north_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void south_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void top_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void bottom_ghost(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
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
    MoL_method: str,
    initial_data_desc: str = "",
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
    boundary_conditions_desc: str = "",
    prefunc: str = "",
    initialize_constant_auxevol: bool = False,
    pre_MoL_step_forward_in_time: str = "",
    post_MoL_step_forward_in_time: str = "",
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:

    output_timestepping_h(
      project_dir=project_dir,
    )

    output_timestepping_cpp(
      project_dir=project_dir,
      MoL_method=MoL_method,
      enable_rfm_precompute=enable_rfm_precompute,
      enable_CurviBCs=True,
      boundary_conditions_desc=boundary_conditions_desc,
    )

    output_timestepping_ci(
      project_dir=project_dir,
      MoL_method=MoL_method,
      pre_MoL_step_forward_in_time=pre_MoL_step_forward_in_time,
      post_MoL_step_forward_in_time=post_MoL_step_forward_in_time,
    )

    register_CFunction_timestepping_malloc()

    register_CFunction_timestepping_free_memory()

    # Register temporary buffers for face data communication to griddata_struct:
    griddata_commondata.register_griddata_commondata(
        __name__,
        "tmpBuffers_struct tmpBuffers",
        "temporary buffer for sending face data to neighbor chares",
    )
