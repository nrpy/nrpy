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

import itertools

def generate_mol_step_forward_code(rk_substep):
    return f"""
    serial{{
        MoL_step_forward_in_time(&commondata, griddata_chare, time_start, {rk_substep});
    }}
    """

def generate_send_neighbor_data_code(which_gf, neighbor_direction):
    return f"""
    serial {{
        send_neighbor_data({which_gf}, {neighbor_direction});
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
                process_ghost({pos_ghost_type}, type_gfs, len_tmpBuffer, tmpBuffer);
            }}
        }}
    }}
    if({this_index_var} > 0) {{
        when {neg_ghost_func}(int type_gfs, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]) {{
            serial {{
                process_ghost({neg_ghost_type}, type_gfs, len_tmpBuffer, tmpBuffer);
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
        thisProxy.diagnostics_ckio_{dimension}_{direction}(m_{dimension}_{direction}->session);
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
    bool contains_gridcenter;
    bool write_diagnostics_this_step;
    const int which_grid_diagnostics = 0;
    Ck::IO::File f_1d_y;
    Ck::IO::File f_1d_z;
    Ck::IO::File f_2d_xy;
    Ck::IO::File f_2d_yz;

    /// Member Functions (private) ///

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
#include "BHaH_defines.h"
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

  if (thisIndex.x == 0 ||
      thisIndex.y == 0 ||
      thisIndex.z == 0 ||
      thisIndex.x == (Nchare0-1) ||
      thisIndex.y == (Nchare1-1) ||
      thisIndex.z == (Nchare2-1)) {
      is_boundarychare = true;
    }

  commondata = inData.commondata;

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
    numerical_grids_chare(&commondata, griddata, griddata_chare, {thisIndex.x, thisIndex.y, thisIndex.z});
  }

  for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
    // Step 2: Initial data are set on y_n_gfs gridfunctions. Allocate storage for them first.
    MoL_malloc_y_n_gfs(&commondata, griddata_chare[grid].params, griddata_chare[grid].gridfuncs);
  }

  // Step 3: Finalize initialization: set up initial data, etc.
  initial_data(&commondata, griddata_chare);

  // Step 4: Allocate storage for non-y_n gridfunctions, needed for the Runge-Kutta-like timestepping
  for(int grid=0; grid<commondata.NUMGRIDS; grid++)
    MoL_malloc_non_y_n_gfs(&commondata, griddata_chare[grid].params, griddata_chare[grid].gridfuncs);
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
// destructor
Timestepping::~Timestepping() {
  // Step 5: Free all allocated memory
  for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
    MoL_free_memory_y_n_gfs(griddata_chare[grid].gridfuncs);
    MoL_free_memory_non_y_n_gfs(griddata_chare[grid].gridfuncs);"""
    if enable_rfm_precompute:
        file_output_str += r"""
    rfm_precompute_free(&commondata, griddata_chare[grid].params, griddata_chare[grid].rfmstruct);"""
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
  }
  free(griddata);
  free(griddata_chare);
}
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
  include "BHaH_defines.h"
  include "BHaH_function_prototypes.h"
  include "CommondataObject.h";
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
        serial {
          if (write_diagnostics_this_step && contains_gridcenter) {
            diagnostics_center(&commondata, griddata_chare);
          }
        }
        // Create sessions for ckio file writing from first chare only
        if (write_diagnostics_this_step && thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
          serial {
            {
              char filename[256];
              sprintf(filename, "out1d-y-conv_factor%.2f-t%08.2f.txt", commondata.convergence_factor, commondata.time);
              Ck::IO::Options opts;
              CkCallback opened_1d_y(CkIndex_Timestepping::ready_1d_y(NULL), thisProxy);
              Ck::IO::open(filename, opened_1d_y, opts);
            }
            {
              char filename[256];
              sprintf(filename, "out1d-z-conv_factor%.2f-t%08.2f.txt", commondata.convergence_factor, commondata.time);
              Ck::IO::Options opts;
              CkCallback opened_1d_z(CkIndex_Timestepping::ready_1d_z(NULL), thisProxy);
              Ck::IO::open(filename, opened_1d_z, opts);
            }
            {
              char filename[256];
              sprintf(filename, "out2d-xy-conv_factor%.2f-t%08.2f.txt", commondata.convergence_factor, commondata.time);
              Ck::IO::Options opts;
              CkCallback opened_2d_xy(CkIndex_Timestepping::ready_2d_xy(NULL), thisProxy);
              Ck::IO::open(filename, opened_2d_xy, opts);
            }
            {
              char filename[256];
              sprintf(filename, "out2d-yz-conv_factor%.2f-t%08.2f.txt", commondata.convergence_factor, commondata.time);
              Ck::IO::Options opts;
              CkCallback opened_2d_yz(CkIndex_Timestepping::ready_2d_yz(NULL), thisProxy);
              Ck::IO::open(filename, opened_2d_yz, opts);
            }
          }
"""
    # Generate code for 1d y diagnostics
    file_output_str += generate_diagnostics_code("1d", "y", "griddata[which_grid_diagnostics].diagnosticptoffsetstruct.num_output_quantities + 1", "griddata[which_grid_diagnostics].diagnosticptoffsetstruct.tot_num_diagnostic_1d_y_pts")

    # Generate code for 1d z diagnostics
    file_output_str += generate_diagnostics_code("1d", "z", "griddata[which_grid_diagnostics].diagnosticptoffsetstruct.num_output_quantities + 1", "griddata[which_grid_diagnostics].diagnosticptoffsetstruct.tot_num_diagnostic_1d_z_pts")

    # Generate code for 2d xy diagnostics
    file_output_str += generate_diagnostics_code("2d", "xy", "griddata[which_grid_diagnostics].diagnosticptoffsetstruct.num_output_quantities + 2", "griddata[which_grid_diagnostics].diagnosticptoffsetstruct.tot_num_diagnostic_2d_xy_pts")

    # Generate code for 2d yz diagnostics
    file_output_str += generate_diagnostics_code("2d", "yz", "griddata[which_grid_diagnostics].diagnosticptoffsetstruct.num_output_quantities + 2", "griddata[which_grid_diagnostics].diagnosticptoffsetstruct.tot_num_diagnostic_2d_yz_pts")

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
    for rk_substep, axis in itertools.product([f"RK_SUBSTEP_K{k}" for k in range(1, 5)], ['x', 'y', 'z']):
        if axis == 'x':
            pos_ghost_type = x_pos_ghost_type
            neg_ghost_type = x_neg_ghost_type
            nchare_var = "Nchare0"
        elif axis == 'y':
            pos_ghost_type = y_pos_ghost_type
            neg_ghost_type = y_neg_ghost_type
            nchare_var = "Nchare1"
        elif axis == 'z':
            pos_ghost_type = z_pos_ghost_type
            neg_ghost_type = z_neg_ghost_type
            nchare_var = "Nchare2"

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

        file_output_str += generate_mol_step_forward_code(rk_substep)
        file_output_str += generate_send_neighbor_data_code(which_gf, 'EAST_WEST')
        file_output_str += generate_ghost_code(axis, 0, pos_ghost_type, neg_ghost_type, nchare_var)
        file_output_str += generate_send_neighbor_data_code(which_gf, 'NORTH_SOUTH')
        file_output_str += generate_ghost_code(axis, 1, pos_ghost_type, neg_ghost_type, nchare_var)
        file_output_str += generate_send_neighbor_data_code(which_gf, 'TOP_BOTTOM')
        file_output_str += generate_ghost_code(axis, 2, pos_ghost_type, neg_ghost_type, nchare_var)

    file_output_str += r"""
        """


    if post_MoL_step_forward_in_time != "":
        file_output_str += post_MoL_step_forward_in_time
    else:
        file_output_str += "  // (nothing here; specify by setting post_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"

    file_output_str += r"""
        serial {
          // Adding dt to commondata->time many times will induce roundoff error,
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

def output_timestepping_h_cpp_ci(
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


    BHaH_defines_h.register_BHaH_defines(
    __name__,
    f"""
    #define K_ODD 0
    #define K_EVEN 1
    #define Y_N 2
    #define EAST_WEST 0
    #define NORTH_SOUTH 1
    #define TOP_BOTTOM 2
    #define {x_pos_ghost_type} 1
    #define {x_neg_ghost_type} 2
    #define {y_pos_ghost_type} 3
    #define {y_neg_ghost_type} 4
    #define {z_pos_ghost_type} 5
    #define {z_neg_ghost_type} 6
    #define RK_SUBSTEP_K1 1
    #define RK_SUBSTEP_K2 2
    #define RK_SUBSTEP_K3 3
    #define RK_SUBSTEP_K4 4
    """,
    )


