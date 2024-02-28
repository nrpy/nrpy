"""

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from typing import List, Tuple, Optional, Dict
import nrpy.c_function as cfc
import sys
from pathlib import Path
import nrpy.params as par
import nrpy.grid as gri
from nrpy.infrastructures.BHaH import griddata_commondata
from nrpy.helpers.generic import clang_format


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

    file_output_str = """#ifndef __TIMESTEPPING_H__
#define __TIMESTEPPING_H__

#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "timestepping.decl.h"

const int K_ODD = 0;
const int K_EVEN = 1;
const int Y_N = 2;

const int EAST_WEST_DIR = 0;
const int NORTH_SOUTH_DIR = 1;
const int TOP_BOTTOM_DIR = 2;

const int EAST_GHOST = 0;
const int WEST_GHOST = 1;
const int NORTH_GHOST = 2:
const int SOUTH_GHOST = 3;
const int TOP_GHOST = 4;
const int BOTTOM_GHOST = 5;

const int RK_SUBSTEP_K1 = 1;
const int RK_SUBSTEP_K2 = 2;
const int RK_SUBSTEP_K3 = 3;
const int RK_SUBSTEP_K4 = 4;

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
    Ck::IO::File f_1d_y;
    Ck::IO::File f_1d_z;
    Ck::IO::File f_2d_xy;
    Ck::IO::File f_2d_yz;
    bool bwrite_diagnostics;

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
    pre_MoL_step_forward_in_time: str = "",
    post_MoL_step_forward_in_time: str = "",
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


# ~ def output_timestepping_ci(
    # ~ MoL_method: str,
    # ~ initial_data_desc: str = "",
    # ~ enable_rfm_precompute: bool = False,
    # ~ enable_CurviBCs: bool = False,
    # ~ boundary_conditions_desc: str = "",
    # ~ prefunc: str = "",
    # ~ initialize_constant_auxevol: bool = False,
    # ~ pre_MoL_step_forward_in_time: str = "",
    # ~ post_MoL_step_forward_in_time: str = "",
    # ~ clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
# ~ ) -> None:
    # ~ """
    # ~ Generate timestepping.cpp

    # ~ :param MoL_method: Method of Lines algorithm used to step forward in time.
    # ~ :param initial_data_desc: Description for initial data, default is an empty string.
    # ~ :param enable_rfm_precompute: Enable rfm precomputation, default is False.
    # ~ :param enable_CurviBCs: Enable CurviBCs, default is False.
    # ~ :param boundary_conditions_desc: Description of the boundary conditions, default is an empty string.
    # ~ :param prefunc: String that appears before main(). DO NOT populate this, except when debugging, default is an empty string.
    # ~ :param initialize_constant_auxevol: If set to True, `initialize_constant_auxevol` function will be called during the simulation initialization phase to set these constants. Default is False.
    # ~ :param pre_MoL_step_forward_in_time: Code for handling pre-right-hand-side operations, default is an empty string.
    # ~ :param post_MoL_step_forward_in_time: Code for handling post-right-hand-side operations, default is an empty string.
    # ~ :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".
    # ~ """
    # ~ initial_data_desc += " "
    # ~ # Make sure all required C functions are registered
    # ~ missing_functions: List[Tuple[str, str]] = []
    # ~ for func_tuple in [
        # ~ ("params_struct_set_to_default", "CodeParameters.py"),
        # ~ (
            # ~ "numerical_grids_and_timestep",
            # ~ "e.g., numerical_grids_and_timestep.py or user defined",
        # ~ ),
        # ~ ("MoL_malloc_y_n_gfs", "MoL.py"),
        # ~ ("MoL_malloc_non_y_n_gfs", "MoL.py"),
        # ~ ("initial_data", "initial_data.py"),
        # ~ ("MoL_step_forward_in_time", "MoL.py"),
        # ~ ("diagnostics", "diagnostics.py"),
        # ~ ("MoL_free_memory_y_n_gfs", "MoL.py"),
        # ~ ("MoL_free_memory_non_y_n_gfs", "MoL.py"),
    # ~ ]:
        # ~ if func_tuple[0] not in cfc.CFunction_dict:
            # ~ missing_functions += [func_tuple]
    # ~ if missing_functions:
        # ~ error_msg = "Error: These functions are required and are not registered.\n"
        # ~ for func_tuple in missing_functions:
            # ~ error_msg += (
                # ~ f'  {func_tuple[0]}, registered by function within "{func_tuple[1]}"\n'
            # ~ )
        # ~ raise ValueError(error_msg)


    # ~ project_Path = Path(project_dir)
    # ~ project_Path.mkdir(parents=True, exist_ok=True)

    # ~ file_output_str = """#include "BHaH_defines.h"
# ~ #include "BHaH_defines.h"
# ~ #include "BHaH_function_prototypes.h"
# ~ #include "timestepping.h"
# ~ #include "main.h"

# ~ extern /* readonly */ CProxy_Main mainProxy;

# ~ Step 1.c: Allocate NUMGRIDS griddata arrays, each containing data specific to an individual grid.
# ~ Step 1.d: Set each CodeParameter in griddata.params to default.
# ~ """
    # ~ if enable_CurviBCs:
        # ~ file_output_str += "Step 1.e: Set non-parfile parameters related to numerical grid, then set up numerical grids and CFL-limited timestep.\n"
    # ~ if enable_rfm_precompute:
        # ~ file_output_str += "Step 1.f: Set up boundary condition struct (bcstruct).\n"
    # ~ file_output_str += f"""Step 2: Initial data are set on y_n_gfs gridfunctions. Allocate storage for them first.
# ~ Step 3: Finalize initialization: set up {initial_data_desc}initial data, etc.
# ~ Step 4: Allocate storage for non-y_n gridfunctions, needed for the Runge-Kutta-like timestepping.
# ~ """
    # ~ if initialize_constant_auxevol:
        # ~ file_output_str += "Step 4.a: Set AUXEVOL gridfunctions that will never change in time."
    # ~ file_output_str += f"""
# ~ Step 5: MAIN SIMULATION LOOP
# ~ - Step 5.a: Output diagnostics.
# ~ - Step 5.b: Prepare to step forward in time.
# ~ - Step 5.c: Step forward in time using Method of Lines with {MoL_method} algorithm, applying {boundary_conditions_desc} boundary conditions.
# ~ - Step 5.d: Finish up step in time.
# ~ Step 6: Free all allocated memory."""    
    # ~ file_output_str += r"""  
    
# ~ Timestepping::Timestepping(CommondataObject &&inData) {

  # ~ if (thisIndex.x == 0 || 
      # ~ thisIndex.y == 0 || 
      # ~ thisIndex.z == 0 || 
      # ~ thisIndex.x == (Nchare0-1) || 
      # ~ thisIndex.y == (Nchare1-1) || 
      # ~ thisIndex.z == (Nchare2-1)) {
      # ~ is_boundarychare = true;
    # ~ }
    
  # ~ commondata = inData.commondata;

  # ~ // Step 1.c: Allocate NUMGRIDS griddata arrays, each containing data specific to an individual grid.
  # ~ griddata = (griddata_struct *restrict)malloc(sizeof(griddata_struct) * commondata.NUMGRIDS);  
  # ~ griddata_chare = (griddata_struct *restrict)malloc(sizeof(griddata_struct) * commondata.NUMGRIDS);

  # ~ // Step 1.d: Set each CodeParameter in griddata.params to default.
  # ~ params_struct_set_to_default(&commondata, griddata);
  # ~ params_struct_set_to_default(&commondata, griddata_chare);

  # ~ // Step 1.e: Set up numerical grids: xx[3], masks, Nxx, dxx, invdxx, bcstruct, rfm_precompute, timestep, etc.
  # ~ {
    # ~ // if calling_for_first_time, then initialize commondata time=nn=t_0=nn_0 = 0
    # ~ const bool calling_for_first_time = true;
    # ~ numerical_grids_and_timestep(&commondata, griddata, calling_for_first_time);
    # ~ numerical_grids_chare(&commondata, griddata, griddata_chare, {thisIndex.x, thisIndex.y, thisIndex.z});
  # ~ }

  # ~ for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
    # ~ // Step 2: Initial data are set on y_n_gfs gridfunctions. Allocate storage for them first.
    # ~ MoL_malloc_y_n_gfs(&commondata, griddata_chare[grid].params, griddata_chare[grid].gridfuncs);
  # ~ }

  # ~ // Step 3: Finalize initialization: set up initial data, etc.
  # ~ initial_data(&commondata, griddata_chare);

  # ~ // Step 4: Allocate storage for non-y_n gridfunctions, needed for the Runge-Kutta-like timestepping
  # ~ for(int grid=0; grid<commondata.NUMGRIDS; grid++)
    # ~ MoL_malloc_non_y_n_gfs(&commondata, griddata_chare[grid].params, griddata_chare[grid].gridfuncs);
# ~ """
    # ~ if initialize_constant_auxevol:
        # ~ file_output_str += """
  # ~ // Step 4.a: Set AUXEVOL gridfunctions that will never change in time.
  # ~ initialize_constant_auxevol(&commondata, griddata_chare);
# ~ """
    # ~ body += """
# ~ // Step 5: MAIN SIMULATION LOOP
# ~ while(commondata.time < commondata.t_final) { // Main loop to progress forward in time.
  # ~ // Step 5.a: Main loop, part 1: Output diagnostics
  # ~ diagnostics(&commondata, griddata);

  # ~ // Step 5.b: Main loop, part 2 (pre_MoL_step_forward_in_time): Prepare to step forward in time
# ~ """
    # ~ if pre_MoL_step_forward_in_time != "":
        # ~ body += pre_MoL_step_forward_in_time
    # ~ else:
        # ~ body += "  // (nothing here; specify by setting pre_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"
    # ~ body += f"""
  # ~ // Step 5.c: Main loop, part 3: Step forward in time using Method of Lines with {MoL_method} algorithm,
  # ~ //           applying {boundary_conditions_desc} boundary conditions.
  # ~ MoL_step_forward_in_time(&commondata, griddata);

  # ~ // Step 5.d: Main loop, part 4 (post_MoL_step_forward_in_time): Finish up step in time
# ~ """
    # ~ if post_MoL_step_forward_in_time != "":
        # ~ body += post_MoL_step_forward_in_time
    # ~ else:
        # ~ body += "  // (nothing here; specify by setting post_MoL_step_forward_in_time string in register_CFunction_main_c().)\n"
    # ~ body += r"""
# ~ } // End main loop to progress forward in time.

# ~ // Step 5: Free all allocated memory
# ~ for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
  # ~ MoL_free_memory_y_n_gfs(&griddata[grid].gridfuncs);
  # ~ MoL_free_memory_non_y_n_gfs(&griddata[grid].gridfuncs);"""
    # ~ if enable_rfm_precompute:
        # ~ body += r"""
  # ~ rfm_precompute_free(&commondata, &griddata[grid].params, &griddata[grid].rfmstruct);"""
    # ~ if enable_CurviBCs:
        # ~ body += r"""
  # ~ free(griddata[grid].bcstruct.inner_bc_array);
  # ~ for(int ng=0;ng<NGHOSTS*3;ng++) free(griddata[grid].bcstruct.pure_outer_bc_array[ng]);
# ~ """
    # ~ body += r"""
  # ~ for(int i=0;i<3;i++) free(griddata[grid].xx[i]);
# ~ }
# ~ free(griddata);
# ~ return 0;
# ~ """

    # ~ timestepping_cpp_file = project_Path / "timestepping.cpp"
    # ~ with timestepping_cpp_file.open("w", encoding="utf-8") as file:
        # ~ file.write(
            # ~ clang_format(file_output_str, clang_format_options=clang_format_options)
        # ~ )



