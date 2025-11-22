#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "main.h"
#include "timestepping.h"
//NEW
#include <vector>
#include <cstring>
extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Timestepping timesteppingArray;
/*
*Step 1.c: Allocate NUMGRIDS griddata arrays, each containing data specific to an individual grid.
*Step 1.d: Set each CodeParameter in griddata.params to default.
*Step 1.e: Set non-parfile parameters related to numerical grid, then set up numerical grids and CFL-limited timestep.
*Step 1.f: Set up boundary condition struct (bcstruct).
Step 2: Initial data are set on y_n_gfs gridfunctions. Allocate storage for them first.
*Step 3: Finalize initialization: set up  initial data, etc.
*Step 4: Allocate storage for non-y_n gridfunctions, needed for the Runge-Kutta-like timestepping.
*/

Timestepping::Timestepping(CommondataObject &&inData) {

  CkPrintf("Timestepping chare %d,%d,%d created on PE %d\n", thisIndex.x, thisIndex.y, thisIndex.z, CkMyPe());

  commondata = inData.commondata;

  if (thisIndex.x == 0 || thisIndex.y == 0 || thisIndex.z == 0 || thisIndex.x == (commondata.Nchare0 - 1) ||
      thisIndex.y == (commondata.Nchare1 - 1) || thisIndex.z == (commondata.Nchare2 - 1)) {
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
  for (int grid = 0; grid < commondata.NUMGRIDS; grid++) {
    // Define data needed for syncing gfs across chares
    MoL_sync_data_defines(&griddata_chare[grid].gridfuncs);
  }

  // Allocate storage for gridfunctions on each grid.
  for (int grid = 0; grid < commondata.NUMGRIDS; grid++) {
    const int Nxx_plus_2NGHOSTS_tot = (griddata[grid].params.Nxx_plus_2NGHOSTS0 * //
                                       griddata[grid].params.Nxx_plus_2NGHOSTS1 * //
                                       griddata[grid].params.Nxx_plus_2NGHOSTS2);

    BHAH_MALLOC(griddata_chare[grid].gridfuncs.y_n_gfs, sizeof(REAL) * Nxx_plus_2NGHOSTS_tot * NUM_EVOL_GFS);

    MoL_malloc_intermediate_stage_gfs(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].gridfuncs);

    if (NUM_AUXEVOL_GFS > 0) {
      BHAH_MALLOC(griddata_chare[grid].gridfuncs.auxevol_gfs, sizeof(REAL) * Nxx_plus_2NGHOSTS_tot * NUM_AUXEVOL_GFS);
    } // END IF NUM_AUXEVOL_GFS > 0

  } // END LOOP over grids

  // Initialize y n and non-y n gfs to NAN to avoid uninitialized memory errors
  for (int grid = 0; grid < commondata.NUMGRIDS; grid++)
    initialize_yn_and_non_yn_gfs_to_nan(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].gridfuncs);

  // Allocate storage for temporary buffers, needed for communicating face data
  for (int grid = 0; grid < commondata.NUMGRIDS; grid++)
    timestepping_malloc_tmpBuffer(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].gridfuncs,
                                  &griddata_chare[grid].nonlocalinnerbcstruct, &griddata_chare[grid].tmpBuffers);
}

// migration constructor
Timestepping::Timestepping(CkMigrateMessage *msg) : CBase_Timestepping(msg) {}

// destructor
Timestepping::~Timestepping() {

  griddata_free(&commondata, griddata, free_non_y_n_gfs_and_core_griddata_pointers);

  for (int grid = 0; grid < commondata.NUMGRIDS; grid++) {
    timestepping_free_memory_tmpBuffer(&griddata_chare[grid].nonlocalinnerbcstruct, &griddata_chare[grid].tmpBuffers);
    free(griddata_chare[grid].bcstruct.inner_bc_array);
    free(griddata_chare[grid].bcstruct.inner_bc_array_nonlocal);
    for (int ng = 0; ng < NGHOSTS * 3; ng++) {
      free(griddata_chare[grid].bcstruct.pure_outer_bc_array[ng]);
    }
    for (int i = 0; i < 3; i++) {
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
  const int *gfs_to_sync = nullptr;

  switch (type_gfs) {
  case Y_NPLUS1_RUNNING_TOTAL_GFS:
    gfs = griddata_chare[grid].gridfuncs.y_nplus1_running_total_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    break;
  case K_ODD_GFS:
    gfs = griddata_chare[grid].gridfuncs.k_odd_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    break;
  case K_EVEN_GFS:
    gfs = griddata_chare[grid].gridfuncs.k_even_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    break;
  case Y_N_GFS:
  case Y_N_GFS_INITIALDATA_PART1:
  case Y_N_GFS_INITIALDATA_PART2:
    gfs = griddata_chare[grid].gridfuncs.y_n_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    break;
  case AUXEVOL_GFS:
    gfs = griddata_chare[grid].gridfuncs.auxevol_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.auxevol_gfs_to_sync;
    break;
  case DIAGNOSTIC_OUTPUT_GFS:
    gfs = griddata_chare[grid].gridfuncs.diagnostic_output_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_aux_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.aux_gfs_to_sync;
    break;

  default:
    break;
  }

  switch (dir) {
  case EAST_WEST:
    // send to west
    if (thisIndex.x > 0) {
      for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
        int i0 = 2 * NGHOSTS - 1;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
            for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
              tmpBuffer_EW[IDXFACES0(which_gf, which_inner, i1, i2)] = gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)];
            }
          }
          i0--;
        }
      }
      thisProxy[CkArrayIndex3D(thisIndex.x - 1, thisIndex.y, thisIndex.z)].east_ghost(
          type_gfs, NUM_GFS * NGHOSTS * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2, tmpBuffer_EW);
    }
    // send to east
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
      thisProxy[CkArrayIndex3D(thisIndex.x + 1, thisIndex.y, thisIndex.z)].west_ghost(
          type_gfs, NUM_GFS * NGHOSTS * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2, tmpBuffer_EW);
    }
    break;
  case NORTH_SOUTH:
    // send to south
    if (thisIndex.y > 0) {
      for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
        int i1 = 2 * NGHOSTS - 1;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
              tmpBuffer_NS[IDXFACES1(which_gf, which_inner, i0, i2)] = gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)];
            }
          }
          i1--;
        }
      }
      thisProxy[CkArrayIndex3D(thisIndex.x, thisIndex.y - 1, thisIndex.z)].north_ghost(
          type_gfs, NUM_GFS * NGHOSTS * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS2, tmpBuffer_NS);
    }
    // send to north
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
      thisProxy[CkArrayIndex3D(thisIndex.x, thisIndex.y + 1, thisIndex.z)].south_ghost(
          type_gfs, NUM_GFS * NGHOSTS * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS2, tmpBuffer_NS);
    }
    break;
  case TOP_BOTTOM:
    // send to bottom
    if (thisIndex.z > 0) {
      for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
        int i2 = 2 * NGHOSTS - 1;
        for (int which_inner = 0; which_inner < NGHOSTS; which_inner++) {
          for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
              tmpBuffer_TB[IDXFACES2(which_gf, which_inner, i0, i1)] = gfs[IDX4(gfs_to_sync[which_gf], i0, i1, i2)];
            }
          }
          i2--;
        }
      }
      thisProxy[CkArrayIndex3D(thisIndex.x, thisIndex.y, thisIndex.z - 1)].top_ghost(
          type_gfs, NUM_GFS * NGHOSTS * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1, tmpBuffer_TB);
    }
    // send to top
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
      thisProxy[CkArrayIndex3D(thisIndex.x, thisIndex.y, thisIndex.z + 1)].bottom_ghost(
          type_gfs, NUM_GFS * NGHOSTS * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1, tmpBuffer_TB);
    }
    break;
  default:
    return;
  }
}

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
  const int *gfs_to_sync = nullptr;

  switch (type_gfs) {
  case Y_NPLUS1_RUNNING_TOTAL_GFS:
    gfs = griddata_chare[grid].gridfuncs.y_nplus1_running_total_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    break;
  case K_ODD_GFS:
    gfs = griddata_chare[grid].gridfuncs.k_odd_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    break;
  case K_EVEN_GFS:
    gfs = griddata_chare[grid].gridfuncs.k_even_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    break;
  case Y_N_GFS:
  case Y_N_GFS_INITIALDATA_PART1:
  case Y_N_GFS_INITIALDATA_PART2:
    gfs = griddata_chare[grid].gridfuncs.y_n_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    break;
  case AUXEVOL_GFS:
    gfs = griddata_chare[grid].gridfuncs.auxevol_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.auxevol_gfs_to_sync;
    break;
  case DIAGNOSTIC_OUTPUT_GFS:
    gfs = griddata_chare[grid].gridfuncs.diagnostic_output_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_aux_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.aux_gfs_to_sync;
    break;

  default:
    break;
  }

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

//OLD
//~ void Timestepping::contribute_localsums_for_residualH(REAL localsums_for_residualH[2]) {
  //~ std::vector<double> outdoubles(2);
  //~ outdoubles[0] = localsums_for_residualH[0];
  //~ outdoubles[1] = localsums_for_residualH[1];
  //~ CkCallback cb(CkIndex_Timestepping::report_sums_for_residualH(NULL), thisProxy);
  //~ contribute(outdoubles, CkReduction::sum_double, cb);
//~ }

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
  const int pt_at_outer_boundary_i2 = Nxx_plus_2NGHOSTS2 / 2;
  const int pt_at_outer_boundary_idx3 = IDX3(pt_at_outer_boundary_i0, pt_at_outer_boundary_i1, pt_at_outer_boundary_i2);
  const int idx3_this_chare = IDX3_OF_CHARE(thisIndex.x, thisIndex.y, thisIndex.z);

  if (griddata_chare[grid].charecommstruct.globalidx3pt_to_chareidx3[pt_at_outer_boundary_idx3] == idx3_this_chare) {
    const int locali0 = MAP_GLOBAL_TO_LOCAL_IDX0(thisIndex.x, pt_at_outer_boundary_i0, Nxx0chare);
    const int locali1 = MAP_GLOBAL_TO_LOCAL_IDX1(thisIndex.y, pt_at_outer_boundary_i1, Nxx1chare);
    const int locali2 = MAP_GLOBAL_TO_LOCAL_IDX2(thisIndex.z, pt_at_outer_boundary_i2, Nxx2chare);
    const REAL wavespeed_at_outer_boundary = griddata_chare[grid].gridfuncs.auxevol_gfs[IDX4GENERAL(
        VARIABLE_WAVESPEEDGF, locali0, locali1, locali2, Nxx_plus_2NGHOSTS0chare, Nxx_plus_2NGHOSTS1chare, Nxx_plus_2NGHOSTS2chare)];
    thisProxy.receiv_wavespeed_at_outer_boundary(wavespeed_at_outer_boundary);
  }
}

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
    int num_srcpts = num_srcpts_each_chare[src_chare_id];
    thisProxy[CkArrayIndex3D(src_chare_index0, src_chare_index1, src_chare_index2)].receiv_nonlocalinnerbc_idx3srcpt_tosend(
        idx3_of_thischare, num_srcpts, globalidx3_srcpts[src_chare_id]);
  }
}

void Timestepping::process_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, const int *restrict globalidx3_srcpts) {
  // Unpack griddata_chare[grid].nonlocalinnerbcstruct
  int grid = 0;
  const int *restrict idx3chare_to_dst_chare_id = griddata_chare[grid].nonlocalinnerbcstruct.idx3chare_to_dst_chare_id;
  int **restrict globalidx3_srcpts_tosend = griddata_chare[grid].nonlocalinnerbcstruct.globalidx3_srcpts_tosend;

  const int dst_chare_id_val = idx3chare_to_dst_chare_id[idx3_of_sendingchare];
  for (int count = 0; count < num_srcpts; count++) {
    globalidx3_srcpts_tosend[dst_chare_id_val][count] = globalidx3_srcpts[count];
  }
}

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
  const int *gfs_to_sync = nullptr;

  switch (type_gfs) {
  case Y_NPLUS1_RUNNING_TOTAL_GFS:
    gfs = griddata_chare[grid].gridfuncs.y_nplus1_running_total_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    break;
  case K_ODD_GFS:
    gfs = griddata_chare[grid].gridfuncs.k_odd_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    break;
  case K_EVEN_GFS:
    gfs = griddata_chare[grid].gridfuncs.k_even_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    break;
  case Y_N_GFS:
  case Y_N_GFS_INITIALDATA_PART1:
  case Y_N_GFS_INITIALDATA_PART2:
    gfs = griddata_chare[grid].gridfuncs.y_n_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    break;
  case AUXEVOL_GFS:
    gfs = griddata_chare[grid].gridfuncs.auxevol_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.auxevol_gfs_to_sync;
    break;
  case DIAGNOSTIC_OUTPUT_GFS:
    gfs = griddata_chare[grid].gridfuncs.diagnostic_output_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_aux_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.aux_gfs_to_sync;
    break;

  default:
    break;
  }

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
    REVERSE_IDX3GENERAL(idx3_of_dst_chares[which_dst_chare], Nchare0, Nchare1, dst_chare_index0, dst_chare_index1, dst_chare_index2);
    switch (type_gfs) {
    case Y_NPLUS1_RUNNING_TOTAL_GFS:
      thisProxy[CkArrayIndex3D(dst_chare_index0, dst_chare_index1, dst_chare_index2)].receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(
          idx3_this_chare, type_gfs, NUM_GFS * num_srcpts_tosend_each_chare[which_dst_chare], tmpBuffer_innerbc_send);
      break;
    case K_ODD_GFS:
      thisProxy[CkArrayIndex3D(dst_chare_index0, dst_chare_index1, dst_chare_index2)].receiv_nonlocalinnerbc_data_k_odd_gfs(
          idx3_this_chare, type_gfs, NUM_GFS * num_srcpts_tosend_each_chare[which_dst_chare], tmpBuffer_innerbc_send);
      break;
    case K_EVEN_GFS:
      thisProxy[CkArrayIndex3D(dst_chare_index0, dst_chare_index1, dst_chare_index2)].receiv_nonlocalinnerbc_data_k_even_gfs(
          idx3_this_chare, type_gfs, NUM_GFS * num_srcpts_tosend_each_chare[which_dst_chare], tmpBuffer_innerbc_send);
      break;
    case Y_N_GFS:
      thisProxy[CkArrayIndex3D(dst_chare_index0, dst_chare_index1, dst_chare_index2)].receiv_nonlocalinnerbc_data_y_n_gfs(
          idx3_this_chare, type_gfs, NUM_GFS * num_srcpts_tosend_each_chare[which_dst_chare], tmpBuffer_innerbc_send);
      break;
    case AUXEVOL_GFS:
      thisProxy[CkArrayIndex3D(dst_chare_index0, dst_chare_index1, dst_chare_index2)].receiv_nonlocalinnerbc_data_auxevol_gfs(
          idx3_this_chare, type_gfs, NUM_GFS * num_srcpts_tosend_each_chare[which_dst_chare], tmpBuffer_innerbc_send);
      break;
    case DIAGNOSTIC_OUTPUT_GFS:
      thisProxy[CkArrayIndex3D(dst_chare_index0, dst_chare_index1, dst_chare_index2)].receiv_nonlocalinnerbc_data_diagnostic_output_gfs(
          idx3_this_chare, type_gfs, NUM_GFS * num_srcpts_tosend_each_chare[which_dst_chare], tmpBuffer_innerbc_send);
      break;
    case Y_N_GFS_INITIALDATA_PART1:
      thisProxy[CkArrayIndex3D(dst_chare_index0, dst_chare_index1, dst_chare_index2)].receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(
          idx3_this_chare, type_gfs, NUM_GFS * num_srcpts_tosend_each_chare[which_dst_chare], tmpBuffer_innerbc_send);
      break;
    case Y_N_GFS_INITIALDATA_PART2:
      thisProxy[CkArrayIndex3D(dst_chare_index0, dst_chare_index1, dst_chare_index2)].receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(
          idx3_this_chare, type_gfs, NUM_GFS * num_srcpts_tosend_each_chare[which_dst_chare], tmpBuffer_innerbc_send);
      break;

    default:
      break;
    }
  }
}

void Timestepping::set_tmpBuffer_innerbc_receiv(const int src_chare_idx3, const int len_tmpBuffer, const REAL *restrict vals, const int grid) {
  const int src_chare_id = griddata_chare[grid].nonlocalinnerbcstruct.idx3chare_to_src_chare_id[src_chare_idx3];
  REAL *restrict tmpBuffer_innerbc_receiv = griddata_chare[grid].tmpBuffers.tmpBuffer_innerbc_receiv[src_chare_id];

  for (int i = 0; i < len_tmpBuffer; i++) {
    tmpBuffer_innerbc_receiv[i] = vals[i];
  }
}

void Timestepping::process_nonlocalinnerbc(const int type_gfs, const int grid) {
  REAL *restrict gfs = nullptr;
  int NUM_GFS;
  const int *gfs_to_sync = nullptr;
  const int8_t *gf_parity_types = nullptr;

  switch (type_gfs) {
  case Y_NPLUS1_RUNNING_TOTAL_GFS:
    gfs = griddata_chare[grid].gridfuncs.y_nplus1_running_total_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    gf_parity_types = evol_gf_parity;
    break;
  case K_ODD_GFS:
    gfs = griddata_chare[grid].gridfuncs.k_odd_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    gf_parity_types = evol_gf_parity;
    break;
  case K_EVEN_GFS:
    gfs = griddata_chare[grid].gridfuncs.k_even_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    gf_parity_types = evol_gf_parity;
    break;
  case Y_N_GFS:
  case Y_N_GFS_INITIALDATA_PART1:
  case Y_N_GFS_INITIALDATA_PART2:
    gfs = griddata_chare[grid].gridfuncs.y_n_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_evol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.evol_gfs_to_sync;
    gf_parity_types = evol_gf_parity;
    break;
  case AUXEVOL_GFS:
    gfs = griddata_chare[grid].gridfuncs.auxevol_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.auxevol_gfs_to_sync;
    gf_parity_types = auxevol_gf_parity;
    break;
  case DIAGNOSTIC_OUTPUT_GFS:
    gfs = griddata_chare[grid].gridfuncs.diagnostic_output_gfs;
    NUM_GFS = griddata_chare[grid].gridfuncs.num_aux_gfs_to_sync;
    gfs_to_sync = griddata_chare[grid].gridfuncs.aux_gfs_to_sync;
    gf_parity_types = aux_gf_parity;
    break;

  default:
    break;
  }

  apply_bcs_inner_only_nonlocal(&commondata, &griddata_chare[grid].params, &griddata_chare[grid].bcstruct,
                                &griddata_chare[grid].nonlocalinnerbcstruct, NUM_GFS, gfs, gfs_to_sync, gf_parity_types,
                                griddata_chare[grid].tmpBuffers.tmpBuffer_innerbc_receiv);
}


//NEW
//~ void Timestepping::contribute_integration_sums(const int which_grid) {
  //~ diags_integration_results_t &integration_results =
      //~ griddata[which_grid].diagnosticstruct.integration_results;

  //~ const int recipe_count = integration_results.num_recipe_results;
  //~ if (recipe_count == 0) {
    //~ std::vector<double> empty_payload(0);
    //~ CkCallback cb(CkReductionTarget(Timestepping, receive_integration_sums), timesteppingProxy[0]);
    //~ contribute(empty_payload, CkReduction::sum_double, cb);
    //~ return;
  //~ }
  //~ const int integrand_count = integration_results.recipe_results[0].num_integrand_results;

  //~ // Total doubles to reduce per chare: R * (1 + M)
  //~ const int doubles_per_recipe = 1 + integrand_count; // 1 volume + M integrals
  //~ const int total_doubles      = recipe_count * doubles_per_recipe;

  //~ // Flat layout per recipe: [proper_volume, integral(0), ..., integral(M-1)]
  //~ std::vector<double> reduction_values(total_doubles);
  //~ int write_index = 0;

  //~ for (int r = 0; r < recipe_count; ++r) {
    //~ const diags_integration_recipe_result_t &recipe_result = integration_results.recipe_results[r];
    //~ reduction_values[write_index++] = recipe_result.proper_volume;

    //~ for (int k = 0; k < integrand_count; ++k) {
      //~ const diags_integration_integrand_result_t &integrand_result = recipe_result.integrand_results[k];
      //~ reduction_values[write_index++] = integrand_result.integral;
    //~ }
  //~ }

  //~ CkCallback cb(CkReductionTarget(Timestepping, receive_integration_sums), timesteppingProxy[0]);
  //~ contribute(reduction_values, CkReduction::sum_double, cb);
//~ }


//NEW
//~ void Timestepping::receive_integration_sums(CkReductionMsg* msg) {
  //~ const int which_grid = 0; // set appropriately for your code path
  //~ diags_integration_results_t &integration_results =
      //~ griddata[which_grid].diagnosticstruct.integration_results;

  //~ const int recipe_count = integration_results.num_recipe_results;
  //~ if (recipe_count == 0) { delete msg; return; }

  //~ const int integrand_count = integration_results.recipe_results[0].num_integrand_results;

  //~ const int doubles_per_recipe = 1 + integrand_count;
  //~ const int expected_total_doubles = recipe_count * doubles_per_recipe;

  //~ const int received_bytes   = msg->getSize();
  //~ const int received_doubles = received_bytes / (int)sizeof(double);

  //~ if (received_doubles != expected_total_doubles) {
    //~ CkAbort("receive_integration_sums: payload size mismatch.");
  //~ }

  //~ std::vector<double> sums(received_doubles);
  //~ std::memcpy(sums.data(), msg->getData(), received_bytes);

  //~ // Scatter global sums back into the container
  //~ int read_index = 0;
  //~ for (int r = 0; r < recipe_count; ++r) {
    //~ diags_integration_recipe_result_t &recipe_result = integration_results.recipe_results[r];

    //~ recipe_result.proper_volume = sums[read_index++];

    //~ for (int k = 0; k < integrand_count; ++k) {
      //~ diags_integration_integrand_result_t &integrand_result = recipe_result.integrand_results[k];
      //~ integrand_result.integral = sums[read_index++];
    //~ }
  //~ }

  //~ delete msg;

  //~ thisProxy[0,0,0].integration_results_set();

//~ }


#include "timestepping.def.h"
