#ifndef __TIMESTEPPING_H__
#define __TIMESTEPPING_H__
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "superB/superB_pup_function_prototypes.h"
#include "timestepping.decl.h"
#include "diagnostics/diagnostics_nearest_common.h"


class Timestepping : public CBase_Timestepping {
Timestepping_SDAG_CODE

    private :
    /// Member Variables (Object State) ///
    commondata_struct commondata;
  griddata_struct *griddata;
  griddata_struct *griddata_chare;
  bool is_boundarychare;
  REAL time_start;
  // bool contains_gridcenter;
  const int grid = 0;
  const int which_grid_diagnostics = 0;
  const bool free_non_y_n_gfs_and_core_griddata_pointers = true;
  bool write_diagnostics_this_step;
  Ck::IO::File f_1d_y;
  Ck::IO::File f_1d_z;
  Ck::IO::File f_2d_xy;
  Ck::IO::File f_2d_yz;
  int count_filewritten = 0;
  const int expected_count_filewritten = 4;
  int iter = 0;
  int type_gfs_nonlocal_innerbc;
  //NEW
  REAL *diagnostic_gfs[MAXNUMGRIDS];
  int NUM_RECIPES;
  //~ diags_integration_recipe_t recipes[DIAGS_INTEGRATION_MAX_RECIPES];
  //~ diags_integration_results_t integration_results;

  /// Member Functions (private) ///
  void send_neighbor_data(const int type_gfs, const int dir, const int grid);
  void process_ghost(const int type_ghost, const int type_gfs, const int len_tmpBuffer, const REAL *restrict tmpBuffer, const int grid);
  void send_nonlocalinnerbc_idx3srcpts_toreceiv();
  void process_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, const int *restrict globalidx3_srcpts);
  void send_nonlocalinnerbc_data(const int type_gfs, const int grid);
  void set_tmpBuffer_innerbc_receiv(const int src_chare_idx3, const int len_tmpBuffer, const REAL *restrict vals, const int grid);
  void process_nonlocalinnerbc(const int type_gfs, const int grid);
  void contribute_localsums_for_residualH(REAL localsums_for_residualH[2]);
  void send_wavespeed_at_outer_boundary(const int grid);
  //NEW
  //~ void contribute_integration_sums(const int which_grid);


public:
  /// Constructors ///
  Timestepping(CommondataObject &&inData);
  Timestepping(CkMigrateMessage *msg);
  /// Destructor ///
  ~Timestepping();
};

#endif //__TIMESTEPPING_H__
