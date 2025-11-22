#include "BHaH_defines.h"

/**
 * Apply BCs to nonlocal inner boundary points only,
 * using data stored in bcstruct->inner_bc_array.
 * These structs are set in bcstruct_set_up().
 * Inner boundary points map to either the grid
 * interior ("pure inner") or to pure outer
 * boundary points ("inner maps to outer").
 */
void apply_bcs_inner_only_nonlocal(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                   const bc_struct *restrict bcstruct, const nonlocalinnerbc_struct *restrict nonlocalinnerbcstruct,
                                   const int NUM_GFS, REAL *restrict gfs, const int *gfs_to_sync, const int8_t *gf_parity_types,
                                   REAL **restrict tmpBuffer_innerbc_receiv) {
#include "set_CodeParameters.h"
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

  // Unpack nonlocalinnerbcstruct
  const int tot_num_src_chares = nonlocalinnerbcstruct->tot_num_src_chares;
  const int *restrict num_srcpts_each_chare = nonlocalinnerbcstruct->num_srcpts_each_chare;
  int **restrict map_srcchare_and_srcpt_id_to_linear_id = nonlocalinnerbcstruct->map_srcchare_and_srcpt_id_to_linear_id;

  for (int which_chare = 0; which_chare < tot_num_src_chares; which_chare++) {
    const REAL *restrict tmpBuffer = tmpBuffer_innerbc_receiv[which_chare];
    for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
      for (int which_srcpt = 0; which_srcpt < num_srcpts_each_chare[which_chare]; which_srcpt++) {
        const int linear_id = map_srcchare_and_srcpt_id_to_linear_id[which_chare][which_srcpt];
        const int dstpt = bcstruct->inner_bc_array_nonlocal[linear_id].dstpt;
        const int idx2 = IDX2NONLOCALINNERBC(which_gf, which_srcpt, num_srcpts_each_chare[which_chare]);
        gfs[IDX4pt(gfs_to_sync[which_gf], dstpt)] =
            bcstruct->inner_bc_array_nonlocal[linear_id].parity[gf_parity_types[gfs_to_sync[which_gf]]] * tmpBuffer[idx2];
      }
    }
  }
} // END FUNCTION apply_bcs_inner_only_nonlocal
