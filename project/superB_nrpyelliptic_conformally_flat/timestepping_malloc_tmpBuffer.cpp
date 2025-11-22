#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Allocate memory for temporary buffers used to communicate face data
 */
void timestepping_malloc_tmpBuffer(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                   const MoL_gridfunctions_struct *restrict gridfuncs, const nonlocalinnerbc_struct *restrict nonlocalinnerbcstruct,
                                   tmpBuffers_struct *restrict tmpBuffers) {
#include "set_CodeParameters.h"
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
} // END FUNCTION timestepping_malloc_tmpBuffer
