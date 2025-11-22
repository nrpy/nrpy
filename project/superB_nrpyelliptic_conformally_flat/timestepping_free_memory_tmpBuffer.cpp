#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Free memory for temporary buffers used to communicate face data
 */
void timestepping_free_memory_tmpBuffer(const nonlocalinnerbc_struct *restrict nonlocalinnerbcstruct, tmpBuffers_struct *restrict tmpBuffers) {
  if (tmpBuffers->tmpBuffer_EW != NULL)
    free(tmpBuffers->tmpBuffer_EW);
  if (tmpBuffers->tmpBuffer_NS != NULL)
    free(tmpBuffers->tmpBuffer_NS);
  if (tmpBuffers->tmpBuffer_TB != NULL)
    free(tmpBuffers->tmpBuffer_TB);

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
} // END FUNCTION timestepping_free_memory_tmpBuffer
