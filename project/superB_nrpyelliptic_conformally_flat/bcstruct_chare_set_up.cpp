#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Setup bcstruct_chare from bcstruct
 */
void bcstruct_chare_set_up(const commondata_struct *restrict commondata, const params_struct *restrict params,
                           const params_struct *restrict params_chare, const charecomm_struct *restrict charecommstruct, REAL *restrict xx[3],
                           const bc_struct *restrict bcstruct, bc_struct *restrict bcstruct_chare,
                           nonlocalinnerbc_struct *restrict nonlocalinnerbcstruct, const int chare_index[3]) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    bcstruct_chare_set_up__rfm__SinhSpherical(commondata, params, params_chare, charecommstruct, xx, bcstruct, bcstruct_chare, nonlocalinnerbcstruct,
                                              chare_index);
    break;
  default:
    fprintf(stderr, "ERROR in bcstruct_chare_set_up(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION bcstruct_chare_set_up
