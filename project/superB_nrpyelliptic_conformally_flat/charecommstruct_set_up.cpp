#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Setup charecommstruct
 */
void charecommstruct_set_up(const commondata_struct *restrict commondata, const params_struct *restrict params,
                            const params_struct *restrict params_chare, charecomm_struct *charecommstruct, const int thischareindex[3]) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    printf("testtttttttttttttttttttttttttt\n");
    charecommstruct_set_up__rfm__SinhSpherical(commondata, params, params_chare, charecommstruct, thischareindex);
    break;
  default:
    fprintf(stderr, "ERROR in charecommstruct_set_up(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION charecommstruct_set_up
