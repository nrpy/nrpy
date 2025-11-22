#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Set RHSs for hyperbolic relaxation equation.
 */
void rhs_eval(const commondata_struct *restrict commondata, const params_struct *restrict params, const rfm_struct *restrict rfmstruct,
              const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict rhs_gfs) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    rhs_eval__rfm__SinhSpherical(commondata, params, rfmstruct, auxevol_gfs, in_gfs, rhs_gfs);
    break;
  default:
    fprintf(stderr, "ERROR in rhs_eval(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION rhs_eval
