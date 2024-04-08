#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * This function is responsible for applying boundary conditions (BCs) to both pure outer and inner
 * boundary points. In the first step, it parallelizes the task using OpenMP and starts by applying BCs to
 * the outer boundary points layer-by-layer, prioritizing the faces in the order x0, x1, x2. The second step
 * applies BCs to the inner boundary points, which may map either to the grid interior or to the outer boundary.
 *
 */
void apply_bcs_outerradiation_and_inner(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                        const bc_struct *restrict bcstruct, REAL *restrict xx[3], const REAL custom_wavespeed[NUM_EVOL_GFS],
                                        const REAL custom_f_infinity[NUM_EVOL_GFS], REAL *restrict gfs, REAL *restrict rhs_gfs) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    apply_bcs_outerradiation_and_inner__rfm__SinhSpherical(commondata, params, bcstruct, xx, custom_wavespeed, custom_f_infinity, gfs, rhs_gfs);
    break;
  default:
    fprintf(stderr, "ERROR in apply_bcs_outerradiation_and_inner(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}
