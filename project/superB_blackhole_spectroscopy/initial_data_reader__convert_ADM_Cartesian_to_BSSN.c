#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Read ADM data in the Cartesian basis, and output rescaled BSSN data in the SinhSpherical basis
 */
void initial_data_reader__convert_ADM_Cartesian_to_BSSN(
    const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], bc_struct *restrict bcstruct,
    MoL_gridfunctions_struct *restrict gridfuncs, ID_persist_struct *restrict ID_persist,
    void ID_function(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
                     const ID_persist_struct *restrict ID_persist, initial_data_struct *restrict initial_data)) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    initial_data_reader__convert_ADM_Cartesian_to_BSSN__rfm__SinhSpherical(commondata, params, xx, bcstruct, gridfuncs, ID_persist, ID_function);
    break;
  default:
    fprintf(stderr, "ERROR in initial_data_reader__convert_ADM_Cartesian_to_BSSN(): CoordSystem hash = %d not #define'd!\n",
            params->CoordSystem_hash);
    exit(1);
  }
}
