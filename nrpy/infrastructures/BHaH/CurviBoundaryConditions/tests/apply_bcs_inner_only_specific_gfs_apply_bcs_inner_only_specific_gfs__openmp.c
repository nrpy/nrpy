#include "BHaH_defines.h"
#include "intrinsics/simd_intrinsics.h"

static inline REAL apply_parity_branchless(const REAL v, const int8_t p) {
#ifdef __cplusplus
  static_assert(sizeof(REAL) == sizeof(uint32_t) || sizeof(REAL) == sizeof(uint64_t), "REAL must be float or double");
#else
  _Static_assert(sizeof(REAL) == sizeof(uint32_t) || sizeof(REAL) == sizeof(uint64_t), "REAL must be float or double");
#endif

  if (sizeof(REAL) == sizeof(uint64_t)) {
    uint64_t bits;
    REAL out;
    memcpy(&bits, &v, sizeof(bits));
    bits ^= (uint64_t)(p < 0) << 63;
    memcpy(&out, &bits, sizeof(out));
    return out;
  } else {
    uint32_t bits;
    REAL out;
    memcpy(&bits, &v, sizeof(bits));
    bits ^= (uint32_t)(p < 0) << 31;
    memcpy(&out, &bits, sizeof(out));
    return out;
  } // END IF 64 bits vs 32 bits
} // END FUNCTION apply_parity_branchless

/**
 * WARNING: CPU ONLY.
 *
 * Apply BCs to specific grid functions at inner boundary points only,
 * using data stored in bcstruct->inner_bc_array.
 * These structs are set in bcstruct_set_up().
 * Inner boundary points map to either the grid
 * interior ("pure inner") or to pure outer
 * boundary points ("inner maps to outer").
 */
void apply_bcs_inner_only_specific_gfs(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                       const bc_struct *restrict bcstruct, REAL *restrict gfs, const int num_gfs, const int8_t *gf_parities,
                                       const int *gfs_to_sync) {
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;
  const innerpt_bc_struct *restrict inner_bc_array = bcstruct->inner_bc_array;
  const int num_inner_boundary_points = bc_info->num_inner_boundary_points;
  // Needed for IDX macros
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

#pragma omp parallel for schedule(static)
  for (int which_gf = 0; which_gf < num_gfs; ++which_gf) {
    const int gf = gfs_to_sync[which_gf];
    const int parity_idx = gf_parities[gf];
    REAL *restrict gf_data = &gfs[IDX4pt(gf, 0)];

    for (int pt = 0; pt < num_inner_boundary_points; ++pt) {
      const innerpt_bc_struct *restrict bc = &inner_bc_array[pt];
      const REAL v = gf_data[bc->srcpt];
      const int8_t p = bc->parity[parity_idx];
      gf_data[bc->dstpt] = apply_parity_branchless(v, p);
    } // END LOOP over inner boundary points
  } // END LOOP over specific gridfunctions
} // END FUNCTION apply_bcs_inner_only_specific_gfs
