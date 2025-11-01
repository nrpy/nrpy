#include "BHaH_defines.h"

/**
 * Kernel: apply_bcs_inner_only_specific_gfs_host.
 * Apply BCs to inner boundary points only for specified GFs.
 */
static void apply_bcs_inner_only_specific_gfs_host(const params_struct *restrict params, const int num_inner_boundary_points,
                                                   const innerpt_bc_struct *restrict inner_bc_array, REAL *restrict gfs, const int num_gfs,
                                                   const int *restrict gfs_to_sync) {
  // Needed for IDX macros
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;

  // collapse(2) results in a nice speedup here, esp in 2D. Two_BHs_collide goes from
  //    5550 M/hr to 7264 M/hr on a Ryzen 9 5950X running on all 16 cores with core affinity.
#pragma omp parallel for collapse(2) // spawn threads and distribute across them
  for (int which_gf = 0; which_gf < num_gfs; which_gf++) {
    for (int pt = 0; pt < num_inner_boundary_points; pt++) {
      const int dstpt = inner_bc_array[pt].dstpt;
      const int srcpt = inner_bc_array[pt].srcpt;
      gfs[IDX4pt(gfs_to_sync[which_gf], dstpt)] =
          inner_bc_array[pt].parity[aux_gf_parity[gfs_to_sync[which_gf]]] * gfs[IDX4pt(gfs_to_sync[which_gf], srcpt)];
    } // END for(int pt=0;pt<num_inner_pts;pt++)
  } // END for(int which_gf=0;which_gf<num_gfs;which_gf++)
} // END FUNCTION apply_bcs_inner_only_specific_gfs_host

/**
 * Apply BCs to specific grid functions at inner boundary points only,
 * using data stored in bcstruct->inner_bc_array.
 * These structs are set in bcstruct_set_up().
 * Inner boundary points map to either the grid
 * interior ("pure inner") or to pure outer
 * boundary points ("inner maps to outer").
 */
void apply_bcs_inner_only_specific_gfs(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                       const bc_struct *restrict bcstruct, REAL *restrict gfs, const int num_gfs, const int *gfs_to_sync) {
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;
  const innerpt_bc_struct *restrict inner_bc_array = bcstruct->inner_bc_array;
  const int num_inner_boundary_points = bc_info->num_inner_boundary_points;
  apply_bcs_inner_only_specific_gfs_host(params, num_inner_boundary_points, inner_bc_array, gfs, num_gfs, gfs_to_sync);
} // END FUNCTION apply_bcs_inner_only_specific_gfs
