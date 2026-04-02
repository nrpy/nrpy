#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Apply copy boundary conditions to scalar primitives and outflow boundary conditions to the velocity.
 */
void apply_copy_and_outflow_bcs__rfm__Spherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                const ghl_parameters *restrict ghl_params, const bc_struct *restrict bcstruct, REAL *restrict xx[3],
                                                REAL *restrict in_gfs, REAL *restrict auxevol_gfs) {
#include "set_CodeParameters.h"
  // Unpack bc_info from bcstruct.
  const bc_info_struct *bc_info = &bcstruct->bc_info;

  const int NUM_PRIM_GFs = 5;
  const int prims_gfs[5] = {RHOBGF, PGF, RESCALEDVU0GF, RESCALEDVU1GF, RESCALEDVU2GF};

////////////////////////////////////////////////////////
// Step 1 of 2: Apply outer-boundary conditions first.
//              We sweep from the innermost ghost zone
//              outward so edge and corner values are
//              available as soon as they are needed.
#pragma omp parallel
  {
    for (int which_gz = 0; which_gz < NGHOSTS; which_gz++)
      for (int dirn = 0; dirn < 3; dirn++) {
        if (bc_info->num_pure_outer_boundary_points[which_gz][dirn] > 0) {
#pragma omp for
          for (int idx2d = 0; idx2d < bc_info->num_pure_outer_boundary_points[which_gz][dirn]; idx2d++) {
            const short i0 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i0;
            const short i1 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i1;
            const short i2 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i2;
            const short FACEX0 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX0;
            const short FACEX1 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX1;
            const short FACEX2 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX2;
            const int idx = IDX3(i0, i1, i2);
            const int idx_offset1 = IDX3(i0 + 1 * FACEX0, i1 + 1 * FACEX1, i2 + 1 * FACEX2);
            for (int which_gf = 0; which_gf < NUM_PRIM_GFs; which_gf++) {
              auxevol_gfs[IDX4pt(prims_gfs[which_gf], idx)] = auxevol_gfs[IDX4pt(prims_gfs[which_gf], idx_offset1)];
            }

            // Update all velocity components together before recomputing u^t.

            const REAL xx0 = xx[0][i0];
            const REAL xx1 = xx[1][i1];
            const REAL xx2 = xx[2][i2];
            /*
             * NRPy-Generated GF Access/FD Code, Step 1 of 2:
             * Read gridfunction(s) from main memory and compute FD stencils as needed.
             */
            const REAL cf = in_gfs[IDX4(CFGF, i0, i1, i2)];
            const REAL hDD00 = in_gfs[IDX4(HDD00GF, i0, i1, i2)];
            const REAL hDD01 = in_gfs[IDX4(HDD01GF, i0, i1, i2)];
            const REAL hDD02 = in_gfs[IDX4(HDD02GF, i0, i1, i2)];

            /*
             * NRPy-Generated GF Access/FD Code, Step 2 of 2:
             * Evaluate SymPy expressions and write to main memory.
             */
            const REAL FDPart3tmp0 = hDD00 + 1;
            const REAL FDPart3tmp1 = ((cf) * (cf));
            const REAL FDPart3tmp2 = (1.0 / (FDPart3tmp1));
            const REAL FDPart3tmp3 = sqrt(FDPart3tmp0 * FDPart3tmp2);
            const REAL FDPart3tmp5 = FDPart3tmp0 * rescaledvU0 + hDD01 * rescaledvU1 + hDD02 * rescaledvU2;
            const REAL FDPart3tmp6 = FDPart3tmp1 * FDPart3tmp3 * TINYDOUBLE - FDPart3tmp5;
            auxevol_gfs[IDX4(RESCALEDVU0GF, i0, i1, i2)] =
                (FDPart3tmp0 * FDPart3tmp6 * rescaledvU0 -
                 1.0 / 2.0 * FDPart3tmp5 * (FDPart3tmp1 * FDPart3tmp3 * fabs(FDPart3tmp2 * FDPart3tmp5 / FDPart3tmp3) - FDPart3tmp5)) /
                (FDPart3tmp0 * FDPart3tmp6);
            auxevol_gfs[IDX4(RESCALEDVU1GF, i0, i1, i2)] = rescaledvU1;
            auxevol_gfs[IDX4(RESCALEDVU2GF, i0, i1, i2)] = rescaledvU2;

            compute_up_index_velocity_time_component_pointwise(
                commondata, params, &commondata->ghl_params, in_gfs[IDX4pt(ALPHAGF, idx)], in_gfs[IDX4pt(VETU0GF, idx)], in_gfs[IDX4pt(VETU1GF, idx)],
                in_gfs[IDX4pt(VETU2GF, idx)], in_gfs[IDX4pt(HDD00GF, idx)], in_gfs[IDX4pt(HDD01GF, idx)], in_gfs[IDX4pt(HDD02GF, idx)],
                in_gfs[IDX4pt(HDD11GF, idx)], in_gfs[IDX4pt(HDD12GF, idx)], in_gfs[IDX4pt(HDD22GF, idx)], in_gfs[IDX4pt(CFGF, idx)],
                &auxevol_gfs[IDX4pt(RESCALEDVU0GF, idx)], &auxevol_gfs[IDX4pt(RESCALEDVU1GF, idx)], &auxevol_gfs[IDX4pt(RESCALEDVU2GF, idx)],
                &auxevol_gfs[IDX4pt(U4UTGF, idx)]);
          }
        }
      }
  }

  ///////////////////////////////////////////////////////
  // Step 2 of 2: Populate inner-boundary points from
  //              their mapped source points once the
  //              outer boundary data are available.

  // collapse(2) improves throughput here, especially in 2D.
#pragma omp parallel for collapse(2)
  for (int which_gf = 0; which_gf < NUM_PRIM_GFs; which_gf++) {
    for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) {
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt;
      const int srcpt = bcstruct->inner_bc_array[pt].srcpt;
      auxevol_gfs[IDX4pt(prims_gfs[which_gf], dstpt)] =
          bcstruct->inner_bc_array[pt].parity[auxevol_gf_parity[prims_gfs[which_gf]]] * auxevol_gfs[IDX4pt(prims_gfs[which_gf], srcpt)];
    }
  }
} // END FUNCTION apply_copy_and_outflow_bcs__rfm__Spherical
