#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Compute 1st derivative finite-difference derivative with arbitrary upwind
 */
static inline REAL FD1_arbitrary_upwind_x0_dirn(const params_struct *restrict params, const REAL *restrict gf, const int i0, const int i1,
                                                const int i2, const int offset) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;

  switch (offset) {
  case 0: {
    static const REAL Rational_decl__1_2 = 1.0 / 2.0;
    return (-Rational_decl__1_2 * gf[IDX3(i0 - 1, i1, i2)] + Rational_decl__1_2 * gf[IDX3(i0 + 1, i1, i2)]) * invdxx0;
  }
  case 1: {
    static const REAL Rational_decl__3_2 = 3.0 / 2.0;
    static const REAL Rational_decl__2_1 = 2;
    static const REAL Rational_decl__1_2 = 1.0 / 2.0;
    return (-Rational_decl__3_2 * gf[IDX3(i0, i1, i2)] + Rational_decl__2_1 * gf[IDX3(i0 + 1, i1, i2)] -
            Rational_decl__1_2 * gf[IDX3(i0 + 2, i1, i2)]) *
           invdxx0;
  }
  case -1: {
    static const REAL Rational_decl__1_2 = 1.0 / 2.0;
    static const REAL Rational_decl__2_1 = 2;
    static const REAL Rational_decl__3_2 = 3.0 / 2.0;
    return (+Rational_decl__1_2 * gf[IDX3(i0 - 2, i1, i2)] - Rational_decl__2_1 * gf[IDX3(i0 - 1, i1, i2)] +
            Rational_decl__3_2 * gf[IDX3(i0, i1, i2)]) *
           invdxx0;
  }
  }
  return 0.0 / 0.0; // poison output if offset computed incorrectly
} // END FUNCTION: FD1_arbitrary_upwind_x0_dirn
/**
 * Compute 1st derivative finite-difference derivative with arbitrary upwind
 */
static inline REAL FD1_arbitrary_upwind_x1_dirn(const params_struct *restrict params, const REAL *restrict gf, const int i0, const int i1,
                                                const int i2, const int offset) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;

  switch (offset) {
  case 0: {
    static const REAL Rational_decl__1_2 = 1.0 / 2.0;
    return (-Rational_decl__1_2 * gf[IDX3(i0, i1 - 1, i2)] + Rational_decl__1_2 * gf[IDX3(i0, i1 + 1, i2)]) * invdxx1;
  }
  case 1: {
    static const REAL Rational_decl__3_2 = 3.0 / 2.0;
    static const REAL Rational_decl__2_1 = 2;
    static const REAL Rational_decl__1_2 = 1.0 / 2.0;
    return (-Rational_decl__3_2 * gf[IDX3(i0, i1, i2)] + Rational_decl__2_1 * gf[IDX3(i0, i1 + 1, i2)] -
            Rational_decl__1_2 * gf[IDX3(i0, i1 + 2, i2)]) *
           invdxx1;
  }
  case -1: {
    static const REAL Rational_decl__1_2 = 1.0 / 2.0;
    static const REAL Rational_decl__2_1 = 2;
    static const REAL Rational_decl__3_2 = 3.0 / 2.0;
    return (+Rational_decl__1_2 * gf[IDX3(i0, i1 - 2, i2)] - Rational_decl__2_1 * gf[IDX3(i0, i1 - 1, i2)] +
            Rational_decl__3_2 * gf[IDX3(i0, i1, i2)]) *
           invdxx1;
  }
  }
  return 0.0 / 0.0; // poison output if offset computed incorrectly
} // END FUNCTION: FD1_arbitrary_upwind_x1_dirn
/**
 * Compute 1st derivative finite-difference derivative with arbitrary upwind
 */
static inline REAL FD1_arbitrary_upwind_x2_dirn(const params_struct *restrict params, const REAL *restrict gf, const int i0, const int i1,
                                                const int i2, const int offset) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;

  switch (offset) {
  case 0: {
    static const REAL Rational_decl__1_2 = 1.0 / 2.0;
    return (-Rational_decl__1_2 * gf[IDX3(i0, i1, i2 - 1)] + Rational_decl__1_2 * gf[IDX3(i0, i1, i2 + 1)]) * invdxx2;
  }
  case 1: {
    static const REAL Rational_decl__3_2 = 3.0 / 2.0;
    static const REAL Rational_decl__2_1 = 2;
    static const REAL Rational_decl__1_2 = 1.0 / 2.0;
    return (-Rational_decl__3_2 * gf[IDX3(i0, i1, i2)] + Rational_decl__2_1 * gf[IDX3(i0, i1, i2 + 1)] -
            Rational_decl__1_2 * gf[IDX3(i0, i1, i2 + 2)]) *
           invdxx2;
  }
  case -1: {
    static const REAL Rational_decl__1_2 = 1.0 / 2.0;
    static const REAL Rational_decl__2_1 = 2;
    static const REAL Rational_decl__3_2 = 3.0 / 2.0;
    return (+Rational_decl__1_2 * gf[IDX3(i0, i1, i2 - 2)] - Rational_decl__2_1 * gf[IDX3(i0, i1, i2 - 1)] +
            Rational_decl__3_2 * gf[IDX3(i0, i1, i2)]) *
           invdxx2;
  }
  }
  return 0.0 / 0.0; // poison output if offset computed incorrectly
} // END FUNCTION: FD1_arbitrary_upwind_x2_dirn
/**
 * Compute \partial_r f
 */
static inline REAL compute_partial_r_f(const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict gfs, const int which_gf,
                                       const int dest_i0, const int dest_i1, const int dest_i2, const int FACEi0, const int FACEi1, const int FACEi2,
                                       const REAL partial_x0_partial_r, const REAL partial_x1_partial_r, const REAL partial_x2_partial_r) {
  ///////////////////////////////////////////////////////////

  // FD1_stencil_radius = radiation_BC_fd_order/2 = 1
  const int FD1_stencil_radius = 1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;

  const int ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

  ///////////////////////////////////////////////////////////
  // Next we'll compute partial_xi f, using a maximally-centered stencil.
  //   The {{i0,i1,i2}}_offset parameters set the offset of the maximally-centered
  //   stencil, such that an offset=0 implies a centered stencil.

  // CHECK: Nxx_plus_2NGHOSTS0=10; FD1_stencil_radius=2. Then Nxx_plus_2NGHOSTS0-FD1_stencil_radius-1 = 7
  //  if dest_i0 = 9, we get i0_offset=7-9=-2, so the (4th order) deriv
  //  stencil is: -4,-3,-2,-1,0

  // CHECK: if FD1_stencil_radius=2 and dest_i0 = 1, we get i0_offset = FD1_stencil_radius-dest_i0 = 1,
  //  so the (4th order) deriv stencil is: -1,0,1,2,3

  // CHECK: if FD1_stencil_radius=2 and dest_i0 = 0, we get i0_offset = FD1_stencil_radius-1 = 2,
  //  so the (4th order) deriv stencil is: 0,1,2,3,4
  int i0_offset = FACEi0; // Shift stencil away from the face we're updating.
  // Next adjust i0_offset so that FD stencil never goes out of bounds.
  if (dest_i0 < FD1_stencil_radius)
    i0_offset = FD1_stencil_radius - dest_i0;
  else if (dest_i0 > (Nxx_plus_2NGHOSTS0 - FD1_stencil_radius - 1))
    i0_offset = (Nxx_plus_2NGHOSTS0 - FD1_stencil_radius - 1) - dest_i0;
  const REAL partial_x0_f = FD1_arbitrary_upwind_x0_dirn(params, &gfs[which_gf * ntot], dest_i0, dest_i1, dest_i2, i0_offset);
  int i1_offset = FACEi1; // Shift stencil away from the face we're updating.
  // Next adjust i1_offset so that FD stencil never goes out of bounds.
  if (dest_i1 < FD1_stencil_radius)
    i1_offset = FD1_stencil_radius - dest_i1;
  else if (dest_i1 > (Nxx_plus_2NGHOSTS1 - FD1_stencil_radius - 1))
    i1_offset = (Nxx_plus_2NGHOSTS1 - FD1_stencil_radius - 1) - dest_i1;
  const REAL partial_x1_f = FD1_arbitrary_upwind_x1_dirn(params, &gfs[which_gf * ntot], dest_i0, dest_i1, dest_i2, i1_offset);
  int i2_offset = FACEi2; // Shift stencil away from the face we're updating.
  // Next adjust i2_offset so that FD stencil never goes out of bounds.
  if (dest_i2 < FD1_stencil_radius)
    i2_offset = FD1_stencil_radius - dest_i2;
  else if (dest_i2 > (Nxx_plus_2NGHOSTS2 - FD1_stencil_radius - 1))
    i2_offset = (Nxx_plus_2NGHOSTS2 - FD1_stencil_radius - 1) - dest_i2;
  const REAL partial_x2_f = FD1_arbitrary_upwind_x2_dirn(params, &gfs[which_gf * ntot], dest_i0, dest_i1, dest_i2, i2_offset);
  return partial_x0_partial_r * partial_x0_f + partial_x1_partial_r * partial_x1_f + partial_x2_partial_r * partial_x2_f;
} // END FUNCTION: compute_partial_r_f

/**
 * *** Apply radiation BCs to all outer boundaries. ***
 */
static inline REAL radiation_bcs(const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict gfs, const REAL partial_t_f_int,
                                 const int which_gf, const REAL gf_wavespeed, const REAL gf_f_infinity, const int dest_i0, const int dest_i1,
                                 const int dest_i2, const short FACEi0, const short FACEi1, const short FACEi2, const REAL r,
                                 const REAL partial_x0_partial_r, const REAL partial_x1_partial_r, const REAL partial_x2_partial_r, const REAL r_int,
                                 const REAL partial_x0_partial_r_int, const REAL partial_x1_partial_r_int, const REAL partial_x2_partial_r_int) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;

  // Nearest "interior" neighbor of this gridpoint, based on current face
  const int dest_i0_int = dest_i0 + 1 * FACEi0, dest_i1_int = dest_i1 + 1 * FACEi1, dest_i2_int = dest_i2 + 1 * FACEi2;
  const REAL partial_r_f = compute_partial_r_f(params, xx, gfs, which_gf, dest_i0, dest_i1, dest_i2, FACEi0, FACEi1, FACEi2, partial_x0_partial_r,
                                               partial_x1_partial_r, partial_x2_partial_r);
  const REAL partial_r_f_int = compute_partial_r_f(params, xx, gfs, which_gf, dest_i0_int, dest_i1_int, dest_i2_int, FACEi0, FACEi1, FACEi2,
                                                   partial_x0_partial_r_int, partial_x1_partial_r_int, partial_x2_partial_r_int);

  const int idx3 = IDX3(dest_i0, dest_i1, dest_i2);
  const int idx3_int = IDX3(dest_i0_int, dest_i1_int, dest_i2_int);

  const REAL c = gf_wavespeed;
  const REAL f_infinity = gf_f_infinity;
  const REAL f = gfs[IDX4pt(which_gf, idx3)];
  const REAL f_int = gfs[IDX4pt(which_gf, idx3_int)];
  const REAL partial_t_f_int_outgoing_wave = -c * (partial_r_f_int + (f_int - f_infinity) / r_int);

  const REAL k = r_int * r_int * r_int * (partial_t_f_int - partial_t_f_int_outgoing_wave);

  const REAL rinv = 1.0 / r;
  const REAL partial_t_f_outgoing_wave = -c * (partial_r_f + (f - f_infinity) * rinv);

  return partial_t_f_outgoing_wave + k * rinv * rinv * rinv;
} // END FUNCTION: radiation_bcs
/**
 * Kernel: apply_bcs_inner_accumulate_host.
 * Prepare or complete inner-boundary RHS accumulation in two RK arrays.
 */
static void apply_bcs_inner_accumulate_host(const params_struct *restrict params, const int num_inner_boundary_points,
                                            const innerpt_bc_struct *restrict inner_bc_array, REAL *restrict stage_gfs, REAL *restrict accum_gfs,
                                            const REAL stage_dt, const REAL accum_dt, const int prepare) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;

#pragma omp parallel for
  for (int pt = 0; pt < num_inner_boundary_points; pt++) {

    const innerpt_bc_struct *restrict bc = &inner_bc_array[pt];
    const REAL sign = prepare ? -1.0 : 1.0;
    for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
      const int dst = IDX4pt(which_gf, bc->dstpt);
      const int src = IDX4pt(which_gf, bc->srcpt);
      const REAL weight = sign * bc->parity[evol_gf_parity[which_gf]];
      if (stage_dt != 0.0)
        stage_gfs[dst] += weight * stage_gfs[src];
      if (accum_dt != 0.0)
        accum_gfs[dst] += weight * accum_gfs[src];
    }
  }
} // END FUNCTION: apply_bcs_inner_accumulate_host
/**
 * Kernel: apply_bcs_pure_accumulate_host.
 * Prepare or complete radiation-boundary RHS accumulation in two RK arrays.
 */
static void apply_bcs_pure_accumulate_host(const params_struct *restrict params, const int num_pure_outer_boundary_points,
                                           const outerpt_bc_struct *restrict pure_outer_bc_array, const REAL *restrict in_gfs,
                                           REAL *restrict stage_gfs, REAL *restrict accum_gfs, REAL *restrict x0, REAL *restrict x1,
                                           REAL *restrict x2, const REAL stage_dt, const REAL accum_dt, const int prepare,
                                           const REAL *custom_wavespeed, const REAL *custom_f_infinity) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;

#pragma omp parallel for
  for (int pt = 0; pt < num_pure_outer_boundary_points; pt++) {

    const outerpt_bc_struct *restrict bc = &pure_outer_bc_array[pt];
    const int dstpt = IDX3(bc->i0, bc->i1, bc->i2);
    const int srcpt = IDX3(bc->i0 + bc->FACEX0, bc->i1 + bc->FACEX1, bc->i2 + bc->FACEX2);
    const REAL rinv = 1.0 / bc->r;
    const REAL weight = bc->r_int * bc->r_int * bc->r_int * rinv * rinv * rinv;
    const REAL sign = prepare ? -1.0 : 1.0;
    REAL *xx[3] = {x0, x1, x2};
    for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
      const int dst = IDX4pt(which_gf, dstpt);
      const int src = IDX4pt(which_gf, srcpt);
      if (stage_dt != 0.0)
        stage_gfs[dst] += sign * weight * stage_gfs[src];
      if (accum_dt != 0.0)
        accum_gfs[dst] += sign * weight * accum_gfs[src];
      if (!prepare) {
        const REAL beta =
            radiation_bcs(params, xx, in_gfs, 0.0, which_gf, custom_wavespeed[which_gf], custom_f_infinity[which_gf], bc->i0, bc->i1, bc->i2,
                          bc->FACEX0, bc->FACEX1, bc->FACEX2, bc->r, bc->partial_x0_partial_r, bc->partial_x1_partial_r, bc->partial_x2_partial_r,
                          bc->r_int, bc->partial_x0_partial_r_int, bc->partial_x1_partial_r_int, bc->partial_x2_partial_r_int);
        if (stage_dt != 0.0)
          stage_gfs[dst] += stage_dt * beta;
        if (accum_dt != 0.0)
          accum_gfs[dst] += accum_dt * beta;
      }
    }
  }
} // END FUNCTION: apply_bcs_pure_accumulate_host

/**
 * Apply radiation and inner RHS boundary conditions while accumulating into two RK arrays.
 */
void apply_bcs_outerradiation_and_inner_accumulate__rfm__Cartesian(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                                   const bc_struct *restrict bcstruct, REAL *restrict xx[3],
                                                                   const REAL custom_wavespeed[NUM_EVOL_GFS],
                                                                   const REAL custom_f_infinity[NUM_EVOL_GFS], const REAL *restrict in_gfs,
                                                                   REAL *restrict stage_gfs, REAL *restrict accum_gfs, const REAL stage_dt,
                                                                   const REAL accum_dt, const int prepare) {
  const bc_info_struct *bc_info = &bcstruct->bc_info;
  const int num_inner_boundary_points = bc_info->num_inner_boundary_points;
  const innerpt_bc_struct *restrict inner_bc_array = bcstruct->inner_bc_array;
  REAL *restrict x0 = xx[0];
  REAL *restrict x1 = xx[1];
  REAL *restrict x2 = xx[2];

  // Inner sources may be pure outer points: subtract them before changing outer values.
  if (prepare) {
    apply_bcs_inner_accumulate_host(params, num_inner_boundary_points, inner_bc_array, stage_gfs, accum_gfs, stage_dt, accum_dt, prepare);
  }

  // Reverse the complete layer/face dependency order only during preparation.
  for (int layer = 0; layer < NGHOSTS; layer++) {
    const int which_gz = prepare ? NGHOSTS - 1 - layer : layer;
    for (int face = 0; face < 3; face++) {
      const int dirn = prepare ? 2 - face : face;
      const int num_pure_outer_boundary_points = bc_info->num_pure_outer_boundary_points[which_gz][dirn];
      if (num_pure_outer_boundary_points > 0) {
        const int gz_idx = dirn + 3 * which_gz;
        const outerpt_bc_struct *restrict pure_outer_bc_array = bcstruct->pure_outer_bc_array[gz_idx];
        apply_bcs_pure_accumulate_host(params, num_pure_outer_boundary_points, pure_outer_bc_array, in_gfs, stage_gfs, accum_gfs, x0, x1, x2,
                                       stage_dt, accum_dt, prepare, custom_wavespeed, custom_f_infinity);
      }
    }
  }

  // All outer sources must have their final accumulated values before inner completion.
  if (!prepare) {
    apply_bcs_inner_accumulate_host(params, num_inner_boundary_points, inner_bc_array, stage_gfs, accum_gfs, stage_dt, accum_dt, prepare);
  }
} // END FUNCTION: apply_bcs_outerradiation_and_inner_accumulate__rfm__Cartesian
