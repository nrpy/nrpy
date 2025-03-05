#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/**
 * Compute 1st derivative finite-difference derivative with arbitrary upwind
 */
__device__ static inline REAL FD1_arbitrary_upwind_x0_dirn(const size_t streamid, const REAL *restrict gf, const int i0, const int i1, const int i2,
                                                           const int offset) {
  MAYBE_UNUSED int const Nxx_plus_2NGHOSTS0 = d_params[streamid].Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED int const Nxx_plus_2NGHOSTS1 = d_params[streamid].Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED int const Nxx_plus_2NGHOSTS2 = d_params[streamid].Nxx_plus_2NGHOSTS2;
  REAL const invdxx0 = d_params[streamid].invdxx0;
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
} // END FUNCTION FD1_arbitrary_upwind_x0_dirn
/**
 * Compute 1st derivative finite-difference derivative with arbitrary upwind
 */
__device__ static inline REAL FD1_arbitrary_upwind_x2_dirn(const size_t streamid, const REAL *restrict gf, const int i0, const int i1, const int i2,
                                                           const int offset) {
  MAYBE_UNUSED int const Nxx_plus_2NGHOSTS0 = d_params[streamid].Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED int const Nxx_plus_2NGHOSTS1 = d_params[streamid].Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED int const Nxx_plus_2NGHOSTS2 = d_params[streamid].Nxx_plus_2NGHOSTS2;
  REAL const invdxx2 = d_params[streamid].invdxx2;
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
} // END FUNCTION FD1_arbitrary_upwind_x2_dirn
/**
 * Compute r(xx0,xx1,xx2) and partial_r x^i.
 */
__device__ static inline void r_and_partial_xi_partial_r_derivs(const size_t streamid, const REAL xx0, const REAL xx1, const REAL xx2, REAL *r,
                                                                REAL *partial_x0_partial_r, REAL *partial_x1_partial_r, REAL *partial_x2_partial_r) {
  const REAL AMPLRHO = d_params[streamid].AMPLRHO;
  const REAL AMPLZ = d_params[streamid].AMPLZ;
  const REAL SINHWRHO = d_params[streamid].SINHWRHO;
  const REAL SINHWZ = d_params[streamid].SINHWZ;

  const REAL tmp0 = (1.0 / (SINHWRHO));
  const REAL tmp6 = (1.0 / (SINHWZ));
  const REAL tmp5 = ((AMPLRHO) * (AMPLRHO)) / ((exp(tmp0) - exp(-tmp0)) * (exp(tmp0) - exp(-tmp0)));
  const REAL tmp7 = exp(tmp6) - exp(-tmp6);
  const REAL tmp2 = exp(tmp0 * xx0);
  const REAL tmp3 = exp(-tmp0 * xx0);
  const REAL tmp9 = exp(tmp6 * xx2);
  const REAL tmp10 = exp(-tmp6 * xx2);
  const REAL tmp24 = (1.0 / (tmp7));
  const REAL tmp4 = tmp2 - tmp3;
  const REAL tmp11 = -tmp10 + tmp9;
  const REAL tmp21 = tmp4 * tmp5 * (2 * tmp0 * tmp2 + 2 * tmp0 * tmp3);
  const REAL tmp13 = ((AMPLZ) * (AMPLZ)) * ((tmp11) * (tmp11)) / ((tmp7) * (tmp7));
  const REAL tmp18 = ((AMPLZ) * (AMPLZ) * (AMPLZ)) * ((tmp11) * (tmp11)) * (2 * tmp10 * tmp6 + 2 * tmp6 * tmp9) / ((tmp7) * (tmp7) * (tmp7));
  const REAL tmp14 = tmp13 + ((tmp4) * (tmp4)) * tmp5;
  const REAL tmp15 = sqrt(tmp14);
  const REAL tmp19 = (1.0 / sqrt(-tmp13 / tmp14 + 1));
  const REAL tmp25 = (1.0 / 2.0) / pow(tmp14, 3.0 / 2.0);
  const REAL tmp23 = (1.0 / (tmp15));
  const REAL tmp26 = tmp19 * (AMPLZ * tmp23 * tmp24 * (tmp10 * tmp6 + tmp6 * tmp9) - tmp18 * tmp25);
  const REAL tmp27 = (1.0 / ((1.0 / 2.0) * tmp21 * tmp23 * tmp26 + (1.0 / 4.0) * tmp18 * tmp19 * tmp21 / ((tmp14) * (tmp14))));
  *r = tmp15;
  *partial_x0_partial_r = tmp26 * tmp27;
  *partial_x1_partial_r = 0;
  *partial_x2_partial_r = AMPLZ * tmp11 * tmp19 * tmp21 * tmp24 * tmp25 * tmp27;
} // END FUNCTION r_and_partial_xi_partial_r_derivs
/**
 * Compute \partial_r f
 */
__device__ static inline REAL compute_partial_r_f(const size_t streamid, REAL *restrict xx[3], const REAL *restrict gfs, const int which_gf,
                                                  const int dest_i0, const int dest_i1, const int dest_i2, const int FACEi0, const int FACEi1,
                                                  const int FACEi2, const REAL partial_x0_partial_r, const REAL partial_x1_partial_r,
                                                  const REAL partial_x2_partial_r) {
  ///////////////////////////////////////////////////////////

  // FD1_stencil_radius = radiation_BC_fd_order/2 = 1
  const int FD1_stencil_radius = 1;
  int const Nxx_plus_2NGHOSTS0 = d_params[streamid].Nxx_plus_2NGHOSTS0;
  int const Nxx_plus_2NGHOSTS1 = d_params[streamid].Nxx_plus_2NGHOSTS1;
  int const Nxx_plus_2NGHOSTS2 = d_params[streamid].Nxx_plus_2NGHOSTS2;
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
  const REAL partial_x0_f = FD1_arbitrary_upwind_x0_dirn(streamid, &gfs[which_gf * ntot], dest_i0, dest_i1, dest_i2, i0_offset);
  const REAL partial_x1_f = 0.0;
  int i2_offset = FACEi2; // Shift stencil away from the face we're updating.
  // Next adjust i2_offset so that FD stencil never goes out of bounds.
  if (dest_i2 < FD1_stencil_radius)
    i2_offset = FD1_stencil_radius - dest_i2;
  else if (dest_i2 > (Nxx_plus_2NGHOSTS2 - FD1_stencil_radius - 1))
    i2_offset = (Nxx_plus_2NGHOSTS2 - FD1_stencil_radius - 1) - dest_i2;
  const REAL partial_x2_f = FD1_arbitrary_upwind_x2_dirn(streamid, &gfs[which_gf * ntot], dest_i0, dest_i1, dest_i2, i2_offset);
  return partial_x0_partial_r * partial_x0_f + partial_x1_partial_r * partial_x1_f + partial_x2_partial_r * partial_x2_f;
} // END FUNCTION compute_partial_r_f

/**
 * *** Apply radiation BCs to all outer boundaries. ***
 *
 */
__device__ static inline REAL radiation_bcs(const size_t streamid, REAL *restrict xx[3], const REAL *restrict gfs, REAL *restrict gfs_rhss,
                                            const int which_gf, const REAL gf_wavespeed, const REAL gf_f_infinity, const int dest_i0,
                                            const int dest_i1, const int dest_i2, const short FACEi0, const short FACEi1, const short FACEi2) {
  int const Nxx_plus_2NGHOSTS0 = d_params[streamid].Nxx_plus_2NGHOSTS0;
  int const Nxx_plus_2NGHOSTS1 = d_params[streamid].Nxx_plus_2NGHOSTS1;
  int const Nxx_plus_2NGHOSTS2 = d_params[streamid].Nxx_plus_2NGHOSTS2;
  // Nearest "interior" neighbor of this gridpoint, based on current face
  const int dest_i0_int = dest_i0 + 1 * FACEi0, dest_i1_int = dest_i1 + 1 * FACEi1, dest_i2_int = dest_i2 + 1 * FACEi2;
  REAL r, partial_x0_partial_r, partial_x1_partial_r, partial_x2_partial_r;
  REAL r_int, partial_x0_partial_r_int, partial_x1_partial_r_int, partial_x2_partial_r_int;
  r_and_partial_xi_partial_r_derivs(streamid, xx[0][dest_i0], xx[1][dest_i1], xx[2][dest_i2], &r, &partial_x0_partial_r, &partial_x1_partial_r,
                                    &partial_x2_partial_r);
  r_and_partial_xi_partial_r_derivs(streamid, xx[0][dest_i0_int], xx[1][dest_i1_int], xx[2][dest_i2_int], &r_int, &partial_x0_partial_r_int,
                                    &partial_x1_partial_r_int, &partial_x2_partial_r_int);
  const REAL partial_r_f = compute_partial_r_f(streamid, xx, gfs, which_gf, dest_i0, dest_i1, dest_i2, FACEi0, FACEi1, FACEi2, partial_x0_partial_r,
                                               partial_x1_partial_r, partial_x2_partial_r);
  const REAL partial_r_f_int = compute_partial_r_f(streamid, xx, gfs, which_gf, dest_i0_int, dest_i1_int, dest_i2_int, FACEi0, FACEi1, FACEi2,
                                                   partial_x0_partial_r_int, partial_x1_partial_r_int, partial_x2_partial_r_int);

  const int idx3 = IDX3(dest_i0, dest_i1, dest_i2);
  const int idx3_int = IDX3(dest_i0_int, dest_i1_int, dest_i2_int);

  const REAL partial_t_f_int = gfs_rhss[IDX4pt(which_gf, idx3_int)];

  const REAL c = gf_wavespeed;
  const REAL f_infinity = gf_f_infinity;
  const REAL f = gfs[IDX4pt(which_gf, idx3)];
  const REAL f_int = gfs[IDX4pt(which_gf, idx3_int)];
  const REAL partial_t_f_int_outgoing_wave = -c * (partial_r_f_int + (f_int - f_infinity) / r_int);

  const REAL k = r_int * r_int * r_int * (partial_t_f_int - partial_t_f_int_outgoing_wave);

  const REAL rinv = 1.0 / r;
  const REAL partial_t_f_outgoing_wave = -c * (partial_r_f + (f - f_infinity) * rinv);

  return partial_t_f_outgoing_wave + k * rinv * rinv * rinv;
} // END FUNCTION radiation_bcs
/**
 * Kernel: apply_bcs_pure_only_gpu.
 * Apply BCs to pure points.
 */
__global__ static void apply_bcs_pure_only_gpu(const size_t streamid, const int num_pure_outer_boundary_points, const int which_gz, const int dirn,
                                               const outerpt_bc_struct *restrict pure_outer_bc_array, REAL *restrict gfs, REAL *restrict rhs_gfs,
                                               REAL *restrict x0, REAL *restrict x1, REAL *restrict x2) {
  int const Nxx_plus_2NGHOSTS0 = d_params[streamid].Nxx_plus_2NGHOSTS0;
  int const Nxx_plus_2NGHOSTS1 = d_params[streamid].Nxx_plus_2NGHOSTS1;
  int const Nxx_plus_2NGHOSTS2 = d_params[streamid].Nxx_plus_2NGHOSTS2;

  // Thread indices
  // Global data index - expecting a 1D dataset
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;

  // Thread strides
  const int stride0 = blockDim.x * gridDim.x;
  for (int idx2d = tid0; idx2d < num_pure_outer_boundary_points; idx2d += stride0) {
    const short i0 = pure_outer_bc_array[idx2d].i0;
    const short i1 = pure_outer_bc_array[idx2d].i1;
    const short i2 = pure_outer_bc_array[idx2d].i2;
    const short FACEX0 = pure_outer_bc_array[idx2d].FACEX0;
    const short FACEX1 = pure_outer_bc_array[idx2d].FACEX1;
    const short FACEX2 = pure_outer_bc_array[idx2d].FACEX2;
    const int idx3 = IDX3(i0, i1, i2);
    REAL *xx[3] = {x0, x1, x2};
    for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
      // *** Apply radiation BCs to all outer boundary points. ***
      rhs_gfs[IDX4pt(which_gf, idx3)] = radiation_bcs(streamid, xx, gfs, rhs_gfs, which_gf, d_gridfunctions_wavespeed[which_gf],
                                                      d_gridfunctions_f_infinity[which_gf], i0, i1, i2, FACEX0, FACEX1, FACEX2);
    }
  }
} // END FUNCTION apply_bcs_pure_only_gpu
/**
 * Kernel: apply_bcs_pure_only.
 * Apply BCs to pure boundary points
 */
static void apply_bcs_pure_only(const params_struct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict *xx, REAL *restrict gfs,
                                REAL *restrict rhs_gfs, const REAL *custom_wavespeed, const REAL *custom_f_infinity) {

  const bc_info_struct *bc_info = &bcstruct->bc_info;
  REAL *restrict x0 = xx[0];
  REAL *restrict x1 = xx[1];
  REAL *restrict x2 = xx[2];
  for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
    for (int dirn = 0; dirn < 3; dirn++) {
      if (bc_info->num_pure_outer_boundary_points[which_gz][dirn] > 0) {
        size_t gz_idx = dirn + (3 * which_gz);
        const outerpt_bc_struct *restrict pure_outer_bc_array = bcstruct->pure_outer_bc_array[gz_idx];
        int num_pure_outer_boundary_points = bc_info->num_pure_outer_boundary_points[which_gz][dirn];
        {

          const size_t threads_in_x_dir = 32;
          const size_t threads_in_y_dir = 1;
          const size_t threads_in_z_dir = 1;
          dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
          dim3 blocks_per_grid((num_pure_outer_boundary_points + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
          size_t sm = 0;
          size_t streamid = params->grid_idx % NUM_STREAMS;
          apply_bcs_pure_only_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, num_pure_outer_boundary_points, which_gz,
                                                                                                 dirn, pure_outer_bc_array, gfs, rhs_gfs, x0, x1, x2);
          cudaCheckErrors(cudaKernel, "apply_bcs_pure_only_gpu failure");
        }
      }
    }
  }
} // END FUNCTION apply_bcs_pure_only

/**
 * This function is responsible for applying boundary conditions (BCs) to both pure outer and inner
 * boundary points. In the first step, it parallelizes the task using OpenMP and starts by applying BCs to
 * the outer boundary points layer-by-layer, prioritizing the faces in the order x0, x1, x2. The second step
 * applies BCs to the inner boundary points, which may map either to the grid interior or to the outer boundary.
 *
 */
void apply_bcs_outerradiation_and_inner__rfm__SinhCylindrical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                              const bc_struct *restrict bcstruct, REAL *restrict xx[3],
                                                              const REAL custom_wavespeed[NUM_EVOL_GFS], const REAL custom_f_infinity[NUM_EVOL_GFS],
                                                              REAL *restrict gfs, REAL *restrict rhs_gfs) {
#include "../set_CodeParameters.h"

  // Update device constants
  cudaMemcpyToSymbol(d_gridfunctions_wavespeed, custom_wavespeed, NUM_EVOL_GFS * sizeof(REAL));
  cudaCheckErrors(copy, "Copy to d_gridfunctions_wavespeed failed");
  cudaMemcpyToSymbol(d_gridfunctions_f_infinity, custom_f_infinity, NUM_EVOL_GFS * sizeof(REAL));
  cudaCheckErrors(copy, "Copy to d_gridfunctions_f_infinity failed");

  ////////////////////////////////////////////////////////
  // STEP 1 of 2: Apply BCs to pure outer boundary points.
  //              By "pure" we mean that these points are
  //              on the outer boundary and not also on
  //              an inner boundary.
  //              Here we fill in the innermost ghost zone
  //              layer first and move outward. At each
  //              layer, we fill in +/- x0 faces first,
  //              then +/- x1 faces, finally +/- x2 faces,
  //              filling in the edges as we go.
  apply_bcs_pure_only(params, bcstruct, xx, gfs, rhs_gfs, custom_wavespeed, custom_f_infinity);

  ///////////////////////////////////////////////////////
  // STEP 2 of 2: Apply BCs to inner boundary points.
  //              These map to either the grid interior
  //              ("pure inner") or to pure outer boundary
  //              points ("inner maps to outer"). Those
  //              that map to outer require that outer be
  //              populated first; hence this being
  //              STEP 2 OF 2.
  apply_bcs_inner_only(commondata, params, bcstruct, rhs_gfs); // <- apply inner BCs to RHS gfs only
} // END FUNCTION apply_bcs_outerradiation_and_inner__rfm__SinhCylindrical
