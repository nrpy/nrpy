/*
 * This file is part of the nrpy.infrastructures.superB package.
 * It contains definitions and macros for the superB infrastructure.
 *
 * Author: Nishita Jadoo
 *         njadoo **at** uidaho **dot* edu
 */

/* superB:  File  "superB.h"*/
#ifndef __SUPERB_H__
#define __SUPERB_H__

#include "ckio.h"
#include "pup.h"
#ifndef REAL
#define REAL double // for the REALs below
#endif

#define restrict __restrict__
#define MAX3(A, B, C)                                                                                                                                    \
  ({                                                                                                                                                     \
    __typeof__(A) _a = (A);                                                                                                                              \
    __typeof__(B) _b = (B);                                                                                                                              \
    __typeof__(C) _c = (C);                                                                                                                              \
    __typeof__(A) _max = (_a > _b) ? _a : _b;                                                                                                           \
    _max > _c ? _max : _c;                                                                                                                               \
  })
#define IDX2GENERAL(i, j, Ni) ((i) + (Ni) * (j))
#define REVERSE_IDX2GENERAL(index, Ni, i, j) \
{ \
    j = (index) / (Ni); \
    i = (index) % (Ni); \
}
#define IDX3_OF_CHARE(i, j, k) ((i) + Nchare0 * ((j) + Nchare1 * ((k))))
#define IDX3GENERAL(i, j, k, Ni, Nj) ((i) + (Ni) * ((j) + (Nj) * (k)))
#define REVERSE_IDX3GENERAL(index, Ni, Nj, i, j, k) \
{ \
    k = (index) / ((Ni) * (Nj)); \
    int temp = (index) % ((Ni) * (Nj)); \
    j = temp / (Ni); \
    i = temp % (Ni); \
}
#define IDX4GENERAL(g, i, j, k, Ni, Nj, Nk) ((i) + Ni * ((j) + Nj * ((k) + Nk * (g))))
#define MAP_LOCAL_TO_GLOBAL_IDX0(chareidx0, local_idx0, Nxx0chare) ((chareidx0 * Nxx0chare) + local_idx0)
#define MAP_LOCAL_TO_GLOBAL_IDX1(chareidx1, local_idx1, Nxx1chare) ((chareidx1 * Nxx1chare) + local_idx1)
#define MAP_LOCAL_TO_GLOBAL_IDX2(chareidx2, local_idx2, Nxx2chare) ((chareidx2 * Nxx2chare) + local_idx2)

#define MAP_GLOBAL_TO_LOCAL_IDX0(chareidx0, global_idx0, Nxx0chare) (global_idx0 - (chareidx0 * Nxx0chare))// Assumes gridpoint lies within local grid of chare
#define MAP_GLOBAL_TO_LOCAL_IDX1(chareidx1, global_idx1, Nxx1chare) (global_idx1 - (chareidx1 * Nxx1chare))// Assumes gridpoint lies within local grid of chare
#define MAP_GLOBAL_TO_LOCAL_IDX2(chareidx2, global_idx2, Nxx2chare) (global_idx2 - (chareidx2 * Nxx2chare))// Assumes gridpoint point lies within local grid of chare

#define IDXFACES0(g, inner, j, k) ((j) + Nxx_plus_2NGHOSTS1 * ((k) + Nxx_plus_2NGHOSTS2 * ((inner) + NGHOSTS * (g))))
#define IDXFACES1(g, inner, i, k) ((i) + Nxx_plus_2NGHOSTS0 * ((k) + Nxx_plus_2NGHOSTS2 * ((inner) + NGHOSTS * (g))))
#define IDXFACES2(g, inner, i, j) ((i) + Nxx_plus_2NGHOSTS0 * ((j) + Nxx_plus_2NGHOSTS1 * ((inner) + NGHOSTS * (g))))

#define IDX2NONLOCALINNERBC(g, idx, Nidx) ((idx) + Nidx * (g))

#define IDX4PSI4(R, l, m, r_i, Nl, Nm, Nr_i) ((r_i) + (Nr_i) * ((m) + (Nm) * ((l) + (Nl) * (R))))

#define OUTPUT_0D 0
#define OUTPUT_1D_Y 1
#define OUTPUT_1D_Z 2
#define OUTPUT_2D_XY 3
#define OUTPUT_2D_YZ 4
#define OUTPUT_PSI4 5
#define OUTPUT_RESIDUAL 6
#define OUTPUT_L2NORM_BSSN_CONSTRAINTS 7

#define EAST_WEST 0
#define NORTH_SOUTH 1
#define TOP_BOTTOM 2
#define EAST_GHOST 1
#define WEST_GHOST 2
#define NORTH_GHOST 3
#define SOUTH_GHOST 4
#define TOP_GHOST 5
#define BOTTOM_GHOST 6

#define MOL_PRE_RK_UPDATE 0
#define MOL_RK_UPDATE 1
#define MOL_POST_RK_UPDATE_APPLY_BCS 2
#define MOL_POST_RK_UPDATE 3

#define INITIALDATA_BIN_ONE 0
#define INITIALDATA_APPLYBCS_INNERONLY 1
#define INITIALDATA_BIN_TWO 2
#define INITIALDATA_APPLYBCS_OUTEREXTRAPANDINNER 3

typedef struct __charecomm_struct__ {
  int *restrict globalidx3pt_to_chareidx3;    // which chare is evolving or applying bcs to grid point
  int *restrict globalidx3pt_to_localidx3pt;  // local index of grid point on chare that is evolving or setting bcs for gridpoint
  int *restrict localidx3pt_to_globalidx3pt;  // local to this chare
} charecomm_struct;

typedef struct __diagnostic_struct__ {
  int num_output_quantities;
  int tot_num_diagnostic_1d_y_pts;
  int tot_num_diagnostic_1d_z_pts;
  int tot_num_diagnostic_2d_xy_pts;
  int tot_num_diagnostic_2d_yz_pts;
  int num_diagnostic_1d_y_pts;
  int num_diagnostic_1d_z_pts;
  int num_diagnostic_2d_xy_pts;
  int num_diagnostic_2d_yz_pts;
  int *restrict localidx3_diagnostic_1d_y_pt;
  int *restrict locali0_diagnostic_1d_y_pt;
  int *restrict locali1_diagnostic_1d_y_pt;
  int *restrict locali2_diagnostic_1d_y_pt;
  int *restrict offset_diagnostic_1d_y_pt;
  int *restrict localidx3_diagnostic_1d_z_pt;
  int *restrict locali0_diagnostic_1d_z_pt;
  int *restrict locali1_diagnostic_1d_z_pt;
  int *restrict locali2_diagnostic_1d_z_pt;
  int *restrict offset_diagnostic_1d_z_pt;
  int *restrict localidx3_diagnostic_2d_xy_pt;
  int *restrict locali0_diagnostic_2d_xy_pt;
  int *restrict locali1_diagnostic_2d_xy_pt;
  int *restrict locali2_diagnostic_2d_xy_pt;
  int *restrict offset_diagnostic_2d_xy_pt;
  int *restrict localidx3_diagnostic_2d_yz_pt;
  int *restrict locali0_diagnostic_2d_yz_pt;
  int *restrict locali1_diagnostic_2d_yz_pt;
  int *restrict locali2_diagnostic_2d_yz_pt;
  int *restrict offset_diagnostic_2d_yz_pt;
  char filename_1d_y[256];
  char filename_1d_z[256];
  char filename_2d_xy[256];
  char filename_2d_yz[256];
  // psi4:
  int num_of_R_exts_chare;
  int psi4_spinweightm2_sph_harmonics_max_l;
  int length_localsums_for_psi4_decomp;
  REAL *restrict list_of_R_exts_chare;
  REAL *restrict localsums_for_psi4_decomp;
  REAL *restrict globalsums_for_psi4_decomp;
  // psi4 cylindrical-like coords only:
  int tot_N_shell_pts_chare;
  REAL dtheta;
  int *restrict N_shell_pts_chare; // of shape int [num_of_R_exts_chare]
  int *restrict N_theta_shell_chare; // of shape int [num_of_R_exts_chare]
  REAL ***restrict xx_shell_chare; // of shape [num_of_R_exts_chare][N_shell_pts_chare][3]
  REAL **restrict theta_shell_chare; // of shape [num_of_R_exts_chare][N_theta_shell_chare]
} diagnostic_struct;

typedef struct __tmpBuffers_struct__ {
  REAL *restrict tmpBuffer_EW;
  REAL *restrict tmpBuffer_NS;
  REAL *restrict tmpBuffer_TB;
  REAL **restrict tmpBuffer_innerbc_send;
  REAL **restrict tmpBuffer_innerbc_receiv;
} tmpBuffers_struct;

typedef struct __nonlocalinnerbc_struct__ {
  // variables for this chare having the dst pt but not the src pt
  int tot_num_src_chares;
  int *restrict idx3_of_src_chares;
  int *restrict idx3chare_to_src_chare_id;
  int *restrict num_srcpts_each_chare;
  int **restrict map_srcchare_and_srcpt_id_to_linear_id;
  int **restrict globalidx3_srcpts; // of shape [tot_num_src_chares][num_srcpts_each_chare]
  // variables for this chare having the src pt but not the dst pt
  int tot_num_dst_chares;
  int *restrict idx3_of_dst_chares;
  int *restrict idx3chare_to_dst_chare_id;
  int *restrict num_srcpts_tosend_each_chare;
  int **restrict globalidx3_srcpts_tosend; // of shape [tot_num_dst_chares][num_srcpts_tosend_each_chare]
} nonlocalinnerbc_struct;

#endif // #ifndef __SUPERB_H__
