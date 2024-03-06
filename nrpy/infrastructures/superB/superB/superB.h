/*

*/

/* superB:  File  "superB.h"*/
#ifndef __SUPERB_H__
#define __SUPERB_H__

#define restrict __restrict__

#include "ckio.h"

#define IDX3_OF_CHARE(i, j, k) ((i) + Nchare0 * ((j) + Nchare1 * ((k))))
#define IDX3GENERAL(i, j, k, Ni, Nj) ((i) + (Ni) * ((j) + (Nj) * (k)))
#define REVERSE_IDX3GENERAL(index, Ni, Nj) \
{ \
	int k = (index) % (Nj); \
	int temp = (index) / (Nj); \
	int j = temp % (Ni); \
	int i = temp / (Ni); \
	i, j, k; \
}

#define OUTPUT_0D 0
#define OUTPUT_1D_Y 1
#define OUTPUT_1D_Z 2
#define OUTPUT_2D_XY 3
#define OUTPUT_2D_YZ 4

#define K_ODD 0
#define K_EVEN 1
#define Y_N 2
#define EAST_WEST 0
#define NORTH_SOUTH 1
#define TOP_BOTTOM 2
#define EAST_GHOST 1
#define WEST_GHOST 2
#define NORTH_GHOST 3
#define SOUTH_GHOST 4
#define TOP_GHOST 5
#define BOTTOM_GHOST 6
#define RK_SUBSTEP_K1 1
#define RK_SUBSTEP_K2 2
#define RK_SUBSTEP_K3 3
#define RK_SUBSTEP_K4 4
#define IDXFACES0(g, inner, j, k) ((j) + Nxx_plus_2NGHOSTS1 * ((k) + Nxx_plus_2NGHOSTS2 * ((inner) + NGHOSTS * (g))))
#define IDXFACES1(g, inner, i, k) ((i) + Nxx_plus_2NGHOSTS0 * ((k) + Nxx_plus_2NGHOSTS2 * ((inner) + NGHOSTS * (g))))
#define IDXFACES2(g, inner, i, j) ((i) + Nxx_plus_2NGHOSTS0 * ((j) + Nxx_plus_2NGHOSTS1 * ((inner) + NGHOSTS * (g))))

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
} diagnostic_struct;

typedef struct __tmpBuffers_struct__ {
  REAL *restrict tmpBuffer_EW;
  REAL *restrict tmpBuffer_NS;
  REAL *restrict tmpBuffer_TB;
} tmpBuffers_struct;

#endif // #ifndef __SUPERB_H__
