#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * At each coordinate point (x0,x1,x2) situated at grid index (i0,i1,i2):
 * Step 1: Set up inner boundary structs bcstruct->inner_bc_array[].
 * Recall that at each inner boundary point we must set innerpt_bc_struct:
 * typedef struct __innerpt_bc_struct__ {
 * int dstpt;  // dstpt is the 3D grid index IDX3(i0,i1,i2) of the inner boundary point (i0,i1,i2)
 * int srcpt;  // srcpt is the 3D grid index (a la IDX3) to which the inner boundary point maps
 * int8_t parity[10];  // parity[10] is a calculation of dot products for the 10 independent parity types
 * } innerpt_bc_struct;
 * At each ghostzone (i.e., each point within NGHOSTS points from grid boundary):
 * Call EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt().
 * This function converts the curvilinear coordinate (x0,x1,x2) to the corresponding
 * Cartesian coordinate (x,y,z), then finds the grid point
 * (i0_inbounds,i1_inbounds,i2_inbounds) in the grid interior or outer boundary
 * corresponding to this Cartesian coordinate (x,y,z).
 * If (i0,i1,i2) *is not* the same as (i0_inbounds,i1_inbounds,i2_inbounds),
 * then we are at an inner boundary point. We must set
 * Set bcstruct->inner_bc_array for this point, which requires we specify
 * both (i0_inbounds,i1_inbounds,i2_inbounds) [just found!] and parity
 * conditions for this gridpoint. The latter is found & specified within the
 * function set_parity_for_inner_boundary_single_pt().
 * If (i0,i1,i2) *is* the same as (i0_inbounds,i1_inbounds,i2_inbounds),
 * then we are at an outer boundary point. Take care of outer BCs in Step 2.
 * Step 2: Set up outer boundary structs bcstruct->outer_bc_array[which_gz][face][idx2d]:
 * Recall that at each inner boundary point we must set outerpt_bc_struct:
 * typedef struct __outerpt_bc_struct__ {
 * short i0,i1,i2;  // the outer boundary point grid index (i0,i1,i2), on the 3D grid
 * int8_t FACEX0,FACEX1,FACEX2;  // 1-byte integers that store
 * //                               FACEX0,FACEX1,FACEX2 = +1, 0, 0 if on the i0=i0min face,
 * //                               FACEX0,FACEX1,FACEX2 = -1, 0, 0 if on the i0=i0max face,
 * //                               FACEX0,FACEX1,FACEX2 =  0,+1, 0 if on the i1=i2min face,
 * //                               FACEX0,FACEX1,FACEX2 =  0,-1, 0 if on the i1=i1max face,
 * //                               FACEX0,FACEX1,FACEX2 =  0, 0,+1 if on the i2=i2min face, or
 * //                               FACEX0,FACEX1,FACEX2 =  0, 0,-1 if on the i2=i2max face,
 * } outerpt_bc_struct;
 * Outer boundary points are filled from the inside out, two faces at a time.
 * E.g., consider a Cartesian coordinate grid that has 14 points in each direction,
 * including the ghostzones, with NGHOSTS=2.
 * We first fill in the lower x0 face with (i0=1,i1={2,11},i2={2,11}). We fill these
 * points in first, since they will in general (at least in the case of extrapolation
 * outer BCs) depend on e.g., i0=2 and i0=3 points.
 * Simultaneously we can fill in the upper x0 face with (i0=12,i1={2,11},i2={2,11}),
 * since these points depend only on e.g., i0=11 and i0=10 (again assuming extrap. BCs).
 * Next we can fill in the lower x1 face: (i0={1,12},i1=2,i2={2,11}). Notice these
 * depend on i0 min and max faces being filled. The remaining pattern goes like this:
 * Upper x1 face: (i0={1,12},i1=12,i2={2,11})
 * Lower x2 face: (i0={1,12},i1={1,12},i2=1)
 * Upper x2 face: (i0={1,12},i1={1,12},i2=12)
 * Lower x0 face: (i0=0,i1={1,12},i2={1,12})
 * Upper x0 face: (i0=13,i1={1,12},i2={1,12})
 * Lower x1 face: (i0={0,13},i1=0,i2={2,11})
 * Upper x1 face: (i0={0,13},i1=13,i2={2,11})
 * Lower x2 face: (i0={0,13},i1={0,13},i2=0)
 * Upper x2 face: (i0={0,13},i1={0,13},i2=13)
 * Note that we allocate a outerpt_bc_struct at *all* boundary points,
 * regardless of whether the point is an outer or inner point. However
 * the struct is set only at outer boundary points. This is slightly
 * wasteful, but only in memory, not in CPU.
 */
void bcstruct_set_up(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],
                     bc_struct *restrict bcstruct) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    bcstruct_set_up__rfm__SinhSpherical(commondata, params, xx, bcstruct);
    break;
  default:
    fprintf(stderr, "ERROR in bcstruct_set_up(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}
