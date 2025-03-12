#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/**
 * Map a reference metric grid point index (i0, i1, i2) in a ghost zone to an interior point index (i0, i1, i2)':
 * (i0, i1, i2) -> (i0, i1, i2)',
 * if it is an inner boundary point. If the grid point maps to itself; i.e.,
 * (i0, i1, i2) -> (i0, i1, i2),
 * it should have been marked as an outer boundary point. This process involves the following double-map:
 * (x0, x1, x2) -> (Cartx, Carty, Cartz) -> (x0, x1, x2)'
 * However, the second map from Cartesian to the reference metric does not always have a closed-form expression,
 * and this simple algorithm will fail. To address this, we exploit the fact that an arbitrary reference metric and
 * its "eigencoordinate" share the exact same index mapping:
 * (i0, i1, i2) -> (i0, i1, i2)'
 * Therefore, while the mapping
 * (Cartx, Carty, Cartz) -> (x0, x1, x2)'
 * may not be closed-form for the chosen CoordSystem, the eigencoordinate mapping
 * (Cartx, Carty, Cartz) -> (x0, x1, x2)'
 * will be, for all reference metrics in NRPy.
 *
 * Definition of Eigencoordinate:
 *
 * A coordinate system's "eigencoordinate" is the simplest member of its family:
 * - Spherical-like: Spherical
 * - Cylindrical-like: Cylindrical
 * - Cartesian-like: Cartesian
 * - SymTP-like: SymTP
 *
 * Key Steps:
 * 1. Convert to Cartesian Coordinates:
 * - Transform (x0, x1, x2) from eigencoordinates to Cartesian (Cartx, Carty, Cartz).
 * 2. Map Back to Eigencoordinates:
 * - Convert (Cartx, Carty, Cartz) back to eigencoordinates (x0', x1', x2').
 * 3. Sanity Check and Data Handling:
 * - If (x0, x1, x2) != (x0', x1', x2'), the point is an inner boundary point.
 * - For example, in Spherical coordinates, a negative radius becomes positive.
 * - On a cell-centered grid, the mapped point lies within the grid interior.
 * - Update the data at (i0, i1, i2) by copying from (i0_inbounds, i1_inbounds, i2_inbounds).
 * - Apply a sign change (+1 or -1) if the data represents tensors or vectors.
 *
 */
static void EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(const commondata_struct *restrict commondata,
                                                                      const params_struct *restrict params, REAL *restrict xx[3], const int i0,
                                                                      const int i1, const int i2, REAL x0x1x2_inbounds[3], int i0i1i2_inbounds[3]) {
#include "../set_CodeParameters.h"

  // Step 1: Convert the (curvilinear) coordinate (x0,x1,x2) to Cartesian coordinates:
  //         (x0,x1,x2) -> (Cartx,Carty,Cartz)
  //         Find the Cartesian coordinate that (x0,x1,x2)
  //         maps to, assuming (x0,x1,x2) is the eigen-
  //         coordinate. Note that we assume (x0,x1,x2)
  //         has the same grid boundaries in both the
  //         original coordinate and the eigencoordinate.
  DOUBLE xCart[3]; // where (x,y,z) is output
  {
    // xx_to_Cart for EigenCoordinate Spherical (orig coord = SinhSpherical):
    DOUBLE xx0 = xx[0][i0];
    DOUBLE xx1 = xx[1][i1];
    DOUBLE xx2 = xx[2][i2];
    /*
     *  Original SymPy expressions:
     *  "[xCart[0] = xx0*sin(xx1)*cos(xx2)]"
     *  "[xCart[1] = xx0*sin(xx1)*sin(xx2)]"
     *  "[xCart[2] = xx0*cos(xx1)]"
     */
    {
      const REAL tmp0 = xx0 * sin(xx1);
      xCart[0] = tmp0 * cos(xx2);
      xCart[1] = tmp0 * sin(xx2);
      xCart[2] = xx0 * cos(xx1);
    }
  }
  DOUBLE Cartx = xCart[0];
  DOUBLE Carty = xCart[1];
  DOUBLE Cartz = xCart[2];

  // Step 2: Convert Cartesian coordinates (Cartx, Carty, Cartz) to eigencoordinates (x0, x1, x2)
  //         and determine the corresponding gridpoint indices (i0, i1, i2).
  //         - For cell-centered grids, (x0, x1, x2) will exactly (within roundoff error) align
  //           with a numerical grid point.
  //         - Identify the in-bounds indices (i0_inbounds, i1_inbounds, i2_inbounds) that
  //           correspond to the eigencoordinates.
  //         - Check the location of (i0_inbounds, i1_inbounds, i2_inbounds):
  //             - If it resides within a ghost zone:
  //                 - It must be equal to (i0, i1, i2), indicating that the point is an outer boundary point.
  //             - Otherwise:
  //                 - The indices are within the grid interior.
  //                 - Replace the data at (i0, i1, i2) with the data from (i0_inbounds, i1_inbounds, i2_inbounds),
  //                   multiplied by the appropriate parity condition (+1 or -1).
  DOUBLE Cart_to_xx0_inbounds, Cart_to_xx1_inbounds, Cart_to_xx2_inbounds;

  // Cart_to_xx for EigenCoordinate Spherical (orig coord = SinhSpherical)
  /*
   *  Original SymPy expressions:
   *  "[Cart_to_xx0_inbounds = sqrt(Cartx**2 + Carty**2 + Cartz**2)]"
   *  "[Cart_to_xx1_inbounds = acos(Cartz/sqrt(Cartx**2 + Carty**2 + Cartz**2))]"
   *  "[Cart_to_xx2_inbounds = atan2(Carty, Cartx)]"
   */
  {
    const REAL tmp0 = sqrt(((Cartx) * (Cartx)) + ((Carty) * (Carty)) + ((Cartz) * (Cartz)));
    Cart_to_xx0_inbounds = tmp0;
    Cart_to_xx1_inbounds = acos(Cartz / tmp0);
    Cart_to_xx2_inbounds = atan2(Carty, Cartx);
  }

  // Next compute xxmin[i]. By definition,
  //    xx[i][j] = xxmin[i] + ((DOUBLE)(j-NGHOSTS) + (1.0/2.0))*dxxi;
  // -> xxmin[i] = xx[i][0] - ((DOUBLE)(0-NGHOSTS) + (1.0/2.0))*dxxi
  const DOUBLE xxmin[3] = {xx[0][0] - ((DOUBLE)(0 - NGHOSTS) + (1.0 / 2.0)) * dxx0, xx[1][0] - ((DOUBLE)(0 - NGHOSTS) + (1.0 / 2.0)) * dxx1,
                           xx[2][0] - ((DOUBLE)(0 - NGHOSTS) + (1.0 / 2.0)) * dxx2};

  // Finally compute i{0,1,2}_inbounds (add 0.5 to account for rounding down)
  const int i0_inbounds = (int)((Cart_to_xx0_inbounds - xxmin[0] - (1.0 / 2.0) * dxx0 + ((DOUBLE)NGHOSTS) * dxx0) / dxx0 + 0.5);
  const int i1_inbounds = (int)((Cart_to_xx1_inbounds - xxmin[1] - (1.0 / 2.0) * dxx1 + ((DOUBLE)NGHOSTS) * dxx1) / dxx1 + 0.5);
  const int i2_inbounds = (int)((Cart_to_xx2_inbounds - xxmin[2] - (1.0 / 2.0) * dxx2 + ((DOUBLE)NGHOSTS) * dxx2) / dxx2 + 0.5);

  // Step 3: Sanity Check
  //         - Convert eigencoordinates x0(i0_inbounds), x1(i1_inbounds), x2(i2_inbounds)
  //           to Cartesian coordinates (Cartx, Carty, Cartz).
  //         - Verify that the converted coordinates match the original mapping:
  //               (Cartx, Carty, Cartz) == (Cartx(x0(i0)), Carty(x1(i1)), Cartz(x2(i2))).
  //         - If the coordinates do not match, trigger an error.

  // Step 3.a: Compute {x,y,z}Cart_from_xx, as a function of i0,i1,i2
  DOUBLE xCart_from_xx, yCart_from_xx, zCart_from_xx;
  {
    // xx_to_Cart for Coordinate SinhSpherical):
    DOUBLE xx0 = xx[0][i0];
    DOUBLE xx1 = xx[1][i1];
    DOUBLE xx2 = xx[2][i2];
    /*
     *  Original SymPy expressions:
     *  "[xCart_from_xx = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))*sin(xx1)*cos(xx2)/(exp(1/SINHW) - exp(-1/SINHW))]"
     *  "[yCart_from_xx = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))*sin(xx1)*sin(xx2)/(exp(1/SINHW) - exp(-1/SINHW))]"
     *  "[zCart_from_xx = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))*cos(xx1)/(exp(1/SINHW) - exp(-1/SINHW))]"
     */
    const REAL tmp0 = (1.0 / (SINHW));
    const REAL tmp2 = AMPL * (exp(tmp0 * xx0) - exp(-tmp0 * xx0)) / (exp(tmp0) - exp(-tmp0));
    const REAL tmp3 = tmp2 * sin(xx1);
    xCart_from_xx = tmp3 * cos(xx2);
    yCart_from_xx = tmp3 * sin(xx2);
    zCart_from_xx = tmp2 * cos(xx1);
  }

  // Step 3.b: Compute {x,y,z}Cart_from_xx_inbounds, as a
  //           function of i0_inbounds,i1_inbounds,i2_inbounds
  DOUBLE xCart_from_xx_inbounds, yCart_from_xx_inbounds, zCart_from_xx_inbounds;
  {
    // xx_to_Cart_inbounds for Coordinate SinhSpherical):
    DOUBLE xx0 = xx[0][i0_inbounds];
    DOUBLE xx1 = xx[1][i1_inbounds];
    DOUBLE xx2 = xx[2][i2_inbounds];
    /*
     *  Original SymPy expressions:
     *  "[xCart_from_xx_inbounds = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))*sin(xx1)*cos(xx2)/(exp(1/SINHW) - exp(-1/SINHW))]"
     *  "[yCart_from_xx_inbounds = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))*sin(xx1)*sin(xx2)/(exp(1/SINHW) - exp(-1/SINHW))]"
     *  "[zCart_from_xx_inbounds = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))*cos(xx1)/(exp(1/SINHW) - exp(-1/SINHW))]"
     */
    const REAL tmp0 = (1.0 / (SINHW));
    const REAL tmp2 = AMPL * (exp(tmp0 * xx0) - exp(-tmp0 * xx0)) / (exp(tmp0) - exp(-tmp0));
    const REAL tmp3 = tmp2 * sin(xx1);
    xCart_from_xx_inbounds = tmp3 * cos(xx2);
    yCart_from_xx_inbounds = tmp3 * sin(xx2);
    zCart_from_xx_inbounds = tmp2 * cos(xx1);
  }

  // Step 3.c: Compare xCart_from_xx to xCart_from_xx_inbounds;
  //           they should be identical!!!
#define EPS_REL 1e-8
  const DOUBLE norm_factor = sqrt(xCart_from_xx * xCart_from_xx + yCart_from_xx * yCart_from_xx + zCart_from_xx * zCart_from_xx) + 1e-15;
  if (fabs((DOUBLE)(xCart_from_xx - xCart_from_xx_inbounds)) > EPS_REL * norm_factor ||
      fabs((DOUBLE)(yCart_from_xx - yCart_from_xx_inbounds)) > EPS_REL * norm_factor ||
      fabs((DOUBLE)(zCart_from_xx - zCart_from_xx_inbounds)) > EPS_REL * norm_factor) {
    fprintf(stderr,
            "Error in SinhSpherical coordinate system: Inner boundary point does not map to grid interior point: ( %.15e %.15e %.15e ) != ( %.15e "
            "%.15e %.15e ) | xx: %e %e %e -> %e %e %e | %d %d %d\n",
            (DOUBLE)xCart_from_xx, (DOUBLE)yCart_from_xx, (DOUBLE)zCart_from_xx, (DOUBLE)xCart_from_xx_inbounds, (DOUBLE)yCart_from_xx_inbounds,
            (DOUBLE)zCart_from_xx_inbounds, xx[0][i0], xx[1][i1], xx[2][i2], xx[0][i0_inbounds], xx[1][i1_inbounds], xx[2][i2_inbounds],
            Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2);
    exit(1);
  }
#undef EPS_REL

  // Step 4: Set output arrays.
  x0x1x2_inbounds[0] = xx[0][i0_inbounds];
  x0x1x2_inbounds[1] = xx[1][i1_inbounds];
  x0x1x2_inbounds[2] = xx[2][i2_inbounds];
  i0i1i2_inbounds[0] = i0_inbounds;
  i0i1i2_inbounds[1] = i1_inbounds;
  i0i1i2_inbounds[2] = i2_inbounds;
} // END FUNCTION EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt
/**
 * set_parity_for_inner_boundary_single_pt():
 * Given (x0,x1,x2)=(xx0,xx1,xx2) and
 * (x0,x1,x2)'=(x0x1x2_inbounds[0],x0x1x2_inbounds[1],x0x1x2_inbounds[2])
 * (see description of
 * EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt()
 * above for more details), here we compute the parity conditions
 * for all 10 tensor types supported by NRPy+.
 */
static void set_parity_for_inner_boundary_single_pt(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                    const REAL xx0, const REAL xx1, const REAL xx2, const REAL x0x1x2_inbounds[3], const int idx,
                                                    innerpt_bc_struct *restrict innerpt_bc_arr) {
#include "../set_CodeParameters.h"

#define EPS_REL 1e-8

  const DOUBLE xx0_inbounds = x0x1x2_inbounds[0];
  const DOUBLE xx1_inbounds = x0x1x2_inbounds[1];
  const DOUBLE xx2_inbounds = x0x1x2_inbounds[2];

  DOUBLE REAL_parity_array[10];
  {
    // Evaluate dot products needed for setting parity
    //     conditions at a given point (xx0,xx1,xx2),
    //     using C code generated by NRPy+
    /*
    NRPy+ Curvilinear Boundary Conditions: Unit vector dot products for all
         ten parity conditions, in given coordinate system.
         Needed for automatically determining sign of tensor across coordinate boundary.
    Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb
    */
    /*
     *  Original SymPy expressions:
     *  "[REAL_parity_array[0] = 1]"
     *  "[REAL_parity_array[1] = sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds) + sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds) +
     * cos(xx1)*cos(xx1_inbounds)]"
     *  "[REAL_parity_array[2] = sin(xx1)*sin(xx1_inbounds) + sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds) +
     * cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)]"
     *  "[REAL_parity_array[3] = sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds)]"
     *  "[REAL_parity_array[4] = (sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds) + sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds) +
     * cos(xx1)*cos(xx1_inbounds))**2]"
     *  "[REAL_parity_array[5] = (sin(xx1)*sin(xx1_inbounds) + sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds) +
     * cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds))*(sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds) +
     * sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds) + cos(xx1)*cos(xx1_inbounds))]"
     *  "[REAL_parity_array[6] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))*(sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds) +
     * sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds) + cos(xx1)*cos(xx1_inbounds))]"
     *  "[REAL_parity_array[7] = (sin(xx1)*sin(xx1_inbounds) + sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds) +
     * cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds))**2]"
     *  "[REAL_parity_array[8] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))*(sin(xx1)*sin(xx1_inbounds) +
     * sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds) + cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds))]"
     *  "[REAL_parity_array[9] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))**2]"
     */
    {
      const REAL tmp0 = cos(xx1) * cos(xx1_inbounds);
      const REAL tmp1 = sin(xx1) * sin(xx1_inbounds);
      const REAL tmp2 = sin(xx2) * sin(xx2_inbounds);
      const REAL tmp3 = cos(xx2) * cos(xx2_inbounds);
      const REAL tmp4 = tmp0 + tmp1 * tmp2 + tmp1 * tmp3;
      const REAL tmp5 = tmp0 * tmp2 + tmp0 * tmp3 + tmp1;
      const REAL tmp6 = tmp2 + tmp3;
      REAL_parity_array[0] = 1;
      REAL_parity_array[1] = tmp4;
      REAL_parity_array[2] = tmp5;
      REAL_parity_array[3] = tmp6;
      REAL_parity_array[4] = ((tmp4) * (tmp4));
      REAL_parity_array[5] = tmp4 * tmp5;
      REAL_parity_array[6] = tmp4 * tmp6;
      REAL_parity_array[7] = ((tmp5) * (tmp5));
      REAL_parity_array[8] = tmp5 * tmp6;
      REAL_parity_array[9] = ((tmp6) * (tmp6));
    }
  }

  // Next perform sanity check on parity array output: should be +1 or -1 to within 8 significant digits:
  for (int whichparity = 0; whichparity < 10; whichparity++) {
    if (fabs(REAL_parity_array[whichparity]) < 1 - EPS_REL || fabs(REAL_parity_array[whichparity]) > 1 + EPS_REL) {
      fprintf(stderr, "Error at point (%e %e %e), which maps to (%e %e %e).\n", xx0, xx1, xx2, xx0_inbounds, xx1_inbounds, xx2_inbounds);
      fprintf(stderr, "Parity evaluated to %e , which is not within 8 significant digits of +1 or -1.\n", REAL_parity_array[whichparity]);
      exit(1);
    }
    for (int parity = 0; parity < 10; parity++) {
      innerpt_bc_arr[idx].parity[parity] = 1;
      if (REAL_parity_array[parity] < 0)
        innerpt_bc_arr[idx].parity[parity] = -1;
    }
  } // END for(int whichparity=0;whichparity<10;whichparity++)
#undef EPS_REL
} // END FUNCTION set_parity_for_inner_boundary_single_pt

/**
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
void bcstruct_set_up__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],
                                         bc_struct *restrict bcstruct) {
#include "../set_CodeParameters.h"

  ////////////////////////////////////////
  // STEP 1: SET UP INNER BOUNDARY STRUCTS
  {
    // First count the number of inner points.
    int num_inner = 0;
    LOOP_OMP("omp parallel for reduction(+:num_inner)", i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
      const int i0i1i2[3] = {i0, i1, i2};
      if (!IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, NGHOSTS)) {
        REAL x0x1x2_inbounds[3];
        int i0i1i2_inbounds[3];
        EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, params, xx, i0, i1, i2, x0x1x2_inbounds, i0i1i2_inbounds);
        if (i0 == i0i1i2_inbounds[0] && i1 == i0i1i2_inbounds[1] && i2 == i0i1i2_inbounds[2]) {
          // this is a pure outer boundary point.
        } else {
          // this is an inner boundary point, which maps either
          //  to the grid interior or to an outer boundary point
          num_inner++;
        }
      }
    }
    // Store num_inner to bc_info:
    bcstruct->bc_info.num_inner_boundary_points = num_inner;

    // Next allocate memory for inner_boundary_points:
    bcstruct->inner_bc_array = (innerpt_bc_struct *restrict)malloc(sizeof(innerpt_bc_struct) * num_inner);
  }

  // Then set inner_bc_array:
  {
    int which_inner = 0;
    LOOP_NOOMP(i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
      const int i0i1i2[3] = {i0, i1, i2};
      if (!IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, NGHOSTS)) {
        REAL x0x1x2_inbounds[3];
        int i0i1i2_inbounds[3];
        EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, params, xx, i0, i1, i2, x0x1x2_inbounds, i0i1i2_inbounds);
        if (i0 == i0i1i2_inbounds[0] && i1 == i0i1i2_inbounds[1] && i2 == i0i1i2_inbounds[2]) {
          // this is a pure outer boundary point.
        } else {
          bcstruct->inner_bc_array[which_inner].dstpt = IDX3(i0, i1, i2);
          bcstruct->inner_bc_array[which_inner].srcpt = IDX3(i0i1i2_inbounds[0], i0i1i2_inbounds[1], i0i1i2_inbounds[2]);
          // printf("%d / %d\n",which_inner, bc_info->num_inner_boundary_points);
          set_parity_for_inner_boundary_single_pt(commondata, params, xx[0][i0], xx[1][i1], xx[2][i2], x0x1x2_inbounds, which_inner,
                                                  bcstruct->inner_bc_array);

          which_inner++;
        }
      }
    }
  }

  ////////////////////////////////////////
  // STEP 2: SET UP OUTER BOUNDARY STRUCTS
  // First set up loop bounds for outer boundary condition updates,
  //   store to bc_info->bc_loop_bounds[which_gz][face][]. Also
  //   allocate memory for outer_bc_array[which_gz][face][]:
  int imin[3] = {NGHOSTS, NGHOSTS, NGHOSTS};
  int imax[3] = {Nxx_plus_2NGHOSTS0 - NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS};
  for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
    const int x0min_face_range[6] = {imin[0] - 1, imin[0], imin[1], imax[1], imin[2], imax[2]};
    imin[0]--;
    const int x0max_face_range[6] = {imax[0], imax[0] + 1, imin[1], imax[1], imin[2], imax[2]};
    imax[0]++;
    const int x1min_face_range[6] = {imin[0], imax[0], imin[1] - 1, imin[1], imin[2], imax[2]};
    imin[1]--;
    const int x1max_face_range[6] = {imin[0], imax[0], imax[1], imax[1] + 1, imin[2], imax[2]};
    imax[1]++;
    const int x2min_face_range[6] = {imin[0], imax[0], imin[1], imax[1], imin[2] - 1, imin[2]};
    imin[2]--;
    const int x2max_face_range[6] = {imin[0], imax[0], imin[1], imax[1], imax[2], imax[2] + 1};
    imax[2]++;

    int face = 0;
    ////////////////////////
    // x0min and x0max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x0min and x0max faces have exactly the same size.
    //                   Also, note that face/2 --v   offsets this factor of 2 ------------------------------------------v
    bcstruct->pure_outer_bc_array[3 * which_gz + face / 2] = (outerpt_bc_struct *restrict)malloc(
        sizeof(outerpt_bc_struct) * 2 *
        ((x0min_face_range[1] - x0min_face_range[0]) * (x0min_face_range[3] - x0min_face_range[2]) * (x0min_face_range[5] - x0min_face_range[4])));
    // x0min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x0min_face_range[i];
    }
    face++;
    // x0max face: Set loop bounds & allocate memory for outer_bc_array:
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x0max_face_range[i];
    }
    face++;
    ////////////////////////

    ////////////////////////
    // x1min and x1max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x1min and x1max faces have exactly the same size.
    //                   Also, note that face/2 --v   offsets this factor of 2 ------------------------------------------v
    bcstruct->pure_outer_bc_array[3 * which_gz + face / 2] = (outerpt_bc_struct *restrict)malloc(
        sizeof(outerpt_bc_struct) * 2 *
        ((x1min_face_range[1] - x1min_face_range[0]) * (x1min_face_range[3] - x1min_face_range[2]) * (x1min_face_range[5] - x1min_face_range[4])));
    // x1min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x1min_face_range[i];
    }
    face++;
    // x1max face: Set loop bounds & allocate memory for outer_bc_array:
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x1max_face_range[i];
    }
    face++;
    ////////////////////////

    ////////////////////////
    // x2min and x2max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x2min and x2max faces have exactly the same size.
    //                   Also, note that face/2 --v   offsets this factor of 2 ------------------------------------------v
    bcstruct->pure_outer_bc_array[3 * which_gz + face / 2] = (outerpt_bc_struct *restrict)malloc(
        sizeof(outerpt_bc_struct) * 2 *
        ((x2min_face_range[1] - x2min_face_range[0]) * (x2min_face_range[3] - x2min_face_range[2]) * (x2min_face_range[5] - x2min_face_range[4])));
    // x2min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x2min_face_range[i];
    }
    face++;
    // x2max face: Set loop bounds & allocate memory for outer_bc_array:
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x2max_face_range[i];
    }
    face++;
    ////////////////////////
  }

  for (int which_gz = 0; which_gz < NGHOSTS; which_gz++)
    for (int dirn = 0; dirn < 3; dirn++) {
      int idx2d = 0;
      // LOWER FACE: dirn=0 -> x0min; dirn=1 -> x1min; dirn=2 -> x2min
      {
        const int face = dirn * 2;
#define IDX2D_BCS(i0, i0min, i0max, i1, i1min, i1max, i2, i2min, i2max)                                                                              \
  (((i0) - (i0min)) + ((i0max) - (i0min)) * (((i1) - (i1min)) + ((i1max) - (i1min)) * ((i2) - (i2min))))
        const int FACEX0 = (face == 0) - (face == 1); // +1 if face==0 (x0min) ; -1 if face==1 (x0max). Otherwise 0.
        const int FACEX1 = (face == 2) - (face == 3); // +1 if face==2 (x1min) ; -1 if face==3 (x1max). Otherwise 0.
        const int FACEX2 = (face == 4) - (face == 5); // +1 if face==4 (x2min) ; -1 if face==5 (x2max). Otherwise 0.
        LOOP_NOOMP(i0, bcstruct->bc_info.bc_loop_bounds[which_gz][face][0], bcstruct->bc_info.bc_loop_bounds[which_gz][face][1], i1,
                   bcstruct->bc_info.bc_loop_bounds[which_gz][face][2], bcstruct->bc_info.bc_loop_bounds[which_gz][face][3], i2,
                   bcstruct->bc_info.bc_loop_bounds[which_gz][face][4], bcstruct->bc_info.bc_loop_bounds[which_gz][face][5]) {
          REAL x0x1x2_inbounds[3];
          int i0i1i2_inbounds[3];
          EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, params, xx, i0, i1, i2, x0x1x2_inbounds, i0i1i2_inbounds);
          if (i0 == i0i1i2_inbounds[0] && i1 == i0i1i2_inbounds[1] && i2 == i0i1i2_inbounds[2]) {
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i0 = i0;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i1 = i1;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i2 = i2;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX0 = FACEX0;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX1 = FACEX1;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX2 = FACEX2;
            idx2d++;
          }
        }
      }
      // UPPER FACE: dirn=0 -> x0max; dirn=1 -> x1max; dirn=2 -> x2max
      {
        const int face = dirn * 2 + 1;
        const int FACEX0 = (face == 0) - (face == 1); // +1 if face==0 ; -1 if face==1. Otherwise 0.
        const int FACEX1 = (face == 2) - (face == 3); // +1 if face==2 ; -1 if face==3. Otherwise 0.
        const int FACEX2 = (face == 4) - (face == 5); // +1 if face==4 ; -1 if face==5. Otherwise 0.
        LOOP_NOOMP(i0, bcstruct->bc_info.bc_loop_bounds[which_gz][face][0], bcstruct->bc_info.bc_loop_bounds[which_gz][face][1], i1,
                   bcstruct->bc_info.bc_loop_bounds[which_gz][face][2], bcstruct->bc_info.bc_loop_bounds[which_gz][face][3], i2,
                   bcstruct->bc_info.bc_loop_bounds[which_gz][face][4], bcstruct->bc_info.bc_loop_bounds[which_gz][face][5]) {
          REAL x0x1x2_inbounds[3];
          int i0i1i2_inbounds[3];
          EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, params, xx, i0, i1, i2, x0x1x2_inbounds, i0i1i2_inbounds);
          if (i0 == i0i1i2_inbounds[0] && i1 == i0i1i2_inbounds[1] && i2 == i0i1i2_inbounds[2]) {
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i0 = i0;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i1 = i1;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i2 = i2;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX0 = FACEX0;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX1 = FACEX1;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX2 = FACEX2;
            idx2d++;
          }
        }
      }
      bcstruct->bc_info.num_pure_outer_boundary_points[which_gz][dirn] = idx2d;
    }
} // END FUNCTION bcstruct_set_up__rfm__SinhSpherical
