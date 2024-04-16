"""
Module providing functions for setting up Curvilinear boundary conditions.

This is documented in Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
"""

# Step P1: Import needed NRPy+ core modules:
from typing import List, Tuple
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import sympy.codegen.ast as sp_ast
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.params as par  # NRPy+: Parameter interface
import nrpy.grid as gri  # NRPy+: Functions having to do with numerical grids
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support
import nrpy.finite_difference as fin  # NRPy+: Finite-difference module
from nrpy.infrastructures.BHaH import griddata_commondata
from nrpy.infrastructures.BHaH import BHaH_defines_h
from nrpy.validate_expressions.validate_expressions import check_zero
import nrpy.infrastructures.BHaH.CurviBoundaryConditions.base_CurviBoundaryConditions as base_cbc_classes

_ = par.CodeParameter(
    "char[50]", __name__, "outer_bc_type", "radiation", commondata=True
)

core_modules_list = [
    "general",
    "nrpy.infrastructures.BHaH.diagnostics.progress_indicator",
    "commondata_struct",
    "params_struct",
    "finite_difference",
    "reference_metric",
    "nrpy.infrastructures.BHaH.CurviBoundaryConditions.openmp.CurviBoundaryConditions",
    "nrpy.infrastructures.BHaH.MoLtimestepping.MoL",
    "nrpy.infrastructures.BHaH.interpolation.interpolation",
    "grid",  # griddata struct depends upon other core modules
]


# bcstruct_set_up():
#      This function is documented in desc= and body= fields below.
class register_CFunction_bcstruct_set_up(
    base_cbc_classes.base_register_CFunction_bcstruct_set_up
):
    def __init__(self, CoordSystem: str, fp_type: str = "double") -> None:
        """
        Register C function for setting up bcstruct.

        This function prescribes how inner and outer boundary points on the
        computational grid are filled, based on the given coordinate system (CoordSystem).

        :param CoordSystem: The coordinate system for which to set up boundary conditions.
        :param fp_type: Floating point type, e.g., "double".
        """
        super().__init__(CoordSystem, fp_type=fp_type)

        self.body = r"""
  ////////////////////////////////////////
  // STEP 1: SET UP INNER BOUNDARY STRUCTS
  {
    // First count the number of inner points.
    int num_inner = 0;
    LOOP_OMP("omp parallel for reduction(+:num_inner)",
             i0,0,Nxx_plus_2NGHOSTS0,  i1,0,Nxx_plus_2NGHOSTS1,  i2,0,Nxx_plus_2NGHOSTS2) {
      const int i0i1i2[3] = { i0,i1,i2 };
      if(!IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2, NGHOSTS)) {
        REAL x0x1x2_inbounds[3];
        int i0i1i2_inbounds[3];
        EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
        if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
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
    bcstruct->inner_bc_array = (innerpt_bc_struct *restrict)malloc( sizeof(innerpt_bc_struct)*num_inner );
  }

  // Then set inner_bc_array:
  {
    int which_inner = 0;
    LOOP_NOOMP(i0,0,Nxx_plus_2NGHOSTS0,  i1,0,Nxx_plus_2NGHOSTS1,  i2,0,Nxx_plus_2NGHOSTS2) {
      const int i0i1i2[3] = { i0,i1,i2 };
      if(!IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2, NGHOSTS)) {
        REAL x0x1x2_inbounds[3];
        int i0i1i2_inbounds[3];
        EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
        if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
          // this is a pure outer boundary point.
        } else {
          bcstruct->inner_bc_array[which_inner].dstpt = IDX3(i0,i1,i2);
          bcstruct->inner_bc_array[which_inner].srcpt = IDX3(i0i1i2_inbounds[0],i0i1i2_inbounds[1],i0i1i2_inbounds[2]);
          //printf("%d / %d\n",which_inner, bc_info->num_inner_boundary_points);
          set_parity_for_inner_boundary_single_pt(commondata, params, xx[0][i0],xx[1][i1],xx[2][i2],
                                                  x0x1x2_inbounds, which_inner, bcstruct->inner_bc_array);

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
  int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
  int imax[3] = { Nxx_plus_2NGHOSTS0-NGHOSTS, Nxx_plus_2NGHOSTS1-NGHOSTS, Nxx_plus_2NGHOSTS2-NGHOSTS };
  for(int which_gz=0;which_gz<NGHOSTS;which_gz++) {
    const int x0min_face_range[6] = { imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2] };  imin[0]--;
    const int x0max_face_range[6] = { imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2] };  imax[0]++;
    const int x1min_face_range[6] = { imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2] };  imin[1]--;
    const int x1max_face_range[6] = { imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2] };  imax[1]++;
    const int x2min_face_range[6] = { imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2] };  imin[2]--;
    const int x2max_face_range[6] = { imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1 };  imax[2]++;

    int face=0;
    ////////////////////////
    // x0min and x0max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x0min and x0max faces have exactly the same size.
    //                   Also, note that face/2 --v   offsets this factor of 2 ------------------------------------------v
    bcstruct->pure_outer_bc_array[3*which_gz + face/2] = (outerpt_bc_struct *restrict)malloc(sizeof(outerpt_bc_struct) * 2 *
                                                                                             ((x0min_face_range[1]-x0min_face_range[0]) *
                                                                                              (x0min_face_range[3]-x0min_face_range[2]) *
                                                                                              (x0min_face_range[5]-x0min_face_range[4])));
    // x0min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x0min_face_range[i]; }
    face++;
    // x0max face: Set loop bounds & allocate memory for outer_bc_array:
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x0max_face_range[i]; }
    face++;
    ////////////////////////

    ////////////////////////
    // x1min and x1max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x1min and x1max faces have exactly the same size.
    //                   Also, note that face/2 --v   offsets this factor of 2 ------------------------------------------v
    bcstruct->pure_outer_bc_array[3*which_gz + face/2] = (outerpt_bc_struct *restrict)malloc(sizeof(outerpt_bc_struct) * 2 *
                                                                                             ((x1min_face_range[1]-x1min_face_range[0]) *
                                                                                              (x1min_face_range[3]-x1min_face_range[2]) *
                                                                                              (x1min_face_range[5]-x1min_face_range[4])));
    // x1min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x1min_face_range[i]; }
    face++;
    // x1max face: Set loop bounds & allocate memory for outer_bc_array:
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x1max_face_range[i]; }
    face++;
    ////////////////////////


    ////////////////////////
    // x2min and x2max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x2min and x2max faces have exactly the same size.
    //                   Also, note that face/2 --v   offsets this factor of 2 ------------------------------------------v
    bcstruct->pure_outer_bc_array[3*which_gz + face/2] = (outerpt_bc_struct *restrict)malloc(sizeof(outerpt_bc_struct) * 2 *
                                                                                             ((x2min_face_range[1]-x2min_face_range[0]) *
                                                                                              (x2min_face_range[3]-x2min_face_range[2]) *
                                                                                              (x2min_face_range[5]-x2min_face_range[4])));
    // x2min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x2min_face_range[i]; }
    face++;
    // x2max face: Set loop bounds & allocate memory for outer_bc_array:
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x2max_face_range[i]; }
    face++;
    ////////////////////////
  }

  for(int which_gz=0;which_gz<NGHOSTS;which_gz++) for(int dirn=0;dirn<3;dirn++) {
      int idx2d = 0;
      // LOWER FACE: dirn=0 -> x0min; dirn=1 -> x1min; dirn=2 -> x2min
      {
        const int face = dirn*2;
#define IDX2D_BCS(i0,i0min,i0max, i1,i1min,i1max ,i2,i2min,i2max)       \
        ( ((i0)-(i0min)) + ((i0max)-(i0min)) * ( ((i1)-(i1min)) + ((i1max)-(i1min)) * ((i2)-(i2min)) ) )
        const int FACEX0=(face==0) - (face==1); // +1 if face==0 (x0min) ; -1 if face==1 (x0max). Otherwise 0.
        const int FACEX1=(face==2) - (face==3); // +1 if face==2 (x1min) ; -1 if face==3 (x1max). Otherwise 0.
        const int FACEX2=(face==4) - (face==5); // +1 if face==4 (x2min) ; -1 if face==5 (x2max). Otherwise 0.
        LOOP_NOOMP(i0,bcstruct->bc_info.bc_loop_bounds[which_gz][face][0],bcstruct->bc_info.bc_loop_bounds[which_gz][face][1],
                   i1,bcstruct->bc_info.bc_loop_bounds[which_gz][face][2],bcstruct->bc_info.bc_loop_bounds[which_gz][face][3],
                   i2,bcstruct->bc_info.bc_loop_bounds[which_gz][face][4],bcstruct->bc_info.bc_loop_bounds[which_gz][face][5]) {
          REAL x0x1x2_inbounds[3];
          int i0i1i2_inbounds[3];
          EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
          if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0 = i0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1 = i1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2 = i2;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0 = FACEX0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1 = FACEX1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2 = FACEX2;
            idx2d++;
          }
        }
      }
      // UPPER FACE: dirn=0 -> x0max; dirn=1 -> x1max; dirn=2 -> x2max
      {
        const int face = dirn*2+1;
        const int FACEX0=(face==0) - (face==1); // +1 if face==0 ; -1 if face==1. Otherwise 0.
        const int FACEX1=(face==2) - (face==3); // +1 if face==2 ; -1 if face==3. Otherwise 0.
        const int FACEX2=(face==4) - (face==5); // +1 if face==4 ; -1 if face==5. Otherwise 0.
        LOOP_NOOMP(i0,bcstruct->bc_info.bc_loop_bounds[which_gz][face][0],bcstruct->bc_info.bc_loop_bounds[which_gz][face][1],
                   i1,bcstruct->bc_info.bc_loop_bounds[which_gz][face][2],bcstruct->bc_info.bc_loop_bounds[which_gz][face][3],
                   i2,bcstruct->bc_info.bc_loop_bounds[which_gz][face][4],bcstruct->bc_info.bc_loop_bounds[which_gz][face][5]) {
          REAL x0x1x2_inbounds[3];
          int i0i1i2_inbounds[3];
          EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
          if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0 = i0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1 = i1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2 = i2;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0 = FACEX0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1 = FACEX1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2 = FACEX2;
            idx2d++;
          }
        }
      }
      bcstruct->bc_info.num_pure_outer_boundary_points[which_gz][dirn] = idx2d;
    }
"""
        cfc.register_CFunction(
            includes=self.includes,
            prefunc=self.prefunc,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            CoordSystem_for_wrapper_func=self.CoordSystem,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=True,
            body=self.body,
        )


###############################
## apply_bcs_inner_only(): Apply inner boundary conditions.
##  Function is documented below in desc= and body=.
class register_CFunction_apply_bcs_inner_only(
    base_cbc_classes.base_register_CFunction_apply_bcs_inner_only
):

    def __init__(self) -> None:
        """
        Register C function for filling inner boundary points on the computational grid,
        as prescribed by bcstruct.
        """
        super().__init__()

        self.body = r"""
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

  // collapse(2) results in a nice speedup here, esp in 2D. Two_BHs_collide goes from
  //    5550 M/hr to 7264 M/hr on a Ryzen 9 5950X running on all 16 cores with core affinity.
#pragma omp parallel for collapse(2)  // spawn threads and distribute across them
  for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
    for(int pt=0;pt<bc_info->num_inner_boundary_points;pt++) {
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt;
      const int srcpt = bcstruct->inner_bc_array[pt].srcpt;
      gfs[IDX4pt(which_gf, dstpt)] = bcstruct->inner_bc_array[pt].parity[evol_gf_parity[which_gf]] * gfs[IDX4pt(which_gf, srcpt)];
    } // END for(int pt=0;pt<num_inner_pts;pt++)
  } // END for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++)
"""
        cfc.register_CFunction(
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=True,
            body=self.body,
        )


###############################
## apply_bcs_outerextrap_and_inner(): Apply extrapolation outer boundary conditions.
##  Function is documented below in desc= and body=.
class register_CFunction_apply_bcs_outerextrap_and_inner(
    base_cbc_classes.base_register_CFunction_apply_bcs_outerextrap_and_inner
):

    def __init__(self) -> None:
        """Register C function for filling boundary points with extrapolation and prescribed bcstruct."""
        super().__init__()
        self.body = r"""
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

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
  // Spawn N OpenMP threads, either across all cores, or according to e.g., taskset.
#pragma omp parallel
  {
    for(int which_gz=0;which_gz<NGHOSTS;which_gz++) for(int dirn=0;dirn<3;dirn++) {
        // This option results in about 1.6% slower runtime for SW curvilinear at 64x24x24 on 8-core Ryzen 9 4900HS
        //#pragma omp for collapse(2)
        //for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++) {
        //  {
        // Don't spawn a thread if there are no boundary points to fill; results in a nice little speedup.
        if(bc_info->num_pure_outer_boundary_points[which_gz][dirn] > 0) {
#pragma omp for  // threads have been spawned; here we distribute across them
          for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++) {
            const short i0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0;
            const short i1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1;
            const short i2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2;
            const short FACEX0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0;
            const short FACEX1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1;
            const short FACEX2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2;
            const int idx_offset0 = IDX3(i0,i1,i2);
            const int idx_offset1 = IDX3(i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2);
            const int idx_offset2 = IDX3(i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2);
            const int idx_offset3 = IDX3(i0+3*FACEX0,i1+3*FACEX1,i2+3*FACEX2);
            for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
              // *** Apply 2nd-order polynomial extrapolation BCs to all outer boundary points. ***
              gfs[IDX4pt(which_gf, idx_offset0)] =
                +3.0*gfs[IDX4pt(which_gf, idx_offset1)]
                -3.0*gfs[IDX4pt(which_gf, idx_offset2)]
                +1.0*gfs[IDX4pt(which_gf, idx_offset3)];
            }
          }
        }
      }
  }

  ///////////////////////////////////////////////////////
  // STEP 2 of 2: Apply BCs to inner boundary points.
  //              These map to either the grid interior
  //              ("pure inner") or to pure outer boundary
  //              points ("inner maps to outer"). Those
  //              that map to outer require that outer be
  //              populated first; hence this being
  //              STEP 2 OF 2.
  apply_bcs_inner_only(commondata, params, bcstruct, gfs);
"""
        cfc.register_CFunction(
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=True,
            body=self.body,
        )


###############################
## RADIATION (NewRad-like) BOUNDARY CONDITIONS.
##  Functions are fully documented in nrpytutorial's
##   Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb,
##   as well as below, in desc= and body=.
# r_and_partial_xi_partial_r_derivs(): Compute r(x0,x1,x2) and dx^i / dr
def setup_Cfunction_r_and_partial_xi_partial_r_derivs(
    CoordSystem: str, fp_type: str = "double"
) -> str:
    """
    Generate C code to compute the radial coordinate r(x0, x1, x2) and its derivatives.

    Compute the radial coordinate r(x0, x1, x2) and its partial derivatives
    partial x^i / partial r for a given coordinate system.

    :param CoordSystem: The coordinate system for which to compute r and its derivatives.
    :param fp_type: Floating point type, e.g., "double".
    :return: A string containing the generated C code for the function.
    """
    desc = "Compute r(xx0,xx1,xx2) and partial_r x^i."
    cfunc_type = "static inline void"
    name = "r_and_partial_xi_partial_r_derivs"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    const REAL xx0,const REAL xx1,const REAL xx2,    REAL *r,
    REAL *partial_x0_partial_r,REAL *partial_x1_partial_r,REAL *partial_x2_partial_r"""
    rfm = refmetric.reference_metric[CoordSystem]
    # sp.simplify(expr) is too slow here for SinhCylindrical
    body = ccg.c_codegen(
        [
            rfm.xxSph[0],
            rfm.Jac_dUrfm_dDSphUD[0][0],
            rfm.Jac_dUrfm_dDSphUD[1][0],
            rfm.Jac_dUrfm_dDSphUD[2][0],
        ],
        [
            "*r",
            "*partial_x0_partial_r",
            "*partial_x1_partial_r",
            "*partial_x2_partial_r",
        ],
        verbose=False,
        include_braces=False,
        fp_type=fp_type,
    )

    cf = cfc.CFunction(
        subdirectory=CoordSystem,
        includes=[],
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cf.full_function


# partial_r f term: generate finite-difference coefficients
#   for \partial_i f with arbitrary upwinding:
def get_arb_offset_FD_coeffs_indices(
    FDORDER: int, offset: int, deriv: int
) -> Tuple[List[float], List[int]]:
    """
    Generate finite-difference coefficients for partial derivatives with arbitrary upwinding.

    :param FDORDER: Order of the finite difference
    :param offset: Offset for upwinding
    :param deriv: Order of the derivative (e.g., 1 for 1st derivative)

    :return: A tuple containing the list of coefficients and the list of indices

    Doctest:
    >>> get_arb_offset_FD_coeffs_indices(3, 0, 1)
    ([1/24, -9/8, 9/8, -1/24], [-1, 0, 1, 2])
    """
    # deriv = 1 <-- 1st derivative
    Minv = fin.setup_FD_matrix__return_inverse(FDORDER + 1, offset)
    indices = []
    coeffs = []
    for i in range(FDORDER + 1):
        indices.append(i - int(FDORDER / 2) + offset)
        coeffs.append(Minv[i, deriv])
    return coeffs, indices


# partial_r f term: FD1_arbitrary_upwind(): C function to evaluate
#   partial_i f with arbitrary upwinding
def setup_Cfunction_FD1_arbitrary_upwind(
    dirn: int,
    radiation_BC_fd_order: int = -1,
    fp_type: str = "double",
) -> str:
    """
    Set up the C function for computing the 1st derivative finite-difference.

    Supports arbitrary upwind for a given direction and order.

    :param dirn: Direction in which to compute the derivative.
    :param radiation_BC_fd_order: Finite difference order for radiation boundary condition.
                                  If -1, will use default finite difference order.
    :param fp_type: Floating point type, e.g., "double".
    :return: The full C function as a string.
    """
    default_FDORDER = par.parval_from_str("fd_order")
    if radiation_BC_fd_order == -1:
        radiation_BC_fd_order = default_FDORDER

    par.set_parval_from_str("fd_order", radiation_BC_fd_order)

    includes: List[str] = []
    desc = "Compute 1st derivative finite-difference derivative with arbitrary upwind"
    cfunc_type = "static inline REAL"
    name = f"FD1_arbitrary_upwind_x{dirn}_dirn"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
const REAL *restrict gf,  const int i0,const int i1,const int i2, const int offset"""
    body = "switch(offset) {\n"

    tmp_list: List[int] = []
    fp_ccg_type = ccg.fp_type_to_sympy_type[fp_type]
    sp_type_alias = {sp_ast.real: fp_ccg_type}
    for offset in range(
        0, int(radiation_BC_fd_order // 2) + 1
    ):  # Use // for integer division
        tmp_list.append(offset)
        if offset > 0:
            tmp_list.append(-offset)

    for offset in tmp_list:
        body += f"case {offset}:\n"
        body += "  return ("
        coeffs, indices = get_arb_offset_FD_coeffs_indices(
            radiation_BC_fd_order, offset, 1
        )

        for i, coeff in enumerate(coeffs):
            if coeff == 0:
                continue
            offset_str: str = str(indices[i])
            if i > 0:
                body += "          "
            if offset_str == "0":
                body += f"+{sp.ccode(coeff, type_aliases=sp_type_alias)}*gf[IDX3(i0,i1,i2)]\n"
            else:
                if dirn == 0:
                    body += f"+{sp.ccode(coeff, type_aliases=sp_type_alias)}*gf[IDX3(i0+{offset_str},i1,i2)]\n"
                elif dirn == 1:
                    body += f"+{sp.ccode(coeff, type_aliases=sp_type_alias)}*gf[IDX3(i0,i1+{offset_str},i2)]\n"
                elif dirn == 2:
                    body += f"+{sp.ccode(coeff, type_aliases=sp_type_alias)}*gf[IDX3(i0,i1,i2+{offset_str})]\n"

        body = body[:-1].replace("+-", "-") + f") * invdxx{dirn};\n"

    body += """}
return 0.0 / 0.0;  // poison output if offset computed incorrectly
"""

    cf = cfc.CFunction(
        subdirectory="one_subdirectory_down",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )

    par.set_parval_from_str("fd_order", default_FDORDER)

    return cf.full_function


# partial_r f term: Numerically evaluate partial_r f,
#   calling functions defined above.
def setup_Cfunction_compute_partial_r_f(
    CoordSystem: str, radiation_BC_fd_order: int = -1
) -> str:
    """
    Set up a C function for computing the partial derivative of f with respect to r.

    :param CoordSystem: Coordinate system to be used for the computation
    :param radiation_BC_fd_order: Order of finite difference for radiation boundary conditions, default is -1
    :return: A C function for computing the partial derivative
    """
    desc = "Compute \\partial_r f"
    cfunc_type = "static inline REAL"
    name = "compute_partial_r_f"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
REAL *restrict xx[3], const REAL *restrict gfs,
const int which_gf, const int dest_i0,const int dest_i1,const int dest_i2,
const int FACEi0,const int FACEi1,const int FACEi2,
const REAL partial_x0_partial_r, const REAL partial_x1_partial_r, const REAL partial_x2_partial_r"""
    rfm = refmetric.reference_metric[CoordSystem]

    default_FDORDER = par.parval_from_str("fd_order")
    if radiation_BC_fd_order == -1:
        radiation_BC_fd_order = default_FDORDER

    FD1_stencil_radius = int(radiation_BC_fd_order / 2)

    body = f"""  ///////////////////////////////////////////////////////////

  // FD1_stencil_radius = radiation_BC_fd_order/2 = {FD1_stencil_radius}
  const int FD1_stencil_radius = {FD1_stencil_radius};

  const int ntot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

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
"""
    for i in range(3):
        si = str(i)
        if check_zero(rfm.Jac_dUrfm_dDSphUD[i][0]):
            body += f"  const REAL partial_x{si}_f=0.0;\n"
        else:
            body += (
                f"  int i{si}_offset = FACEi{si};  // Shift stencil away from the face we're updating.\n"
                f"  // Next adjust i{si}_offset so that FD stencil never goes out of bounds.\n"
                f"  if(dest_i{si} < FD1_stencil_radius) i{si}_offset = FD1_stencil_radius-dest_i{si};\n"
                f"  else if(dest_i{si} > (Nxx_plus_2NGHOSTS{si}-FD1_stencil_radius-1)) i{si}_offset = (Nxx_plus_2NGHOSTS{si}-FD1_stencil_radius-1) - dest_i{si};\n"
                f"  const REAL partial_x{si}_f=FD1_arbitrary_upwind_x{si}_dirn(commondata, params,&gfs[which_gf*ntot],dest_i0,dest_i1,dest_i2,i{si}_offset);\n"
            )
    body += "  return partial_x0_partial_r*partial_x0_f + partial_x1_partial_r*partial_x1_f + partial_x2_partial_r*partial_x2_f;\n"

    cf = cfc.CFunction(
        subdirectory=CoordSystem,
        includes=[],
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cf.full_function


# radiation_bcs(): Put it all together, for a single outer boundary point.
def setup_Cfunction_radiation_bcs(
    CoordSystem: str,
    radiation_BC_fd_order: int = -1,
    fp_type: str = "double",
) -> str:
    """
    Generate C code to apply radiation boundary conditions in a given coordinate system.

    :param CoordSystem: The coordinate system to use.
    :param radiation_BC_fd_order: Finite differencing order to use. Default is -1.
    :param fp_type: Floating point type, e.g., "double".
    :return: A string containing the generated C code for the function.
    """
    includes: List[str] = []
    prefunc = ""
    rfm = refmetric.reference_metric[CoordSystem]
    for i in range(3):
        # Do not generate FD1_arbitrary_upwind_xj_dirn() if the symbolic expression for dxj/dr == 0!
        if not check_zero(rfm.Jac_dUrfm_dDSphUD[i][0]):
            prefunc += setup_Cfunction_FD1_arbitrary_upwind(
                dirn=i,
                radiation_BC_fd_order=radiation_BC_fd_order,
                fp_type=fp_type,
            )
    prefunc += setup_Cfunction_r_and_partial_xi_partial_r_derivs(
        CoordSystem=CoordSystem,
        fp_type=fp_type,
    )
    prefunc += setup_Cfunction_compute_partial_r_f(
        CoordSystem=CoordSystem, radiation_BC_fd_order=radiation_BC_fd_order
    )
    desc = r"""*** Apply radiation BCs to all outer boundaries. ***
"""
    cfunc_type = "static inline REAL"
    name = "radiation_bcs"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    const bc_struct *restrict bcstruct,REAL *restrict xx[3],
    const REAL *restrict gfs, REAL *restrict gfs_rhss,
    const int which_gf, const REAL gf_wavespeed, const REAL gf_f_infinity,
    const int dest_i0,const int dest_i1,const int dest_i2,
    const short FACEi0,const short FACEi1,const short FACEi2"""
    body = r"""// Nearest "interior" neighbor of this gridpoint, based on current face
const int dest_i0_int=dest_i0+1*FACEi0, dest_i1_int=dest_i1+1*FACEi1, dest_i2_int=dest_i2+1*FACEi2;
REAL r, partial_x0_partial_r,partial_x1_partial_r,partial_x2_partial_r;
REAL r_int, partial_x0_partial_r_int,partial_x1_partial_r_int,partial_x2_partial_r_int;
r_and_partial_xi_partial_r_derivs(commondata, params,xx[0][dest_i0],xx[1][dest_i1],xx[2][dest_i2],
                                  &r, &partial_x0_partial_r, &partial_x1_partial_r,  &partial_x2_partial_r);
r_and_partial_xi_partial_r_derivs(commondata, params, xx[0][dest_i0_int], xx[1][dest_i1_int], xx[2][dest_i2_int],
                                  &r_int, &partial_x0_partial_r_int, &partial_x1_partial_r_int, &partial_x2_partial_r_int);
const REAL partial_r_f     = compute_partial_r_f(commondata, params,xx,gfs, which_gf,dest_i0,    dest_i1,    dest_i2,
                                                 FACEi0,FACEi1,FACEi2,
                                                 partial_x0_partial_r    ,partial_x1_partial_r    ,partial_x2_partial_r);
const REAL partial_r_f_int = compute_partial_r_f(commondata, params,xx,gfs, which_gf,dest_i0_int,dest_i1_int,dest_i2_int,
                                                 FACEi0,FACEi1,FACEi2,
                                                 partial_x0_partial_r_int,partial_x1_partial_r_int,partial_x2_partial_r_int);

const int idx3 = IDX3(dest_i0,dest_i1,dest_i2);
const int idx3_int = IDX3(dest_i0_int,dest_i1_int,dest_i2_int);

const REAL partial_t_f_int = gfs_rhss[IDX4pt(which_gf, idx3_int)];

const REAL c = gf_wavespeed;
const REAL f_infinity = gf_f_infinity;
const REAL f     = gfs[IDX4pt(which_gf, idx3)];
const REAL f_int = gfs[IDX4pt(which_gf, idx3_int)];
const REAL partial_t_f_int_outgoing_wave = -c * (partial_r_f_int + (f_int - f_infinity) / r_int);

const REAL k = r_int*r_int*r_int * (partial_t_f_int - partial_t_f_int_outgoing_wave);

const REAL rinv = 1.0 / r;
const REAL partial_t_f_outgoing_wave = -c * (partial_r_f + (f - f_infinity) * rinv);

return partial_t_f_outgoing_wave + k * rinv*rinv*rinv;
"""

    cf = cfc.CFunction(
        subdirectory=CoordSystem,
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cf.full_function


# apply_bcs_outerradiation_and_inner():
#   Apply radiation BCs at outer boundary points, and
#   inner boundary conditions at inner boundary points.
class register_CFunction_apply_bcs_outerradiation_and_inner(
    base_cbc_classes.base_register_CFunction_apply_bcs_outerradiation_and_inner
):
    def __init__(
        self,
        CoordSystem: str,
        radiation_BC_fd_order: int = 2,
        fp_type: str = "double",
    ) -> None:
        """
        Register a C function to apply boundary conditions to both pure outer and inner boundary points.

        :param CoordSystem: The coordinate system to use.
        :param radiation_BC_fd_order: Finite differencing order for the radiation boundary conditions. Default is 2.
        :param fp_type: Floating point type, e.g., "double".
        """
        super().__init__(
            CoordSystem,
            radiation_BC_fd_order=radiation_BC_fd_order,
            fp_type=fp_type,
        )
        self.body = r"""
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

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
  // Spawn N OpenMP threads, either across all cores, or according to e.g., taskset.
#pragma omp parallel
  {
    for(int which_gz=0;which_gz<NGHOSTS;which_gz++) for(int dirn=0;dirn<3;dirn++) {
        // This option results in about 1.6% slower runtime for SW curvilinear at 64x24x24 on 8-core Ryzen 9 4900HS
        //#pragma omp for collapse(2)
        //for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++) {
        //  {
        // Don't spawn a thread if there are no boundary points to fill; results in a nice little speedup.
        if(bc_info->num_pure_outer_boundary_points[which_gz][dirn] > 0) {
#pragma omp for  // threads have been spawned; here we distribute across them
          for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++) {
            const short i0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0;
            const short i1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1;
            const short i2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2;
            const short FACEX0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0;
            const short FACEX1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1;
            const short FACEX2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2;
            const int idx3 = IDX3(i0,i1,i2);
            for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
              // *** Apply radiation BCs to all outer boundary points. ***
              rhs_gfs[IDX4pt(which_gf, idx3)] = radiation_bcs(commondata, params, bcstruct, xx, gfs, rhs_gfs, which_gf,
                                                               custom_wavespeed[which_gf], custom_f_infinity[which_gf],
                                                               i0,i1,i2, FACEX0,FACEX1,FACEX2);
            }
          }
        }
      }
  }

  ///////////////////////////////////////////////////////
  // STEP 2 of 2: Apply BCs to inner boundary points.
  //              These map to either the grid interior
  //              ("pure inner") or to pure outer boundary
  //              points ("inner maps to outer"). Those
  //              that map to outer require that outer be
  //              populated first; hence this being
  //              STEP 2 OF 2.
  apply_bcs_inner_only(commondata, params, bcstruct, rhs_gfs); // <- apply inner BCs to RHS gfs only
"""
        cfc.register_CFunction(
            includes=self.includes,
            prefunc=self.prefunc,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            CoordSystem_for_wrapper_func=self.CoordSystem,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=True,
            body=self.body,
        )


class CurviBoundaryConditions_register_C_functions(
    base_cbc_classes.base_CurviBoundaryConditions_register_C_functions
):

    def __init__(
        self,
        list_of_CoordSystems: List[str],
        radiation_BC_fd_order: int = 2,
        set_parity_on_aux: bool = False,
        set_parity_on_auxevol: bool = False,
        fp_type: str = "double",
    ) -> None:
        """
        Register various C functions responsible for handling boundary conditions.

        :param list_of_CoordSystems: List of coordinate systems to use.
        :param radiation_BC_fd_order: Finite differencing order for the radiation boundary conditions. Default is 2.
        :param set_parity_on_aux: If True, set parity on auxiliary grid functions.
        :param set_parity_on_auxevol: If True, set parity on auxiliary evolution grid functions.
        :param fp_type: Floating point type, e.g., "double".
        """
        super().__init__(
            list_of_CoordSystems,
            radiation_BC_fd_order=radiation_BC_fd_order,
            set_parity_on_aux=set_parity_on_aux,
            set_parity_on_auxevol=set_parity_on_auxevol,
            fp_type=fp_type,
        )
        # self.post_register_BHAH_header(self)
        for CoordSystem in self.list_of_CoordSystems:
            # Register C function to set up the boundary condition struct.
            register_CFunction_bcstruct_set_up(
                CoordSystem=CoordSystem, fp_type=self.fp_type
            )

            # Register C function to apply boundary conditions to both pure outer and inner boundary points.
            register_CFunction_apply_bcs_outerradiation_and_inner(
                CoordSystem=CoordSystem,
                radiation_BC_fd_order=self.radiation_BC_fd_order,
                fp_type=self.fp_type,
            )

        # Register C function to apply boundary conditions to inner-only boundary points.
        register_CFunction_apply_bcs_inner_only()

        # Register C function to apply boundary conditions to outer-extrapolated and inner boundary points.
        register_CFunction_apply_bcs_outerextrap_and_inner()

        griddata_commondata.register_griddata_commondata(
            __name__,
            "bc_struct bcstruct",
            "all data needed to perform boundary conditions in curvilinear coordinates",
        )

        BHaH_defines_h.register_BHaH_defines(__name__, self.CBC_BHd_str)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
