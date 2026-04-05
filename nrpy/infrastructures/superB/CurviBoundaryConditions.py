"""
C function for setting up the boundary condition struct for the local chare grid using the boundary condition struct for the global grid for Curvilinear boundary conditions.

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
         Nishita Jadoo
         njadoo **at** uidaho **dot* edu
"""

# Step P1: Import needed NRPy core modules:
from typing import Set

import nrpy.c_function as cfc
from nrpy.infrastructures import BHaH


def register_CFunction_apply_bcs_inner_only_nonlocal() -> None:
    """Register C function for filling inner boundary points on the computational grid, as prescribed by bcstruct."""
    includes = ["BHaH_defines.h"]
    desc = r"""
Apply BCs to nonlocal inner boundary points only,
using data stored in bcstruct->inner_bc_array.
These structs are set in bcstruct_set_up().
Inner boundary points map to either the grid
interior ("pure inner") or to pure outer
boundary points ("inner maps to outer").
"""
    cfunc_type = "void"
    name = "apply_bcs_inner_only_nonlocal"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct, const nonlocalinnerbc_struct *restrict nonlocalinnerbcstruct, const int NUM_GFS, REAL *restrict gfs, const int *gfs_to_sync, const int8_t* gf_parity_types, REAL **restrict tmpBuffer_innerbc_receiv"
    body = r"""
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

  // Unpack nonlocalinnerbcstruct
  const int tot_num_src_chares = nonlocalinnerbcstruct->tot_num_src_chares;
  const int *restrict num_srcpts_each_chare = nonlocalinnerbcstruct->num_srcpts_each_chare;
  int **restrict map_srcchare_and_srcpt_id_to_linear_id = nonlocalinnerbcstruct->map_srcchare_and_srcpt_id_to_linear_id;

  for (int which_chare = 0; which_chare < tot_num_src_chares; which_chare++) {
    const REAL *restrict tmpBuffer = tmpBuffer_innerbc_receiv[which_chare];
    for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
      for (int which_srcpt = 0; which_srcpt < num_srcpts_each_chare[which_chare]; which_srcpt++) {
        const int linear_id = map_srcchare_and_srcpt_id_to_linear_id[which_chare][which_srcpt];
        const int dstpt = bcstruct->inner_bc_array_nonlocal[linear_id].dstpt;
        const int idx2 = IDX2NONLOCALINNERBC(which_gf, which_srcpt, num_srcpts_each_chare[which_chare]);
        gfs[IDX4pt(gfs_to_sync[which_gf], dstpt)] = bcstruct->inner_bc_array_nonlocal[linear_id].parity[gf_parity_types[gfs_to_sync[which_gf]]] * tmpBuffer[idx2];
      }
    }
  }
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )


def register_CFunction_bcstruct_chare_set_up(CoordSystem: str) -> None:
    """
    Register C function for setting up bc_struct for the local chare grid.

    :param CoordSystem: The coordinate system for which to set up boundary conditions.
    """
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = r"""Setup bcstruct_chare from bcstruct"""
    cfunc_type = "void"
    name = "bcstruct_chare_set_up"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const params_struct *restrict params_chare, const charecomm_struct *restrict charecommstruct, REAL *restrict xx[3], const bc_struct *restrict bcstruct, bc_struct *restrict bcstruct_chare, nonlocalinnerbc_struct *restrict nonlocalinnerbcstruct, const int chare_index[3]"
    body = r"""
  const int Nchare0 = commondata->Nchare0;
  const int Nchare1 = commondata->Nchare1;
  const int Nchare2 = commondata->Nchare2;
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  const int Nxx0 = params->Nxx0;
  const int Nxx1 = params->Nxx1;
  const int Nxx2 = params->Nxx2;
  const int Nxx_plus_2NGHOSTS0chare = params_chare->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1chare = params_chare->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2chare = params_chare->Nxx_plus_2NGHOSTS2;
  const int Nxx0chare = params_chare->Nxx0;
  const int Nxx1chare = params_chare->Nxx1;
  const int Nxx2chare = params_chare->Nxx2;
  const int tot_num_chares = Nchare0*Nchare1*Nchare2;
  const int idx3_this_chare = IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2]);

  // Unpack globalidx3pt_to_chareidx and globalidx3pt_to_localidx3pt from charecommstruct
  const int *restrict globalidx3pt_to_chareidx3 = charecommstruct->globalidx3pt_to_chareidx3;
  const int *restrict globalidx3pt_to_localidx3pt = charecommstruct->globalidx3pt_to_localidx3pt;

  ////////////////////////////////////////
  // STEP 1: SET UP INNER BOUNDARY STRUCTS
  {
    // Count inner bc pts that are on this chare's grid with a local or nonlocal scrpt
    int num_inner_chare = 0;
    int num_inner_chare_nonlocal = 0;
    for (int pt = 0; pt < bcstruct->bc_info.num_inner_boundary_points; pt++) {
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt;
      const int srcpt = bcstruct->inner_bc_array[pt].srcpt;
      if (globalidx3pt_to_chareidx3[dstpt] == idx3_this_chare) {
        if (globalidx3pt_to_chareidx3[srcpt] == idx3_this_chare) {
          num_inner_chare++;
        } else {
          num_inner_chare_nonlocal++;
        }
      }
    }
    // Store num_inner to bc_info:
    bcstruct_chare->bc_info.num_inner_boundary_points = num_inner_chare;
    bcstruct_chare->bc_info.num_inner_boundary_points_nonlocal = num_inner_chare_nonlocal;

    // Next allocate memory for inner_boundary_points:
    bcstruct_chare->inner_bc_array = (innerpt_bc_struct *restrict)malloc(sizeof(innerpt_bc_struct) * num_inner_chare);
    bcstruct_chare->inner_bc_array_nonlocal = (innerpt_bc_struct *restrict)malloc(sizeof(innerpt_bc_struct) * num_inner_chare_nonlocal);

  // Then set inner_bc_array and count nonlocal srcpts and dstpts
    int which_inner_chare = 0;
    int which_inner_chare_nonlocal = 0;
    int count_nonlocalinnerbc_srcpt_eachchare[tot_num_chares]= {0};
    int count_nonlocalinnerbc_dstpt_eachchare[tot_num_chares]= {0};
    for (int pt = 0; pt < bcstruct->bc_info.num_inner_boundary_points; pt++) {
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt;
      const int srcpt = bcstruct->inner_bc_array[pt].srcpt;
      // if dstpt is in local grid
      if (globalidx3pt_to_chareidx3[dstpt] == idx3_this_chare) {
        // if srcpt is in local grid
        if (globalidx3pt_to_chareidx3[srcpt] == idx3_this_chare) {
          bcstruct_chare->inner_bc_array[which_inner_chare].dstpt = globalidx3pt_to_localidx3pt[dstpt]; // store the local idx3 of dstpt
          bcstruct_chare->inner_bc_array[which_inner_chare].srcpt = globalidx3pt_to_localidx3pt[srcpt]; // store the local idx3 of srcpt
          // Also copy parity
          for (int i = 0; i < 10; ++i) {
            bcstruct_chare->inner_bc_array[which_inner_chare].parity[i] = bcstruct->inner_bc_array[pt].parity[i];
          }
          which_inner_chare++;
        // if srcpt is not in local grid
        } else {
          count_nonlocalinnerbc_srcpt_eachchare[globalidx3pt_to_chareidx3[srcpt]]++;
        }
      // if dstpt is not in local grid
      } else {
        // if srcpt is in local grid
        if (globalidx3pt_to_chareidx3[srcpt] == idx3_this_chare) {
          count_nonlocalinnerbc_dstpt_eachchare[globalidx3pt_to_chareidx3[dstpt]]++;
        }
      }
    }

    // Next set up nonlocalinnerbcstruct

    // dstpt is in local grid, srcpt is not in local grid, count number of src chares
    int count_src_chares = 0;
    for (int chareidx3 = 0; chareidx3 < tot_num_chares; chareidx3++) {
      if (count_nonlocalinnerbc_srcpt_eachchare[chareidx3] > 0) {
        count_src_chares++;
      }
    }
    // Set value
    nonlocalinnerbcstruct->tot_num_src_chares = count_src_chares;
    // Allocate memory
    nonlocalinnerbcstruct->idx3_of_src_chares = (int *restrict)malloc(sizeof(int) * count_src_chares);
    nonlocalinnerbcstruct->num_srcpts_each_chare = (int *restrict)malloc(sizeof(int) * count_src_chares);
    nonlocalinnerbcstruct->idx3chare_to_src_chare_id = (int *restrict)malloc(sizeof(int) * tot_num_chares);
    // Initialize idx3chare_to_src_chare_id to default value of -1
    for (int i = 0; i < tot_num_chares; i++) {
      nonlocalinnerbcstruct->idx3chare_to_src_chare_id[i] = -1;
    }
    // Set arrays in nonlocalinnerbcstruct
    count_src_chares = 0;
    for (int chareidx3 = 0; chareidx3 < tot_num_chares; chareidx3++) {
      if (count_nonlocalinnerbc_srcpt_eachchare[chareidx3] > 0) {
        nonlocalinnerbcstruct->idx3_of_src_chares[count_src_chares] = chareidx3;
        nonlocalinnerbcstruct->num_srcpts_each_chare[count_src_chares] = count_nonlocalinnerbc_srcpt_eachchare[chareidx3];
        nonlocalinnerbcstruct->idx3chare_to_src_chare_id[chareidx3] = count_src_chares;
        count_src_chares++;
      }
    }
    // Map src chare id and src pt id to a single id
    nonlocalinnerbcstruct-> map_srcchare_and_srcpt_id_to_linear_id = (int **)malloc(nonlocalinnerbcstruct->tot_num_src_chares * sizeof(int *));
    int current_index = 0;
    for (int src_chare_id = 0; src_chare_id < nonlocalinnerbcstruct->tot_num_src_chares; src_chare_id++) {
        nonlocalinnerbcstruct-> map_srcchare_and_srcpt_id_to_linear_id[src_chare_id] = (int *restrict)malloc(sizeof(int) * nonlocalinnerbcstruct->num_srcpts_each_chare[src_chare_id]);
        for (int srcpt_id = 0; srcpt_id < nonlocalinnerbcstruct->num_srcpts_each_chare[src_chare_id]; srcpt_id++) {
            nonlocalinnerbcstruct-> map_srcchare_and_srcpt_id_to_linear_id[src_chare_id][srcpt_id] = current_index;
            current_index++;
        }
    }
    // Allocate memory for globalidx3_srcpts
    nonlocalinnerbcstruct->globalidx3_srcpts = (int **)malloc(nonlocalinnerbcstruct->tot_num_src_chares * sizeof(int *));
    for (int src_chare_id = 0; src_chare_id < nonlocalinnerbcstruct->tot_num_src_chares; src_chare_id++) {
      nonlocalinnerbcstruct->globalidx3_srcpts[src_chare_id] = (int *restrict)malloc(sizeof(int) * nonlocalinnerbcstruct->num_srcpts_each_chare[src_chare_id]);
    }
    // Set globalidx3_srcpts and bcstruct_chare->inner_bc_array_nonlocal
    int count_nonlocalinnerbc_srcpt_each_src_chare[nonlocalinnerbcstruct->tot_num_src_chares] = {0};
    for (int pt = 0; pt < bcstruct->bc_info.num_inner_boundary_points; pt++) {
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt;
      const int srcpt = bcstruct->inner_bc_array[pt].srcpt;
      // if dstpt is in local grid
      if (globalidx3pt_to_chareidx3[dstpt] == idx3_this_chare) {
        int chareidx3 = globalidx3pt_to_chareidx3[srcpt];
        // if srcpt is not in local grid
        if ( chareidx3 != idx3_this_chare) {
          const int src_chare_id = nonlocalinnerbcstruct->idx3chare_to_src_chare_id[chareidx3];
          nonlocalinnerbcstruct->globalidx3_srcpts[src_chare_id][count_nonlocalinnerbc_srcpt_each_src_chare[src_chare_id]] =  srcpt;
          const int linear_id = nonlocalinnerbcstruct-> map_srcchare_and_srcpt_id_to_linear_id[src_chare_id][count_nonlocalinnerbc_srcpt_each_src_chare[src_chare_id]];
          bcstruct_chare->inner_bc_array_nonlocal[linear_id].dstpt = globalidx3pt_to_localidx3pt[dstpt]; // store the local idx3 of dstpt
          // Also copy parity
          for (int i = 0; i < 10; ++i) {
            bcstruct_chare->inner_bc_array_nonlocal[linear_id].parity[i] = bcstruct->inner_bc_array[pt].parity[i];
          }
          count_nonlocalinnerbc_srcpt_each_src_chare[src_chare_id]++;
        }
      }
    }
    // Next: dstpt is not in local grid, srcpt is in local grid, count number of dst chares
    int count_dst_chares = 0;
    for (int chareidx3 = 0; chareidx3 < tot_num_chares; chareidx3++) {
      if (count_nonlocalinnerbc_dstpt_eachchare[chareidx3] > 0) {
        count_dst_chares++;
      }
    }
    // Set value
    nonlocalinnerbcstruct->tot_num_dst_chares = count_dst_chares;
    // Allocate memory
    nonlocalinnerbcstruct->idx3chare_to_dst_chare_id = (int *restrict)malloc(sizeof(int) * tot_num_chares);
    nonlocalinnerbcstruct->num_srcpts_tosend_each_chare = (int *restrict)malloc(sizeof(int) * count_dst_chares);
    nonlocalinnerbcstruct->idx3_of_dst_chares = (int *restrict)malloc(sizeof(int) * count_dst_chares);
    // Initialize idx3chare_to_dst_chare_id to default value of -1
    for (int i = 0; i < tot_num_chares; i++) {
      nonlocalinnerbcstruct->idx3chare_to_dst_chare_id[i] = -1;
    }
    // Set arrays in nonlocalinnerbcstruct
    count_dst_chares = 0;
    for (int chareidx3 = 0; chareidx3 < tot_num_chares; chareidx3++) {
      if (count_nonlocalinnerbc_dstpt_eachchare[chareidx3] > 0) {
        nonlocalinnerbcstruct->idx3_of_dst_chares[count_dst_chares] = chareidx3;
        nonlocalinnerbcstruct->num_srcpts_tosend_each_chare[count_dst_chares] = count_nonlocalinnerbc_dstpt_eachchare[chareidx3];
        nonlocalinnerbcstruct->idx3chare_to_dst_chare_id[chareidx3] = count_dst_chares;
        count_dst_chares++;
      }
    }
    // Allocate memory for globalidx3_srcpts_tosend
    nonlocalinnerbcstruct->globalidx3_srcpts_tosend = (int **)malloc(nonlocalinnerbcstruct->tot_num_dst_chares * sizeof(int *));
    for (int dst_chare_id = 0; dst_chare_id < nonlocalinnerbcstruct->tot_num_dst_chares; dst_chare_id++) {
      nonlocalinnerbcstruct->globalidx3_srcpts_tosend[dst_chare_id] = (int *restrict)malloc(sizeof(int) * nonlocalinnerbcstruct->num_srcpts_tosend_each_chare[dst_chare_id]);
    }
  }

  ////////////////////////////////////////
  // STEP 2: SET UP OUTER BOUNDARY STRUCTS
  for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
    for (int dirn = 0; dirn < 3; dirn++) {
      int num_pure_outer_boundary_points_chare = 0;
      if (bcstruct->bc_info.num_pure_outer_boundary_points[which_gz][dirn] > 0) {
        for (int idx2d = 0; idx2d < bcstruct->bc_info.num_pure_outer_boundary_points[which_gz][dirn]; idx2d++) {
          const short i0 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i0;
          const short i1 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i1;
          const short i2 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i2;
          const int globalidx3 =  IDX3GENERAL(i0, i1, i2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1);
          if (charecommstruct->globalidx3pt_to_chareidx3[globalidx3] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])){
            num_pure_outer_boundary_points_chare++;
          }
        }
      }
      bcstruct_chare->bc_info.num_pure_outer_boundary_points[which_gz][dirn] = num_pure_outer_boundary_points_chare;
      bcstruct_chare->pure_outer_bc_array[dirn + (3 * which_gz)] = (outerpt_bc_struct *restrict)malloc(num_pure_outer_boundary_points_chare * sizeof(outerpt_bc_struct));
    }
  }

  for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
    for (int dirn = 0; dirn < 3; dirn++) {
      int which_idx2d_chare = 0;
      if (bcstruct->bc_info.num_pure_outer_boundary_points[which_gz][dirn] > 0) {
        for (int idx2d = 0; idx2d < bcstruct->bc_info.num_pure_outer_boundary_points[which_gz][dirn]; idx2d++) {
          const short i0 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i0;
          const short i1 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i1;
          const short i2 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i2;
          const int globalidx3 =  IDX3GENERAL(i0, i1, i2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1);
          if (charecommstruct->globalidx3pt_to_chareidx3[globalidx3] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])){
            bcstruct_chare->pure_outer_bc_array[dirn + (3 * which_gz)][which_idx2d_chare].i0 = MAP_GLOBAL_TO_LOCAL_IDX0(chare_index[0], i0, Nxx0chare);
            bcstruct_chare->pure_outer_bc_array[dirn + (3 * which_gz)][which_idx2d_chare].i1 = MAP_GLOBAL_TO_LOCAL_IDX1(chare_index[1], i1, Nxx1chare);
            bcstruct_chare->pure_outer_bc_array[dirn + (3 * which_gz)][which_idx2d_chare].i2 = MAP_GLOBAL_TO_LOCAL_IDX2(chare_index[2], i2, Nxx2chare);
            const short FACEX0 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX0;
            const short FACEX1 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX1;
            const short FACEX2 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX2;
            bcstruct_chare->pure_outer_bc_array[dirn + (3 * which_gz)][which_idx2d_chare].FACEX0 = FACEX0;
            bcstruct_chare->pure_outer_bc_array[dirn + (3 * which_gz)][which_idx2d_chare].FACEX1 = FACEX1;
            bcstruct_chare->pure_outer_bc_array[dirn + (3 * which_gz)][which_idx2d_chare].FACEX2 = FACEX2;
            which_idx2d_chare++;
          }
        }
      }
    }
  }
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


def CurviBoundaryConditions_register_C_functions(
    set_of_CoordSystems: Set[str],
    radiation_BC_fd_order: int = 2,
    set_parity_on_aux: bool = False,
) -> None:
    """
    Register various C functions responsible for handling boundary conditions.

    :param set_of_CoordSystems: Set of coordinate systems to use.
    :param radiation_BC_fd_order: Finite differencing order for the radiation boundary conditions. Default is 2.
    :param set_parity_on_aux: If True, set parity on auxiliary grid functions.
    """
    for CoordSystem in set_of_CoordSystems:
        # Register C function to set up the boundary condition struct.
        BHaH.CurviBoundaryConditions.bcstruct_set_up.register_CFunction_bcstruct_set_up(
            CoordSystem=CoordSystem
        )

        # Register C function to apply boundary conditions to both pure outer and inner boundary points.
        BHaH.CurviBoundaryConditions.apply_bcs_outerradiation_and_inner.register_CFunction_apply_bcs_outerradiation_and_inner(
            CoordSystem=CoordSystem,
            radiation_BC_fd_order=radiation_BC_fd_order,
        )

    # Register C function to apply boundary conditions to inner-only boundary points.
    BHaH.CurviBoundaryConditions.apply_bcs_inner_only.register_CFunction_apply_bcs_inner_only()

    # Register C function to apply boundary conditions to inner-only boundary points for specific gfs.
    BHaH.CurviBoundaryConditions.apply_bcs_inner_only_specific_gfs.register_CFunction_apply_bcs_inner_only_specific_gfs()

    # Register C function to apply boundary conditions to nonlocal inner-only boundary points.
    register_CFunction_apply_bcs_inner_only_nonlocal()

    # Register C function to apply boundary conditions to outer-extrapolated and inner boundary points.
    BHaH.CurviBoundaryConditions.apply_bcs_outerextrap_and_inner.register_CFunction_apply_bcs_outerextrap_and_inner()

    # Register C function to apply boundary conditions to outer-extrapolated and inner boundary points for specific gfs.
    # Was once needed for psi4, now unused.
    # CurviBoundaryConditions.apply_bcs_outerextrap_and_inner.register_CFunction_apply_bcs_outerextrap_and_inner_specific_gfs()

    # Register bcstruct's contribution to griddata_struct:
    BHaH.griddata_commondata.register_griddata_commondata(
        __name__,
        "bc_struct bcstruct",
        "all data needed to perform boundary conditions in curvilinear coordinates",
    )

    # Register bcstruct's contribution to BHaH_defines.h:
    CBC_BHd_str = r"""
// NRPy Curvilinear Boundary Conditions: Core data structures
// Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

typedef struct __innerpt_bc_struct__ {
  int dstpt;  // dstpt is the 3D grid index IDX3S(i0,i1,i2) of the inner boundary point (i0,i1,i2)
  int srcpt;  // srcpt is the 3D grid index (a la IDX3S) to which the inner boundary point maps
  int8_t parity[10];  // parity[10] is a calculation of dot products for the 10 independent parity types
} innerpt_bc_struct;

typedef struct __outerpt_bc_struct__ {
  short i0,i1,i2;  // the outer boundary point grid index (i0,i1,i2), on the 3D grid
  int8_t FACEX0,FACEX1,FACEX2;  // 1-byte integers that store
  //                               FACEX0,FACEX1,FACEX2 = +1, 0, 0 if on the i0=i0min face,
  //                               FACEX0,FACEX1,FACEX2 = -1, 0, 0 if on the i0=i0max face,
  //                               FACEX0,FACEX1,FACEX2 =  0,+1, 0 if on the i1=i2min face,
  //                               FACEX0,FACEX1,FACEX2 =  0,-1, 0 if on the i1=i1max face,
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,+1 if on the i2=i2min face, or
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,-1 if on the i2=i2max face,
} outerpt_bc_struct;

typedef struct __bc_info_struct__ {
  int num_inner_boundary_points;  // stores total number of inner boundary points
  int num_inner_boundary_points_nonlocal;  // stores total number of inner boundary points that lie on another's grid
  int num_pure_outer_boundary_points[NGHOSTS][3];  // stores number of outer boundary points on each
  //                                                  ghostzone level and direction (update min and
  //                                                  max faces simultaneously on multiple cores)
  int bc_loop_bounds[NGHOSTS][6][6];  // stores outer boundary loop bounds. Unused after bcstruct_set_up()
} bc_info_struct;

typedef struct __bc_struct__ {
  innerpt_bc_struct *inner_bc_array;  // information needed for updating each inner boundary point
  innerpt_bc_struct *inner_bc_array_nonlocal;  // information needed for updating each nonlocal inner boundary point
  outerpt_bc_struct *pure_outer_bc_array[NGHOSTS*3]; // information needed for updating each outer
  //                                                             boundary point
  bc_info_struct bc_info;  // stores number of inner and outer boundary points, needed for setting loop
  //                          bounds and parallelizing over as many boundary points as possible.
} bc_struct;
"""
    # inter-chare communication assumes auxevol parity types are set
    CBC_BHd_str += BHaH.CurviBoundaryConditions.BHaH_defines.BHaH_defines_set_gridfunction_defines_with_parity_types(
        set_parity_on_aux=set_parity_on_aux,
        set_parity_on_auxevol=True,
        verbose=True,
    )
    # ~ BHaH_defines_h.register_BHaH_defines(__name__, CBC_BHd_str)
    BHaH.BHaH_defines_h.register_BHaH_defines(
        "nrpy.infrastructures.BHaH.CurviBoundaryConditions.BHaH_defines",
        CBC_BHd_str,
    )

    for CoordSystem in set_of_CoordSystems:
        # Register C function to set up the boundary condition struct for local chare grid.
        register_CFunction_bcstruct_chare_set_up(CoordSystem=CoordSystem)

    # Register temporary buffers for face data communication to griddata_struct:
    BHaH.griddata_commondata.register_griddata_commondata(
        __name__,
        "nonlocalinnerbc_struct nonlocalinnerbcstruct",
        "for communication of non-local inner boundary data across chares",
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
