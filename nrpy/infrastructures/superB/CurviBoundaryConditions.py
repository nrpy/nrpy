"""

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
         Nishita Jadoo
         njadoo **at** uidaho **dot* edu
"""

# Step P1: Import needed NRPy+ core modules:
from typing import List, Dict, Tuple, Union
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
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


def register_CFunction_bcstruct_set_up(CoordSystem: str) -> None:
    """
    Register C function for setting up bcstruct_chare.

    :param CoordSystem: The coordinate system for which to set up boundary conditions.
    """
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = r"""Setup bcstruct_chare from bcstruct"""
    c_type = "void"
    name = "bcstruct_chare_set_up"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const params_struct *restrict params_chare, REAL *restrict xx[3], const bc_struct *restrict bcstruct, bc_struct *restrict bcstruct_chare, const int chare_index[3]"
    body = r"""
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

  ////////////////////////////////////////
  // STEP 1: SET UP INNER BOUNDARY STRUCTS
  {
    // First count the number of inner points that lie inside the chare grid.
    int num_inner_chare = 0;
    for (int pt = 0; pt < bcstruct->bc_info.num_inner_boundary_points; pt++) {
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt;
      if (charecommstruct->globalidx3pt_to_chareidx3[dstpt] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])){
        num_inner_chare++;
      }
    }
    // Store num_inner to bc_info:
    bcstruct_chare->bc_info.num_inner_boundary_points = num_inner_chare;

    // Next allocate memory for inner_boundary_points:
    bcstruct_chare->inner_bc_array = (innerpt_bc_struct *restrict)malloc(sizeof(innerpt_bc_struct) * num_inner_chare);
  }

  // Then set inner_bc_array:
  {
    int  which_inner_chare = 0;
    for (int pt = 0; pt < bcstruct->bc_info.num_inner_boundary_points; pt++) {
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt;
      const int srcpt = bcstruct->inner_bc_array[pt].srcpt;
      // if dstpt is part of local grid
      if (charecommstruct->globalidx3pt_to_chareidx3[dstpt] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])) {
        bcstruct_chare->inner_bc_array[which_inner_chare].dstpt = charecommstruct->globalidx3pt_to_localidx3pt[dstpt];// store the local idx3 of dstpt
        // Now check srcpt is in local grid
        // If yes, save local index of srcpt
        if (charecommstruct->globalidx3pt_to_chareidx3[srcpt] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])) {
          bcstruct_chare->inner_bc_array[which_inner_chare].srcpt = charecommstruct->global_to_local_idx3[srcpt];
        } else {
          printf("Error: dst pt is in chare's grid but not src pt\n");
        }
        // Also copy parity
        for (int i = 0; i < 10; ++i) {
          bcstruct_chare->inner_bc_array[which_inner_chare].parity[i] = bcstruct->inner_bc_array[pt].parity[i];
        }
        which_inner_chare++;
      }
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
          const int globalidx3 =  IDX3GENERAL(i0, i1, i2, Nxx0, Nxx1);
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
          const int globalidx3 =  IDX3GENERAL(i0, i1, i2, Nxx0, Nxx1);
          if (charecommstruct->globalidx3pt_to_chareidx3[globalidx3] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])){
            // convert i0, etc to local i0
            const int localidx3 = charecommstruct->globalidx3pt_to_localidx3pt[globalidx3];
            int locali0, locali1, locali2;
            REVERSE_IDX3GENERAL(localidx3, Nxx0chare, Nxx1chare);
            bcstruct_chare->pure_outer_bc_array[dirn + (3 * which_gz)][which_idx2d_chare].i0 = locali0;
            bcstruct_chare->pure_outer_bc_array[dirn + (3 * which_gz)][which_idx2d_chare].i1 = locali1;
            bcstruct_chare->pure_outer_bc_array[dirn + (3 * which_gz)][which_idx2d_chare].i2 = locali2;
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
        c_type=c_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


def CurviBoundaryConditions_register_C_functions(
    list_of_CoordSystems: List[str],
    radiation_BC_fd_order: int = 2,
    set_parity_on_aux: bool = False,
    set_parity_on_auxevol: bool = False,
) -> None:
    """
    Register various C functions responsible for handling boundary conditions.

    :param list_of_CoordSystems: List of coordinate systems to use.
    :param radiation_BC_fd_order: Finite differencing order for the radiation boundary conditions. Default is 2.
    :param set_parity_on_aux: If True, set parity on auxiliary grid functions.
    :param set_parity_on_auxevol: If True, set parity on auxiliary evolution grid functions.
    :return: None
    """
    for CoordSystem in list_of_CoordSystems:
        # Register C function to set up the boundary condition struct.
        register_CFunction_bcstruct_set_up(CoordSystem=CoordSystem)

    # Register bcstruct's contribution to griddata_struct:
    griddata_commondata.register_griddata_commondata(
        __name__,
        "bc_struct bcstruct_chare",
        "all data needed to perform boundary conditions in curvilinear coordinates in a chare",
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
