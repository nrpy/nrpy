"""
Register function and CodeParameters for managing and interpolating external input grid data on a 3D spherical source grid.

Specifically,
- Define and register external input grid parameters.
- Register gridfunctions for the external input data (gamma_{ij} and K_{ij} tensors).
- Set up and initialize the external input grid, including memory allocation and boundary conditions.
- Perform coordinate transformations from Cartesian to spherical basis and apply necessary rescaling.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import sys
from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Tuple, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.BHaHAHA.bcstruct_set_up as locCBC
import nrpy.params as par
from nrpy.infrastructures.BHaH import BHaH_defines_h, griddata_commondata


def register_CFunction_numgrid__interp_src_set_up() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for reading original source metric data.

    :return: None if in registration phase, else the updated NRPy environment.

    DocTests:
    >>> env = register_CFunction_numgrid__interp_src_set_up()
    Grid function "src_aDD00" with rank 2 has parity type 4.
    Grid function "src_aDD01" with rank 2 has parity type 5.
    Grid function "src_aDD02" with rank 2 has parity type 6.
    Grid function "src_aDD11" with rank 2 has parity type 7.
    Grid function "src_aDD12" with rank 2 has parity type 8.
    Grid function "src_aDD22" with rank 2 has parity type 9.
    Grid function "src_hDD00" with rank 2 has parity type 4.
    Grid function "src_hDD01" with rank 2 has parity type 5.
    Grid function "src_hDD02" with rank 2 has parity type 6.
    Grid function "src_hDD11" with rank 2 has parity type 7.
    Grid function "src_hDD12" with rank 2 has parity type 8.
    Grid function "src_hDD22" with rank 2 has parity type 9.
    Grid function "src_partial_D_hDD000" with rank 3 has parity type 10.
    Grid function "src_partial_D_hDD001" with rank 3 has parity type 11.
    Grid function "src_partial_D_hDD002" with rank 3 has parity type 12.
    Grid function "src_partial_D_hDD011" with rank 3 has parity type 13.
    Grid function "src_partial_D_hDD012" with rank 3 has parity type 14.
    Grid function "src_partial_D_hDD022" with rank 3 has parity type 15.
    Grid function "src_partial_D_hDD100" with rank 3 has parity type 16.
    Grid function "src_partial_D_hDD101" with rank 3 has parity type 17.
    Grid function "src_partial_D_hDD102" with rank 3 has parity type 18.
    Grid function "src_partial_D_hDD111" with rank 3 has parity type 19.
    Grid function "src_partial_D_hDD112" with rank 3 has parity type 20.
    Grid function "src_partial_D_hDD122" with rank 3 has parity type 21.
    Grid function "src_partial_D_hDD200" with rank 3 has parity type 22.
    Grid function "src_partial_D_hDD201" with rank 3 has parity type 23.
    Grid function "src_partial_D_hDD202" with rank 3 has parity type 24.
    Grid function "src_partial_D_hDD211" with rank 3 has parity type 25.
    Grid function "src_partial_D_hDD212" with rank 3 has parity type 26.
    Grid function "src_partial_D_hDD222" with rank 3 has parity type 27.
    Grid function "src_partial_D_WW0" with rank 1 has parity type 1.
    Grid function "src_partial_D_WW1" with rank 1 has parity type 2.
    Grid function "src_partial_D_WW2" with rank 1 has parity type 3.
    Grid function "src_trK" with rank 0 has parity type 0.
    Grid function "src_WW" with rank 0 has parity type 0.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 1: Construct a list of interp src gridfunctions; i.e., gridfunctions that act as the source for 1D interpolations.
    # List of gridfunction names with their corresponding ranks
    list_of_interp_src_gf_names_ranks: List[Tuple[str, int]] = []
    # Scalar gridfunctions
    list_of_interp_src_gf_names_ranks.append(("src_WW", 0))
    list_of_interp_src_gf_names_ranks.append(("src_trK", 0))
    # Rank-2 gridfunctions (hDD and aDD)
    for i in range(3):
        for j in range(i, 3):
            list_of_interp_src_gf_names_ranks.append((f"src_hDD{i}{j}", 2))
            list_of_interp_src_gf_names_ranks.append((f"src_aDD{i}{j}", 2))
    # Rank-1 gridfunctions (partial_D_WW)
    for i in range(3):
        list_of_interp_src_gf_names_ranks.append((f"src_partial_D_WW{i}", 1))
    # Rank-3 gridfunctions (partial_D_hDD)
    for k in range(3):
        for i in range(3):
            for j in range(i, 3):
                list_of_interp_src_gf_names_ranks.append(
                    (f"src_partial_D_hDD{k}{i}{j}", 3)
                )
    # Sort the list by gridfunction name
    list_of_interp_src_gf_names_ranks.sort(key=lambda x: x[0].upper())

    # Step 2: Register CodeParameters specific to this grid.
    # fmt: off
    for i in range(3):
        _ = par.CodeParameter("int", __name__, f"interp_src_Nxx{i}", 128, commondata=True, add_to_parfile=True)
        _ = par.CodeParameter("int", __name__, f"interp_src_Nxx_plus_2NGHOSTS{i}", 128, commondata=True,
                              add_to_parfile=True)
        _ = par.CodeParameter("REAL", __name__, f"interp_src_dxx{i}", 128, commondata=True, add_to_parfile=True)
        _ = par.CodeParameter("REAL", __name__, f"interp_src_invdxx{i}", 128, commondata=True, add_to_parfile=True)
    # fmt: on

    # Step 3: Register contributions to BHaH_defines.h and commondata.
    BHaH_defines_contrib = f"""
#define NUM_INTERP_SRC_GFS {len(list_of_interp_src_gf_names_ranks)} // Number of interp_src grid functions
enum {{
"""
    for name, _ in list_of_interp_src_gf_names_ranks:
        BHaH_defines_contrib += f"    {name.upper()}GF,\n"
    BHaH_defines_contrib += "};\n"

    # Finally, define parity types for interp_src gridfunctions.
    BHaH_defines_contrib += (
        locCBC.BHaH_defines_set_gridfunction_defines_with_parity_types(
            grid_name="interp_src",
            list_of_gf_names_ranks=list_of_interp_src_gf_names_ranks,
            verbose=True,
        )
    )

    BHaH_defines_h.register_BHaH_defines(__name__, BHaH_defines_contrib)

    griddata_commondata.register_griddata_commondata(
        __name__,
        "REAL *restrict interp_src_r_theta_phi[3]",
        "Source grid coordinates",
        is_commondata=True,
    )
    griddata_commondata.register_griddata_commondata(
        __name__,
        "REAL *restrict interp_src_gfs",
        f"{len(list_of_interp_src_gf_names_ranks)} 3D volume-filling gridfunctions with same angular sampling as evolved grid, but same radial sampling as input_gfs. GFs include: h_ij, h_ij,k, a_ij, trK, W, and W_,k. In RESCALED SPHERICAL basis.",
        is_commondata=True,
    )

    # Step 4: Register numgrid__interp_src_set_up().
    prefunc = ""
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Initializes the interp_src numerical grid, i.e., the source grid for 1D radial-spoke
interpolations during the hyperbolic relaxation.

This function sets up the numerical grid used as the "interpolation source" for metric data
on the evolved grids. It configures grid parameters, allocates memory for grid functions
and coordinate arrays, performs interpolation from external input, applies boundary conditions,
and computes necessary spatial derivatives.

@param commondata Pointer to the common data structure containing simulation parameters and data.
@param Nx_evol_grid Array specifying the number of grid points in each dimension for the evolved grid.
@return Returns BHAHAHA_SUCCESS on successful setup, or an error code if memory allocation fails."""
    cfunc_type = "int"
    name = "numgrid__interp_src_set_up"
    params = "commondata_struct *restrict commondata, const int Nx_evol_grid[3]"
    body = r"""
  int i0_min_shift = 0;
  if (commondata->bhahaha_params_and_data->r_min_external_input == 0)
    i0_min_shift = NGHOSTS;

  // Step 1: Configure grid parameters for the interpolation source.
  {
    // Align the radial grid with external input data.
    commondata->interp_src_Nxx0 = commondata->external_input_Nxx0;
    commondata->interp_src_Nxx1 = Nx_evol_grid[1];
    commondata->interp_src_Nxx2 = Nx_evol_grid[2];

    // Calculate grid sizes including ghost zones.
    commondata->interp_src_Nxx_plus_2NGHOSTS0 = commondata->interp_src_Nxx0 + 2 * NGHOSTS;
    commondata->interp_src_Nxx_plus_2NGHOSTS1 = commondata->interp_src_Nxx1 + 2 * NGHOSTS;
    commondata->interp_src_Nxx_plus_2NGHOSTS2 = commondata->interp_src_Nxx2 + 2 * NGHOSTS;

    // Set grid spacing based on external input and predefined angular ranges.
    commondata->interp_src_dxx0 = commondata->external_input_dxx0;
    const REAL xxmin1 = 0.0, xxmax1 = M_PI;
    const REAL xxmin2 = -M_PI, xxmax2 = M_PI;
    commondata->interp_src_dxx1 = (xxmax1 - xxmin1) / ((REAL)commondata->interp_src_Nxx1);
    commondata->interp_src_dxx2 = (xxmax2 - xxmin2) / ((REAL)commondata->interp_src_Nxx2);

    // Precompute inverse grid spacings for efficiency in derivative calculations.
    commondata->interp_src_invdxx0 = 1.0 / commondata->interp_src_dxx0;
    commondata->interp_src_invdxx1 = 1.0 / commondata->interp_src_dxx1;
    commondata->interp_src_invdxx2 = 1.0 / commondata->interp_src_dxx2;

    // Allocate memory for interpolation source grid functions.
    commondata->interp_src_gfs = malloc(sizeof(REAL) * commondata->interp_src_Nxx_plus_2NGHOSTS0 * commondata->interp_src_Nxx_plus_2NGHOSTS1 *
                                        commondata->interp_src_Nxx_plus_2NGHOSTS2 * NUM_INTERP_SRC_GFS);
    if (commondata->interp_src_gfs == NULL) {
      // Memory allocation failed for grid functions.
      return NUMGRID_INTERP_MALLOC_ERROR_GFS;
    }
  } // END STEP 1: Configure grid parameters for the interpolation source.

  // Step 2: Initialize coordinate arrays for the interpolation source grid.
  {
    // Step 2.a: Allocate memory for radial, theta, and phi coordinate arrays.
    commondata->interp_src_r_theta_phi[0] = (REAL *)malloc(sizeof(REAL) * commondata->interp_src_Nxx_plus_2NGHOSTS0);
    commondata->interp_src_r_theta_phi[1] = (REAL *)malloc(sizeof(REAL) * commondata->interp_src_Nxx_plus_2NGHOSTS1);
    commondata->interp_src_r_theta_phi[2] = (REAL *)malloc(sizeof(REAL) * commondata->interp_src_Nxx_plus_2NGHOSTS2);
    if (commondata->interp_src_r_theta_phi[0] == NULL || commondata->interp_src_r_theta_phi[1] == NULL ||
        commondata->interp_src_r_theta_phi[2] == NULL) {
      // Free previously allocated grid functions before exiting due to memory allocation failure.
      free(commondata->interp_src_gfs);
      return NUMGRID_INTERP_MALLOC_ERROR_RTHETAPHI;
    } // END IF memory allocation for coordinate arrays failed

    // Step 2.b: Populate coordinate arrays for a uniform, cell-centered spherical grid.
    const REAL xxmin1 = 0.0;
    const REAL xxmin2 = -M_PI;

    // Initialize radial coordinates by copying from external input.
    for (int j = 0; j < commondata->interp_src_Nxx_plus_2NGHOSTS0; j++)
      commondata->interp_src_r_theta_phi[0][j] = commondata->external_input_r_theta_phi[0][j];

    // Initialize theta coordinates with cell-centered values.
    for (int j = 0; j < commondata->interp_src_Nxx_plus_2NGHOSTS1; j++)
      commondata->interp_src_r_theta_phi[1][j] = xxmin1 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * commondata->interp_src_dxx1;

    // Initialize phi coordinates with cell-centered values.
    for (int j = 0; j < commondata->interp_src_Nxx_plus_2NGHOSTS2; j++)
      commondata->interp_src_r_theta_phi[2][j] = xxmin2 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * commondata->interp_src_dxx2;
  } // END STEP 2: Initialize coordinate arrays for the interpolation source grid.

  // Step 2.c: Extract grid sizes for use in indexing macros.
  const int Nxx_plus_2NGHOSTS0 = commondata->interp_src_Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = commondata->interp_src_Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = commondata->interp_src_Nxx_plus_2NGHOSTS2;

  // Step 3: Perform interpolation from external input to the interpolation source grid.
  // This involves multiple 2D interpolations corresponding to the radial grid and ghost zones.
  bah_interpolation_2d_external_input_to_interp_src_grid(commondata);

  // Step 4: Transfer interpolated data from external grid functions to interpolation source grid functions.
  {
    const REAL *restrict r_theta_phi[3] = {commondata->interp_src_r_theta_phi[0], commondata->interp_src_r_theta_phi[1],
                                           commondata->interp_src_r_theta_phi[2]};
    REAL *restrict in_gfs = commondata->interp_src_gfs;

#pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
      const MAYBE_UNUSED REAL xx2 = r_theta_phi[2][i2];
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
        const MAYBE_UNUSED REAL xx1 = r_theta_phi[1][i1];
        for (int i0 = i0_min_shift; i0 < Nxx_plus_2NGHOSTS0; i0++) {
          const MAYBE_UNUSED REAL xx0 = r_theta_phi[0][i0];
          // We perform this transformation in place; data read in will be written to the same points.
          const REAL external_Sph_W = in_gfs[IDX4(EXTERNAL_SPHERICAL_WWGF, i0, i1, i2)];
          const REAL external_Sph_trK = in_gfs[IDX4(EXTERNAL_SPHERICAL_TRKGF, i0, i1, i2)];
"""
    for i in range(3):
        for j in range(i, 3):
            body += f"const REAL external_Sph_hDD{i}{j} = in_gfs[IDX4(EXTERNAL_SPHERICAL_HDD{i}{j}GF, i0, i1, i2)];\n"
    for i in range(3):
        for j in range(i, 3):
            body += f"const REAL external_Sph_aDD{i}{j} = in_gfs[IDX4(EXTERNAL_SPHERICAL_ADD{i}{j}GF, i0, i1, i2)];\n"
    body += """
            in_gfs[IDX4(SRC_WWGF, i0, i1, i2)] = external_Sph_W;
            in_gfs[IDX4(SRC_TRKGF, i0, i1, i2)] = external_Sph_trK;
"""
    for i in range(3):
        for j in range(i, 3):
            body += (
                f"in_gfs[IDX4(SRC_HDD{i}{j}GF, i0, i1, i2)] = external_Sph_hDD{i}{j};\n"
            )
    for i in range(3):
        for j in range(i, 3):
            body += (
                f"in_gfs[IDX4(SRC_ADD{i}{j}GF, i0, i1, i2)] = external_Sph_aDD{i}{j};\n"
            )
    body += """
        } // END LOOP over i0
      } // END LOOP over i1
    } // END LOOP over i2
  } // END STEP 4: Transfer interpolated data to interpolation source grid functions.

  // Step 5: Initialize boundary condition structure for the interpolation source grid.
  bc_struct interp_src_bcstruct;
  {
    // Assign grid spacing and sizes to the boundary condition structure.
    commondata->bcstruct_dxx0 = commondata->interp_src_dxx0;
    commondata->bcstruct_dxx1 = commondata->interp_src_dxx1;
    commondata->bcstruct_dxx2 = commondata->interp_src_dxx2;
    commondata->bcstruct_Nxx_plus_2NGHOSTS0 = commondata->interp_src_Nxx_plus_2NGHOSTS0;
    commondata->bcstruct_Nxx_plus_2NGHOSTS1 = commondata->interp_src_Nxx_plus_2NGHOSTS1;
    commondata->bcstruct_Nxx_plus_2NGHOSTS2 = commondata->interp_src_Nxx_plus_2NGHOSTS2;

    // Set up boundary conditions based on the initialized grid.
    bah_bcstruct_set_up(commondata, commondata->interp_src_r_theta_phi, &interp_src_bcstruct);
  } // END STEP 5: Initialize boundary condition structure.

  // Step 6: Apply inner boundary conditions to specific grid functions to ensure smoothness.
  {
    // Step 6.a: Access boundary condition information from the boundary condition structure.
    const bc_info_struct *restrict bc_info = &interp_src_bcstruct.bc_info;

    // Step 6.b: Iterate over relevant grid functions and apply inner boundary conditions.
#pragma omp parallel
    for (int which_gf = 0; which_gf < NUM_INTERP_SRC_GFS; which_gf++) {
      switch (which_gf) {
      case SRC_WWGF:
      case SRC_HDD00GF:
      case SRC_HDD01GF:
      case SRC_HDD02GF:
      case SRC_HDD11GF:
      case SRC_HDD12GF:
      case SRC_HDD22GF: {
#pragma omp for
        for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) {
          const int dstpt = interp_src_bcstruct.inner_bc_array[pt].dstpt;
          const int srcpt = interp_src_bcstruct.inner_bc_array[pt].srcpt;

          // Apply boundary condition by copying and adjusting with parity.
          commondata->interp_src_gfs[IDX4pt(which_gf, dstpt)] =
              interp_src_bcstruct.inner_bc_array[pt].parity[interp_src_gf_parity[which_gf]] * commondata->interp_src_gfs[IDX4pt(which_gf, srcpt)];
        } // END LOOP over inner boundary points
        break;
      }
      default:
        // No boundary conditions needed for other grid functions.
        break;
      } // END SWITCH
    } // END LOOP over gridfunctions
  } // END STEP 6: Apply inner boundary conditions to specific grid functions.

  // Step 7: Compute spatial derivatives of h_{ij} within the interior of the interpolation source grid.
  bah_hDD_dD_and_W_dD_in_interp_src_grid_interior(commondata);

  // Step 8: Calculate radial derivatives at the outer boundaries using upwinding for stability.
  // If r_min is non-zero, apply the same procedure at the inner radial boundary.
  bah_apply_bcs_r_maxmin_partial_r_hDD_upwinding(commondata, commondata->interp_src_r_theta_phi, commondata->interp_src_gfs,
                                                 commondata->bhahaha_params_and_data->r_min_external_input != 0);

  // Step 9: Enforce boundary conditions on all interpolation source grid functions.
  {
    // Step 9.a: Access boundary condition information.
    const bc_info_struct *restrict bc_info = &interp_src_bcstruct.bc_info;

    // Step 9.b: Apply boundary conditions across all grid functions and boundary points.
#pragma omp parallel for collapse(2)
    for (int which_gf = 0; which_gf < NUM_INTERP_SRC_GFS; which_gf++) {
      for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) {
        const int dstpt = interp_src_bcstruct.inner_bc_array[pt].dstpt;
        const int srcpt = interp_src_bcstruct.inner_bc_array[pt].srcpt;

        // Apply boundary condition with parity correction for derivative calculations.
        commondata->interp_src_gfs[IDX4pt(which_gf, dstpt)] =
            interp_src_bcstruct.inner_bc_array[pt].parity[interp_src_gf_parity[which_gf]] * commondata->interp_src_gfs[IDX4pt(which_gf, srcpt)];
      } // END LOOP over inner boundary points
    } // END LOOP over gridfunctions
  } // END STEP 9: Enforce boundary conditions on all interpolation source grid functions.

  // Step 10: Release allocated memory for boundary condition structures.
  {
    free(interp_src_bcstruct.inner_bc_array);
    for (int ng = 0; ng < NGHOSTS * 3; ng++)
      free(interp_src_bcstruct.pure_outer_bc_array[ng]);
  } // END STEP 10: Free allocated memory for boundary condition structures.

  return BHAHAHA_SUCCESS;
"""
    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
