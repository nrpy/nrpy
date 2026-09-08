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

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, List, Tuple, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.infrastructures import BHaH


def register_CFunction_numgrid__interp_src_set_up() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for reading original source metric data.

    :return: None if in registration phase, else the updated NRPy environment.
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
    BHaH_defines_contrib += BHaH.BHaHAHA.bcstruct_set_up.BHaH_defines_set_gridfunction_defines_with_parity_types(
        grid_name="interp_src",
        list_of_gf_names_ranks=list_of_interp_src_gf_names_ranks,
        verbose=True,
    )

    interp_src_gf_index: Dict[str, int] = {
        name: idx for idx, (name, _) in enumerate(list_of_interp_src_gf_names_ranks)
    }
    interp_src_deriv_dst_dirn: List[int] = []
    interp_src_deriv_src_gf: List[List[int]] = []
    for name, _ in list_of_interp_src_gf_names_ranks:
        if name.startswith("src_partial_D_WW"):
            interp_src_deriv_dst_dirn.append(int(name[-1]))
            interp_src_deriv_src_gf.append(
                [
                    interp_src_gf_index[f"src_partial_D_WW{src_dirn}"]
                    for src_dirn in range(3)
                ]
            )
        elif name.startswith("src_partial_D_hDD"):
            interp_src_deriv_dst_dirn.append(int(name[-3]))
            interp_src_deriv_src_gf.append(
                [
                    interp_src_gf_index[
                        f"src_partial_D_hDD{src_dirn}{name[-2]}{name[-1]}"
                    ]
                    for src_dirn in range(3)
                ]
            )
        else:
            interp_src_deriv_dst_dirn.append(-1)
            interp_src_deriv_src_gf.append([-1, -1, -1])

    BHaH_defines_contrib += """// DERIVATIVE METADATA FOR INTERP_SRC GRID GRIDFUNCTIONS.
// interp_src_gf_parity stores base-field parity for stored coordinate derivatives.
// Values 0..2 in interp_src_gf_deriv_dst_dirn select the destination derivative direction;
// -1 marks non-derivative fields. interp_src_gf_deriv_src_gf maps each derivative
// family to the sibling source gridfunction for each mapped source derivative direction.
"""
    deriv_dst_dirn_values = ", ".join(map(str, interp_src_deriv_dst_dirn))
    deriv_src_gf_rows = ",\n  ".join(
        "{ " + ", ".join(map(str, row)) + " }" for row in interp_src_deriv_src_gf
    )
    BHaH_defines_contrib += rf"""static const int8_t interp_src_gf_deriv_dst_dirn[{len(list_of_interp_src_gf_names_ranks)}] = {{ {deriv_dst_dirn_values} }};
static const int16_t interp_src_gf_deriv_src_gf[{len(list_of_interp_src_gf_names_ranks)}][3] = {{
  {deriv_src_gf_rows}
}};
"""

    BHaH.BHaH_defines_h.register_BHaH_defines(__name__, BHaH_defines_contrib)

    BHaH.griddata_commondata.register_griddata_commondata(
        __name__,
        "REAL *restrict interp_src_r_theta_phi[3]",
        "Source grid coordinates",
        is_commondata=True,
    )
    BHaH.griddata_commondata.register_griddata_commondata(
        __name__,
        "REAL *restrict interp_src_gfs",
        f"{len(list_of_interp_src_gf_names_ranks)} 3D volume-filling gridfunctions with same angular sampling as evolved grid, but same radial sampling as input_gfs. GFs include: h_ij, h_ij,k, a_ij, trK, W, and W_,k. In RESCALED SPHERICAL basis.",
        is_commondata=True,
    )

    # Step 4: Register numgrid__interp_src_set_up().
    prefunc = r"""
/**
* This function is a substep of bah_numgrid__interp_src_set_up which Initializes the interp_src numerical grid, i.e., the source grid for 1D radial-spoke
* interpolations during the hyperbolic relaxation.
*
* This function initializes coordinate arrays and performs interpolation from external input
*
* @param commondata Pointer to the common data structure containing simulation parameters and data.
*/
#ifdef __CUDACC__
__global__
#endif
void initialize_and_interpolate(commondata_struct *restrict commondata)
{
#ifdef __CUDACC__
  // Set up cooperative group
  namespace cg = cooperative_groups;
  cg::grid_group gpu_grid = cg::this_grid();
#endif

  // Step 1: Extract grid sizes for use in indexing macros.
  const int Nxx_plus_2NGHOSTS0 = commondata->interp_src_Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = commondata->interp_src_Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = commondata->interp_src_Nxx_plus_2NGHOSTS2;

  {
    // Step 2: Populate coordinate arrays for a uniform, cell-centered spherical grid.
    const REAL xxmin1 = 0.0;
    const REAL xxmin2 = -M_PI;
  
    // Initialize radial coordinates by copying from external input.
    PARALLEL_1D_LOOP(j, 0, Nxx_plus_2NGHOSTS0) {
      commondata->interp_src_r_theta_phi[0][j] = commondata->external_input_r_theta_phi[0][j];
    } END_PARALLEL_1D_LOOP
    // Initialize theta coordinates with cell-centered values.
    PARALLEL_1D_LOOP(j, 0, Nxx_plus_2NGHOSTS1) {
      commondata->interp_src_r_theta_phi[1][j] = xxmin1 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * commondata->interp_src_dxx1;
    } END_PARALLEL_1D_LOOP
    // Initialize phi coordinates with cell-centered values.
    PARALLEL_1D_LOOP(j, 0, Nxx_plus_2NGHOSTS2) {
      commondata->interp_src_r_theta_phi[2][j] = xxmin2 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * commondata->interp_src_dxx2;
    } END_PARALLEL_1D_LOOP
  } // END STEP 2: Initialize coordinate arrays for the interpolation source grid.
#ifdef __CUDACC__
  gpu_grid.sync();
#endif

  // Step 3: Interpolate external data to interpolation source grid.
  bah_interpolation_2d_external_input_to_interp_src_grid(commondata);
#ifdef __CUDACC__
  gpu_grid.sync();
#endif
  if (commondata->error_flag != BHAHAHA_SUCCESS)
    return;
  // Step 4: Transfer interpolated data from external grid functions to interpolation source grid functions.
  {
     
    const REAL *restrict r_theta_phi[3] = {commondata->interp_src_r_theta_phi[0], commondata->interp_src_r_theta_phi[1],
                                           commondata->interp_src_r_theta_phi[2]};
    
    REAL *restrict in_gfs = commondata->interp_src_gfs;

  int i0_min_shift = 0;
  if (commondata->bhahaha_params_and_data->r_min_external_input == 0)
    i0_min_shift = NGHOSTS;

#ifdef __CUDACC__
    PARALLEL_LOOP(i0, i0_min_shift, Nxx_plus_2NGHOSTS0, i1, NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS, i2, NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS) {
      MAYBE_UNUSED const REAL xx2 = r_theta_phi[2][i2];
      MAYBE_UNUSED const REAL xx1 = r_theta_phi[1][i1];
      MAYBE_UNUSED const REAL xx0 = r_theta_phi[0][i0];
#else
#pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
      MAYBE_UNUSED const REAL xx2 = r_theta_phi[2][i2];
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
        MAYBE_UNUSED const REAL xx1 = r_theta_phi[1][i1];
        for (int i0 = i0_min_shift; i0 < Nxx_plus_2NGHOSTS0; i0++) {
          MAYBE_UNUSED const REAL xx0 = r_theta_phi[0][i0];
#endif
          // We perform this transformation in place; data read in will be written to the same points.
          const REAL external_Sph_W = in_gfs[IDX4(EXTERNAL_SPHERICAL_WWGF, i0, i1, i2)];
          const REAL external_Sph_trK = in_gfs[IDX4(EXTERNAL_SPHERICAL_TRKGF, i0, i1, i2)];
"""
    for i in range(3):
        for j in range(i, 3):
            prefunc += f"const REAL external_Sph_hDD{i}{j} = in_gfs[IDX4(EXTERNAL_SPHERICAL_HDD{i}{j}GF, i0, i1, i2)];\n"
    for i in range(3):
        for j in range(i, 3):
            prefunc += f"const REAL external_Sph_aDD{i}{j} = in_gfs[IDX4(EXTERNAL_SPHERICAL_ADD{i}{j}GF, i0, i1, i2)];\n"
    prefunc += """
            in_gfs[IDX4(SRC_WWGF, i0, i1, i2)] = external_Sph_W;
            in_gfs[IDX4(SRC_TRKGF, i0, i1, i2)] = external_Sph_trK;
"""
    for i in range(3):
        for j in range(i, 3):
            prefunc += (
                f"in_gfs[IDX4(SRC_HDD{i}{j}GF, i0, i1, i2)] = external_Sph_hDD{i}{j};\n"
            )
    for i in range(3):
        for j in range(i, 3):
            prefunc += (
                f"in_gfs[IDX4(SRC_ADD{i}{j}GF, i0, i1, i2)] = external_Sph_aDD{i}{j};\n"
            )
    prefunc += """
#ifndef __CUDACC__
        } // END LOOP over i0
      } // END LOOP over i1
    } // END LOOP over i2
#else
    } END_PARALLEL_LOOP
#endif
  } // END STEP 4: Transfer interpolated data to interpolation source grid functions.
}


#ifdef __CUDACC__
__constant__ int8_t c_interp_src_gf_parity[35];
/**
* This function is a substep of bah_numgrid__interp_src_set_up, which initializes the interp_src numerical grid, i.e., the source grid for 1D radial-spoke
* interpolations during the hyperbolic relaxation.
*
* This function applies boundary conditions, and computes necessary spatial derivatives.
*
* @param commondata Pointer to the common data structure containing simulation parameters and data.
* @param interp_sr_bcstruct Pointer to the bcstruct for the interp_src numerical grid.
*/

__global__
#endif
void apply_bcs_interp_src(commondata_struct *restrict commondata, bc_struct *restrict interp_src_bcstruct)
{
#ifdef __CUDACC__
  // Set interpolation source grid function parity array from device constant memory
  int8_t *interp_src_gf_parity = c_interp_src_gf_parity;

  //Set up cooperative group
  namespace cg = cooperative_groups;
  cg::grid_group gpu_grid = cg::this_grid();
#endif

  // Extract grid sizes for use in indexing macros.
  const int Nxx_plus_2NGHOSTS0 = commondata->interp_src_Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = commondata->interp_src_Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = commondata->interp_src_Nxx_plus_2NGHOSTS2;

  // Apply inner boundary conditions to specific grid functions to ensure smoothness.
  {
    // Step 1.a: Access boundary condition information from the boundary condition structure.
    const bc_info_struct *restrict bc_info = &interp_src_bcstruct->bc_info;

    // Step 1.b: Iterate over relevant grid functions and apply inner boundary conditions.
#ifndef __CUDACC__
#pragma omp parallel
#endif
    for (int which_gf = 0; which_gf < NUM_INTERP_SRC_GFS; which_gf++) {
      switch (which_gf) {
      case SRC_WWGF:
      case SRC_HDD00GF:
      case SRC_HDD01GF:
      case SRC_HDD02GF:
      case SRC_HDD11GF:
      case SRC_HDD12GF:
      case SRC_HDD22GF: {
#ifdef __CUDACC__
        PARALLEL_1D_LOOP(pt, 0, bc_info->num_inner_boundary_points)
#else
#pragma omp for
        for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) 
#endif          
        {
          const int dstpt = interp_src_bcstruct->inner_bc_array[pt].dstpt;
          const int srcpt = interp_src_bcstruct->inner_bc_array[pt].srcpt;

          // Apply boundary condition by copying and adjusting with parity.
          commondata->interp_src_gfs[IDX4pt(which_gf, dstpt)] =
              interp_src_bcstruct->inner_bc_array[pt].parity[interp_src_gf_parity[which_gf]] * commondata->interp_src_gfs[IDX4pt(which_gf, srcpt)];
        } // END LOOP over inner boundary points
        #ifdef __CUDACC__
        END_PARALLEL_1D_LOOP // END LOOP over inner boundary points
        #endif
        break;
      }
      default:
        // No boundary conditions needed for other grid functions.
        break;
      } // END SWITCH
    } // END LOOP over gridfunctions
  } // END STEP 1: Apply inner boundary conditions to specific grid functions.
  
  #ifdef __CUDACC__
  gpu_grid.sync();
  #endif

  // Step 2: Compute spatial derivatives of h_{ij} within the interior of the interpolation source grid.
  bah_hDD_dD_and_W_dD_in_interp_src_grid_interior(commondata);

  // Step 3: Calculate radial derivatives at the outer boundaries using upwinding for stability.
  // If r_min is non-zero, apply the same procedure at the inner radial boundary.
  bah_apply_bcs_r_maxmin_partial_r_hDD_upwinding(commondata, commondata->interp_src_r_theta_phi, commondata->interp_src_gfs,
                                                 commondata->bhahaha_params_and_data->r_min_external_input != 0);

  // Step 4: Enforce boundary conditions on all interpolation source grid functions.
  {
    // Step 4.a: Access boundary condition information.
    const bc_info_struct *restrict bc_info = &interp_src_bcstruct->bc_info;

    // Step 4.b: Apply boundary conditions across all grid functions and boundary points.
#ifdef __CUDACC__
    PARALLEL_2D_LOOP(pt, 0, bc_info->num_inner_boundary_points, which_gf, 0, NUM_INTERP_SRC_GFS) {
#else
#pragma omp parallel for collapse(2)
    for (int which_gf = 0; which_gf < NUM_INTERP_SRC_GFS; which_gf++) {
      for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) {
#endif        
        const int dstpt = interp_src_bcstruct->inner_bc_array[pt].dstpt;
        const int srcpt = interp_src_bcstruct->inner_bc_array[pt].srcpt;

        // Apply boundary condition with parity correction for derivative calculations.
        commondata->interp_src_gfs[IDX4pt(which_gf, dstpt)] =
            interp_src_bcstruct->inner_bc_array[pt].parity[interp_src_gf_parity[which_gf]] * commondata->interp_src_gfs[IDX4pt(which_gf, srcpt)];
#ifndef __CUDACC__
      } // END LOOP over inner boundary points
    } // END LOOP over gridfunctions
#else
    } END_PARALLEL_2D_LOOP // END LOOP over inner boundary points and gridfunctions
#endif
  } // END STEP 4: Enforce boundary conditions on all interpolation source grid functions.
} // END FUNCTION apply_bcs_interp_src

"""
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Initializes the interp_src numerical grid, i.e., the source grid for 1D radial-spoke
interpolations during the hyperbolic relaxation.

This function sets up the numerical grid used as the "interpolation source" for metric data
on the evolved grids. It configures grid parameters, allocates memory for grid functions
and coordinate arrays, performs interpolation from external input, applies boundary conditions,
and computes necessary spatial derivatives.

@param[in,out] commondata Pointer to the common data structure containing simulation parameters and data.
@param[in] Nx_evol_grid Array specifying the number of grid points in each dimension for the evolved grid.
@return Returns BHAHAHA_SUCCESS on successful setup, or an error code if memory allocation fails."""
    cfunc_type = "int"
    name = "numgrid__interp_src_set_up"
    params = "commondata_struct *restrict commondata, const int Nx_evol_grid[3]"
    body = r"""
#ifdef __CUDACC__
  //Allocate a device side copy of commondata
  commondata_struct *d_commondata = NULL;
  gpuErrchk( cudaMalloc((void**)&d_commondata, sizeof(commondata_struct)) );
#endif

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
  
  } // END STEP 1: Configure grid parameters for the interpolation source.

  // Step 2: Allocate interpolation source grid functions and coordinate arrays for the interpolation source grid.
  {
    // Step 2.a: Allocate memory for interpolation source grid functions.
#ifdef __CUDACC__
    cudaMalloc((void**)&commondata->interp_src_gfs, sizeof(REAL) * commondata->interp_src_Nxx_plus_2NGHOSTS0 * commondata->interp_src_Nxx_plus_2NGHOSTS1 *
                                        commondata->interp_src_Nxx_plus_2NGHOSTS2 * NUM_INTERP_SRC_GFS);
#else
    commondata->interp_src_gfs = (double*)malloc(sizeof(REAL) * commondata->interp_src_Nxx_plus_2NGHOSTS0 * commondata->interp_src_Nxx_plus_2NGHOSTS1 *
                                        commondata->interp_src_Nxx_plus_2NGHOSTS2 * NUM_INTERP_SRC_GFS);
#endif
    if (commondata->interp_src_gfs == NULL) {
      // Memory allocation failed for grid functions.
      commondata->error_flag = NUMGRID_INTERP_MALLOC_ERROR_GFS;
      return;
    }

    // Step 2.b: Allocate memory for radial, theta, and phi coordinate arrays.
#ifdef __CUDACC__
    cudaMalloc((void**)&commondata->interp_src_r_theta_phi[0], sizeof(REAL) * commondata->interp_src_Nxx_plus_2NGHOSTS0);
    cudaMalloc((void**)&commondata->interp_src_r_theta_phi[1], sizeof(REAL) * commondata->interp_src_Nxx_plus_2NGHOSTS1);
    cudaMalloc((void**)&commondata->interp_src_r_theta_phi[2], sizeof(REAL) * commondata->interp_src_Nxx_plus_2NGHOSTS2);
#else
    commondata->interp_src_r_theta_phi[0] = (REAL *)malloc(sizeof(REAL) * commondata->interp_src_Nxx_plus_2NGHOSTS0);
    commondata->interp_src_r_theta_phi[1] = (REAL *)malloc(sizeof(REAL) * commondata->interp_src_Nxx_plus_2NGHOSTS1);
    commondata->interp_src_r_theta_phi[2] = (REAL *)malloc(sizeof(REAL) * commondata->interp_src_Nxx_plus_2NGHOSTS2);
#endif
    if (commondata->interp_src_r_theta_phi[0] == NULL || commondata->interp_src_r_theta_phi[1] == NULL ||
        commondata->interp_src_r_theta_phi[2] == NULL) {
      // Free previously allocated grid functions before exiting due to memory allocation failure.
      FREE(commondata->interp_src_gfs);
      commondata->error_flag = NUMGRID_INTERP_MALLOC_ERROR_RTHETAPHI;
      return;
  } // END IF memory allocation for coordinate arrays failed
} // END STEP 2: Allocate interpolation source grid functions and coordinate arrays for the interpolation source grid.

#ifdef __CUDACC__
  // Update device side commondata
  gpuErrchk( cudaMemcpy(d_commondata, commondata, sizeof(commondata_struct), cudaMemcpyHostToDevice) );
#endif

  // Step 3: Initialize coordinate arrays and Perform interpolation from external input to the interpolation source grid.
  // This involves multiple 2D interpolations corresponding to the radial grid and ghost zones.
  {
#ifdef __CUDACC__
    void *Args[] = {&d_commondata};
    COOPERATIVE_KERNEL(initialize_and_interpolate, Args);
    gpuErrchk( cudaMemcpy(commondata, d_commondata, sizeof(commondata_struct), cudaMemcpyDeviceToHost) );
#else
    initialize_and_interpolate(commondata);
#endif 
    if (commondata->error_flag != BHAHAHA_SUCCESS)
      return;
  } 
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
#ifdef __CUDACC__
    // Update device side commondata 
    gpuErrchk( cudaMemcpy(d_commondata, commondata, sizeof(commondata_struct), cudaMemcpyHostToDevice) );
#endif

    // Set up boundary conditions based on the initialized grid.
#ifdef __CUDACC__
    bah_bcstruct_set_up(d_commondata, NULL, commondata->interp_src_r_theta_phi, &interp_src_bcstruct);
    gpuErrchk( cudaMemcpy(commondata, d_commondata, sizeof(commondata_struct), cudaMemcpyDeviceToHost) );
#else
    bah_bcstruct_set_up(commondata, NULL, commondata->interp_src_r_theta_phi, &interp_src_bcstruct);
#endif
    if (commondata->error_flag != BHAHAHA_SUCCESS)
      return;
  } // END STEP 5: Initialize boundary condition structure.


  // Step 6. Appply boundary conditions and compute necessary spacial derivatives
  {
#ifdef __CUDACC__
    // Create a device side copy of interp_src_bcstruct
    bc_struct *d_interp_src_bcstruct = NULL;
    cudaMalloc((void**)&d_interp_src_bcstruct, sizeof(bc_struct));
    cudaMemcpy(d_interp_src_bcstruct, &interp_src_bcstruct, sizeof(bc_struct), cudaMemcpyHostToDevice);

    // Set interpolation source grid function parity in device constant memory
    cudaMemcpyToSymbol(c_interp_src_gf_parity, &interp_src_gf_parity, 35*sizeof(int8_t));
#endif

    // Apply boundary conditions and compute necessary spacial derivatives
#ifdef __CUDACC__
    void *Args[] = {&d_commondata, &d_interp_src_bcstruct};
    COOPERATIVE_KERNEL(apply_bcs_interp_src, Args);
#else
    apply_bcs_interp_src(commondata, &interp_src_bcstruct);
#endif


    // Release allocated memory for boundary condition structures.
#ifdef __CUDACC__
    gpuErrchk( cudaFree(d_interp_src_bcstruct) );
#endif
    FREE(interp_src_bcstruct.inner_bc_array);
    for (int ng = 0; ng < NGHOSTS * 3; ng++)
      FREE(interp_src_bcstruct.pure_outer_bc_array[ng]);
  } //END STEP 6: Appply boundary conditions and compute necessary spacial derivatives

  //Copy then free device side copy of commondata
#ifdef __CUDACC__
  gpuErrchk( cudaFree(d_commondata) );
#endif

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
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
