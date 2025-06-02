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

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.jacobians as jac
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.infrastructures.BHaH.BHaHAHA.bcstruct_set_up as locCBC
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.infrastructures.BHaH import BHaH_defines_h, griddata_commondata


def register_CFunction_numgrid__external_input_set_up() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register the C function for reading original source metric data.

    :return: None if in registration phase, else the updated NRPy environment.

    DocTests:
    >>> env = register_CFunction_numgrid__external_input_set_up()
    Grid function "external_spherical_aDD00" with rank 2 has parity type 4.
    Grid function "external_spherical_aDD01" with rank 2 has parity type 5.
    Grid function "external_spherical_aDD02" with rank 2 has parity type 6.
    Grid function "external_spherical_aDD11" with rank 2 has parity type 7.
    Grid function "external_spherical_aDD12" with rank 2 has parity type 8.
    Grid function "external_spherical_aDD22" with rank 2 has parity type 9.
    Grid function "external_spherical_hDD00" with rank 2 has parity type 4.
    Grid function "external_spherical_hDD01" with rank 2 has parity type 5.
    Grid function "external_spherical_hDD02" with rank 2 has parity type 6.
    Grid function "external_spherical_hDD11" with rank 2 has parity type 7.
    Grid function "external_spherical_hDD12" with rank 2 has parity type 8.
    Grid function "external_spherical_hDD22" with rank 2 has parity type 9.
    Grid function "external_spherical_trK" with rank 0 has parity type 0.
    Grid function "external_spherical_WW" with rank 0 has parity type 0.
    Setting up reference_metric[Spherical]...
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 1: Construct a list of external input gridfunctions.
    # List of gridfunction names with their corresponding ranks
    list_of_external_input_gf_names_ranks: List[Tuple[str, int]] = []
    list_of_external_input_gf_names_ranks.append(("external_spherical_WW", 0))
    list_of_external_input_gf_names_ranks.append(("external_spherical_trK", 0))
    for i in range(3):
        for j in range(i, 3):
            list_of_external_input_gf_names_ranks.append(
                (f"external_spherical_hDD{i}{j}", 2)
            )
            list_of_external_input_gf_names_ranks.append(
                (f"external_spherical_aDD{i}{j}", 2)
            )
    # Sort the list by gridfunction name
    list_of_external_input_gf_names_ranks.sort(key=lambda x: x[0].upper())

    # Step 2: Register CodeParameters.
    # fmt: off
    for i in range(3):
        _ = par.CodeParameter("int", __name__, f"external_input_Nxx{i}", 128, commondata=True, add_to_parfile=True)
        _ = par.CodeParameter("int", __name__, f"external_input_Nxx_plus_2NGHOSTS{i}", 128, commondata=True,
                              add_to_parfile=True)
        _ = par.CodeParameter("REAL", __name__, f"external_input_dxx{i}", 128, commondata=True, add_to_parfile=True)
        _ = par.CodeParameter("REAL", __name__, f"external_input_invdxx{i}", 128, commondata=True, add_to_parfile=True)
    # fmt: on

    # Step 3: Register contributions to BHaH_defines.h and commondata.
    BHaH_defines_contrib = f"""
#define NUM_EXT_INPUT_CONFORMAL_GFS {len(list_of_external_input_gf_names_ranks)} // Number of external input grid functions
enum {{
"""
    for name, _ in list_of_external_input_gf_names_ranks:
        BHaH_defines_contrib += f"    {name.upper()}GF,\n"
    BHaH_defines_contrib += "};\n"
    BHaH_defines_contrib += (
        locCBC.BHaH_defines_set_gridfunction_defines_with_parity_types(
            grid_name="external_input",
            list_of_gf_names_ranks=list_of_external_input_gf_names_ranks,
            verbose=True,
        )
    )
    BHaH_defines_h.register_BHaH_defines(__name__, BHaH_defines_contrib)

    griddata_commondata.register_griddata_commondata(
        __name__,
        "REAL *restrict external_input_gfs_Cart_basis_no_gzs",
        f"{len(list_of_external_input_gf_names_ranks)} gridfunctions provided by external source, including gamma_ij and K_ij, in Cartesian basis, with no ghostzones.",
        is_commondata=True,
    )
    griddata_commondata.register_griddata_commondata(
        __name__,
        "REAL *restrict external_input_gfs",
        f"{len(list_of_external_input_gf_names_ranks)} gridfunctions provided by external source, including gamma_ij and K_ij, in spherical rescaled basis, with ghostzones.",
        is_commondata=True,
    )

    griddata_commondata.register_griddata_commondata(
        __name__,
        "REAL *restrict external_input_r_theta_phi[3]",
        "Three 1D arrays storing uniform (r, theta, phi) coordinates.",
        is_commondata=True,
    )

    # Step 4: Register numgrid__external_input_set_up().
    prefunc = r"""
// Indexing macros
#define EX_IDX4(g, i, j, k)                                                                                                                          \
  ((i) + commondata->external_input_Nxx_plus_2NGHOSTS0 *                                                                                             \
             ((j) + commondata->external_input_Nxx_plus_2NGHOSTS1 * ((k) + commondata->external_input_Nxx_plus_2NGHOSTS2 * (g))))
#define EX_NOGZ_IDX4(g, i, j, k)                                                                                                                     \
  ((i) + bhahaha_params_and_data->Nr_external_input * ((j) + commondata->external_input_Nxx1 * ((k) + commondata->external_input_Nxx2 * (g))))
"""
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r"""Initializes and processes the external input metric data (gamma_{ij}, K_{ij}).

This function performs the following steps:
1. Unpacks input parameters from the common data structure.
2. Adds ghost zones to the external input arrays to facilitate boundary condition application.
3. Allocates memory for the external input gridfunctions with ghost zones and assigns it to the common data structure.
4. Transfers metric data from arrays without ghost zones into the newly allocated arrays with ghost zones.
5. Sets up coordinate arrays for a uniform, cell-centered spherical grid.
6. Transforms the metric components (gamma_{ij}, K_{ij}) from Cartesian to spherical coordinates, including necessary rescaling.
7. Sets up boundary condition structures and applies inner boundary conditions, including parity corrections for all gridfunctions.

@param commondata - Pointer to the common data structure containing simulation parameters and data.
@param n_resolutions - Number of angular resolutions.
@param Ntheta - Array containing the number of theta points for each resolution.
@param Nphi - Array containing the number of phi points for each resolution.

@return BHAHAHA_SUCCESS on successful setup, or an error code indicating the failure reason.
"""
    cfunc_type = "int"
    name = "numgrid__external_input_set_up"
    params = "commondata_struct *restrict commondata, const int n_resolutions, const int *restrict Ntheta, const int *restrict Nphi"
    body = r"""
  // Step 1: Unpack input parameters from the common data structure.
  const bhahaha_params_and_data_struct *restrict bhahaha_params_and_data = commondata->bhahaha_params_and_data;

  // Calculate the number of interior (non-ghost) radial points by subtracting ghost zones.
  // Nr from external includes r ~ r_max NGHOSTS.
  commondata->external_input_Nxx0 = bhahaha_params_and_data->Nr_external_input - NGHOSTS;
  if (bhahaha_params_and_data->r_min_external_input > 0) {
    commondata->external_input_Nxx0 = bhahaha_params_and_data->Nr_external_input - 2 * NGHOSTS;
  }
  int i0_min_shift = 0;
  if (bhahaha_params_and_data->r_min_external_input == 0)
    i0_min_shift = NGHOSTS;

  // Set fixed angular resolutions for theta and phi directions.
  {
    const int max_resolution_i = bhahaha_params_and_data->num_resolutions_multigrid-1;
    commondata->external_input_Nxx1 = bhahaha_params_and_data->Ntheta_array_multigrid[max_resolution_i];
    commondata->external_input_Nxx2 = bhahaha_params_and_data->Nphi_array_multigrid[max_resolution_i];
  }

  // Step 1.a: Calculate grid spacing in each coordinate direction based on the simulation domain and resolution.
  // x_i = min_i + (j + 0.5) * dx_i, where dx_i = (max_i - min_i) / N_i

  commondata->external_input_dxx0 = bhahaha_params_and_data->dr_external_input;
  commondata->external_input_dxx1 = M_PI / ((REAL)commondata->external_input_Nxx1);
  commondata->external_input_dxx2 = 2 * M_PI / ((REAL)commondata->external_input_Nxx2);

  // Precompute inverse grid spacings for performance optimization in calculations.
  commondata->external_input_invdxx0 = 1.0 / commondata->external_input_dxx0;
  commondata->external_input_invdxx1 = 1.0 / commondata->external_input_dxx1;
  commondata->external_input_invdxx2 = 1.0 / commondata->external_input_dxx2;

  // Step 2: Add ghost zones to the external input data.
  // Ghost zones are added to so that inner boundary conditions may be applied; 2 * NGHOSTS in each angular direction and NGHOSTS in the radial
  // direction.
  commondata->external_input_Nxx_plus_2NGHOSTS0 = commondata->external_input_Nxx0 + 2 * NGHOSTS;
  commondata->external_input_Nxx_plus_2NGHOSTS1 = commondata->external_input_Nxx1 + 2 * NGHOSTS;
  commondata->external_input_Nxx_plus_2NGHOSTS2 = commondata->external_input_Nxx2 + 2 * NGHOSTS;

  // Calculate the total number of grid points including ghost zones.
  const int total_elements_incl_gzs =
      commondata->external_input_Nxx_plus_2NGHOSTS0 * commondata->external_input_Nxx_plus_2NGHOSTS1 * commondata->external_input_Nxx_plus_2NGHOSTS2;

  // Pointers to the external input gridfunctions without ghost zones.
  REAL *restrict external_input_gfs_no_gzs = commondata->external_input_gfs_Cart_basis_no_gzs;

  // Allocate memory for the external input gridfunctions with ghost zones.
  REAL *restrict external_input_gfs = (REAL *)malloc(NUM_EXT_INPUT_CONFORMAL_GFS * total_elements_incl_gzs * sizeof(REAL));
  if (external_input_gfs == NULL) {
    return NUMGRID_EXTERN_MALLOC_ERROR_GFS;
  } // END IF memory allocation for external_input_gfs failed

  // Step 3: Assign the allocated array to commondata for use outside this function.
  commondata->external_input_gfs = external_input_gfs;

  // Step 4: Transfer data from the no-ghost zones array to the array with ghost zones.
  // This involves copying metric data (gamma_{ij} and K_{ij}) from the Cartesian basis into the newly allocated arrays with ghost zones.
  LOOP_OMP("omp parallel for",                        //
           i0, 0, bhahaha_params_and_data->Nr_external_input, //
           i1, 0, commondata->external_input_Nxx1,    //
           i2, 0, commondata->external_input_Nxx2) {
    for (int gf = 0; gf < NUM_EXT_INPUT_CARTESIAN_GFS; gf++) {
      external_input_gfs[EX_IDX4(gf, i0 + i0_min_shift, i1 + NGHOSTS, i2 + NGHOSTS)] = external_input_gfs_no_gzs[EX_NOGZ_IDX4(gf, i0, i1, i2)];
    }
  } // END LOOP: iterating through the external input grid points

  // Step 5: Set up coordinate arrays for a uniform, cell-centered spherical grid.
  {
    const int Nxx_plus_2NGHOSTS0 = commondata->external_input_Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = commondata->external_input_Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = commondata->external_input_Nxx_plus_2NGHOSTS2;

    // Step 5.a: Allocate memory for coordinate arrays in radial, theta, and phi directions.
    commondata->external_input_r_theta_phi[0] = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS0);
    commondata->external_input_r_theta_phi[1] = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS1);
    commondata->external_input_r_theta_phi[2] = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS2);
    if (commondata->external_input_r_theta_phi[0] == NULL || commondata->external_input_r_theta_phi[1] == NULL ||
        commondata->external_input_r_theta_phi[2] == NULL) {
      free(external_input_gfs);
      return NUMGRID_EXTERN_MALLOC_ERROR_RTHETAPHI;
    } // END IF memory allocation for external_input_r_theta_phi arrays failed

    // Step 5.b: Initialize coordinate arrays for a uniform, cell-centered spherical grid.
    // The coordinates are centered within each cell by adding 0.5 to the index before scaling.
    const REAL xxmin0 = bhahaha_params_and_data->r_min_external_input;
    const REAL xxmin1 = 0.0;
    const REAL xxmin2 = -M_PI;

    for (int j = 0; j < Nxx_plus_2NGHOSTS0; j++)
      commondata->external_input_r_theta_phi[0][j] = xxmin0 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * commondata->external_input_dxx0;
    for (int j = 0; j < Nxx_plus_2NGHOSTS1; j++)
      commondata->external_input_r_theta_phi[1][j] = xxmin1 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * commondata->external_input_dxx1;
    for (int j = 0; j < Nxx_plus_2NGHOSTS2; j++)
      commondata->external_input_r_theta_phi[2][j] = xxmin2 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * commondata->external_input_dxx2;
  } // END BLOCK: setting up coordinate arrays

  // Step 6: Transform the metric components (gamma_{ij}, K_{ij}) from Cartesian to spherical coordinates,
  // including necessary rescaling.
  {
    const int Nxx_plus_2NGHOSTS0 = commondata->external_input_Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = commondata->external_input_Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = commondata->external_input_Nxx_plus_2NGHOSTS2;

    // Step 6.a: Extract coordinate arrays for easier access during transformation.
    REAL *restrict external_input_r_theta_phi[3];
    for (int ww = 0; ww < 3; ww++)
      external_input_r_theta_phi[ww] = commondata->external_input_r_theta_phi[ww];

    // Step 6.b: Metric components: basis transform from Cartesian to Spherical & convert ADM->rescaled BSSN.
#pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < commondata->external_input_Nxx2 + NGHOSTS; i2++) {
      const REAL xx2 = external_input_r_theta_phi[2][i2];
      for (int i1 = NGHOSTS; i1 < commondata->external_input_Nxx1 + NGHOSTS; i1++) {
        const REAL xx1 = external_input_r_theta_phi[1][i1];
        // Include all valid points, including those near r ~ r_max.
        for (int i0 = i0_min_shift; i0 < commondata->external_input_Nxx_plus_2NGHOSTS0; i0++) {
          const REAL xx0 = external_input_r_theta_phi[0][i0];

          // Read Cartesian metric components at the current grid point.
"""
    labels = ["X", "Y", "Z"]
    for prefix in ["gammaDD", "KDD"]:
        for i in range(3):
            for j in range(i, 3):
                body += f"const REAL Cart_{prefix}{i}{j} = external_input_gfs[IDX4(INTERP_{prefix.upper()}{labels[i]}{labels[j]}GF, i0, i1, i2)];\n"
    # Cartesian -> Spherical basis transform
    Cart_gammaDD = ixp.declarerank2("Cart_gammaDD", symmetry="sym01")
    Sph_gammaDD = jac.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
        "Spherical", Cart_gammaDD
    )
    Cart_KDD = ixp.declarerank2("Cart_KDD", symmetry="sym01")
    Sph_KDD = jac.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
        "Spherical", Cart_KDD
    )
    rfm = refmetric.reference_metric["Spherical"]

    hDD = ixp.zerorank2()
    gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(Sph_gammaDD)  # _gammaUU unused.
    gammabarDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            gammabarDD[i][j] = (rfm.detgammahat / gammaDET) ** (
                sp.Rational(1, 3)
            ) * Sph_gammaDD[i][j]
            # gammabar_{ij}  = h_{ij}*ReDD[i][j] + gammahat_{ij}
            # -> h_{ij} = ( gammabar_{ij}  -   gammahat_{ij} ) / ReDD[i][j]
            hDD[i][j] = (gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[i][j]
    # Derivation of W:
    # gamma_{ij} = ADM "physical" 3-metric.
    # gammabar_{ij} = BSSN conformal 3-metric
    # gammahat_{ij} = reference metric.
    # gamma_{ij} = e^{4 phi} gammabar_{ij}
    # det(gamma_{ij}) = det(e^{4 phi} gammabar_{ij})
    #                 = e^{12 phi} det(gammabar_{ij})
    #                 = e^{12 phi} det(gammahat_{ij})  # We are allowed to choose det(gammabar) = det(gammahat)
    # -> e^{12 phi} = det(gamma_{ij}) / det(gammahat_{ij})
    # -> W = e^{-2 phi} = (e^{12 phi})^{-1/6} = (det(gamma_{ij}) / det(gammahat_{ij}))^{-1/6}
    # We are allowed to choose det(gammabar) = det(gammahat)
    # Finally, note that sometimes gammaDET < 0 after interpolation. This happens
    #    only at the puncture, which should be far from the horizon provided dr is small enough.
    #    To correct this and avoid NaNs in our evolutions (caused by [negative number]^{-1/6}, we
    #    simply take the absolute value of gammaDET/gammahatDET when computing W.
    W = (sp.Abs(gammaDET / rfm.detgammahat)) ** (-sp.Rational(1, 6))

    # K = gamma^{ij} K_{ij}, and
    # \bar{A}_{ij} &= (\frac{\bar{gamma}}{gamma})^{1/3}*(K_{ij} - \frac{1}{3}*gamma_{ij}*K)
    trK = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            trK += gammaUU[i][j] * Sph_KDD[i][j]
    AbarDD = ixp.zerorank2()
    aDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            AbarDD[i][j] = (rfm.detgammahat / gammaDET) ** (sp.Rational(1, 3)) * (
                Sph_KDD[i][j] - sp.Rational(1, 3) * Sph_gammaDD[i][j] * trK
            )
            aDD[i][j] = AbarDD[i][j] / rfm.ReDD[i][j]

    external_input_gfs_no_derivs_names: List[str] = []
    external_input_gfs_no_derivs_exprs: List[sp.Expr] = []
    external_input_gfs_no_derivs_names += [
        "external_input_gfs[IDX4(EXTERNAL_SPHERICAL_WWGF, i0, i1, i2)]"
    ]
    external_input_gfs_no_derivs_exprs += [W]

    external_input_gfs_no_derivs_names += [
        "external_input_gfs[IDX4(EXTERNAL_SPHERICAL_TRKGF, i0, i1, i2)]"
    ]
    external_input_gfs_no_derivs_exprs += [trK]
    for i in range(3):
        for j in range(i, 3):
            external_input_gfs_no_derivs_names += [
                f"external_input_gfs[IDX4(EXTERNAL_SPHERICAL_HDD{i}{j}GF, i0, i1, i2)]"
            ]
            external_input_gfs_no_derivs_exprs += [hDD[i][j]]
    for i in range(3):
        for j in range(i, 3):
            external_input_gfs_no_derivs_names += [
                f"external_input_gfs[IDX4(EXTERNAL_SPHERICAL_ADD{i}{j}GF, i0, i1, i2)]"
            ]
            external_input_gfs_no_derivs_exprs += [aDD[i][j]]

    # Input gf names & expressions into c_codegen, to complete the conversion.
    body += ccg.c_codegen(
        external_input_gfs_no_derivs_exprs,
        external_input_gfs_no_derivs_names,
        verbose=False,
        include_braces=False,
    )
    body += """
        } // END LOOP over i0
      } // END LOOP over i1
    } // END LOOP over i2
  } // END BLOCK: transformation and rescaling

  // Step 7: Set up boundary condition structures and apply inner boundary conditions.
  {
    const int Nxx_plus_2NGHOSTS0 = commondata->external_input_Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = commondata->external_input_Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = commondata->external_input_Nxx_plus_2NGHOSTS2;

    // Assign grid spacings and sizes to the boundary condition structure within commondata.
    commondata->bcstruct_dxx0 = commondata->external_input_dxx0;
    commondata->bcstruct_dxx1 = commondata->external_input_dxx1;
    commondata->bcstruct_dxx2 = commondata->external_input_dxx2;

    commondata->bcstruct_Nxx_plus_2NGHOSTS0 = Nxx_plus_2NGHOSTS0;
    commondata->bcstruct_Nxx_plus_2NGHOSTS1 = Nxx_plus_2NGHOSTS1;
    commondata->bcstruct_Nxx_plus_2NGHOSTS2 = Nxx_plus_2NGHOSTS2;

    // Initialize the boundary condition structure for external input data.
    bc_struct external_input_bcstruct;
    bah_bcstruct_set_up(commondata, commondata->external_input_r_theta_phi, &external_input_bcstruct);

    // Step 7.a: Unpack boundary condition information from the boundary condition structure.
    const bc_info_struct *restrict bc_info = &external_input_bcstruct.bc_info;

    // Step 7.b: Apply inner boundary conditions to all gridfunctions.
    // This involves copying values from source points to destination points with parity corrections.
#pragma omp parallel
    for (int which_gf = 0; which_gf < NUM_EXT_INPUT_CONFORMAL_GFS; which_gf++) {
#pragma omp for
      for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) {
        const int dstpt = external_input_bcstruct.inner_bc_array[pt].dstpt;
        const int srcpt = external_input_bcstruct.inner_bc_array[pt].srcpt;
        // Apply the boundary condition by copying values from the source point to the destination point,
        // applying the appropriate parity correction for the gridfunction.
        commondata->external_input_gfs[IDX4pt(which_gf, dstpt)] =
            external_input_bcstruct.inner_bc_array[pt].parity[external_input_gf_parity[which_gf]] *
            commondata->external_input_gfs[IDX4pt(which_gf, srcpt)];
      } // END LOOP over inner boundary points
    } // END LOOP over gridfunctions

    // Step 7.c: Free allocated memory for boundary condition structures to prevent memory leaks.
    free(external_input_bcstruct.inner_bc_array);
    for (int ng = 0; ng < NGHOSTS * 3; ng++)
      free(external_input_bcstruct.pure_outer_bc_array[ng]);
  } // END BLOCK: applying boundary conditions

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
