"""
Register C function for finding horizons with BHaHAHA.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, List, Tuple, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.jacobians as jac
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.infrastructures.BHaH import griddata_commondata


def register_CFunction_bhahaha_find_horizons(
    BHaHAHA_subdir: str,
    CoordSystem: str,
    max_horizons: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for general-purpose 3D Lagrange interpolation.

    :param BHaHAHA_subdir: Subdirectory where BHaHAHA is located.
    :param CoordSystem: CoordSystem of project, where horizon finding will take place.
    :param max_horizons: Maximum number of horizons to search for.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If EvolvedConformalFactor_cf set to unsupported value.

    >>> env = register_CFunction_bhahaha_find_horizons()
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    griddata_commondata.register_griddata_commondata(
        __name__,
        f"bhahaha_params_and_data_struct bhahaha_params_and_data[{max_horizons}]",
        "BHaHAHA parameters and data, including previous horizon data",
        is_commondata=True,
    )
    griddata_commondata.register_griddata_commondata(
        __name__,
        f"bhahaha_diagnostics_struct bhahaha_diagnostics[{max_horizons}]",
        "BHaHAHA diagnostics",
        is_commondata=True,
    )

    # Register commondata parameters
    # Dictionary for REAL type parameters
    real_params = {
        "bah_max_search_radius": ("Maximum radius to search for horizons", 1.5),
        "bah_cfl_factor": ("CFL factor for timestep. 0.98 recommended.", 0.98),
        "bah_Theta_Linf_times_M_tolerance": (
            "Convergence criterion: Linf norm of residual times M_scale, 1e-2 recommended",
            1e-2,
        ),
        "bah_Theta_L2_times_M_tolerance": (
            "Convergence criterion: L2 norm of residual times M_scale, 2e-5 recommended",
            2e-5,
        ),
        "bah_eta_damping_times_M": (
            "Exponential damping parameter, 3.5 recommended",
            7.0,
        ),
        "bah_M_scale": ("Mass scale for horizon", 1.0),
        "bah_KO_strength": (
            "Kreiss-Oliger dissipation parameter, 0.0 strongly recommended",
            0.0,
        ),
        "bah_initial_grid_x_center": ("Initial x-center coordinates", 0.0),
        "bah_initial_grid_y_center": ("Initial y-center coordinates", 0.0),
        "bah_initial_grid_z_center": ("Initial z-center coordinates", 0.0),
    }
    for key, value in real_params.items():
        _ = par.register_CodeParameter(
            cparam_type=f"REAL[{max_horizons}]",
            module="BHaHAHA",
            name=key,
            defaultvalue=value[1],
            commondata=True,
            add_to_parfile=True,
            add_to_set_CodeParameters_h=False,
            add_to_glb_code_params_dict=True,
            description=value[0],
        )

    # Dictionary for INT type parameters
    int_params = {
        "bah_Nr_interp_max": (
            "Maximum number of radial grid points used in 3D interpolation",
            48,
        ),
        "bah_verbosity_level": ("Verbosity level", 2),
        "bah_max_iterations": (
            "Maximum iterations before giving up. 10000 recommended",
            10000,
        ),
        "bah_max_num_horizons": (
            "Maximum number of horizons to find",
            max_horizons,
        ),
        "bah_num_resolutions_multigrid": (
            "Number of resolutions for low-to-high multigrid pass",
            3,
        ),
        "bah_BBH_mode_enable": (
            "BBH mode? Enable=1 ; Disable=0",
            0,
        ),
        "bah_BBH_mode_common_horizon_idx": (
            "BBH mode: Horizon index (zero-offset) for common (merged BH) horizon",
            2,
        ),
        "bah_enable_eta_varying_alg_for_precision_common_horizon": (
            "Enable varying-eta prescription for precision common-horizon finding",
            0,
        ),
    }
    for key, value in int_params.items():
        _ = par.register_CodeParameter(
            cparam_type="int",
            module="BHaHAHA",
            name=key,
            defaultvalue=value[1],
            commondata=True,
            add_to_parfile=True,
            add_to_set_CodeParameters_h=False,
            add_to_glb_code_params_dict=True,
            description=value[0],
        )

    # Dictionary for INT array type parameters
    int_array_params: Dict[str, Tuple[str, List[int]]] = {
        "bah_Ntheta_array_multigrid": (
            "Array of theta resolutions for multigrid cycle",
            [8, 16, 32],
        ),
        "bah_Nphi_array_multigrid": (
            "Array of phi resolutions for multigrid cycle",
            [16, 32, 64],
        ),
        "bah_BBH_mode_inspiral_BH_idxs": (
            "BBH mode: Horizon indices (zero-offset) for inspiralling BHs",
            [0, 1],
        ),
    }
    for key, intvalue in int_array_params.items():
        _ = par.register_CodeParameter(
            cparam_type=f"int[{len(intvalue[1])}]",
            module="BHaHAHA",
            name=key,
            defaultvalue=intvalue[1],
            commondata=True,
            add_to_parfile=True,
            add_to_set_CodeParameters_h=False,
            add_to_glb_code_params_dict=True,
            description=intvalue[0],
        )

    _ = par.register_CodeParameter(
        cparam_type=f"int[{max_horizons}]",
        module="BHaHAHA",
        name="bah_BBH_mode_horizon_active",
        defaultvalue=[1] * max_horizons,
        commondata=True,
        add_to_parfile=False,
        add_to_set_CodeParameters_h=False,
        add_to_glb_code_params_dict=True,
        description="BBH mode: Record of which horizons are active.",
    )

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        f"{BHaHAHA_subdir}/BHaHAHA.h",
    ]
    prefunc = r"""
#define BHAHAHA_NUM_INTERP_GFS 14

// BSSN gridfunctions input into the interpolator. Must be in the same order as the enum list below.
const int bhahaha_gf_interp_indices[BHAHAHA_NUM_INTERP_GFS] = {
    ADD00GF, ADD01GF, ADD02GF, ADD11GF, ADD12GF, ADD22GF, // Traceless, rescaled extrinsic curvature components.
    CFGF,                                                 // Conformal factor.
    HDD00GF, HDD01GF, HDD02GF, HDD11GF, HDD12GF, HDD22GF, // Rescaled conformal 3-metric components.
    TRKGF                                                 // Trace of extrinsic curvature.
};

// Rescaled BSSN gridfunctions output by the interpolator. Must be in the same order as the list above.
enum {
  INTERP_ADD00GF,
  INTERP_ADD01GF,
  INTERP_ADD02GF,
  INTERP_ADD11GF,
  INTERP_ADD12GF,
  INTERP_ADD22GF,
  INTERP_CFGF,
  INTERP_HDD00GF,
  INTERP_HDD01GF,
  INTERP_HDD02GF,
  INTERP_HDD11GF,
  INTERP_HDD12GF,
  INTERP_HDD22GF,
  INTERP_TRKGF
};

// Transformed gridfunctions in the ADM formalism after basis transformation.
enum {
  INTERP_GAMMADD00GF,
  INTERP_GAMMADD01GF,
  INTERP_GAMMADD02GF,
  INTERP_GAMMADD11GF,
  INTERP_GAMMADD12GF,
  INTERP_GAMMADD22GF,
  INTERP_KDD00GF,
  INTERP_KDD01GF,
  INTERP_KDD02GF,
  INTERP_KDD11GF,
  INTERP_KDD12GF,
  INTERP_KDD22GF
};

/**
 * Initializes BHaHAHA parameters and data structures for each horizon.
 *
 * This function allocates memory for previous horizon data and initializes various parameters
 * required for the BHaHAHA horizon finding algorithm. It iterates over all horizons specified
 * in the common data structure and sets up initial values for time, position, and search parameters.
 *
 * @param commondata Pointer to the common data structure containing simulation parameters.
 *
 * @note Modifies the commondata structure by initializing its bhahaha_params_and_data array.
 */
void initialize_bhahaha_params_and_data_struct(commondata_struct *restrict commondata) {
  for (int which_horizon = 0; which_horizon < commondata->bah_max_num_horizons; which_horizon++) {
    bhahaha_params_and_data_struct *restrict bhahaha_params_and_data = &commondata->bhahaha_params_and_data[which_horizon];

    //------------------------
    // Horizon-tracking & interpolation setup
    //------------------------
    // Allocate storage for previous horizons' data.
    const int Ntheta_max = commondata->bah_Ntheta_array_multigrid[commondata->bah_num_resolutions_multigrid - 1];
    const int Nphi_max = commondata->bah_Nphi_array_multigrid[commondata->bah_num_resolutions_multigrid - 1];
    bhahaha_params_and_data->prev_horizon_m1 = malloc(Ntheta_max * Nphi_max * sizeof(REAL));
    bhahaha_params_and_data->prev_horizon_m2 = malloc(Ntheta_max * Nphi_max * sizeof(REAL));
    bhahaha_params_and_data->prev_horizon_m3 = malloc(Ntheta_max * Nphi_max * sizeof(REAL));
    // t_m? = -1.0 indicates to the quadratic extrapolation algorithm that this horizon has not yet been found.
    bhahaha_params_and_data->t_m1 = bhahaha_params_and_data->t_m2 = bhahaha_params_and_data->t_m3 = -1.0;
    bhahaha_params_and_data->x_center_m1 = commondata->bah_initial_grid_x_center[which_horizon];
    bhahaha_params_and_data->y_center_m1 = commondata->bah_initial_grid_y_center[which_horizon];
    bhahaha_params_and_data->z_center_m1 = commondata->bah_initial_grid_z_center[which_horizon];
    bhahaha_params_and_data->r_min_m1 = 0.0;
    bhahaha_params_and_data->r_max_m1 = commondata->bah_max_search_radius[which_horizon];
    //------------------------
    // Parameters common to all horizons:
    //------------------------
    // Initialize evolution parameters required for horizon finding.
    bhahaha_params_and_data->enable_eta_varying_alg_for_precision_common_horizon =
        commondata->bah_enable_eta_varying_alg_for_precision_common_horizon;
    bhahaha_params_and_data->use_fixed_radius_guess_on_full_sphere = 1;
    bhahaha_params_and_data->verbosity_level = commondata->bah_verbosity_level;
    bhahaha_params_and_data->max_iterations = commondata->bah_max_iterations;
    //------------------------
    // Per-horizon parameters:
    //------------------------
    bhahaha_params_and_data->cfl_factor = commondata->bah_cfl_factor[which_horizon];
    bhahaha_params_and_data->Theta_L2_times_M_tolerance = commondata->bah_Theta_L2_times_M_tolerance[which_horizon];
    bhahaha_params_and_data->Theta_Linf_times_M_tolerance = commondata->bah_Theta_Linf_times_M_tolerance[which_horizon];
    bhahaha_params_and_data->M_scale = commondata->bah_M_scale[which_horizon];
    bhahaha_params_and_data->eta_damping_times_M = commondata->bah_eta_damping_times_M[which_horizon];
    bhahaha_params_and_data->KO_strength = commondata->bah_KO_strength[which_horizon];

    // Set multigrid cycle resolutions.
    bhahaha_params_and_data->num_resolutions_multigrid = commondata->bah_num_resolutions_multigrid;
    memcpy(bhahaha_params_and_data->Ntheta_array_multigrid, commondata->bah_Ntheta_array_multigrid, sizeof(int) * commondata->bah_num_resolutions_multigrid);
    memcpy(bhahaha_params_and_data->Nphi_array_multigrid, commondata->bah_Nphi_array_multigrid, sizeof(int) * commondata->bah_num_resolutions_multigrid);
  } // END LOOP: Initialize parameters for each horizon.
} // END FUNCTION: initialize_bhahaha_params_and_data_struct()

/**
 * Interpolates BSSN gridfunctions onto a spherical grid and transforms them into Cartesian ADM variables.
 *
 * This function performs the following operations:
 * 1. Allocates memory for the input metric data array.
 * 2. Sets up a spherical grid with specified radial and angular resolutions.
 * 3. Converts grid points from spherical to reference-metric coordinates.
 * 4. Interpolates the BSSN gridfunctions onto the spherical grid points.
 * 5. Transforms the interpolated data from BSSN to ADM variables in the Cartesian basis.
 * 6. Stores the transformed metric data in the output array for further processing.
 *
 * @param commondata       Pointer to common simulation data structures.
 * @param params           Pointer to simulation parameters.
 * @param xx               Array of pointers to reference-metric coordinate arrays.
 * @param y_n_gfs          Pointer to input gridfunctions for interpolation.
 * @param bhahaha_params   Pointer to BHaHAHA parameters and data structures.
 * @param x_center         x-coordinate of the spherical grid center in Cartesian coordinates.
 * @param y_center         y-coordinate of the spherical grid center in Cartesian coordinates.
 * @param z_center         z-coordinate of the spherical grid center in Cartesian coordinates.
 * @param radii            Array of radial coordinates for the spherical grid.
 * @param input_metric_data Output pointer to store the interpolated and transformed metric data.
 *
 */
static void BHaHAHA_interpolate_metric_data(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],
                                            const REAL *restrict y_n_gfs, bhahaha_params_and_data_struct *bhahaha_params, const REAL x_center,
                                            const REAL y_center, const REAL z_center, const REAL radii[commondata->bah_Nr_interp_max],
                                            REAL *restrict input_metric_data) {

  // Define spherical grid resolution parameters.
  const int Ntheta_max = commondata->bah_Ntheta_array_multigrid[commondata->bah_num_resolutions_multigrid - 1];
  const int Nphi_max = commondata->bah_Nphi_array_multigrid[commondata->bah_num_resolutions_multigrid - 1];
  const int Ntheta = Ntheta_max;               // Number of theta divisions.
  const int Nphi = Nphi_max;                   // Number of phi divisions.
  const REAL dtheta = M_PI / ((REAL)Ntheta);   // Theta grid spacing.
  const REAL dphi = 2.0 * M_PI / ((REAL)Nphi); // Phi grid spacing.

  // Calculate the total number of grid points in the spherical grid.
  const int total_elements = bhahaha_params->Nr_external_input * Ntheta * Nphi;

  // Allocate memory for destination points on the spherical grid.
  // Each grid point stores Cartesian coordinates (x0, x1, x2) in reference-metric coordinates.
  REAL(*dst_x0x1x2)[3] = (REAL(*)[3])malloc(total_elements * 3 * sizeof(REAL));

  // Verify successful memory allocation for destination coordinates.
  if (dst_x0x1x2 == NULL) {
    fprintf(stderr, "Failed to allocate memory for dst_x0x1x2.\n");
    exit(EXIT_FAILURE);
  }

// Macro to compute a linear index for spherical grid coordinates.
#define IDX3_SPH_INTERP(ir, itheta, iphi) ((ir) + bhahaha_params->Nr_external_input * ((itheta) + Ntheta * (iphi)))

// Populate dst_x0x1x2 with reference-metric coordinates corresponding to each spherical grid point.
#pragma omp parallel for
  for (int iphi = 0; iphi < Nphi; iphi++) {
    const REAL phi = -M_PI + ((REAL)iphi + 0.5) * dphi;
    const REAL sinphi = sin(phi);
    const REAL cosphi = cos(phi);
    for (int itheta = 0; itheta < Ntheta; itheta++) {
      const REAL theta = ((REAL)itheta + 0.5) * dtheta;
      const REAL sintheta = sin(theta);
      const REAL costheta = cos(theta);
      for (int ir = 0; ir < bhahaha_params->Nr_external_input; ir++) {
        const REAL r = radii[ir];

        const int idx3 = IDX3_SPH_INTERP(ir, itheta, iphi);
        const REAL xCart[3] = {
            x_center + r * sintheta * cosphi, // x-coordinate of spherical grid point.
            y_center + r * sintheta * sinphi, // y-coordinate of spherical grid point.
            z_center + r * costheta           // z-coordinate of spherical grid point.
        };
        int Cart_to_i0i1i2_not_stored_to_save_memory[3];

        // Convert Cartesian coordinates to local grid reference-metric coordinates and store them.
        Cart_to_xx_and_nearest_i0i1i2(params, xCart, dst_x0x1x2[idx3], Cart_to_i0i1i2_not_stored_to_save_memory);
      } // END LOOP: Radial coordinate (ir).
    } // END LOOP: Theta coordinate (itheta).
  } // END LOOP: Phi coordinate (iphi).

  // Initialize source gridfunction pointers based on predefined interpolation indices.
  const REAL *restrict src_gf_ptrs[BHAHAHA_NUM_INTERP_GFS];
  {
    const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

    for (int idx = 0; idx < BHAHAHA_NUM_INTERP_GFS; idx++)
      src_gf_ptrs[idx] = &y_n_gfs[IDX4pt(bhahaha_gf_interp_indices[idx], 0)];
  }

  // Allocate memory for destination gridfunctions after interpolation.
  REAL *dst_data_ptrs[BHAHAHA_NUM_INTERP_GFS];
  for (int i = 0; i < BHAHAHA_NUM_INTERP_GFS; i++)
    dst_data_ptrs[i] = (REAL *)malloc(total_elements * sizeof(REAL));

  // Perform interpolation from source gridfunctions to spherical grid points.
  interpolation_3d_general__uniform_src_grid((NGHOSTS - 1), params->dxx0, params->dxx1, params->dxx2, params->Nxx_plus_2NGHOSTS0,
                                             params->Nxx_plus_2NGHOSTS1, params->Nxx_plus_2NGHOSTS2, BHAHAHA_NUM_INTERP_GFS, xx, src_gf_ptrs,
                                             total_elements, dst_x0x1x2, dst_data_ptrs);

  // Execute basis transformation from the reference metric (rfm) to the Cartesian basis and convert variables from BSSN to ADM formalism.
  {
#include "set_CodeParameters.h" // Incorporate additional code parameters required for transformation.

// Perform the basis transformation and variable conversion at each grid point.
#pragma omp parallel for
    for (int iphi = 0; iphi < Nphi; iphi++)
      for (int itheta = 0; itheta < Ntheta; itheta++)
        for (int ir = 0; ir < bhahaha_params->Nr_external_input; ir++) {
          const int idx3 = IDX3_SPH_INTERP(ir, itheta, iphi);
          const REAL xx0 = dst_x0x1x2[idx3][0];
          const REAL xx1 = dst_x0x1x2[idx3][1];
          const REAL xx2 = dst_x0x1x2[idx3][2];
          const REAL cf = dst_data_ptrs[INTERP_CFGF][idx3];
          const REAL trK = dst_data_ptrs[INTERP_TRKGF][idx3];
"""
    defines_list: List[str] = []
    for i in range(3):
        for j in range(i, 3):
            defines_list += [
                f"const REAL rfm_hDD{i}{j} = dst_data_ptrs[INTERP_HDD{i}{j}GF][idx3];\n"
            ]
            defines_list += [
                f"const REAL rfm_aDD{i}{j} = dst_data_ptrs[INTERP_ADD{i}{j}GF][idx3];\n"
            ]
    prefunc += "".join(sorted(defines_list, key=str.casefold))

    rfm = refmetric.reference_metric[CoordSystem]
    rfm_aDD = ixp.declarerank2("rfm_aDD", symmetry="sym01")
    rfm_hDD = ixp.declarerank2("rfm_hDD", symmetry="sym01")
    rfm_gammabarDD = ixp.zerorank2()
    rfm_AbarDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            rfm_gammabarDD[i][j] = rfm_hDD[i][j] * rfm.ReDD[i][j] + rfm.ghatDD[i][j]
            rfm_AbarDD[i][j] = rfm_aDD[i][j] * rfm.ReDD[i][j]
    Cart_gammabarDD = jac.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(
        CoordSystem, rfm_gammabarDD
    )
    Cart_AbarDD = jac.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(
        CoordSystem, rfm_AbarDD
    )
    exp4phi = sp.sympify(0)
    cf = sp.Symbol("cf", real=True)
    EvolvedConformalFactor_cf = par.parval_from_str("EvolvedConformalFactor_cf")
    if EvolvedConformalFactor_cf == "phi":
        exp4phi = sp.exp(4 * cf)
    elif EvolvedConformalFactor_cf == "chi":
        exp4phi = 1 / cf
    elif EvolvedConformalFactor_cf == "W":
        exp4phi = 1 / cf**2
    else:
        raise ValueError(
            f"Error EvolvedConformalFactor_cf type = {EvolvedConformalFactor_cf} unknown."
        )
    expr_list: List[sp.Expr] = []
    name_list: List[str] = []
    Cart_gammaDD = ixp.zerorank2()
    labels = ["X", "Y", "Z"]
    for i in range(3):
        for j in range(i, 3):
            Cart_gammaDD[i][j] = exp4phi * Cart_gammabarDD[i][j]
            expr_list += [Cart_gammaDD[i][j]]
            name_list += [
                f"dst_data_ptrs[INTERP_GAMMADD{labels[i]}{labels[j]}GF][idx3]"
            ]
    Cart_KDD = ixp.zerorank2()
    trK = sp.Symbol("trK", real=True)
    for i in range(3):
        for j in range(i, 3):
            Cart_KDD[i][j] = (
                exp4phi * Cart_AbarDD[i][j]
                + sp.Rational(1, 3) * Cart_gammaDD[i][j] * trK
            )
            expr_list += [Cart_KDD[i][j]]
            name_list += [f"dst_data_ptrs[INTERP_KDD{labels[i]}{labels[j]}GF][idx3]"]
    prefunc += ccg.c_codegen(expr_list, name_list, include_braces=True, verbose=False)

    prefunc += r"""
        } // END LOOP: Over all interpolated grid points.
  }

// Assign interpolated and transformed data to the output metric data array.
#pragma omp parallel for
  for (int gf = 0; gf < 12; gf++)
    for (int iphi = 0; iphi < Nphi; iphi++)
      for (int itheta = 0; itheta < Ntheta; itheta++) {
        for (int ir = 0; ir < bhahaha_params->Nr_external_input; ir++)
          input_metric_data[gf * total_elements + IDX3_SPH_INTERP(ir, itheta, iphi)] = dst_data_ptrs[gf][IDX3_SPH_INTERP(ir, itheta, iphi)];
      } // END LOOP: Assign transformed metric data to output array.

  // Release allocated memory for destination coordinates to prevent memory leaks.
  free(dst_x0x1x2);

  // Free memory allocated for destination gridfunctions after transferring data.
  for (int gf = 0; gf < BHAHAHA_NUM_INTERP_GFS; gf++)
    free(dst_data_ptrs[gf]);
} // END FUNCTION: BHaHAHA_interpolate_metric_data
"""

    desc = r"""Finds black hole apparent horizons using BHaHAHA, the BlackHoles@Home Apparent Horizon Algorithm.

This function iterates through each horizon to:
1. Initialize horizon parameters and diagnostic structures on the first iteration.
2. Perform interpolation of metric data onto a spherical grid.
3. Execute the BHaHAHA horizon finding algorithm.
4. Output diagnostic information upon successful horizon detection.
5. Adjust search parameters and retry if horizon finding fails.

@param commondata Pointer to the common data structure containing simulation parameters and state.
@param griddata   Pointer to the grid data structure containing grid parameters and gridfunctions.

@note Updates the commondata structure with found horizons and diagnostic information.
"""

    cfunc_type = "void"
    name = "bhahaha_find_horizons"
    params = (
        """commondata_struct *restrict commondata, griddata_struct *restrict griddata"""
    )

    body = r"""
  // Validate the number of horizons to ensure it is within acceptable bounds.
  if (commondata->bah_max_num_horizons <= 0 || commondata->bah_max_num_horizons > 100) {
    fprintf(stderr, "Invalid 'max_num_horizons': %d. Must be between 1 and 100.\n", commondata->bah_max_num_horizons);
    return;
  } // END IF: Validate 'max_num_horizons'.

  if (commondata->nn == 0) {
    // Initialize BHaHAHA parameters for all horizons on the first iteration.
    initialize_bhahaha_params_and_data_struct(commondata);
    if (commondata->bah_BBH_mode_enable) {
      if (commondata->bah_max_num_horizons != 3) {
        fprintf(stderr, "Error: bah_BBH_mode_enable requires bah_max_num_horizons==3\n");
        exit(1);
      }
      // At iteration 0, activate inspiralling BHs when BBH mode enabled;
      //   deactivate common horizon.
      const int inspiral_bh1 = commondata->bah_BBH_mode_inspiral_BH_idxs[0];
      const int inspiral_bh2 = commondata->bah_BBH_mode_inspiral_BH_idxs[1];
      commondata->bah_BBH_mode_horizon_active[inspiral_bh1] = commondata->bah_BBH_mode_horizon_active[inspiral_bh2] = 1;
      commondata->bah_BBH_mode_horizon_active[commondata->bah_BBH_mode_common_horizon_idx] = 0;
    } // END IF bah_BBH_mode_enable
  } // END IF: Initialization on first iteration.

  const REAL currtime = commondata->time;                     // Current simulation time.
  const REAL currdt = commondata->dt;                         // Current time step size.
  const REAL outevery = commondata->diagnostics_output_every; // Interval for diagnostics output.

  // Determine if the current time aligns with a diagnostics output interval.
  // This ensures 'currtime' is within half a time step of a multiple of 'outevery'.
  if (fabs(round(currtime / outevery) * outevery - currtime) < 0.5 * currdt) {

    // Iterate over each horizon to perform initialization and horizon finding.
    for (int which_horizon = 0; which_horizon < commondata->bah_max_num_horizons; which_horizon++) {
      // In BBH mode, skip horizons that aren't active.
      if (commondata->bah_BBH_mode_enable && !commondata->bah_BBH_mode_horizon_active[which_horizon])
        continue;

      // Set parameters for the current horizon.
      bhahaha_params_and_data_struct *restrict bhahaha_params_and_data = &commondata->bhahaha_params_and_data[which_horizon];

      // which_horizon & num_horizons are used for printing output within BHaHAHA (verbosity_level>0).
      bhahaha_params_and_data->which_horizon = which_horizon + 1;
      bhahaha_params_and_data->num_horizons = commondata->bah_max_num_horizons;
      bhahaha_params_and_data->iteration_external_input = commondata->nn;
      bhahaha_params_and_data->time_external_input = commondata->time;

      // Predict {x,y,z}_center, r_{min,max}, and then set radii[] for BHaHAHA's "external_input" spherical grid.
      REAL radii[commondata->bah_Nr_interp_max];
      REAL x_center, y_center, z_center, r_min, r_max;
      {
        // Predict {x,y,z}_center and r_{min,max} at the current time using previous data.
        bah_xyz_center_r_minmax(bhahaha_params_and_data, &x_center, &y_center, &z_center, &r_min, &r_max);
        if (bhahaha_params_and_data->use_fixed_radius_guess_on_full_sphere) {
          r_min = 0.0;
          r_max = commondata->bah_max_search_radius[which_horizon];
        }

        // Set radii[], a cell-centered radial grid for interpolation.
        bah_radial_grid_cell_centered_set_up(commondata->bah_Nr_interp_max, commondata->bah_max_search_radius[which_horizon], r_min, r_max,
                                             &bhahaha_params_and_data->Nr_external_input, &bhahaha_params_and_data->r_min_external_input,
                                             &bhahaha_params_and_data->dr_external_input, radii);
      } // END predict {x,y,z}_center, r_{min,max}, and then set radii[]

      // Allocate memory for input metric data, accounting for gridfunctions and radial points.
      const int Ntheta_max = commondata->bah_Ntheta_array_multigrid[commondata->bah_num_resolutions_multigrid - 1];
      const int Nphi_max = commondata->bah_Nphi_array_multigrid[commondata->bah_num_resolutions_multigrid - 1];
      REAL *restrict input_metric_data =
          malloc(sizeof(REAL) * BHAHAHA_NUM_INTERP_GFS * bhahaha_params_and_data->Nr_external_input * Ntheta_max * Nphi_max);
      if (input_metric_data == NULL) {
        fprintf(stderr, "Error: Memory allocation for input_metric_data failed.\n");
        exit(EXIT_FAILURE);
      }

      // Interpolate metric data from the NR code's grid(s) to the spherical BHaHAHA "external_input" grid.
      BHaHAHA_interpolate_metric_data(commondata, &griddata[0].params, griddata[0].xx, griddata[0].gridfuncs.y_n_gfs, bhahaha_params_and_data,
                                      x_center, y_center, z_center, radii, input_metric_data);
      bhahaha_params_and_data->input_metric_data = input_metric_data;

      // Initialize diagnostics structure to capture horizon finding results.
      bhahaha_diagnostics_struct *restrict bhahaha_diags = &commondata->bhahaha_diagnostics[which_horizon];

      printf("Horizon %d / %d : (x,y,z)_{relative to evolution grid origin} = (%.3f, %.3f, %.3f) | Nr_interp = %d\n", which_horizon + 1,
             commondata->bah_max_num_horizons, x_center, y_center, z_center, bhahaha_params_and_data->Nr_external_input);

      // Find horizons with BHaHAHA, the BlackHoles@Home Apparent Horizon Algorithm.
      const int return_code = bah_find_horizon(bhahaha_params_and_data, bhahaha_diags);

      // Release allocated memory for interpolated metric data, as it's not needed anymore.
      free(input_metric_data);

      if (return_code == BHAHAHA_SUCCESS) {
        // Output diagnostic information for the successfully found horizon. BHaHAHA input {x,y,z}_center
        //   is needed for 3D output into h.t*.gp, as h is with respect to this location.
        bah_diagnostics_file_output(bhahaha_diags, bhahaha_params_and_data, commondata->bah_max_num_horizons, x_center, y_center, z_center, ".");

        // Next time use this horizon as input into the next initial guess.
        bhahaha_params_and_data->use_fixed_radius_guess_on_full_sphere = 0;
      } else {
        // Horizon finding failed. Log the error and adjust search parameters.
        fprintf(stderr, "Horizon finding failed for horizon %d with error code %d:\n  %s\n", which_horizon, return_code,
                bah_error_message((bhahaha_error_codes)return_code));
        bhahaha_params_and_data->use_fixed_radius_guess_on_full_sphere = 1;
      } // END IF: Check horizon finding success.
    } // END LOOP: Over each horizon.

    if (commondata->bah_BBH_mode_enable) {
      const int inspiral_bh1 = commondata->bah_BBH_mode_inspiral_BH_idxs[0];
      const int inspiral_bh2 = commondata->bah_BBH_mode_inspiral_BH_idxs[1];
      const int common_horizon = commondata->bah_BBH_mode_common_horizon_idx;
      const bool found_bh1 = (commondata->bhahaha_params_and_data[inspiral_bh1].use_fixed_radius_guess_on_full_sphere == 0);
      const bool found_bh2 = (commondata->bhahaha_params_and_data[inspiral_bh2].use_fixed_radius_guess_on_full_sphere == 0);
      const bool found_common_horizon = (commondata->bhahaha_params_and_data[common_horizon].use_fixed_radius_guess_on_full_sphere == 0);
      const bool common_horizon_active = commondata->bah_BBH_mode_horizon_active[common_horizon];
      const bool bh1_active = (commondata->bah_BBH_mode_horizon_active[inspiral_bh1]);
      const bool bh2_active = (commondata->bah_BBH_mode_horizon_active[inspiral_bh2]);

      // BBH mode: Deactivate individual horizons if they are active & common horizon has been found.
      if (common_horizon_active && found_common_horizon && bh1_active && bh2_active) {
        commondata->bah_BBH_mode_horizon_active[inspiral_bh1] = commondata->bah_BBH_mode_horizon_active[inspiral_bh2] = 0;
        if (commondata->bhahaha_params_and_data->verbosity_level > 0)
          fprintf(stderr, "BHaHAHA: Deactivating individual horizon search.\n");
      } // END IF: Deactivate individual horizons if common horizon found.

      // BBH mode: Activate common horizon search if individual horizons close enough
      if (found_bh1 && found_bh2 && !common_horizon_active) {
        const REAL x_bh1 = commondata->bhahaha_params_and_data[inspiral_bh1].x_center_m1;
        const REAL x_bh2 = commondata->bhahaha_params_and_data[inspiral_bh2].x_center_m1;
        const REAL y_bh1 = commondata->bhahaha_params_and_data[inspiral_bh1].y_center_m1;
        const REAL y_bh2 = commondata->bhahaha_params_and_data[inspiral_bh2].y_center_m1;
        const REAL z_bh1 = commondata->bhahaha_params_and_data[inspiral_bh1].z_center_m1;
        const REAL z_bh2 = commondata->bhahaha_params_and_data[inspiral_bh2].z_center_m1;
        const REAL rmax_bh1 = commondata->bhahaha_params_and_data[inspiral_bh1].r_max_m1;
        const REAL rmax_bh2 = commondata->bhahaha_params_and_data[inspiral_bh2].r_max_m1;
        const REAL dist_between_horizons = sqrt(SQR(x_bh1 - x_bh2) + SQR(y_bh1 - y_bh2) + SQR(z_bh1 - z_bh2));
        if (dist_between_horizons <= 1.5 * (rmax_bh1 + rmax_bh2)) {
          const REAL M_scale_bh1 = commondata->bhahaha_params_and_data[inspiral_bh1].M_scale;
          const REAL M_scale_bh2 = commondata->bhahaha_params_and_data[inspiral_bh2].M_scale;
          commondata->bah_BBH_mode_horizon_active[common_horizon] = 1;
          // Set the common horizon centroid guess to be the center of mass.
          commondata->bhahaha_params_and_data[common_horizon].x_center_m1 = (M_scale_bh1 * x_bh1 + M_scale_bh2 * x_bh2) / (M_scale_bh1 + M_scale_bh2);
          commondata->bhahaha_params_and_data[common_horizon].y_center_m1 = (M_scale_bh1 * y_bh1 + M_scale_bh2 * y_bh2) / (M_scale_bh1 + M_scale_bh2);
          commondata->bhahaha_params_and_data[common_horizon].z_center_m1 = (M_scale_bh1 * z_bh1 + M_scale_bh2 * z_bh2) / (M_scale_bh1 + M_scale_bh2);

          if (commondata->bhahaha_params_and_data->verbosity_level > 0)
            fprintf(stderr, "BHaHAHA: Activating common horizon search.\n");
        } // END IF: Individual horizons close enough to warrant activating search for common horizon.
      } // END IF: Activate common horizon search
    } // END IF: BBH mode enabled

  } // END IF: It's time to search for horizons.
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
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
