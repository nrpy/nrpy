"""
Generate Psi4 decomposition on spherical like grids using spin-weighted spherical harmonics.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

import nrpy.c_function as cfc
import nrpy.params as par


def register_CFunction_psi4_spinweightm2_decomposition() -> None:
    """
    Register the C function `psi4_spinweightm2_decomposition`.

    The registered C function:
    - Sets parameters for psi4 extraction radii and angular resolution.
    - Interpolates psi4 from the 3D grid onto spherical shells.
    - Decomposes psi4 into spin-weight -2 spherical-harmonic modes.
    - Writes R_ext * psi4_{l,m} time series to disk.
    """
    par.register_CodeParameters(
        "int",
        __name__,
        [
            "num_theta_points_on_shell_for_psi4_interp",
            "num_phi_points_on_shell_for_psi4_interp",
        ],
        [32, 64],
        add_to_parfile=True,
        commondata=True,
    )
    par.register_CodeParameter(
        cparam_type="int",
        module=__name__,
        name="num_psi4_extraction_radii",
        defaultvalue=6,
        description="Number of radii at which psi4 is extracted.",
        add_to_parfile=True,
        commondata=True,
        add_to_set_CodeParameters_h=False,
    )

    # Register list_of_psi4_extraction_radii: the array containing extraction radii values.
    par.register_CodeParameter(
        cparam_type="REAL[6]",
        module=__name__,
        name="list_of_psi4_extraction_radii",
        defaultvalue=[15, 30, 45, 60, 75, 90],
        description="List of radii at which psi4 is extracted. Must set num_psi4_extraction_radii consistently.",
        add_to_parfile=True,
        commondata=True,
        add_to_set_CodeParameters_h=False,
        assumption="RealPositive",
    )

    prefunc = r"""
#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

// Needed for interpolation error codes (copied from interpolation_3d_general__uniform_src_grid.c)
MAYBE_UNUSED static enum {
  INTERP_SUCCESS,
  INTERP3D_GENERAL_NULL_PTRS,
  INTERP3D_GENERAL_INTERP_ORDER_GT_NXX123,
  INTERP3D_GENERAL_HORIZON_OUT_OF_BOUNDS
} error_codes;

#define IDX2GENERAL(i, j, Ni) ((i) + (Ni) * (j))

static void lowlevel_decompose_psi4_into_swm2_modes(const REAL dtheta, const REAL dphi, const int swm2sh_maximum_l_mode_to_compute,
                                                    const REAL curr_time, const REAL R_ext, const int N_theta, const int N_phi,
                                                    const REAL *restrict theta_array, const REAL *restrict phi_array,
                                                    const REAL *restrict psi4r_at_R_ext, const REAL *restrict psi4i_at_R_ext);

void psi4_spinweightm2_shell_init(const commondata_struct *restrict commondata, psi4_shell_angular_grid_t *restrict shell) {
  shell->N_theta = commondata->num_theta_points_on_shell_for_psi4_interp;
  shell->N_phi = commondata->num_phi_points_on_shell_for_psi4_interp;
  shell->num_pts = shell->N_theta * shell->N_phi;
  shell->dtheta = M_PI / shell->N_theta;
  shell->dphi = (2.0 * M_PI) / shell->N_phi;

  shell->theta_array = (REAL *)malloc(shell->N_theta * sizeof(REAL));
  shell->sin_theta_array = (REAL *)malloc(shell->N_theta * sizeof(REAL));
  shell->cos_theta_array = (REAL *)malloc(shell->N_theta * sizeof(REAL));
  shell->phi_array = (REAL *)malloc(shell->N_phi * sizeof(REAL));
  shell->sin_phi_array = (REAL *)malloc(shell->N_phi * sizeof(REAL));
  shell->cos_phi_array = (REAL *)malloc(shell->N_phi * sizeof(REAL));

  for (int i = 0; i < shell->N_theta; i++) {
    shell->theta_array[i] = ((REAL)i + 0.5) * shell->dtheta;
    shell->sin_theta_array[i] = sin(shell->theta_array[i]);
    shell->cos_theta_array[i] = cos(shell->theta_array[i]);
  }
  for (int i = 0; i < shell->N_phi; i++) {
    shell->phi_array[i] = -M_PI + ((REAL)i + 0.5) * shell->dphi;
    shell->sin_phi_array[i] = sin(shell->phi_array[i]);
    shell->cos_phi_array[i] = cos(shell->phi_array[i]);
  }
}

void psi4_spinweightm2_shell_free(psi4_shell_angular_grid_t *restrict shell) {
  free(shell->theta_array);
  free(shell->sin_theta_array);
  free(shell->cos_theta_array);
  free(shell->phi_array);
  free(shell->sin_phi_array);
  free(shell->cos_phi_array);
  shell->theta_array = NULL;
  shell->sin_theta_array = NULL;
  shell->cos_theta_array = NULL;
  shell->phi_array = NULL;
  shell->sin_phi_array = NULL;
  shell->cos_phi_array = NULL;
}

void psi4_spinweightm2_shell_fill_points(const params_struct *restrict params, const psi4_shell_angular_grid_t *restrict shell, const REAL R_ext,
                                         REAL (*dst_pts)[3], int *all_points_interior) {
  int all_interior = 1;
  for (int i_ph = 0; i_ph < shell->N_phi; i_ph++) {
    const REAL sin_phi = shell->sin_phi_array[i_ph];
    const REAL cos_phi = shell->cos_phi_array[i_ph];
    for (int i_th = 0; i_th < shell->N_theta; i_th++) {
      const REAL sin_theta = shell->sin_theta_array[i_th];
      const REAL cos_theta = shell->cos_theta_array[i_th];
      const int idx = IDX2GENERAL(i_th, i_ph, shell->N_theta);

      const REAL xCart[3] = {R_ext * sin_theta * cos_phi, R_ext * sin_theta * sin_phi, R_ext * cos_theta};
      int i0i1i2[3];
      REAL dst_xx[3];
      Cart_to_xx_and_nearest_i0i1i2(params, xCart, dst_xx, i0i1i2);

      dst_pts[idx][0] = dst_xx[0];
      dst_pts[idx][1] = dst_xx[1];
      dst_pts[idx][2] = dst_xx[2];

      if (all_points_interior && !IS_IN_GRID_INTERIOR(i0i1i2, params->Nxx_plus_2NGHOSTS0, params->Nxx_plus_2NGHOSTS1,
                                                      params->Nxx_plus_2NGHOSTS2, NGHOSTS)) {
        all_interior = 0;
      }
    }
  }
  if (all_points_interior) {
    *all_points_interior = all_interior;
  }
}

void psi4_spinweightm2_decompose_shell(const commondata_struct *restrict commondata, const psi4_shell_angular_grid_t *restrict shell, const REAL curr_time,
                                       const REAL R_ext, const REAL *restrict psi4r_at_R_ext, const REAL *restrict psi4i_at_R_ext) {
  lowlevel_decompose_psi4_into_swm2_modes(shell->dtheta, shell->dphi, commondata->swm2sh_maximum_l_mode_to_compute, curr_time, R_ext, shell->N_theta,
                                         shell->N_phi, shell->theta_array, shell->phi_array, psi4r_at_R_ext, psi4i_at_R_ext);
}

/**
 * @brief Decomposes psi4 on a spherical shell into spin-weight -2 spherical harmonic modes and appends them to output files.
 *
 * @details This function performs a numerical integral of psi4 multiplied by the complex conjugate of the spin-weight -2
 *          spherical harmonics to compute the mode coefficients psi4_{l,m}. The integration is performed over a
 *          spherical shell at a given extraction radius.
 *
 *          For each l-mode from 2 to swm2sh_maximum_l_mode_to_compute, this function appends one line of data to a
 *          file named "Rpsi4_l<l>-r<R_ext>.txt".
 *
 *          At the first call (curr_time == 0), it creates the files and writes a header describing the columns.
 *          On all calls, it appends a row containing the retarded time (t - R_ext) followed by the real and
 *          imaginary parts of (R_ext * psi4_{l,m}) for each m from -l to l.
 *
 *          The numerical integration is parallelized using OpenMP.
 *
 * @param[in] dtheta The angular spacing in the theta direction.
 * @param[in] dphi The angular spacing in the phi direction.
 * @param[in] swm2sh_maximum_l_mode_to_compute The maximum l-mode to compute for the decomposition.
 * @param[in] curr_time The current simulation time.
 * @param[in] R_ext The radius of the spherical shell where psi4 is sampled.
 * @param[in] N_theta The number of points in the theta direction on the shell.
 * @param[in] N_phi The number of points in the phi direction on the shell.
 * @param[in] theta_array Array of theta coordinates for the points on the shell.
 * @param[in] phi_array Array of phi coordinates for the points on the shell.
 * @param[in] psi4r_at_R_ext Pointer to an array containing the real part of psi4 on the shell.
 * @param[in] psi4i_at_R_ext Pointer to an array containing the imaginary part of psi4 on the shell.
 *
 * @return void
 */
static void lowlevel_decompose_psi4_into_swm2_modes(const REAL dtheta, const REAL dphi, const int swm2sh_maximum_l_mode_to_compute,
                                                    const REAL curr_time, const REAL R_ext, const int N_theta, const int N_phi,
                                                    const REAL *restrict theta_array, const REAL *restrict phi_array,
                                                    const REAL *restrict psi4r_at_R_ext, const REAL *restrict psi4i_at_R_ext) {
  char filename[100];
  // Initialize files at t=0
  if (curr_time == 0) {
    for (int l = 2; l <= swm2sh_maximum_l_mode_to_compute; l++) {
      sprintf(filename, "Rpsi4_l%d-r%06.1f.txt", l, (double)R_ext);
      FILE *outfile = fopen(filename, "w");
      if (!outfile) {
        perror("Error opening psi4 output file");
        continue;
      } // END IF file open error
      fprintf(outfile, "# column 1: t-R_ext = [retarded time]\n");
      int col = 2;
      for (int m = -l; m <= l; m++) {
        fprintf(outfile, "# column %d: Re(psi4_{l=%d,m=%d}) * R_ext\n", col++, l, m);
        fprintf(outfile, "# column %d: Im(psi4_{l=%d,m=%d}) * R_ext\n", col++, l, m);
      } // END LOOP over m modes
      fclose(outfile);
    } // END LOOP over l modes
  } // END IF t==0 header output

  // Perform decomposition
  for (int l = 2; l <= swm2sh_maximum_l_mode_to_compute; l++) {
    sprintf(filename, "Rpsi4_l%d-r%06.1f.txt", l, (double)R_ext);
    FILE *outfile = fopen(filename, "a");
    if (!outfile) {
      perror("Error opening psi4 output file for appending");
      continue;
    } // END IF file open error
    fprintf(outfile, "%e", (double)(curr_time - R_ext));
    for (int m = -l; m <= l; m++) {
      REAL psi4r_l_m = 0.0, psi4i_l_m = 0.0;
#pragma omp parallel for collapse(2) reduction(+ : psi4r_l_m, psi4i_l_m)
      for (int i_ph = 0; i_ph < N_phi; i_ph++) {
        for (int i_th = 0; i_th < N_theta; i_th++) {
          const REAL ph = phi_array[i_ph];
          const REAL th = theta_array[i_th];
          const int idx = IDX2GENERAL(i_th, i_ph, N_theta);
          REAL ReY, ImY;
          spin_weight_minus2_sph_harmonics(l, m, th, ph, &ReY, &ImY);
          const REAL psi4r = psi4r_at_R_ext[idx];
          const REAL psi4i = psi4i_at_R_ext[idx];
          const REAL dOmega = sin(th) * dtheta * dphi;
          // Integral: psi4 * conj(Y) * sin(theta) * dtheta * dphi
          psi4r_l_m += (psi4r * ReY + psi4i * ImY) * dOmega;
          psi4i_l_m += (psi4i * ReY - psi4r * ImY) * dOmega;
        } // END LOOP over theta
      } // END LOOP over phi
      fprintf(outfile, " %.15e %.15e", (double)(R_ext * psi4r_l_m), (double)(R_ext * psi4i_l_m));
    } // END LOOP over m modes
    fprintf(outfile, "\n");
    fclose(outfile);
  } // END LOOP over l modes
} // END FUNCTION lowlevel_decompose_psi4_into_swm2_modes()
"""

    desc = """
 * @brief Manages the extraction and decomposition of the Weyl scalar psi4 at multiple radii.
 *
 * @details This function orchestrates the process of computing the spin-weight -2 spherical harmonic modes of psi4.
 *          It iterates through a list of specified extraction radii. For each radius, it constructs a spherical
 *          shell of interpolation points.
 *
 *          A crucial check is performed to ensure that all points on the spherical shell lie within the grid's
 *          interior (i.e., not in the ghost zones). If any point is outside, the interpolation for that shell is
 *          skipped, and the resulting modes are effectively zero.
 *
 *          If the shell is fully within the interior, the function performs a batch interpolation of the real and
 *          imaginary components of psi4 from the 3D grid onto the shell.
 *
 *          Finally, it calls the low-level function `lowlevel_decompose_psi4_into_swm2_modes` to perform the
 *          numerical integration and write the mode coefficients to disk.
 *
 * @param[in] commondata Pointer to the commondata_struct, containing global simulation parameters like time,
 *                       and psi4 extraction settings (radii, number of angular points).
 * @param[in] griddata Pointer to the griddata_struct, containing grid-specific data like parameters and
 *                     coordinate arrays. This function uses data from grid 0.
 * @param[in] gridfuncs_diags Array of pointers to the diagnostic gridfunctions for each grid. This function
 *                            reads psi4 data from gridfuncs_diags[0].
 *
 * @return void
"""
    name = "psi4_spinweightm2_decomposition"
    params = r"""commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                                     const REAL *restrict gridfuncs_diags[MAXNUMGRIDS]"""
    body = r"""
  // FIXME: Extract grid 0 data only
  const params_struct *restrict params = &griddata[0].params;
  const REAL *restrict diagnostic_gfs = gridfuncs_diags[0];
  REAL *restrict xx[3] = {griddata[0].xx[0], griddata[0].xx[1], griddata[0].xx[2]};

  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  const int num_radii = commondata->num_psi4_extraction_radii;
  const REAL *restrict radii = commondata->list_of_psi4_extraction_radii;
  psi4_shell_angular_grid_t shell;
  psi4_spinweightm2_shell_init(commondata, &shell);

  const REAL *src_gf_ptrs[2] = {&diagnostic_gfs[IDX4(DIAG_PSI4_REGF, 0, 0, 0)], &diagnostic_gfs[IDX4(DIAG_PSI4_IMGF, 0, 0, 0)]};

  for (int which_R_ext = 0; which_R_ext < num_radii; which_R_ext++) {
    const REAL R_ext = radii[which_R_ext];

    // Allocate arrays for full shell
    REAL(*all_dst_pts)[3] = (REAL(*)[3])malloc(sizeof(REAL) * shell.num_pts * 3);
    REAL *psi4r_at_R_ext = (REAL *)calloc(shell.num_pts, sizeof(REAL));
    REAL *psi4i_at_R_ext = (REAL *)calloc(shell.num_pts, sizeof(REAL));

    int all_points_interior = 1;
    psi4_spinweightm2_shell_fill_points(params, &shell, R_ext, all_dst_pts, &all_points_interior);

    // Only interpolate if ALL points on the sphere are in the interior
    if (all_points_interior) {
      // Batch interpolate all points at once
      REAL *dst_data[2] = {psi4r_at_R_ext, psi4i_at_R_ext};
      int error_code = interpolation_3d_general__uniform_src_grid(NGHOSTS, params->dxx0, params->dxx1, params->dxx2, params->Nxx_plus_2NGHOSTS0,
                                                                  params->Nxx_plus_2NGHOSTS1, params->Nxx_plus_2NGHOSTS2, 2, xx, src_gf_ptrs, shell.num_pts,
                                                                  (const REAL(*)[3])all_dst_pts, dst_data);

      if (error_code != 0)
        fprintf(stderr, "Interpolation error code %d for R_ext = %e\n", error_code, (double)R_ext);
    } // END IF all_points_interior
    // If not all points in interior, arrays remain zero from calloc and all l,m modes will be zero

    // Perform decomposition (with zeros for all modes if sphere not entirely in interior)
    psi4_spinweightm2_decompose_shell(commondata, &shell, commondata->time, R_ext, psi4r_at_R_ext, psi4i_at_R_ext);

    free(all_dst_pts);
    free(psi4r_at_R_ext);
    free(psi4i_at_R_ext);
  } // END LOOP over extraction radii

  psi4_spinweightm2_shell_free(&shell);
"""

    cfc.register_CFunction(
        includes=[
            "BHaH_defines.h",
            "BHaH_function_prototypes.h",
            "diagnostics/diagnostic_gfs.h",
        ],
        prefunc=prefunc,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )
