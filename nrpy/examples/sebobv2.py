"""
Set up a complete C code project for the second version of the Spinning Effective-to-Backwards One Body (SEBOBv2) model.

This is a work in progress. The final SEBOBv2 contains
the following improvements over the current SEBOBv1:

[]1. Higher multipolar modes
[*]    1.a. Inspiral modes (2,2), (2,1), (3,3), (3,2), (4,4), (4,3), (5,5); Currently computes all inspiral modes but only outputs (2,2) for now
[]    1.b. Merger-ringdown modes (2,2), (2,1), (3,3), (3,2), (4,4), (4,3), (5,5)
[]    1.c. NQCs for higher modes
[]    1.d. Computing mode-mixing coefficients for higher modes
[x]    1.e. Transformation from co-precessing to inertial frame
[]2. Improved inspiral dynamics, including precession
[x]    2.a. Decoupled PN spin evolution equations
[]    2.b. Decoupled Quasi-precessing EOB Hamiltonian
[]    2.c. Generic spin EOB Hamiltonian
[x]    2.d. Post-adiabatic dynamics
[]3. Improved BOBv2 merger-ringdown
[x]    3.a. News-to-strain conversion
[]    3.b. Improved handling of precessing merger modes
[]    3.c. Analytical multipolar strain peak times
[]    3.d. Memory modes

By default, this example calculates aligned-spin (2,2) IMR modes using
SEOBNRv5 and BOBv2. Setting ``use_projected_attachment = True`` in the
parameter file enables projected-spin attachment, with ``chi1_z`` and
``chi2_z`` as the authoritative initial spin projections. Projected attachment
does not accept scalar ``chi1`` and ``chi2`` command-line overrides.

Authors:
        Anuj Kankani
        aak00009 **at** mix **dot** wvu **dot** edu
        Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Suchindram Dasgupta
        sd00113 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import os
import pkgutil

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import shutil
from pathlib import Path

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.infrastructures import BHaH

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "sebobv2"
use_projected_attachment_default = False

enable_parallel_codegen = True

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("enable_parallel_codegen", enable_parallel_codegen)

# Development flags (NOT command-line-tunable).
# Output the commondata struct to commondata.bin for low-level debugging.
output_commondata_flag = False
# Print the final IMR waveform to stdout. Set to False for performance checks.
output_waveform_flag = True
# Enable the optional inspiral-only coprecessing-rotation sandbox.
validate_sandbox_flag = False
# Write sandbox diagnostics. This only has an effect when validate_sandbox_flag is also True.
enable_sandbox_diagnostics_flag = False

#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]


def register_CFunction_main_c(
    output_waveform: bool = True,
    output_commondata: bool = True,
    validate_sandbox: bool = False,
    enable_sandbox_diagnostics: bool = False,
) -> None:
    """
    Generate a simplified C main() function for computing the SEBOB waveform.

    :param output_waveform: Flag to enable/disable printing the waveform
    :param output_commondata: Flag to enable/disable outputting commondata to a binary file.
    :param validate_sandbox: Whether to run the isolated coprecessing-rotation sandbox.
    :param enable_sandbox_diagnostics: Whether to write sandbox diagnostic sidecar files.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """-={ main() function }=-
Main function for computing the SEBOBv2 waveform.
"""
    cfunc_type = "int"
    name = "main"
    params = "int argc, const char *argv[]"
    body = r"""  commondata_struct commondata; // commondata contains parameters common to all grids.
// Step 1: Set commondata parameters to defaults.
commondata_struct_set_to_default(&commondata);
// Step 2: Apply parameter-file and command-line inputs.
cmdline_input_and_parfile_parser(&commondata, argc, argv);
if (commondata.use_projected_attachment) {
  if (argc > 2) {
    fprintf(stderr,
        "Error: projected attachment requires vector spin inputs from a parameter file; "
        "scalar chi1 and chi2 command-line overrides are not supported.\n");
    return EXIT_FAILURE;
  } // END IF: projected mode has CLI inputs
  commondata.chi1 = commondata.chi1_z;
  commondata.chi2 = commondata.chi2_z;
} // END IF: projected attachment input normalization
"""

    if validate_sandbox:
        body += r"""
if (fabs(commondata.chi1 - commondata.chi1_z) > 1e-14 ||
    fabs(commondata.chi2 - commondata.chi2_z) > 1e-14) {
  fprintf(stderr,
      "Warning: the optional coprecessing-rotation sandbox uses vector spin "
      "inputs while aligned attachment uses scalar chi1 and chi2. For "
      "self-consistent aligned sandbox validation, set chi1=chi1_z and "
      "chi2=chi2_z.\n");
  fflush(stderr);
} // END IF: sandbox spin representations differ
"""

    body += r"""
// Step 3: Initialize masses and model coefficients.
SEOBNRv5_quasi_precessing_spin_coefficients(&commondata);
// Step 4: Evolve the spin dynamics.
SEOBNRv5_quasi_precessing_spin_dynamics(&commondata);
// Step 5: Compute conservative initial conditions.
SEOBNRv5_aligned_spin_initial_conditions_conservative(&commondata);
// Step 6: Generate the inspiral trajectory.
SEOBNRv5_aligned_spin_pa_integration(&commondata);
// Step 7: Calculate special amplitude coefficients.
SEOBNRv5_aligned_spin_special_coefficients(&commondata);
// Step 8: Generate the inspiral waveform.
SEOBNRv5_aligned_spin_waveform_from_dynamics(&commondata);
"""

    if validate_sandbox:
        sandbox_diagnostics_flag = "1" if enable_sandbox_diagnostics else "0"
        body += r"""
// Optional validation sandbox: inspiral-only coprecessing rotations.
// This block exercises the J->P Euler-angle construction of Sec. III A,
// Eqs. (15)-(17), and the direct polarization rotation of Sec. III C,
// Eqs. (26)-(28), for fixed observer angles. It does not implement the
// merger-ringdown angle extension of Eqs. (18)-(20), and it does not alter
// the IMR waveform pipeline.
{
  const size_t n_low = commondata.nsteps_low;
  const size_t n_fine = commondata.nsteps_fine;
  const size_t n_insp = n_low + n_fine;
  const int enable_sandbox_diagnostics = """
        body += sandbox_diagnostics_flag
        body += r""";
  const REAL sandbox_iota = 0.9;
  const REAL sandbox_varphi_0 = 0.3;

  REAL *real_buffers = (REAL *)calloc(2 * n_insp, sizeof(REAL));
  COMPLEX *complex_buffers = (COMPLEX *)calloc(7 * n_insp, sizeof(COMPLEX));

  if (real_buffers != NULL && complex_buffers != NULL) {
    // Partition sandbox work buffers.
    REAL *h_plus_I  = real_buffers + 0 * n_insp;
    REAL *h_cross_I = real_buffers + 1 * n_insp;

    COMPLEX *hP_22 = complex_buffers + 0 * n_insp;
    COMPLEX *hP_21 = complex_buffers + 1 * n_insp;
    COMPLEX *hP_33 = complex_buffers + 2 * n_insp;
    COMPLEX *hP_32 = complex_buffers + 3 * n_insp;
    COMPLEX *hP_44 = complex_buffers + 4 * n_insp;
    COMPLEX *hP_43 = complex_buffers + 5 * n_insp;
    COMPLEX *hP_55 = complex_buffers + 6 * n_insp;

    // Splice and unpack all seven positive-m coprecessing inspiral modes.
    for (size_t i = 0; i < n_low; i++) {
      hP_22[i] = commondata.waveform_low[IDX_WF(i, STRAIN22)];
      hP_21[i] = commondata.waveform_low[IDX_WF(i, STRAIN21)];
      hP_33[i] = commondata.waveform_low[IDX_WF(i, STRAIN33)];
      hP_32[i] = commondata.waveform_low[IDX_WF(i, STRAIN32)];
      hP_44[i] = commondata.waveform_low[IDX_WF(i, STRAIN44)];
      hP_43[i] = commondata.waveform_low[IDX_WF(i, STRAIN43)];
      hP_55[i] = commondata.waveform_low[IDX_WF(i, STRAIN55)];
    } // END LOOP: low-frequency inspiral modes
    for (size_t i = 0; i < n_fine; i++) {
      size_t dest_idx = i + n_low;
      hP_22[dest_idx] = commondata.waveform_fine[IDX_WF(i, STRAIN22)];
      hP_21[dest_idx] = commondata.waveform_fine[IDX_WF(i, STRAIN21)];
      hP_33[dest_idx] = commondata.waveform_fine[IDX_WF(i, STRAIN33)];
      hP_32[dest_idx] = commondata.waveform_fine[IDX_WF(i, STRAIN32)];
      hP_44[dest_idx] = commondata.waveform_fine[IDX_WF(i, STRAIN44)];
      hP_43[dest_idx] = commondata.waveform_fine[IDX_WF(i, STRAIN43)];
      hP_55[dest_idx] = commondata.waveform_fine[IDX_WF(i, STRAIN55)];
    } // END LOOP: fine-frequency inspiral modes

    // Generate physical J->P Euler angles from the precessing dynamics.
    SEOBNRv5_coprecessing_angles(&commondata);

    // Execute the physical coprecessing-to-observer rotation.
    SEOBNRv5_coprecessing_rotations(
        n_insp, sandbox_iota, sandbox_varphi_0,
        commondata.J_f_x, commondata.J_f_y, commondata.J_f_z,
        commondata.alpha_JP, commondata.beta_JP, commondata.gamma_JP,
        hP_22, hP_21, hP_33, hP_32, hP_44, hP_43, hP_55,
        h_plus_I, h_cross_I
    );

    // Optionally dump raw modes, physical Euler angles, and rotated polarizations.
    FILE *fp = enable_sandbox_diagnostics ? fopen("validation_waveform.txt", "w") : NULL;
    if (fp != NULL) {
      fprintf(fp, "# Physical coprecessing-rotation validation: iota=%.15e varphi_0=%.15e J_f=(%.15e, %.15e, %.15e)\n",
          sandbox_iota, sandbox_varphi_0, commondata.J_f_x, commondata.J_f_y, commondata.J_f_z);
      fprintf(fp, "# Time | alpha_JP | beta_JP | gamma_JP | Re(hP_22) | Im(hP_22) | Re(hP_21) | Im(hP_21) | Re(hP_33) | Im(hP_33) | Re(hP_32) | Im(hP_32) | Re(hP_44) | Im(hP_44) | Re(hP_43) | Im(hP_43) | Re(hP_55) | Im(hP_55) | h_plus_I | h_cross_I\n");
      for (size_t i = 0; i < n_low; i++) {
        fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
            commondata.dynamics_low[IDX(i, TIME)],
            commondata.alpha_JP[i], commondata.beta_JP[i], commondata.gamma_JP[i],
            creal(hP_22[i]), cimag(hP_22[i]),
            creal(hP_21[i]), cimag(hP_21[i]),
            creal(hP_33[i]), cimag(hP_33[i]),
            creal(hP_32[i]), cimag(hP_32[i]),
            creal(hP_44[i]), cimag(hP_44[i]),
            creal(hP_43[i]), cimag(hP_43[i]),
            creal(hP_55[i]), cimag(hP_55[i]),
            h_plus_I[i], h_cross_I[i]);
      } // END LOOP: low-frequency validation samples
      for (size_t i = 0; i < n_fine; i++) {
        size_t dest_idx = i + n_low;
        fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
            commondata.dynamics_fine[IDX(i, TIME)],
            commondata.alpha_JP[dest_idx], commondata.beta_JP[dest_idx], commondata.gamma_JP[dest_idx],
            creal(hP_22[dest_idx]), cimag(hP_22[dest_idx]),
            creal(hP_21[dest_idx]), cimag(hP_21[dest_idx]),
            creal(hP_33[dest_idx]), cimag(hP_33[dest_idx]),
            creal(hP_32[dest_idx]), cimag(hP_32[dest_idx]),
            creal(hP_44[dest_idx]), cimag(hP_44[dest_idx]),
            creal(hP_43[dest_idx]), cimag(hP_43[dest_idx]),
            creal(hP_55[dest_idx]), cimag(hP_55[dest_idx]),
            h_plus_I[dest_idx], h_cross_I[dest_idx]);
      } // END LOOP: fine-frequency validation samples
      fclose(fp);
    } // END IF: diagnostic file opened
    else if (enable_sandbox_diagnostics) {
      fprintf(stderr, "Warning: Could not open validation_waveform.txt for writing.\n");
      fflush(stderr);
    } // END ELSE IF: diagnostic file open failed
  } // END IF: sandbox buffers allocated
  else {
    if (enable_sandbox_diagnostics) {
      fprintf(stderr, "Warning: Validation memory allocation failed. Skipping sandbox diagnostics.\n");
      fflush(stderr);
    } // END IF: diagnostic allocation failure
  } // END ELSE: sandbox work-buffer allocation failed

  // Free sandbox work buffers immediately.
  free(real_buffers);
  free(complex_buffers);
} // END BLOCK: optional inspiral-only coprecessing-rotation sandbox
"""
    body += r"""
// Step 9: Compute and apply NQC corrections.
SEBOBv2_NQC_corrections(&commondata);
// Step 10: Compute the IMR waveform.
SEBOBv2_IMR_waveform(&commondata);
"""

    if output_waveform:
        body += r"""
// Step 11: Print the IMR waveform.
for (size_t i = 0; i < commondata.nsteps_IMR; i++) {
    printf("%.15e %.15e %.15e\n", creal(commondata.waveform_IMR[IDX_WF(i,TIME)])
    , creal(commondata.waveform_IMR[IDX_WF(i,STRAIN22)]), cimag(commondata.waveform_IMR[IDX_WF(i,STRAIN22)]));
} // END LOOP: IMR waveform samples
"""

    if output_commondata:
        body += r"""
commondata_io(&commondata, "commondata.bin");
"""
    body += r"""
if (commondata.chi1_lnhat.spline != NULL) gsl_spline_free(commondata.chi1_lnhat.spline);
if (commondata.chi1_lnhat.acc != NULL) gsl_interp_accel_free(commondata.chi1_lnhat.acc);
if (commondata.chi2_lnhat.spline != NULL) gsl_spline_free(commondata.chi2_lnhat.spline);
if (commondata.chi2_lnhat.acc != NULL) gsl_interp_accel_free(commondata.chi2_lnhat.acc);
if (commondata.chi1_l.spline != NULL) gsl_spline_free(commondata.chi1_l.spline);
if (commondata.chi1_l.acc != NULL) gsl_interp_accel_free(commondata.chi1_l.acc);
if (commondata.chi2_l.spline != NULL) gsl_spline_free(commondata.chi2_l.spline);
if (commondata.chi2_l.acc != NULL) gsl_interp_accel_free(commondata.chi2_l.acc);
if (commondata.chi1_x_spline.spline != NULL) gsl_spline_free(commondata.chi1_x_spline.spline);
if (commondata.chi1_x_spline.acc != NULL) gsl_interp_accel_free(commondata.chi1_x_spline.acc);
if (commondata.chi1_y_spline.spline != NULL) gsl_spline_free(commondata.chi1_y_spline.spline);
if (commondata.chi1_y_spline.acc != NULL) gsl_interp_accel_free(commondata.chi1_y_spline.acc);
if (commondata.chi1_z_spline.spline != NULL) gsl_spline_free(commondata.chi1_z_spline.spline);
if (commondata.chi1_z_spline.acc != NULL) gsl_interp_accel_free(commondata.chi1_z_spline.acc);
if (commondata.chi2_x_spline.spline != NULL) gsl_spline_free(commondata.chi2_x_spline.spline);
if (commondata.chi2_x_spline.acc != NULL) gsl_interp_accel_free(commondata.chi2_x_spline.acc);
if (commondata.chi2_y_spline.spline != NULL) gsl_spline_free(commondata.chi2_y_spline.spline);
if (commondata.chi2_y_spline.acc != NULL) gsl_interp_accel_free(commondata.chi2_y_spline.acc);
if (commondata.chi2_z_spline.spline != NULL) gsl_spline_free(commondata.chi2_z_spline.spline);
if (commondata.chi2_z_spline.acc != NULL) gsl_interp_accel_free(commondata.chi2_z_spline.acc);
if (commondata.lnhat_x.spline != NULL) gsl_spline_free(commondata.lnhat_x.spline);
if (commondata.lnhat_x.acc != NULL) gsl_interp_accel_free(commondata.lnhat_x.acc);
if (commondata.lnhat_y.spline != NULL) gsl_spline_free(commondata.lnhat_y.spline);
if (commondata.lnhat_y.acc != NULL) gsl_interp_accel_free(commondata.lnhat_y.acc);
if (commondata.lnhat_z.spline != NULL) gsl_spline_free(commondata.lnhat_z.spline);
if (commondata.lnhat_z.acc != NULL) gsl_interp_accel_free(commondata.lnhat_z.acc);
if (commondata.L_x.spline != NULL) gsl_spline_free(commondata.L_x.spline);
if (commondata.L_x.acc != NULL) gsl_interp_accel_free(commondata.L_x.acc);
if (commondata.L_y.spline != NULL) gsl_spline_free(commondata.L_y.spline);
if (commondata.L_y.acc != NULL) gsl_interp_accel_free(commondata.L_y.acc);
if (commondata.L_z.spline != NULL) gsl_spline_free(commondata.L_z.spline);
if (commondata.L_z.acc != NULL) gsl_interp_accel_free(commondata.L_z.acc);
free(commondata.dynamics_low);
free(commondata.dynamics_fine);
free(commondata.dynamics_raw);
free(commondata.waveform_low);
free(commondata.waveform_fine);
free(commondata.waveform_inspiral);
free(commondata.waveform_IMR);
free(commondata.alpha_JP);
free(commondata.beta_JP);
free(commondata.gamma_JP);
return 0;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )


# register utilities needed by the waveform code
BHaH.seobnr.utils.commondata_io.register_CFunction_commondata_io()
BHaH.seobnr.utils.handle_gsl_return_status.register_CFunction_handle_gsl_return_status()
BHaH.seobnr.utils.SEOBNRv5_aligned_spin_unwrap.register_CFunction_SEOBNRv5_aligned_spin_unwrap()
BHaH.seobnr.utils.dy_dx.register_CFunction_dy_dx()
BHaH.seobnr.utils.finite_difference_stencil.register_CFunction_finite_difference_stencil()
BHaH.seobnr.utils.root_finding_1d.register_CFunction_root_finding_1d()
BHaH.seobnr.utils.root_finding_multidimensional.register_CFunction_root_finding_multidimensional()
BHaH.seobnr.utils.integration_stencil.register_CFunction_integration_stencil()
BHaH.seobnr.utils.cumulative_integration.register_CFunction_cumulative_integration()

# register SEOBNRv5 coefficients
BHaH.seobnr.SEOBNRv5_aligned_spin_coefficients.register_CFunction_SEOBNRv5_evaluate_l2m2_qnm()
BHaH.seobnr.SEOBNRv5_quasi_precessing_spin_coefficients.register_CFunction_SEOBNRv5_quasi_precessing_spin_coefficients(
    use_projected_attachment_default
)

# register h_NR fits
BHaH.seobnr.SEOBNRv5_aligned_spin_hNR_fits_at_t_attach.register_Cfunction_SEOBNRv5_aligned_spin_hNR_fits_at_t_attach()
BHaH.seobnr.SEOBNRv5_aligned_spin_omegaNR_fits_at_t_attach.register_Cfunction_SEOBNRv5_aligned_spin_omegaNR_fits_at_t_attach()

# register initial condition routines
BHaH.seobnr.initial_conditions.SEOBNRv5_aligned_spin_multidimensional_root_wrapper.register_CFunction_SEOBNRv5_multidimensional_root_wrapper()
BHaH.seobnr.initial_conditions.SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit.register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit()
BHaH.seobnr.initial_conditions.SEOBNRv5_aligned_spin_initial_conditions_conservative_nodf.register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_conservative_nodf()
BHaH.seobnr.initial_conditions.SEOBNRv5_aligned_spin_radial_momentum_condition.register_CFunction_SEOBNRv5_aligned_spin_radial_momentum_condition()
BHaH.seobnr.initial_conditions.SEOBNRv5_aligned_spin_initial_conditions_dissipative.register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_dissipative()

# register quasi-precessing spin initial conditions
BHaH.seobnr.dynamics.SEOBNRv5_quasi_precessing_spin_angular_momentum.register_CFunction_SEOBNRv5_quasi_precessing_spin_angular_momentum()
BHaH.seobnr.dynamics.SEOBNRv5_quasi_precessing_spin_dynamics.register_CFunction_SEOBNRv5_quasi_precessing_spin_dynamics()
BHaH.seobnr.dynamics.SEOBNRv5_quasi_precessing_spin_equations.register_CFunction_SEOBNRv5_quasi_precessing_spin_equations()

# register PA integration routines
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_pa_integration.register_CFunction_SEOBNRv5_aligned_spin_pa_integration()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_pphi0_equation.register_CFunction_SEOBNRv5_pphi0_equation()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_pphi_equation.register_CFunction_SEOBNRv5_pphi_equation()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_pr_equation.register_CFunction_SEOBNRv5_pr_equation()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_phi_equation.register_CFunction_SEOBNRv5_phi_equation()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_t_equation.register_CFunction_SEOBNRv5_t_equation()

# register trajectory integration and processing routines
BHaH.seobnr.dynamics.eval_abs_deriv.register_CFunction_eval_abs_deriv()
BHaH.seobnr.dynamics.find_local_minimum_index.register_CFunction_find_local_minimum_index()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_augments.register_CFunction_SEOBNRv5_aligned_spin_augments()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_interpolate_dynamics.register_CFunction_SEOBNRv5_aligned_spin_interpolate_dynamics()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_iterative_refinement.register_CFunction_SEOBNRv5_aligned_spin_iterative_refinement()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_right_hand_sides.register_CFunction_SEOBNRv5_aligned_spin_right_hand_sides()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_ode_integration.register_CFunction_SEOBNRv5_aligned_spin_ode_integration()

# register inspiral waveform routines
BHaH.seobnr.inspiral_waveform.SEOBNRv5_aligned_spin_gamma_wrapper.register_CFunction_SEOBNRv5_aligned_spin_gamma_wrapper()
BHaH.seobnr.inspiral_waveform.SEOBNRv5_aligned_spin_interpolate_modes.register_CFunction_SEOBNRv5_aligned_spin_interpolate_modes()
BHaH.seobnr.inspiral_waveform.SEOBNRv5_aligned_spin_waveform_from_dynamics_higher_mode.register_CFunction_SEOBNRv5_aligned_spin_waveform_from_dynamics()
BHaH.seobnr.inspiral_waveform.SEOBNRv5_aligned_spin_waveform_higher_mode.register_CFunction_SEOBNRv5_aligned_spin_waveform()
BHaH.seobnr.inspiral_waveform.SEOBNRv5_aligned_spin_special_amplitude_coefficients.register_Cfunction_SEOBNRv5_aligned_spin_special_amplitude_coefficients_rholm()
BHaH.seobnr.inspiral_waveform.SEOBNRv5_aligned_spin_special_amplitude_coefficients.register_Cfunction_SEOBNRv5_aligned_spin_special_amplitude_coefficients()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_flux.register_CFunction_SEOBNRv5_aligned_spin_flux()
BHaH.seobnr.SEOBNRv5_coprecessing_angles.register_CFunction_SEOBNRv5_coprecessing_angles(
    enable_forensic_diagnostics=enable_sandbox_diagnostics_flag,
)
BHaH.seobnr.SEOBNRv5_coprecessing_rotations.register_CFunction_SEOBNRv5_coprecessing_rotations()

# Register the SEBOBv2-specific BOB peak-time parameter.
par.register_CodeParameters(
    "REAL",
    __name__,
    ["t_p_BOB"],
    commondata=True,
    add_to_parfile=False,
    add_to_set_CodeParameters_h=False,
)

if __name__ == "__main__":
    # Register some functions/code parameters based on input flags
    BHaH.seobnr.nqc_corrections.SEBOBv2_NQC_corrections.register_CFunction_SEBOBv2_NQC_corrections()
    BHaH.seobnr.nqc_corrections.BOB_v2_NQC_rhs.register_CFunction_BOB_v2_NQC_rhs()
    BHaH.seobnr.nqc_corrections.BOB_aligned_spin_NQC_corrections_higher_modes.register_CFunction_BOB_aligned_spin_NQC_corrections_higher_modes(
        False
    )
    BHaH.seobnr.nqc_corrections.BOB_aligned_spin_NQC_rhs_higher_modes.register_CFunction_BOB_aligned_spin_NQC_rhs_HM()

    # set up merger-ringdown routines based on input flags
    BHaH.seobnr.merger_waveform.BOB_v2_setup_peak_attachment.register_CFunction_BOB_v2_setup_peak_attachment()
    BHaH.seobnr.merger_waveform.BOB_v2_waveform.register_CFunction_BOB_v2_waveform()
    BHaH.seobnr.merger_waveform.BOB_v2_waveform_from_times.register_CFunction_BOB_v2_waveform_from_times()
    # register IMR waveform generation routine
    BHaH.seobnr.SEBOBv2_IMR_waveform.register_CFunction_SEBOBv2_IMR_waveform()
    pcg.do_parallel_codegen()
#########################################################
# STEP 3: Generate header files, register C functions and
#         command line parameters, set up boundary conditions,
#         and create a Makefile for this project.
#         Project is output to project/[project_name]/
BHaH.CodeParameters.write_CodeParameters_h_files(
    set_commondata_only=True, project_dir=project_dir
)
BHaH.CodeParameters.register_CFunctions_params_commondata_struct_set_to_default()
BHaH.cmdline_input_and_parfiles.generate_default_parfile(
    project_dir=project_dir, project_name=project_name
)
BHaH.cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
    project_name=project_name,
    cmdline_inputs=["mass_ratio", "chi1", "chi2", "initial_omega", "total_mass", "dt"],
)
# Load the header file using pkgutil
data_bytes = pkgutil.get_data("nrpy.infrastructures.BHaH.seobnr", "spline_struct.h")
if data_bytes is None:
    raise FileNotFoundError("spline_struct.h not found via pkgutil.get_data")

spline_struct_h = data_bytes.decode("utf-8")

# Write the updated content to the output file
with Path(project_dir, "spline_struct.h").open("w", encoding="utf-8") as output_file:
    output_file.write(spline_struct_h)

additional_includes = [
    str(Path("gsl") / Path("gsl_vector.h")),
    str(Path("gsl") / Path("gsl_multiroots.h")),
    str(Path("gsl") / Path("gsl_errno.h")),
    str(Path("gsl") / Path("gsl_roots.h")),
    str(Path("gsl") / Path("gsl_matrix.h")),
    str(Path("gsl") / Path("gsl_odeiv2.h")),
    str(Path("gsl") / Path("gsl_spline.h")),
    str(Path("gsl") / Path("gsl_interp.h")),
    str(Path("gsl") / Path("gsl_sf_gamma.h")),
    str(Path("gsl") / Path("gsl_linalg.h")),
    "complex.h",
    "spline_struct.h",
]
BHaH.BHaH_defines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    additional_includes=additional_includes,
    enable_rfm_precompute=False,
    supplemental_defines_dict={"SEOBNR": """
#include<complex.h>
#define COMPLEX double complex
#define NUMVARS_SPIN 10
#define LN_X 0
#define LN_Y 1
#define LN_Z 2
#define CHI1_X 3
#define CHI1_Y 4
#define CHI1_Z 5
#define CHI2_X 6
#define CHI2_Y 7
#define CHI2_Z 8
#define OMEGA_PN 9
#define IDX_SPIN(idx,var) ((idx)*NUMVARS_SPIN + (var))
#define NUMVARS 8
#define TIME 0
#define R 1
#define PHI 2
#define PRSTAR 3
#define PPHI 4
#define H 5
#define OMEGA 6
#define OMEGA_CIRC 7
#define IDX(idx, var) ((idx)*NUMVARS + (var))
#define NUMVARS_COEFFICIENTS 3
#define RHO21 0
#define RHO43 1
#define RHO55 2
#define IDX_COEFFICIENTS(idx, var) ((idx)*NUMVARS_COEFFICIENTS + (var))
#define NUMVARS_HNRFITS 7
#define HNR22 0
#define HNR21 1
#define HNR33 2
#define HNR32 3
#define HNR43 4
#define HNR44 5
#define HNR55 6
#define IDX_HNRFITS(idx, var) ((idx)*NUMVARS_HNRFITS + (var))
#define NUMMODES 8
#define STRAIN22 1
#define STRAIN21 2
#define STRAIN33 3
#define STRAIN32 4
#define STRAIN44 5
#define STRAIN43 6
#define STRAIN55 7
#define NUMMODESSTORED 2 // process 2,2 mode for now
#define STRAIN 1
#define IDX_WF(idx,var) ((idx)*NUMMODES + (var))
"""},
)
register_CFunction_main_c(
    output_waveform_flag,
    output_commondata_flag,
    validate_sandbox_flag,
    enable_sandbox_diagnostics_flag,
)

addl_cflags = ["$(shell gsl-config --cflags)"]
BHaH.Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir,
    project_name=project_name,
    exec_or_library_name=project_name,
    addl_CFLAGS=addl_cflags,
    addl_libraries=["$(shell gsl-config --libs)"],
)

print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
)
print(f"    Parameter file can be found in {project_name}.par")
