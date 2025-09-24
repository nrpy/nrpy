"""
Set up a complete C code project for the second version of the Spinning Effective-to-Backwards One Body (SEBOBv2) model.

This is a work in progress. The final SEBOBv2 contains
the following improvements over the current SEBOBv1:

[]1. Higher multipolar modes
[]    1.a. Inspiral modes (2,2), (2,1), (3,3), (3,2), (4,4), (4,3), (5,5)
[]    1.b. Merger-ringdown modes (2,2), (2,1), (3,3), (3,2), (4,4), (4,3), (5,5)
[]    1.c. NQCs for higher modes
[]    1.d. Computing mode-mixing coefficients for higher modes
[]    1.e. Transformation from co-precessing to inertial frame
[]2. Improved inspiral dynamics, including precession
[]    2.a. Decoupled Quasi-precessing EOB Hamiltonian
[]    2.b. Decoupled PN spin evolution equations
[]    2.c. Generic spin EOB Hamiltonian
[]    2.d. Post-adiabatic dynamics
[]3. Improved BOBv2 merger-ringdown
[x]    3.a. News-to-strain conversion
[]    3.b. Improved handling of precessing merger modes
[]    3.c. Analytical multipolar strain peak times
[]    3.d. Memory modes

Currently, this examples calculates:
Aligned-spin (2,2) IMR modes using SEOBNRv5 and BOBv2.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import os

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

parallel_codegen_enable = True

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)

# Development flags (NOT command-line-tunable)
# Flag to compute the SEOBNR waveform in the frequency domain.
frequency_domain_flag = False
# Flag to precompute the waveform coefficients. Only works for aligned spins.
# Once we have precession in place, we can make this more tunable
precompute_waveform_coefficients_flag = False
# Flag to output the commondata struct to a file.
output_commondata_flag = False
# Flag to output the SEOBNRv5 waveform using a print statement like lalsimulation does.
# (set to False for performance checks)
output_waveform_flag = True
# SEOBNRv5 uses an iterative refinement routine to find the location of the peak of the orbital frequency.
# This is not done in a robust manner and disabling it does not impact accuracy.
# This flag helps enable/disable the routine at the codegen level so that we can assess the impact on performance and accuracy.
perform_iterative_refinement = True
# Removed all command line flags since we do not have multiple usage options currently.

#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]


def register_CFunction_main_c(
    output_waveform: bool = True,
    frequency_domain: bool = False,
    precompute_waveform_coefficients: bool = False,
    output_commondata: bool = False,
) -> None:
    """
    Generate a simplified C main() function for computing the SEBOB waveform.

    :param output_waveform: Flag to enable/disable printing the waveform
    :param frequency_domain: Flag to enable/disable FFT to get a frequency domain waveform
    :param precompute_waveform_coefficients: Flag to enable/disable precomputing the waveform coefficients
    :param output_commondata: Flag to enable/disable outputting the commondata struct to a binary file
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """-={ main() function }=-
Step 1.a: Set each commondata CodeParameter to default.
Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
Step 2: Compute Hamiltonian and derivatives."""
    cfunc_type = "int"
    name = "main"
    params = "int argc, const char *argv[]"
    body = r"""  commondata_struct commondata; // commondata contains parameters common to all grids.
//Step 0: Initialize a loop parameter for outputs
// Step 1.a: Set each commondata CodeParameter to default.
commondata_struct_set_to_default(&commondata);
// Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
cmdline_input_and_parfile_parser(&commondata, argc, argv);
// Step 1.c: Overwrite default values of m1, m2, a6, and dSO.
SEOBNRv5_aligned_spin_coefficients(&commondata);
"""
    if precompute_waveform_coefficients:
        body += r"""// Step 1.d: Set the waveform coefficients
SEOBNRv5_aligned_spin_waveform_coefficients(&commondata);
"""
    body += r"""// Step 2.a: Compute SEOBNRv5 conservative initial conditions.
SEOBNRv5_aligned_spin_initial_conditions_conservative(&commondata);
// Step 2.b: Print out the conservative initial conditions.
//printf("r = %.15e\n",commondata.r);
//printf("pphi = %.15e\n",commondata.pphi);
// Step 3.a: Compute SEOBNRv5 dissipative initial conditions.
SEOBNRv5_aligned_spin_initial_conditions_dissipative(&commondata);
// Step 3.b: Print out the dissipative initial conditions.
//printf("prstar = %.15e\n",commondata.prstar);
// Step 4: Run the ODE integration.
SEOBNRv5_aligned_spin_ode_integration(&commondata);
// Step 5. Generate the waveform.
SEOBNRv5_aligned_spin_waveform_from_dynamics(&commondata);
// Step 6. Compute and apply the NQC corrections
SEBOBv2_NQC_corrections(&commondata);
// Step 7.a Compute the IMR waveform
SEBOBv2_IMR_waveform(&commondata);
"""
    if frequency_domain:
        body += r"""
// Step 7.b Compute the FFT-ed IMR waveform
// Specify wisdom file
const char *wisdom_file = "fftw_wisdom.dat";
SEOBNRv5_aligned_spin_FD_waveform(wisdom_file, &commondata);
"""
    if output_waveform:
        if frequency_domain:
            body += r"""
// Step 6.b: Print the resulting waveform.
for (size_t i = 0; i < commondata.nsteps_IMR_FD; i++) {
    printf("%.15e %.15e %.15e\n", creal(commondata.waveform_IMR_FD[IDX_WF(i,FREQ)])
    , creal(commondata.waveform_IMR_FD[IDX_WF(i,STRAIN)]), cimag(commondata.waveform_IMR_FD[IDX_WF(i,STRAIN)]));
}
"""
        else:
            body += r"""
// Step 6.b: Print the resulting waveform.
for (size_t i = 0; i < commondata.nsteps_IMR; i++) {
    printf("%.15e %.15e %.15e\n", creal(commondata.waveform_IMR[IDX_WF(i,TIME)])
    , creal(commondata.waveform_IMR[IDX_WF(i,STRAIN)]), cimag(commondata.waveform_IMR[IDX_WF(i,STRAIN)]));
}
"""
    if output_commondata:
        body += r"""
commondata_io(&commondata, "commondata.bin");
"""
    body += r"""
free(commondata.dynamics_low);
free(commondata.dynamics_fine);
free(commondata.dynamics_raw);
free(commondata.waveform_low);
free(commondata.waveform_fine);
free(commondata.waveform_inspiral);
free(commondata.waveform_IMR);
"""
    if frequency_domain:
        body += r"""
free(commondata.waveform_IMR_FD);
"""
    body += r"""
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

# register SEOBNRv5 coefficients
BHaH.seobnr.SEOBNRv5_aligned_spin_coefficients.register_CFunction_SEOBNRv5_aligned_spin_coefficients()

# register initial condition routines
BHaH.seobnr.initial_conditions.SEOBNRv5_aligned_spin_multidimensional_root_wrapper.register_CFunction_SEOBNRv5_multidimensional_root_wrapper()
BHaH.seobnr.initial_conditions.SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit.register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit()
# SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_dRHS.register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_dRHS()
# SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_RHSdRHS.register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_RHSdRHS()
# SEOBNRv5_aligned_spin_initial_conditions_conservative.register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_conservative()
BHaH.seobnr.initial_conditions.SEOBNRv5_aligned_spin_initial_conditions_conservative_nodf.register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_conservative_nodf()
BHaH.seobnr.initial_conditions.SEOBNRv5_aligned_spin_radial_momentum_condition.register_CFunction_SEOBNRv5_aligned_spin_radial_momentum_condition()
BHaH.seobnr.initial_conditions.SEOBNRv5_aligned_spin_initial_conditions_dissipative.register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_dissipative()

# register trajectory integration and processing routines
BHaH.seobnr.dynamics.eval_abs_deriv.register_CFunction_eval_abs_deriv()
BHaH.seobnr.dynamics.find_local_minimum_index.register_CFunction_find_local_minimum_index()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_argrelmin.register_CFunction_SEOBNRv5_aligned_spin_argrelmin()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_augments.register_CFunction_SEOBNRv5_aligned_spin_augments()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_interpolate_dynamics.register_CFunction_SEOBNRv5_aligned_spin_interpolate_dynamics()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_iterative_refinement.register_CFunction_SEOBNRv5_aligned_spin_iterative_refinement()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_right_hand_sides.register_CFunction_SEOBNRv5_aligned_spin_right_hand_sides()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_ode_integration.register_CFunction_SEOBNRv5_aligned_spin_ode_integration(
    perform_iterative_refinement
)

# register inspiral waveform routines
BHaH.seobnr.inspiral_waveform.SEOBNRv5_aligned_spin_gamma_wrapper.register_CFunction_SEOBNRv5_aligned_spin_gamma_wrapper()
BHaH.seobnr.inspiral_waveform.SEOBNRv5_aligned_spin_interpolate_modes.register_CFunction_SEOBNRv5_aligned_spin_interpolate_modes()
BHaH.seobnr.inspiral_waveform.SEOBNRv5_aligned_spin_waveform_from_dynamics.register_CFunction_SEOBNRv5_aligned_spin_waveform_from_dynamics()
if precompute_waveform_coefficients_flag:
    BHaH.seobnr.inspiral_waveform_precomputed.SEOBNRv5_aligned_spin_waveform_coefficients.register_CFunction_SEOBNRv5_aligned_spin_waveform_coefficients()
    BHaH.seobnr.inspiral_waveform_precomputed.SEOBNRv5_aligned_spin_waveform_precomputed.register_CFunction_SEOBNRv5_aligned_spin_waveform()
    BHaH.seobnr.inspiral_waveform_precomputed.SEOBNRv5_aligned_spin_flux_precomputed.register_CFunction_SEOBNRv5_aligned_spin_flux()
else:
    BHaH.seobnr.inspiral_waveform.SEOBNRv5_aligned_spin_waveform.register_CFunction_SEOBNRv5_aligned_spin_waveform()
    BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_flux.register_CFunction_SEOBNRv5_aligned_spin_flux()

# register additional commondata parameters needed for SEBOBv2 (but not needed for SEOBNR)
par.register_CodeParameters(
    "REAL",
    __name__,
    ["t_p_BOB", "Omega_0_BOB"],
    commondata=True,
    add_to_parfile=False,
    add_to_set_CodeParameters_h=False,
)

if __name__ == "__main__":
    # retaining this print statement as we will want to add usage options (aligned or precessing, etc) in the future along with help statements
    #    print(
    #        """Generating a compileable C project to calculate gravitational waveforms using the SEOBNRv5 and BOB model!
    # To learn more about usage options, run: python nrpy/example/seobnrv5_aligned_spin_inspiral.py -h
    # """
    #    )

    # Register some functions/code parameters based on input flags
    if frequency_domain_flag:
        par.register_CodeParameter(
            "double complex *restrict",
            __name__,
            "waveform_IMR_FD",
            commondata=True,
            add_to_parfile=False,
            add_to_set_CodeParameters_h=False,
        )
        par.register_CodeParameter(
            "size_t",
            __name__,
            "nsteps_IMR_FD",
            commondata=True,
            add_to_parfile=False,
            add_to_set_CodeParameters_h=False,
        )
        BHaH.seobnr.fft_utils.SEOBNRv5_aligned_spin_FD_waveform.register_CFunction_SEOBNRv5_aligned_spin_FD_waveform()
        BHaH.seobnr.fft_utils.SEOBNRv5_aligned_spin_process_waveform.register_CFunction_SEOBNRv5_aligned_spin_process_waveform()

    BHaH.seobnr.nqc_corrections.SEBOBv2_NQC_corrections.register_CFunction_SEBOBv2_NQC_corrections()
    BHaH.seobnr.nqc_corrections.BOB_v2_NQC_rhs.register_CFunction_BOB_v2_NQC_rhs()

    # set up merger-ringdown routines based on input flags
    BHaH.seobnr.merger_waveform.BOB_v2_find_tp_Omega0.register_CFunction_BOB_v2_find_tp_Omega0()
    BHaH.seobnr.merger_waveform.BOB_v2_peak_strain_conditions.register_CFunction_BOB_v2_peak_strain_conditions()
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
]
if frequency_domain_flag:
    additional_includes.append(str(Path("fftw3.h")))
BHaH.BHaH_defines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    additional_includes=additional_includes,
    enable_intrinsics=False,
    enable_rfm_precompute=False,
    supplemental_defines_dict={
        "SEOBNR": """
#include<complex.h>
#define COMPLEX double complex
#define NUMVARS 8
#define TIME 0
#define FREQ 0
#define R 1
#define PHI 2
#define PRSTAR 3
#define PPHI 4
#define H 5
#define OMEGA 6
#define OMEGA_CIRC 7
#define IDX(idx, var) ((idx)*NUMVARS + (var))
#define NUMMODES 2
#define STRAIN 1
#define IDX_WF(idx,var) ((idx)*NUMMODES + (var))
typedef struct {
  gsl_spline *spline;
  gsl_interp_accel *acc;
} spline_data;
"""
    },
)
register_CFunction_main_c(
    output_waveform_flag,
    frequency_domain_flag,
    precompute_waveform_coefficients_flag,
    output_commondata_flag,
)

addl_cflags = ["$(shell gsl-config --cflags)"]
if frequency_domain_flag:
    addl_cflags.append("-lfftw3 -lm")
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
