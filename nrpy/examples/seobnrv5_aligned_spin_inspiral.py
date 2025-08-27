"""
Set up a complete C code project for setting up the Spinning Effective-to-Backwards One Body (SEBOB) model with SEOBNRv5 and BOB.

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
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.cmdline_input_and_parfiles as cmdpar
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile
import nrpy.infrastructures.BHaH.seobnr.BOB_C_codegen_library as BOB_CCL
import nrpy.infrastructures.BHaH.seobnr.SEOBNR_BOB_C_waveform_codegen_library as seobnr_wf_CCL
import nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library as seobnr_CCL
import nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library_precomputed as seobnr_CCL_precomp
import nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_dynamics_codegen_library as seobnr_dyn_CCL
import nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_gsl_routines_library as seobnr_gsl
import nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_initial_conditions_codegen_library as seobnr_ic_CCL
import nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_merger_codegen_library as seobnr_mr_CCL
import nrpy.params as par

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "seobnrv5_aligned_spin_inspiral"

parallel_codegen_enable = True

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)

# Development flags (NOT user-tunable)
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


# User-tunable flags (tuned through choice of approximant, defaults to all BOBs)
# Flag to use numerical relativity fits to evaluate Non Quasi-Circular corrections to the inspiral.
numerical_relativity_nqc_flag = False
# Flag to compute SEOBNRv5's phenomological fit for the merger-ringdown waveform.
seobnv5_merger_ringdown_flag = False

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
SEOBNRv5_aligned_spin_NQC_corrections(&commondata);
// Step 7.a Compute the IMR waveform
SEOBNRv5_aligned_spin_IMR_waveform(&commondata);
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


seobnr_CCL.register_CFunction_commondata_io()
seobnr_gsl.register_CFunction_handle_gsl_return_status()
seobnr_gsl.register_CFunction_SEOBNRv5_multidimensional_root_wrapper()
seobnr_CCL.register_CFunction_SEOBNRv5_aligned_spin_right_hand_sides()
seobnr_ic_CCL.register_CFunction_SEOBNRv5_aligned_spin_coefficients()
seobnr_ic_CCL.register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit()
# seobnr_ic_CCL.register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_dRHS()
# seobnr_ic_CCL.register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_RHSdRHS()
seobnr_ic_CCL.register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_conservative_nodf()
seobnr_ic_CCL.register_CFunction_SEOBNRv5_aligned_spin_radial_momentum_condition()
seobnr_ic_CCL.register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_dissipative()
seobnr_dyn_CCL.register_CFunction_SEOBNRv5_aligned_spin_augments()
seobnr_dyn_CCL.register_CFunction_SEOBNRv5_aligned_spin_argrelmin()
seobnr_dyn_CCL.register_CFunction_eval_abs_deriv()
seobnr_dyn_CCL.register_CFunction_find_local_minimum_index()
seobnr_dyn_CCL.register_CFunction_SEOBNRv5_aligned_spin_iterative_refinement()
seobnr_dyn_CCL.register_CFunction_SEOBNRv5_aligned_spin_intepolate_dynamics()
seobnr_gsl.register_CFunction_SEOBNRv5_aligned_spin_gamma_wrapper()
seobnr_CCL.register_CFunction_SEOBNRv5_aligned_spin_waveform_from_dynamics()
if precompute_waveform_coefficients_flag:
    seobnr_CCL_precomp.register_CFunction_SEOBNRv5_aligned_spin_waveform_coefficients()
    seobnr_CCL_precomp.register_CFunction_SEOBNRv5_aligned_spin_waveform()
    seobnr_CCL_precomp.register_CFunction_SEOBNRv5_aligned_spin_flux()
else:
    seobnr_CCL.register_CFunction_SEOBNRv5_aligned_spin_waveform()
    seobnr_CCL.register_CFunction_SEOBNRv5_aligned_spin_flux()
seobnr_wf_CCL.register_CFunction_SEOBNRv5_aligned_spin_unwrap()
seobnr_wf_CCL.register_CFunction_SEOBNRv5_aligned_spin_interpolate_modes()

if __name__ == "__main__":
    print(
        """Generating a compileable C project to calculate gravitational waveforms using the SEOBNRv5 and BOB model! 
To learn more about usage options, run: python nrpy/example/seobnrv5_aligned_spin_inspiral.py -h
"""
    )
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate a compileable C project to calculate gravitational waveforms using the SEOBNRv5 and BOB model!"
    )

    parser.add_argument(
        "-seobnrv5_bob",
        action="store_true",
        help="SEOBNRv5 model with BOB-informed NQC corrections and merger-ringdown (default)",
    )
    parser.add_argument(
        "-seobnrv5_nrnqc_bob",
        action="store_true",
        help="SEOBNRv5 model with native NQC corrections and BOB-informed merger-ringdown",
    )
    parser.add_argument(
        "-seobnrv5_nrpy", action="store_true", help="native SEOBNRv5 model"
    )
    args = parser.parse_args()

    if not args.seobnrv5_bob and not args.seobnrv5_nrnqc_bob and not args.seobnrv5_nrpy:
        print("Defaulting to seobnrv5_bob.")
        args.seobnrv5_bob = True
    if args.seobnrv5_nrnqc_bob:
        use_numerical_relativity_nqc_flag = True
    if args.seobnrv5_nrpy:
        use_numerical_relativity_nqc_flag = True
        use_seobnrv5_merger_ringdown_flag = True

    seobnr_dyn_CCL.register_CFunction_SEOBNRv5_aligned_spin_ode_integration(
        perform_iterative_refinement
    )
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
        seobnr_CCL.register_CFunction_SEOBNRv5_aligned_spin_FD_waveform()
        seobnr_gsl.register_CFunction_SEOBNRv5_aligned_spin_process_waveform()

    seobnr_wf_CCL.register_CFunction_SEOBNRv5_NQC_corrections(
        numerical_relativity_nqc_flag
    )
    seobnr_wf_CCL.register_CFunction_SEOBNRv5_aligned_spin_IMR_waveform(
        seobnv5_merger_ringdown_flag
    )
    if numerical_relativity_nqc_flag:
        seobnr_mr_CCL.register_CFunction_SEOBNRv5_aligned_spin_NQC_rhs()
    else:
        BOB_CCL.register_CFunction_BOB_aligned_spin_NQC_rhs()
    if seobnv5_merger_ringdown_flag:
        seobnr_mr_CCL.register_CFunction_SEOBNRv5_aligned_spin_merger_waveform()
        seobnr_mr_CCL.register_CFunction_SEOBNRv5_aligned_spin_merger_waveform_from_times()
    else:
        BOB_CCL.register_CFunction_BOB_aligned_spin_waveform()
        BOB_CCL.register_CFunction_BOB_aligned_spin_waveform_from_times()

    pcg.do_parallel_codegen()
#########################################################
# STEP 3: Generate header files, register C functions and
#         command line parameters, set up boundary conditions,
#         and create a Makefile for this project.
#         Project is output to project/[project_name]/
CPs.write_CodeParameters_h_files(set_commondata_only=True, project_dir=project_dir)
CPs.register_CFunctions_params_commondata_struct_set_to_default()
cmdpar.generate_default_parfile(project_dir=project_dir, project_name=project_name)
cmdpar.register_CFunction_cmdline_input_and_parfile_parser(
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
Bdefines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    additional_includes=additional_includes,
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
    enable_intrinsics=False,
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
Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
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
