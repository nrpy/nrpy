"""
Set up a complete C code project for setting 3.5PN quasicircular momenta for binary black holes, using NRPyPN.

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
import sys
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
import nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_dynamics_codegen_library as seobnr_dyn_CCL
import nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_gsl_routines_library as seobnr_gsl
import nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_initial_conditions_codegen_library as seobnr_ic_CCL
import nrpy.params as par

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "seobnrv5_aligned_spin_inspiral"

parallel_codegen_enable = True

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)


# Define a function to print instructions in case there is a help flag
def print_help() -> None:
    """Print a list of usage options for running seobnrv5_aligned_spin_inspiral."""
    print(
        """
Generate a compileable C project to calculate gravitational waveforms using the SEOBNRv5 and BOB model!
Usage: python nrpy/example/seobnrv5_aligned_spin_inspiral.py [options] [arguments]

Options:
  -h, --help, -?, --?           Show this help message and exit.
  --output_waveform=True        Enable/disable printing the waveforms (defaults to True)
  --frequency_domain=True       Returns waveform in the frequency domain (defaults to False)

Example:
  python nrpy/example/seobnrv5_aligned_spin_inspiral.py --output_waveform=True --frequency_domain=True
          Generates a C project that calculates the FFT-ed frequency domain waveform using SEOBNRv5 and BOB
          and prints the real and imaginary part of the strain as a frequency series.
    """
    )


# Set flags to default values

# Flag to compute the SEOBNR waveform in the frequency domain.
frequency_domain_flag = False

# Flag to output the SEOBNRv5 waveform using a print statement like lalsimulation does.
output_waveform_flag = True


#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]


def register_CFunction_main_c(
    output_waveform: bool = True, frequency_domain: bool = False
) -> None:
    """
    Generate a simplified C main() function for computing SEOBNRv5 Hamiltonian and its derivatives.

    :param output_waveform: Flag to enable/disable printing the waveform
    :param frequency_domain: Flag to enable/disable FFT to get a frequency domain waveform
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
// Step 2.a: Compute SEOBNRv5 conservative initial conditions.
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


seobnr_gsl.register_CFunction_handle_gsl_return_status()
seobnr_gsl.register_CFunction_SEOBNRv5_multidimensional_root_wrapper()
seobnr_CCL.register_CFunction_SEOBNRv5_aligned_spin_flux()
seobnr_CCL.register_CFunction_SEOBNRv5_aligned_spin_right_hand_sides()
seobnr_ic_CCL.register_CFunction_SEOBNRv5_aligned_spin_coefficients()
seobnr_ic_CCL.register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit()
seobnr_ic_CCL.register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_dRHS()
seobnr_ic_CCL.register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_RHSdRHS()
seobnr_ic_CCL.register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_conservative()
seobnr_ic_CCL.register_CFunction_SEOBNRv5_aligned_spin_radial_momentum_condition()
seobnr_ic_CCL.register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_dissipative()
seobnr_dyn_CCL.register_CFunction_SEOBNRv5_aligned_spin_augments()
seobnr_dyn_CCL.register_CFunction_SEOBNRv5_aligned_spin_argrelmin()
seobnr_dyn_CCL.register_CFunction_SEOBNRv5_aligned_spin_iterative_refinement()
seobnr_dyn_CCL.register_CFunction_SEOBNRv5_aligned_spin_intepolate_dynamics()
seobnr_gsl.register_CFunction_SEOBNRv5_aligned_spin_gamma_wrapper()
seobnr_CCL.register_CFunction_SEOBNRv5_aligned_spin_waveform()
seobnr_CCL.register_CFunction_SEOBNRv5_aligned_spin_waveform_from_dynamics()
seobnr_wf_CCL.register_CFunction_SEOBNRv5_aligned_spin_unwrap()
seobnr_wf_CCL.register_CFunction_SEOBNRv5_aligned_spin_interpolate_modes()
seobnr_wf_CCL.register_CFunction_SEOBNRv5_NQC_corrections()
BOB_CCL.register_CFunction_BOB_aligned_spin_waveform()
BOB_CCL.register_CFunction_BOB_aligned_spin_waveform_from_times()
BOB_CCL.register_CFunction_BOB_aligned_spin_NQC_rhs()
seobnr_wf_CCL.register_CFunction_SEOBNRv5_aligned_spin_IMR_waveform()

if __name__ == "__main__":
    # SEOBNRv5 uses an iterative refinement routine to find the location of the peak of the orbital frequency.
    # This is not done in a robust manner and disabling it does not impact accuracy.
    # This flag helps enable/disable the routine at the codegen level so that we can assess the impact on performance and accuracy.
    perform_iterative_refinement = False
    seobnr_dyn_CCL.register_CFunction_SEOBNRv5_aligned_spin_ode_integration(
        perform_iterative_refinement
    )

    # If running without arguments, output the ideal waveform code and inform about running help.
    if len(sys.argv) == 1:
        print(
            """
Generating a compileable C project to calculate gravitational waveforms using the SEOBNRv5 and BOB model!
To learn more about usage options, run:
python nrpy/example/seobnrv5_aligned_spin_inspiral.py -h
or 
python nrpy/example/seobnrv5_aligned_spin_inspiral.py --help
"""
        )

    # Check for help flags
    if any(flag in sys.argv for flag in ["-h", "--help", "-?", "--?"]):
        print_help()
        sys.exit(0)

    # Check for additional flags
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            if arg.startswith("--output_waveform="):
                output_waveform_flag = arg.split("=")[1].lower() == "true"
            elif arg.startswith("--frequency_domain="):
                frequency_domain_flag = arg.split("=")[1].lower() == "true"

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
"""
    },
    enable_intrinsics=False,
)
register_CFunction_main_c(output_waveform_flag, frequency_domain_flag)

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
