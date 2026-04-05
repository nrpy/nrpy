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
import nrpy.params as par
from nrpy.infrastructures import BHaH

par.set_parval_from_str("Infrastructure", "BHaH")
enable_parallel_codegen = True
par.set_parval_from_str("enable_parallel_codegen", enable_parallel_codegen)

# Development flags (NOT command-line-tunable)
# Flag to output the commondata struct to a file.
output_commondata_flag = True
# Flag to output the SEOBNRv5 waveform using a print statement like lalsimulation does.
# (set to False for performance checks)
output_waveform_flag = True


# Command-line-tunable flags (tuned through choice of approximant, defaults to all BOBs)
# Flag to use numerical relativity fits to evaluate Non Quasi-Circular corrections to the inspiral.
numerical_relativity_nqc_flag = False
# Flag to compute SEOBNRv5's phenomenological fit for the merger-ringdown waveform.
seobnrv5_merger_ringdown_flag = False

#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]


def register_CFunction_main_c(
    output_waveform: bool = True,
    output_commondata: bool = False,
) -> None:
    """
    Generate a simplified C main() function for computing the SEBOB waveform.

    :param output_waveform: Flag to enable/disable printing the waveform
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
    if output_waveform:
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
BHaH.seobnr.utils.root_finding_1d.register_CFunction_root_finding_1d()
BHaH.seobnr.utils.root_finding_multidimensional.register_CFunction_root_finding_multidimensional()


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
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_augments.register_CFunction_SEOBNRv5_aligned_spin_augments()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_interpolate_dynamics.register_CFunction_SEOBNRv5_aligned_spin_interpolate_dynamics()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_iterative_refinement.register_CFunction_SEOBNRv5_aligned_spin_iterative_refinement()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_right_hand_sides.register_CFunction_SEOBNRv5_aligned_spin_right_hand_sides()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_ode_integration.register_CFunction_SEOBNRv5_aligned_spin_ode_integration()

# register inspiral waveform routines
BHaH.seobnr.inspiral_waveform.SEOBNRv5_aligned_spin_gamma_wrapper.register_CFunction_SEOBNRv5_aligned_spin_gamma_wrapper()
BHaH.seobnr.inspiral_waveform.SEOBNRv5_aligned_spin_interpolate_modes.register_CFunction_SEOBNRv5_aligned_spin_interpolate_modes()
BHaH.seobnr.inspiral_waveform.SEOBNRv5_aligned_spin_waveform_from_dynamics.register_CFunction_SEOBNRv5_aligned_spin_waveform_from_dynamics()
BHaH.seobnr.inspiral_waveform.SEOBNRv5_aligned_spin_waveform.register_CFunction_SEOBNRv5_aligned_spin_waveform()
BHaH.seobnr.dynamics.SEOBNRv5_aligned_spin_flux.register_CFunction_SEOBNRv5_aligned_spin_flux()

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
    parser.add_argument(
        "-calibration_no_spin",
        action="store_true",
        help="Set up par file and Hamiltonian coefficients for non-spinning calibration",
    )
    parser.add_argument(
        "-calibration_spin",
        action="store_true",
        help="Set up par file and Hamiltonian coefficients for spin-dependent calibration",
    )
    args = parser.parse_args()
    # The SEOBNRv5 calibration process is done in two steps:
    # 1. Calibration of the non-spinning coefficients
    # 2. Calibration of the spin-dependent coefficients
    # Therefore, the C code can only be generated for one of the above calibration options.
    if not args.seobnrv5_bob and not args.seobnrv5_nrnqc_bob and not args.seobnrv5_nrpy:
        print("Defaulting to seobnrv5_bob.")
        args.seobnrv5_bob = True
    if args.calibration_no_spin and args.calibration_spin:
        raise ValueError(
            "calibration_no_spin and calibration_spin cannot both be True."
        )
    # Code-generation-time parameters:
    project_name = "seobnrv5_nrpy"
    if args.seobnrv5_bob:
        project_name = "seobnrv5_bob"
    if args.seobnrv5_nrnqc_bob:
        project_name = "seobnrv5_nrnqc_bob"
    if args.calibration_no_spin:
        project_name = f"{project_name}_calibration_no_spin"
    if args.calibration_spin:
        project_name = f"{project_name}_calibration_spin"
    project_dir = os.path.join("project", project_name)

    # First clean the project directory, if it exists.
    shutil.rmtree(project_dir, ignore_errors=True)

    # register SEOBNRv5 coefficients
    BHaH.seobnr.SEOBNRv5_aligned_spin_coefficients.register_CFunction_SEOBNRv5_aligned_spin_coefficients(
        args.calibration_no_spin, args.calibration_spin
    )

    if args.seobnrv5_nrnqc_bob:
        numerical_relativity_nqc_flag = True
    if args.seobnrv5_nrpy:
        print("Using native SEOBNRv5 model.")
        numerical_relativity_nqc_flag = True
        seobnrv5_merger_ringdown_flag = True
    # Register some functions/code parameters based on input flags
    # set up NQC correction routines based on input flags
    BHaH.seobnr.nqc_corrections.SEOBNRv5_aligned_spin_NQC_corrections.register_CFunction_SEOBNRv5_aligned_spin_NQC_corrections(
        numerical_relativity_nqc_flag
    )
    if numerical_relativity_nqc_flag:
        BHaH.seobnr.nqc_corrections.SEOBNRv5_aligned_spin_NQC_rhs.register_CFunction_SEOBNRv5_aligned_spin_NQC_rhs()
    else:
        BHaH.seobnr.nqc_corrections.BOB_aligned_spin_NQC_rhs.register_CFunction_BOB_aligned_spin_NQC_rhs()

    # set up merger-ringdown routines based on input flags
    if seobnrv5_merger_ringdown_flag:
        BHaH.seobnr.merger_waveform.SEOBNRv5_aligned_spin_merger_waveform.register_CFunction_SEOBNRv5_aligned_spin_merger_waveform()
        BHaH.seobnr.merger_waveform.SEOBNRv5_aligned_spin_merger_waveform_from_times.register_CFunction_SEOBNRv5_aligned_spin_merger_waveform_from_times()
    else:
        BHaH.seobnr.merger_waveform.BOB_aligned_spin_waveform.register_CFunction_BOB_aligned_spin_waveform()
        BHaH.seobnr.merger_waveform.BOB_aligned_spin_waveform_from_times.register_CFunction_BOB_aligned_spin_waveform_from_times()

    # register IMR waveform generation routine
    BHaH.seobnr.SEOBNRv5_aligned_spin_IMR_waveform.register_CFunction_SEOBNRv5_aligned_spin_IMR_waveform(
        seobnrv5_merger_ringdown_flag
    )
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
BHaH.BHaH_defines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    additional_includes=additional_includes,
    enable_rfm_precompute=False,
    supplemental_defines_dict={"SEOBNR": """
#include<complex.h>
#define COMPLEX double complex
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
#define NUMMODES 2
#define STRAIN 1
#define IDX_WF(idx,var) ((idx)*NUMMODES + (var))
typedef struct {
  gsl_spline *spline;
  gsl_interp_accel *acc;
} spline_data;
"""},
)
register_CFunction_main_c(
    output_waveform_flag,
    output_commondata_flag,
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
