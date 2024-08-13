"""
Set up a complete C code project for setting 3.5PN quasicircular momenta for binary black holes, using NRPyPN.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import os

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import shutil

import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.cmdline_input_and_parfiles as cmdpar
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile
import nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library as seobnr_CCL
import nrpy.params as par

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "seobnrv5_aligned_spin_inspiral"

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)


#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]
def register_CFunction_main_c() -> None:
    """Generate a simplified C main() function for computing SEOBNRv5 Hamiltonian and its derivatives."""
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """-={ main() function }=-
Step 1.a: Set each commondata CodeParameter to default.
Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
Step 2: Compute Hamiltonian and derivatives."""
    cfunc_type = "int"
    name = "main"
    params = "int argc, const char *argv[]"
    body = r"""  commondata_struct commondata; // commondata contains parameters common to all grids.

// Step 1.a: Set each commondata CodeParameter to default.
commondata_struct_set_to_default(&commondata);
// Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
cmdline_input_and_parfile_parser(&commondata, argc, argv);

// Step 2: Compute SEOBNRv5 Hamiltonian.
SEOBNRv5_aligned_spin_Hamiltonian(&commondata);

// Step 3: Compute SEOBNRv5 Hamiltonian AND derivatives.
SEOBNRv5_aligned_spin_Hamiltonian_and_derivs(&commondata);

// Step 4: Compute SEOBNRv5 Hamiltonian's circular derivatives.
SEOBNRv5_aligned_spin_Hamiltonian_circular_derivs(&commondata);

// Step 5: Compute SEOBNRv5's gravitational wave flux.
SEOBNRv5_aligned_spin_flux(&commondata);
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


seobnr_CCL.register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian()
seobnr_CCL.register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_and_derivs()
seobnr_CCL.register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_derivs()
seobnr_CCL.register_CFunction_SEOBNRv5_aligned_spin_flux()
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
    cmdline_inputs=[
        "initial_separation",
        "mass_Msun",
        "mass_ratio",
        "chi1",
        "chi2",
    ],
)
Bdefines_h.output_BHaH_defines_h(project_dir=project_dir, enable_simd=False)
register_CFunction_main_c()

Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir,
    project_name=project_name,
    exec_or_library_name=project_name,
)

print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
)
print(f"    Parameter file can be found in {project_name}.par")
