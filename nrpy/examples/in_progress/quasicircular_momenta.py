"""
Sets up a complete C code project for solving the wave equation
  in Cartesian coordinates, on a cell-centered Cartesian
  grid.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import shutil
import os
from nrpypn.NRPyPN_shortcuts import (
    m1,
    m2,
    chi1U,
    chi2U,
    r,
    n12U,
    n21U,
    S1U,
    S2U,
    p1U,
    p2U,
    num_eval,
)
from nrpypn.PN_p_t import PN_p_t
from nrpypn.PN_p_r import PN_p_r

import nrpy.params as par
import nrpy.c_function as cfc
import nrpy.c_codegen as ccg

import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.main_c as main
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile
import nrpy.infrastructures.BHaH.cmdline_input_and_parfiles as cmdpar

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "quasicircular_momenta"

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]


#########################################################
# STEP 3: Generate header files, register C functions and
#         command line parameters, set up boundary conditions,
#         and create a Makefile for this project.
#         Project is output to project/[project_name]/
CPs.write_CodeParameters_h_files(project_dir=project_dir)
CPs.register_CFunctions_params_commondata_struct_set_to_default()
cmdpar.generate_default_parfile(project_dir=project_dir, project_name=project_name)
cmdpar.register_CFunction_cmdline_input_and_parfile_parser(
    project_name=project_name, cmdline_inputs=["convergence_factor"]
)
Bdefines_h.output_BHaH_defines_h(
    project_dir=project_dir, enable_simd=False, MoL_method="disabled"
)
main.register_CFunction_main_c(
    MoL_method="disabled",
    initial_data_desc="Quasi-circular initial data for",
)

Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir, project_name=project_name, exec_name=project_name
)
print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
)
print(f"    Parameter file can be found in {project_name}.par")

# print(cfc.CFunction_dict["initial_data"].full_function)
# print(cfc.CFunction_dict["rhs_eval"].full_function)
# print(cfc.CFunction_dict["apply_bcs"].full_function)
# print(cfc.CFunction_dict["parameter_file_read_and_parse"].full_function)
# print(cfc.CFunction_dict["main"].full_function)
