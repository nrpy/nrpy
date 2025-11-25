"""
Set up a complete JAX project for setting up the Spinning Effective-to-Backwards One Body (SEBOB) model with SEOBNRv5 and BOB.
The model corresponds to the published version in Mahesh et. al. 2025 (https://arxiv.org/abs/2508.20418)

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

import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.infrastructures import JAX

par.set_parval_from_str("Infrastructure", "JAX")

# Code-generation-time parameters:
project_name = "sebobv1_jax"

enable_parallel_codegen = True

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("enable_parallel_codegen", enable_parallel_codegen)

#########################################################
# STEP 2: Declare core JAX functions & register each to
#         pyfc.py_function_dict["function_name"]


# register SEOBNRv5 coefficients
JAX.SEOBNRv5_aligned_spin_coefficients.register_PyFunction_SEOBNRv5_aligned_spin_coefficients()

if __name__ == "__main__":
    print(
        """Generating a JAX project to calculate gravitational waveforms using the SEOBNRv5 and BOB model! 
Currently, this example only sets up the SEOBNRv5 aligned spin coefficients.
"""
    )
    pcg.do_parallel_codegen()
    # Write the generated JAX functions to the project directory
    JAX.jax_project_generator.output_PyFunction_files_and_construct_project(
        project_dir=project_dir,
        project_name=project_name,
        lib_function_prefix="",
    )
    print(f"Finished! Now go into project/{project_name}")
