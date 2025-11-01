# nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_inner_only.py
"""
Generates C functions for applying inner boundary conditions on a curvilinear grid.

The generated functions apply parity conditions to grid functions at inner boundary
points (e.g., for black hole excision). The boundary mappings are pre-computed
and stored in a `bc_struct`. This module provides routines to apply these
conditions to either all evolved grid functions or a specific list of auxiliary
grid functions.

This process is documented in the NRPy tutorial:
Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
         Samuel D. Tootle
         sdtootle **at** gmail **dot* com
"""
import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par
from nrpy.infrastructures import BHaH


def register_CFunction_apply_bcs_inner_only_specific_gfs() -> None:
    """
    Register C function for filling inner boundary points for specific gfs on the computational grid.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    cfc.CFunction_dict.clear()
    ...    register_CFunction_apply_bcs_inner_only_specific_gfs()
    ...    generated_str = cfc.CFunction_dict[f'apply_bcs_inner_only_specific_gfs'].full_function
    ...    validation_desc = f"apply_bcs_inner_only_specific_gfs__{parallelization}"
    ...    validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    includes = ["BHaH_defines.h"]
    desc = r"""
Apply BCs to specific grid functions at inner boundary points only,
using data stored in bcstruct->inner_bc_array.
These structs are set in bcstruct_set_up().
Inner boundary points map to either the grid
interior ("pure inner") or to pure outer
boundary points ("inner maps to outer").
"""
    cfunc_type = "void"
    name = "apply_bcs_inner_only_specific_gfs"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs, const int num_gfs, const int *gf_parities, const int *gfs_to_sync"
    parallelization = par.parval_from_str("parallelization")

    kernel_body = BHaH.CurviBoundaryConditions.apply_bcs_inner_only.generate_apply_bcs_inner_only__kernel_body(
        loop_bounds="num_gfs",
        parity_ary="gf_parities",
        gf_index="gfs_to_sync[which_gf]",
    )

    comments = "Apply BCs to inner boundary points only for specified GFs."
    # Prepare the argument dicts
    arg_dict_cuda = {
        "num_inner_boundary_points": "const int",
        "inner_bc_array": "const innerpt_bc_struct *restrict",
        "gfs": "REAL *restrict",
        "num_gfs": "const int",
        "gf_parities": "const int *restrict",
        "gfs_to_sync": "const int *restrict",
    }
    arg_dict_host = {
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }
    prefunc, new_body = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body,
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=comments,
        launch_dict={
            "blocks_per_grid": [
                "MAX(1U, (num_inner_boundary_points + threads_in_x_dir - 1) / threads_in_x_dir)"
            ],
            "stream": "params->grid_idx % NUM_STREAMS",
        },
        thread_tiling_macro_suffix="CURVIBC_INNER",
    )
    kernel_launch_body = """
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;
  const innerpt_bc_struct *restrict inner_bc_array = bcstruct->inner_bc_array;
  const int num_inner_boundary_points = bc_info->num_inner_boundary_points;
"""
    if parallelization in ["cuda"]:
        kernel_launch_body += """
        // Allocate device memory for gfs_to_sync
        int *gfs_to_sync_device;
        BHAH_MALLOC_DEVICE(gfs_to_sync_device, num_gfs * sizeof(int));
        cudaMemcpy(gfs_to_sync_device, gfs_to_sync, num_gfs * sizeof(int), cudaMemcpyHostToDevice);
        """
        new_body = new_body.replace("gfs_to_sync", "gfs_to_sync_device")
    kernel_launch_body += new_body

    if parallelization in ["cuda"]:
        kernel_launch_body += """
        BHAH_DEVICE_SYNC();
        // Free device memory for gfs_to_sync
        BHAH_FREE_DEVICE(gfs_to_sync_device);
        """

    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=kernel_launch_body,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
