# nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_outerextrap_and_inner.py
"""
Generates C code for applying boundary conditions on curvilinear grids.

This module builds and registers C functions that apply boundary conditions by
first setting outer boundary ghost zones using 2nd-order extrapolation, and
then setting the inner boundary ghost zones.

This two-step approach ensures that inner boundary points that map to the
outer boundary are handled correctly. Functions are generated to handle all
evolved gridfunctions, as well as specific auxiliary gridfunctions.

The method is documented in the NRPy tutorial:
"Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb".

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
         Samuel D. Tootle
         sdtootle **at** gmail **dot* com
"""

import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par
from nrpy import c_function


def generate_apply_bcs_outerextrap_and_inner_only__kernel_body(
    loop_bounds: str = "NUM_EVOL_GFS", gf_index: str = "which_gf"
) -> str:
    """
    Generate the kernel body for the apply_bcs_inner_only function.

    This function is used to apply boundary conditions to inner boundary points
    on the computational grid.

    :param loop_bounds: The number of loop bounds (default is NUM_EVOL_GFS).
    :param gf_index: The storage for the index of the grid function (default is which_gf).

    :return: The kernel body as a string.
    """
    parallelization = par.parval_from_str("parallelization")

    # Specify kernel body
    kernel_body = f"{parallel_utils.get_loop_parameters(parallelization)}\n"

    kernel_body += (
        """
for (int idx2d = tid0; idx2d < num_pure_outer_boundary_points; idx2d+=stride0) {"""
        if parallelization == "cuda"
        else """
#pragma omp parallel for
    for (int idx2d = 0; idx2d < num_pure_outer_boundary_points; idx2d++) {"""
    )
    kernel_body += f"""
const short i0 = pure_outer_bc_array[idx2d].i0;
const short i1 = pure_outer_bc_array[idx2d].i1;
const short i2 = pure_outer_bc_array[idx2d].i2;
const short FACEX0 = pure_outer_bc_array[idx2d].FACEX0;
const short FACEX1 = pure_outer_bc_array[idx2d].FACEX1;
const short FACEX2 = pure_outer_bc_array[idx2d].FACEX2;
const int idx_offset0 = IDX3(i0, i1, i2);
const int idx_offset1 = IDX3(i0 + 1 * FACEX0, i1 + 1 * FACEX1, i2 + 1 * FACEX2);
const int idx_offset2 = IDX3(i0 + 2 * FACEX0, i1 + 2 * FACEX1, i2 + 2 * FACEX2);
const int idx_offset3 = IDX3(i0 + 3 * FACEX0, i1 + 3 * FACEX1, i2 + 3 * FACEX2);
for (int which_gf = 0; which_gf < {loop_bounds}; which_gf++) {{
    // *** Apply 2nd-order polynomial extrapolation BCs to all outer boundary points. ***
    gfs[IDX4pt({gf_index}, idx_offset0)] =
        + 3.0 * gfs[IDX4pt({gf_index}, idx_offset1)]
        - 3.0 * gfs[IDX4pt({gf_index}, idx_offset2)]
        + 1.0 * gfs[IDX4pt({gf_index}, idx_offset3)];
}}
}}
"""

    return kernel_body


def generate_prefunc__apply_bcs_outerextrap_and_inner_only() -> str:
    """
    Generate the prefunction string for apply_bcs_outerextrap_and_inner.

    This requires a function that will launch the device kernel as well
    as the device kernel itself.

    :returns: prefunc string
    """
    # Header details for function that will launch the GPU kernel
    desc = "Apply BCs to pure boundary points"
    params = "const params_struct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs"
    name = "apply_bcs_outerextrap_and_inner_only"
    cfunc_type = "static void"
    return_str = ""
    prefunc = ""
    parallelization = par.parval_from_str("parallelization")

    kernel_body = generate_apply_bcs_outerextrap_and_inner_only__kernel_body()

    # Generate compute Kernel
    comments = "Apply extrapolation BCs to pure points."
    # Prepare the argument dicts
    arg_dict_cuda = {
        "num_pure_outer_boundary_points": "const int",
        "which_gz": "const int",
        "dirn": "const int",
        "pure_outer_bc_array": "const outerpt_bc_struct *restrict",
        "gfs": "REAL *restrict",
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
                "MAX(1U, (num_pure_outer_boundary_points + threads_in_x_dir -1) / threads_in_x_dir)"
            ],
            "stream": "default",
        },
        thread_tiling_macro_suffix="CURVIBC_EXTRAP",
    )
    # Specify the function body for launching the kernel
    kernel_launch_body = f"""
const bc_info_struct *bc_info = &bcstruct->bc_info;
for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {{
for (int dirn = 0; dirn < 3; dirn++) {{
    if (bc_info->num_pure_outer_boundary_points[which_gz][dirn] > 0) {{
    size_t gz_idx = dirn + (3 * which_gz);
    const outerpt_bc_struct *restrict pure_outer_bc_array = bcstruct->pure_outer_bc_array[gz_idx];
    int num_pure_outer_boundary_points = bc_info->num_pure_outer_boundary_points[which_gz][dirn];
    {new_body}
    }}
  }}
}}
"""
    # Generate the Launch kernel CFunction
    kernel_launch_CFunction = c_function.CFunction(
        includes=[],
        desc=desc,
        cfunc_type=cfunc_type,
        name=f"{name}__launcher",
        params=params,
        body=kernel_launch_body,
    )

    # Append Launch kernel to prefunc
    return_str = prefunc + kernel_launch_CFunction.full_function
    return return_str


def generate_prefunc__apply_bcs_outerextrap_and_inner_only_specific_gfs() -> str:
    """
    Generate the prefunction string for apply_bcs_outerextrap_and_inner.

    This requires a function that will launch the device kernel as well
    as the device kernel itself.

    :returns: prefunc string
    """
    # Header details for function that will launch the GPU kernel
    desc = "Apply BCs to pure boundary points"
    params = "const params_struct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs, const int num_gfs, const int *gfs_to_sync"
    name = "apply_bcs_outerextrap_and_inner_only_specific_gfs"
    cfunc_type = "static void"
    return_str = ""
    prefunc = ""
    parallelization = par.parval_from_str("parallelization")

    kernel_body = generate_apply_bcs_outerextrap_and_inner_only__kernel_body(
        loop_bounds="num_gfs", gf_index="gfs_to_sync[which_gf]"
    )

    # Generate compute Kernel
    comments = "Apply extrapolation BCs to pure points."
    # Prepare the argument dicts
    arg_dict_cuda = {
        "num_pure_outer_boundary_points": "const int",
        "which_gz": "const int",
        "dirn": "const int",
        "pure_outer_bc_array": "const outerpt_bc_struct *restrict",
        "gfs": "REAL *restrict",
        "num_gfs": "const int",
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
                "MAX(1U, (num_pure_outer_boundary_points + threads_in_x_dir -1) / threads_in_x_dir)"
            ],
            "stream": "default",
        },
        thread_tiling_macro_suffix="CURVIBC_EXTRAP",
    )
    # Specify the function body for launching the kernel
    kernel_launch_body = "const bc_info_struct *bc_info = &bcstruct->bc_info;"
    if parallelization in ["cuda"]:
        kernel_launch_body += """
        // Allocate device memory for gfs_to_sync
        int *gfs_to_sync_device;
        BHAH_MALLOC_DEVICE(gfs_to_sync_device, num_gfs * sizeof(int));
        cudaMemcpy(gfs_to_sync_device, gfs_to_sync, num_gfs * sizeof(int), cudaMemcpyHostToDevice);
        """
        new_body = new_body.replace("gfs_to_sync", "gfs_to_sync_device")
    kernel_launch_body += f"""
for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {{
for (int dirn = 0; dirn < 3; dirn++) {{
    if (bc_info->num_pure_outer_boundary_points[which_gz][dirn] > 0) {{
    size_t gz_idx = dirn + (3 * which_gz);
    const outerpt_bc_struct *restrict pure_outer_bc_array = bcstruct->pure_outer_bc_array[gz_idx];
    int num_pure_outer_boundary_points = bc_info->num_pure_outer_boundary_points[which_gz][dirn];
    {new_body}
    }}
  }}
}}
"""
    # Generate the Launch kernel CFunction
    kernel_launch_CFunction = c_function.CFunction(
        includes=[],
        desc=desc,
        cfunc_type=cfunc_type,
        name=f"{name}__launcher",
        params=params,
        body=kernel_launch_body,
    )

    # Append Launch kernel to prefunc
    return_str = prefunc + kernel_launch_CFunction.full_function
    return return_str


###############################
## apply_bcs_outerextrap_and_inner(): Apply extrapolation outer boundary conditions.
##  Function is documented below in desc= and body=.
def register_CFunction_apply_bcs_outerextrap_and_inner() -> None:
    """
    Register C function for filling boundary points with extrapolation and prescribed bcstruct.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name = "apply_bcs_outerextrap_and_inner"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    c_function.CFunction_dict.clear()
    ...    register_CFunction_apply_bcs_outerextrap_and_inner()
    ...    generated_str = c_function.CFunction_dict[f'{name}'].full_function
    ...    validation_desc = f"{name}__{parallelization}"
    ...    validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r"""Suppose the outer boundary point is at the i0=max(i0) face. Then we fit known data at i0-3, i0-2, and i0-1
  to the unique quadratic polynomial that passes through those points, and fill the data at
  i0 with the value implied from the polynomial.
As derived in nrpytutorial's Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb,
  the coefficients must be f_{i0} = f_{i0-3} - 3 f_{i0-2} + 3 f_{i0-1}.
  To check these coefficients are correct, consider
  * f(x0 = constant. Then f_{i0} = f_{i0-3} <- CHECK!
  * f(x) = x. WOLOG suppose x0=0. Then f_{i0} = (-3dx) - 3(-2dx) + 3(-dx) = + dx(-3+6-3) = 0 <- CHECK!
  * f(x) = x^2. WOLOG suppose x0=0. Then f_{i0} = (-3dx)^2 - 3(-2dx)^2 + 3(-dx)^2 = + dx^2(9-12+3) = 0 <- CHECK!"""
    cfunc_type = "void"
    name = "apply_bcs_outerextrap_and_inner"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs"
    body = r"""
  ////////////////////////////////////////////////////////
  // STEP 1 of 2: Apply BCs to pure outer boundary points.
  //              By "pure" we mean that these points are
  //              on the outer boundary and not also on
  //              an inner boundary.
  //              Here we fill in the innermost ghost zone
  //              layer first and move outward. At each
  //              layer, we fill in +/- x0 faces first,
  //              then +/- x1 faces, finally +/- x2 faces,
  //              filling in the edges as we go.
  apply_bcs_outerextrap_and_inner_only__launcher(params, bcstruct, gfs);

  ///////////////////////////////////////////////////////
  // STEP 2 of 2: Apply BCs to inner boundary points.
  //              These map to either the grid interior
  //              ("pure inner") or to pure outer boundary
  //              points ("inner maps to outer"). Those
  //              that map to outer require that outer be
  //              populated first; hence this being
  //              STEP 2 OF 2.
  apply_bcs_inner_only(commondata, params, bcstruct, gfs);
"""
    prefunc = generate_prefunc__apply_bcs_outerextrap_and_inner_only()
    c_function.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )


###############################
## apply_bcs_outerextrap_and_inner(): Apply extrapolation outer boundary conditions.
##  Function is documented below in desc= and body=.
def register_CFunction_apply_bcs_outerextrap_and_inner_specific_gfs() -> None:
    """
    Register C function for filling boundary points with extrapolation and prescribed bcstruct for specific aux gridfunctions.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name = "apply_bcs_outerextrap_and_inner_specific_gfs"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    c_function.CFunction_dict.clear()
    ...    register_CFunction_apply_bcs_outerextrap_and_inner_specific_gfs()
    ...    generated_str = c_function.CFunction_dict[f'{name}'].full_function
    ...    validation_desc = f"{name}__{parallelization}"
    ...    validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r"""*Suppose the outer boundary point is at the i0=max(i0) face. Then we fit known data at i0-3, i0-2, and i0-1
*  to the unique quadratic polynomial that passes through those points, and fill the data at
*  i0 with the value implied from the polynomial.
*As derived in nrpytutorial's Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb,
*  the coefficients must be f_{i0} = f_{i0-3} - 3 f_{i0-2} + 3 f_{i0-1}.
*  To check these coefficients are correct, consider
*  * f(x0 = constant. Then f_{i0} = f_{i0-3} <- CHECK!
*  * f(x) = x. WOLOG suppose x0=0. Then f_{i0} = (-3dx) - 3(-2dx) + 3(-dx) = + dx(-3+6-3) = 0 <- CHECK!
*  * f(x) = x^2. WOLOG suppose x0=0. Then f_{i0} = (-3dx)^2 - 3(-2dx)^2 + 3(-dx)^2 = + dx^2(9-12+3) = 0 <- CHECK!"""
    cfunc_type = "void"
    name = "apply_bcs_outerextrap_and_inner_specific_gfs"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs, const int num_gfs, const int *gfs_to_sync"
    body = r"""
  ////////////////////////////////////////////////////////
  // STEP 1 of 2: Apply BCs to pure outer boundary points.
  //              By "pure" we mean that these points are
  //              on the outer boundary and not also on
  //              an inner boundary.
  //              Here we fill in the innermost ghost zone
  //              layer first and move outward. At each
  //              layer, we fill in +/- x0 faces first,
  //              then +/- x1 faces, finally +/- x2 faces,
  //              filling in the edges as we go.
  apply_bcs_outerextrap_and_inner_only_specific_gfs__launcher(params, bcstruct, gfs, num_gfs, gfs_to_sync);

  ///////////////////////////////////////////////////////
  // STEP 2 of 2: Apply BCs to inner boundary points.
  //              These map to either the grid interior
  //              ("pure inner") or to pure outer boundary
  //              points ("inner maps to outer"). Those
  //              that map to outer require that outer be
  //              populated first; hence this being
  //              STEP 2 OF 2.
  apply_bcs_inner_only_specific_gfs(commondata, params, bcstruct, gfs, num_gfs, gfs_to_sync);
"""
    prefunc = generate_prefunc__apply_bcs_outerextrap_and_inner_only_specific_gfs()
    c_function.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
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
