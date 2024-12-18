"""
Register and handle error codes for BHaHAHA, including error message definitions.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List, Tuple

import nrpy.c_function as cfc
from nrpy.infrastructures.BHaH import griddata_commondata

# Define a list of (str, str) tuples representing error code macro names and messages:
error_code_msg_tuples_list: List[Tuple[str, str]] = []

# fmt: off
error_code_msg_tuples_list += [("BHAHAHA_SUCCESS", "Success.")]
error_code_msg_tuples_list += [("BHAHAHA_UNKNOWN_ERROR", "BHaHAHA: unknown error... developer did not add appropriate error handling (grumble grumble).")]
error_code_msg_tuples_list += [("FIND_HORIZON_GETTIMEOFDAY_BROKEN", "bah_find_horizon(): gettimeofday() returned an error exit code, indicating a non-POSIX-compatible or broken system. Sorry.")]
error_code_msg_tuples_list += [("FIND_HORIZON_WRONG_ANGULAR_RESOLUTIONS", "bah_find_horizon(): Did not find Ntheta[n_resolutions-1]=32 & Nphi[n_resolutions-1]=64. Please set them correctly.")]
error_code_msg_tuples_list += [("FIND_HORIZON_MAX_ITERATIONS_EXCEEDED", "bah_find_horizon(): maximum iterations exceeded. Set verbosity to 2, or increase max_iterations.")]
error_code_msg_tuples_list += [("FIND_HORIZON_HORIZON_TOO_SMALL", "bah_find_horizon(): Horizon radius < 3 * dr of BHaHAHA input grid. Set verbosity to 2, and/or try a smaller/higher-resolution grid.")]
error_code_msg_tuples_list += [("BCSTRUCT_EIGENCOORD_FAILURE", "bah_bcstruct_set_up(): problem setting up Eigen-coordinates.")]
error_code_msg_tuples_list += [("BCSTRUCT_SET_PARITY_ERROR", "bah_bcstruct_set_up(): problem computing parity conditions.")]
error_code_msg_tuples_list += [("INITIAL_DATA_MALLOC_ERROR", "bah_initial_data(): Failed to allocate memory to coarse_to_fine or dst_pts[][].")]
error_code_msg_tuples_list += [("NUMGRID_EXTERN_MALLOC_ERROR_GFS", "bah_numgrid__external_input_set_up(): Failed to allocate memory to external_input_gfs.")]
error_code_msg_tuples_list += [("NUMGRID_EXTERN_MALLOC_ERROR_RTHETAPHI", "bah_numgrid__external_input_set_up(): Failed to allocate memory to commondata->external_input_r_theta_phi[][].")]
error_code_msg_tuples_list += [("NUMGRID_INTERP_MALLOC_ERROR_GFS", "bah_numgrid__interp_src_set_up(): Failed to allocate memory to commondata->interp_src_gfs.")]
error_code_msg_tuples_list += [("NUMGRID_INTERP_MALLOC_ERROR_RTHETAPHI", "bah_numgrid__interp_src_set_up(): Failed to allocate memory to commondata->interp_src_r_theta_phi.")]
error_code_msg_tuples_list += [("INTERP1D_NULL_PTRS", "bah_interpolation_1d_radial_spokes_on_3d_src_grid(): Found NULL pointer(s); one or more input arrays not malloc'ed.")]
error_code_msg_tuples_list += [("INTERP1D_INTERP_ORDER_GT_NXX_PLUS_2NINTERPGHOSTS0", "bah_interpolation_1d_radial_spokes_on_3d_src_grid(): Interpolation order > Nxx0 + 2*NinterpGHOSTS.")]
error_code_msg_tuples_list += [("INTERP1D_HORIZON_TOO_LARGE", "bah_interpolation_1d_radial_spokes_on_3d_src_grid(): Horizon extends beyond input grid. Set verbosity to 2, and/or try a larger grid.")]
error_code_msg_tuples_list += [("INTERP1D_HORIZON_TOO_SMALL", "bah_interpolation_1d_radial_spokes_on_3d_src_grid(): Horizon radius dropped below r_min_external_input. Set verbosity to 2, and/or try a larger/higher-resolution grid.")]
error_code_msg_tuples_list += [("INTERP2D_EXT_TO_INTERPSRC_NULL_PTRS", "bah_interpolation_2d_external_input_to_interp_src_grid(): Found NULL pointer(s); one or more input arrays not malloc'ed.")]
error_code_msg_tuples_list += [("INTERP2D_EXT_TO_INTERPSRC_INTERP_ORDER_GT_NXX_PLUS_2NINTERPGHOSTS12", "bah_interpolation_2d_external_input_to_interp_src_grid(): Interpolation order > Nxx{1,2} + 2*NinterpGHOSTS.")]
error_code_msg_tuples_list += [("INTERP2D_EXT_TO_INTERPSRC_HORIZON_OUT_OF_BOUNDS", "bah_interpolation_2d_external_input_to_interp_src_grid(): Horizon extends beyond angular directions theta,phi. Should never happen...")]
error_code_msg_tuples_list += [("INTERP2D_GENERAL_NULL_PTRS", "bah_interpolation_2d_general__uniform_src_grid(): Found NULL pointer(s); one or more input arrays not malloc'ed.")]
error_code_msg_tuples_list += [("INTERP2D_GENERAL_INTERP_ORDER_GT_NXX_PLUS_2NGHOSTS12", "bah_interpolation_2d_general__uniform_src_grid(): Interpolation order > Nxx{1,2} + 2*NGHOSTS.")]
error_code_msg_tuples_list += [("INTERP2D_GENERAL_HORIZON_OUT_OF_BOUNDS", "bah_interpolation_2d_general__uniform_src_grid(): Horizon extends beyond angular directions theta,phi. Should never happen...")]
error_code_msg_tuples_list += [("DIAG_PROPER_CIRCUM_MALLOC_ERROR", "diagnostics_proper_circumferences(): One or more malloc's failed.")]
# fmt: on


def register_CFunction_error_message() -> None:
    """
    Register and handle C function error codes for BHaHAHA.

    DocTests:
        >>> register_CFunction_error_message()
    """
    griddata_commondata.register_griddata_commondata(
        __name__,
        "int error_flag",
        "Enables subroutines to pass error flags to parent routines.",
        is_commondata=True,
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
 Function: bah_error_handling()

 Description:
  - Driver function for BHaHAHA error reporting, including when horizon not found!
  - This function interprets error messages from throughout BHaHAHA.

 Parameter:
  - error_code - error code from BHaHAHA.

 Returns:
  - Error message string.
"""
    cfunc_type = "const char *"
    name = "error_message"
    params = "const bhahaha_error_codes error_code"
    body = "  switch (error_code) {\n"
    for item in error_code_msg_tuples_list:
        body += f"    case {item[0]}:\n"
        body += f'        return "{item[1]}";\n'
    body += """  }
    fprintf(stderr, "BHaHAHA error code %d not defined!\\n", error_code);
    return NULL;
"""
    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
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
