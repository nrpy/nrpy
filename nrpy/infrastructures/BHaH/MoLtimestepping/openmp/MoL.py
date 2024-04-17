"""
Module for producing C codes related to MoL timestepping within the BHaH infrastructure.
This includes implementation details and functions for allocating and deallocating the necessary memory.
This modular is specifically focused on utilizing OpenMP parallelization when generating code

Authors: Brandon Clark
         Zachariah B. Etienne (maintainer)
         zachetie **at** gmail **dot** com
         Samuel D. Tootle
         sdtootle **at** gmail **dot** com
"""

from typing import List, Union, Dict, Tuple
import os  # Standard Python module for multiplatform OS-level functions
import sympy as sp  # Import SymPy, a computer algebra system written entirely in Python
import nrpy.params as par  # NRPy+: Parameter interface
import nrpy.c_function as cfc
from nrpy.infrastructures.BHaH.MoLtimestepping import base_MoL
from nrpy.infrastructures.BHaH import griddata_commondata
from nrpy.infrastructures.BHaH import BHaH_defines_h

# fmt: off
_ = par.CodeParameter("int", __name__, "nn_0", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("int", __name__, "nn", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("REAL", __name__, "CFL_FACTOR", 0.5, commondata=True)
_ = par.CodeParameter("REAL", __name__, "dt", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("REAL", __name__, "t_0", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("REAL", __name__, "time", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("REAL", __name__, "t_final", 10.0, commondata=True)
# fmt: on

# Update core_modules to use correct key for ordering
for i, key in enumerate(BHaH_defines_h.core_modules_list):
    if "nrpy.infrastructures.BHaH.MoLtimestepping" in key:
        BHaH_defines_h.core_modules_list[i] = (
            "nrpy.infrastructures.BHaH.MoLtimestepping.openmp.MoL"
        )


class register_CFunction_MoL_malloc(base_MoL.base_register_CFunction_MoL_malloc):
    """
    Register MoL_malloc_y_n_gfs() and MoL_malloc_non_y_n_gfs(), allocating memory for the gridfunctions indicated.

    :param Butcher_dict: Dictionary of Butcher tables for the MoL method.
    :param MoL_method: Method for the Method of Lines.
    :param which_gfs: Specifies which gridfunctions to consider ("y_n_gfs" or "non_y_n_gfs").

    :raises ValueError: If the which_gfs parameter is neither "y_n_gfs" nor "non_y_n_gfs".

    Doctest: FIXME
    # >>> register_CFunction_MoL_malloc("Euler", "y_n_gfs")
    """

    def __init__(
        self,
        Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
        MoL_method: str,
        which_gfs: str,
    ) -> None:

        super().__init__(Butcher_dict, MoL_method, which_gfs)

        # Generate the body of the function

        for gridfunctions in self.gridfunctions_list:
            num_gfs = (
                "NUM_EVOL_GFS" if gridfunctions != "auxevol_gfs" else "NUM_AUXEVOL_GFS"
            )
            # Don't malloc a zero-sized array.
            if num_gfs == "NUM_AUXEVOL_GFS":
                self.body += "  if(NUM_AUXEVOL_GFS > 0) "
            self.body += (
                f"gridfuncs->{gridfunctions} = (REAL *restrict)malloc(sizeof(REAL) * {num_gfs} * "
                "Nxx_plus_2NGHOSTS_tot);\n"
            )

        self.body += f"\ngridfuncs->diagnostic_output_gfs  = gridfuncs->{self.diagnostic_gridfunctions_point_to};\n"
        self.body += f"gridfuncs->diagnostic_output_gfs2 = gridfuncs->{self.diagnostic_gridfunctions2_point_to};\n"

        cfc.register_CFunction(
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=True,
            body=self.body,
        )


########################################################################################################################
# EXAMPLE
# ODE: y' = f(t,y), y(t_0) = y_0
# Starting at time t_n with solution having value y_n and trying to update to y_nplus1 with timestep dt

# Example of scheme for RK4 with k_1, k_2, k_3, k_4 (Using non-diagonal algorithm) Notice this requires storage of
# y_n, y_nplus1, k_1 through k_4

# k_1      = dt*f(t_n, y_n)
# k_2      = dt*f(t_n + 1/2*dt, y_n + 1/2*k_1)
# k_3      = dt*f(t_n + 1/2*dt, y_n + 1/2*k_2)
# k_4      = dt*f(t_n + dt, y_n + k_3)
# y_nplus1 = y_n + 1/3k_1 + 1/6k_2 + 1/6k_3 + 1/3k_4

# Example of scheme RK4 using only k_odd and k_even (Diagonal algroithm) Notice that this only requires storage


# k_odd     = dt*f(t_n, y_n)
# y_nplus1  = 1/3*k_odd
# k_even    = dt*f(t_n + 1/2*dt, y_n + 1/2*k_odd)
# y_nplus1 += 1/6*k_even
# k_odd     = dt*f(t_n + 1/2*dt, y_n + 1/2*k_even)
# y_nplus1 += 1/6*k_odd
# k_even    = dt*f(t_n + dt, y_n + k_odd)
# y_nplus1 += 1/3*k_even
########################################################################################################################
class register_CFunction_MoL_step_forward_in_time(
    base_MoL.base_register_CFunction_MoL_step_forward_in_time
):
    """
    Register MoL_step_forward_in_time() C function, which is the core driver for time evolution in BHaH codes.

    :param Butcher_dict: A dictionary containing the Butcher tables for various RK-like methods.
    :param MoL_method: The method of lines (MoL) used for time-stepping.
    :param rhs_string: Right-hand side string of the C code.
    :param post_rhs_string: Input string for post-RHS phase in the C code.
    :param post_post_rhs_string: String to be used after the post-RHS phase.
    :param enable_rfm_precompute: Flag to enable reference metric functionality.
    :param enable_curviBCs: Flag to enable curvilinear boundary conditions.
    :param enable_simd: Flag to enable SIMD functionality.
    :param fp_type: Floating point type, e.g., "double".

    Doctest:
    # FIXME
    """

    def __init__(
        self,
        Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
        MoL_method: str,
        rhs_string: str = "",
        post_rhs_string: str = "",
        post_post_rhs_string: str = "",
        enable_rfm_precompute: bool = False,
        enable_curviBCs: bool = False,
        enable_simd: bool = False,
        fp_type: str = "double",
    ) -> None:

        super().__init__(
            Butcher_dict,
            MoL_method,
            rhs_string=rhs_string,
            post_rhs_string=post_rhs_string,
            post_post_rhs_string=post_post_rhs_string,
            enable_rfm_precompute=enable_rfm_precompute,
            enable_curviBCs=enable_curviBCs,
            enable_simd=enable_simd,
            fp_type=fp_type,
        )
        if enable_simd:
            self.includes += [os.path.join("SIMD", "SIMD_intrinsics.h")]
        self.setup_gf_aliases()
        self.generate_RK_steps()
        self.register_final_code()


# register_CFunction_MoL_free_memory() registers
#           MoL_free_memory_y_n_gfs() and
#           MoL_free_memory_non_y_n_gfs(), which free memory for
#           the indicated sets of gridfunctions
class register_CFunction_MoL_free_memory(
    base_MoL.base_register_CFunction_MoL_free_memory
):
    """
    Free memory for the specified Method of Lines (MoL) gridfunctions, given an MoL_method.

    :param Butcher_dict: Dictionary containing Butcher tableau for MoL methods.
    :param MoL_method: The Method of Lines method.
    :param which_gfs: The gridfunctions to be freed, either 'y_n_gfs' or 'non_y_n_gfs'.
    """

    def __init__(
        self,
        Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
        MoL_method: str,
        which_gfs: str,
    ) -> None:

        super().__init__(Butcher_dict, MoL_method, which_gfs)
        for gridfunction in self.gridfunctions_list:
            # Don't free a zero-sized array.
            if gridfunction == "auxevol_gfs":
                self.body += (
                    f"  if(NUM_AUXEVOL_GFS > 0) free(gridfuncs->{gridfunction});"
                )
            else:
                self.body += f"  free(gridfuncs->{gridfunction});"
        cfc.register_CFunction(
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            body=self.body,
        )


# Register all the CFunctions and NRPy basic defines
class register_CFunctions(base_MoL.base_register_CFunctions):
    r"""
    Register all MoL C functions and NRPy basic defines.

    :param MoL_method: The method to be used for MoL. Default is 'RK4'.
    :param rhs_string: RHS function call as string. Default is "rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, RK_INPUT_GFS, RK_OUTPUT_GFS);"
    :param post_rhs_string: Post-RHS function call as string. Default is "apply_bcs(Nxx, Nxx_plus_2NGHOSTS, RK_OUTPUT_GFS);"
    :param post_post_rhs_string: Post-post-RHS function call as string. Default is an empty string.
    :param enable_rfm_precompute: Enable reference metric support. Default is False.
    :param enable_curviBCs: Enable curvilinear boundary conditions. Default is False.
    :param enable_simd: Enable Single Instruction, Multiple Data (SIMD). Default is False.
    :param register_MoL_step_forward_in_time: Whether to register the MoL step forward function. Default is True.
    :param fp_type: Floating point type, e.g., "double".

    Doctests:
    >>> from nrpy.helpers.generic import compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> cfc.CFunction_dict.clear()
    >>> _ = register_CFunctions()
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4CNMBHNdABGaScZHDxOiAHcc/CXr2duHb+UiUyv83OVALtvJ+o7uK/PoSVGe7rPuvil8asOnzsX/43MrS1ReBRFZgQkawI+ZAWFEbmoUaaqYpfoO3aZEtuq8ccmII6pcpTvJSXJAdqsfFHNdhI+ukUE4ucIC0Aud6Kiq7SYm+yXe7fWO2X7YgBtJVshZo9Zw7wYSXfL/ysxDXrqrCktTcTM8HtIxOwos9sw2I89sGN1LrstSEtNhf0HqVtkNCCKpMerpP8FgrhYWA1dwTTbAvMguFRUcdNJGzS3g1ZaJGX7G5jYySoHU0pJ2d+pFJESKVStj0gkUxweXpmt8+fNlgI6NElWpm2V5ABWA9tfHnNBSFmrepJYtJ0N/6zPMXrc1lpYPN0WdRnqFyjJ19zwYCvkDovokVVSj/SO+3B8YVQ6xPkqE5lIqVZikjXHN4sEwTDLVlHppA3M0DKuYrcuDDgUSx/Up2AvFXszeBvWDXEGaVKd9YAeIBHMo980fgOQZGaCpuAl/ZT6RsUc9Z6ci1Rg8FTzBPG+1/NwpppD7Oaamj+kkAqift7YRzBHigxeKAy0talLQzeSZMzJcVNMm7172U/nK2Cbl5zWU4fgSt80zSayfFq/iIHX5vETgZUcGTzRh6WdCcuLt/GmjD7kKsk0qfDo9JGeXvpPDlmmgQLwwnsCl/jCYSFtv0Uhd9QeESXjyDWEUChr/jh/YW+ySr0vi1kSumG0w7IXNP/B35qR62ioxCMmAxUMFE86RYIWqlxcaUx1ESF4rnZ4Sn7Gqs31nt4WLHETh+GDb4cLqwKDSjx9UvYxP74Q3erRrArgpaRN2slvlvuc7Pmm06H/Su6yWcpYujcNR85IF6n99Jy3n43BKGd5IfdjbOdxgoSBBsdm2O9+D9RCYX4qKj3d7LSrafhx/KtyLQlGin1wxB3VWDDpqb+N1TfPEiLYH4e53hMH/DpNOlUz1fEE0MI9MlztNoJ/ThVkMY5zX+gaQXqEAErIBggd7JgEF1V9A+G12tDIu693TondIu19evcc2Xg+tDiSLyjERKCgbyFbtbI2AFm6evfLhRSPm/KIfJbZWPaZ/BDKXn/vavhTYMsr7zdmz0uOassSZdqc8QbLDiyDjcjOXGmROVaytkzptGob+Dg/wEmRsvkGUsEegv0qjROxGygIr8Mj1Ykh66f4edb5IQxt+hV1p2qVsdNDChxwkTdU6MTQyzZEXZ/itnE4BQ8yvwNS+BvNetkfnnGIiVRwEq6gn7TVjarI2fKnqO4AC8N1tWvTJQxDTpWXECv6z9T4p7BZggk4KebwN6CgiSNseoFUELlMPpAu6IgeeorOzO5B9uYiLdnSHKmMk0nqbxYBu7aWXC4HUyCq5rMv4V8cfSAywfdf8+nayTOEUuW0myS+aTKqVaWT8EouiMzAY2cLFFQKQp97iyKLVpRlIK67JX9TyFIwpGAclLEm8kGcfURxMEYjKtyC5kr99Tsm3Ln92ZggPaq8zi9jjSwtj+4t0227tAMQdYem3AABRp+aKG1DUjAABjwnNRgAACNeT0bHEZ/sCAAAAAARZWg==")
    >>> returned_string = cfc.CFunction_dict["MoL_step_forward_in_time"].full_function
    >>> if returned_string != expected_string:
    ...    compressed_str = compress_string_to_base64(returned_string)
    ...    error_message = "Trusted MoL_step_forward_in_time.full_function string changed!\n Here's the diff:\n"
    ...    error_message += "Here's the diff:\n" + diff_strings(expected_string, returned_string) + "\n"
    ...    raise ValueError(error_message + f"base64-encoded output: {compressed_str}")
    >>> sorted(cfc.CFunction_dict.keys())
    ['MoL_free_memory_non_y_n_gfs', 'MoL_free_memory_y_n_gfs', 'MoL_malloc_non_y_n_gfs', 'MoL_malloc_y_n_gfs', 'MoL_step_forward_in_time']
    >>> print(cfc.CFunction_dict["MoL_free_memory_non_y_n_gfs"].full_function)
    #include "BHaH_defines.h"
    #include "BHaH_function_prototypes.h"
    /*
     * Method of Lines (MoL) for "RK4" method: Free memory for "non_y_n_gfs" gridfunctions
     * - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep
     * - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method
     *
     */
    void MoL_free_memory_non_y_n_gfs(MoL_gridfunctions_struct *restrict gridfuncs) {
      free(gridfuncs->y_nplus1_running_total_gfs);
      free(gridfuncs->k_odd_gfs);
      free(gridfuncs->k_even_gfs);
      if (NUM_AUXEVOL_GFS > 0)
        free(gridfuncs->auxevol_gfs);
    }
    <BLANKLINE>
    >>> print(cfc.CFunction_dict["MoL_malloc_non_y_n_gfs"].full_function)
    #include "BHaH_defines.h"
    #include "BHaH_function_prototypes.h"
    /*
     * Method of Lines (MoL) for "RK4" method: Allocate memory for "non_y_n_gfs" gridfunctions
     * - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep
     * - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method
     */
    void MoL_malloc_non_y_n_gfs(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                MoL_gridfunctions_struct *restrict gridfuncs) {
    #include "set_CodeParameters.h"
      const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
      gridfuncs->y_nplus1_running_total_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
      gridfuncs->k_odd_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
      gridfuncs->k_even_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
      if (NUM_AUXEVOL_GFS > 0)
        gridfuncs->auxevol_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_AUXEVOL_GFS * Nxx_plus_2NGHOSTS_tot);
    <BLANKLINE>
      gridfuncs->diagnostic_output_gfs = gridfuncs->y_nplus1_running_total_gfs;
      gridfuncs->diagnostic_output_gfs2 = gridfuncs->k_odd_gfs;
    }
    <BLANKLINE>
    """

    def __init__(
        self,
        MoL_method: str = "RK4",
        rhs_string: str = "rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, RK_INPUT_GFS, RK_OUTPUT_GFS);",
        post_rhs_string: str = "apply_bcs(Nxx, Nxx_plus_2NGHOSTS, RK_OUTPUT_GFS);",
        post_post_rhs_string: str = "",
        enable_rfm_precompute: bool = False,
        enable_curviBCs: bool = False,
        enable_simd: bool = False,
        register_MoL_step_forward_in_time: bool = True,
        fp_type: str = "double",
    ) -> None:
        super().__init__(
            MoL_method=MoL_method,
            rhs_string=rhs_string,
            post_rhs_string=post_rhs_string,
            post_post_rhs_string=post_post_rhs_string,
            enable_rfm_precompute=enable_rfm_precompute,
            enable_curviBCs=enable_curviBCs,
            register_MoL_step_forward_in_time=register_MoL_step_forward_in_time,
            fp_type=fp_type,
        )
        for which_gfs in ["y_n_gfs", "non_y_n_gfs"]:
            register_CFunction_MoL_malloc(self.Butcher_dict, MoL_method, which_gfs)
            register_CFunction_MoL_free_memory(self.Butcher_dict, MoL_method, which_gfs)
        if register_MoL_step_forward_in_time:
            register_CFunction_MoL_step_forward_in_time(
                self.Butcher_dict,
                self.MoL_method,
                self.rhs_string,
                post_rhs_string=post_rhs_string,
                post_post_rhs_string=post_post_rhs_string,
                enable_rfm_precompute=self.enable_rfm_precompute,
                enable_curviBCs=self.enable_curviBCs,
                enable_simd=enable_simd,
                fp_type=self.fp_type,
            )

        griddata_commondata.register_griddata_commondata(
            __name__, "MoL_gridfunctions_struct gridfuncs", "MoL gridfunctions"
        )

        # Add OpenMP specific Loop
        self.BHaH_MoL_body += """
#define LOOP_ALL_GFS_GPS(ii) \
_Pragma("omp parallel for") \
  for(int (ii)=0;(ii)<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;(ii)++)
"""
        BHaH_defines_h.register_BHaH_defines(__name__, self.BHaH_MoL_body)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
