# nrpy/infrastructures/BHaH/MoLtimestepping/register_all.py
"""
Registers all C code components for Method of Lines (MoL) time integration.

This module provides a "one-stop shop" for setting up the time-stepping
framework within the NRPy BHaH infrastructure. It uses the Method of Lines
approach, where spatial derivatives are discretized first, yielding a large
system of ordinary differential equations (ODEs) that are then solved
numerically.

The main function, `register_CFunctions`, generates and registers all
necessary C components, including:
- The core `MoL_step_forward_in_time` function, which evolves gridfunctions
  by one timestep using a specified algorithm (e.g., Runge-Kutta).
- Memory allocation and deallocation routines for the temporary storage
  (e.g., Runge-Kutta stages) required by the chosen integrator.
- The `MoL_gridfunctions_struct` to manage all MoL-related data.
- C preprocessor definitions related to the MoL algorithm.

Authors: Zachariah B. Etienne (lead maintainer)
         zachetie **at** gmail **dot* com
         Samuel Tootle (GPU support in NRPy2)
         Brandon Clark (original, NRPy1 version)
"""
import nrpy.params as par
from nrpy.infrastructures.BHaH import MoLtimestepping, griddata_commondata

# fmt: off
_ = par.CodeParameter("int", __name__, "nn_0", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("int", __name__, "nn", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("REAL", __name__, "CFL_FACTOR", 0.5, commondata=True)
_ = par.CodeParameter("REAL", __name__, "dt", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("REAL", __name__, "t_0", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("REAL", __name__, "time", add_to_parfile=False, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("REAL", __name__, "t_final", 10.0, commondata=True)


# fmt: on


def register_CFunctions(
    MoL_method: str = "RK4",
    rhs_string: str = "rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, RK_INPUT_GFS, RK_OUTPUT_GFS);",
    post_rhs_string: str = "apply_bcs(Nxx, Nxx_plus_2NGHOSTS, RK_OUTPUT_GFS);",
    post_post_rhs_string: str = "",
    enable_rfm_precompute: bool = False,
    enable_curviBCs: bool = False,
    enable_intrinsics: bool = False,
    register_MoL_step_forward_in_time: bool = True,
    rational_const_alias: str = "const",
) -> None:
    r"""
    Register all MoL C functions and NRPy basic defines.

    :param MoL_method: The method to be used for MoL. Default is 'RK4'.
    :param rhs_string: RHS function call as string. Default is "rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, RK_INPUT_GFS, RK_OUTPUT_GFS);"
    :param post_rhs_string: Post-RHS function call as string. Default is "apply_bcs(Nxx, Nxx_plus_2NGHOSTS, RK_OUTPUT_GFS);"
    :param post_post_rhs_string: Post-post-RHS function call as string. Default is an empty string.
    :param enable_rfm_precompute: Enable reference metric support. Default is False.
    :param enable_curviBCs: Enable curvilinear boundary conditions. Default is False.
    :param enable_intrinsics: A flag to specify if hardware instructions should be used.
    :param register_MoL_step_forward_in_time: Whether to register the MoL step-forward function. Default is True.
    :param rational_const_alias: Overload const specifier for rational definitions

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.c_function as cfc
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunctions()
    >>> validate_strings(cfc.CFunction_dict["MoL_step_forward_in_time"].full_function, "MoL_step_forward_in_time")
    >>> sorted(cfc.CFunction_dict.keys())
    ['MoL_free_memory_non_y_n_gfs', 'MoL_free_memory_y_n_gfs', 'MoL_malloc_non_y_n_gfs', 'MoL_malloc_y_n_gfs', 'MoL_step_forward_in_time']
    >>> print(cfc.CFunction_dict["MoL_free_memory_non_y_n_gfs"].full_function) # doctest: +ELLIPSIS
    #include "BHaH_defines.h"
    /**
     * Method of Lines (MoL) for "RK4" method: Free memory for "non_y_n_gfs" gridfunctions
     *    - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep
     *    - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method
     *
     */
    void MoL_free_memory_non_y_n_gfs(MoL_gridfunctions_struct *restrict gridfuncs) {
      BHAH_FREE(gridfuncs->y_nplus1_running_total_gfs);
      BHAH_FREE(gridfuncs->k_odd_gfs);
      BHAH_FREE(gridfuncs->k_even_gfs);
      if (NUM_AUXEVOL_GFS > 0)
        BHAH_FREE(gridfuncs->auxevol_gfs);
    } // END FUNCTION MoL_free_memory_non_y_n_gfs
    <BLANKLINE>
    >>> print(cfc.CFunction_dict["MoL_malloc_non_y_n_gfs"].full_function) # doctest: +ELLIPSIS
    #include "BHaH_defines.h"
    /**
     * Method of Lines (MoL) for "RK4" method: Allocate memory for "non_y_n_gfs" gridfunctions
     *    - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep
     *    - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method
     */
    void MoL_malloc_non_y_n_gfs(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                MoL_gridfunctions_struct *restrict gridfuncs) {
      const int Nxx_plus_2NGHOSTS_tot = params->Nxx_plus_2NGHOSTS0 * params->Nxx_plus_2NGHOSTS1 * params->Nxx_plus_2NGHOSTS2;
      BHAH_MALLOC(gridfuncs->y_nplus1_running_total_gfs, sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
      BHAH_MALLOC(gridfuncs->k_odd_gfs, sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
      BHAH_MALLOC(gridfuncs->k_even_gfs, sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
      if (NUM_AUXEVOL_GFS > 0)
        BHAH_MALLOC(gridfuncs->auxevol_gfs, sizeof(REAL) * NUM_AUXEVOL_GFS * Nxx_plus_2NGHOSTS_tot);
    <BLANKLINE>
      gridfuncs->diagnostic_output_gfs = gridfuncs->y_nplus1_running_total_gfs;
      gridfuncs->diagnostic_output_gfs2 = gridfuncs->k_odd_gfs;
    } // END FUNCTION MoL_malloc_non_y_n_gfs
    <BLANKLINE>
    """
    MoLtimestepping.rk_substep.check_supported_parallelization("register_CFunctions")

    Butcher_dict = MoLtimestepping.rk_butcher_table_dictionary.generate_Butcher_tables()

    # Step 1: Build all memory alloc and free:
    for which_gfs in ["y_n_gfs", "non_y_n_gfs"]:
        MoLtimestepping.allocators.register_CFunction_MoL_malloc(
            Butcher_dict, MoL_method, which_gfs
        )
        MoLtimestepping.allocators.register_CFunction_MoL_free_memory(
            Butcher_dict, MoL_method, which_gfs
        )

    # Step 2: Possibly register the main stepping function:
    if register_MoL_step_forward_in_time:
        MoLtimestepping.step_forward.register_CFunction_MoL_step_forward_in_time(
            Butcher_dict,
            MoL_method,
            rhs_string,
            post_rhs_string,
            post_post_rhs_string,
            enable_rfm_precompute=enable_rfm_precompute,
            enable_curviBCs=enable_curviBCs,
            enable_intrinsics=enable_intrinsics,
            rational_const_alias=rational_const_alias,
        )

    # Step 3: Register the struct in BHaH_defines_h:
    griddata_commondata.register_griddata_commondata(
        __name__, "MoL_gridfunctions_struct gridfuncs", "MoL gridfunctions"
    )

    MoLtimestepping.BHaH_defines.register_BHaH_defines_h(Butcher_dict, MoL_method)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
