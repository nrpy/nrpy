# nrpy/infrastructures/BHaH/MoLtimestepping/MoL_register_all.py
"""
A one-stop function to register *all* MoL C-functions in the correct order.
This includes memory allocation, free, and the main step-forward function.

Authors: Zachariah B. Etienne (lead maintainer)
         zachetie **at** gmail **dot* com
         Samuel Tootle (GPU support in NRPy2)
         Brandon Clark (original, NRPy1 version)
"""

import nrpy.params as par
from nrpy.infrastructures.BHaH import BHaH_defines_h, griddata_commondata
from nrpy.infrastructures.BHaH.MoLtimestepping.MoL_allocators import (
    register_CFunction_MoL_free_memory,
    register_CFunction_MoL_malloc,
)
from nrpy.infrastructures.BHaH.MoLtimestepping.MoL_gridfunction_names import (
    generate_gridfunction_names,
)
from nrpy.infrastructures.BHaH.MoLtimestepping.MoL_rk_substep import (
    check_supported_parallelization,
)
from nrpy.infrastructures.BHaH.MoLtimestepping.MoL_step_forward import (
    register_CFunction_MoL_step_forward_in_time,
)
from nrpy.infrastructures.BHaH.MoLtimestepping.RK_Butcher_Table_Dictionary import (
    generate_Butcher_tables,
)

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
    parallelization: str = "openmp",
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
    :param parallelization: Parameter to specify parallelization (openmp or cuda).
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
    #include "BHaH_function_prototypes.h"
    /**
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
    } // END FUNCTION MoL_free_memory_non_y_n_gfs
    <BLANKLINE>
    >>> print(cfc.CFunction_dict["MoL_malloc_non_y_n_gfs"].full_function) # doctest: +ELLIPSIS
    #include "BHaH_defines.h"
    #include "BHaH_function_prototypes.h"
    /**
     * Method of Lines (MoL) for "RK4" method: Allocate memory for "non_y_n_gfs" gridfunctions
     * - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep
     * - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method
     */
    void MoL_malloc_non_y_n_gfs(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                MoL_gridfunctions_struct *restrict gridfuncs) {
      const int Nxx_plus_2NGHOSTS_tot = params->Nxx_plus_2NGHOSTS0 * params->Nxx_plus_2NGHOSTS1 * params->Nxx_plus_2NGHOSTS2;
      gridfuncs->y_nplus1_running_total_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
      gridfuncs->k_odd_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
      gridfuncs->k_even_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
      if (NUM_AUXEVOL_GFS > 0)
        gridfuncs->auxevol_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_AUXEVOL_GFS * Nxx_plus_2NGHOSTS_tot);
    <BLANKLINE>
      gridfuncs->diagnostic_output_gfs = gridfuncs->y_nplus1_running_total_gfs;
      gridfuncs->diagnostic_output_gfs2 = gridfuncs->k_odd_gfs;
    } // END FUNCTION MoL_malloc_non_y_n_gfs
    <BLANKLINE>
    """
    check_supported_parallelization("register_CFunctions", parallelization)

    Butcher_dict = generate_Butcher_tables()

    # Step 1: Build all memory alloc and free:
    for which_gfs in ["y_n_gfs", "non_y_n_gfs"]:
        register_CFunction_MoL_malloc(
            Butcher_dict, MoL_method, which_gfs, parallelization=parallelization
        )
        register_CFunction_MoL_free_memory(
            Butcher_dict, MoL_method, which_gfs, parallelization=parallelization
        )

    # Step 2: Possibly register the main stepping function:
    if register_MoL_step_forward_in_time:
        register_CFunction_MoL_step_forward_in_time(
            Butcher_dict,
            MoL_method,
            rhs_string,
            post_rhs_string,
            post_post_rhs_string,
            enable_rfm_precompute=enable_rfm_precompute,
            enable_curviBCs=enable_curviBCs,
            enable_intrinsics=enable_intrinsics,
            parallelization=parallelization,
            rational_const_alias=rational_const_alias,
        )

    # Step 3: Register the struct in BHaH_defines_h:
    griddata_commondata.register_griddata_commondata(
        __name__, "MoL_gridfunctions_struct gridfuncs", "MoL gridfunctions"
    )

    (
        y_n_gridfunctions,
        non_y_n_gridfunctions_list,
        _diag_pt,
        _diag_pt2,
    ) = generate_gridfunction_names(Butcher_dict, MoL_method=MoL_method)

    BHaH_MoL_body: str = (
        "typedef struct __MoL_gridfunctions_struct__ {\n"
        + f"REAL *restrict {y_n_gridfunctions};\n"
        + "".join(f"REAL *restrict {gfs};\n" for gfs in non_y_n_gridfunctions_list)
        + r"""REAL *restrict diagnostic_output_gfs;
REAL *restrict diagnostic_output_gfs2;
} MoL_gridfunctions_struct;
"""
    )
    if parallelization != "openmp":
        BHaH_MoL_body = BHaH_MoL_body.replace("REAL *restrict ", "REAL * ")

    BHaH_defines_h.register_BHaH_defines(__name__, BHaH_MoL_body)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
