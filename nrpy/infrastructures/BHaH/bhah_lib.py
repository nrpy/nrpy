# nrpy/infrastructures/BHaH/bhah_lib.py
"""
Functions needed to generate bhah_lib, a NRPy2 library that interfaces with the unstructured mesh hydrodynamics code MANGA.

Author: Leonardo Rosa Werneck
        wernecklr **at** gmail **dot* com
"""

from typing import Dict

from nrpy.c_function import register_CFunction


def register_CFunction_bhah_initialize() -> None:
    """Register an interface function for initializing bhah_lib."""
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = "BlackHoles@Home library setup"
    cfunc_type = "BHaH_struct *"
    name = "bhah_initialize"
    params = ""
    include_CodeParameters_h = False
    body = """

// Step 1.a: Allocate memory for the BH@H struct
BHaH_struct *bhah_struct = (BHaH_struct *)malloc(sizeof(BHaH_struct));
commondata_struct *commondata = (commondata_struct *)malloc(sizeof(commondata_struct));

// Step 1.b: Set each commondata CodeParameter to default.
commondata_struct_set_to_default(commondata);

// Step 1.c: Allocate NUMGRIDS griddata arrays, each containing
//           data specific to an individual grid.
const int n_grids = commondata->NUMGRIDS;
griddata_struct *griddata = (griddata_struct *)malloc(sizeof(griddata_struct) * n_grids);

// Step 1.d: Set each CodeParameter in griddata.params to default.
params_struct_set_to_default(commondata, griddata);

// Step 1.e: Set up numerical grids: xx[3], masks, Nxx, dxx, invdxx, bcstruct,
//           rfm_precompute, timestep, etc. With calling_for_first_time = true,
//           commondata time is set to zero (as is the iteration, etc).
const bool calling_for_first_time = true;
numerical_grids_and_timestep(commondata, griddata, calling_for_first_time);

// Step 1.f: Allocate memory for all gridfunctions
for(int grid = 0; grid < n_grids; grid++) {
  MoL_malloc_y_n_gfs(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
  MoL_malloc_non_y_n_gfs(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
}

// Step 1.g: Set up initial data
initial_data(commondata, griddata);

// Step 1.h: Set pointers to griddata and commondata
bhah_struct->griddata = griddata;
bhah_struct->commondata = commondata;

return bhah_struct;
"""
    register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
    )


def register_CFunction_bhah_finalize() -> None:
    """Register an interface function for finalizing bhah_lib."""
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = "Finalize the BlackHoles@Home library"
    cfunc_type = "void"
    name = "bhah_finalize"
    params = "BHaH_struct *bhah_struct"
    include_CodeParameters_h = False
    body = """
    commondata_struct *commondata = bhah_struct->commondata;
    griddata_struct *griddata = bhah_struct->griddata;
    for(int grid = 0; grid < commondata->NUMGRIDS; grid++) {
      const bool free_auxevol_gfs_if_exist = true;
      MoL_free_memory_y_n_gfs(&griddata[grid].gridfuncs);
      MoL_free_memory_non_y_n_gfs(&griddata[grid].gridfuncs, free_auxevol_gfs_if_exist);
      rfm_precompute_free(commondata, &griddata[grid].params, griddata[grid].rfmstruct);
    }
    free(bhah_struct->commondata);
    free(bhah_struct->griddata);
    free(bhah_struct);
"""

    register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
    )


def register_CFunction_bhah_diagnostics() -> None:
    """Register an interface function for diagnostics."""
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = "Output BH@H diagnostics"
    cfunc_type = "void"
    name = "bhah_diagnostics"
    params = "BHaH_struct *bhah_struct"
    include_CodeParameters_h = False
    body = "diagnostics(bhah_struct->commondata, bhah_struct->griddata);"

    register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
    )


def register_CFunction_bhah_evolve() -> None:
    """Register an interface function for time evolution."""
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = "Perform spacetime evolution in BlackHoles@Home"
    cfunc_type = "void"
    name = "bhah_evolve"
    params = "BHaH_struct *bhah_struct"
    include_CodeParameters_h = False
    body = """
    commondata_struct *commondata = bhah_struct->commondata;
    griddata_struct *griddata = bhah_struct->griddata;
    while(commondata->time < commondata->t_final) {
      bhah_diagnostics(bhah_struct);
      MoL_step_forward_in_time(commondata, griddata);
    }
"""

    register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
    )


def register_CFunctions_bhah_lib() -> None:
    """Register all C functions required by the bhah_lib library."""
    register_CFunction_bhah_initialize()
    register_CFunction_bhah_evolve()
    register_CFunction_bhah_diagnostics()
    register_CFunction_bhah_finalize()


def supplemental_defines_dict_bhah_lib() -> Dict[str, str]:
    """
    Return the supplemental dictionary required by bhah_lib.

    :return: Dictionary of Infrastructure : Struct.
    """
    return {"BHaH Lib": """
typedef struct BHaH_struct {
  commondata_struct *commondata;
  griddata_struct *griddata;
} BHaH_struct;"""}
