"""
Registration of boundary condition functions.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
import nrpy.c_function as cfc
import nrpy.grid as gri


def register_CFunction_specify_Driver_BoundaryConditions(thorn_name: str) -> None:
    """
    Register C functions for specifying boundary conditions within a given thorn.

    This function generates and registers C code that sets the boundary conditions
    for both auxiliary and evolved grid functions within a specific thorn. The
    boundary conditions are specified based on the `Driver` method in the
    Einstein Toolkit.

    :param thorn_name: The name of the thorn for which to specify boundary conditions.
    :return: None
    """
    includes = [
        "stdio.h",
        "cctk.h",
        "cctk_Arguments.h",
        "cctk_Parameters.h",
        "cctk_Faces.h",
        "util_Table.h",
    ]

    desc = """
EVOL variables are set to use `none` boundary conditions, as these are set via NewRad.

Since we choose NewRad boundary conditions, we must register all
evolved gridfunctions to have boundary type "none". This is because
NewRad directly modifies the RHSs.

This code is based on Kranc's McLachlan/ML_BSSN/src/Boundaries.cc code.
"""

    c_type = "void"
    name = f"{thorn_name}_specify_Driver_BoundaryConditions"
    params = "CCTK_ARGUMENTS"

    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
DECLARE_CCTK_PARAMETERS;
CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
"""
    for gfname, gf in sorted(gri.glb_gridfcs_dict.items()):
        if gf.group == "EVOL":
            body += f"""
ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "{thorn_name}::{gfname}GF", "none");
if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for {thorn_name}::{gfname}GF!");
"""

    ET_schedule_bins_entries = [
        (
            "Driver_BoundarySelect",
            """
schedule FUNC_NAME in Driver_BoundarySelect
{
  LANG: C
  OPTIONS: LEVEL
} "Register boundary conditions in PreSync bin Driver_BoundarySelect."
""",
        ),
        (
            "MoL_PostStep",
            """
schedule FUNC_NAME in MoL_PostStep
{
  LANG: C
  OPTIONS: LEVEL
  SYNC: evol_variables
} "Dummy function to force AMR+interprocessor synchronization"
""",
        ),
    ]
    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=ET_schedule_bins_entries,
    )


def register_CFunction_specify_NewRad_BoundaryConditions_parameters(
    thorn_name: str,
) -> None:
    """
    Set up NewRad boundary conditions for the given thorn.

    This function specifies NewRad boundary conditions for evolved gridfunctions
    for a thorn in the Einstein Toolkit.

    :param thorn_name: The name of the NRPy+ generated thorn for which to set up NewRad boundary conditions.
    :return: None
    """
    includes = ["math.h", "cctk.h", "cctk_Arguments.h", "cctk_Parameters.h"]
    desc = """
Set up NewRad boundary conditions.
   As explained in lean_public/LeanBSSNMoL/src/calc_bssn_rhs.F90,
   the function NewRad_Apply takes the following arguments:
   NewRad_Apply(cctkGH, var, rhs, var0, v0, radpower),
     which implement the boundary condition:
       var  =  var_at_infinite_r + u(r-var_char_speed*t)/r^var_radpower
  Obviously for var_radpower>0, var_at_infinite_r is the value of
    the variable at r->infinity. var_char_speed is the propagation
    speed at the outer boundary, and var_radpower is the radial
    falloff rate.
"""
    c_type = "void"
    name = f"{thorn_name}_specify_NewRad_BoundaryConditions_parameters"
    params = "CCTK_ARGUMENTS"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
  DECLARE_CCTK_PARAMETERS;
"""
    for gfname, gf in gri.glb_gridfcs_dict.items():
        if gf.group == "EVOL" and isinstance(gf, gri.ETLegacyGridFunction):
            var_radpower = "1.0"
            body += f"  NewRad_Apply(cctkGH, {gfname}GF, {gfname}_rhsGF, {gf.f_infinity}, {gf.wavespeed}, {var_radpower});\n"

    ET_schedule_bins_entries = [
        (
            "MoL_CalcRHS",
            """
schedule FUNC_NAME in MoL_CalcRHS after BaikalETK_RHS
{
  LANG: C
  READS: evol_variables(everywhere)
  WRITES: evol_variables_rhs(boundary)
} "NewRad boundary conditions, scheduled right after RHS eval."
""",
        ),
        (
            "MoL_PseudoEvolution",
            """
# This schedule call is not required for PreSync but remains in the schedule for backward compatibility.
schedule GROUP ApplyBCs as WaveToy_auxgfs_ApplyBCs in MoL_PseudoEvolution after specify_BoundaryConditions
{
} "Apply boundary conditions"
""",
        ),
    ]
    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=ET_schedule_bins_entries,
    )


def register_CFunctions(thorn_name: str) -> None:
    """
    Register C functions related to boundary conditions for the given thorn.

    :param thorn_name: The name of the thorn for which to register C functions related to boundary conditions.
    :return: None
    """
    register_CFunction_specify_Driver_BoundaryConditions(thorn_name=thorn_name)
    register_CFunction_specify_NewRad_BoundaryConditions_parameters(
        thorn_name=thorn_name
    )
