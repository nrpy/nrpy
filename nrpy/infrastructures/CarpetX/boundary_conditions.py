"""
Registration of boundary condition functions.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_function as cfc
import nrpy.grid as gri


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
    includes = [
        "math.h",
        "cctk.h",
        "cctk_Arguments.h",
        "cctk_Parameters.h",
        "newradx.hxx",
    ]
    desc = """
  Set up NewRad boundary conditions.
  The function NewRadX_Apply takes the following arguments:
  NewRadX_Apply(cctkGH, var, rhs, var0, v0, radpower),
  which implement the boundary condition:
    var  =  var_at_infinite_r + u(r-var_char_speed*t)/r^var_radpower
  Obviously for var_radpower>0, var_at_infinite_r is the value of
  the variable at r->infinity. var_char_speed is the propagation
  speed at the outer boundary, and var_radpower is the radial
  falloff rate."""
    c_type = 'extern "C" void'
    name = f"{thorn_name}_specify_NewRad_BoundaryConditions_parameters"
    params = "CCTK_ARGUMENTS"
    body = f"""  DECLARE_CCTK_ARGUMENTSX_{name};
  DECLARE_CCTK_PARAMETERS;

  using namespace NewRadX;

"""
    for gfname, gf in gri.glb_gridfcs_dict.items():
        if gf.group == "EVOL" and isinstance(gf, gri.CarpetXGridFunction):
            var_radpower = "1.0"
            body += f"  NewRadX_Apply(cctkGH, {gfname}GF, {gfname}_rhsGF, {gf.f_infinity}, {gf.wavespeed}, {var_radpower});\n"

    ET_schedule_bins_entries = [
        (
            "ODESolvers_RHS",
            f"""
schedule FUNC_NAME in ODESolvers_RHS after {thorn_name}_RHS
{{
  LANG: C
  READS:  evol_variables(everywhere)
  WRITES: evol_variables_rhs(boundary)
}} "NewRad boundary conditions, scheduled right after RHS eval."
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
    register_CFunction_specify_NewRad_BoundaryConditions_parameters(
        thorn_name=thorn_name
    )
