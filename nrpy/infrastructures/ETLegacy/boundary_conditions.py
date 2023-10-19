import nrpy.c_function as cfc


def add_to_Cfunction_dict_specify_NewRad_BoundaryConditions_parameters(
    thorn_name: str,
) -> None:
    """
    Set up NewRad boundary conditions for the given thorn.

    This function specifies NewRad boundary conditions for evolved gridfunctions
    for a thorn in the Einstein Toolkit.

    :param thorn_name: The name of the NRPy+ generated thorn for which to set up NewRad boundary conditions.
    :return: None
    :raises: None
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
    name = f"specify_NewRad_BoundaryConditions_parameters_{thorn_name}"
    params = "CCTK_ARGUMENTS"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
  DECLARE_CCTK_PARAMETERS;
"""

    for gf in evol_gfs_list:
        var_at_infinite_r = "0.0"
        var_char_speed = "1.0"
        for reg_gf in gri.glb_gridfcs_list:
            if reg_gf.name + "GF" == gf:
                var_at_infinite_r = str(reg_gf.f_infinity)
                var_char_speed = str(reg_gf.wavespeed)
        var_radpower = "1.0"
        if "aDD" in gf or "trK" in gf:  # consistent with Lean code.
            var_radpower = "2.0"

        body += f"  NewRad_Apply(cctkGH, {gf}, {gf.replace('GF', '')}_rhsGF, {var_at_infinite_r}, {var_char_speed}, {var_radpower});\n"

    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
        enableCparameters=False,
    )
