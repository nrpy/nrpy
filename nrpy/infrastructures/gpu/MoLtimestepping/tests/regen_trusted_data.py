import nrpy.c_function as cfc, json

def regen_base_register_CFunction_MoL_step_forward_in_time():
    """Regenerate trusted data stored in JSON file for base_register_CFunction_MoL_step_forward_in_time"""
    from nrpy.infrastructures.gpu.MoLtimestepping.base_MoL import base_register_CFunction_MoL_step_forward_in_time, MoL_Functions_dict
    from nrpy.infrastructures.BHaH.MoLtimestepping.RK_Butcher_Table_Dictionary import (
        generate_Butcher_tables,
    )
    from nrpy.helpers.generic import compress_string_to_base64

    Butcher_dict = generate_Butcher_tables()

    expected_str_dict=dict()
    for k, Butcher_tuple in Butcher_dict.items():
        Butcher = Butcher_tuple[0]
        # Ignore adaptive Butcher tables
        if Butcher[-1][0] != "":
            continue

        # Reset stored C functions
        cfc.CFunction_dict.clear()
        MoL_Functions_dict.clear()

        MoLclass = base_register_CFunction_MoL_step_forward_in_time(Butcher_dict, k, rhs_string="rhs_eval(commondata, params, rfmstruct,  auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);", post_rhs_string="""if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
    apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);""")
        MoLclass.setup_gf_aliases()
        MoLclass.generate_RK_steps()
        MoLclass.register()
        expected_str_dict[k] = compress_string_to_base64(
            cfc.CFunction_dict["MoL_step_forward_in_time"].full_function
        )
    with open("nrpy/infrastructures/gpu/MoLtimestepping/tests/DOCTEST-base_register_CFunction_MoL_step_forward_in_time.json",'w') as f:
        f.write(json.dumps(expected_str_dict, sort_keys=True, indent=4))

def regen__cuda__register_CFunction_MoL_step_forward_in_time():
    """Regenerate trusted data stored in JSON file for base_register_CFunction_MoL_step_forward_in_time"""
    from nrpy.infrastructures.gpu.MoLtimestepping.base_MoL import MoL_Functions_dict
    from nrpy.infrastructures.gpu.MoLtimestepping.cuda.MoL import register_CFunction_MoL_step_forward_in_time
    from nrpy.infrastructures.BHaH.MoLtimestepping.RK_Butcher_Table_Dictionary import (
        generate_Butcher_tables,
    )
    from nrpy.helpers.generic import compress_string_to_base64

    Butcher_dict = generate_Butcher_tables()

    expected_str_dict=dict()
    for k, Butcher_tuple in Butcher_dict.items():
        Butcher = Butcher_tuple[0]
        # Ignore adaptive Butcher tables
        if Butcher[-1][0] != "":
            continue

        # Reset stored C functions
        cfc.CFunction_dict.clear()
        MoL_Functions_dict.clear()

        rhs_string = "rhs_eval(commondata, params, rfmstruct,  auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);"
        post_rhs_string=(
            "if (strncmp(commondata->outer_bc_type, \"extrapolation\", 50) == 0)\n"
            "  apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);"
        )

        _MoLclass = register_CFunction_MoL_step_forward_in_time(
            Butcher_dict,
            k,
            rhs_string=rhs_string,
            post_rhs_string=post_rhs_string
        )

        expected_str_dict[k] = compress_string_to_base64(
            cfc.CFunction_dict["MoL_step_forward_in_time"].full_function
        )
    with open("nrpy/infrastructures/gpu/MoLtimestepping/tests/DOCTEST-cuda__register_CFunction_MoL_step_forward_in_time.json",'w') as f:
        f.write(json.dumps(expected_str_dict, sort_keys=True, indent=4))

if __name__ == "__main__":
    regen_base_register_CFunction_MoL_step_forward_in_time()
    regen__cuda__register_CFunction_MoL_step_forward_in_time()