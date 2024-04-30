import nrpy.c_function as cfc, json

def regen__openmp__register_CFunction_read_checkpoint():
    """Regenerate trusted data stored in JSON file for register_CFunction_read_checkpoint"""
    from nrpy.infrastructures.BHaH.checkpoints.openmp.checkpointing import register_CFunction_read_checkpoint
    from nrpy.helpers.generic import compress_string_to_base64

    expected_str_dict=dict()
    
    k = "default"
    _ = register_CFunction_read_checkpoint()
    expected_str_dict[k] = compress_string_to_base64(
        cfc.CFunction_dict["read_checkpoint"].full_function
    )
    with open("nrpy/infrastructures/BHaH/checkpoints/tests/DOCTEST-openmp__register_CFunction_read_checkpoint.json",'w') as f:
        f.write(json.dumps(expected_str_dict, sort_keys=True, indent=4))

def regen__openmp__register_CFunction_write_checkpoint():
    """Regenerate trusted data stored in JSON file for register_CFunction_write_checkpoint"""
    from nrpy.infrastructures.BHaH.checkpoints.openmp.checkpointing import register_CFunction_write_checkpoint
    from nrpy.helpers.generic import compress_string_to_base64

    expected_str_dict=dict()
    
    k = "default"
    _ = register_CFunction_write_checkpoint()
    expected_str_dict[k] = compress_string_to_base64(
        cfc.CFunction_dict["write_checkpoint"].full_function
    )
    with open("nrpy/infrastructures/BHaH/checkpoints/tests/DOCTEST-openmp__register_CFunction_write_checkpoint.json",'w') as f:
        f.write(json.dumps(expected_str_dict, sort_keys=True, indent=4))

if __name__ == "__main__":
    regen__openmp__register_CFunction_read_checkpoint()
    regen__openmp__register_CFunction_write_checkpoint()