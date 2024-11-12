import json
import nrpy.helpers.gpu_kernels.cuda_utilities as gputils

def regen__gpu_kernel__CUDA_reductions__trusted_data():
    """Regenerate trusted data stored in JSON file for gpu_kernel.py"""
    from nrpy.helpers.generic import compress_string_to_base64

    expected_str_dict={}
    for fp_type in ['float', 'double']:
        expected_sub_str_dict={}
        
        for reduction_type, _reduction_func in gputils.implemented_reduction_dict.items():
            reduction = gputils.CUDA_reductions(reduction_type=reduction_type, fp_type=fp_type)
            reduction.generate_CFunction()

            # MoLclass.register_final_code()
            expected_sub_str_dict[reduction_type] = compress_string_to_base64(
                reduction.CFunction.full_function
            )
        expected_str_dict[fp_type] = expected_sub_str_dict
    with open("nrpy/tests/DOCTEST-gpu_kernel__CUDA_reductions.json",'w') as f:
        f.write(json.dumps(expected_str_dict, sort_keys=True, indent=4))
        
if __name__ == "__main__":
    regen__gpu_kernel__CUDA_reductions__trusted_data()
