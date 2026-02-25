"""
Generates the C placeholder for the external numerical interpolation engine.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc

def placeholder_interpolation_engine(spacetime_name: str, PARTICLE: str) -> None:
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = f"""@brief Placeholder for the external batch-processing interpolation engine.
    Adapted for flattened SoA arrays. Passes the global batch capacity stride
    to underlying analytic workers with dual-architecture offloading.
    """

    name = f"placeholder_interpolation_engine_{spacetime_name}"
    params = """const commondata_struct *restrict commondata,
                const int num_photons,
                const int *restrict req_photon_ids,
                const double *restrict req_pos,
                double *restrict metric_g4DD,
                double *restrict conn_GammaUDD"""

    metric_worker_func = f"g4DD_metric_{spacetime_name}"
    conn_worker_func = f"connections_{spacetime_name}"

    body = f"""
    (void)req_photon_ids;
    
    #ifdef USE_GPU
        #pragma omp target teams distribute parallel for \
                    map(to: commondata[0:1]) \
                    is_device_ptr(req_photon_ids, req_pos, metric_g4DD, conn_GammaUDD)
    #else
        #pragma omp parallel for
    #endif
    for (int batch_id = 0; batch_id < num_photons; ++batch_id) {{
    
        {metric_worker_func}(commondata, req_pos, metric_g4DD, BUNDLE_CAPACITY, batch_id);

        if (conn_GammaUDD != NULL) {{
            {conn_worker_func}(commondata, req_pos, conn_GammaUDD, BUNDLE_CAPACITY, batch_id);
        }}
    }}
    """

    cfc.register_CFunction(
        includes=includes, desc=desc, name=name, params=params, body=body,
    )

if __name__ == "__main__":
    import logging
    import os
    import sys
    sys.path.append(os.getcwd())
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("TestPlaceholderInterpolationEngine")
    SPACETIME = "KerrSchild_Cartesian"
    PARTICLE = "photon"
    logger.info("Test: Generating Placeholder Interpolation Engine C-code...")
    try:
        placeholder_interpolation_engine(SPACETIME, PARTICLE)
        for func_name, c_function in cfc.CFunction_dict.items():
            filename = f"{func_name}.c"
            with open(filename, "w", encoding="utf-8") as f:
                f.write(c_function.full_function)
    except Exception as e:
        logger.error(" -> FAIL: placeholder_interpolation test failed: %s", e)
        sys.exit(1)