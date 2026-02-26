import nrpy.c_function as cfc

def calculate_and_fill_blueprint_data_universal() -> None:
    name = "calculate_and_fill_blueprint_data_universal"
    cfunc_type = "void" # Pointer-based update
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "@brief Updates a photon's final state in the buffer without overwriting initial crossing data."

    params = """const PhotonStateSoA *restrict all_photons, 
                const long int num_rays, 
                const long int photon_idx,
                blueprint_data_t *restrict result""" 
    body = r"""
    // result->y_w, z_w, y_s, and z_s already contain the crossing points from the integrator!
    result->termination_type = all_photons->status[photon_idx];

    // --- Termination Dispatch ---
    // NOTE: We intentionally DO NOT re-call handle_source_plane_intersection here 
    // to prevent overwriting the valid data stored during the active integration loop.

    if (all_photons->status[photon_idx] == TERMINATION_TYPE_CELESTIAL_SPHERE) {
        const double x = all_photons->f[IDX_GLOBAL(1, photon_idx, num_rays)];
        const double y = all_photons->f[IDX_GLOBAL(2, photon_idx, num_rays)];
        const double z = all_photons->f[IDX_GLOBAL(3, photon_idx, num_rays)];
        const double r = sqrt((x*x) + (y*y) + (z*z));
        if (r > 1e-9) {
            result->final_theta = acos(z / r);
            result->final_phi = atan2(y, x);
        }
    }
    """
    prefunc = """
    #ifdef USE_GPU
    #pragma omp declare target
    #endif
    """
    
    postfunc = """
    #ifdef USE_GPU
    #pragma omp end declare target
    #endif
    """
    
    # Step 6: Register the C function
    cfc.register_CFunction(
        prefunc=prefunc,      
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
        postfunc=postfunc  
    )