"""
Generates the C function to finalize photon blueprint data.

Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
This module extracts the final state from the SoA and populates the 
results buffer, specifically handling the celestial sphere mapping.
"""

import nrpy.c_function as cfc

def calculate_and_fill_blueprint_data_universal() -> None:
    """
    Generate the C function for universal blueprint data finalization.

    This function is called after the integration loop finishes to ensure 
    all termination conditions (Escape, Failure, etc.) are correctly logged.
    """
    # 1. Define C-Function metadata
    name = "calculate_and_fill_blueprint_data_universal"
    cfunc_type = "void"
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r"""@brief Finalizes the blueprint data for a single photon trajectory.
    
    Algorithm:
    1. Synchronizes the termination status from the master SoA.
    2. For escaped photons, projects the final 3D position onto a celestial sphere $(\theta, \phi)$."""

    params = """const PhotonStateSoA *restrict all_photons, const long int num_rays, 
                const long int photon_idx, blueprint_data_t *restrict result""" 
    
    # 2. Build the C body
    body = r"""
    // result->y_w, z_w, y_s, and z_s are persistent and were populated during event detection.
    result->termination_type = all_photons->status[photon_idx]; // Synchronize final exit status.

    // === Termination Dispatch: Celestial Sphere (Escape) ===
    if (all_photons->status[photon_idx] == TERMINATION_TYPE_CELESTIAL_SPHERE) {
        // Unpack final coordinates from the master state vector f.
        const double x_final = all_photons->f[IDX_GLOBAL(1, photon_idx, num_rays)]; // Final $x$ at escape.
        const double y_final = all_photons->f[IDX_GLOBAL(2, photon_idx, num_rays)]; // Final $y$ at escape.
        const double z_final = all_photons->f[IDX_GLOBAL(3, photon_idx, num_rays)]; // Final $z$ at escape.
        
        const double r_final = sqrt(x_final*x_final + y_final*y_final + z_final*z_final); // Final radial distance.
        
        if (r_final > 1e-9) {
            result->final_theta = acos(z_final / r_final); // Polar angle $\theta$ relative to the $z$-axis.
            result->final_phi = atan2(y_final, x_final); // Azimuthal angle $\phi$ in the $xy$-plane.
        }
    }
    """
    
    prefunc = "#ifdef USE_GPU\n#pragma omp declare target\n#endif"
    postfunc = "#ifdef USE_GPU\n#pragma omp end declare target\n#endif"
    
    cfc.register_CFunction(
        prefunc=prefunc, includes=includes, desc=desc, cfunc_type=cfunc_type, 
        name=name, params=params, body=body, postfunc=postfunc  
    )