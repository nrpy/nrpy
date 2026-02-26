"""
Register C function for event finding.
Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
Uses LDG (Read-Only Data Cache) for geometry and direct SoA indexing.
"""
import nrpy.c_function as cfc

def find_event_time_and_state() -> None:
    includes = ["BHaH_defines.h", "<math.h>"]
    desc = "@brief Portable high-performance second-order root-finding."
    name = "find_event_time_and_state"
    
    params = """PhotonStateSoA *restrict all_photons, 
                const long int num_rays, 
                const long int photon_idx, 
                const double *restrict normal, 
                const double dist, 
                const event_type_t event_type"""

    body = r"""
    #define GET_COMP(soa_ptr, comp) (soa_ptr[IDX_GLOBAL(comp, photon_idx, num_rays)])
    #define PLANE_EVAL(soa_ptr) (GET_COMP(soa_ptr, 1)*normal[0] + GET_COMP(soa_ptr, 2)*normal[1] + GET_COMP(soa_ptr, 3)*normal[2] - dist)

    const double f0 = PLANE_EVAL(all_photons->f_p_p);
    const double f1 = PLANE_EVAL(all_photons->f_p);
    const double f2 = PLANE_EVAL(all_photons->f);

    const double t0 = all_photons->affine_param_p_p[photon_idx];
    const double t1 = all_photons->affine_param_p[photon_idx];
    const double t2 = all_photons->affine_param[photon_idx];

    double t_linear;
    if ( (f1 * f2 <= 0.0 || fabs(f1) < 1e-12) && fabs(f2 - f1) > 1e-15 ) { 
        t_linear = (f2 * t1 - f1 * t2) / (f2 - f1);
    } else if ( (f0 * f1 <= 0.0 || fabs(f0) < 1e-12) && fabs(f1 - f0) > 1e-15 ) { 
        t_linear = (f1 * t0 - f0 * t1) / (f1 - f0);
    } else {
        t_linear = t1;
    }

    const double h0 = t1 - t0;
    const double h1 = t2 - t1;
    double lambda_event = t_linear;

    if (fabs(h0) > 1e-15 && fabs(h1) > 1e-15) {
        const double delta0 = (f1 - f0) / h0;
        const double delta1 = (f2 - f1) / h1;
        const double a = (delta1 - delta0) / (h1 + h0);
        const double b = a * h1 + delta1;
        const double discriminant = b*b - 4.0 * a * f2;

        if (discriminant >= 0.0 && fabs(a) > 1e-16) {
            double denom = (b >= 0.0) ? (b + sqrt(discriminant)) : (b - sqrt(discriminant));
            if (fabs(denom) > 1e-16) {
                double t_quad = t2 - (2.0 * f2 / denom);
                double t_min = (t0 < t2) ? t0 : t2;
                double t_max = (t0 < t2) ? t2 : t0;
                if (t_quad >= t_min && t_quad <= t_max) lambda_event = t_quad;
            }
        }
    }

    const double t = lambda_event;
    double L0, L1, L2;
    if (fabs(h0) < 1e-15 || fabs(h1) < 1e-15) {
        L0 = 0.0; L1 = (t2 - t) / (t2 - t1); L2 = (t - t1) / (t2 - t1);
    } else {
        L0 = ((t - t1) * (t - t2)) / ((t0 - t1) * (t0 - t2));
        L1 = ((t - t0) * (t - t2)) / ((t1 - t0) * (t1 - t2));
        L2 = ((t - t0) * (t - t1)) / ((t2 - t0) * (t2 - t1));
    }

    double *intersect_ptr = (event_type == SOURCE_EVENT) ? all_photons->source_event_f_intersect : all_photons->window_event_f_intersect;
    
    if (event_type == SOURCE_EVENT) {
        all_photons->source_event_lambda[photon_idx] = lambda_event;
        all_photons->source_event_found[photon_idx] = true;
    } else {
        all_photons->window_event_lambda[photon_idx] = lambda_event;
        all_photons->window_event_found[photon_idx] = true;
    }

    for (int i = 0; i < 9; i++) {
        intersect_ptr[IDX_GLOBAL(i, photon_idx, num_rays)] = GET_COMP(all_photons->f_p_p, i) * L0 + 
                                                            GET_COMP(all_photons->f_p,   i) * L1 + 
                                                            GET_COMP(all_photons->f,     i) * L2;
    }

    #undef GET_COMP
    #undef PLANE_EVAL
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
        name=name,
        params=params,
        body=body,
        postfunc=postfunc  
    )