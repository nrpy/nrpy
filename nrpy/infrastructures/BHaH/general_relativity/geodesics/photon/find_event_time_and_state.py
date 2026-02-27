"""
Registers the high-precision event-finding C kernel.

Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
Utilizes second-order quadratic interpolation for root-finding and
Lagrange interpolation for state reconstruction at plane intersections.
"""

import nrpy.c_function as cfc


def find_event_time_and_state() -> None:
    """
    Register the find_event_time_and_state C function.

    This function resolves the exact coordinate intersection when a plane
    crossing is detected. It ensures sub-step accuracy for photon mapping.

    >>> find_event_time_and_state()
    """
    # 1. Define C-Function metadata
    includes = ["BHaH_defines.h", "math.h"]
    desc = """@brief Portable high-performance second-order root-finding.
    Detailed algorithm: Uses position data from the current and two previous
    integration steps to construct a quadratic model of the trajectory relative
    to the target plane. The intersection is solved via the quadratic formula
    and the full state is reconstructed via Lagrange polynomials."""
    name = "find_event_time_and_state"

    params = """PhotonStateSoA *restrict all_photons,
                const long int num_rays,
                const long int photon_idx,
                const double *restrict normal,
                const double dist,
                const event_type_t event_type"""

    # 2. Build the C body with internal descriptive comments
    body = r"""
    // --- Helper Macros for SoA Access ---
    #define GET_COMP(soa_ptr, comp) (soa_ptr[IDX_GLOBAL(comp, photon_idx, num_rays)])
    #define PLANE_EVAL(soa_ptr) (GET_COMP(soa_ptr, 1)*normal[0] + GET_COMP(soa_ptr, 2)*normal[1] + GET_COMP(soa_ptr, 3)*normal[2] - dist)

    // --- Preamble: Unpack Temporal and Geometric State History ---
    // f2: current state, f1: previous state, f0: two steps prior.
    const double f0 = PLANE_EVAL(all_photons->f_p_p); // Plane evaluation at step n-2.
    const double f1 = PLANE_EVAL(all_photons->f_p);   // Plane evaluation at step n-1.
    const double f2 = PLANE_EVAL(all_photons->f);     // Plane evaluation at step n.

    const double t0 = all_photons->affine_param_p_p[photon_idx]; // Affine parameter lambda at step n-2.
    const double t1 = all_photons->affine_param_p[photon_idx];   // Affine parameter lambda at step n-1.
    const double t2 = all_photons->affine_param[photon_idx];     // Affine parameter lambda at step n.

    // --- Step 1: Linear Seed Calculation ---
    // Provides a fallback value if the quadratic interpolation fails or is unnecessary.
    double t_linear;
    if ( (f1 * f2 <= 0.0 || fabs(f1) < 1e-12) && fabs(f2 - f1) > 1e-15 ) {
        t_linear = (f2 * t1 - f1 * t2) / (f2 - f1);
    } else if ( (f0 * f1 <= 0.0 || fabs(f0) < 1e-12) && fabs(f1 - f0) > 1e-15 ) {
        t_linear = (f1 * t0 - f0 * t1) / (f1 - f0);
    } else {
        t_linear = t1;
    }

    // --- Step 2: Quadratic Interpolation (Second-Order Accuracy) ---
    const double h0 = t1 - t0; // Interval size between step n-2 and n-1.
    const double h1 = t2 - t1; // Interval size between step n-1 and n.
    double lambda_event = t_linear; // The resulting affine parameter for the crossing.

    if (fabs(h0) > 1e-15 && fabs(h1) > 1e-15) {
        // Newton form divided differences for the quadratic model: f(t) = a*(t-t2)^2 + b*(t-t2) + f2
        const double delta0 = (f1 - f0) / h0;
        const double delta1 = (f2 - f1) / h1;
        const double a = (delta1 - delta0) / (h1 + h0); // The second-order coefficient (acceleration relative to plane).
        const double b = a * h1 + delta1;               // The first-order coefficient (velocity relative to plane).
        const double discriminant = b*b - 4.0 * a * f2; // Discriminant for the quadratic intersection.

        if (discriminant >= 0.0 && fabs(a) > 1e-16) {
            // Use the stable form of the quadratic formula to find the root closest to the current step.
            double denom = (b >= 0.0) ? (b + sqrt(discriminant)) : (b - sqrt(discriminant));
            if (fabs(denom) > 1e-16) {
                double t_quad = t2 - (2.0 * f2 / denom);
                double t_min = (t0 < t2) ? t0 : t2;
                double t_max = (t0 < t2) ? t2 : t0;
                // Bound check: ensure the quadratic root lies within the historical integration window.
                if (t_quad >= t_min && t_quad <= t_max) lambda_event = t_quad;
            }
        }
    }

    // --- Step 3: Lagrange State Reconstruction ---
    // Reconstruct the full 9-component state vector at the exact lambda_event.
    const double t = lambda_event;
    double L0, L1, L2; // Lagrange basis polynomials.
    if (fabs(h0) < 1e-15 || fabs(h1) < 1e-15) {
        // Fallback to linear weights if the steps are degenerate.
        L0 = 0.0; L1 = (t2 - t) / (t2 - t1); L2 = (t - t1) / (t2 - t1);
    } else {
        L0 = ((t - t1) * (t - t2)) / ((t0 - t1) * (t0 - t2));
        L1 = ((t - t0) * (t - t2)) / ((t1 - t0) * (t1 - t2));
        L2 = ((t - t0) * (t - t1)) / ((t2 - t0) * (t2 - t1));
    }

    // Assign results to the appropriate intersection buffer (Source vs Window).
    double *intersect_ptr = (event_type == SOURCE_EVENT) ? all_photons->source_event_f_intersect : all_photons->window_event_f_intersect;

    if (event_type == SOURCE_EVENT) {
        all_photons->source_event_lambda[photon_idx] = lambda_event;
        all_photons->source_event_found[photon_idx] = true;
    } else {
        all_photons->window_event_lambda[photon_idx] = lambda_event;
        all_photons->window_event_found[photon_idx] = true;
    }

    // Interpolate all 9 components: {t, x, y, z, p_t, p_x, p_y, p_z, lambda}
    for (int i = 0; i < 9; i++) {
        intersect_ptr[IDX_GLOBAL(i, photon_idx, num_rays)] = GET_COMP(all_photons->f_p_p, i) * L0 +
                                                            GET_COMP(all_photons->f_p,   i) * L1 +
                                                            GET_COMP(all_photons->f,     i) * L2;
    }

    #undef GET_COMP
    #undef PLANE_EVAL
    """

    prefunc = "#ifdef USE_GPU\n#pragma omp declare target\n#endif"
    postfunc = "#ifdef USE_GPU\n#pragma omp end declare target\n#endif"

    # 3. Register the C function
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
        postfunc=postfunc,
    )
