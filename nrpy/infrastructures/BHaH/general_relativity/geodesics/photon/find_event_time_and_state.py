# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/find_event_time_and_state.py
"""
Defines the event-finding C kernel.

This module resolves exact coordinate intersections when a boundary plane crossing is
detected by evaluating historical trajectory points against the target plane equation.
It accepts the thread-local state arrays and integration parameters from the current and
two previous integration steps, along with the geometric unit normal and orthogonal
distance of the target plane. The algorithm executes a fallback linear seed followed
by second-order quadratic interpolation to resolve the exact event parameter. The full
tensor state is then reconstructed at the boundary using three-point Lagrange
polynomials. Mapping this logic directly to local memory space minimizes global data
transfers and maintains all intermediate state reconstructions within local hardware
threads.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc


def find_event_time_and_state() -> None:
    """Define find_event_time_and_state C function configuration."""
    includes = ["BHaH_defines.h"]
    desc = r""" Portable high-performance second-order root-finding.

    @param f_local The thread-local state array for step $f^\mu_{n}$.
    @param f_p_local The thread-local state array for step $f^\mu_{n-1}$.
    @param f_p_p_local The thread-local state array for step $f^\mu_{n-2}$.
    @param integration_param The current integration parameter.
    @param integration_param_prev The preceding integration parameter.
    @param integration_param_pre_prev The integration parameter two steps earlier.
    @param normal The geometric unit normal vector $n_i$ of the target plane.
    @param dist The scalar distance $d$ from the origin to the target plane.
    @param event_integration_param Pointer to the interpolated event parameter.
    @param event_f_intersect The thread-local array where the reconstructed state $f^\mu$ is stored.

    Detailed algorithm: Uses position data $x^i$ from the current and two previous
    integration steps to construct a quadratic model of the trajectory relative
    to the target plane. The event parameter is solved via the quadratic formula
    and the full state $f^\mu$ is reconstructed via Lagrange polynomials.
    Mapping logic directly to thread-local registers minimizes global memory fetches
    and keeps all intermediates within hardware register limits."""
    cfunc_type = "BHAH_HD_INLINE void"
    name = "find_event_time_and_state"
    params = (
        "const double *restrict f_local, "
        "const double *restrict f_p_local, "
        "const double *restrict f_p_p_local, "
        "const double integration_param, "
        "const double integration_param_prev, "
        "const double integration_param_pre_prev, "
        "const double *restrict normal, "
        "const double dist, "
        "double *restrict event_integration_param, "
        "double *restrict event_f_intersect"
    )
    include_CodeParameters_h = False
    body = r"""
    //==========================================
    // TEMPORAL AND GEOMETRIC STATE UNPACKING
    //==========================================
    // Evaluates the geometric plane equation $n_i x^i - d$ at the historical steps.
    // Executing this via thread-local registers minimizes instruction latency.

    const double f0 = f_p_p_local[1]*normal[0] + f_p_p_local[2]*normal[1] + f_p_p_local[3]*normal[2] - dist; // Plane evaluation at step $n-2$.
    const double f1 = f_p_local[1]*normal[0] + f_p_local[2]*normal[1] + f_p_local[3]*normal[2] - dist;     // Plane evaluation at step $n-1$.
    const double f2 = f_local[1]*normal[0] + f_local[2]*normal[1] + f_local[3]*normal[2] - dist;           // Plane evaluation at step $n$.

    //==========================================
    // STEP 1: LINEAR SEED CALCULATION
    //==========================================
    // Provides a fallback value if the quadratic interpolation fails or is unnecessary.
    // Uses linear interpolation $\frac{y_2 x_1 - y_1 x_2}{y_2 - y_1}$ for numerical stability.

    double integration_param_linear; // Fallback linear interpolation for the event parameter.
    if ( (f1 * f2 <= 0.0 || fabs(f1) < 1e-12) && fabs(f2 - f1) > 1e-15 ) {
        integration_param_linear =
            (f2 * integration_param_prev - f1 * integration_param) / (f2 - f1);
    } else if ( (f0 * f1 <= 0.0 || fabs(f0) < 1e-12) && fabs(f1 - f0) > 1e-15 ) {
        integration_param_linear =
            (f1 * integration_param_pre_prev - f0 * integration_param_prev) /
            (f1 - f0);
    } else {
        integration_param_linear = integration_param_prev;
    } // END ELSE: fallback to previous parameter

    //==========================================
    // STEP 2: QUADRATIC INTERPOLATION (SECOND-ORDER ACCURACY)
    //==========================================
    const double h0 = integration_param_prev - integration_param_pre_prev;
    const double h1 = integration_param - integration_param_prev;
    double interpolated_event_integration_param = integration_param_linear;

    if (fabs(h0) > 1e-15 && fabs(h1) > 1e-15) {
        // Newton form divided differences for the quadratic model: $f(t) = a(t-t_2)^2 + b(t-t_2) + f_2$.
        const double delta0 = (f1 - f0) / h0; // The first divided difference for interval $n-2$ to $n-1$.
        const double delta1 = (f2 - f1) / h1; // The first divided difference for interval $n-1$ to $n$.
        const double a = (delta1 - delta0) / (h1 + h0); // The second-order coefficient representing relative acceleration.
        const double b = a * h1 + delta1;               // The first-order coefficient representing relative velocity.
        const double discriminant = b*b - 4.0 * a * f2; // Discriminant $b^2 - 4ac$ for the quadratic intersection.

        if (discriminant >= 0.0 && fabs(a) > 1e-16) {
            // Use the stable form of the quadratic formula $t = t_2 - \frac{2c}{-b \pm \sqrt{b^2 - 4ac}}$ to minimize truncation error.
            double denom = (b >= 0.0) ? (b + sqrt(discriminant)) : (b - sqrt(discriminant)); // Stable denominator for the quadratic root.
            if (fabs(denom) > 1e-16) {
                const double integration_param_quad =
                    integration_param - (2.0 * f2 / denom);
                const double integration_param_min =
                    (integration_param_pre_prev < integration_param)
                        ? integration_param_pre_prev
                        : integration_param;
                const double integration_param_max =
                    (integration_param_pre_prev < integration_param)
                        ? integration_param
                        : integration_param_pre_prev;
                if (integration_param_quad >= integration_param_min &&
                    integration_param_quad <= integration_param_max) {
                    interpolated_event_integration_param = integration_param_quad;
                } // END IF: quadratic root lies inside the interpolation interval
            } // END IF: stable denominator for quadratic root
        } // END IF: discriminant >= 0 and a > 1e-16
    } // END IF: intervals are not degenerate

    //==========================================
    // STEP 3: LAGRANGE STATE RECONSTRUCTION
    //==========================================
    // Reconstruct the full 9-component state vector at the event parameter.
    const double event_param = interpolated_event_integration_param;
    double L0, L1, L2; // Lagrange basis polynomials $L_i(t)$ mapped to the thread-local state.
    if (fabs(h0) < 1e-15 || fabs(h1) < 1e-15) {
        // Fallback to linear weights $L_i$ if the intervals are degenerate.
        L0 = 0.0;
        L1 = (integration_param - event_param) /
             (integration_param - integration_param_prev);
        L2 = (event_param - integration_param_prev) /
             (integration_param - integration_param_prev);
    } else {
        L0 = ((event_param - integration_param_prev) *
              (event_param - integration_param)) /
             ((integration_param_pre_prev - integration_param_prev) *
              (integration_param_pre_prev - integration_param));
        L1 = ((event_param - integration_param_pre_prev) *
              (event_param - integration_param)) /
             ((integration_param_prev - integration_param_pre_prev) *
              (integration_param_prev - integration_param));
        L2 = ((event_param - integration_param_pre_prev) *
              (event_param - integration_param_prev)) /
             ((integration_param - integration_param_pre_prev) *
              (integration_param - integration_param_prev));
    } // END ELSE: compute full quadratic Lagrange weights

    // Assign the interpolated event parameter.
    *event_integration_param = interpolated_event_integration_param;

    // Interpolate all nine state components.
    for (int i = 0; i < 9; i++) { // Loop index $i$ over the 9 state vector components.
        event_f_intersect[i] = f_p_p_local[i] * L0 +
                               f_p_local[i]   * L1 +
                               f_local[i]     * L2;
    } // END LOOP: for i over 9 state vector components
    """

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
