"""
This module contains the high-precision event-finding C kernel.

This module resolves exact coordinate intersections when a plane
crossing is detected. It utilizes second-order quadratic interpolation
for root-finding and Lagrange interpolation for state reconstruction at
the intersection boundaries. Executing strictly within thread-local registers
bypasses global memory fetches, preserving the sm_86 architecture limit of 255 registers per thread.
Author: Dalton J. Moone."""

import nrpy.c_function as cfc


def find_event_time_and_state() -> None:
    """
    This function defines the find_event_time_and_state C function configuration.

    :raises SystemError: If C function registration fails within the NRPy+ pipeline.
    """
    includes = ["BHaH_defines.h", "<math.h>"]
    desc = r"""@brief Portable high-performance second-order root-finding.

    @param f_local The thread-local state array for step $f^\mu_{n}$.
    @param f_p_local The thread-local state array for step $f^\mu_{n-1}$.
    @param f_p_p_local The thread-local state array for step $f^\mu_{n-2}$.
    @param lam The current affine parameter $\lambda_{n}$.
    @param lam_p The history affine parameter $\lambda_{n-1}$.
    @param lam_p_p The history affine parameter $\lambda_{n-2}$.
    @param normal The geometric unit normal vector $n_i$ of the target plane.
    @param dist The scalar distance $d$ from the origin to the target plane.
    @param event_lambda Pointer to the local affine parameter $\lambda$ to be updated.
    @param event_f_intersect The thread-local array where the reconstructed state $f^\mu$ is stored.

    Detailed algorithm: Uses position data $x^i$ from the current and two previous
    integration steps to construct a quadratic model of the trajectory relative 
    to the target plane. The intersection $\lambda$ is solved via the quadratic formula 
    and the full state $f^\mu$ is reconstructed via Lagrange polynomials.
    Mapping logic directly to thread-local registers preserves the strict sm_86 architecture limits by
    bypassing global memory fetches and keeping all intermediates within the 255
    registers per thread."""
    cfunc_type = "static BHAH_HD_INLINE void"
    name = "find_event_time_and_state"
    params = (
        "const double *restrict f_local, "
        "const double *restrict f_p_local, "
        "const double *restrict f_p_p_local, "
        "const double lam, "
        "const double lam_p, "
        "const double lam_p_p, "
        "const double *restrict normal, "
        "const double dist, "
        "double *restrict event_lambda, "
        "double *restrict event_f_intersect"
    )
    include_CodeParameters_h = False
    body = r"""
    // --- TEMPORAL AND GEOMETRIC STATE UNPACKING ---
    // Evaluates the geometric plane equation $n_i x^i - d$ at the historical steps.
    // Executing this via thread-local registers minimizes instruction latency.

    const double f0 = f_p_p_local[1]*normal[0] + f_p_p_local[2]*normal[1] + f_p_p_local[3]*normal[2] - dist; // Plane evaluation at step $n-2$.
    const double f1 = f_p_local[1]*normal[0] + f_p_local[2]*normal[1] + f_p_local[3]*normal[2] - dist;     // Plane evaluation at step $n-1$.
    const double f2 = f_local[1]*normal[0] + f_local[2]*normal[1] + f_local[3]*normal[2] - dist;           // Plane evaluation at step $n$.

    // --- STEP 1: LINEAR SEED CALCULATION ---
    // Provides a fallback value if the quadratic interpolation fails or is unnecessary.
    // Uses linear interpolation $\frac{y_2 x_1 - y_1 x_2}{y_2 - y_1}$ for numerical stability.

    double t_linear; // Fallback linear interpolation for the affine parameter $\lambda$.
    if ( (f1 * f2 <= 0.0 || fabs(f1) < 1e-12) && fabs(f2 - f1) > 1e-15 ) {
        t_linear = (f2 * lam_p - f1 * lam) / (f2 - f1);
    } else if ( (f0 * f1 <= 0.0 || fabs(f0) < 1e-12) && fabs(f1 - f0) > 1e-15 ) {
        t_linear = (f1 * lam_p_p - f0 * lam_p) / (f1 - f0);
    } else {
        t_linear = lam_p;
    }

    // --- STEP 2: QUADRATIC INTERPOLATION (SECOND-ORDER ACCURACY) ---
    const double h0 = lam_p - lam_p_p; // Interval size $\Delta\lambda$ between step $n-2$ and $n-1$.
    const double h1 = lam - lam_p;     // Interval size $\Delta\lambda$ between step $n-1$ and $n$.
    double lambda_event = t_linear;    // The resulting affine parameter $\lambda$ for the crossing.

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
                double t_quad = lam - (2.0 * f2 / denom); // The computed quadratic root $\lambda_{root}$.
                double t_min = (lam_p_p < lam) ? lam_p_p : lam; // Minimum bound of the historical integration window.
                double t_max = (lam_p_p < lam) ? lam : lam_p_p; // Maximum bound of the historical integration window.
                // Bound check: ensure the quadratic root $\lambda_{root}$ lies within the historical integration window.
                if (t_quad >= t_min && t_quad <= t_max) lambda_event = t_quad;
            }
        }
    }

    // --- STEP 3: LAGRANGE STATE RECONSTRUCTION ---
    // Reconstruct the full 9-component state vector $f^\mu$ at the exact $\lambda_{event}$.
    const double t = lambda_event; // Local copy of the event parameter $\lambda_{event}$.
    double L0, L1, L2; // Lagrange basis polynomials $L_i(t)$ mapped to the thread-local state.
    if (fabs(h0) < 1e-15 || fabs(h1) < 1e-15) {
        // Fallback to linear weights $L_i$ if the intervals are degenerate.
        L0 = 0.0; L1 = (lam - t) / (lam - lam_p); L2 = (t - lam_p) / (lam - lam_p);
    } else {
        L0 = ((t - lam_p) * (t - lam)) / ((lam_p_p - lam_p) * (lam_p_p - lam));
        L1 = ((t - lam_p_p) * (t - lam)) / ((lam_p - lam_p_p) * (lam_p - lam));
        L2 = ((t - lam_p_p) * (t - lam_p)) / ((lam - lam_p_p) * (lam - lam_p));
    }

    // Assign the final computed event parameter $\lambda$.
    *event_lambda = lambda_event;

    // Interpolate all 9 components $f^\mu$: $\{t, x, y, z, p_t, p_x, p_y, p_z, \lambda\}$.
    for (int i = 0; i < 9; i++) { // Loop index $i$ over the 9 state vector components.
        event_f_intersect[i] = f_p_p_local[i] * L0 +
                               f_p_local[i]   * L1 +
                               f_local[i]     * L2;
    }
    """

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body
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