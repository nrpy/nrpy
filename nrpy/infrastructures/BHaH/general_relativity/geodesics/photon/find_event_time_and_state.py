# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/find_event_time_and_state.py
"""
Generate centered interpolation for photon plane crossings.

The generated C function reconstructs all nine photon state components from
accepted states. The integration parameter is affine parameter for geodesic
EOM and coordinate time for normalized EOM.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc


def find_event_time_and_state(maximum_degree: int) -> None:
    """
    Register centered interpolation over an accepted plane-crossing step.

    :param maximum_degree: Largest generated polynomial degree.
    :raises ValueError: If maximum_degree is below three.
    """
    if maximum_degree < 3:
        raise ValueError("Plane interpolation degree must be at least three")
    evaluator = r"""
BHAH_HD_INLINE double evaluate_plane_interpolation_scalar(
    const double *restrict nodes,
    const double *restrict weights,
    const double *restrict values,
    const int point_count,
    const double argument) {
    for (int node = 0; node < point_count; ++node) {
        if (argument == nodes[node]) return values[node];
    } // END LOOP: return exact interpolation node
    double numerator = 0.0;
    double denominator = 0.0;
    for (int node = 0; node < point_count; ++node) {
        const double term = weights[node] / (argument - nodes[node]);
        numerator += term * values[node];
        denominator += term;
    } // END LOOP: accumulate barycentric terms
    return numerator / denominator;
} // END FUNCTION: evaluate_plane_interpolation_scalar
"""
    # Try the requested degree first, then each lower even degree. Degree two
    # remains available when a physical stop leaves only one post-crossing state.
    degrees = list(range(2, maximum_degree, 2))
    degrees.append(maximum_degree)
    body = r"""
    // History index zero is the newest accepted state. The crossing step is
    // between indices steps_past_plane - 1 and steps_past_plane.
    *actual_degree = 0;
    if (steps_past_plane < 1 || steps_past_plane > MAXIMUM_DEGREE ||
        !isfinite(dist) || !isfinite(normal[0]) ||
        !isfinite(normal[1]) || !isfinite(normal[2])) return -1;
""".replace("MAXIMUM_DEGREE", str(maximum_degree))
    for degree in reversed(degrees):
        points_after = (degree + 1) // 2
        # Even degrees have two equally centered contiguous stencils. Try the
        # older-side stencil first, then the newer-side stencil if history repeats.
        maximum_shift = 1 if degree % 2 == 0 else 0
        body += f"""
        // Degree {degree} uses {degree + 1} contiguous accepted states.
        for (int stencil_shift = 0; stencil_shift <= {maximum_shift}; ++stencil_shift) {{
            if (steps_past_plane < {points_after} + stencil_shift) continue;
            const int candidate_first =
                steps_past_plane - {points_after} - stencil_shift;
            bool usable = candidate_first >= 0 &&
                          candidate_first + {degree} <= {maximum_degree};
            if (usable) {{
                const double first_param = integration_params[candidate_first];
                const double second_param = integration_params[candidate_first + 1];
                const bool increasing = second_param > first_param;
                const bool decreasing = second_param < first_param;
                usable = isfinite(first_param) && isfinite(second_param) &&
                         (increasing || decreasing);
                for (int node = 1; usable && node <= {degree}; ++node) {{
                    const double previous = integration_params[candidate_first + node - 1];
                    const double current = integration_params[candidate_first + node];
                    usable = isfinite(current) &&
                             ((increasing && current > previous) ||
                              (decreasing && current < previous));
                }} // END LOOP: check integration parameter ordering
                for (int node = 0; usable && node <= {degree}; ++node) {{
                    for (int component = 0; component < 9; ++component) {{
                        if (!isfinite(state_history[9 * (candidate_first + node) + component])) {{
                            usable = false;
                            break;
                        }} // END IF: nonfinite history component
                    }} // END LOOP: inspect photon components
                }} // END LOOP: inspect accepted states
            }} // END IF: candidate history usable
            if (usable) {{
                const int point_count = {degree + 1};
                const int first_state_index = candidate_first;
                const double origin = integration_params[steps_past_plane];
                double scale = 0.0;
                for (int node = 0; node < point_count; ++node) {{
                    const double distance = fabs(
                        integration_params[first_state_index + node] - origin);
                    if (distance > scale) scale = distance;
                }} // END LOOP: find parameter scale
                if (!(scale > 0.0) || !isfinite(scale)) return -1;

                double nodes[{degree + 1}];
                double weights[{degree + 1}];
                double plane_values[{degree + 1}];
                for (int node = 0; node < point_count; ++node) {{
                    const int history_index = first_state_index + node;
                    nodes[node] = (integration_params[history_index] - origin) / scale;
                    const double *restrict state = &state_history[9 * history_index];
                    plane_values[node] = state[1] * normal[0] +
                                         state[2] * normal[1] +
                                         state[3] * normal[2] - dist;
                    if (!isfinite(plane_values[node])) return -1;
                }} // END LOOP: evaluate plane distances

                double largest_weight = 0.0;
                for (int node = 0; node < point_count; ++node) {{
                    double weight = 1.0;
                    for (int other = 0; other < point_count; ++other) {{
                        if (node != other) weight /= nodes[node] - nodes[other];
                    }} // END LOOP: calculate barycentric weight
                    if (!isfinite(weight)) return -1;
                    weights[node] = weight;
                    if (fabs(weight) > largest_weight) largest_weight = fabs(weight);
                }} // END LOOP: calculate interpolation weights
                if (!(largest_weight > 0.0)) return -1;
                for (int node = 0; node < point_count; ++node) {{
                    weights[node] /= largest_weight;
                }} // END LOOP: normalize barycentric weights

                const int after_index = steps_past_plane - 1 - first_state_index;
                const int before_index = after_index + 1;
                double left = nodes[after_index];
                double right = nodes[before_index];
                double left_value = plane_values[after_index];
                double right_value = plane_values[before_index];
                if (!((left_value <= 0.0 && right_value >= 0.0) ||
                      (left_value >= 0.0 && right_value <= 0.0))) return -1;

                for (int iteration = 0; iteration < 80; ++iteration) {{
                    const double middle = left + 0.5 * (right - left);
                    if (middle == left || middle == right) break;
                    const double middle_value = evaluate_plane_interpolation_scalar(
                        nodes, weights, plane_values, point_count, middle);
                    if (!isfinite(middle_value)) return -1;
                    if (middle_value == 0.0) {{
                        left = middle;
                        right = middle;
                        break;
                    }} // END IF: exact plane root
                    if ((left_value < 0.0 && middle_value < 0.0) ||
                        (left_value > 0.0 && middle_value > 0.0)) {{
                        left = middle;
                        left_value = middle_value;
                    }} // END IF: root remains left
                    else {{
                        right = middle;
                        right_value = middle_value;
                    }} // END ELSE: root remains right
                }} // END LOOP: refine bracketed plane root
                const double root = fabs(left_value) < fabs(right_value)
                    ? left : right;
                *event_integration_param = origin + scale * root;
                if (!isfinite(*event_integration_param)) return -1;

                double component_values[{degree + 1}];
                for (int component = 0; component < 9; ++component) {{
                    for (int node = 0; node < point_count; ++node) {{
                        component_values[node] = state_history[
                            9 * (first_state_index + node) + component];
                    }} // END LOOP: read photon component history
                    event_f_intersect[component] = evaluate_plane_interpolation_scalar(
                        nodes, weights, component_values, point_count, root);
                    if (!isfinite(event_f_intersect[component])) return -1;
                }} // END LOOP: reconstruct photon components
                *actual_degree = {degree};
                return 1;
            }} // END IF: interpolate usable history
        }} // END LOOP: select centered stencil
    """
    body += "\n    return 0; // No generated degree has distinct accepted states.\n"
    cfc.register_CFunction(
        prefunc=evaluator,
        includes=["BHaH_defines.h", "math.h"],
        desc=r"""Interpolate the full photon state at a bracketed plane crossing.

        The history is newest first. For geodesic EOM, the integration
        parameter is affine parameter and state component zero is coordinate
        time. For normalized EOM, the integration parameter is coordinate
        time and state component zero is affine parameter. Return 1 on
        success, 0 if no generated degree has usable accepted states, or -1
        if the plane root cannot be reconstructed.

        @param[in] state_history Nine-component accepted states, newest first.
        @param[in] integration_params Accepted integration parameters, newest first.
        @param steps_past_plane Number of accepted states after crossing.
        @param[in] normal Unit normal of the crossing plane.
        @param dist Plane offset from the Cartesian origin.
        @param[out] event_integration_param Integration parameter at the crossing.
        @param[out] event_f_intersect Interpolated nine-component photon state.
        @param[out] actual_degree Polynomial degree selected from generated cases.
        """,
        cfunc_type="BHAH_HD_INLINE int",
        name="find_event_time_and_state_centered",
        params=(
            "const double *restrict state_history, "
            "const double *restrict integration_params, "
            "const int steps_past_plane, "
            "const double *restrict normal, "
            "const double dist, "
            "double *restrict event_integration_param, "
            "double *restrict event_f_intersect, "
            "int *restrict actual_degree"
        ),
        include_CodeParameters_h=False,
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
