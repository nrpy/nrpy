# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/rkf45_numerical_control.py
r"""
Defines the CPU/OpenMP helper that applies numerical RKF45 step control.

This module registers a host-only C helper for the normalized-photon
numerical-spacetime RKF45 pipeline. It combines the existing embedded
Cash-Karp truncation estimate with the normalized-photon constraint
``\gamma^{ij} \Pi_i \Pi_j = 1`` evaluated at the interpolated 5th-order
candidate state. The helper then uses the larger of the two normalized error
measures to accept or reject each ray independently, update the next step size,
and preserve the numerical time-window contract used by slot-based metric
mapping.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

from typing import List, Union

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.params as par


def rkf45_numerical_control(
    norm_expr: sp.Expr, enable_numerical_time_window_step_cap: bool = True
) -> None:
    r"""
    Register the CPU/OpenMP helper that applies RKF45 numerical step control.

    This helper reads the anchored start state, the 5th-order candidate state,
    the metric interpolated at that candidate state, and the stored RKF45 stage
    derivatives. It computes the standard embedded Cash-Karp truncation error,
    evaluates the normalized-photon constraint residual at the candidate state,
    and uses the larger normalized error to accept or reject each ray.

    :param norm_expr: SymPy expression for the normalized-photon constraint
        ``\gamma^{ij} \Pi_i \Pi_j``.
    :param enable_numerical_time_window_step_cap: Whether to preserve the
        slot-based numerical time-window step-size cap driven by
        ``rkf45_max_delta_t``.
    :raises ValueError: If ``norm_expr`` is not a SymPy expression.
    :raises ValueError: If the configured parallelization mode is not OpenMP.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> import nrpy.equations.general_relativity.geodesics.geodesics as geo
    >>> import os
    >>> import tempfile
    >>> cfc.CFunction_dict.clear()
    >>> generic_geodesic_equations = geo.GeodesicEquations.__new__(geo.GeodesicEquations)
    >>> norm_expr = generic_geodesic_equations.normalization_constraint_photon_normalized()
    >>> with tempfile.TemporaryDirectory(dir=os.getcwd()) as temp_dir:
    ...     old_cache_home = os.environ.get("XDG_CACHE_HOME")
    ...     _ = os.environ.__setitem__("XDG_CACHE_HOME", temp_dir)
    ...     rkf45_numerical_control(norm_expr, enable_numerical_time_window_step_cap=True)
    ...     generated = cfc.CFunction_dict["rkf45_numerical_control"].full_function
    ...     if old_cache_home is None:
    ...         _ = os.environ.pop("XDG_CACHE_HOME", None)
    ...     else:
    ...         _ = os.environ.__setitem__("XDG_CACHE_HOME", old_cache_home)
    >>> "#pragma omp parallel for" in generated
    True
    >>> "constraint_err = fabs(C_5th - 1.0) / constraint_tol;" in generated
    True
    >>> "slot_position" in generated
    True
    """
    if not isinstance(norm_expr, sp.Expr):
        raise ValueError("norm_expr must be a valid SymPy expression.")

    parallelization = par.parval_from_str("parallelization")
    if parallelization != "openmp":
        raise ValueError(
            "rkf45_numerical_control currently supports only "
            "parallelization='openmp'."
        )

    if (
        enable_numerical_time_window_step_cap
        and "time_window_manager_numerical"
        not in par.glb_extras_dict.get("BHaH_defines", {})
    ):
        # Import lazily to avoid a package-level circular import between the
        # photon and interpolation aggregation modules.
        from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.time_window_manager_numerical import (  # pylint: disable=import-outside-toplevel
            time_window_manager_numerical,
        )

        time_window_manager_numerical()

    # Step 1: Register runtime parameters consumed by the numerical control helper.
    real_param_names: List[str] = [
        "rkf45_error_tolerance",
        "rkf45_absolute_error_tolerance",
        "rkf45_log_energy_tolerance",
        "rkf45_h_min",
        "rkf45_h_max",
        "rkf45_constraint_tolerance",
    ]
    real_param_defaults: List[Union[str, int, float]] = [
        1e-8,
        1e-8,
        1e-8,
        1e-10,
        10.0,
        1e-8,
    ]
    par.register_CodeParameters(
        "REAL",
        __name__,
        real_param_names,
        real_param_defaults,
        commondata=True,
        add_to_parfile=True,
    )
    par.register_CodeParameter(
        "int", __name__, "rkf45_max_retries", 10, commondata=True, add_to_parfile=True
    )

    # Step 2: Lower the symbolic normalized-photon constraint to C code.
    used_symbol_names = {str(sym) for sym in norm_expr.free_symbols}
    constraint_math = ccg.c_codegen(
        [norm_expr],
        ["C_5th"],
        enable_cse=True,
        verbose=False,
        include_braces=False,
    )

    hydration_lines = [
        "    //==========================================",
        "    // CANDIDATE-STATE HYDRATION FOR CONSTRAINT EVALUATION",
        "    //==========================================",
    ]

    if "u" in used_symbol_names:
        hydration_lines.append(
            "    const double u = d_f_5th_bundle[IDX_F(4, i)]; // Normalized photon energy variable $u$ at the candidate state."
        )
        hydration_lines.append(
            "    (void)u; // Suppress unused variable warnings when symbolic simplification removes $u$."
        )

    for idx in range(3):
        symbol_name = f"PiD{idx}"
        if symbol_name in used_symbol_names:
            hydration_lines.append(
                f"    const double {symbol_name} = d_f_5th_bundle[IDX_F({idx + 5}, i)]; // Candidate normalized covariant momentum component $\\Pi_{{{idx}}}$."
            )
            hydration_lines.append(
                f"    (void){symbol_name}; // Suppress unused variable warnings for candidate $\\Pi_{{{idx}}}$."
            )

    metric_idx = 0
    for i in range(4):
        for j in range(i, 4):
            comp_name = f"metric_g4DD{i}{j}"
            if comp_name in used_symbol_names:
                hydration_lines.append(
                    f"    const double {comp_name} = d_metric_5th_bundle[IDX_METRIC({metric_idx}, i)]; // Metric tensor component $g_{{{i}{j}}}$ at the candidate state."
                )
                hydration_lines.append(
                    f"    (void){comp_name}; // Suppress unused variable warnings for $g_{{{i}{j}}}$."
                )
            metric_idx += 1
    hydration_str = "\n".join(hydration_lines)

    # Step 3: Emit helper code above the registered C function.
    prefunc = r"""
static inline int rkf45_checked_floor_to_long(
    const double value,
    long int *out) {
  if (out == NULL || !isfinite(value) ||
      value < (double)LONG_MIN || value > (double)LONG_MAX) {
    return 1;
  } // END IF: floor() source value was not representable as long int
  *out = (long int)floor(value);
  return 0;
} // END FUNCTION: rkf45_checked_floor_to_long
"""
    postfunc = ""

    direct_step_cap = ""
    if enable_numerical_time_window_step_cap:
        direct_step_cap = r"""
    if (commondata->rkf45_max_delta_t > 0.0 &&
        h_new_abs > commondata->rkf45_max_delta_t) {
      h_new_abs = commondata->rkf45_max_delta_t;
      h_new = h_sign * h_new_abs;
    } // END IF: next-step magnitude exceeded rkf45_max_delta_t
"""

    accepted_step_cap = ""
    if enable_numerical_time_window_step_cap:
        accepted_step_cap = r"""
      //==========================================
      // NUMERICAL TIME-WINDOW ACCEPTED-STEP CAP
      //==========================================
      // Preserve the slot-based mapped-window contract used by the numerical
      // interpolation pipeline after the direct rkf45_max_delta_t cap has
      // already been enforced for both accepted and rejected next steps.
      {
        const double t_accepted = d_f_5th_bundle[IDX_F(0, i)];
        const double slot_position =
            (t_accepted - commondata->slot_manager_t_min) /
            commondata->slot_manager_delta_t;
        long int slot_idx = 0L;

        if (rkf45_checked_floor_to_long(slot_position, &slot_idx) != 0) {
          h_new_abs = commondata->rkf45_h_min;
          h_new = h_sign * h_new_abs;
          d_status[i] = TERMINATION_TYPE_FAILURE;
        } else {
          const double slot_lower =
              commondata->slot_manager_t_min +
              (double)slot_idx * commondata->slot_manager_delta_t;
          const double delta_to_lower_slot_edge =
              fmax(0.0, t_accepted - slot_lower);
          const double allowed_delta_t =
              delta_to_lower_slot_edge + commondata->rkf45_max_delta_t;

          if (allowed_delta_t > 0.0 && h_new_abs > allowed_delta_t) {
            h_new_abs = allowed_delta_t;
          } // END IF: next-step magnitude exceeded the mapped window lookahead
          h_new = h_sign * h_new_abs;
        } // END ELSE: slot position was representable as long int
      } // END BLOCK: numerical time-window accepted-step cap
"""

    # Step 4: Define the registered C function metadata.
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<limits.h>",
        "<math.h>",
        "<stdbool.h>",
    ]
    desc = r"""Apply RKF45 numerical error control for one active photon chunk.

    This helper combines the embedded Cash-Karp truncation estimate with the
    normalized-photon constraint residual evaluated using the metric
    interpolated at the 5th-order candidate state. Each ray is then accepted
    or rejected independently, preserving the persistent state on rejection.

    @param[in] commondata Global runtime parameters.
    @param[in,out] d_f_persistent Persistent photon state bundle updated only on acceptance.
    @param[in] d_f_start Anchored start-state bundle $f_n$.
    @param[in] d_f_5th_bundle Scratch bundle containing the 5th-order candidate state.
    @param[in] d_metric_5th_bundle Metric bundle interpolated at the 5th-order candidate state.
    @param[in] d_k_bundle Flattened RKF45 stage-derivative bundle.
    @param[in,out] d_h Per-ray step sizes.
    @param[in,out] d_status Per-ray termination status values.
    @param[in,out] d_affine Per-ray coordinate-time tracker used by the numerical interpolation pipeline.
    @param[in,out] d_retries Per-ray adaptive-step rejection counters.
    @param[in] chunk_size Number of active rays in the current chunk.
    @param[in] stream_idx Unused placeholder kept for interface compatibility with other photon helpers.

    @note When enabled, the accepted-step cap preserves the mapped numerical
    time-window contract driven by rkf45_max_delta_t and slot-manager metadata.
    """
    name = "rkf45_numerical_control"
    params = (
        "const commondata_struct *restrict commondata, "
        "double *restrict d_f_persistent, "
        "const double *restrict d_f_start, "
        "const double *restrict d_f_5th_bundle, "
        "const double *restrict d_metric_5th_bundle, "
        "const double *restrict d_k_bundle, "
        "double *restrict d_h, "
        "termination_type_t *restrict d_status, "
        "double *restrict d_affine, "
        "int *restrict d_retries, "
        "const long int chunk_size, "
        "const int stream_idx"
    )

    # Step 5: Build the CPU/OpenMP C body.
    body = f"""
  (void)stream_idx;

  //==========================================
  // MACRO DEFINITIONS FOR BUNDLE ACCESS
  //==========================================
  // IDX_F maps a component to the flattened state bundle using SoA layout.
  #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
  // IDX_K maps a stage/component pair into the flattened RKF45 derivative bundle.
  #define IDX_K(s, c, ray_id) ((s) * 9 * BUNDLE_CAPACITY + (c) * BUNDLE_CAPACITY + (ray_id))
  // IDX_METRIC maps a component to the flattened symmetric metric bundle.
  #define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))

  //==========================================
  // OPENMP LOOP ARCHITECTURE
  //==========================================
  // Distribute photon rays across available CPU threads for parallel evaluation.
  #pragma omp parallel for
  for (long int i = 0; i < chunk_size; i++) {{
    // Step 1: Load tolerances, step-size state, and retry state.
    const double rtol = commondata->rkf45_error_tolerance;
    const double atol = commondata->rkf45_absolute_error_tolerance;
    const double log_energy_tol = commondata->rkf45_log_energy_tolerance;
    const double constraint_tol =
        fmax(fabs(commondata->rkf45_constraint_tolerance), 1.0e-300);
    const double h_local = d_h[i];
    const double h_sign = (h_local < 0.0) ? -1.0 : 1.0;
    const int retries = d_retries[i];

    // Step 2: Compute the embedded Cash-Karp error exactly as in the
    // original RKF45 finalize logic, but using the precomputed 5th-order
    // candidate bundle for state scaling.
    double err_embedded = 0.0;
    for (int comp = 0; comp < 9; ++comp) {{
      const double f_n = d_f_start[IDX_F(comp, i)];
      const double f_5th_val = d_f_5th_bundle[IDX_F(comp, i)];
      const double k0 = d_k_bundle[IDX_K(0, comp, i)];
      // Bundle slot 1 carries zero weight in the embedded Cash-Karp estimator.
      const double k2 = d_k_bundle[IDX_K(2, comp, i)];
      const double k3 = d_k_bundle[IDX_K(3, comp, i)];
      const double k4 = d_k_bundle[IDX_K(4, comp, i)];
      const double k5 = d_k_bundle[IDX_K(5, comp, i)];

      if (!isfinite(f_5th_val)) {{
        err_embedded = 1.0e30;
      }} // END IF: candidate-state component was not finite

      double err_val = (1.0 / 360.0) * k0;
      err_val += (-128.0 / 4275.0) * k2;
      err_val += (-2197.0 / 75240.0) * k3;
      err_val += (1.0 / 50.0) * k4;
      err_val += (2.0 / 55.0) * k5;
      const double err_abs = fabs(h_local * err_val);

      if (comp > 0 && comp < 8) {{
        double scale = 0.0;
        const double state_mag = fmax(fabs(f_n), fabs(f_5th_val));

        if (comp <= 3) {{
          scale = atol + rtol * state_mag;
        }} else if (comp == 4) {{
          scale = log_energy_tol;
        }} else {{
          scale = atol + rtol * fmax(1.0, state_mag);
        }} // END ELSE: candidate normalized-momentum component scaling

        const double current_err = err_abs / scale;
        if (!isfinite(current_err)) {{
          err_embedded = 1.0e30;
        }} else {{
          err_embedded = fmax(err_embedded, current_err);
        }} // END ELSE: embedded-error accumulation
      }} // END IF: exclude slots 0 and 8 from the embedded error norm
    }} // END LOOP: for comp over 9 state-vector components

    // Step 3: Check the interpolated metric and evaluate the normalized-photon
    // constraint directly at the 5th-order candidate state.
    bool metric_data_finite = true;
    for (int metric_comp = 0; metric_comp < 10; ++metric_comp) {{
      const double metric_val = d_metric_5th_bundle[IDX_METRIC(metric_comp, i)];
      if (!isfinite(metric_val)) {{
        metric_data_finite = false;
        break;
      }} // END IF: one candidate-state metric component was not finite
    }} // END LOOP: for metric_comp over 10 metric components

{hydration_str}

    double C_5th = NAN;
    double constraint_err = 1.0e30;
    if (metric_data_finite) {{
      {constraint_math}
      if (isfinite(C_5th)) {{
        constraint_err = fabs(C_5th - 1.0) / constraint_tol;
      }} // END IF: constraint evaluation produced a finite scalar
    }} // END IF: all candidate-state metric components were finite

    // Step 4: Combine the embedded RK error and normalized-photon constraint
    // residual into one acceptance norm.
    double err_total = fmax(err_embedded, constraint_err);
    if (!isfinite(err_total)) {{
      err_total = 1.0e30;
    }} // END IF: combined acceptance norm became non-finite

    // Step 5: Update the next-step size using the standard RKF45 safety rule.
    const double safety = 0.9;
    const double factor = (err_total > 1.0e-15) ? pow(1.0 / err_total, 0.2) : 2.0;
    double h_new_abs = safety * fabs(h_local) * factor;
    h_new_abs = fmax(h_new_abs, commondata->rkf45_h_min);
    h_new_abs = fmin(h_new_abs, commondata->rkf45_h_max);
    double h_new = h_sign * h_new_abs;

{direct_step_cap}

    if (err_total <= 1.0) {{
      // Step 6: Accept the candidate state, synchronize the numerical
      // interpolation history tracker, and prepare the next step.
      for (int comp = 0; comp < 9; ++comp) {{
        d_f_persistent[IDX_F(comp, i)] = d_f_5th_bundle[IDX_F(comp, i)];
      }} // END LOOP: for comp over 9 accepted state-vector components

      d_affine[i] = d_f_5th_bundle[IDX_F(0, i)];
      d_retries[i] = 0;
      d_status[i] = ACTIVE;

{accepted_step_cap}
      d_h[i] = h_new;
    }} else {{
      // Step 7: Reject the attempted step for this ray only, preserving the
      // persistent state while shrinking the retry step size.
      const int new_retries = retries + 1;
      d_retries[i] = new_retries;

      if (new_retries > commondata->rkf45_max_retries) {{
        d_status[i] = FAILURE_RKF45_REJECTION_LIMIT;
      }} else {{
        d_status[i] = REJECTED;
      }} // END ELSE: rejection remained recoverable

      d_h[i] = h_new;
    }} // END ELSE: rejected candidate state
  }} // END LOOP: for i over chunk_size rays

  //==========================================
  // MACRO CLEANUP
  //==========================================
  #undef IDX_F
  #undef IDX_K
  #undef IDX_METRIC
"""

    # Step 6: Register the host-only numerical control helper.
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=False,
        postfunc=postfunc,
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
