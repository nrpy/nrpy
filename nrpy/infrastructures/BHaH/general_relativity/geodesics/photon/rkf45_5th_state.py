# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/rkf45_5th_state.py
"""
Defines the CPU/OpenMP helper that builds the RKF45 5th-order state.

This module registers a host-only C helper for the normalized-photon RKF45
pipeline. After all six RKF45 stage derivatives have been accumulated into the
flattened stage bundle, this helper evaluates the Cash-Karp 5th-order solution
component-by-component and writes the resulting candidate state into a scratch
bundle. The helper is intentionally limited to candidate construction so later
adaptive-control logic can combine the embedded truncation estimate with
additional physics-aware acceptance criteria.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.params as par


def rkf45_5th_state() -> None:
    r"""
    Register the CPU/OpenMP helper that builds the RKF45 5th-order state.

    This helper reads the anchored start state ``f_start``, the six RKF45 stage
    derivative vectors, and the per-ray step sizes ``h``. It then evaluates the
    Cash-Karp 5th-order update and writes the candidate state into a dedicated
    scratch bundle without modifying the persistent trajectory state.

    :raises ValueError: If the configured parallelization mode is not OpenMP.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> import os
    >>> import tempfile
    >>> cfc.CFunction_dict.clear()
    >>> with tempfile.TemporaryDirectory(dir=os.getcwd()) as temp_dir:
    ...     old_cache_home = os.environ.get("XDG_CACHE_HOME")
    ...     _ = os.environ.__setitem__("XDG_CACHE_HOME", temp_dir)
    ...     rkf45_5th_state()
    ...     generated = cfc.CFunction_dict["rkf45_5th_state"].full_function
    ...     if old_cache_home is None:
    ...         _ = os.environ.pop("XDG_CACHE_HOME", None)
    ...     else:
    ...         _ = os.environ.__setitem__("XDG_CACHE_HOME", old_cache_home)
    >>> "#pragma omp parallel for" in generated
    True
    >>> "d_f_5th_bundle[IDX_F(comp, i)] = f_5th_val;" in generated
    True
    """
    parallelization = par.parval_from_str("parallelization")
    if parallelization != "openmp":
        raise ValueError(
            "rkf45_5th_state currently supports only parallelization='openmp'."
        )

    # Step 1: Define C function metadata.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = ""
    postfunc = ""
    desc = r"""Build the RKF45 5th-order state for one active photon chunk.

    The helper evaluates the Cash-Karp 5th-order solution

    $f_{n+1}^{(5)} = f_n + h \sum_s b_s^{(5)} k_s$

    component-by-component for the normalized photon state bundle and writes the
    result into a scratch output bundle, leaving the persistent state unchanged.

    @param[in] d_f_start Pointer to the anchored start-state bundle $f_n$.
    @param[in] d_k_bundle Pointer to the flattened RKF45 stage-derivative bundle.
    @param[in] d_h Pointer to the per-ray RKF45 step sizes.
    @param[out] d_f_5th_bundle Pointer to the output scratch bundle storing $f_{n+1}^{(5)}$.
    @param[in] chunk_size Number of active rays in the current chunk.
    @param[in] stream_idx Unused placeholder kept for interface compatibility with other photon helpers.
    """
    name = "rkf45_5th_state"
    params = (
        "const double *restrict d_f_start, "
        "const double *restrict d_k_bundle, "
        "const double *restrict d_h, "
        "double *restrict d_f_5th_bundle, "
        "const long int chunk_size, "
        "const int stream_idx"
    )

    # Step 2: Build the OpenMP C body.
    body = r"""
  (void)stream_idx;

  //==========================================
  // MACRO DEFINITIONS FOR BUNDLE ACCESS
  //==========================================
  // IDX_F maps a component to the flattened state bundle using SoA layout.
  #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
  // IDX_K maps a stage/component pair into the flattened RKF45 derivative bundle.
  #define IDX_K(s, c, ray_id) ((s) * 9 * BUNDLE_CAPACITY + (c) * BUNDLE_CAPACITY + (ray_id))

  //==========================================
  // OPENMP LOOP ARCHITECTURE
  //==========================================
  // Distribute photon rays across available CPU threads for parallel evaluation.
  #pragma omp parallel for
  for (long int i = 0; i < chunk_size; i++) {
    // Step 1: Load the signed RKF45 step size for this ray.
    const double h = d_h[i];

    // Step 2: Build the 5th-order candidate state component-by-component.
    for (int comp = 0; comp < 9; ++comp) {
      const double f_n = d_f_start[IDX_F(comp, i)];
      const double k0 = d_k_bundle[IDX_K(0, comp, i)];
      // Bundle slot 1 carries zero weight in the Cash-Karp 5th-order solution.
      const double k2 = d_k_bundle[IDX_K(2, comp, i)];
      const double k3 = d_k_bundle[IDX_K(3, comp, i)];
      const double k4 = d_k_bundle[IDX_K(4, comp, i)];
      const double k5 = d_k_bundle[IDX_K(5, comp, i)];

      const double update_val =
          (16.0 / 135.0) * k0 +
          (6656.0 / 12825.0) * k2 +
          (28561.0 / 56430.0) * k3 +
          (-9.0 / 50.0) * k4 +
          (2.0 / 55.0) * k5;
      const double f_5th_val = f_n + h * update_val;

      d_f_5th_bundle[IDX_F(comp, i)] = f_5th_val;
    } // END LOOP: for comp over 9 state-vector components
  } // END LOOP: for i over chunk_size rays

  //==========================================
  // MACRO CLEANUP
  //==========================================
  #undef IDX_F
  #undef IDX_K
"""

    # Step 3: Register the host-only C helper.
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
