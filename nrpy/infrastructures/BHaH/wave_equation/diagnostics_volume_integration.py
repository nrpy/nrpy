"""
Register and emit the C `diagnostics_volume_integration()` routine for recipe-driven volume integrals.

This module exposes a single entry point,
`register_CFunction_diagnostics_volume_integration()`, which generates and registers
the C function `diagnostics_volume_integration()`. The generated routine constructs a
small "recipe book" of integration definitions, allows user edits in a clearly marked
block, and then dispatches to an executor that performs domain selections, integrand
evaluations, numerical quadrature, and file I/O.

Recipe model:
  - A recipe is a named collection of:
      * Zero or more "rules" that include or exclude geometric regions.
        The emitted template demonstrates spherical rules using
        `diags_integration_sphere_rule_t` with fields:
          .center_xyz = {x0, x1, x2}, .radius, .exclude_inside (0 include, 1 exclude)
      * One or more "integrands" declared via
        `diags_integration_integrand_spec_t`:
          .gf_index (enum selecting a diagnostic gridfunction)
          .is_squared (0 for integral of f dV, 1 for integral of f^2 dV, enabling RMS/L2)
  - The third function parameter, `gridfuncs_diags[grid]`, is a caller-provided,
    gf-major pointer to REAL diagnostic gridfunctions available to the integrator.

The emitted template defines three example recipes:
  - "whole_domain_L2_norms": L2-like integrals (is_squared=1) of relative errors over
    the entire domain (no rules).
  - "shell_integrals": plain integrals (is_squared=0) of exact fields inside a
    spherical shell Rmin <= r <= Rmax, plus a unit integrand to measure total volume.
  - "outside_ball_RMS": RMS-like integrals (is_squared=1) of numerical fields outside
    a small excision sphere.

During code generation:
  - The routine is registered under the `diagnostics/` subdirectory.
  - It includes headers `BHaH_defines.h`, `BHaH_function_prototypes.h`, and
    `diagnostic_gfs.h`.

At runtime, `diagnostics_volume_integration()` initializes a fixed-size recipe array,
invokes a user-edit block to define recipes, and then calls
`diags_integration_execute_recipes(...)` to compute and write results. Memory for
diagnostic gridfunctions is owned by the caller; this routine does not allocate
or free those buffers.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_diagnostics_volume_integration() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Generate and register the C `diagnostics_volume_integration()` routine that builds and executes integration recipes.

    This routine:
      1) Checks for the parallel codegen registration phase. If active, it registers
           the function call and returns. Otherwise, it proceeds with C function generation.
      2) Constructs C function metadata (name, params, includes) and a body that:
           - Allocates a fixed-size array `diags_integration_recipe_t recipes[...]`.
           - Calls `diags_integration_initialize_recipes(recipes)` to zero/init state.
           - Enters a clearly marked USER-EDIT block where recipes are defined.
           - Calls `diags_integration_execute_recipes(commondata, griddata, gridfuncs_diags, recipes, NUM_RECIPES)`
             to perform the actual computations and I/O.
      3) Registers the function under `diagnostics/` and returns the updated NRPy environment.

    Parameters and ownership:
      - None (Python side). The generated C routine takes:
          * `commondata_struct *restrict commondata`
          * `griddata_struct *restrict griddata`
          * `const REAL *restrict gridfuncs_diags[MAXNUMGRIDS]`
        The caller owns and initializes `gridfuncs_diags[grid]` and any result buffers
        written by the executor; this generator does not allocate or free memory.

    :returns: The updated NRPy environment after registration.

    Doctests: TBD
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # --- C Function Registration ---
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "diagnostic_gfs.h",
        "diagnostics_volume_integration_helpers.h",
    ]

    desc = r"""
 * @brief Build user-defined diagnostic integration recipes and execute them.
 *
 * @details
 * This routine constructs a small "recipe book" describing spatial domains and
 * integrands, then dispatches them to the integration executor. Each recipe is a
 * named collection of:
 *
 *  - Rules: geometric include/exclude predicates. The template uses spherical
 *    rules defined by `diags_integration_sphere_rule_t` with fields:
 *      .center_xyz = {x0, x1, x2}, .radius, .exclude_inside (0 include, 1 exclude).
 *  - Integrands: selected by gridfunction enum index and an `is_squared` flag,
 *    using `diags_integration_integrand_spec_t`:
 *      .gf_index = DIAG_* token, .is_squared = 0 or 1.
 *    Setting `.is_squared = 1` enables L2/RMS-style integrals of f^2; otherwise
 *    plain integrals of f are computed.
 *
 * After defining recipes, the routine calls:
 *   diags_integration_execute_recipes(commondata, griddata, gridfuncs_diags, recipes, NUM_RECIPES);
 * which performs all domain checks, integral accumulation, and file I/O for each recipe.
 *
 * @param[in,out] commondata Shared simulation metadata and runtime context.
 * @param[in,out] griddata   Grid-specific data (coordinates, parameters, gridfuncs).
 * @param[in]     gridfuncs_diags Array of length MAXNUMGRIDS; entry `gridfuncs_diags[grid]`
 *                                points to REAL diagnostic gridfunction data available
 *                                to integrand evaluators.
 *
 * @pre
 *  - `gridfuncs_diags[grid]` is non-null and points to valid diagnostic gridfunctions.
 *  - Enum tokens referenced by `.gf_index` (for example DIAG_RELERROR_UUGF) are defined.
 *  - Integration helper APIs declared in included headers are available.
 *
 * @post
 *  - All defined recipes are processed; integrals are computed and written via the executor.
 *  - No memory is allocated or freed by this routine.
 *
 * @warning
 *  - The "USER-EDIT" block is intended for end-user customization. Ensure that
 *    `NUM_RECIPES` never exceeds `DIAGS_INTEGRATION_MAX_RECIPES`.
 *
 * @return void
 """
    cfunc_type = "void"
    name = "diagnostics_volume_integration"
    params = """commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                const REAL *restrict gridfuncs_diags[MAXNUMGRIDS]"""

    body = r"""
  // Initialize a recipe book.
  diags_integration_recipe_t recipes[DIAGS_INTEGRATION_MAX_RECIPES];
  diags_integration_initialize_recipes(recipes);
  int NUM_RECIPES = 0;

  // ========================= USER-EDIT: Define recipes =========================
  // Each recipe has its own set of rules and its own list of integrands.
  {
    // Recipe 0: Integrals over the whole domain
    if (NUM_RECIPES < DIAGS_INTEGRATION_MAX_RECIPES) {
      recipes[NUM_RECIPES].name = "whole_domain";
      recipes[NUM_RECIPES].num_rules = 0; // No rules means the whole domain is used.
      // Define the integrands for this recipe:

      // Important: is_squared=1 enables computation of L2 norm & RMS; RMS_f = sqrt(int f^2 dV / int dV)
      recipes[NUM_RECIPES].integrands[0] = (diags_integration_integrand_spec_t){.gf_index = DIAG_RELERROR_UUGF, .is_squared = 1};
      recipes[NUM_RECIPES].integrands[1] = (diags_integration_integrand_spec_t){.gf_index = DIAG_RELERROR_VVGF, .is_squared = 1};
      recipes[NUM_RECIPES].num_integrands = 2;

      NUM_RECIPES++;
    } // END IF define Recipe 0

    //  Recipe 1: RMS of numerical solution inside a sphere of radius 80
    if (NUM_RECIPES < DIAGS_INTEGRATION_MAX_RECIPES) {
      const REAL R_outer = 80;
      recipes[NUM_RECIPES].name = "sphere_R_80";
      recipes[NUM_RECIPES].num_rules = 1;

      // Important: exclude_inside=0 implies outer is excluded.
      recipes[NUM_RECIPES].rules[0] = (diags_integration_sphere_rule_t){.center_xyz = {0, 0, 0}, .radius = R_outer, .exclude_inside = 0};

      // Important: is_squared=1 enables computation of L2 norm & RMS; RMS_f = sqrt(int f^2 dV / int dV)
      recipes[NUM_RECIPES].integrands[0] = (diags_integration_integrand_spec_t){.gf_index = DIAG_RELERROR_UUGF, .is_squared = 1};
      recipes[NUM_RECIPES].integrands[1] = (diags_integration_integrand_spec_t){.gf_index = DIAG_RELERROR_VVGF, .is_squared = 1};
      recipes[NUM_RECIPES].num_integrands = 2;

      NUM_RECIPES++;
    } // END IF define Recipe 1

    // --- You can continue to add more recipes here, using any combination of rules and integrands ---
  } // END USER-EDIT recipes block

  // ========================= End of user edits =========================

  // Declare a container to hold the integration results
  diags_integration_results_t integration_results;
  diags_integration_execute_recipes(commondata, griddata, gridfuncs_diags, recipes, NUM_RECIPES, &integration_results);

  // ========================= USER-EDIT: Capture results =========================
  // Execute all defined recipes and capture results
  {
    REAL rms_relerror_uugf = NAN;
    diags_integration_get_rms(&integration_results, "whole_domain", DIAG_RELERROR_UUGF, &rms_relerror_uugf);
    printf("%e %e | time rms error of u\n", commondata->time, rms_relerror_uugf);
  } // END USER-EDIT capture results block
  // ========================= End of user edits =========================
"""
    cfc.register_CFunction(
        subdirectory="diagnostics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
