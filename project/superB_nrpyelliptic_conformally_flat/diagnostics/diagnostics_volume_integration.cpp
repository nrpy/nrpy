#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "diagnostic_gfs.h"
#include "diagnostics_volume_integration_helpers.h"

/**
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
 *
 */
void diagnostics_volume_integration(commondata_struct *restrict commondata, griddata_struct *restrict griddata, griddata_struct *restrict griddata_chare,
                                    const REAL *restrict gridfuncs_diags[MAXNUMGRIDS],
                                    const int chare_index[3],
                                    const int which_diagnostics_part) {
  // Initialize a recipe book.
  //~ diags_integration_recipe_t recipes[DIAGS_INTEGRATION_MAX_RECIPES];



  int which_grid = 0;
  diags_integration_recipe_t *const recipes = griddata_chare[which_grid].diagnosticstruct.recipes;
  diags_integration_results_t *const integration_results = griddata_chare[which_grid].diagnosticstruct.integration_results;

  switch (which_diagnostics_part) {

    case DIAGNOSTICS_VOLUME_EXECUTE_RECIPE_FOR_CHARE_GRID: {

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
          recipes[NUM_RECIPES].integrands[0] = (diags_integration_integrand_spec_t){.gf_index = DIAG_RESIDUAL, .is_squared = 1};
          recipes[NUM_RECIPES].num_integrands = 1;

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
          recipes[NUM_RECIPES].integrands[0] = (diags_integration_integrand_spec_t){.gf_index = DIAG_RESIDUAL, .is_squared = 1};
          recipes[NUM_RECIPES].num_integrands = 1;

          NUM_RECIPES++;
        } // END IF define Recipe 1

        // --- You can continue to add more recipes here, using any combination of rules and integrands ---
      } // END USER-EDIT recipes block

      // ========================= End of user edits =========================

      griddata[which_grid].diagnosticstruct.NUM_RECIPES = NUM_RECIPES;

      //~ diags_integration_execute_recipes(commondata, griddata, gridfuncs_diags, recipes, NUM_RECIPES, integration_results, chare_index, which_diagnostics_part);

      break;
    }

    case DIAGNOSTICS_VOLUME_WRITE: {

      int which_grid = 0;
      diags_integration_recipe_t *const recipes = griddata_chare[which_grid].diagnosticstruct.recipes;
      diags_integration_results_t *const integration_results = griddata_chare[which_grid].diagnosticstruct.integration_results;
      int NUM_RECIPES = griddata[which_grid].diagnosticstruct.NUM_RECIPES;

      //~ diags_integration_execute_recipes(commondata, griddata, gridfuncs_diags, recipes, NUM_RECIPES, integration_results, which_diagnostics_part);

      // ========================= USER-EDIT: Capture results =========================
      // Execute all defined recipes and capture results
      {
        //~ REAL rms_residual_r_80 = NAN;
        //~ diags_integration_get_rms(&integration_results, "sphere_R_80", DIAG_RESIDUAL, &rms_residual_r_80);
        //~ commondata->log10_current_residual = log10(rms_residual_r_80);
      } // END USER-EDIT capture results block
      // ========================= End of user edits =========================


      break;
    }

  }


} // END FUNCTION diagnostics_volume_integration
