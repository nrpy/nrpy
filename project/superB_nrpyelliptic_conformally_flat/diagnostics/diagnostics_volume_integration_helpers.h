/**
 * @file diagnostics_volume_integration_helpers.h
 * @brief Provides a recipe-based interface to compute domain integrals of diagnostic gridfunctions.
 * @details This header provides a recipe-based interface to compute domain integrals
 * of diagnostic gridfunctions over spatial regions defined by spherical rules,
 * with powerful query capabilities for extracting results without dependency on
 * recipe or integrand ordering.
 *
 * @section usage Usage
 *
 * 1. Define recipes with spherical rules and integrands
 * 2. Call diags_integration_execute_recipes() with a results container
 * 3. Query results using convenience functions or generic query interface
 *
 * @code
 * diags_integration_results_t results;
 * diags_integration_execute_recipes(..., &results);
 *
 * REAL rms;
 * if (diags_integration_get_rms(&results, "whole_domain", DIAG_RESIDUAL, &rms)) {
 *   printf("RMS = %e\n", rms);
 * }
 * @endcode
 *
 * @section features Key Features
 * - Type-safe extraction via enums
 * - Order-independent queries by name
 * - Self-contained results with metadata
 * - Backward compatible (optional results parameter)
 * - Comprehensive error handling
 */

#ifndef DIAGNOSTICS_VOLUME_INTEGRATION_HELPERS_H
#define DIAGNOSTICS_VOLUME_INTEGRATION_HELPERS_H

#include "BHaH_defines.h"
#include "diagnostic_gfs.h" // diagnostic_gf_names[] and enum


//NEW
#include "diagnostics_volume_integration_helpers_defines.h"
//~ // ======================== Limits (stack-only) ========================
//~ // Upper bounds for small, fixed-size containers used in this module.
//~ enum {
  //~ DIAGS_INTEGRATION_MAX_RECIPES = 8,    // Maximum number of recipes processed at once.
  //~ DIAGS_INTEGRATION_MAX_RULES = 4,      // Maximum number of spherical rules per recipe.
  //~ DIAGS_INTEGRATION_MAX_INTEGRANDS = 16 // Maximum number of integrands per recipe.
//~ };

//~ // Special GF index indicating the unit integrand (proper volume element).
//~ #define diags_integration_UNIT_INTEGRAND_GFINDEX (-1)

//~ // ======================== Rule / Recipe types ========================

//~ // Per-integrand specification.
//~ typedef struct {
  //~ int gf_index;   // Diagnostic GF enum value, or diags_integration_UNIT_INTEGRAND_GFINDEX (-1) for the unit integrand.
  //~ int is_squared; // If nonzero, integrate f^2 and also print L2 and RMS.
//~ } diags_integration_integrand_spec_t;

//~ // Spherical include/exclude rule in Cartesian space.
//~ typedef struct {
  //~ REAL center_xyz[3]; // Sphere center (x,y,z) in Cartesian coordinates.
  //~ REAL radius;        // Sphere radius (>= 0).
  //~ int exclude_inside; // If 1, exclude r <= radius (inside); if 0, exclude r > radius (outside).
//~ } diags_integration_sphere_rule_t;

//~ // Describes an integration domain and the integrands to measure.
//~ typedef struct {
  //~ // Rules that define the integration domain
  //~ int num_rules;                                                      // Number of active rules (0..DIAGS_INTEGRATION_MAX_RULES).
  //~ diags_integration_sphere_rule_t rules[DIAGS_INTEGRATION_MAX_RULES]; // Rule array (only first num_rules are used).

  //~ // Integrands to compute over that domain
  //~ int num_integrands;                                                              // Number of integrands (0..DIAGS_INTEGRATION_MAX_INTEGRANDS).
  //~ diags_integration_integrand_spec_t integrands[DIAGS_INTEGRATION_MAX_INTEGRANDS]; // Integrand specifications.

  //~ const char *name; // Mandatory recipe name; used to construct the output filename.
//~ } diags_integration_recipe_t;

//~ // ======================== Result query types ========================

//~ // Specifies the type of value to extract from integration results.
//~ typedef enum {
  //~ DIAGS_EXTRACT_VOLUME,   // Remaining proper volume after rule application
  //~ DIAGS_EXTRACT_INTEGRAL, // Raw integral: integral of f dV or integral of f^2 dV
  //~ DIAGS_EXTRACT_L2,       // L2 norm: sqrt(integral of f^2 dV), requires is_squared=1
  //~ DIAGS_EXTRACT_RMS       // RMS: sqrt(integral of f^2 dV / V), requires is_squared=1
//~ } diags_integration_extraction_type_t;

//~ // Results for a single integrand within a recipe.
//~ typedef struct {
  //~ int gf_index;        // Diagnostic GF enum value (or diags_integration_UNIT_INTEGRAND_GFINDEX)
  //~ const char *gf_name; // Human-readable name (pointer to diagnostic_gf_names[] or "one")
  //~ int is_squared;      // Whether this integrand was squared (0 or 1)
  //~ REAL integral;       // integral of f dV if not squared, integral of f^2 dV if squared
  //~ REAL L2;             // sqrt(max(integral of f^2 dV, 0)); only valid if is_squared=1
  //~ REAL RMS;            // sqrt(max(integral of f^2 dV / V, 0)); only valid if is_squared=1
  //~ int L2_valid;        // 1 if L2 value is meaningful, 0 otherwise
  //~ int RMS_valid;       // 1 if RMS value is meaningful, 0 otherwise
//~ } diags_integration_integrand_result_t;

//~ // Results for one recipe execution.
//~ typedef struct {
  //~ const char *recipe_name;   // Name of the recipe (pointer to user's string)
  //~ REAL proper_volume;        // Remaining volume after applying all rules
  //~ int num_integrand_results; // Number of valid integrand results
  //~ diags_integration_integrand_result_t integrand_results[DIAGS_INTEGRATION_MAX_INTEGRANDS];
//~ } diags_integration_recipe_result_t;

//~ // Container for all recipe results from one execution.
//~ typedef struct {
  //~ int num_recipe_results; // Number of recipes executed
  //~ diags_integration_recipe_result_t recipe_results[DIAGS_INTEGRATION_MAX_RECIPES];
//~ } diags_integration_results_t;

// ======================== Small helpers (inline) ========================

/**
 * @brief Tests whether a point is excluded by the first `prefix_rule_count` rules.
 *
 * @param[in] xx  Point x-coordinate in Cartesian space.
 * @param[in] yy  Point y-coordinate in Cartesian space.
 * @param[in] zz  Point z-coordinate in Cartesian space.
 * @param[in] prefix_rule_count Number of rules from the beginning of @p rules to apply.
 * @param[in] rules Array of spherical include/exclude rules.
 * @return `true` if any applied rule excludes the point; `false` otherwise.
 */
static inline bool diags_integration_point_is_excluded_prefix(const REAL xx, const REAL yy, const REAL zz, int prefix_rule_count,
                                                              const diags_integration_sphere_rule_t *restrict rules) {
  for (int rule = 0; rule < prefix_rule_count; rule++) {
    const REAL dx = xx - rules[rule].center_xyz[0];
    const REAL dy = yy - rules[rule].center_xyz[1];
    const REAL dz = zz - rules[rule].center_xyz[2];
    const REAL rr = sqrt(dx * dx + dy * dy + dz * dz);
    if (rules[rule].exclude_inside) {
      if (rr <= rules[rule].radius)
        return true;
    } // END IF for exclude_inside
    else {
      if (rr > rules[rule].radius)
        return true;
    } // END ELSE for exclude_inside
  } // END LOOP over rules
  return false;
} // END FUNCTION diags_integration_point_is_excluded_prefix

// ======================== File I/O helpers ========================

/**
 * @brief Opens (or appends) the per-recipe output file.
 *
 * @details Opens in write mode (`"w"`) when `nn == 0`, otherwise appends (`"a"`).
 * On failure, prints an error and exits the process.
 *
 * @param[in] nn        Global iteration counter (from `commondata_struct`).
 * @param[in] filename  Path to the file to open.
 * @return File pointer opened for text output (never `NULL`; exits on error).
 */
static inline FILE *diags_integration_open_output_file(const int nn, const char *filename) {
  FILE *file_ptr = fopen(filename, (nn == 0) ? "w" : "a");
  if (file_ptr == NULL) {
    fprintf(stderr, "diagnostics_integration: could not open %s\n", filename);
    exit(1);
  } // END IF file pointer is NULL
  return file_ptr;
} // END FUNCTION diags_integration_open_output_file

/**
 * @brief Write a per-recipe header describing column mappings, rules, and integrands.
 *
 * @details Every header line begins with '#', so plotting and parsing tools can ignore or
 * consume it as metadata. The header documents:
 *   (A) The recipe (index & name).
 *   (B) Each spherical rule with index and an auto-generated rule name.
 *   (C) Exact 1-based column indices and contiguous batches for each integrand.
 *   (D) A per-prefix guide that embeds a "prefix:" tag into each column label so
 *       rows (which differ by rules_applied = k) can be mapped to correct labels.
 *   (E) Two compact, machine-friendly lines:
 *         - "# COLUMNS:" : explicit tokens for all columns, with a "prefix:" tag.
 *         - "# time ..." : legacy-style one-liner for drop-in compatibility (also with "prefix:" tags).
 *
 * @param[in] file_ptr      Open file pointer to which the header is written.
 * @param[in] recipe_index  Index of the recipe in the current execution batch.
 * @param[in] recipe        Pointer to the recipe descriptor (name, rules, integrands).
 * @return void
 */
static inline void diags_integration_write_header(FILE *file_ptr, int recipe_index, const diags_integration_recipe_t *restrict recipe) {
  // Defensive extraction
  const char *recipe_name = (recipe && recipe->name) ? recipe->name : "(null)";
  const int num_rules = (recipe) ? recipe->num_rules : 0;
  const int num_specs = (recipe) ? recipe->num_integrands : 0;

  // ---------------------------------------------------------------------------
  // Section A: Recipe summary
  // ---------------------------------------------------------------------------
  fprintf(file_ptr, "# diagnostics_integration: column mapping & rule documentation\n");
  fprintf(file_ptr, "# recipe_index : %d\n", recipe_index);
  fprintf(file_ptr, "# recipe_name  : %s\n", recipe_name);
  fprintf(file_ptr, "# row semantics: only the *full* rules-applied row is emitted (k = num_rules)\n");

  // ---------------------------------------------------------------------------
  // Section B: Spherical rules (index + auto-generated rule name + geometry)
  // ---------------------------------------------------------------------------
  if (num_rules > 0) {
    for (int rule = 0; rule < num_rules; rule++) {
      const diags_integration_sphere_rule_t *ruleptr = &recipe->rules[rule];
      const int excl_in = ruleptr->exclude_inside;

      // Auto-generated rule name:
      //   exclude_inside==1 => EXCLUDE inside (KEEP OUTSIDE)  : "KeepOutside_r>R_at(x,y,z)"
      //   exclude_inside==0 => EXCLUDE outside (KEEP INSIDE)  : "KeepInside_r<=R_at(x,y,z)"
      char rulename[256];
      if (excl_in) {
        snprintf(rulename, sizeof(rulename), "KeepOutside_r>%.17g_at(%.17g,%.17g,%.17g)", (double)ruleptr->radius, (double)ruleptr->center_xyz[0],
                 (double)ruleptr->center_xyz[1], (double)ruleptr->center_xyz[2]);
      } else {
        snprintf(rulename, sizeof(rulename), "KeepInside_r<=%.17g_at(%.17g,%.17g,%.17g)", (double)ruleptr->radius, (double)ruleptr->center_xyz[0],
                 (double)ruleptr->center_xyz[1], (double)ruleptr->center_xyz[2]);
      } // END IF-ELSE for rule name generation

      fprintf(file_ptr, "#   rule[%d] name=\"%s\" : center=(%.17g, %.17g, %.17g), radius=%.17g, action=%s\n", rule, rulename,
              (double)ruleptr->center_xyz[0], (double)ruleptr->center_xyz[1], (double)ruleptr->center_xyz[2], (double)ruleptr->radius,
              excl_in ? "exclude INSIDE (keep OUTSIDE: r > radius)" : "exclude OUTSIDE (keep INSIDE: r <= radius)");
    } // END LOOP over rules
  } else {
    fprintf(file_ptr, "#   (no rules; the entire computational domain is used)\n");
  } // END ELSE for num_rules > 0

  // ---------------------------------------------------------------------------
  // Section C: Column schema (1-based indices; contiguous integrand batches)
  // ---------------------------------------------------------------------------
  int col = 1;
  fprintf(file_ptr, "#   %3d: time                    (current simulation time; commondata->time)\n", col++);
  fprintf(file_ptr, "#   %3d: remaining_volume        (proper volume after applying the prefix rules)\n", col++);

  // Prepare compact one-liners (without recipe_idx and rules_applied):
  char compact_columns[8192];
  char legacy_columns[8192];
  snprintf(compact_columns, sizeof(compact_columns), "# COLUMNS: time remaining_volume");
  snprintf(legacy_columns, sizeof(legacy_columns), "# time remaining_volume");

  const int num_specs_local = num_specs;
  if (num_specs_local > 0) {
    for (int spec = 0; spec < num_specs_local; spec++) {
      const diags_integration_integrand_spec_t *specptr = &recipe->integrands[spec];
      const char *nm = (specptr->gf_index == diags_integration_UNIT_INTEGRAND_GFINDEX) ? "one" : diagnostic_gf_names[specptr->gf_index];

      const int col_INT = col++;
      if (specptr->is_squared) {
        // INT[name^2] for clarity that it is integral of (name)^2 dV
        fprintf(file_ptr, "#   %3d: INT[prefix:%s^2]      (integral of (%s)^2 dV)\n", col_INT, nm, nm);
      } else {
        fprintf(file_ptr, "#   %3d: INT[prefix:%s]        (integral of %s dV)\n", col_INT, nm, nm);
      } // END IF-ELSE for squared integral printing

      // Extend compact/legacy one-liners
      {
        size_t used = strlen(compact_columns);
        if (specptr->is_squared)
          snprintf(compact_columns + used, sizeof(compact_columns) - used, " INT[prefix:%s^2]", nm);
        else
          snprintf(compact_columns + used, sizeof(compact_columns) - used, " INT[prefix:%s]", nm);

        used = strlen(legacy_columns);
        if (specptr->is_squared)
          snprintf(legacy_columns + used, sizeof(legacy_columns) - used, " INT[prefix:%s^2]", nm);
        else
          snprintf(legacy_columns + used, sizeof(legacy_columns) - used, " INT[prefix:%s]", nm);
      } // END SCOPE for one-liners

      if (specptr->is_squared) {
        const int col_L2 = col++;
        const int col_RMS = col++;
        fprintf(file_ptr, "#   %3d: L2[prefix:%s]         (sqrt( max( INT[prefix:%s^2], 0 ) ))\n", col_L2, nm, nm);
        fprintf(file_ptr, "#   %3d: RMS[prefix:%s]        ((remaining_volume>0) ? sqrt( max( INT[prefix:%s^2], 0 ) / remaining_volume ) : 0)\n",
                col_RMS, nm, nm);

        // Extend compact/legacy one-liners
        size_t used = strlen(compact_columns);
        snprintf(compact_columns + used, sizeof(compact_columns) - used, " L2[prefix:%s] RMS[prefix:%s]", nm, nm);
        used = strlen(legacy_columns);
        snprintf(legacy_columns + used, sizeof(legacy_columns) - used, " L2[prefix:%s] RMS[prefix:%s]", nm, nm);
      } // END IF for is_squared extras
    } // END LOOP over integrands
  } else {
    fprintf(file_ptr, "#   (no integrands defined for this recipe)\n");
  } // END ELSE for num_specs_local > 0

  // ---------------------------------------------------------------------------
  // Section D: Prefix mapping guide (embeds the actual prefix tag per k)
  // ---------------------------------------------------------------------------
  fprintf(file_ptr, "# --- Prefix mapping (how to interpret \"prefix:\" per row) ----------------------\n");
  fprintf(file_ptr, "#   For a row with rules_applied = k, substitute:\n");
  fprintf(file_ptr, "#     k = 0 : prefix = p0\n");
  if (num_rules > 0) {
    for (int k0 = 1; k0 <= num_rules; k0++) {
      const diags_integration_sphere_rule_t *ruleptr = &recipe->rules[k0 - 1];
      const int excl_in = ruleptr->exclude_inside;
      char rulename[256];
      if (excl_in) {
        snprintf(rulename, sizeof(rulename), "KeepOutside_r>%.17g_at(%.17g,%.17g,%.17g)", (double)ruleptr->radius, (double)ruleptr->center_xyz[0],
                 (double)ruleptr->center_xyz[1], (double)ruleptr->center_xyz[2]);
      } else {
        snprintf(rulename, sizeof(rulename), "KeepInside_r<=%.17g_at(%.17g,%.17g,%.17g)", (double)ruleptr->radius, (double)ruleptr->center_xyz[0],
                 (double)ruleptr->center_xyz[1], (double)ruleptr->center_xyz[2]);
      } // END IF-ELSE for rule name generation
      fprintf(file_ptr, "#     k = %d : prefix = p%d|rule[%d]:%s\n", k0, k0, k0 - 1, rulename);
    } // END LOOP for prefix mapping
  } // END IF for prefix mapping

  // ---------------------------------------------------------------------------
  // Section E: Compact machine-friendly lines
  // ---------------------------------------------------------------------------
  fprintf(file_ptr, "%s\n", compact_columns); // "# COLUMNS: ..."
  fprintf(file_ptr, "%s\n", legacy_columns);  // "# time ..."

  // Footer delimiter lines with dashes are intentionally omitted.
} // END FUNCTION diags_integration_write_header

/**
 * @brief Writes one data row for a specific recipe prefix at a given time.
 *
 * @param[in] file_ptr      Open file pointer.
 * @param[in] time          Simulation time (printed with `%.17g`).
 * @param[in] volume        Remaining proper volume after applying the prefix rules.
 * @param[in] specs         Integrand specifications for column labeling and `L2`/`RMS` logic.
 * @param[in] integrals     Per-integrand integrals (`integral of f dV` or `integral of f^2 dV` depending on `is_squared`).
 * @param[in] num_specs     Number of entries in @p specs and @p integrals.
 * @return void
 */
static inline void diags_integration_write_row(FILE *file_ptr, double time, double volume, const diags_integration_integrand_spec_t *restrict specs,
                                               const REAL *restrict integrals, int num_specs) {
  fprintf(file_ptr, "%.15e %.15e", time, volume);
  for (int spec = 0; spec < num_specs; spec++) {
    const double integralvalue = (double)integrals[spec]; // integral of f dV or integral of f^2 dV depending on is_squared
    fprintf(file_ptr, " %.15e", integralvalue);
    if (specs[spec].is_squared) {
      const double L2 = sqrt(fmax(integralvalue, 0.0));                                  // ||f||_2
      const double RMS = (volume > 0.0) ? sqrt(fmax(integralvalue, 0.0) / volume) : 0.0; // sqrt(mean(f^2))
      fprintf(file_ptr, " %.15e %.15e", L2, RMS);
    } // END IF is_squared
  } // END LOOP over specs
  fprintf(file_ptr, "\n");
} // END FUNCTION diags_integration_write_row

/**
 * @brief Apply a prefix of spherical rules and integrate requested integrands.
 *
 * @details For each grid and each in-domain cell, this function:
 *  1) Computes the proper volume element dV via sqrt_detgammahat_d3xx_volume_element.
 *  2) Converts curvilinear coordinates xx to Cartesian xCart for rule tests.
 *  3) Optionally skips points excluded by the rule prefix.
 *  4) Accumulates proper volume and per-integrand integrals (squaring f if requested).
 *
 * OpenMP parallelizes the outermost 3D loop with array-section reductions.
 *
 * @param[in]  commondata         Global/common simulation metadata (e.g., NUMGRIDS, time).
 * @param[in]  griddata           Per-grid geometry, coordinates, and masks.
 * @param[in]  diagnostic_gfs_major     For each grid function index, major-order diagnostic GF array.
 * @param[in]  rules              Array of spherical rules.
 * @param[in]  prefix_rule_count  Number of leading rules from @p rules to apply as an exclusion prefix.
 * @param[in]  integrands         Array of integrand specifications (gf indices and square flags).
 * @param[in]  num_integrands     Number of entries in @p integrands.
 * @param[out] proper_volume_out  Remaining proper volume after applying the rules.
 * @param[out] integrals_out      Output array (size >= num_integrands) receiving integral of (f or f^2) dV.
 * @return void
 */
static void diagnostics_integration_apply_rules(const commondata_struct *restrict commondata, const griddata_struct *restrict griddata,
                                                const REAL *restrict diagnostic_gfs_major[MAXNUMGRIDS],
                                                const diags_integration_sphere_rule_t *restrict rules, const int prefix_rule_count,
                                                const diags_integration_integrand_spec_t *restrict integrands, const int num_integrands,
                                                REAL *restrict proper_volume_out, REAL *restrict integrals_out) {
  const int NUMGRIDS = commondata->NUMGRIDS;

  // Accumulators reduced across all threads (and across all grids via outer loop)
  REAL volume_sum = 0.0;
  REAL integrals_sum[DIAGS_INTEGRATION_MAX_INTEGRANDS];
  for (int k0 = 0; k0 < num_integrands; k0++)
    integrals_sum[k0] = 0.0;

  // Precompute per-integrand metadata once (GF indices and flags).
  // NOTE: gf_index == diags_integration_UNIT_INTEGRAND_GFINDEX (-1) means "unit" integrand (value 1).
  int gf_index_arr[DIAGS_INTEGRATION_MAX_INTEGRANDS];
  unsigned char is_unit[DIAGS_INTEGRATION_MAX_INTEGRANDS];
  unsigned char is_squared[DIAGS_INTEGRATION_MAX_INTEGRANDS];
  for (int k0 = 0; k0 < num_integrands; k0++) {
    const int gf_index = integrands[k0].gf_index;
    const unsigned char unit = (gf_index == diags_integration_UNIT_INTEGRAND_GFINDEX);
    gf_index_arr[k0] = gf_index;
    is_unit[k0] = unit;
    is_squared[k0] = (unsigned char)(integrands[k0].is_squared ? 1 : 0);
  }

  for (int grid = 0; grid < NUMGRIDS; grid++) {
    const params_struct *restrict params = &griddata[grid].params;
    SET_NXX_PLUS_2NGHOSTS_VARS(grid);

    // const int8_t *restrict mask = griddata[grid].mask;
    // const int8_t curr_grid = (int8_t)grid;

    // Current grid's diagnostic GF-major slab (must be non-null).
    const REAL *restrict diag = diagnostic_gfs_major[grid];
#ifdef DEBUG
    if (diag == NULL) {
      fprintf(stderr, "diagnostics_integration_apply_rules: diagnostic_gfs_major[%d] is NULL\n", grid);
      exit(1);
    } // END IF diag is NULL
#endif

    // Parallelize the 3D sweep with reductions
#pragma omp parallel for schedule(static) reduction(+ : volume_sum, integrals_sum[ : DIAGS_INTEGRATION_MAX_INTEGRANDS])
    for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
      const REAL xx2 = griddata[grid].xx[2][i2];
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
        const REAL xx1 = griddata[grid].xx[1][i1];
        for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
          const int idx3 = IDX3(i0, i1, i2);
          const REAL xx0 = griddata[grid].xx[0][i0];

          // Skip cells not owned by this grid
          // if (blah)
          //   continue;

          // Build coordinates once
          const REAL xx[3] = {xx0, xx1, xx2};
          REAL xCart[3];
          xx_to_Cart(params, xx, xCart);

          // Early exclusion by prefix rules
          if (prefix_rule_count > 0 && diags_integration_point_is_excluded_prefix(xCart[0], xCart[1], xCart[2], prefix_rule_count, rules)) {
            continue;
          } // END IF point is excluded

          // Volume element; skip if zero
          REAL dV;
          sqrt_detgammahat_d3xx_volume_element(params, xx0, xx1, xx2, &dV);
          if (dV == (REAL)0.0)
            continue;

          volume_sum += dV;

          // Accumulate all requested integrands for this cell
          for (int k0 = 0; k0 < num_integrands; k0++) {
            REAL val;
            if (is_unit[k0]) {
              val = (REAL)1.0;
            } else {
              const int gf_index = gf_index_arr[k0];
#ifdef DEBUG
              if (gf_index < 0 || gf_index >= TOTAL_NUM_DIAG_GFS) {
                fprintf(stderr, "diagnostics_integration_apply_rules: gf_index %d out of range [0,%d)\n", gf_index, (int)TOTAL_NUM_DIAG_GFS);
                exit(1);
              } // END IF gf_index is out of range
#endif
              val = diag[IDX4pt(gf_index, idx3)];
            } // END IF-ELSE for unit integrand

            if (is_squared[k0])
              val *= val;
            integrals_sum[k0] += val * dV;
          } // END LOOP over integrands
        } // END LOOP over i0
      } // END LOOP over i1
    } // END LOOP over i2
  } // END LOOP over grids

  *proper_volume_out = volume_sum;
  for (int k0 = 0; k0 < num_integrands; k0++)
    integrals_out[k0] = integrals_sum[k0];
} // END FUNCTION diagnostics_integration_apply_rules

/**
 * @brief Initializes an array of recipes to a default empty state.
 *
 * @param[out] recipes Array of size #DIAGS_INTEGRATION_MAX_RECIPES to initialize.
 *
 * @post For each entry: `name=NULL`, `num_rules=0`, `num_integrands=0`,
 *       rule centers set to `(0,0,0)`, `radius=0`, and `exclude_inside=0`.
 * @return void
 */
static inline void diags_integration_initialize_recipes(diags_integration_recipe_t recipes[DIAGS_INTEGRATION_MAX_RECIPES]) {
  for (int recipe = 0; recipe < DIAGS_INTEGRATION_MAX_RECIPES; recipe++) {
    recipes[recipe].name = NULL;
    recipes[recipe].num_rules = 0;
    recipes[recipe].num_integrands = 0;
    for (int rule = 0; rule < DIAGS_INTEGRATION_MAX_RULES; rule++) {
      recipes[recipe].rules[rule].center_xyz[0] = 0.0;
      recipes[recipe].rules[rule].center_xyz[1] = 0.0;
      recipes[recipe].rules[rule].center_xyz[2] = 0.0;
      recipes[recipe].rules[rule].radius = 0.0;
      recipes[recipe].rules[rule].exclude_inside = 0;
    } // END LOOP over rules
  } // END LOOP over recipes
} // END FUNCTION diags_integration_initialize_recipes

// ======================== Result query functions ========================

/**
 * @brief Find a recipe result by name.
 *
 * @param[in]  results           Pointer to the overall results structure.
 * @param[in]  recipe_name       Name of the recipe to find.
 * @param[out] recipe_result_out Pointer to receive the recipe result pointer.
 * @return 1 if found, 0 otherwise.
 *
 * @note This function performs a linear search. For small numbers of recipes
 *       (< 10), this is more efficient than hash table overhead.
 */
static inline int diags_integration_get_recipe_result(const diags_integration_results_t *restrict results, const char *restrict recipe_name,
                                                      const diags_integration_recipe_result_t **recipe_result_out) {

  // Defensive checks
  if (!results || !recipe_name || !recipe_result_out) {
    return 0;
  } // END IF for defensive checks

  for (int i0 = 0; i0 < results->num_recipe_results; i0++) {
    if (results->recipe_results[i0].recipe_name && strcmp(results->recipe_results[i0].recipe_name, recipe_name) == 0) {
      *recipe_result_out = &results->recipe_results[i0];
      return 1;
    } // END IF recipe found
  } // END LOOP over recipe results
  return 0;
} // END FUNCTION diags_integration_get_recipe_result

/**
 * @brief Find an integrand result within a recipe by GF index.
 *
 * @param[in]  recipe_result        Pointer to a recipe result.
 * @param[in]  gf_index             Diagnostic GF index to find.
 * @param[out] integrand_result_out Pointer to receive the integrand result pointer.
 * @return 1 if found, 0 otherwise.
 */
static inline int diags_integration_get_integrand_result(const diags_integration_recipe_result_t *restrict recipe_result, int gf_index,
                                                         const diags_integration_integrand_result_t **integrand_result_out) {

  if (!recipe_result || !integrand_result_out) {
    return 0;
  } // END IF for defensive checks

  for (int i0 = 0; i0 < recipe_result->num_integrand_results; i0++) {
    if (recipe_result->integrand_results[i0].gf_index == gf_index) {
      *integrand_result_out = &recipe_result->integrand_results[i0];
      return 1;
    } // END IF integrand found
  } // END LOOP over integrand results
  return 0;
} // END FUNCTION diags_integration_get_integrand_result

/**
 * @brief Extract a specific calculated value from the results.
 *
 * @param[in]  results          Pointer to the overall results structure.
 * @param[in]  recipe_name      Name of the recipe.
 * @param[in]  gf_index         Diagnostic GF index.
 * @param[in]  extraction_type  Type of value to extract (enum).
 * @param[out] value_out        Pointer to receive the extracted value.
 * @return 1 if successful, 0 if not found or invalid request.
 *
 * @note If extracting L2 or RMS, the integrand must have is_squared=1.
 *       The function will return 0 if this constraint is violated.
 */
static inline int diags_integration_query(const diags_integration_results_t *restrict results, const char *restrict recipe_name, int gf_index,
                                          diags_integration_extraction_type_t extraction_type, REAL *restrict value_out) {

  if (!results || !recipe_name || !value_out) {
    return 0;
  } // END IF for defensive checks

  // Step 1: Find recipe
  const diags_integration_recipe_result_t *recipe_result;
  if (!diags_integration_get_recipe_result(results, recipe_name, &recipe_result)) {
    return 0; // Recipe not found
  } // END IF recipe not found

  // Step 2: Handle VOLUME extraction (doesn't need integrand)
  if (extraction_type == DIAGS_EXTRACT_VOLUME) {
    *value_out = recipe_result->proper_volume;
    return 1;
  } // END IF extracting volume

  // Step 3: Find integrand
  const diags_integration_integrand_result_t *integrand_result;
  if (!diags_integration_get_integrand_result(recipe_result, gf_index, &integrand_result)) {
    return 0; // Integrand not found
  } // END IF integrand not found

  // Step 4: Extract requested quantity
  switch (extraction_type) {
  case DIAGS_EXTRACT_INTEGRAL:
    *value_out = integrand_result->integral;
    return 1;

  case DIAGS_EXTRACT_L2:
    if (!integrand_result->L2_valid) {
      return 0; // L2 only valid for squared integrands
    } // END IF L2 is not valid
    *value_out = integrand_result->L2;
    return 1;

  case DIAGS_EXTRACT_RMS:
    if (!integrand_result->RMS_valid) {
      return 0; // RMS only valid for squared integrands
    } // END IF RMS is not valid
    *value_out = integrand_result->RMS;
    return 1;

  case DIAGS_EXTRACT_VOLUME:
    // Already handled above
    return 0;

  default:
    return 0; // Unknown extraction type
  } // END SWITCH on extraction_type
} // END FUNCTION diags_integration_query

// Helpful diags_integration_query() macros
#define diags_integration_get_rms(results, recipe_name, gf_index, out_ptr)                                                                           \
  diags_integration_query((results), (recipe_name), (gf_index), DIAGS_EXTRACT_RMS, (out_ptr))
#define diags_integration_get_integral(results, recipe_name, gf_index, out_ptr)                                                                      \
  diags_integration_query((results), (recipe_name), (gf_index), DIAGS_EXTRACT_INTEGRAL, (out_ptr))
#define diags_integration_get_l2(results, recipe_name, gf_index, out_ptr)                                                                            \
  diags_integration_query((results), (recipe_name), (gf_index), DIAGS_EXTRACT_L2, (out_ptr))
#define diags_integration_get_volume(results, recipe_name, out_ptr)                                                                                  \
  diags_integration_query((results), (recipe_name), -1, DIAGS_EXTRACT_VOLUME, (out_ptr))

/**
 * @brief Check if a specific result exists without retrieving it.
 *
 * @param[in] results     Pointer to results structure.
 * @param[in] recipe_name Name of the recipe.
 * @param[in] gf_index    Diagnostic GF index.
 * @return 1 if the integrand exists in the recipe, 0 otherwise.
 */
static inline int diags_integration_result_exists(const diags_integration_results_t *restrict results, const char *restrict recipe_name,
                                                  int gf_index) {

  const diags_integration_recipe_result_t *recipe_result;
  if (!diags_integration_get_recipe_result(results, recipe_name, &recipe_result)) {
    return 0;
  } // END IF recipe not found

  const diags_integration_integrand_result_t *integrand_result;
  return diags_integration_get_integrand_result(recipe_result, gf_index, &integrand_result);
} // END FUNCTION diags_integration_result_exists

/**
 * @brief Print all available results for debugging.
 *
 * @param[in] results Pointer to results structure.
 * @param[in] file_ptr      File pointer (use stdout for console output).
 * @return void
 */
static inline void diags_integration_print_all_results(const diags_integration_results_t *restrict results, FILE *restrict file_ptr) {

  if (!results || !file_ptr) {
    return;
  } // END IF for defensive checks

  fprintf(file_ptr, "\n=== Diagnostics Integration Results ===\n");
  fprintf(file_ptr, "Total recipes: %d\n\n", results->num_recipe_results);

  for (int i0 = 0; i0 < results->num_recipe_results; i0++) {
    const diags_integration_recipe_result_t *recipe = &results->recipe_results[i0];
    fprintf(file_ptr, "Recipe[%d]: %s\n", i0, recipe->recipe_name);
    fprintf(file_ptr, "  Volume: %.15e\n", (double)recipe->proper_volume);
    fprintf(file_ptr, "  Integrands: %d\n", recipe->num_integrand_results);

    for (int j0 = 0; j0 < recipe->num_integrand_results; j0++) {
      const diags_integration_integrand_result_t *integ = &recipe->integrand_results[j0];
      fprintf(file_ptr, "    [%d] %s (gf_index=%d, squared=%d):\n", j0, integ->gf_name, integ->gf_index, integ->is_squared);
      fprintf(file_ptr, "        INT  = %.15e\n", (double)integ->integral);
      if (integ->L2_valid) {
        fprintf(file_ptr, "        L2   = %.15e\n", (double)integ->L2);
      } // END IF L2 is valid
      if (integ->RMS_valid) {
        fprintf(file_ptr, "        RMS  = %.15e\n", (double)integ->RMS);
      } // END IF RMS is valid
    } // END LOOP over integrand results
    fprintf(file_ptr, "\n");
  } // END LOOP over recipe results
  fprintf(file_ptr, "======================================\n\n");
} // END FUNCTION diags_integration_print_all_results

/**
 * @brief Executes a list of recipes, writes results to files, and optionally populates a results structure.
 *
 * @details For each recipe:
 *   - Validates the `name` and `num_integrands`.
 *   - Opens (or appends) `out3d-integrals-conv_factor<conv_factor>-<name>.txt`.
 *   - Writes a header at the first timestep (`nn == 0`).
 *   - Computes integrals over the domain defined by all rules.
 *   - Writes one row to the output file.
 *   - If `results_out` is provided, populates it with all computed quantities.
 *
 * @param[in,out] commondata             Global simulation metadata (reads `nn`, `time`).
 * @param[in,out] griddata               Per-grid geometry, coordinates, and masks.
 * @param[in]     diagnostic_gfs_major   For each grid, major-order diagnostic GF array.
 * @param[in]     recipes                Array of recipe descriptors.
 * @param[in]     num_recipes            Number of recipes to execute.
 * @param[out]    results_out            (Optional) Pointer to a results structure. If NULL, results
 *                                       are only written to files. If non-NULL, structure is populated
 *                                       with all computed values for programmatic access.
 * @return void
 */
//~ static inline void diags_integration_execute_recipes(commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                                                     //~ const REAL *restrict diagnostic_gfs_major[MAXNUMGRIDS],
                                                     //~ const diags_integration_recipe_t *recipes, int num_recipes,
                                                     //~ diags_integration_results_t *restrict results_out, const int chare_index[3],const int which_diagnostics_part) {

  //~ switch (which_diagnostics_part) {

    //~ case DIAGNOSTICS_VOLUME_EXECUTE_RECIPE_FOR_CHARE_GRID: {

      //~ // Initialize results if provided
      //~ if (results_out) {
        //~ results_out->num_recipe_results = 0;
      //~ } // END IF results_out is provided

      //~ for (int recipe = 0; recipe < num_recipes; recipe++) {
        //~ if (recipes[recipe].name == NULL || recipes[recipe].name[0] == '\0') {
          //~ fprintf(stderr, "diags_integration_execute_recipes: recipe %d is missing a name, skipping.\n", recipe);
          //~ continue;
        //~ } // END IF recipe name is invalid
        //~ if (recipes[recipe].num_integrands <= 0 || recipes[recipe].num_integrands > DIAGS_INTEGRATION_MAX_INTEGRANDS) {
          //~ fprintf(stderr, "diags_integration_execute_recipes: recipe '%s' has an invalid number of integrands (%d), skipping.\n", recipes[recipe].name,
                  //~ recipes[recipe].num_integrands);
          //~ continue;
        //~ } // END IF integrand count is invalid


        //~ // Apply the *full* set of rules (prefix = num_rules)
        //~ const int prefix = recipes[recipe].num_rules;
        //~ REAL selected_volume = 0.0;
        //~ REAL integrals[DIAGS_INTEGRATION_MAX_INTEGRANDS] = {0.0};

        //~ diagnostics_integration_apply_rules(commondata, griddata, diagnostic_gfs_major, recipes[recipe].rules, prefix, recipes[recipe].integrands,
                                            //~ recipes[recipe].num_integrands, &selected_volume, integrals);

      //~ } // END LOOP over recipes
      //~ break;
    //~ }

    //~ case DIAGNOSTICS_VOLUME_WRITE: {

      //~ for (int recipe = 0; recipe < num_recipes; recipe++) {
        //~ if (recipes[recipe].name == NULL || recipes[recipe].name[0] == '\0') {
          //~ fprintf(stderr, "diags_integration_execute_recipes: recipe %d is missing a name, skipping.\n", recipe);
          //~ continue;
        //~ } // END IF recipe name is invalid
        //~ if (recipes[recipe].num_integrands <= 0 || recipes[recipe].num_integrands > DIAGS_INTEGRATION_MAX_INTEGRANDS) {
          //~ fprintf(stderr, "diags_integration_execute_recipes: recipe '%s' has an invalid number of integrands (%d), skipping.\n", recipes[recipe].name,
                  //~ recipes[recipe].num_integrands);
          //~ continue;
        //~ } // END IF integrand count is invalid


        //~ // Open a unique output file for this recipe
        //~ char filename[256];
        //~ snprintf(filename, sizeof(filename), "out3d-integrals-conv_factor%.2f-%s.txt", commondata->convergence_factor, recipes[recipe].name);
        //~ FILE *file_ptr = diags_integration_open_output_file(commondata->nn, filename);

        //~ // Write header once at the beginning of the simulation
        //~ if (commondata->nn == 0) {
          //~ diags_integration_write_header(file_ptr, recipe, &recipes[recipe]);
        //~ } // END IF first timestep



        //~ // Write the data row to the file

        //~ REAL selected_volume = result->proper_volume;
        //~ diags_integration_write_row(file_ptr, (double)commondata->time, (double)selected_volume, recipes[recipe].integrands, integrals,
                                    //~ recipes[recipe].num_integrands);

        //~ // Populate results structure if requested
        //~ if (results_out) {
          //~ // Bounds check
          //~ if (results_out->num_recipe_results >= DIAGS_INTEGRATION_MAX_RECIPES) {
            //~ fprintf(stderr,
                    //~ "Warning: Maximum number of recipe results (%d) exceeded. "
                    //~ "Recipe '%s' results will not be stored.\n",
                    //~ DIAGS_INTEGRATION_MAX_RECIPES, recipes[recipe].name);
          //~ } else {
            //~ diags_integration_recipe_result_t *result = &results_out->recipe_results[results_out->num_recipe_results];

            //~ // Store recipe-level data
            //~ result->recipe_name = recipes[recipe].name;
            //~ result->proper_volume = selected_volume;
            //~ result->num_integrand_results = recipes[recipe].num_integrands;

            //~ // Store integrand-level data
            //~ for (int spec = 0; spec < recipes[recipe].num_integrands; spec++) {
              //~ diags_integration_integrand_result_t *int_result = &result->integrand_results[spec];

              //~ const diags_integration_integrand_spec_t *spec_ptr = &recipes[recipe].integrands[spec];

              //~ // Store metadata
              //~ int_result->gf_index = spec_ptr->gf_index;
              //~ int_result->is_squared = spec_ptr->is_squared;

              //~ // Store name (pointer into diagnostic_gf_names or literal)
              //~ if (spec_ptr->gf_index == diags_integration_UNIT_INTEGRAND_GFINDEX) {
                //~ int_result->gf_name = "one";
              //~ } else {
                //~ int_result->gf_name = diagnostic_gf_names[spec_ptr->gf_index];
              //~ } // END IF unit integrand gf_index

              //~ // Store integral
              //~ int_result->integral = integrals[spec];

              //~ // Compute and store L2/RMS if applicable
              //~ if (spec_ptr->is_squared) {
                //~ const double integral_value = (double)integrals[spec];
                //~ int_result->L2 = sqrt(fmax(integral_value, 0.0));
                //~ int_result->RMS = (selected_volume > 0.0) ? sqrt(fmax(integral_value, 0.0) / (double)selected_volume) : 0.0;
                //~ int_result->L2_valid = 1;
                //~ int_result->RMS_valid = 1;
              //~ } else {
                //~ int_result->L2 = 0.0;
                //~ int_result->RMS = 0.0;
                //~ int_result->L2_valid = 0;
                //~ int_result->RMS_valid = 0;
              //~ } // END IF integrand is squared
            //~ } // END LOOP over integrands

            //~ results_out->num_recipe_results++;
          //~ } // END ELSE for bounds check
        //~ } // END IF results_out
        //~ fclose(file_ptr);
      //~ } // END LOOP over recipes
      //~ break;
    //~ }
  //~ }
//~ } // END FUNCTION diags_integration_execute_recipes

#endif // DIAGNOSTICS_VOLUME_INTEGRATION_HELPERS_H
