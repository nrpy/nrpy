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

#ifndef DIAGNOSTICS_VOLUME_INTEGRATION_HELPERS_DEFINES_H
#define DIAGNOSTICS_VOLUME_INTEGRATION_HELPERS_DEFINES_H


#ifndef REAL
#define REAL double // for the REALs below
#endif



#include "diagnostic_gfs.h" // diagnostic_gf_names[] and enum

// ======================== Limits (stack-only) ========================
// Upper bounds for small, fixed-size containers used in this module.
enum {
  DIAGS_INTEGRATION_MAX_RECIPES = 8,    // Maximum number of recipes processed at once.
  DIAGS_INTEGRATION_MAX_RULES = 4,      // Maximum number of spherical rules per recipe.
  DIAGS_INTEGRATION_MAX_INTEGRANDS = 16 // Maximum number of integrands per recipe.
};

// Special GF index indicating the unit integrand (proper volume element).
#define diags_integration_UNIT_INTEGRAND_GFINDEX (-1)

// ======================== Rule / Recipe types ========================

// Per-integrand specification.
typedef struct {
  int gf_index;   // Diagnostic GF enum value, or diags_integration_UNIT_INTEGRAND_GFINDEX (-1) for the unit integrand.
  int is_squared; // If nonzero, integrate f^2 and also print L2 and RMS.
} diags_integration_integrand_spec_t;

// Spherical include/exclude rule in Cartesian space.
typedef struct {
  REAL center_xyz[3]; // Sphere center (x,y,z) in Cartesian coordinates.
  REAL radius;        // Sphere radius (>= 0).
  int exclude_inside; // If 1, exclude r <= radius (inside); if 0, exclude r > radius (outside).
} diags_integration_sphere_rule_t;

// Describes an integration domain and the integrands to measure.
typedef struct {
  // Rules that define the integration domain
  int num_rules;                                                      // Number of active rules (0..DIAGS_INTEGRATION_MAX_RULES).
  diags_integration_sphere_rule_t rules[DIAGS_INTEGRATION_MAX_RULES]; // Rule array (only first num_rules are used).

  // Integrands to compute over that domain
  int num_integrands;                                                              // Number of integrands (0..DIAGS_INTEGRATION_MAX_INTEGRANDS).
  diags_integration_integrand_spec_t integrands[DIAGS_INTEGRATION_MAX_INTEGRANDS]; // Integrand specifications.

  const char *name; // Mandatory recipe name; used to construct the output filename.
} diags_integration_recipe_t;

// ======================== Result query types ========================

// Specifies the type of value to extract from integration results.
typedef enum {
  DIAGS_EXTRACT_VOLUME,   // Remaining proper volume after rule application
  DIAGS_EXTRACT_INTEGRAL, // Raw integral: integral of f dV or integral of f^2 dV
  DIAGS_EXTRACT_L2,       // L2 norm: sqrt(integral of f^2 dV), requires is_squared=1
  DIAGS_EXTRACT_RMS       // RMS: sqrt(integral of f^2 dV / V), requires is_squared=1
} diags_integration_extraction_type_t;

// Results for a single integrand within a recipe.
typedef struct {
  int gf_index;        // Diagnostic GF enum value (or diags_integration_UNIT_INTEGRAND_GFINDEX)
  const char *gf_name; // Human-readable name (pointer to diagnostic_gf_names[] or "one")
  int is_squared;      // Whether this integrand was squared (0 or 1)
  REAL integral;       // integral of f dV if not squared, integral of f^2 dV if squared
  REAL L2;             // sqrt(max(integral of f^2 dV, 0)); only valid if is_squared=1
  REAL RMS;            // sqrt(max(integral of f^2 dV / V, 0)); only valid if is_squared=1
  int L2_valid;        // 1 if L2 value is meaningful, 0 otherwise
  int RMS_valid;       // 1 if RMS value is meaningful, 0 otherwise
} diags_integration_integrand_result_t;

// Results for one recipe execution.
typedef struct {
  const char *recipe_name;   // Name of the recipe (pointer to user's string)
  REAL proper_volume;        // Remaining volume after applying all rules
  int num_integrand_results; // Number of valid integrand results
  diags_integration_integrand_result_t integrand_results[DIAGS_INTEGRATION_MAX_INTEGRANDS];
} diags_integration_recipe_result_t;

// Container for all recipe results from one execution.
typedef struct {
  int num_recipe_results; // Number of recipes executed
  diags_integration_recipe_result_t recipe_results[DIAGS_INTEGRATION_MAX_RECIPES];
} diags_integration_results_t;


#endif // DIAGNOSTICS_VOLUME_INTEGRATION_HELPERS_H
