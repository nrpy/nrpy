#ifndef __BHAH_DEFINES_H__
#define __BHAH_DEFINES_H__
// BHaH core header file, automatically generated from output_BHaH_defines_h within BHaH_defines_h.py,
//    DO NOT EDIT THIS FILE BY HAND.

// ----------------------------
// Basic definitions for module
// general:
// ----------------------------
#include <ctype.h>   // Character type functions, such as isdigit, isalpha, etc.
#include <errno.h>   // Error number definitions
#include <math.h>    // Transcendental functions, etc.
#include <stdbool.h> // bool-typed variables
#include <stdint.h>  // int8_t-typed variables
#include <stdio.h>   // Basic input/output functions, such as *printf, fopen, fwrite, etc.
#include <stdlib.h>  // malloc/free, etc.
#include <string.h>  // String handling functions, such as strlen, strcmp, etc.
#include <time.h>    // Time-related functions and types, such as time(), clock(),
// output_BHaH_defines_h(...,enable_intrinsics=True) was called so we intrinsics headers:
#include "intrinsics/simd_intrinsics.h"
#define REAL double
#define DOUBLE double

// These macros for MIN(), MAX(), and SQR() ensure that if the arguments inside
//   are a function/complex expression, the function/expression is evaluated
//   *only once* per argument. See https://lwn.net/Articles/983965/ for details.
// They are improvements over the original implementations:
// #define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
// #define MAX(A, B) ( ((A) > (B)) ? (A) : (B) )
// #define SQR(A) ((A) * (A))
#define MIN(A, B)                                                                                                                                    \
  ({                                                                                                                                                 \
    __typeof__(A) _a = (A);                                                                                                                          \
    __typeof__(B) _b = (B);                                                                                                                          \
    _a < _b ? _a : _b;                                                                                                                               \
  })
#define MAX(A, B)                                                                                                                                    \
  ({                                                                                                                                                 \
    __typeof__(A) _a = (A);                                                                                                                          \
    __typeof__(B) _b = (B);                                                                                                                          \
    _a > _b ? _a : _b;                                                                                                                               \
  })
#define SQR(A)                                                                                                                                       \
  ({                                                                                                                                                 \
    __typeof__(A) _a = (A);                                                                                                                          \
    _a *_a;                                                                                                                                          \
  })
#ifndef MAYBE_UNUSED
#if __cplusplus >= 201703L
#define MAYBE_UNUSED [[maybe_unused]]
#elif defined(__GNUC__) || defined(__clang__) || defined(__NVCC__)
#define MAYBE_UNUSED __attribute__((unused))
#else
#define MAYBE_UNUSED
#endif // END check for GCC, Clang, or NVCC
#endif // END MAYBE_UNUSED
// START: CodeParameters declared as #define.
#ifndef MAXNUMGRIDS
#define MAXNUMGRIDS 15 // nrpy.grid
#endif
// END: CodeParameters declared as #define.

// ----------------------------
// Basic definitions for module
// commondata_struct:
// ----------------------------
typedef struct __commondata_struct__ {
  REAL CFL_FACTOR; // (nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all)
  REAL dt;         // (nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all)
  REAL t_0;        // (nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all)
  REAL t_final;    // (nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all)
  REAL time;       // (nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all)
  int NUMGRIDS;    // (nrpy.grid)
  int nn;          // (nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all)
  int nn_0;        // (nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all)
} commondata_struct;

// ----------------------------
// Basic definitions for module
// params_struct:
// ----------------------------
typedef struct __params_struct__ {
  REAL Cart_originx; // (nrpy.grid)
  REAL Cart_originy; // (nrpy.grid)
  REAL Cart_originz; // (nrpy.grid)
  bool grid_rotates; // (nrpy.grid)
} params_struct;

// ----------------------------
// Basic definitions for module
// finite_difference:
// ----------------------------

// Set the number of ghost zones
// Note that upwinding in e.g., BSSN requires that NGHOSTS = fd_order/2 + 1 <- Notice the +1.
#define NGHOSTS 2

// Declare NO_INLINE macro, used in FD functions. GCC v10+ compilations hang on complex RHS expressions (like BSSN) without this.
#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define NO_INLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define NO_INLINE __declspec(noinline)
#else
#define NO_INLINE // Fallback for unknown compilers
#endif

// ----------------------------
// Basic definitions for module
// nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all:
// ----------------------------
typedef struct __MoL_gridfunctions_struct__ {
  REAL *restrict y_n_gfs;
  REAL *restrict y_nplus1_running_total_gfs;
  REAL *restrict k_odd_gfs;
  REAL *restrict k_even_gfs;
  REAL *restrict auxevol_gfs;
  REAL *restrict diagnostic_output_gfs;
  REAL *restrict diagnostic_output_gfs2;
} MoL_gridfunctions_struct;

// ----------------------------
// Basic definitions for module
// grid:
// ----------------------------

// EVOL VARIABLES:
#define NUM_EVOL_GFS 0

// AUX VARIABLES:
#define NUM_AUX_GFS 0

// AUXEVOL VARIABLES:
#define NUM_AUXEVOL_GFS 0

// ----------------------------
// Indexing macros
// ----------------------------
// IDX4: Converts 4D grid indices (gf, i, j, k) into a 1D array index using the strides
//       Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, and Nxx_plus_2NGHOSTS2. This macro assumes
//       that the "i" index varies fastest in memory.
#define IDX4(gf, i, j, k) ((i) + Nxx_plus_2NGHOSTS0 * ((j) + Nxx_plus_2NGHOSTS1 * ((k) + Nxx_plus_2NGHOSTS2 * (gf))))
// IDX4P: Similar to IDX4, but retrieves grid dimensions from the provided parameter structure
//        "params" instead of using global variables.
#define IDX4P(params, gf, i, j, k)                                                                                                                   \
  ((i) + (params)->Nxx_plus_2NGHOSTS0 * ((j) + (params)->Nxx_plus_2NGHOSTS1 * ((k) + (params)->Nxx_plus_2NGHOSTS2 * (gf))))
// IDX4pt: Computes the 1D index offset for a given grid function index (gf) based on an existing index (idx)
//         by using the total number of elements in one grid function, defined as the product of the grid strides.
#define IDX4pt(gf, idx) ((idx) + (Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2) * (gf))
// IDX4ptP: Similar to IDX4pt, but retrieves grid dimensions from the provided parameter structure
//        "params" instead of using global variables.
#define IDX4Ppt(params, gf, idx) ((idx) + ((params)->Nxx_plus_2NGHOSTS0 * (params)->Nxx_plus_2NGHOSTS1 * (params)->Nxx_plus_2NGHOSTS2) * (gf))
// IDX3: Converts 3D grid indices (i, j, k) into a 1D array index using the strides Nxx_plus_2NGHOSTS0
//       and Nxx_plus_2NGHOSTS1. Like IDX4, this macro assumes the "i" index varies fastest.
#define IDX3(i, j, k) ((i) + Nxx_plus_2NGHOSTS0 * ((j) + Nxx_plus_2NGHOSTS1 * ((k))))
// IDX3P: Similar to IDX3, but retrieves grid dimensions from the provided parameter structure "params".
#define IDX3P(params, i, j, k) ((i) + (params)->Nxx_plus_2NGHOSTS0 * ((j) + (params)->Nxx_plus_2NGHOSTS1 * ((k))))
// END: Indexing macros

// ----------------------------
// Loop-related macros
// ----------------------------
// SET_NXX_PLUS_2NGHOSTS_VARS: Declares local constants for the grid dimensions (including ghost zones) by extracting
// the values from griddata[whichgrid].params.
#define SET_NXX_PLUS_2NGHOSTS_VARS(whichgrid)                                                                                                        \
  const int Nxx_plus_2NGHOSTS0 = griddata[whichgrid].params.Nxx_plus_2NGHOSTS0;                                                                      \
  const int Nxx_plus_2NGHOSTS1 = griddata[whichgrid].params.Nxx_plus_2NGHOSTS1;                                                                      \
  const int Nxx_plus_2NGHOSTS2 = griddata[whichgrid].params.Nxx_plus_2NGHOSTS2;
// LOOP_REGION: Iterates over a 3D region defined by the inclusive lower bounds (i0min, i1min, i2min)
// and exclusive upper bounds (i0max, i1max, i2max) for each dimension.
#define LOOP_REGION(i0min, i0max, i1min, i1max, i2min, i2max)                                                                                        \
  for (int i2 = i2min; i2 < i2max; i2++)                                                                                                             \
    for (int i1 = i1min; i1 < i1max; i1++)                                                                                                           \
      for (int i0 = i0min; i0 < i0max; i0++)
// LOOP_OMP: Similar to LOOP_REGION but inserts an OpenMP pragma (via __OMP_PRAGMA__) for parallelization.
#define LOOP_OMP(__OMP_PRAGMA__, i0, i0min, i0max, i1, i1min, i1max, i2, i2min, i2max)                                                               \
  _Pragma(__OMP_PRAGMA__) for (int(i2) = (i2min); (i2) < (i2max); (i2)++) for (int(i1) = (i1min); (i1) < (i1max);                                    \
                                                                               (i1)++) for (int(i0) = (i0min); (i0) < (i0max); (i0)++)
// LOOP_NOOMP: A non-parallel version of the 3D loop, identical in structure to LOOP_REGION.
#define LOOP_NOOMP(i0, i0min, i0max, i1, i1min, i1max, i2, i2min, i2max)                                                                             \
  for (int(i2) = (i2min); (i2) < (i2max); (i2)++)                                                                                                    \
    for (int(i1) = (i1min); (i1) < (i1max); (i1)++)                                                                                                  \
      for (int(i0) = (i0min); (i0) < (i0max); (i0)++)
// LOOP_BREAKOUT: Forces an exit from the nested loops by setting the loop indices to their maximum values and executing a break.
#define LOOP_BREAKOUT(i0, i1, i2, i0max, i1max, i2max)                                                                                               \
  {                                                                                                                                                  \
    i0 = (i0max);                                                                                                                                    \
    i1 = (i1max);                                                                                                                                    \
    i2 = (i2max);                                                                                                                                    \
    break;                                                                                                                                           \
  }
// IS_IN_GRID_INTERIOR: Checks whether the provided 3D index array (i0i1i2) lies within the grid interior,
// defined as the region excluding NG ghost cells on each boundary.
#define IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, NG)                                                  \
  (i0i1i2[0] >= (NG) && i0i1i2[0] < (Nxx_plus_2NGHOSTS0) - (NG) && i0i1i2[1] >= (NG) && i0i1i2[1] < (Nxx_plus_2NGHOSTS1) - (NG) &&                   \
   i0i1i2[2] >= (NG) && i0i1i2[2] < (Nxx_plus_2NGHOSTS2) - (NG))

// ----------------------------
// Define griddata struct
// ----------------------------
typedef struct __griddata__ {
  // griddata_struct stores data needed on each grid
  // xx[3] stores the uniform grid coordinates.
  REAL *restrict xx[3];
  // NRPy+ MODULE: nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all
  MoL_gridfunctions_struct gridfuncs; // <- MoL gridfunctions
  // NRPy+ MODULE: params
  params_struct params; // <- BHaH parameters, generated from NRPy+'s CodeParameters
} griddata_struct;

#define BHAH_FREE(a)                                                                                                                                 \
  do {                                                                                                                                               \
    {                                                                                                                                                \
      if (a) {                                                                                                                                       \
        {                                                                                                                                            \
          free((void *)(a));                                                                                                                         \
          (a) = NULL;                                                                                                                                \
        }                                                                                                                                            \
      }                                                                                                                                              \
    }                                                                                                                                                \
  } while (0);
#endif

#define BHAH_MALLOC___PtrMember(a, b, sz)                                                                                                            \
  do {                                                                                                                                               \
    if (a) {                                                                                                                                         \
      a->b = malloc(sz);                                                                                                                             \
    }                                                                                                                                                \
  } while (0);
#define BHAH_FREE___PtrMember(a, b)                                                                                                                  \
  do {                                                                                                                                               \
    if (a) {                                                                                                                                         \
      BHAH_FREE(a->b);                                                                                                                               \
    }                                                                                                                                                \
  } while (0);
