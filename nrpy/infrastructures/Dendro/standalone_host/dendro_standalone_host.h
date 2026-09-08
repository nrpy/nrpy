// Standalone-host Dendro-GR host API for the NRPy-generated Dendro solver.
//
// This header provides the minimal host types needed to compile the
// NRPy-generated Dendro CFunctions in a standalone-host translation unit.  It
// declares just enough of the pinned Dendrolib API for the generated solver to
// compile and run; it does not implement Dendro.
#ifndef DENDRO_STANDALONE_HOST_H
#define DENDRO_STANDALONE_HOST_H

#include <cstddef>
#include <cstdint>

#include "block_geometry.h"

// --- Host lifecycle declarations (standalone-host build only) ---------------
// These declarations support the standalone test vehicle only. Real builds
// include the pinned Dendrolib headers and use the separate real context branch.
// The generated CMake records the Dendrolib pin; the validation page records
// which host capabilities have been exercised.

namespace standalone_host {

// Unzipped EVOL vector: one padded variable-major, x-fastest block array per
// component, in generated NRPy order.  In the standalone-host build unzip and
// zip are identity copies.
struct DVector {
  DendroScalar** comp;  // [num_components] -> [block][vol]
  unsigned num_blocks;
  unsigned num_components;
};  // END STRUCT: DVector

// Minimal timestep-context stub: the generated context owns the vectors and
// the local block list; the host integrator (here: one-stage Euler in the
// generated context) advances state.  LTS and real Dendro integration are
// out of scope for the standalone-host build.
struct Ctx {
  DVector in;
  DVector rhs;
  DVector out;
  // Diagnostic vector: recomputed from the evolved state,
  // never checkpoint state, so it is a separate vector from `in`/`out`.
  DVector diag;
  StandaloneHostMesh mesh;
};  // END STRUCT: Ctx

// clang-format off
}  // END NAMESPACE: standalone_host
// clang-format on

#endif  // DENDRO_STANDALONE_HOST_H
