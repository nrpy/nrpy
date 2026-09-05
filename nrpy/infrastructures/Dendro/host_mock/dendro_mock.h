// Mock Dendro-GR host API for NRPy Dendro generated-code compile tests.
//
// This header provides the minimal host types needed to compile the
// NRPy-generated Dendro CFunctions in a mock translation unit. It is a test
// double for the pinned Dendrolib API; it does not implement Dendro.
#ifndef DENDRO_MOCK_H
#define DENDRO_MOCK_H

#include <cstddef>
#include <cstdint>

// The generated scalar contract header (<solver_stem>_types.h) owns the
// DendroScalar definition under the DENDRO_SCALAR_DEFINED guard; the mock
// host respects the same guard so a translation unit can include both
// without redefining the alias.  In the mock vehicle the alias is always
// double.
#ifndef DENDRO_SCALAR_DEFINED
#define DENDRO_SCALAR_DEFINED
using DendroScalar = double;
#endif

// Normalized block geometry: interior extents, padding,
// the padded origin, the per-axis spacing, the per-component base offset,
// and the boundary flags.
struct BlockGeometry {
  unsigned nx;
  unsigned ny;
  unsigned nz;
  unsigned padding;
  std::size_t component_offset;
  DendroScalar pmin_padded[3];
  DendroScalar dx[3];
  std::uint32_t boundary_flags;
};  // END STRUCT: BlockGeometry

namespace mock {
// The mock world holds a small fixed number of blocks; the bound is named so
// callers can check `num_blocks` against it rather than trusting an unchecked
// array size (the real ot::Block list is dynamic and arrives with the pinned
// Dendrolib gates).
inline constexpr unsigned MAX_MOCK_BLOCKS = 2;
}  // END NAMESPACE: mock

// A tiny world of a few blocks so the mock can exercise the NRPy block loop.
struct MockWorld {
  BlockGeometry geom[mock::MAX_MOCK_BLOCKS];
  unsigned num_blocks;
};  // END STRUCT: MockWorld

// --- Mock-vehicle host lifecycle stubs (dendro_mock only) ---------------
// Dendrolib is UNPINNED in this repository (dendrolib_pin.json); the real
// ot::DVector / ts::Ctx APIs are test doubles below.  They exist only so the
// generated solver context can be built and the Minkowski lifecycle run in
// the mock vehicle; the real-host signatures stay frozen until the I0-1
// capability proof lands; see dendrolib_pin.json and
// dendrolib_capabilities.json.

namespace mock {

// Unzipped EVOL vector: one padded variable-major, x-fastest block array per
// component, in generated NRPy order (mock unzip/zip are identity copies).
struct DVector {
  DendroScalar** comp;  // [num_components] -> [block][vol]
  unsigned num_blocks;
  unsigned num_components;
};  // END STRUCT: DVector

// Minimal timestep-context stub: the generated context owns the vectors and
// the local block list; the host integrator (here: one-stage Euler in the
// generated context) advances state.  LTS and real Dendro integration are
// out of scope for the mock vehicle.
struct Ctx {
  DVector in;
  DVector rhs;
  DVector out;
  // Diagnostic vector: recomputed from the evolved state,
  // never checkpoint state, so it is a separate vector from `in`/`out`.
  DVector diag;
  MockWorld world;
};  // END STRUCT: Ctx

}  // END NAMESPACE: mock

#endif  // DENDRO_MOCK_H
