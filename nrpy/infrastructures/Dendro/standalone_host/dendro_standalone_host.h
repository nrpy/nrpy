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

// The generated scalar contract header (<solver_stem>_types.h) owns the
// DendroScalar definition under the DENDRO_SCALAR_DEFINED guard, and this
// header respects the same guard so a translation unit can include both without
// redefining the alias.  In the standalone-host build the alias is always
// double.
#ifndef DENDRO_SCALAR_DEFINED
#define DENDRO_SCALAR_DEFINED
using DendroScalar = double;
#endif

// Normalized block geometry: the padded per-axis extents (nx, ny and nz each
// count the ghost points on both sides, so they size the allocation and set
// the strides; a point loop subtracts the padding to reach the interior), the
// padding, the padded origin, the per-axis spacing, the per-component base
// offset, and the boundary flags.
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

namespace standalone_host {
// The standalone-host mesh holds a small fixed number of blocks, and the bound
// is named so callers can check `num_blocks` against it rather than trusting an
// unchecked array size.  The real ot::Block list is dynamic.
inline constexpr unsigned MAX_STANDALONE_HOST_BLOCKS = 2;
// The padded extent per axis is bounded for the same reason: the standalone
// host allocates one dense array per component per block -- three evolved
// vectors, the diagnostic vector, and the drift snapshot -- so an absurd extent
// is refused with a message instead of reaching std::bad_alloc.  512 padded
// points per axis is 1.1 GB for one component of one block at double precision,
// so the bound refuses the absurd rather than promising that every accepted
// extent fits in memory.
inline constexpr unsigned MAX_STANDALONE_HOST_EXTENT = 512;
// clang-format off
}  // END NAMESPACE: standalone_host
// clang-format on

// A mesh of a few blocks, enough for the NRPy block loop to run.  The real
// ot::Mesh block list is dynamic and arrives with the pinned Dendrolib gates.
struct StandaloneHostMesh {
  BlockGeometry geom[standalone_host::MAX_STANDALONE_HOST_BLOCKS];
  unsigned num_blocks;
};  // END STRUCT: StandaloneHostMesh

// --- Host lifecycle declarations (standalone-host build only) ---------------
// The generated solver is not yet built against the real host, so the
// ot::DVector and ts::Ctx APIs are declared below only as far as the generated
// solver uses them.  They exist so the generated solver context can be built
// and the Minkowski lifecycle run in the standalone-host build; the real-host
// signatures stay frozen until the generated solver is built against a real
// Dendro-GR checkout.  The generated CMakeLists.txt records the Dendrolib
// commit these assumptions were proven against, and the Dendro validation page
// in the NRPy knowledge base carries the proven axes.

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
