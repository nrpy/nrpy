// Shared geometry for generated block kernels and host adapters.
#ifndef NRPY_DENDRO_BLOCK_GEOMETRY_H
#define NRPY_DENDRO_BLOCK_GEOMETRY_H
#include <cstddef>
#include <cstdint>
// The generated scalar contract header (<solver_stem>_types.h) owns the
// DendroScalar definition under the DENDRO_SCALAR_DEFINED guard, and this
// header respects the same guard so a translation unit can include both without
// redefining the alias.  In the standalone-host build the alias is always
// double.
#if !defined(DENDRO_SCALAR_DEFINED) && !defined(DendroScalar)
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
}; // END STRUCT: BlockGeometry

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
}; // END STRUCT: StandaloneHostMesh

#endif // NRPY_DENDRO_BLOCK_GEOMETRY_H
