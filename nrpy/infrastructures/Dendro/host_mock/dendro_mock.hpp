// Mock Dendro-GR host API for NRPy Dendro generated-code compile tests.
//
// This header provides the minimal host types needed to compile the
// NRPy-generated Dendro CFunctions in a mock translation unit. It is a test
// double for the pinned Dendrolib API; it does not implement Dendro.
#pragma once

#include <cstddef>
#include <cstdint>

// The generated scalar contract (fccz4_types.hpp) owns the DendroScalar
// definition under the DENDRO_SCALAR_DEFINED guard; the mock host respects
// the same guard so a translation unit can include both without redefining
// the alias.  In the mock vehicle the alias is always double.
#ifndef DENDRO_SCALAR_DEFINED
#define DENDRO_SCALAR_DEFINED
using DendroScalar = double;
#endif

// Normalized block geometry (whitepaper 7.4): interior extents, padding,
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
};

namespace mock {
// The mock world holds a small fixed number of blocks; the bound is named so
// callers can check `num_blocks` against it rather than trusting an unchecked
// array size (the real ot::Block list is dynamic and arrives with the pinned
// Dendrolib gates).
inline constexpr unsigned MAX_MOCK_BLOCKS = 2;
}  // namespace mock

// A tiny world of a few blocks so the mock can exercise the NRPy block loop.
struct MockWorld {
    BlockGeometry geom[mock::MAX_MOCK_BLOCKS];
    unsigned num_blocks;
};

// --- PR 7 mock-vehicle host lifecycle stubs (dendro_mock only) -----------
// Dendrolib is UNPINNED in this repository (dendrolib_pin.json); the real
// ot::DVector / ts::Ctx APIs are test doubles below.  They exist only so the
// generated FCCZ4GR context can be built and the Minkowski lifecycle run in
// the mock vehicle; the real-host signatures stay frozen until the I0-1
// capability proof lands (whitepaper sections 13.2/13.4/13.5).

namespace mock {

// Unzipped EVOL vector: one padded variable-major, x-fastest block array per
// component, in generated NRPy order (mock unzip/zip are identity copies).
struct DVector {
    DendroScalar** comp;  // [num_components] -> [block][vol]
    unsigned num_blocks;
    unsigned num_components;
};

// Minimal timestep-context stub: the generated context owns the vectors and
// the local block list; the host integrator (here: one-stage Euler in the
// generated context) advances state.  LTS and real Dendro integration are
// out of scope for the mock vehicle.
struct Ctx {
    DVector in;
    DVector rhs;
    DVector out;
    // Diagnostic vector (whitepaper 14.6): recomputed from the evolved state,
    // never checkpoint state, so it is a separate vector from `in`/`out`.
    DVector diag;
    MockWorld world;
};

}  // namespace mock