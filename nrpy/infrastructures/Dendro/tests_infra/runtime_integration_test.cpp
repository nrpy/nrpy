// Compile against a generated real context, not standalone host declarations.
// RUNTIME_HEADER and RUNTIME_NAMESPACE select the generated formulation.
//
// Author: Zachariah B. Etienne
//         zachetie **at** gmail **dot* com
#include RUNTIME_HEADER
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

#include "meshUtils.h"
#include "octUtils.h"
namespace app = RUNTIME_NAMESPACE;

namespace allocation_measurement {
std::atomic<bool> enabled{false};
std::atomic<unsigned long long> count{0};
std::atomic<unsigned long long> bytes{0};

void record(std::size_t size) noexcept {
    if (enabled.load(std::memory_order_relaxed)) {
        count.fetch_add(1, std::memory_order_relaxed);
        bytes.fetch_add(size, std::memory_order_relaxed);
    }  // END IF: allocation measurement enabled
}  // END FUNCTION: record
// clang-format off
}  // END NAMESPACE: allocation measurement
// clang-format on

void *operator new(std::size_t size) {
    allocation_measurement::record(size);
    if (void *memory = std::malloc(size == 0 ? 1 : size)) return memory;
    throw std::bad_alloc();
}  // END FUNCTION: operator new

void *operator new[](std::size_t size) {
    allocation_measurement::record(size);
    if (void *memory = std::malloc(size == 0 ? 1 : size)) return memory;
    throw std::bad_alloc();
}  // END FUNCTION: operator new[]

void operator delete(void *memory) noexcept {
    std::free(memory);
}  // END FUNCTION: operator delete
void operator delete[](void *memory) noexcept {
    std::free(memory);
}  // END FUNCTION: operator delete[]
void operator delete(void *memory, std::size_t) noexcept {
    std::free(memory);
}  // END FUNCTION: operator delete
void operator delete[](void *memory, std::size_t) noexcept {
    std::free(memory);
}  // END FUNCTION: operator delete[]

void *operator new(std::size_t size, std::align_val_t alignment) {
    allocation_measurement::record(size);
    void *memory = nullptr;
    if (posix_memalign(&memory, static_cast<std::size_t>(alignment),
                       size == 0 ? 1 : size) == 0)
        return memory;
    throw std::bad_alloc();
}  // END FUNCTION: operator new

void *operator new[](std::size_t size, std::align_val_t alignment) {
    return operator new(size, alignment);
}  // END FUNCTION: operator new[]

void operator delete(void *memory, std::align_val_t) noexcept {
    std::free(memory);
}  // END FUNCTION: operator delete
void operator delete[](void *memory, std::align_val_t) noexcept {
    std::free(memory);
}  // END FUNCTION: operator delete[]
void operator delete(void *memory, std::size_t, std::align_val_t) noexcept {
    std::free(memory);
}  // END FUNCTION: operator delete
void operator delete[](void *memory, std::size_t,
                       std::align_val_t) noexcept {
    std::free(memory);
}  // END FUNCTION: operator delete[]

namespace {
double field(unsigned f, double x, double y, double z) {
    return 101.0 * (f + 1) + (f + 2) * x - (2 * f + 3) * y + (3 * f + 5) * z;
}  // END FUNCTION: field

/**
 * Compare generated callbacks on one local Dendro block.
 *
 * 1. Selects a local block with a nonzero component offset when available.
 * 2. Compares whole-vector, unzipped, and flat component-major RHS values.
 * 3. Compares flat algebraic projection with whole-vector projection.
 *
 * @param[in,out] context Generated context whose state and work vectors change.
 * @param[in] mesh Borrowed mesh that supplies local blocks and vector layouts.
 * @param[in] minimum Physical domain minimum used to construct block geometry.
 * @param[in] maximum Physical domain maximum used to construct block geometry.
 * @param[in] fault Optional invalid input selected by the negative tests.
 * @param rank MPI rank used for root-only diagnostics.
 *
 * @note The caller retains ownership of the context and mesh.
 * @warning Invalid input aborts MPI_COMM_WORLD; failed comparisons throw
 * std::runtime_error.
 */
void qualify_block_callbacks(app::Ctx &context, ot::Mesh &mesh,
                             const Point &minimum, const Point &maximum,
                             const std::string &fault, int rank) {
    const unsigned dof                 = app::generated::NUM_EVOL_GFS;
    constexpr double sentinel          = 9.87654321e200;
    unsigned long long callback_points = 0, callback_offsets = 0,
                       preserved_points = 0, projection_points = 0,
                       block_rhs_allocations = 0,
                       block_rhs_allocation_bytes = 0;
    double flat_rhs_error = 0.0, whole_rhs_error = 0.0, projection_error = 0.0,
           whole_rhs_scale = 0.0;
    context.initialize();
    if (mesh.isActive()) {
        // A constant nonzero shift-driver field gives d(B^0)/dt=-eta*B^0.
        // This prevents a zero Minkowski RHS from hiding component-routing or
        // input-selection errors in the whole-vector/blockwise comparison.
        context.params.eta = 1.0;
        const unsigned b0 =
            static_cast<unsigned>(app::generated::EvolVar::betU0);
        const unsigned state_stride = mesh.getDegOfFreedom();
        for (unsigned i = mesh.getNodeLocalBegin(); i < mesh.getNodeLocalEnd();
             ++i)
            context.state.get_vec_ptr()[std::size_t(b0) * state_stride + i] =
                0.01;
        const auto &blocks = mesh.getLocalBlockList();
        if (blocks.empty())
            throw std::runtime_error("active rank has no local blocks");
        unsigned selected = 0;
        for (unsigned id = 0; id < blocks.size(); ++id) {
            const auto candidate =
                app::block_geometry(mesh, blocks[id], minimum, maximum);
            if (candidate.component_offset != 0) {
                selected = id;
                break;
            }  // END IF: nonzero component offset found
        }  // END LOOP: for id seeking nonzero offset
        auto geometry =
            app::block_geometry(mesh, blocks[selected], minimum, maximum);
        if (geometry.component_offset != 0) ++callback_offsets;
        double block_time        = 0.0;
        const std::size_t stride = mesh.getDegOfFreedomUnZip();
        const std::size_t volume =
            std::size_t(geometry.nx) * geometry.ny * geometry.nz;

        // Retain the selected block result from the whole-vector callback.
        // rhs() owns exchange, exterior values, all-block traversal, and zip;
        // rhs_blkwise() below starts from the unzipped input it prepared.
        app::DVec whole_output;
        whole_output.create_vector(&mesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                                   ot::DVEC_LOC::HOST, dof, true);
        context.rhs(&context.state, &whole_output, 1, block_time);
        std::vector<double *> whole_pointers(dof);
        context.unzipped_rhs.to_2d(whole_pointers.data());
        std::vector<double> whole_block(std::size_t(dof) * volume);
        for (unsigned f = 0; f < dof; ++f)
            for (std::size_t cell = 0; cell < volume; ++cell)
                whole_block[std::size_t(f) * volume + cell] =
                    whole_pointers[f][geometry.component_offset + cell];
        whole_output.destroy_vector();
        std::fill_n(context.unzipped_rhs.get_vec_ptr(),
                    context.unzipped_rhs.get_size(), sentinel);
        if (fault == "blockwise_null_ids")
            context.rhs_blkwise(context.unzipped, context.unzipped_rhs, nullptr,
                                1, &block_time);
        if (fault == "blockwise_bad_id") {
            const unsigned bad_id = blocks.size();
            context.rhs_blkwise(context.unzipped, context.unzipped_rhs, &bad_id,
                                1, &block_time);
        }  // END IF: inject out-of-range block id
        if (fault == "blockwise_bad_dof") {
            app::DVec wrong_dof;
            wrong_dof.create_vector(&mesh,
                                    ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,
                                    ot::DVEC_LOC::HOST, dof - 1, true);
            context.rhs_blkwise(wrong_dof, context.unzipped_rhs, &selected, 1,
                                &block_time);
            wrong_dof.destroy_vector();
        }  // END IF: inject blockwise field-count fault
        const unsigned selected_ids[1] = {selected};
        context.rhs_blkwise(context.unzipped, context.unzipped_rhs,
                            selected_ids, 1, &block_time);

        // 1: selected block interior. 2: x padding of an interior row of the
        // selected block, which a SIMD kernel may write when the row has fewer
        // interior points than SIMD_WIDTH.
        std::vector<unsigned char> selected_interior(stride, 0);
        for (unsigned k = geometry.padding; k < geometry.nz - geometry.padding;
             ++k)
            for (unsigned j = geometry.padding;
                 j < geometry.ny - geometry.padding; ++j)
                for (unsigned i = 0; i < geometry.nx; ++i) {
                    const std::size_t cell =
                        geometry.component_offset + i +
                        std::size_t(geometry.nx) *
                            (j + std::size_t(geometry.ny) * k);
                    const bool interior =
                        i >= geometry.padding &&
                        i < geometry.nx - geometry.padding;
                    selected_interior[cell] = interior ? 1 : 2;
                    if (interior) ++callback_points;
                }  // END LOOP: for i over padded x
        std::vector<double *> input(dof), output(dof);
        context.unzipped.to_2d(input.data());
        context.unzipped_rhs.to_2d(output.data());
        for (unsigned f = 0; f < dof; ++f)
            for (std::size_t cell = 0; cell < stride; ++cell) {
                if (selected_interior[cell] == 2) continue;
                if (selected_interior[cell] == 1) {
                    whole_rhs_scale = std::max(
                        whole_rhs_scale,
                        std::abs(whole_block[std::size_t(f) * volume + cell -
                                             geometry.component_offset]));
                    if (!std::isfinite(output[f][cell]) ||
                        output[f][cell] == sentinel)
                        throw std::runtime_error(
                            "selected block interior was not evaluated");
                    whole_rhs_error = std::max(
                        whole_rhs_error,
                        std::abs(output[f][cell] -
                                 whole_block[std::size_t(f) * volume + cell -
                                             geometry.component_offset]));
                }  // END IF: selected block interior
                else {
                    ++preserved_points;
                    if (output[f][cell] != sentinel)
                        throw std::runtime_error(
                            "blockwise RHS wrote outside selected interior");
                }  // END ELSE: point outside selected block
            }  // END LOOP: for cell over unzipped points

        std::vector<double> flat_input(std::size_t(dof) * volume);
        std::vector<double> flat_output(std::size_t(dof) * volume, sentinel);
        for (unsigned f = 0; f < dof; ++f)
            for (std::size_t cell = 0; cell < volume; ++cell)
                flat_input[std::size_t(f) * volume + cell] =
                    input[f][geometry.component_offset + cell];
        if (fault == "block_null")
            context.rhs_blk(nullptr, flat_output.data(), dof, selected,
                            block_time);
        if (fault == "block_bad_dof")
            context.rhs_blk(flat_input.data(), flat_output.data(), dof - 1,
                            selected, block_time);
        if (fault == "block_bad_id")
            context.rhs_blk(flat_input.data(), flat_output.data(), dof,
                            blocks.size(), block_time);
        allocation_measurement::count.store(0, std::memory_order_relaxed);
        allocation_measurement::bytes.store(0, std::memory_order_relaxed);
        allocation_measurement::enabled.store(true, std::memory_order_relaxed);
        context.rhs_blk(flat_input.data(), flat_output.data(), dof, selected,
                        block_time);
        allocation_measurement::enabled.store(false, std::memory_order_relaxed);
        block_rhs_allocations =
            allocation_measurement::count.load(std::memory_order_relaxed);
        block_rhs_allocation_bytes =
            allocation_measurement::bytes.load(std::memory_order_relaxed);
        const unsigned long long mesh_storage_bytes =
            static_cast<unsigned long long>(flat_input.size() * sizeof(double));
        if (block_rhs_allocation_bytes >= mesh_storage_bytes)
            throw std::runtime_error("flat block RHS allocated mesh-sized storage");
        for (unsigned f = 0; f < dof; ++f)
            for (std::size_t cell = 0; cell < volume; ++cell)
                flat_rhs_error = std::max(
                    flat_rhs_error,
                    std::abs(flat_output[std::size_t(f) * volume + cell] -
                             output[f][geometry.component_offset + cell]));

        const std::vector<double> before_hooks = flat_input;
        context.pre_stage_blk(flat_input.data(), dof, selected, block_time);
        context.post_stage_blk(flat_input.data(), dof, selected, block_time);
        context.pre_timestep_blk(flat_input.data(), dof, selected, block_time);
        if (std::memcmp(flat_input.data(), before_hooks.data(),
                        flat_input.size() * sizeof(double)) != 0)
            throw std::runtime_error("no-op block hook changed its input");

        // Give both projection paths the same nontrivial positive-definite
        // conformal metric and nonzero conformal extrinsic curvature.
        context.initialize();
        const unsigned h00 =
            static_cast<unsigned>(app::generated::EvolVar::hDD00);
        const unsigned a00 =
            static_cast<unsigned>(app::generated::EvolVar::aDD00);
        const unsigned zipped_stride = mesh.getDegOfFreedom();
        for (unsigned i = mesh.getNodeLocalBegin(); i < mesh.getNodeLocalEnd();
             ++i) {
            context.state.get_vec_ptr()[std::size_t(h00) * zipped_stride + i] +=
                0.125;
            context.state.get_vec_ptr()[std::size_t(a00) * zipped_stride + i] +=
                0.03125;
        }  // END LOOP: for i over owned nodes
        context.unzip(context.state, context.unzipped, 1);
        context.unzipped.to_2d(input.data());
        std::vector<double> flat_projection(std::size_t(dof) * volume);
        for (unsigned f = 0; f < dof; ++f)
            for (std::size_t cell = 0; cell < volume; ++cell)
                flat_projection[std::size_t(f) * volume + cell] =
                    input[f][geometry.component_offset + cell];
        if (fault == "projection_null")
            context.post_timestep_blk(nullptr, dof, selected, block_time);
        if (fault == "projection_bad_dof")
            context.post_timestep_blk(flat_projection.data(), dof - 1, selected,
                                      block_time);
        if (fault == "projection_bad_id")
            context.post_timestep_blk(flat_projection.data(), dof,
                                      blocks.size(), block_time);
        context.post_timestep_blk(flat_projection.data(), dof, selected,
                                  block_time);
        context.post_timestep(context.state);
        context.unzip(context.state, context.unzipped, 1);
        context.unzipped.to_2d(input.data());
        for (unsigned f = 0; f < dof; ++f)
            for (std::size_t cell = 0; cell < volume; ++cell) {
                projection_error = std::max(
                    projection_error,
                    std::abs(flat_projection[std::size_t(f) * volume + cell] -
                             input[f][geometry.component_offset + cell]));
                ++projection_points;
            }  // END LOOP: for cell over selected block
    }  // END IF: qualify active-rank callbacks

    unsigned long long totals[6] = {},
                       local[6]  = {callback_points, callback_offsets,
                                    preserved_points, projection_points,
                                    block_rhs_allocations,
                                    block_rhs_allocation_bytes};
    MPI_Allreduce(local, totals, 6, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
    double local_errors[4] = {flat_rhs_error, whole_rhs_error, projection_error,
                              whole_rhs_scale},
           errors[4]       = {};
    MPI_Allreduce(local_errors, errors, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    const bool ok = totals[0] > 0 && totals[1] > 0 && totals[2] > 0 &&
                    totals[3] > 0 && errors[0] < 1e-12 && errors[1] < 1e-12 &&
                    errors[2] < 1e-12 && errors[3] > 1e-6;
    if (rank == 0)
        std::printf(
            "REAL_CALLBACKS %s interior=%llu nonzero_offsets=%llu "
            "preserved=%llu projected=%llu flat_rhs_error=%.17g "
            "whole_rhs_error=%.17g projection_error=%.17g "
            "whole_rhs_scale=%.17g block_rhs_allocations=%llu "
            "block_rhs_allocation_bytes=%llu\n",
            ok ? "PASS" : "FAIL", totals[0], totals[1], totals[2], totals[3],
            errors[0], errors[1], errors[2], errors[3], totals[4], totals[5]);
    if (!ok) throw std::runtime_error("block callback qualification failed");
}  // END FUNCTION: qualify_block_callbacks
// clang-format off
} // END NAMESPACE: independent field oracle
// clang-format on
/**
 * Qualify real block geometry, distributed transport, runtime parameters, and
 * failures.
 *
 * @param argc Number of command-line arguments.
 * @param[in,out] argv Argument vector parsed by MPI and this entry point.
 * @return 0 on success, 1 if MPI_Abort unexpectedly returns after a failure.
 *
 * @note Detected failures abort MPI_COMM_WORLD, including inactive ranks.
 */
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0, ranks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try {
        const std::string fault = argc > 1 ? argv[1] : "";
        m_uiMaxDepth            = 8;
        _InitializeHcurve(m_uiDim);
        std::function<double(double, double, double)> refine =
            [](double x, double y, double z) {
                return std::exp(-(x * x + y * y + z * z) / 0.5);
            };  // END LAMBDA: choose initial octree refinement
        std::vector<ot::TreeNode> octree;
        const unsigned order = app::generated::FD_ORDER;
        const unsigned mesh_order = fault == "element_order" ? order + 2 : order;
        function2Octree(refine, octree, 5, 1e-3, mesh_order, MPI_COMM_WORLD);
        std::unique_ptr<ot::Mesh> mesh(
            ot::createMesh(octree.data(), octree.size(), mesh_order,
                           MPI_COMM_WORLD, 0, ot::SM_TYPE::FDM));
        const double lo[3] = {-1, -2, -4}, hi[3] = {3, 2, 4};
        const Point minimum(lo[0], lo[1], lo[2]), maximum(hi[0], hi[1], hi[2]);
        mesh->setDomainBounds(minimum, maximum);
        {
            app::Ctx context(mesh.get(), minimum, maximum, 0.001);
            if (fault == "nonfinite") {
                context.initialize();
                if (rank == ranks - 1 && mesh->isActive())
                    context.state.get_vec_ptr()[mesh->getNodeLocalBegin()] =
                        std::numeric_limits<double>::quiet_NaN();
                context.max_rhs();
            }  // END IF: inject rank local NaN
            else {
                const unsigned dof = app::generated::NUM_EVOL_GFS;
                std::function<void(double, double, double, double *)> fill =
                    [dof](double x, double y, double z, double *out) {
                        for (unsigned f = 0; f < dof; ++f)
                            out[f] = field(f, x, y, z);
                    };  // END LAMBDA: initialize distinct nodal components
                std::unique_ptr<double[]> expected(
                    mesh->createCGVector<double>(fill, dof));
                unsigned long long halos = 0, remote = 0, blocks = 0,
                                   offset_blocks = 0;
                double error                     = 0;
                if (mesh->isActive()) {
                    blocks                = mesh->getLocalBlockList().size();
                    const unsigned stride = mesh->getDegOfFreedom();
                    std::copy_n(expected.get(), context.state.get_size(),
                                context.state.get_vec_ptr());
                    // Poison every non-owned zipped node. A missing MPI
                    // exchange cannot pass by reading zero-filled or
                    // analytically prefilled ghosts.
                    for (unsigned f = 0; f < dof; ++f)
                        for (unsigned i = 0; i < stride; ++i)
                            if (i < mesh->getNodeLocalBegin() ||
                                i >= mesh->getNodeLocalEnd())
                                context.state
                                    .get_vec_ptr()[std::size_t(f) * stride +
                                                   i] =
                                    std::numeric_limits<double>::quiet_NaN();
                    std::fill_n(context.unzipped.get_vec_ptr(),
                                context.unzipped.get_size(),
                                std::numeric_limits<double>::quiet_NaN());
                    context.unzip(context.state, context.unzipped, 1);
                    // Dendrolib may reserve unused non-owned slots. Only its
                    // receive scatter map names nodes promised by the MPI
                    // exchange.
                    for (const unsigned i : mesh->getRecvNodeSM()) {
                        for (unsigned f = 0; f < dof; ++f)
                            if (!std::isfinite(
                                    context.state
                                        .get_vec_ptr()[std::size_t(f) * stride +
                                                       i]))
                                throw std::runtime_error(
                                    "unfilled distributed zipped ghost");
                        ++remote;
                    }  // END LOOP: for i over received ghosts
                    std::vector<double *> pointers(dof);
                    context.unzipped.to_2d(pointers.data());
                    for (const auto &b : mesh->getLocalBlockList()) {
                        auto g =
                            app::block_geometry(*mesh, b, minimum, maximum);
                        if (g.component_offset) ++offset_blocks;
                        if (fault == "offset" && rank == ranks - 1)
                            g.component_offset = 0;
                        // Expected coordinates derive independently from the
                        // raw octree.
                        const auto node      = b.getBlockNode();
                        const unsigned bflag = b.getBlkNodeFlag();
                        const double base[3] = {double(node.minX()),
                                                double(node.minY()),
                                                double(node.minZ())};
                        const double grid_step =
                            double(1u
                                   << (m_uiMaxDepth - b.getRegularGridLev())) /
                            order;
                        for (unsigned k = 0; k < g.nz; ++k)
                            for (unsigned j = 0; j < g.ny; ++j)
                                for (unsigned i = 0; i < g.nx; ++i) {
                                    const unsigned index[3] = {i, j, k};
                                    double x[3];
                                    for (unsigned a = 0; a < 3; ++a) {
                                        x[a] =
                                            lo[a] +
                                            (hi[a] - lo[a]) *
                                                (base[a] + (double(index[a]) -
                                                            b.get1DPadWidth()) *
                                                               grid_step) /
                                                double(1u << m_uiMaxDepth);
                                    }  // END LOOP: for a over coordinate axes
                                    const bool halo = i < g.padding ||
                                                      j < g.padding ||
                                                      k < g.padding ||
                                                      i >= g.nx - g.padding ||
                                                      j >= g.ny - g.padding ||
                                                      k >= g.nz - g.padding;
                                    const bool exterior =
                                        ((bflag & (1u << OCT_DIR_LEFT)) &&
                                         i < g.padding) ||
                                        ((bflag & (1u << OCT_DIR_RIGHT)) &&
                                         i >= g.nx - g.padding) ||
                                        ((bflag & (1u << OCT_DIR_DOWN)) &&
                                         j < g.padding) ||
                                        ((bflag & (1u << OCT_DIR_UP)) &&
                                         j >= g.ny - g.padding) ||
                                        ((bflag & (1u << OCT_DIR_BACK)) &&
                                         k < g.padding) ||
                                        ((bflag & (1u << OCT_DIR_FRONT)) &&
                                         k >= g.nz - g.padding);
                                    if (exterior) continue;
                                    if (halo) ++halos;
                                    const std::size_t cell =
                                        g.component_offset + i +
                                        std::size_t(g.nx) *
                                            (j + std::size_t(g.ny) * k);
                                    for (unsigned f = 0; f < dof; ++f) {
                                        if (fault == "halo" &&
                                            rank == ranks - 1 && halo)
                                            pointers[f][cell] = 0;
                                        if (!std::isfinite(pointers[f][cell]))
                                            throw std::runtime_error(
                                                "unfilled block halo or "
                                                "interior");
                                        error = std::max(
                                            error, std::abs(pointers[f][cell] -
                                                            field(f, x[0], x[1],
                                                                  x[2])));
                                        // Nonconstant zip oracle: change values
                                        // in the block array, then compare
                                        // owned zipped nodes with independent
                                        // nodal data.
                                        pointers[f][cell] =
                                            2 * field(f, x[0], x[1], x[2]) + 7;
                                    }  // END LOOP: for f over field components
                                }  // END LOOP: for i over padded x
                    }  // END LOOP: for b over local blocks
                    context.zip(context.unzipped, context.state);
                    for (unsigned f = 0; f < dof; ++f)
                        for (unsigned i = mesh->getNodeLocalBegin();
                             i < mesh->getNodeLocalEnd(); ++i) {
                            const std::size_t cell =
                                std::size_t(f) * stride + i;
                            if (!std::isfinite(
                                    context.state.get_vec_ptr()[cell]))
                                throw std::runtime_error(
                                    "nonfinite zipped result");
                            error = std::max(
                                error,
                                std::abs(context.state.get_vec_ptr()[cell] -
                                         (2 * expected[cell] + 7)));
                        }  // END LOOP: for i over owned nodes
                }  // END IF: verify active rank transport
                unsigned long long total_halos = 0, total_remote = 0,
                                   max_blocks = 0, total_offsets = 0;
                MPI_Allreduce(&halos, &total_halos, 1, MPI_UNSIGNED_LONG_LONG,
                              MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&remote, &total_remote, 1, MPI_UNSIGNED_LONG_LONG,
                              MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&blocks, &max_blocks, 1, MPI_UNSIGNED_LONG_LONG,
                              MPI_MAX, MPI_COMM_WORLD);
                MPI_Allreduce(&offset_blocks, &total_offsets, 1,
                              MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
                double global_error = 0;
                MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX,
                              MPI_COMM_WORLD);
                int active = mesh->isActive(), active_ranks = 0;
                MPI_Allreduce(&active, &active_ranks, 1, MPI_INT, MPI_SUM,
                              MPI_COMM_WORLD);
                const bool distributed_ok =
                    ranks == 1 || (active_ranks >= 2 && total_remote > 0);
                const bool ok = active_ranks >= 1 && max_blocks >= 2 &&
                                total_offsets > 0 && total_halos > 0 &&
                                distributed_ok && global_error < 1e-9;
                if (rank == 0)
                    std::printf(
                        "REAL_TRANSPORT %s active_ranks=%d "
                        "max_local_blocks=%llu "
                        "offset_blocks=%llu halos=%llu remote_nodes=%llu "
                        "error=%.17g\n",
                        ok ? "PASS" : "FAIL", active_ranks, max_blocks,
                        total_offsets, total_halos, total_remote, global_error);
                if (!ok)
                    throw std::runtime_error("transport qualification failed");
                qualify_block_callbacks(context, *mesh, minimum, maximum, fault,
                                        rank);
                // Exercise the real RHS callback with a parameter-dependent
                // nonzero state. For the shipped shift gauge, d(B^0)/dt
                // contains -eta*B^0.
                context.initialize();
                const unsigned b0 =
                    static_cast<unsigned>(app::generated::EvolVar::betU0);
                const unsigned stride = mesh->getDegOfFreedom();
                if (mesh->isActive())
                    for (unsigned i = mesh->getNodeLocalBegin();
                         i < mesh->getNodeLocalEnd(); ++i)
                        context.state
                            .get_vec_ptr()[std::size_t(b0) * stride + i] = 0.01;
                app::DVec first, second;
                for (auto *vector : {&first, &second})
                    vector->create_vector(mesh.get(),
                                          ot::DVEC_TYPE::OCT_SHARED_NODES,
                                          ot::DVEC_LOC::HOST, dof, true);
                context.params.eta = 0.0;
                context.rhs(&context.state, &first, 1, 0.0);
                context.params.eta = 1.0;
                context.rhs(&context.state, &second, 1, 0.0);
                double parameter_error = 0.0, global_parameter_error = 0.0;
                if (mesh->isActive())
                    for (unsigned f = 0; f < dof; ++f)
                        for (unsigned i = mesh->getNodeLocalBegin();
                             i < mesh->getNodeLocalEnd(); ++i) {
                            const auto cell = std::size_t(f) * stride + i;
                            const double difference =
                                second.get_vec_ptr()[cell] -
                                first.get_vec_ptr()[cell];
                            if (!std::isfinite(difference))
                                throw std::runtime_error(
                                    "nonfinite parameter response");
                            parameter_error = std::max(
                                parameter_error,
                                std::abs(difference - (f == b0 ? -0.01 : 0.0)));
                        }  // END LOOP: for i over owned nodes
                first.destroy_vector();
                second.destroy_vector();
                MPI_Allreduce(&parameter_error, &global_parameter_error, 1,
                              MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                if (global_parameter_error > 1e-12)
                    throw std::runtime_error("runtime eta response mismatch");
                if (rank == 0)
                    std::printf("REAL_PARAMETER PASS error=%.17g\n",
                                global_parameter_error);
            }  // END ELSE: qualify transport and parameters
        }  // END BLOCK: own context before mesh destruction
    }  // END TRY: qualify real host integration
    catch (const std::exception &e) {
        std::fprintf(stderr, "rank %d: %s\n", rank, e.what());
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }  // END CATCH: abort all parent ranks
    MPI_Finalize();
    return 0;
}  // END FUNCTION: main
