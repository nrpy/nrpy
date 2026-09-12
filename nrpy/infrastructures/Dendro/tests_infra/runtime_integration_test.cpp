// Compile against a generated real context, not standalone host declarations.
// RUNTIME_HEADER and RUNTIME_NAMESPACE select the generated formulation.
#include RUNTIME_HEADER
#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include "meshUtils.h"
#include "octUtils.h"
namespace app = RUNTIME_NAMESPACE;
namespace {
double field(unsigned f, double x, double y, double z) {
    return 101.0 * (f + 1) + (f + 2) * x - (2 * f + 3) * y + (3 * f + 5) * z;
}  // END FUNCTION: field
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
        const unsigned order = 2 * app::generated::REQUIRED_PADDING;
        function2Octree(refine, octree, 5, 1e-3, order, MPI_COMM_WORLD);
        std::unique_ptr<ot::Mesh> mesh(
            ot::createMesh(octree.data(), octree.size(), order, MPI_COMM_WORLD,
                           0, ot::SM_TYPE::FDM));
        const double lo[3] = {-1, -2, -4}, hi[3] = {3, 2, 4};
        const Point minimum(lo[0], lo[1], lo[2]), maximum(hi[0], hi[1], hi[2]);
        mesh->setDomainBounds(minimum, maximum);
        {
            app::Ctx context(mesh.get(), minimum, maximum, 0.001);
            if (fault == "nonfinite") {
                context.initialize();
                if (rank == 1 && mesh->isActive())
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
                    }  // END LOOP: verify received ghost nodes
                    std::vector<double *> pointers(dof);
                    context.unzipped.to_2d(pointers.data());
                    for (const auto &b : mesh->getLocalBlockList()) {
                        auto g =
                            app::block_geometry(*mesh, b, minimum, maximum);
                        if (g.component_offset) ++offset_blocks;
                        if (fault == "offset" && rank == 1)
                            g.component_offset = 0;
                        // Expected coordinates derive independently from the
                        // raw octree.
                        const auto node      = b.getBlockNode();
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
                                    bool inside = true;
                                    for (unsigned a = 0; a < 3; ++a) {
                                        x[a] =
                                            lo[a] +
                                            (hi[a] - lo[a]) *
                                                (base[a] + (double(index[a]) -
                                                            b.get1DPadWidth()) *
                                                               grid_step) /
                                                double(1u << m_uiMaxDepth);
                                        inside &= x[a] >= lo[a] - 1e-12 &&
                                                  x[a] <= hi[a] + 1e-12;
                                        const double actual =
                                            g.pmin_padded[a] +
                                            index[a] * g.dx[a];
                                        error = std::max(
                                            error, std::abs(actual - x[a]));
                                    }  // END LOOP: recompute physical padded
                                       // coordinates
                                    if (!inside) continue;
                                    const bool halo = i < g.padding ||
                                                      j < g.padding ||
                                                      k < g.padding ||
                                                      i >= g.nx - g.padding ||
                                                      j >= g.ny - g.padding ||
                                                      k >= g.nz - g.padding;
                                    if (halo) ++halos;
                                    const std::size_t cell =
                                        g.component_offset + i +
                                        std::size_t(g.nx) *
                                            (j + std::size_t(g.ny) * k);
                                    for (unsigned f = 0; f < dof; ++f) {
                                        if (fault == "halo" && rank == 1 &&
                                            halo)
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
                                    }  // END LOOP: compare and transform
                                       // components
                                }  // END LOOP: inspect padded block points
                    }  // END LOOP: verify all local blocks
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
                        }  // END LOOP: check nonconstant zipped values
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
                const bool ok = ranks >= 2 && active_ranks >= 2 &&
                                max_blocks >= 2 && total_offsets > 0 &&
                                total_halos > 0 && total_remote > 0 &&
                                global_error < 1e-9;
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
                        }  // END LOOP: compare parameter response
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
