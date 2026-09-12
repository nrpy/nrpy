// Dendrolib capability mini-tests for the NRPy Dendro infrastructure.
//
// Each check proves one of the host-contract axes recorded on the Dendro
// validation page in the NRPy knowledge base.  Every check is
// written so that a host that does not honour the assumed contract fails it:
// the expected values are recomputed here from the block record and the
// physical domain, never read back from the same call under test.
#include <mpi.h>

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "block.h"
#include "dendro.h"
#include "mesh.h"
#include "meshUtils.h"
#include "octUtils.h"

namespace {

// Physical domain the octree is mapped onto; deliberately asymmetric and
// anisotropic so an axis swap or a dropped origin shift cannot pass.
const double kDomainMin[3] = {-1.0, -2.0, -4.0};
const double kDomainMax[3] = {3.0, 2.0, 4.0};

int g_failures = 0;

// Fault injection: every checker below must be demonstrated able to fail, or
// its passing result proves nothing.  CAPTEST_INJECT names one defect to
// introduce; the matching axis is then expected to report FAILED.
const char *injected() {
    const char *v = std::getenv("CAPTEST_INJECT");
    return v == nullptr ? "" : v;
}  // END FUNCTION: injected
bool inject(const char *what) { return std::string(injected()) == what; }

void report(const char *axis, bool ok, const std::string &detail) {
    if (!ok) ++g_failures;
    std::printf("%-16s %-6s %s\n", axis, ok ? "PROVEN" : "FAILED",
                detail.c_str());
}  // END FUNCTION: report

// Grid coordinate -> physical coordinate, the mapping Dendro-GR's
// GRIDX_TO_X performs, written out independently here.
double grid_to_phys(double xg, unsigned axis) {
    const double span = kDomainMax[axis] - kDomainMin[axis];
    return (span / static_cast<double>(1u << m_uiMaxDepth)) * xg +
           kDomainMin[axis];
}  // END FUNCTION: grid_to_phys

// Three independent linear fields.  Linear data is reproduced exactly by the
// unzip interpolation, so a mismatch anywhere -- interior or padding -- is a
// layout, origin, or ghost-validity defect rather than interpolation error.
double field(unsigned v, double x, double y, double z) {
    switch (v) {
        case 0:
            return 1.0 + 2.0 * x - 3.0 * y + 5.0 * z;
        case 1:
            return -7.0 + 11.0 * x + 0.5 * y - 13.0 * z;
        default:
            return 17.0 - 19.0 * x + 23.0 * y + 29.0 * z;
    }  // END SWITCH: field component
}  // END FUNCTION: field

// One element order: build a uniform mesh, unzip three fields, and check
// every axis against values recomputed from the block record.
bool run_order(unsigned eleOrder, unsigned level, MPI_Comm comm, bool verbose) {
    int rank = 0, npes = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);
    const bool faultRank = rank == (npes > 1 ? 1 : 0);
    // A wavelet-adaptive octree, so the block decomposition yields several
    // blocks whose padding is filled from neighbours rather than one block
    // whose padding is entirely outside the domain.
    std::vector<ot::TreeNode> octree;
    std::function<double(double, double, double)> refine =
        [](double x, double y, double z) {
            const double r2 = x * x + y * y + z * z;
            return std::exp(-r2 / 0.5);
        };
    function2Octree(refine, octree, level + 2, 1e-3, eleOrder, comm);
    ot::Mesh *mesh = ot::createMesh(octree.data(), octree.size(), eleOrder,
                                    comm, 0, ot::SM_TYPE::FDM);
    const int meshFailed = mesh == nullptr || (faultRank && inject("mesh"));
    int setupFailed = 0;
    MPI_Allreduce(&meshFailed, &setupFailed, 1, MPI_INT, MPI_MAX, comm);
    if (setupFailed) {
        if (meshFailed) report("mesh_setup", false, "mesh creation failed");
        // A partial mesh cannot be destroyed collectively: a missing mesh
        // cannot participate in its destructor's communicator free. Reclaim
        // it at process teardown after every rank reports qualification failure.
        return false;
    }  // END IF: mesh setup failed

    const Point pmin(kDomainMin[0], kDomainMin[1], kDomainMin[2]);
    const Point pmax(kDomainMax[0], kDomainMax[1], kDomainMax[2]);
    mesh->setDomainBounds(pmin, pmax);

    const unsigned dof = 3;
    const unsigned unzipSz = mesh->getDegOfFreedomUnZip();

    // createCGVector hands the callback physical coordinates, so the same
    // field function serves for the expected values below.
    std::function<void(double, double, double, double *)> fill =
        [](double x, double y, double z, double *out) {
            for (unsigned v = 0; v < 3; ++v) out[v] = field(v, x, y, z);
        };
    double *zipped = mesh->createCGVector<double>(fill, dof);
    double *unzipped = mesh->createUnZippedVector<double>(dof);
    if (faultRank && inject("zipped")) {
        delete[] zipped;
        zipped = nullptr;
    }  // END IF: inject zipped allocation failure
    if (faultRank && inject("unzipped")) {
        delete[] unzipped;
        unzipped = nullptr;
    }  // END IF: inject unzipped allocation failure
    const int vectorFailed = mesh->isActive() &&
                             (zipped == nullptr || unzipped == nullptr);
    MPI_Allreduce(&vectorFailed, &setupFailed, 1, MPI_INT, MPI_MAX, comm);
    if (setupFailed) {
        if (vectorFailed) report("vector_setup", false, "vector allocation failed");
        delete[] zipped;
        delete[] unzipped;
        delete mesh;
        return false;
    }  // END IF: vector allocation failed
    if (!mesh->isActive()) {
        delete mesh;
        return true;
    }  // END IF: no local work here
    for (std::size_t i = 0; i < static_cast<std::size_t>(unzipSz) * dof; ++i)
        unzipped[i] = 0.0;

    mesh->readFromGhostBegin(zipped, dof);
    mesh->readFromGhostEnd(zipped, dof);
    mesh->unzip(zipped, unzipped, dof);
    if (inject("ghost")) {
        // Blank the halo the host just filled.
        const std::vector<ot::Block> &bl = mesh->getLocalBlockList();
        for (std::size_t b = 0; b < bl.size(); ++b) {
            const unsigned pw = bl[b].get1DPadWidth();
            const unsigned sx = bl[b].getAllocationSzX();
            const unsigned sy = bl[b].getAllocationSzY();
            const unsigned sz2 = bl[b].getAllocationSzZ();
            for (unsigned v = 0; v < dof; ++v)
                for (unsigned k = 0; k < sz2; ++k)
                    for (unsigned j = 0; j < sy; ++j)
                        for (unsigned i = 0; i < sx; ++i)
                            if (i < pw || j < pw || k < pw || i >= sx - pw ||
                                j >= sy - pw || k >= sz2 - pw)
                                unzipped[static_cast<std::size_t>(v) * unzipSz +
                                         bl[b].getOffset() + i +
                                         static_cast<std::size_t>(sx) *
                                             (j + static_cast<std::size_t>(sy) *
                                                      k)] = 0.0;
        }  // END LOOP: for b over local blocks
    }  // END IF: halo injection requested

    const std::vector<ot::Block> &blkList = mesh->getLocalBlockList();
    const unsigned numBlocks = blkList.size();

    bool dimsOk = true, layoutOk = true, offsetsOk = true, ghostOk = true;
    bool originOk = true, padOk = true, mapOk = true;
    double worstInterior = 0.0, worstPadded = 0.0;
    unsigned badPadded = 0, badInterior = 0, checkedPadded = 0, badCorner = 0;
    double worstCorner = 0.0;

    for (unsigned blk = 0; blk < numBlocks; ++blk) {
        const ot::Block &b = blkList[blk];
        const unsigned pw = b.get1DPadWidth();
        const unsigned sz[3] = {b.getAllocationSzX(), b.getAllocationSzY(),
                                b.getAllocationSzZ()};
        const std::size_t off = b.getOffset();
        const std::size_t vol =
            static_cast<std::size_t>(sz[0]) * sz[1] * sz[2];

        if (pw != (eleOrder >> 1u) + (inject("padding") ? 1u : 0u))
            padOk = false;

        const unsigned elems1D =
            1u << (b.getRegularGridLev() - b.getBlockNode().getLevel());
        const unsigned expect1D =
            eleOrder * elems1D + 1 + 2 * pw + (inject("dimensions") ? 1u : 0u);
        if (sz[0] != expect1D || sz[1] != expect1D || sz[2] != expect1D)
            dimsOk = false;

        if (off + vol + (inject("offsets") ? unzipSz : 0u) > unzipSz)
            offsetsOk = false;
        for (unsigned other = 0; other < numBlocks; ++other) {
            if (other == blk) continue;
            const ot::Block &o = blkList[other];
            const std::size_t ovol =
                static_cast<std::size_t>(o.getAllocationSzX()) *
                o.getAllocationSzY() * o.getAllocationSzZ();
            if (off < o.getOffset() + ovol && o.getOffset() < off + vol)
                offsetsOk = false;
        }  // END LOOP: for other over local blocks

        const double dx[3] = {b.computeDx(pmin, pmax), b.computeDy(pmin, pmax),
                              b.computeDz(pmin, pmax)};

        // Cross-check the affine grid-to-physical map against the mesh's own
        // conversion before relying on it for the origin claim.
        Point checkPt;
        mesh->octCoordToDomainCoord(
            Point(static_cast<double>(b.getBlockNode().minX()),
                  static_cast<double>(b.getBlockNode().minY()),
                  static_cast<double>(b.getBlockNode().minZ())),
            checkPt);
        const double mine[3] = {grid_to_phys(b.getBlockNode().minX(), 0),
                                grid_to_phys(b.getBlockNode().minY(), 1),
                                grid_to_phys(b.getBlockNode().minZ(), 2)};
        if (std::fabs(checkPt.x() - mine[0]) > 1e-12 ||
            std::fabs(checkPt.y() - mine[1]) > 1e-12 ||
            std::fabs(checkPt.z() - mine[2]) > 1e-12)
            mapOk = false;

        const double padShift = inject("origin") ? 0.0 : 1.0;
        const double org[3] = {mine[0] - padShift * pw * dx[0],
                               mine[1] - padShift * pw * dx[1],
                               mine[2] - padShift * pw * dx[2]};

        for (unsigned v = 0; v < dof; ++v) {
            const std::size_t cornerIdx =
                static_cast<std::size_t>(v) * unzipSz + off + pw +
                static_cast<std::size_t>(sz[0]) *
                    (pw + static_cast<std::size_t>(sz[1]) * pw);
            // The known-bad input for this probe is the classic off-by-one-
            // padding error: treating padded index zero as the block corner.
            const double cornerShift = inject("origin") ? 1.0 : 0.0;
            const double cornerWant =
                field(v, mine[0] + cornerShift * pw * dx[0],
                      mine[1] + cornerShift * pw * dx[1],
                      mine[2] + cornerShift * pw * dx[2]);
            const double cornerErr = std::fabs(unzipped[cornerIdx] - cornerWant) /
                                     (1.0 + std::fabs(cornerWant));
            if (cornerErr > 1e-9) ++badCorner;
            if (cornerErr > worstCorner) worstCorner = cornerErr;

            for (unsigned k = 0; k < sz[2]; ++k) {
                for (unsigned j = 0; j < sz[1]; ++j) {
                    for (unsigned i = 0; i < sz[0]; ++i) {
                        // Transposing i and k turns the x-fastest claim into
                        // a z-fastest one; on a cube of equal extents that is
                        // exactly the defect the layout axis must catch.
                        const unsigned ii = inject("layout") ? k : i;
                        const unsigned kk = inject("layout") ? i : k;
                        const std::size_t idx =
                            static_cast<std::size_t>(v) * unzipSz + off + ii +
                            static_cast<std::size_t>(sz[0]) *
                                (j + static_cast<std::size_t>(sz[1]) * kk);
                        const double want =
                            field(v, org[0] + i * dx[0], org[1] + j * dx[1],
                                  org[2] + k * dx[2]);
                        const double err = std::fabs(unzipped[idx] - want) /
                                           (1.0 + std::fabs(want));
                        const bool padded =
                            (i < pw || j < pw || k < pw || i >= sz[0] - pw ||
                             j >= sz[1] - pw || k >= sz[2] - pw);
                        if (padded) {
                            // Padding outside the domain has no neighbour to
                            // supply it; only in-domain halo is a ghost claim.
                            const double px = org[0] + i * dx[0];
                            const double py = org[1] + j * dx[1];
                            const double pz = org[2] + k * dx[2];
                            const double eps = 1e-9;
                            if (px < kDomainMin[0] - eps ||
                                px > kDomainMax[0] + eps ||
                                py < kDomainMin[1] - eps ||
                                py > kDomainMax[1] + eps ||
                                pz < kDomainMin[2] - eps ||
                                pz > kDomainMax[2] + eps)
                                continue;
                            ++checkedPadded;
                            if (err > worstPadded) worstPadded = err;
                            if (err > 1e-9) ++badPadded;
                        } else {
                            if (err > worstInterior) worstInterior = err;
                            if (err > 1e-9) ++badInterior;
                        }  // END ELSE: interior point
                    }  // END LOOP: for i over padded x
                }  // END LOOP: for j over padded y
            }  // END LOOP: for k over padded z
        }  // END LOOP: for v over unzipped components
    }  // END LOOP: for blk over local blocks

    if (badInterior != 0) layoutOk = false;
    if (badCorner != 0) originOk = false;
    // A vacuous pass is not a proof: the axis fails if nothing was checked.
    if (badPadded != 0 || checkedPadded == 0) ghostOk = false;

    char detail[256];
    std::snprintf(detail, sizeof(detail), "eleOrder=%u pad=%u blocks=%u",
                  eleOrder, eleOrder >> 1u, numBlocks);
    report("padding", padOk, detail);
    report("dimensions", dimsOk, detail);
    report("offsets", offsetsOk, detail);
    report("grid_to_phys", mapOk, detail);
    std::snprintf(detail, sizeof(detail),
                  "eleOrder=%u interior_bad=%u worst=%.2e", eleOrder,
                  badInterior, worstInterior);
    report("layout", layoutOk, detail);
    std::snprintf(detail, sizeof(detail),
                  "eleOrder=%u corner_bad=%u worst=%.2e", eleOrder, badCorner,
                  worstCorner);
    report("origin", originOk, detail);
    std::snprintf(detail, sizeof(detail),
                  "eleOrder=%u halo_checked=%u halo_bad=%u worst=%.2e",
                  eleOrder, checkedPadded, badPadded, worstPadded);
    report("ghost_validity", ghostOk, detail);

    if (verbose)
        std::printf("  eleOrder=%u blocks=%u unzipSz=%u\n", eleOrder,
                    numBlocks, unzipSz);

    delete[] zipped;
    delete[] unzipped;
    delete mesh;
    return true;
}  // END FUNCTION: run_order

// clang-format off
}  // END NAMESPACE: internal linkage
// clang-format on

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0, npes = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    m_uiMaxDepth = 8;
    // Dendro's space-filling-curve tables must be initialized before any
    // octree construction; without this TreeNode ordering is undefined.
    _InitializeHcurve(3);

    if (rank == 0) {
        std::printf("Dendrolib capability mini-tests, %d rank(s)\n", npes);
        // Axis: scalar ABI.  The generated solver's DendroScalar contract.
        const bool scalarOk = sizeof(DendroScalar) == (inject("scalar") ? 4u : 8u) &&
                              std::is_same<DendroScalar, double>::value;
        char detail[128];
        std::snprintf(detail, sizeof(detail), "sizeof=%zu is_double=%d",
                      sizeof(DendroScalar),
                      static_cast<int>(std::is_same<DendroScalar, double>::value));
        report("scalar_abi", scalarOk, detail);
    }  // END IF: rank 0 reports scalar ABI

    // CAPTEST_ORDERS overrides the element orders under test; padding is
    // half the element order, so order 10 is the padding-5 probe that an
    // eighth-order finite-difference profile would need.
    std::vector<unsigned> orders = {2, 4, 6, 8};
    if (const char *env = std::getenv("CAPTEST_ORDERS"); rank == 0 && env != nullptr) {
        orders.clear();
        for (const char *s = env; *s != '\0';) {
            char *end = nullptr;
            errno = 0;
            const unsigned long v = std::strtoul(s, &end, 10);
            if (*s < '0' || *s > '9' || end == s || errno == ERANGE ||
                v < 2 || v % 2 != 0 || v > std::numeric_limits<unsigned>::max() ||
                (*end != '\0' && (*end != ',' || end[1] == '\0'))) {
                orders.clear();
                break;
            }  // END IF: invalid element order list
            orders.push_back(static_cast<unsigned>(v));
            s = (*end == '\0') ? end : end + 1;
        }  // END LOOP: for s over CAPTEST_ORDERS
        if (orders.empty())
            report("orders", false, "CAPTEST_ORDERS requires comma-separated positive even orders");
    }  // END IF: element orders overridden
    // Rank zero owns the request, so ranks cannot enter different order loops.
    unsigned orderCount = static_cast<unsigned>(orders.size());
    MPI_Bcast(&orderCount, 1, MPI_UNSIGNED, 0, comm);
    orders.resize(orderCount);
    MPI_Bcast(orders.data(), orderCount, MPI_UNSIGNED, 0, comm);
    for (unsigned o : orders) {
        if (!run_order(o, 3, comm, rank == 0)) {
            ++g_failures;
            if (rank == 0) std::printf("eleOrder=%u setup failed\n", o);
            break;
        }  // END IF: order did not complete
        MPI_Barrier(comm);
    }  // END LOOP: for o over element orders

    int total = 0;
    MPI_Allreduce(&g_failures, &total, 1, MPI_INT, MPI_SUM, comm);
    if (rank == 0)
        std::printf("%s\n", total == 0 ? "CAPABILITY_TESTS_OK"
                                       : "CAPABILITY_TESTS_FAILED");
    MPI_Finalize();
    return total == 0 ? 0 : 1;
}  // END FUNCTION: main
