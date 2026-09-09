# nrpy/infrastructures/Dendro/solver_context.py
"""
Emit standalone and real Dendro runtime contexts for a generated solver.

The context allocates the generated-count vectors, runs the startup checks and
invokes the registered generated CFunctions.  It carries no formulation
content: every count, name and constant comes from a generated header, and
every kernel name is read back from the Dendro role registry.

SCOPE OF EVIDENCE.  This vehicle establishes lifecycle plumbing -- allocation
with generated counts, the call path, the block and flat-block entry points
agreeing, per-rank decomposition, and an application-supplied fixed point.  It cannot
detect a uniform sign or scale error in the finite-difference coefficients:
such a kernel approximates a different continuum operator and still converges
at the requested order.  Pointwise correctness is established against an
independent evaluator, never here.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.grid as gri
from nrpy.infrastructures.Dendro import CFunction_roles as roles
from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.header_guards import header_guard

BANNER = generated_file_banner()

_HEADER = """// The umbrella header carries the host declarations behind their guard, the
// generated types, constants, parameters and state, and the CFunction
// prototypes this context calls, so the guard lives in exactly one place.
#include "$STEM_defines.h"

namespace $NAMESPACE {

class Ctx {
 public:
  ~Ctx();
  // Allocate this rank's standalone host mesh (n_blocks padded blocks of the given
  // extent and spacing) and the generated-count EVOL vectors, then run the
  // registered parameter CFunctions.  `rank` shifts this rank's subdomain so
  // that ranks own disjoint physical blocks rather than duplicating one.
  int initialize_mesh(int n_blocks, int extent, double dx, int rank,
                       const char* parfile_path);
  // Startup checks: scalar/state contracts, vector component counts, per-block
  // padding, FD order, parameter validation.
  int startup_checks(int rank);
$STANDALONE_APPLICATION_DECLARATIONS
  // Generated all-block RHS.
  int rhs_eval_all_blocks();
  // RHS magnitude at the current state, maximised over this rank's interior.
  double max_interior_rhs();
  // Largest difference between the per-block entry point and the LTS
  // flat-block adapter on the same state (one numerical body).
  double flat_adapter_max_difference();
  // Host-owned explicit-Euler lifecycle step: u += dt * f(u).  The host
  // integrator owns this loop exactly as the real Dendro integrator would; it
  // evaluates no formulation -- every derivative and every equation term comes
  // from the registered generated CFunctions.
  int euler_step(double dt);
  // Host reduction: maximum |field| over the interior of every local block,
  // over `ncomp` components.  The count is a parameter so the same reduction
  // serves the evolved and the diagnostic vectors without naming either.
  double max_interior_value(const $SCALAR* const* fields, unsigned ncomp);
  // Resolve host-supplied exact NRPy names (output or refinement selection).
  // An unknown name is fatal: this prints every valid generated name and
  // returns nonzero.  No name is written down here.
  int select_variables(const char* const* names, unsigned count);

  // Host mesh and EVOL vectors (in / rhs / out), standalone-host lifecycle.
  // Value-initialized: standalone_host::DVector is an aggregate with no default member
  // initializers, and ~Ctx frees all four vectors unconditionally.  When
  // initialize_mesh rejects its inputs it returns before assigning them, so
  // without the braces the destructor would free indeterminate pointers.
  standalone_host::Ctx host{};
  // Generated runtime parameter table, owned by the context.
  $NAMESPACE::generated::params_struct params;
$STANDALONE_APPLICATION_MEMBERS
};  // END CLASS: Ctx

$STANDALONE_APPLICATION_FREE_DECLARATIONS

// clang-format off
}  // END NAMESPACE: $NAMESPACE
// clang-format on
"""

_SOURCE = """// Host-owned lifecycle only.  Loops here are host reductions and integrator
// sweeps; they carry their own index names (`bx`, `by`, `bz`, `cell`) because
// `i0`/`i1`/`i2`/`blk_id` are reserved for NRPy-emitted numerical loops, which
// never appear in this emitter.

#include "$STEMCtx.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <optional>
#include <string_view>
#include <vector>

namespace $NAMESPACE {

namespace {

void swap_vectors(standalone_host::DVector& a, standalone_host::DVector& b) {
  $SCALAR** tmp_comp = a.comp;
  unsigned tmp_b = a.num_blocks;
  unsigned tmp_c = a.num_components;
  a.comp = b.comp;
  a.num_blocks = b.num_blocks;
  a.num_components = b.num_components;
  b.comp = tmp_comp;
  b.num_blocks = tmp_b;
  b.num_components = tmp_c;
}  // END FUNCTION: swap_vectors

void free_vector(standalone_host::DVector& vec) {
  if (vec.comp == nullptr) return;
  for (unsigned f = 0; f < vec.num_components; ++f) delete[] vec.comp[f];
  delete[] vec.comp;
  vec.comp = nullptr;
}  // END FUNCTION: free_vector

// clang-format off
}  // END NAMESPACE: internal linkage
// clang-format on

Ctx::~Ctx() {
  free_vector(host.in);
  free_vector(host.rhs);
  free_vector(host.out);
  free_vector(host.diag);
$STANDALONE_APPLICATION_DESTRUCTOR
}  // END FUNCTION: Ctx::~Ctx

int Ctx::initialize_mesh(int n_blocks, int extent, double dx, int rank,
                          const char* parfile_path) {
  const int pad = static_cast<int>($NAMESPACE::generated::REQUIRED_PADDING);
  if (n_blocks < 1 || n_blocks > static_cast<int>(standalone_host::MAX_STANDALONE_HOST_BLOCKS) ||
      extent < 2 * pad + 1 ||
      extent > static_cast<int>(standalone_host::MAX_STANDALONE_HOST_EXTENT) ||
      !(dx > 0.0) || !std::isfinite(dx)) {
    std::fprintf(stderr,
                 "ERROR: standalone-host mesh needs 1..%d blocks, extent in "
                 "[2*padding+1 (%d), %d] and a finite dx > 0; got n_blocks=%d "
                 "extent=%d dx=%g\\n",
                 static_cast<int>(standalone_host::MAX_STANDALONE_HOST_BLOCKS), 2 * pad + 1,
                 static_cast<int>(standalone_host::MAX_STANDALONE_HOST_EXTENT),
                 n_blocks, extent, dx);
    return 1;
  }  // END IF: mesh request rejected
  host.mesh.num_blocks = static_cast<unsigned>(n_blocks);
  // Each rank owns a disjoint physical subdomain: rank r starts where rank
  // r-1's blocks end.  Without this every rank would evolve an identical
  // private copy and a multi-rank run would prove nothing beyond a repeated
  // serial run.
  const double rank_origin =
      static_cast<double>(rank) * static_cast<double>(n_blocks) * extent * dx;
  for (int b = 0; b < n_blocks; ++b) {
    BlockGeometry& g = host.mesh.geom[b];
    g.nx = g.ny = g.nz = static_cast<unsigned>(extent);
    g.padding = static_cast<unsigned>(pad);
    g.component_offset =
        static_cast<std::size_t>(b) * extent * extent * extent;
    g.pmin_padded[0] = rank_origin + static_cast<double>(b) * extent * dx;
    g.pmin_padded[1] = 0.0;
    g.pmin_padded[2] = 0.0;
    g.dx[0] = g.dx[1] = g.dx[2] = dx;
    g.boundary_flags = 0;
  }  // END LOOP: for b over local blocks
  const unsigned ncomp = $NAMESPACE::generated::NUM_EVOL_GFS;
  const std::size_t vol = static_cast<std::size_t>(extent) * extent * extent;
  standalone_host::DVector* vectors[] = {&host.in, &host.rhs, &host.out, &host.diag};
  // The diagnostic vector carries the generated DIAG count, which is a
  // different cardinality from the evolved vectors.
  const unsigned counts[] = {ncomp, ncomp, ncomp,
                             $NAMESPACE::generated::NUM_DIAG_GFS};
  for (unsigned v = 0; v < sizeof(vectors) / sizeof(vectors[0]); ++v) {
    standalone_host::DVector* vec = vectors[v];
    vec->num_blocks = static_cast<unsigned>(n_blocks);
    vec->num_components = counts[v];
    vec->comp = new $SCALAR*[counts[v]];
    for (unsigned f = 0; f < counts[v]; ++f) {
      vec->comp[f] = new $SCALAR[static_cast<std::size_t>(n_blocks) * vol];
      std::memset(vec->comp[f], 0,
                  sizeof($SCALAR) * static_cast<std::size_t>(n_blocks) * vol);
    }  // END LOOP: for f over vector components
  }  // END LOOP: for v over host vectors
  // The registered parameter CFunctions own the parameter lifecycle.
  $PARAMS_STRUCT_SET_TO_DEFAULT(params);
$STANDALONE_APPLICATION_POST_MESH
  if (parfile_path != nullptr) {
    std::fprintf(stderr,
                 "ERROR: a parameter file was supplied (%s), but parameter-file "
                 "parsing belongs to the Dendro-GR host, which reads its own "
                 "TOML and calls the generated setters.  The standalone build "
                 "runs on the compiled defaults; rerun without -t.\\n",
                 parfile_path);
    return 1;
  }  // END IF: a parameter file was supplied
  return 0;
}  // END FUNCTION: Ctx::initialize_mesh

int Ctx::startup_checks(int rank) {
  // The generated scalar/state/padding/FD contracts hold for every local block
  // and vector.
  if (host.in.num_components != $NAMESPACE::generated::NUM_EVOL_GFS ||
      host.rhs.num_components != $NAMESPACE::generated::NUM_EVOL_GFS ||
      host.out.num_components != $NAMESPACE::generated::NUM_EVOL_GFS ||
      host.diag.num_components != $NAMESPACE::generated::NUM_DIAG_GFS) {
    std::fprintf(stderr, "ERROR: vector component count mismatch\\n");
    return 1;
  }  // END IF: component count mismatch
  for (unsigned b = 0; b < host.mesh.num_blocks; ++b) {
    if (host.mesh.geom[b].padding <
        $NAMESPACE::generated::REQUIRED_PADDING) {
      std::fprintf(stderr,
                   "ERROR: block %u padding below generated minimum\\n", b);
      return 1;
    }  // END IF: block padding below minimum
  }  // END LOOP: for b over local blocks
$STANDALONE_APPLICATION_STARTUP_CHECKS
  if (!$VALIDATE(params)) {
    std::fprintf(stderr, "ERROR: parameter validation failed\\n");
    return 1;
  }  // END IF: parameter validation failed
  // Rank-0 only, as every other diagnostic in the entry point is: an N-rank
  // run printing N identical parameter tables buries the gate output.
  if (rank == 0) $PRINT_EFFECTIVE(params);
  return 0;
}  // END FUNCTION: Ctx::startup_checks

$STANDALONE_APPLICATION_INITIALIZATION

int Ctx::rhs_eval_all_blocks() {
  $RHS_EVAL(host.mesh, host.in.comp, host.rhs.comp$RHS_EVAL_TAIL);
  return 0;
}  // END FUNCTION: Ctx::rhs_eval_all_blocks

double Ctx::max_interior_rhs() {
  rhs_eval_all_blocks();
  return max_interior_value(host.rhs.comp,
                            $NAMESPACE::generated::NUM_EVOL_GFS);
}  // END FUNCTION: Ctx::max_interior_rhs

double Ctx::flat_adapter_max_difference() {
  // The LTS flat-block adapter must invoke the same numerical body as the
  // per-block entry point.  Evaluate both on this rank's first block and
  // compare pointwise.
  const unsigned ncomp = $NAMESPACE::generated::NUM_EVOL_GFS;
  const BlockGeometry& g = host.mesh.geom[0];
  const std::size_t vol = static_cast<std::size_t>(g.nx) * g.ny * g.nz;
  std::vector<$SCALAR> flat_in(static_cast<std::size_t>(ncomp) * vol, 0.0);
  std::vector<$SCALAR> flat_rhs(static_cast<std::size_t>(ncomp) * vol, 0.0);
  for (unsigned f = 0; f < ncomp; ++f) {
    std::memcpy(&flat_in[static_cast<std::size_t>(f) * vol],
                host.in.comp[f] + g.component_offset,
                sizeof($SCALAR) * vol);
  }  // END LOOP: for f over evolved components
  BlockGeometry flat_geom = g;
  flat_geom.component_offset = 0;
  $RHS_EVAL_FLAT_BLOCK(flat_geom, flat_in.data(), flat_rhs.data()$RHS_EVAL_TAIL);
  std::vector<const $SCALAR*> block_in(ncomp);
  std::vector<$SCALAR*> block_rhs(ncomp);
  std::vector<std::vector<$SCALAR>> storage(
      ncomp, std::vector<$SCALAR>(vol, 0.0));
  for (unsigned f = 0; f < ncomp; ++f) {
    block_in[f] = host.in.comp[f] + g.component_offset;
    block_rhs[f] = storage[f].data();
  }  // END LOOP: for f over evolved components
  $RHS_EVAL_BLOCK(flat_geom, block_in.data(), block_rhs.data()$RHS_EVAL_TAIL);
  double worst = 0.0;
  for (unsigned f = 0; f < ncomp; ++f) {
    for (unsigned bz = g.padding; bz < g.nz - g.padding; ++bz) {
      for (unsigned by = g.padding; by < g.ny - g.padding; ++by) {
        for (unsigned bx = g.padding; bx < g.nx - g.padding; ++bx) {
          const std::size_t cell = bx + g.nx * (by + g.ny * bz);
          const double d = std::fabs(
              static_cast<double>(block_rhs[f][cell]) -
              static_cast<double>(
                  flat_rhs[static_cast<std::size_t>(f) * vol + cell]));
          if (d > worst) worst = d;
        }  // END LOOP: for bx over interior x
      }  // END LOOP: for by over interior y
    }  // END LOOP: for bz over interior z
  }  // END LOOP: for f over evolved components
  return worst;
}  // END FUNCTION: Ctx::flat_adapter_max_difference

int Ctx::euler_step(double dt) {
  rhs_eval_all_blocks();
  // Host integrator: u_out = u_in + dt * rhs.  No formulation content here:
  // every derivative and equation term was computed by the registered
  // generated CFunctions.
  const unsigned ncomp = $NAMESPACE::generated::NUM_EVOL_GFS;
  const unsigned nb = host.in.num_blocks;
  const unsigned extent = host.mesh.geom[0].nx;
  const std::size_t vol = static_cast<std::size_t>(extent) * extent * extent;
  const std::size_t total = static_cast<std::size_t>(nb) * vol;
  for (unsigned f = 0; f < ncomp; ++f) {
    const $SCALAR* u = host.in.comp[f];
    const $SCALAR* f_rhs = host.rhs.comp[f];
    $SCALAR* u_new = host.out.comp[f];
    for (std::size_t cell = 0; cell < total; ++cell) {
      u_new[cell] = u[cell] + static_cast<$SCALAR>(dt) * f_rhs[cell];
    }
  }  // END LOOP: for f over evolved components
  swap_vectors(host.in, host.out);
  return 0;
}  // END FUNCTION: Ctx::euler_step

double Ctx::max_interior_value(const $SCALAR* const* fields, unsigned ncomp) {
  double worst = 0.0;
  for (unsigned b = 0; b < host.mesh.num_blocks; ++b) {
    const BlockGeometry& g = host.mesh.geom[b];
    const std::size_t base = g.component_offset;
    for (unsigned f = 0; f < ncomp; ++f) {
      const $SCALAR* v = fields[f];
      for (unsigned bz = g.padding; bz < g.nz - g.padding; ++bz) {
        for (unsigned by = g.padding; by < g.ny - g.padding; ++by) {
          for (unsigned bx = g.padding; bx < g.nx - g.padding; ++bx) {
            const double a = std::fabs(static_cast<double>(
                v[base + bx + g.nx * (by + g.ny * bz)]));
            if (a > worst) worst = a;
          }  // END LOOP: for bx over interior x
        }  // END LOOP: for by over interior y
      }  // END LOOP: for bz over interior z
    }  // END LOOP: for f over vector components
  }  // END LOOP: for b over local blocks
  return worst;
}  // END FUNCTION: Ctx::max_interior_value

$STANDALONE_APPLICATION_AFTER_RHS

int Ctx::select_variables(const char* const* names, unsigned count) {
  // Output and refinement candidates are selected by exact NRPy name.
  // Matching is case-sensitive, and an unknown name is fatal.
  int rc = 0;
  for (unsigned k = 0; k < count; ++k) {
    const std::optional<$NAMESPACE::generated::VariableRef> ref =
        $NAMESPACE::generated::find_variable(names[k]);
    if (!ref.has_value()) {
      std::fprintf(stderr, "ERROR: unknown variable name '%s'\\n", names[k]);
      rc = 1;
      continue;
    }  // END IF: name did not resolve
    std::printf("SELECTED %s group=%u index=%u\\n", names[k],
                static_cast<unsigned>(ref->group), ref->index);
  }  // END LOOP: for k over selected names
  if (rc != 0) {
    std::fprintf(stderr, "Valid generated variable names:\\n");
    for (const $NAMESPACE::generated::VariableRef::Group group :
         $NAMESPACE::generated::VARIABLE_GROUPS) {
      const unsigned n = $NAMESPACE::generated::variable_count(group);
      for (unsigned index = 0; index < n; ++index) {
        const std::string_view name =
            $NAMESPACE::generated::variable_name(group, index);
        std::fprintf(stderr, "  %.*s\\n", static_cast<int>(name.size()),
                     name.data());
      }  // END LOOP: for index over group names
    }  // END LOOP: for group over generated groups
  }  // END IF: an unknown name was supplied
  return rc;
}  // END FUNCTION: Ctx::select_variables

$STANDALONE_APPLICATION_FREE_DEFINITIONS

// clang-format off
}  // END NAMESPACE: $NAMESPACE
// clang-format on
"""


_REAL_HEADER = r"""#include "$STEM_defines.h"
#include "ctx.h"
#include "meshUtils.h"
#include <vector>
namespace $NAMESPACE {
using DVec = ot::DVector<DendroScalar, unsigned int>;
/**
 * Normalize a real block allocation, component offset, and physical padded origin.
 *
 * @param[in] mesh Mesh providing the unzip stride and coordinate transformation.
 * @param[in] block Local block whose padded allocation is described.
 * @param[in] minimum Physical domain minimum used to compute spacing.
 * @param[in] maximum Physical domain maximum used to compute spacing.
 * @return Validated geometry with an offset relative to each component base.
 *
 * @note Throws std::runtime_error for invalid allocation, padding, or geometry.
 */
BlockGeometry block_geometry(const ot::Mesh& mesh, const ot::Block& block,
                             const Point& minimum, const Point& maximum);
// The fixed-mesh context owns storage. DVector itself is a shallow handle.
class Ctx : public ts::Ctx<Ctx, DendroScalar, unsigned int> {
 public:
  generated::params_struct params{};
  DVec state, unzipped, unzipped_rhs, diagnostics;
  /**
   * Allocate owned vectors and communication buffers for a fixed mesh.
   *
   * @param[in,out] mesh Borrowed mesh used for storage and communication setup.
   * @param[in] minimum Physical domain minimum, matching the mesh domain bounds.
   * @param[in] maximum Physical domain maximum, matching the mesh domain bounds.
   * @param dt Positive finite timestep used by the host time integrator.
   *
   * @note The caller retains mesh ownership and must outlive this context.
   */
  Ctx(ot::Mesh* mesh, const Point& minimum, const Point& maximum, double dt);
  ~Ctx();
  Ctx(const Ctx&) = delete;
  Ctx& operator=(const Ctx&) = delete;
  DVec& get_evolution_vars() { return state; }
$REAL_APPLICATION_DECLARATIONS
  /**
   * Exchange halos, evaluate generated block RHS kernels, and zip the result.
   *
   * @param[in,out] in One packed evolution vector; unzip updates its ghost nodes.
   * @param[out] out One packed vector receiving the zipped RHS.
   * @param count Number of packed evolution vectors; must equal one.
   * @param time Host stage time; unused by this autonomous generated profile.
   * @return 0 on success or an inactive rank; invalid data aborts MPI_COMM_WORLD.
   */
  int rhs(DVec* in, DVec* out, unsigned int count, DendroScalar time);
 private:
  /**
   * Prescribe application values outside the physical domain only.
   *
   * @note Updates unzipped exterior points; preserves all in-domain halo values.
   */
  void fill_exterior();
  void require_finite(DVec& value);
}; // END CLASS: fixed mesh Dendro context
// clang-format off
}  // END NAMESPACE: $NAMESPACE
// clang-format on
"""

_REAL_SOURCE = r"""#include "$STEMCtx.h"
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <stdexcept>
namespace $NAMESPACE {
namespace {
[[noreturn]] void fail(const char* reason) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::fprintf(stderr, "rank %d: %s\n", rank, reason);
  // ETS ignores callback return codes. Abort the parent communicator, including
  // inactive ranks, so one failed rank cannot strand peers in a halo exchange.
  MPI_Abort(MPI_COMM_WORLD, 1);
  std::abort();
} // END FUNCTION: abort all parent ranks
double maximum(double local) {
  double global = 0.0;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return global;
} // END FUNCTION: reduce world maximum
/**
 * Measure finite generated values over block interiors.
 *
 * @param[in] mesh Mesh providing local block interiors and component offsets.
 * @param[in] values Unzipped generated values inspected without modification.
 * @return Local maximum absolute value, or zero on an inactive rank.
 *
 * @note Nonfinite interior values abort MPI_COMM_WORLD.
 */
double interior_max(const ot::Mesh& mesh, DVec& values) {
  double result = 0.0;
  if (!mesh.isActive()) return result;
  std::vector<DendroScalar*> pointers(values.get_dof());
  values.to_2d(pointers.data());
  for (const auto& b : mesh.getLocalBlockList()) {
    const unsigned p = b.get1DPadWidth();
    const unsigned nx = b.getAllocationSzX(), ny = b.getAllocationSzY();
    for (unsigned f = 0; f < values.get_dof(); ++f)
      for (unsigned k = p; k < b.getAllocationSzZ() - p; ++k)
        for (unsigned j = p; j < ny - p; ++j)
          for (unsigned i = p; i < nx - p; ++i) {
            const auto v = pointers[f][b.getOffset() + i + std::size_t(nx) * (j + std::size_t(ny) * k)];
            if (!std::isfinite(v)) fail("nonfinite generated output");
            result = std::max(result, std::abs(v));
          } // END LOOP: inspect interior values
  } // END LOOP: visit local blocks
  return result;
} // END FUNCTION: measure interior maximum
// clang-format off
}  // END NAMESPACE: internal linkage
// clang-format on
BlockGeometry block_geometry(const ot::Mesh& mesh, const ot::Block& block,
                             const Point& minimum, const Point& maximum) {
  BlockGeometry g{};
  g.nx = block.getAllocationSzX();
  g.ny = block.getAllocationSzY();
  g.nz = block.getAllocationSzZ();
  g.padding = block.get1DPadWidth();
  g.component_offset = block.getOffset();
  g.dx[0] = block.computeDx(minimum, maximum);
  g.dx[1] = block.computeDy(minimum, maximum);
  g.dx[2] = block.computeDz(minimum, maximum);
  Point origin;
  const auto node = block.getBlockNode();
  mesh.octCoordToDomainCoord(Point(node.minX(), node.minY(), node.minZ()), origin);
  g.pmin_padded[0] = origin.x() - g.padding * g.dx[0];
  g.pmin_padded[1] = origin.y() - g.padding * g.dx[1];
  g.pmin_padded[2] = origin.z() - g.padding * g.dx[2];
  g.boundary_flags = block.getBlkNodeFlag();
  const std::size_t volume = std::size_t(g.nx) * g.ny * g.nz;
  if (g.padding < generated::REQUIRED_PADDING || g.nx <= 2*g.padding ||
      g.ny <= 2*g.padding || g.nz <= 2*g.padding ||
      g.component_offset > mesh.getDegOfFreedomUnZip() ||
      volume > mesh.getDegOfFreedomUnZip() - g.component_offset)
    throw std::runtime_error("invalid Dendro block allocation or padding");
  for (unsigned a = 0; a < 3; ++a)
    if (!(g.dx[a] > 0.0) || !std::isfinite(g.dx[a]) || !std::isfinite(g.pmin_padded[a]))
      throw std::runtime_error("invalid physical block geometry");
  return g;
} // END FUNCTION: normalize real block geometry
Ctx::Ctx(ot::Mesh* mesh, const Point& minimum, const Point& maximum, double dt) {
  set_mesh(mesh);
  m_uiElementOrder = mesh->getElementOrder();
  m_uiMinPt = minimum;
  m_uiMaxPt = maximum;
  m_uiTinfo = {};
  m_uiTinfo._m_uiTh = dt;
  $PARAMS_STRUCT_SET_TO_DEFAULT(params);
  state.create_vector(mesh, ot::DVEC_TYPE::OCT_SHARED_NODES, ot::DVEC_LOC::HOST, generated::NUM_EVOL_GFS, true);
  unzipped.create_vector(mesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST, generated::NUM_EVOL_GFS, true);
  unzipped_rhs.create_vector(mesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST, generated::NUM_EVOL_GFS, true);
  diagnostics.create_vector(mesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST, generated::NUM_DIAG_GFS, true);
  if (mesh->isActive()) {
    for (DVec* v : {&state, &unzipped, &unzipped_rhs, &diagnostics})
      std::fill_n(v->get_vec_ptr(), v->get_size(), 0.0);
    for (const auto& block : mesh->getLocalBlockList())
      block_geometry(*mesh, block, minimum, maximum);
  } // END IF: initialize active storage
  ot::alloc_mpi_ctx<DendroScalar>(mesh, m_mpi_ctx, generated::NUM_EVOL_GFS, 1);
} // END FUNCTION: construct fixed mesh context
Ctx::~Ctx() {
  ot::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, generated::NUM_EVOL_GFS, 1);
  for (DVec* v : {&state, &unzipped, &unzipped_rhs, &diagnostics}) v->destroy_vector();
} // END FUNCTION: release owned host storage
void Ctx::require_finite(DVec& value) {
  if (!m_uiMesh->isActive()) return;
  std::vector<DendroScalar*> pointers(value.get_dof());
  value.to_2d(pointers.data());
  for (const auto* ptr : pointers)
    for (unsigned i = m_uiMesh->getNodeLocalBegin(); i < m_uiMesh->getNodeLocalEnd(); ++i)
      if (!std::isfinite(ptr[i])) fail("nonfinite evolved state");
} // END FUNCTION: validate owned evolved values
$REAL_APPLICATION_INITIALIZATION
void Ctx::fill_exterior() {
  // Only points outside the physical domain are prescribed. Interior halos
  // remain the output of real unzip. The application supplies field values;
  // coordinate construction, classification, traversal, and writes stay here.
  std::vector<DendroScalar*> pointers(generated::NUM_EVOL_GFS);
  unzipped.to_2d(pointers.data());
  std::vector<DendroScalar> flat(generated::NUM_EVOL_GFS);
$APPLICATION_EXTERIOR_VALUES
  const double low[3] = {m_uiMinPt.x(), m_uiMinPt.y(), m_uiMinPt.z()};
  const double high[3] = {m_uiMaxPt.x(), m_uiMaxPt.y(), m_uiMaxPt.z()};
  for (const auto& b : m_uiMesh->getLocalBlockList()) {
    const auto g = block_geometry(*m_uiMesh, b, m_uiMinPt, m_uiMaxPt);
    for (unsigned k = 0; k < g.nz; ++k)
      for (unsigned j = 0; j < g.ny; ++j)
        for (unsigned i = 0; i < g.nx; ++i) {
          const double x[3] = {g.pmin_padded[0] + i*g.dx[0], g.pmin_padded[1] + j*g.dx[1], g.pmin_padded[2] + k*g.dx[2]};
          bool exterior = false;
          for (unsigned a = 0; a < 3; ++a) exterior |= x[a] < low[a] - 1e-12*g.dx[a] || x[a] > high[a] + 1e-12*g.dx[a];
          if (exterior)
            for (unsigned f = 0; f < flat.size(); ++f)
              pointers[f][g.component_offset + i + std::size_t(g.nx)*(j + std::size_t(g.ny)*k)] = flat[f];
        } // END LOOP: prescribe exterior points
  } // END LOOP: visit exterior block padding
} // END FUNCTION: fill application exterior
int Ctx::rhs(DVec* in, DVec* out, unsigned int count, DendroScalar) {
  if (count != 1) fail("RHS requires one packed evolution vector");
  if (!m_uiMesh->isActive()) return 0;
  require_finite(*in);
  unzip(*in, unzipped, 1);
  fill_exterior();
  std::fill_n(unzipped_rhs.get_vec_ptr(), unzipped_rhs.get_size(), 0.0);
  std::vector<DendroScalar*> input(generated::NUM_EVOL_GFS), output(generated::NUM_EVOL_GFS);
  unzipped.to_2d(input.data());
  unzipped_rhs.to_2d(output.data());
  for (const auto& b : m_uiMesh->getLocalBlockList()) {
    const auto g = block_geometry(*m_uiMesh, b, m_uiMinPt, m_uiMaxPt);
    $RHS_EVAL_BLOCK(g, input.data(), output.data()$RHS_EVAL_BLOCK_TAIL);
  } // END LOOP: evaluate generated block RHS
  interior_max(*m_uiMesh, unzipped_rhs);
  zip(unzipped_rhs, *out);
  require_finite(*out);
  return 0;
} // END FUNCTION: unzip evaluate and zip
$REAL_APPLICATION_AFTER_RHS
// clang-format off
}  // END NAMESPACE: $NAMESPACE
// clang-format on
"""


def _codeparameter_tail(cfunction_name: str, table: str) -> str:
    """
    Render the trailing CodeParameter arguments one host call must forward.

    :param cfunction_name: Registered CFunction the host calls.
    :param table: C++ expression naming the generated parameter table.
    :return: ``", table.a, table.b"``, or the empty string when the signature
        declares no CodeParameter.
    """
    names = roles.CFunction_codeparameters(cfunction_name)
    if not names:
        return ""
    return ", " + ", ".join(f"{table}.{name}" for name in names)


def substitute_solver_identifiers(
    text: str,
    solver_stem: str,
    solver_namespace: str,
) -> str:
    """
    Substitute the solver identifiers and registered kernel names into one body.

    Every kernel name is read back from the Dendro role registry, so a renamed
    kernel cannot leave a stale call behind in the host adapter.

    :param text: The emitter body carrying ``$`` placeholders.
    :param solver_stem: Lowercase formulation stem for emitted file names.
    :param solver_namespace: Solver namespace.
    :return: The substituted text.

    Doctests:
    >>> substitute_solver_identifiers("namespace $NAMESPACE { using T = $SCALAR; }", "wave", "wave")
    'namespace wave { using T = DendroScalar; }'
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> _saved_cfuncs = dict(cfc.CFunction_dict)
    >>> _saved_dendro = par.glb_extras_dict.get("Dendro")
    >>> cfc.CFunction_dict.clear()
    >>> par.glb_extras_dict.pop("Dendro", None)
    >>> cfc.register_CFunction(desc="rhs", name="wave_rhs", params="", body="(void)0;")
    >>> roles.set_CFunction_role("wave_rhs", "rhs_eval")
    >>> roles.set_CFunction_codeparameters(
    ...     "wave_rhs", ("amplitude", "num_steps", "enable_filter")
    ... )
    >>> _codeparameter_tail("wave_rhs", "runtime")
    ', runtime.amplitude, runtime.num_steps, runtime.enable_filter'
    >>> substitute_solver_identifiers("call $RHS_EVAL$RHS_EVAL_TAIL;", "wave", "wave")
    'call wave_rhs, params.amplitude, params.num_steps, params.enable_filter;'
    >>> cfc.CFunction_dict.clear(); cfc.CFunction_dict.update(_saved_cfuncs)
    >>> _ = par.glb_extras_dict.pop("Dendro", None)
    >>> _ = par.glb_extras_dict.setdefault("Dendro", _saved_dendro) if _saved_dendro is not None else None
    """
    replacements = (
        (
            "$PARAMS_STRUCT_SET_TO_DEFAULT",
            f"{solver_stem}_params_struct_set_to_default",
        ),
        ("$VALIDATE", f"{solver_stem}_params_validate"),
        ("$PRINT_EFFECTIVE", f"{solver_stem}_params_print_effective"),
        ("$NAMESPACE", solver_namespace),
        ("$STEM", solver_stem),
        ("$SCALAR", gri.DENDRO_SCALAR_TYPE),
    )
    for placeholder, value in replacements:
        text = text.replace(placeholder, value)

    role_tokens = (
        ("$RHS_EVAL_FLAT_BLOCK", "rhs_eval_flat_block", False),
        ("$RHS_EVAL_BLOCK_TAIL", "rhs_eval_block", True),
        ("$RHS_EVAL_BLOCK", "rhs_eval_block", False),
        ("$RHS_EVAL_TAIL", "rhs_eval", True),
        ("$RHS_EVAL", "rhs_eval", False),
    )
    for placeholder, role, is_tail in role_tokens:
        if placeholder not in text:
            continue
        function_name = roles.CFunction_name_for_role(role)
        value = (
            _codeparameter_tail(function_name, "params") if is_tail else function_name
        )
        text = text.replace(placeholder, value)
    return text


def output_solver_context_h(
    solver_stem: str,
    solver_namespace: str,
    standalone_application_declarations: str,
    standalone_application_members: str,
    standalone_application_free_declarations: str,
    real_application_declarations: str,
) -> str:
    r"""
    Emit the context header with explicit host selection.

    :param solver_stem: Lowercase formulation stem for emitted file names.
    :param solver_namespace: Solver namespace.
    :param standalone_application_declarations: Application methods for the
        standalone context public interface.
    :param standalone_application_members: Application-owned standalone state.
    :param standalone_application_free_declarations: Application free-function
        declarations in the solver namespace.
    :param real_application_declarations: Real-host application methods and
        state, inserted before the context's private generic helpers.
    :return: The complete C++ header text.
    :raises ValueError: If an application declaration remains unresolved.

    Doctests:
    >>> _HEADER.count("// clang-format off") == _HEADER.count("}  // END NAMESPACE:")
    True
    >>> _HEADER.count("}  // END NAMESPACE:")
    1
    >>> _SOURCE.count("// clang-format off") == _SOURCE.count("}  // END NAMESPACE:")
    True
    >>> _SOURCE.count("}  // END NAMESPACE:")
    2
    >>> unrelated = output_solver_context_h(
    ...     "wave", "wave", "  int initialize_scalar();", "", "", "  int initialize();"
    ... )
    >>> "initialize_scalar" in unrelated and "detgtrazero" not in unrelated
    True
    """
    opening, closing = header_guard(f"{solver_stem}Ctx.h")
    body = (
        "#if defined(NRPY_DENDRO_STANDALONE_HOST)\n"
        + _HEADER
        + "\n#else\n"
        + _REAL_HEADER
        + "\n#endif\n"
    )
    for token, value in (
        ("$STANDALONE_APPLICATION_DECLARATIONS", standalone_application_declarations),
        ("$STANDALONE_APPLICATION_MEMBERS", standalone_application_members),
        (
            "$STANDALONE_APPLICATION_FREE_DECLARATIONS",
            standalone_application_free_declarations,
        ),
        ("$REAL_APPLICATION_DECLARATIONS", real_application_declarations),
    ):
        body = body.replace(token, value)
    unresolved = tuple(
        token
        for token in (
            "$STANDALONE_APPLICATION_DECLARATIONS",
            "$STANDALONE_APPLICATION_MEMBERS",
            "$STANDALONE_APPLICATION_FREE_DECLARATIONS",
            "$REAL_APPLICATION_DECLARATIONS",
        )
        if token in body
    )
    if unresolved:
        raise ValueError(f"Application declarations were not resolved: {unresolved}")
    return (
        BANNER
        + f"{opening}\n\n"
        + substitute_solver_identifiers(
            body,
            solver_stem,
            solver_namespace,
        )
        + f"\n{closing}\n"
    )


def output_solver_context_cpp(
    solver_stem: str,
    solver_namespace: str,
    standalone_application_destructor: str,
    standalone_application_post_mesh: str,
    standalone_application_startup_checks: str,
    standalone_application_initialization: str,
    standalone_application_after_rhs: str,
    standalone_application_free_definitions: str,
    real_application_initialization: str,
    real_application_exterior_values: str,
    real_application_after_rhs: str,
) -> str:
    r"""
    Emit the context implementation with explicit host selection.

    :param solver_stem: Lowercase formulation stem for emitted file names.
    :param solver_namespace: Solver namespace.
    :param standalone_application_destructor: Cleanup of application-owned
        standalone context state.
    :param standalone_application_post_mesh: Application setup after generic
        standalone mesh and parameter initialization.
    :param standalone_application_startup_checks: Additional application
        qualification checks.
    :param standalone_application_initialization: Application initialization
        method definitions before the generic RHS methods.
    :param standalone_application_after_rhs: Application diagnostics and
        accepted-step policy definitions after generic reductions.
    :param standalone_application_free_definitions: Application free-function
        definitions emitted before the standalone namespace closes.
    :param real_application_initialization: Application initialization method
        definitions inserted after generic storage setup.
    :param real_application_exterior_values: Application statements that fill
        the existing ``flat`` value vector before generic exterior traversal.
    :param real_application_after_rhs: Application hook and diagnostic method
        definitions inserted after the generic RHS callback.
    :return: The complete C++ source text.
    :raises ValueError: If an application insertion remains unresolved.

    An unrelated scalar/vector application crosses every generic assembly
    boundary without importing GR policy.  The live registries and sidecars are
    restored by object identity afterward.

    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.infrastructures.Dendro import main_cpp as generic_main
    >>> from nrpy.infrastructures.Dendro import state_h, types_h
    >>> _saved_infrastructure = par.parval_from_str("Infrastructure")
    >>> _saved_fd_order = par.parval_from_str("fd_order")
    >>> _saved_fields = dict(gri.glb_gridfcs_dict)
    >>> _saved_parameters = dict(par.glb_code_params_dict)
    >>> _saved_cfuncs = dict(cfc.CFunction_dict)
    >>> _saved_extras = dict(par.glb_extras_dict)
    >>> _saved_dendro_present = "Dendro" in par.glb_extras_dict
    >>> _saved_dendro = par.glb_extras_dict.pop("Dendro", None)
    >>> _saved_roles_present = _saved_dendro is not None and "CFunction_roles" in _saved_dendro
    >>> _saved_roles = None if _saved_dendro is None else _saved_dendro.get("CFunction_roles")
    >>> _saved_codeparameters_present = _saved_dendro is not None and "CFunction_codeparameters" in _saved_dendro
    >>> _saved_codeparameters = None if _saved_dendro is None else _saved_dendro.get("CFunction_codeparameters")
    >>> try:
    ...     gri.glb_gridfcs_dict.clear()
    ...     cfc.CFunction_dict.clear()
    ...     par.set_parval_from_str("Infrastructure", "Dendro")
    ...     _ = gri.register_gridfunctions("wave_scalar", group="EVOL")
    ...     _ = gri.register_gridfunctions_for_single_rank1(
    ...         "waveU", dimension=2, group="EVOL"
    ...     )
    ...     roles.set_upwind_control_fields(("waveU0", "waveU1"))
    ...     _function_specs = (
    ...         ("wave_rhs", "rhs_eval",
    ...          "const StandaloneHostMesh& mesh, const DendroScalar* const* in_gfs, DendroScalar* const* rhs_gfs",
    ...          "(void)mesh; rhs_gfs[0][0] = in_gfs[0][0] + in_gfs[1][0];"),
    ...         ("wave_rhs_block", "rhs_eval_block",
    ...          "const BlockGeometry& geom, const DendroScalar* const* in_gfs, DendroScalar* const* rhs_gfs",
    ...          "const auto p = geom.component_offset; rhs_gfs[0][p] = in_gfs[0][p] + in_gfs[1][p];"),
    ...         ("wave_rhs_flat", "rhs_eval_flat_block",
    ...          "const BlockGeometry& geom, const DendroScalar* in_gfs, DendroScalar* rhs_gfs",
    ...          "const auto v = geom.nx*geom.ny*geom.nz; rhs_gfs[0] = in_gfs[0] + in_gfs[v];"),
    ...     )
    ...     for _name, _role, _params, _body in _function_specs:
    ...         cfc.register_CFunction(
    ...             desc="wave fixture", name=_name, params=_params, body=_body
    ...         )
    ...         roles.set_CFunction_role(_name, _role)
    ...     _types = types_h.output_types_h(
    ...         "wave", "wave", "struct wave_status { int accepted = 0; };"
    ...     )
    ...     _state = state_h.output_state_h("wave", "wave")
    ...     _header = output_solver_context_h(
    ...         "wave", "wave",
    ...         "  int initialize_scalar_vector();\n  double max_wave_rhs();\n  void apply_wave_boundary();",
    ...         "  double energy = 0.0;", "", "  double max_wave_rhs();"
    ...     )
    ...     _source = output_solver_context_cpp(
    ...         "wave", "wave", "  energy = 0.0;", "  energy = dx;", "  if (energy < 0.0) return 1;",
    ...         "int Ctx::initialize_scalar_vector() { energy = 1.0; return 0; }\nvoid Ctx::apply_wave_boundary() { energy += 1.0; }",
    ...         "double Ctx::max_wave_rhs() { return max_interior_rhs(); }", "", "",
    ...         "  flat[0] = 1.0; flat[1] = 2.0; flat[2] = 3.0;",
    ...         "double Ctx::max_wave_rhs() { return interior_max(*m_uiMesh, unzipped_rhs); }"
    ...     )
    ...     _main = generic_main.output_main_cpp(
    ...         "wave", "wave", "waveSolver", "  if (ctx.initialize_scalar_vector()) return 1;",
    ...         "  const double wave_norm = ctx.max_wave_rhs();", "    ctx.apply_wave_boundary();",
    ...         "  if (!std::isfinite(wave_norm)) return 1;", "4", "0.25*dx",
    ...         "      if (!std::isfinite(context.max_wave_rhs())) return 1;"
    ...     )
    ...     assert all(name in _state for name in ("wave_scalar", "waveU0", "waveU1"))
    ...     assert "struct wave_status" in _types
    ...     assert _header.count("double max_wave_rhs();") == 2
    ...     assert _source.count("double Ctx::max_wave_rhs()") == 2
    ...     assert "int Ctx::initialize_scalar_vector()" in _source
    ...     assert "void Ctx::apply_wave_boundary()" in _source
    ...     assert "ctx.initialize_scalar_vector()" in _main
    ...     assert "ctx.apply_wave_boundary()" in _main
    ...     assert "context.max_wave_rhs()" in _main
    ...     assert "in_gfs[0][p] + in_gfs[1][p]" in cfc.CFunction_dict["wave_rhs_block"].body
    ...     assert all(token in _source for token in ("wave_rhs", "wave_rhs_block", "wave_rhs_flat"))
    ...     assert "flat[0] = 1.0" in _source
    ...     assert all(token not in (_types + _state + _header + _source + _main)
    ...                for token in ("minkowski", "detgtrazero", "constraints_eval"))
    ... finally:
    ...     gri.glb_gridfcs_dict.clear(); gri.glb_gridfcs_dict.update(_saved_fields)
    ...     cfc.CFunction_dict.clear(); cfc.CFunction_dict.update(_saved_cfuncs)
    ...     _ = par.glb_extras_dict.pop("Dendro", None)
    ...     _ = par.glb_extras_dict.setdefault("Dendro", _saved_dendro) if _saved_dendro is not None else None
    ...     par.set_parval_from_str("Infrastructure", _saved_infrastructure)
    >>> set(gri.glb_gridfcs_dict) == set(_saved_fields)
    True
    >>> all(gri.glb_gridfcs_dict[name] is value for name, value in _saved_fields.items())
    True
    >>> set(par.glb_code_params_dict) == set(_saved_parameters)
    True
    >>> all(par.glb_code_params_dict[name] is value for name, value in _saved_parameters.items())
    True
    >>> set(cfc.CFunction_dict) == set(_saved_cfuncs)
    True
    >>> all(cfc.CFunction_dict[name] is value for name, value in _saved_cfuncs.items())
    True
    >>> set(par.glb_extras_dict) == set(_saved_extras)
    True
    >>> all(par.glb_extras_dict[name] is value for name, value in _saved_extras.items())
    True
    >>> ("Dendro" in par.glb_extras_dict) == _saved_dendro_present
    True
    >>> par.glb_extras_dict.get("Dendro") is _saved_dendro
    True
    >>> (_saved_dendro is not None and "CFunction_roles" in _saved_dendro) == _saved_roles_present
    True
    >>> _saved_dendro is None or _saved_dendro.get("CFunction_roles") is _saved_roles
    True
    >>> (_saved_dendro is not None and "CFunction_codeparameters" in _saved_dendro) == _saved_codeparameters_present
    True
    >>> _saved_dendro is None or _saved_dendro.get("CFunction_codeparameters") is _saved_codeparameters
    True
    >>> par.parval_from_str("Infrastructure") == _saved_infrastructure
    True
    >>> par.parval_from_str("fd_order") == _saved_fd_order
    True
    """
    text = (
        "#if defined(NRPY_DENDRO_STANDALONE_HOST)\n"
        + _SOURCE
        + "\n#else\n"
        + _REAL_SOURCE
        + "\n#endif\n"
    )
    for token, value in (
        ("$STANDALONE_APPLICATION_DESTRUCTOR", standalone_application_destructor),
        ("$STANDALONE_APPLICATION_POST_MESH", standalone_application_post_mesh),
        (
            "$STANDALONE_APPLICATION_STARTUP_CHECKS",
            standalone_application_startup_checks,
        ),
        (
            "$STANDALONE_APPLICATION_INITIALIZATION",
            standalone_application_initialization,
        ),
        ("$STANDALONE_APPLICATION_AFTER_RHS", standalone_application_after_rhs),
        (
            "$STANDALONE_APPLICATION_FREE_DEFINITIONS",
            standalone_application_free_definitions,
        ),
        ("$REAL_APPLICATION_INITIALIZATION", real_application_initialization),
        ("$APPLICATION_EXTERIOR_VALUES", real_application_exterior_values),
        ("$REAL_APPLICATION_AFTER_RHS", real_application_after_rhs),
    ):
        text = text.replace(token, value)
    unresolved = tuple(
        token
        for token in (
            "$STANDALONE_APPLICATION_DESTRUCTOR",
            "$STANDALONE_APPLICATION_POST_MESH",
            "$STANDALONE_APPLICATION_STARTUP_CHECKS",
            "$STANDALONE_APPLICATION_INITIALIZATION",
            "$STANDALONE_APPLICATION_AFTER_RHS",
            "$STANDALONE_APPLICATION_FREE_DEFINITIONS",
            "$REAL_APPLICATION_INITIALIZATION",
            "$APPLICATION_EXTERIOR_VALUES",
            "$REAL_APPLICATION_AFTER_RHS",
        )
        if token in text
    )
    if unresolved:
        raise ValueError(f"Application insertions were not resolved: {unresolved}")
    return BANNER + substitute_solver_identifiers(
        text,
        solver_stem,
        solver_namespace,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
