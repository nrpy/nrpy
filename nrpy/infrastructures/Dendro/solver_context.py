# nrpy/infrastructures/Dendro/solver_context.py
"""
Emit the standalone-host runtime context for a generated Dendro solver.

The context allocates the generated-count vectors, runs the startup checks and
invokes the registered generated CFunctions.  It carries no formulation
content: every count, name and constant comes from a generated header, and
every kernel name is read back from the Dendro role registry.

SCOPE OF EVIDENCE.  This vehicle establishes lifecycle plumbing -- allocation
with generated counts, the call path, the block and flat-block entry points
agreeing, per-rank decomposition, and the Minkowski fixed point.  It cannot
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
  // Fill the state through the generated initial-data CFunction.
  int minkowski_initial_data();
  // Add a smooth, analytic, spatially varying perturbation to every evolved
  // component over the whole padded block, so the interior RHS exercises the
  // generated derivative stencils.  No field is named: the same continuum
  // profile is applied to every component of the generated state.
  int perturb_state();
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
  // Run the generated constraint enforcement over every local block.  Returns
  // nonzero when the generated status reports a refused point; the enforcement
  // itself never calls exit().
  int enforce_detgbar_equals_detghat_trAzero_all_blocks();
  // Maximum |diagnostic| over every DIAG component and local block, after
  // recomputing the generated constraint diagnostics.  No diagnostic is named:
  // the reduction runs over the generated count.
  double max_constraint_violation();
  // Resolve host-supplied exact NRPy names (output or refinement selection).
  // An unknown name is fatal: this prints every valid generated name and
  // returns nonzero.  No name is written down here.
  int select_variables(const char* const* names, unsigned count);
  // Copy the current state (bitwise snapshot for the drift gate).
  int snapshot_state();
  // Maximum |u - u_snapshot| over all components and blocks.
  double max_drift_from_snapshot();

  // Host mesh and EVOL vectors (in / rhs / out), standalone-host lifecycle.
  // Value-initialized: standalone_host::DVector is an aggregate with no default member
  // initializers, and ~Ctx frees all four vectors unconditionally.  When
  // initialize_mesh rejects its inputs it returns before assigning them, so
  // without the braces the destructor would free indeterminate pointers.
  standalone_host::Ctx host{};
  // Generated runtime parameter table, owned by the context.
  $NAMESPACE::generated::params_struct params;
  // Bitwise snapshot of the state taken before evolution.
  $SCALAR** u0 = nullptr;
  // Status of the most recent enforcement pass.
  $NAMESPACE::generated::detgtrazero_status_struct last_detgtrazero_status;
  // Number of completed enforcement passes.  A flat state looks the same
  // whether or not it was projected, so a count the host can check is the
  // evidence that the hooks ran as configured.
  unsigned long long detgtrazero_passes = 0;
};  // END CLASS: Ctx

// Observed convergence order of the generated RHS under grid refinement: the
// same smooth analytic state is evaluated at spacings h, h/2 and h/4 on a
// single block, and the Richardson ratio of the successive differences at
// coincident physical points is returned.  A kernel whose stencils are absent,
// mis-shaped, or applied at the wrong offsets does not exhibit the requested
// order.  Returns a negative value when the differences are too small to form
// a ratio.
double observed_convergence_order(
    double base_dx, const $NAMESPACE::generated::params_struct& params);

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
  if (u0 != nullptr) {
    for (unsigned f = 0; f < $NAMESPACE::generated::NUM_EVOL_GFS; ++f) {
      delete[] u0[f];
    }
    delete[] u0;
    u0 = nullptr;
  }  // END IF: snapshot buffer allocated
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
  // One default the generator cannot know: the perturbation wavelength is a
  // length, so it is meaningful only against the grid the host just built.
  // 0.618 of a block is deliberately incommensurate with the per-rank stride
  // (rank * n_blocks * extent * dx); a wavelength equal to the stride would
  // make every rank see an identical field, and the multi-rank gates would
  // prove nothing beyond a repeated serial run.
  params.smooth_perturbation_wavelength = 0.618 * extent * dx;
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
  const unsigned fd_order = $NAMESPACE::generated::FD_ORDER;
  if (fd_order != 2 && fd_order != 4 && fd_order != 6) {
    std::fprintf(stderr,
                 "ERROR: unsupported generated FD order %u; the qualified set "
                 "is 2, 4, 6\\n",
                 fd_order);
    return 1;
  }  // END IF: unsupported generated FD order
  if (!$VALIDATE(params)) {
    std::fprintf(stderr, "ERROR: parameter validation failed\\n");
    return 1;
  }  // END IF: parameter validation failed
  // Rank-0 only, as every other diagnostic in the entry point is: an N-rank
  // run printing N identical parameter tables buries the gate output.
  if (rank == 0) $PRINT_EFFECTIVE(params);
  return 0;
}  // END FUNCTION: Ctx::startup_checks

int Ctx::minkowski_initial_data() {
  $MINKOWSKI_INITIAL_DATA(host.mesh, host.in.comp);
  return 0;
}  // END FUNCTION: Ctx::minkowski_initial_data

int Ctx::perturb_state() {
  // The analytic profile and the loop that applies it are NRPy-authored and
  // registered; the host only supplies the state and the two scalars.  Those
  // two come from params, so the values printed as effective are the values
  // the run perturbs with.
  $SMOOTH_PERTURBATION(host.mesh, host.in.comp,
                static_cast<$SCALAR>(params.smooth_perturbation_amplitude),
                static_cast<$SCALAR>(params.smooth_perturbation_wavelength));
  return 0;
}  // END FUNCTION: Ctx::perturb_state

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

int Ctx::snapshot_state() {
  const unsigned ncomp = $NAMESPACE::generated::NUM_EVOL_GFS;
  const unsigned nb = host.in.num_blocks;
  const unsigned extent = host.mesh.geom[0].nx;
  const std::size_t total =
      static_cast<std::size_t>(nb) * extent * extent * extent;
  if (u0 != nullptr) {
    for (unsigned f = 0; f < ncomp; ++f) delete[] u0[f];
    delete[] u0;
    u0 = nullptr;
  }  // END IF: snapshot buffer allocated
  u0 = new $SCALAR*[ncomp];
  for (unsigned f = 0; f < ncomp; ++f) {
    u0[f] = new $SCALAR[total];
    std::memcpy(u0[f], host.in.comp[f], sizeof($SCALAR) * total);
  }  // END LOOP: for f over evolved components
  return 0;
}  // END FUNCTION: Ctx::snapshot_state

double Ctx::max_drift_from_snapshot() {
  const unsigned ncomp = $NAMESPACE::generated::NUM_EVOL_GFS;
  const unsigned nb = host.in.num_blocks;
  const unsigned extent = host.mesh.geom[0].nx;
  const std::size_t total =
      static_cast<std::size_t>(nb) * extent * extent * extent;
  double worst = 0.0;
  for (unsigned f = 0; f < ncomp; ++f) {
    for (std::size_t cell = 0; cell < total; ++cell) {
      const double d = std::fabs(static_cast<double>(host.in.comp[f][cell]) -
                                 static_cast<double>(u0[f][cell]));
      if (d > worst) worst = d;
    }  // END LOOP: for cell over local-block cells
  }  // END LOOP: for f over evolved components
  return worst;
}  // END FUNCTION: Ctx::max_drift_from_snapshot

int Ctx::enforce_detgbar_equals_detghat_trAzero_all_blocks() {
  // The enforcement is scheduled after initial-data construction and after
  // every accepted timestep.  The generated kernel never calls exit(); it
  // reports a structured status, and the host decides.
  last_detgtrazero_status = $NAMESPACE::generated::detgtrazero_status_struct{};
  $ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO(host.mesh, host.in.comp, &last_detgtrazero_status);
  ++detgtrazero_passes;
  if (last_detgtrazero_status.failed_points != 0) {
    std::fprintf(stderr,
                 "ERROR: constraint enforcement refused %llu point(s); first at index "
                 "%lld, first nonfinite field index %d\\n",
                 last_detgtrazero_status.failed_points,
                 last_detgtrazero_status.first_failing_index,
                 last_detgtrazero_status.first_failing_field);
    return 1;
  }  // END IF: enforcement refused a point
  return 0;
}  // END FUNCTION: Ctx::enforce_detgbar_equals_detghat_trAzero_all_blocks

double Ctx::max_constraint_violation() {
  // Recompute the diagnostics from the current evolved state and reduce over
  // every generated DIAG component.  No diagnostic is named here: the count
  // and the ordering are generated.
  $CONSTRAINTS_EVAL(host.mesh, host.in.comp, host.diag.comp);
  return max_interior_value(host.diag.comp,
                            $NAMESPACE::generated::NUM_DIAG_GFS);
}  // END FUNCTION: Ctx::max_constraint_violation

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

namespace {

// Evaluate the generated RHS on one block of the given resolution, filled with
// the same continuum profile, and return the interior values sampled at the
// coincident physical point (the block centre) for every component.
void sample_rhs_at_centre(
    int extent, double dx, int sample_index,
    const $NAMESPACE::generated::params_struct& params,
    std::vector<double>& out) {
  const unsigned ncomp = $NAMESPACE::generated::NUM_EVOL_GFS;
  const std::size_t vol = static_cast<std::size_t>(extent) * extent * extent;
  std::vector<std::vector<$SCALAR>> state(
      ncomp, std::vector<$SCALAR>(vol, 0.0));
  std::vector<std::vector<$SCALAR>> rhs(
      ncomp, std::vector<$SCALAR>(vol, 0.0));
  BlockGeometry g;
  g.nx = g.ny = g.nz = static_cast<unsigned>(extent);
  g.padding = $NAMESPACE::generated::REQUIRED_PADDING;
  g.component_offset = 0;
  g.pmin_padded[0] = g.pmin_padded[1] = g.pmin_padded[2] = 0.0;
  g.dx[0] = g.dx[1] = g.dx[2] = dx;
  g.boundary_flags = 0;
  std::vector<$SCALAR*> state_ptr(ncomp);
  std::vector<const $SCALAR*> state_cptr(ncomp);
  std::vector<$SCALAR*> rhs_ptr(ncomp);
  for (unsigned f = 0; f < ncomp; ++f) {
    state_ptr[f] = state[f].data();
    state_cptr[f] = state[f].data();
    rhs_ptr[f] = rhs[f].data();
  }  // END LOOP: for f over evolved components
  $MINKOWSKI_INITIAL_DATA_BLOCK(g, state_ptr.data());
  $SMOOTH_PERTURBATION_BLOCK(
      g, state_ptr.data(),
      static_cast<$SCALAR>(params.smooth_perturbation_amplitude),
      static_cast<$SCALAR>(params.smooth_perturbation_wavelength));
  $RHS_EVAL_BLOCK(g, state_cptr.data(), rhs_ptr.data()$RHS_EVAL_BLOCK_TAIL);
  // The caller passes the index of the shared physical point: at spacing
  // dx/2^k that point is index sample_index*2^k, so the three grids sample the
  // same location and the Richardson ratio is meaningful.
  const unsigned sample = static_cast<unsigned>(sample_index);
  out.assign(ncomp, 0.0);
  for (unsigned f = 0; f < ncomp; ++f) {
    out[f] = static_cast<double>(
        rhs[f][sample + g.nx * (sample + g.ny * sample)]);
  }  // END LOOP: for f over evolved components
}  // END FUNCTION: sample_rhs_at_centre

double max_abs_difference(const std::vector<double>& a,
                          const std::vector<double>& b) {
  double worst = 0.0;
  for (std::size_t k = 0; k < a.size() && k < b.size(); ++k) {
    const double d = std::fabs(a[k] - b[k]);
    if (d > worst) worst = d;
  }  // END LOOP: for k over paired values
  return worst;
}  // END FUNCTION: max_abs_difference

// clang-format off
}  // END NAMESPACE: internal linkage
// clang-format on

double observed_convergence_order(
    double base_dx,
    const $NAMESPACE::generated::params_struct& params) {
  const int pad = static_cast<int>($NAMESPACE::generated::REQUIRED_PADDING);
  const int coarse = 2 * pad + 9;  // odd, so the centre is a grid point
  const int centre = (coarse - 1) / 2;
  std::vector<double> r_h;
  std::vector<double> r_h2;
  std::vector<double> r_h4;
  sample_rhs_at_centre(coarse, base_dx, centre, params,
                       r_h);
  sample_rhs_at_centre(2 * coarse - 1, base_dx / 2.0, 2 * centre, params, r_h2);
  sample_rhs_at_centre(4 * coarse - 3, base_dx / 4.0, 4 * centre, params, r_h4);
  const double d1 = max_abs_difference(r_h, r_h2);
  const double d2 = max_abs_difference(r_h2, r_h4);
  if (!(d1 > 0.0) || !(d2 > 0.0)) return -1.0;
  return std::log2(d1 / d2);
}  // END FUNCTION: observed_convergence_order

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
    """
    rhs_block = roles.CFunction_name_for_role("rhs_eval_block")
    replacements = (
        (
            "$PARAMS_STRUCT_SET_TO_DEFAULT",
            f"{solver_stem}_params_struct_set_to_default",
        ),
        ("$VALIDATE", f"{solver_stem}_params_validate"),
        ("$PRINT_EFFECTIVE", f"{solver_stem}_params_print_effective"),
        (
            "$MINKOWSKI_INITIAL_DATA_BLOCK",
            roles.CFunction_name_for_role("minkowski_initial_data_block"),
        ),
        (
            "$MINKOWSKI_INITIAL_DATA",
            roles.CFunction_name_for_role("minkowski_initial_data"),
        ),
        (
            "$SMOOTH_PERTURBATION_BLOCK",
            roles.CFunction_name_for_role("smooth_perturbation_block"),
        ),
        ("$SMOOTH_PERTURBATION", roles.CFunction_name_for_role("smooth_perturbation")),
        ("$RHS_EVAL_FLAT_BLOCK", roles.CFunction_name_for_role("rhs_eval_flat_block")),
        ("$RHS_EVAL_BLOCK_TAIL", _codeparameter_tail(rhs_block, "params")),
        ("$RHS_EVAL_BLOCK", rhs_block),
        (
            "$RHS_EVAL_TAIL",
            _codeparameter_tail(roles.CFunction_name_for_role("rhs_eval"), "params"),
        ),
        ("$RHS_EVAL", roles.CFunction_name_for_role("rhs_eval")),
        (
            "$ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_BLOCK",
            roles.CFunction_name_for_role(
                "enforce_detgbar_equals_detghat_trAzero_block"
            ),
        ),
        (
            "$ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO",
            roles.CFunction_name_for_role("enforce_detgbar_equals_detghat_trAzero"),
        ),
        (
            "$CONSTRAINTS_EVAL_BLOCK",
            roles.CFunction_name_for_role("constraints_eval_block"),
        ),
        ("$CONSTRAINTS_EVAL", roles.CFunction_name_for_role("constraints_eval")),
        ("$NAMESPACE", solver_namespace),
        ("$STEM", solver_stem),
        ("$SCALAR", gri.DENDRO_SCALAR_TYPE),
    )
    for placeholder, value in replacements:
        text = text.replace(placeholder, value)
    return text


def output_solver_context_h(solver_stem: str, solver_namespace: str) -> str:
    """
    Emit the standalone-host context header.

    :param solver_stem: Lowercase formulation stem for emitted file names.
    :param solver_namespace: Solver namespace.
    :return: The complete C++ header text.

    Doctests:
    >>> from nrpy.infrastructures.Dendro.clang_format_guards import (
    ...     unguarded_end_namespace_markers,
    ... )
    >>> unguarded_end_namespace_markers(_HEADER)
    []
    >>> _HEADER.count("}  // END NAMESPACE:")
    1
    >>> unguarded_end_namespace_markers(_SOURCE)
    []
    >>> _SOURCE.count("}  // END NAMESPACE:")
    3
    """
    opening, closing = header_guard(f"{solver_stem}Ctx.h")
    body = substitute_solver_identifiers(_HEADER, solver_stem, solver_namespace)
    return BANNER + f"{opening}\n\n" + body + f"\n{closing}\n"


def output_solver_context_cpp(solver_stem: str, solver_namespace: str) -> str:
    """
    Emit the standalone-host context implementation.

    :param solver_stem: Lowercase formulation stem for emitted file names.
    :param solver_namespace: Solver namespace.
    :return: The complete C++ source text.
    """
    return BANNER + substitute_solver_identifiers(
        _SOURCE, solver_stem, solver_namespace
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
