# nrpy/infrastructures/Dendro/general_relativity/solver_context.py
"""
Assemble GR application policy into the generic Dendro context.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Tuple

from nrpy.infrastructures.Dendro import CFunction_roles as roles
from nrpy.infrastructures.Dendro import solver_context as generic_context

_APPLICATION_REPLACEMENTS: Tuple[Tuple[str, str], ...] = (
    (
        "$ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_BLOCK",
        "enforce_detgbar_equals_detghat_trAzero_block",
    ),
    (
        "$ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO",
        "enforce_detgbar_equals_detghat_trAzero",
    ),
    ("$MINKOWSKI_INITIAL_DATA_BLOCK", "minkowski_initial_data_block"),
    ("$SMOOTH_PERTURBATION_BLOCK", "smooth_perturbation_block"),
    ("$CONSTRAINTS_EVAL_BLOCK", "constraints_eval_block"),
    ("$MINKOWSKI_INITIAL_DATA", "minkowski_initial_data"),
    ("$SMOOTH_PERTURBATION", "smooth_perturbation"),
    ("$CONSTRAINTS_EVAL", "constraints_eval"),
)


def substitute_application_identifiers(text: str) -> str:
    """
    Resolve only GR roles whose exact tokens occur in ``text``.

    :param text: Generated text that may request registered GR roles.
    :return: Text with every requested GR role resolved.

    """
    for token, role in _APPLICATION_REPLACEMENTS:
        if token in text:
            text = text.replace(token, roles.CFunction_name_for_role(role))
    return text


def output_solver_context_h(solver_stem: str, solver_namespace: str) -> str:
    """
    Emit the generic context with registered GR operations resolved.

    :param solver_stem: Lowercase formulation stem used in emitted names.
    :param solver_namespace: Namespace containing the generated solver.
    :return: Complete generated context header.
    """
    standalone_declarations = """  // Fill the state through generated GR initial data.
  int minkowski_initial_data();
  // Add the smooth GR qualification perturbation to every evolved component.
  int perturb_state();
  /**
   * Enforce the conformal determinant and trace-free constraints on all blocks.
   *
   * @return 0 when every point is projected; 1 if any point is rejected.
   */
  int enforce_detgbar_equals_detghat_trAzero_all_blocks();
  double max_constraint_violation();
  /**
   * Replace the owned reference snapshot with the current evolved state.
   *
   * @return 0 after allocating and copying the snapshot.
   */
  int snapshot_state();
  /**
   * Measure drift between the evolved state and the owned snapshot.
   *
   * @return Maximum absolute component drift over every local point.
   */
  double max_drift_from_snapshot();"""
    standalone_members = """  $SCALAR** u0 = nullptr;
  $NAMESPACE::generated::detgtrazero_status_struct last_detgtrazero_status;
  unsigned long long detgtrazero_passes = 0;"""
    standalone_free = """/**
 * Estimate the generated RHS convergence order at one shared physical point.
 *
 * @param base_dx Coarsest grid spacing.
 * @param[in] params Generated parameters forwarded to the RHS kernel.
 * @return Observed base-two convergence order, or -1 if differences vanish.
 */
double observed_convergence_order(
    double base_dx, const $NAMESPACE::generated::params_struct& params);"""
    real_declarations = """  unsigned long long projection_passes = 0;
  unsigned long long dendrogr_initial_data_calls = 0;
  double projection_residual = 0.0;
  /**
   * Initialize the real-host state from Dendro-GR Minkowski data.
   *
   * @return 0 on success or an inactive rank; invalid data aborts MPI_COMM_WORLD.
   */
  int initialize();
  int pre_stage(DVec&) { return 0; }
  int post_stage(DVec&) { return 0; }
  int pre_timestep(DVec&) { return 0; }
  /**
   * Project algebraic constraints after a real-host timestep.
   *
   * @param[in,out] input Packed evolved state updated by the projection.
   * @return 0 on success or an inactive rank; invalid data aborts MPI_COMM_WORLD.
   */
  int post_timestep(DVec& input);
  /**
   * Evaluate and reduce the generated constraint diagnostics.
   *
   * @return World maximum absolute constraint value.
   */
  double max_constraints();
  double max_rhs();
  /**
   * Compare the current real-host state with its initial state.
   *
   * @return World maximum absolute evolved-state drift.
   */
  double max_drift();
 private:
  std::vector<DendroScalar> initial;
 public:"""
    return substitute_application_identifiers(
        generic_context.output_solver_context_h(
            solver_stem,
            solver_namespace,
            standalone_declarations,
            standalone_members,
            standalone_free,
            real_declarations,
        )
    )


def output_solver_context_cpp(
    solver_stem: str, solver_namespace: str, enable_fCCZ4: bool = False
) -> str:
    """
    Emit the generic context implementation with GR operations resolved.

    :param solver_stem: Lowercase formulation stem used in emitted names.
    :param solver_namespace: Namespace containing the generated solver.
    :param enable_fCCZ4: Whether the generated state extends BSSN with Theta.
    :return: Complete generated context implementation.
    """
    standalone_destructor = """  if (u0 != nullptr) {
    for (unsigned f = 0; f < $NAMESPACE::generated::NUM_EVOL_GFS; ++f)
      delete[] u0[f];
    delete[] u0;
    u0 = nullptr;
  }  // END IF: snapshot buffer allocated"""
    standalone_post_mesh = """  // The perturbation wavelength is a length and
  // is fixed only after the standalone mesh spacing is known.
  params.smooth_perturbation_wavelength = 0.618 * extent * dx;"""
    standalone_startup = """  const unsigned fd_order = $NAMESPACE::generated::FD_ORDER;
  if (fd_order != 2 && fd_order != 4 && fd_order != 6) {
    std::fprintf(stderr,
                 "ERROR: unsupported generated FD order %u; the qualified set "
                 "is 2, 4, 6\\n", fd_order);
    return 1;
  }  // END IF: unsupported generated FD order"""
    standalone_initialization = """int Ctx::minkowski_initial_data() {
  $MINKOWSKI_INITIAL_DATA(host.mesh, host.in.comp);
  return 0;
}  // END FUNCTION: Ctx::minkowski_initial_data

int Ctx::perturb_state() {
  $SMOOTH_PERTURBATION(host.mesh, host.in.comp,
      static_cast<$SCALAR>(params.smooth_perturbation_amplitude),
      static_cast<$SCALAR>(params.smooth_perturbation_wavelength));
  return 0;
}  // END FUNCTION: Ctx::perturb_state"""
    standalone_after_rhs = """int Ctx::snapshot_state() {
  const unsigned ncomp = $NAMESPACE::generated::NUM_EVOL_GFS;
  const unsigned nb = host.in.num_blocks;
  const unsigned extent = host.mesh.geom[0].nx;
  const std::size_t total = static_cast<std::size_t>(nb)*extent*extent*extent;
  if (u0 != nullptr) {
    for (unsigned f = 0; f < ncomp; ++f) delete[] u0[f];
    delete[] u0;
    u0 = nullptr;
  }  // END IF: replace snapshot buffer
  u0 = new $SCALAR*[ncomp];
  for (unsigned f = 0; f < ncomp; ++f) {
    u0[f] = new $SCALAR[total];
    std::memcpy(u0[f], host.in.comp[f], sizeof($SCALAR)*total);
  }  // END LOOP: copy evolved snapshot
  return 0;
}  // END FUNCTION: Ctx::snapshot_state
double Ctx::max_drift_from_snapshot() {
  const unsigned ncomp = $NAMESPACE::generated::NUM_EVOL_GFS;
  const unsigned nb = host.in.num_blocks;
  const unsigned extent = host.mesh.geom[0].nx;
  const std::size_t total = static_cast<std::size_t>(nb)*extent*extent*extent;
  double worst = 0.0;
  for (unsigned f = 0; f < ncomp; ++f)
    for (std::size_t cell = 0; cell < total; ++cell)
      worst = std::max(worst, std::fabs(static_cast<double>(host.in.comp[f][cell]) -
                                        static_cast<double>(u0[f][cell])));
  return worst;
}  // END FUNCTION: Ctx::max_drift_from_snapshot
int Ctx::enforce_detgbar_equals_detghat_trAzero_all_blocks() {
  last_detgtrazero_status = $NAMESPACE::generated::detgtrazero_status_struct{};
  $ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO(
      host.mesh, host.in.comp, &last_detgtrazero_status);
  ++detgtrazero_passes;
  if (last_detgtrazero_status.failed_points != 0) {
    std::fprintf(stderr,
        "ERROR: constraint enforcement refused %llu point(s); first at index "
        "%lld, first nonfinite field index %d\\n",
        last_detgtrazero_status.failed_points,
        last_detgtrazero_status.first_failing_index,
        last_detgtrazero_status.first_failing_field);
    return 1;
  }  // END IF: constraint enforcement failed
  return 0;
}  // END FUNCTION: Ctx::enforce_detgbar_equals_detghat_trAzero_all_blocks
double Ctx::max_constraint_violation() {
  $CONSTRAINTS_EVAL(host.mesh, host.in.comp, host.diag.comp);
  return max_interior_value(host.diag.comp,
                            $NAMESPACE::generated::NUM_DIAG_GFS);
}  // END FUNCTION: Ctx::max_constraint_violation"""
    standalone_free_definitions = """namespace {

/**
 * Evaluate the generated RHS on one block at a shared physical point.
 *
 * @param extent Number of points on each block axis.
 * @param dx Uniform grid spacing.
 * @param sample_index Index of the sampled point on each axis.
 * @param[in] params Generated parameters forwarded to the RHS kernel.
 * @param[out] out Sampled RHS values in generated evolved-variable order.
 */
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
  block_geometry_struct g;
  g.nx = g.ny = g.nz = static_cast<unsigned>(extent);
  g.padding = $NAMESPACE::generated::REQUIRED_PADDING;
  g.component_offset = 0;
  g.pmin_padded[0] = g.pmin_padded[1] = g.pmin_padded[2] = 0.0;
  g.dx[0] = g.dx[1] = g.dx[2] = dx;
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
    const double difference = std::fabs(a[k] - b[k]);
    if (difference > worst) worst = difference;
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
  const int coarse = 2 * pad + 9;
  const int centre = (coarse - 1) / 2;
  std::vector<double> r_h;
  std::vector<double> r_h2;
  std::vector<double> r_h4;
  sample_rhs_at_centre(coarse, base_dx, centre, params, r_h);
  sample_rhs_at_centre(2 * coarse - 1, base_dx / 2.0, 2 * centre, params, r_h2);
  sample_rhs_at_centre(4 * coarse - 3, base_dx / 4.0, 4 * centre, params, r_h4);
  const double d1 = max_abs_difference(r_h, r_h2);
  const double d2 = max_abs_difference(r_h2, r_h4);
  if (!(d1 > 0.0) || !(d2 > 0.0)) return -1.0;
  return std::log2(d1 / d2);
}  // END FUNCTION: observed_convergence_order"""
    cf_from_dendrogr_chi = (
        "native[::bssn::VAR::U_CHI]"
        if enable_fCCZ4
        else "std::sqrt(native[::bssn::VAR::U_CHI])"
    )
    initialization = f"""namespace {{

/**
 * Adapt one Dendro-GR Minkowski point to the generated state layout.
 *
 * Dendro-GR supplies the physical values; this function changes only ordering and
 * storage.
 * Fields absent from Dendro-GR's BSSN state start at generated asymptotic values.
 *
 * @param x Physical x coordinate.
 * @param y Physical y coordinate.
 * @param z Physical z coordinate.
 * @param[out] output Generated component arrays receiving the adapted state.
 * @param cell Flat destination index in every component array.
 */
void dendrogr_minkowski_point(double x, double y, double z,
                              DendroScalar* const* output,
                              std::size_t cell) {{
  std::array<double, ::bssn::BSSN_NUM_VARS> native{{}};
  ::bssn::minkowskiInitialData(x, y, z, native.data());
  for (unsigned f = 0; f < generated::NUM_EVOL_GFS; ++f)
    output[f][cell] = generated::EVOL_GF_F_INFINITY[f];
  const auto store = [&](generated::EvolVar target, ::bssn::VAR source,
                         double offset = 0.0) {{
    output[generated::to_index(target)][cell] = native[source] + offset;
  }};
  store(generated::EvolVar::aDD00, ::bssn::VAR::U_SYMAT0);
  store(generated::EvolVar::aDD01, ::bssn::VAR::U_SYMAT1);
  store(generated::EvolVar::aDD02, ::bssn::VAR::U_SYMAT2);
  store(generated::EvolVar::aDD11, ::bssn::VAR::U_SYMAT3);
  store(generated::EvolVar::aDD12, ::bssn::VAR::U_SYMAT4);
  store(generated::EvolVar::aDD22, ::bssn::VAR::U_SYMAT5);
  store(generated::EvolVar::alpha, ::bssn::VAR::U_ALPHA);
  store(generated::EvolVar::betU0, ::bssn::VAR::U_B0);
  store(generated::EvolVar::betU1, ::bssn::VAR::U_B1);
  store(generated::EvolVar::betU2, ::bssn::VAR::U_B2);
  output[generated::to_index(generated::EvolVar::cf)][cell] =
      {cf_from_dendrogr_chi};
  store(generated::EvolVar::hDD00, ::bssn::VAR::U_SYMGT0, -1.0);
  store(generated::EvolVar::hDD01, ::bssn::VAR::U_SYMGT1);
  store(generated::EvolVar::hDD02, ::bssn::VAR::U_SYMGT2);
  store(generated::EvolVar::hDD11, ::bssn::VAR::U_SYMGT3, -1.0);
  store(generated::EvolVar::hDD12, ::bssn::VAR::U_SYMGT4);
  store(generated::EvolVar::hDD22, ::bssn::VAR::U_SYMGT5, -1.0);
  store(generated::EvolVar::lambdaU0, ::bssn::VAR::U_GT0);
  store(generated::EvolVar::lambdaU1, ::bssn::VAR::U_GT1);
  store(generated::EvolVar::lambdaU2, ::bssn::VAR::U_GT2);
  store(generated::EvolVar::trK, ::bssn::VAR::U_K);
  store(generated::EvolVar::vetU0, ::bssn::VAR::U_BETA0);
  store(generated::EvolVar::vetU1, ::bssn::VAR::U_BETA1);
  store(generated::EvolVar::vetU2, ::bssn::VAR::U_BETA2);
}} // END FUNCTION: dendrogr_minkowski_point

// clang-format off
}}  // END NAMESPACE: internal linkage
// clang-format on

int Ctx::initialize() {{
  if (!$VALIDATE(params)) fail("invalid runtime parameters");
  if (!m_uiMesh->isActive()) return 0;
  std::vector<DendroScalar*> pointers(generated::NUM_EVOL_GFS);
  unzipped.to_2d(pointers.data());
  for (const auto& b : m_uiMesh->getLocalBlockList()) {{
    const auto g = block_geometry(*m_uiMesh, b, m_uiMinPt, m_uiMaxPt);
    for (unsigned k = 0; k < g.nz; ++k)
      for (unsigned j = 0; j < g.ny; ++j)
        for (unsigned i = 0; i < g.nx; ++i) {{
          const std::size_t cell =
              g.component_offset + i + std::size_t(g.nx) *
                                       (j + std::size_t(g.ny) * k);
          dendrogr_minkowski_point(
              g.pmin_padded[0] + i * g.dx[0],
              g.pmin_padded[1] + j * g.dx[1],
              g.pmin_padded[2] + k * g.dx[2], pointers.data(), cell);
          ++dendrogr_initial_data_calls;
        }} // END LOOP: initialize one padded Dendro point
  }} // END LOOP: initialize local blocks
  zip(unzipped, state);
  post_timestep(state);
  initial.assign(state.get_vec_ptr(), state.get_vec_ptr() + state.get_size());
  return 0;
}} // END FUNCTION: Ctx::initialize"""
    exterior_values = """// Exterior values use the same Dendro-GR initial-data
  // routine as the evolved initial state.
  std::vector<DendroScalar*> fp(flat.size());
  for (unsigned f = 0; f < flat.size(); ++f) fp[f] = &flat[f];
  dendrogr_minkowski_point(0.0, 0.0, 0.0, fp.data(), 0);"""
    after_rhs = """int Ctx::post_timestep(DVec& input) {
  if (!m_uiMesh->isActive()) return 0;
  require_finite(input);
  unzip(input, unzipped, 1);
  std::vector<DendroScalar*> pointers(generated::NUM_EVOL_GFS);
  unzipped.to_2d(pointers.data());
  for (const auto& b : m_uiMesh->getLocalBlockList()) {
    const auto g = block_geometry(*m_uiMesh, b, m_uiMinPt, m_uiMaxPt);
    generated::detgtrazero_status_struct status{};
    $ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_BLOCK(g, pointers.data(), &status);
    if (status.failed_points || status.nonfinite_points) fail("projection failed");
    projection_residual = std::max(projection_residual,
        std::max(status.max_abs_det_minus_one, status.max_abs_trace_residual));
  } // END LOOP: project local blocks
  zip(unzipped, input);
  ++projection_passes;
  return 0;
} // END FUNCTION: project evolved stage state
double Ctx::max_constraints() {
  double local = 0.0;
  if (m_uiMesh->isActive()) {
    unzip(state, unzipped, 1);
    fill_exterior();
    std::vector<DendroScalar*> input(generated::NUM_EVOL_GFS), output(generated::NUM_DIAG_GFS);
    unzipped.to_2d(input.data());
    diagnostics.to_2d(output.data());
    for (const auto& b : m_uiMesh->getLocalBlockList()) {
      const auto g = block_geometry(*m_uiMesh, b, m_uiMinPt, m_uiMaxPt);
      $CONSTRAINTS_EVAL_BLOCK(g, input.data(), output.data());
    } // END LOOP: evaluate block constraints
    local = interior_max(*m_uiMesh, diagnostics);
  } // END IF: evaluate active rank diagnostics
  return maximum(local);
} // END FUNCTION: reduce constraint maximum
double Ctx::max_rhs() {
  DVec output;
  output.create_vector(m_uiMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                       ot::DVEC_LOC::HOST, generated::NUM_EVOL_GFS, true);
  rhs(&state, &output, 1, m_uiTinfo._m_uiT);
  const double local = interior_max(*m_uiMesh, unzipped_rhs);
  output.destroy_vector();
  return maximum(local);
} // END FUNCTION: measure current RHS
double Ctx::max_drift() {
  double local = 0.0;
  require_finite(state);
  if (m_uiMesh->isActive()) {
    const unsigned stride = m_uiMesh->getDegOfFreedom();
    for (unsigned f = 0; f < generated::NUM_EVOL_GFS; ++f)
      for (unsigned i = m_uiMesh->getNodeLocalBegin();
           i < m_uiMesh->getNodeLocalEnd(); ++i)
        local = std::max(local, std::abs(
            state.get_vec_ptr()[std::size_t(f)*stride+i] -
            initial[std::size_t(f)*stride+i]));
  } // END IF: compare owned evolved values
  return maximum(local);
} // END FUNCTION: measure evolved state drift"""
    return substitute_application_identifiers(
        generic_context.output_solver_context_cpp(
            solver_stem,
            solver_namespace,
            standalone_destructor,
            standalone_post_mesh,
            standalone_startup,
            standalone_initialization,
            standalone_after_rhs,
            standalone_free_definitions,
            initialization,
            exterior_values,
            after_rhs,
            '#include "grUtils.h"\n#include <array>',
        )
    )


if __name__ == "__main__":
    import doctest

    raise SystemExit(doctest.testmod().failed)
