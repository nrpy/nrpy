# nrpy/infrastructures/Dendro/general_relativity/solver_context.py
"""Assemble GR application policy into the generic Dendro context."""

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

    Doctests:
    >>> substitute_application_identifiers("identity")
    'identity'
    >>> from nrpy.infrastructures.Dendro import CFunction_roles as roles
    >>> _saved = dict(roles._CFunction_roles())
    >>> roles._CFunction_roles().clear()
    >>> try:
    ...     substitute_application_identifiers("$CONSTRAINTS_EVAL")
    ... except ValueError as exc:
    ...     "constraints_eval" in str(exc)
    True
    >>> roles._CFunction_roles().update(_saved)
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
  int enforce_detgbar_equals_detghat_trAzero_all_blocks();
  double max_constraint_violation();
  int snapshot_state();
  double max_drift_from_snapshot();"""
    standalone_members = """  $SCALAR** u0 = nullptr;
  $NAMESPACE::generated::detgtrazero_status_struct last_detgtrazero_status;
  unsigned long long detgtrazero_passes = 0;"""
    standalone_free = """double observed_convergence_order(
    double base_dx, const $NAMESPACE::generated::params_struct& params);"""
    real_declarations = """  unsigned long long projection_passes = 0;
  double projection_residual = 0.0;
  int initialize();
  int pre_stage(DVec&) { return 0; }
  int post_stage(DVec&) { return 0; }
  int pre_timestep(DVec&) { return 0; }
  int post_timestep(DVec& input);
  double max_constraints();
  double max_rhs();
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


def output_solver_context_cpp(solver_stem: str, solver_namespace: str) -> str:
    """
    Emit the generic context implementation with GR operations resolved.

    :param solver_stem: Lowercase formulation stem used in emitted names.
    :param solver_namespace: Namespace containing the generated solver.
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
  }
  u0 = new $SCALAR*[ncomp];
  for (unsigned f = 0; f < ncomp; ++f) {
    u0[f] = new $SCALAR[total];
    std::memcpy(u0[f], host.in.comp[f], sizeof($SCALAR)*total);
  }
  return 0;
}
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
}
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
  }
  return 0;
}
double Ctx::max_constraint_violation() {
  $CONSTRAINTS_EVAL(host.mesh, host.in.comp, host.diag.comp);
  return max_interior_value(host.diag.comp,
                            $NAMESPACE::generated::NUM_DIAG_GFS);
}"""
    standalone_free_definitions = """namespace {

// Evaluate the generated RHS on one block at a shared physical point.
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
    initialization = """int Ctx::initialize() {
  if (!$VALIDATE(params)) fail("invalid runtime parameters");
  if (!m_uiMesh->isActive()) return 0;
  std::vector<DendroScalar*> pointers(generated::NUM_EVOL_GFS);
  unzipped.to_2d(pointers.data());
  for (const auto& b : m_uiMesh->getLocalBlockList()) {
    const auto g = block_geometry(*m_uiMesh, b, m_uiMinPt, m_uiMaxPt);
    $MINKOWSKI_INITIAL_DATA_BLOCK(g, pointers.data());
  } // END LOOP: initialize local blocks
  zip(unzipped, state);
  post_timestep(state);
  initial.assign(state.get_vec_ptr(), state.get_vec_ptr() + state.get_size());
  return 0;
} // END FUNCTION: initialize Minkowski state"""
    exterior_values = """BlockGeometry one{};
  one.nx = one.ny = one.nz = 1;
  one.dx[0] = one.dx[1] = one.dx[2] = 1.0;
  std::vector<DendroScalar*> fp(flat.size());
  for (unsigned f = 0; f < flat.size(); ++f) fp[f] = &flat[f];
  $MINKOWSKI_INITIAL_DATA_BLOCK(one, fp.data());"""
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
        )
    )


if __name__ == "__main__":
    import doctest

    raise SystemExit(doctest.testmod().failed)
