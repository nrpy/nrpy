# nrpy/infrastructures/Dendro/general_relativity/self_tests_cpp.py
"""
Assemble GR scientific sections into the generated Dendro self tests.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import math
from typing import Any, Dict, List, NamedTuple, Tuple

import sympy as sp
from mpmath import mp  # type: ignore[import-untyped]

import nrpy.grid as gri
import nrpy.params as par
from nrpy.finite_difference import (
    compute_fdcoeffs_fdstencl,
    extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars,
)
from nrpy.infrastructures.Dendro import CFunction_roles as roles
from nrpy.infrastructures.Dendro import gridfunction_name_decorations as gf_names
from nrpy.infrastructures.Dendro import self_tests_cpp as generic_tests
from nrpy.infrastructures.Dendro import solver_context as generic_context
from nrpy.infrastructures.Dendro.general_relativity.rhs_eval import RHSBuild
from nrpy.infrastructures.Dendro.general_relativity.solver_context import (
    substitute_application_identifiers,
)
from nrpy.validate_expressions.validate_expressions import (
    inject_mpfs_into_cse_expression,
)

SampleKey = Tuple[str, Tuple[int, int, int]]
SampleMap = Dict[SampleKey, float]


class ReferenceValue(NamedTuple):
    """
    Record one high-precision result and its pre-runtime roundoff allowance.

    :param value: Reference rounded once to the generated double type.
    :param bound: Componentwise double-evaluation error allowance.
    :param precision_delta: Difference between the 80- and 100-digit results.
    :param operation_count: Operations in the expression's actual CSE DAG.
    :param evaluation_scale: Largest stencil sum, input, or CSE temporary.
    """

    value: float
    bound: float
    precision_delta: float
    operation_count: int
    evaluation_scale: float


# Scientific names are explicit application input to CMake.  The two fixture
# sections exercise formulation-neutral lowering through the same executable.
def test_sections() -> Tuple[str, ...]:
    """
    Return meaningful test sections for the configured FD order.

    :return: Ordered CTest section names for the current generation profile.
    """
    sections: Tuple[str, ...] = (
        "state",
        "params",
        "padding",
        "offsets",
        "upwind",
        "rhs",
        "init",
        "names",
        "detgtrazero",
        "constraints",
    )
    if int(par.parval_from_str("fd_order")) == 4:
        sections += ("address_values", "parameter_forwarding", "gr_nonflat_reference")
    else:
        sections += ("address_values", "parameter_forwarding")
    return sections


_GR_SCIENTIFIC_TESTS = r"""
struct GRTestBlock {
  explicit GRTestBlock(unsigned blocks = 1)
      : extent(2 * $NAMESPACE::generated::REQUIRED_PADDING + 3),
        vol(static_cast<std::size_t>(extent) * extent * extent),
        state(NUM_EVOL_GFS, std::vector<$SCALAR>(vol * blocks, $SCALAR{0})),
        rhs(NUM_EVOL_GFS, std::vector<$SCALAR>(vol * blocks, $SCALAR{0})) {
    geometry.nx = geometry.ny = geometry.nz = extent;
    geometry.padding = $NAMESPACE::generated::REQUIRED_PADDING;
    geometry.component_offset = 0;
    geometry.pmin_padded[0] = geometry.pmin_padded[1] =
        geometry.pmin_padded[2] = 0.0;
    geometry.dx[0] = geometry.dx[1] = geometry.dx[2] = 0.1;
  }  // END FUNCTION: GRTestBlock::GRTestBlock
  std::size_t index(unsigned i, unsigned j, unsigned k) const {
    return geometry.component_offset + i + extent * (j + extent * k);
  }
  std::vector<$SCALAR*> state_pointers() {
    std::vector<$SCALAR*> result(NUM_EVOL_GFS);
    for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
      result[field] = state[field].data();
    return result;
  }  // END FUNCTION: GRTestBlock::state_pointers
  std::vector<const $SCALAR*> const_state_pointers() {
    std::vector<const $SCALAR*> result(NUM_EVOL_GFS);
    for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
      result[field] = state[field].data();
    return result;
  }  // END FUNCTION: GRTestBlock::const_state_pointers
  std::vector<$SCALAR*> rhs_pointers() {
    std::vector<$SCALAR*> result(NUM_EVOL_GFS);
    for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
      result[field] = rhs[field].data();
    return result;
  }  // END FUNCTION: GRTestBlock::rhs_pointers
  unsigned extent;
  std::size_t vol;
  std::vector<std::vector<$SCALAR>> state;
  std::vector<std::vector<$SCALAR>> rhs;
  BlockGeometry geometry{};
};  // END STRUCT: GRTestBlock

double gr_interior_rhs_max(const GRTestBlock& block) {
  double worst = 0.0;
  const unsigned padding = block.geometry.padding;
  for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
    for (unsigned k = padding; k < block.extent - padding; ++k)
      for (unsigned j = padding; j < block.extent - padding; ++j)
        for (unsigned i = padding; i < block.extent - padding; ++i)
          worst = std::max(worst, std::fabs(static_cast<double>(
              block.rhs[field][block.index(i, j, k)])));
  return worst;
}  // END FUNCTION: gr_interior_rhs_max

int test_padding() {
  const unsigned centred_radius = $NAMESPACE::generated::FD_ORDER / 2;
  return $NAMESPACE::generated::REQUIRED_PADDING >= centred_radius ? 0 : 1;
}  // END FUNCTION: test_padding

int test_offsets() {
  GRTestBlock block(2);
  block.geometry.component_offset = block.vol;
  const $SCALAR sentinel = static_cast<$SCALAR>(-7.5);
  for (auto& component : block.state)
    std::fill(component.begin(), component.end(), sentinel);
  std::vector<$SCALAR*> pointers = block.state_pointers();
  $MINKOWSKI_INITIAL_DATA_BLOCK(block.geometry, pointers.data());
  for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
    for (std::size_t cell = 0; cell < block.vol; ++cell) {
      if (block.state[field][cell] != sentinel) return 1;
      if (block.state[field][block.vol + cell] !=
          $NAMESPACE::generated::EVOL_GF_F_INFINITY[field]) return 2;
    }  // END LOOP: check component offsets
  return 0;
}  // END FUNCTION: test_offsets

int test_upwind() {
  if ($NAMESPACE::generated::NUM_UPWIND_CONTROL_GFS < 3) return 1;
  $NAMESPACE::generated::params_struct params;
  $PARAMS_STRUCT_SET_TO_DEFAULT(params);
  for (unsigned axis = 0; axis < 3; ++axis) {
    double sensitivity[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
    for (int sign_index = 0; sign_index < 2; ++sign_index)
      for (int side = 0; side < 2; ++side) {
        GRTestBlock block;
        std::vector<$SCALAR*> state = block.state_pointers();
        $MINKOWSKI_INITIAL_DATA_BLOCK(block.geometry, state.data());
        for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
          for (unsigned k = 0; k < block.extent; ++k)
            for (unsigned j = 0; j < block.extent; ++j)
              for (unsigned i = 0; i < block.extent; ++i)
                block.state[field][block.index(i, j, k)] +=
                    static_cast<$SCALAR>(1e-3 * (i + j + k));
        for (unsigned control = 0;
             control < $NAMESPACE::generated::NUM_UPWIND_CONTROL_GFS;
             ++control) {
          const unsigned field =
              $NAMESPACE::generated::EVOL_UPWIND_CONTROL_INDICES[control];
          const double value = control == axis
                                   ? (sign_index == 0 ? 0.5 : -0.5)
                                   : 0.5;
          std::fill(block.state[field].begin(), block.state[field].end(),
                    static_cast<$SCALAR>(value));
        }  // END LOOP: set upwind controls
        const unsigned p = block.geometry.padding;
        const unsigned r = $NAMESPACE::generated::REQUIRED_PADDING;
        const unsigned moved = side == 0 ? p + r : p - r;
        std::size_t moved_cell = block.index(moved, p, p);
        if (axis == 1) moved_cell = block.index(p, moved, p);
        if (axis == 2) moved_cell = block.index(p, p, moved);
        const std::size_t probe = block.index(p, p, p);
        std::vector<const $SCALAR*> input = block.const_state_pointers();
        std::vector<$SCALAR*> output = block.rhs_pointers();
        $RHS_EVAL_BLOCK(block.geometry, input.data(), output.data()$RHS_EVAL_BLOCK_TAIL);
        std::vector<double> before(NUM_EVOL_GFS);
        for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
          before[field] = block.rhs[field][probe];
        for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
          block.state[field][moved_cell] += static_cast<$SCALAR>(0.25);
        $RHS_EVAL_BLOCK(block.geometry, input.data(), output.data()$RHS_EVAL_BLOCK_TAIL);
        for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
          sensitivity[sign_index][side] = std::max(
              sensitivity[sign_index][side],
              std::fabs(static_cast<double>(block.rhs[field][probe]) -
                        before[field]));
      }  // END LOOP: test shifted stencil side
    const double scale = std::max(
        std::max(sensitivity[0][0], sensitivity[0][1]),
        std::max(sensitivity[1][0], sensitivity[1][1]));
    const double margin = 1e-6 * scale;
    if (!(scale > 0.0) ||
        !(sensitivity[0][0] - sensitivity[0][1] > margin) ||
        !(sensitivity[1][1] - sensitivity[1][0] > margin)) return 2;
  }  // END LOOP: test upwind axes
  return 0;
}  // END FUNCTION: test_upwind

int test_rhs() {
  $NAMESPACE::generated::params_struct params;
  $PARAMS_STRUCT_SET_TO_DEFAULT(params);
  GRTestBlock block;
  std::vector<$SCALAR*> state = block.state_pointers();
  $MINKOWSKI_INITIAL_DATA_BLOCK(block.geometry, state.data());
  std::vector<const $SCALAR*> input = block.const_state_pointers();
  std::vector<$SCALAR*> output = block.rhs_pointers();
  $RHS_EVAL_BLOCK(block.geometry, input.data(), output.data()$RHS_EVAL_BLOCK_TAIL);
  return gr_interior_rhs_max(block) <= 1e-13 ? 0 : 1;
}  // END FUNCTION: test_rhs

int test_init() {
  GRTestBlock block;
  std::vector<$SCALAR*> state = block.state_pointers();
  $MINKOWSKI_INITIAL_DATA_BLOCK(block.geometry, state.data());
  double norm = 0.0;
  for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
    for (std::size_t cell = 0; cell < block.vol; ++cell) {
      const double value = block.state[field][cell];
      norm += std::fabs(value);
      if (std::fabs(value - $NAMESPACE::generated::EVOL_GF_F_INFINITY[field]) >
          1e-15) return 1;
    }  // END LOOP: check initialized state
  return norm > 0.0 ? 0 : 2;
}  // END FUNCTION: test_init

int test_detgtrazero() {
  GRTestBlock block;
  std::vector<$SCALAR*> state = block.state_pointers();
  $MINKOWSKI_INITIAL_DATA_BLOCK(block.geometry, state.data());
  const auto flat = block.state;
  $NAMESPACE::generated::detgtrazero_status_struct flat_status{};
  $ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_BLOCK(block.geometry, state.data(), &flat_status);
  if (flat_status.failed_points || block.state != flat ||
      flat_status.max_abs_det_minus_one > 5e-13 ||
      flat_status.max_abs_trace_residual > 5e-13) return 1;
  $SMOOTH_PERTURBATION_BLOCK(block.geometry, state.data(),
                             static_cast<$SCALAR>(1e-2),
                             static_cast<$SCALAR>(1.0));
  $NAMESPACE::generated::detgtrazero_status_struct first_status{};
  $ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_BLOCK(block.geometry, state.data(), &first_status);
  if (first_status.failed_points || first_status.projected_points == 0) return 2;
  double state_scale = 1.0;
  for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
    for (std::size_t cell = 0; cell < block.vol; ++cell)
      state_scale = std::max(state_scale, std::fabs(static_cast<double>(
          block.state[field][cell])));
  const auto projected = block.state;
  $NAMESPACE::generated::detgtrazero_status_struct second_status{};
  $ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_BLOCK(block.geometry, state.data(), &second_status);
  double movement = 0.0;
  for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
    for (std::size_t cell = 0; cell < block.vol; ++cell)
      movement = std::max(movement, std::fabs(static_cast<double>(
          block.state[field][cell] - projected[field][cell])));
  return !second_status.failed_points && movement <= 5e-13 * state_scale &&
                 second_status.max_abs_det_minus_one <= 5e-13 &&
                 second_status.max_abs_trace_residual <= 5e-13 * state_scale
             ? 0
             : 3;
}  // END FUNCTION: test_detgtrazero

int test_constraints() {
  if ($NAMESPACE::generated::NUM_DIAG_GFS == 0) return 1;
  GRTestBlock block;
  std::vector<$SCALAR*> state = block.state_pointers();
  $MINKOWSKI_INITIAL_DATA_BLOCK(block.geometry, state.data());
  std::vector<std::vector<$SCALAR>> diagnostics(
      $NAMESPACE::generated::NUM_DIAG_GFS,
      std::vector<$SCALAR>(block.vol, $SCALAR{0}));
  std::vector<$SCALAR*> output(diagnostics.size());
  for (unsigned field = 0; field < output.size(); ++field)
    output[field] = diagnostics[field].data();
  std::vector<const $SCALAR*> input = block.const_state_pointers();
  $CONSTRAINTS_EVAL_BLOCK(block.geometry, input.data(), output.data());
  const unsigned padding = block.geometry.padding;
  double worst = 0.0;
  for (unsigned field = 0; field < diagnostics.size(); ++field)
    for (unsigned k = padding; k < block.extent - padding; ++k)
      for (unsigned j = padding; j < block.extent - padding; ++j)
        for (unsigned i = padding; i < block.extent - padding; ++i)
          worst = std::max(worst, std::fabs(static_cast<double>(
              diagnostics[field][block.index(i, j, k)])));
  return worst <= 1e-13 ? 0 : 2;
}  // END FUNCTION: test_constraints
"""

_GR_DISPATCH = r"""  if (std::strcmp(section, "padding") == 0) return test_padding();
  if (std::strcmp(section, "offsets") == 0) return test_offsets();
  if (std::strcmp(section, "upwind") == 0) return test_upwind();
  if (std::strcmp(section, "rhs") == 0) return test_rhs();
  if (std::strcmp(section, "init") == 0) return test_init();
  if (std::strcmp(section, "detgtrazero") == 0) return test_detgtrazero();
  if (std::strcmp(section, "constraints") == 0) return test_constraints();
"""

_GR_ALL = (
    "+ test_padding() + test_offsets() + test_upwind() + test_rhs() "
    "+ test_init() + test_detgtrazero() + test_constraints()"
)


def _coordinates(
    point: Tuple[int, int, int], spacings: Tuple[float, float, float]
) -> Tuple[float, float, float]:
    """
    Map an integer fixture point to physical coordinates.

    :param point: Integer grid indices.
    :param spacings: Grid spacing in each coordinate direction.
    :return: Physical coordinates about the fixed fixture centre.
    """
    centres = (6, 7, 8)
    return (
        (point[0] - centres[0]) * spacings[0],
        (point[1] - centres[1]) * spacings[1],
        (point[2] - centres[2]) * spacings[2],
    )


def _field_value(
    name: str,
    field_index: int,
    x: float,
    y: float,
    z: float,
    upwind_controls: Tuple[str, ...],
) -> float:
    """
    Return deterministic, algebraically admissible GR block data.

    :param name: Scientific evolved-field name.
    :param field_index: Canonical evolved component index.
    :param x: Physical x coordinate.
    :param y: Physical y coordinate.
    :param z: Physical z coordinate.
    :param upwind_controls: Scientific names controlling upwind selection.
    :return: Once-rounded field value at the requested coordinate.
    """
    if name in upwind_controls:
        value = (field_index + 1) * (x + x**6)
        return float(value)
    if name == "nrpy_reference_quadratic":
        return float(x * x + x * y)
    if name == "nrpy_reference_degree6":
        return float(x**6)
    q = 1.0e-3 * (x + 0.5 * y * y - 0.25 * z + 0.1 * x * y)
    amplitude = 2.0e-4 * (1.0 + x * x + y + 0.2 * z * z)
    if name == "hDD00":
        value = math.exp(2.0 * q) - 1.0
    elif name in ("hDD11", "hDD22"):
        value = math.exp(-q) - 1.0
    elif name in ("hDD01", "hDD02", "hDD12"):
        value = 0.0
    elif name == "aDD00":
        value = 2.0 * amplitude * math.exp(2.0 * q)
    elif name in ("aDD11", "aDD22"):
        value = -amplitude * math.exp(-q)
    elif name in ("aDD01", "aDD02", "aDD12"):
        value = 0.1 * amplitude * (x + y - z)
    else:
        polynomial = (
            0.1
            + x
            + 0.3 * y * y
            - 0.2 * z
            + 0.05 * x * y
            + 0.02 * x**6
            + 0.015 * y**6
            + 0.01 * z**6
        )
        asymptotic = float(gri.glb_gridfcs_dict[name].f_infinity)
        value = asymptotic + 1.0e-4 * (field_index + 1) * polynomial
    return float(value)


def _sample_value(
    samples: SampleMap,
    name: str,
    field_index: int,
    point: Tuple[int, int, int],
    spacings: Tuple[float, float, float],
    upwind_controls: Tuple[str, ...],
) -> float:
    """
    Return and retain the one binary64 sample shared by both oracle paths.

    :param samples: Mutable exact-sample cache.
    :param name: Scientific evolved-field name.
    :param field_index: Canonical evolved component index.
    :param point: Integer grid indices of the sample.
    :param spacings: Grid spacing in each coordinate direction.
    :param upwind_controls: Scientific names controlling upwind selection.
    :return: The cached binary64 value.
    """
    key = (name, point)
    if key not in samples:
        samples[key] = _field_value(
            name,
            field_index,
            *_coordinates(point, spacings),
            upwind_controls,
        )
    return samples[key]


def _exact_mpf(value: Any) -> Any:
    """
    Convert a binary64, integer, or exact SymPy rational to ``mpf`` exactly.

    :param value: Numeric value to convert.
    :return: Exact multiprecision representation at the active precision.
    """
    if isinstance(value, float):
        numerator, denominator = value.as_integer_ratio()
        return mp.mpf(numerator) / mp.mpf(denominator)
    if isinstance(value, sp.Rational):
        return mp.mpf(int(value.p)) / mp.mpf(int(value.q))
    return mp.mpf(value)


def _evaluate_reference(
    expression: sp.Expr,
    point: Tuple[int, int, int],
    spacings: Tuple[float, float, float],
    evol_order: Tuple[str, ...],
    upwind_controls: Tuple[str, ...],
    fd_order: int,
    samples: SampleMap,
) -> ReferenceValue:
    """
    Evaluate an RHS with shared CSE and explicit block substitutions.

    :param expression: Canonical scientific RHS expression.
    :param point: Integer grid indices at which to evaluate.
    :param spacings: Grid spacing in each coordinate direction.
    :param evol_order: Canonical evolved-field order.
    :param upwind_controls: Scientific names controlling upwind selection.
    :param fd_order: Finite-difference order.
    :param samples: Binary64 samples shared with the generated executable.
    :return: Stable reference and a scale-derived binary64 error allowance.
    :raises ValueError: If evaluation is complex or precision is unstable.
    """
    replaced, reduced = sp.cse(expression, order="none")
    operation_count = max(
        1,
        sum(int(sp.count_ops(rhs)) for _lhs, rhs in replaced)
        + sum(int(sp.count_ops(rhs)) for rhs in reduced),
    )
    values: List[Any] = []
    scales: List[Any] = []
    previous_dps = mp.dps
    try:
        for precision in (80, 100):
            mp.dps = precision
            environment: Dict[sp.Symbol, object] = {}
            stencil_scale = mp.mpf(0)
            field_indices = {name: index for index, name in enumerate(evol_order)}
            coordinates = _coordinates(point, spacings)
            for symbol in expression.free_symbols:
                name = str(symbol)
                if name in field_indices:
                    value = _exact_mpf(
                        _sample_value(
                            samples,
                            name,
                            field_indices[name],
                            point,
                            spacings,
                            upwind_controls,
                        )
                    )
                elif name in par.glb_code_params_dict:
                    value = _exact_mpf(par.glb_code_params_dict[name].defaultvalue)
                elif name in ("xx0", "xx1", "xx2"):
                    value = _exact_mpf(coordinates[int(name[-1])])
                elif any(
                    family in name for family in ("_dD", "_dupD", "_ddnD", "_dKOD")
                ):
                    bases, operators = (
                        extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars(
                            [symbol]
                        )
                    )
                    base_name, operator = bases[0], operators[0]
                    selected_operator = operator
                    if operator.startswith("dupD"):
                        axis = int(operator[-1])
                        control_name = upwind_controls[axis]
                        control_index = evol_order.index(control_name)
                        control = _sample_value(
                            samples,
                            control_name,
                            control_index,
                            point,
                            spacings,
                            upwind_controls,
                        )
                        selected_operator = (
                            operator
                            if control > 0.0
                            else operator.replace("dupD", "ddnD")
                        )
                    coefficients, offsets = compute_fdcoeffs_fdstencl(
                        selected_operator, fd_order
                    )
                    value = mp.mpf(0)
                    absolute_sum = mp.mpf(0)
                    for coefficient, offset in zip(coefficients, offsets):
                        shifted = (
                            point[0] + offset[0],
                            point[1] + offset[1],
                            point[2] + offset[2],
                        )
                        term = _exact_mpf(coefficient) * _exact_mpf(
                            _sample_value(
                                samples,
                                base_name,
                                field_indices[base_name],
                                shifted,
                                spacings,
                                upwind_controls,
                            )
                        )
                        value += term
                        absolute_sum += abs(term)
                    axes = [int(digit) for digit in operator if digit.isdigit()]
                    scale = mp.mpf(1)
                    for axis in axes:
                        scale /= _exact_mpf(spacings[axis])
                    value *= scale
                    absolute_sum *= abs(scale)
                    stencil_scale = max(stencil_scale, absolute_sum)
                else:
                    raise ValueError(f"Unclassified GR reference symbol {name!r}.")
                environment[symbol] = value
            result = inject_mpfs_into_cse_expression(environment, replaced, reduced)
            if hasattr(result, "imag") and result.imag != 0:
                raise ValueError(
                    "GR fixed-block reference unexpectedly became complex."
                )
            real_result = result.real if hasattr(result, "real") else result
            values.append(mp.mpf(real_result))
            environment_scale = max(
                (abs(mp.mpf(value)) for value in environment.values()),
                default=mp.mpf(0),
            )
            scales.append(
                max(mp.mpf(1), stencil_scale, environment_scale, abs(values[-1]))
            )
    finally:
        mp.dps = previous_dps
    precision_delta = abs(values[1] - values[0])
    # The lower-precision pass has 80 decimal digits.  Five guard digits and
    # the actual CSE operation count bound accumulated reference arithmetic;
    # this is independent of generated binary64 results.
    stability_bound = 32 * operation_count * mp.power(10, -75) * scales[1]
    if precision_delta > stability_bound:
        raise ValueError(
            "GR fixed-block reference did not stabilize between 80 and 100 digits."
        )
    # Standard gamma_n forward-error scaling for the actual per-component CSE
    # DAG.  The scale includes every substituted value, high-precision CSE
    # temporary, and the absolute stencil sums after inverse-spacing factors.
    # A fixed factor of 32, selected before kernel execution, covers different
    # valid CSE association and coefficient rounding in emitted binary64 C++.
    double_epsilon = mp.power(2, -52)
    gamma_n = operation_count * double_epsilon / (1 - operation_count * double_epsilon)
    roundoff_bound = (
        32 * gamma_n * scales[1] + 8 * precision_delta + double_epsilon * abs(values[1])
    )
    return ReferenceValue(
        float(values[1]),
        float(roundoff_bound),
        float(precision_delta),
        operation_count,
        float(scales[1]),
    )


def output_self_test_artifacts(
    solver_stem: str,
    solver_namespace: str,
    rhs_build: RHSBuild,
    enable_ko: bool,
) -> Dict[str, str]:
    """
    Return the GR test source and fixture companion headers.

    :param solver_stem: Lowercase formulation stem used in emitted paths.
    :param solver_namespace: Namespace containing the production solver.
    :param rhs_build: Canonical scientific RHS expressions and field order.
    :param enable_ko: Whether this generation profile includes KO dissipation.
    :return: Solver-root-relative paths mapped to complete artifact text.
    :raises ValueError: If configuration or reference validation is invalid.
    """
    fd_order = int(par.parval_from_str("fd_order"))
    if fd_order != 4:
        artifacts = generic_tests.output_self_test_artifacts(
            solver_stem,
            solver_namespace,
            application_test_functions=_GR_SCIENTIFIC_TESTS,
            application_dispatch=_GR_DISPATCH,
            application_all=_GR_ALL,
        )
        source_path = f"tests/{solver_stem}_self_tests.cpp"
        artifacts[source_path] = substitute_application_identifiers(
            artifacts[source_path]
        )
        return artifacts
    evol_order = tuple(rhs_build.evol_order)
    upwind_controls = tuple(rhs_build.upwind_control_fields)
    if len(upwind_controls) < 3:
        raise ValueError("GR nonflat reference requires three upwind controls.")
    spacings = (0.125, 0.25, 0.5)
    analytic_evol_order = (
        "nrpy_reference_quadratic",
        "nrpy_reference_degree6",
        "nrpy_reference_control0",
        "nrpy_reference_control1",
        "nrpy_reference_control2",
    )
    analytic_controls = analytic_evol_order[2:]
    analytic_expression = (
        sp.Symbol("nrpy_reference_quadratic_dD0")
        + sp.Symbol("nrpy_reference_quadratic_dDD01")
        + sp.Symbol("nrpy_reference_degree6_dupD0")
        + sp.Symbol("nrpy_reference_degree6_dKOD0")
    )
    analytic_samples: SampleMap = {}
    for analytic_point in ((3, 3, 3), (6, 7, 8), (9, 11, 13)):
        x, y, _z = _coordinates(analytic_point, spacings)
        control = _sample_value(
            analytic_samples,
            analytic_controls[0],
            2,
            analytic_point,
            spacings,
            analytic_controls,
        )
        upwind_operator = "dupD0" if control > 0.0 else "ddnD0"
        upwind_coefficients, upwind_offsets = compute_fdcoeffs_fdstencl(
            upwind_operator, 4
        )
        ko_coefficients, ko_offsets = compute_fdcoeffs_fdstencl("dKOD0", 4)
        upwind = mp.mpf(0)
        ko = mp.mpf(0)
        for coefficient, offset in zip(upwind_coefficients, upwind_offsets):
            shifted = (
                analytic_point[0] + offset[0],
                analytic_point[1],
                analytic_point[2],
            )
            sample = _sample_value(
                analytic_samples,
                analytic_evol_order[1],
                1,
                shifted,
                spacings,
                analytic_controls,
            )
            upwind += _exact_mpf(coefficient) * _exact_mpf(sample)
        for coefficient, offset in zip(ko_coefficients, ko_offsets):
            shifted = (
                analytic_point[0] + offset[0],
                analytic_point[1],
                analytic_point[2],
            )
            sample = _sample_value(
                analytic_samples,
                analytic_evol_order[1],
                1,
                shifted,
                spacings,
                analytic_controls,
            )
            ko += _exact_mpf(coefficient) * _exact_mpf(sample)
        upwind /= _exact_mpf(spacings[0])
        ko /= _exact_mpf(spacings[0])
        analytic = _exact_mpf(2.0 * x + y) + 1 + upwind + ko
        reference = _evaluate_reference(
            analytic_expression,
            analytic_point,
            spacings,
            analytic_evol_order,
            analytic_controls,
            4,
            analytic_samples,
        )
        if abs(reference.value - float(analytic)) > reference.bound:
            raise ValueError(
                "Complete FD4 reference pipeline failed an analytic identity."
            )
        if abs(ko) <= _exact_mpf(1.0e-12):
            raise ValueError("Independent FD4 reference has a zero KO discriminator.")
    points = ((3, 3, 3), (6, 7, 8), (9, 11, 13))
    first_control = upwind_controls[0]
    first_control_index = evol_order.index(first_control)
    signs = []
    for point in points:
        coordinates = _coordinates(point, spacings)
        control = _field_value(
            first_control, first_control_index, *coordinates, upwind_controls
        )
        signs.append(1 if control > 0.0 else -1 if control < 0.0 else 0)
    if tuple(signs) != (-1, 0, 1):
        raise ValueError(
            "GR reference points do not cover negative, zero, positive upwind controls."
        )
    expression_by_field = {
        gf_names.rhs_symbol_to_gridfunction_name(name): expression
        for name, expression in rhs_build.rhs_by_symbol_name.items()
    }
    samples: SampleMap = {}
    references: List[ReferenceValue] = []
    reference_by_field_point: Dict[Tuple[str, Tuple[int, int, int]], ReferenceValue] = (
        {}
    )
    for point in points:
        for name in evol_order:
            reference = _evaluate_reference(
                expression_by_field[name],
                point,
                spacings,
                evol_order,
                upwind_controls,
                fd_order,
                samples,
            )
            references.append(reference)
            reference_by_field_point[(name, point)] = reference
    ko_discriminators = 0
    minimum_ko_ratio = math.inf
    if enable_ko:
        family_discriminated: Dict[str, bool] = {}
        family_best_ratio: Dict[str, float] = {}
        for name, expression in expression_by_field.items():
            ko_symbols = {
                symbol for symbol in expression.free_symbols if "_dKOD" in str(symbol)
            }
            if not ko_symbols:
                continue
            tensor_family = gf_names.tensor_family_of(name)
            family_name = tensor_family[0] if tensor_family is not None else name
            family_discriminated.setdefault(family_name, False)
            field_discriminated = False
            without_ko = expression.xreplace(
                {symbol: sp.Integer(0) for symbol in ko_symbols}
            )
            for point in points:
                off_reference = _evaluate_reference(
                    without_ko,
                    point,
                    spacings,
                    evol_order,
                    upwind_controls,
                    fd_order,
                    samples,
                )
                on_reference = reference_by_field_point[(name, point)]
                difference = abs(on_reference.value - off_reference.value)
                combined_bound = on_reference.bound + off_reference.bound
                family_best_ratio[family_name] = max(
                    family_best_ratio.get(family_name, 0.0),
                    difference / combined_bound,
                )
                if difference > 8.0 * combined_bound:
                    field_discriminated = True
                    family_discriminated[family_name] = True
                    minimum_ko_ratio = min(
                        minimum_ko_ratio, difference / combined_bound
                    )
            ko_discriminators += int(field_discriminated)
        resolved_families = sorted(
            family
            for family, discriminated in family_discriminated.items()
            if discriminated
        )
        # The fixture deliberately designates several robust families instead
        # of claiming sensitivity where a tiny tensor KO term lies below the
        # independently derived full-RHS roundoff bound.
        if ko_discriminators < 3 or len(resolved_families) < 3:
            raise ValueError(
                "KO-on profile needs at least three component and family "
                "discriminators; "
                f"resolved={resolved_families}, ratios={family_best_ratio}."
            )
    expected_values = ",\n      ".join(
        f"{reference.value:.17g}" for reference in references
    )
    expected_bounds = ",\n      ".join(
        f"{reference.bound:.17g}" for reference in references
    )
    maximum_operation_count = max(reference.operation_count for reference in references)
    maximum_evaluation_scale = max(
        reference.evaluation_scale for reference in references
    )
    maximum_precision_delta = max(reference.precision_delta for reference in references)
    block_name = roles.CFunction_name_for_role("rhs_eval_block")
    flat_name = roles.CFunction_name_for_role("rhs_eval_flat_block")
    block_tail = generic_context._codeparameter_tail(block_name, "params")
    value_lines = [
        "double gr_reference_value(unsigned f, double x, double y, double z) {",
        "  const double q = 1.0e-3*(x + 0.5*y*y - 0.25*z + 0.1*x*y);",
        "  const double amplitude = 2.0e-4*(1.0 + x*x + y + 0.2*z*z);",
    ]
    for index, name in enumerate(evol_order):
        if name in upwind_controls:
            value = f"{index + 1}.0*(x + std::pow(x, 6))"
        elif name == "hDD00":
            value = "std::exp(2.0*q) - 1.0"
        elif name in ("hDD11", "hDD22"):
            value = "std::exp(-q) - 1.0"
        elif name in ("hDD01", "hDD02", "hDD12"):
            value = "0.0"
        elif name == "aDD00":
            value = "2.0*amplitude*std::exp(2.0*q)"
        elif name in ("aDD11", "aDD22"):
            value = "-amplitude*std::exp(-q)"
        elif name in ("aDD01", "aDD02", "aDD12"):
            value = "0.1*amplitude*(x+y-z)"
        else:
            asymptotic = float(gri.glb_gridfcs_dict[name].f_infinity)
            value = (
                f"{asymptotic:.17g} + 1.0e-4*{index + 1}.0*"
                "(0.1+x+0.3*y*y-0.2*z+0.05*x*y+0.02*std::pow(x,6)"
                "+0.015*std::pow(y,6)+0.01*std::pow(z,6))"
            )
        value_lines.append(f"  if (f == {index}u) return {value};  // {name}")
    value_lines.extend(("  return 0.0;", "}  // END FUNCTION: gr_reference_value"))
    value_function = "\n".join(value_lines)
    field_indices = {name: index for index, name in enumerate(evol_order)}
    exact_sample_entries = []
    for (name, point), sample_value in sorted(
        samples.items(), key=lambda item: (field_indices[item[0][0]], item[0][1])
    ):
        exact_sample_entries.append(
            "    {"
            f"{field_indices[name]}u,{point[0]}u,{point[1]}u,{point[2]}u,"
            f"{sample_value.hex()}"
            f"}},  // {name}"
        )
    exact_samples = "\n".join(exact_sample_entries)
    formulation = "fCCZ4/chi" if "Theta_fCCZ4" in evol_order else "BSSN/W"
    reference_cpp = f"""
{value_function}

int test_gr_nonflat_reference() {{
  constexpr unsigned nx=13, ny=15, nz=17, pad=3;
  constexpr std::size_t offset=11;
  constexpr std::size_t vol=static_cast<std::size_t>(nx)*ny*nz;
  constexpr double sentinel=-54321.25;
  const double dx[3]={{0.125,0.25,0.5}};
  BlockGeometry geom{{}}; geom.nx=nx; geom.ny=ny; geom.nz=nz;
  geom.padding=pad; geom.component_offset=offset;
  geom.dx[0]=dx[0]; geom.dx[1]=dx[1]; geom.dx[2]=dx[2];
  std::vector<std::vector<DendroScalar>> input(NUM_EVOL_GFS,
      std::vector<DendroScalar>(offset+vol+17,sentinel));
  std::vector<std::vector<DendroScalar>> output(NUM_EVOL_GFS,
      std::vector<DendroScalar>(offset+vol+17,sentinel));
  auto index=[](unsigned a,unsigned b,unsigned c){{
    return static_cast<std::size_t>(a)+nx*(static_cast<std::size_t>(b)+ny*c);
  }};
  for (unsigned c=0;c<nz;++c) for (unsigned b=0;b<ny;++b)
    for (unsigned a=0;a<nx;++a) {{
      const double x=(static_cast<int>(a)-6)*dx[0];
      const double y=(static_cast<int>(b)-7)*dx[1];
      const double z=(static_cast<int>(c)-8)*dx[2];
      for (unsigned f=0;f<NUM_EVOL_GFS;++f)
        input[f][offset+index(a,b,c)]=gr_reference_value(f,x,y,z);
    }}  // END LOOP: populate reference input
  struct ExactSample {{unsigned f,a,b,c; double value;}};  // END STRUCT: ExactSample
  const ExactSample exact_samples[]={{
{exact_samples}
  }};
  // Every input used by a checked stencil is one hexadecimal binary64 value
  // emitted from the same cache used by the multiprecision reference.
  for(const ExactSample& sample:exact_samples)
    input[sample.f][offset+index(sample.a,sample.b,sample.c)]=sample.value;
  for(const ExactSample& sample:exact_samples)
    if(input[sample.f][offset+index(sample.a,sample.b,sample.c)]!=sample.value) {{
      std::fprintf(stderr,"FAIL: exact input sample f=%u cell=(%u,%u,%u)\\n",
          sample.f,sample.a,sample.b,sample.c); return 1;
    }}  // END IF: exact input mismatch
  const auto before=input;
  std::vector<const DendroScalar*> in(NUM_EVOL_GFS);
  std::vector<DendroScalar*> out(NUM_EVOL_GFS);
  for (unsigned f=0;f<NUM_EVOL_GFS;++f) {{in[f]=input[f].data();out[f]=output[f].data();}}  // END LOOP: bind block pointers
  $NAMESPACE::generated::params_struct params;
  $PARAMS_STRUCT_SET_TO_DEFAULT(params);
  {block_name}(geom,in.data(),out.data(){block_tail});
  const double expected[]={{
      {expected_values}
  }};
  const double expected_bounds[]={{
      {expected_bounds}
  }};
  const unsigned points[][3]={{{{3,3,3}},{{6,7,8}},{{9,11,13}}}};
  // Per-component bounds were fixed before execution from each actual CSE DAG,
  // its exact scaled stencil sums, all inputs and CSE temporaries. Across this
  // profile: max operations={maximum_operation_count}, max scale={maximum_evaluation_scale:.17g},
  // max 80/100-digit delta={maximum_precision_delta:.17g}; KO discriminators={ko_discriminators},
  // families={','.join(resolved_families) if enable_ko else 'none'},
  // minimum KO difference/bound ratio={minimum_ko_ratio if enable_ko else 0.0:.17g}.
  double actual_norm=0.0;
  for (unsigned p=0;p<3;++p) for (unsigned f=0;f<NUM_EVOL_GFS;++f) {{
    const double actual=output[f][offset+index(points[p][0],points[p][1],points[p][2])];
    const double reference=expected[p*NUM_EVOL_GFS+f];
    const double bound=expected_bounds[p*NUM_EVOL_GFS+f];
    actual_norm=std::max(actual_norm,std::fabs(actual));
    if (!std::isfinite(actual) || std::fabs(actual-reference)>bound) {{
      std::fprintf(stderr,"FAIL: {formulation} FD4 KO={'on' if enable_ko else 'off'} "
          "component %u point %u actual %.17g reference %.17g bound %.17g\\n",
          f,p,actual,reference,bound); return 1;
    }}  // END IF: block reference mismatch
  }}  // END LOOP: check block references
  if (!(actual_norm>1.0e-8) || input!=before) {{
    std::fprintf(stderr,"FAIL: zero output norm or modified {formulation} input\\n");
    return 2;
  }}  // END IF: invalid block output
  for(unsigned f=0;f<NUM_EVOL_GFS;++f) for(std::size_t cell=0;
      cell<output[f].size();++cell) {{
    const bool interior=cell>=offset && cell<offset+vol &&
        ((cell-offset)%nx)>=pad && ((cell-offset)%nx)<nx-pad &&
        (((cell-offset)/nx)%ny)>=pad && (((cell-offset)/nx)%ny)<ny-pad &&
        ((cell-offset)/(nx*ny))>=pad && ((cell-offset)/(nx*ny))<nz-pad;
    if(!interior && output[f][cell]!=sentinel) {{
      std::fprintf(stderr,"FAIL: block sentinel f=%u cell=%zu value=%.17g\\n",
          f,cell,static_cast<double>(output[f][cell])); return 3;
    }}  // END IF: block sentinel changed
  }}  // END LOOP: check block sentinels
  std::vector<DendroScalar> flat_in(offset+NUM_EVOL_GFS*vol+17,sentinel);
  std::vector<DendroScalar> flat_out(offset+NUM_EVOL_GFS*vol+17,sentinel);
  for(unsigned f=0;f<NUM_EVOL_GFS;++f) for(std::size_t cell=0;cell<vol;++cell)
    flat_in[offset+static_cast<std::size_t>(f)*vol+cell]=input[f][offset+cell];
  {flat_name}(geom,flat_in.data(),flat_out.data(){block_tail});
  for(unsigned p=0;p<3;++p) for(unsigned f=0;f<NUM_EVOL_GFS;++f) {{
    const std::size_t cell=index(points[p][0],points[p][1],points[p][2]);
    const double actual=flat_out[offset+static_cast<std::size_t>(f)*vol+cell];
    const double reference=expected[p*NUM_EVOL_GFS+f];
    const double bound=expected_bounds[p*NUM_EVOL_GFS+f];
    if (!std::isfinite(actual) || std::fabs(actual-reference)>bound) {{
      std::fprintf(stderr,"FAIL: {formulation} flat FD4 KO={'on' if enable_ko else 'off'} "
          "component %u point %u actual %.17g reference %.17g bound %.17g\\n",
          f,p,actual,reference,bound); return 4;
    }}  // END IF: flat reference mismatch
  }}  // END LOOP: check flat references
  for(std::size_t cell=0;cell<flat_out.size();++cell) {{
    bool writable=false;
    if(cell>=offset && cell<offset+NUM_EVOL_GFS*vol) {{
      const std::size_t local=(cell-offset)%vol;
      writable=(local%nx)>=pad && (local%nx)<nx-pad &&
          ((local/nx)%ny)>=pad && ((local/nx)%ny)<ny-pad &&
          (local/(nx*ny))>=pad && (local/(nx*ny))<nz-pad;
    }}  // END IF: identify writable flat cell
    if(!writable && flat_out[cell]!=sentinel) {{
      std::fprintf(stderr,"FAIL: flat sentinel cell=%zu value=%.17g\\n",
          cell,static_cast<double>(flat_out[cell])); return 5;
    }}  // END IF: flat sentinel changed
  }}  // END LOOP: check flat sentinels
  return 0;
}}  // END FUNCTION: test_gr_nonflat_reference
"""
    dispatch = (
        '  if (std::strcmp(section, "gr_nonflat_reference") == 0) '
        "return test_gr_nonflat_reference();\n"
    )
    artifacts = generic_tests.output_self_test_artifacts(
        solver_stem,
        solver_namespace,
        application_test_functions=_GR_SCIENTIFIC_TESTS + reference_cpp,
        application_dispatch=_GR_DISPATCH + dispatch,
        application_all=_GR_ALL + " + test_gr_nonflat_reference()",
    )
    source_path = f"tests/{solver_stem}_self_tests.cpp"
    artifacts[source_path] = substitute_application_identifiers(artifacts[source_path])
    return artifacts
