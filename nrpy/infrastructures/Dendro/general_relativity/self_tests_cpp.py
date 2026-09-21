# nrpy/infrastructures/Dendro/general_relativity/self_tests_cpp.py
"""
Assemble GR scientific sections into the generated Dendro self tests.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import math
from typing import Dict, List, Mapping, NamedTuple, Tuple, Union, cast

import sympy as sp
from mpmath import mp, mpf  # type: ignore[import-untyped]

import nrpy.grid as gri
import nrpy.params as par
from nrpy.finite_difference import (
    extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars,
)
from nrpy.infrastructures.Dendro import CFunction_roles as roles
from nrpy.infrastructures.Dendro import gridfunction_name_decorations as gf_names
from nrpy.infrastructures.Dendro import self_tests_cpp as generic_tests
from nrpy.infrastructures.Dendro import solver_context as generic_context
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


# Scientific names are explicit application input to CMake.  The two test
# sections exercise formulation-neutral lowering through the same executable.
def test_sections() -> Tuple[str, ...]:
    """
    Return meaningful test sections for the configured FD order.

    :return: Ordered CTest section names for the current generation profile.
    """
    sections: Tuple[str, ...] = (
        "state",
        "params",
        "offsets",
        "stencil_reach",
        "rhs",
        "init",
        "names",
        "detgtrazero",
    )
    return sections + ("address_values", "parameter_forwarding", "gr_nonflat_reference")


_GR_SCIENTIFIC_TESTS = r"""
struct GRTestBlock {
  /**
   * Allocate padded arrays for generated state.
   *
   * @param blocks Number of component-offset-sized blocks to allocate.
   * @param padding Padding width on each block face.
   */
  explicit GRTestBlock(
      unsigned blocks = 1,
      unsigned padding = $NAMESPACE::generated::REQUIRED_PADDING)
      : extent(2 * padding + 3),
        vol(static_cast<std::size_t>(extent) * extent * extent),
        state(NUM_EVOL_GFS, std::vector<$SCALAR>(vol * blocks, $SCALAR{0})),
        rhs(NUM_EVOL_GFS, std::vector<$SCALAR>(vol * blocks, $SCALAR{0})) {
    geometry.nx = geometry.ny = geometry.nz = extent;
    geometry.padding = padding;
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
  block_geometry_struct geometry{};
};  // END STRUCT: GRTestBlock

/**
 * Measure the maximum generated RHS magnitude over a test-grid interior.
 *
 * @param[in] block Test state whose RHS arrays are inspected.
 * @return Maximum absolute interior RHS value.
 */
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

/**
 * Check that generated initial data honors component offsets.
 *
 * @return 0 on success; 1 for a modified prefix or 2 for bad initial data.
 */
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

/**
 * Check centered sensitivity and declared stencil reach on every axis.
 *
 * @return 0 on success; 1 for a sensitivity or reach failure.
 */
int test_stencil_reach() {
  $NAMESPACE::generated::params_struct params;
  $PARAMS_STRUCT_SET_TO_DEFAULT(params);
  const unsigned reach = $NAMESPACE::generated::REQUIRED_PADDING;
  const unsigned safe_padding = reach + 1;
  for (unsigned axis = 0; axis < 3; ++axis) {
    double sensitivity[2] = {0.0, 0.0};
    double outside_sensitivity = 0.0;
    for (int side = 0; side < 2; ++side) {
        GRTestBlock block(1, safe_padding);
        std::vector<$SCALAR*> state = block.state_pointers();
        $MINKOWSKI_INITIAL_DATA_BLOCK(block.geometry, state.data());
        for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
          for (unsigned k = 0; k < block.extent; ++k)
            for (unsigned j = 0; j < block.extent; ++j)
              for (unsigned i = 0; i < block.extent; ++i)
                block.state[field][block.index(i, j, k)] +=
                    static_cast<$SCALAR>(1e-3 * (i + j + k));
        const unsigned p = safe_padding + 1;
        const unsigned moved = side == 0 ? p + reach : p - reach;
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
        std::vector<$SCALAR> moved_values(NUM_EVOL_GFS);
        for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
          moved_values[field] = block.state[field][moved_cell];
        for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
          block.state[field][moved_cell] =
              moved_values[field] + static_cast<$SCALAR>(0.25);
        $RHS_EVAL_BLOCK(block.geometry, input.data(), output.data()$RHS_EVAL_BLOCK_TAIL);
        for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
          sensitivity[side] = std::max(
              sensitivity[side],
              std::fabs(static_cast<double>(block.rhs[field][probe]) -
                        before[field]));
        for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
          block.state[field][moved_cell] = moved_values[field];
        const unsigned outside =
            side == 0 ? p + reach + 1 : p - reach - 1;
        std::size_t outside_cell = block.index(outside, p, p);
        if (axis == 1) outside_cell = block.index(p, outside, p);
        if (axis == 2) outside_cell = block.index(p, p, outside);
        for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
          block.state[field][outside_cell] += static_cast<$SCALAR>(0.25);
        $RHS_EVAL_BLOCK(block.geometry, input.data(), output.data()$RHS_EVAL_BLOCK_TAIL);
        for (unsigned field = 0; field < NUM_EVOL_GFS; ++field)
          outside_sensitivity = std::max(
              outside_sensitivity,
              std::fabs(static_cast<double>(block.rhs[field][probe]) -
                        before[field]));
    }  // END LOOP: test both centered stencil sides
    if (!(sensitivity[0] > 0.0) || !(sensitivity[1] > 0.0) ||
        outside_sensitivity != 0.0) return 1;
  }  // END LOOP: test centered axes
  return 0;
}  // END FUNCTION: test_stencil_reach

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

/**
 * Check generated Minkowski initial data against asymptotic values.
 *
 * @return 0 on success; 1 for a value mismatch or 2 for an empty state norm.
 */
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

/**
 * Check algebraic-constraint projection, movement, and idempotence.
 *
 * @return 0 on success; 1 for changed flat data, 2 for failed first projection,
 * or 3 for a non-idempotent or inaccurate second projection.
 */
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

"""

_GR_DISPATCH = r"""  if (std::strcmp(section, "offsets") == 0) return test_offsets();
  if (std::strcmp(section, "stencil_reach") == 0) return test_stencil_reach();
  if (std::strcmp(section, "rhs") == 0) return test_rhs();
  if (std::strcmp(section, "init") == 0) return test_init();
  if (std::strcmp(section, "detgtrazero") == 0) return test_detgtrazero();
"""

_GR_ALL = (
    "+ test_offsets() + test_stencil_reach() + test_rhs() "
    "+ test_init() + test_detgtrazero()"
)


def _coordinates(
    point: Tuple[int, int, int], spacings: Tuple[float, float, float]
) -> Tuple[float, float, float]:
    """
    Map an integer test point to physical coordinates.

    :param point: Integer grid indices.
    :param spacings: Grid spacing in each coordinate direction.
    :return: Physical coordinates about the fixed test centre.
    """
    centres = (8, 9, 10)
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
) -> float:
    """
    Return deterministic, algebraically admissible GR block data.

    :param name: Scientific evolved-field name.
    :param field_index: Canonical evolved component index.
    :param x: Physical x coordinate.
    :param y: Physical y coordinate.
    :param z: Physical z coordinate.
    :return: Once-rounded field value at the requested coordinate.
    """
    if name == "nrpy_reference_quadratic":
        return float(x * x + x * y)
    if name == "nrpy_reference_degree8":
        return float(x**8)
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
            + 0.02 * x**8
            + 0.015 * y**8
            + 0.01 * z**8
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
) -> float:
    """
    Return and retain the one binary64 sample shared by both oracle paths.

    :param samples: Mutable exact-sample cache.
    :param name: Scientific evolved-field name.
    :param field_index: Canonical evolved component index.
    :param point: Integer grid indices of the sample.
    :param spacings: Grid spacing in each coordinate direction.
    :return: The cached binary64 value.
    """
    key = (name, point)
    if key not in samples:
        samples[key] = _field_value(
            name,
            field_index,
            *_coordinates(point, spacings),
        )
    return samples[key]


def _exact_mpf(value: Union[float, int, str, sp.Rational, mpf]) -> mpf:
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
    fd_order: int,
    ko_fd_order: int,
    samples: SampleMap,
) -> ReferenceValue:
    """
    Evaluate an RHS with shared CSE and explicit block substitutions.

    :param expression: Canonical scientific RHS expression.
    :param point: Integer grid indices at which to evaluate.
    :param spacings: Grid spacing in each coordinate direction.
    :param evol_order: Canonical evolved-field order.
    :param fd_order: Finite-difference order.
    :param ko_fd_order: Base order of the KO difference.
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
    values: List[mpf] = []
    scales: List[mpf] = []
    previous_dps = mp.dps
    try:
        for precision in (80, 100):
            mp.dps = precision
            environment: Dict[sp.Symbol, object] = {}
            stencil_scale = mp.mpf(0)
            field_indices = {name: index for index, name in enumerate(evol_order)}
            coordinates = _coordinates(point, spacings)
            for free_symbol in expression.free_symbols:
                symbol = cast(sp.Symbol, free_symbol)
                name = str(symbol)
                if name in field_indices:
                    value = _exact_mpf(
                        _sample_value(
                            samples,
                            name,
                            field_indices[name],
                            point,
                            spacings,
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
                    if operator.startswith(("dupD", "ddnD")):
                        raise ValueError(
                            f"Directional operator {operator!r} remains in centered GR reference."
                        )
                    coefficients, offsets = generic_tests._independent_stencil(
                        operator, fd_order, ko_fd_order
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
    rhs_by_symbol_name: Mapping[str, sp.Expr],
    diagnostics_by_name: Mapping[str, sp.Expr],
    *,
    fd_order: int,
    ko_fd_order: int,
    enable_ko: bool,
) -> Dict[str, str]:
    """
    Return the GR test source and companion headers.

    :param solver_stem: Lowercase formulation stem used in emitted paths.
    :param solver_namespace: Namespace containing the production solver.
    :param rhs_by_symbol_name: Scientific RHS expressions keyed by output name.
    :param diagnostics_by_name: Diagnostic expressions keyed by gridfunction.
    :param fd_order: Centered finite-difference order.
    :param ko_fd_order: Base order supplied to the KO difference operator.
    :param enable_ko: Whether the RHS contains Kreiss-Oliger dissipation.
    :return: Solver-root-relative paths mapped to complete file contents.
    :raises ValueError: If configuration or reference validation is invalid.
    """
    if fd_order not in (4, 6, 8) or ko_fd_order != fd_order - 2:
        raise ValueError("GR reference requires a Dendro FD4/6/8 profile.")
    evol_order = roles.registered_evol_order()
    spacings = (0.125, 0.25, 0.5)
    analytic_evol_order = (
        "nrpy_reference_quadratic",
        "nrpy_reference_degree8",
    )
    analytic_expression = (
        sp.Symbol("nrpy_reference_quadratic_dD0")
        + sp.Symbol("nrpy_reference_quadratic_dDD01")
        + sp.Symbol("nrpy_reference_degree8_dD0")
        + sp.Symbol("nrpy_reference_degree8_dKOD0")
    )
    analytic_samples: SampleMap = {}
    for analytic_point in ((4, 4, 4), (8, 9, 10), (12, 14, 16)):
        x, y, _z = _coordinates(analytic_point, spacings)
        centered_coefficients, centered_offsets = generic_tests._independent_stencil(
            "dD0", fd_order, ko_fd_order
        )
        ko_coefficients, ko_offsets = generic_tests._independent_stencil(
            "dKOD0", fd_order, ko_fd_order
        )
        centered = mp.mpf(0)
        ko = mp.mpf(0)
        for coefficient, offset in zip(centered_coefficients, centered_offsets):
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
            )
            centered += _exact_mpf(coefficient) * _exact_mpf(sample)
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
            )
            ko += _exact_mpf(coefficient) * _exact_mpf(sample)
        centered /= _exact_mpf(spacings[0])
        ko /= _exact_mpf(spacings[0])
        analytic = _exact_mpf(2.0 * x + y) + 1 + centered + ko
        reference = _evaluate_reference(
            analytic_expression,
            analytic_point,
            spacings,
            analytic_evol_order,
            fd_order,
            ko_fd_order,
            analytic_samples,
        )
        if abs(reference.value - float(analytic)) > reference.bound:
            raise ValueError(
                f"Complete FD{fd_order} reference calculation failed an analytic identity."
            )
        if abs(ko) <= _exact_mpf(1.0e-12):
            raise ValueError(
                f"Independent FD{fd_order} reference has a zero KO discriminator."
            )
    points = ((4, 4, 4), (8, 9, 10), (12, 14, 16))
    expression_by_field = {
        gf_names.rhs_symbol_to_gridfunction_name(name): expression
        for name, expression in rhs_by_symbol_name.items()
    }
    diag_order = tuple(gri.GridFunction.gridfunction_lists()[2])
    diagnostic_expressions = dict(diagnostics_by_name)
    if set(diag_order) != set(diagnostic_expressions):
        raise ValueError(
            "GR nonflat reference diagnostic order does not match the "
            f"registered expressions: order={diag_order}, "
            f"expressions={tuple(diagnostic_expressions)}."
        )
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
                fd_order,
                ko_fd_order,
                samples,
            )
            references.append(reference)
            reference_by_field_point[(name, point)] = reference
    diagnostic_references: List[ReferenceValue] = []
    for point in points:
        for name in diag_order:
            diagnostic_references.append(
                _evaluate_reference(
                    diagnostic_expressions[name],
                    point,
                    spacings,
                    evol_order,
                    fd_order,
                    ko_fd_order,
                    samples,
                )
            )
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
                    fd_order,
                    ko_fd_order,
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
        # The test deliberately designates several robust families instead
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
    expected_diagnostic_values = ",\n      ".join(
        f"{reference.value:.17g}" for reference in diagnostic_references
    )
    expected_diagnostic_bounds = ",\n      ".join(
        f"{reference.bound:.17g}" for reference in diagnostic_references
    )
    all_references = references + diagnostic_references
    maximum_operation_count = max(
        reference.operation_count for reference in all_references
    )
    maximum_evaluation_scale = max(
        reference.evaluation_scale for reference in all_references
    )
    maximum_precision_delta = max(
        reference.precision_delta for reference in all_references
    )
    block_name = roles.CFunction_name_for_role("rhs_eval_block")
    flat_name = roles.CFunction_name_for_role("rhs_eval_flat_block")
    constraints_name = roles.CFunction_name_for_role("constraints_eval_block")
    block_tail = generic_context.codeparameter_tail(block_name, "params")
    value_lines = [
        "/**",
        " * Evaluate one analytic nonflat reference-state component.",
        " *",
        " * @param f Generated evolved-component index.",
        " * @param x Physical x coordinate.",
        " * @param y Physical y coordinate.",
        " * @param z Physical z coordinate.",
        " * @return Analytic component value, or zero for an unknown index.",
        " */",
        "double gr_reference_value(unsigned f, double x, double y, double z) {",
        "  const double q = 1.0e-3*(x + 0.5*y*y - 0.25*z + 0.1*x*y);",
        "  const double amplitude = 2.0e-4*(1.0 + x*x + y + 0.2*z*z);",
    ]
    for index, name in enumerate(evol_order):
        if name == "hDD00":
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
                "(0.1+x+0.3*y*y-0.2*z+0.05*x*y+0.02*std::pow(x,8)"
                "+0.015*std::pow(y,8)+0.01*std::pow(z,8))"
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

/**
 * Check nonflat block, diagnostic, and flat adapters against frozen oracles.
 *
 * @return 0 on success; 1 for reference/input encoding failure, 2 for mutated
 * input, 3 or 7 for sentinel corruption, 4 for flat-reference mismatch, 5 for
 * flat sentinel corruption, or 6 for diagnostic-reference mismatch.
 */
int test_gr_nonflat_reference() {{
  constexpr unsigned nx=17, ny=19, nz=21;
  constexpr unsigned pad=$NAMESPACE::generated::REQUIRED_PADDING;
  constexpr std::size_t offset=11;
  constexpr std::size_t vol=static_cast<std::size_t>(nx)*ny*nz;
  constexpr double sentinel=-54321.25;
  const double dx[3]={{0.125,0.25,0.5}};
  block_geometry_struct geom{{}}; geom.nx=nx; geom.ny=ny; geom.nz=nz;
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
      const double x=(static_cast<int>(a)-8)*dx[0];
      const double y=(static_cast<int>(b)-9)*dx[1];
      const double z=(static_cast<int>(c)-10)*dx[2];
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
  const double expected_diagnostics[]={{
      {expected_diagnostic_values}
  }};
  const double expected_diagnostic_bounds[]={{
      {expected_diagnostic_bounds}
  }};
  const unsigned points[][3]={{{{4,4,4}},{{8,9,10}},{{12,14,16}}}};
  // Per-component bounds were fixed before execution from each actual CSE DAG,
  // its exact scaled stencil sums, all inputs and CSE temporaries. Across this
  // profile: max operations={maximum_operation_count}, max scale={maximum_evaluation_scale:.17g},
  // max 80/100-digit delta={maximum_precision_delta:.17g}; KO discriminators={ko_discriminators},
  // families={','.join(resolved_families) if enable_ko else 'none'},
  // minimum KO difference/bound ratio={minimum_ko_ratio if enable_ko else 0.0:.17g}.
  for (unsigned p=0;p<3;++p) for (unsigned f=0;f<NUM_EVOL_GFS;++f) {{
    const double actual=output[f][offset+index(points[p][0],points[p][1],points[p][2])];
    const double reference=expected[p*NUM_EVOL_GFS+f];
    const double bound=expected_bounds[p*NUM_EVOL_GFS+f];
    if (!std::isfinite(actual) || std::fabs(actual-reference)>bound) {{
      std::fprintf(stderr,"FAIL: {formulation} FD{fd_order} KO={'on' if enable_ko else 'off'} "
          "component %u point %u actual %.17g reference %.17g bound %.17g\\n",
          f,p,actual,reference,bound); return 1;
    }}  // END IF: block reference mismatch
  }}  // END LOOP: check block references
  if (input!=before) {{
    std::fprintf(stderr,"FAIL: modified {formulation} input\\n");
    return 2;
  }}  // END IF: invalid block output
  std::vector<std::vector<DendroScalar>> diagnostics(
      $NAMESPACE::generated::NUM_DIAG_GFS,
      std::vector<DendroScalar>(offset+vol+17,sentinel));
  std::vector<DendroScalar*> diagnostic_out(
      $NAMESPACE::generated::NUM_DIAG_GFS);
  for(unsigned f=0;f<$NAMESPACE::generated::NUM_DIAG_GFS;++f)
    diagnostic_out[f]=diagnostics[f].data();
  {constraints_name}(geom,in.data(),diagnostic_out.data());
  for(unsigned p=0;p<3;++p)
    for(unsigned f=0;f<$NAMESPACE::generated::NUM_DIAG_GFS;++f) {{
    const double actual=diagnostics[f][offset+index(
        points[p][0],points[p][1],points[p][2])];
    const double reference=expected_diagnostics[
        p*$NAMESPACE::generated::NUM_DIAG_GFS+f];
    const double bound=expected_diagnostic_bounds[
        p*$NAMESPACE::generated::NUM_DIAG_GFS+f];
    if(!std::isfinite(actual) || std::fabs(actual-reference)>bound) {{
      std::fprintf(stderr,"FAIL: {formulation} FD{fd_order} diagnostic component %u "
          "point %u actual %.17g reference %.17g bound %.17g\\n",
          f,p,actual,reference,bound); return 6;
    }}  // END IF: diagnostic reference mismatch
  }}  // END LOOP: check diagnostic references
  for(unsigned f=0;f<$NAMESPACE::generated::NUM_DIAG_GFS;++f)
    for(std::size_t cell=0;
      cell<diagnostics[f].size();++cell) {{
    const bool interior=cell>=offset && cell<offset+vol &&
        ((cell-offset)%nx)>=pad && ((cell-offset)%nx)<nx-pad &&
        (((cell-offset)/nx)%ny)>=pad && (((cell-offset)/nx)%ny)<ny-pad &&
        ((cell-offset)/(nx*ny))>=pad && ((cell-offset)/(nx*ny))<nz-pad;
    if(!interior && diagnostics[f][cell]!=sentinel) {{
      std::fprintf(stderr,"FAIL: diagnostic sentinel f=%u cell=%zu value=%.17g\\n",
          f,cell,static_cast<double>(diagnostics[f][cell])); return 7;
    }}  // END IF: diagnostic sentinel changed
  }}  // END LOOP: check diagnostic sentinels
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
      std::fprintf(stderr,"FAIL: {formulation} flat FD{fd_order} KO={'on' if enable_ko else 'off'} "
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
