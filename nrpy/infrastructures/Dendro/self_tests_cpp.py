# nrpy/infrastructures/Dendro/self_tests_cpp.py
"""
Emit the generated project's self-test executable.

Every test is driven from generated headers and registered CFunctions; the
standalone host types stand in for the pinned Dendrolib API.  Production field
names and asymptotic values remain application-owned.  The local two-field
fixture is temporary, formulation-neutral, and removed from the registries after
its generated artifacts are assembled.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Dict, List, Tuple

import sympy as sp

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.finite_difference import compute_fdcoeffs_fdstencl
from nrpy.infrastructures.Dendro import CFunction_roles as roles
from nrpy.infrastructures.Dendro import CodeParameters, block_kernel_helpers, types_h
from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.solver_context import (
    _codeparameter_tail,
    substitute_solver_identifiers,
)

BANNER = generated_file_banner()

# The CTest case names, in the order the runner dispatches them.  The tests
# CMake registers one case per entry, so the two cannot drift.
SECTIONS: Tuple[str, ...] = (
    "state",
    "params",
    "names",
)

_TESTS = """// Usage: $STEM_self_tests {state|params|names|application-section|all}
// Exits 0 on success.  Each section runs only what it names, so a failure
// localises to one gate.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "$STEM_defines.h"
#include "$STEM_fixture_parameters.h"

$FIXTURE_FUNCTIONS

namespace {

using $NAMESPACE::generated::NUM_EVOL_GFS;

int test_state() {
  // The generated state registry must be internally complete: one name per
  // enum member, all names distinct (a duplicate would alias two evolved
  // fields onto one identity), and one metadata entry per field.
  if (NUM_EVOL_GFS == 0) return 1;
  if ($NAMESPACE::generated::to_index($NAMESPACE::generated::EvolVar::END) !=
      NUM_EVOL_GFS) {
    return 1;
  }
  if ($NAMESPACE::generated::EVOL_GF_NAMES.size() != NUM_EVOL_GFS) return 1;
  if ($NAMESPACE::generated::EVOL_GF_RANK.size() != NUM_EVOL_GFS) return 1;
  if ($NAMESPACE::generated::EVOL_GF_F_INFINITY.size() != NUM_EVOL_GFS)
    return 1;
  if ($NAMESPACE::generated::EVOL_GF_WAVESPEED.size() != NUM_EVOL_GFS)
    return 1;
  for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
    if ($NAMESPACE::generated::EVOL_GF_NAMES[f].empty()) return 1;
    for (unsigned g = f + 1; g < NUM_EVOL_GFS; ++g) {
      if ($NAMESPACE::generated::EVOL_GF_NAMES[f] ==
          $NAMESPACE::generated::EVOL_GF_NAMES[g]) {
        return 1;
      }
    }  // END LOOP: for g over later components
    if (!($NAMESPACE::generated::EVOL_GF_WAVESPEED[f] > 0.0)) return 1;
  }  // END LOOP: for f over evolved components
  return 0;
}  // END FUNCTION: test_state

int test_params() {
  // The registered parameter CFunctions set, validate, and print the
  // generated defaults.
  $NAMESPACE::generated::params_struct params;
  $PARAMS_STRUCT_SET_TO_DEFAULT(params);
  if (!$VALIDATE(params)) return 1;
  $PRINT_EFFECTIVE(params);
  return 0;
}  // END FUNCTION: test_params

int test_names() {
  // Every generated name resolves to its own group and index, matching is
  // case-sensitive, and an unknown name resolves to nothing.
  for (const $NAMESPACE::generated::VariableRef::Group group :
       $NAMESPACE::generated::VARIABLE_GROUPS) {
    const unsigned count = $NAMESPACE::generated::variable_count(group);
    for (unsigned index = 0; index < count; ++index) {
      const std::string_view name =
          $NAMESPACE::generated::variable_name(group, index);
      if (name.empty()) return 1;
      const std::optional<$NAMESPACE::generated::VariableRef> ref =
          $NAMESPACE::generated::find_variable(name);
      if (!ref.has_value()) return 2;
      if (ref->group != group || ref->index != index) return 3;
      // A case-flipped spelling must not resolve unless it is itself a
      // generated name: NRPy variance suffixes are case-sensitive.
      std::string flipped(name);
      bool changed = false;
      for (char& c : flipped) {
        if (c >= 'a' && c <= 'z') {
          c = static_cast<char>(c - 'a' + 'A');
          changed = true;
        } else if (c >= 'A' && c <= 'Z') {  // END IF: lowercase letter flipped
          c = static_cast<char>(c - 'A' + 'a');
          changed = true;
        }  // END ELSE IF: uppercase letter flipped
      }  // END LOOP: for c over name characters
      if (changed && flipped != std::string(name)) {
        bool flipped_is_generated = false;
        for (const $NAMESPACE::generated::VariableRef::Group other :
             $NAMESPACE::generated::VARIABLE_GROUPS) {
          const unsigned n = $NAMESPACE::generated::variable_count(other);
          for (unsigned k = 0; k < n; ++k) {
            if ($NAMESPACE::generated::variable_name(other, k) == flipped) {
              flipped_is_generated = true;
            }  // END IF: flipped spelling is itself generated
          }  // END LOOP: for k over other names
        }  // END LOOP: for other over generated groups
        const bool resolved =
            $NAMESPACE::generated::find_variable(flipped).has_value();
        if (resolved != flipped_is_generated) return 4;
      }  // END IF: case flip changed spelling
    }  // END LOOP: for index over group names
  }  // END LOOP: for group over generated groups
  if ($NAMESPACE::generated::find_variable("no_such_generated_variable")
          .has_value()) {
    return 5;
  }
  if ($NAMESPACE::generated::find_variable("").has_value()) return 6;
  return 0;
}  // END FUNCTION: test_names

$ADDITIONAL_TEST_FUNCTIONS

int run_section(const char* section) {
  if (std::strcmp(section, "state") == 0) return test_state();
  if (std::strcmp(section, "params") == 0) return test_params();
  if (std::strcmp(section, "names") == 0) return test_names();
$ADDITIONAL_DISPATCH
  if (std::strcmp(section, "all") == 0) {
    const int rc = test_state() + test_params() + test_names()
                   $ADDITIONAL_ALL;
    return rc == 0 ? 0 : 1;
  }  // END IF: section is all
  std::fprintf(stderr, "unknown section %s\\n", section);
  return 2;
}  // END FUNCTION: run_section

// clang-format off
}  // END NAMESPACE: internal linkage
// clang-format on

int main(int argc, char* argv[]) {
  const char* section = (argc > 1) ? argv[1] : "all";
  const int rc = run_section(section);
  if (rc != 0) std::fprintf(stderr, "FAIL: section %s rc=%d\\n", section, rc);
  return rc;
}  // END FUNCTION: main
"""


def output_self_tests_cpp(
    solver_stem: str,
    solver_namespace: str,
    fixture_functions: str,
    additional_test_functions: str,
    additional_dispatch: str,
    additional_all: str,
) -> str:
    r"""
    Emit the generated project's self-test source.

    :param solver_stem: Lowercase formulation stem for emitted file names.
    :param solver_namespace: Solver namespace.
    :param fixture_functions: Complete generated fixture functions.
    :param additional_test_functions: Application-owned test functions.
    :param additional_dispatch: Application section dispatch statements.
    :param additional_all: Application terms in the aggregate runner.
    :return: The complete C++ source text.

    Doctests:
    >>> _TESTS.count("// clang-format off") == _TESTS.count("}  // END NAMESPACE:")
    True
    >>> _TESTS.count("}  // END NAMESPACE:")
    1
    """
    text = _TESTS
    replacements = (
        ("$FIXTURE_FUNCTIONS", fixture_functions),
        ("$ADDITIONAL_TEST_FUNCTIONS", additional_test_functions),
        ("$ADDITIONAL_DISPATCH", additional_dispatch),
        ("$ADDITIONAL_ALL", additional_all),
    )
    for token, value in replacements:
        text = text.replace(token, value)
    return BANNER + substitute_solver_identifiers(text, solver_stem, solver_namespace)


def _stencil_terms(name: str, operator: str) -> str:
    """
    Emit one coefficient/offset table from the shared FD owner.

    :param name: C++ table name.
    :param operator: Canonical finite-difference operator name.
    :return: Complete C++ table declaration.
    """
    coefficients, stencil = compute_fdcoeffs_fdstencl(operator, 4)
    entries = ",\n".join(
        "    {"
        + f"{float(coefficient):.17g}, {offset[0]}, {offset[1]}, {offset[2]}"
        + "}"
        for coefficient, offset in zip(coefficients, stencil)
    )
    return f"static constexpr StencilTerm {name}[] = {{\n{entries}\n}};\n"


_FIXTURE_CPP = """
struct StencilTerm {{ double coefficient; int di; int dj; int dk; }};
{tables}

double fixture_u(double x, double, double) {{ return x + x*x*x*x*x*x; }}
double fixture_v(double x, double y, double z) {{
  return 2.0 + 0.3*x*x + 0.1*y*y*y*y*y*y + 0.2*z + 0.05*x*y;
}}

template <std::size_t N, class Sample>
double fixture_stencil(const StencilTerm (&terms)[N], Sample sample,
                       unsigned a, unsigned b, unsigned c, double scale) {{
  double value = 0.0;
  for (const StencilTerm& term : terms) {{
    value += term.coefficient * sample(
        static_cast<int>(a) + term.di, static_cast<int>(b) + term.dj,
        static_cast<int>(c) + term.dk);
  }}  // END LOOP: stencil terms
  return value * scale;
}}  // END FUNCTION: fixture_stencil

int test_address_values(bool alternate_parameters = false) {{
  constexpr unsigned nx = 13, ny = 15, nz = 17, pad = 3;
  constexpr std::size_t offset = 11;
  constexpr std::size_t vol = static_cast<std::size_t>(nx) * ny * nz;
  constexpr double sentinel = -9876.25;
  const double dx[3] = {{0.125, 0.25, 0.5}};
  block_geometry_struct geom{{}};
  geom.nx = nx; geom.ny = ny; geom.nz = nz; geom.padding = pad;
  geom.component_offset = offset;
  geom.dx[0] = dx[0]; geom.dx[1] = dx[1]; geom.dx[2] = dx[2];
  std::vector<std::vector<DendroScalar>> input(
      2, std::vector<DendroScalar>(offset + vol + 13, sentinel));
  std::vector<std::vector<DendroScalar>> output(
      2, std::vector<DendroScalar>(offset + vol + 13, sentinel));
  auto index = [](int a, int b, int c) {{
    return static_cast<std::size_t>(a) + nx *
           (static_cast<std::size_t>(b) + ny * static_cast<std::size_t>(c));
  }};
  for (unsigned c = 0; c < nz; ++c) for (unsigned b = 0; b < ny; ++b)
    for (unsigned a = 0; a < nx; ++a) {{
      const double x = (static_cast<int>(a) - 6) * dx[0];
      const double y = (static_cast<int>(b) - 7) * dx[1];
      const double z = (static_cast<int>(c) - 8) * dx[2];
      input[0][offset + index(a,b,c)] = fixture_u(x,y,z);
      input[1][offset + index(a,b,c)] = fixture_v(x,y,z);
    }}  // END LOOP: input grid points
  const auto input_before = input;
  const DendroScalar* in_ptr[] = {{input[0].data(), input[1].data()}};
  DendroScalar* out_ptr[] = {{output[0].data(), output[1].data()}};
  $NAMESPACE_fixture::generated::params_struct fixture_params{{}};
  fixture_params.nrpy_fixture_common = alternate_parameters ? -0.25 : 0.375;
  fixture_params.nrpy_fixture_count = alternate_parameters ? 5 : 3;
  fixture_params.nrpy_fixture_real = alternate_parameters ? 0.75 : 1.25;
  fixture_params.nrpy_fixture_toggle = !alternate_parameters;
  {block_name}(geom, in_ptr, out_ptr{call_tail});

  auto sample_u = [&](int a, int b, int c) {{
    return static_cast<double>(input[0][offset + index(a,b,c)]);
  }};
  auto sample_v = [&](int a, int b, int c) {{
    return static_cast<double>(input[1][offset + index(a,b,c)]);
  }};
  const unsigned points[][3] = {{{{3,3,3}}, {{6,7,8}}, {{9,11,13}}}};
  for (const auto& point : points) {{
    const unsigned a = point[0], b = point[1], c = point[2];
    const double control = sample_u(a,b,c);
    const double u_d0 = fixture_stencil(D1, sample_u, a,b,c, 1.0/dx[0]);
    const double u_up = fixture_stencil(control > 0.0 ? UP : DOWN,
                                         sample_u, a,b,c, 1.0/dx[0]);
    const double u_ko = fixture_stencil(KO, sample_u, a,b,c, 1.0/dx[0]);
    const double v_mix = fixture_stencil(
        MIXED, sample_v, a,b,c, 1.0/(dx[0]*dx[1]));
    const double v_d2 = fixture_stencil(D1, [&](int da,int db,int dc) {{
      return sample_v(dc,db,da);
    }}, c,b,a, 1.0/dx[2]);
    const double v_up1 = fixture_stencil(control > 0.0 ? UP : DOWN,
        [&](int da,int db,int dc) {{ return sample_v(db,da,dc); }},
        b,a,c, 1.0/dx[1]);
    const double v_ko1 = fixture_stencil(KO,
        [&](int da,int db,int dc) {{ return sample_v(db,da,dc); }},
        b,a,c, 1.0/dx[1]);
    const double expected_u = fixture_params.nrpy_fixture_real*sample_u(a,b,c)
        + fixture_params.nrpy_fixture_count*u_d0
        + fixture_params.nrpy_fixture_common*v_mix
        + (fixture_params.nrpy_fixture_toggle ? sample_v(a,b,c)
                                              : sample_u(a,b,c))
        + u_up + u_ko;
    const double expected_v = fixture_params.nrpy_fixture_real*sample_v(a,b,c)
        - fixture_params.nrpy_fixture_count*v_d2
        + fixture_params.nrpy_fixture_common*v_mix
        + (fixture_params.nrpy_fixture_toggle ? sample_u(a,b,c)
                                              : sample_v(a,b,c))
        + v_up1 + v_ko1;
    const double actual_u = output[0][offset + index(a,b,c)];
    const double actual_v = output[1][offset + index(a,b,c)];
    const double bound_u = 2e-12 * std::max(1.0, std::fabs(expected_u));
    const double bound_v = 2e-12 * std::max(1.0, std::fabs(expected_v));
    if (std::fabs(actual_u-expected_u) > bound_u ||
        std::fabs(actual_v-expected_v) > bound_v) {{
      std::fprintf(stderr, "FAIL: address/value point (%u,%u,%u): "
          "u %.17g expected %.17g; v %.17g expected %.17g\\n",
          a,b,c,actual_u,expected_u,actual_v,expected_v);
      return 1;
    }}  // END IF: reference mismatch
  }}  // END LOOP: reference points
  if (input != input_before) return 2;
  for (unsigned f = 0; f < 2; ++f) {{
    for (std::size_t cell = 0; cell < output[f].size(); ++cell) {{
      const bool interior = cell >= offset && cell < offset + vol &&
          ((cell-offset) % nx) >= pad && ((cell-offset) % nx) < nx-pad &&
          (((cell-offset)/nx) % ny) >= pad &&
          (((cell-offset)/nx) % ny) < ny-pad &&
          ((cell-offset)/(nx*ny)) >= pad &&
          ((cell-offset)/(nx*ny)) < nz-pad;
      if (!interior && output[f][cell] != sentinel) return 3;
    }}  // END LOOP: output cells
  }}  // END LOOP: output fields

  std::vector<DendroScalar> flat_input(offset + 2*vol + 13, sentinel);
  std::vector<DendroScalar> flat_output(offset + 2*vol + 13, sentinel);
  for (unsigned f = 0; f < 2; ++f)
    for (std::size_t cell = 0; cell < vol; ++cell)
      flat_input[offset + static_cast<std::size_t>(f)*vol + cell] =
          input[f][offset + cell];
  {flat_name}(geom, flat_input.data(), flat_output.data(){call_tail});
  for (const auto& point : points) {{
    const std::size_t cell = index(point[0],point[1],point[2]);
    for (unsigned f = 0; f < 2; ++f) {{
      if (flat_output[offset + static_cast<std::size_t>(f)*vol + cell] !=
          output[f][offset + cell]) return 4;
    }}  // END LOOP: flat output fields
  }}  // END LOOP: flat reference points
  for (std::size_t cell = 0; cell < flat_output.size(); ++cell) {{
    bool writable = false;
    if (cell >= offset && cell < offset + 2*vol) {{
      const std::size_t local = (cell-offset) % vol;
      writable = (local % nx) >= pad && (local % nx) < nx-pad &&
          ((local/nx) % ny) >= pad && ((local/nx) % ny) < ny-pad &&
          (local/(nx*ny)) >= pad && (local/(nx*ny)) < nz-pad;
    }}  // END IF: candidate writable cell
    if (!writable && flat_output[cell] != sentinel) return 5;
  }}  // END LOOP: flat sentinel cells
  return 0;
}}  // END FUNCTION: test_address_values

int test_parameter_forwarding() {{
  using Signature = void (*)(const block_geometry_struct&, const DendroScalar* const*,
      DendroScalar* const*, const DendroScalar, const int,
      const DendroScalar, const bool);
  static_assert(std::is_same_v<decltype(&{block_name}), Signature>,
                "fixture parameter types/order changed");
  static_assert(std::is_same_v<decltype(
      $NAMESPACE_fixture::generated::params_struct::nrpy_fixture_count), int>);
  static_assert(std::is_same_v<decltype(
      $NAMESPACE_fixture::generated::params_struct::nrpy_fixture_toggle), bool>);
  $NAMESPACE_fixture::generated::params_struct params{{}};
  params.nrpy_fixture_unused = 19.0;
  if (params.nrpy_fixture_unused != 19.0) return 1;
  // A second complete execution changes every forwarded value and flips the
  // Boolean branch.  The same independent stencil oracle checks the result.
  return test_address_values(true);
}}  // END FUNCTION: test_parameter_forwarding
"""


def output_self_test_artifacts(
    solver_stem: str,
    solver_namespace: str,
    application_test_functions: str,
    application_dispatch: str,
    application_all: str,
) -> Dict[str, str]:
    """
    Emit the self-test source and isolated tiny-fixture headers.

    Temporary registrations are collision checked and removed individually;
    production registry objects are never cleared or replaced.

    :param solver_stem: Lowercase formulation stem used in emitted paths.
    :param solver_namespace: Namespace containing the production solver.
    :param application_test_functions: Application-owned C++ test functions.
    :param application_dispatch: Application section dispatch statements.
    :param application_all: Application terms in the aggregate runner.
    :return: Solver-root-relative paths mapped to complete artifact text.
    :raises ValueError: If fixture names collide or fixture contracts disagree.

    Success and a failure after both fixture functions are registered restore
    every pre-existing registry and sidecar object.

    >>> _saved_infrastructure = par.parval_from_str("Infrastructure")
    >>> _saved_parallelization = par.parval_from_str("parallelization")
    >>> _saved_fp_type = par.parval_from_str("fp_type")
    >>> _saved_fd_order = par.parval_from_str("fd_order")
    >>> _saved_fields = dict(gri.glb_gridfcs_dict)
    >>> _saved_parameters = dict(par.glb_code_params_dict)
    >>> _saved_functions = dict(cfc.CFunction_dict)
    >>> _saved_extras = dict(par.glb_extras_dict)
    >>> _saved_dendro_present = "Dendro" in par.glb_extras_dict
    >>> _saved_dendro = par.glb_extras_dict.get("Dendro")
    >>> _saved_codeparameters_present = _saved_dendro is not None and "CFunction_codeparameters" in _saved_dendro
    >>> _saved_codeparameters = None if _saved_dendro is None else _saved_dendro.get("CFunction_codeparameters")
    >>> _saved_roles_present = _saved_dendro is not None and "CFunction_roles" in _saved_dendro
    >>> _saved_roles = None if _saved_dendro is None else _saved_dendro.get("CFunction_roles")
    >>> try:
    ...     par.set_parval_from_str("Infrastructure", "Dendro")
    ...     par.set_parval_from_str("parallelization", "none")
    ...     par.set_parval_from_str("fp_type", "double")
    ...     _artifacts = output_self_test_artifacts("probe", "probe", "", "", "")
    ...     assert sorted(_artifacts) == ["tests/probe_fixture_parameters.h", "tests/probe_fixture_types.h", "tests/probe_self_tests.cpp"]
    ...     assert set(gri.glb_gridfcs_dict) == set(_saved_fields)
    ...     assert all(gri.glb_gridfcs_dict[name] is value for name, value in _saved_fields.items())
    ...     assert set(par.glb_code_params_dict) == set(_saved_parameters)
    ...     assert all(par.glb_code_params_dict[name] is value for name, value in _saved_parameters.items())
    ...     assert set(cfc.CFunction_dict) == set(_saved_functions)
    ...     assert all(cfc.CFunction_dict[name] is value for name, value in _saved_functions.items())
    ...     assert set(par.glb_extras_dict) == set(_saved_extras)
    ...     assert all(par.glb_extras_dict[name] is value for name, value in _saved_extras.items())
    ...     assert ("Dendro" in par.glb_extras_dict) == _saved_dendro_present
    ...     assert par.glb_extras_dict.get("Dendro") is _saved_dendro
    ...     assert (_saved_dendro is not None and "CFunction_codeparameters" in _saved_dendro) == _saved_codeparameters_present
    ...     assert (_saved_dendro is None or _saved_dendro.get("CFunction_codeparameters") is _saved_codeparameters)
    ...     assert (_saved_dendro is not None and "CFunction_roles" in _saved_dendro) == _saved_roles_present
    ...     assert (_saved_dendro is None or _saved_dendro.get("CFunction_roles") is _saved_roles)
    ...     assert par.parval_from_str("fd_order") == _saved_fd_order
    ...     _owner_globals = output_self_test_artifacts.__globals__
    ...     _original_output = _owner_globals["output_self_tests_cpp"]
    ...     def _fail_after_registration(*_args, **_kwargs):
    ...         raise RuntimeError("injected post-registration failure")
    ...     try:
    ...         _owner_globals["output_self_tests_cpp"] = _fail_after_registration
    ...         try:
    ...             output_self_test_artifacts("probe", "probe", "", "", "")
    ...         except RuntimeError as error:
    ...             assert str(error) == "injected post-registration failure"
    ...         else:
    ...             raise AssertionError("injected fixture failure was not observed")
    ...     finally:
    ...         _owner_globals["output_self_tests_cpp"] = _original_output
    ...     assert set(gri.glb_gridfcs_dict) == set(_saved_fields)
    ...     assert all(gri.glb_gridfcs_dict[name] is value for name, value in _saved_fields.items())
    ...     assert set(par.glb_code_params_dict) == set(_saved_parameters)
    ...     assert all(par.glb_code_params_dict[name] is value for name, value in _saved_parameters.items())
    ...     assert set(cfc.CFunction_dict) == set(_saved_functions)
    ...     assert all(cfc.CFunction_dict[name] is value for name, value in _saved_functions.items())
    ...     assert set(par.glb_extras_dict) == set(_saved_extras)
    ...     assert all(par.glb_extras_dict[name] is value for name, value in _saved_extras.items())
    ...     assert ("Dendro" in par.glb_extras_dict) == _saved_dendro_present
    ...     assert par.glb_extras_dict.get("Dendro") is _saved_dendro
    ...     assert (_saved_dendro is not None and "CFunction_codeparameters" in _saved_dendro) == _saved_codeparameters_present
    ...     assert (_saved_dendro is None or _saved_dendro.get("CFunction_codeparameters") is _saved_codeparameters)
    ...     assert (_saved_dendro is not None and "CFunction_roles" in _saved_dendro) == _saved_roles_present
    ...     assert (_saved_dendro is None or _saved_dendro.get("CFunction_roles") is _saved_roles)
    ...     assert par.parval_from_str("fd_order") == _saved_fd_order
    ... finally:
    ...     par.set_parval_from_str("Infrastructure", _saved_infrastructure)
    ...     par.set_parval_from_str("parallelization", _saved_parallelization)
    ...     par.set_parval_from_str("fp_type", _saved_fp_type)
    >>> par.parval_from_str("Infrastructure") == _saved_infrastructure
    True
    >>> par.parval_from_str("parallelization") == _saved_parallelization
    True
    >>> par.parval_from_str("fp_type") == _saved_fp_type
    True
    >>> par.parval_from_str("fd_order") == _saved_fd_order
    True
    """
    field_names = ("nrpy_fixture_u", "nrpy_fixture_v")
    parameter_specs = (
        ("REAL", "nrpy_fixture_common", 0.375, True),
        ("int", "nrpy_fixture_count", 3, False),
        ("REAL", "nrpy_fixture_real", 1.25, False),
        ("bool", "nrpy_fixture_toggle", True, False),
        ("REAL", "nrpy_fixture_unused", 9.5, False),
    )
    function_names = (
        f"{solver_stem}_fixture_rhs_block",
        f"{solver_stem}_fixture_rhs_flat_block",
    )
    collisions: List[str] = [
        name for name in field_names if name in gri.glb_gridfcs_dict
    ]
    collisions += [
        name
        for _kind, name, _default, _common in parameter_specs
        if name in par.glb_code_params_dict
    ]
    collisions += [name for name in function_names if name in cfc.CFunction_dict]
    if collisions:
        raise ValueError(
            f"Self-test fixture registration collision: {sorted(collisions)}"
        )
    added_fields: List[str] = []
    added_parameters: List[str] = []
    added_functions: List[str] = []
    saved_fd_order = par.parval_from_str("fd_order")
    dendro_extras_existed = "Dendro" in par.glb_extras_dict
    dendro_extras = par.glb_extras_dict.get("Dendro", {})
    codeparameter_sidecar_existed = "CFunction_codeparameters" in dendro_extras
    role_sidecar_existed = "CFunction_roles" in dendro_extras
    try:
        # The formulation-neutral fixture is deliberately FD4 in every
        # application build; production FD order is restored in ``finally``.
        par.set_parval_from_str("fd_order", 4)
        u, v = gri.register_gridfunctions(list(field_names), group="EVOL")
        added_fields.extend(field_names)
        parameter_symbols = {}
        for kind, name, default, commondata in parameter_specs:
            parameter_symbols[name] = par.register_CodeParameter(
                kind, __name__, name, default, commondata=commondata
            )
            added_parameters.append(name)
        u_dD = ixp.declarerank1("nrpy_fixture_u_dD")
        u_dupD = ixp.declarerank1("nrpy_fixture_u_dupD")
        u_dKOD = ixp.declarerank1("nrpy_fixture_u_dKOD")
        v_dD = ixp.declarerank1("nrpy_fixture_v_dD")
        v_dupD = ixp.declarerank1("nrpy_fixture_v_dupD")
        v_dKOD = ixp.declarerank1("nrpy_fixture_v_dKOD")
        v_dDD = ixp.declarerank2("nrpy_fixture_v_dDD")
        toggle = parameter_symbols["nrpy_fixture_toggle"]
        conditional_u = sp.Piecewise((v, sp.Eq(toggle, 1)), (u, True))
        conditional_v = sp.Piecewise((u, sp.Eq(toggle, 1)), (v, True))
        expressions = (
            parameter_symbols["nrpy_fixture_real"] * u
            + parameter_symbols["nrpy_fixture_count"] * u_dD[0]
            + parameter_symbols["nrpy_fixture_common"] * v_dDD[0][1]
            + conditional_u
            + u_dupD[0]
            + u_dKOD[0],
            parameter_symbols["nrpy_fixture_real"] * v
            - parameter_symbols["nrpy_fixture_count"] * v_dD[2]
            + parameter_symbols["nrpy_fixture_common"] * v_dDD[0][1]
            + conditional_v
            + v_dupD[1]
            + v_dKOD[1],
        )
        kernel = c_codegen(
            list(expressions),
            [f"rhs_{name}[pp]" for name in field_names],
            enable_fd_codegen=True,
            enable_fd_functions=False,
            enable_simd=False,
            fp_type=str(par.parval_from_str("fp_type")),
            fp_type_alias=gri.DENDRO_SCALAR_TYPE,
            mem_alloc_style="210",
            rational_const_alias="static const",
            verbose=False,
            upwind_control_vec=[u, u, u],
        )
        used = block_kernel_helpers.used_codeparameters(expressions)
        declarations = block_kernel_helpers.cparam_declarations(used)
        block_body = (
            block_kernel_helpers.block_pointer_bindings(
                field_names, gri.DENDRO_SCALAR_TYPE
            )
            + "\n"
            + block_kernel_helpers.point_loop(kernel, "3")
        )
        params = (
            f"const block_geometry_struct& geom, const {gri.DENDRO_SCALAR_TYPE}* const* in_gfs, "
            f"{gri.DENDRO_SCALAR_TYPE}* const* rhs_gfs, {declarations}"
        )
        cfc.register_CFunction(
            includes=[f"{solver_stem}_fixture_parameters.h", "block_geometry.h"],
            desc="Tiny executable Dendro address and parameter fixture.",
            cfunc_type="void",
            name=function_names[0],
            params=params,
            body=block_body,
        )
        added_functions.append(function_names[0])
        roles.set_CFunction_codeparameters(function_names[0], used)
        flat_body = block_kernel_helpers.flat_block_pointer_bindings(
            field_names, gri.DENDRO_SCALAR_TYPE
        ) + (
            f"\nconst {gri.DENDRO_SCALAR_TYPE}* const in_gfs_call[] = "
            "{in_nrpy_fixture_u, in_nrpy_fixture_v};\n"
            f"{gri.DENDRO_SCALAR_TYPE}* rhs_gfs_call[] = "
            "{rhs_nrpy_fixture_u, rhs_nrpy_fixture_v};\n"
            f"{function_names[0]}(geom, in_gfs_call, rhs_gfs_call, "
            + block_kernel_helpers.cparam_arguments(used)
            + ");"
        )
        flat_params = (
            f"const block_geometry_struct& geom, const {gri.DENDRO_SCALAR_TYPE}* in_gfs_flat, "
            f"{gri.DENDRO_SCALAR_TYPE}* rhs_gfs_flat, {declarations}"
        )
        cfc.register_CFunction(
            includes=[f"{solver_stem}_fixture_parameters.h", "block_geometry.h"],
            desc="Flat-layout adapter for the tiny executable fixture.",
            cfunc_type="void",
            name=function_names[1],
            params=flat_params,
            body=flat_body,
        )
        added_functions.append(function_names[1])
        roles.set_CFunction_codeparameters(function_names[1], used)
        function_source = "\n".join(
            cfc.CFunction_dict[name].full_function for name in function_names
        )
        block_tail = _codeparameter_tail(function_names[0], "fixture_params")
        flat_tail = _codeparameter_tail(function_names[1], "fixture_params")
        if block_tail != flat_tail:
            raise ValueError(
                "Fixture block and flat adapters must forward identical parameters."
            )
        tables = "".join(
            (
                _stencil_terms("D1", "dD0"),
                _stencil_terms("UP", "dupD0"),
                _stencil_terms("DOWN", "ddnD0"),
                _stencil_terms("KO", "dKOD0"),
                _stencil_terms("MIXED", "dDD01"),
            )
        )
        fixture_tests = _FIXTURE_CPP.format(
            tables=tables,
            block_name=function_names[0],
            call_tail=block_tail,
            flat_name=function_names[1],
        ).replace("$NAMESPACE", solver_namespace)
        source = output_self_tests_cpp(
            solver_stem,
            solver_namespace,
            function_source,
            fixture_tests + application_test_functions,
            '  if (std::strcmp(section, "address_values") == 0) return test_address_values();\n'
            '  if (std::strcmp(section, "parameter_forwarding") == 0) return test_parameter_forwarding();\n'
            + application_dispatch,
            "+ test_address_values() + test_parameter_forwarding() " + application_all,
        )
        fixture_stem = f"{solver_stem}_fixture"
        fixture_namespace = f"{solver_namespace}_fixture"
        return {
            f"tests/{solver_stem}_self_tests.cpp": source,
            f"tests/{solver_stem}_fixture_types.h": types_h.output_types_h(
                fixture_stem, fixture_namespace, ""
            ),
            f"tests/{solver_stem}_fixture_parameters.h": CodeParameters.output_parameters_h(
                fixture_stem, fixture_namespace
            ),
        }
    finally:
        par.set_parval_from_str("fd_order", saved_fd_order)
        extras = par.glb_extras_dict.get("Dendro", {})
        codeparameters = extras.get("CFunction_codeparameters", {})
        role_table = extras.get("CFunction_roles", {})
        for name in added_functions:
            cfc.CFunction_dict.pop(name, None)
            codeparameters.pop(name, None)
            role_table.pop(name, None)
        for name in added_parameters:
            par.glb_code_params_dict.pop(name, None)
        for name in added_fields:
            gri.glb_gridfcs_dict.pop(name, None)
        if not codeparameter_sidecar_existed and not codeparameters:
            extras.pop("CFunction_codeparameters", None)
        if not role_sidecar_existed and not role_table:
            extras.pop("CFunction_roles", None)
        if not dendro_extras_existed and not extras:
            par.glb_extras_dict.pop("Dendro", None)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
