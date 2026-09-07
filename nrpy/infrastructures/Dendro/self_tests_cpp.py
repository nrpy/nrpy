# nrpy/infrastructures/Dendro/self_tests_cpp.py
"""
Emit the generated project's self-test executable.

Every test is driven from the generated headers and the registered CFunctions;
the standalone host types stand in for the pinned Dendrolib API.  No field name and
no asymptotic value is written here: the state is established entirely by the
generated initial-data CFunction and only invariants are asserted.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Tuple

from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.solver_context import (
    substitute_solver_identifiers,
)

BANNER = generated_file_banner()

# The CTest case names, in the order the runner dispatches them.  The tests
# CMake registers one case per entry, so the two cannot drift.
SECTIONS: Tuple[str, ...] = (
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

_TESTS = """// Usage: $STEM_self_tests {state|params|padding|offsets|upwind|rhs|init|names|
//                            detgtrazero|constraints|all}
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

namespace {

using $NAMESPACE::generated::NUM_EVOL_GFS;

// One padded single block, sized from the generated padding constant.
struct TestBlock {
  explicit TestBlock(unsigned blocks = 1)
      : extent(2 * $NAMESPACE::generated::REQUIRED_PADDING + 3),
        vol(static_cast<std::size_t>(extent) * extent * extent),
        store(NUM_EVOL_GFS,
              std::vector<$SCALAR>(vol * blocks, $SCALAR{0})),
        rhs(NUM_EVOL_GFS,
            std::vector<$SCALAR>(vol * blocks, $SCALAR{0})) {
    geom.nx = geom.ny = geom.nz = extent;
    geom.padding = $NAMESPACE::generated::REQUIRED_PADDING;
    geom.component_offset = 0;
    geom.pmin_padded[0] = geom.pmin_padded[1] = geom.pmin_padded[2] = 0.0;
    geom.dx[0] = geom.dx[1] = geom.dx[2] = 0.1;
    geom.boundary_flags = 0;
  }  // END FUNCTION: TestBlock constructor

  std::vector<$SCALAR*> state_pointers() {
    std::vector<$SCALAR*> out(NUM_EVOL_GFS);
    for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) out[f] = store[f].data();
    return out;
  }  // END FUNCTION: state_pointers
  std::vector<const $SCALAR*> const_state_pointers() {
    std::vector<const $SCALAR*> out(NUM_EVOL_GFS);
    for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) out[f] = store[f].data();
    return out;
  }  // END FUNCTION: const_state_pointers
  std::vector<$SCALAR*> rhs_pointers() {
    std::vector<$SCALAR*> out(NUM_EVOL_GFS);
    for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) out[f] = rhs[f].data();
    return out;
  }  // END FUNCTION: rhs_pointers
  std::size_t index(unsigned a, unsigned b, unsigned c) const {
    return geom.component_offset + a + extent * (b + extent * c);
  }

  unsigned extent;
  std::size_t vol;
  std::vector<std::vector<$SCALAR>> store;
  std::vector<std::vector<$SCALAR>> rhs;
  BlockGeometry geom;
};  // END STRUCT: TestBlock

double max_abs_interior(const TestBlock& block) {
  double worst = 0.0;
  const unsigned pad = block.geom.padding;
  for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
    for (unsigned c = pad; c < block.extent - pad; ++c) {
      for (unsigned b = pad; b < block.extent - pad; ++b) {
        for (unsigned a = pad; a < block.extent - pad; ++a) {
          worst = std::max(worst, std::fabs(static_cast<double>(
                                      block.rhs[f][block.index(a, b, c)])));
        }  // END LOOP: for a over interior x
      }  // END LOOP: for b over interior y
    }  // END LOOP: for c over interior z
  }  // END LOOP: for f over evolved components
  return worst;
}  // END FUNCTION: max_abs_interior

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

int test_padding() {
  // The padding must cover the centred radius, plus one more point when the
  // kernel upwinds or carries Kreiss-Oliger dissipation, because both of those
  // families reach one point further.  The four constants come from different
  // registries -- the padding from the reach of the emitted derivative
  // operators, the order and the dissipation switch from the generation
  // parameters, the control-field count from the state registry -- so this
  // compares them rather than restating one.  Re-deriving the reach here would
  // be a second stencil model.
  const unsigned needed =
      $NAMESPACE::generated::FD_ORDER / 2 +
      (($NAMESPACE::generated::NUM_UPWIND_CONTROL_GFS > 0 ||
        $NAMESPACE::generated::KO_ENABLED)
           ? 1u
           : 0u);
  if ($NAMESPACE::generated::REQUIRED_PADDING < needed) {
    std::fprintf(stderr,
                 "FAIL: REQUIRED_PADDING %u is below the reach %u implied by "
                 "FD_ORDER %u, NUM_UPWIND_CONTROL_GFS %u and KO_ENABLED %d\\n",
                 $NAMESPACE::generated::REQUIRED_PADDING, needed,
                 $NAMESPACE::generated::FD_ORDER,
                 $NAMESPACE::generated::NUM_UPWIND_CONTROL_GFS,
                 static_cast<int>($NAMESPACE::generated::KO_ENABLED));
    return 1;
  }  // END IF: padding below stencil reach
  return 0;
}  // END FUNCTION: test_padding

int test_offsets() {
  // The generated bindings add `geom.component_offset`, so a writer aimed at
  // component offset V must fill exactly the second block of a two-block
  // allocation and leave the first untouched.  This exercises the emitted
  // pointer arithmetic rather than restating it.
  //
  // The whole allocation is pre-filled with a sentinel first.  Summing the
  // untouched region is not enough on its own: the initial-data fill writes
  // the registered asymptotic value, which is zero for every field but the
  // lapse and the conformal factor, so a mis-aimed writer would write zeros
  // over an already-zero region and stay invisible for most components.
  // Against a sentinel every component is discriminating, and each is checked
  // on its own rather than through a summed norm.
  TestBlock block(2);
  block.geom.component_offset = block.vol;
  const $SCALAR sentinel = static_cast<$SCALAR>(-7.5);
  for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
    for (std::size_t cell = 0; cell < 2 * block.vol; ++cell) {
      block.store[f][cell] = sentinel;
    }  // END LOOP: for cell over both blocks
  }  // END LOOP: for f over evolved components
  std::vector<$SCALAR*> out = block.state_pointers();
  $MINKOWSKI_INITIAL_DATA_BLOCK(block.geom, out.data());
  for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
    for (std::size_t cell = 0; cell < block.vol; ++cell) {
      if (block.store[f][cell] != sentinel) {
        std::fprintf(stderr,
                     "FAIL: component %u was written outside its own block, "
                     "so the emitted binding dropped geom.component_offset\\n",
                     f);
        return 1;
      }  // END IF: wrote outside its component
      if (block.store[f][block.vol + cell] != $NAMESPACE::generated::EVOL_GF_F_INFINITY[f]) {
        std::fprintf(stderr,
                     "FAIL: component %u was not filled at its component "
                     "offset\\n",
                     f);
        return 1;
      }  // END IF: not filled at the offset
    }  // END LOOP: for cell over the block
  }  // END LOOP: for f over evolved components

  // The same contract binds the right-hand side: its in_ and rhs_ bindings
  // add the same per-component base, so with a nonzero component offset the
  // kernel must leave the first block's right-hand side untouched.  Without
  // this the initial-data writer is the only binding any gate covers.
  //
  // The first block's state is filled with a spatially *varying* decoy the
  // initial-data fill never writes.  Minkowski is a fixed point, so with
  // correct input bindings the right-hand side is identically zero; an input
  // binding that dropped the offset would read the decoy and produce a nonzero
  // result.  The decoy has to vary in space: a constant decoy has vanishing
  // derivatives, and a constant perturbation of a single field leaves the
  // right-hand side zero, so most bindings would stay invisible.
  TestBlock rhs_block(2);
  rhs_block.geom.component_offset = rhs_block.vol;
  for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
    for (std::size_t cell = 0; cell < 2 * rhs_block.vol; ++cell) {
      rhs_block.store[f][cell] = static_cast<$SCALAR>(
          0.05 * static_cast<double>((f + 1) * (cell % 13) + 1));
      rhs_block.rhs[f][cell] = sentinel;
    }  // END LOOP: for cell over both blocks
  }  // END LOOP: for f over evolved components
  std::vector<$SCALAR*> rhs_state = rhs_block.state_pointers();
  $MINKOWSKI_INITIAL_DATA_BLOCK(rhs_block.geom, rhs_state.data());
  $NAMESPACE::generated::params_struct params;
  $PARAMS_STRUCT_SET_TO_DEFAULT(params);
  std::vector<const $SCALAR*> rhs_in = rhs_block.const_state_pointers();
  std::vector<$SCALAR*> rhs_out = rhs_block.rhs_pointers();
  $RHS_EVAL_BLOCK(rhs_block.geom, rhs_in.data(), rhs_out.data()$RHS_EVAL_BLOCK_TAIL);
  for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
    for (std::size_t cell = 0; cell < rhs_block.vol; ++cell) {
      if (rhs_block.rhs[f][cell] != sentinel) {
        std::fprintf(stderr,
                     "FAIL: the right-hand side wrote component %u outside "
                     "its own block, so an emitted output binding dropped "
                     "geom.component_offset\\n",
                     f);
        return 1;
      }  // END IF: RHS wrote outside its component
    }  // END LOOP: for cell over the block
  }  // END LOOP: for f over evolved components
  const unsigned pad = rhs_block.geom.padding;
  for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
    for (unsigned c = pad; c < rhs_block.extent - pad; ++c) {
      for (unsigned b = pad; b < rhs_block.extent - pad; ++b) {
        for (unsigned a = pad; a < rhs_block.extent - pad; ++a) {
          const std::size_t at = rhs_block.vol + rhs_block.index(a, b, c) -
                                 rhs_block.geom.component_offset;
          if (rhs_block.rhs[f][at] != static_cast<$SCALAR>(0)) {
            std::fprintf(stderr,
                         "FAIL: the Minkowski right-hand side of component %u "
                         "is nonzero at a component offset, so an emitted "
                         "input binding dropped geom.component_offset\\n",
                         f);
            return 1;
          }  // END IF: RHS nonzero at offset
        }  // END LOOP: for a over interior x
      }  // END LOOP: for b over interior y
    }  // END LOOP: for c over interior z
  }  // END LOOP: for f over evolved components

  // The flat-block adapter must apply the component offset exactly once.  It
  // forwards its geometry unchanged to the per-block kernel, which applies the
  // offset itself, so an adapter that also applied it would address
  // `base + f*vol + 2*offset`.  Nothing else reaches this: the lifecycle's
  // flat-versus-block comparison zeroes the offset on both sides, so without
  // this check a reintroduced double application passes every gate.
  //
  // Minkowski data is staged at the component offset and the rest of the
  // allocation carries the sentinel, so a correct binding reads a fixed point
  // and must leave an identically zero interior.  Any binding that shifts the
  // per-component base fails, in either direction: the shift leaves at least
  // one asserted window unwritten, holding the sentinel, and moves the rest
  // onto a neighbouring component's slot.  Under the doubled offset this guards
  // against, the first component's window is the unwritten one, and each
  // component reads its neighbour's staged data rather than its own.  The allocation is sized for the doubled reach as well, so a
  // shifted binding fails this check rather than running past the end.
  const std::size_t flat_offset = rhs_block.vol;
  const std::size_t flat_span =
      static_cast<std::size_t>(NUM_EVOL_GFS + 2) * rhs_block.vol;
  std::vector<$SCALAR> flat_in(flat_span, sentinel);
  std::vector<$SCALAR> flat_rhs(flat_span, sentinel);
  for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
    for (std::size_t cell = 0; cell < rhs_block.vol; ++cell) {
      flat_in[flat_offset + static_cast<std::size_t>(f) * rhs_block.vol + cell] =
          rhs_block.store[f][rhs_block.vol + cell];
    }  // END LOOP: for cell over the block
  }  // END LOOP: for f over evolved components
  $RHS_EVAL_FLAT_BLOCK(rhs_block.geom, flat_in.data(), flat_rhs.data()$RHS_EVAL_BLOCK_TAIL);
  for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
    for (unsigned c = pad; c < rhs_block.extent - pad; ++c) {
      for (unsigned b = pad; b < rhs_block.extent - pad; ++b) {
        for (unsigned a = pad; a < rhs_block.extent - pad; ++a) {
          const std::size_t at =
              flat_offset + static_cast<std::size_t>(f) * rhs_block.vol + a +
              rhs_block.extent * (b + rhs_block.extent * c);
          if (flat_rhs[at] != static_cast<$SCALAR>(0)) {
            std::fprintf(stderr,
                         "FAIL: the flat-block adapter's Minkowski right-hand "
                         "side of component %u is wrong at a component offset, "
                         "so it did not apply geom.component_offset exactly "
                         "once\\n",
                         f);
            return 1;
          }  // END IF: flat RHS nonzero at offset
        }  // END LOOP: for a over interior x
      }  // END LOOP: for b over interior y
    }  // END LOOP: for c over interior z
  }  // END LOOP: for f over evolved components
  return 0;
}  // END FUNCTION: test_offsets

int test_upwind() {
  // The generated kernel selects between the forward- and backward-shifted
  // stencils using the registered upwind control fields, whose indices are
  // generated (no field is named here).
  //
  // The probe is two-sided and per axis.  Holding one control sign fixed, it
  // perturbs the cell REQUIRED_PADDING points ahead on the axis, then the
  // cell the same distance behind, and compares how much each moves the
  // interior right-hand side.  Within one control sign the algebraic
  // contribution of the control value and the symmetric Kreiss-Oliger reach
  // are identical ahead and behind, so both cancel in the difference and only
  // the shifted stencil survives.  A positive control must lean on the cell
  // ahead and a negative control on the cell behind.
  //
  // A one-sided or magnitude-only comparison is not sufficient: reversing the
  // selection sense merely exchanges the two sensitivities, and comparing
  // across control signs also sees the control's algebraic contribution.
  if ($NAMESPACE::generated::NUM_UPWIND_CONTROL_GFS == 0) {
    std::fprintf(stderr,
                 "FAIL: no upwind control fields are generated, so stencil "
                 "selection cannot be exercised\\n");
    return 1;
  }  // END IF: no upwind control fields
  $NAMESPACE::generated::params_struct params;
  $PARAMS_STRUCT_SET_TO_DEFAULT(params);
  const unsigned axes =
      ($NAMESPACE::generated::NUM_UPWIND_CONTROL_GFS < 3u)
          ? $NAMESPACE::generated::NUM_UPWIND_CONTROL_GFS
          : 3u;
  for (unsigned axis = 0; axis < axes; ++axis) {
    // sens[sign_index][side]: side 0 is the cell ahead, side 1 behind.
    double sens[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
    for (int sign_index = 0; sign_index < 2; ++sign_index) {
      const double sign = (sign_index == 0) ? 1.0 : -1.0;
      for (int side = 0; side < 2; ++side) {
        TestBlock block;
        std::vector<$SCALAR*> out = block.state_pointers();
        $MINKOWSKI_INITIAL_DATA_BLOCK(block.geom, out.data());
        // A ramp along every axis makes each first derivative nonzero.
        for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
          for (unsigned c = 0; c < block.extent; ++c) {
            for (unsigned b = 0; b < block.extent; ++b) {
              for (unsigned a = 0; a < block.extent; ++a) {
                block.store[f][block.index(a, b, c)] +=
                    static_cast<$SCALAR>(1e-3 * static_cast<double>(a + b + c));
              }  // END LOOP: for a over padded x
            }  // END LOOP: for b over padded y
          }  // END LOOP: for c over padded z
        }  // END LOOP: for f over evolved components
        // Only the control for this axis changes sign; the others stay
        // positive, so a permuted or partly dead selection is exposed too.
        for (unsigned k = 0;
             k < $NAMESPACE::generated::NUM_UPWIND_CONTROL_GFS; ++k) {
          const unsigned f =
              $NAMESPACE::generated::EVOL_UPWIND_CONTROL_INDICES[k];
          const double value = (k == axis) ? sign * 0.5 : 0.5;
          for (std::size_t cell = 0; cell < block.vol; ++cell) {
            block.store[f][cell] = static_cast<$SCALAR>(value);
          }  // END LOOP: for cell over the block
        }  // END LOOP: for k over control fields
        const unsigned p = block.geom.padding;
        const unsigned r = $NAMESPACE::generated::REQUIRED_PADDING;
        const unsigned q = (side == 0) ? p + r : p - r;
        const std::size_t probe = block.index(p, p, p);
        std::size_t moved_cell = block.index(q, p, p);
        if (axis == 1) moved_cell = block.index(p, q, p);
        if (axis == 2) moved_cell = block.index(p, p, q);
        std::vector<const $SCALAR*> in = block.const_state_pointers();
        std::vector<$SCALAR*> rhs = block.rhs_pointers();
        $RHS_EVAL_BLOCK(block.geom, in.data(), rhs.data()$RHS_EVAL_BLOCK_TAIL);
        std::vector<double> before(NUM_EVOL_GFS, 0.0);
        for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
          before[f] = static_cast<double>(block.rhs[f][probe]);
        }  // END LOOP: for f over evolved components
        for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
          block.store[f][moved_cell] += static_cast<$SCALAR>(0.25);
        }  // END LOOP: for f over evolved components
        $RHS_EVAL_BLOCK(block.geom, in.data(), rhs.data()$RHS_EVAL_BLOCK_TAIL);
        for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
          const double moved =
              std::fabs(static_cast<double>(block.rhs[f][probe]) - before[f]);
          if (moved > sens[sign_index][side]) sens[sign_index][side] = moved;
        }  // END LOOP: for f over evolved components
      }  // END LOOP: for side ahead and behind
    }  // END LOOP: for sign_index over control signs
    const double scale = std::fmax(std::fmax(sens[0][0], sens[0][1]),
                                   std::fmax(sens[1][0], sens[1][1]));
    const double margin = 1.0e-6 * scale;
    if (!(scale > 0.0) || !(sens[0][0] - sens[0][1] > margin) ||
        !(sens[1][1] - sens[1][0] > margin)) {
      std::fprintf(stderr,
                   "FAIL: axis %u upwind stencil is not selected in the "
                   "direction the control demands (positive control: ahead "
                   "%.17g vs behind %.17g; negative control: ahead %.17g vs "
                   "behind %.17g)\\n",
                   axis, sens[0][0], sens[0][1], sens[1][0], sens[1][1]);
      return 1;
    }  // END IF: wrong stencil selected
  }  // END LOOP: for axis over upwinded directions
  return 0;
}  // END FUNCTION: test_upwind

int test_rhs() {
  // The Minkowski state is a fixed point: the generated initial-data CFunction
  // fills every EVOL field to its registered asymptotic value and the
  // generated direct-FD RHS must then vanish at every interior cell.
  $NAMESPACE::generated::params_struct params;
  $PARAMS_STRUCT_SET_TO_DEFAULT(params);
  TestBlock block;
  std::vector<$SCALAR*> out = block.state_pointers();
  $MINKOWSKI_INITIAL_DATA_BLOCK(block.geom, out.data());
  std::vector<const $SCALAR*> in = block.const_state_pointers();
  std::vector<$SCALAR*> rhs = block.rhs_pointers();
  $RHS_EVAL_BLOCK(block.geom, in.data(), rhs.data()$RHS_EVAL_BLOCK_TAIL);
  return max_abs_interior(block) <= 1e-13 ? 0 : 2;
}  // END FUNCTION: test_rhs

int test_init() {
  // The generated Minkowski fill writes every EVOL field to its registered
  // asymptotic value at every cell: the result must be non-trivial (some
  // component is nonzero) and must match the generated metadata, so the writer
  // and the state registry agree without naming a field.
  TestBlock block;
  std::vector<$SCALAR*> out = block.state_pointers();
  $MINKOWSKI_INITIAL_DATA_BLOCK(block.geom, out.data());
  double norm = 0.0;
  for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
    for (std::size_t cell = 0; cell < block.vol; ++cell) {
      const double written = static_cast<double>(block.store[f][cell]);
      norm += std::fabs(written);
      if (std::fabs(written - static_cast<double>(
                                  $NAMESPACE::generated::EVOL_GF_F_INFINITY[f])) >
          1e-15) {
        return 2;  // fill disagrees with the generated asymptotic metadata
      }
    }  // END LOOP: for cell over the block
  }  // END LOOP: for f over evolved components
  if (norm <= 0.0) return 3;  // the fill wrote only zeros (not Minkowski)
  return 0;
}  // END FUNCTION: test_init

// One padded block of DIAG storage, sized from the generated DIAG count.
struct DiagBlock {
  explicit DiagBlock(std::size_t vol)
      : store($NAMESPACE::generated::NUM_DIAG_GFS,
              std::vector<$SCALAR>(vol, $SCALAR{0})) {}
  std::vector<$SCALAR*> pointers() {
    std::vector<$SCALAR*> out($NAMESPACE::generated::NUM_DIAG_GFS);
    for (unsigned f = 0; f < $NAMESPACE::generated::NUM_DIAG_GFS; ++f) {
      out[f] = store[f].data();
    }  // END LOOP: for f over diagnostic components
    return out;
  }  // END FUNCTION: DiagBlock::pointers
  std::vector<std::vector<$SCALAR>> store;
};  // END STRUCT: diagnostic block storage

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

int test_detgtrazero() {
  // A flat state is unchanged, a perturbed state is projected onto the
  // determinant and trace constraints, and enforcement is idempotent.  The
  // perturbation is applied through the generated CFunction to every
  // component, so no field is named and no asymptotic value is written here.
  TestBlock block;
  std::vector<$SCALAR*> state = block.state_pointers();
  $MINKOWSKI_INITIAL_DATA_BLOCK(block.geom, state.data());
  std::vector<std::vector<$SCALAR>> flat = block.store;
  $NAMESPACE::generated::detgtrazero_status_struct flat_status;
  $ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_BLOCK(block.geom, state.data(), &flat_status);
  if (flat_status.failed_points != 0) return 1;
  if (block.store != flat) return 2;  // flat state must be untouched
  if (flat_status.max_abs_det_minus_one > 5e-13) return 3;
  if (flat_status.max_abs_trace_residual > 5e-13) return 4;

  $SMOOTH_PERTURBATION_BLOCK(block.geom, state.data(), static_cast<$SCALAR>(1e-2),
                      static_cast<$SCALAR>(1.0));
  $NAMESPACE::generated::detgtrazero_status_struct first;
  $ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_BLOCK(block.geom, state.data(), &first);
  if (first.failed_points != 0) return 5;
  if (first.projected_points == 0) return 6;
  $NAMESPACE::generated::detgtrazero_status_struct second;
  $ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_BLOCK(block.geom, state.data(), &second);
  if (second.failed_points != 0) return 7;
  // After one enforcement pass the constraints hold to the stated tolerance.
  if (second.max_abs_det_minus_one > 5e-13) return 8;
  // Tolerance scale: max(1, largest |evolved component|) over the block.
  double state_scale = 1.0;
  for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
    for (std::size_t cell = 0; cell < block.vol; ++cell) {
      state_scale = std::max(
          state_scale, std::fabs(static_cast<double>(block.store[f][cell])));
    }  // END LOOP: for cell over the block
  }  // END LOOP: for f over evolved components
  if (second.max_abs_trace_residual > 5e-13 * state_scale) return 9;
  // Idempotence: a further pass over an already-projected state moves no
  // value by more than the tolerance.  Bitwise equality is not the claim --
  // the kernel recomputes a cube root and an inverse, so the last ulp may
  // move.
  std::vector<std::vector<$SCALAR>> projected = block.store;
  $NAMESPACE::generated::detgtrazero_status_struct third;
  $ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_BLOCK(block.geom, state.data(), &third);
  if (third.failed_points != 0) return 10;
  double moved = 0.0;
  for (unsigned f = 0; f < NUM_EVOL_GFS; ++f) {
    for (std::size_t cell = 0; cell < block.vol; ++cell) {
      moved = std::max(moved,
                       std::fabs(static_cast<double>(block.store[f][cell]) -
                                 static_cast<double>(projected[f][cell])));
    }  // END LOOP: for cell over the block
  }  // END LOOP: for f over evolved components
  if (moved > 5e-13 * state_scale) return 11;
  if (third.max_abs_det_minus_one > 5e-13) return 12;
  if (third.max_abs_trace_residual > 5e-13 * state_scale) return 13;
  return 0;
}  // END FUNCTION: test_detgtrazero

int test_constraints() {
  // The constraint diagnostics of the Minkowski solution vanish.  The
  // reduction runs over the generated DIAG count, so no diagnostic is named.
  // Pointwise diagnostic correctness is not established here.  What is
  // pinned elsewhere is the DIAG registration and the exact write set of the
  // lowered kernel, in the diagnostic owner's doctests.  The equations layer's
  // trusted dictionaries pin a different construction profile than this
  // solver lowers, so none of these expressions is pinned there.
  if ($NAMESPACE::generated::NUM_DIAG_GFS == 0) return 1;
  TestBlock block;
  std::vector<$SCALAR*> state = block.state_pointers();
  $MINKOWSKI_INITIAL_DATA_BLOCK(block.geom, state.data());
  DiagBlock diag(block.vol);
  std::vector<$SCALAR*> diag_ptr = diag.pointers();
  std::vector<const $SCALAR*> in = block.const_state_pointers();
  $CONSTRAINTS_EVAL_BLOCK(block.geom, in.data(), diag_ptr.data());
  const unsigned pad = block.geom.padding;
  double worst = 0.0;
  for (unsigned f = 0; f < $NAMESPACE::generated::NUM_DIAG_GFS; ++f) {
    for (unsigned c = pad; c < block.extent - pad; ++c) {
      for (unsigned b = pad; b < block.extent - pad; ++b) {
        for (unsigned a = pad; a < block.extent - pad; ++a) {
          worst = std::max(worst, std::fabs(static_cast<double>(
                                      diag.store[f][block.index(a, b, c)])));
        }  // END LOOP: for a over interior x
      }  // END LOOP: for b over interior y
    }  // END LOOP: for c over interior z
  }  // END LOOP: for f over diagnostic components
  return worst <= 1e-13 ? 0 : 2;
}  // END FUNCTION: test_constraints

int run_section(const char* section) {
  if (std::strcmp(section, "state") == 0) return test_state();
  if (std::strcmp(section, "params") == 0) return test_params();
  if (std::strcmp(section, "padding") == 0) return test_padding();
  if (std::strcmp(section, "offsets") == 0) return test_offsets();
  if (std::strcmp(section, "upwind") == 0) return test_upwind();
  if (std::strcmp(section, "rhs") == 0) return test_rhs();
  if (std::strcmp(section, "init") == 0) return test_init();
  if (std::strcmp(section, "names") == 0) return test_names();
  if (std::strcmp(section, "detgtrazero") == 0) return test_detgtrazero();
  if (std::strcmp(section, "constraints") == 0) return test_constraints();
  if (std::strcmp(section, "all") == 0) {
    const int rc = test_state() + test_params() + test_padding() +
                   test_offsets() + test_upwind() + test_rhs() + test_init() +
                   test_names() + test_detgtrazero() + test_constraints();
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


def output_self_tests_cpp(solver_stem: str, solver_namespace: str) -> str:
    """
    Emit the generated project's self-test source.

    :param solver_stem: Lowercase formulation stem for emitted file names.
    :param solver_namespace: Solver namespace.
    :return: The complete C++ source text.

    Doctests:
    >>> from nrpy.infrastructures.Dendro.clang_format_guards import (
    ...     unguarded_end_namespace_markers,
    ... )
    >>> unguarded_end_namespace_markers(_TESTS)
    []
    >>> _TESTS.count("}  // END NAMESPACE:")
    1
    """
    return BANNER + substitute_solver_identifiers(_TESTS, solver_stem, solver_namespace)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
