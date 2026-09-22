"""
Emit the Dendro-GR host adapter for the generated BSSN right-hand side.

The adapter preserves Dendro-GR's parameter reader, TwoPunctures initial data,
Berger--Oliger mesh evolution, boundary conditions, diagnostics, and output.
Only the legacy BSSN right-hand-side call is replaced by the generated NRPy
block kernel.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Dict

from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner

_DENDRO_GR_CMAKE = r"""option(NRPY_BSSN_BUILD_DENDRO_GR_DRIVER
       "Build the Dendro-GR BSSN application with the generated NRPy RHS"
       OFF)

if(NRPY_BSSN_BUILD_DENDRO_GR_DRIVER)
  if(NOT TARGET bssn_common)
    message(FATAL_ERROR
      "nrpy_bssnSolver requires the Dendro-GR BSSN_GR directory to be added before nrpy_bssn")
  endif()
  if(WITH_CUDA OR BSSN_ENABLE_CUDA)
    message(FATAL_ERROR "nrpy_bssnSolver is qualified only for the CPU Dendro-GR path")
  endif()

  set(_nrpy_bssn_dendro_gr_dir "${CMAKE_SOURCE_DIR}/BSSN_GR")
  foreach(_nrpy_bssn_source
      "${_nrpy_bssn_dendro_gr_dir}/src/bssngr_main.cpp"
      "${_nrpy_bssn_dendro_gr_dir}/src/bssnAEH.cpp"
      "${_nrpy_bssn_dendro_gr_dir}/src/bssnCtx.cpp")
    if(NOT EXISTS "${_nrpy_bssn_source}")
      message(FATAL_ERROR "nrpy_bssnSolver cannot find ${_nrpy_bssn_source}")
    endif()
  endforeach()

  add_executable(nrpy_bssnSolver
    "${_nrpy_bssn_dendro_gr_dir}/src/bssngr_main.cpp"
    "${_nrpy_bssn_dendro_gr_dir}/src/bssnAEH.cpp"
    "${_nrpy_bssn_dendro_gr_dir}/src/bssnCtx.cpp"
    "${BSSN_MODULE_ROOT}/src/bssn_dendro_gr_adapter.cpp"
  )
  set_source_files_properties(
    "${_nrpy_bssn_dendro_gr_dir}/src/bssnCtx.cpp"
    PROPERTIES COMPILE_DEFINITIONS "bssnRHS=nrpy_bssnRHS"
  )
  target_compile_features(nrpy_bssnSolver PRIVATE cxx_std_17)
  target_compile_options(nrpy_bssnSolver PRIVATE -Wall -fext-numeric-literals)
  target_include_directories(nrpy_bssnSolver PRIVATE
    "${_nrpy_bssn_dendro_gr_dir}/include"
    "${_nrpy_bssn_dendro_gr_dir}/scripts"
    "${_nrpy_bssn_dendro_gr_dir}/src"
  )
  target_link_libraries(nrpy_bssnSolver PRIVATE
    nrpy_bssn_dendro
    bssn_common
    dendro5
    MPI::MPI_CXX
    GSL::gsl
    toml11::toml11
    m
    dendro_git_version_and_date
  )
endif()
"""


_ADAPTER = r"""#include "bssn_defines.h"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <type_traits>
#include <vector>

#include "derivs.h"
#include "grDef.h"
#include "parameters.h"
#include "rhs.h"

namespace {

[[noreturn]] void fail_host_profile(const char* reason) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::fprintf(stderr, "rank %d: NRPy BSSN host adapter: %s\n", rank, reason);
  MPI_Abort(MPI_COMM_WORLD, 1);
  std::abort();
}  // END FUNCTION: fail_host_profile

// clang-format off
}  // END NAMESPACE: internal linkage
// clang-format on

/**
 * Evaluate the generated NRPy BSSN RHS through Dendro-GR's production driver.
 *
 * Input and output retain Dendro-GR's component order and full conformal metric.
 * The adapter maps those pointers to NRPy's BSSN order and subtracts the flat
 * Cartesian metric only for the three diagonal hDD input arrays. Dendro-GR
 * retains ownership of all input, output, mesh, AMR, and communication storage.
 *
 * @param[out] legacy_rhs Dendro-GR unzipped RHS component pointers.
 * @param[in] legacy_state Dendro-GR unzipped evolved-state component pointers.
 * @param[in] blocks Local Dendro blocks, including padding and boundary flags.
 * @param num_blocks Number of local blocks.
 * @param time Current stage time; the generated autonomous RHS does not use it.
 * @param[in] legacy_constraints Dendro-GR constraint arrays; unused by the
 *     generated vacuum BSSN equations.
 */
void nrpy_bssnRHS(double** legacy_rhs, const double** legacy_state,
                  const ot::Block* blocks, unsigned int num_blocks,
                  const double time, const double** legacy_constraints) {
  static_assert(std::is_same<DendroScalar, double>::value,
                "Dendro-GR and generated BSSN scalar types must match");
  static_assert(bssn::BSSN_NUM_VARS ==
                nrpy::bssn::generated::NUM_EVOL_GFS,
                "Dendro-GR and generated BSSN field counts must match");
  if (legacy_rhs == nullptr || legacy_state == nullptr ||
      (num_blocks != 0 && blocks == nullptr))
    fail_host_profile("received a null RHS, state, or block pointer");
  if (bssn::BSSN_ELE_ORDER != nrpy::bssn::generated::FD_ORDER)
    fail_host_profile("BSSN_ELE_ORDER does not match generated FD_ORDER");
  if (bssn::BSSN_PADDING_WIDTH != nrpy::bssn::generated::REQUIRED_PADDING)
    fail_host_profile(
        "BSSN_PADDING_WIDTH does not match generated REQUIRED_PADDING");
$CAKO_CHECK

  nrpy::bssn::generated::params_struct params{};
  bssn_params_struct_set_to_default(params);
  params.eta = static_cast<DendroScalar>(bssn::ETA_CONST);
$KO_PARAMETERS
  if (!bssn_params_validate(params))
    fail_host_profile("Dendro-GR values produced invalid generated parameters");

  using nrpy::bssn::generated::EvolVar;
  using nrpy::bssn::generated::to_index;
  constexpr std::array<unsigned, nrpy::bssn::generated::NUM_EVOL_GFS>
      legacy_field_for_generated = {
          bssn::VAR::U_SYMAT0, bssn::VAR::U_SYMAT1, bssn::VAR::U_SYMAT2,
          bssn::VAR::U_SYMAT3, bssn::VAR::U_SYMAT4, bssn::VAR::U_SYMAT5,
          bssn::VAR::U_ALPHA,  bssn::VAR::U_BETA0, bssn::VAR::U_BETA1,
          bssn::VAR::U_BETA2, bssn::VAR::U_CHI,   bssn::VAR::U_SYMGT0,
          bssn::VAR::U_SYMGT1, bssn::VAR::U_SYMGT2, bssn::VAR::U_SYMGT3,
          bssn::VAR::U_SYMGT4, bssn::VAR::U_SYMGT5, bssn::VAR::U_GT0,
          bssn::VAR::U_GT1, bssn::VAR::U_GT2, bssn::VAR::U_K,
          bssn::VAR::U_B0, bssn::VAR::U_B1, bssn::VAR::U_B2};
  constexpr std::array<double, bssn::BSSN_NUM_VARS> falloff = {
      1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
  constexpr std::array<double, bssn::BSSN_NUM_VARS> asymptotic = {
      1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
      0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const Point domain_minimum(bssn::BSSN_COMPD_MIN[0],
                             bssn::BSSN_COMPD_MIN[1],
                             bssn::BSSN_COMPD_MIN[2]);
  const Point domain_maximum(bssn::BSSN_COMPD_MAX[0],
                             bssn::BSSN_COMPD_MAX[1],
                             bssn::BSSN_COMPD_MAX[2]);
  std::vector<DendroScalar> diagonal_h;
  std::array<std::vector<double>, 3> boundary_gradient;

  for (unsigned int block_id = 0; block_id < num_blocks; ++block_id) {
    const ot::Block& block = blocks[block_id];
    block_geometry_struct geometry{};
    geometry.nx = block.getAllocationSzX();
    geometry.ny = block.getAllocationSzY();
    geometry.nz = block.getAllocationSzZ();
    geometry.padding = block.get1DPadWidth();
    geometry.component_offset = 0;
    geometry.dx[0] = block.computeDx(domain_minimum, domain_maximum);
    geometry.dx[1] = block.computeDy(domain_minimum, domain_maximum);
    geometry.dx[2] = block.computeDz(domain_minimum, domain_maximum);
    const auto node = block.getBlockNode();
    geometry.pmin_padded[0] =
        GRIDX_TO_X(node.minX()) - geometry.padding * geometry.dx[0];
    geometry.pmin_padded[1] =
        GRIDY_TO_Y(node.minY()) - geometry.padding * geometry.dx[1];
    geometry.pmin_padded[2] =
        GRIDZ_TO_Z(node.minZ()) - geometry.padding * geometry.dx[2];
    if (geometry.padding != nrpy::bssn::generated::REQUIRED_PADDING ||
        geometry.nx <= 2 * geometry.padding ||
        geometry.ny <= 2 * geometry.padding ||
        geometry.nz <= 2 * geometry.padding)
      fail_host_profile("received incompatible Dendro block geometry");
    const std::size_t nx = geometry.nx;
    const std::size_t ny = geometry.ny;
    const std::size_t nz = geometry.nz;
    if (ny > std::numeric_limits<std::size_t>::max() / nx)
      fail_host_profile("Dendro block volume overflows size_t");
    const std::size_t plane = nx * ny;
    if (nz > std::numeric_limits<std::size_t>::max() / plane)
      fail_host_profile("Dendro block volume overflows size_t");
    const std::size_t volume = plane * nz;
    if (volume > std::numeric_limits<std::size_t>::max() / 3)
      fail_host_profile("Dendro diagonal storage overflows size_t");
    const std::size_t offset = block.getOffset();

    diagonal_h.resize(3 * volume);
    for (std::size_t cell = 0; cell < volume; ++cell) {
      diagonal_h[cell] = legacy_state[bssn::VAR::U_SYMGT0][offset + cell] - 1.0;
      diagonal_h[volume + cell] =
          legacy_state[bssn::VAR::U_SYMGT3][offset + cell] - 1.0;
      diagonal_h[2 * volume + cell] =
          legacy_state[bssn::VAR::U_SYMGT5][offset + cell] - 1.0;
    }  // END LOOP: convert diagonal conformal metric

    std::array<const DendroScalar*,
               nrpy::bssn::generated::NUM_EVOL_GFS> input{};
    std::array<DendroScalar*, nrpy::bssn::generated::NUM_EVOL_GFS> output{};
    for (unsigned generated_field = 0;
         generated_field < nrpy::bssn::generated::NUM_EVOL_GFS;
         ++generated_field) {
      const unsigned legacy_field = legacy_field_for_generated[generated_field];
      input[generated_field] = legacy_state[legacy_field] + offset;
      output[generated_field] = legacy_rhs[legacy_field] + offset;
    }  // END LOOP: map Dendro-GR component order

    input[to_index(EvolVar::hDD00)] = diagonal_h.data();
    input[to_index(EvolVar::hDD11)] = diagonal_h.data() + volume;
    input[to_index(EvolVar::hDD22)] = diagonal_h.data() + 2 * volume;

    bssn_rhs_eval_block(geometry, input.data(), output.data()$RHS_ARGUMENTS);

    const unsigned int boundary_flags = block.getBlkNodeFlag();
    if (boundary_flags != 0) {
      const unsigned int size[3] = {geometry.nx, geometry.ny, geometry.nz};
      const double pmin[3] = {geometry.pmin_padded[0],
                              geometry.pmin_padded[1],
                              geometry.pmin_padded[2]};
      const double pmax[3] = {
          pmin[0] + (geometry.nx - 1) * geometry.dx[0],
          pmin[1] + (geometry.ny - 1) * geometry.dx[1],
          pmin[2] + (geometry.nz - 1) * geometry.dx[2]};
      for (auto& component : boundary_gradient)
        component.resize(volume);
      for (unsigned int field = 0; field < bssn::BSSN_NUM_VARS; ++field) {
        const double* const field_state = legacy_state[field] + offset;
        double* const field_rhs = legacy_rhs[field] + offset;
        deriv_x(boundary_gradient[0].data(), field_state, geometry.dx[0], size,
                boundary_flags);
        deriv_y(boundary_gradient[1].data(), field_state, geometry.dx[1], size,
                boundary_flags);
        deriv_z(boundary_gradient[2].data(), field_state, geometry.dx[2], size,
                boundary_flags);
        bssn_bcs(field_rhs, field_state, boundary_gradient[0].data(),
                 boundary_gradient[1].data(), boundary_gradient[2].data(), pmin, pmax,
                 falloff[field], asymptotic[field], size, boundary_flags);
$BOUNDARY_KO
      }  // END LOOP: impose physical boundaries
    }  // END IF: block touches physical boundary
  }  // END LOOP: evaluate local Dendro blocks
  (void)time;
  (void)legacy_constraints;
}  // END FUNCTION: nrpy_bssnRHS
"""


def output_bssn_host_files(enable_ko: bool) -> Dict[str, str]:
    """
    Emit source and CMake files for the Dendro-GR BSSN application adapter.

    :param enable_ko: Whether the generated BSSN kernel contains KO dissipation.
    :return: Solver-root-relative generated file contents.

    Doctests:
    >>> files = output_bssn_host_files(enable_ko=True)
    >>> sorted(files)
    ['generated/cmake/dendro_gr_host.cmake', 'src/bssn_dendro_gr_adapter.cpp']
    >>> adapter = files['src/bssn_dendro_gr_adapter.cpp']
    >>> all(name in adapter for name in ('U_ALPHA', 'U_CHI', 'U_SYMGT0', 'U_SYMAT5'))
    True
    >>> 'KreissOliger_strength_gauge' in adapter
    True
    >>> 'nrpy_bssnSolver' in files['generated/cmake/dendro_gr_host.cmake']
    True
    """
    if enable_ko:
        cako_check = (
            "  if (bssn::BSSN_CAKO_ENABLED)\n"
            '    fail_host_profile("conformal-factor-scaled KO is not supported by this generated profile");'
        )
        ko_parameters = (
            "  params.KreissOliger_strength_gauge =\n"
            "      static_cast<DendroScalar>(bssn::KO_DISS_SIGMA);\n"
            "  params.KreissOliger_strength_nongauge =\n"
            "      static_cast<DendroScalar>(bssn::KO_DISS_SIGMA);"
        )
        rhs_arguments = (
            ", params.KreissOliger_strength_gauge, "
            "params.KreissOliger_strength_nongauge, params.eta"
        )
        boundary_ko = r"""        ko_deriv_x(boundary_gradient[0].data(), field_state, geometry.dx[0], size,
                   boundary_flags);
        ko_deriv_y(boundary_gradient[1].data(), field_state, geometry.dx[1], size,
                   boundary_flags);
        ko_deriv_z(boundary_gradient[2].data(), field_state, geometry.dx[2], size,
                   boundary_flags);
        for (unsigned int k = geometry.padding;
             k < geometry.nz - geometry.padding; ++k)
          for (unsigned int j = geometry.padding;
               j < geometry.ny - geometry.padding; ++j)
            for (unsigned int i = geometry.padding;
                 i < geometry.nx - geometry.padding; ++i) {
              const bool physical_boundary =
                  ((boundary_flags & (1u << OCT_DIR_LEFT)) &&
                   i == geometry.padding) ||
                  ((boundary_flags & (1u << OCT_DIR_RIGHT)) &&
                   i == geometry.nx - geometry.padding - 1) ||
                  ((boundary_flags & (1u << OCT_DIR_DOWN)) &&
                   j == geometry.padding) ||
                  ((boundary_flags & (1u << OCT_DIR_UP)) &&
                   j == geometry.ny - geometry.padding - 1) ||
                  ((boundary_flags & (1u << OCT_DIR_BACK)) &&
                   k == geometry.padding) ||
                  ((boundary_flags & (1u << OCT_DIR_FRONT)) &&
                   k == geometry.nz - geometry.padding - 1);
              if (physical_boundary) {
                const std::size_t cell =
                    i + static_cast<std::size_t>(geometry.nx) *
                            (j + static_cast<std::size_t>(geometry.ny) * k);
                field_rhs[cell] +=
                    bssn::KO_DISS_SIGMA *
                    (boundary_gradient[0][cell] +
                     boundary_gradient[1][cell] +
                     boundary_gradient[2][cell]);
              }  // END IF: point lies on physical boundary
            }  // END LOOP: add KO at physical boundary"""
    else:
        cako_check = (
            "  if (bssn::KO_DISS_SIGMA != 0.0)\n"
            '    fail_host_profile("KO_DISS_SIGMA must be zero for a no-KO generated profile");'
        )
        ko_parameters = ""
        rhs_arguments = ", params.eta"
        boundary_ko = ""

    adapter = generated_file_banner() + _ADAPTER.replace(
        "$CAKO_CHECK", cako_check
    ).replace("$KO_PARAMETERS", ko_parameters).replace(
        "$RHS_ARGUMENTS", rhs_arguments
    ).replace(
        "$BOUNDARY_KO", boundary_ko
    )
    return {
        "generated/cmake/dendro_gr_host.cmake": generated_file_banner("#")
        + _DENDRO_GR_CMAKE,
        "src/bssn_dendro_gr_adapter.cpp": adapter,
    }


if __name__ == "__main__":
    import doctest

    raise SystemExit(doctest.testmod().failed)
