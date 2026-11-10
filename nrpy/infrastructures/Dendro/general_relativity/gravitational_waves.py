"""
Generate Dendro spin-weight minus-two wave decomposition.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Union, cast

import sympy as sp

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.equations.special_functions import spin_weighted_spherical_harmonics


def register_CFunction_gravitational_waves(
    solver_stem: str,
    *,
    maximum_l_mode_generated: int = 8,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register Lebedev Psi4 mode decomposition.

    The generated function interpolates zipped Psi4 fields onto coordinate
    spheres, computes ``Psi4_lm`` without an extraction-radius factor, and
    accumulates the separate real and imaginary Psi4 L2 norms.

    :param solver_stem: Lowercase solver name used by the generated header.
    :param maximum_l_mode_generated: Largest compiled spin-weighted mode.
    :return: Updated NRPy environment, or ``None`` during task collection.
    :raises ValueError: If the infrastructure or mode range is invalid.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("Infrastructure must be 'Dendro' for wave extraction.")
    if maximum_l_mode_generated < 2:
        raise ValueError("Psi4 extraction requires maximum_l_mode_generated >= 2.")

    scalar_type = gri.DENDRO_SCALAR_TYPE
    phi = sp.Symbol("phi", real=True)
    harmonic_theta = sp.Symbol("harmonic_theta", real=True)
    harmonic_cases: List[str] = []
    for ell in range(2, maximum_l_mode_generated + 1):
        for mode in range(-ell, ell + 1):
            harmonic = spin_weighted_spherical_harmonics.Y(
                -2, ell, mode, harmonic_theta, phi
            )
            assignments = c_codegen(
                [sp.re(harmonic), sp.im(harmonic)],
                ["harmonic_real", "harmonic_imag"],
                include_braces=False,
                enable_simd=False,
                fp_type=str(par.parval_from_str("fp_type")),
                fp_type_alias=scalar_type,
                cse_sorting="none",
                verbose=False,
            )
            harmonic_cases.append(
                f"case {ell * ell + ell + mode}: {{\n{assignments}    break;\n}}  // END BLOCK: harmonic ell={ell} m={mode}"
            )
    switch_body = "\n".join(harmonic_cases)
    decomposition = rf"""if (mesh == nullptr || extraction_radii == nullptr ||
    modes_real == nullptr || modes_imag == nullptr ||
    l2_real == nullptr || l2_imag == nullptr) {{
    throw std::invalid_argument("gravitational_waves received a null array");
}}  // END IF: null input array pointers
if (maximum_l < 2 || maximum_l > {maximum_l_mode_generated}) {{
    throw std::invalid_argument("requested Psi4 mode was not generated");
}}  // END IF: maximum_l outside generated range
const unsigned required_modes = (maximum_l + 1) * (maximum_l + 1);
if (mode_stride < required_modes) {{
    throw std::invalid_argument("Psi4 mode stride is too small");
}}  // END IF: mode stride too small
if (num_radii != 0 &&
    mode_stride > std::numeric_limits<unsigned>::max() / num_radii) {{
    throw std::overflow_error("Psi4 mode array size overflows unsigned");
}}  // END IF: output size overflows unsigned
const unsigned output_size = num_radii * mode_stride;
if (output_size > static_cast<unsigned>(std::numeric_limits<int>::max())) {{
    throw std::overflow_error("Psi4 MPI reduction count exceeds INT_MAX");
}}  // END IF: MPI count exceeds INT_MAX
std::fill(modes_real, modes_real + output_size, 0.0);
std::fill(modes_imag, modes_imag + output_size, 0.0);
std::fill(l2_real, l2_real + num_radii, 0.0);
std::fill(l2_imag, l2_imag + num_radii, 0.0);
if (!mesh->isActive() || num_radii == 0) {{
    return;
}}  // END IF: inactive rank or no radii
if (psi4_real_zipped == nullptr || psi4_imag_zipped == nullptr) {{
    throw std::invalid_argument("gravitational_waves received null Psi4 data");
}}  // END IF: null Psi4 field data
constexpr unsigned num_points = LEBEDEV_025_NUM_PTS;
const Point grid_limits[2] = {{grid_min, grid_max}};
const Point domain_limits[2] = {{domain_min, domain_max}};
std::vector<{scalar_type}> local_real(output_size, 0.0);
std::vector<{scalar_type}> local_imag(output_size, 0.0);
std::vector<{scalar_type}> local_l2_real(num_radii, 0.0);
std::vector<{scalar_type}> local_l2_imag(num_radii, 0.0);
std::vector<{scalar_type}> coordinates(3 * num_points);
std::vector<{scalar_type}> shell_real(num_points);
std::vector<{scalar_type}> shell_imag(num_points);
std::vector<unsigned> valid_real;
std::vector<unsigned> valid_imag;
for (unsigned radius_index = 0; radius_index < num_radii; ++radius_index) {{
    const {scalar_type} radius = extraction_radii[radius_index];
    for (unsigned point = 0; point < num_points; ++point) {{
        const {scalar_type} theta = LEBEDEV_025_THETA[point];
        const {scalar_type} phi = LEBEDEV_025_PHI[point];
        coordinates[3 * point] = extraction_center.x() +
            radius * std::sin(theta) * std::cos(phi);
        coordinates[3 * point + 1] = extraction_center.y() +
            radius * std::sin(theta) * std::sin(phi);
        coordinates[3 * point + 2] = extraction_center.z() +
            radius * std::cos(theta);
    }}  // END LOOP: for point over Lebedev points
    valid_real.clear();
    valid_imag.clear();
    ot::da::interpolateToCoords(
        mesh, psi4_real_zipped, coordinates.data(), 3 * num_points,
        grid_limits, domain_limits, shell_real.data(), valid_real);
    ot::da::interpolateToCoords(
        mesh, psi4_imag_zipped, coordinates.data(), 3 * num_points,
        grid_limits, domain_limits, shell_imag.data(), valid_imag);
    if (valid_real != valid_imag) {{
        throw std::logic_error("Psi4 fields have different interpolation owners");
    }}  // END IF: mismatched Psi4 interpolation owners
    for (const unsigned valid_index : valid_real) {{
        if (!std::isfinite(shell_real[valid_index]) ||
            !std::isfinite(shell_imag[valid_index])) {{
            std::ostringstream message;
            message << "non-finite interpolated Psi4 at extraction radius "
                    << radius << ", point " << valid_index << " ("
                    << coordinates[3 * valid_index] << ", "
                    << coordinates[3 * valid_index + 1] << ", "
                    << coordinates[3 * valid_index + 2] << "): ("
                    << shell_real[valid_index] << ", "
                    << shell_imag[valid_index] << ")";
            std::cerr << "rank " << mesh->getMPIRank() << ": "
                      << message.str() << std::endl;
            throw std::runtime_error(message.str());
        }}  // END IF: non-finite interpolated Psi4 value
        local_l2_real[radius_index] +=
            shell_real[valid_index] * shell_real[valid_index];
        local_l2_imag[radius_index] +=
            shell_imag[valid_index] * shell_imag[valid_index];
    }}  // END LOOP: for valid_index over owned points
    for (unsigned ell = 2; ell <= maximum_l; ++ell) {{
        for (int mode = -static_cast<int>(ell);
             mode <= static_cast<int>(ell); ++mode) {{
            const unsigned mode_index = static_cast<unsigned>(
                static_cast<int>(ell * ell + ell) + mode);
            {scalar_type} integral_real = 0.0;
            {scalar_type} integral_imag = 0.0;
            for (unsigned valid_index : valid_real) {{
          const {scalar_type} theta = LEBEDEV_025_THETA[valid_index];
          const {scalar_type} phi = LEBEDEV_025_PHI[valid_index];
          const {scalar_type} pole_offset =
              32.0 * std::numeric_limits<{scalar_type}>::epsilon();
          const {scalar_type} harmonic_theta = std::clamp(
              theta, pole_offset, M_PI - pole_offset);
                {scalar_type} harmonic_real = 0.0;
                {scalar_type} harmonic_imag = 0.0;
                switch (mode_index) {{
{switch_body}
                default:
                    throw std::logic_error(
                        "missing generated spin-weighted harmonic");
                }}  // END SWITCH: harmonic by mode_index
                const {scalar_type} weight = LEBEDEV_025_WEIGHT[valid_index];
                integral_real += weight *
                    (shell_real[valid_index] * harmonic_real +
                     shell_imag[valid_index] * harmonic_imag);
                integral_imag += weight *
                    (shell_imag[valid_index] * harmonic_real -
                     shell_real[valid_index] * harmonic_imag);
            }}  // END LOOP: for valid_index over Lebedev quadrature
            const unsigned output_index = radius_index * mode_stride + mode_index;
            local_real[output_index] = 4.0 * M_PI * integral_real;
            local_imag[output_index] = 4.0 * M_PI * integral_imag;
        }}  // END LOOP: for mode over -ell..ell
    }}  // END LOOP: for ell over 2..maximum_l
}}  // END LOOP: for radius_index over extraction radii
par::Mpi_Allreduce(
    local_real.data(), modes_real, static_cast<int>(output_size), MPI_SUM,
    mesh->getMPICommunicator());
par::Mpi_Allreduce(
    local_imag.data(), modes_imag, static_cast<int>(output_size), MPI_SUM,
    mesh->getMPICommunicator());
par::Mpi_Allreduce(
    local_l2_real.data(), l2_real, static_cast<int>(num_radii), MPI_SUM,
    mesh->getMPICommunicator());
par::Mpi_Allreduce(
    local_l2_imag.data(), l2_imag, static_cast<int>(num_radii), MPI_SUM,
    mesh->getMPICommunicator());"""
    desc = """Interpolate Psi4 to extraction spheres and project it onto spin-weight -2 harmonics.

@param[in] mesh Dendro mesh that holds Psi4.
@param[in] psi4_real_zipped Zipped real part of Psi4.
@param[in] psi4_imag_zipped Zipped imaginary part of Psi4.
@param[in] extraction_radii Radii of the extraction spheres.
@param num_radii Number of extraction spheres.
@param maximum_l Largest l mode computed.
@param mode_stride Number of mode slots per extraction radius in the output arrays.
@param[in] extraction_center Center of the extraction spheres.
@param[in] grid_min Lower corner of the octree coordinates.
@param[in] grid_max Upper corner of the octree coordinates.
@param[in] domain_min Lower corner of the physical domain.
@param[in] domain_max Upper corner of the physical domain.
@param[out] modes_real Real part of each (l, m) mode at each radius, summed over ranks.
@param[out] modes_imag Imaginary part of each (l, m) mode at each radius, summed over ranks.
@param[out] l2_real Sum of squared real Psi4 samples at each radius, over ranks.
@param[out] l2_imag Sum of squared imaginary Psi4 samples at each radius, over ranks.
"""
    cfunc_type = "void"
    name = "gravitational_waves"
    params = (
        f"const ot::Mesh* mesh, const {scalar_type}* psi4_real_zipped, "
        f"const {scalar_type}* psi4_imag_zipped, "
        f"const {scalar_type}* extraction_radii, unsigned num_radii, "
        "unsigned maximum_l, unsigned mode_stride, "
        "const Point& extraction_center, const Point& grid_min, "
        "const Point& grid_max, const Point& domain_min, "
        f"const Point& domain_max, {scalar_type}* modes_real, "
        f"{scalar_type}* modes_imag, {scalar_type}* l2_real, "
        f"{scalar_type}* l2_imag"
    )
    body = decomposition
    cfc.register_CFunction(
        subdirectory="generated/src/gravitational_waves",
        includes=[
            f"{solver_stem}_defines.h",
            "daUtils.h",
            "lebedev.h",
            "parUtils.h",
            "<algorithm>",
            "<cmath>",
            "<iostream>",
            "<limits>",
            "<sstream>",
            "<stdexcept>",
            "<vector>",
        ],
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )
    return pcg.NRPyEnv()
