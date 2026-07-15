"""
Register azimuthal-symmetry spatial Lagrange interpolation.

This module emits the spatial stage of the numerical-spacetime interpolation
pipeline used by the geodesic integrators. The generated C API evaluates one
spatial position ``(x, y, z)`` against caller-supplied mapped time-slice
payloads. It writes the normally interpolated Cartesian-basis ``g4DD`` values
and temporarily repurposes the existing 40-component ``Gamma4UDD`` output
buffer to store first Cartesian derivatives of the metric.

The derivative bundle is serialized in symmetric metric-component order,
with derivative direction fastest. For metric pairs
``(00, 01, 02, 03, 11, 12, 13, 22, 23, 33)``, each four-value group is
``(d/dt, d/dx, d/dy, d/dz)``. This spatial helper sets every time derivative
to zero and fills the three Cartesian spatial derivatives.

The target spatial position is converted from Cartesian to the dataset's
native coordinates with the generated BHaH inverse map. Metric values and the
two non-azimuthal native derivatives are obtained from one uniform 2D
Lagrange reconstruction. No interpolation in ``phi`` is performed. Instead,
the metric is rotated from the selected stored reference-phi plane to the
target azimuth, and the native phi derivative of its fixed Cartesian
components is obtained analytically by differentiating that rotation.

Finally, the native derivative index is transformed to Cartesian with the
inverse coordinate Jacobian generated directly from
``reference_metric[CoordSystem].Jac_dUrfm_dDCartUD``. Thus nonlinear mappings
such as ``SinhCylindricalv2n2`` inherit their stretching factors from the
reference-metric definition rather than from duplicated handwritten formulas.

The existing payload stride remains unchanged at 53 doubles per point:
three coordinates, ten metric values, and forty stored Christoffel values.
The forty stored Christoffel values are intentionally ignored by this helper.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, List, Tuple, Union, cast

import sympy as sp

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.BHaH.xx_tofrom_Cart import (
    register_CFunction__Cart_to_xx_and_nearest_i0i1i2,
)
from nrpy.reference_metric import reference_metric

AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_DEFINES = r"""
// Constants for azimuthal-symmetry spatial Lagrange interpolation.
#define AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI 3.14159265358979323846264338327950288
#define AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT 10
#define AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT 40
#define AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_METRIC_DERIVATIVE_COMPONENT_COUNT \
  AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT
#define AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_TENSOR_COMPONENT_COUNT \
  (AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT + AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT)
#define AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_RECORD_COMPONENT_COUNT 53
#define AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_MAX_HALF_WIDTH 32
// Status codes returned by the azimuthal-symmetry spatial Lagrange helper.
typedef enum {
  AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS = 0,
  AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET = 1,
  AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL = 2
} azimuthal_symmetry_spatial_lagrange_interp_status; // END ENUM: azimuthal_symmetry_spatial_lagrange_interp_status

// Reusable parsed context for many spatial interpolation calls.
typedef struct {
  double stored_phi_samples[2];
} azimuthal_symmetry_spatial_lagrange_context_struct; // END STRUCT: azimuthal_symmetry_spatial_lagrange_context_struct
"""

_METRIC_COMPONENT_ORDER: Tuple[Tuple[int, int], ...] = (
    (0, 0),
    (0, 1),
    (0, 2),
    (0, 3),
    (1, 1),
    (1, 2),
    (1, 3),
    (2, 2),
    (2, 3),
    (3, 3),
)
_METRIC_COMPONENT_INDEX: Dict[Tuple[int, int], int] = {
    component: idx for idx, component in enumerate(_METRIC_COMPONENT_ORDER)
}


def register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation(
    CoordSystem: str,
    enable_simd: bool,
    project_dir: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the azimuthal-symmetry spatial Lagrange interpolation helper.

    This is the first stage of the numerical-spacetime interpolation pipeline.
    For one requested photon position and one set of already-mapped numerical
    time slices, it computes the interpolated Cartesian-basis metric and its
    Cartesian spatial derivatives on each supplied slice. The existing
    40-component connection output is used as a temporary metric-derivative
    bundle with all time-derivative entries set to zero.

    :param CoordSystem: Coordinate system used by the source dataset; must be
        `"Spherical"`, `"SinhSpherical"`, `"Cylindrical"`,
        `"SinhCylindrical"`, or `"SinhCylindricalv2n2"`.
    :param enable_simd: Whether SIMD helper headers are already available.
    :param project_dir: Destination project directory for copied headers.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If `CoordSystem` is not supported.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import os
    >>> import tempfile
    >>> import nrpy.c_function as cfc
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> cfc.CFunction_dict.clear()
    >>> with tempfile.TemporaryDirectory(dir=os.getcwd()) as project_dir:
    ...     old_cache_home = os.environ.get("XDG_CACHE_HOME")
    ...     _ = os.environ.__setitem__("XDG_CACHE_HOME", project_dir)
    ...     with contextlib.redirect_stdout(io.StringIO()):
    ...         _ = register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation(
    ...             "SinhCylindricalv2n2", enable_simd=True, project_dir=project_dir
    ...         )
    ...         generated = clang_format(
    ...             cfc.CFunction_dict[
    ...                 "azimuthal_symmetry_spatial_lagrange_interpolation__rfm__SinhCylindricalv2n2"
    ...             ].full_function
    ...         )
    ...         _ = validate_strings(
    ...             generated,
    ...             "azimuthal_symmetry_spatial_lagrange_interpolation__rfm__SinhCylindricalv2n2",
    ...             file_ext="c",
    ...         )
    ...     if old_cache_home is None:
    ...         _ = os.environ.pop("XDG_CACHE_HOME", None)
    ...     else:
    ...         _ = os.environ.__setitem__("XDG_CACHE_HOME", old_cache_home)
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    _ = par.register_CodeParameter(
        "int",
        __name__,
        "numerical_spacetime_spatial_interp_order",
        2,
        commondata=True,
        add_to_parfile=True,
        description=(
            "Centered spatial Lagrange interpolation half-width n; the generated "
            "helper uses 2*n+1 stencil points in each interpolation direction."
        ),
    )
    Bdefines_h.register_BHaH_defines(
        "azimuthal_symmetry_spatial_lagrange_interpolation",
        AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_DEFINES,
    )

    # Step 1: Validate high-level code-generation assumptions and determine
    # which native coordinates define the 2D interpolation grid.
    if CoordSystem in ("Spherical", "SinhSpherical"):
        phi_dim, interp_dim0, interp_dim1 = 2, 0, 1
    elif CoordSystem in (
        "Cylindrical",
        "SinhCylindrical",
        "SinhCylindricalv2n2",
    ):
        phi_dim, interp_dim0, interp_dim1 = 1, 0, 2
    else:
        raise ValueError(
            "azimuthal_symmetry_spatial_lagrange_interpolation is currently "
            "implemented only for CoordSystem in ('Spherical', "
            "'SinhSpherical', 'Cylindrical', 'SinhCylindrical', "
            "'SinhCylindricalv2n2'); found "
            f"'{CoordSystem}'."
        )

    # Step 2: Ensure required helper headers are copied into the project.
    if not enable_simd:
        copy_files(
            package="nrpy.helpers",
            filenames_list=["simd_intrinsics.h"],
            project_dir=project_dir,
            subdirectory="intrinsics",
        )
    copy_files(
        package="nrpy.infrastructures.BHaH.interpolation",
        filenames_list=[
            "differentiate_interpolation_lagrange_uniform.h",
            "interpolation_lagrange_uniform.h",
        ],
        project_dir=project_dir,
        subdirectory="./",
    )
    register_CFunction__Cart_to_xx_and_nearest_i0i1i2(CoordSystem)
    cart_to_xx_name = f"Cart_to_xx_and_nearest_i0i1i2__rfm__{CoordSystem}"
    storage_exprs = {
        interp_dim0: "interp_0_storage_stencil[u]",
        interp_dim1: "interp_1_storage_stencil[v]",
        phi_dim: "phi_plane_storage_index",
    }

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<stdint.h>",
        "differentiate_interpolation_lagrange_uniform.h",
        "interpolation_lagrange_uniform.h",
    ]
    # Step 3: Generate the inverse reference-metric Jacobian and pre-expand
    # metric rotations so the emitted C helper stays linear and reviewable.
    inverse_jacobian_assignments = _build_inverse_jacobian_c_code(CoordSystem)
    direct_metric_assignments = "\n".join(
        _emit_wrapped_assignment(
            f"g4dd_rot[{idx}]",
            _build_metric_rotation_terms(mu, nu, "g4dd_ref"),
        )
        for idx, (mu, nu) in enumerate(_METRIC_COMPONENT_ORDER)
    )
    direct_interp_0_derivative_assignments = "\n".join(
        _emit_wrapped_assignment(
            f"g4dd_native_derivatives_rot[{interp_dim0}][{idx}]",
            _build_metric_rotation_terms(mu, nu, "g4dd_interp_0_ref"),
        )
        for idx, (mu, nu) in enumerate(_METRIC_COMPONENT_ORDER)
    )
    direct_phi_derivative_assignments = "\n".join(
        _emit_wrapped_assignment(
            f"g4dd_native_derivatives_rot[{phi_dim}][{idx}]",
            _build_metric_rotation_derivative_terms(mu, nu, "g4dd_ref"),
        )
        for idx, (mu, nu) in enumerate(_METRIC_COMPONENT_ORDER)
    )
    direct_interp_1_derivative_assignments = "\n".join(
        _emit_wrapped_assignment(
            f"g4dd_native_derivatives_rot[{interp_dim1}][{idx}]",
            _build_metric_rotation_terms(mu, nu, "g4dd_interp_1_ref"),
        )
        for idx, (mu, nu) in enumerate(_METRIC_COMPONENT_ORDER)
    )

    prefunc = (
        r"""
/**
 * Map one full-payload storage-grid index triplet to the serialized point order.
 *
 * @param[in] params Generated BHaH parameter struct.
 * @param i0_storage Storage-grid x0 index, including ghost zones.
 * @param i1_storage Storage-grid x1 index, including ghost zones.
 * @param i2_storage Storage-grid x2 index, including ghost zones.
 * @param[out] point_index Serialized zero-based point-record index.
 * @return Interpolation status code.
 */
static int azimuthal_symmetry_spatial_lagrange_point_index_from_full_payload_indices(
    const params_struct *restrict params,
    const int i0_storage,
    const int i1_storage,
    const int i2_storage,
    uint64_t *restrict point_index) {
  const int payload_i0_count = params->Nxx_plus_2NGHOSTS0;
  const int payload_i1_count = params->Nxx_plus_2NGHOSTS1;
  const int payload_i2_count = params->Nxx_plus_2NGHOSTS2;

  if (i0_storage < 0 || i0_storage >= payload_i0_count ||
      i1_storage < 0 || i1_storage >= payload_i1_count ||
      i2_storage < 0 || i2_storage >= payload_i2_count)
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;

  *point_index =
      (uint64_t)i0_storage +
      (uint64_t)payload_i0_count *
          ((uint64_t)i1_storage +
           (uint64_t)payload_i1_count * (uint64_t)i2_storage);
  return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS;
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_point_index_from_full_payload_indices

/**
 * Evaluate the inverse reference-metric coordinate Jacobian at one native point.
 *
 * The generated entries satisfy
 * `inverse_jacobian[native_direction][cartesian_direction] =
 * d xx^native_direction / d xCart^cartesian_direction`.
 *
 * @param[in] params Generated BHaH parameter struct.
 * @param[in] xx Native target coordinates.
 * @param[out] inverse_jacobian Native-with-respect-to-Cartesian Jacobian.
 */
static void azimuthal_symmetry_spatial_lagrange_inverse_jacobian(
    const params_struct *restrict params,
    const REAL xx[3],
    REAL inverse_jacobian[3][3]) {
  const REAL xx0 = xx[0];
  const REAL xx1 = xx[1];
  const REAL xx2 = xx[2];
{inverse_jacobian_assignments}
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_inverse_jacobian

/**
 * Accumulate metric values and two native derivatives from one payload record.
 *
 * The record begins with the ten Cartesian-basis metric components. Its
 * following forty stored Christoffel components are intentionally ignored.
 *
 * @param value_weight Two-dimensional Lagrange value weight.
 * @param interp_0_derivative_weight Derivative weight in the first native interpolation direction.
 * @param interp_1_derivative_weight Derivative weight in the second native interpolation direction.
 * @param[in] tensor_record Serialized payload record beginning at the metric fields.
 * @param[in,out] g4dd_ref Interpolated metric on the selected stored phi plane.
 * @param[in,out] g4dd_interp_0_ref First native interpolation-direction derivative.
 * @param[in,out] g4dd_interp_1_ref Second native interpolation-direction derivative.
 */
static void azimuthal_symmetry_spatial_lagrange_accumulate_metric_record_direct(
    const REAL value_weight,
    const REAL interp_0_derivative_weight,
    const REAL interp_1_derivative_weight,
    const double *restrict tensor_record,
    REAL g4dd_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
    REAL g4dd_interp_0_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
    REAL g4dd_interp_1_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT]) {
  for (int comp = 0;
       comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT;
       comp++) {
    const REAL metric_value = (REAL)tensor_record[comp];
    g4dd_ref[comp] += value_weight * metric_value;
    g4dd_interp_0_ref[comp] += interp_0_derivative_weight * metric_value;
    g4dd_interp_1_ref[comp] += interp_1_derivative_weight * metric_value;
  } // END LOOP: for comp over serialized metric components
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_accumulate_metric_record_direct

/**
 * Rotate the metric and native metric derivatives about the z axis.
 *
 * The two differentiated Lagrange reconstructions rotate as ordinary metric
 * tensors. The native phi derivative is generated analytically by
 * differentiating the same z-axis rotation used for the metric values. No
 * interpolation in phi is performed.
 *
 * @param delta_phi Active rotation angle from stored plane to target azimuth.
 * @param[in] g4dd_ref Metric components on the selected stored phi plane.
 * @param[in] g4dd_interp_0_ref First interpolated native derivative on that plane.
 * @param[in] g4dd_interp_1_ref Second interpolated native derivative on that plane.
 * @param[out] g4dd_rot Metric components at the target azimuth.
 * @param[out] g4dd_native_derivatives_rot Native derivatives at the target azimuth.
 */
static void azimuthal_symmetry_spatial_lagrange_rotate_metric_and_derivatives_about_z(
    const REAL delta_phi,
    const REAL g4dd_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
    const REAL g4dd_interp_0_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
    const REAL g4dd_interp_1_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
    REAL g4dd_rot[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
    REAL g4dd_native_derivatives_rot[3][AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT]) {
  const REAL cos_delta = cos(delta_phi);
  const REAL sin_delta = sin(delta_phi);

  for (int native_direction = 0; native_direction < 3; native_direction++)
    for (int comp = 0;
         comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT;
         comp++)
      g4dd_native_derivatives_rot[native_direction][comp] = 0.0;

{direct_metric_assignments}
{direct_interp_0_derivative_assignments}
{direct_phi_derivative_assignments}
{direct_interp_1_derivative_assignments}
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_rotate_metric_and_derivatives_about_z
""".replace("{inverse_jacobian_assignments}", inverse_jacobian_assignments)
        .replace("{direct_metric_assignments}", direct_metric_assignments)
        .replace(
            "{direct_interp_0_derivative_assignments}",
            direct_interp_0_derivative_assignments,
        )
        .replace(
            "{direct_phi_derivative_assignments}", direct_phi_derivative_assignments
        )
        .replace(
            "{direct_interp_1_derivative_assignments}",
            direct_interp_1_derivative_assignments,
        )
    )

    cfunc_type = "int"
    name = "azimuthal_symmetry_spatial_lagrange_interpolation"
    params = """const azimuthal_symmetry_spatial_lagrange_context_struct *restrict context,
                const commondata_struct *restrict commondata,
                const params_struct *restrict params,
                const REAL x, const REAL y, const REAL z,
                const int num_target_slices,
                const double *const *restrict slice_payloads,
                REAL *restrict g4dd_out,
                REAL *restrict gamma4udd_out"""
    # Step 4: Use placeholder replacement for generated symbols so the raw C
    # body can keep ordinary braces instead of escaped f-string braces.
    body = (
        r"""
  // Step 1: Validate pointers, then convert the target Cartesian point to
  // native coordinates.
  if (context == NULL || commondata == NULL || params == NULL ||
      slice_payloads == NULL || g4dd_out == NULL || gamma4udd_out == NULL ||
      num_target_slices <= 0)
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;

  const REAL xCart[3] = {x, y, z};
  REAL xx_target[3];
  int center_idx[3];
  {cart_to_xx_name}(params, xCart, xx_target, center_idx);

  const REAL target_interp_0 = xx_target[{interp_dim0}];
  const REAL target_interp_1 = xx_target[{interp_dim1}];
  REAL target_phi = xx_target[{phi_dim}];
  const REAL rho_sq = x * x + y * y;
  const REAL origin_epsilon = 1.0e-14;
  const REAL axis_rho_epsilon = 1.0e-14;
  const int n_interp_ghosts = commondata->numerical_spacetime_spatial_interp_order;
  const long int interp_order_long = 2L * (long int)n_interp_ghosts + 1L;
  const REAL phi0 = (REAL)context->stored_phi_samples[0];
  const REAL phi1 = (REAL)context->stored_phi_samples[1];
  REAL phi_delta = phi1 - phi0;

  if (!isfinite((double)target_interp_0) ||
      !isfinite((double)target_interp_1) || !isfinite((double)target_phi))
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;
  if (target_phi >= (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI - (REAL)1.0e-12 &&
      target_phi <= (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI + (REAL)1.0e-12)
    target_phi -= (REAL)(2.0 * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI);
  if (target_phi < (REAL)(-AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI - 1.0e-12) ||
      target_phi > (REAL)(AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI + 1.0e-12))
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;

  // The azimuthal rotation and inverse cylindrical/spherical Jacobians are
  // singular on the symmetry axis, so retain the existing axis rejection.
  if (target_interp_0 <= origin_epsilon ||
      rho_sq <= axis_rho_epsilon * axis_rho_epsilon)
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;
  while (phi_delta <= (REAL)(-AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI))
    phi_delta += (REAL)(2.0 * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI);
  while (phi_delta > (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI)
    phi_delta -= (REAL)(2.0 * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI);
  if (params->Nxx{phi_dim} != 2 || !isfinite((double)phi0) ||
      !isfinite((double)phi1) ||
      fabs((double)(fabs((double)phi_delta) -
                    AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI)) > 1.0e-12)
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;
  if (n_interp_ghosts < 0 ||
      n_interp_ghosts > AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_MAX_HALF_WIDTH ||
      interp_order_long > (long int)params->Nxx_plus_2NGHOSTS{interp_dim0} ||
      interp_order_long > (long int)params->Nxx_plus_2NGHOSTS{interp_dim1})
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;

  // Step 2: Select one stored phi plane and construct value and derivative
  // coefficients on the two-dimensional native stencil.
  const int interp_order = (int)interp_order_long;
  const int phi_plane = target_phi < (REAL)0.0 ? 0 : 1;
  const int phi_plane_storage_index = NGHOSTS + phi_plane;
  const REAL phi_ref = (REAL)context->stored_phi_samples[phi_plane];

  REAL inv_denom[interp_order];
  REAL src_interp_0_stencil[interp_order];
  REAL src_interp_1_stencil[interp_order];
  REAL diffs_interp_0[interp_order];
  REAL diffs_interp_1[interp_order];
  REAL coeff_interp_0[interp_order];
  REAL coeff_interp_1[interp_order];
  REAL derivative_coeff_interp_0[interp_order];
  REAL derivative_coeff_interp_1[interp_order];
  uint64_t mapped_point_index[interp_order][interp_order];
  int interp_0_storage_stencil[interp_order];
  int interp_1_storage_stencil[interp_order];
  const REAL normalization_2d =
      pow((REAL)(params->dxx{interp_dim0} * params->dxx{interp_dim1}),
          -(interp_order - 1));

  for (int u = 0; u < interp_order; u++) {
    const int interp_0_raw = center_idx[{interp_dim0}] + (u - n_interp_ghosts);
    interp_0_storage_stencil[u] = interp_0_raw;
    src_interp_0_stencil[u] =
        (REAL)(params->xxmin{interp_dim0} +
               (((interp_0_raw - NGHOSTS) + 0.5) * params->dxx{interp_dim0}));
  } // END LOOP: for u over first interpolation-dimension stencil nodes
  for (int v = 0; v < interp_order; v++) {
    const int interp_1_raw = center_idx[{interp_dim1}] + (v - n_interp_ghosts);
    interp_1_storage_stencil[v] = interp_1_raw;
    src_interp_1_stencil[v] =
        (REAL)(params->xxmin{interp_dim1} +
               (((interp_1_raw - NGHOSTS) + 0.5) * params->dxx{interp_dim1}));
  } // END LOOP: for v over second interpolation-dimension stencil nodes

  compute_inv_denom(interp_order, inv_denom);
  compute_diffs_xi(
      interp_order, target_interp_0, src_interp_0_stencil, diffs_interp_0);
  compute_diffs_xi(
      interp_order, target_interp_1, src_interp_1_stencil, diffs_interp_1);
  compute_lagrange_basis_coeffs_xi(
      interp_order, inv_denom, diffs_interp_0, coeff_interp_0);
  compute_lagrange_basis_coeffs_xi(
      interp_order, inv_denom, diffs_interp_1, coeff_interp_1);
  compute_lagrange_basis_derivative_coeffs_xi(
      interp_order, inv_denom, diffs_interp_0, derivative_coeff_interp_0);
  compute_lagrange_basis_derivative_coeffs_xi(
      interp_order, inv_denom, diffs_interp_1, derivative_coeff_interp_1);

  // Step 3: Convert each raw storage-grid stencil node to one payload index.
  for (int v = 0; v < interp_order; v++) {
    for (int u = 0; u < interp_order; u++) {
      const int i0_storage = {i0_storage_expr};
      const int i1_storage = {i1_storage_expr};
      const int i2_storage = {i2_storage_expr};
      const int index_status =
          azimuthal_symmetry_spatial_lagrange_point_index_from_full_payload_indices(
              params, i0_storage, i1_storage, i2_storage,
              &mapped_point_index[v][u]);
      if (index_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS)
        return index_status;
    } // END LOOP: for u over first interpolation-dimension payload indices
  } // END LOOP: for v over second interpolation-dimension payload indices

  // Step 4: Evaluate the native-with-respect-to-Cartesian Jacobian once for
  // this target position. It is common to every requested time slice.
  REAL inverse_jacobian[3][3];
  azimuthal_symmetry_spatial_lagrange_inverse_jacobian(
      params, xx_target, inverse_jacobian);
  for (int native_direction = 0; native_direction < 3; native_direction++)
    for (int cartesian_direction = 0; cartesian_direction < 3;
         cartesian_direction++)
      if (!isfinite((double)inverse_jacobian[native_direction][cartesian_direction]))
        return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;

  // Step 5: Interpolate and differentiate each requested time slice.
  for (int which_slice = 0; which_slice < num_target_slices; which_slice++) {
    const double *restrict slice_payload = slice_payloads[which_slice];
    REAL g4dd_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT] = {0};
    REAL g4dd_interp_0_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT] = {0};
    REAL g4dd_interp_1_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT] = {0};
    REAL g4dd_rot[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT];
    REAL g4dd_native_derivatives_rot[3][AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT];

    if (slice_payload == NULL)
      return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;

    for (int v = 0; v < interp_order; v++) {
      for (int u = 0; u < interp_order; u++) {
        const double *restrict tensor_record =
            slice_payload +
            mapped_point_index[v][u] *
                (uint64_t)
                    AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_RECORD_COMPONENT_COUNT +
            3ULL;
        const REAL value_weight =
            normalization_2d * coeff_interp_1[v] * coeff_interp_0[u];
        const REAL interp_0_derivative_weight =
            normalization_2d * coeff_interp_1[v] *
            derivative_coeff_interp_0[u];
        const REAL interp_1_derivative_weight =
            normalization_2d * derivative_coeff_interp_1[v] *
            coeff_interp_0[u];

        azimuthal_symmetry_spatial_lagrange_accumulate_metric_record_direct(
            value_weight, interp_0_derivative_weight,
            interp_1_derivative_weight, tensor_record, g4dd_ref,
            g4dd_interp_0_ref, g4dd_interp_1_ref);
      } // END LOOP: for u over first interpolation-dimension stencil values
    } // END LOOP: for v over second interpolation-dimension stencil values

    // Step 5.a: Rotate the metric and differentiated native reconstruction to
    // the target phi. The native phi derivative comes from d(rotation)/dphi.
    azimuthal_symmetry_spatial_lagrange_rotate_metric_and_derivatives_about_z(
        target_phi - phi_ref, g4dd_ref, g4dd_interp_0_ref,
        g4dd_interp_1_ref, g4dd_rot, g4dd_native_derivatives_rot);

    REAL *restrict g4dd_slice_out =
        &g4dd_out[which_slice *
                  AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT];
    REAL *restrict metric_derivative_slice_out =
        &gamma4udd_out[
            which_slice *
            AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_METRIC_DERIVATIVE_COMPONENT_COUNT];

    // Step 5.b: Transform only the derivative index from native coordinates to
    // Cartesian coordinates. Pack four derivatives per symmetric metric pair:
    // (dt, dx, dy, dz), with dt deliberately set to zero in this spatial stage.
    for (int metric_component = 0;
         metric_component < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT;
         metric_component++) {
      g4dd_slice_out[metric_component] = g4dd_rot[metric_component];
      metric_derivative_slice_out[4 * metric_component] = 0.0;
      for (int cartesian_direction = 0; cartesian_direction < 3;
           cartesian_direction++) {
        REAL cartesian_derivative = 0.0;
        for (int native_direction = 0; native_direction < 3;
             native_direction++)
          cartesian_derivative +=
              inverse_jacobian[native_direction][cartesian_direction] *
              g4dd_native_derivatives_rot[native_direction][metric_component];
        metric_derivative_slice_out[
            4 * metric_component + cartesian_direction + 1] =
            cartesian_derivative;
      } // END LOOP: for cartesian_direction over x, y, z
    } // END LOOP: for metric_component over symmetric g4DD components
  } // END LOOP: for which_slice over requested slice indices

  return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS;
""".replace("{cart_to_xx_name}", cart_to_xx_name)
        .replace("{phi_dim}", str(phi_dim))
        .replace("{interp_dim0}", str(interp_dim0))
        .replace("{interp_dim1}", str(interp_dim1))
        .replace("{i0_storage_expr}", storage_exprs[0])
        .replace("{i1_storage_expr}", storage_exprs[1])
        .replace("{i2_storage_expr}", storage_exprs[2])
    )

    desc = r"""Interpolate the Cartesian metric and its spatial derivatives.

The helper retains the existing numerical payload and output ABI for a quick
metric-derivative test. It reads only the ten metric values from each 53-double
point record and ignores the following forty stored Christoffel values. It
interpolates the metric and analytically differentiates the same native 2D
Lagrange reconstruction in both non-phi directions.

Axisymmetry is handled without phi interpolation. Metric values and the two
interpolated native derivatives are rotated from one stored reference-phi plane
to the target azimuth. The native phi derivative of the fixed Cartesian metric
components is obtained by analytically differentiating that rotation. The
native derivative index is then transformed to Cartesian with the inverse
Jacobian generated from `reference_metric[CoordSystem].Jac_dUrfm_dDCartUD`.

`gamma4udd_out` is temporarily repurposed as a 40-component metric derivative
bundle. For symmetric metric pairs `(00, 01, 02, 03, 11, 12, 13, 22, 23, 33)`,
each four-value group is `(d/dt, d/dx, d/dy, d/dz)`. This spatial helper sets
all ten time derivatives to zero.

@param[in] context Trusted spatial context.
@param[in] commondata Common interpolation parameters.
@param[in] params Generated BHaH grid and reference-metric parameters.
@param x Cartesian x coordinate.
@param y Cartesian y coordinate.
@param z Cartesian z coordinate.
@param num_target_slices Number of mapped slice payload pointers.
@param[in] slice_payloads Mapped ghost-zone-inclusive slice payload pointers.
@param[out] g4dd_out Flat metric output, ten values per slice.
@param[out] gamma4udd_out Temporary metric derivative output, forty values per slice.
@return Interpolation status code.
"""

    cfc.register_CFunction(
        subdirectory=CoordSystem,
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()


# These Python-side helpers generate the inverse reference-metric Jacobian and
# pre-expand metric rotations so the emitted C stays reviewable.


def _indent_c_code(code: str, spaces: int) -> str:
    """Indent generated C code by a fixed number of spaces."""
    prefix = " " * spaces
    return "\n".join(prefix + line if line else line for line in code.splitlines())


def _build_inverse_jacobian_c_code(CoordSystem: str) -> str:
    """
    Generate C assignments for `d xx^A / d xCart^i` from the reference metric.

    Native coordinate symbols remain local aliases `xx0`, `xx1`, and `xx2`.
    All other free symbols are registered reference-metric CodeParameters and
    are emitted as members of `params`.
    """
    rfm = reference_metric[CoordSystem]
    expressions: List[sp.Expr] = []
    output_names: List[str] = []
    for native_direction in range(3):
        for cartesian_direction in range(3):
            expressions.append(
                sp.sympify(
                    rfm.Jac_dUrfm_dDCartUD[native_direction][cartesian_direction]
                )
            )
            output_names.append(
                f"inverse_jacobian[{native_direction}][{cartesian_direction}]"
            )

    substitutions: Dict[sp.Basic, sp.Basic] = {}
    for expression in expressions:
        for symbol in expression.free_symbols:
            symbol_name = str(symbol)
            if symbol_name not in ("xx0", "xx1", "xx2"):
                substitutions[symbol] = sp.Symbol(f"params->{symbol_name}")
    expressions = [expression.xreplace(substitutions) for expression in expressions]

    generated = c_codegen(
        expressions,
        output_names,
        include_braces=False,
        verbose=False,
    ).rstrip()
    return _indent_c_code(generated, 2)


def _rotation_source_terms(output_index: int) -> List[Tuple[int, int, int, int]]:
    """Return sparse source-index terms for one active z-axis rotation row."""
    if output_index == 0:
        return [(0, 1, 0, 0)]
    if output_index == 1:
        return [(1, 1, 1, 0), (2, -1, 0, 1)]
    if output_index == 2:
        return [(1, 1, 0, 1), (2, 1, 1, 0)]
    if output_index == 3:
        return [(3, 1, 0, 0)]
    raise ValueError(f"Unsupported rotated tensor index: {output_index}")


def _format_scaled_source_term(
    source_expr: str,
    coefficient: int,
    cos_power: int,
    sin_power: int,
) -> str:
    """Format one signed trigonometric monomial times a source component."""
    factors: List[str] = []
    coefficient_abs = abs(coefficient)
    if coefficient_abs != 1:
        factors.append(f"{coefficient_abs}.0")
    factors.extend(["cos_delta"] * cos_power)
    factors.extend(["sin_delta"] * sin_power)
    factors.append(source_expr)
    term = " * ".join(factors)
    return f"-{term}" if coefficient < 0 else term


def _emit_wrapped_assignment(lhs: str, terms: List[str]) -> str:
    """Emit one wrapped C assignment from signed monomial terms."""
    if not terms:
        raise ValueError("Cannot emit an assignment from an empty term list.")
    lines = [f"  {lhs} = {terms[0]}"]
    for term in terms[1:]:
        if term.startswith("-"):
            lines.append(f"      - {term[1:]}")
        else:
            lines.append(f"      + {term}")
    lines[-1] += ";"
    return "\n".join(lines)


def _metric_rotation_term_map(mu: int, nu: int) -> Dict[Tuple[int, int, int], int]:
    """Collect polynomial rotation coefficients for one metric component."""
    term_map: Dict[Tuple[int, int, int], int] = {}
    for source_mu, sign_mu, cos_mu, sin_mu in _rotation_source_terms(mu):
        for source_nu, sign_nu, cos_nu, sin_nu in _rotation_source_terms(nu):
            metric_mu, metric_nu = source_mu, source_nu
            if metric_mu > metric_nu:
                metric_mu, metric_nu = metric_nu, metric_mu
            source_idx = _METRIC_COMPONENT_INDEX[(metric_mu, metric_nu)]
            term_key = (source_idx, cos_mu + cos_nu, sin_mu + sin_nu)
            term_map[term_key] = term_map.get(term_key, 0) + sign_mu * sign_nu
    return {key: coefficient for key, coefficient in term_map.items() if coefficient}


def _format_metric_term_map(
    term_map: Dict[Tuple[int, int, int], int], source_array: str
) -> List[str]:
    """Format a collected metric rotation term map for emitted C code."""
    terms: List[str] = []
    for source_idx, cos_power, sin_power in sorted(term_map):
        coefficient = term_map[(source_idx, cos_power, sin_power)]
        terms.append(
            _format_scaled_source_term(
                f"{source_array}[{source_idx}]",
                coefficient,
                cos_power,
                sin_power,
            )
        )
    return terms or ["0.0"]


def _build_metric_rotation_terms(mu: int, nu: int, source_array: str) -> List[str]:
    """Build terms for one rotated serialized metric component."""
    return _format_metric_term_map(_metric_rotation_term_map(mu, nu), source_array)


def _build_metric_rotation_derivative_terms(
    mu: int, nu: int, source_array: str
) -> List[str]:
    """Differentiate one rotated metric component with respect to azimuth."""
    derivative_map: Dict[Tuple[int, int, int], int] = {}
    for (source_idx, cos_power, sin_power), coefficient in _metric_rotation_term_map(
        mu, nu
    ).items():
        if cos_power > 0:
            key = (source_idx, cos_power - 1, sin_power + 1)
            derivative_map[key] = derivative_map.get(key, 0) - coefficient * cos_power
        if sin_power > 0:
            key = (source_idx, cos_power + 1, sin_power - 1)
            derivative_map[key] = derivative_map.get(key, 0) + coefficient * sin_power
    derivative_map = {
        key: coefficient for key, coefficient in derivative_map.items() if coefficient
    }
    return _format_metric_term_map(derivative_map, source_array)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
