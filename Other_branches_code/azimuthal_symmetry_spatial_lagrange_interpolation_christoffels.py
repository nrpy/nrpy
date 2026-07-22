"""
Register azimuthal-symmetry spatial Lagrange interpolation.

This module emits the spatial stage of the numerical-spacetime interpolation
pipeline used by the geodesic integrators. The generated C API evaluates one
spatial position `(x, y, z)` against caller-supplied mapped time-slice
payloads and flat output buffers for the interpolated Cartesian-basis `g4DD`
and `Gamma4UDD` components. Higher-level code applies this helper once per
time slice in the active temporal stencil, then passes the resulting per-slice
tensors to the temporal interpolation stage.

The target spatial position is converted from Cartesian to the dataset's
native coordinates with the generated BHaH inverse map. The interpolation
coefficients are built on the raw uniform two-dimensional native stencil
required by `interpolation_lagrange_uniform.h`. However, payload reads now use the full
ghost-zone-inclusive logical/storage grid directly. The stored ghost-zone point
records are treated as authoritative Cartesian-basis tensor data at those
ghost-zone grid locations, so this helper no longer remaps raw stencil nodes
back into the stored interior through explicit spherical reflections.

This helper intentionally performs no interpolation in `phi`. Its contract is
dataset-specific: the combined numerical container stores exactly two native
phi planes, and this helper interpolates only on the uniform two-dimensional
native grid. It selects one stored phi plane, performs the weighted native
interpolation there, and finally rotates the interpolated Cartesian-basis
tensors to the target azimuth. In other words, axisymmetry is enforced through
tensor rotation rather than through a separate interpolation in `phi`. This
means the caller must provide metadata for exactly those two stored phi planes,
and the native spatial grid spacing in the two interpolation directions must be
constant across the stored payload.

This module does not implement JSON parsing or mmap ownership for the combined
metadata. The caller must parse the combined file elsewhere, map the required
time-window payloads, and pass already-mapped slice payload pointers to the
interpolation routine.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, List, Tuple, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.BHaH.xx_tofrom_Cart import (
    register_CFunction__Cart_to_xx_and_nearest_i0i1i2,
)

AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_DEFINES = r"""
// Constants for azimuthal-symmetry spatial Lagrange interpolation.
#define AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI 3.14159265358979323846264338327950288
#define AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT 10
#define AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT 40
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
_GAMMA_COMPONENT_ORDER: Tuple[Tuple[int, int, int], ...] = tuple(
    (alpha, mu, nu) for alpha in range(4) for mu in range(4) for nu in range(mu, 4)
)
_METRIC_COMPONENT_INDEX: Dict[Tuple[int, int], int] = {
    component: idx for idx, component in enumerate(_METRIC_COMPONENT_ORDER)
}
_GAMMA_COMPONENT_INDEX: Dict[Tuple[int, int, int], int] = {
    component: idx for idx, component in enumerate(_GAMMA_COMPONENT_ORDER)
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
    time slices, it computes spatially interpolated Cartesian-basis tensors on
    each supplied slice. The companion temporal interpolation helper then
    consumes those per-slice tensors and interpolates them in physical time.

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
    ...             "Spherical", enable_simd=True, project_dir=project_dir
    ...         )
    ...         generated = clang_format(
    ...             cfc.CFunction_dict[
    ...                 "azimuthal_symmetry_spatial_lagrange_interpolation__rfm__Spherical"
    ...             ].full_function
    ...         )
    ...         _ = validate_strings(
    ...             generated,
    ...             "azimuthal_symmetry_spatial_lagrange_interpolation__rfm__Spherical",
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
        filenames_list=["interpolation_lagrange_uniform.h"],
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
        "<stdint.h>",
        "interpolation_lagrange_uniform.h",
    ]
    # Step 3: Precompute wrapped direct-rotation assignments once in Python so
    # the emitted C helper stays linear and reviewable.
    direct_metric_assignments = "\n".join(
        _emit_wrapped_assignment(
            f"g4dd_rot[{idx}]",
            _build_metric_rotation_terms(mu, nu),
        )
        for idx, (mu, nu) in enumerate(_METRIC_COMPONENT_ORDER)
    )
    direct_gamma_assignments = "\n".join(
        _emit_wrapped_assignment(
            f"gamma_rot[{idx}]",
            _build_gamma_rotation_terms(alpha, mu, nu),
        )
        for idx, (alpha, mu, nu) in enumerate(_GAMMA_COMPONENT_ORDER)
    )
    reflected_metric_signs = ", ".join(
        str(_pi_rotation_axis_sign(mu) * _pi_rotation_axis_sign(nu))
        for mu, nu in _METRIC_COMPONENT_ORDER
    )
    reflected_gamma_signs = ", ".join(
        str(
            _pi_rotation_axis_sign(alpha)
            * _pi_rotation_axis_sign(mu)
            * _pi_rotation_axis_sign(nu)
        )
        for alpha, mu, nu in _GAMMA_COMPONENT_ORDER
    )

    prefunc = r"""
/**
 * Map one full-payload storage-grid index triplet to the serialized point order.
 *
 * The combined numerical payload stores one point record for every storage-grid
 * location, including ghost zones, in i2-major, i0-fast order. For the current
 * writer contract the raw stencil indices are already storage indices, so this
 * helper only validates bounds and converts them to the serialized payload
 * point-record index.
 *
 * @param[in] params Generated BHaH parameter struct.
 * @param i0_storage Storage-grid x0 index, including ghost zones.
 * @param i1_storage Storage-grid x1 index, including ghost zones.
 * @param i2_storage Storage-grid x2 index, including ghost zones.
 * @param[out] point_index Serialized zero-based point-record index.
 * @return `AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS` if the storage
 * indices are inside the full payload, otherwise
 * `AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL`.
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
 * Accumulate one serialized tensor record onto the selected stored phi plane.
 *
 * The ghost-zone-aware payload already stores Cartesian-basis tensor records at
 * the storage-grid location that the raw stencil node requested. After the
 * direct payload read no additional per-node axisymmetry transform is needed
 * before weighted accumulation on the selected stored phi plane.
 *
 * @param weight Weighted 2D Lagrange coefficient for this stencil node.
 * @param[in] tensor_record Serialized payload record beginning at the metric fields.
 * @param[in,out] tensor_ref Spatial interpolation accumulator on the
 * selected common stored phi plane.
 */
static void azimuthal_symmetry_spatial_lagrange_accumulate_tensor_record_direct(
    const REAL weight,
    const double *restrict tensor_record,
    REAL tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_TENSOR_COMPONENT_COUNT]) {
  const double *restrict gamma_record =
      tensor_record + AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT;

  for (int comp = 0;
       comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT;
       comp++)
    tensor_ref[comp] += weight * (REAL)tensor_record[comp];
  for (int comp = 0;
       comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT;
       comp++)
    tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT + comp] +=
        weight * (REAL)gamma_record[comp];
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_accumulate_tensor_record_direct

/**
 * Rotate Cartesian-basis metric and Christoffel components about the z axis.
 *
 * The payload stores tensors on one of exactly two native reference azimuthal
 * planes. This helper does not interpolate in `phi`; instead, axisymmetry
 * recovers the target azimuth by an active spatial rotation through
 * `delta_phi` after the two-dimensional native interpolation has completed.
 * This keeps the interpolation strictly two-dimensional in the stored data
 * while still returning tensors at the requested azimuth, embedded in 4D so
 * that:
 *
 *   x'^i = Lambda^i_j x^j
 *   g'_{mu nu} = (Lambda^T)^alpha_mu (Lambda^T)^beta_nu g_{alpha beta}
 *   Gamma'^alpha_{mu nu} =
 *       Lambda^alpha_beta (Lambda^T)^gamma_mu (Lambda^T)^delta_nu
 *       Gamma^beta_{gamma delta}
 *
 * The time and z directions are unchanged, while the x-y block is the usual
 * planar rotation matrix.
 *
 * @param delta_phi Active rotation angle from stored plane to target azimuth.
 * @param[in] g4dd_ref Upper-triangular metric components on the stored plane.
 * @param[in] gamma_ref Serialized Christoffel components on the stored plane.
 * @param[out] g4dd_rot Rotated upper-triangular metric components.
 * @param[out] gamma_rot Rotated serialized Christoffel components.
 */
static void azimuthal_symmetry_spatial_lagrange_rotate_about_z(
    const REAL delta_phi,
    const REAL g4dd_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
    const REAL gamma_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT],
    REAL g4dd_rot[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
    REAL gamma_rot[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT]) {
  const REAL cos_delta = cos(delta_phi);
  const REAL sin_delta = sin(delta_phi);

  // Step 1: Rotate only the serialized outputs required by the caller.
{direct_metric_assignments}
{direct_gamma_assignments}
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_rotate_about_z
""".replace(
        "{direct_metric_assignments}",
        direct_metric_assignments,
    ).replace(
        "{direct_gamma_assignments}",
        direct_gamma_assignments,
    )
    if CoordSystem in ("Cylindrical", "SinhCylindrical", "SinhCylindricalv2n2"):
        prefunc += r"""

/**
 * Accumulate one serialized tensor record after a rho reflection in cylindrical
 * coordinates.
 *
 * Lower-rho stencil nodes are mirrored into the stored rho >= 0 region, which
 * shifts the native azimuth by pi. In Cartesian components that half turn about
 * the z axis reduces to fixed sign flips on the serialized tensor entries.
 *
 * @param weight Weighted 2D Lagrange coefficient for this stencil node.
 * @param[in] tensor_record Serialized payload record beginning at the metric fields.
 * @param[in,out] tensor_ref Spatial interpolation accumulator on the
 * selected common stored phi plane.
 */
static void azimuthal_symmetry_spatial_lagrange_accumulate_tensor_record_reflected_pi(
    const REAL weight,
    const double *restrict tensor_record,
    REAL tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_TENSOR_COMPONENT_COUNT]) {
  const double *restrict gamma_record =
      tensor_record + AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT;
  static const int8_t metric_sign[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT] = {
      {reflected_metric_signs}};
  static const int8_t gamma_sign[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT] = {
      {reflected_gamma_signs}};

  for (int comp = 0;
       comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT;
       comp++)
    tensor_ref[comp] +=
        weight * (REAL)metric_sign[comp] * (REAL)tensor_record[comp];
  for (int comp = 0;
       comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT;
       comp++)
    tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT + comp] +=
        weight * (REAL)gamma_sign[comp] * (REAL)gamma_record[comp];
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_accumulate_tensor_record_reflected_pi
""".replace(
            "{reflected_metric_signs}",
            reflected_metric_signs,
        ).replace(
            "{reflected_gamma_signs}",
            reflected_gamma_signs,
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
    payload_index_setup = (
        """      const int reflect_rho = reflected_rho_u[u];
      const int i0_storage =
          reflect_rho ? (2 * NGHOSTS - 1 - interp_0_storage_stencil[u])
                      : interp_0_storage_stencil[u];
      const int i1_storage =
          reflect_rho ? phi_plane_storage_index_reflected
                      : phi_plane_storage_index;
      const int i2_storage = interp_1_storage_stencil[v];"""
        if CoordSystem in ("Cylindrical", "SinhCylindrical", "SinhCylindricalv2n2")
        else f"""      const int i0_storage = {storage_exprs[0]};
      const int i1_storage = {storage_exprs[1]};
      const int i2_storage = {storage_exprs[2]};"""
    )
    accumulation_call = (
        """        if (reflected_rho_u[u]) {
          azimuthal_symmetry_spatial_lagrange_accumulate_tensor_record_reflected_pi(
              weight, tensor_record, tensor_ref);
        } else {
          azimuthal_symmetry_spatial_lagrange_accumulate_tensor_record_direct(
              weight, tensor_record, tensor_ref);
        } // END ELSE: rho reflection was not needed for this stencil node"""
        if CoordSystem in ("Cylindrical", "SinhCylindrical", "SinhCylindricalv2n2")
        else """        azimuthal_symmetry_spatial_lagrange_accumulate_tensor_record_direct(
            weight, tensor_record, tensor_ref);"""
    )
    # Step 4: Use placeholder replacement for the inverse-map symbol so the
    # raw C body can keep ordinary braces instead of escaped f-string braces.
    body = (
        r"""
  // Step 1: Bind the trusted parsed context.
  const REAL xCart[3] = {x, y, z};
  REAL xx_target[3];
  int center_idx[3];

  // Step 2: Convert the target Cartesian point to native coordinates.
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
  int phi_plane;
  int phi_plane_storage_index;
  REAL phi_ref = 0.0;
  REAL phi_delta = phi1 - phi0;

  if (!isfinite((double)target_interp_0) || !isfinite((double)target_interp_1) ||
      !isfinite((double)target_phi))
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;
  if (target_phi >= (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI - (REAL)1.0e-12 &&
      target_phi <= (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI + (REAL)1.0e-12)
    target_phi -= (REAL)(2.0 * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI);
  if (target_phi < (REAL)(-AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI - 1.0e-12) ||
      target_phi > (REAL)(AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI + 1.0e-12))
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;

  // Reject the origin/axis because the azimuth is degenerate and this helper
  // recovers phi through rotation, not through a separate interpolation.
  if (target_interp_0 <= origin_epsilon ||
      rho_sq <= axis_rho_epsilon * axis_rho_epsilon)
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;
  while (phi_delta <= (REAL)(-AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI))
    phi_delta += (REAL)(2.0 * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI);
  while (phi_delta > (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI)
    phi_delta -= (REAL)(2.0 * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI);
  if (params->Nxx{phi_dim} != 2 || !isfinite((double)phi0) || !isfinite((double)phi1) ||
      fabs((double)(fabs((double)phi_delta) -
                    AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI)) > 1.0e-12)
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;
  if (n_interp_ghosts < 0 ||
      n_interp_ghosts > AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_MAX_HALF_WIDTH ||
      interp_order_long > (long int)params->Nxx_plus_2NGHOSTS{interp_dim0} ||
      interp_order_long > (long int)params->Nxx_plus_2NGHOSTS{interp_dim1})
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;

  // Step 3: Select one of the two stored phi planes and build 2D interpolation
  // coefficients. No interpolation in phi is performed anywhere in this helper;
  // phi is handled later by rotating tensors from the stored plane to the
  // target azimuth.
  const int interp_order = (int)interp_order_long;
  phi_plane = target_phi < (REAL)0.0 ? 0 : 1;
  phi_plane_storage_index = NGHOSTS + phi_plane;
  phi_ref = (REAL)context->stored_phi_samples[phi_plane];{reflected_phi_setup}
  REAL inv_denom[interp_order];
  REAL src_interp_0_stencil[interp_order];
  REAL src_interp_1_stencil[interp_order];
  REAL diffs_interp_0[interp_order];
  REAL diffs_interp_1[interp_order];
  REAL coeff_interp_0[interp_order];
  REAL coeff_interp_1[interp_order];
  uint64_t mapped_point_index[interp_order][interp_order];
  int interp_0_storage_stencil[interp_order];
  int interp_1_storage_stencil[interp_order];{reflected_node_storage}
  const REAL normalization_2d =
      pow((REAL)(params->dxx{interp_dim0} * params->dxx{interp_dim1}),
          -(interp_order - 1));

  for (int u = 0; u < interp_order; u++) {
    const int interp_0_raw = center_idx[{interp_dim0}] + (u - n_interp_ghosts);
    interp_0_storage_stencil[u] = interp_0_raw;{reflected_node_setup}
    // Keep the polynomial nodes on raw ghost-zone-inclusive storage-grid
    // indices.
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

  // Step 4: Convert each raw storage-grid stencil node to one payload point index.
  for (int v = 0; v < interp_order; v++) {
    for (int u = 0; u < interp_order; u++) {
{payload_index_setup}
      const int index_status =
          azimuthal_symmetry_spatial_lagrange_point_index_from_full_payload_indices(
              params,
              i0_storage,
              i1_storage,
              i2_storage,
              &mapped_point_index[v][u]);
      if (index_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS)
        return index_status;
    } // END LOOP: for u over first interpolation-dimension stencil nodes during payload-index precompute
  } // END LOOP: for v over second interpolation-dimension stencil nodes during payload-index precompute

  // Step 5: Interpolate each requested time slice independently.
  for (int which_slice = 0; which_slice < num_target_slices; which_slice++) {
    const double *restrict slice_payload = slice_payloads[which_slice];
    REAL tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_TENSOR_COMPONENT_COUNT] = {0};

    // Step 5.a: Read ghost-zone-aware payload records and directly accumulate
    // the weighted stencil nodes for this slice.
    for (int v = 0; v < interp_order; v++) {
      for (int u = 0; u < interp_order; u++) {
        const double *restrict tensor_record =
            slice_payload +
            mapped_point_index[v][u] *
                (uint64_t)
                    AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_RECORD_COMPONENT_COUNT +
            3ULL;
        const REAL weight =
            normalization_2d * coeff_interp_1[v] * coeff_interp_0[u];

{accumulation_call}
      } // END LOOP: for u over first interpolation-dimension stencil nodes
    } // END LOOP: for v over second interpolation-dimension stencil nodes

    // Step 5.b: Rotate the interpolated Cartesian-basis tensors to the target azimuth.
    azimuthal_symmetry_spatial_lagrange_rotate_about_z(
        target_phi - phi_ref,
        &tensor_ref[0],
        &tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
        &g4dd_out[which_slice * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
        &gamma4udd_out[which_slice * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT]);
  } // END LOOP: for which_slice over requested slice indices

  return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS;
""".replace("{cart_to_xx_name}", cart_to_xx_name)
        .replace("{phi_dim}", str(phi_dim))
        .replace("{interp_dim0}", str(interp_dim0))
        .replace("{interp_dim1}", str(interp_dim1))
        .replace(
            "{reflected_phi_setup}",
            (
                "\n  const int phi_plane_storage_index_reflected = NGHOSTS + "
                "(1 - phi_plane);"
                if CoordSystem
                in ("Cylindrical", "SinhCylindrical", "SinhCylindricalv2n2")
                else "\n"
            ),
        )
        .replace(
            "{reflected_node_storage}",
            (
                "\n  int reflected_rho_u[interp_order];"
                if CoordSystem
                in ("Cylindrical", "SinhCylindrical", "SinhCylindricalv2n2")
                else ""
            ),
        )
        .replace(
            "{reflected_node_setup}",
            (
                "\n    reflected_rho_u[u] = interp_0_raw < NGHOSTS;"
                if CoordSystem
                in ("Cylindrical", "SinhCylindrical", "SinhCylindricalv2n2")
                else ""
            ),
        )
        .replace(
            "{payload_index_setup}",
            payload_index_setup,
        )
        .replace(
            "{accumulation_call}",
            accumulation_call,
        )
    )

    desc = r"""Interpolate Cartesian geodesic tensors at one spatial position.

The caller supplies a trusted spatial context and already-mapped time-slice
payload pointers. The helper builds one native two-dimensional stencil, converts
each raw storage-grid stencil node to the corresponding full-payload point
record index once, reads tensor components directly from mapped payload memory
for each requested slice, interpolates on the selected stored phi plane,
rotates to the target azimuth, and writes flat per-slice outputs.

This helper assumes the stored payload has constant native grid spacing in
the two non-phi interpolation directions, and that axisymmetry is represented
by exactly two stored phi planes. It therefore performs no interpolation in
`phi`: it interpolates only on the uniform two-dimensional native stencil. The
payload is expected to include
ghost-zone point records, and those records are treated as authoritative
Cartesian-basis tensor values at their storage-grid locations. The final
interpolated Cartesian-basis tensors are then rotated from the selected stored
reference plane to the target azimuth.

@param[in] context Trusted spatial context.
@param[in] commondata Common runtime parameters.
@param[in] params Generated BHaH grid parameters.
@param x Cartesian x coordinate.
@param y Cartesian y coordinate.
@param z Cartesian z coordinate.
@param num_target_slices Number of mapped slice payload pointers.
@param[in] slice_payloads Mapped slice payload pointers.
@param[out] g4dd_out Flat metric output.
@param[out] gamma4udd_out Flat Christoffel output.
@return `AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS` on success,
`AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET` for invalid target
coordinates, or
`AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL` for
unsupported stencil or payload geometry.

@note Each `slice_payloads` entry must point to the beginning of one mapped
ghost-zone-inclusive 3D-grid payload and remain valid for the duration of this
call. The spatial stencil half-width is read from
`commondata->numerical_spacetime_spatial_interp_order`; the actual number of
Lagrange nodes per interpolation direction is `2*n+1`.
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


# These Python-side helpers pre-expand wrapped tensor-rotation assignments so
# the emitted C stays reviewable without changing the formulas.


def _rotation_source_terms(output_index: int) -> List[Tuple[int, int, int, int]]:
    """
    Return sparse source-index coefficients for one rotated Cartesian index.

    Each returned tuple is `(source_index, sign, cos_power, sin_power)`, using
    the active z-axis rotation convention emitted into the generated C helper.

    :param output_index: Destination tensor index in `(t, x, y, z)` ordering.
    :return: Sparse source-index contributions for the requested output index.
    :raises ValueError: If `output_index` falls outside the 4D tensor range.
    """
    if output_index == 0:
        return [(0, 1, 0, 0)]
    if output_index == 1:
        return [(1, 1, 1, 0), (2, -1, 0, 1)]
    if output_index == 2:
        return [(1, 1, 0, 1), (2, 1, 1, 0)]
    if output_index == 3:
        return [(3, 1, 0, 0)]
    raise ValueError(f"Unsupported rotated tensor index: {output_index}")


def _pi_rotation_axis_sign(output_index: int) -> int:
    """
    Return the z-axis half-turn sign for one Cartesian tensor index.

    :param output_index: Tensor index in `(t, x, y, z)` ordering.
    :return: `+1` for unchanged axes and `-1` for the x-y plane directions.
    :raises ValueError: If `output_index` falls outside the 4D tensor range.
    """
    if output_index in (0, 3):
        return 1
    if output_index in (1, 2):
        return -1
    raise ValueError(f"Unsupported pi-rotation tensor index: {output_index}")


def _format_scaled_source_term(
    source_expr: str,
    coefficient: int,
    cos_power: int,
    sin_power: int,
) -> str:
    """
    Format one signed polynomial source term for emitted C code.

    :param source_expr: Source array expression, such as `g4dd_ref[4]`.
    :param coefficient: Integer prefactor after symbolic term collection.
    :param cos_power: Power of `cos_delta`.
    :param sin_power: Power of `sin_delta`.
    :return: C expression for one signed monomial times one source component.
    """
    factors: List[str] = []
    coefficient_abs = abs(coefficient)
    if coefficient_abs != 1:
        factors.append(f"{coefficient_abs}.0")
    factors.extend(["cos_delta"] * cos_power)
    factors.extend(["sin_delta"] * sin_power)
    factors.append(source_expr)
    term = " * ".join(factors)
    if coefficient < 0:
        return f"-{term}"
    return term


def _emit_wrapped_assignment(lhs: str, terms: List[str]) -> str:
    """
    Emit one wrapped C assignment from signed monomial terms.

    :param lhs: Left-hand-side C lvalue.
    :param terms: Signed monomial terms in emitted order.
    :return: One wrapped C assignment string.
    :raises ValueError: If `terms` is empty.
    """
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


def _build_metric_rotation_terms(mu: int, nu: int) -> List[str]:
    """
    Build signed monomial terms for one rotated serialized metric component.

    :param mu: First destination metric index.
    :param nu: Second destination metric index.
    :return: Signed monomial terms for the rotated serialized metric component.
    """
    term_map: Dict[Tuple[int, int, int], int] = {}
    for source_mu, sign_mu, cos_mu, sin_mu in _rotation_source_terms(mu):
        for source_nu, sign_nu, cos_nu, sin_nu in _rotation_source_terms(nu):
            metric_mu, metric_nu = source_mu, source_nu
            if metric_mu > metric_nu:
                metric_mu, metric_nu = metric_nu, metric_mu
            source_idx = _METRIC_COMPONENT_INDEX[(metric_mu, metric_nu)]
            term_key = (source_idx, cos_mu + cos_nu, sin_mu + sin_nu)
            term_map[term_key] = term_map.get(term_key, 0) + sign_mu * sign_nu
    terms: List[str] = []
    for source_idx, cos_power, sin_power in sorted(term_map):
        coefficient = term_map[(source_idx, cos_power, sin_power)]
        if coefficient == 0:
            continue
        terms.append(
            _format_scaled_source_term(
                f"g4dd_ref[{source_idx}]",
                coefficient,
                cos_power,
                sin_power,
            )
        )
    return terms


def _build_gamma_rotation_terms(alpha: int, mu: int, nu: int) -> List[str]:
    """
    Build signed monomial terms for one rotated serialized Christoffel component.

    :param alpha: Destination upper Christoffel index.
    :param mu: First destination lower Christoffel index.
    :param nu: Second destination lower Christoffel index.
    :return: Signed monomial terms for the rotated serialized Christoffel
        component.
    """
    term_map: Dict[Tuple[int, int, int], int] = {}
    for source_alpha, sign_alpha, cos_alpha, sin_alpha in _rotation_source_terms(alpha):
        for source_mu, sign_mu, cos_mu, sin_mu in _rotation_source_terms(mu):
            for source_nu, sign_nu, cos_nu, sin_nu in _rotation_source_terms(nu):
                gamma_mu, gamma_nu = source_mu, source_nu
                if gamma_mu > gamma_nu:
                    gamma_mu, gamma_nu = gamma_nu, gamma_mu
                source_idx = _GAMMA_COMPONENT_INDEX[(source_alpha, gamma_mu, gamma_nu)]
                term_key = (
                    source_idx,
                    cos_alpha + cos_mu + cos_nu,
                    sin_alpha + sin_mu + sin_nu,
                )
                term_map[term_key] = (
                    term_map.get(term_key, 0) + sign_alpha * sign_mu * sign_nu
                )
    terms: List[str] = []
    for source_idx, cos_power, sin_power in sorted(term_map):
        coefficient = term_map[(source_idx, cos_power, sin_power)]
        if coefficient == 0:
            continue
        terms.append(
            _format_scaled_source_term(
                f"gamma_ref[{source_idx}]",
                coefficient,
                cos_power,
                sin_power,
            )
        )
    return terms


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
