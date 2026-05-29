"""
Register a C function for dataset-specific raytracing spatial interpolation.

This module emits a reader-side spatial interpolation helper for the current
`two_blackholes_collide` workflow. The generated C API evaluates one photon
position `(x, y, z)` against a caller-specified list of slice indices and flat
output buffers for the interpolated Cartesian-basis `g4DD` and `Gamma4UDD`
components.

The target photon position is still converted from Cartesian to native
spherical coordinates with the generated BHaH inverse map. In contrast, raw
stencil nodes are now remapped with explicit dataset-specific spherical
reflection rules. This keeps the interpolation nodes themselves on the raw
uniform `(r, theta)` stencil required by
`interpolation_lagrange_uniform.h`, while remapping only the payload-read
indices back into the stored interior.

This module does not implement JSON parsing for the combined metadata. The
caller must parse the combined file elsewhere and populate the lightweight
runtime context consumed by the interpolation routine.

Author: Dalton Moone
        daltonmoone **at** gmail **dot** com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.BHaH.xx_tofrom_Cart import (
    register_CFunction__Cart_to_xx_and_nearest_i0i1i2,
)

RAYTRACING_SPATIAL_INTERP_DEFINES = r"""
// Dataset-specific constants for raytracing spatial interpolation.
#define RAYTRACING_SPATIAL_INTERP_PI 3.14159265358979323846264338327950288
#define RAYTRACING_RT_G4_COMPONENT_COUNT 10
#define RAYTRACING_RT_GAMMA_COMPONENT_COUNT 40
#define RAYTRACING_RT_TENSOR_COMPONENT_COUNT \
  (RAYTRACING_RT_G4_COMPONENT_COUNT + RAYTRACING_RT_GAMMA_COMPONENT_COUNT)
#define RAYTRACING_RT_RECORD_COMPONENT_COUNT 53
#define RAYTRACING_RT_POINT_RECORD_BYTES \
  (RAYTRACING_RT_RECORD_COMPONENT_COUNT * sizeof(double))
// Status codes returned by the raytracing spatial interpolation helper.
typedef enum {
  RAYTRACING_SPATIAL_INTERP_SUCCESS = 0,
  RAYTRACING_SPATIAL_INTERP_INVALID_METADATA = 1,
  RAYTRACING_SPATIAL_INTERP_INVALID_TARGET = 2,
  RAYTRACING_SPATIAL_INTERP_UNSUPPORTED_STENCIL = 3,
  RAYTRACING_SPATIAL_INTERP_IO_ERROR = 4,
  RAYTRACING_SPATIAL_INTERP_INVALID_SLICE_INDEX = 5
} raytracing_spatial_interp_status; // END ENUM: raytracing_spatial_interp_status

// Reusable parsed context for many photon spatial interpolation calls.
typedef struct {
  const uint64_t *restrict payload_offsets;
  uint64_t num_time_slices;
  double stored_phi_samples[2];
} raytracing_spatial_context_struct; // END STRUCT: raytracing_spatial_context_struct
"""
Bdefines_h.register_BHaH_defines(
    "raytracing_spatial_interpolation", RAYTRACING_SPATIAL_INTERP_DEFINES
)


def register_CFunction_raytracing_spatial_interpolation(
    CoordSystem: str,
    enable_simd: bool,
    project_dir: str,
) -> None:
    """
    Register the dataset-specific raytracing spatial interpolation helper.

    :param CoordSystem: Coordinate system used by the source dataset.
    :param enable_simd: Whether SIMD helper headers are already available.
    :param project_dir: Destination project directory for copied headers.
    :return: None.
    :raises ValueError: If `CoordSystem` is not the dataset-specific value.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import os
    >>> import tempfile
    >>> import nrpy.c_function as cfc
    >>> cfc.CFunction_dict.clear()
    >>> with tempfile.TemporaryDirectory() as project_dir:
    ...     old_cache_home = os.environ.get("XDG_CACHE_HOME")
    ...     _ = os.environ.__setitem__("XDG_CACHE_HOME", project_dir)
    ...     with contextlib.redirect_stdout(io.StringIO()):
    ...         _ = register_CFunction_raytracing_spatial_interpolation(
    ...             "Spherical", enable_simd=True, project_dir=project_dir
    ...         )
    ...     if old_cache_home is None:
    ...         _ = os.environ.pop("XDG_CACHE_HOME", None)
    ...     else:
    ...         _ = os.environ.__setitem__("XDG_CACHE_HOME", old_cache_home)
    ...     generated = cfc.CFunction_dict[
    ...         "raytracing_spatial_interpolation__rfm__Spherical"
    ...     ].full_function
    ...     (
    ...         "phi_toggle ^= 1;" in generated
    ...         and "xx_to_Cart(params" not in generated
    ...         and "target_phi - phi_ref" in generated
    ...     )
    True
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return

    # Step 1: Validate high-level code-generation assumptions.
    if CoordSystem != "Spherical":
        raise ValueError(
            "raytracing_spatial_interpolation is currently implemented only "
            f"for CoordSystem='Spherical'; found '{CoordSystem}'."
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

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "sys/types.h",
        "interpolation_lagrange_uniform.h",
    ]

    prefunc = r"""
/**
 * Map one raw extended-grid stencil node back into the stored payload.
 *
 * The interpolation coefficients are always built on the raw uniform stencil.
 * This helper only recovers the in-domain payload index for one raw node by
 * applying the dataset-specific spherical reflection rules directly.
 *
 * @param[in] params Generated BHaH parameter struct.
 * @param[in] r_ext Raw stencil radial coordinate, possibly outside the interior.
 * @param[in] theta_ext Raw stencil polar coordinate, possibly outside the interior.
 * @param[in] i2_base Logical phi index for the selected stored reference plane.
 * @param[out] i0i1i2_map In-domain logical payload index recovered by reflection.
 * @return Status code indicating whether the mapped node lands inside the payload.
 */
static int raytracing_spatial_map_extended_node(
    const params_struct *restrict params,
    const REAL r_ext,
    const REAL theta_ext,
    const int i2_base,
    int i0i1i2_map[3]) {
  const int payload_i0_start = NGHOSTS;
  const int payload_i0_end = NGHOSTS + params->Nxx0;
  const int payload_i1_start = NGHOSTS;
  const int payload_i1_end = NGHOSTS + params->Nxx1;
  const int payload_i2_start = NGHOSTS;
  const int payload_i2_end = NGHOSTS + params->Nxx2;
  const REAL theta_tol = (REAL)1.0e-12;
  const REAL r_endpoint_tol =
      (REAL)(1.0e-12 * fmax(1.0, fabs((double)params->dxx0)));
  const REAL theta_endpoint_tol =
      (REAL)(1.0e-12 * fmax(1.0, fabs((double)params->dxx1)));
  const REAL theta_max_supported =
      (REAL)(params->xxmin1 +
             (((payload_i1_end - 1 - NGHOSTS) + 0.5) * params->dxx1));
  const REAL r_max_supported =
      (REAL)(params->xxmin0 +
             (((payload_i0_end - 1 - NGHOSTS) + 0.5) * params->dxx0));
  REAL r_map = r_ext;
  REAL theta_map = theta_ext;
  REAL theta_index_map = theta_ext;
  int phi_toggle = 0;
  int i0_map;
  int i1_map;
  int i2_map;

  // Step 1: Reflect negative-radius raw nodes through the origin.
  if (r_map < (REAL)0.0) {
    r_map = -r_map;
    theta_map = (REAL)RAYTRACING_SPATIAL_INTERP_PI - theta_map;
    phi_toggle ^= 1;
  } // END IF: raw stencil node crossed the inner radial boundary

  // Step 2: Fold theta back into the physical [0, pi] range.
  if (theta_map < (REAL)0.0) {
    theta_map = -theta_map;
    phi_toggle ^= 1;
  } else if (theta_map > (REAL)RAYTRACING_SPATIAL_INTERP_PI) {
    theta_map = (REAL)(2.0 * RAYTRACING_SPATIAL_INTERP_PI) - theta_map;
    phi_toggle ^= 1;
  } // END IF: raw stencil node crossed a polar boundary

  // Step 3: Only accept theta values that now land inside the physical range.
  if (theta_map < -theta_tol ||
      theta_map > (REAL)RAYTRACING_SPATIAL_INTERP_PI + theta_tol) {
    return RAYTRACING_SPATIAL_INTERP_UNSUPPORTED_STENCIL;
  } // END IF: reflected theta remained materially out of range
  if (theta_map < (REAL)0.0) {
    theta_map = (REAL)0.0;
  } else if (theta_map > (REAL)RAYTRACING_SPATIAL_INTERP_PI) {
    theta_map = (REAL)RAYTRACING_SPATIAL_INTERP_PI;
  } // END IF: theta exceeded the physical range only by roundoff

  // Step 4: Repair only small endpoint roundoff before converting to logical indices.
  if (r_map > r_max_supported && fabs((double)(r_map - r_max_supported)) <= r_endpoint_tol) {
    r_map = r_max_supported;
  } // END IF: radial coordinate landed on the last stored cell center to roundoff
  theta_index_map = theta_map;
  if (theta_index_map > theta_max_supported &&
      theta_index_map <=
          theta_max_supported + (REAL)(0.5 * params->dxx1) + theta_endpoint_tol) {
    theta_index_map = theta_max_supported;
  } // END IF: theta landed in the last physical half-cell and should read the endpoint payload cell

  // Step 5: Convert the reflected coordinates and phi toggle to payload indices.
  i0_map = (int)((r_map - params->xxmin0) / params->dxx0 + (REAL)NGHOSTS);
  i1_map = (int)((theta_index_map - params->xxmin1) / params->dxx1 + (REAL)NGHOSTS);
  if (phi_toggle == 0) {
    i2_map = i2_base;
  } else {
    i2_map =
        payload_i2_start + payload_i2_end - 1 - (i2_base - payload_i2_start);
  } // END IF: apply the pi shift through the stored two-plane layout

  // Step 6: Confirm that the remapped logical index lies inside the stored payload.
  if (i0_map < payload_i0_start || i0_map >= payload_i0_end ||
      i1_map < payload_i1_start || i1_map >= payload_i1_end ||
      i2_map < payload_i2_start || i2_map >= payload_i2_end) {
    return RAYTRACING_SPATIAL_INTERP_UNSUPPORTED_STENCIL;
  } // END IF: reflected stencil node could not be supported by the payload
  i0i1i2_map[0] = i0_map;
  i0i1i2_map[1] = i1_map;
  i0i1i2_map[2] = i2_map;
  return RAYTRACING_SPATIAL_INTERP_SUCCESS;
} // END FUNCTION: raytracing_spatial_map_extended_node

/**
 * Rotate Cartesian-basis metric and Christoffel components about the z axis.
 *
 * The payload stores tensors on one reference azimuthal plane. Axisymmetry
 * recovers the target azimuth by an active spatial rotation through
 * `delta_phi`, embedded in 4D so that:
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
static void raytracing_spatial_rotate_about_z(
    const REAL delta_phi,
    const REAL g4dd_ref[RAYTRACING_RT_G4_COMPONENT_COUNT],
    const REAL gamma_ref[RAYTRACING_RT_GAMMA_COMPONENT_COUNT],
    REAL g4dd_rot[RAYTRACING_RT_G4_COMPONENT_COUNT],
    REAL gamma_rot[RAYTRACING_RT_GAMMA_COMPONENT_COUNT]) {
  REAL g_ref[4][4] = {{0}};
  REAL g_dst[4][4] = {{0}};
  REAL Gamma_ref[4][4][4] = {{{0}}};
  REAL Gamma_dst[4][4][4] = {{{0}}};
  REAL Lambda[4][4] = {{0}};
  REAL LambdaT[4][4] = {{0}};
  const REAL cosD = cos(delta_phi);
  const REAL sinD = sin(delta_phi);
  int gamma_index = 0;

  // Step 1: Unpack the stored upper-triangular metric components.
  g_ref[0][0] = g4dd_ref[0];
  g_ref[0][1] = g_ref[1][0] = g4dd_ref[1];
  g_ref[0][2] = g_ref[2][0] = g4dd_ref[2];
  g_ref[0][3] = g_ref[3][0] = g4dd_ref[3];
  g_ref[1][1] = g4dd_ref[4];
  g_ref[1][2] = g_ref[2][1] = g4dd_ref[5];
  g_ref[1][3] = g_ref[3][1] = g4dd_ref[6];
  g_ref[2][2] = g4dd_ref[7];
  g_ref[2][3] = g_ref[3][2] = g4dd_ref[8];
  g_ref[3][3] = g4dd_ref[9];

  // Step 2: Unpack the serialized Christoffels, mirroring lower-index symmetry.
  for (int alpha = 0; alpha < 4; alpha++) {
    for (int mu = 0; mu < 4; mu++) {
      for (int nu = mu; nu < 4; nu++) {
        Gamma_ref[alpha][mu][nu] = gamma_ref[gamma_index];
        Gamma_ref[alpha][nu][mu] = gamma_ref[gamma_index];
        gamma_index++;
      }
    }
  }

  // Step 3: Build the active z-axis rotation matrix in Cartesian components.
  Lambda[0][0] = 1.0;
  Lambda[1][1] = cosD;
  Lambda[1][2] = -sinD;
  Lambda[2][1] = sinD;
  Lambda[2][2] = cosD;
  Lambda[3][3] = 1.0;

  // Step 4: For pure rotations, the inverse acting on lower indices is Lambda^T.
  for (int a = 0; a < 4; a++) {
    for (int b = 0; b < 4; b++) {
      LambdaT[a][b] = Lambda[b][a];
    }
  }

  // Step 5: Rotate the covariant metric with the lower-index transform.
  for (int mu = 0; mu < 4; mu++) {
    for (int nu = 0; nu < 4; nu++) {
      REAL sum = 0.0;
      for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 4; b++) {
          sum += LambdaT[a][mu] * LambdaT[b][nu] * g_ref[a][b];
        }
      }
      g_dst[mu][nu] = sum;
    }
  }

  // Step 6: Rotate the connection with one upper and two lower tensor factors.
  for (int alpha = 0; alpha < 4; alpha++) {
    for (int mu = 0; mu < 4; mu++) {
      for (int nu = 0; nu < 4; nu++) {
        REAL sum = 0.0;
        for (int b = 0; b < 4; b++) {
          for (int g = 0; g < 4; g++) {
            for (int d = 0; d < 4; d++) {
              sum += Lambda[alpha][b] * LambdaT[g][mu] * LambdaT[d][nu] *
                     Gamma_ref[b][g][d];
            }
          }
        }
        Gamma_dst[alpha][mu][nu] = sum;
      }
    }
  }

  // Step 7: Repack into the writer's serialized ordering.
  g4dd_rot[0] = g_dst[0][0];
  g4dd_rot[1] = g_dst[0][1];
  g4dd_rot[2] = g_dst[0][2];
  g4dd_rot[3] = g_dst[0][3];
  g4dd_rot[4] = g_dst[1][1];
  g4dd_rot[5] = g_dst[1][2];
  g4dd_rot[6] = g_dst[1][3];
  g4dd_rot[7] = g_dst[2][2];
  g4dd_rot[8] = g_dst[2][3];
  g4dd_rot[9] = g_dst[3][3];

  gamma_index = 0;
  for (int alpha = 0; alpha < 4; alpha++) {
    for (int mu = 0; mu < 4; mu++) {
      for (int nu = mu; nu < 4; nu++) {
        gamma_rot[gamma_index] = Gamma_dst[alpha][mu][nu];
        gamma_index++;
      }
    }
  }
} // END FUNCTION: raytracing_spatial_rotate_about_z
"""

    cfunc_type = "int"
    name = "raytracing_spatial_interpolation"
    params = """FILE *restrict fp,
                const raytracing_spatial_context_struct *restrict context,
                const params_struct *restrict params,
                const REAL x, const REAL y, const REAL z,
                const int n_interp_ghosts,
                const int num_target_slices,
                const int *restrict slice_indices,
                REAL *restrict g4dd_out,
                REAL *restrict gamma4udd_out"""
    body = r"""
  // Step 1: Bind the trusted parsed context.
  const REAL xCart[3] = {x, y, z};
  REAL xx_target[3];
  int center_idx[3];

  if (context == NULL || fp == NULL || params == NULL || slice_indices == NULL ||
      g4dd_out == NULL || gamma4udd_out == NULL || context->payload_offsets == NULL) {
    return RAYTRACING_SPATIAL_INTERP_INVALID_METADATA;
  }
  if (num_target_slices <= 0) {
    return RAYTRACING_SPATIAL_INTERP_INVALID_METADATA;
  }
  if (n_interp_ghosts < 1) {
    return RAYTRACING_SPATIAL_INTERP_UNSUPPORTED_STENCIL;
  }
  const uint64_t *restrict payload_offsets = context->payload_offsets;
  if (!isfinite((double)x) || !isfinite((double)y) || !isfinite((double)z)) {
    return RAYTRACING_SPATIAL_INTERP_INVALID_TARGET;
  }

  // Step 2: Convert the target Cartesian point to native spherical coordinates.
  Cart_to_xx_and_nearest_i0i1i2(params, xCart, xx_target, center_idx);

  const REAL target_r = xx_target[0];
  const REAL target_theta = xx_target[1];
  REAL target_phi = xx_target[2];
  const REAL rho_sq = x * x + y * y;
  const REAL origin_epsilon = 1.0e-14;
  const REAL axis_rho_epsilon = 1.0e-14;
  const int interp_order = 2 * n_interp_ghosts + 1;
  const REAL r_max_supported =
      (REAL)(params->xxmin0 + (((params->Nxx0 - 1) + 0.5) * params->dxx0));
  const REAL r_support_tol =
      (REAL)(1.0e-12 * fmax(1.0, fabs((double)r_max_supported)));
  int phi_plane;
  int i2_base;
  REAL phi_ref = 0.0;

  if (!isfinite((double)target_r) || !isfinite((double)target_theta) ||
      !isfinite((double)target_phi)) {
    return RAYTRACING_SPATIAL_INTERP_INVALID_TARGET;
  }
  if (target_phi >= (REAL)RAYTRACING_SPATIAL_INTERP_PI - (REAL)1.0e-12 &&
      target_phi <= (REAL)RAYTRACING_SPATIAL_INTERP_PI + (REAL)1.0e-12) {
    target_phi -= (REAL)(2.0 * RAYTRACING_SPATIAL_INTERP_PI);
  }
  if (target_phi < (REAL)(-RAYTRACING_SPATIAL_INTERP_PI - 1.0e-12) ||
      target_phi > (REAL)(RAYTRACING_SPATIAL_INTERP_PI + 1.0e-12)) {
    return RAYTRACING_SPATIAL_INTERP_INVALID_TARGET;
  }

  if (target_r <= origin_epsilon || rho_sq <= axis_rho_epsilon * axis_rho_epsilon) {
    return RAYTRACING_SPATIAL_INTERP_INVALID_TARGET;
  }
  if (interp_order > params->Nxx0 || interp_order > params->Nxx1) {
    return RAYTRACING_SPATIAL_INTERP_UNSUPPORTED_STENCIL;
  }

  // Step 3: Select the stored phi plane and build 2D interpolation coefficients.
  phi_plane = target_phi < (REAL)0.0 ? 0 : 1;
  i2_base = NGHOSTS + phi_plane;
  phi_ref = (REAL)context->stored_phi_samples[phi_plane];

  REAL inv_denom[interp_order];
  REAL src_r_stencil[interp_order];
  REAL src_theta_stencil[interp_order];
  REAL diffs_r[interp_order];
  REAL diffs_theta[interp_order];
  REAL coeff_r[interp_order];
  REAL coeff_theta[interp_order];
  uint64_t mapped_point_index[interp_order][interp_order];
  const REAL normalization_2d =
      pow((REAL)(params->dxx0 * params->dxx1), -(interp_order - 1));

  for (int u = 0; u < interp_order; u++) {
    const int i0_raw = center_idx[0] + (u - n_interp_ghosts);
    // Keep the polynomial nodes on the raw extended grid, even if the later
    // payload reads are recovered back into the interior.
    src_r_stencil[u] =
        (REAL)(params->xxmin0 + (((i0_raw - NGHOSTS) + 0.5) * params->dxx0));
  } // END LOOP: for u over radial stencil nodes
  for (int v = 0; v < interp_order; v++) {
    const int i1_raw = center_idx[1] + (v - n_interp_ghosts);
    src_theta_stencil[v] =
        (REAL)(params->xxmin1 + (((i1_raw - NGHOSTS) + 0.5) * params->dxx1));
  } // END LOOP: for v over theta stencil nodes

  if (src_r_stencil[interp_order - 1] > r_max_supported + r_support_tol) {
    return RAYTRACING_SPATIAL_INTERP_UNSUPPORTED_STENCIL;
  } // END IF: largest raw radial stencil node exceeded the supported payload radius

  compute_inv_denom(interp_order, inv_denom);
  compute_diffs_xi(interp_order, target_r, src_r_stencil, diffs_r);
  compute_diffs_xi(interp_order, target_theta, src_theta_stencil, diffs_theta);
  compute_lagrange_basis_coeffs_xi(interp_order, inv_denom, diffs_r, coeff_r);
  compute_lagrange_basis_coeffs_xi(
      interp_order, inv_denom, diffs_theta, coeff_theta);

  // Step 4: Map each raw stencil node once for this photon target.
  for (int v = 0; v < interp_order; v++) {
    const REAL theta_ext = src_theta_stencil[v];

    for (int u = 0; u < interp_order; u++) {
      const REAL r_ext = src_r_stencil[u];
      int i0i1i2_map[3] = {-1, -1, -1};

      const int map_status = raytracing_spatial_map_extended_node(
          params, r_ext, theta_ext, i2_base, i0i1i2_map);
      if (map_status != RAYTRACING_SPATIAL_INTERP_SUCCESS) {
        return map_status;
      }
      const uint64_t j0 = (uint64_t)(i0i1i2_map[0] - NGHOSTS);
      const uint64_t j1 = (uint64_t)(i0i1i2_map[1] - NGHOSTS);
      const uint64_t j2 = (uint64_t)(i0i1i2_map[2] - NGHOSTS);
      mapped_point_index[v][u] =
          j0 + (uint64_t)params->Nxx0 * (j1 + (uint64_t)params->Nxx1 * j2);
    } // END LOOP: for u over radial stencil nodes during remap precompute
  } // END LOOP: for v over theta stencil nodes during remap precompute

  // Step 5: Interpolate each requested time slice independently.
  for (int which_slice = 0; which_slice < num_target_slices; which_slice++) {
    const int slice_index = slice_indices[which_slice];
    if (slice_index < 0 || (uint64_t)slice_index >= context->num_time_slices) {
      return RAYTRACING_SPATIAL_INTERP_INVALID_SLICE_INDEX;
    }
    const uint64_t payload_offset = payload_offsets[slice_index];
    REAL tensor_ref[RAYTRACING_RT_TENSOR_COMPONENT_COUNT] = {0};

    // Step 5.a: Read and accumulate the weighted stencil nodes for this slice.
    for (int v = 0; v < interp_order; v++) {
      for (int u = 0; u < interp_order; u++) {
        const uint64_t record_offset =
            payload_offset +
            mapped_point_index[v][u] * (uint64_t)RAYTRACING_RT_POINT_RECORD_BYTES;
        const uint64_t tensor_offset =
            record_offset + 3ULL * (uint64_t)sizeof(double);
        const REAL weight = normalization_2d * coeff_theta[v] * coeff_r[u];
        double tensor_record[RAYTRACING_RT_TENSOR_COMPONENT_COUNT];

        if (fseeko(fp, (off_t)tensor_offset, SEEK_SET) != 0) {
          return RAYTRACING_SPATIAL_INTERP_IO_ERROR;
        }
        if (fread(tensor_record, sizeof(double), RAYTRACING_RT_TENSOR_COMPONENT_COUNT, fp) !=
            RAYTRACING_RT_TENSOR_COMPONENT_COUNT) {
          return RAYTRACING_SPATIAL_INTERP_IO_ERROR;
        }
        for (int comp = 0; comp < RAYTRACING_RT_TENSOR_COMPONENT_COUNT; comp++) {
          tensor_ref[comp] += weight * (REAL)tensor_record[comp];
        } // END LOOP: for comp over metric and Christoffel payload components
      } // END LOOP: for u over radial stencil nodes
    } // END LOOP: for v over theta stencil nodes

    // Step 5.b: Rotate the interpolated Cartesian-basis tensors to the target azimuth.
    raytracing_spatial_rotate_about_z(
        target_phi - phi_ref,
        &tensor_ref[0],
        &tensor_ref[RAYTRACING_RT_G4_COMPONENT_COUNT],
        &g4dd_out[which_slice * RAYTRACING_RT_G4_COMPONENT_COUNT],
        &gamma4udd_out[which_slice * RAYTRACING_RT_GAMMA_COMPONENT_COUNT]);
  } // END LOOP: for which_slice over requested slice indices

  return RAYTRACING_SPATIAL_INTERP_SUCCESS;
"""

    desc = r"""Interpolate Cartesian raytracing tensors at one photon position.

The caller supplies a trusted spatial context and requested time-slice indices.
The helper builds one native `(r, theta)` stencil, remaps the payload-read
indices once with explicit spherical reflections, reads tensor components for
each requested slice, interpolates them, rotates to the target azimuth, and
writes flat per-slice outputs.

@param[in,out] fp Open combined-file stream.
@param[in] context Trusted spatial context.
@param[in] params Generated BHaH grid parameters.
@param x Cartesian x coordinate.
@param y Cartesian y coordinate.
@param z Cartesian z coordinate.
@param n_interp_ghosts Interpolation half-width.
@param num_target_slices Number of requested time slices.
@param[in] slice_indices Requested slice indices.
@param[out] g4dd_out Flat metric output.
@param[out] gamma4udd_out Flat Christoffel output.
@return Status code indicating success or the interpolation failure reason.

@note Callers must serialize access if sharing one `FILE *`.
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


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
