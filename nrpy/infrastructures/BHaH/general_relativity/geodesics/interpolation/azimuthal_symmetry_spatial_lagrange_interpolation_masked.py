"""
Generate spatial Lagrange interpolation with excluded native grid nodes.

Retained rho nodes are interpolated independently in each z row, followed by
interpolation over the usable z rows. The ordinary spatial generator supplies
the same Cartesian tensor rotations and derivative transformations for both
variants. Node counts and weights are calculated at runtime.

Author: Dalton J. Moone
"""

from typing import Union

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg

MASKED_SPATIAL_WEIGHTS_C_CODE = r"""
/**
 * Calculate nonuniform Lagrange values and derivatives, including exact nodes.
 *
 * @param count Number of distinct retained nodes, at least two.
 * @param target Native interpolation coordinate.
 * @param[in] nodes Retained native coordinates in ascending order.
 * @param[out] values Lagrange basis values.
 * @param[out] derivatives Native-coordinate basis derivatives.
 */
static void spatial_lagrange_masked_basis(
    const int count, const REAL target, const REAL nodes[count],
    REAL values[count], REAL derivatives[count]) {
  const REAL scale = nodes[count - 1] - nodes[0];
  REAL weights[count];
  int exact_node = -1;
  for (int i = 0; i < count; ++i) {
    weights[i] = 1.0;
    for (int j = 0; j < count; ++j)
      if (i != j)
        weights[i] /= (nodes[i] - nodes[j]) / scale;
    values[i] = 0.0;
    derivatives[i] = 0.0;
    if (target == nodes[i])
      exact_node = i;
  } // END LOOP: form barycentric weights
  if (exact_node >= 0) {
    values[exact_node] = 1.0;
    for (int i = 0; i < count; ++i) {
      if (i == exact_node)
        continue;
      derivatives[i] = weights[i] /
          (weights[exact_node] * (nodes[exact_node] - nodes[i]));
      derivatives[exact_node] -= derivatives[i];
    } // END LOOP: differentiate at exact node
    return;
  } // END IF: target equals retained node
  REAL denominator = 0.0;
  for (int i = 0; i < count; ++i) {
    const REAL delta = (target - nodes[i]) / scale;
    values[i] = weights[i] / delta;
    denominator += values[i];
  } // END LOOP: accumulate barycentric denominator
  for (int i = 0; i < count; ++i)
    values[i] /= denominator;
  // Pairwise differences avoid subtracting two nearly equal reciprocal
  // distances when the target is very close to a retained node.
  for (int i = 0; i < count; ++i) {
    REAL slope = 0.0;
    for (int j = 0; j < count; ++j)
      if (i != j)
        slope += values[j] * (nodes[j] - nodes[i]) /
            ((target - nodes[j]) * (target - nodes[i]));
    derivatives[i] = values[i] * slope;
  } // END LOOP: normalize basis and derivatives
} // END FUNCTION: nonuniform Lagrange basis

/**
 * Form rho-then-z weights using the retained nodes of each spatial row.
 *
 * @param count Original number of nodes per spatial direction.
 * @param target0 Native rho interpolation coordinate.
 * @param target2 Native z interpolation coordinate.
 * @param[in] nodes0 Original native rho nodes.
 * @param[in] nodes2 Original native z nodes.
 * @param[in] mask Row-major exclusion bytes, one for each original node.
 * @param[out] value Two-dimensional interpolation weights.
 * @param[out] derivative0 Native rho derivative weights.
 * @param[out] derivative2 Native z derivative weights.
 * @return SUCCESS or UNSUPPORTED_STENCIL when fewer than two rows remain.
 */
static int spatial_lagrange_masked_weights(
    const int count, const REAL target0, const REAL target2,
    const REAL nodes0[count], const REAL nodes2[count],
    const unsigned char *restrict mask,
    REAL value[count][count], REAL derivative0[count][count],
    REAL derivative2[count][count]) {
  int usable_count = 0;
  int usable_rows[count];
  REAL usable_z[count];
  for (int v = 0; v < count; ++v) {
    int retained_count = 0;
    int retained_indices[count];
    REAL retained_nodes[count];
    for (int u = 0; u < count; ++u) {
      value[v][u] = derivative0[v][u] = derivative2[v][u] = 0.0;
      if (mask[v * count + u] == 0) {
        retained_indices[retained_count] = u;
        retained_nodes[retained_count++] = nodes0[u];
      } // END IF: retain rho node
    } // END LOOP: gather retained rho nodes
    if (retained_count < 2)
      continue;
    REAL row_values[retained_count], row_derivatives[retained_count];
    spatial_lagrange_masked_basis(
        retained_count, target0, retained_nodes, row_values, row_derivatives);
    for (int k = 0; k < retained_count; ++k) {
      const int u = retained_indices[k];
      value[v][u] = row_values[k];
      derivative0[v][u] = row_derivatives[k];
    } // END LOOP: scatter row basis weights
    usable_rows[usable_count] = v;
    usable_z[usable_count++] = nodes2[v];
  } // END LOOP: interpolate retained rho rows
  if (usable_count < 2)
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;
  REAL z_values[usable_count], z_derivatives[usable_count];
  spatial_lagrange_masked_basis(
      usable_count, target2, usable_z, z_values, z_derivatives);
  for (int j = 0; j < usable_count; ++j) {
    const int v = usable_rows[j];
    for (int u = 0; u < count; ++u) {
      derivative2[v][u] = value[v][u] * z_derivatives[j];
      derivative0[v][u] *= z_values[j];
      value[v][u] *= z_values[j];
    } // END LOOP: combine rho and z weights
  } // END LOOP: combine usable z rows
  return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS;
} // END FUNCTION: masked spatial polynomial weights
"""


def register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation_masked(
    CoordSystem: str,
    interpolation_method: str = "g4DD",
    enable_simd: bool = False,
    project_dir: str = ".",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the masked variant of the Cartesian spatial interpolation function.

    :param CoordSystem: Numerical dataset coordinate system.
    :param interpolation_method: Stored metric or connection components.
    :param enable_simd: Whether SIMD helper headers are already available.
    :param project_dir: Generated project directory for copied headers.
    :return: None during registration, otherwise the updated NRPy environment.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import os
    >>> import tempfile
    >>> from unittest.mock import patch
    >>> with tempfile.TemporaryDirectory() as directory, patch.dict(os.environ, {"XDG_CACHE_HOME": directory}):
    ...     cfc.CFunction_dict.clear()
    ...     with contextlib.redirect_stdout(io.StringIO()):
    ...         _ = register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation_masked(
    ...             "SinhCylindricalv2n2", enable_simd=True, project_dir=directory
    ...         )
    ...     masked = cfc.CFunction_dict[
    ...         "azimuthal_symmetry_spatial_lagrange_interpolation_masked__rfm__SinhCylindricalv2n2"
    ...     ]
    >>> masked.params.endswith("const unsigned char *restrict spatial_mask")
    True
    """
    from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.azimuthal_symmetry_spatial_lagrange_interpolation import (  # pylint: disable=import-outside-toplevel
        register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation,
    )

    if (
        f"azimuthal_symmetry_spatial_lagrange_interpolation__rfm__{CoordSystem}"
        not in cfc.CFunction_dict
    ):
        register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation(
            CoordSystem,
            interpolation_method=interpolation_method,
            enable_simd=enable_simd,
            project_dir=project_dir,
        )
    return register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation(
        CoordSystem,
        interpolation_method=interpolation_method,
        enable_simd=enable_simd,
        project_dir=project_dir,
        masked_weight_code=MASKED_SPATIAL_WEIGHTS_C_CODE,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
