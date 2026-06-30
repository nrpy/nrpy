"""
Register a C function to perform numerical integrations over the apparent horizon surface for spin and mass diagnostics.
Needs: Integrands built in equations/general_relativity/bhahaha/SpECTRESpinEstimate.py

Derivation and AKV convention summary:

- The scalar potentials z_alpha are the AKV representatives. With the chosen
  surface orientation eps^{theta phi} = +1/sqrt(q), a potential zeta generates
  the tangential vector field phi^A = eps^{AB} D_B zeta. Reversing the
  orientation reverses eps^{AB}, Omega, and the reported spin components.
- The spin function used by the integrands is
  Omega = eps^{AB} D_A X_B, where X_B = e_B^i K_ij s^j. On a closed surface,
  integration by parts gives the AKV angular momentum
  J[zeta] = (1/(8*pi)) int phi^A X_A dA
          = (1/(8*pi)) int zeta Omega dA
  for this sign convention.
- The runtime potential solve discretizes the symmetric scalar AKV weak form
  K(eta,z) = int (Delta eta)(Delta z) dA
           - int R D_A eta D^A z dA,
  M(eta,z) = int D_A eta D^A z dA, with K z = lambda M z and constants
  removed by an area-weighted mean-zero constraint.
- Normalization is part of the convention, not a post-processing detail:
  rescaling z_alpha rescales S_alpha. This file therefore normalizes each
  centered potential with int z_alpha^2 dA = A^3/(48*pi^2) before evaluating
  S_alpha = (1/(8*pi)) int z_alpha Omega dA.

Author: Ralston Graves
        ralstonkgraves **at** gmail **dot** com
"""

from typing import Union

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.equations.general_relativity.bhahaha.SpECTRESpinEstimate import (
    SpECTRESpinEstimate,
)
from nrpy.helpers.generic import clang_format
from nrpy.infrastructures.BHaH.CurviBoundaryConditions.apply_bcs_inner_only import (
    APPLY_PARITY_BRANCHLESS_PREFUNC,
)

_SPECTRE_SPIN_SCRATCH_GFS = (
    "SE_qDD00",
    "SE_qDD01",
    "SE_qDD11",
    "SE_XD0",
    "SE_XD1",
    "zU0",
    "zU1",
    "zU2",
)

_SPECTRE_SPIN_COORD_COVARIANT_GFS = (
    "SE_qDD00",
    "SE_qDD01",
    "SE_qDD11",
    "SE_XD0",
    "SE_XD1",
)


def register_CFunction_diagnostics_spectre_spin(
    CoordSystem: str = "Spherical",
    enable_rfm_precompute: bool = False,
    enable_fd_functions: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register a C function that computes the SpECTRE-style spin vector.

    :param CoordSystem: The coordinate system to use, defaults to "Spherical".
    :param enable_rfm_precompute: Whether to enable RFM precompute.
    :param enable_fd_functions: Whether to enable finite-difference functions.
    :return: An NRPyEnv_type object if registration is successful, otherwise None.
    :raises ValueError: If a precompute gridfunction has an unsupported rank.

    """
    if par.parval_from_str("fp_type") != "double":
        raise ValueError(
            "SpECTRE spin diagnostics require fp_type=double because the "
            "runtime eigensolver calls double-precision PRIMME (dprimme)."
        )

    if pcg.pcg_registration_phase():
        pcg.register_func_call(
            f"{__name__}.{register_CFunction_diagnostics_spectre_spin.__name__}",
            locals(),
        )
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "akv_primme.h"]
    desc = r"""
    Compute the SpECTRE-style dimensionless spin vector diagnostic and store it in the diagnostics struct.
    """
    cfunc_type = "int"
    cfunc_name = "diagnostics_spectre_spin"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    # Step 1: Get an instance of the symbolic calculator.
    spin_calc = SpECTRESpinEstimate[
        CoordSystem + ("_rfm_precompute" if enable_rfm_precompute else "")
    ]

    # Step 2: Temporarily register spin scratch symbols for FD codegen only.
    saved_spectre_spin_gfs = {
        gf_name: gri.glb_gridfcs_dict.get(gf_name)
        for gf_name in _SPECTRE_SPIN_SCRATCH_GFS
    }
    try:
        for gf_name in _SPECTRE_SPIN_SCRATCH_GFS:
            gri.glb_gridfcs_dict.pop(gf_name, None)

        _ = gri.register_gridfunctions_for_single_rank2(
            "SE_qDD",
            symmetry="sym01",
            dimension=2,
            group="AUX",
            gf_array_name="spectre_spin_gfs",
        )
        _ = gri.register_gridfunctions_for_single_rank1(
            "SE_XD",
            dimension=2,
            group="AUX",
            gf_array_name="spectre_spin_gfs",
        )
        _ = gri.register_gridfunctions_for_single_rank1(
            "zU",
            dimension=3,
            group="AUX",
            gf_array_name="spectre_spin_gfs",
        )

        gf_assignments = spin_calc.get_gridfunction_assignments(
            include_flux_density=False
        )
        gf_names = [str(sym) for sym in gf_assignments.keys()]

        lhss_precompute = [
            gri.BHaHGridFunction.access_gf(gf_name, gf_array_name="spectre_spin_gfs")
            for gf_name in gf_names
        ]
        rhss_precompute = list(gf_assignments.values())
        gf_macros = [f"{gf_name.upper()}GF" for gf_name in gf_names]

        parity_conditions_rank2 = {
            (0, 0): 4,
            (0, 1): 5,
            (0, 2): 6,
            (1, 1): 7,
            (1, 2): 8,
            (2, 2): 9,
        }
        parity_entries = []
        for gf_name, gf_macro in zip(gf_names, gf_macros):
            gf = gri.glb_gridfcs_dict[gf_name]
            if gf.name in _SPECTRE_SPIN_COORD_COVARIANT_GFS:
                # These surface coordinate-basis covariant fields are transformed
                # with bc->deriv_jacobian in the custom ghost-fill helper below.
                parity_value = 0
            elif gf.rank == 0:
                parity_value = 0
            elif gf.rank == 1:
                parity_value = int(gf.name[-1]) + 1
            elif gf.rank == 2:
                parity_value = parity_conditions_rank2[
                    (int(gf.name[-2]), int(gf.name[-1]))
                ]
            else:
                raise ValueError(
                    f"Unsupported spin-diagnostic precompute gridfunction rank: {gf.name}, rank={gf.rank}"
                )
            parity_entries.append(f"  [{gf_macro}] = {parity_value},")
        parity_table_entries = "\n".join(parity_entries)
        selected_precompute_gfs = ", ".join(gf_macros)
        scratch_gf_defines = "\n".join(
            f"#define {gf_name.upper()}GF {idx}"
            for idx, gf_name in enumerate(_SPECTRE_SPIN_SCRATCH_GFS)
        )

        # Step 3: Generate C code for intermediate surface fields.
        precompute_c_code = ccg.c_codegen(
            rhss_precompute,
            lhss_precompute,
            enable_fd_codegen=True,
            enable_fd_functions=enable_fd_functions,
        )

        # Step 4: Retrieve the dictionary of all per-point integrands.
        integrands_dict = spin_calc.get_public_integrands()

        # Step 5: Extract the symbolic expressions we need to integrate.
        # According to the SpECTRESpinEstimate documentation, we need to integrate:
        # - 1 (for the Area A)
        # - x^i (for the centroid XU)
        # - R (for R0)
        # - x^i * R (for XRU)
        # - Omega (for O0)
        # - x^i * Omega (for XOU)
        # - z_alpha * Omega (for ZOU)
        # - |Omega| (for Oabs, used in the near-zero policy)

        # Note: The 'integrand' is the quantity 'f' in ∮ f dA.
        # The differential area element is dA = sqrt(q) * weights * dθ * dφ.
        # We will pass sqrt(q) to c_codegen as 'area_density' and handle the
        # weights and coordinate steps in the C code loop.

        area_density = integrands_dict["area_density"]

        # List of all symbolic quantities to be evaluated inside the loop
        integrand_c_vars = [
            "const REAL area_density",
            "const REAL A_integrand",
            "const REAL XU0_integrand",
            "const REAL XU1_integrand",
            "const REAL XU2_integrand",
            "const REAL R0_integrand",
            "const REAL XRU0_integrand",
            "const REAL XRU1_integrand",
            "const REAL XRU2_integrand",
            "const REAL O0_integrand",
            "const REAL XOU0_integrand",
            "const REAL XOU1_integrand",
            "const REAL XOU2_integrand",
            "const REAL ZOU0_integrand",
            "const REAL ZOU1_integrand",
            "const REAL ZOU2_integrand",
            "const REAL Oabs_integrand",
        ]

        sympy_expressions = [
            area_density,
            integrands_dict["area_integrand"],  # This is just 1.
            *integrands_dict["measurement_frame_xU"],
            integrands_dict["ricci_scalar"],
            *integrands_dict["x_times_R_integrand"],
            integrands_dict["spin_function"],
            *integrands_dict["xOmega_momentU"],
            *integrands_dict["zOmegaU"],
            integrands_dict["abs_omega_integrand"],
        ]

        ricci_c_code = ccg.c_codegen(
            [area_density, integrands_dict["ricci_scalar"]],
            ["const REAL spin_area_density", "const REAL spin_ricci_scalar"],
            enable_fd_codegen=True,
            enable_fd_functions=enable_fd_functions,
        )
        fd_order = int(par.parval_from_str("finite_difference::fd_order"))
        supported_fd_orders = (2, 4, 6, 8)
        if fd_order not in supported_fd_orders:
            raise ValueError(
                f"SpECTRE spin-potential solve supports centered fd_order in {list(supported_fd_orders)}, got {fd_order}"
            )

        fd_radius = fd_order // 2
        fd_width = fd_order + 1
        fd_matrix_inverse = fin.setup_FD_matrix__return_inverse(fd_width, 0)
        fd_first_coeffs = ", ".join(
            f"{float(fd_matrix_inverse[i, 1]):.17e}" for i in range(fd_width)
        )
        fd_second_coeffs = ", ".join(
            f"{float(2 * fd_matrix_inverse[i, 2]):.17e}" for i in range(fd_width)
        )
        max_row_nnz = 2 * fd_width**2 + 4 * fd_width

        prefunc = APPLY_PARITY_BRANCHLESS_PREFUNC + rf"""
#define NUM_SPECTRE_SPIN_SCRATCH_GFS {len(_SPECTRE_SPIN_SCRATCH_GFS)}
{scratch_gf_defines}

static const int8_t spectre_spin_scratch_gf_parity[NUM_SPECTRE_SPIN_SCRATCH_GFS] = {{
{parity_table_entries}
}};

/**
 * SE_XD and SE_qDD are surface coordinate-basis covariant components on the
 * horizon angular grid, not ambient unit-basis tensor components. Therefore
 * they must be transformed with bc->deriv_jacobian rather than the usual
 * unit-vector parity table.
 *
 * @param[in] gfs Spin diagnostic scratch storage.
 * @param srcpt Source point index in flattened grid storage.
 * @param src_coord Source angular coordinate index.
 * @param spectre_spin_npoints Number of points in one scratch gridfunction.
 * @return Requested source covector component, or zero for unsupported indices.
 */
static REAL spectre_spin_src_XD(const REAL *restrict gfs, const int srcpt, const int src_coord,
                                const size_t spectre_spin_npoints) {{
  if (src_coord == 1)
    return gfs[(size_t)srcpt + spectre_spin_npoints * (size_t)SE_XD0GF];
  if (src_coord == 2)
    return gfs[(size_t)srcpt + spectre_spin_npoints * (size_t)SE_XD1GF];
  return 0.0;
}} // END FUNCTION: spectre_spin_src_XD

/**
 * Transform one surface covector component across an inner boundary map.
 *
 * @param[in] gfs Spin diagnostic scratch storage.
 * @param[in] bc Inner-boundary source/destination metadata.
 * @param dst_surface_A Destination surface-coordinate component.
 * @param spectre_spin_npoints Number of points in one scratch gridfunction.
 * @return Transformed destination covector component.
 */
static REAL spectre_spin_transform_XD(const REAL *restrict gfs, const innerpt_bc_struct *restrict bc,
                                      const int dst_surface_A, const size_t spectre_spin_npoints) {{
  const int dst_coord = dst_surface_A + 1; // 0 -> theta, 1 -> phi

  REAL val = 0.0;
  for (int src_coord = 1; src_coord <= 2; src_coord++) {{
    const int8_t J = bc->deriv_jacobian[dst_coord][src_coord];
    if (J != 0)
      val += (REAL)J * spectre_spin_src_XD(gfs, bc->srcpt, src_coord, spectre_spin_npoints);
  }} // END LOOP: for src_coord over surface source coordinates

  return val;
}} // END FUNCTION: spectre_spin_transform_XD

/**
 * Read one symmetric surface metric component from a source point.
 *
 * @param[in] gfs Spin diagnostic scratch storage.
 * @param srcpt Source point index in flattened grid storage.
 * @param src_coord_A First source angular coordinate index.
 * @param src_coord_B Second source angular coordinate index.
 * @param spectre_spin_npoints Number of points in one scratch gridfunction.
 * @return Requested source metric component, or zero for unsupported indices.
 */
static REAL spectre_spin_src_qDD(const REAL *restrict gfs, const int srcpt, const int src_coord_A,
                                 const int src_coord_B, const size_t spectre_spin_npoints) {{
  if (src_coord_A == 1 && src_coord_B == 1)
    return gfs[(size_t)srcpt + spectre_spin_npoints * (size_t)SE_QDD00GF];
  if ((src_coord_A == 1 && src_coord_B == 2) || (src_coord_A == 2 && src_coord_B == 1))
    return gfs[(size_t)srcpt + spectre_spin_npoints * (size_t)SE_QDD01GF];
  if (src_coord_A == 2 && src_coord_B == 2)
    return gfs[(size_t)srcpt + spectre_spin_npoints * (size_t)SE_QDD11GF];
  return 0.0;
}} // END FUNCTION: spectre_spin_src_qDD

/**
 * Transform one symmetric surface metric component across an inner boundary map.
 *
 * @param[in] gfs Spin diagnostic scratch storage.
 * @param[in] bc Inner-boundary source/destination metadata.
 * @param dst_surface_A First destination surface-coordinate component.
 * @param dst_surface_B Second destination surface-coordinate component.
 * @param spectre_spin_npoints Number of points in one scratch gridfunction.
 * @return Transformed destination metric component.
 */
static REAL spectre_spin_transform_qDD(const REAL *restrict gfs, const innerpt_bc_struct *restrict bc,
                                       const int dst_surface_A, const int dst_surface_B,
                                       const size_t spectre_spin_npoints) {{
  const int dst_coord_A = dst_surface_A + 1; // 0 -> theta, 1 -> phi
  const int dst_coord_B = dst_surface_B + 1;

  REAL val = 0.0;
  for (int src_coord_C = 1; src_coord_C <= 2; src_coord_C++) {{
    const int8_t JA = bc->deriv_jacobian[dst_coord_A][src_coord_C];
    if (JA == 0)
      continue;

    for (int src_coord_D = 1; src_coord_D <= 2; src_coord_D++) {{
      const int8_t JB = bc->deriv_jacobian[dst_coord_B][src_coord_D];
      if (JB != 0)
        val += (REAL)JA * (REAL)JB * spectre_spin_src_qDD(gfs, bc->srcpt, src_coord_C, src_coord_D, spectre_spin_npoints);
    }} // END LOOP: for src_coord_D over surface source coordinates
  }} // END LOOP: for src_coord_C over surface source coordinates

  return val;
}} // END FUNCTION: spectre_spin_transform_qDD

/**
 * Apply inner boundary conditions to selected spin scratch gridfunctions.
 *
 * The spin diagnostic precomputes SE_qDD and SE_XD on physical angular points
 * before generated finite-difference code differentiates those fields. This
 * helper fills their angular ghost zones for every active radial horizon slab.
 *
 * @param[in] bcstruct Boundary metadata for inner curvilinear points.
 * @param[in,out] spectre_spin_gfs Spin diagnostic scratch storage.
 * @param Nxx0 Number of active radial horizon slabs.
 * @param Nxx_plus_2NGHOSTS0 Radial grid size including ghost zones.
 * @param Nxx_plus_2NGHOSTS1 Theta grid size including ghost zones.
 * @param Nxx_plus_2NGHOSTS2 Phi grid size including ghost zones.
 * @param[in] which_gfs Selected spin scratch gridfunction indices.
 * @param num_gfs Number of selected gridfunctions.
 * @param[in] scratch_gf_parity Parity table indexed by spin scratch gridfunction.
 */
static void apply_inner_bc_for_selected_spectre_spin_gfs(const bc_struct *restrict bcstruct, REAL *restrict spectre_spin_gfs, const int Nxx0,
                                                         const int Nxx_plus_2NGHOSTS0, const int Nxx_plus_2NGHOSTS1,
                                                         const int Nxx_plus_2NGHOSTS2, const int *restrict which_gfs, const int num_gfs,
                                                         const int8_t *restrict scratch_gf_parity) {{
  const bc_info_struct *bc_info = &bcstruct->bc_info;
  const size_t spectre_spin_npoints =
      (size_t)Nxx_plus_2NGHOSTS0 * (size_t)Nxx_plus_2NGHOSTS1 * (size_t)Nxx_plus_2NGHOSTS2;

#pragma omp parallel for collapse(2)
  for (int gf_idx = 0; gf_idx < num_gfs; gf_idx++) {{
    for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) {{
      const int which_gf = which_gfs[gf_idx];
      const innerpt_bc_struct *restrict bc = &bcstruct->inner_bc_array[pt];
      const int dstpt = bc->dstpt;
      const int dst_i0 = dstpt % Nxx_plus_2NGHOSTS0;

      if (dst_i0 < NGHOSTS || dst_i0 >= NGHOSTS + Nxx0)
        continue;

      switch (which_gf) {{
      case SE_XD0GF:
        spectre_spin_gfs[IDX4pt(which_gf, dstpt)] = spectre_spin_transform_XD(spectre_spin_gfs, bc, 0, spectre_spin_npoints);
        break;
      case SE_XD1GF:
        spectre_spin_gfs[IDX4pt(which_gf, dstpt)] = spectre_spin_transform_XD(spectre_spin_gfs, bc, 1, spectre_spin_npoints);
        break;
      case SE_QDD00GF:
        spectre_spin_gfs[IDX4pt(which_gf, dstpt)] = spectre_spin_transform_qDD(spectre_spin_gfs, bc, 0, 0, spectre_spin_npoints);
        break;
      case SE_QDD01GF:
        spectre_spin_gfs[IDX4pt(which_gf, dstpt)] = spectre_spin_transform_qDD(spectre_spin_gfs, bc, 0, 1, spectre_spin_npoints);
        break;
      case SE_QDD11GF:
        spectre_spin_gfs[IDX4pt(which_gf, dstpt)] = spectre_spin_transform_qDD(spectre_spin_gfs, bc, 1, 1, spectre_spin_npoints);
        break;
      default: {{
        const int8_t p = bc->parity[scratch_gf_parity[which_gf]];
        spectre_spin_gfs[IDX4pt(which_gf, dstpt)] = apply_parity_branchless(spectre_spin_gfs[IDX4pt(which_gf, bc->srcpt)], p);
        break;
      }} // END DEFAULT: apply scalar parity to selected scratch gridfunction
      }} // END SWITCH: selected spin scratch gridfunction
    }} // END LOOP: for pt over inner boundary points
  }} // END LOOP: for selected spin scratch gridfunctions
}} // END FUNCTION: apply_inner_bc_for_selected_spectre_spin_gfs

/**
 * Check selected scratch gridfunctions for finite values over an index box.
 *
 * @param[in] spectre_spin_gfs Spin diagnostic scratch storage.
 * @param Nxx0 Number of active radial horizon slabs.
 * @param Nxx_plus_2NGHOSTS0 Radial grid size including ghost zones.
 * @param Nxx_plus_2NGHOSTS1 Theta grid size including ghost zones.
 * @param Nxx_plus_2NGHOSTS2 Phi grid size including ghost zones.
 * @param[in] which_gfs Selected spin scratch gridfunction indices.
 * @param num_gfs Number of selected gridfunctions.
 * @param i1_min Minimum theta index to check.
 * @param i1_max One-past-maximum theta index to check.
 * @param i2_min Minimum phi index to check.
 * @param i2_max One-past-maximum phi index to check.
 * @param[in] stage Diagnostic stage name for warning output.
 * @return BHAHAHA_SUCCESS, or a geometry error if a value is non-finite.
 */
static int spectre_spin_check_finite_scratch_gfs(const REAL *restrict spectre_spin_gfs, const int Nxx0, const int Nxx_plus_2NGHOSTS0,
                                                 const int Nxx_plus_2NGHOSTS1, const int Nxx_plus_2NGHOSTS2,
                                                 const int *restrict which_gfs, const int num_gfs, const int i1_min,
                                                 const int i1_max, const int i2_min, const int i2_max, const char *restrict stage) {{
  if (NGHOSTS < 0 || NGHOSTS + Nxx0 > Nxx_plus_2NGHOSTS0 ||
      i1_min < 0 || i2_min < 0 || i1_max > Nxx_plus_2NGHOSTS1 || i2_max > Nxx_plus_2NGHOSTS2)
    return DIAG_SPECTRE_SPIN_POTENTIAL_GEOMETRY_ERROR;

  for (int gf_idx = 0; gf_idx < num_gfs; gf_idx++) {{
    const int gf = which_gfs[gf_idx];
    for (int i2 = i2_min; i2 < i2_max; i2++) {{
      for (int i1 = i1_min; i1 < i1_max; i1++) {{
        for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {{
          const REAL val = spectre_spin_gfs[IDX4(gf, i0, i1, i2)];
          if (!isfinite(val)) {{
            fprintf(stderr,
                    "WARNING: SpECTRE spin scratch check failed after %s: "
                    "gf=%d i0=%d i1=%d i2=%d value=%+.17e\n",
                    stage, gf, i0, i1, i2, (double)val);
            return DIAG_SPECTRE_SPIN_POTENTIAL_GEOMETRY_ERROR;
          }} // END IF: scratch-gridfunction value is not finite
        }} // END LOOP: i0
      }} // END LOOP: i1
    }} // END LOOP: i2
  }} // END LOOP: gf_idx
  return BHAHAHA_SUCCESS;
}} // END FUNCTION: spectre_spin_check_finite_scratch_gfs
"""
        prefunc += (
            r"""
#include "akv_primme.h"

#ifndef PRIMME_VERSION_MAJOR
#define BHAHAHA_PRIMME_USES_LEGACY_MATVEC 1
#elif PRIMME_VERSION_MAJOR < 3
#define BHAHAHA_PRIMME_USES_LEGACY_MATVEC 1
#endif

#define SPECTRE_SPIN_FD_RADIUS @FD_RADIUS@
#define SPECTRE_SPIN_FD_WIDTH @FD_WIDTH@
#define SPECTRE_SPIN_MAX_ROW_NNZ @MAX_ROW_NNZ@

static const REAL spectre_spin_fd_first[SPECTRE_SPIN_FD_WIDTH] = {@FD_FIRST@};
static const REAL spectre_spin_fd_second[SPECTRE_SPIN_FD_WIDTH] = {@FD_SECOND@};

typedef struct {
  int col;
  REAL val;
} spectre_spin_row_entry;

typedef struct {
  int n;
  spectre_spin_row_entry e[SPECTRE_SPIN_MAX_ROW_NNZ];
} spectre_spin_sparse_row;

typedef struct {
  int row;
  int col;
  REAL val;
} spectre_spin_triplet;

typedef struct {
  int rows;
  int cols;
  int nnz;
  int capacity;
  spectre_spin_triplet *restrict entries;
} spectre_spin_triplet_builder;

typedef struct {
  int rows;
  int cols;
  int nnz;
  int *restrict rowptr;
  int *restrict colind;
  REAL *restrict vals;
} spectre_spin_csr_matrix;

typedef struct {
  const spectre_spin_csr_matrix *restrict K;
  const spectre_spin_csr_matrix *restrict M;
  const REAL *restrict mu;
  const int *restrict red_to_full;
  const int *restrict full_to_red;
  int nfull;
  int nred;
  int anchor;
  REAL mu_anchor;
  REAL *restrict full_x;
  REAL *restrict full_y;
} spectre_spin_primme_ctx;

static void spectre_spin_row_clear(spectre_spin_sparse_row *restrict row) {
  row->n = 0;
} // END FUNCTION: spectre_spin_row_clear

/**
 * Add or accumulate one sparse-row matrix entry.
 *
 * @param[in,out] row Sparse row being assembled.
 * @param col Column index to add.
 * @param val Entry value to add.
 * @return BHAHAHA_SUCCESS, or an error code if the row is full.
 */
static int spectre_spin_row_add(spectre_spin_sparse_row *restrict row, const int col, const REAL val) {
  if (val == 0.0)
    return BHAHAHA_SUCCESS;
  for (int i = 0; i < row->n; i++) {
    if (row->e[i].col == col) {
      row->e[i].val += val;
      return BHAHAHA_SUCCESS;
    } // END IF: sparse-row entry already exists
  } // END LOOP: for i over sparse-row entries
  if (row->n >= SPECTRE_SPIN_MAX_ROW_NNZ)
    return DIAG_SPECTRE_SPIN_POTENTIAL_GEOMETRY_ERROR;
  row->e[row->n].col = col;
  row->e[row->n].val = val;
  row->n++;
  return BHAHAHA_SUCCESS;
} // END FUNCTION: spectre_spin_row_add

/**
 * Remove zero or non-finite entries from a sparse row in place.
 *
 * @param[in,out] row Sparse row to prune.
 */
static void spectre_spin_row_prune(spectre_spin_sparse_row *restrict row) {
  int out = 0;
  for (int i = 0; i < row->n; i++) {
    if (row->e[i].val != 0.0 && isfinite(row->e[i].val)) {
      row->e[out++] = row->e[i];
    } // END IF: row entry is nonzero and finite
  } // END LOOP: for i over sparse-row entries
  row->n = out;
} // END FUNCTION: spectre_spin_row_prune

/**
 * Initialize a triplet sparse-matrix builder.
 *
 * @param[in,out] builder Triplet builder to initialize.
 * @param rows Number of matrix rows.
 * @param cols Number of matrix columns.
 * @param initial_capacity Requested initial entry capacity.
 * @return BHAHAHA_SUCCESS, or an allocation error code.
 */
static int spectre_spin_builder_init(spectre_spin_triplet_builder *restrict builder, const int rows, const int cols,
                                     const int initial_capacity) {
  builder->rows = rows;
  builder->cols = cols;
  builder->nnz = 0;
  builder->capacity = initial_capacity > 0 ? initial_capacity : 1024;
  builder->entries = (spectre_spin_triplet *)malloc((size_t)builder->capacity * sizeof(spectre_spin_triplet));
  return builder->entries == NULL ? DIAG_SPECTRE_SPIN_POTENTIAL_MALLOC_ERROR : BHAHAHA_SUCCESS;
} // END FUNCTION: spectre_spin_builder_init

/**
 * Free storage owned by a triplet sparse-matrix builder.
 *
 * @param[in,out] builder Triplet builder to release.
 */
static void spectre_spin_builder_free(spectre_spin_triplet_builder *restrict builder) {
  free(builder->entries);
  builder->entries = NULL;
  builder->nnz = 0;
  builder->capacity = 0;
} // END FUNCTION: spectre_spin_builder_free

/**
 * Append one finite matrix entry to a triplet sparse-matrix builder.
 *
 * @param[in,out] builder Triplet builder being assembled.
 * @param row Row index for the entry.
 * @param col Column index for the entry.
 * @param val Entry value.
 * @return BHAHAHA_SUCCESS, or an allocation/geometry error code.
 */
static int spectre_spin_builder_add(spectre_spin_triplet_builder *restrict builder, const int row, const int col, const REAL val) {
  if (val == 0.0)
    return BHAHAHA_SUCCESS;
  if (row < 0 || row >= builder->rows || col < 0 || col >= builder->cols || !isfinite(val))
    return DIAG_SPECTRE_SPIN_POTENTIAL_GEOMETRY_ERROR;
  if (builder->nnz == builder->capacity) {
    const int new_capacity = builder->capacity < 1048576 ? 2 * builder->capacity : builder->capacity + 1048576;
    spectre_spin_triplet *restrict new_entries =
        (spectre_spin_triplet *)realloc(builder->entries, (size_t)new_capacity * sizeof(spectre_spin_triplet));
    if (new_entries == NULL)
      return DIAG_SPECTRE_SPIN_POTENTIAL_MALLOC_ERROR;
    builder->entries = new_entries;
    builder->capacity = new_capacity;
  } // END IF: triplet builder storage is full
  builder->entries[builder->nnz].row = row;
  builder->entries[builder->nnz].col = col;
  builder->entries[builder->nnz].val = val;
  builder->nnz++;
  return BHAHAHA_SUCCESS;
} // END FUNCTION: spectre_spin_builder_add

static int spectre_spin_triplet_cmp(const void *a, const void *b) {
  const spectre_spin_triplet *ta = (const spectre_spin_triplet *)a;
  const spectre_spin_triplet *tb = (const spectre_spin_triplet *)b;
  if (ta->row != tb->row)
    return ta->row < tb->row ? -1 : 1;
  if (ta->col != tb->col)
    return ta->col < tb->col ? -1 : 1;
  return 0;
} // END FUNCTION: spectre_spin_triplet_cmp

/**
 * Convert assembled triplets to compressed sparse row storage.
 *
 * Duplicate triplet entries are sorted, summed, and zero/non-finite sums are
 * discarded before the CSR arrays are allocated and populated.
 *
 * @param[in,out] builder Triplet builder containing entries to compress.
 * @param[in,out] csr Output CSR matrix.
 * @return BHAHAHA_SUCCESS, or an allocation error code.
 */
static int spectre_spin_builder_to_csr(spectre_spin_triplet_builder *restrict builder, spectre_spin_csr_matrix *restrict csr) {
  csr->rows = builder->rows;
  csr->cols = builder->cols;
  csr->nnz = 0;
  csr->rowptr = NULL;
  csr->colind = NULL;
  csr->vals = NULL;

  qsort(builder->entries, (size_t)builder->nnz, sizeof(spectre_spin_triplet), spectre_spin_triplet_cmp);

  int compressed_nnz = 0;
  for (int i = 0; i < builder->nnz;) {
    const int row = builder->entries[i].row;
    const int col = builder->entries[i].col;
    REAL val = 0.0;
    while (i < builder->nnz && builder->entries[i].row == row && builder->entries[i].col == col) {
      val += builder->entries[i].val;
      i++;
    } // END WHILE: accumulate duplicate triplet entries
    if (val != 0.0 && isfinite(val)) {
      builder->entries[compressed_nnz].row = row;
      builder->entries[compressed_nnz].col = col;
      builder->entries[compressed_nnz].val = val;
      compressed_nnz++;
    } // END IF: compressed triplet entry is nonzero and finite
  } // END LOOP: for i over sorted triplet entries

  csr->rowptr = (int *)calloc((size_t)csr->rows + 1, sizeof(int));
  csr->colind = (int *)malloc((size_t)compressed_nnz * sizeof(int));
  csr->vals = (REAL *)malloc((size_t)compressed_nnz * sizeof(REAL));
  if (csr->rowptr == NULL || csr->colind == NULL || csr->vals == NULL) {
    free(csr->rowptr);
    free(csr->colind);
    free(csr->vals);
    csr->rowptr = NULL;
    csr->colind = NULL;
    csr->vals = NULL;
    return DIAG_SPECTRE_SPIN_POTENTIAL_MALLOC_ERROR;
  } // END IF: CSR array allocation failed

  for (int i = 0; i < compressed_nnz; i++)
    csr->rowptr[builder->entries[i].row + 1]++;
  for (int row = 0; row < csr->rows; row++)
    csr->rowptr[row + 1] += csr->rowptr[row];
  int *restrict cursor = (int *)malloc((size_t)csr->rows * sizeof(int));
  if (cursor == NULL) {
    free(csr->rowptr);
    free(csr->colind);
    free(csr->vals);
    csr->rowptr = NULL;
    csr->colind = NULL;
    csr->vals = NULL;
    return DIAG_SPECTRE_SPIN_POTENTIAL_MALLOC_ERROR;
  } // END IF: CSR row cursor allocation failed
  for (int row = 0; row < csr->rows; row++)
    cursor[row] = csr->rowptr[row];
  for (int i = 0; i < compressed_nnz; i++) {
    const int row = builder->entries[i].row;
    const int dst = cursor[row]++;
    csr->colind[dst] = builder->entries[i].col;
    csr->vals[dst] = builder->entries[i].val;
  } // END LOOP: for i over compressed triplet entries
  free(cursor);
  csr->nnz = compressed_nnz;
  return BHAHAHA_SUCCESS;
} // END FUNCTION: spectre_spin_builder_to_csr

/**
 * Free storage owned by a CSR matrix.
 *
 * @param[in,out] csr CSR matrix to release.
 */
static void spectre_spin_csr_free(spectre_spin_csr_matrix *restrict csr) {
  free(csr->rowptr);
  free(csr->colind);
  free(csr->vals);
  csr->rowptr = NULL;
  csr->colind = NULL;
  csr->vals = NULL;
  csr->nnz = 0;
} // END FUNCTION: spectre_spin_csr_free

/**
 * Apply a CSR matrix to a vector.
 *
 * @param[in] csr Matrix in compressed sparse row storage.
 * @param[in] x Input vector.
 * @param[out] y Output vector.
 */
static void spectre_spin_csr_matvec(const spectre_spin_csr_matrix *restrict csr, const REAL *restrict x, REAL *restrict y) {
  for (int row = 0; row < csr->rows; row++) {
    REAL sum = 0.0;
    for (int jj = csr->rowptr[row]; jj < csr->rowptr[row + 1]; jj++)
      sum += csr->vals[jj] * x[csr->colind[jj]];
    y[row] = sum;
  } // END LOOP: for row over CSR matrix rows
} // END FUNCTION: spectre_spin_csr_matvec

/**
 * Seed reduced-space eigenvectors from centered coordinate functions.
 *
 * @param N Number of full-space grid points.
 * @param Nred Number of reduced-space grid points.
 * @param[in] red_to_full Map from reduced-space indices to full-space indices.
 * @param[in] x_ref Reference coordinate functions.
 * @param[in] x_centroid Area centroid.
 * @param[out] evecs_red Initial reduced-space eigenvector guesses.
 */
static void spectre_spin_seed_coordinate_reduced(const int N, const int Nred, const int *restrict red_to_full, const REAL *restrict x_ref,
                                                 const REAL x_centroid[3], double *restrict evecs_red) {
  for (int mode = 0; mode < 3; mode++) {
    for (int r = 0; r < Nred; r++) {
      const int full = red_to_full[r];
      evecs_red[(size_t)mode * (size_t)Nred + r] = (double)(x_ref[(size_t)mode * (size_t)N + full] - x_centroid[mode]);
    } // END LOOP: for r over reduced-space grid points
  } // END LOOP: for mode over coordinate seed modes
} // END FUNCTION: spectre_spin_seed_coordinate_reduced

/**
 * Expand a reduced mean-zero vector into full-space coordinates.
 *
 * @param[in] ctx PRIMME matrix-vector context.
 * @param[in] xred Reduced-space vector.
 * @param[out] xfull Full-space vector satisfying the weighted mean constraint.
 */
static void spectre_spin_expand_reduced(const spectre_spin_primme_ctx *restrict ctx, const double *restrict xred, REAL *restrict xfull) {
  for (int i = 0; i < ctx->nfull; i++)
    xfull[i] = 0.0;
  REAL anchor_sum = 0.0;
  for (int r = 0; r < ctx->nred; r++) {
    const int full = ctx->red_to_full[r];
    xfull[full] = (REAL)xred[r];
    anchor_sum += ctx->mu[full] * xfull[full];
  } // END LOOP: for r over reduced-space entries
  xfull[ctx->anchor] = -anchor_sum / ctx->mu_anchor;
} // END FUNCTION: spectre_spin_expand_reduced

/**
 * Project a full-space operator result back into reduced coordinates.
 *
 * @param[in] ctx PRIMME matrix-vector context.
 * @param[in] yfull Full-space vector.
 * @param[out] yred Reduced-space projected vector.
 */
static void spectre_spin_project_reduced(const spectre_spin_primme_ctx *restrict ctx, const REAL *restrict yfull, double *restrict yred) {
  const REAL y_anchor = yfull[ctx->anchor];
  for (int r = 0; r < ctx->nred; r++) {
    const int full = ctx->red_to_full[r];
    yred[r] = (double)(yfull[full] - (ctx->mu[full] / ctx->mu_anchor) * y_anchor);
  } // END LOOP: for r over reduced-space entries
} // END FUNCTION: spectre_spin_project_reduced

/**
 * Apply one reduced operator to a block of PRIMME vectors.
 *
 * @param[in] ctx PRIMME matrix-vector context.
 * @param[in] A Full-space CSR operator to apply.
 * @param[in] x Reduced-space input block.
 * @param ldx Leading dimension of the input block.
 * @param[out] y Reduced-space output block.
 * @param ldy Leading dimension of the output block.
 * @param block_size Number of vectors in the block.
 */
static void spectre_spin_apply_reduced_operator(const spectre_spin_primme_ctx *restrict ctx, const spectre_spin_csr_matrix *restrict A,
                                                const double *restrict x, const PRIMME_INT ldx, double *restrict y, const PRIMME_INT ldy,
                                                const int block_size) {
  for (int block = 0; block < block_size; block++) {
    const double *restrict xcol = x + (size_t)block * (size_t)ldx;
    double *restrict ycol = y + (size_t)block * (size_t)ldy;
    spectre_spin_expand_reduced(ctx, xcol, ctx->full_x);
    spectre_spin_csr_matvec(A, ctx->full_x, ctx->full_y);
    spectre_spin_project_reduced(ctx, ctx->full_y, ycol);
  } // END LOOP: for block over PRIMME block vectors
} // END FUNCTION: spectre_spin_apply_reduced_operator

#ifdef BHAHAHA_PRIMME_USES_LEGACY_MATVEC
static void spectre_spin_primme_K_matvec(void *x, void *y, int *blockSize, primme_params *primme, int *err) {
  spectre_spin_primme_ctx *restrict ctx = (spectre_spin_primme_ctx *)primme->matrix;
  spectre_spin_apply_reduced_operator(ctx, ctx->K, (const double *)x, (PRIMME_INT)ctx->nred, (double *)y, (PRIMME_INT)ctx->nred, *blockSize);
  *err = 0;
} // END FUNCTION: spectre_spin_primme_K_matvec

static void spectre_spin_primme_M_matvec(void *x, void *y, int *blockSize, primme_params *primme, int *err) {
  spectre_spin_primme_ctx *restrict ctx = (spectre_spin_primme_ctx *)primme->massMatrix;
  spectre_spin_apply_reduced_operator(ctx, ctx->M, (const double *)x, (PRIMME_INT)ctx->nred, (double *)y, (PRIMME_INT)ctx->nred, *blockSize);
  *err = 0;
} // END FUNCTION: spectre_spin_primme_M_matvec
#else
static void spectre_spin_primme_K_matvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme,
                                         int *err) {
  spectre_spin_primme_ctx *restrict ctx = (spectre_spin_primme_ctx *)primme->matrix;
  spectre_spin_apply_reduced_operator(ctx, ctx->K, (const double *)x, *ldx, (double *)y, *ldy, *blockSize);
  *err = 0;
} // END FUNCTION: spectre_spin_primme_K_matvec

static void spectre_spin_primme_M_matvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme,
                                         int *err) {
  spectre_spin_primme_ctx *restrict ctx = (spectre_spin_primme_ctx *)primme->massMatrix;
  spectre_spin_apply_reduced_operator(ctx, ctx->M, (const double *)x, *ldx, (double *)y, *ldy, *blockSize);
  *err = 0;
} // END FUNCTION: spectre_spin_primme_M_matvec
#endif

static int spectre_spin_active_index(const int j1, const int j2, const int Ntheta) {
  return j1 + Ntheta * j2;
} // END FUNCTION: spectre_spin_active_index

/**
 * Reflect a scalar angular-grid index through the spherical-polar axis.
 *
 * @param[in,out] j1 Theta-like angular index.
 * @param[in,out] j2 Phi-like angular index.
 * @param Ntheta Number of theta points.
 * @param Nphi Number of phi points.
 */
static void spectre_spin_reflect_scalar_index(int *restrict j1, int *restrict j2, const int Ntheta, const int Nphi) {
  while (*j1 < 0 || *j1 >= Ntheta) {
    if (*j1 < 0)
      *j1 = -*j1 - 1;
    else
      *j1 = 2 * Ntheta - *j1 - 1;
    *j2 += Nphi / 2;
  } // END WHILE: theta index lies outside the active angular domain
  *j2 %= Nphi;
  if (*j2 < 0)
    *j2 += Nphi;
} // END FUNCTION: spectre_spin_reflect_scalar_index

/**
 * Build one finite-difference derivative row for a scalar on the horizon grid.
 *
 * @param[in,out] row Sparse row receiving the finite-difference stencil.
 * @param j1 Theta-like angular point.
 * @param j2 Phi-like angular point.
 * @param Ntheta Number of theta points.
 * @param Nphi Number of phi points.
 * @param deriv_type Derivative selector.
 * @param invdtheta Inverse theta spacing.
 * @param invdphi Inverse phi spacing.
 * @return BHAHAHA_SUCCESS, or an error code from sparse-row assembly.
 */
static int spectre_spin_build_scalar_derivative_row(spectre_spin_sparse_row *restrict row, const int j1, const int j2, const int Ntheta,
                                                    const int Nphi, const int deriv_type, const REAL invdtheta, const REAL invdphi) {
  spectre_spin_row_clear(row);
  if (deriv_type == 0 || deriv_type == 2) {
    const REAL scale = deriv_type == 0 ? invdtheta : invdtheta * invdtheta;
    const REAL *restrict coeffs = deriv_type == 0 ? spectre_spin_fd_first : spectre_spin_fd_second;
    for (int s = -SPECTRE_SPIN_FD_RADIUS; s <= SPECTRE_SPIN_FD_RADIUS; s++) {
      int jj1 = j1 + s;
      int jj2 = j2;
      spectre_spin_reflect_scalar_index(&jj1, &jj2, Ntheta, Nphi);
      const int col = spectre_spin_active_index(jj1, jj2, Ntheta);
      const int status = spectre_spin_row_add(row, col, coeffs[s + SPECTRE_SPIN_FD_RADIUS] * scale);
      if (status != BHAHAHA_SUCCESS)
        return status;
    } // END LOOP: for s over theta derivative stencil offsets
  } else if (deriv_type == 1 || deriv_type == 3) {
    const REAL scale = deriv_type == 1 ? invdphi : invdphi * invdphi;
    const REAL *restrict coeffs = deriv_type == 1 ? spectre_spin_fd_first : spectre_spin_fd_second;
    for (int s = -SPECTRE_SPIN_FD_RADIUS; s <= SPECTRE_SPIN_FD_RADIUS; s++) {
      int jj1 = j1;
      int jj2 = j2 + s;
      spectre_spin_reflect_scalar_index(&jj1, &jj2, Ntheta, Nphi);
      const int col = spectre_spin_active_index(jj1, jj2, Ntheta);
      const int status = spectre_spin_row_add(row, col, coeffs[s + SPECTRE_SPIN_FD_RADIUS] * scale);
      if (status != BHAHAHA_SUCCESS)
        return status;
    } // END LOOP: for s over phi derivative stencil offsets
  } else {
    const REAL scale = invdtheta * invdphi;
    for (int s1 = -SPECTRE_SPIN_FD_RADIUS; s1 <= SPECTRE_SPIN_FD_RADIUS; s1++) {
      for (int s2 = -SPECTRE_SPIN_FD_RADIUS; s2 <= SPECTRE_SPIN_FD_RADIUS; s2++) {
        int jj1 = j1 + s1;
        int jj2 = j2 + s2;
        spectre_spin_reflect_scalar_index(&jj1, &jj2, Ntheta, Nphi);
        const int col = spectre_spin_active_index(jj1, jj2, Ntheta);
        const REAL val = spectre_spin_fd_first[s1 + SPECTRE_SPIN_FD_RADIUS] * spectre_spin_fd_first[s2 + SPECTRE_SPIN_FD_RADIUS] * scale;
        const int status = spectre_spin_row_add(row, col, val);
        if (status != BHAHAHA_SUCCESS)
          return status;
      } // END LOOP: for s2 over phi mixed-derivative stencil offsets
    } // END LOOP: for s1 over theta mixed-derivative stencil offsets
  } // END ELSE: build mixed theta-phi derivative row
  spectre_spin_row_prune(row);
  return BHAHAHA_SUCCESS;
} // END FUNCTION: spectre_spin_build_scalar_derivative_row

/**
 * Compute one angular derivative of a scratch metric gridfunction.
 *
 * @param[in] gfs Spin scratch gridfunction storage.
 * @param gf Gridfunction index to differentiate.
 * @param i0 Radial grid index.
 * @param i1 Theta grid index.
 * @param i2 Phi grid index.
 * @param dir Angular derivative direction.
 * @param Nxx_plus_2NGHOSTS0 Radial grid size including ghost zones.
 * @param Nxx_plus_2NGHOSTS1 Theta grid size including ghost zones.
 * @param Nxx_plus_2NGHOSTS2 Phi grid size including ghost zones.
 * @param invdtheta Inverse theta spacing.
 * @param invdphi Inverse phi spacing.
 * @return Angular derivative of the requested gridfunction.
 */
static REAL spectre_spin_metric_deriv(const REAL *restrict gfs, const int gf, const int i0, const int i1, const int i2,
                                      const int dir, const int Nxx_plus_2NGHOSTS0, const int Nxx_plus_2NGHOSTS1,
                                      const int Nxx_plus_2NGHOSTS2, const REAL invdtheta, const REAL invdphi) {
  REAL deriv = 0.0;
  if (dir == 0) {
    for (int s = -SPECTRE_SPIN_FD_RADIUS; s <= SPECTRE_SPIN_FD_RADIUS; s++) {
      const int idx = i0 + Nxx_plus_2NGHOSTS0 * ((i1 + s) + Nxx_plus_2NGHOSTS1 * (i2 + Nxx_plus_2NGHOSTS2 * gf));
      deriv += spectre_spin_fd_first[s + SPECTRE_SPIN_FD_RADIUS] * gfs[idx];
    } // END LOOP: for s over theta derivative stencil offsets
    deriv *= invdtheta;
  } else {
    for (int s = -SPECTRE_SPIN_FD_RADIUS; s <= SPECTRE_SPIN_FD_RADIUS; s++) {
      const int idx = i0 + Nxx_plus_2NGHOSTS0 * (i1 + Nxx_plus_2NGHOSTS1 * ((i2 + s) + Nxx_plus_2NGHOSTS2 * gf));
      deriv += spectre_spin_fd_first[s + SPECTRE_SPIN_FD_RADIUS] * gfs[idx];
    } // END LOOP: for s over phi derivative stencil offsets
    deriv *= invdphi;
  } // END ELSE: differentiate in phi direction
  return deriv;
} // END FUNCTION: spectre_spin_metric_deriv

/**
 * Add a scaled sparse outer product to a triplet matrix builder.
 *
 * @param[in,out] builder Triplet builder receiving the outer product.
 * @param[in] a Left sparse row.
 * @param[in] b Right sparse row.
 * @param factor Scalar multiplier.
 * @return BHAHAHA_SUCCESS, or an error code from triplet assembly.
 */
static int spectre_spin_add_outer(spectre_spin_triplet_builder *restrict builder, const spectre_spin_sparse_row *restrict a,
                                  const spectre_spin_sparse_row *restrict b, const REAL factor) {
  if (factor == 0.0)
    return BHAHAHA_SUCCESS;
  for (int ia = 0; ia < a->n; ia++) {
    for (int ib = 0; ib < b->n; ib++) {
      const int status = spectre_spin_builder_add(builder, a->e[ia].col, b->e[ib].col, factor * a->e[ia].val * b->e[ib].val);
      if (status != BHAHAHA_SUCCESS)
        return status;
    } // END LOOP: for ib over right sparse-row entries
  } // END LOOP: for ia over left sparse-row entries
  return BHAHAHA_SUCCESS;
} // END FUNCTION: spectre_spin_add_outer

/**
 * Add one q^{AB} grad_A eta grad_B z bilinear form contribution.
 *
 * @param[in,out] builder Triplet builder receiving the contribution.
 * @param[in] dtheta Theta derivative row.
 * @param[in] dphi Phi derivative row.
 * @param factor Scalar quadrature and physical multiplier.
 * @param q00 Inverse surface metric theta-theta component.
 * @param q01 Inverse surface metric theta-phi component.
 * @param q11 Inverse surface metric phi-phi component.
 * @return BHAHAHA_SUCCESS, or an error code from triplet assembly.
 */
static int spectre_spin_add_gradient_form(spectre_spin_triplet_builder *restrict builder, const spectre_spin_sparse_row *restrict dtheta,
                                          const spectre_spin_sparse_row *restrict dphi, const REAL factor, const REAL q00,
                                          const REAL q01, const REAL q11) {
  int status = spectre_spin_add_outer(builder, dtheta, dtheta, factor * q00);
  if (status != BHAHAHA_SUCCESS)
    return status;
  status = spectre_spin_add_outer(builder, dtheta, dphi, factor * q01);
  if (status != BHAHAHA_SUCCESS)
    return status;
  status = spectre_spin_add_outer(builder, dphi, dtheta, factor * q01);
  if (status != BHAHAHA_SUCCESS)
    return status;
  return spectre_spin_add_outer(builder, dphi, dphi, factor * q11);
} // END FUNCTION: spectre_spin_add_gradient_form

/**
 * Compute eigenpairs of a symmetric 3x3 matrix using Jacobi rotations.
 *
 * @param[in] input Symmetric input matrix.
 * @param[out] evals Eigenvalues in ascending order.
 * @param[out] evecs Corresponding eigenvectors.
 * @return BHAHAHA_SUCCESS, or a normalization error if an eigenvalue is not finite.
 */
static int spectre_spin_jacobi_eigen_3x3(const REAL input[3][3], REAL evals[3], REAL evecs[3][3]) {
  REAL a[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      a[i][j] = input[i][j];
      evecs[i][j] = i == j ? 1.0 : 0.0;
    } // END LOOP: for j over matrix columns
  } // END LOOP: for i over matrix rows
  for (int iter = 0; iter < 64; iter++) {
    int p = 0;
    int q = 1;
    REAL max_offdiag = fabs(a[0][1]);
    if (fabs(a[0][2]) > max_offdiag) {
      p = 0;
      q = 2;
      max_offdiag = fabs(a[0][2]);
    } // END IF: a[0][2] is the largest candidate off-diagonal entry
    if (fabs(a[1][2]) > max_offdiag) {
      p = 1;
      q = 2;
      max_offdiag = fabs(a[1][2]);
    } // END IF: a[1][2] is the largest candidate off-diagonal entry
    if (max_offdiag < 1.0e-13)
      break;
    const REAL tau = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
    const REAL t = copysign(1.0 / (fabs(tau) + sqrt(1.0 + tau * tau)), tau);
    const REAL c = 1.0 / sqrt(1.0 + t * t);
    const REAL s = t * c;
    const REAL app = a[p][p];
    const REAL aqq = a[q][q];
    const REAL apq = a[p][q];
    a[p][p] = app - t * apq;
    a[q][q] = aqq + t * apq;
    a[p][q] = 0.0;
    a[q][p] = 0.0;
    for (int k = 0; k < 3; k++) {
      if (k != p && k != q) {
        const REAL akp = a[k][p];
        const REAL akq = a[k][q];
        a[k][p] = c * akp - s * akq;
        a[p][k] = a[k][p];
        a[k][q] = s * akp + c * akq;
        a[q][k] = a[k][q];
      } // END IF: k is outside the active Jacobi pivot pair
      const REAL vkp = evecs[k][p];
      const REAL vkq = evecs[k][q];
      evecs[k][p] = c * vkp - s * vkq;
      evecs[k][q] = s * vkp + c * vkq;
    } // END LOOP: for k over Jacobi-rotation rows
  } // END LOOP: for iter over Jacobi sweeps
  for (int i = 0; i < 3; i++) {
    evals[i] = a[i][i];
    if (!isfinite(evals[i]))
      return DIAG_SPECTRE_SPIN_POTENTIAL_NORMALIZATION_ERROR;
  } // END LOOP: for i over eigenvalues
  for (int i = 0; i < 2; i++) {
    for (int j = i + 1; j < 3; j++) {
      if (evals[j] < evals[i]) {
        const REAL tmp_eval = evals[i];
        evals[i] = evals[j];
        evals[j] = tmp_eval;
        for (int k = 0; k < 3; k++) {
          const REAL tmp_vec = evecs[k][i];
          evecs[k][i] = evecs[k][j];
          evecs[k][j] = tmp_vec;
        } // END LOOP: for k over eigenvector components
      } // END IF: eigenvalues are out of order
    } // END LOOP: for j over later eigenvalues
  } // END LOOP: for i over eigenvalues to sort
  return BHAHAHA_SUCCESS;
} // END FUNCTION: spectre_spin_jacobi_eigen_3x3

static REAL spectre_spin_det3(const REAL A[3][3]) {
  return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
         A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
         A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
} // END FUNCTION: spectre_spin_det3

/**
 * Compute the orthogonal Procrustes rotation from a 3x3 correlation matrix.
 *
 * @param[in] C Correlation matrix between numerical and reference modes.
 * @param[out] O Orientation-preserving orthogonal alignment matrix.
 * @return BHAHAHA_SUCCESS, or an error code if the alignment is singular.
 */
static int spectre_spin_procrustes(const REAL C[3][3], REAL O[3][3]) {
  REAL CtC[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        CtC[i][j] += C[k][i] * C[k][j];

  REAL evals[3];
  REAL V[3][3];
  int status = spectre_spin_jacobi_eigen_3x3(CtC, evals, V);
  if (status != BHAHAHA_SUCCESS)
    return status;
  for (int i = 0; i < 3; i++) {
    if (!(evals[i] > 0.0) || !isfinite(evals[i]))
      return DIAG_SPECTRE_SPIN_POTENTIAL_NORMALIZATION_ERROR;
  } // END LOOP: for i over Procrustes singular values

  REAL invsqrt[3] = {1.0 / sqrt(evals[0]), 1.0 / sqrt(evals[1]), 1.0 / sqrt(evals[2])};
  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      O[a][b] = 0.0;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          O[a][b] += C[a][i] * V[i][j] * invsqrt[j] * V[b][j];
    } // END LOOP: for b over output matrix columns
  } // END LOOP: for a over output matrix rows
  if (spectre_spin_det3(O) < 0.0) {
    invsqrt[0] = -invsqrt[0];
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        O[a][b] = 0.0;
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            O[a][b] += C[a][i] * V[i][j] * invsqrt[j] * V[b][j];
      } // END LOOP: for b over orientation-corrected output matrix columns
    } // END LOOP: for a over orientation-corrected output matrix rows
  } // END IF: Procrustes rotation needs orientation correction
  return BHAHAHA_SUCCESS;
} // END FUNCTION: spectre_spin_procrustes

/**
 * Compute normalized scalar spin-potential modes z_alpha for the SpECTRE-style
 * spin diagnostic.
 *
 * This diagnostic owns ZU0GF, ZU1GF, and ZU2GF. It solves the Owen/Beetle
 * scalar AKV generalized eigenproblem in symmetric weak form:
 *
 *   K(eta,z) = int (Delta eta)(Delta z) dA - int R grad eta . grad z dA
 *   M(eta,z) = int grad eta . grad z dA
 *   K z = lambda M z
 *
 * Constants are removed by an area-weighted mean-zero reduced space. The three
 * lowest modes are aligned to measurement-frame reference potentials and
 * normalized with the scalar-potential convention:
 * int (z_alpha - <z_alpha>)^2 dA = A^3 / (48*pi^2).
 */
static int bah_compute_spectre_spin_potentials(commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                                               const REAL *restrict auxevol_gfs, REAL *restrict spectre_spin_gfs) {
  const int grid = 0;
  const params_struct *restrict params = &griddata[grid].params;
  const REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];
#include "set_CodeParameters.h"

#ifdef NDEBUG
  const REAL z_init = 0.0;
#else
  const REAL z_init = (REAL)NAN;
#endif

#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
    for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
      for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
        spectre_spin_gfs[IDX4(ZU0GF, i0, i1, i2)] = z_init;
        spectre_spin_gfs[IDX4(ZU1GF, i0, i1, i2)] = z_init;
        spectre_spin_gfs[IDX4(ZU2GF, i0, i1, i2)] = z_init;
      }
    }
  }

  if (Nxx0 != 1)
    return DIAG_SPECTRE_SPIN_POTENTIAL_GEOMETRY_ERROR;
  if (Nxx2 % 2 != 0)
    return DIAG_SPECTRE_SPIN_POTENTIAL_GEOMETRY_ERROR;

  const int Ntheta = Nxx1;
  const int Nphi = Nxx2;
  const int N = Ntheta * Nphi;
  const int Nred = N - 1;
  if (Nred < 3)
    return DIAG_SPECTRE_SPIN_POTENTIAL_GEOMETRY_ERROR;

  const REAL *restrict weights;
  int weight_stencil_size;
  bah_diagnostics_integration_weights(Nxx1, Nxx2, &weights, &weight_stencil_size);

  REAL *restrict mu = (REAL *)malloc((size_t)N * sizeof(REAL));
  REAL *restrict sqrtq = (REAL *)malloc((size_t)N * sizeof(REAL));
  REAL *restrict ricci = (REAL *)malloc((size_t)N * sizeof(REAL));
  REAL *restrict qUU00 = (REAL *)malloc((size_t)N * sizeof(REAL));
  REAL *restrict qUU01 = (REAL *)malloc((size_t)N * sizeof(REAL));
  REAL *restrict qUU11 = (REAL *)malloc((size_t)N * sizeof(REAL));
  REAL *restrict x_ref = (REAL *)malloc((size_t)3 * (size_t)N * sizeof(REAL));
  int *restrict red_to_full = (int *)malloc((size_t)Nred * sizeof(int));
  int *restrict full_to_red = (int *)malloc((size_t)N * sizeof(int));
  double *restrict evals = (double *)malloc(3 * sizeof(double));
  double *restrict evecs_red = (double *)malloc((size_t)Nred * 3 * sizeof(double));
  double *restrict resnorms = (double *)malloc(3 * sizeof(double));
  REAL *restrict evecs_full = (REAL *)malloc((size_t)3 * (size_t)N * sizeof(REAL));
  REAL *restrict modes = (REAL *)malloc((size_t)3 * (size_t)N * sizeof(REAL));
  if (mu == NULL || sqrtq == NULL || ricci == NULL || qUU00 == NULL || qUU01 == NULL || qUU11 == NULL || x_ref == NULL ||
      red_to_full == NULL || full_to_red == NULL || evals == NULL || evecs_red == NULL || resnorms == NULL || evecs_full == NULL ||
      modes == NULL) {
    free(mu);
    free(sqrtq);
    free(ricci);
    free(qUU00);
    free(qUU01);
    free(qUU11);
    free(x_ref);
    free(red_to_full);
    free(full_to_red);
    free(evals);
    free(evecs_red);
    free(resnorms);
    free(evecs_full);
    free(modes);
    return DIAG_SPECTRE_SPIN_POTENTIAL_MALLOC_ERROR;
  }

  const int i0 = NGHOSTS;
  const REAL invdtheta = 1.0 / dxx1;
  const REAL invdphi = 1.0 / dxx2;
  const REAL surface_weight = dxx1 * dxx2;
  REAL area = 0.0;
  REAL x_centroid[3] = {0.0, 0.0, 0.0};
  int anchor = 0;
  REAL mu_anchor = 0.0;

  for (int j2 = 0; j2 < Nphi; j2++) {
    const int i2 = NGHOSTS + j2;
    const REAL weight2 = weights[j2 % weight_stencil_size];
    const REAL phi = xx[2][i2];
    const REAL cos_phi = cos(phi);
    const REAL sin_phi = sin(phi);
    for (int j1 = 0; j1 < Ntheta; j1++) {
      const int i1 = NGHOSTS + j1;
      const int p = spectre_spin_active_index(j1, j2, Ntheta);
      const REAL weight1 = weights[j1 % weight_stencil_size];
      const REAL theta = xx[1][i1];
      const REAL sin_theta = sin(theta);
      const REAL cos_theta = cos(theta);

@RICCI_CODE@

      const REAL q00 = spectre_spin_gfs[IDX4(SE_QDD00GF, i0, i1, i2)];
      const REAL q01 = spectre_spin_gfs[IDX4(SE_QDD01GF, i0, i1, i2)];
      const REAL q11 = spectre_spin_gfs[IDX4(SE_QDD11GF, i0, i1, i2)];
      const REAL detq = q00 * q11 - q01 * q01;
      if (!(detq > 0.0) || !isfinite(detq) || !(spin_area_density > 0.0) || !isfinite(spin_area_density) ||
          !isfinite(spin_ricci_scalar)) {
        free(mu);
        free(sqrtq);
        free(ricci);
        free(qUU00);
        free(qUU01);
        free(qUU11);
        free(x_ref);
        free(red_to_full);
        free(full_to_red);
        free(evals);
        free(evecs_red);
        free(resnorms);
        free(evecs_full);
        free(modes);
        return DIAG_SPECTRE_SPIN_POTENTIAL_GEOMETRY_ERROR;
      }
      sqrtq[p] = spin_area_density;
      ricci[p] = spin_ricci_scalar;
      qUU00[p] = q11 / detq;
      qUU01[p] = -q01 / detq;
      qUU11[p] = q00 / detq;
      mu[p] = spin_area_density * weight1 * weight2 * surface_weight;
      if (!(mu[p] > 0.0) || !isfinite(mu[p])) {
        free(mu);
        free(sqrtq);
        free(ricci);
        free(qUU00);
        free(qUU01);
        free(qUU11);
        free(x_ref);
        free(red_to_full);
        free(full_to_red);
        free(evals);
        free(evecs_red);
        free(resnorms);
        free(evecs_full);
        free(modes);
        return DIAG_SPECTRE_SPIN_POTENTIAL_GEOMETRY_ERROR;
      }
      const REAL hh = in_gfs[IDX4(HHGF, i0, i1, i2)];
      x_ref[0 * N + p] = hh * sin_theta * cos_phi;
      x_ref[1 * N + p] = hh * sin_theta * sin_phi;
      x_ref[2 * N + p] = hh * cos_theta;
      area += mu[p];
      for (int a = 0; a < 3; a++)
        x_centroid[a] += x_ref[a * N + p] * mu[p];
      if (mu[p] > mu_anchor) {
        mu_anchor = mu[p];
        anchor = p;
      }
    }
  }
  if (!(area > 0.0) || !(mu_anchor > 0.0) || !isfinite(area) || !isfinite(mu_anchor)) {
    free(mu);
    free(sqrtq);
    free(ricci);
    free(qUU00);
    free(qUU01);
    free(qUU11);
    free(x_ref);
    free(red_to_full);
    free(full_to_red);
    free(evals);
    free(evecs_red);
    free(resnorms);
    free(evecs_full);
    free(modes);
    return DIAG_SPECTRE_SPIN_POTENTIAL_GEOMETRY_ERROR;
  }
  for (int a = 0; a < 3; a++)
    x_centroid[a] /= area;

  for (int p = 0, r = 0; p < N; p++) {
    full_to_red[p] = -1;
    if (p != anchor) {
      full_to_red[p] = r;
      red_to_full[r] = p;
      r++;
    }
  }

  spectre_spin_triplet_builder K_builder;
  spectre_spin_triplet_builder M_builder;
  const int initial_capacity = N * 512;
  int status = spectre_spin_builder_init(&K_builder, N, N, initial_capacity);
  if (status == BHAHAHA_SUCCESS)
    status = spectre_spin_builder_init(&M_builder, N, N, initial_capacity);
  if (status != BHAHAHA_SUCCESS) {
    if (K_builder.entries != NULL)
      spectre_spin_builder_free(&K_builder);
    free(mu);
    free(sqrtq);
    free(ricci);
    free(qUU00);
    free(qUU01);
    free(qUU11);
    free(x_ref);
    free(red_to_full);
    free(full_to_red);
    free(evals);
    free(evecs_red);
    free(resnorms);
    free(evecs_full);
    free(modes);
    return status;
  }

  spectre_spin_sparse_row dtheta_row, dphi_row, dthetatheta_row, dthetaphi_row, dphiphi_row, lap_row;
  for (int j2 = 0; j2 < Nphi; j2++) {
    const int i2 = NGHOSTS + j2;
    for (int j1 = 0; j1 < Ntheta; j1++) {
      const int i1 = NGHOSTS + j1;
      const int p = spectre_spin_active_index(j1, j2, Ntheta);
      status = spectre_spin_build_scalar_derivative_row(&dtheta_row, j1, j2, Ntheta, Nphi, 0, invdtheta, invdphi);
      if (status == BHAHAHA_SUCCESS)
        status = spectre_spin_build_scalar_derivative_row(&dphi_row, j1, j2, Ntheta, Nphi, 1, invdtheta, invdphi);
      if (status == BHAHAHA_SUCCESS)
        status = spectre_spin_build_scalar_derivative_row(&dthetatheta_row, j1, j2, Ntheta, Nphi, 2, invdtheta, invdphi);
      if (status == BHAHAHA_SUCCESS)
        status = spectre_spin_build_scalar_derivative_row(&dphiphi_row, j1, j2, Ntheta, Nphi, 3, invdtheta, invdphi);
      if (status == BHAHAHA_SUCCESS)
        status = spectre_spin_build_scalar_derivative_row(&dthetaphi_row, j1, j2, Ntheta, Nphi, 4, invdtheta, invdphi);
      if (status != BHAHAHA_SUCCESS)
        break;

      const REAL dq00_dtheta = spectre_spin_metric_deriv(spectre_spin_gfs, SE_QDD00GF, i0, i1, i2, 0, Nxx_plus_2NGHOSTS0,
                                                         Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, invdtheta, invdphi);
      const REAL dq01_dtheta = spectre_spin_metric_deriv(spectre_spin_gfs, SE_QDD01GF, i0, i1, i2, 0, Nxx_plus_2NGHOSTS0,
                                                         Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, invdtheta, invdphi);
      const REAL dq11_dtheta = spectre_spin_metric_deriv(spectre_spin_gfs, SE_QDD11GF, i0, i1, i2, 0, Nxx_plus_2NGHOSTS0,
                                                         Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, invdtheta, invdphi);
      const REAL dq00_dphi = spectre_spin_metric_deriv(spectre_spin_gfs, SE_QDD00GF, i0, i1, i2, 1, Nxx_plus_2NGHOSTS0,
                                                       Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, invdtheta, invdphi);
      const REAL dq01_dphi = spectre_spin_metric_deriv(spectre_spin_gfs, SE_QDD01GF, i0, i1, i2, 1, Nxx_plus_2NGHOSTS0,
                                                       Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, invdtheta, invdphi);
      const REAL dq11_dphi = spectre_spin_metric_deriv(spectre_spin_gfs, SE_QDD11GF, i0, i1, i2, 1, Nxx_plus_2NGHOSTS0,
                                                       Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, invdtheta, invdphi);

      const REAL d_q[2][2][2] = {
          {{dq00_dtheta, dq00_dphi}, {dq01_dtheta, dq01_dphi}},
          {{dq01_dtheta, dq01_dphi}, {dq11_dtheta, dq11_dphi}}};
      const REAL qinv[2][2] = {{qUU00[p], qUU01[p]}, {qUU01[p], qUU11[p]}};
      REAL Gamma[2][2][2] = {{{0.0, 0.0}, {0.0, 0.0}}, {{0.0, 0.0}, {0.0, 0.0}}};
      for (int C = 0; C < 2; C++) {
        for (int A = 0; A < 2; A++) {
          for (int B = 0; B < 2; B++) {
            for (int D = 0; D < 2; D++)
              Gamma[C][A][B] += 0.5 * qinv[C][D] * (d_q[D][B][A] + d_q[A][D][B] - d_q[A][B][D]);
          }
        }
      }
      const REAL gradtheta_coeff = -(qUU00[p] * Gamma[0][0][0] + 2.0 * qUU01[p] * Gamma[0][0][1] + qUU11[p] * Gamma[0][1][1]);
      const REAL gradphi_coeff = -(qUU00[p] * Gamma[1][0][0] + 2.0 * qUU01[p] * Gamma[1][0][1] + qUU11[p] * Gamma[1][1][1]);

      spectre_spin_row_clear(&lap_row);
      for (int n = 0; n < dthetatheta_row.n; n++)
        status = spectre_spin_row_add(&lap_row, dthetatheta_row.e[n].col, qUU00[p] * dthetatheta_row.e[n].val);
      for (int n = 0; status == BHAHAHA_SUCCESS && n < dthetaphi_row.n; n++)
        status = spectre_spin_row_add(&lap_row, dthetaphi_row.e[n].col, 2.0 * qUU01[p] * dthetaphi_row.e[n].val);
      for (int n = 0; status == BHAHAHA_SUCCESS && n < dphiphi_row.n; n++)
        status = spectre_spin_row_add(&lap_row, dphiphi_row.e[n].col, qUU11[p] * dphiphi_row.e[n].val);
      for (int n = 0; status == BHAHAHA_SUCCESS && n < dtheta_row.n; n++)
        status = spectre_spin_row_add(&lap_row, dtheta_row.e[n].col, gradtheta_coeff * dtheta_row.e[n].val);
      for (int n = 0; status == BHAHAHA_SUCCESS && n < dphi_row.n; n++)
        status = spectre_spin_row_add(&lap_row, dphi_row.e[n].col, gradphi_coeff * dphi_row.e[n].val);
      if (status != BHAHAHA_SUCCESS)
        break;
      spectre_spin_row_prune(&lap_row);

      status = spectre_spin_add_gradient_form(&M_builder, &dtheta_row, &dphi_row, mu[p], qUU00[p], qUU01[p], qUU11[p]);
      if (status == BHAHAHA_SUCCESS)
        status = spectre_spin_add_outer(&K_builder, &lap_row, &lap_row, mu[p]);
      if (status == BHAHAHA_SUCCESS)
        status = spectre_spin_add_gradient_form(&K_builder, &dtheta_row, &dphi_row, -mu[p] * ricci[p], qUU00[p], qUU01[p], qUU11[p]);
      if (status != BHAHAHA_SUCCESS)
        break;
    }
    if (status != BHAHAHA_SUCCESS)
      break;
  }
  if (status != BHAHAHA_SUCCESS) {
    spectre_spin_builder_free(&K_builder);
    spectre_spin_builder_free(&M_builder);
    free(mu);
    free(sqrtq);
    free(ricci);
    free(qUU00);
    free(qUU01);
    free(qUU11);
    free(x_ref);
    free(red_to_full);
    free(full_to_red);
    free(evals);
    free(evecs_red);
    free(resnorms);
    free(evecs_full);
    free(modes);
    return status;
  }

  spectre_spin_csr_matrix K_csr, M_csr;
  status = spectre_spin_builder_to_csr(&K_builder, &K_csr);
  if (status == BHAHAHA_SUCCESS)
    status = spectre_spin_builder_to_csr(&M_builder, &M_csr);
  spectre_spin_builder_free(&K_builder);
  spectre_spin_builder_free(&M_builder);
  if (status != BHAHAHA_SUCCESS) {
    if (K_csr.rowptr != NULL)
      spectre_spin_csr_free(&K_csr);
    free(mu);
    free(sqrtq);
    free(ricci);
    free(qUU00);
    free(qUU01);
    free(qUU11);
    free(x_ref);
    free(red_to_full);
    free(full_to_red);
    free(evals);
    free(evecs_red);
    free(resnorms);
    free(evecs_full);
    free(modes);
    return status;
  }

  REAL *restrict full_x = (REAL *)malloc((size_t)N * sizeof(REAL));
  REAL *restrict full_y = (REAL *)malloc((size_t)N * sizeof(REAL));
  if (full_x == NULL || full_y == NULL) {
    spectre_spin_csr_free(&K_csr);
    spectre_spin_csr_free(&M_csr);
    free(full_x);
    free(full_y);
    free(mu);
    free(sqrtq);
    free(ricci);
    free(qUU00);
    free(qUU01);
    free(qUU11);
    free(x_ref);
    free(red_to_full);
    free(full_to_red);
    free(evals);
    free(evecs_red);
    free(resnorms);
    free(evecs_full);
    free(modes);
    return DIAG_SPECTRE_SPIN_POTENTIAL_MALLOC_ERROR;
  }

  spectre_spin_primme_ctx ctx = {
      .K = &K_csr,
      .M = &M_csr,
      .mu = mu,
      .red_to_full = red_to_full,
      .full_to_red = full_to_red,
      .nfull = N,
      .nred = Nred,
      .anchor = anchor,
      .mu_anchor = mu_anchor,
      .full_x = full_x,
      .full_y = full_y};

  primme_params primme;
  primme_initialize(&primme);
  primme.n = Nred;
  primme.nLocal = Nred;
  primme.numEvals = 3;

  
  // Target the smallest-magnitude eigenvalues near zero.
  static double target_shift = 0.0;
  primme.target = primme_closest_abs;
  primme.numTargetShifts = 1;
  primme.targetShifts = &target_shift;
  
  primme.matrixMatvec = spectre_spin_primme_K_matvec;
  primme.massMatrixMatvec = spectre_spin_primme_M_matvec;
  primme.matrixMatvec_type = primme_op_double;
  primme.massMatrixMatvec_type = primme_op_double;
  primme.ldOPs = Nred;
  primme.ldevecs = Nred;
  primme.matrix = &ctx;
  primme.massMatrix = &ctx;
  primme.eps = 1.0e-6;
  primme.maxMatvecs = 50000;
  primme.maxOuterIterations = 5000;
  primme.printLevel = 1;
  primme_set_method(PRIMME_DEFAULT_MIN_TIME, &primme);
  primme.initSize = 3;
  primme.maxBasisSize = 60;
  primme.minRestartSize = 20;
  primme.maxBlockSize = 4;
  spectre_spin_seed_coordinate_reduced(N, Nred, red_to_full, x_ref, x_centroid, evecs_red);

  const int primme_status = dprimme(evals, evecs_red, resnorms, &primme);
  primme_free(&primme);
  if (primme_status != 0) {
    spectre_spin_csr_free(&K_csr);
    spectre_spin_csr_free(&M_csr);
    free(full_x);
    free(full_y);
    free(mu);
    free(sqrtq);
    free(ricci);
    free(qUU00);
    free(qUU01);
    free(qUU11);
    free(x_ref);
    free(red_to_full);
    free(full_to_red);
    free(evals);
    free(evecs_red);
    free(resnorms);
    free(evecs_full);
    free(modes);
    return DIAG_SPECTRE_SPIN_POTENTIAL_PRIMME_ERROR;
  }
  for (int a = 0; a < 3; a++) {
    if (!isfinite(evals[a]) || !isfinite(resnorms[a])) {
      spectre_spin_csr_free(&K_csr);
      spectre_spin_csr_free(&M_csr);
      free(full_x);
      free(full_y);
      free(mu);
      free(sqrtq);
      free(ricci);
      free(qUU00);
      free(qUU01);
      free(qUU11);
      free(x_ref);
      free(red_to_full);
      free(full_to_red);
      free(evals);
      free(evecs_red);
      free(resnorms);
      free(evecs_full);
      free(modes);
      return DIAG_SPECTRE_SPIN_POTENTIAL_PRIMME_ERROR;
    }
  }

  for (int mode = 0; mode < 3; mode++) {
    spectre_spin_expand_reduced(&ctx, &evecs_red[(size_t)mode * (size_t)Nred], &evecs_full[(size_t)mode * (size_t)N]);
  }

  const REAL area_radius = sqrt(area / (4.0 * M_PI));
  for (int p = 0; p < N; p++)
    for (int a = 0; a < 3; a++)
      x_ref[a * N + p] = area_radius * (x_ref[a * N + p] - x_centroid[a]);

  REAL C[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  for (int a = 0; a < 3; a++) {
    spectre_spin_csr_matvec(&M_csr, &x_ref[a * N], full_y);
    for (int mode = 0; mode < 3; mode++) {
      REAL overlap = 0.0;
      for (int p = 0; p < N; p++)
        overlap += full_y[p] * evecs_full[mode * N + p];
      C[a][mode] = overlap;
    }
  }
  REAL O[3][3];
  status = spectre_spin_procrustes(C, O);
  if (status != BHAHAHA_SUCCESS) {
    spectre_spin_csr_free(&K_csr);
    spectre_spin_csr_free(&M_csr);
    free(full_x);
    free(full_y);
    free(mu);
    free(sqrtq);
    free(ricci);
    free(qUU00);
    free(qUU01);
    free(qUU11);
    free(x_ref);
    free(red_to_full);
    free(full_to_red);
    free(evals);
    free(evecs_red);
    free(resnorms);
    free(evecs_full);
    free(modes);
    return status;
  }

  for (int a = 0; a < 3; a++) {
    for (int p = 0; p < N; p++) {
      modes[a * N + p] = 0.0;
      for (int mode = 0; mode < 3; mode++)
        modes[a * N + p] += O[a][mode] * evecs_full[mode * N + p];
    }
  }

  const REAL target_potential_norm = area * area * area / (48.0 * M_PI * M_PI);
  if (!(target_potential_norm > 0.0) || !isfinite(target_potential_norm))
    status = DIAG_SPECTRE_SPIN_POTENTIAL_NORMALIZATION_ERROR;
  for (int a = 0; status == BHAHAHA_SUCCESS && a < 3; a++) {
    REAL mean = 0.0;
    for (int p = 0; p < N; p++)
      mean += mu[p] * modes[a * N + p];
    mean /= area;
    REAL norm = 0.0;
    for (int p = 0; p < N; p++) {
      const REAL centered = modes[a * N + p] - mean;
      norm += mu[p] * centered * centered;
    }
    if (!(norm > 0.0) || !isfinite(norm)) {
      status = DIAG_SPECTRE_SPIN_POTENTIAL_NORMALIZATION_ERROR;
      break;
    }
    const REAL scale = sqrt(target_potential_norm / norm);
    if (!isfinite(scale)) {
      status = DIAG_SPECTRE_SPIN_POTENTIAL_NORMALIZATION_ERROR;
      break;
    }
    for (int p = 0; p < N; p++)
      modes[a * N + p] = scale * (modes[a * N + p] - mean);
  }

  if (status == BHAHAHA_SUCCESS) {
#pragma omp parallel for
    for (int j2 = 0; j2 < Nphi; j2++) {
      for (int j1 = 0; j1 < Ntheta; j1++) {
        const int p = spectre_spin_active_index(j1, j2, Ntheta);
        const int ii1 = NGHOSTS + j1;
        const int ii2 = NGHOSTS + j2;
        for (int ii0 = NGHOSTS; ii0 < NGHOSTS + Nxx0; ii0++) {
          spectre_spin_gfs[IDX4(ZU0GF, ii0, ii1, ii2)] = modes[0 * N + p];
          spectre_spin_gfs[IDX4(ZU1GF, ii0, ii1, ii2)] = modes[1 * N + p];
          spectre_spin_gfs[IDX4(ZU2GF, ii0, ii1, ii2)] = modes[2 * N + p];
        }
      }
    }

    int finite_error = 0;
#pragma omp parallel for reduction(| : finite_error)
    for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
      for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
        for (int ii0 = NGHOSTS; ii0 < NGHOSTS + Nxx0; ii0++) {
          finite_error |= !isfinite(spectre_spin_gfs[IDX4(ZU0GF, ii0, i1, i2)]);
          finite_error |= !isfinite(spectre_spin_gfs[IDX4(ZU1GF, ii0, i1, i2)]);
          finite_error |= !isfinite(spectre_spin_gfs[IDX4(ZU2GF, ii0, i1, i2)]);
        }
      }
    }
    if (finite_error)
      status = DIAG_SPECTRE_SPIN_POTENTIAL_NORMALIZATION_ERROR;
  }

  spectre_spin_csr_free(&K_csr);
  spectre_spin_csr_free(&M_csr);
  free(full_x);
  free(full_y);
  free(mu);
  free(sqrtq);
  free(ricci);
  free(qUU00);
  free(qUU01);
  free(qUU11);
  free(x_ref);
  free(red_to_full);
  free(full_to_red);
  free(evals);
  free(evecs_red);
  free(resnorms);
  free(evecs_full);
  free(modes);
  return status;
} // END FUNCTION: bah_compute_spectre_spin_potentials
""".replace("@FD_RADIUS@", str(fd_radius))
            .replace("@FD_WIDTH@", str(fd_width))
            .replace("@MAX_ROW_NNZ@", str(max_row_nnz))
            .replace("@FD_FIRST@", fd_first_coeffs)
            .replace("@FD_SECOND@", fd_second_coeffs)
            .replace("@RICCI_CODE@", ricci_c_code.rstrip())
        )

        # Step 6: Construct the body of the C function.
        body = r"""
    const int grid=0;
    const params_struct *restrict params = &griddata[grid].params;
    REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    const REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs; // for hh and its time-derivs
    REAL *restrict xx[3];
    for(int ww=0;ww<3;ww++) xx[ww] = griddata[grid].xx[ww];
#include "set_CodeParameters.h"

    const size_t spectre_spin_npoints =
        (size_t)Nxx_plus_2NGHOSTS0 * (size_t)Nxx_plus_2NGHOSTS1 * (size_t)Nxx_plus_2NGHOSTS2;
    REAL *restrict spectre_spin_gfs =
        (REAL *)malloc((size_t)NUM_SPECTRE_SPIN_SCRATCH_GFS * spectre_spin_npoints * sizeof(REAL));
    if (spectre_spin_gfs == NULL)
        return DIAG_SPECTRE_SPIN_POTENTIAL_MALLOC_ERROR;
    for (size_t idx = 0; idx < (size_t)NUM_SPECTRE_SPIN_SCRATCH_GFS * spectre_spin_npoints; idx++)
        spectre_spin_gfs[idx] = (REAL)NAN;

    // This diagnostic owns private ZU0GF, ZU1GF, and ZU2GF scratch slots. After
    // surface precompute and ghost-zone fill, bah_compute_spectre_spin_potentials()
    // writes a normalized spin-potential basis before z_alpha * Omega is evaluated.
    
    // Initialize all RunSums accumulators to zero.
    REAL A_sum = 0.0;
    REAL XU_sum[3] = {0.0, 0.0, 0.0};
    REAL R0_sum = 0.0;
    REAL XRU_sum[3] = {0.0, 0.0, 0.0};
    REAL O0_sum = 0.0;
    REAL XOU_sum[3] = {0.0, 0.0, 0.0};
    REAL ZOU_sum[3] = {0.0, 0.0, 0.0};
    REAL Oabs_sum = 0.0;
    
    // Set integration weights (e.g., for Simpson's rule).
    // This external function provides the 1D weights array.
    const REAL *restrict weights;
    int weight_stencil_size;
    bah_diagnostics_integration_weights(Nxx1, Nxx2, &weights, &weight_stencil_size);

    // Precompute SE_qDD and SE_XD only where generated finite-difference
    // stencils are valid.
#pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
        const REAL xx2 = xx[2][i2];
        (void)xx2;
        for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
            const REAL xx1 = xx[1][i1];
            (void)xx1;
            for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
"""
        body += precompute_c_code
        body += rf"""
            }} // END LOOP: for i0 over active radial horizon-grid slabs
        }} // END LOOP: for i1 over physical theta horizon-grid points
    }} // END LOOP: for i2 over physical phi horizon-grid points

    {{
        const int spectre_spin_precompute_gfs[{len(gf_macros)}] = {{{selected_precompute_gfs}}};
        int spectre_spin_scratch_status = spectre_spin_check_finite_scratch_gfs(
            spectre_spin_gfs, Nxx0, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1,
            Nxx_plus_2NGHOSTS2, spectre_spin_precompute_gfs, {len(gf_macros)},
            NGHOSTS, NGHOSTS + Nxx1, NGHOSTS, NGHOSTS + Nxx2, "physical precompute");
        if (spectre_spin_scratch_status != BHAHAHA_SUCCESS) {{
            free(spectre_spin_gfs);
            return spectre_spin_scratch_status;
        }}

        apply_inner_bc_for_selected_spectre_spin_gfs(
            &griddata[grid].bcstruct, spectre_spin_gfs, Nxx0, Nxx_plus_2NGHOSTS0,
            Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2,
            spectre_spin_precompute_gfs, {len(gf_macros)}, spectre_spin_scratch_gf_parity);

        spectre_spin_scratch_status = spectre_spin_check_finite_scratch_gfs(
            spectre_spin_gfs, Nxx0, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1,
            Nxx_plus_2NGHOSTS2, spectre_spin_precompute_gfs, {len(gf_macros)},
            NGHOSTS - {fd_order // 2}, NGHOSTS + Nxx1 + {fd_order // 2},
            NGHOSTS - {fd_order // 2}, NGHOSTS + Nxx2 + {fd_order // 2}, "ghost-zone fill");
        if (spectre_spin_scratch_status != BHAHAHA_SUCCESS) {{
            free(spectre_spin_gfs);
            return spectre_spin_scratch_status;
        }}
    }} // END BLOCK: fill SE_qDD and SE_XD ghost zones before differentiating them

    const int spin_potential_status =
        bah_compute_spectre_spin_potentials(commondata, griddata, auxevol_gfs, spectre_spin_gfs);
    if (spin_potential_status != BHAHAHA_SUCCESS) {{
        free(spectre_spin_gfs);
        return spin_potential_status;
    }}

    {{
        const int spectre_spin_z_gfs[3] = {{ZU0GF, ZU1GF, ZU2GF}};
        const int spectre_spin_z_status = spectre_spin_check_finite_scratch_gfs(
            spectre_spin_gfs, Nxx0, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1,
            Nxx_plus_2NGHOSTS2, spectre_spin_z_gfs, 3,
            NGHOSTS, NGHOSTS + Nxx1, NGHOSTS, NGHOSTS + Nxx2, "spin-potential solve");
        if (spectre_spin_z_status != BHAHAHA_SUCCESS) {{
            free(spectre_spin_gfs);
            return spectre_spin_z_status;
        }}
    }}

#pragma omp parallel
{{
    // Private accumulators for each thread
    REAL A_sum_private = 0.0;
    REAL XU_sum_private[3] = {{0.0, 0.0, 0.0}};
    REAL R0_sum_private = 0.0;
    REAL XRU_sum_private[3] = {{0.0, 0.0, 0.0}};
    REAL O0_sum_private = 0.0;
    REAL XOU_sum_private[3] = {{0.0, 0.0, 0.0}};
    REAL ZOU_sum_private[3] = {{0.0, 0.0, 0.0}};
    REAL Oabs_sum_private = 0.0;

#pragma omp for
    for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {{
        const REAL weight2 = weights[(i2 - NGHOSTS) % weight_stencil_size];
        const REAL xx2 = xx[2][i2];
        (void)xx2;
        for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {{
            const REAL weight1 = weights[(i1 - NGHOSTS) % weight_stencil_size];
            const REAL xx1 = xx[1][i1];
            (void)xx1;
            for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {{
"""
        # Step 7: Generate C code for all integrands and the area density.
        # enable_fd_codegen=True tells c_codegen to automatically handle all
        # finite difference derivatives of gridfunctions.
        body += ccg.c_codegen(
            sympy_expressions,
            integrand_c_vars,
            enable_fd_codegen=True,
            enable_fd_functions=enable_fd_functions,
        )

        body += r"""
                // The differential area element, excluding coordinate steps (dθ, dφ)
                const REAL dA_unscaled = area_density * weight1 * weight2;
                
                // Accumulate into thread-private variables
                A_sum_private += A_integrand * dA_unscaled;
                XU_sum_private[0] += XU0_integrand * dA_unscaled;
                XU_sum_private[1] += XU1_integrand * dA_unscaled;
                XU_sum_private[2] += XU2_integrand * dA_unscaled;
                XRU_sum_private[0] += XRU0_integrand * dA_unscaled;
                XRU_sum_private[1] += XRU1_integrand * dA_unscaled;
                XRU_sum_private[2] += XRU2_integrand * dA_unscaled;
                XOU_sum_private[0] += XOU0_integrand * dA_unscaled;
                XOU_sum_private[1] += XOU1_integrand * dA_unscaled;
                XOU_sum_private[2] += XOU2_integrand * dA_unscaled;
                ZOU_sum_private[0] += ZOU0_integrand * dA_unscaled;
                ZOU_sum_private[1] += ZOU1_integrand * dA_unscaled;
                ZOU_sum_private[2] += ZOU2_integrand * dA_unscaled;
                R0_sum_private   += R0_integrand * dA_unscaled;
                O0_sum_private   += O0_integrand * dA_unscaled;
                Oabs_sum_private += Oabs_integrand * dA_unscaled;
            } // END LOOP: for i0 over active radial horizon-grid slabs
        } // END LOOP: for i1 over physical theta horizon-grid points
    } // END LOOP: for i2 over physical phi horizon-grid points

    // Use a critical section for the final reduction from private to shared sums
    #pragma omp critical
    {
        A_sum += A_sum_private;
        for(int i=0; i<3; ++i) {
            XU_sum[i]  += XU_sum_private[i];
            XRU_sum[i] += XRU_sum_private[i];
            XOU_sum[i] += XOU_sum_private[i];
            ZOU_sum[i] += ZOU_sum_private[i];
        }
        R0_sum   += R0_sum_private;
        O0_sum   += O0_sum_private;
        Oabs_sum += Oabs_sum_private;
    } // END OMP CRITICAL: update shared spin-diagnostic sums
} // END OMP PARALLEL: integrate spin diagnostic

// Step 8: Compute the dimensionless spin vector from the integrated quantities.
bhahaha_diagnostics_struct *restrict bhahaha_diags = commondata->bhahaha_diagnostics;

const REAL surface_weight = dxx1 * dxx2;
const REAL spin_norm_tolerance = 1.0e-14;
const REAL A = A_sum * surface_weight;
const REAL XU[3] = {
    XU_sum[0] * surface_weight,
    XU_sum[1] * surface_weight,
    XU_sum[2] * surface_weight};
const REAL R0 = R0_sum * surface_weight;
const REAL XRU[3] = {
    XRU_sum[0] * surface_weight,
    XRU_sum[1] * surface_weight,
    XRU_sum[2] * surface_weight};
const REAL O0 = O0_sum * surface_weight;
const REAL XOU[3] = {
    XOU_sum[0] * surface_weight,
    XOU_sum[1] * surface_weight,
    XOU_sum[2] * surface_weight};
const REAL ZOU[3] = {
    ZOU_sum[0] * surface_weight,
    ZOU_sum[1] * surface_weight,
    ZOU_sum[2] * surface_weight};
const REAL Oabs = Oabs_sum * surface_weight;

REAL S_U[3] = {0.0, 0.0, 0.0};
REAL chi_U[3] = {0.0, 0.0, 0.0};

if (fabs(A) > spin_norm_tolerance) {
    REAL x0U[3], xRcorrU[3], IU[3], SalphaU[3];

    for (int i = 0; i < 3; i++) {
        x0U[i] = XU[i] / A;
        xRcorrU[i] = (XRU[i] - x0U[i] * R0) / (8.0 * M_PI);
        IU[i] = XOU[i] - (x0U[i] + xRcorrU[i]) * O0;
        SalphaU[i] = ZOU[i] / (8.0 * M_PI);
    } // END LOOP: for i over spatial spin-vector components

    const REAL normI = sqrt(IU[0] * IU[0] + IU[1] * IU[1] + IU[2] * IU[2]);
    const REAL Salpha_norm = sqrt(SalphaU[0] * SalphaU[0] + SalphaU[1] * SalphaU[1] + SalphaU[2] * SalphaU[2]);
    const REAL S = Salpha_norm;

    if (Oabs > spin_norm_tolerance && normI > spin_norm_tolerance) {
        for (int i = 0; i < 3; i++)
            S_U[i] = S * IU[i] / normI;
    } else if (Salpha_norm > spin_norm_tolerance) {
        for (int i = 0; i < 3; i++)
            S_U[i] = S * SalphaU[i] / Salpha_norm;
    } // END ELSE IF: use z_alpha fallback direction when nominal direction is unsafe

    const REAL M_irr_squared = A / (16.0 * M_PI);
    if (M_irr_squared > 0.0 && isfinite(M_irr_squared)) {
        const REAL M_horizon_squared = M_irr_squared + S * S / (4.0 * M_irr_squared);
        if (M_horizon_squared > 0.0 && isfinite(M_horizon_squared)) {
            for (int i = 0; i < 3; i++)
                chi_U[i] = S_U[i] / M_horizon_squared;
        } // END IF: Christodoulou mass squared is safe

    } // END IF: irreducible mass squared is safe
} // END IF: area is safe for centroid reduction

bhahaha_diags->spin_chi_x_spectre = chi_U[0];
bhahaha_diags->spin_chi_y_spectre = chi_U[1];
bhahaha_diags->spin_chi_z_spectre = chi_U[2];

free(spectre_spin_gfs);
return BHAHAHA_SUCCESS;
"""
    finally:
        for gf_name in _SPECTRE_SPIN_SCRATCH_GFS:
            saved_gf = saved_spectre_spin_gfs[gf_name]
            if saved_gf is None:
                gri.glb_gridfcs_dict.pop(gf_name, None)
            else:
                gri.glb_gridfcs_dict[gf_name] = saved_gf
    formatted_body = clang_format(body)

    # Format and register the C function using the standard helper.
    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=cfunc_name,
        params=params,
        include_CodeParameters_h=False,
        body=formatted_body,
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
