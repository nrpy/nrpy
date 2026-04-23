"""
C function registration for GRoovy TOV initial data.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Optional, Union, cast

import sympy as sp
from typing_extensions import Literal

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.reference_metric as refmetric

TOVMode = Literal[
    "hybrid",
    "tabulated_temperature",
    "tabulated_constant_entropy",
]


def _build_apply_grhayl_limits_code(declare_grhayl_structs: bool) -> str:
    """
    Build the common GRHayL limiting and writeback block.

    :param declare_grhayl_structs: Whether to declare the primitive, metric, and
        conservative structs in this block.
    :return: C code that rebuilds, limits, and writes back GRHayL variables.
    """
    declarations = ""
    if declare_grhayl_structs:
        declarations = r"""
      ghl_primitive_quantities prims;
      ghl_metric_quantities metric;
      ghl_conservative_quantities cons;
"""

    return declarations + r"""
      // Rebuild GRHayL structs from the TOV-seeded gridfunction data.
      basis_transform_rfm_basis_to_Cartesian(
          commondata, params, &prims, &cons, &metric, i0, i1, i2, xx,
          auxevol_gfs, in_gfs);

      bool speed_limited;
      ghl_enforce_primitive_limits_and_compute_u0(
          &commondata->ghl_params, &commondata->eos, &metric, &prims,
          &speed_limited);

      // Store the limited primitive variables and reconstructed conservatives.
      basis_transform_Cartesian_to_rfm_basis(
          commondata, params, &prims, &cons, i0, i1, i2, xx,
          auxevol_gfs, in_gfs);

      auxevol_gfs[IDX4(U4UTGF, i0, i1, i2)] = prims.u0;
"""


def _build_constant_entropy_prefunc() -> str:
    """
    Build the helper C routine used to read a constant-entropy EOS slice.

    :return: C helper-function code emitted into the `prefunc` region.
    """
    return r"""
/**
 * Read a fixed-entropy, three-column EOS slice from disk.
 *
 * @param[in] filename Path to the slice file.
 * @param[in] num_rows Number of rows expected in the file.
 * @return Dynamically allocated array with shape [num_rows][3], or NULL on failure.
 */
static double (*eos_slice_file(const char *filename, int num_rows))[3] {
  FILE *file_ptr = NULL;
  double (*data_array)[3] = malloc(num_rows * sizeof(double[3]));
  char line_buffer[256];

  if (num_rows <= 0) {
    fprintf(stderr, "Error: Number of rows must be positive.\n");
    return NULL;
  } // END IF: invalid row count
  if (data_array == NULL) {
    fprintf(stderr, "Error: Could not allocate memory for %d rows.\n", num_rows);
    return NULL;
  } // END IF: allocation failure

  file_ptr = fopen(filename, "r");
  if (file_ptr == NULL) {
    fprintf(stderr, "Error opening file '%s': %s\n", filename, strerror(errno));
    free(data_array);
    return NULL;
  } // END IF: failed file open

  for (int i = 0; i < num_rows; ++i) {
    if (fgets(line_buffer, sizeof(line_buffer), file_ptr) == NULL) {
      fprintf(stderr, "Error: Unexpected EOF or read error at row %d.\n", i + 1);
      free(data_array);
      fclose(file_ptr);
      return NULL;
    } // END IF: failed line read

    if (sscanf(
            line_buffer, "%lf %lf %lf",
            &data_array[i][0], &data_array[i][1], &data_array[i][2]) != 3) {
      fprintf(stderr, "Error: Could not parse 3 doubles on line %d.\n", i + 1);
      free(data_array);
      fclose(file_ptr);
      return NULL;
    } // END IF: failed line parse
  } // END LOOP: for i over EOS slice rows

  fclose(file_ptr);
  return data_array;
} // END FUNCTION: eos_slice_file
"""


def register_CFunction_TOV_initial_data(
    grhayl_setup_str: str,
    CoordSystem: str,
    OMP_collapse: int,
    enable_GoldenKernels: bool = False,
    tov_mode: TOVMode = "hybrid",
    tov_temperature: Optional[float] = None,
    tov_entropy: Optional[float] = None,
    tov_data_filename: str = "TOVdata.txt",
    entropy_slice_filename: str = "constant_ent_Beq_slice.txt",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register a routine that fills GRoovy gridfunctions with TOV initial data.

    This master routine supports the baseline hybrid-EOS TOV data as well as
    the tabulated fixed-temperature and tabulated fixed-entropy variants.

    :param grhayl_setup_str: C code used to initialize the GRHayL EOS structs.
    :param CoordSystem: The coordinate system.
    :param OMP_collapse: Level of OpenMP loop collapsing.
    :param enable_GoldenKernels: Whether to enable Golden Kernel generation.
    :param tov_mode: TOV initial-data mode. Supported values are `"hybrid"`,
        `"tabulated_temperature"`, and `"tabulated_constant_entropy"`.
    :param tov_temperature: Fixed temperature used by the tabulated-temperature
        branch.
    :param tov_entropy: Fixed entropy used by the tabulated constant-entropy
        branch.
    :param tov_data_filename: Path to the TOV radial profile file.
    :param entropy_slice_filename: Path to the constant-entropy tabulated EOS
        slice file.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If the requested TOV mode is unsupported or if its
        required parameters are inconsistent.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> mode_configs = [
    ...     (
    ...         "hybrid",
    ...         None,
    ...         None,
    ...         "tov_initial_data__openmp__Spherical__hybrid",
    ...     ),
    ...     (
    ...         "tabulated_temperature",
    ...         0.3,
    ...         None,
    ...         "tov_initial_data__openmp__Spherical__tabulated_temperature",
    ...     ),
    ...     (
    ...         "tabulated_constant_entropy",
    ...         None,
    ...         0.3,
    ...         "tov_initial_data__openmp__Spherical__tabulated_constant_entropy",
    ...     ),
    ... ]
    >>> for tov_mode, tov_temperature, tov_entropy, validation_desc in mode_configs:
    ...     cfc.CFunction_dict.clear()
    ...     with contextlib.redirect_stdout(io.StringIO()):
    ...         _ = register_CFunction_TOV_initial_data(
    ...             grhayl_setup_str="",
    ...             CoordSystem="Spherical",
    ...             OMP_collapse=1,
    ...             tov_mode=tov_mode,
    ...             tov_temperature=tov_temperature,
    ...             tov_entropy=tov_entropy,
    ...         )
    ...     generated_str = clang_format(cfc.CFunction_dict["initial_data"].full_function)
    ...     _ = validate_strings(generated_str, validation_desc, file_ext="c")
    """
    # Step 1: Validate the requested TOV mode.
    if tov_mode not in (
        "hybrid",
        "tabulated_temperature",
        "tabulated_constant_entropy",
    ):
        raise ValueError(f"Unsupported tov_mode={tov_mode}.")
    if tov_mode == "hybrid":
        if tov_temperature is not None or tov_entropy is not None:
            raise ValueError(
                "hybrid TOV mode does not accept tov_temperature or tov_entropy."
            )
    elif tov_mode == "tabulated_temperature":
        if tov_temperature is None or tov_entropy is not None:
            raise ValueError(
                "tabulated_temperature mode requires tov_temperature and does "
                "not accept tov_entropy."
            )
    elif tov_entropy is None or tov_temperature is not None:
        raise ValueError(
            "tabulated_constant_entropy mode requires tov_entropy and does not "
            "accept tov_temperature."
        )

    # Step 2: Register the function with NRPy's parallel codegen infrastructure.
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 3: Set up the C-function metadata.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Spherically symmetric TOV initial data"
    if tov_mode == "tabulated_temperature":
        desc += " with tabulated EOS data at fixed temperature"
    elif tov_mode == "tabulated_constant_entropy":
        desc += " with tabulated EOS data at fixed entropy"
    cfunc_type = "void"
    name = "initial_data"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    tov_data_filename_escaped = tov_data_filename.replace("\\", "\\\\").replace(
        '"', '\\"'
    )
    entropy_slice_filename_escaped = entropy_slice_filename.replace(
        "\\", "\\\\"
    ).replace('"', '\\"')

    # Step 4: Build the common TOV interpolation code.
    rfm = refmetric.reference_metric[CoordSystem]
    press, rho_baryon = sp.symbols("press rho_baryon", real=True)
    rescaledvU = ixp.zerorank1()

    hybrid_output_vars: List[sp.Expr] = [
        rho_baryon,
        press,
        rescaledvU[0],
        rescaledvU[1],
        rescaledvU[2],
    ]
    hybrid_grid_access: List[str] = [
        gri.BHaHGridFunction.access_gf("rhob", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("P", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledvU0", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledvU1", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledvU2", gf_array_name="auxevol_gfs"),
    ]

    common_loop_prefix = (
        r"""      // Interpolate the TOV matter profile at this grid point.
      REAL rbar, rho, rho_baryon, press, M, expnu, exp4phi;
"""
        + ccg.c_codegen(
            rfm.xxSph[0],
            "rbar",
            include_braces=False,
            enable_simd=False,
            enable_GoldenKernels=enable_GoldenKernels,
        )
        + r"""
      TOV_interpolate_1D(
          rbar, &ID_persist, &rho, &rho_baryon, &press, &M, &expnu, &exp4phi);
"""
    )

    # Step 5: Build any mode-specific helper code emitted before the main loop.
    prefunc = ""
    mode_setup = ""
    if tov_mode == "tabulated_temperature":
        assert tov_temperature is not None
        mode_setup = f"""
  const double beq_temperature = {tov_temperature};
  ghl_tabulated_compute_Ye_of_rho_beq_constant_T(
      beq_temperature, &commondata->eos);
"""
    elif tov_mode == "tabulated_constant_entropy":
        assert tov_entropy is not None
        prefunc = _build_constant_entropy_prefunc()
        mode_setup = f"""
  const char *eos_slice_filename = "{entropy_slice_filename_escaped}";
  double (*constant_ent_Beq_slice)[3] =
      eos_slice_file(eos_slice_filename, commondata->eos.N_rho);

  commondata->eos.lp_of_lr =
      (double *)malloc(sizeof(double) * commondata->eos.N_rho);
  commondata->eos.le_of_lr =
      (double *)malloc(sizeof(double) * commondata->eos.N_rho);
  commondata->eos.Ye_of_lr =
      (double *)malloc(sizeof(double) * commondata->eos.N_rho);

  for (int i = 0; i < commondata->eos.N_rho; ++i) {{
    commondata->eos.Ye_of_lr[i] = constant_ent_Beq_slice[i][0];
    commondata->eos.lp_of_lr[i] = constant_ent_Beq_slice[i][1];
    commondata->eos.le_of_lr[i] = constant_ent_Beq_slice[i][2];
  }} // END LOOP: for i over EOS entropy-slice rows

  free(constant_ent_Beq_slice);
  const double beq_entropy = {tov_entropy};
"""

    # Step 6: Build the main function body and import the TOV spacetime.
    body = grhayl_setup_str + mode_setup + f"""
  ID_persist_struct ID_persist;
  TOV_read_data_file_set_ID_persist("{tov_data_filename_escaped}", &ID_persist);

  for(int grid=0; grid<commondata->NUMGRIDS; grid++) {{
    // Unpack the current grid and its gridfunctions.
    params_struct *restrict params = &griddata[grid].params;
    bc_struct *restrict bcstruct = &griddata[grid].bcstruct;
#include "set_CodeParameters.h"
    REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++) {{
      xx[ww] = griddata[grid].xx[ww];
    }} // END LOOP: for ww over coordinate directions
    REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
    REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;

    // Import the TOV spacetime data into the BSSN gridfunctions.
    initial_data_reader__convert_ADM_Spherical_to_BSSN(
        commondata, params, xx, bcstruct, &griddata[grid].gridfuncs,
        &ID_persist, TOV_ID_function);
"""

    # Step 7: Build the pointwise matter-seeding code for the requested mode.
    if tov_mode == "hybrid":
        loop_body = common_loop_prefix + ccg.c_codegen(
            hybrid_output_vars,
            hybrid_grid_access,
            enable_simd=False,
            enable_GoldenKernels=enable_GoldenKernels,
        )
        loop_body += _build_apply_grhayl_limits_code(declare_grhayl_structs=True)
    elif tov_mode == "tabulated_temperature":
        loop_body = common_loop_prefix + r"""
      // Seed the tabulated primitive variables at fixed temperature.
      ghl_primitive_quantities prims;
      ghl_metric_quantities metric;
      ghl_conservative_quantities cons;
      double dummy0, dummy1, dummy2, dummy3, dummy4;

      prims.rho = rho_baryon;
      prims.press = press;
      prims.vU[0] = 0.0;
      prims.vU[1] = 0.0;
      prims.vU[2] = 0.0;

      const double rhoL = prims.rho;
      if (rhoL <= 1.01 * commondata->eos.rho_atm) {
        prims.rho = commondata->eos.rho_atm;
        prims.press = commondata->eos.press_atm;
        prims.eps = commondata->eos.eps_atm;
        prims.Y_e = commondata->eos.Y_e_atm;
        prims.temperature = commondata->eos.T_atm;
      } else {
        const double tempL = beq_temperature;
        const double YeL =
            ghl_tabulated_compute_Ye_from_rho(&commondata->eos, rhoL);

        double pressL, epsL;
        ghl_tabulated_compute_P_eps_from_T(
            &commondata->eos, rhoL, YeL, tempL, &pressL, &epsL);

        prims.press = pressL;
        prims.eps = epsL;
        prims.Y_e = YeL;
        prims.temperature = tempL;
      } // END IF: atmosphere reset versus tabulated-temperature seed

      ghl_return_primitives(
          &prims,
          &auxevol_gfs[IDX4(RHOBGF, i0, i1, i2)],
          &auxevol_gfs[IDX4(PGF, i0, i1, i2)],
          &dummy0,
          &auxevol_gfs[IDX4(RESCALEDVU0GF, i0, i1, i2)],
          &auxevol_gfs[IDX4(RESCALEDVU1GF, i0, i1, i2)],
          &auxevol_gfs[IDX4(RESCALEDVU2GF, i0, i1, i2)],
          &dummy1,
          &dummy2,
          &dummy3,
          &dummy4,
          &auxevol_gfs[IDX4(YEGF, i0, i1, i2)],
          &auxevol_gfs[IDX4(TEMPERATUREGF, i0, i1, i2)]);
"""
        loop_body += _build_apply_grhayl_limits_code(declare_grhayl_structs=False)
    else:
        loop_body = common_loop_prefix + r"""
      // Seed the tabulated primitive variables at fixed entropy.
      ghl_primitive_quantities prims;
      ghl_metric_quantities metric;
      ghl_conservative_quantities cons;
      double dummy0, dummy1, dummy2, dummy3, dummy4;

      prims.rho = rho_baryon;
      prims.press = press;
      prims.vU[0] = 0.0;
      prims.vU[1] = 0.0;
      prims.vU[2] = 0.0;

      const double rhoL = prims.rho;
      if (rhoL <= 1.01 * commondata->eos.rho_atm) {
        prims.rho = commondata->eos.rho_atm;
        prims.press = commondata->eos.press_atm;
        prims.eps = commondata->eos.eps_atm;
        prims.Y_e = commondata->eos.Y_e_atm;
        prims.temperature = commondata->eos.T_atm;
      } else {
        prims.Y_e = ghl_tabulated_compute_Ye_from_rho(&commondata->eos, rhoL);

        // Supply a temperature guess before inverting for the fixed-entropy state.
        prims.temperature = 1.0;
        ghl_tabulated_compute_P_T_from_S(
            &commondata->eos,
            prims.rho,
            prims.Y_e,
            beq_entropy,
            &prims.press,
            &prims.temperature);
      } // END IF: atmosphere reset versus tabulated-entropy seed

      ghl_return_primitives(
          &prims,
          &auxevol_gfs[IDX4(RHOBGF, i0, i1, i2)],
          &auxevol_gfs[IDX4(PGF, i0, i1, i2)],
          &dummy0,
          &auxevol_gfs[IDX4(RESCALEDVU0GF, i0, i1, i2)],
          &auxevol_gfs[IDX4(RESCALEDVU1GF, i0, i1, i2)],
          &auxevol_gfs[IDX4(RESCALEDVU2GF, i0, i1, i2)],
          &dummy1,
          &dummy2,
          &dummy3,
          &dummy4,
          &auxevol_gfs[IDX4(YEGF, i0, i1, i2)],
          &auxevol_gfs[IDX4(TEMPERATUREGF, i0, i1, i2)]);
"""
        loop_body += _build_apply_grhayl_limits_code(declare_grhayl_structs=False)

    body += lp.simple_loop(
        loop_body=loop_body,
        loop_region="all points",
        enable_intrinsics=False,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=False,
        read_xxs=True,
        OMP_collapse=OMP_collapse,
    )

    # Step 8: Finalize cleanup for TOV profile and optional EOS-slice data.
    body += r"""
  } // END LOOP: for grid over grids

  // Release the interpolated TOV profile storage after initial-data setup.
  free(ID_persist.r_Schw_arr);
  free(ID_persist.rho_arr);
  free(ID_persist.rho_baryon_arr);
  free(ID_persist.P_arr);
  free(ID_persist.M_arr);
  free(ID_persist.expnu_arr);
  free(ID_persist.exp4phi_arr);
  free(ID_persist.rbar_arr);
"""
    if tov_mode == "tabulated_constant_entropy":
        body += r"""
  free(commondata->eos.Ye_of_lr);
  commondata->eos.Ye_of_lr = NULL;
  free(commondata->eos.lp_of_lr);
  commondata->eos.lp_of_lr = NULL;
  free(commondata->eos.le_of_lr);
  commondata->eos.le_of_lr = NULL;
"""

    # Step 9: Register the generated C function.
    cfc.register_CFunction(
        include_CodeParameters_h=False,
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
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
