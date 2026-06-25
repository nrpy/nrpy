"""
Function for TOV initial data, for use by GRoovy.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot** com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri  # NRPy+: Functionality for handling numerical grids
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.infrastructures.BHaH.simple_loop as lp


def register_CFunction_hybrid_EoS_TOV_initial_data(
    grhayl_setup_str: str,
    CoordSystem: str,
    OMP_collapse: int,
    enable_GoldenKernels: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function to initialize grid functions with TOV initial data.

    :param grhayl_setup_str: GRHayL EOS and parameter setup code.
    :param CoordSystem: The coordinate system.
    :param OMP_collapse: Level of OpenMP loop collapsing.
    :param enable_GoldenKernels: Enable Golden Kernels.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Initialize GRoovy TOV initial data from TOVola and GRHayL.

@param[in,out] commondata Global simulation data, including GRHayL parameters and EOS.
@param[in,out] griddata Per-grid coordinates, parameters, and gridfunction storage."""
    cfunc_type = "void"
    # Has to be named this, for generic call in main.c
    name = "initial_data"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    body = grhayl_setup_str
    body += r"""

// TOVola CodeParameters are already populated from defaults and the parfile, and placed within commondata.
ID_persist_struct ID_persist;
TOVola_solve(commondata, &ID_persist);

for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
  bc_struct *restrict bcstruct = &griddata[grid].bcstruct;

  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];
  REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;

  initial_data_reader__convert_ADM_Spherical_to_BSSN(
      commondata, params, (const REAL *restrict *)xx, bcstruct,
      &griddata[grid].gridfuncs, &ID_persist, TOVola_interp);
"""

    loop_body = r"""  // Interpolate through TOVola's public pointwise interface.
  const REAL xOrig[3] = {xx0, xx1, xx2};
  REAL xCart[3];
  xx_to_Cart(params, xOrig, xCart);

  initial_data_struct initial_data;
  TOVola_interp(commondata, params, xCart, &ID_persist, &initial_data);

  // For a static perfect fluid, T^{rr} gamma_{rr} = P. Recover rho_b from the hybrid EOS.
  // In the TOVola_interp function call we set: 
  // initial_data->gammaSphorCartDD00 = exp4phi_val;
  // initial_data->T4SphorCartUU11 = P_val / exp4phi_val;
  // Thus, we can get P_val using:
  const REAL P_val = initial_data.T4SphorCartUU11 * initial_data.gammaSphorCartDD00;
  const REAL rho_baryon_val =
      ghl_hybrid_compute_rho_cold_from_P_cold(&commondata->eos, NRPYMAX(P_val, 0.0));

"""

    rho_baryon_val, P_val = sp.symbols("rho_baryon_val P_val", real=True)

    rescaledvU = ixp.zerorank1()

    output_vars = [
        rho_baryon_val,
        P_val,
        rescaledvU[0],
        rescaledvU[1],
        rescaledvU[2],
    ]

    vars_grid_access = [
        gri.BHaHGridFunction.access_gf("rhob", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("P", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledvU0", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledvU1", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledvU2", gf_array_name="auxevol_gfs"),
    ]

    apply_grhayl_limits_str = r"""

	ghl_primitive_quantities prims = {0};
	ghl_metric_quantities metric = {0};
	ghl_conservative_quantities cons = {0};

	basis_transform_rfm_basis_to_Cartesian(
        commondata, params, &prims, &cons, &metric,
        i0, i1, i2, xx, auxevol_gfs, in_gfs);

	bool speed_limited = false;
	ghl_error_codes_t error = ghl_enforce_primitive_limits_and_compute_u0(
			&commondata->ghl_params, &commondata->eos, &metric, &prims, &speed_limited);
	if (error != ghl_success) {
	fprintf(stderr,
            "GRHayL primitive limiting failed in TOV initial data at (%d, %d, %d).\n",
            i0, i1, i2);
	exit(EXIT_FAILURE);
	} // END IF: primitive limiting failed for TOV initial data

	basis_transform_Cartesian_to_rfm_basis(commondata, params, &prims, &cons,
																					i0, i1, i2, xx, auxevol_gfs, in_gfs);

	auxevol_gfs[IDX4(U4UTGF, i0, i1, i2)] = prims.u0;
"""

    body += lp.simple_loop(
        loop_body=loop_body
        + ccg.c_codegen(
            output_vars,
            vars_grid_access,
            enable_simd=False,
            enable_GoldenKernels=enable_GoldenKernels,
        )
        + apply_grhayl_limits_str,
        loop_region="all points",
        enable_intrinsics=False,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=False,
        read_xxs=True,
        OMP_collapse=OMP_collapse,
    )

    body += r"""
	} // END LOOP: for grid over all numerical grids

	free(ID_persist.r_Schw_arr);
	free(ID_persist.rho_energy_arr);
	free(ID_persist.rho_baryon_arr);
	free(ID_persist.P_arr);
	free(ID_persist.M_arr);
	free(ID_persist.expnu_arr);
	free(ID_persist.exp4phi_arr);
	free(ID_persist.r_iso_arr);
   """
    cfc.register_CFunction(
        include_CodeParameters_h=False,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )
    return pcg.NRPyEnv()
