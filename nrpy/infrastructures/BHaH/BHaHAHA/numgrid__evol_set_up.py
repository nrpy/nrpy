"""
Register function for setting up the evolved numerical grid for h(theta, phi).

Set up evolved numerical grids for use with BHaHAHA.
Original code stolen from nrpy/infrastructures/BHaH/numerical_grids_and_timestep.py

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List, Tuple

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par
from nrpy.infrastructures import BHaH


def register_CFunction_numgrid__evol_set_up() -> None:
    """Register a C function to set up evolved numerical grid for use with BHaHAHA."""
    # Step 1: Register gridfunctions on this grid.
    # Register EVOL gridfunctions
    _, __ = gri.register_gridfunctions(["hh", "vv"], group="EVOL")

    # Register AUXEVOL gridfunctions
    _ = gri.register_gridfunctions_for_single_rank2(
        "hDD", symmetry="sym01", group="AUXEVOL", gf_array_name="auxevol_gfs"
    )
    _ = gri.register_gridfunctions("WW", group="AUXEVOL", gf_array_name="auxevol_gfs")[
        0
    ]
    _ = gri.register_gridfunctions_for_single_rank1(
        "partial_D_WW", group="AUXEVOL", gf_array_name="auxevol_gfs"
    )
    _ = gri.register_gridfunctions_for_single_rankN(
        "partial_D_hDD",
        rank=3,
        symmetry="sym12",
        group="AUXEVOL",
        gf_array_name="auxevol_gfs",
    )
    _ = gri.register_gridfunctions("trK", group="AUXEVOL", gf_array_name="auxevol_gfs")[
        0
    ]
    _ = gri.register_gridfunctions_for_single_rank2(
        "aDD", symmetry="sym01", group="AUXEVOL", gf_array_name="auxevol_gfs"
    )

    # Step 2: Register (curvilinear boundary condition) contributions to BHaH_defines.h and griddata.
    list_of_evol_gf_names_ranks: List[Tuple[str, int]] = []
    # Scalar gridfunctions
    list_of_evol_gf_names_ranks.append(("hh", 0))
    list_of_evol_gf_names_ranks.append(("vv", 0))
    BHaH_defines_contrib = BHaH.BHaHAHA.bcstruct_set_up.BHaH_defines_set_gridfunction_defines_with_parity_types(
        grid_name="evol",
        list_of_gf_names_ranks=list_of_evol_gf_names_ranks,
        verbose=True,
    )
    BHaH.BHaH_defines_h.register_BHaH_defines(__name__, BHaH_defines_contrib)

    # Carrying over previously found horizon to the current horizon evolution.
    BHaH.griddata_commondata.register_griddata_commondata(
        __name__,
        "bhahaha_params_and_data_struct *restrict bhahaha_params_and_data",
        "input parameters and data set by the external code",
        is_commondata=True,
    )
    BHaH.griddata_commondata.register_griddata_commondata(
        __name__,
        "REAL *restrict coarse_horizon",
        "most recently found coarse-resolution horizon h(theta, phi) on the evolution grid",
        is_commondata=True,
    )
    BHaH.griddata_commondata.register_griddata_commondata(
        __name__,
        "int use_coarse_horizon",
        "prolongate h from coarse to finer, in initial_data(); 1=yes, 0=no",
        is_commondata=True,
    )
    for i in range(1, 3):
        BHaH.griddata_commondata.register_griddata_commondata(
            __name__,
            f"REAL coarse_horizon_dxx{i}",
            f"coarse horizon dxx{i}",
            is_commondata=True,
        )
        BHaH.griddata_commondata.register_griddata_commondata(
            __name__,
            f"int coarse_horizon_Nxx_plus_2NGHOSTS{i}",
            f"coarse horizon Nxx_plus_2NGHOSTS{i}",
            is_commondata=True,
        )
    BHaH.griddata_commondata.register_griddata_commondata(
        __name__,
        "REAL *restrict coarse_horizon_r_theta_phi[3]",
        "coarse horizon r_theta_phi",
        is_commondata=True,
    )
    # By NRPy/BHaH convention, bcstruct goes in griddata for the evolved grid.
    BHaH.griddata_commondata.register_griddata_commondata(
        __name__,
        "bc_struct bcstruct",
        "all data needed to apply boundary conditions in curvilinear coordinates",
    )

    # Step 3: Register CodeParameters specific to this grid.
    # fmt: off
    for i in range(3):
        _ = par.CodeParameter("int", __name__, f"Nxx{i}", 64)
    for i in range(3):
        _ = par.CodeParameter("int", __name__, f"Nxx_plus_2NGHOSTS{i}", add_to_parfile=False,
                              add_to_set_CodeParameters_h=True)
        _ = par.CodeParameter("REAL", __name__, f"dxx{i}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
        _ = par.CodeParameter("REAL", __name__, f"invdxx{i}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
    # fmt: on

    # Step 4: Register numgrid__evol_set_up().
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Set up numerical grid for BHaHAHA 2D evolution grids; Nxx0 = Nr = 1."
    cfunc_type = "void"
    name = "numgrid__evol_set_up"
    params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata, const int Nx_evol_grid[3]"
    body = r"""
  // Step 1: Set default parameters in griddata.params
  bah_params_struct_set_to_default(commondata, griddata);

  // Step 2: Set number of grids to 1 for BHaHAHA
  commondata->NUMGRIDS = 1;

  // Step 3: Initialize grid parameters
  const int grid = 0;
  params_struct *restrict params = &griddata[grid].params;

  snprintf(params->CoordSystemName, 100, "Spherical"); // Must be set, or valgrind will complain about reading this in set_CodeParameters.h
  params->grid_physical_size = 1.0; // Unused, since h sets the actual radius

  // Set grid sizes from Nx_evol_grid
  params->Nxx0 = Nx_evol_grid[0];
  params->Nxx1 = Nx_evol_grid[1];
  params->Nxx2 = Nx_evol_grid[2];

  // Include ghost zones
  params->Nxx_plus_2NGHOSTS0 = params->Nxx0 + 2 * NGHOSTS;
  params->Nxx_plus_2NGHOSTS1 = params->Nxx1 + 2 * NGHOSTS;
  params->Nxx_plus_2NGHOSTS2 = params->Nxx2 + 2 * NGHOSTS;

  // Step 4: Set grid boundaries
  params->RMAX = 1.0; // Completely arbitrary; must be set to some reasonable-sized number > 0, just so that BC set up doesn't error out.
  const REAL xxmin0 = 0.0;
  const REAL xxmax0 = params->RMAX;
  const REAL xxmin1 = 0.0;
  const REAL xxmax1 = M_PI;
  const REAL xxmin2 = -M_PI;
  const REAL xxmax2 = M_PI;

  // Step 5: Compute dxx and invdxx
  params->dxx0 = (xxmax0 - xxmin0) / ((REAL)params->Nxx0);
  params->dxx1 = (xxmax1 - xxmin1) / ((REAL)params->Nxx1);
  params->dxx2 = (xxmax2 - xxmin2) / ((REAL)params->Nxx2);

  params->invdxx0 = 1.0 / params->dxx0;
  params->invdxx1 = 1.0 / params->dxx1;
  params->invdxx2 = 1.0 / params->dxx2;

  // Arrays for simplifying loops
  const int Nxx_plus_2NGHOSTS[3] = {params->Nxx_plus_2NGHOSTS0, params->Nxx_plus_2NGHOSTS1, params->Nxx_plus_2NGHOSTS2};
  const REAL xxmin[3] = {xxmin0, xxmin1, xxmin2};
  const REAL dxx[3] = {params->dxx0, params->dxx1, params->dxx2};

  // Step 6: Allocate and initialize cell-centered grid coordinate arrays xx[0], xx[1], xx[2]
  for (int dir = 0; dir < 3; dir++) {
    griddata[grid].xx[dir] = (REAL *restrict)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS[dir]);
    for (int i = 0; i < Nxx_plus_2NGHOSTS[dir]; i++)
      griddata[grid].xx[dir][i] = xxmin[dir] + ((REAL)(i - NGHOSTS) + (1.0 / 2.0)) * dxx[dir];
  }

  // Step 7: Allocate and define reference-metric precompute lookup arrays.
  griddata[grid].rfmstruct = (rfm_struct *)malloc(sizeof(rfm_struct));
  bah_rfm_precompute_malloc(commondata, params, griddata[grid].rfmstruct);
  bah_rfm_precompute_defines(commondata, params, griddata[grid].rfmstruct, griddata[grid].xx);

  // Step 8: Set up bcstruct, for setting inner boundary conditions (i.e., BCs in theta & phi ghost zones, like theta < 0).
  {
    commondata->bcstruct_dxx0 = params->dxx0;
    commondata->bcstruct_dxx1 = params->dxx1;
    commondata->bcstruct_dxx2 = params->dxx2;
    commondata->bcstruct_Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
    commondata->bcstruct_Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
    commondata->bcstruct_Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
    bah_bcstruct_set_up(commondata, griddata[grid].xx, &griddata[grid].bcstruct);
  }

  // Step 9: Initialize time-stepping parameters
  commondata->nn = 0;
  commondata->nn_0 = 0;
  commondata->t_0 = 0.0;
  commondata->time = 0.0;
"""
    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
