"""
Register numerical_grids_and_timestep() C function, as well as functions called by this one.

These functions set up numerical grids for use within the BHaH infrastructure using
OpenMP parallelization

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Samuel D. Tootle
        sdtootle **at** gmail **dot* com        
"""

from typing import Dict, List
import nrpy.c_function as cfc
import nrpy.params as par
import nrpy.infrastructures.BHaH.grid_management.base_numerical_grids_and_timestep as base_gm_classes

# fmt: off
for i in range(3):
    _ = par.CodeParameter("int", __name__, f"Nxx_plus_2NGHOSTS{i}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("int", __name__, f"Nxx{i}", 64)
    # reference_metric sets xxmin and xxmax below.
    _ = par.CodeParameter("REAL", __name__, f"xxmin{i}", -10.0, add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("REAL", __name__, f"xxmax{i}", 10.0, add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("REAL", __name__, f"invdxx{i}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("REAL", __name__, f"dxx{i}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
_ = par.CodeParameter("REAL", __name__, "convergence_factor", 1.0, commondata=True)
_ = par.CodeParameter("int", __name__, "CoordSystem_hash", commondata=False, add_to_parfile=False)
_ = par.CodeParameter("int", __name__, "grid_idx", commondata=False, add_to_parfile=False)
_ = par.CodeParameter("char[200]", __name__, "gridding_choice", "independent grid(s)", commondata=True, add_to_parfile=True)
# fmt: on


class register_CFunction_numerical_grid_params_Nxx_dxx_xx(
    base_gm_classes.base_register_CFunction_numerical_grid_params_Nxx_dxx_xx
):
    """
    Register a C function to Set up a cell-centered grid of size grid_physical_size.
    Set params: Nxx, Nxx_plus_2NGHOSTS, dxx, invdxx, and xx.

    :param CoordSystem: The coordinate system used for the simulation.
    :param grid_physical_size: The physical size of the grid.
    :param Nxx_dict: A dictionary that maps coordinate systems to lists containing the number of grid points along each direction.

    :return: None.
    :raises ValueError: If CoordSystem is not in Nxx_dict.
    """

    def __init__(
        self,
        CoordSystem: str,
        grid_physical_size: float,
        Nxx_dict: Dict[str, List[int]],
    ) -> None:
        super().__init__(CoordSystem, grid_physical_size, Nxx_dict)

        self.body+="""
    // Set up cell-centered Cartesian coordinate grid, centered at the origin.
    xx[0] = (REAL *restrict)malloc(sizeof(REAL)*params->Nxx_plus_2NGHOSTS0);
    xx[1] = (REAL *restrict)malloc(sizeof(REAL)*params->Nxx_plus_2NGHOSTS1);
    xx[2] = (REAL *restrict)malloc(sizeof(REAL)*params->Nxx_plus_2NGHOSTS2);
    for (int j = 0; j < params->Nxx_plus_2NGHOSTS0; j++) xx[0][j] = params->xxmin0 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * params->dxx0;
    for (int j = 0; j < params->Nxx_plus_2NGHOSTS1; j++) xx[1][j] = params->xxmin1 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * params->dxx1;
    for (int j = 0; j < params->Nxx_plus_2NGHOSTS2; j++) xx[2][j] = params->xxmin2 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * params->dxx2;
    """
        cfc.register_CFunction(
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            CoordSystem_for_wrapper_func=CoordSystem,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=True,  # keep this False or regret having to debug the mess.
            body=self.body,
        )


class register_CFunction_cfl_limited_timestep(
    base_gm_classes.base_register_CFunction_cfl_limited_timestep
):
    """
    Register a C function to find the CFL-limited timestep dt on a numerical grid.

    The timestep is determined by the relation dt = CFL_FACTOR * ds_min, where ds_min
    is the minimum spacing between neighboring gridpoints on a numerical grid.

    :param CoordSystem: The coordinate system used for the simulation.
    :param fp_type: Floating point type, e.g., "double".
    :return: None.
    """

    def __init__(self, CoordSystem: str, fp_type: str = "double") -> None:
        super().__init__(CoordSystem, fp_type=fp_type)
        # could be replaced by simple loop?
        self.body = r"""
REAL ds_min = 1e38;
#pragma omp parallel for reduction(min:ds_min)
LOOP_NOOMP(i0, 0, Nxx_plus_2NGHOSTS0,
           i1, 0, Nxx_plus_2NGHOSTS1,
           i2, 0, Nxx_plus_2NGHOSTS2) {
    const REAL xx0 = xx[0][i0];
    const REAL xx1 = xx[1][i1];
    const REAL xx2 = xx[2][i2];
    REAL dsmin0, dsmin1, dsmin2;
"""
        self.body += self.min_body

        cfc.register_CFunction(
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            CoordSystem_for_wrapper_func=self.CoordSystem,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=True,
            body=self.body,
        )


class register_CFunction_numerical_grids_and_timestep(
    base_gm_classes.base_register_CFunction_numerical_grids_and_timestep
):
    """
    Register a C function to set up all numerical grids and timestep.

    The function configures the numerical grids based on given parameters, specifically
    focusing on the usage of reference metric precomputations and curvilinear boundary
    conditions.

    :param list_of_CoordSystems: List of CoordSystems
    :param enable_rfm_precompute: Whether to enable reference metric precomputation (default: False).
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions (default: False).
    :return: None.
    """

    def __init__(
        self,
        list_of_CoordSystems: List[str],
        enable_rfm_precompute: bool = False,
        enable_CurviBCs: bool = False,
    ) -> None:
        super().__init__(
            list_of_CoordSystems,
            enable_rfm_precompute=enable_rfm_precompute,
            enable_CurviBCs=enable_CurviBCs,
        )
        self.body = r"""
    // Step 1.a: Set each CodeParameter in griddata.params to default, for MAXNUMGRIDS grids.
    params_struct_set_to_default(commondata, griddata);"""
        self.body += r"""
      if(strncmp(commondata->gridding_choice, "independent grid(s)", 200) == 0) {
        // Independent grids
        bool grid_is_resized=false;
        int Nx[3] = { -1, -1, -1 };


        // Step 1.b: For each grid, set Nxx & Nxx_plus_2NGHOSTS, as well as dxx, invdxx, & xx based on grid_physical_size
        int grid=0;
    """
        for CoordSystem in self.list_of_CoordSystems:
            self.body += (
                f"griddata[grid].params.CoordSystem_hash = {CoordSystem.upper()};\n"
            )
            self.body += "numerical_grid_params_Nxx_dxx_xx(commondata, &griddata[grid].params, griddata[grid].xx, Nx, grid_is_resized);\n"
            self.body += "grid++;\n\n"
        self.body += r"""}

// Step 1.c: Allocate memory for and define reference-metric precomputation lookup tables
"""
        if self.enable_rfm_precompute:
            self.body += r"""
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  rfm_precompute_malloc(commondata, &griddata[grid].params, &griddata[grid].rfmstruct);
  rfm_precompute_defines(commondata, &griddata[grid].params, &griddata[grid].rfmstruct, griddata[grid].xx);
}
"""
        else:
            self.body += "// (reference-metric precomputation disabled)\n"
        self.body += (
            "\n// Step 1.d: Set up curvilinear boundary condition struct (bcstruct)\n"
        )

        if self.enable_CurviBCs:
            self.body += r"""
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  bcstruct_set_up(commondata, &griddata[grid].params, griddata[grid].xx, &griddata[grid].bcstruct);
}
"""
        else:
            self.body += "// (curvilinear boundary conditions bcstruct disabled)\n"

        self.body += r"""
// Step 1.e: Set timestep based on minimum spacing between neighboring gridpoints.
commondata->dt = 1e30;
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  cfl_limited_timestep(commondata, &griddata[grid].params, griddata[grid].xx, &griddata[grid].bcstruct);
}

// Step 1.f: Initialize timestepping parameters to zero if this is the first time this function is called.
if(calling_for_first_time) {
  commondata->nn = 0;
  commondata->nn_0 = 0;
  commondata->t_0 = 0.0;
  commondata->time = 0.0;
}
"""
        cfc.register_CFunction(
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=False,
            body=self.body,
        )


def register_CFunctions(
    list_of_CoordSystems: List[str],
    grid_physical_size: float,
    Nxx_dict: Dict[str, List[int]],
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
    fp_type: str = "double",
) -> None:
    """
    Register C functions related to coordinate systems and grid parameters.

    :param list_of_CoordSystems: List of CoordSystems
    :param grid_physical_size: Physical size of the grid.
    :param Nxx_dict: Dictionary containing number of grid points.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions.
    :param fp_type: Floating point type, e.g., "double".
    """
    for CoordSystem in list_of_CoordSystems:
        register_CFunction_numerical_grid_params_Nxx_dxx_xx(
            CoordSystem=CoordSystem,
            grid_physical_size=grid_physical_size,
            Nxx_dict=Nxx_dict,
        )
        register_CFunction_cfl_limited_timestep(
            CoordSystem=CoordSystem, fp_type=fp_type
        )
    register_CFunction_numerical_grids_and_timestep(
        list_of_CoordSystems=list_of_CoordSystems,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_CurviBCs=enable_CurviBCs,
    )
