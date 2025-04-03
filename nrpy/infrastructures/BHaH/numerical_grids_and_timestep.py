"""
Register numerical_grids_and_timestep() C function, as well as functions called by this one.

These functions set up numerical grids for use within the BHaH infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, List, Set, Tuple, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)

# fmt: off
for j in range(3):
    _ = par.CodeParameter("int", __name__, f"Nxx_plus_2NGHOSTS{j}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("int", __name__, f"Nxx{j}", 64)
    # reference_metric sets xxmin and xxmax below.
    _ = par.CodeParameter("REAL", __name__, f"xxmin{j}", -10.0, add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("REAL", __name__, f"xxmax{j}", 10.0, add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("REAL", __name__, f"invdxx{j}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("REAL", __name__, f"dxx{j}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
_ = par.CodeParameter("REAL", __name__, "convergence_factor", 1.0, commondata=True)
_ = par.CodeParameter("int", __name__, "CoordSystem_hash", commondata=False, add_to_parfile=False)
_ = par.CodeParameter("int", __name__, "grid_idx", 0, commondata=False, add_to_parfile=False)
_ = par.CodeParameter("char[100]", __name__, "gridname", commondata=False, add_to_parfile=False)
# fmt: on


def register_CFunction_numerical_grid_params_Nxx_dxx_xx(
    CoordSystem: str, Nxx_dict: Dict[str, List[int]]
) -> None:
    """
    Register a C function to set up a cell-centered grid of size grid_physical_size.
       Set params: Nxx, Nxx_plus_2NGHOSTS, dxx, invdxx, and xx.

    :param CoordSystem: The coordinate system used for the simulation.
    :param Nxx_dict: A dictionary that maps coordinate systems to lists containing the number of grid points along each direction.

    :raises ValueError: If CoordSystem is not in Nxx_dict.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import unittest_CoordSystems
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name = "numerical_grid_params_Nxx_dxx_xx"
    >>> Nxx_dict = { k : [72, 12, 8] for k in unittest_CoordSystems }
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in unittest_CoordSystems:
    ...       cfc.CFunction_dict.clear()
    ...       register_CFunction_numerical_grid_params_Nxx_dxx_xx(CoordSystem, Nxx_dict)
    ...       generated_str = cfc.CFunction_dict[f'{name}__rfm__{CoordSystem}'].full_function
    ...       validation_desc = f"{name}__{parallelization}__{CoordSystem}"
    ...       validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    if CoordSystem not in Nxx_dict:
        raise ValueError(
            f"{CoordSystem} is not in Nxx_dict = {Nxx_dict}. Please add it."
        )
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"""Initializes a cell-centered grid in {CoordSystem} coordinates based on physical dimensions (grid_physical_size).

Inputs:
- Nx[] inputs: Specifies new grid dimensions, if needed.
- params.convergence_factor (set to 1.0 by default): Factor by which grid resolution is increased; set to 1.0 by default.
- set_xxmin_xxmax_to_defaults: Whether to set xxmin[3], xxmax[3] to default values set in reference_metric.py.

Parameter outputs:
- Nxx: Number of grid points in each direction.
- Nxx_plus_2NGHOSTS: Total grid points including ghost zones.
- dxx: Grid spacing.
- invdxx: Inverse of grid spacing.

Grid setup output:
- xx: Coordinate values for each (cell-centered) grid point.
"""
    parallelization = par.parval_from_str("parallelization")
    cfunc_type = "void"
    name = "numerical_grid_params_Nxx_dxx_xx"
    params = "const commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3], const int Nx[3], const bool set_xxmin_xxmax_to_defaults".replace(
        "REAL *restrict xx[3]",
        "REAL * xx[3]" if parallelization in ["cuda"] else "REAL *restrict xx[3]",
    )
    body = "// Set default values for the grid resolution in each dimension.\n"
    for dirn in range(3):
        body += f"params->Nxx{dirn} = {Nxx_dict[CoordSystem][dirn]};\n"
    body += """
// If all components of Nx[] are set to valid values (i.e., not -1), override the default values with Nx[].
if( Nx[0]!=-1 && Nx[1]!=-1 && Nx[2]!=-1 ) {
"""
    for dirn in range(3):
        body += f"params->Nxx{dirn} = Nx[{dirn}];\n"
    body += f"""}}
snprintf(params->CoordSystemName, 100, "{CoordSystem}");

// Resize grid by convergence_factor; used for convergence testing.
{{
  // convergence_factor does not increase resolution across an axis of symmetry (Nxx == 2):
  if(params->Nxx0 != 2) params->Nxx0 *= commondata->convergence_factor;
  if(params->Nxx1 != 2) params->Nxx1 *= commondata->convergence_factor;
  if(params->Nxx2 != 2) params->Nxx2 *= commondata->convergence_factor;
}}

// Set the full grid size; including the ghostzones (of width NGHOSTS) on the boundaries.
params->Nxx_plus_2NGHOSTS0 = params->Nxx0 + 2*NGHOSTS;
params->Nxx_plus_2NGHOSTS1 = params->Nxx1 + 2*NGHOSTS;
params->Nxx_plus_2NGHOSTS2 = params->Nxx2 + 2*NGHOSTS;

"""
    rfm = refmetric.reference_metric[CoordSystem]
    # Set grid_physical_size & grid_hole_radius
    body += """{
#include "../set_CodeParameters.h"
// Set grid size to a function of grid_physical_size (grid_physical_size set in set_CodeParameters.h above):
"""
    for key, value in rfm.grid_physical_size_dict.items():
        body += f"params->{key} = {value};\n"
    body += "}\n"

    # Set grid_hole_radius
    if "Holey" in CoordSystem or "Wedge" in CoordSystem:
        body += """{
#include "../set_CodeParameters.h"
// Set grid hole radius to a function of grid_hole_radius (grid_hole_radius set in set_CodeParameters.h above):
"""
        for key, value in rfm.grid_hole_radius_dict.items():
            body += f"params->{key} = {value};\n"
        body += "}\n"

    # Set minimum and maximum values of xx[][] for each grid.
    body += """if (set_xxmin_xxmax_to_defaults) {
#include "../set_CodeParameters.h"
// Set {xxmin[], xxmax[]} to default values, which could be functions of other rfm params (set in set_CodeParameters.h above):
"""
    for minmax in ["min", "max"]:
        for dirn in range(3):
            rfm_value = rfm.xxmin[dirn] if minmax == "min" else rfm.xxmax[dirn]
            body += f"params->xx{minmax}{dirn} = {rfm_value};\n"
    body += "}\n"

    # Set quantities that depend on Nxx and {xxmin, xxmax}, then set up coordinate arrays xx[3][Nxxi].
    body += """
// Set quantities that depend on Nxx and {xxmin, xxmax}: dxx, invdxx.
params->dxx0 = (params->xxmax0 - params->xxmin0) / ((REAL)params->Nxx0);
params->dxx1 = (params->xxmax1 - params->xxmin1) / ((REAL)params->Nxx1);
params->dxx2 = (params->xxmax2 - params->xxmin2) / ((REAL)params->Nxx2);

params->invdxx0 = ((REAL)params->Nxx0) / (params->xxmax0 - params->xxmin0);
params->invdxx1 = ((REAL)params->Nxx1) / (params->xxmax1 - params->xxmin1);
params->invdxx2 = ((REAL)params->Nxx2) / (params->xxmax2 - params->xxmin2);

// Set up uniform, cell-centered, topologically Cartesian numerical grid,
//   centered at (xxmin[i] + xxmax[i])/2 in direction i, and store
//   {xx[0], xx[1], xx[2]} arrays.
BHAH_MALLOC(xx[0], sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
BHAH_MALLOC(xx[1], sizeof(REAL) * params->Nxx_plus_2NGHOSTS1);
BHAH_MALLOC(xx[2], sizeof(REAL) * params->Nxx_plus_2NGHOSTS2);
""".replace(
        "BHAH_MALLOC",
        "BHAH_MALLOC_DEVICE" if parallelization not in ["openmp"] else "BHAH_MALLOC",
    )

    if parallelization in ["cuda"]:
        body += (
            "cpyHosttoDevice_params__constant(params, params->grid_idx % NUM_STREAMS);"
        )

    prefunc = (
        """
        #define SET_XX_CELL_CENTERED_COORDS(COORD_DIR) \
        const int index  = blockIdx.x * blockDim.x + threadIdx.x; \
        const int stride = blockDim.x * gridDim.x; \
        static constexpr REAL onehalf = 1.0 / 2.0; \
        for (int j = index; j < d_params[streamid].Nxx_plus_2NGHOSTS##COORD_DIR; j+=stride) { \
            xx##COORD_DIR[j] = d_params[streamid].xxmin##COORD_DIR + ((REAL)(j - NGHOSTS) + onehalf) * d_params[streamid].dxx##COORD_DIR;\
        }
        """
        if parallelization == "cuda"
        else """
        #define SET_XX_CELL_CENTERED_COORDS(COORD_DIR) \
        static const REAL onehalf = 1.0 / 2.0; \
        for (int j = 0; j < params->Nxx_plus_2NGHOSTS##COORD_DIR; j+=1) { \
            xx##COORD_DIR[j] = params->xxmin##COORD_DIR + ((REAL)(j - NGHOSTS) + onehalf) * params->dxx##COORD_DIR;\
        }
        """
    )
    for i in range(3):
        kernel_body = f"SET_XX_CELL_CENTERED_COORDS({i});"
        # Kernel name
        kernel_name: str = f"initialize_grid_xx{i}"
        # Comments
        comments = f"Kernel to compute xx{i} coordinates."

        # Prepare the argument dicts
        arg_dict_cuda = {
            f"xx{i}": "REAL *restrict",
        }
        arg_dict_host = {
            "params": "const params_struct *restrict",
            **arg_dict_cuda,
        }

        # We replicate the original approach of computing param_streamid in the body for CUDA
        # so define_param_streamid=True, which forces the helper to insert that line.
        new_prefunc, new_body = parallel_utils.generate_kernel_and_launch_code(
            kernel_name,
            kernel_body,
            arg_dict_cuda,
            arg_dict_host,
            par.parval_from_str("parallelization"),
            cfunc_type="static void",
            comments=comments,
            launchblock_with_braces=True,
            launch_dict={
                "blocks_per_grid": [
                    f"(params->Nxx_plus_2NGHOSTS{i} + threads_in_x_dir - 1) / threads_in_x_dir"
                ],
                "stream": "",
            },
        )
        prefunc += new_prefunc
        body += new_body.replace(f"xx{i})", f"xx[{i}])")

    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,  # keep this False or regret having to debug the mess.
        body=body,
    )


def register_CFunction_ds_min_radial_like_dirns_single_pt(
    CoordSystem: str,
) -> None:
    """
    Register a C function to find the minimum grid spacing ds_min_radial_like_dirns, the ds_min in radial-like directions only.

    ds_min_radial_like_dirns is the minimum spacing between neighboring gridpoints on a numerical grid.

    :param CoordSystem: The coordinate system of the numerical grid.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Examining only radial-like (non-angular) directions at a given point on a numerical grid, find the minimum grid spacing ds_min."
    cfunc_type = "void"
    name = "ds_min_radial_like_dirns_single_pt"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xx0, const REAL xx1, const REAL xx2, REAL *restrict ds_min_radial_like_dirns"
    rfm = refmetric.reference_metric[CoordSystem]
    # These are set in CodeParameters.h
    dxx = sp.symbols("dxx0 dxx1 dxx2", real=True)
    body = "MAYBE_UNUSED REAL ds0=1e38, ds1=1e38, ds2=1e38;\n"
    ds_expr_list: List[sp.Expr] = []
    ds_str_list: List[str] = []
    for dirn in rfm.radial_like_dirns:
        ds_expr_list += [sp.Abs(rfm.scalefactor_orthog[dirn] * dxx[dirn])]
        ds_str_list += [f"ds{dirn}"]
    body += ccg.c_codegen(
        ds_expr_list,
        ds_str_list,
        include_braces=False,
    )
    body += "*ds_min_radial_like_dirns = MIN(ds0, MIN(ds1, ds2));\n"
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )


def register_CFunction_ds_min_single_pt(
    CoordSystem: str,
) -> None:
    """
    Register a C function to find the minimum grid spacing ds_min.

    ds_min is the minimum spacing between neighboring gridpoints on a numerical grid.

    :param CoordSystem: The coordinate system of the numerical grid.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import unittest_CoordSystems
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name = "ds_min_single_pt"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in unittest_CoordSystems:
    ...       cfc.CFunction_dict.clear()
    ...       register_CFunction_ds_min_single_pt(CoordSystem)
    ...       generated_str = cfc.CFunction_dict[f'{name}__rfm__{CoordSystem}'].full_function
    ...       validation_desc = f"{name}__{parallelization}__{CoordSystem}"
    ...       validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    Setting up reference_metric[SinhSymTP]...
    Setting up reference_metric[HoleySinhSpherical]...
    Setting up reference_metric[Cartesian]...
    Setting up reference_metric[SinhCylindricalv2n2]...
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Examining all three directions at a given point on a numerical grid, find the minimum grid spacing ds_min."
    cfunc_type = "void"
    name = "ds_min_single_pt"
    params = "const params_struct *restrict params, const REAL xx0, const REAL xx1, const REAL xx2, REAL *restrict ds_min"
    rfm = refmetric.reference_metric[CoordSystem]
    # These are set in CodeParameters.h
    dxx0, dxx1, dxx2 = sp.symbols("dxx0 dxx1 dxx2", real=True)
    expr_list = [
        sp.Abs(rfm.scalefactor_orthog[0] * dxx0),
        sp.Abs(rfm.scalefactor_orthog[1] * dxx1),
        sp.Abs(rfm.scalefactor_orthog[2] * dxx2),
    ]
    body = ccg.c_codegen(
        expr_list,
        ["const REAL ds0", "const REAL ds1", "const REAL ds2"],
        include_braces=False,
    )
    body += "*ds_min = MIN(ds0, MIN(ds1, ds2));\n"

    parallelization = par.parval_from_str("parallelization")
    param_symbols, _ = get_params_commondata_symbols_from_expr_list(expr_list)
    params_definitions = generate_definition_header(
        param_symbols,
        enable_intrinsics=False,
        var_access=parallel_utils.get_params_access("openmp"),
    )

    kernel_body = f"{params_definitions}\n{body}"
    cfunc_decorators = "__host__ __device__" if parallelization == "cuda" else ""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=kernel_body,
        cfunc_decorators=cfunc_decorators,
    )


def generate_grid_minimum_gridspacing_prefunc() -> Tuple[str, str]:
    """
    Compute global minimum of grid spacing.

    The timestep is determined by the relation dt = CFL_FACTOR * ds_min, where ds_min
    is the minimum spacing between neighboring gridpoints on a numerical grid.

    :return: Tuple of prefunc and function_launch code strings.
    """
    comments = "Kernel to find minimum grid spacing."
    cfunc_type = "static void"
    # Kernel name
    name = "compute_ds_min"
    params = "params_struct *restrict params"
    parallelization = par.parval_from_str("parallelization")
    for i in range(3):
        params += f", REAL *restrict xx{i}"
    loop_body = (
        r"""REAL local_ds_min;
    ds_min_single_pt(params, xx0[i0], xx1[i1], xx2[i2], &local_ds_min);
    ds_min = MIN(ds_min, local_ds_min);
    """
        if parallelization not in ["cuda"]
        else r"""ds_min_single_pt(params, xx0[i0], xx1[i1], xx2[i2], &ds_min[IDX3(i0, i1, i2)]);
    """
    )
    body = parallel_utils.get_loop_parameters(parallelization)
    OMP_custom_pragma: str = ""
    if parallelization != "cuda":
        body += "REAL ds_min = 1e38;\n"
        OMP_custom_pragma = "#pragma omp parallel for reduction(min:ds_min)"
    body += rf"""{lp.simple_loop(
        loop_body,
        loop_region="all points",
        OMP_custom_pragma=OMP_custom_pragma)}
"""
    if parallelization != "cuda":
        body += "*ds_min_result = ds_min;\n"

    # Prepare the argument dicts
    arg_dict_cuda = {
        "params": "const params_struct *restrict",
        **{f"xx{i}": "REAL *restrict" for i in range(3)},
        "ds_min": "REAL *restrict",
    }
    arg_dict_host = {
        "params": "const params_struct *restrict",
        **{f"xx{i}": "REAL *restrict" for i in range(3)},
        "ds_min_result": "REAL *restrict",
    }

    prefunc, function_launch = parallel_utils.generate_kernel_and_launch_code(
        name,
        body,
        arg_dict_cuda,
        arg_dict_host,
        parallelization,
        cfunc_type=cfunc_type,
        comments=comments,
        launchblock_with_braces=False,
    )

    for i in range(3):
        function_launch = function_launch.replace(f"xx{i}", f"xx[{i}]")
    # Fix call for device code
    function_launch = function_launch.replace("ds_min)", "ds_min_device)")
    # Fix call for host code
    function_launch = function_launch.replace("ds_min_result)", "ds_min)")
    return prefunc, function_launch


def register_CFunction_cfl_limited_timestep() -> None:
    """
    Register a C function to find the CFL-limited timestep dt on a numerical grid.

    The timestep is determined by the relation dt = CFL_FACTOR * ds_min, where ds_min
    is the minimum spacing between neighboring gridpoints on a numerical grid.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name = "cfl_limited_timestep"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    cfc.CFunction_dict.clear()
    ...    register_CFunction_cfl_limited_timestep()
    ...    generated_str = cfc.CFunction_dict[f'{name}'].full_function
    ...    validation_desc = f"{name}__{parallelization}"
    ...    validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Compute minimum timestep dt = CFL_FACTOR * ds_min."
    cfunc_type = "void"
    name = "cfl_limited_timestep"
    params = "commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3]"
    body = ""
    ds_min_prefunc, ds_min_launch = generate_grid_minimum_gridspacing_prefunc()

    if par.parval_from_str("parallelization") in ["cuda"]:
        body += rf"""
  // Allocate memory for ds_min on the device
  const int Nxx_tot = params->Nxx_plus_2NGHOSTS0 * params->Nxx_plus_2NGHOSTS1 * params->Nxx_plus_2NGHOSTS2;
  REAL *ds_min_device;
  BHAH_MALLOC_DEVICE(ds_min_device, sizeof(REAL) * Nxx_tot);
  params_struct *device_params;
  BHAH_MALLOC_DEVICE(device_params, sizeof(params_struct));
  BHAH_MEMCPY_HOST_TO_DEVICE(device_params, params, sizeof(params_struct));
  {ds_min_launch.replace("params,", "device_params,")}
  REAL ds_min = find_global__minimum(ds_min_device, Nxx_tot);
  """
    else:
        body += rf"""REAL ds_min = 1e38;
        {ds_min_launch.replace(", ds_min", ", &ds_min")}
    """
    body += "commondata->dt = MIN(commondata->dt, ds_min * commondata->CFL_FACTOR);\n"

    if par.parval_from_str("parallelization") in ["cuda"]:
        body += "BHAH_FREE_DEVICE(ds_min_device);\n"
    cfc.register_CFunction(
        prefunc=ds_min_prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


def register_CFunction_numerical_grids_and_timestep(
    set_of_CoordSystems: Set[str],
    list_of_grid_physical_sizes: List[float],
    gridding_approach: str = "independent grid(s)",
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
    enable_set_cfl_timestep: bool = True,
) -> None:
    """
    Register a C function to set up all numerical grids and timestep.

    The function configures the numerical grids based on given parameters, specifically
    focusing on the usage of reference metric precomputations and curvilinear boundary
    conditions.

    :param set_of_CoordSystems: Set of CoordSystems
    :param list_of_grid_physical_sizes: List of grid_physical_size for each CoordSystem; needed for Independent grids.
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :param enable_rfm_precompute: Whether to enable reference metric precomputation (default: False).
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions (default: False).
    :param enable_set_cfl_timestep: Whether to enable computation of dt, the CFL timestep. A custom version can be implemented later.

    :raises ValueError: If invalid gridding_approach selected.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import unittest_CoordSystems
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name = "numerical_grids_and_timestep"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    cfc.CFunction_dict.clear()
    ...    register_CFunction_numerical_grids_and_timestep({"Spherical"}, [10.0], gridding_approach="independent grid(s)", enable_rfm_precompute=True, enable_CurviBCs=True)
    ...    generated_str = cfc.CFunction_dict[f'{name}'].full_function
    ...    validation_desc = f"{name}__{parallelization}"
    ...    validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Set up numerical grids and timestep."
    cfunc_type = "void"
    name = "numerical_grids_and_timestep"
    parallelization = par.parval_from_str("parallelization")
    params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata, bool calling_for_first_time".replace(
        "griddata,",
        (
            "griddata, griddata_struct *restrict griddata_host,"
            if parallelization in ["cuda"]
            else "griddata,"
        ),
    )
    body = r"""
  // Step 1.a: Set each CodeParameter in griddata.params to default, for MAXNUMGRIDS grids.
  if(calling_for_first_time)
    params_struct_set_to_default(commondata, griddata);"""
    body += rf"""
  // Step 1.b: Set commondata->NUMGRIDS to number of CoordSystems we have
  commondata->NUMGRIDS = {len(set_of_CoordSystems)};
"""
    if gridding_approach == "independent grid(s)":
        body += """
  {
    // Independent grids
    int Nx[3] = { -1, -1, -1 };

    // Step 1.c: For each grid, set Nxx & Nxx_plus_2NGHOSTS, as well as dxx, invdxx, & xx based on grid_physical_size
    const bool set_xxmin_xxmax_to_defaults = true;
    int grid=0;
"""
        for which_CoordSystem, CoordSystem in enumerate(set_of_CoordSystems):
            body += f"""// In multipatch, gridname is a helpful alias indicating position of the patch. E.g., "lower {CoordSystem} patch"
    snprintf(griddata[grid].params.gridname, 100, "grid_{CoordSystem}");
"""
            body += (
                f"  griddata[grid].params.CoordSystem_hash = {CoordSystem.upper()};\n"
            )
            body += f"  griddata[grid].params.grid_physical_size = {list_of_grid_physical_sizes[which_CoordSystem]};\n"
            body += "  numerical_grid_params_Nxx_dxx_xx(commondata, &griddata[grid].params, griddata[grid].xx, Nx, set_xxmin_xxmax_to_defaults);\n"
            body += (
                "  memcpy(&griddata_host[grid].params, &griddata[grid].params, sizeof(params_struct));\n"
                if parallelization in ["cuda"]
                else ""
            )
            body += "  grid++;\n\n"
        body += "}\n"
    elif gridding_approach == "multipatch":
        # fmt: off
        _ = par.CodeParameter("char[200]", __name__, "multipatch_choice", "", commondata=True, add_to_parfile=True)
        # fmt: on
        unit_vector_dict = {"x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1]}
        for dirn in ["x", "y", "z"]:
            # Direction of unit vectors relative to original, accounting for accumulation of regrids."
            _ = par.register_CodeParameter(
                "REAL[3]",
                __name__,
                f"cumulatively_rotated_{dirn}hatU",
                unit_vector_dict[
                    dirn
                ],  # Set below in C code when calling_for_first_time.
                commondata=False,
                add_to_parfile=False,
                add_to_set_CodeParameters_h=False,
            )
        body += """
  // Step 1.c: Multipatch grid structures are set up algorithmically.
  multipatch_grids_set_up(commondata, griddata);
"""
    else:
        raise ValueError(
            f"""gridding_approach == "{gridding_approach}" not supported.
        Supported approaches include: "independent grid(s)" (default) and "multipatch"."""
        )
    body += "\n// Step 1.d: Allocate memory for and define reference-metric precomputation lookup tables\n"
    if enable_rfm_precompute:
        body += r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  BHAH_MALLOC(griddata[grid].rfmstruct, sizeof(rfm_struct))
"""
        if parallelization in ["cuda"]:
            body = body.replace("BHAH_MALLOC", "BHAH_MALLOC_DEVICE")
            body += "\ncpyHosttoDevice_params__constant(&griddata[grid].params, griddata[grid].params.grid_idx % NUM_STREAMS);"
        body += r"""
    rfm_precompute_malloc(commondata, &griddata[grid].params, griddata[grid].rfmstruct);
    rfm_precompute_defines(commondata, &griddata[grid].params, griddata[grid].rfmstruct, griddata[grid].xx);
    }
"""
    else:
        body += "// (reference-metric precomputation disabled)\n"

    if parallelization in ["cuda"]:
        body += """
  cpyDevicetoHost__grid(commondata, griddata_host, griddata);
  BHAH_DEVICE_SYNC();
"""
    body += "\n// Step 1.e: Set up curvilinear boundary condition struct (bcstruct)\n"
    if enable_CurviBCs:
        body += (
            r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  bcstruct_set_up(commondata, &griddata[grid].params, griddata_host[grid].xx, &griddata_host[grid].bcstruct, &griddata[grid].bcstruct);
}
"""
            if parallelization in ["cuda"]
            else r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  bcstruct_set_up(commondata, &griddata[grid].params, griddata[grid].xx, &griddata[grid].bcstruct);
}
"""
        )
    else:
        body += "// (curvilinear boundary conditions bcstruct disabled)\n"
    if enable_set_cfl_timestep:
        sync_params = (
            "cpyHosttoDevice_params__constant(&griddata[grid].params, griddata[grid].params.grid_idx % NUM_STREAMS);"
            if parallelization in ["cuda"]
            else ""
        )
        body += rf"""
// Step 1.f: Set timestep based on minimum spacing between neighboring gridpoints.
commondata->dt = 1e30;
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {{
  {sync_params}
  cfl_limited_timestep(commondata, &griddata[grid].params, griddata[grid].xx);
}}"""
    body += r"""
// Step 1.g: Initialize timestepping parameters to zero if this is the first time this function is called.
if(calling_for_first_time) {
  commondata->nn = 0;
  commondata->nn_0 = 0;
  commondata->t_0 = 0.0;
  commondata->time = 0.0;
}

// Step 1.h: Set grid_idx for each grid.
for(int grid=0;grid<commondata->NUMGRIDS;grid++) {
   griddata[grid].params.grid_idx = grid;
}
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


def register_CFunctions(
    set_of_CoordSystems: Set[str],
    list_of_grid_physical_sizes: List[float],
    Nxx_dict: Dict[str, List[int]],
    gridding_approach: str = "independent grid(s)",
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
    enable_set_cfl_timestep: bool = True,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C functions related to coordinate systems and grid parameters.

    :param set_of_CoordSystems: Set of CoordSystems
    :param list_of_grid_physical_sizes: List of grid_physical_size for each CoordSystem; needed for Independent grids.
    :param Nxx_dict: Dictionary containing number of grid points.
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions.
    :param enable_set_cfl_timestep: Whether to enable computation of dt, the CFL timestep. A custom version can be implemented later.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    # Parallel codegen support.
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    for CoordSystem in set_of_CoordSystems:
        register_CFunction_numerical_grid_params_Nxx_dxx_xx(
            CoordSystem=CoordSystem,
            Nxx_dict=Nxx_dict,
        )
    if enable_set_cfl_timestep:
        register_CFunction_cfl_limited_timestep()
    register_CFunction_numerical_grids_and_timestep(
        set_of_CoordSystems=set_of_CoordSystems,
        list_of_grid_physical_sizes=list_of_grid_physical_sizes,
        gridding_approach=gridding_approach,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_CurviBCs=enable_CurviBCs,
        enable_set_cfl_timestep=enable_set_cfl_timestep,
    )

    if gridding_approach == "multipatch" or enable_set_cfl_timestep:
        for CoordSystem in set_of_CoordSystems:
            register_CFunction_ds_min_single_pt(CoordSystem)
    if gridding_approach == "multipatch":
        # Register regrid & masking functions
        for CoordSystem in set_of_CoordSystems:
            register_CFunction_ds_min_radial_like_dirns_single_pt(CoordSystem)

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
