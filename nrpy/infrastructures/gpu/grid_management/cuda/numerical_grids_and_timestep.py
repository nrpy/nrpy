"""
Register numerical_grids_and_timestep() C function, as well as functions called by this one.

These functions set up numerical grids for use within the BHaH infrastructure using
CUDA parallelization

Author: Samuel D. Tootle
        sdtootle **at** gmail **dot* com
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Dict, List

import nrpy.helpers.gpu_kernels.kernel_base as gputils
import nrpy.infrastructures.gpu.grid_management.base_numerical_grids_and_timestep as base_gm_classes
import nrpy.infrastructures.gpu.loop_utilities.cuda.simple_loop as lp
import nrpy.params as par

# fmt: off
for idx in range(3):
    _ = par.CodeParameter("int", __name__, f"Nxx_plus_2NGHOSTS{idx}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("int", __name__, f"Nxx{idx}", 64)
    # reference_metric sets xxmin and xxmax below.
    _ = par.CodeParameter("REAL", __name__, f"xxmin{idx}", -10.0, add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("REAL", __name__, f"xxmax{idx}", 10.0, add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("REAL", __name__, f"invdxx{idx}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("REAL", __name__, f"dxx{idx}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
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
    :param Nxx_dict: A dictionary that maps coordinate systems to lists containing the number of grid points along each direction.

    :return: None.
    :raises ValueError: If CoordSystem is not in Nxx_dict.
    """

    def __init__(
        self,
        CoordSystem: str,
        Nxx_dict: Dict[str, List[int]],
    ) -> None:
        super().__init__(CoordSystem, Nxx_dict)

        self.params.replace("REAL *restrict xx[3]", "REAL * xx[3]")
        self.prefunc = ""
        self.body += """
    // Allocate device storage
    cudaMalloc(&xx[0], sizeof(REAL) * Nxx_plus_2NGHOSTS0);
    cudaCheckErrors(malloc, "Malloc failed");
    cudaMalloc(&xx[1], sizeof(REAL) * Nxx_plus_2NGHOSTS1);
    cudaCheckErrors(malloc, "Malloc failed");
    cudaMalloc(&xx[2], sizeof(REAL) * Nxx_plus_2NGHOSTS2);
    cudaCheckErrors(malloc, "Malloc failed");

    cpyHosttoDevice_params__constant(params);

    dim3 block_threads, grid_blocks;
    auto set_grid_block = [&block_threads, &grid_blocks](auto Nx) {
        size_t threads_in_x_dir = 32;
        block_threads = dim3(threads_in_x_dir, 1, 1);
        grid_blocks = dim3((Nx + threads_in_x_dir - 1)/threads_in_x_dir, 1, 1);
    };

    size_t streamid = params->grid_idx % nstreams;
    set_grid_block(Nxx_plus_2NGHOSTS0);
    initialize_grid_xx0_gpu<<<grid_blocks, block_threads, 0, streams[streamid]>>>(xx[0]);
    cudaCheckErrors(initialize_grid_xx0_gpu, "kernel failed");

    streamid = (params->grid_idx + 1) % nstreams;
    set_grid_block(Nxx_plus_2NGHOSTS1);
    initialize_grid_xx1_gpu<<<grid_blocks, block_threads, 0, streams[streamid]>>>(xx[1]);
    cudaCheckErrors(initialize_grid_xx1_gpu, "kernel failed");

    streamid = (params->grid_idx + 2) % nstreams;
    set_grid_block(Nxx_plus_2NGHOSTS2);
    initialize_grid_xx2_gpu<<<grid_blocks, block_threads, 0, streams[streamid]>>>(xx[2]);
    cudaCheckErrors(initialize_grid_xx2_gpu, "kernel failed");
    """
        for i in range(3):
            kernel_body = f"""
  const int index  = blockIdx.x * blockDim.x + threadIdx.x;
  const int stride = blockDim.x * gridDim.x;

  REAL const xxmin{i} = d_params.xxmin{i};

  REAL const dxx{i} = d_params.dxx{i};

  REAL const Nxx_plus_2NGHOSTS{i} = d_params.Nxx_plus_2NGHOSTS{i};

  static constexpr REAL onehalf = 1.0 / 2.0;

  for (int j = index; j < Nxx_plus_2NGHOSTS{i}; j+=stride)
    xx{i}[j] = xxmin{i} + ((REAL)(j - NGHOSTS) + onehalf) * dxx{i};
"""
            xx0_kernel = gputils.GPU_Kernel(
                kernel_body,
                {f"xx{i}": "REAL *restrict"},
                f"initialize_grid_xx{i}_gpu",
                launch_dict={
                    "blocks_per_grid": [],
                    "threads_per_block": ["64"],
                    "stream": "default",
                },
            )
            self.prefunc += xx0_kernel.CFunction.full_function

        self.register()


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
const int Nxx_tot = (Nxx_plus_2NGHOSTS0)*(Nxx_plus_2NGHOSTS1)*(Nxx_plus_2NGHOSTS2);
  REAL *ds_min;
  REAL *restrict x0 = xx[0];
  REAL *restrict x1 = xx[1];
  REAL *restrict x2 = xx[2];

  // We only loop over a single GF array length
  cudaMalloc(&ds_min,sizeof(REAL) * Nxx_tot);
  cudaCheckErrors(cudaMalloc, "cudaMalloc failure"); // error checking
"""
        # lp_body = "REAL ds_min = 1e38;\n"
        lp_body = "REAL dsmin0, dsmin1, dsmin2;\n" + self.min_body_compute
        lp_body += """
  int idx = IDX3(i0,i1,i2);
  ds_min[idx] = MIN(dsmin0, MIN(dsmin1, dsmin2));
"""
        self.loop_body = ""
        for param_sym in self.unique_symbols:
            self.loop_body += f"const REAL {param_sym} = d_params.{param_sym};\n"
        self.loop_body += lp.simple_loop(
            loop_body=lp_body,
            read_xxs=True,
            loop_region="all points",
            fp_type=self.fp_type,
            CoordSystem=self.CoordSystem,
        ).full_loop_body

        # Put loop_body into a device kernel
        self.device_kernel = gputils.GPU_Kernel(
            self.loop_body,
            {
                "x0": "const REAL *restrict",
                "x1": "const REAL *restrict",
                "x2": "const REAL *restrict",
                "ds_min": "REAL *restrict",
            },
            "compute_ds_min__gpu",
            launch_dict={
                "blocks_per_grid": [],
                "threads_per_block": ["64"],
                "stream": "default",
            },
            fp_type=self.fp_type,
            comments="GPU Kernel to compute local ds_min per grid point.",
        )
        self.body += f"{self.device_kernel.launch_block}"
        self.body += f"{self.device_kernel.c_function_call()}"
        self.body += """
  REAL ds_min__global = find_global__minimum(ds_min, Nxx_tot);

  commondata->dt = MIN(commondata->dt, ds_min__global * commondata->CFL_FACTOR);
  cudaFree(ds_min);
"""

        self.prefunc = self.device_kernel.CFunction.full_function
        self.register()


class register_CFunction_numerical_grids_and_timestep(
    base_gm_classes.base_register_CFunction_numerical_grids_and_timestep
):
    """
    Register a C function to set up all numerical grids and timestep.

    The function configures the numerical grids based on given parameters, specifically
    focusing on the usage of reference metric precomputations and curvilinear boundary
    conditions.

    :param list_of_CoordSystems: List of CoordSystems
    :param list_of_grid_physical_sizes: List of grid_physical_size for each CoordSystem; needed for Independent grids.
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :param enable_rfm_precompute: Whether to enable reference metric precomputation (default: False).
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions (default: False).

    :raises ValueError: If invalid gridding_approach selected.
    """

    def __init__(
        self,
        list_of_CoordSystems: List[str],
        list_of_grid_physical_sizes: List[float],
        gridding_approach: str = "independent grid(s)",
        enable_rfm_precompute: bool = False,
        enable_CurviBCs: bool = False,
    ) -> None:
        super().__init__(
            list_of_CoordSystems,
            list_of_grid_physical_sizes,
            enable_rfm_precompute=enable_rfm_precompute,
            enable_CurviBCs=enable_CurviBCs,
        )
        self.params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata, "
        self.params += (
            "griddata_struct *restrict griddata_host, bool calling_for_first_time"
        )

        if self.enable_rfm_precompute:
            self.body += r"""
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  rfm_precompute_malloc(commondata, &griddata[grid].params, &griddata[grid].rfmstruct);
  cpyHosttoDevice_params__constant(&griddata[grid].params);
  rfm_precompute_defines(commondata, &griddata[grid].params, &griddata[grid].rfmstruct, griddata[grid].xx);
}
  cpyDevicetoHost__grid(commondata, griddata_host, griddata);
  cudaDeviceSynchronize();
"""
        else:
            self.body += "// (reference-metric precomputation disabled)\n"
        self.body += (
            "\n// Step 1.e: Set up curvilinear boundary condition struct (bcstruct)\n"
        )

        if self.enable_CurviBCs:
            self.body += r"""
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  cpyHosttoDevice_params__constant(&griddata[grid].params);
  bcstruct_set_up(commondata, &griddata[grid].params, griddata_host[grid].xx, &griddata[grid].bcstruct);
}
"""
        else:
            self.body += "// (curvilinear boundary conditions bcstruct disabled)\n"

        self.body += r"""
// Step 1.e: Set timestep based on minimum spacing between neighboring gridpoints.
commondata->dt = 1e30;
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  cpyHosttoDevice_params__constant(&griddata[grid].params);
  cfl_limited_timestep(commondata, &griddata[grid].params, griddata[grid].xx);
}

// Step 1.f: Initialize timestepping parameters to zero if this is the first time this function is called.
if(calling_for_first_time) {
  commondata->nn = 0;
  commondata->nn_0 = 0;
  commondata->t_0 = 0.0;
  commondata->time = 0.0;
}
"""
        self.register()


def register_CFunctions(
    list_of_CoordSystems: List[str],
    list_of_grid_physical_sizes: List[float],
    Nxx_dict: Dict[str, List[int]],
    gridding_approach: str = "independent grid(s)",
    enable_rfm_precompute: bool = False,
    enable_CurviBCs: bool = False,
    fp_type: str = "double",
) -> None:
    """
    Register C functions related to coordinate systems and grid parameters.

    :param list_of_CoordSystems: List of CoordSystems
    :param list_of_grid_physical_sizes: List of grid_physical_size for each CoordSystem; needed for Independent grids.
    :param Nxx_dict: Dictionary containing number of grid points.
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_CurviBCs: Whether to enable curvilinear boundary conditions.
    :param fp_type: Floating point type, e.g., "double".
    """
    for CoordSystem in list_of_CoordSystems:
        register_CFunction_numerical_grid_params_Nxx_dxx_xx(
            CoordSystem=CoordSystem,
            Nxx_dict=Nxx_dict,
        )
        register_CFunction_cfl_limited_timestep(
            CoordSystem=CoordSystem, fp_type=fp_type
        )
    register_CFunction_numerical_grids_and_timestep(
        list_of_CoordSystems=list_of_CoordSystems,
        list_of_grid_physical_sizes=list_of_grid_physical_sizes,
        gridding_approach=gridding_approach,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_CurviBCs=enable_CurviBCs,
    )

    if gridding_approach == "multipatch":
        # Register regrid & masking functions
        pass
