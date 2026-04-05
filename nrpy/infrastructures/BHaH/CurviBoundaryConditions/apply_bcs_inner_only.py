# nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_inner_only.py
"""
Generates C functions for applying inner boundary conditions on a curvilinear grid.

The generated functions apply parity conditions to grid functions at inner boundary
points (e.g., for black hole excision). The boundary mappings are pre-computed
and stored in a `bc_struct`. This module provides routines to apply these
conditions to either all evolved grid functions or a specific list of auxiliary
grid functions.

This process is documented in the NRPy tutorial:
Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
         Samuel D. Tootle
         sdtootle **at** gmail **dot* com
"""

import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par

APPLY_PARITY_BRANCHLESS_PREFUNC = r"""
static inline REAL apply_parity_branchless(const REAL v, const int8_t p) {
#ifdef __cplusplus
  static_assert(sizeof(REAL) == sizeof(uint32_t) || sizeof(REAL) == sizeof(uint64_t), "REAL must be float or double");
#else
  _Static_assert(sizeof(REAL) == sizeof(uint32_t) || sizeof(REAL) == sizeof(uint64_t), "REAL must be float or double");
#endif

  if (sizeof(REAL) == sizeof(uint64_t)) {
    uint64_t bits;
    REAL out;
    memcpy(&bits, &v, sizeof(bits));
    bits ^= (uint64_t)(p < 0) << 63;
    memcpy(&out, &bits, sizeof(out));
    return out;
  } else {
    uint32_t bits;
    REAL out;
    memcpy(&bits, &v, sizeof(bits));
    bits ^= (uint32_t)(p < 0) << 31;
    memcpy(&out, &bits, sizeof(out));
    return out;
  } // END IF 64 bits vs 32 bits
} // END FUNCTION apply_parity_branchless
"""


def generate_apply_bcs_inner_only__kernel_body(
    loop_bounds: str = "NUM_EVOL_GFS",
    parity_ary: str = "evol_gf_parity",
    gf_index: str = "which_gf",
) -> str:
    """
    Generate the kernel body for the apply_bcs_inner_only function.

    This function is used to apply boundary conditions to inner boundary points
    on the computational grid.

    :param loop_bounds: The number of loop bounds (default is NUM_EVOL_GFS).
    :param parity_ary: The name of the parity array (default is evol_gf_parity).
    :param gf_index: The storage for the index of the grid function (default is which_gf).

    :return: The kernel body as a string.
    """
    parallelization = par.parval_from_str("parallelization")

    # Specify kernel body
    kernel_body = "// Needed for IDX macros\n"
    kernel_body += parallel_utils.get_loop_parameters(parallelization)
    kernel_body += (
        f"""
for(int which_gf=0;which_gf<{loop_bounds};which_gf++) {{
for (int pt = tid0; pt < num_inner_boundary_points; pt+=stride0) {{"""
        if parallelization == "cuda"
        else f"""
#pragma omp parallel for schedule(static)
  for(int which_gf=0;which_gf<{loop_bounds};++which_gf) {{
    const int parity_idx = {parity_ary}[{gf_index}];
    REAL *restrict gf = &gfs[IDX4pt({gf_index}, 0)];

    for(int pt=0;pt<num_inner_boundary_points;++pt) {{"""
    )
    kernel_body += (
        f"""
      const int dstpt = inner_bc_array[pt].dstpt;
      const int srcpt = inner_bc_array[pt].srcpt;
      gfs[IDX4pt({gf_index}, dstpt)] = inner_bc_array[pt].parity[{parity_ary}[{gf_index}]] * gfs[IDX4pt({gf_index}, srcpt)];
    }} // END for(int pt=0;pt<num_inner_pts;pt++)
  }} // END for(int which_gf=0;which_gf<{loop_bounds};which_gf++)
""".replace(f"{parity_ary}[{gf_index}]", f"d_{parity_ary}[{gf_index}]")
        if parallelization == "cuda"
        else """
      const innerpt_bc_struct *restrict bc = &inner_bc_array[pt];
      const REAL v = gf[bc->srcpt];
      const int8_t p = bc->parity[parity_idx];
      gf[bc->dstpt] = apply_parity_branchless(v, p);
    } // END for(int pt=0;pt<num_inner_pts;pt++)
  } // END for(int which_gf=0;which_gf<loop_bounds;which_gf++)
""".replace("loop_bounds", loop_bounds)
    )
    return kernel_body


###############################
## apply_bcs_inner_only(): Apply inner boundary conditions.
##  Function is documented below in desc= and body=.
def register_CFunction_apply_bcs_inner_only() -> None:
    """
    Register C function for filling inner boundary points on the computational grid, as prescribed by bcstruct.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    cfc.CFunction_dict.clear()
    ...    register_CFunction_apply_bcs_inner_only()
    ...    generated_str = cfc.CFunction_dict[f'apply_bcs_inner_only'].full_function
    ...    validation_desc = f"apply_bcs_inner_only__{parallelization}"
    ...    validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    includes = ["BHaH_defines.h"]
    desc = r"""
Apply BCs to inner boundary points only,
using data stored in bcstruct->inner_bc_array.
These structs are set in bcstruct_set_up().
Inner boundary points map to either the grid
interior ("pure inner") or to pure outer
boundary points ("inner maps to outer").
"""
    cfunc_type = "void"
    name = "apply_bcs_inner_only"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs"
    parallelization = par.parval_from_str("parallelization")

    kernel_body = generate_apply_bcs_inner_only__kernel_body()

    comments = "Apply BCs to inner boundary points only."
    # Prepare the argument dicts
    arg_dict_cuda = {
        "num_inner_boundary_points": "const int",
        "inner_bc_array": "const innerpt_bc_struct *restrict",
        "gfs": "REAL *restrict",
    }
    arg_dict_host = {
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }
    prefunc, new_body = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body,
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=comments,
        launch_dict={
            "blocks_per_grid": [
                "MAX(1U, (num_inner_boundary_points + threads_in_x_dir - 1) / threads_in_x_dir)"
            ],
            "stream": "params->grid_idx % NUM_STREAMS",
        },
        thread_tiling_macro_suffix="CURVIBC_INNER",
    )
    if parallelization != "cuda":
        prefunc = APPLY_PARITY_BRANCHLESS_PREFUNC + prefunc
    kernel_launch_body = rf"""
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;
  const innerpt_bc_struct *restrict inner_bc_array = bcstruct->inner_bc_array;
  const int num_inner_boundary_points = bc_info->num_inner_boundary_points;
  {new_body}"""
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=kernel_launch_body,
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
