"""
C function for hyperbolic relaxation diagnostics.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.helpers.expression_utils import get_params_commondata_symbols_from_expr_list


# Define function to compute the l^2 of a gridfunction
def register_CFunction_log10_L2norm_gf(
    CoordSystem: str,
) -> None:
    """
    Register function to compute l2-norm of a gridfunction assuming a single grid.

    Note that parallel codegen is disabled for this function, as it sometimes causes a
    multiprocess race condition on Python 3.6.7

    :param CoordSystem: the rfm coordinate system.

    Doctest:
    >>> import nrpy.c_function as cfc
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import unittest_CoordSystems
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name="log10_L2norm_gf"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in unittest_CoordSystems:
    ...       cfc.CFunction_dict.clear()
    ...       _ = register_CFunction_log10_L2norm_gf(CoordSystem)
    ...       generated_str = cfc.CFunction_dict[f'{name}__rfm__{CoordSystem}'].full_function
    ...       validation_desc = f"{name}__{parallelization}__{CoordSystem}".replace(" ", "_")
    ...       validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    Setting up reference_metric[SinhSymTP]...
    Setting up reference_metric[HoleySinhSpherical]...
    Setting up reference_metric[Cartesian]...
    Setting up reference_metric[SinhCylindricalv2n2]...
    """
    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Compute l2-norm of a gridfunction assuming a single grid."
    cfunc_type = "void"
    name = "log10_L2norm_gf"
    params = """commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3],
                const REAL integration_radius, const int gf_index, REAL * l2norm, const REAL *restrict in_gfs"""
    params += ", REAL *restrict aux_gfs" if parallelization in ["cuda"] else ""

    rfm = refmetric.reference_metric[CoordSystem]

    fp_type = par.parval_from_str("fp_type")

    # Enforce at least double precision since calcuations
    # can go beyond single precision numerical limits for some
    # coordinate systems
    fp_type_alias = "DOUBLE" if fp_type == "float" else "REAL"
    ccg_fp_type = "double" if fp_type == "float" else fp_type
    expr_list = [
        rfm.xxSph[0],
        rfm.detgammahat,
    ]

    # Define the norm calculations within the loop
    reduction_loop_body = ccg.c_codegen(
        expr_list,
        [
            "const DOUBLE r",
            "const DOUBLE sqrtdetgamma",
        ],
        include_braces=False,
        fp_type=ccg_fp_type,
        fp_type_alias=fp_type_alias,
    )

    reduction_loop_body += r"""
if(r < integration_radius) {
  const DOUBLE gf_of_x = in_gfs[IDX4(gf_index, i0, i1, i2)];
  const DOUBLE dV = sqrtdetgamma * dxx0 * dxx1 * dxx2;
"""
    reduction_loop_body += (
        r"""
  aux_gfs[IDX4(L2_SQUARED_DVGF, i0, i1, i2)] = gf_of_x * gf_of_x * dV;
  aux_gfs[IDX4(L2_DVGF, i0, i1, i2)] = dV;
} // END if(r < integration_radius)
"""
        if parallelization in ["cuda"]
        else r"""
  squared_sum += gf_of_x * gf_of_x * dV;
  volume_sum  += dV;
} // END if(r < integration_radius)
"""
    )
    reduction_loop_body = reduction_loop_body.replace("REAL", fp_type_alias)

    OMP_custom_pragma = (
        r"#pragma omp parallel for reduction(+:squared_sum,volume_sum)"
        if parallelization not in ["cuda"]
        else ""
    )

    # Generate the loop for the reduction_loop_body
    loop_body = lp.simple_loop(
        loop_body="\n" + reduction_loop_body,
        read_xxs=True,
        loop_region="interior",
        OMP_custom_pragma=OMP_custom_pragma,
    )

    # Device code computes the local L2 quantities which are reduced
    # in a separate algorithm, find_global__sum.
    # Host code uses OpenMP reduction #pragma to compute the
    # L2 quantities and the global norms in a single kernel
    comments = (
        "Kernel to compute L2 quantities pointwise (not summed)."
        if parallelization in ["cuda"]
        else "Kernel to compute L2 quantities pointwise (summed)."
    )
    # Prepare the argument dictionaries
    arg_dict_cuda = {
        "x0": "const REAL *restrict",
        "x1": "const REAL *restrict",
        "x2": "const REAL *restrict",
        "in_gfs": "const REAL *restrict",
        "aux_gfs": "REAL *restrict",
        "integration_radius": "const REAL",
        "gf_index": "const int",
    }
    arg_dict_host = {
        "params": "const params_struct *restrict",
        "x0": "const REAL *restrict",
        "x1": "const REAL *restrict",
        "x2": "const REAL *restrict",
        "in_gfs": "const REAL *restrict",
        "squared_sum_final": "REAL *restrict",
        "volume_sum_final": "REAL *restrict",
        "integration_radius": "const REAL",
        "gf_index": "const int",
    }
    loop_params = parallel_utils.get_loop_parameters(parallelization)
    params_symbols, _ = get_params_commondata_symbols_from_expr_list(
        expr_list, exclude=[f"xx{i}" for i in range(3)]
    )
    loop_params += "// Load necessary parameters from params_struct\n"
    # We have to manually add dxx{i} here since they are not in the SymPy
    # expression list, but, instead, are manually used in calculations
    # in reduction_loop_body above
    for param in params_symbols + [f"dxx{i}" for i in range(3)]:
        loop_params += f"const REAL {param} = {parallel_utils.get_params_access(parallelization)}{param};\n"
    loop_params += (
        "\n"
        if parallelization in ["cuda"]
        else r"""
  DOUBLE squared_sum = 0.0;
  DOUBLE volume_sum  = 0.0;
        """
    )

    kernel_body = f"{loop_params}\n{loop_body}"
    for i in range(3):
        kernel_body = kernel_body.replace(f"xx[{i}]", f"x{i}")

    if parallelization not in ["cuda"]:
        kernel_body += r"""*squared_sum_final = squared_sum;
*volume_sum_final = volume_sum;
        """

    prefunc, new_body = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body,
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=comments,
        launch_dict={
            "blocks_per_grid": [],
            "threads_per_block": ["32", "NGHOSTS"],
            "stream": "default",
        },
        thread_tiling_macro_suffix="NELL_GRIDL2",
    )

    # Define launch kernel body
    body = r"""
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
  REAL *restrict x0 = xx[0];
  REAL *restrict x1 = xx[1];
  REAL *restrict x2 = xx[2];
"""

    body += (
        r"""
  // Since we're performing sums, make sure arrays are zero'd
  cudaMemset(aux_gfs, 0, sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
"""
        if parallelization in ["cuda"]
        else r"""
  // Set summation variables to compute l2-norm
  DOUBLE squared_sum = 0.0;
  DOUBLE volume_sum  = 0.0;
"""
    )

    body += f"{new_body}\n"
    body += (
        r"""
  // Set summation variables to compute l2-norm
  REAL squared_sum = find_global__sum(&aux_gfs[IDX4(L2_SQUARED_DVGF, 0, 0, 0)], Nxx_plus_2NGHOSTS_tot);
  REAL volume_sum = find_global__sum(&aux_gfs[IDX4(L2_DVGF, 0, 0, 0)], Nxx_plus_2NGHOSTS_tot);
  // Compute and output the log of the l2-norm.
  REAL local_norm = log10(1e-16 + sqrt(squared_sum / volume_sum));  // 1e-16 + ... avoids log10(0)
"""
        if parallelization in ["cuda"]
        else ""
    )

    # For host reduction code.
    body = body.replace("squared_sum_final", "&squared_sum").replace(
        "volume_sum_final", "&volume_sum"
    )

    body += r"""
  // Compute and output the log of the l2-norm.
  *l2norm = log10(1e-16 + sqrt(squared_sum / volume_sum));  // 1e-16 + ... avoids log10(0)
"""

    cfc.register_CFunction(
        subdirectory="diagnostics",
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
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
