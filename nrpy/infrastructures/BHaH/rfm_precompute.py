"""
Set up basic functions and loop insertions for precomputed reference metric infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Any, Dict, List, Tuple

import sympy as sp
import sympy.codegen.ast as sp_ast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.helpers.expression_utils import get_unique_expression_symbols_as_strings
from nrpy.helpers.generic import superfast_uniq
from nrpy.helpers.gpu.utilities import generate_kernel_and_launch_code
from nrpy.infrastructures.BHaH.BHaH_defines_h import register_BHaH_defines


class ReferenceMetricPrecompute:
    """
    Base class for reference metric precomputation.

    This class stores contributions to BHaH_defines.h, as well as functions for memory allocation,
    definition, and freeing of rfm precomputation data. It also provides strings for reading rfm
    precompute quantities within loops for both SIMD-ized and non-SIMD loops.
    """

    def __init__(self, CoordSystem: str, parallelization: str = "openmp"):
        rfm = refmetric.reference_metric[CoordSystem + "_rfm_precompute"]
        # Step 7: Construct needed C code for declaring rfmstruct, allocating storage for
        #         rfmstruct arrays, defining each element in each array, reading the
        #         rfmstruct data from memory (both with and without SIMD enabled), and
        #         freeing allocated memory for the rfmstruct arrays.

        # BHaH_defines_list: String that contains the body of the rfmstruct struct.
        self.BHaH_defines_list: List[str] = []
        self.rfm_struct__define = ""
        # rfmstruct stores pointers to (so far) 1D arrays. The rfm_struct__malloc string allocates space for the arrays.
        self.rfm_struct__malloc = ""
        self.rfm_struct__freemem = ""

        # readvr_str reads the arrays from memory as needed
        self.readvr_str = ["", "", ""]
        self.readvr_intrinsics_outer_str = ["", "", ""]
        self.readvr_intrinsics_inner_str = ["", "", ""]

        # Sort freevars_uniq_vals and freevars_uniq_xx_indep, according to alphabetized freevars_uniq_xx_indep.
        #    Without this step, the ordering of elements in rfmstruct would be random, and would change each time
        #    this function was called.
        freevars_uniq_xx_indep: List[sp.Expr] = []
        freevars_uniq_vals: List[sp.Expr] = []
        if len(rfm.freevars_uniq_xx_indep) > 0:
            freevars_uniq_xx_indep, freevars_uniq_vals = (
                list(x)
                for x in zip(
                    *sorted(
                        zip(rfm.freevars_uniq_xx_indep, rfm.freevars_uniq_vals), key=str
                    )
                )
            )

        # Tease out how many variables each function in freevars_uniq_vals
        Nxx_plus_2NGHOSTS = [
            sp.Symbol("Nxx_plus_2NGHOSTS0", real=True),
            sp.Symbol("Nxx_plus_2NGHOSTS1", real=True),
            sp.Symbol("Nxx_plus_2NGHOSTS2", real=True),
        ]
        which_freevar: int = 0
        fp_ccg_type = ccg.fp_type_to_sympy_type[par.parval_from_str("fp_type")]
        sp_type_alias = {sp_ast.real: fp_ccg_type}
        self.rfm_struct__define_kernel_dict: Dict[sp.Expr, Any] = {}
        for expr in freevars_uniq_vals:
            if "_of_xx" in str(freevars_uniq_xx_indep[which_freevar]):
                frees = list(expr.free_symbols)
                frees_uniq = superfast_uniq(frees)
                xx_list: List[sp.Basic] = []
                malloc_size: int = 1
                for i in range(3):
                    if rfm.xx[i] in frees_uniq:
                        xx_list.append(rfm.xx[i])
                        malloc_size *= Nxx_plus_2NGHOSTS[i]

                self.BHaH_defines_list += [
                    f"REAL *restrict {freevars_uniq_xx_indep[which_freevar]};"
                ]
                self.rfm_struct__malloc += f"rfmstruct->{freevars_uniq_xx_indep[which_freevar]} = (REAL *)malloc(sizeof(REAL)*params->{malloc_size});\n"
                self.rfm_struct__freemem += (
                    f"free(rfmstruct->{freevars_uniq_xx_indep[which_freevar]});\n"
                )

                output_define_and_readvr = False
                for dirn in range(3):
                    if (
                        (rfm.xx[dirn] in frees_uniq)
                        and not (rfm.xx[(dirn + 1) % 3] in frees_uniq)
                        and not (rfm.xx[(dirn + 2) % 3] in frees_uniq)
                    ):
                        defines_kernel_body = (
                            "// Kernel thread/stride setup\n"
                            "const int tid0 = threadIdx.x + blockIdx.x*blockDim.x;\n"
                            "const int stride0 = blockDim.x * gridDim.x;\n\n"
                            f"for(int i{dirn}=tid0;i{dirn}<d_params[streamid].Nxx_plus_2NGHOSTS{dirn};i{dirn}+=stride0) {{\n"
                            if parallelization == "cuda"
                            else f"for(int i{dirn}=0;i{dirn}<params->Nxx_plus_2NGHOSTS{dirn};i{dirn}++) {{\n"
                        )
                        defines_kernel_body += (
                            f"  const REAL xx{dirn} = x{dirn}[i{dirn}];\n"
                            f"  rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i{dirn}] = {sp.ccode(freevars_uniq_vals[which_freevar], type_aliases=sp_type_alias)};\n"
                            "}"
                        )
                        # This is needed by register_CFunctions_rfm_precompute
                        self.rfm_struct__define_kernel_dict[
                            freevars_uniq_xx_indep[which_freevar]
                        ] = {
                            "body": defines_kernel_body,
                            "expr": freevars_uniq_vals[which_freevar],
                            "coord": [f"x{dirn}"],
                        }
                        self.readvr_str[
                            dirn
                        ] += f"MAYBE_UNUSED const REAL {freevars_uniq_xx_indep[which_freevar]} = rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i{dirn}];\n"
                        self.readvr_intrinsics_outer_str[
                            dirn
                        ] += f"const double NOSIMD{freevars_uniq_xx_indep[which_freevar]} = rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i{dirn}]; "
                        self.readvr_intrinsics_outer_str[
                            dirn
                        ] += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {freevars_uniq_xx_indep[which_freevar]} = ConstSIMD(NOSIMD{freevars_uniq_xx_indep[which_freevar]});\n"
                        self.readvr_intrinsics_inner_str[
                            dirn
                        ] += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {freevars_uniq_xx_indep[which_freevar]} = ReadSIMD(&rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i{dirn}]);\n"
                        output_define_and_readvr = True

                if (
                    (not output_define_and_readvr)
                    and (rfm.xx[0] in frees_uniq)
                    and (rfm.xx[1] in frees_uniq)
                ):
                    defines_kernel_body = (
                        "const int Nxx_plus_2NGHOSTS0 = d_params[streamid].Nxx_plus_2NGHOSTS0;\n"
                        "const int Nxx_plus_2NGHOSTS1 = d_params[streamid].Nxx_plus_2NGHOSTS1;\n\n"
                        "// Kernel thread/stride setup\n"
                        "const int tid0 = threadIdx.x + blockIdx.x*blockDim.x;\n"
                        "const int tid1 = blockIdx.y * blockDim.y + threadIdx.y;\n\n"
                        "const int stride0 = blockDim.x * gridDim.x;\n"
                        "const int stride1 = blockDim.y * gridDim.y;\n\n"
                        "for(int i1=tid1;i1<Nxx_plus_2NGHOSTS1;i1+=stride1) for(int i0=tid0;i0<Nxx_plus_2NGHOSTS0;i0+=stride0) {"
                        if parallelization == "cuda"
                        else "const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;\n"
                        "const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;\n\n"
                        "for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {"
                    )
                    defines_kernel_body += f"""
                  const REAL xx0 = x0[i0];
                  const REAL xx1 = x1[i1];
                  rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i0 + Nxx_plus_2NGHOSTS0*i1] = {sp.ccode(freevars_uniq_vals[which_freevar], type_aliases=sp_type_alias)};
                }}\n\n"""
                    # This is needed by register_CFunctions_rfm_precompute
                    self.rfm_struct__define_kernel_dict[
                        freevars_uniq_xx_indep[which_freevar]
                    ] = {
                        "body": defines_kernel_body,
                        "expr": freevars_uniq_vals[which_freevar],
                        "coord": ["x0", "x1"],
                    }
                    self.readvr_str[
                        0
                    ] += f"MAYBE_UNUSED const REAL {freevars_uniq_xx_indep[which_freevar]} = rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i0 + Nxx_plus_2NGHOSTS0*i1];\n"
                    self.readvr_intrinsics_outer_str[
                        0
                    ] += f"const double NOSIMD{freevars_uniq_xx_indep[which_freevar]} = rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i0 + Nxx_plus_2NGHOSTS0*i1]; "
                    self.readvr_intrinsics_outer_str[
                        0
                    ] += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {freevars_uniq_xx_indep[which_freevar]} = ConstSIMD(NOSIMD{freevars_uniq_xx_indep[which_freevar]});\n"
                    self.readvr_intrinsics_inner_str[
                        0
                    ] += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {freevars_uniq_xx_indep[which_freevar]} = ReadSIMD(&rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i0 + Nxx_plus_2NGHOSTS0*i1]);\n"
                    output_define_and_readvr = True

                if not output_define_and_readvr:
                    raise RuntimeError(
                        f"ERROR: Could not figure out the (xx0,xx1,xx2) dependency within the expression for {freevars_uniq_xx_indep[which_freevar]}: {freevars_uniq_vals[which_freevar]}"
                    )
            if parallelization == "cuda":
                self.readvr_intrinsics_outer_str = [
                    s.replace("SIMD", "CUDA") for s in self.readvr_intrinsics_outer_str
                ]
                self.readvr_intrinsics_inner_str = [
                    s.replace("SIMD", "CUDA") for s in self.readvr_intrinsics_inner_str
                ]

            which_freevar += 1


def generate_rfmprecompute_defines(
    rfm_precompute: ReferenceMetricPrecompute, parallelization: str = "openmp"
) -> Tuple[str, str]:
    """
    Generate the body and prefunctions for the rfm precompute defines.

    :param rfm_precompute: ReferenceMetricPrecompute object.
    :param parallelization: Parallelization method to use.
    :return: Prefunction and body strings.
    """
    prefunc = ""
    body = ""
    func_name = "defines"
    kernel_dicts = rfm_precompute.rfm_struct__define_kernel_dict

    for i, (key_sym, kernel_dict) in enumerate(kernel_dicts.items()):
        # Gather unique symbols, excluding xx0, xx1, xx2
        unique_symbols = get_unique_expression_symbols_as_strings(
            kernel_dict["expr"], exclude=[f"xx{j}" for j in range(3)]
        )

        # Build the kernel body, pulling param access
        kernel_body = "// Temporary parameters\n"
        for sym in unique_symbols:
            kernel_body += f"const REAL {sym} = params->{sym};\n"
        kernel_body += kernel_dict["body"]

        # Kernel name
        kernel_name = f"rfm_precompute_{func_name}__{key_sym}"
        # Comments
        comments = f"Kernel to precompute metric quantity {key_sym}."

        # Prepare the argument dicts
        arg_dict_cuda = {
            "rfmstruct": "rfm_struct *restrict",
        }
        arg_dict_host = {
            "params": "const params_struct *restrict",
            "rfmstruct": "rfm_struct *restrict",
        }
        # Add coordinate arrays from kernel_dict["coord"]
        for coord_var in kernel_dict["coord"]:
            arg_dict_cuda[coord_var] = "const REAL *restrict"
            arg_dict_host[coord_var] = "const REAL *restrict"

        # We replicate the original approach of computing param_streamid in the body for CUDA
        # so define_param_streamid=True, which forces the helper to insert that line.
        new_prefunc, new_body = generate_kernel_and_launch_code(
            kernel_name,
            kernel_body,
            arg_dict_cuda,
            arg_dict_host,
            parallelization,
            cfunc_type="static void",
            comments=comments,
        )
        prefunc += new_prefunc
        body += new_body

    return prefunc, body


def generate_rfmprecompute_malloc(
    rfm_precompute: ReferenceMetricPrecompute, parallelization: str = "openmp"
) -> Tuple[str, str]:
    """
    Generate the body and prefunctions for allocating the rfmstruct arrays.

    :param rfm_precompute: ReferenceMetricPrecompute object.
    :param parallelization: Parallelization method to use.
    :return: Prefunction and body strings.
    """
    # Build the kernel body
    kernel_body = "// Temporary parameters\n"
    for i in range(3):
        kernel_body += (
            f"MAYBE_UNUSED const int Nxx_plus_2NGHOSTS{i} = "
            f"params->Nxx_plus_2NGHOSTS{i};\n"
        )
    kernel_body += rfm_precompute.rfm_struct__malloc

    kernel_name = "rfm_precompute_malloc__allocate"
    comments = "Kernel to allocate rfmstruct arrays."

    # Arg dict
    arg_dict_cuda = {"rfmstruct": "rfm_struct *restrict"}
    arg_dict_host = {
        "params": "const params_struct *restrict",
        "rfmstruct": "rfm_struct *restrict",
    }

    prefunc, body = generate_kernel_and_launch_code(
        kernel_name,
        kernel_body,
        arg_dict_cuda,
        arg_dict_host,
        parallelization,
        comments=comments,
        # For "cudaMalloc", only a single thread is needed
        launch_dict={
            "blocks_per_grid": ["1"],
            "threads_per_block": ["1"],
            "stream": "default",
        },
    )
    return prefunc, body


def generate_rfmprecompute_free(
    rfm_precompute: ReferenceMetricPrecompute, parallelization: str = "openmp"
) -> Tuple[str, str]:
    """
    Generate the body and prefunctions for deallocating the rfmstruct arrays.

    :param rfm_precompute: ReferenceMetricPrecompute object.
    :param parallelization: Parallelization method to use.
    :return: Prefunction and body strings.
    """
    # Build the kernel body
    kernel_body = "// Temporary parameters\n"
    kernel_body += rfm_precompute.rfm_struct__freemem

    kernel_name = "rfm_precompute_free__deallocate"
    comments = "Kernel to deallocate rfmstruct arrays."

    # Arg dict
    arg_dict_cuda = {"rfmstruct": "rfm_struct *restrict"}
    arg_dict_host = {"rfmstruct": "rfm_struct *restrict"}

    prefunc, body = generate_kernel_and_launch_code(
        kernel_name,
        kernel_body,
        arg_dict_cuda,
        arg_dict_host,
        parallelization,
        comments=comments,
        # For "cudaMalloc", only a single thread is needed
        launch_dict={
            "blocks_per_grid": ["1"],
            "threads_per_block": ["1"],
        },
    )
    return prefunc, body


def register_CFunctions_rfm_precompute(
    list_of_CoordSystems: List[str],
    parallelization: str = "openmp",
) -> None:
    """
    Register C functions for reference metric precomputed lookup arrays.

    :param list_of_CoordSystems: List of coordinate systems to register the C functions.
    :param parallelization: Parallelization method to use.

    Doctest:
    >>> import nrpy.c_function as cfc
    >>> from nrpy.infrastructures.BHaH import rfm_precompute
    >>> from nrpy.reference_metric import supported_CoordSystems
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.params as par
    >>> par.set_parval_from_str("fp_type", "float")
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> for parallelization in supported_Parallelizations:
    ...    for CoordSystem in supported_CoordSystems:
    ...       cfc.CFunction_dict.clear()
    ...       rfm_precompute.register_CFunctions_rfm_precompute([CoordSystem], parallelization=parallelization) # doctest: +SKIP
    ...       for rfm_base_function in ["malloc", "defines", "free"]:
    ...          generated_str = cfc.CFunction_dict[f'rfm_precompute_{rfm_base_function}__rfm__{CoordSystem}'].full_function
    ...          validation_desc = f"{rfm_base_function}__{parallelization}__{CoordSystem}".replace(" ", "_")
    ...          validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    combined_BHaH_defines_list = []
    for CoordSystem in list_of_CoordSystems:
        rfm_precompute = ReferenceMetricPrecompute(
            CoordSystem, parallelization=parallelization
        )

        includes = ["BHaH_defines.h"]
        cfunc_type = "void"

        body = ""
        for i in range(3):
            body += f"MAYBE_UNUSED const REAL *restrict x{i} = xx[{i}];\n"
            body += (
                f"MAYBE_UNUSED const int Nxx_plus_2NGHOSTS{i} = params->Nxx_plus_2NGHOSTS{i};\n"
                if parallelization == "cuda"
                else ""
            )

        defines_prefunc, defines_body = generate_rfmprecompute_defines(
            rfm_precompute, parallelization=parallelization
        )
        malloc_prefunc, malloc_body = generate_rfmprecompute_malloc(
            rfm_precompute, parallelization=parallelization
        )
        free_prefunc, free_body = generate_rfmprecompute_free(
            rfm_precompute, parallelization=parallelization
        )

        rfm_precompute.rfm_struct__define = body + defines_body
        rfm_precompute.rfm_struct__malloc = malloc_body
        rfm_precompute.rfm_struct__freemem = free_body

        prefunc_dict = {
            "malloc": malloc_prefunc,
            "defines": defines_prefunc,
            "free": free_prefunc,
        }
        for func in [
            ("malloc", rfm_precompute.rfm_struct__malloc),
            ("defines", rfm_precompute.rfm_struct__define),
            ("free", rfm_precompute.rfm_struct__freemem),
        ]:

            desc = f"rfm_precompute_{func[0]}: reference metric precomputed lookup arrays: {func[0]}"
            name = "rfm_precompute_" + func[0]
            params = "const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct"
            include_CodeParameters_h = False
            if func[0] == "defines":
                params += ", REAL *restrict xx[3]"

            body = " "
            body += func[1]

            defines_list = [
                s.replace("restrict", "") if parallelization == "cuda" else s
                for s in rfm_precompute.BHaH_defines_list
            ]
            prefunc = prefunc_dict[func[0]] if func[0] in prefunc_dict else ""
            combined_BHaH_defines_list.extend(defines_list)
            cfc.register_CFunction(
                prefunc=prefunc,
                includes=includes,
                desc=desc,
                cfunc_type=cfunc_type,
                CoordSystem_for_wrapper_func=CoordSystem,
                name=name,
                params=params,
                include_CodeParameters_h=include_CodeParameters_h,
                body=body,
            )

    BHaH_defines = "typedef struct __rfmstruct__ {\n"
    BHaH_defines += "\n".join(sorted(superfast_uniq(combined_BHaH_defines_list)))
    BHaH_defines += "\n} rfm_struct;\n"
    register_BHaH_defines("reference_metric", BHaH_defines)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
