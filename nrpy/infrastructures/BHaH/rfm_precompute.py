"""
Set up basic functions and loop insertions for precomputed reference metric infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Any, Dict, List, Set, Tuple

import sympy as sp
import sympy.codegen.ast as sp_ast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.helpers.expression_utils import get_unique_expression_symbols_as_strings
from nrpy.helpers.generic import superfast_uniq
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code
from nrpy.infrastructures.BHaH.BHaH_defines_h import register_BHaH_defines


class ReferenceMetricPrecompute:
    """
    Base class for reference metric precomputation.

    This class stores contributions to BHaH_defines.h, as well as functions for memory allocation,
    definition, and freeing of rfm precomputation data. It also provides strings for reading rfm
    precompute quantities within loops for both SIMD-ized and non-SIMD loops.
    """

    def __init__(self, CoordSystem: str) -> None:
        """
        Initialize and orchestrate the generation of C code for precomputation.

        This constructor acts as a high-level coordinator, delegating tasks to specialized
        helper methods for configuration, initialization, and C code generation.

        :param CoordSystem: The coordinate system to use for precomputation.
        """
        # Step 7: Construct needed C code for declaring rfmstruct, allocating storage for
        #         rfmstruct arrays, defining each element in each array, reading the
        #         rfmstruct data from memory (both with and without SIMD enabled), and
        #         freeing allocated memory for the rfmstruct arrays.

        self._initialize_configuration(CoordSystem)
        self._initialize_properties()

        # Sort expressions to ensure deterministic code generation.
        sorted_expressions = self._get_sorted_precomputed_expressions()

        for symbol, expr in sorted_expressions:
            # We only process expressions that depend on grid coordinates ("_of_xx").
            if "_of_xx" in str(symbol):
                self._process_expression(symbol, expr)

        # Final adjustments for CUDA if necessary.
        if self.is_cuda:
            self._use_cuda_intrinsics()

        # Finalize the 'defines' C function body.
        self._finalize_rfm_struct_define()

    def _initialize_configuration(self, CoordSystem: str) -> None:
        """
        Set up basic configuration properties based on global parameters.

        :param CoordSystem: The coordinate system to use for precomputation.
        """
        self.N_DIMS = 3
        self.parallelization = par.parval_from_str("parallelization")
        self.is_cuda = self.parallelization == "cuda"
        self.rfm = refmetric.reference_metric[CoordSystem + "_rfm_precompute"]
        fp_type = par.parval_from_str("fp_type")
        self.sp_type_alias = {sp_ast.real: ccg.fp_type_to_sympy_type[fp_type]}
        self.nxx_plus_2NGHOSTS_syms = [
            sp.Symbol(f"Nxx_plus_2NGHOSTS{i}", real=True) for i in range(self.N_DIMS)
        ]

    def _initialize_properties(self) -> None:
        """Initialize storage for C code snippets and data structures."""
        # BHaH_defines_list: String that contains the body of the rfmstruct struct.
        self.BHaH_defines_list: List[str] = []
        # rfmstruct stores pointers to (so far) 1D arrays.
        self.rfm_struct__malloc = ""
        self.rfm_struct__freemem = ""
        # The rfm_struct__define string will be populated later, after generating kernels.
        self.rfm_struct__define = ""
        self.rfm_struct__define_kernel_dict: Dict[sp.Expr, Any] = {}

        # readvr_str reads the arrays from memory as needed
        self.readvr_str = [""] * self.N_DIMS
        self.readvr_intrinsics_outer_str = [""] * self.N_DIMS
        self.readvr_intrinsics_inner_str = [""] * self.N_DIMS

    def _get_sorted_precomputed_expressions(
        self,
    ) -> List[Tuple[sp.Expr, sp.Expr]]:
        """
        Sort precomputed expressions to ensure deterministic code generation.

        Sorting is based on the string representation of the symbol name.
        Without this step, the ordering of elements in rfmstruct would be random.

        :return: A list of (symbol, expression) tuples, sorted by symbol name.
        """
        if not self.rfm.freevars_uniq_xx_indep:
            return []

        zipped_pairs = zip(self.rfm.freevars_uniq_xx_indep, self.rfm.freevars_uniq_vals)
        return sorted(zipped_pairs, key=lambda pair: str(pair[0]))

    def _process_expression(self, symbol: sp.Expr, expr: sp.Expr) -> None:
        """
        Process a single precomputed expression, dispatching to the correct handler.

        This method analyzes the expression's coordinate dependencies and calls the
        appropriate function to generate the C code for it.

        :param symbol: The SymPy symbol for the precomputed quantity.
        :param expr: The SymPy expression for the precomputed quantity.
        :raises RuntimeError: If the coordinate dependency is not supported.
        """
        free_symbols = superfast_uniq(list(expr.free_symbols))
        dependencies = [self.rfm.xx[i] in free_symbols for i in range(self.N_DIMS)]

        self._add_memory_management_code(symbol, dependencies)

        # Check for supported dependency patterns (1D or 2D in xx0, xx1).
        is_1d = sum(dependencies) == 1
        is_2d_01 = dependencies[0] and dependencies[1] and not dependencies[2]

        if is_1d:
            dep_dir = dependencies.index(True)
            self._handle_1d_dependency(symbol, expr, dep_dir)
        elif is_2d_01:
            self._handle_2d_dependency(symbol, expr)
        else:
            raise RuntimeError(
                f"ERROR: Could not figure out the (xx0,xx1,xx2) dependency "
                f"within the expression for {symbol}: {expr}"
            )

    def _add_memory_management_code(
        self, symbol: sp.Expr, dependencies: List[bool]
    ) -> None:
        """
        Generate and store C code for memory allocation and deallocation.

        :param symbol: The SymPy symbol for the precomputed quantity.
        :param dependencies: A boolean list indicating coordinate dependencies.
        """
        malloc_size_expr = sp.sympify(1)
        for i, is_dep in enumerate(dependencies):
            if is_dep:
                malloc_size_expr *= self.nxx_plus_2NGHOSTS_syms[i]

        malloc_macro = (
            "BHAH_MALLOC_DEVICE__PtrMember"
            if self.is_cuda
            else "BHAH_MALLOC__PtrMember"
        )
        free_macro = (
            "BHAH_FREE_DEVICE__PtrMember" if self.is_cuda else "BHAH_FREE__PtrMember"
        )

        self.BHaH_defines_list.append(f"REAL *restrict {symbol};")
        self.rfm_struct__malloc += f"{malloc_macro}(rfmstruct, {symbol}, sizeof(REAL)*params->{malloc_size_expr});"
        self.rfm_struct__freemem += f"{free_macro}(rfmstruct, {symbol});"

    def _handle_1d_dependency(
        self, symbol: sp.Expr, expr: sp.Expr, dep_dir: int
    ) -> None:
        """
        Handle expressions with a 1D coordinate dependency.

        :param symbol: The SymPy symbol for the precomputed quantity.
        :param expr: The SymPy expression for the precomputed quantity.
        :param dep_dir: The direction (0, 1, or 2) of the 1D dependency.
        """
        if self.is_cuda:
            kernel_body = (
                "// Kernel thread/stride setup\n"
                "const int tid0 = threadIdx.x + blockIdx.x*blockDim.x;\n"
                "const int stride0 = blockDim.x * gridDim.x;\n\n"
                f"for(int i{dep_dir}=tid0; i{dep_dir}<d_params[streamid].Nxx_plus_2NGHOSTS{dep_dir}; i{dep_dir}+=stride0) {{\n"
            )
        else:
            kernel_body = f"for(int i{dep_dir}=0; i{dep_dir}<params->Nxx_plus_2NGHOSTS{dep_dir}; i{dep_dir}++) {{\n"

        kernel_body += (
            f"  const REAL xx{dep_dir} = x{dep_dir}[i{dep_dir}];\n"
            f"  rfmstruct->{symbol}[i{dep_dir}] = {sp.ccode(expr, type_aliases=self.sp_type_alias)};\n"
            "}"
        )
        self.rfm_struct__define_kernel_dict[symbol] = {
            "body": kernel_body,
            "expr": expr,
            "coord": [f"x{dep_dir}"],
        }

        # Reader strings
        self.readvr_str[
            dep_dir
        ] += f"MAYBE_UNUSED const REAL {symbol} = rfmstruct->{symbol}[i{dep_dir}];\n"
        self.readvr_intrinsics_outer_str[
            dep_dir
        ] += f"const double NOSIMD{symbol} = rfmstruct->{symbol}[i{dep_dir}]; "
        self.readvr_intrinsics_outer_str[
            dep_dir
        ] += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol} = ConstSIMD(NOSIMD{symbol});\n"
        self.readvr_intrinsics_inner_str[
            dep_dir
        ] += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol} = ReadSIMD(&rfmstruct->{symbol}[i{dep_dir}]);\n"

    def _handle_2d_dependency(self, symbol: sp.Expr, expr: sp.Expr) -> None:
        """
        Handle expressions with a 2D coordinate dependency on xx0 and xx1.

        :param symbol: The SymPy symbol for the precomputed quantity.
        :param expr: The SymPy expression for the precomputed quantity.
        """
        mem_idx = "i0 + Nxx_plus_2NGHOSTS0*i1"

        if self.is_cuda:
            kernel_body = (
                "const int Nxx_plus_2NGHOSTS0 = d_params[streamid].Nxx_plus_2NGHOSTS0;\n"
                "const int Nxx_plus_2NGHOSTS1 = d_params[streamid].Nxx_plus_2NGHOSTS1;\n\n"
                "// Kernel thread/stride setup\n"
                "const int tid0 = threadIdx.x + blockIdx.x*blockDim.x;\n"
                "const int tid1 = blockIdx.y * blockDim.y + threadIdx.y;\n\n"
                "const int stride0 = blockDim.x * gridDim.x;\n"
                "const int stride1 = blockDim.y * gridDim.y;\n\n"
                "for(int i1=tid1;i1<Nxx_plus_2NGHOSTS1;i1+=stride1) for(int i0=tid0;i0<Nxx_plus_2NGHOSTS0;i0+=stride0) {"
            )
        else:
            kernel_body = (
                "const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;\n"
                "const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;\n\n"
                "for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {"
            )
        kernel_body += f"""
          const REAL xx0 = x0[i0];
          const REAL xx1 = x1[i1];
          rfmstruct->{symbol}[{mem_idx}] = {sp.ccode(expr, type_aliases=self.sp_type_alias)};
        }}\n\n"""

        self.rfm_struct__define_kernel_dict[symbol] = {
            "body": kernel_body,
            "expr": expr,
            "coord": ["x0", "x1"],
        }
        # For 2D dependencies, reader strings are added to the dirn=0 component.
        self.readvr_str[
            0
        ] += f"MAYBE_UNUSED const REAL {symbol} = rfmstruct->{symbol}[{mem_idx}];\n"
        self.readvr_intrinsics_outer_str[
            0
        ] += f"const double NOSIMD{symbol} = rfmstruct->{symbol}[{mem_idx}]; "
        self.readvr_intrinsics_outer_str[
            0
        ] += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol} = ConstSIMD(NOSIMD{symbol});\n"
        self.readvr_intrinsics_inner_str[
            0
        ] += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol} = ReadSIMD(&rfmstruct->{symbol}[{mem_idx}]);\n"

    def _use_cuda_intrinsics(self) -> None:
        """Replace SIMD intrinsic calls with CUDA equivalents."""
        self.readvr_intrinsics_outer_str = [
            s.replace("SIMD", "CUDA") for s in self.readvr_intrinsics_outer_str
        ]
        self.readvr_intrinsics_inner_str = [
            s.replace("SIMD", "CUDA") for s in self.readvr_intrinsics_inner_str
        ]

    def _finalize_rfm_struct_define(self) -> None:
        """Construct the final C code for the 'defines' function body."""
        body_prefix = ""
        for i in range(3):
            body_prefix += f"MAYBE_UNUSED const REAL *restrict x{i} = xx[{i}];\n"
            if self.is_cuda:
                body_prefix += f"MAYBE_UNUSED const int Nxx_plus_2NGHOSTS{i} = params->Nxx_plus_2NGHOSTS{i};\n"

        _, defines_body = generate_rfmprecompute_defines(self)
        self.rfm_struct__define = body_prefix + defines_body


def generate_rfmprecompute_defines(
    rfm_precompute: ReferenceMetricPrecompute,
) -> Tuple[str, str]:
    """
    Generate the body and prefunctions for the rfm precompute defines.

    :param rfm_precompute: ReferenceMetricPrecompute object.
    :return: A tuple containing the C code for pre-functions and the main function body.
    """
    prefuncs = ""
    bodies = ""
    func_name = "defines"
    kernel_definitions = rfm_precompute.rfm_struct__define_kernel_dict

    for symbol, kernel_info in kernel_definitions.items():
        # Gather unique symbols from the expression, excluding coordinate symbols.
        unique_symbols = get_unique_expression_symbols_as_strings(
            kernel_info["expr"], exclude=[f"xx{j}" for j in range(3)]
        )

        # Build the kernel body, accessing necessary parameters from the params struct.
        kernel_body = "// Temporary parameters\n"
        for sym_str in unique_symbols:
            kernel_body += f"const REAL {sym_str} = params->{sym_str};\n"
        kernel_body += kernel_info["body"]

        # Prepare for C kernel generation.
        kernel_name = f"rfm_precompute_{func_name}__{symbol}"
        comments = f"Kernel to precompute metric quantity {symbol}."

        # Define arguments for CUDA device and host functions.
        arg_dict_cuda = {"rfmstruct": "rfm_struct *restrict"}
        arg_dict_host = {
            "params": "const params_struct *restrict",
            "rfmstruct": "rfm_struct *restrict",
        }
        for coord_var in kernel_info["coord"]:
            arg_dict_cuda[coord_var] = "const REAL *restrict"
            arg_dict_host[coord_var] = "const REAL *restrict"

        # Generate and append the C code for the kernel and its launcher.
        new_prefunc, new_body = generate_kernel_and_launch_code(
            kernel_name,
            kernel_body,
            arg_dict_cuda,
            arg_dict_host,
            par.parval_from_str("parallelization"),
            cfunc_type="static void",
            comments=comments,
            launchblock_with_braces=True,
            thread_tiling_macro_suffix="DEFAULT",
        )
        prefuncs += new_prefunc
        bodies += new_body

    return prefuncs, bodies


def register_CFunctions_rfm_precompute(
    set_of_CoordSystems: Set[str],
) -> None:
    """
    Register C functions for reference metric precomputed lookup arrays.

    :param set_of_CoordSystems: Set of coordinate systems to register the C functions.

    Doctest:
    >>> import nrpy.c_function as cfc
    >>> from nrpy.infrastructures.BHaH import rfm_precompute
    >>> from nrpy.reference_metric import unittest_CoordSystems
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.params as par
    >>> par.set_parval_from_str("fp_type", "float")
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in unittest_CoordSystems:
    ...       cfc.CFunction_dict.clear()
    ...       rfm_precompute.register_CFunctions_rfm_precompute({CoordSystem})
    ...       for rfm_base_function in ["malloc", "defines", "free"]:
    ...          generated_str = cfc.CFunction_dict[f'rfm_precompute_{rfm_base_function}__rfm__{CoordSystem}'].full_function
    ...          validation_desc = f"{rfm_base_function}__{parallelization}__{CoordSystem}".replace(" ", "_")
    ...          validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    Setting up reference_metric[SinhSymTP_rfm_precompute]...
    Setting up reference_metric[HoleySinhSpherical_rfm_precompute]...
    Setting up reference_metric[Cartesian_rfm_precompute]...
    Setting up reference_metric[SinhCylindricalv2n2_rfm_precompute]...
    """
    combined_BHaH_defines_list = []
    parallelization = par.parval_from_str("parallelization")
    is_cuda = parallelization == "cuda"

    for CoordSystem in set_of_CoordSystems:
        rfm_precompute = ReferenceMetricPrecompute(CoordSystem)

        # In CUDA, 'restrict' is a keyword and cannot be used with pointer-to-pointer members
        # like those in the rfm_struct.
        defines_list = [
            s.replace("restrict", "") if is_cuda else s
            for s in rfm_precompute.BHaH_defines_list
        ]
        combined_BHaH_defines_list.extend(defines_list)

        # Common setup for all registered C functions
        includes = ["BHaH_defines.h"]
        cfunc_type = "void"
        base_params = "const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct"

        defines_prefunc, _ = generate_rfmprecompute_defines(rfm_precompute)

        prefunc_dict = {"defines": defines_prefunc}
        c_functions_to_register = {
            "malloc": rfm_precompute.rfm_struct__malloc,
            "defines": rfm_precompute.rfm_struct__define,
            "free": rfm_precompute.rfm_struct__freemem,
        }

        for func_name, func_body in c_functions_to_register.items():
            desc = f"rfm_precompute_{func_name}: reference metric precomputed lookup arrays: {func_name}"
            name = f"rfm_precompute_{func_name}"
            params = base_params
            if func_name == "defines":
                params += ", REAL *restrict xx[3]"

            # Ensure the body is never an empty string to satisfy the CFunction constructor
            # and the doctest. A single space is a valid, empty C function body.
            final_body = func_body if func_body else " "

            cfc.register_CFunction(
                prefunc=prefunc_dict.get(func_name, ""),
                includes=includes,
                desc=desc,
                cfunc_type=cfunc_type,
                CoordSystem_for_wrapper_func=CoordSystem,
                name=name,
                params=params,
                include_CodeParameters_h=False,
                body=final_body,
            )

    # After processing all coordinate systems, define the rfm_struct
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
