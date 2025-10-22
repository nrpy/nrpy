# nrpy/infrastructures/BHaH/rfm_precompute.py
"""
C function registration for reference metric precomputation.

This module constructs and registers C functions that allocate, populate,
and free memory for grid-dependent reference metric quantities.
The generated C functions feature a runtime dispatch to support
both CPU and GPU (CUDA) execution.

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
from nrpy.infrastructures import BHaH


class ReferenceMetricPrecompute:
    """
    Manages the generation of C code for precomputing reference metric quantities.

    This class orchestrates the construction of C functions for allocating, defining,
    and freeing arrays of precomputed reference metric values. It handles both CPU
    and GPU (CUDA) code paths with a runtime dispatch mechanism.

    The residency contract for the generated code is as follows:
      - When params->is_host!=0, rfmstruct and xx[] are HOST pointers.
      - When params->is_host==0, rfmstruct and xx[] are DEVICE-resident or UVM pointers.

    Host-side SIMD reader strings are generated, which are rewritten to CUDA variants
    on CUDA-enabled builds.
    """

    def __init__(self, CoordSystem: str) -> None:
        """
        Initialize the ReferenceMetricPrecompute class.

        This constructor sets up the configuration, initializes internal properties,
        processes all precomputed expressions for the given coordinate system,
        and finalizes the generated C code snippets.

        :param CoordSystem: The coordinate system for which to generate precomputation code.
        """
        self._initialize_configuration(CoordSystem)
        self._initialize_properties()

        for symbol, expr in self._get_sorted_precomputed_expressions():
            if "_of_xx" in str(symbol):
                self._process_expression(symbol, expr)

        # Rewrite host-only reader intrinsics to CUDA variants when targeting CUDA.
        self._apply_intrinsics_mode()

        self._finalize_rfm_struct_define()

    def _initialize_configuration(self, CoordSystem: str) -> None:
        """
        Set up initial configuration based on the coordinate system and NRPy parameters.

        This method configures dimensions, detects the parallelization mode (e.g., CUDA),
        loads the appropriate reference metric data, and sets up symbolic parameters.

        :param CoordSystem: The name of the coordinate system.
        """
        self.N_DIMS = 3
        # Detect CUDA target from NRPy parameter "parallelization"
        self.parallelization_mode = par.parval_from_str("parallelization")
        self.is_cuda = self.parallelization_mode == "cuda"

        self.rfm = refmetric.reference_metric[CoordSystem + "_rfm_precompute"]
        fp_type = par.parval_from_str("fp_type")
        self.sp_type_alias = {sp_ast.real: ccg.fp_type_to_sympy_type[fp_type]}
        self.nxx_plus_2NGHOSTS_syms = [
            sp.Symbol(f"Nxx_plus_2NGHOSTS{dim_index}", real=True)
            for dim_index in range(self.N_DIMS)
        ]
        # Ensure params->is_host exists in params_struct (no-op if already registered)
        _ = par.register_CodeParameter(
            "bool", __name__, "is_host", True, commondata=False, add_to_parfile=False
        )

    def _initialize_properties(self) -> None:
        """
        Initialize internal lists and dictionaries for storing generated code and metadata.

        This method resets all containers used to accumulate C code for struct definitions,
        memory allocation/deallocation, and reader strings.
        """
        self.BHaH_defines_list: List[str] = []
        # Build malloc/free bodies from per-symbol host/device lines.
        self.rfm_struct__malloc_host_lines: List[str] = []
        self.rfm_struct__malloc_device_lines: List[str] = []
        self.rfm_struct__free_host_lines: List[str] = []
        self.rfm_struct__free_device_lines: List[str] = []

        self.rfm_struct__define = ""
        self.rfm_struct__define_kernel_dict: Dict[sp.Expr, Any] = {}

        # Host-only readers (later possibly rewritten from SIMD->CUDA)
        self.readvr_str = ["/* Host-only readers */\n"] * self.N_DIMS
        self.readvr_intrinsics_outer_str = [
            "/* Host-only SIMD readers */\n"
        ] * self.N_DIMS
        self.readvr_intrinsics_inner_str = [
            "/* Host-only SIMD readers */\n"
        ] * self.N_DIMS

    def _get_sorted_precomputed_expressions(self) -> List[Tuple[sp.Expr, sp.Expr]]:
        """
        Retrieve and sort the precomputed reference metric expressions.

        This method extracts the coordinate-dependent free variables and their corresponding
        expressions from the reference metric object, then sorts them alphabetically.

        :return: A sorted list of (symbol, expression) tuples.
        """
        if not self.rfm.freevars_uniq_xx_indep:
            return []
        zipped_pairs = zip(self.rfm.freevars_uniq_xx_indep, self.rfm.freevars_uniq_vals)
        return sorted(zipped_pairs, key=lambda pair: str(pair[0]))

    def _process_expression(self, symbol: sp.Expr, expr: sp.Expr) -> None:
        """
        Process a single symbolic expression to generate corresponding C code.

        This function determines the coordinate dependencies of the expression (1D or 2D),
        generates memory management code, and dispatches to the appropriate handler
        to generate the computational logic.

        :param symbol: The SymPy symbol for the precomputed quantity.
        :param expr: The SymPy expression defining the precomputed quantity.
        :raises RuntimeError: If the expression has an unsupported dependency pattern.
        """
        free_symbols_set = list(expr.free_symbols)
        dependencies = [self.rfm.xx[i] in free_symbols_set for i in range(self.N_DIMS)]
        self._add_memory_management_code(symbol, dependencies)

        has_1d_dep = sum(dependencies) == 1
        has_2d_dep_xx0_xx1 = dependencies[0] and dependencies[1] and not dependencies[2]

        if has_1d_dep:
            dep_axis_index = dependencies.index(True)
            self._handle_1d_dependency(symbol, expr, dep_axis_index)
        elif has_2d_dep_xx0_xx1:
            self._handle_2d_dependency(symbol, expr)
        else:
            raise RuntimeError(f"Unsupported dependency for {symbol}: {expr}")

    def _add_memory_management_code(
        self, symbol: sp.Expr, dependencies: List[bool]
    ) -> None:
        """
        Generate C code lines for memory allocation and deallocation.

        This method constructs C code snippets for both host (`malloc`/`free`) and device
        (`cudaMalloc`/`cudaFree`) memory management for a given precomputed quantity.

        :param symbol: The symbol representing the array to be allocated.
        :param dependencies: A list of booleans indicating coordinate dependencies (xx0, xx1, xx2).
        """
        # No restrict on struct members; project provides its own 'restrict' define.
        self.BHaH_defines_list.append(f"REAL * {symbol};")
        size_terms = [
            f"(size_t)params->Nxx_plus_2NGHOSTS{dim_index}"
            for dim_index, is_dependent in enumerate(dependencies)
            if is_dependent
        ]
        size_expr = " * ".join(size_terms) if size_terms else "(size_t)1"
        # Host allocation/free use existing BHaH macros:
        self.rfm_struct__malloc_host_lines.append(
            f"BHAH_MALLOC__PtrMember(rfmstruct, {symbol}, sizeof(REAL) * ({size_expr}));"
        )
        self.rfm_struct__free_host_lines.append(
            f"BHAH_FREE__PtrMember(rfmstruct, {symbol});"
        )
        # Device allocation/free are defined via local (prefunc) CUDA helpers:
        self.rfm_struct__malloc_device_lines.append(
            f"BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, {symbol}, sizeof(REAL) * ({size_expr}));"
        )
        self.rfm_struct__free_device_lines.append(
            f"BHAH_FREE_DEVICE__PtrMember(rfmstruct, {symbol});"
        )

    def _handle_1d_dependency(
        self, symbol: sp.Expr, expr: sp.Expr, dep_axis_index: int
    ) -> None:
        """
        Generate code for a quantity with a 1D coordinate dependency.

        This method populates the kernel dictionary with metadata for the 1D computation
        and generates the corresponding host-only reader strings for scalar and SIMD access.

        :param symbol: The SymPy symbol for the precomputed quantity.
        :param expr: The SymPy expression defining the quantity.
        :param dep_axis_index: The index (0, 1, or 2) of the coordinate dependency.
        """
        self.rfm_struct__define_kernel_dict[symbol] = {
            "kind": "1d",
            "dep_axis_index": dep_axis_index,
            "expr_cc": sp.ccode(expr, type_aliases=self.sp_type_alias),
            "coord": [f"x{dep_axis_index}"],
            "unique_params": get_unique_expression_symbols_as_strings(
                expr, exclude=[f"xx{j}" for j in range(self.N_DIMS)]
            ),
        }
        # Host-only readers assume enclosing loops define i{axis} variables.
        self.readvr_str[
            dep_axis_index
        ] += f"MAYBE_UNUSED const REAL {symbol} = rfmstruct->{symbol}[i{dep_axis_index}];\n"
        self.readvr_intrinsics_outer_str[dep_axis_index] += (
            f"const double NOSIMD{symbol} = rfmstruct->{symbol}[i{dep_axis_index}]; "
            f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol} = ConstSIMD(NOSIMD{symbol});\n"
        )
        self.readvr_intrinsics_inner_str[
            dep_axis_index
        ] += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol} = ReadSIMD(&rfmstruct->{symbol}[i{dep_axis_index}]);\n"

    def _handle_2d_dependency(self, symbol: sp.Expr, expr: sp.Expr) -> None:
        """
        Generate code for a quantity with a 2D coordinate dependency (xx0, xx1).

        This method populates the kernel dictionary with metadata for the 2D computation
        and generates the corresponding host-only reader strings.

        :param symbol: The SymPy symbol for the precomputed quantity.
        :param expr: The SymPy expression defining the quantity.
        """
        self.rfm_struct__define_kernel_dict[symbol] = {
            "kind": "2d",
            "expr_cc": sp.ccode(expr, type_aliases=self.sp_type_alias),
            "coord": ["x0", "x1"],
            "unique_params": get_unique_expression_symbols_as_strings(
                expr, exclude=[f"xx{j}" for j in range(self.N_DIMS)]
            ),
        }
        mem_index_expr = "(i0 + (size_t)N0 * (size_t)i1)"  # host-only readers
        self.readvr_str[
            0
        ] += f"MAYBE_UNUSED const REAL {symbol} = rfmstruct->{symbol}[{mem_index_expr}];\n"
        self.readvr_intrinsics_outer_str[0] += (
            f"const double NOSIMD{symbol} = rfmstruct->{symbol}[{mem_index_expr}]; "
            f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol} = ConstSIMD(NOSIMD{symbol});\n"
        )
        self.readvr_intrinsics_inner_str[
            0
        ] += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol} = ReadSIMD(&rfmstruct->{symbol}[{mem_index_expr}]);\n"

    def _apply_intrinsics_mode(self) -> None:
        """
        Rewrite host-side reader intrinsics from SIMD to CUDA variants on CUDA builds.

        This function inspects the parallelization mode and, if it is 'cuda', replaces
        all instances of "SIMD" with "CUDA" in the generated reader strings. These strings
        are only used for host-side code; device kernels perform direct memory loads.
        """
        if not self.is_cuda:
            return
        self.readvr_intrinsics_outer_str = [
            s.replace("SIMD", "CUDA") for s in self.readvr_intrinsics_outer_str
        ]
        self.readvr_intrinsics_inner_str = [
            s.replace("SIMD", "CUDA") for s in self.readvr_intrinsics_inner_str
        ]

    def _finalize_rfm_struct_define(self) -> None:
        """
        Assemble the final C code for the body of the `rfm_precompute_defines` function.

        This method combines boilerplate C code (e.g., aliasing coordinate arrays) with
        the dynamically generated loops and kernel launches for all precomputed quantities.
        """
        body_prefix_lines: List[str] = []
        for dim_index in range(3):
            body_prefix_lines.append(
                f"MAYBE_UNUSED const REAL * restrict x{dim_index} = xx[{dim_index}];"
            )
        _, defines_body = generate_rfmprecompute_defines(self)
        self.rfm_struct__define = "\n".join(body_prefix_lines) + "\n" + defines_body


def _emit_kernel_param_copies(
    unique_params: List[str], kernel_param_pointer: str, indent_spaces: int = 4
) -> str:
    """
    Generate C code to copy parameters from a struct pointer to local variables.

    This helper creates a block of C code that declares and initializes local `const REAL`
    variables from the members of a given struct pointer (e.g., `kparams->param_name`).

    :param unique_params: A list of parameter names (strings) to be copied.
    :param kernel_param_pointer: The name of the C struct pointer (e.g., "kparams").
    :param indent_spaces: The number of spaces to use for indentation.
    :return: A string containing the formatted C code block.
    """
    if not unique_params:
        return ""
    indent = " " * indent_spaces
    lines = [
        f"{indent}const REAL {param_name} = {kernel_param_pointer}->{param_name};"
        for param_name in unique_params
    ]
    return "\n".join(lines) + "\n"


def _emit_host_param_copies(unique_params: List[str], indent_spaces: int = 4) -> str:
    """
    Generate C code to copy parameters from the `params` struct to local variables.

    This helper creates a block of C code that declares and initializes local `const REAL`
    variables from the members of the main `params` struct (e.g., `params->param_name`).

    :param unique_params: A list of parameter names (strings) to be copied.
    :param indent_spaces: The number of spaces to use for indentation.
    :return: A string containing the formatted C code block.
    """
    if not unique_params:
        return ""
    indent = " " * indent_spaces
    lines = [
        f"{indent}const REAL {param_name} = params->{param_name};"
        for param_name in unique_params
    ]
    return "\n".join(lines) + "\n"


def _prefunc_common_guards() -> str:
    """
    Generate preprocessor guards for accelerator availability.

    This function returns a C preprocessor block that defines `BHAH_ASSERT_VALID_ACCELERATOR`.
    In non-CUDA builds, this macro will cause the program to abort if a GPU path is
    erroneously requested. In CUDA builds, it is a no-op.

    :return: A string containing the C preprocessor definitions.
    """
    return (
        "/* Local safety guards */\n"
        "#ifndef __CUDACC__\n"
        '#define BHAH_ASSERT_VALID_ACCELERATOR(p) do { if ((p)->is_host==0) { fprintf(stderr, "Error: GPU path requested but CUDA is not enabled.\\n"); abort(); } } while(0)\n'
        "#else\n"
        "#define BHAH_ASSERT_VALID_ACCELERATOR(p) do { (void)(p); } while(0)\n"
        "#endif\n"
    )


def _prefunc_cuda_alloc_helpers() -> str:
    """
    Generate preprocessor definitions for CUDA memory management macros.

    This function returns a C preprocessor block that defines `BHAH_MALLOC_DEVICE__PtrMember`
    and `BHAH_FREE_DEVICE__PtrMember`. In CUDA builds, these macros wrap `cudaMalloc` and
    `cudaFree` with error checking. In non-CUDA builds, they are defined as no-ops.

    :return: A string containing the C preprocessor definitions.
    """
    return (
        "#ifdef __CUDACC__\n"
        "#  include <cuda_runtime.h>\n"
        "#  ifndef BHAH_MALLOC_DEVICE__PtrMember\n"
        "#    define BHAH_MALLOC_DEVICE__PtrMember(structptr, member, nbytes) \\\n"
        "       do { cudaError_t __e = cudaMalloc((void**)&((structptr)->member), (size_t)(nbytes)); \\\n"
        '            if (__e != cudaSuccess) { fprintf(stderr, "cudaMalloc failed: %s\\n", cudaGetErrorString(__e)); abort(); } } while(0)\n'
        "#  endif\n"
        "#  ifndef BHAH_FREE_DEVICE__PtrMember\n"
        "#    define BHAH_FREE_DEVICE__PtrMember(structptr, member) \\\n"
        "       do { if ((structptr)->member) { cudaError_t __e2 = cudaFree((structptr)->member); (structptr)->member = NULL; \\\n"
        '            if (__e2 != cudaSuccess) { fprintf(stderr, "cudaFree failed: %s\\n", cudaGetErrorString(__e2)); abort(); } } } while(0)\n'
        "#  endif\n"
        "#else\n"
        "#  ifndef BHAH_MALLOC_DEVICE__PtrMember\n"
        "#    define BHAH_MALLOC_DEVICE__PtrMember(structptr, member, nbytes) do { (void)(structptr); (void)(member); (void)(nbytes); } while(0)\n"
        "#  endif\n"
        "#  ifndef BHAH_FREE_DEVICE__PtrMember\n"
        "#    define BHAH_FREE_DEVICE__PtrMember(structptr, member) do { (void)(structptr); (void)(member); } while(0)\n"
        "#  endif\n"
        "#endif\n"
    )


def _prefunc_launch_defaults() -> str:
    """
    Generate preprocessor definitions for default CUDA kernel launch configurations.

    This function returns a C preprocessor block defining helper functions for default
    1D and 2D grid and block dimensions for CUDA kernels. These are only defined in
    CUDA-enabled builds.

    :return: A string containing the C preprocessor definitions.
    """
    return (
        "#ifdef __CUDACC__\n"
        "static inline dim3 BHAH_DEFAULT_1D_BLOCK(void) { return dim3(256,1,1); }\n"
        "static inline dim3 BHAH_DEFAULT_1D_GRID(size_t N) { unsigned gx = (unsigned)((N + 256u - 1u)/256u); return dim3(gx?gx:1u,1u,1u); }\n"
        "static inline dim3 BHAH_DEFAULT_2D_BLOCK(void) { return dim3(32,8,1); }\n"
        "static inline dim3 BHAH_DEFAULT_2D_GRID(size_t N0, size_t N1) { unsigned gx = (unsigned)((N0 + 32u - 1u)/32u); unsigned gy = (unsigned)((N1 + 8u - 1u)/8u); return dim3(gx?gx:1u, gy?gy:1u,1u); }\n"
        "#endif\n"
    )


def generate_rfmprecompute_defines(
    rfm_precompute: ReferenceMetricPrecompute,
) -> Tuple[str, str]:
    """
    Generate the C code for CUDA kernels and the body of the `defines` function.

    This function iterates through the precomputed quantities stored in the
    `rfm_precompute` object and generates the corresponding C code for both the
    CUDA kernels (as a pre-function string) and the main function body, which
    contains the host-side loops and runtime dispatch logic for launching GPU kernels.

    :param rfm_precompute: An instance of the ReferenceMetricPrecompute class containing the processed expressions.
    :raises RuntimeError: If an unknown dependency kind is encountered.
    :return: A tuple containing two strings: the C code for the CUDA kernels (`prefuncs_kernels`)
             and the C code for the main `defines` function body (`body`).
    """
    prefuncs_kernels = ""
    body = ""

    for symbol, info in rfm_precompute.rfm_struct__define_kernel_dict.items():
        symbol_name = str(symbol)
        unique_params = info.get("unique_params", [])

        if info["kind"] == "1d":
            dep_axis_index = info["dep_axis_index"]
            expr_cc = info["expr_cc"]

            # ---- CUDA kernel (prefuncs_kernels) ----
            kernel_param_lines = _emit_kernel_param_copies(
                unique_params, "kparams", indent_spaces=4
            )
            prefuncs_kernels += (
                "#ifdef __CUDACC__\n"
                f"__global__ static void rfm_precompute_defines__{symbol_name}(\n"
                "    const params_struct * restrict kparams,\n"
                "    rfm_struct * restrict d_rfm,\n"
                f"    const REAL * restrict dx{dep_axis_index}) {{\n"
                f"  const size_t N = (size_t)kparams->Nxx_plus_2NGHOSTS{dep_axis_index};\n"
                "  size_t thread_linear_index = (size_t)threadIdx.x + (size_t)blockIdx.x * (size_t)blockDim.x;\n"
                "  size_t thread_linear_stride = (size_t)blockDim.x * (size_t)gridDim.x;\n"
                f"  for (size_t i{dep_axis_index} = thread_linear_index; i{dep_axis_index} < N; "
                f"i{dep_axis_index} += thread_linear_stride) {{\n"
                f"    const REAL xx{dep_axis_index} = dx{dep_axis_index}[i{dep_axis_index}];\n"
                f"{kernel_param_lines}"
                f"    d_rfm->{symbol_name}[i{dep_axis_index}] = {expr_cc};\n"
                "  }\n"
                "}\n"
                "#endif\n\n"
            )

            # ---- Host loop + runtime GPU branch ----
            host_param_lines = _emit_host_param_copies(unique_params, indent_spaces=4)
            body += (
                f"/* {symbol_name}: 1D precompute */\n"
                "if (params->is_host) {\n"
                "  {\n"
                f"{host_param_lines}"
                f"    const size_t N = (size_t)params->Nxx_plus_2NGHOSTS{dep_axis_index};\n"
                f"    for (size_t i{dep_axis_index}=0; i{dep_axis_index}<N; i{dep_axis_index}++) {{\n"
                f"      const REAL xx{dep_axis_index} = x{dep_axis_index}[i{dep_axis_index}];\n"
                f"      rfmstruct->{symbol_name}[i{dep_axis_index}] = {expr_cc};\n"
                "    }\n"
                "  }\n"
                "} else {\n"
                "  BHAH_ASSERT_VALID_ACCELERATOR(params);\n"
                "  IFCUDARUN({\n"
                f"    const size_t N = (size_t)params->Nxx_plus_2NGHOSTS{dep_axis_index};\n"
                "    dim3 block = BHAH_DEFAULT_1D_BLOCK();\n"
                "    dim3 grid  = BHAH_DEFAULT_1D_GRID(N);\n"
                f"    rfm_precompute_defines__{symbol_name}<<<grid, block>>>(params, rfmstruct, x{dep_axis_index});\n"
                "  });\n"
                "}\n\n"
            )

        elif info["kind"] == "2d":
            expr_cc = info["expr_cc"]

            kernel_param_lines = _emit_kernel_param_copies(
                unique_params, "kparams", indent_spaces=6
            )
            prefuncs_kernels += (
                "#ifdef __CUDACC__\n"
                f"__global__ static void rfm_precompute_defines__{symbol_name}(\n"
                "    const params_struct * restrict kparams,\n"
                "    rfm_struct * restrict d_rfm,\n"
                "    const REAL * restrict dx0,\n"
                "    const REAL * restrict dx1) {\n"
                "  const size_t N0 = (size_t)kparams->Nxx_plus_2NGHOSTS0;\n"
                "  const size_t N1 = (size_t)kparams->Nxx_plus_2NGHOSTS1;\n"
                "  size_t idx0 = (size_t)threadIdx.x + (size_t)blockIdx.x * (size_t)blockDim.x;\n"
                "  size_t idx1 = (size_t)threadIdx.y + (size_t)blockIdx.y * (size_t)blockDim.y;\n"
                "  size_t stride0 = (size_t)blockDim.x * (size_t)gridDim.x;\n"
                "  size_t stride1 = (size_t)blockDim.y * (size_t)gridDim.y;\n"
                "  for (size_t i1=idx1; i1<N1; i1+=stride1) {\n"
                "    for (size_t i0=idx0; i0<N0; i0+=stride0) {\n"
                "      const REAL xx0 = dx0[i0];\n"
                "      const REAL xx1 = dx1[i1];\n"
                f"{kernel_param_lines}"
                "      const size_t flat_idx = i0 + (size_t)N0 * i1;\n"
                f"      d_rfm->{symbol_name}[flat_idx] = {expr_cc};\n"
                "    }\n"
                "  }\n"
                "}\n"
                "#endif\n\n"
            )

            host_param_lines = _emit_host_param_copies(unique_params, indent_spaces=4)
            body += (
                f"/* {symbol_name}: 2D (xx0,xx1) precompute */\n"
                "if (params->is_host) {\n"
                "  {\n"
                f"{host_param_lines}"
                "    const size_t N0 = (size_t)params->Nxx_plus_2NGHOSTS0;\n"
                "    const size_t N1 = (size_t)params->Nxx_plus_2NGHOSTS1;\n"
                "    for (size_t i1=0; i1<N1; i1++) {\n"
                "      for (size_t i0=0; i0<N0; i0++) {\n"
                "        const REAL xx0 = x0[i0];\n"
                "        const REAL xx1 = x1[i1];\n"
                "        const size_t flat_idx = i0 + (size_t)N0 * i1;\n"
                f"        rfmstruct->{symbol_name}[flat_idx] = {expr_cc};\n"
                "      }\n"
                "    }\n"
                "  }\n"
                "} else {\n"
                "  BHAH_ASSERT_VALID_ACCELERATOR(params);\n"
                "  IFCUDARUN({\n"
                "    const size_t N0 = (size_t)params->Nxx_plus_2NGHOSTS0;\n"
                "    const size_t N1 = (size_t)params->Nxx_plus_2NGHOSTS1;\n"
                "    dim3 block = BHAH_DEFAULT_2D_BLOCK();\n"
                "    dim3 grid  = BHAH_DEFAULT_2D_GRID(N0, N1);\n"
                f"    rfm_precompute_defines__{symbol_name}<<<grid, block>>>(params, rfmstruct, x0, x1);\n"
                "  });\n"
                "}\n\n"
            )

        else:
            raise RuntimeError("Internal: unknown dependency kind")

    return prefuncs_kernels, body


def register_CFunctions_rfm_precompute(set_of_CoordSystems: Set[str]) -> None:
    """
    Construct and register C functions for allocating, populating, and freeing reference metric precomputed arrays.

    This function iterates through a set of coordinate systems and, for each,
    generates and registers three C functions:
    1. rfm_precompute_malloc: Allocates memory for precomputed arrays.
    2. rfm_precompute_defines: Populates the arrays with coordinate-dependent values.
    3. rfm_precompute_free: Deallocates the memory.

    The generated functions include a runtime dispatch mechanism, allowing them to
    execute on either the CPU (host) or a GPU (device) in CUDA-enabled builds,
    based on the `params->is_host` flag.

    :param set_of_CoordSystems: A set of coordinate system names for which to generate the C functions.

    Doctests:
    TBD
    """
    combined_BHaH_defines_list: List[str] = []

    for CoordSystem in set_of_CoordSystems:
        rfm_precompute = ReferenceMetricPrecompute(CoordSystem)

        combined_BHaH_defines_list.extend(list(rfm_precompute.BHaH_defines_list))

        includes = ["BHaH_defines.h"]
        cfunc_type = "void"
        base_params = (
            "const commondata_struct * restrict commondata, "
            "const params_struct * restrict params, "
            "rfm_struct * restrict rfmstruct"
        )

        kernels_prefunc, _defines_body = generate_rfmprecompute_defines(rfm_precompute)
        # Attach guards and launch defaults exactly once:
        defines_prefunc = (
            _prefunc_common_guards() + _prefunc_launch_defaults() + kernels_prefunc
        )

        # Compose malloc/free bodies from accumulated host/device lines:
        malloc_body = (
            "/* rfm_precompute_malloc: allocate rfmstruct arrays on host or device */\n"
            "if (params->is_host) {\n"
            + "".join(
                f"  {line}\n" for line in rfm_precompute.rfm_struct__malloc_host_lines
            )
            + "} else {\n"
            "  BHAH_ASSERT_VALID_ACCELERATOR(params);\n"
            "  IFCUDARUN({\n"
            + "".join(
                f"    {line}\n"
                for line in rfm_precompute.rfm_struct__malloc_device_lines
            )
            + "  });\n"
            "}\n"
        )

        free_body = (
            "/* rfm_precompute_free: free rfmstruct arrays from host or device */\n"
            "if (params->is_host) {\n"
            + "".join(
                f"  {line}\n" for line in rfm_precompute.rfm_struct__free_host_lines
            )
            + "} else {\n"
            "  BHAH_ASSERT_VALID_ACCELERATOR(params);\n"
            "  IFCUDARUN({\n"
            + "".join(
                f"    {line}\n" for line in rfm_precompute.rfm_struct__free_device_lines
            )
            + "  });\n"
            "}\n"
        )

        # Prefuncs for malloc/free define CUDA alloc helpers locally (as requested).
        malloc_prefunc = _prefunc_common_guards() + _prefunc_cuda_alloc_helpers()
        free_prefunc = _prefunc_common_guards() + _prefunc_cuda_alloc_helpers()

        prefunc_dict = {
            "malloc": malloc_prefunc,
            "defines": defines_prefunc,
            "free": free_prefunc,
        }
        c_functions_to_register = {
            "malloc": malloc_body,
            "defines": rfm_precompute.rfm_struct__define,
            "free": free_body,
        }
        doxygen_descs = {
            "malloc": """
 * @file rfm_precompute_malloc.c
 * @brief Allocates memory for precomputed reference metric quantities.
 *
 * This function allocates memory for arrays within the rfm_struct, which store
 * precomputed quantities of the reference metric that depend on grid coordinates.
 *
 * It adheres to the following residency contract:
 *   - When params->is_host==true, it performs HOST memory allocation for CPU execution.
 *   - When params->is_host==false, it performs DEVICE memory allocation for GPU execution.
 * In CPU-only builds, requests for device allocation will cause an abort.
 *
 * @param[in]  commondata  Pointer to global simulation metadata (unused).
 * @param[in]  params      Pointer to grid parameters, containing is_host for dispatch and grid dimensions.
 * @param[out] rfmstruct   Pointer to the struct where pointers to allocated memory will be stored.
 *
 * @return void.
 """,
            "defines": """
 * @file rfm_precompute_defines.c
 * @brief Computes and populates arrays with precomputed reference metric quantities.
 *
 * This function computes reference metric quantities that depend on grid coordinate
 * values (e.g., sin(xx[1])) and stores them in the pre-allocated arrays
 * within the rfm_struct.
 *
 * It adheres to the following residency contract:
 *   - When params->is_host==true, rfmstruct and xx[] are HOST pointers; a CPU path with loops is executed.
 *   - When params->is_host==false, rfmstruct and xx[] are DEVICE or UVM pointers; a GPU path with CUDA kernels is executed.
 * In CPU-only builds, GPU path requests will cause an abort.
 *
 * @param[in]  commondata  Pointer to global simulation metadata (unused).
 * @param[in]  params      Pointer to grid parameters, including is_host for dispatch and grid dimensions.
 * @param[out] rfmstruct   Pointer to the struct containing the arrays to be populated.
 * @param[in]  xx          Array of pointers to the 1D grid coordinate arrays.
 *
 * @return void.
 """,
            "free": """
 * @file rfm_precompute_free.c
 * @brief Frees memory for precomputed reference metric quantities.
 *
 * This function frees the memory for arrays within the rfm_struct that was
 * previously allocated by rfm_precompute_malloc().
 *
 * It adheres to the following residency contract:
 *   - When params->is_host==true, it performs HOST memory deallocation.
 *   - When params->is_host==false, it performs DEVICE memory deallocation.
 * In CPU-only builds, requests for device deallocation will cause an abort.
 *
 * @param[in]  commondata  Pointer to global simulation metadata (unused).
 * @param[in]  params      Pointer to grid parameters, including is_host for dispatch.
 * @param[out] rfmstruct   Pointer to the struct containing pointers to the memory to be freed.
 *
 * @return void.
 """,
        }

        for func_name, func_body in c_functions_to_register.items():
            function_desc = doxygen_descs[func_name]
            function_name = f"rfm_precompute_{func_name}"
            params_sig = base_params
            if func_name == "defines":
                params_sig += ", REAL * restrict xx[3]"

            residency_contract_comment = (
                "/* Residency contract:\n"
                "   - params->is_host==true: rfmstruct, xx[] are HOST pointers; CPU path executes.\n"
                "   - params->is_host==false: rfmstruct, xx[] are DEVICE or UVM pointers; GPU path executes.\n"
                "   - In CPU-only builds, GPU requests abort via BHAH_ASSERT_VALID_ACCELERATOR.\n"
                "*/\n"
            )

            final_body = (
                (residency_contract_comment + func_body)
                if func_body
                else residency_contract_comment + " "
            )

            cfc.register_CFunction(
                prefunc=prefunc_dict.get(func_name, ""),
                includes=includes,
                desc=function_desc,
                cfunc_type=cfunc_type,
                CoordSystem_for_wrapper_func=CoordSystem,
                name=function_name,
                params=params_sig,
                include_CodeParameters_h=False,
                body=final_body,
            )

    rfm_struct_typedef = "typedef struct __rfmstruct__ {\n"
    rfm_struct_typedef += "\n".join(sorted(superfast_uniq(combined_BHaH_defines_list)))
    rfm_struct_typedef += "\n} rfm_struct;\n"
    BHaH.BHaH_defines_h.register_BHaH_defines("reference_metric", rfm_struct_typedef)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
