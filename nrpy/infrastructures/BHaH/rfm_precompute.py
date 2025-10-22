# nrpy/infrastructures/BHaH/rfm_precompute.py
"""
C function registration for reference metric precomputation.

Generates & registers C functions that allocate, populate, and free
grid-dependent reference-metric arrays. Runtime dispatch supports
both CPU and CUDA, with the host path using plain loops and the
device path launching CUDA kernels.

Design choices to match NRPy/BHaH style:
  * Python structure mirrors the C we generate.
  * Avoid fragile string concatenation; prefer clear, block-like f-strings.
  * CUDA launch config uses classic BHaH macros:
       BHAH_THREADS_IN_X_DIR_DEFAULT / _Y_ / _Z_
     with explicit dim3 math (no bespoke "DEFAULT_*" helpers).

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


# --------------------------------------------------------------------------
# Orchestrator: build per-symbol allocation/compute/free & CUDA kernels
# --------------------------------------------------------------------------
class ReferenceMetricPrecompute:
    """
    Manages the generation of C code for precomputing reference metric quantities.

    This class orchestrates the construction of C functions for allocating, defining,
    and freeing arrays of precomputed reference metric values. It handles both CPU
    and GPU (CUDA) code paths with a runtime dispatch mechanism.

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

    # --------------------------- setup ---------------------------
    def _initialize_configuration(self, CoordSystem: str) -> None:
        """
        Set up initial configuration based on the coordinate system and NRPy parameters.

        This method configures dimensions, detects the parallelization mode (e.g., CUDA),
        loads the appropriate reference metric data, and sets up symbolic parameters.

        :param CoordSystem: The name of the coordinate system.
        """
        self.N_DIMS = 3
        self.parallelization_mode = par.parval_from_str("parallelization")
        self.is_cuda = self.parallelization_mode == "cuda"

        self.rfm = refmetric.reference_metric[CoordSystem + "_rfm_precompute"]
        fp_type = par.parval_from_str("fp_type")
        self.sp_type_alias = {sp_ast.real: ccg.fp_type_to_sympy_type[fp_type]}

        self.nxx_plus_2NGHOSTS_syms = [
            sp.Symbol(f"Nxx_plus_2NGHOSTS{d}", real=True) for d in range(self.N_DIMS)
        ]
        # Ensure params->is_host exists (harmless if already present).
        _ = par.register_CodeParameter(
            "bool", __name__, "is_host", True, commondata=False, add_to_parfile=False
        )

    def _initialize_properties(self) -> None:
        """
        Set up initial configuration based on the coordinate system and NRPy parameters.

        This method configures dimensions, detects the parallelization mode (e.g., CUDA),
        loads the appropriate reference metric data, and sets up symbolic parameters.
        """
        self.BHaH_defines_list: List[str] = []

        # Per-member alloc/free lines for host & device:
        self.rfm_struct__malloc_host_lines: List[str] = []
        self.rfm_struct__malloc_device_lines: List[str] = []
        self.rfm_struct__free_host_lines: List[str] = []
        self.rfm_struct__free_device_lines: List[str] = []

        # Main "defines" (compute) body & per-symbol kernel metadata:
        self.rfm_struct__define = ""
        self.rfm_struct__define_kernel_dict: Dict[sp.Expr, Any] = {}

        # Host-only readers (later rewritten to CUDA variants if needed):
        self.readvr_str = ["/* Host-only readers */\n"] * self.N_DIMS
        self.readvr_intrinsics_outer_str = [
            "/* Host-only SIMD readers */\n"
        ] * self.N_DIMS
        self.readvr_intrinsics_inner_str = [
            "/* Host-only SIMD readers */\n"
        ] * self.N_DIMS

    # --------------------------- discovery ---------------------------
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

    # --------------------------- per-symbol processing ---------------------------
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

    # --------------------------- alloc/free snippets ---------------------------
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
        # No restrict on struct members; BHaH provides its own 'restrict' define.
        self.BHaH_defines_list.append(f"REAL * {symbol};")

        size_terms = [
            f"(size_t)params->Nxx_plus_2NGHOSTS{dim}"
            for dim, dep in enumerate(dependencies)
            if dep
        ]
        size_expr = " * ".join(size_terms) if size_terms else "(size_t)1"

        # Host allocation/free use BHaH macros:
        self.rfm_struct__malloc_host_lines.append(
            f"BHAH_MALLOC__PtrMember(rfmstruct, {symbol}, sizeof(REAL) * ({size_expr}));"
        )
        self.rfm_struct__free_host_lines.append(
            f"BHAH_FREE__PtrMember(rfmstruct, {symbol});"
        )

        # Device allocation/free (local CUDA helpers guarded in prefuncs):
        self.rfm_struct__malloc_device_lines.append(
            f"BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, {symbol}, sizeof(REAL) * ({size_expr}));"
        )
        self.rfm_struct__free_device_lines.append(
            f"BHAH_FREE_DEVICE__PtrMember(rfmstruct, {symbol});"
        )

    # --------------------------- dependency kinds ---------------------------
    # --------------------------- small helper ---------------------------
    def _ccg_rhs(self, expr: sp.Expr) -> str:
        """
        Generate a single C expression (RHS only) for `expr` using NRPy's c_codegen().
        We request one simple assignment, then strip off the LHS and trailing semicolon.
        """
        # Ask c_codegen for a single, simple line: "REAL __rhs = <expr>;"
        code = ccg.c_codegen(
            expr,
            "REAL __rhs",
            include_braces=False,
            verbose=False,
            enable_cse=False,
        )
        # Extract the RHS of the last non-empty line containing an '='
        line = next(
            l
            for l in (s.strip() for s in code.splitlines())
            if ("=" in l and l.endswith(";"))
        )
        return line.split("=", 1)[1].strip().rstrip(";")

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
            "expr_cc": self._ccg_rhs(expr),
            "coord": [f"x{dep_axis_index}"],
            "unique_params": get_unique_expression_symbols_as_strings(
                expr, exclude=[f"xx{j}" for j in range(self.N_DIMS)]
            ),
        }

        # Host-only readers assume outer loops define i{axis}.
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
            "expr_cc": self._ccg_rhs(expr),
            "coord": ["x0", "x1"],
            "unique_params": get_unique_expression_symbols_as_strings(
                expr, exclude=[f"xx{j}" for j in range(self.N_DIMS)]
            ),
        }
        flat = "(i0 + (size_t)N0 * (size_t)i1)"
        self.readvr_str[
            0
        ] += f"MAYBE_UNUSED const REAL {symbol} = rfmstruct->{symbol}[{flat}];\n"
        self.readvr_intrinsics_outer_str[0] += (
            f"const double NOSIMD{symbol} = rfmstruct->{symbol}[{flat}]; "
            f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol} = ConstSIMD(NOSIMD{symbol});\n"
        )
        self.readvr_intrinsics_inner_str[
            0
        ] += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol} = ReadSIMD(&rfmstruct->{symbol}[{flat}]);\n"

    # --------------------------- SIMD<->CUDA reader replacements ---------------------------
    def _apply_intrinsics_mode(self) -> None:
        """
        Replace host-side reader intrinsics from SIMD to CUDA variants on CUDA builds.

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

    # --------------------------- finalize defines body ---------------------------
    def _finalize_rfm_struct_define(self) -> None:
        """
        Assemble the final C code for the body of the `rfm_precompute_defines` function.

        This method combines boilerplate C code (e.g., aliasing coordinate arrays) with
        the dynamically generated loops and kernel launches for all precomputed quantities.
        """
        alias_lines = [
            f"MAYBE_UNUSED const REAL * restrict x{d} = xx[{d}];" for d in range(3)
        ]
        _, defines_body = generate_rfmprecompute_defines(self)
        self.rfm_struct__define = "\n".join(alias_lines) + "\n" + defines_body


# --------------------------------------------------------------------------
# Small C helpers: parameter copies into locals (host/kernels)
# --------------------------------------------------------------------------
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
    return "".join(
        f"{indent}const REAL {p} = {kernel_param_pointer}->{p};\n"
        for p in unique_params
    )


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
    return "".join(f"{indent}const REAL {p} = params->{p};\n" for p in unique_params)


# --------------------------------------------------------------------------
# CUDA alloc helpers (kept local to this translation unit for malloc/free)
# --------------------------------------------------------------------------
def _prefunc_cuda_alloc_helpers() -> str:
    """
    Generate preprocessor definitions for CUDA memory management macros.

    This function returns a C preprocessor block that defines `BHAH_MALLOC_DEVICE__PtrMember`
    and `BHAH_FREE_DEVICE__PtrMember`. In CUDA builds, these macros wrap `cudaMalloc` and
    `cudaFree` with error checking. In non-CUDA builds, they are defined as no-ops.

    :return: A string containing the C preprocessor definitions.
    """
    return r"""#ifdef __CUDACC__
#  include <cuda_runtime.h>
#  ifndef BHAH_MALLOC_DEVICE__PtrMember
#    define BHAH_MALLOC_DEVICE__PtrMember(structptr, member, nbytes) \
       do { cudaError_t __e = cudaMalloc((void**)&((structptr)->member), (size_t)(nbytes)); \
            if (__e != cudaSuccess) { fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(__e)); abort(); } } while(0)
#  endif  // BHAH_MALLOC_DEVICE__PtrMember
#  ifndef BHAH_FREE_DEVICE__PtrMember
#    define BHAH_FREE_DEVICE__PtrMember(structptr, member) \
       do { if ((structptr)->member) { cudaError_t __e2 = cudaFree((structptr)->member); (structptr)->member = NULL; \
            if (__e2 != cudaSuccess) { fprintf(stderr, "cudaFree failed: %s\n", cudaGetErrorString(__e2)); abort(); } } } while(0)
#  endif  // BHAH_FREE_DEVICE__PtrMember
#else  // !__CUDACC__
#  ifndef BHAH_MALLOC_DEVICE__PtrMember
#    define BHAH_MALLOC_DEVICE__PtrMember(structptr, member, nbytes) do { (void)(structptr); (void)(member); (void)(nbytes); } while(0)
#  endif  // BHAH_MALLOC_DEVICE__PtrMember
#  ifndef BHAH_FREE_DEVICE__PtrMember
#    define BHAH_FREE_DEVICE__PtrMember(structptr, member) do { (void)(structptr); (void)(member); } while(0)
#  endif  // BHAH_FREE_DEVICE__PtrMember
#endif  // __CUDACC__
"""


# --------------------------------------------------------------------------
# Launch-config snippets (classic BHaH style; readable & copy/paste-friendly)
# --------------------------------------------------------------------------
def _emit_launch_setup_1d(n_sym_c: str) -> str:
    return f"""
    const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_DEFAULT;
    dim3 threads_per_block(threads_in_x_dir, 1, 1); // 1-D block to avoid duplicate writes
    dim3 blocks_per_grid((({n_sym_c}) + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
"""


def _emit_launch_setup_2d(n0_c: str, n1_c: str) -> str:
    return f"""
    const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_DEFAULT;
    const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_DEFAULT;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, 1);
    dim3 blocks_per_grid((({n0_c}) + threads_in_x_dir - 1) / threads_in_x_dir,
                         (({n1_c}) + threads_in_y_dir - 1) / threads_in_y_dir,
                         1);
"""


# --------------------------------------------------------------------------
# Main code generator for kernels + defines() body
# --------------------------------------------------------------------------
def generate_rfmprecompute_defines(
    rfm_precompute: ReferenceMetricPrecompute,
) -> Tuple[str, str]:
    """
    Generate the C code for CUDA kernels and the body of the `defines` function.

    Iterates through the precomputed quantities stored in `rfm_precompute`,
    emitting (1) CUDA kernels (prefunc string) and (2) the main function body
    with host loops + runtime CUDA branch.

    :param rfm_precompute: An instance of the ReferenceMetricPrecompute class containing the processed expressions.
    :raises RuntimeError: If an unknown dependency kind is encountered.
    :return: A tuple containing two strings: the C code for the CUDA kernels (`prefuncs_kernels`)
             and the C code for the main `defines` function body (`body`).
    """
    prefuncs_kernels = ""
    body = ""

    for symbol, info in rfm_precompute.rfm_struct__define_kernel_dict.items():
        sname = str(symbol)
        unique_params = info.get("unique_params", [])

        if info["kind"] == "1d":
            ax = info["dep_axis_index"]
            expr_cc = info["expr_cc"]
            # ---- CUDA kernel (1D): pass sizes & scalars explicitly; no device deref of params ----
            kernel_param_sig = (
                (", " + ", ".join(f"const REAL {p}" for p in unique_params))
                if unique_params
                else ""
            )

            # ---- CUDA kernel (1D) ----
            prefuncs_kernels += f"""
#ifdef __CUDACC__
__global__ static void rfm_precompute_defines__{sname}(
    const size_t N,
    rfm_struct * restrict d_rfm,
    const REAL * restrict dx{ax}{kernel_param_sig}) {{
  size_t thread_linear_index = (size_t)threadIdx.x + (size_t)blockIdx.x * (size_t)blockDim.x;
  size_t thread_linear_stride = (size_t)blockDim.x * (size_t)gridDim.x;
  for (size_t i{ax} = thread_linear_index; i{ax} < N; i{ax} += thread_linear_stride) {{
    const REAL xx{ax} = dx{ax}[i{ax}];
    d_rfm->{sname}[i{ax}] = {expr_cc};
  }}
}}
#endif // __CUDACC__
"""

            # ---- Host loop + CUDA branch (1D) ----
            host_params = _emit_host_param_copies(unique_params, indent_spaces=4)
            launch_setup = _emit_launch_setup_1d("N")
            kernel_param_vals = (
                (", " + ", ".join(f"params->{p}" for p in unique_params))
                if unique_params
                else ""
            )
            body += f"""
/* {sname}: 1D precompute */
if (params->is_host) {{
  {{
{host_params}    const size_t N = (size_t)params->Nxx_plus_2NGHOSTS{ax};
    for (size_t i{ax}=0; i{ax}<N; i{ax}++) {{
      const REAL xx{ax} = x{ax}[i{ax}];
      rfmstruct->{sname}[i{ax}] = {expr_cc};
    }}
  }}
}} else {{
  IFCUDARUN({{
    const size_t N = (size_t)params->Nxx_plus_2NGHOSTS{ax};
{launch_setup}    size_t sm = 0;
    const size_t streamid = params->grid_idx % NUM_STREAMS;
    rfm_precompute_defines__{sname}<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(N, rfmstruct, x{ax}{kernel_param_vals});
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__{sname} failure");
  }});
}}

"""

        elif info["kind"] == "2d":
            expr_cc = info["expr_cc"]
            # ---- CUDA kernel (2D): pass sizes & scalars explicitly; no device deref of params ----
            kernel_param_sig = (
                (", " + ", ".join(f"const REAL {p}" for p in unique_params))
                if unique_params
                else ""
            )

            # ---- CUDA kernel (2D) ----
            prefuncs_kernels += f"""
#ifdef __CUDACC__
__global__ static void rfm_precompute_defines__{sname}(
    const size_t N0,
    const size_t N1,
    rfm_struct * restrict d_rfm,
    const REAL * restrict dx0,
    const REAL * restrict dx1{kernel_param_sig}) {{
  size_t idx0 = (size_t)threadIdx.x + (size_t)blockIdx.x * (size_t)blockDim.x;
  size_t idx1 = (size_t)threadIdx.y + (size_t)blockIdx.y * (size_t)blockDim.y;
  size_t stride0 = (size_t)blockDim.x * (size_t)gridDim.x;
  size_t stride1 = (size_t)blockDim.y * (size_t)gridDim.y;
  for (size_t i1=idx1; i1<N1; i1+=stride1) {{
    for (size_t i0=idx0; i0<N0; i0+=stride0) {{
      const REAL xx0 = dx0[i0];
      const REAL xx1 = dx1[i1];
      const size_t flat_idx = i0 + (size_t)N0 * i1;
      d_rfm->{sname}[flat_idx] = {expr_cc};
    }}
  }}
}}
#endif // __CUDACC__

"""

            # ---- Host loop + CUDA branch (2D) ----
            host_params = _emit_host_param_copies(unique_params, indent_spaces=4)
            launch_setup = _emit_launch_setup_2d("N0", "N1")
            kernel_param_vals = (
                (", " + ", ".join(f"params->{p}" for p in unique_params))
                if unique_params
                else ""
            )
            body += f"""
/* {sname}: 2D (xx0,xx1) precompute */
if (params->is_host) {{
  {{
{host_params}    const size_t N0 = (size_t)params->Nxx_plus_2NGHOSTS0;
    const size_t N1 = (size_t)params->Nxx_plus_2NGHOSTS1;
    for (size_t i1=0; i1<N1; i1++) {{
      for (size_t i0=0; i0<N0; i0++) {{
        const REAL xx0 = x0[i0];
        const REAL xx1 = x1[i1];
        const size_t flat_idx = i0 + (size_t)N0 * i1;
        rfmstruct->{sname}[flat_idx] = {expr_cc};
      }}
    }}
  }}
}} else {{
  IFCUDARUN({{
    const size_t N0 = (size_t)params->Nxx_plus_2NGHOSTS0;
    const size_t N1 = (size_t)params->Nxx_plus_2NGHOSTS1;
{launch_setup}    size_t sm = 0;
    const size_t streamid = params->grid_idx % NUM_STREAMS;
    rfm_precompute_defines__{sname}<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(N0, N1, rfmstruct, x0, x1{kernel_param_vals});
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__{sname} failure");
  }});
}}
"""

        else:
            raise RuntimeError("Internal: unknown dependency kind")

    return prefuncs_kernels, body


# --------------------------------------------------------------------------
# Public entry: register malloc/defines/free for all CoordSystems
# --------------------------------------------------------------------------
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
        params = (
            "const commondata_struct * restrict commondata, "
            "const params_struct * restrict params, "
            "rfm_struct * restrict rfmstruct"
        )

        kernels_prefunc, _ = generate_rfmprecompute_defines(rfm_precompute)

        # -------- malloc/free bodies --------
        host_malloc_lines = "".join(
            f"  {line}\n" for line in rfm_precompute.rfm_struct__malloc_host_lines
        )
        device_malloc_lines = "".join(
            f"    {line}\n" for line in rfm_precompute.rfm_struct__malloc_device_lines
        )
        malloc_body = f"""
  // rfm_precompute_malloc: allocate rfmstruct arrays on host or device
  if (params->is_host) {{
    {host_malloc_lines}
  }} else {{
    IFCUDARUN({{ {device_malloc_lines} }});
  }} // END IF params->is_host
"""

        host_free_lines = "".join(
            f"  {line}\n" for line in rfm_precompute.rfm_struct__free_host_lines
        )
        device_free_lines = "".join(
            f"    {line}\n" for line in rfm_precompute.rfm_struct__free_device_lines
        )
        free_body = f"""
  // rfm_precompute_free: free rfmstruct arrays from host or device
  if (params->is_host) {{
    {host_free_lines}
  }} else {{
    IFCUDARUN({{ {device_free_lines} }});
  }} // END IF params->is_host
"""

        # -------- prefuncs --------
        # malloc/free need local CUDA alloc helpers (in case translation unit is compiled in host-only translation).
        malloc_prefunc = _prefunc_cuda_alloc_helpers()
        free_prefunc = _prefunc_cuda_alloc_helpers()
        # defines() only needs the CUDA kernels:
        defines_prefunc = kernels_prefunc

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

        desc_dict = {
            "malloc": """
 * @file rfm_precompute_malloc.c
 * @brief Allocates memory for precomputed reference metric quantities.
 *
 * - When params->is_host==true, perform HOST allocation.
 * - When params->is_host==false, perform DEVICE allocation.
 *
 * @param[in]  commondata  Global simulation metadata (unused).
 * @param[in]  params      Grid parameters (is_host, dims).
 * @param[out] rfmstruct   Destination struct for allocated pointers.
 *
 * @return void.
 """,
            "defines": """
 * @file rfm_precompute_defines.c
 * @brief Computes and populates arrays with precomputed reference metric quantities.
 *
 * - params->is_host==true: rfmstruct & xx[] are HOST; run CPU loops.
 * - params->is_host==false: rfmstruct & xx[] are DEVICE/UVM; launch CUDA kernels.
 *
 * @param[in]  commondata  Global simulation metadata (unused).
 * @param[in]  params      Grid parameters and dimensions.
 * @param[out] rfmstruct   Struct containing arrays to populate.
 * @param[in]  xx          Pointers to 1D coordinate arrays.
 *
 * @return void.
 """,
            "free": """
 * @file rfm_precompute_free.c
 * @brief Frees memory for precomputed reference metric quantities.
 *
 * - params->is_host==true: free HOST memory.
 * - params->is_host==false: free DEVICE memory.
 *
 * @param[in]  commondata  Global simulation metadata (unused).
 * @param[in]  params      Grid parameters.
 * @param[out] rfmstruct   Struct whose members will be freed.
 *
 * @return void.
 """,
        }

        for func_name, func_body in c_functions_to_register.items():
            function_desc = desc_dict[func_name]
            function_name = f"rfm_precompute_{func_name}"
            params_sig = params + (
                ", REAL * restrict xx[3]" if func_name == "defines" else ""
            )

            if_is_host_comment = """
// If params->is_host==true: rfmstruct, xx[] are HOST pointers; CPU path executes.
// If params->is_host==false: rfmstruct, xx[] are DEVICE or UVM pointers; GPU path executes.
"""
            final_body = if_is_host_comment + (func_body if func_body else " ")

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

    # rfm_struct typedef
    rfm_struct_typedef = "typedef struct __rfmstruct__ {\n"
    rfm_struct_typedef += "\n".join(sorted(superfast_uniq(combined_BHaH_defines_list)))
    rfm_struct_typedef += "\n} rfm_struct;\n"
    BHaH.BHaH_defines_h.register_BHaH_defines("reference_metric", rfm_struct_typedef)


# --------------------------------------------------------------------------
# doctest harness
# --------------------------------------------------------------------------
if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
