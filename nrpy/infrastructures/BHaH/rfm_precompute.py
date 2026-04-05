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

from typing import List, Set, Tuple

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
# Launch-config snippets (classic BHaH style)
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
# Params helper: signature suffix, call suffix, and host param copies (single source)
# --------------------------------------------------------------------------
def _emit_param_fragments(unique_params: List[str]) -> Tuple[str, str, str]:
    """
    Build kernel parameter fragments and host-side parameter copies.

    :param unique_params: Parameter names to include, in order.
    :return: A tuple of three strings: (kernel_sig, kernel_vals, host_copies), where
        kernel_sig is a kernel signature suffix, kernel_vals is a kernel call suffix,
        and host_copies is a newline-terminated block of C statements declaring local
        host copies of parameters.
    """
    if not unique_params:
        return "", "", ""

    kernel_sig = ", " + ", ".join(f"const REAL {p}" for p in unique_params)
    kernel_vals = ", " + ", ".join(f"params->{p}" for p in unique_params)
    host_copies = (
        "\n".join(f"const REAL {p} = params->{p};" for p in unique_params) + "\n"
    )
    return kernel_sig, kernel_vals, host_copies


# --------------------------------------------------------------------------
# Orchestrator: build per-symbol allocation/compute/free & CUDA kernels
# --------------------------------------------------------------------------
class ReferenceMetricPrecompute:
    """
    Manages the generation of C code for precomputing reference metric quantities.

    - When params->is_host!=0, rfmstruct and xx[] are HOST pointers.
    - When params->is_host==0, rfmstruct and xx[] are DEVICE-resident or UVM pointers.

    Host-side SIMD reader strings are generated, which are rewritten to CUDA variants
    on CUDA-enabled builds.
    """

    def __init__(self, CoordSystem: str) -> None:
        self._initialize_configuration(CoordSystem)
        self._initialize_properties()

        for symbol, expr in self._get_sorted_precomputed_expressions():
            if "_of_xx" in str(symbol):
                self._process_expression(symbol, expr)

        # Rewrite host-only reader intrinsics to CUDA variants when targeting CUDA.
        self._apply_intrinsics_mode()

        # Finalize kernels/body once; reuse everywhere.
        self.kernels_prefunc = "".join(self._kernels_parts)
        defines_body = "".join(self._defines_parts)
        alias_lines = "\n".join(
            [f"MAYBE_UNUSED const REAL * restrict x{d} = xx[{d}];" for d in range(3)]
        )
        self.rfm_struct__define = alias_lines + "\n" + defines_body
        if self.rfm_struct__define and not self.rfm_struct__define.endswith("\n"):
            self.rfm_struct__define += "\n"

    # --------------------------- setup ---------------------------
    def _initialize_configuration(self, CoordSystem: str) -> None:
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
        # Struct member decls for BHaH_defines.h:
        self.BHaH_defines_list: List[str] = []

        # Single source of truth for alloc/free: (member_name, size_expr)
        self.member_specs: List[Tuple[str, str]] = []

        # Main emitted code (single-pass, ordered):
        self._kernels_parts: List[str] = []
        self._defines_parts: List[str] = []

        # Host-only readers (later rewritten to CUDA variants if needed):
        self.readvr_str = [""] * self.N_DIMS
        self.readvr_intrinsics_outer_str = [""] * self.N_DIMS
        self.readvr_intrinsics_inner_str = [""] * self.N_DIMS

        # Final outputs (set once at end of __init__):
        self.kernels_prefunc = ""
        self.rfm_struct__define = ""

    # --------------------------- discovery ---------------------------
    def _get_sorted_precomputed_expressions(self) -> List[Tuple[sp.Expr, sp.Expr]]:
        if not self.rfm.freevars_uniq_xx_indep:
            return []
        zipped_pairs = zip(self.rfm.freevars_uniq_xx_indep, self.rfm.freevars_uniq_vals)
        return sorted(zipped_pairs, key=lambda pair: str(pair[0]))

    # --------------------------- small helper ---------------------------
    def _ccg_rhs(self, expr: sp.Expr) -> str:
        """
        Convert a SymPy expression into a C RHS expression string.

        :param expr: SymPy expression to convert.
        :return: The right-hand-side C expression as a string, without an assignment
            target or a trailing semicolon.
        """
        code = ccg.c_codegen(
            expr,
            "REAL __rhs",
            include_braces=False,
            verbose=False,
            enable_cse=False,
        )
        line = next(
            l
            for l in (s.strip() for s in code.splitlines())
            if ("=" in l and l.endswith(";"))
        )
        return line.split("=", 1)[1].strip().rstrip(";")

    def _size_expr_from_dependencies(self, dependencies: List[bool]) -> str:
        size_terms = [
            f"(size_t)params->Nxx_plus_2NGHOSTS{dim}"
            for dim, dep in enumerate(dependencies)
            if dep
        ]
        return " * ".join(size_terms) if size_terms else "(size_t)1"

    # --------------------------- per-symbol processing (single-pass emit) ---------------------------
    def _process_expression(self, symbol: sp.Expr, expr: sp.Expr) -> None:
        # Dependency detection:
        dependencies = [expr.has(self.rfm.xx[i]) for i in range(self.N_DIMS)]

        # Register struct member + alloc/free spec:
        sname = str(symbol)
        self.BHaH_defines_list.append(f"REAL * {sname};")
        self.member_specs.append(
            (sname, self._size_expr_from_dependencies(dependencies))
        )

        has_1d_dep = sum(dependencies) == 1
        has_2d_dep_xx0_xx1 = dependencies[0] and dependencies[1] and not dependencies[2]

        expr_cc = self._ccg_rhs(expr)
        unique_params = get_unique_expression_symbols_as_strings(
            expr, exclude=[f"xx{j}" for j in range(self.N_DIMS)]
        )
        kernel_param_sig, kernel_param_vals, host_param_copies = _emit_param_fragments(
            unique_params
        )

        if has_1d_dep:
            ax = dependencies.index(True)
            self._emit_1d(
                symbol_name=sname,
                ax=ax,
                expr_cc=expr_cc,
                kernel_param_sig=kernel_param_sig,
                kernel_param_vals=kernel_param_vals,
                host_param_copies=host_param_copies,
            )
            self._emit_readers_1d(symbol_name=sname, ax=ax)

        elif has_2d_dep_xx0_xx1:
            self._emit_2d(
                symbol_name=sname,
                expr_cc=expr_cc,
                kernel_param_sig=kernel_param_sig,
                kernel_param_vals=kernel_param_vals,
                host_param_copies=host_param_copies,
            )
            self._emit_readers_2d(symbol_name=sname)

        else:
            raise RuntimeError(f"Unsupported dependency for {symbol}: {expr}")

    # --------------------------- emitters ---------------------------
    def _emit_1d(
        self,
        symbol_name: str,
        ax: int,
        expr_cc: str,
        kernel_param_sig: str,
        kernel_param_vals: str,
        host_param_copies: str,
    ) -> None:
        # CUDA kernel (1D)
        self._kernels_parts.append(f"""
#ifdef __CUDACC__
__global__ static void rfm_precompute_defines__{symbol_name}(
    const size_t N,
    rfm_struct * restrict d_rfm,
    const REAL * restrict dx{ax}{kernel_param_sig}) {{
  size_t thread_linear_index = (size_t)threadIdx.x + (size_t)blockIdx.x * (size_t)blockDim.x;
  size_t thread_linear_stride = (size_t)blockDim.x * (size_t)gridDim.x;
  for (size_t i{ax} = thread_linear_index; i{ax} < N; i{ax} += thread_linear_stride) {{
    const REAL xx{ax} = dx{ax}[i{ax}];
    d_rfm->{symbol_name}[i{ax}] = {expr_cc};
  }}
}}
#endif // __CUDACC__
""")

        launch_setup = _emit_launch_setup_1d("N")

        # Host loop + CUDA branch (1D)
        self._defines_parts.append(f"""
/* {symbol_name}: 1D precompute */
if (params->is_host) {{
  {{
{host_param_copies}    const size_t N = (size_t)params->Nxx_plus_2NGHOSTS{ax};
    for (size_t i{ax}=0; i{ax}<N; i{ax}++) {{
      const REAL xx{ax} = x{ax}[i{ax}];
      rfmstruct->{symbol_name}[i{ax}] = {expr_cc};
    }}
  }}
}} else {{
  IFCUDARUN({{
    const size_t N = (size_t)params->Nxx_plus_2NGHOSTS{ax};
{launch_setup}    size_t sm = 0;
    const size_t streamid = params->grid_idx % NUM_STREAMS;
    rfm_precompute_defines__{symbol_name}<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(N, rfmstruct, x{ax}{kernel_param_vals});
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__{symbol_name} failure");
  }});
}}

""")

    def _emit_2d(
        self,
        symbol_name: str,
        expr_cc: str,
        kernel_param_sig: str,
        kernel_param_vals: str,
        host_param_copies: str,
    ) -> None:
        # CUDA kernel (2D)
        self._kernels_parts.append(f"""
#ifdef __CUDACC__
__global__ static void rfm_precompute_defines__{symbol_name}(
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
      d_rfm->{symbol_name}[flat_idx] = {expr_cc};
    }}
  }}
}}
#endif // __CUDACC__

""")

        launch_setup = _emit_launch_setup_2d("N0", "N1")

        # Host loop + CUDA branch (2D)
        self._defines_parts.append(f"""
/* {symbol_name}: 2D (xx0,xx1) precompute */
if (params->is_host) {{
  {{
{host_param_copies}    const size_t N0 = (size_t)params->Nxx_plus_2NGHOSTS0;
    const size_t N1 = (size_t)params->Nxx_plus_2NGHOSTS1;
    for (size_t i1=0; i1<N1; i1++) {{
      for (size_t i0=0; i0<N0; i0++) {{
        const REAL xx0 = x0[i0];
        const REAL xx1 = x1[i1];
        const size_t flat_idx = i0 + (size_t)N0 * i1;
        rfmstruct->{symbol_name}[flat_idx] = {expr_cc};
      }}
    }}
  }}
}} else {{
  IFCUDARUN({{
    const size_t N0 = (size_t)params->Nxx_plus_2NGHOSTS0;
    const size_t N1 = (size_t)params->Nxx_plus_2NGHOSTS1;
{launch_setup}    size_t sm = 0;
    const size_t streamid = params->grid_idx % NUM_STREAMS;
    rfm_precompute_defines__{symbol_name}<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(N0, N1, rfmstruct, x0, x1{kernel_param_vals});
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__{symbol_name} failure");
  }});
}}
""")

    # --------------------------- readers ---------------------------
    def _emit_readers_1d(self, symbol_name: str, ax: int) -> None:
        self.readvr_str[
            ax
        ] += f"MAYBE_UNUSED const REAL {symbol_name} = rfmstruct->{symbol_name}[i{ax}];\n"
        self.readvr_intrinsics_outer_str[ax] += (
            f"const double NOSIMD{symbol_name} = rfmstruct->{symbol_name}[i{ax}]; "
            f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol_name} = ConstSIMD(NOSIMD{symbol_name});\n"
        )
        self.readvr_intrinsics_inner_str[
            ax
        ] += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol_name} = ReadSIMD(&rfmstruct->{symbol_name}[i{ax}]);\n"

    def _emit_readers_2d(self, symbol_name: str) -> None:
        flat = "(i0 + (size_t)N0 * (size_t)i1)"
        self.readvr_str[
            0
        ] += f"MAYBE_UNUSED const REAL {symbol_name} = rfmstruct->{symbol_name}[{flat}];\n"
        self.readvr_intrinsics_outer_str[0] += (
            f"const double NOSIMD{symbol_name} = rfmstruct->{symbol_name}[{flat}]; "
            f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol_name} = ConstSIMD(NOSIMD{symbol_name});\n"
        )
        self.readvr_intrinsics_inner_str[
            0
        ] += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol_name} = ReadSIMD(&rfmstruct->{symbol_name}[{flat}]);\n"

    # --------------------------- SIMD<->CUDA reader replacements ---------------------------
    def _apply_intrinsics_mode(self) -> None:
        if not self.is_cuda:
            return
        self.readvr_intrinsics_outer_str = [
            s.replace("SIMD", "CUDA") for s in self.readvr_intrinsics_outer_str
        ]
        self.readvr_intrinsics_inner_str = [
            s.replace("SIMD", "CUDA") for s in self.readvr_intrinsics_inner_str
        ]


# --------------------------------------------------------------------------
# Public entry: register malloc/defines/free for all CoordSystems
# --------------------------------------------------------------------------
def register_CFunctions_rfm_precompute(set_of_CoordSystems: Set[str]) -> None:
    """
    Register C functions for reference-metric precomputation.

    For each coordinate system, this registers `rfm_precompute_malloc`,
    `rfm_precompute_defines`, and `rfm_precompute_free`.

    :param set_of_CoordSystems: Coordinate systems to generate and register functions
        for.
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

        # --- malloc/free bodies from single member-spec list (host/device macro swap) ---
        host_malloc_lines = "\n".join(
            f"BHAH_MALLOC__PtrMember(rfmstruct, {m}, sizeof(REAL) * ({sz}));"
            for m, sz in rfm_precompute.member_specs
        )
        device_malloc_lines = "\n".join(
            f"BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, {m}, sizeof(REAL) * ({sz}));"
            for m, sz in rfm_precompute.member_specs
        )
        malloc_body = f"""
// rfm_precompute_malloc: allocate rfmstruct arrays on host or device
if (params->is_host) {{
{host_malloc_lines}
}} else {{
IFCUDARUN({{
{device_malloc_lines}
}});
}} // END IF params->is_host
"""
        host_free_lines = "\n".join(
            f"BHAH_FREE__PtrMember(rfmstruct, {m});"
            for m, _ in rfm_precompute.member_specs
        )
        device_free_lines = "\n".join(
            f"BHAH_FREE_DEVICE__PtrMember(rfmstruct, {m});"
            for m, _ in rfm_precompute.member_specs
        )
        free_body = f"""
// rfm_precompute_free: free rfmstruct arrays from host or device
if (params->is_host) {{
{host_free_lines}
}} else {{
IFCUDARUN({{
{device_free_lines}
}});
}} // END IF params->is_host
"""

        # --- prefuncs ---
        malloc_prefunc = ""
        free_prefunc = ""
        defines_prefunc = rfm_precompute.kernels_prefunc

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

        if_is_host_comment = """
// If params->is_host==true: rfmstruct, xx[] are HOST pointers; CPU path executes.
// If params->is_host==false: rfmstruct, xx[] are DEVICE or UVM pointers; GPU path executes.
"""

        for func_name, func_body in c_functions_to_register.items():
            function_desc = desc_dict[func_name]
            function_name = f"rfm_precompute_{func_name}"
            params_sig = params + (
                ", REAL * restrict xx[3]" if func_name == "defines" else ""
            )

            final_body = if_is_host_comment + (func_body if func_body else " ")
            if final_body and not final_body.endswith("\n"):
                final_body += "\n"

            prefunc = prefunc_dict.get(func_name, "")
            if prefunc and not prefunc.endswith("\n"):
                prefunc += "\n"

            cfc.register_CFunction(
                prefunc=prefunc,
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
