"""
Register a C function that exports Cartesian raytracing data during evolution.

This module generates a BHaH diagnostics helper that, on each scheduled output
step, recomputes same-slice Ricci/RHS data, evaluates the Cartesian covariant
four-metric and four-Christoffel symbols from the native BSSN evolution state
using the transformed symbolic recipes in
``nrpy.equations.general_relativity.geodesics.geodesics``, and writes a
stable binary payload for later raytracing.

The payload stores the final Cartesian tensor data directly. No later
reconstruction or runtime tensor basis transformation is
required to recover the metric and Christoffels written here.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Tuple, Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)
from nrpy.infrastructures import BHaH


def _validate_fixed_width_string(text: str, field_bytes: int, label: str) -> None:
    """
    Validate that a fixed-width on-disk string fits its serialized field.

    The binary writer stores null-padded strings and rejects any text whose
    encoded length would fill or exceed the field. Checking this at Python
    code-generation time is cheaper than discovering the mismatch during a run.

    :param text: String that will be serialized into the fixed-width field.
    :param field_bytes: Exact field width in bytes.
    :param label: Human-readable label for error reporting.
    :raises ValueError: If the string cannot fit in the fixed-width field.
    """
    if len(text.encode("ascii")) >= field_bytes:
        raise ValueError(
            f"String '{text}' for {label} exceeds fixed-width field size "
            f"{field_bytes}."
        )


def register_CFunction_output_raytracing_data(
    CoordSystem: str,
    enable_rfm_precompute: bool,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct and register the stage-1 Cartesian raytracing binary exporter.

    The generated C helper writes one binary file per diagnostics output with
    filename pattern ``./raytracing_data_t########.bin``. Each file contains
    a fixed-size self-describing header, the physical simulation time, and one
    point record per exported logical-grid point. Each point record stores:

    1. Cartesian coordinates ``(x, y, z)``.
    2. The 10 unique covariant four-metric components ``g4DD`` in a fixed
       upper-triangular order.
    3. The 40 unique four-Christoffel components ``Gamma4UDD`` with the lower
       pair serialized in upper-triangular order.

    In v1, the exporter evaluates only interior points. This keeps the
    finite-difference derivative path consistent with the current BHaH
    infrastructure and avoids writing connection data at points where the
    centered stencil would step beyond the available storage. The header
    records that ghost-zone points are excluded from the payload.

    :param CoordSystem: Coordinate system used by the evolved BSSN state.
    :param enable_rfm_precompute: Whether the generated project uses
        reference-metric precomputation.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If ``CoordSystem`` is not a string.
    :raises ValueError: If ``enable_rfm_precompute`` is not a bool.

    Doctests:
    None.
    """
    if not isinstance(CoordSystem, str):
        raise ValueError(
            f"CoordSystem must be a string, got {type(CoordSystem).__name__}"
        )
    if not isinstance(enable_rfm_precompute, bool):
        raise ValueError(
            "enable_rfm_precompute must be a bool, "
            f"got {type(enable_rfm_precompute).__name__}"
        )

    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Import the BSSN quantities module so the conformal-factor parameter is registered.
    # pylint: disable-next=import-outside-toplevel
    __import__("nrpy.equations.general_relativity.BSSN_quantities")

    # pylint: disable-next=import-outside-toplevel
    from nrpy.equations.general_relativity.geodesics.geodesics import (
        symbolic_christoffel_recipe_from_bssn_grid_basis,
        symbolic_g4DD_recipe_from_bssn_grid_basis,
    )

    # Step 1: Build the transformed symbolic recipes for the exported tensors.
    g4DD = symbolic_g4DD_recipe_from_bssn_grid_basis(
        bssn_coord_system=CoordSystem,
        target_basis="Cartesian",
        enable_bssn_rfm_precompute=enable_rfm_precompute,
    )
    Gamma4UDD = symbolic_christoffel_recipe_from_bssn_grid_basis(
        bssn_coord_system=CoordSystem,
        target_basis="Cartesian",
        enable_bssn_rfm_precompute=enable_rfm_precompute,
    )

    # Step 2: Define the on-disk tensor component ordering.
    metric_components: List[Tuple[Tuple[int, int], str]] = [
        ((0, 0), "g4DD00"),
        ((0, 1), "g4DD01"),
        ((0, 2), "g4DD02"),
        ((0, 3), "g4DD03"),
        ((1, 1), "g4DD11"),
        ((1, 2), "g4DD12"),
        ((1, 3), "g4DD13"),
        ((2, 2), "g4DD22"),
        ((2, 3), "g4DD23"),
        ((3, 3), "g4DD33"),
    ]
    christoffel_components: List[Tuple[Tuple[int, int, int], str]] = []
    for alpha in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                christoffel_components.append(
                    ((alpha, mu, nu), f"Gamma4UDD{alpha}{mu}{nu}")
                )

    # Step 3: Flatten the symbolic expressions into the serialized record layout.
    expr_list = [g4DD[mu][nu] for (mu, nu), _ in metric_components]
    expr_list.extend(
        Gamma4UDD[alpha][mu][nu] for (alpha, mu, nu), _ in christoffel_components
    )

    lhs_list = [f"REAL {name}" for _, name in metric_components]
    lhs_list.extend(f"REAL {name}" for _, name in christoffel_components)

    # Step 4: Ask the shared expression helpers which runtime parameters the
    #         symbolic recipes need, then generate the corresponding C locals.
    params_symbols, commondata_symbols = get_params_commondata_symbols_from_expr_list(
        expr_list
    )
    params_definitions = generate_definition_header(
        params_symbols,
        enable_intrinsics=False,
        var_access="params->",
    )
    commondata_definitions = generate_definition_header(
        commondata_symbols,
        enable_intrinsics=False,
        var_access="commondata->",
    )

    # Step 5: Define the fixed-width header schema.
    record_component_names = ["x", "y", "z"]
    record_component_names.extend(name for _, name in metric_components)
    record_component_names.extend(name for _, name in christoffel_components)
    header_label = "NRPy RT spacetime data"
    record_component_name_bytes = 24
    header_label_bytes = 32
    format_name = "Cartesian g4DD+Gamma4UDD"
    format_name_bytes = 32
    target_basis_name = "Cartesian"
    basis_name_bytes = 16
    coord_name_bytes = 32
    loop_order_name = "i2maj_i0fast"
    loop_order_name_bytes = 16
    metric_component_name_bytes = 16
    gamma_component_name_bytes = 16
    num_metric_components = len(metric_components)
    num_gamma_components = len(christoffel_components)
    num_record_components = len(record_component_names)
    format_version = 1
    point_record_real_count = 3 + num_metric_components + num_gamma_components
    point_record_bytes = 8 * point_record_real_count
    num_u32_header_fields = 20
    num_u64_header_fields = 15
    num_f64_header_fields = 15
    cf_convention_bytes = format_name_bytes
    header_size_bytes = (
        16
        + num_u32_header_fields * 4
        + num_u64_header_fields * 8
        + num_f64_header_fields * 8
        + header_label_bytes
        + format_name_bytes
        + basis_name_bytes
        + coord_name_bytes
        + loop_order_name_bytes
        + cf_convention_bytes
        + num_record_components * record_component_name_bytes
        + num_metric_components * metric_component_name_bytes
        + num_gamma_components * gamma_component_name_bytes
    )
    cf_convention = cast(str, par.parval_from_str("EvolvedConformalFactor_cf"))

    # Step 1: Validate every fixed-width string before generating the C exporter.
    _validate_fixed_width_string(header_label, header_label_bytes, "header_label")
    _validate_fixed_width_string(format_name, format_name_bytes, "format_name")
    _validate_fixed_width_string(target_basis_name, basis_name_bytes, "target_basis")
    _validate_fixed_width_string(CoordSystem, coord_name_bytes, "source_coord_system")
    _validate_fixed_width_string(loop_order_name, loop_order_name_bytes, "loop_order")
    _validate_fixed_width_string(cf_convention, cf_convention_bytes, "cf_convention")
    for component_name in record_component_names:
        _validate_fixed_width_string(
            component_name, record_component_name_bytes, "record_component_names"
        )
    for _, component_name in metric_components:
        _validate_fixed_width_string(
            component_name, metric_component_name_bytes, "metric_component_names"
        )
    for _, component_name in christoffel_components:
        _validate_fixed_width_string(
            component_name, gamma_component_name_bytes, "christoffel_component_names"
        )

    # Step 6: Build the interior-point evaluation kernel for one output record.
    loop_body = f"""
{commondata_definitions}
{params_definitions}
  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;
  MAYBE_UNUSED const REAL xx0 = xx[0][i0];
  MAYBE_UNUSED const REAL xx1 = xx[1][i1];
  MAYBE_UNUSED const REAL xx2 = xx[2][i2];
  const int idx3 = IDX3(i0, i1, i2);
  const REAL alpha_d0 = rhs_gfs[IDX4pt(ALPHAGF, idx3)];
  const REAL vetU_d00 = rhs_gfs[IDX4pt(VETU0GF, idx3)];
  const REAL vetU_d01 = rhs_gfs[IDX4pt(VETU1GF, idx3)];
  MAYBE_UNUSED const REAL vetU_d02 = rhs_gfs[IDX4pt(VETU2GF, idx3)];
  const REAL *restrict in_gfs = y_n_gfs;
  const REAL xOrig[3] = {{xx0, xx1, xx2}};
  REAL xCart[3];
  xx_to_Cart(params, xOrig, xCart);

{ccg.c_codegen(
        expr_list,
        lhs_list,
        enable_fd_codegen=True,
        enable_simd=False,
        enable_fd_functions=False,
        include_braces=False,
    )}
  raytracing_data_write_f64_or_abort(fp, (double)xCart[0], "x");
  raytracing_data_write_f64_or_abort(fp, (double)xCart[1], "y");
  raytracing_data_write_f64_or_abort(fp, (double)xCart[2], "z");
"""
    for _, name in metric_components:
        loop_body += (
            f'  raytracing_data_write_f64_or_abort(fp, (double){name}, "{name}");\n'
        )
    for _, name in christoffel_components:
        loop_body += (
            f'  raytracing_data_write_f64_or_abort(fp, (double){name}, "{name}");\n'
        )

    # Step 7: Wrap the pointwise kernel in the standard BHaH interior loop.
    loop = BHaH.simple_loop.simple_loop(
        loop_body=loop_body,
        loop_region="interior",
        enable_intrinsics=False,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=1,
    )

    # Step 8: Register the generated C exporter and its helper functions.
    includes = [
        "<errno.h>",
        "<float.h>",
        "<math.h>",
        "<stdint.h>",
        "<stdio.h>",
        "<stdlib.h>",
        "<string.h>",
        "<sys/stat.h>",
        "<unistd.h>",
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = """
Export final Cartesian metric and Christoffel data for raytracing.

This function writes at most one binary file per diagnostics output to:

  ./raytracing_data_t########.bin

where the eight-digit index is the diagnostics output index. Each file is
written to a unique temporary sibling path and then installed at the final
path without overwriting an existing file. The payload stores the physical
simulation time and one point record per exported interior logical-grid point.

Each point record contains:
  - Cartesian coordinates x, y, z,
  - the 10 unique covariant four-metric components
      (00, 01, 02, 03, 11, 12, 13, 22, 23, 33),
  - the 40 unique four-Christoffel components serialized by nested
      (alpha, mu, nu) order with nu >= mu.

The Cartesian tensor expressions are evaluated directly from the transformed
symbolic recipes in geodesics.py. The exporter does not perform a second
runtime tensor basis transformation after importing those recipes.

Before evaluating the Christoffels, the exporter refreshes same-slice
auxiliary data using:
  Ricci_eval(params, rfmstruct, y_n_gfs, auxevol_gfs);
  rhs_eval(commondata, params, rfmstruct, auxevol_gfs, y_n_gfs, rhs_gfs);

where rhs_gfs is a temporary EVOL-sized scratch buffer allocated inside this
function. The header explicitly records that only interior logical-grid points
are exported in v1, so ghost-zone points are not serialized.

@param[in] commondata  Global simulation metadata used in output naming and headers
@param[in] griddata  Host-side per-grid runtime state to export
@param time_output_index  Diagnostics output index used in the emitted filename
"""
    cfunc_type = "void"
    name = "output_raytracing_data"
    params = (
        "const commondata_struct *restrict commondata, "
        "const griddata_struct *restrict griddata, "
        "const int time_output_index"
    )

    prefunc = r"""
/**
 * Abort the program with a message.
 *
 * @param[in] message  Error message to print before aborting.
 */
static void raytracing_data_abort_with_message(const char *restrict message) {
  fprintf(stderr, "%s\n", message);
  exit(1);
} // END FUNCTION: raytracing_data_abort_with_message

/**
 * Cast a nonnegative integer to uint32_t, aborting on invalid input.
 *
 * @param value  Integer value to convert.
 * @param[in] label  Label for the error message.
 * @return Safely casted uint32_t value.
 */
static uint32_t raytracing_data_u32_from_nonnegative_int_or_abort(
    const int value,
    const char *restrict label
) {
  if (value < 0) {
    fprintf(stderr,
            "Error: output_raytracing_data received negative %s=%d.\n",
            label,
            value);
    exit(1);
  } // END IF: negative integers cannot be serialized as uint32_t
  if ((uint64_t)value > (uint64_t)UINT32_MAX) {
    fprintf(stderr,
            "Error: output_raytracing_data %s=%d does not fit in uint32_t.\n",
            label,
            value);
    exit(1);
  } // END IF: value exceeds the documented header format
  return (uint32_t)value;
} // END FUNCTION: raytracing_data_u32_from_nonnegative_int_or_abort

/**
 * Add two uint64_t values, aborting on overflow.
 *
 * @param a  First value.
 * @param b  Second value.
 * @param[in] label  Label for the error message.
 * @return The sum a + b.
 */
static uint64_t raytracing_data_add_u64_or_abort(
    const uint64_t a,
    const uint64_t b,
    const char *restrict label
) {
  if (UINT64_MAX - a < b) {
    fprintf(stderr,
            "Error: output_raytracing_data uint64 overflow while computing %s.\n",
            label);
    exit(1);
  } // END IF: addition would overflow
  return a + b;
} // END FUNCTION: raytracing_data_add_u64_or_abort

/**
 * Multiply two uint64_t values, aborting on overflow.
 *
 * @param a  First value.
 * @param b  Second value.
 * @param[in] label  Label for the error message.
 * @return The product a * b.
 */
static uint64_t raytracing_data_mul_u64_or_abort(
    const uint64_t a,
    const uint64_t b,
    const char *restrict label
) {
  if (a != 0ULL && b > UINT64_MAX / a) {
    fprintf(stderr,
            "Error: output_raytracing_data uint64 overflow while computing %s.\n",
            label);
    exit(1);
  } // END IF: multiplication would overflow
  return a * b;
} // END FUNCTION: raytracing_data_mul_u64_or_abort

/**
 * Write elements to a file, aborting on failure.
 *
 * @param[in,out] fp  File pointer.
 * @param[in] data  Data to write.
 * @param element_size  Size of each element.
 * @param count  Number of elements.
 * @param[in] label  Label for the error message.
 */
static void raytracing_data_write_or_abort(
    FILE *restrict fp,
    const void *restrict data,
    const size_t element_size,
    const size_t count,
    const char *restrict label
) {
  if (count != 0 && element_size > SIZE_MAX / count) {
    fprintf(stderr,
            "Error: output_raytracing_data size_t overflow before writing %s.\n",
            label);
    fclose(fp);
    exit(1);
  } // END IF: fwrite byte count would overflow
  if (fwrite(data, element_size, count, fp) != count) {
    fprintf(stderr,
            "Error: output_raytracing_data failed while writing %s.\n",
            label);
    fclose(fp);
    exit(1);
  } // END IF: fwrite failed
} // END FUNCTION: raytracing_data_write_or_abort

/**
 * Write a uint32_t value in little-endian byte order.
 *
 * @param[in,out] fp  File pointer.
 * @param value  Value to write.
 * @param[in] label  Label for the error message.
 */
static void raytracing_data_write_u32_or_abort(
    FILE *restrict fp,
    const uint32_t value,
    const char *restrict label
) {
  uint8_t buffer[4];
  buffer[0] = (uint8_t)(value & 0xffU);
  buffer[1] = (uint8_t)((value >> 8U) & 0xffU);
  buffer[2] = (uint8_t)((value >> 16U) & 0xffU);
  buffer[3] = (uint8_t)((value >> 24U) & 0xffU);
  raytracing_data_write_or_abort(fp, buffer, sizeof(uint8_t), 4, label);
} // END FUNCTION: raytracing_data_write_u32_or_abort

/**
 * Write a uint64_t value in little-endian byte order.
 *
 * @param[in,out] fp  File pointer.
 * @param value  Value to write.
 * @param[in] label  Label for the error message.
 */
static void raytracing_data_write_u64_or_abort(
    FILE *restrict fp,
    const uint64_t value,
    const char *restrict label
) {
  uint8_t buffer[8];
  buffer[0] = (uint8_t)(value & 0xffULL);
  buffer[1] = (uint8_t)((value >> 8U) & 0xffULL);
  buffer[2] = (uint8_t)((value >> 16U) & 0xffULL);
  buffer[3] = (uint8_t)((value >> 24U) & 0xffULL);
  buffer[4] = (uint8_t)((value >> 32U) & 0xffULL);
  buffer[5] = (uint8_t)((value >> 40U) & 0xffULL);
  buffer[6] = (uint8_t)((value >> 48U) & 0xffULL);
  buffer[7] = (uint8_t)((value >> 56U) & 0xffULL);
  raytracing_data_write_or_abort(fp, buffer, sizeof(uint8_t), 8, label);
} // END FUNCTION: raytracing_data_write_u64_or_abort

/**
 * Write a binary64 value in little-endian byte order.
 *
 * @param[in,out] fp  File pointer.
 * @param value  Value to write.
 * @param[in] label  Label for the error message.
 */
static void raytracing_data_write_f64_or_abort(
    FILE *restrict fp,
    const double value,
    const char *restrict label
) {
  uint64_t bits = 0ULL;
  memcpy(&bits, &value, sizeof(bits));
  raytracing_data_write_u64_or_abort(fp, bits, label);
} // END FUNCTION: raytracing_data_write_f64_or_abort

/**
 * Validate the binary64 output assumptions before serialization.
 */
static void raytracing_data_validate_binary64_output_or_abort(void) {
#if FLT_RADIX != 2 || DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024 || DBL_MIN_EXP != -1021
  raytracing_data_abort_with_message(
      "Error: raytracing output requires C double to match IEEE-754 binary64 semantics.");
#endif
  if (sizeof(double) != 8)
    raytracing_data_abort_with_message(
        "Error: raytracing output requires 8-byte C double for binary64 serialization.");
  if (sizeof(REAL) > sizeof(double))
    raytracing_data_abort_with_message(
        "Error: raytracing output refuses to silently down-convert REAL values wider than binary64.");
} // END FUNCTION: raytracing_data_validate_binary64_output_or_abort

/**
 * Write a fixed-width string, null-padded to the requested size.
 *
 * @param[in,out] fp  File pointer.
 * @param[in] text  Text to write.
 * @param field_bytes  Exact field width in bytes.
 * @param[in] label  Label for the error message.
 */
static void raytracing_data_write_fixed_length_string_or_abort(
    FILE *restrict fp,
    const char *restrict text,
    const size_t field_bytes,
    const char *restrict label
) {
  char buffer[128];
  const size_t text_len = strlen(text);
  if (field_bytes > sizeof(buffer))
    raytracing_data_abort_with_message(
        "Error: raytracing header string field exceeds the internal serialization buffer.");
  if (field_bytes == 0 || text_len >= field_bytes) {
    fprintf(stderr,
            "Error: output_raytracing_data string field %s is too long for its fixed header field.\n",
            label);
    fclose(fp);
    exit(1);
  } // END IF: fixed-width string would be truncated
  memset(buffer, 0, sizeof(buffer));
  memcpy(buffer, text, text_len);
  raytracing_data_write_or_abort(fp, buffer, sizeof(char), field_bytes, label);
} // END FUNCTION: raytracing_data_write_fixed_length_string_or_abort

"""

    body = rf"""
  if (commondata->NUMGRIDS != 1)
    raytracing_data_abort_with_message(
        "Error: output_raytracing_data currently requires commondata->NUMGRIDS == 1.");
  if (time_output_index < 0)
    raytracing_data_abort_with_message(
        "Error: output_raytracing_data received a negative output index.");

  raytracing_data_validate_binary64_output_or_abort();

  const params_struct *restrict params = &griddata[0].params;
  const rfm_struct *restrict rfmstruct = griddata[0].rfmstruct;
  const REAL *restrict y_n_gfs = griddata[0].gridfuncs.y_n_gfs;
  REAL *restrict auxevol_gfs = griddata[0].gridfuncs.auxevol_gfs;
  REAL *restrict xx[3] = {{
      griddata[0].xx[0],
      griddata[0].xx[1],
      griddata[0].xx[2],
  }};
  const uint32_t format_version = {format_version}U;
  const uint32_t header_size = {header_size_bytes}U;
  const uint32_t output_index =
      raytracing_data_u32_from_nonnegative_int_or_abort(
          time_output_index, "time_output_index");
  const uint32_t num_grids = (uint32_t)commondata->NUMGRIDS;
  const uint32_t serialized_real_bytes = 8U;
  const uint32_t record_component_count = {num_record_components}U;
  const uint32_t metric_component_count = {num_metric_components}U;
  const uint32_t christoffel_component_count = {num_gamma_components}U;
  const uint32_t point_record_real_count = {point_record_real_count}U;
  const uint32_t point_record_bytes = {point_record_bytes}U;
  const uint32_t payload_includes_ghost_zones = 0U;
  const uint32_t file_is_little_endian = 1U;
  const uint32_t time_variable_is_f64 = 1U;
  const uint32_t reserved_u32 = 0U;
  const uint32_t Nxx0 =
      raytracing_data_u32_from_nonnegative_int_or_abort(params->Nxx0, "Nxx0");
  const uint32_t Nxx1 =
      raytracing_data_u32_from_nonnegative_int_or_abort(params->Nxx1, "Nxx1");
  const uint32_t Nxx2 =
      raytracing_data_u32_from_nonnegative_int_or_abort(params->Nxx2, "Nxx2");
  const uint32_t Nxx_plus_2NGHOSTS0_u32 =
      raytracing_data_u32_from_nonnegative_int_or_abort(
          params->Nxx_plus_2NGHOSTS0, "Nxx_plus_2NGHOSTS0");
  const uint32_t Nxx_plus_2NGHOSTS1_u32 =
      raytracing_data_u32_from_nonnegative_int_or_abort(
          params->Nxx_plus_2NGHOSTS1, "Nxx_plus_2NGHOSTS1");
  const uint32_t Nxx_plus_2NGHOSTS2_u32 =
      raytracing_data_u32_from_nonnegative_int_or_abort(
          params->Nxx_plus_2NGHOSTS2, "Nxx_plus_2NGHOSTS2");
  const uint64_t point_record_count = raytracing_data_mul_u64_or_abort(
      raytracing_data_mul_u64_or_abort((uint64_t)Nxx0, (uint64_t)Nxx1, "Nxx0*Nxx1"),
      (uint64_t)Nxx2,
      "interior point count");
  const uint64_t point_records_bytes = raytracing_data_mul_u64_or_abort(
      point_record_count,
      (uint64_t)point_record_bytes,
      "point-record payload size");
  const uint64_t simulation_time_offset = (uint64_t)header_size;
  const uint64_t point_records_offset = raytracing_data_add_u64_or_abort(
      simulation_time_offset, 8ULL, "point-record payload offset");
  const uint64_t total_file_bytes = raytracing_data_add_u64_or_abort(
      point_records_offset, point_records_bytes, "total file bytes");
  SET_NXX_PLUS_2NGHOSTS_VARS(0);
  const double dxx[3] = {{
      (double)params->dxx0,
      (double)params->dxx1,
      (double)params->dxx2,
  }};
  const double invdxx[3] = {{
      (double)params->invdxx0,
      (double)params->invdxx1,
      (double)params->invdxx2,
  }};
  const double xxmin[3] = {{
      (double)params->xxmin0,
      (double)params->xxmin1,
      (double)params->xxmin2,
  }};
  const double xxmax[3] = {{
      (double)params->xxmax0,
      (double)params->xxmax1,
      (double)params->xxmax2,
  }};
  const double cart_origin[3] = {{
      (double)params->Cart_originx,
      (double)params->Cart_originy,
      (double)params->Cart_originz,
  }};
  const char magic[16] = "NRPYRTDATA4D";

  char final_filename[256];
  char temporary_filename[256];
  snprintf(final_filename, sizeof(final_filename), "./raytracing_data_t%08d.bin", time_output_index);
  snprintf(
      temporary_filename,
      sizeof(temporary_filename),
      "./raytracing_data_t%08d.bin.tmp.XXXXXX",
      time_output_index);

  const int temporary_fd = mkstemp(temporary_filename);
  if (temporary_fd == -1) {{
    fprintf(stderr,
            "Error: output_raytracing_data could not create temporary file %s. errno=%d\n",
            temporary_filename,
            errno);
    exit(1);
  }} // END IF: temporary output file could not be created

  const mode_t current_umask = umask((mode_t)0);
  umask(current_umask);
  if (fchmod(temporary_fd, (mode_t)(0666 & ~current_umask)) != 0) {{
    const int saved_errno = errno;
    const int close_result = close(temporary_fd);
    const int close_errno = close_result == 0 ? 0 : errno;
    const int cleanup_result = remove(temporary_filename);
    const int cleanup_errno = cleanup_result == 0 ? 0 : errno;
    if (close_result != 0 || cleanup_result != 0) {{
      fprintf(stderr,
              "Error: output_raytracing_data could not adjust permissions on %s and cleanup was incomplete. fchmod errno=%d close result=%d close errno=%d cleanup result=%d cleanup errno=%d\n",
              temporary_filename,
              saved_errno,
              close_result,
              close_errno,
              cleanup_result,
              cleanup_errno);
      exit(1);
    }} // END IF: cleanup after fchmod failure failed
    fprintf(stderr,
            "Error: output_raytracing_data could not adjust permissions on %s. errno=%d\n",
            temporary_filename,
            saved_errno);
    exit(1);
  }} // END IF: temporary output file permissions could not be adjusted

  FILE *restrict fp = fdopen(temporary_fd, "wb");
  if (fp == NULL) {{
    const int saved_errno = errno;
    const int close_result = close(temporary_fd);
    const int close_errno = close_result == 0 ? 0 : errno;
    const int cleanup_result = remove(temporary_filename);
    const int cleanup_errno = cleanup_result == 0 ? 0 : errno;
    if (close_result != 0 || cleanup_result != 0) {{
      fprintf(stderr,
              "Error: output_raytracing_data could not open %s for writing and cleanup was incomplete. fdopen errno=%d close result=%d close errno=%d cleanup result=%d cleanup errno=%d\n",
              temporary_filename,
              saved_errno,
              close_result,
              close_errno,
              cleanup_result,
              cleanup_errno);
      exit(1);
    }} // END IF: cleanup after fdopen failure failed
    fprintf(stderr,
            "Error: output_raytracing_data could not open %s for writing. errno=%d\n",
            temporary_filename,
            saved_errno);
    exit(1);
  }} // END IF: fdopen failed

  if (access(final_filename, F_OK) == 0) {{
    fclose(fp);
    if (remove(temporary_filename) != 0) {{
      fprintf(stderr,
              "Error: output_raytracing_data found existing final file %s but could not remove temporary file %s. errno=%d\n",
              final_filename,
              temporary_filename,
              errno);
      exit(1);
    }} // END IF: temporary file cleanup failed after no-overwrite check
    return;
  }} // END IF: final output file already exists

  const int Nxx_plus_2NGHOSTS_tot =
      params->Nxx_plus_2NGHOSTS0 * params->Nxx_plus_2NGHOSTS1 * params->Nxx_plus_2NGHOSTS2;
  REAL *restrict rhs_gfs;
  BHAH_MALLOC(rhs_gfs, NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot * sizeof(REAL));
  if (rhs_gfs == NULL) {{
    fclose(fp);
    remove(temporary_filename);
    raytracing_data_abort_with_message(
        "Error: output_raytracing_data could not allocate temporary rhs_gfs storage.");
  }} // END IF: rhs scratch allocation failed

  Ricci_eval(params, rfmstruct, y_n_gfs, auxevol_gfs);
  rhs_eval(commondata, params, rfmstruct, auxevol_gfs, y_n_gfs, rhs_gfs);
"""
    body += r"""

  // Step 1: Write the fixed-size file header in documented little-endian field order.
  raytracing_data_write_or_abort(fp, magic, sizeof(char), 16, "magic");
  raytracing_data_write_u32_or_abort(fp, format_version, "format_version");
  raytracing_data_write_u32_or_abort(fp, header_size, "header_size");
  raytracing_data_write_u32_or_abort(fp, output_index, "output_index");
  raytracing_data_write_u32_or_abort(fp, num_grids, "num_grids");
  raytracing_data_write_u32_or_abort(fp, serialized_real_bytes, "serialized_real_bytes");
  raytracing_data_write_u32_or_abort(fp, record_component_count, "record_component_count");
  raytracing_data_write_u32_or_abort(fp, metric_component_count, "metric_component_count");
  raytracing_data_write_u32_or_abort(
      fp, christoffel_component_count, "christoffel_component_count");
  raytracing_data_write_u32_or_abort(fp, point_record_real_count, "point_record_real_count");
  raytracing_data_write_u32_or_abort(fp, point_record_bytes, "point_record_bytes");
  raytracing_data_write_u32_or_abort(
      fp, payload_includes_ghost_zones, "payload_includes_ghost_zones");
  raytracing_data_write_u32_or_abort(
      fp, file_is_little_endian, "file_is_little_endian");
  raytracing_data_write_u32_or_abort(
      fp, time_variable_is_f64, "time_variable_is_f64");
  raytracing_data_write_u32_or_abort(fp, reserved_u32, "reserved_u32");
  raytracing_data_write_u32_or_abort(fp, Nxx0, "Nxx0");
  raytracing_data_write_u32_or_abort(fp, Nxx1, "Nxx1");
  raytracing_data_write_u32_or_abort(fp, Nxx2, "Nxx2");
  raytracing_data_write_u32_or_abort(fp, Nxx_plus_2NGHOSTS0_u32, "Nxx_plus_2NGHOSTS0");
  raytracing_data_write_u32_or_abort(fp, Nxx_plus_2NGHOSTS1_u32, "Nxx_plus_2NGHOSTS1");
  raytracing_data_write_u32_or_abort(fp, Nxx_plus_2NGHOSTS2_u32, "Nxx_plus_2NGHOSTS2");
  raytracing_data_write_u64_or_abort(fp, point_record_count, "point_record_count");
  raytracing_data_write_u64_or_abort(fp, simulation_time_offset, "simulation_time_offset");
  raytracing_data_write_u64_or_abort(fp, point_records_offset, "point_records_offset");
  raytracing_data_write_u64_or_abort(fp, point_records_bytes, "point_records_bytes");
  raytracing_data_write_u64_or_abort(fp, total_file_bytes, "total_file_bytes");
  raytracing_data_write_u64_or_abort(fp, (uint64_t)NGHOSTS, "NGHOSTS");
  raytracing_data_write_u64_or_abort(fp, (uint64_t)Nxx0, "payload_i0_count");
  raytracing_data_write_u64_or_abort(fp, (uint64_t)Nxx1, "payload_i1_count");
  raytracing_data_write_u64_or_abort(fp, (uint64_t)Nxx2, "payload_i2_count");
  raytracing_data_write_u64_or_abort(fp, (uint64_t)NGHOSTS, "payload_i0_start");
  raytracing_data_write_u64_or_abort(fp, (uint64_t)NGHOSTS, "payload_i1_start");
  raytracing_data_write_u64_or_abort(fp, (uint64_t)NGHOSTS, "payload_i2_start");
  raytracing_data_write_u64_or_abort(
      fp, (uint64_t)(params->Nxx_plus_2NGHOSTS0 - NGHOSTS), "payload_i0_end");
  raytracing_data_write_u64_or_abort(
      fp, (uint64_t)(params->Nxx_plus_2NGHOSTS1 - NGHOSTS), "payload_i1_end");
  raytracing_data_write_u64_or_abort(
      fp, (uint64_t)(params->Nxx_plus_2NGHOSTS2 - NGHOSTS), "payload_i2_end");

  for (int i = 0; i < 3; i++)
    raytracing_data_write_f64_or_abort(fp, dxx[i], "dxx");
  for (int i = 0; i < 3; i++)
    raytracing_data_write_f64_or_abort(fp, invdxx[i], "invdxx");
  for (int i = 0; i < 3; i++)
    raytracing_data_write_f64_or_abort(fp, xxmin[i], "xxmin");
  for (int i = 0; i < 3; i++)
    raytracing_data_write_f64_or_abort(fp, xxmax[i], "xxmax");
  for (int i = 0; i < 3; i++)
    raytracing_data_write_f64_or_abort(fp, cart_origin[i], "cart_origin");

"""
    body += (
        "  raytracing_data_write_fixed_length_string_or_abort("
        f'fp, "{header_label}", {header_label_bytes}, "header_label");\n'
    )
    body += (
        "  raytracing_data_write_fixed_length_string_or_abort("
        f'fp, "{format_name}", {format_name_bytes}, "format_name");\n'
    )
    body += (
        "  raytracing_data_write_fixed_length_string_or_abort("
        f'fp, "{target_basis_name}", {basis_name_bytes}, "target_basis");\n'
    )
    body += (
        f'  raytracing_data_write_fixed_length_string_or_abort(fp, "{CoordSystem}",'
        f' {coord_name_bytes}, "source_coord_system");\n'
    )
    body += (
        "  raytracing_data_write_fixed_length_string_or_abort("
        f'fp, "{loop_order_name}", {loop_order_name_bytes}, "loop_order");\n'
    )
    body += (
        "  raytracing_data_write_fixed_length_string_or_abort("
        f'fp, "{cf_convention}", {cf_convention_bytes}, "cf_convention");\n'
    )
    for component_name in record_component_names:
        body += (
            "  raytracing_data_write_fixed_length_string_or_abort("
            f'fp, "{component_name}", {record_component_name_bytes}, "record_component_names");\n'
        )
    for _, component_name in metric_components:
        body += (
            "  raytracing_data_write_fixed_length_string_or_abort("
            f'fp, "{component_name}", {metric_component_name_bytes}, "metric_component_names");\n'
        )
    for _, component_name in christoffel_components:
        body += (
            "  raytracing_data_write_fixed_length_string_or_abort("
            f'fp, "{component_name}", {gamma_component_name_bytes}, "christoffel_component_names");\n'
        )
    body += r"""

  // Step 2: Write the physical simulation time.
  raytracing_data_write_f64_or_abort(fp, (double)commondata->time, "simulation_time");

  // Step 3: Evaluate Cartesian coordinates, metric, and Christoffels on the interior grid.
"""
    body += loop
    body += r"""

  // Step 4: Finalize the file and install it without overwriting an existing output.
  BHAH_FREE(rhs_gfs);

  if (fflush(fp) != 0) {
    const int saved_errno = errno;
    fclose(fp);
    remove(temporary_filename);
    fprintf(stderr,
            "Error: output_raytracing_data failed while flushing %s. errno=%d\n",
            temporary_filename,
            saved_errno);
    exit(1);
  } // END IF: fflush failed
  if (fclose(fp) != 0) {
    const int saved_errno = errno;
    remove(temporary_filename);
    fprintf(stderr,
            "Error: output_raytracing_data failed while closing %s. errno=%d\n",
            temporary_filename,
            saved_errno);
    exit(1);
  } // END IF: fclose failed
  if (rename(temporary_filename, final_filename) != 0) {
    const int saved_errno = errno;
    remove(temporary_filename);
    fprintf(stderr,
            "Error: output_raytracing_data could not install %s at %s. errno=%d\n",
            temporary_filename,
            final_filename,
            saved_errno);
    exit(1);
  } // END IF: rename failed
"""

    cfc.register_CFunction(
        subdirectory="diagnostics",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
