"""
Register C functions that export raytracing-connection 4D reconstruction data.

This module generates a C helper that writes exactly one binary file per
diagnostics output time. Each file contains:

1. A fixed-size self-describing binary header.
2. The physical simulation time serialized as little-endian IEEE-754 binary64.
3. The logical coordinate arrays ``xx[0..2]`` serialized as little-endian
   IEEE-754 binary64 values.
4. Only the evolved BSSN gridfunctions needed to reconstruct the spatial BSSN
   state consumed by
   ``nrpy.equations.general_relativity.geodesics.symbolic_christoffel_recipe_from_bssn_grid_basis()``,
   serialized as little-endian IEEE-754 binary64 values.

The fixed-size header is designed so that a later post-processing or
interpolation reader does not need to know the in-memory C layout of
``params_struct`` to interpret the payload layout. It stores explicit byte
offsets, per-array lengths, per-array byte sizes, the canonical on-disk
floating-point format, the evolved conformal-factor convention, the
reference-metric precompute mode, coordinate-system labels, and
interpolation-critical grid metadata such as interior and ghosted extents,
coordinate minima/maxima, grid spacings, and coordinate/gridfunction naming.

The portable payload intentionally does not serialize raw ``params_struct``
bytes. All reader-required grid metadata is written explicitly in the fixed
header. If a future reader requires additional runtime parameters, they should
be added to the header or to a documented field-by-field canonical params block,
not by writing native C struct bytes.

The emitted payload intentionally excludes evolved fields that are not consumed
by the Christoffel construction. In particular, it does not write ``betU`` or
``lambdaU`` gridfunctions, because the current ``BSSN_to_g4Christoffel``
construction depends only on ``alpha``, ``cf``, ``trK``, ``vetU``, ``hDD``,
and ``aDD`` plus spatial coordinates and grid parameters. Spatial derivatives
are not exported directly; the file contains the state needed to reconstruct
them later by interpolation and finite differencing. Likewise, reference-metric
precompute helpers such as ``f0_of_xx0`` and ``f1_of_xx1`` are not serialized
as payload arrays because they are deterministic functions of the coordinate
system, the stored coordinate arrays, and the stored runtime parameters.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Tuple, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par


def _raytracing_connections_required_gridfunctions() -> List[Tuple[str, str]]:
    """
    Return the exact evolved BSSN gridfunctions needed for Christoffel reconstruction.

    The list order defines the on-disk payload order for every exported
    ``raytracing_connections_4d_data_t########.bin`` file.

    :return: Ordered ``(macro_name, human_readable_name)`` pairs for each exported field.
    """
    return [
        ("ADD00GF", "aDD00"),
        ("ADD01GF", "aDD01"),
        ("ADD02GF", "aDD02"),
        ("ADD11GF", "aDD11"),
        ("ADD12GF", "aDD12"),
        ("ADD22GF", "aDD22"),
        ("ALPHAGF", "alpha"),
        ("CFGF", "cf"),
        ("HDD00GF", "hDD00"),
        ("HDD01GF", "hDD01"),
        ("HDD02GF", "hDD02"),
        ("HDD11GF", "hDD11"),
        ("HDD12GF", "hDD12"),
        ("HDD22GF", "hDD22"),
        ("TRKGF", "trK"),
        ("VETU0GF", "vetU0"),
        ("VETU1GF", "vetU1"),
        ("VETU2GF", "vetU2"),
    ]


def output_raytracing_connections_4d_data(
    enable_rfm_precompute: bool,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct and register the raytracing-connections 4D binary exporter.

    The generated C helper writes one binary file per diagnostics output with
    filename pattern ``./raytracing_connections_4d_data_t########.bin``. The
    file is first written to a temporary sibling path and renamed into place
    only after all bytes are successfully flushed. If the final file already
    exists, the exporter returns immediately so restart/checkpoint workflows do
    not overwrite or append duplicate data.

    File contents, in order:
    1. A fixed-size binary header describing the payload layout, byte offsets,
       array lengths, array sizes, canonical floating-point precision, byte
       order, the evolved conformal-factor convention, the reference-metric
       precompute mode, and the interpolation-critical grid metadata duplicated
       from ``griddata[0].params``.
    2. The physical simulation time ``commondata->time`` serialized as
       little-endian IEEE-754 binary64.
    3. The logical coordinate arrays ``xx[0]``, ``xx[1]``, and ``xx[2]``
       serialized as little-endian IEEE-754 binary64.
    4. The exact required BSSN gridfunctions, one full gridfunction array at a
       time, serialized as little-endian IEEE-754 binary64.

    One output file contains the spatial BSSN state and metadata needed for
    later interpolation/finite-difference reconstruction. Time derivatives such
    as ``alpha_d0`` and ``vetU_d0`` are not directly stored in a single file.

    The portable payload intentionally does not serialize raw ``params_struct``
    bytes. All reader-required grid metadata is written explicitly in the fixed
    header. If a future reader requires a complete params snapshot, add a
    documented field-by-field canonical params serializer rather than writing
    native C struct bytes.

    :param enable_rfm_precompute: Whether the generated project reconstructs
        reference-metric helpers from the precompute pathway.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If ``enable_rfm_precompute`` is not a bool.

    Doctests:
    None.
    """
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

    required_gfs = _raytracing_connections_required_gridfunctions()
    num_export_gfs = len(required_gfs)
    gf_name_buffer_size = 24
    coord_name_buffer_size = 8
    type_name_buffer_size = 32
    header_scalar_real_count = 19
    format_version = 5
    header_size_bytes = (
        16
        + 20 * 4
        + (15 + num_export_gfs) * 8
        + 16
        + 16
        + type_name_buffer_size
        + 3 * coord_name_buffer_size
        + 100
        + 100
        + 5 * 128
        + header_scalar_real_count * 8
        + 3
        + 5
        + num_export_gfs * gf_name_buffer_size
    )
    cf_convention = cast(str, par.parval_from_str("EvolvedConformalFactor_cf"))

    includes = [
        "<errno.h>",
        "<float.h>",
        "<stdint.h>",
        "<stdio.h>",
        "<stdlib.h>",
        "<string.h>",
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = """
Export raytracing-connection 4D reconstruction data from the stored BSSN state.

This function writes at most one binary file per diagnostics output to:

  ./raytracing_connections_4d_data_t########.bin

where the eight-digit index is the diagnostics output index, not the raw
timestep counter. Each file is written to a temporary sibling path and then
atomically renamed into place after all bytes are flushed successfully. If
the final file already exists, the function returns without modifying it so
restarted runs do not overwrite or append duplicate output.

The exported payload is intentionally limited to the quantities needed to
reconstruct the spatial BSSN state used later for 4D Christoffel
reconstruction:
  - physical simulation time serialized as little-endian IEEE-754 binary64,
  - duplicated interpolation-critical grid metadata in the fixed header,
  - logical coordinate arrays xx[0], xx[1], xx[2] serialized as little-endian
    IEEE-754 binary64,
  - aDD00, aDD01, aDD02, aDD11, aDD12, aDD22,
  - alpha,
  - cf,
  - hDD00, hDD01, hDD02, hDD11, hDD12, hDD22,
  - trK,
  - vetU0, vetU1, vetU2.

The fixed header records byte offsets and byte counts for every payload
section, the EvolvedConformalFactor_cf convention used to interpret the
stored cf field, and whether reference-metric precompute was enabled when
the project was generated. A later reader can therefore parse the coordinate
arrays and each stored gridfunction deterministically without decoding any
opaque native C struct bytes.

The portable payload intentionally does not serialize raw params_struct bytes.
All reader-required grid metadata is written explicitly in the fixed header.
If a future reader requires a complete params snapshot, add a documented
field-by-field canonical params serializer rather than writing native C struct
bytes.

Fields not used by the current BSSN-to-Christoffel construction, such as
betU and lambdaU, are deliberately omitted. Reference-metric helper arrays
such as f0_of_xx0 are also omitted because they are reconstructed
analytically from CoordSystemName, xx[0..2], and the stored runtime
parameters instead of being treated as independent simulation data.

@param[in] commondata  Global simulation metadata used in the file header and output naming
@param[in] griddata  Host-side per-grid runtime state to export
@param[in] time_output_index  Diagnostics output index used in the emitted filename
"""
    cfunc_type = "void"
    name = "output_raytracing_connections_4d_data"
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
static void raytracing_connections_abort_with_message(const char *restrict message) {
  fprintf(stderr, "%s\n", message);
  exit(1);
} // END FUNCTION: raytracing_connections_abort_with_message

/**
 * Cast a nonnegative integer to uint32_t, aborting if negative or if it overflows.
 *
 * @param[in] value  Integer value to convert.
 * @param[in] label  Label for error message.
 * @return Safely casted uint32_t.
 */
static uint32_t raytracing_connections_u32_from_nonnegative_int_or_abort(
    const int value,
    const char *restrict label
) {
  if (value < 0) {
    fprintf(stderr,
            "Error: output_raytracing_connections_4d_data received negative %s=%d.\n",
            label,
            value);
    exit(1);
  } // END IF: negative value cannot be serialized as uint32_t
  if ((uint64_t)value > (uint64_t)UINT32_MAX) {
    fprintf(stderr,
            "Error: output_raytracing_connections_4d_data %s=%d does not fit in uint32_t.\n",
            label,
            value);
    exit(1);
  } // END IF: value cannot be represented in the documented header format
  return (uint32_t)value;
} // END FUNCTION: raytracing_connections_u32_from_nonnegative_int_or_abort

/**
 * Add two uint64_t values, aborting on overflow.
 *
 * @param[in] a      First value.
 * @param[in] b      Second value.
 * @param[in] label  Label for error message.
 * @return Sum of a and b.
 */
static uint64_t raytracing_connections_add_u64_or_abort(
    const uint64_t a,
    const uint64_t b,
    const char *restrict label
) {
  if (UINT64_MAX - a < b) {
    fprintf(stderr,
            "Error: output_raytracing_connections_4d_data uint64 overflow while computing %s.\n",
            label);
    exit(1);
  } // END IF: addition would overflow
  return a + b;
} // END FUNCTION: raytracing_connections_add_u64_or_abort

/**
 * Multiply two uint64_t values, aborting on overflow.
 *
 * @param[in] a      First value.
 * @param[in] b      Second value.
 * @param[in] label  Label for error message.
 * @return Product of a and b.
 */
static uint64_t raytracing_connections_mul_u64_or_abort(
    const uint64_t a,
    const uint64_t b,
    const char *restrict label
) {
  if (a != 0ULL && b > UINT64_MAX / a) {
    fprintf(stderr,
            "Error: output_raytracing_connections_4d_data uint64 overflow while computing %s.\n",
            label);
    exit(1);
  } // END IF: multiplication would overflow
  return a * b;
} // END FUNCTION: raytracing_connections_mul_u64_or_abort

/**
 * Write elements to a file, aborting on failure.
 *
 * @param[in,out] fp            File pointer.
 * @param[in]     data          Data to write.
 * @param[in]     element_size  Size of each element.
 * @param[in]     count         Number of elements.
 * @param[in]     label         Label for error message.
 */
static void raytracing_connections_write_or_abort(
    FILE *restrict fp,
    const void *restrict data,
    const size_t element_size,
    const size_t count,
    const char *restrict label
) {
  if (count != 0 && element_size > SIZE_MAX / count) {
    fprintf(stderr,
            "Error: output_raytracing_connections_4d_data size_t overflow before writing %s.\n",
            label);
    fclose(fp);
    exit(1);
  } // END IF: fwrite byte count would overflow size_t
  if (fwrite(data, element_size, count, fp) != count) {
    fprintf(stderr,
            "Error: output_raytracing_connections_4d_data failed while writing %s.\n",
            label);
    fclose(fp);
    exit(1);
  } // END IF: fwrite failed
} // END FUNCTION: raytracing_connections_write_or_abort

/**
 * Write a uint32_t value as little-endian, aborting on failure.
 *
 * @param[in,out] fp     File pointer.
 * @param[in]     value  Value to write.
 * @param[in]     label  Label for error message.
 */
static void raytracing_connections_write_u32_or_abort(
    FILE *restrict fp,
    const uint32_t value,
    const char *restrict label
) {
  uint8_t buffer[4];
  buffer[0] = (uint8_t)(value & 0xffU);
  buffer[1] = (uint8_t)((value >> 8U) & 0xffU);
  buffer[2] = (uint8_t)((value >> 16U) & 0xffU);
  buffer[3] = (uint8_t)((value >> 24U) & 0xffU);
  raytracing_connections_write_or_abort(fp, buffer, sizeof(uint8_t), 4, label);
} // END FUNCTION: raytracing_connections_write_u32_or_abort

/**
 * Write a uint64_t value as little-endian, aborting on failure.
 *
 * @param[in,out] fp     File pointer.
 * @param[in]     value  Value to write.
 * @param[in]     label  Label for error message.
 */
static void raytracing_connections_write_u64_or_abort(
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
  raytracing_connections_write_or_abort(fp, buffer, sizeof(uint8_t), 8, label);
} // END FUNCTION: raytracing_connections_write_u64_or_abort

/**
 * Write a double value as little-endian binary64, aborting on failure.
 *
 * @param[in,out] fp     File pointer.
 * @param[in]     value  Value to write.
 * @param[in]     label  Label for error message.
 */
static void raytracing_connections_write_f64_or_abort(
    FILE *restrict fp,
    const double value,
    const char *restrict label
) {
  uint64_t bits = 0ULL;
  memcpy(&bits, &value, sizeof(bits));
  raytracing_connections_write_u64_or_abort(fp, bits, label);
} // END FUNCTION: raytracing_connections_write_f64_or_abort

/**
 * Validate that the current C compiler's double matches IEEE-754 binary64 semantics,
 * and REAL is not wider than double. Abort if unsupported.
 */
static void raytracing_connections_validate_binary64_output_or_abort(void) {
#if FLT_RADIX != 2 || DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024 || DBL_MIN_EXP != -1021
  raytracing_connections_abort_with_message(
      "Error: raytracing output requires C double to match IEEE-754 binary64 semantics.");
#endif
  if (sizeof(double) != 8) {
    raytracing_connections_abort_with_message(
        "Error: raytracing output requires 8-byte C double for binary64 serialization.");
  } // END IF: unsupported double size
  if (sizeof(REAL) > sizeof(double)) {
    raytracing_connections_abort_with_message(
        "Error: raytracing output refuses to silently down-convert REAL values wider than binary64.");
  } // END IF: REAL would lose precision in the documented on-disk format
} // END FUNCTION: raytracing_connections_validate_binary64_output_or_abort

/**
 * Write an array of REAL values as little-endian binary64, aborting on failure.
 *
 * @param[in,out] fp      File pointer.
 * @param[in]     values  Array of REAL values.
 * @param[in]     count   Number of values.
 * @param[in]     label   Label for error message.
 */
static void raytracing_connections_write_real_array_as_f64_or_abort(
    FILE *restrict fp,
    const REAL *restrict values,
    const uint64_t count,
    const char *restrict label
) {
  for (uint64_t idx = 0; idx < count; idx++) {
    raytracing_connections_write_f64_or_abort(fp, (double)values[idx], label);
  } // END LOOP: write REAL array as canonical little-endian binary64
} // END FUNCTION: raytracing_connections_write_real_array_as_f64_or_abort

/**
 * Write a fixed-width string (null-padded), aborting if too long.
 *
 * @param[in,out] fp           File pointer.
 * @param[in]     text         String to write.
 * @param[in]     field_bytes  Exact number of bytes to write.
 * @param[in]     label        Label for error message.
 */
static void raytracing_connections_write_fixed_length_string_or_abort(
    FILE *restrict fp,
    const char *restrict text,
    const size_t field_bytes,
    const char *restrict label
) {
  char buffer[128];
  const size_t text_len = strlen(text);
  if (field_bytes > sizeof(buffer)) {
    raytracing_connections_abort_with_message(
        "Error: raytracing header string field exceeds the internal serialization buffer.");
  } // END IF: unsupported fixed-width string size
  if (field_bytes == 0 || text_len >= field_bytes) {
    fprintf(stderr,
            "Error: output_raytracing_connections_4d_data string field %s is too long for its fixed header field.\n",
            label);
    fclose(fp);
    exit(1);
  } // END IF: fixed-width string would be truncated
  memset(buffer, 0, sizeof(buffer));
  memcpy(buffer, text, text_len);
  raytracing_connections_write_or_abort(fp, buffer, sizeof(char), field_bytes, label);
} // END FUNCTION: raytracing_connections_write_fixed_length_string_or_abort

/**
 * Assert that the file pointer is at the expected header size, aborting if not.
 *
 * @param[in,out] fp                    File pointer.
 * @param[in]     expected_header_size  Expected byte offset.
 */
static void raytracing_connections_assert_header_size_or_abort(
    FILE *restrict fp,
    const uint32_t expected_header_size
) {
  const long current_offset = ftell(fp);
  if (current_offset < 0 || (uint64_t)current_offset != (uint64_t)expected_header_size) {
    fprintf(stderr,
            "Error: output_raytracing_connections_4d_data wrote %ld header bytes, expected %u.\n",
            current_offset,
            expected_header_size);
    fclose(fp);
    exit(1);
  } // END IF: hand-maintained header byte count is inconsistent
} // END FUNCTION: raytracing_connections_assert_header_size_or_abort

/**
 * Check if a file exists and is readable.
 *
 * @param[in] path  File path to check.
 * @return 1 if exists, 0 otherwise.
 */
static int raytracing_connections_file_exists(const char *restrict path) {
  FILE *restrict fp = fopen(path, "rb");
  if (fp != NULL) {
    fclose(fp);
    return 1;
  } // END IF: existing file found
  return 0;
} // END FUNCTION: raytracing_connections_file_exists
"""

    body = rf"""
  if (commondata->NUMGRIDS != 1) {{
    raytracing_connections_abort_with_message(
        "Error: output_raytracing_connections_4d_data currently requires commondata->NUMGRIDS == 1.");
  }} // END IF: unsupported number of grids

  if (time_output_index < 0) {{
    raytracing_connections_abort_with_message(
        "Error: output_raytracing_connections_4d_data received a negative output index.");
  }} // END IF: invalid output index

  raytracing_connections_validate_binary64_output_or_abort();

  const params_struct *restrict params = &griddata[0].params;
  const uint32_t format_version = {format_version}U;
  const uint32_t header_size = {header_size_bytes}U;
  const uint32_t output_index =
      raytracing_connections_u32_from_nonnegative_int_or_abort(
          time_output_index, "time_output_index");
  const uint32_t num_grids = (uint32_t)commondata->NUMGRIDS;
  const uint32_t sizeof_real = 8U;
  const uint32_t sizeof_params = 0U;
  const uint32_t coord_array_count = 3U;
  const uint32_t num_export_gridfunctions = {num_export_gfs}U;
  const uint32_t gridfunction_name_bytes = {gf_name_buffer_size}U;
  const uint32_t coord_name_bytes = {coord_name_buffer_size}U;
  const uint32_t grid_idx =
      raytracing_connections_u32_from_nonnegative_int_or_abort(params->grid_idx, "grid_idx");
  const uint32_t file_is_little_endian = 1U;
  const uint32_t time_variable_is_f64 = 1U;
  const uint32_t reserved_u32 = 0U;
  const uint32_t Nxx0 =
      raytracing_connections_u32_from_nonnegative_int_or_abort(params->Nxx0, "Nxx0");
  const uint32_t Nxx1 =
      raytracing_connections_u32_from_nonnegative_int_or_abort(params->Nxx1, "Nxx1");
  const uint32_t Nxx2 =
      raytracing_connections_u32_from_nonnegative_int_or_abort(params->Nxx2, "Nxx2");
  const uint32_t Nxx_plus_2NGHOSTS0_u32 =
      raytracing_connections_u32_from_nonnegative_int_or_abort(
          params->Nxx_plus_2NGHOSTS0, "Nxx_plus_2NGHOSTS0");
  const uint32_t Nxx_plus_2NGHOSTS1_u32 =
      raytracing_connections_u32_from_nonnegative_int_or_abort(
          params->Nxx_plus_2NGHOSTS1, "Nxx_plus_2NGHOSTS1");
  const uint32_t Nxx_plus_2NGHOSTS2_u32 =
      raytracing_connections_u32_from_nonnegative_int_or_abort(
          params->Nxx_plus_2NGHOSTS2, "Nxx_plus_2NGHOSTS2");
  const uint64_t ntot01 = raytracing_connections_mul_u64_or_abort(
      (uint64_t)Nxx_plus_2NGHOSTS0_u32,
      (uint64_t)Nxx_plus_2NGHOSTS1_u32,
      "Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1");
  const uint64_t ntot = raytracing_connections_mul_u64_or_abort(
      ntot01,
      (uint64_t)Nxx_plus_2NGHOSTS2_u32,
      "total grid point count");
  const uint32_t serialized_real_bytes = sizeof_real;
  const uint64_t simulation_time_bytes = 8ULL;
  const uint64_t params_struct_bytes = 0ULL;
  const uint64_t coord_array_bytes0 = raytracing_connections_mul_u64_or_abort(
      (uint64_t)Nxx_plus_2NGHOSTS0_u32,
      (uint64_t)serialized_real_bytes,
      "xx[0] byte count");
  const uint64_t coord_array_bytes1 = raytracing_connections_mul_u64_or_abort(
      (uint64_t)Nxx_plus_2NGHOSTS1_u32,
      (uint64_t)serialized_real_bytes,
      "xx[1] byte count");
  const uint64_t coord_array_bytes2 = raytracing_connections_mul_u64_or_abort(
      (uint64_t)Nxx_plus_2NGHOSTS2_u32,
      (uint64_t)serialized_real_bytes,
      "xx[2] byte count");
  const uint64_t gridfunction_bytes_per_array = raytracing_connections_mul_u64_or_abort(
      ntot,
      (uint64_t)serialized_real_bytes,
      "gridfunction byte count");

  char final_filename[256];
  char temporary_filename[256];
  snprintf(
      final_filename,
      sizeof(final_filename),
      "./raytracing_connections_4d_data_t%08d.bin",
      time_output_index);
  snprintf(
      temporary_filename,
      sizeof(temporary_filename),
      "./raytracing_connections_4d_data_t%08d.bin.tmp",
      time_output_index);

  if (raytracing_connections_file_exists(final_filename)) {{
    return;
  }} // END IF: output file already exists

  remove(temporary_filename);

  FILE *restrict fp = fopen(temporary_filename, "wb");
  if (fp == NULL) {{
    fprintf(stderr,
            "Error: output_raytracing_connections_4d_data could not open %s for writing. errno=%d\n",
            temporary_filename, errno);
    exit(1);
  }} // END IF: fopen failed

  const uint64_t simulation_time_offset = (uint64_t)header_size;
  const uint64_t params_struct_offset = 0ULL;
  const uint64_t first_payload_offset = raytracing_connections_add_u64_or_abort(
      raytracing_connections_add_u64_or_abort(
          simulation_time_offset, simulation_time_bytes, "simulation-time payload end"),
      params_struct_bytes,
      "first coordinate-array offset");
  const uint64_t coord_array_offsets[3] = {{
      first_payload_offset,
      raytracing_connections_add_u64_or_abort(
          first_payload_offset, coord_array_bytes0, "xx[1] offset"),
      raytracing_connections_add_u64_or_abort(
          raytracing_connections_add_u64_or_abort(
              first_payload_offset, coord_array_bytes0, "temporary xx[2] offset"),
          coord_array_bytes1,
          "xx[2] offset"),
  }};
  const uint64_t gridfunction_payload_offset = raytracing_connections_add_u64_or_abort(
      coord_array_offsets[2], coord_array_bytes2, "gridfunction payload offset");
  const uint64_t total_gridfunction_payload_bytes = raytracing_connections_mul_u64_or_abort(
      (uint64_t){num_export_gfs}, gridfunction_bytes_per_array, "all gridfunction payload bytes");
  const uint64_t total_file_bytes = raytracing_connections_add_u64_or_abort(
      gridfunction_payload_offset, total_gridfunction_payload_bytes, "total file bytes");
  const uint64_t coord_array_lengths[3] = {{
      (uint64_t)Nxx_plus_2NGHOSTS0_u32,
      (uint64_t)Nxx_plus_2NGHOSTS1_u32,
      (uint64_t)Nxx_plus_2NGHOSTS2_u32,
  }};
  const uint64_t coord_array_bytes[3] = {{
      coord_array_bytes0,
      coord_array_bytes1,
      coord_array_bytes2,
  }};
  uint64_t gridfunction_offsets[{num_export_gfs}];
  const double cart_origin[3] = {{
      (double)params->Cart_originx,
      (double)params->Cart_originy,
      (double)params->Cart_originz,
  }};
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
  const double PI = (double)params->PI;
  const double RMAX = (double)params->RMAX;
  const double grid_hole_radius = (double)params->grid_hole_radius;
  const double grid_physical_size = (double)params->grid_physical_size;
  const int8_t enable_rfm_precompute_i8 = {1 if enable_rfm_precompute else 0};
  const int8_t grid_rotates = params->grid_rotates ? 1 : 0;
  const int8_t is_host = params->is_host ? 1 : 0;
  const int8_t reserved_flags[5] = {{0, 0, 0, 0, 0}};
  const char magic[16] = "NRPYRTCONN4D";

  for (size_t gf_index = 0; gf_index < (size_t){num_export_gfs}; gf_index++) {{
    gridfunction_offsets[gf_index] = raytracing_connections_add_u64_or_abort(
        gridfunction_payload_offset,
        raytracing_connections_mul_u64_or_abort(
            (uint64_t)gf_index, gridfunction_bytes_per_array, "gridfunction offset stride"),
        "gridfunction offset");
  }} // END LOOP: assign gridfunction payload offsets
"""
    body += r"""

  // Step 1: Write the fixed-size file header in documented little-endian field order.
  raytracing_connections_write_or_abort(fp, magic, sizeof(char), 16, "magic");
  raytracing_connections_write_u32_or_abort(fp, format_version, "format_version");
  raytracing_connections_write_u32_or_abort(fp, header_size, "header_size");
  raytracing_connections_write_u32_or_abort(fp, output_index, "output_index");
  raytracing_connections_write_u32_or_abort(fp, num_grids, "num_grids");
  raytracing_connections_write_u32_or_abort(fp, sizeof_real, "serialized_real_bytes");
  raytracing_connections_write_u32_or_abort(fp, sizeof_params, "sizeof_params_struct");
  raytracing_connections_write_u32_or_abort(fp, coord_array_count, "coord_array_count");
  raytracing_connections_write_u32_or_abort(
      fp, num_export_gridfunctions, "num_export_gridfunctions");
  raytracing_connections_write_u32_or_abort(
      fp, gridfunction_name_bytes, "gridfunction_name_bytes");
  raytracing_connections_write_u32_or_abort(fp, coord_name_bytes, "coord_name_bytes");
  raytracing_connections_write_u32_or_abort(fp, grid_idx, "grid_idx");
  raytracing_connections_write_u32_or_abort(
      fp, file_is_little_endian, "file_is_little_endian");
  raytracing_connections_write_u32_or_abort(
      fp, time_variable_is_f64, "time_variable_is_f64");
  raytracing_connections_write_u32_or_abort(fp, reserved_u32, "reserved_u32");
  raytracing_connections_write_u32_or_abort(fp, Nxx0, "Nxx0");
  raytracing_connections_write_u32_or_abort(fp, Nxx1, "Nxx1");
  raytracing_connections_write_u32_or_abort(fp, Nxx2, "Nxx2");
  raytracing_connections_write_u32_or_abort(
      fp, Nxx_plus_2NGHOSTS0_u32, "Nxx_plus_2NGHOSTS0");
  raytracing_connections_write_u32_or_abort(
      fp, Nxx_plus_2NGHOSTS1_u32, "Nxx_plus_2NGHOSTS1");
  raytracing_connections_write_u32_or_abort(
      fp, Nxx_plus_2NGHOSTS2_u32, "Nxx_plus_2NGHOSTS2");
  raytracing_connections_write_u64_or_abort(fp, ntot, "ntot");
  raytracing_connections_write_u64_or_abort(
      fp, simulation_time_offset, "simulation_time_offset");
  raytracing_connections_write_u64_or_abort(
      fp, params_struct_offset, "params_struct_offset");
  raytracing_connections_write_u64_or_abort(
      fp, coord_array_offsets[0], "coord_array_offsets[0]");
  raytracing_connections_write_u64_or_abort(
      fp, coord_array_offsets[1], "coord_array_offsets[1]");
  raytracing_connections_write_u64_or_abort(
      fp, coord_array_offsets[2], "coord_array_offsets[2]");
  raytracing_connections_write_u64_or_abort(
      fp, gridfunction_payload_offset, "gridfunction_payload_offset");
  raytracing_connections_write_u64_or_abort(fp, total_file_bytes, "total_file_bytes");
  raytracing_connections_write_u64_or_abort(
      fp, coord_array_lengths[0], "coord_array_lengths[0]");
  raytracing_connections_write_u64_or_abort(
      fp, coord_array_lengths[1], "coord_array_lengths[1]");
  raytracing_connections_write_u64_or_abort(
      fp, coord_array_lengths[2], "coord_array_lengths[2]");
  raytracing_connections_write_u64_or_abort(
      fp, coord_array_bytes[0], "coord_array_bytes[0]");
  raytracing_connections_write_u64_or_abort(
      fp, coord_array_bytes[1], "coord_array_bytes[1]");
  raytracing_connections_write_u64_or_abort(
      fp, coord_array_bytes[2], "coord_array_bytes[2]");
  raytracing_connections_write_u64_or_abort(
      fp, gridfunction_bytes_per_array, "gridfunction_bytes_per_array");
"""
    for gf_index, _ in enumerate(required_gfs):
        body += rf"""
  raytracing_connections_write_u64_or_abort(
      fp, gridfunction_offsets[{gf_index}], "gridfunction_offsets[{gf_index}]");
"""
    body += rf"""
  raytracing_connections_write_fixed_length_string_or_abort(
      fp, "binary64", 16, "real_type");
  raytracing_connections_write_fixed_length_string_or_abort(
      fp, "{cf_convention}", 16, "conformal_factor_cf");
  raytracing_connections_write_fixed_length_string_or_abort(
      fp, "not_serialized", {type_name_buffer_size}, "params_struct_type");
  raytracing_connections_write_fixed_length_string_or_abort(
      fp, "xx0", {coord_name_buffer_size}, "coord_names[0]");
  raytracing_connections_write_fixed_length_string_or_abort(
      fp, "xx1", {coord_name_buffer_size}, "coord_names[1]");
  raytracing_connections_write_fixed_length_string_or_abort(
      fp, "xx2", {coord_name_buffer_size}, "coord_names[2]");
  raytracing_connections_write_fixed_length_string_or_abort(
      fp, params->CoordSystemName, 100, "coord_system_name");
  raytracing_connections_write_fixed_length_string_or_abort(
      fp, params->gridname, 100, "grid_name");
  raytracing_connections_write_fixed_length_string_or_abort(
      fp,
      "Fixed header and payload values are serialized field-by-field in little-endian order.",
      128,
      "header_description");
  raytracing_connections_write_fixed_length_string_or_abort(
      fp,
      "No raw params_struct blob is serialized; required grid metadata is explicit in this header.",
      128,
      "params_payload_description");
  raytracing_connections_write_fixed_length_string_or_abort(
      fp,
      "Three logical coordinate arrays xx0, xx1, xx2 follow as little-endian binary64 values.",
      128,
      "coord_payload_description");
  raytracing_connections_write_fixed_length_string_or_abort(
      fp,
      "Each listed gridfunction follows as little-endian IEEE-754 binary64 values.",
      128,
      "gridfunction_payload_description");
  raytracing_connections_write_fixed_length_string_or_abort(
      fp,
      "Rebuild reference-metric helper functions analytically from CoordSystemName, xx[0..2], and stored parameters.",
      128,
      "reference_metric_reconstruction");
  raytracing_connections_write_f64_or_abort(fp, cart_origin[0], "Cart_origin[0]");
  raytracing_connections_write_f64_or_abort(fp, cart_origin[1], "Cart_origin[1]");
  raytracing_connections_write_f64_or_abort(fp, cart_origin[2], "Cart_origin[2]");
  raytracing_connections_write_f64_or_abort(fp, dxx[0], "dxx[0]");
  raytracing_connections_write_f64_or_abort(fp, dxx[1], "dxx[1]");
  raytracing_connections_write_f64_or_abort(fp, dxx[2], "dxx[2]");
  raytracing_connections_write_f64_or_abort(fp, invdxx[0], "invdxx[0]");
  raytracing_connections_write_f64_or_abort(fp, invdxx[1], "invdxx[1]");
  raytracing_connections_write_f64_or_abort(fp, invdxx[2], "invdxx[2]");
  raytracing_connections_write_f64_or_abort(fp, xxmin[0], "xxmin[0]");
  raytracing_connections_write_f64_or_abort(fp, xxmin[1], "xxmin[1]");
  raytracing_connections_write_f64_or_abort(fp, xxmin[2], "xxmin[2]");
  raytracing_connections_write_f64_or_abort(fp, xxmax[0], "xxmax[0]");
  raytracing_connections_write_f64_or_abort(fp, xxmax[1], "xxmax[1]");
  raytracing_connections_write_f64_or_abort(fp, xxmax[2], "xxmax[2]");
  raytracing_connections_write_f64_or_abort(fp, PI, "PI");
  raytracing_connections_write_f64_or_abort(fp, RMAX, "RMAX");
  raytracing_connections_write_f64_or_abort(fp, grid_hole_radius, "grid_hole_radius");
  raytracing_connections_write_f64_or_abort(fp, grid_physical_size, "grid_physical_size");
  raytracing_connections_write_or_abort(
      fp, &enable_rfm_precompute_i8, sizeof(int8_t), 1, "enable_rfm_precompute");
  raytracing_connections_write_or_abort(fp, &grid_rotates, sizeof(int8_t), 1, "grid_rotates");
  raytracing_connections_write_or_abort(fp, &is_host, sizeof(int8_t), 1, "is_host");
  raytracing_connections_write_or_abort(
      fp, reserved_flags, sizeof(int8_t), 5, "reserved_flags");
"""
    for _, gf_name in required_gfs:
        body += rf"""
  raytracing_connections_write_fixed_length_string_or_abort(
      fp, "{gf_name}", {gf_name_buffer_size}, "{gf_name}");
"""
    body += r"""

  raytracing_connections_assert_header_size_or_abort(fp, header_size);

  // Step 2: Write the physical simulation time in the canonical on-disk format.
  raytracing_connections_write_f64_or_abort(
      fp,
      (double)commondata->time,
      "simulation time");

  // Step 3: Write the logical coordinate arrays as little-endian IEEE-754 binary64.
  raytracing_connections_write_real_array_as_f64_or_abort(
      fp,
      griddata[0].xx[0],
      (uint64_t)Nxx_plus_2NGHOSTS0_u32,
      "xx[0]");
  raytracing_connections_write_real_array_as_f64_or_abort(
      fp,
      griddata[0].xx[1],
      (uint64_t)Nxx_plus_2NGHOSTS1_u32,
      "xx[1]");
  raytracing_connections_write_real_array_as_f64_or_abort(
      fp,
      griddata[0].xx[2],
      (uint64_t)Nxx_plus_2NGHOSTS2_u32,
      "xx[2]");

  // Step 4: Write exactly the BSSN gridfunctions needed for later spatial-state reconstruction,
  // serialized as little-endian IEEE-754 binary64.
"""
    for macro_name, gf_name in required_gfs:
        body += rf"""
  raytracing_connections_write_real_array_as_f64_or_abort(
      fp,
      &griddata[0].gridfuncs.y_n_gfs[ntot * (size_t){macro_name}],
      ntot,
      "{gf_name}");
"""
    body += r"""

  if (fflush(fp) != 0) {
    fprintf(stderr,
            "Error: output_raytracing_connections_4d_data could not flush %s. errno=%d\n",
            temporary_filename, errno);
    fclose(fp);
    exit(1);
  } // END IF: fflush failed

  if (fclose(fp) != 0) {
    fprintf(stderr,
            "Error: output_raytracing_connections_4d_data could not close %s cleanly. errno=%d\n",
            temporary_filename, errno);
    exit(1);
  } // END IF: fclose failed

  if (rename(temporary_filename, final_filename) != 0) {
    fprintf(stderr,
            "Error: output_raytracing_connections_4d_data could not rename %s to %s. errno=%d\n",
            temporary_filename, final_filename, errno);
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


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
