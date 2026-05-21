"""
Register C functions that export raytracing-connection 4D reconstruction data.

This module generates a C helper that writes exactly one binary file per
diagnostics output time. Each file contains:

1. A fixed-size self-describing binary header.
2. The physical simulation time.
3. The full raw ``params_struct`` stored in ``griddata[0].params``.
4. The logical coordinate arrays ``xx[0..2]``.
5. Only the evolved BSSN gridfunctions needed to reconstruct the spatial BSSN
   state consumed by
   ``nrpy.equations.general_relativity.geodesics.symbolic_christoffel_recipe_from_bssn_grid_basis()``.

The fixed-size header is designed so that a later post-processing or
interpolation reader does not need to know the in-memory C layout of
``params_struct`` just to interpret the payload layout. In addition to
recording the raw ``params_struct`` byte count, the header stores explicit
byte offsets, per-array lengths, per-array byte sizes, floating-point format,
the evolved conformal-factor convention, the reference-metric precompute mode,
coordinate-system labels, and interpolation-critical grid metadata such as
interior and ghosted extents, coordinate minima/maxima, grid spacings, and
coordinate/gridfunction naming.

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
    ``raytracing_connections_4D_data_t########.bin`` file.

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


def output_raytracing_connections_4D_data(
    enable_rfm_precompute: bool,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct and register the raytracing-connections 4D binary exporter.

    The generated C helper writes one binary file per diagnostics output with
    filename pattern ``./raytracing_connections_4D_data_t########.bin``. The
    file is first written to a temporary sibling path and renamed into place
    only after all bytes are successfully flushed. If the final file already
    exists, the exporter returns immediately so restart/checkpoint workflows do
    not overwrite or append duplicate data.

    File contents, in order:
    1. A fixed-size binary header describing the payload layout, byte offsets,
       array lengths, array sizes, floating-point precision, byte order, the
       evolved conformal-factor convention, the reference-metric precompute
       mode, and the interpolation-critical grid metadata duplicated from
       ``griddata[0].params``.
    2. The physical simulation time ``commondata->time``.
    3. The raw ``params_struct`` bytes from ``griddata[0].params``.
    4. The logical coordinate arrays ``xx[0]``, ``xx[1]``, and ``xx[2]``.
    5. The exact required BSSN gridfunctions, one full gridfunction array at a
       time.

    One output file contains the spatial BSSN state and metadata needed for
    later interpolation/finite-difference reconstruction. Time derivatives such
    as ``alpha_d0`` and ``vetU_d0`` are not directly stored in a single file.

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
    fp_type = cast(str, par.parval_from_str("fp_type"))
    cf_convention = cast(str, par.parval_from_str("EvolvedConformalFactor_cf"))

    includes = [
        "<errno.h>",
        "<stdint.h>",
        "<stdio.h>",
        "<stdlib.h>",
        "<string.h>",
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = """
 * @file output_raytracing_connections_4D_data.c
 * @brief Export raytracing-connection 4D reconstruction data from the stored BSSN state.
 *
 * This function writes at most one binary file per diagnostics output to:
 *
 *   ./raytracing_connections_4D_data_t########.bin
 *
 * where the eight-digit index is the diagnostics output index, not the raw
 * timestep counter. Each file is written to a temporary sibling path and then
 * atomically renamed into place after all bytes are flushed successfully. If
 * the final file already exists, the function returns without modifying it so
 * restarted runs do not overwrite or append duplicate output.
 *
 * The exported payload is intentionally limited to the quantities needed to
 * reconstruct the spatial BSSN state used later for 4D Christoffel
 * reconstruction:
 *   - physical simulation time,
 *   - full raw params_struct bytes from griddata[0].params,
 *   - duplicated interpolation-critical grid metadata in the fixed header,
 *   - logical coordinate arrays xx[0], xx[1], xx[2],
 *   - aDD00, aDD01, aDD02, aDD11, aDD12, aDD22,
 *   - alpha,
 *   - cf,
 *   - hDD00, hDD01, hDD02, hDD11, hDD12, hDD22,
 *   - trK,
 *   - vetU0, vetU1, vetU2.
 *
 * The fixed header records byte offsets and byte counts for every payload
 * section, the EvolvedConformalFactor_cf convention used to interpret the
 * stored cf field, and whether reference-metric precompute was enabled when
 * the project was generated. A later reader can therefore parse the coordinate
 * arrays and each stored gridfunction deterministically without decoding the
 * opaque raw params_struct bytes first. The raw params_struct payload is still
 * written as the authoritative complete snapshot of griddata[0].params,
 * matching the requirement that all stored grid parameters be present.
 *
 * Fields not used by the current BSSN-to-Christoffel construction, such as
 * betU and lambdaU, are deliberately omitted. Reference-metric helper arrays
 * such as f0_of_xx0 are also omitted because they are reconstructed
 * analytically from CoordSystemName, xx[0..2], and the stored runtime
 * parameters instead of being treated as independent simulation data.
 *
 * @param[in] commondata  Global simulation metadata used in the file header and output naming.
 * @param[in] griddata    Host-side per-grid runtime state to export.
 * @param[in] time_output_index  Diagnostics output index used in the emitted filename.
 *
 * @return void
"""
    cfunc_type = "void"
    name = "output_raytracing_connections_4D_data"
    params = (
        "const commondata_struct *restrict commondata, "
        "const griddata_struct *restrict griddata, "
        "const int time_output_index"
    )

    prefunc = rf"""
typedef struct __raytracing_connections_4D_file_header__ {{
  char magic[16];
  uint32_t format_version;
  uint32_t header_size;
  uint32_t output_index;
  uint32_t num_grids;
  uint32_t sizeof_REAL;
  uint32_t sizeof_params_struct;
  uint32_t coord_array_count;
  uint32_t num_export_gridfunctions;
  uint32_t gridfunction_name_bytes;
  uint32_t coord_name_bytes;
  uint32_t grid_idx;
  uint32_t host_is_little_endian;
  uint32_t time_variable_is_REAL;
  uint32_t reserved_u32;
  uint32_t Nxx0;
  uint32_t Nxx1;
  uint32_t Nxx2;
  uint32_t Nxx_plus_2NGHOSTS0;
  uint32_t Nxx_plus_2NGHOSTS1;
  uint32_t Nxx_plus_2NGHOSTS2;
  uint64_t ntot;
  uint64_t simulation_time_offset;
  uint64_t params_struct_offset;
  uint64_t coord_array_offsets[3];
  uint64_t gridfunction_payload_offset;
  uint64_t total_file_bytes;
  uint64_t coord_array_lengths[3];
  uint64_t coord_array_bytes[3];
  uint64_t gridfunction_bytes_per_array;
  uint64_t gridfunction_offsets[{num_export_gfs}];
  char real_type[16];
  char conformal_factor_cf[16];
  char params_struct_type[{type_name_buffer_size}];
  char coord_names[3][{coord_name_buffer_size}];
  char coord_system_name[100];
  char grid_name[100];
  char header_description[128];
  char params_payload_description[128];
  char coord_payload_description[128];
  char gridfunction_payload_description[128];
  char reference_metric_reconstruction[128];
  REAL Cart_origin[3];
  REAL dxx[3];
  REAL invdxx[3];
  REAL xxmin[3];
  REAL xxmax[3];
  REAL PI;
  REAL RMAX;
  REAL grid_hole_radius;
  REAL grid_physical_size;
  int8_t enable_rfm_precompute;
  int8_t grid_rotates;
  int8_t is_host;
  int8_t reserved_flags[5];
  char gridfunction_names[{num_export_gfs}][{gf_name_buffer_size}];
}} raytracing_connections_4D_file_header;

static int raytracing_connections_host_is_little_endian(void) {{
  const uint16_t one = 1;
  return (*((const uint8_t *)&one) == 1) ? 1 : 0;
}} // END FUNCTION: raytracing_connections_host_is_little_endian

static void raytracing_connections_abort_with_message(const char *restrict message) {{
  fprintf(stderr, "%s\n", message);
  exit(1);
}} // END FUNCTION: raytracing_connections_abort_with_message

static void raytracing_connections_write_or_abort(
    FILE *restrict fp,
    const void *restrict data,
    const size_t element_size,
    const size_t count,
    const char *restrict label
) {{
  if (fwrite(data, element_size, count, fp) != count) {{
    fprintf(stderr, "Error: output_raytracing_connections_4D_data failed while writing %s.\\n", label);
    fclose(fp);
    exit(1);
  }} // END IF: fwrite failed
}} // END FUNCTION: raytracing_connections_write_or_abort

static int raytracing_connections_file_exists(const char *restrict path) {{
  FILE *restrict fp = fopen(path, "rb");
  if (fp != NULL) {{
    fclose(fp);
    return 1;
  }} // END IF: existing file found
  return 0;
}} // END FUNCTION: raytracing_connections_file_exists
"""

    body = rf"""
  if (commondata->NUMGRIDS != 1) {{
    raytracing_connections_abort_with_message(
        "Error: output_raytracing_connections_4D_data currently requires commondata->NUMGRIDS == 1.");
  }} // END IF: unsupported number of grids

  if (time_output_index < 0) {{
    raytracing_connections_abort_with_message(
        "Error: output_raytracing_connections_4D_data received a negative output index.");
  }} // END IF: invalid output index

  const params_struct *restrict params = &griddata[0].params;
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  const uint64_t ntot =
      (size_t)Nxx_plus_2NGHOSTS0 * (size_t)Nxx_plus_2NGHOSTS1 * (size_t)Nxx_plus_2NGHOSTS2;
  const uint64_t coord_array_bytes0 = (uint64_t)Nxx_plus_2NGHOSTS0 * (uint64_t)sizeof(REAL);
  const uint64_t coord_array_bytes1 = (uint64_t)Nxx_plus_2NGHOSTS1 * (uint64_t)sizeof(REAL);
  const uint64_t coord_array_bytes2 = (uint64_t)Nxx_plus_2NGHOSTS2 * (uint64_t)sizeof(REAL);
  const uint64_t gridfunction_bytes_per_array = ntot * (uint64_t)sizeof(REAL);

  char final_filename[256];
  char temporary_filename[256];
  snprintf(
      final_filename,
      sizeof(final_filename),
      "./raytracing_connections_4D_data_t%08d.bin",
      time_output_index);
  snprintf(
      temporary_filename,
      sizeof(temporary_filename),
      "./raytracing_connections_4D_data_t%08d.bin.tmp",
      time_output_index);

  if (raytracing_connections_file_exists(final_filename)) {{
    return;
  }} // END IF: output file already exists

  remove(temporary_filename);

  FILE *restrict fp = fopen(temporary_filename, "wb");
  if (fp == NULL) {{
    fprintf(stderr,
            "Error: output_raytracing_connections_4D_data could not open %s for writing. errno=%d\\n",
            temporary_filename, errno);
    exit(1);
  }} // END IF: fopen failed

  raytracing_connections_4D_file_header header;
  memset(&header, 0, sizeof(raytracing_connections_4D_file_header));
  snprintf(header.magic, sizeof(header.magic), "NRPYRTCONN4D");
  header.format_version = 3U;
  header.header_size = (uint32_t)sizeof(raytracing_connections_4D_file_header);
  header.output_index = (uint32_t)time_output_index;
  header.num_grids = (uint32_t)commondata->NUMGRIDS;
  header.sizeof_REAL = (uint32_t)sizeof(REAL);
  header.sizeof_params_struct = (uint32_t)sizeof(params_struct);
  header.coord_array_count = 3;
  header.num_export_gridfunctions = {num_export_gfs};
  header.gridfunction_name_bytes = {gf_name_buffer_size};
  header.coord_name_bytes = {coord_name_buffer_size};
  header.grid_idx = (uint32_t)params->grid_idx;
  header.host_is_little_endian = (uint32_t)raytracing_connections_host_is_little_endian();
  header.time_variable_is_REAL = 1U;
  header.Nxx0 = (uint32_t)params->Nxx0;
  header.Nxx1 = (uint32_t)params->Nxx1;
  header.Nxx2 = (uint32_t)params->Nxx2;
  header.Nxx_plus_2NGHOSTS0 = (uint32_t)Nxx_plus_2NGHOSTS0;
  header.Nxx_plus_2NGHOSTS1 = (uint32_t)Nxx_plus_2NGHOSTS1;
  header.Nxx_plus_2NGHOSTS2 = (uint32_t)Nxx_plus_2NGHOSTS2;
  header.ntot = ntot;
  header.simulation_time_offset = (uint64_t)sizeof(raytracing_connections_4D_file_header);
  header.params_struct_offset = header.simulation_time_offset + (uint64_t)sizeof(REAL);
  header.coord_array_offsets[0] = header.params_struct_offset + (uint64_t)sizeof(params_struct);
  header.coord_array_offsets[1] = header.coord_array_offsets[0] + coord_array_bytes0;
  header.coord_array_offsets[2] = header.coord_array_offsets[1] + coord_array_bytes1;
  header.gridfunction_payload_offset = header.coord_array_offsets[2] + coord_array_bytes2;
  header.coord_array_lengths[0] = (uint64_t)Nxx_plus_2NGHOSTS0;
  header.coord_array_lengths[1] = (uint64_t)Nxx_plus_2NGHOSTS1;
  header.coord_array_lengths[2] = (uint64_t)Nxx_plus_2NGHOSTS2;
  header.coord_array_bytes[0] = coord_array_bytes0;
  header.coord_array_bytes[1] = coord_array_bytes1;
  header.coord_array_bytes[2] = coord_array_bytes2;
  header.gridfunction_bytes_per_array = gridfunction_bytes_per_array;
  header.total_file_bytes =
      header.gridfunction_payload_offset + (uint64_t){num_export_gfs} * gridfunction_bytes_per_array;
  snprintf(header.real_type, sizeof(header.real_type), "%s", "{fp_type}");
  snprintf(header.conformal_factor_cf, sizeof(header.conformal_factor_cf), "%s", "{cf_convention}");
  snprintf(header.params_struct_type, sizeof(header.params_struct_type), "%s", "params_struct");
  snprintf(header.coord_names[0], sizeof(header.coord_names[0]), "%s", "xx0");
  snprintf(header.coord_names[1], sizeof(header.coord_names[1]), "%s", "xx1");
  snprintf(header.coord_names[2], sizeof(header.coord_names[2]), "%s", "xx2");
  snprintf(header.coord_system_name, sizeof(header.coord_system_name), "%s", params->CoordSystemName);
  snprintf(header.grid_name, sizeof(header.grid_name), "%s", params->gridname);
  snprintf(
      header.header_description,
      sizeof(header.header_description),
      "%s",
      "Fixed header with payload offsets, array lengths, and interpolation metadata.");
  snprintf(
      header.params_payload_description,
      sizeof(header.params_payload_description),
      "%s",
      "Raw params_struct bytes follow immediately after simulation time.");
  snprintf(
      header.coord_payload_description,
      sizeof(header.coord_payload_description),
      "%s",
      "Three logical coordinate arrays xx0, xx1, xx2 follow in this listed order.");
  snprintf(
      header.gridfunction_payload_description,
      sizeof(header.gridfunction_payload_description),
      "%s",
      "Each listed gridfunction follows as one full contiguous array.");
  snprintf(
      header.reference_metric_reconstruction,
      sizeof(header.reference_metric_reconstruction),
      "%s",
      "Rebuild reference-metric helper functions analytically from CoordSystemName, xx[0..2], and stored parameters.");
  header.Cart_origin[0] = params->Cart_originx;
  header.Cart_origin[1] = params->Cart_originy;
  header.Cart_origin[2] = params->Cart_originz;
  header.dxx[0] = params->dxx0;
  header.dxx[1] = params->dxx1;
  header.dxx[2] = params->dxx2;
  header.invdxx[0] = params->invdxx0;
  header.invdxx[1] = params->invdxx1;
  header.invdxx[2] = params->invdxx2;
  header.xxmin[0] = params->xxmin0;
  header.xxmin[1] = params->xxmin1;
  header.xxmin[2] = params->xxmin2;
  header.xxmax[0] = params->xxmax0;
  header.xxmax[1] = params->xxmax1;
  header.xxmax[2] = params->xxmax2;
  header.PI = params->PI;
  header.RMAX = params->RMAX;
  header.grid_hole_radius = params->grid_hole_radius;
  header.grid_physical_size = params->grid_physical_size;
  header.enable_rfm_precompute = {1 if enable_rfm_precompute else 0};
  header.grid_rotates = params->grid_rotates ? 1 : 0;
  header.is_host = params->is_host ? 1 : 0;
"""
    for gf_index, (_, gf_name) in enumerate(required_gfs):
        body += rf"""
  header.gridfunction_offsets[{gf_index}] =
      header.gridfunction_payload_offset + (uint64_t){gf_index} * gridfunction_bytes_per_array;
  snprintf(
      header.gridfunction_names[{gf_index}],
      sizeof(header.gridfunction_names[{gf_index}]),
      "%s",
      "{gf_name}");
"""
    body += """

  // Step 1: Write the fixed-size file header.
  raytracing_connections_write_or_abort(
      fp,
      &header,
      sizeof(raytracing_connections_4D_file_header),
      1,
      "output_raytracing_connections_4D_data header");

  // Step 2: Write the physical simulation time and full raw params_struct bytes.
  raytracing_connections_write_or_abort(fp, &commondata->time, sizeof(REAL), 1, "simulation time");
  raytracing_connections_write_or_abort(fp, params, sizeof(params_struct), 1, "params_struct");

  // Step 3: Write the logical coordinate arrays needed for later interpolation.
  raytracing_connections_write_or_abort(fp, griddata[0].xx[0], sizeof(REAL), (size_t)Nxx_plus_2NGHOSTS0, "xx[0]");
  raytracing_connections_write_or_abort(fp, griddata[0].xx[1], sizeof(REAL), (size_t)Nxx_plus_2NGHOSTS1, "xx[1]");
  raytracing_connections_write_or_abort(fp, griddata[0].xx[2], sizeof(REAL), (size_t)Nxx_plus_2NGHOSTS2, "xx[2]");

  // Step 4: Write exactly the BSSN gridfunctions needed for later spatial-state reconstruction.
"""
    for macro_name, gf_name in required_gfs:
        body += rf"""
  raytracing_connections_write_or_abort(
      fp,
      &griddata[0].gridfuncs.y_n_gfs[ntot * (size_t){macro_name}],
      sizeof(REAL),
      ntot,
      "{gf_name}");
"""
    body += """

  if (fflush(fp) != 0) {
    fprintf(stderr,
            "Error: output_raytracing_connections_4D_data could not flush %s. errno=%d\\n",
            temporary_filename, errno);
    fclose(fp);
    exit(1);
  } // END IF: fflush failed

  if (fclose(fp) != 0) {
    fprintf(stderr,
            "Error: output_raytracing_connections_4D_data could not close %s cleanly. errno=%d\\n",
            temporary_filename, errno);
    exit(1);
  } // END IF: fclose failed

  if (rename(temporary_filename, final_filename) != 0) {
    fprintf(stderr,
            "Error: output_raytracing_connections_4D_data could not rename %s to %s. errno=%d\\n",
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
