"""
Register C functions that export directly stored BHaH runtime state.

This module generates a C helper that writes host-side runtime snapshots for later
offline analysis. The export is intentionally limited to data already stored in
runtime memory, such as commondata metadata, concrete per-grid params metadata,
logical coordinate arrays, evolved gridfunctions, and optional precomputed
reference-metric arrays. The emitted files are organized into an atomic raw-snapshot
directory layout rooted directly under the configured export directory so a later
Python preparation step can normalize and index the dataset for
interpolation-heavy workloads. No finite differencing, interpolation, or tensor
reconstruction is performed here.

Author: Dalton Moone
        daltonmoone **at** gmail **dot** com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Optional, Tuple, Union, cast

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.infrastructures import BHaH


def _runtime_export_c_string_literal(value: str) -> str:
    """
    Escape a Python string for safe insertion into a generated C string literal.

    :param value: String to escape.
    :return: Escaped string content suitable for a C string literal.
    """
    return value.replace("\\", "\\\\").replace('"', '\\"')


def _runtime_export_dependency_axes_from_size_expr(size_expr: str) -> Tuple[int, ...]:
    """
    Infer reference-metric dependency axes from a generated size expression.

    :param size_expr: Generated size expression for one rfmstruct member.
    :return: Tuple of axis indices on which the member depends.
    """
    dependency_axes: List[int] = []
    for axis in range(3):
        if f"Nxx_plus_2NGHOSTS{axis}" in size_expr:
            dependency_axes.append(axis)
    return tuple(dependency_axes)


def _runtime_export_rfm_member_specs(
    CoordSystem: str, enable_rfm_precompute: bool
) -> List[Tuple[str, str, Tuple[int, ...]]]:
    """
    Build the list of rfmstruct member metadata for runtime export.

    :param CoordSystem: Coordinate system whose precomputed reference-metric members should
        be exported.
    :param enable_rfm_precompute: Whether runtime reference-metric precompute is enabled.
    :return: A list of `(member_name, size_expr, dependency_axes)` tuples for host-side
        binary export and metadata generation.
    """
    if not enable_rfm_precompute:
        return []
    rfm_precompute = BHaH.rfm_precompute.ReferenceMetricPrecompute(CoordSystem)
    return [
        (
            member_name,
            size_expr,
            _runtime_export_dependency_axes_from_size_expr(size_expr),
        )
        for member_name, size_expr in rfm_precompute.member_specs
    ]


def register_CFunction_stored_state_export(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    project_name: str = "nrpy_project",
    export_root: str = "stored_state/raw",
    run_id: str = "default_run",
    evol_gf_names: Optional[List[str]] = None,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct and register a C function that exports directly stored runtime state.

    The generated helper writes an atomic raw-snapshot directory tree containing:
    - snapshot-level JSON metadata,
    - concrete per-grid params metadata,
    - the host-side logical coordinate arrays xx[0..2],
    - all evolved gridfunctions in y_n_gfs, and
    - all host-side rfmstruct arrays when reference-metric precompute is enabled.

    It intentionally does not perform finite differencing, interpolation, or tensor
    reconstruction. It writes the values already stored in ``commondata``,
    ``griddata[grid].params``, ``griddata[grid].xx[0..2]``,
    ``griddata[grid].gridfuncs.y_n_gfs``, and optionally
    ``griddata[grid].rfmstruct``. The output is still a raw export of directly
    stored runtime state, but it includes enough machine-readable layout metadata
    for a later preparation/indexing script to validate and reorganize snapshots by
    time, grid, and field.

    :param CoordSystem: Coordinate system used to determine the rfmstruct member list.
        This is used only to determine the compile-time rfmstruct member list.
        Emitted snapshot and grid metadata use the runtime value
        ``params->CoordSystemName`` as the authoritative coordinate-system label.
    :param enable_rfm_precompute: Whether reference-metric precompute is enabled.
    :param project_name: Project name recorded in emitted metadata.
    :param export_root: Root directory under which snapshot_stepXXXXXXXX/ directories
        will be exported directly.
    :param run_id: Run identifier recorded in emitted metadata.
    :param evol_gf_names: Ordered list of EVOL gridfunction names to export. If not
        provided, the current BHaH EVOL registry is queried.
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    None.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    evol_gfs = (
        list(evol_gf_names)
        if evol_gf_names is not None
        else gri.BHaHGridFunction.gridfunction_lists()[0]
    )
    fp_type = cast(str, par.parval_from_str("fp_type"))
    rfm_member_specs = _runtime_export_rfm_member_specs(
        CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
    )
    project_name_c = _runtime_export_c_string_literal(project_name)
    export_root_c = _runtime_export_c_string_literal(export_root)
    run_id_c = _runtime_export_c_string_literal(run_id)

    includes = [
        "<stdarg.h>",
        "<string.h>",
        "<sys/stat.h>",
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = """
 * @file stored_state_export.c
 * @brief Export directly stored host-side runtime state needed for later offline analysis.
 *
 * The function "stored_state_export" writes one atomic raw snapshot per output step
 * under:
 *
 *   <export_root>/snapshot_stepXXXXXXXX/
 *
 * For each output step it stages a temporary snapshot directory, writes JSON
 * metadata plus one subdirectory per grid containing:
 *   1) concrete per-grid params metadata needed to interpret array layout,
 *   2) logical coordinate arrays xx[0], xx[1], and xx[2],
 *   3) all evolved gridfunctions stored in gridfuncs.y_n_gfs, and
 *   4) host-side rfmstruct arrays when reference-metric precompute is enabled.
 *
 * The temporary directory is renamed into place only after all files are written, so
 * the final snapshot directory itself is the completeness signal. No finite
 * differencing, interpolation, tensor reconstruction, or symbolic postprocessing is
 * performed here.
 *
 * @param[in] commondata  Global simulation metadata used in file naming and metadata output.
 * @param[in] griddata    Host-side per-grid runtime state to export.
 *
 * @return void.
"""
    cfunc_type = "void"
    name = "stored_state_export"
    params = (
        "const commondata_struct *restrict commondata, "
        "const griddata_struct *restrict griddata"
    )

    prefunc = r"""
static void stored_state_export_mkdir_or_abort(const char *restrict path) {
  if (mkdir(path, 0777) != 0 && errno != EEXIST) {
    fprintf(stderr, "Error: stored_state_export could not create directory %s. errno=%d\n", path, errno);
    exit(1);
  } // END IF: mkdir failed
} // END FUNCTION: stored_state_export_mkdir_or_abort

static void stored_state_export_mkdirs_or_abort(const char *restrict path) {
  char path_copy[1024];
  size_t length = strlen(path);
  if (length >= sizeof(path_copy)) {
    fprintf(stderr, "Error: stored_state_export path too long: %s\n", path);
    exit(1);
  } // END IF: path length exceeds buffer
  snprintf(path_copy, sizeof(path_copy), "%s", path);
  for (size_t i = 1; path_copy[i] != '\0'; i++) {
    if (path_copy[i] == '/') {
      path_copy[i] = '\0';
      if (strlen(path_copy) > 0) {
        stored_state_export_mkdir_or_abort(path_copy);
      } // END IF: non-empty path prefix
      path_copy[i] = '/';
    } // END IF: path separator found
  } // END LOOP: for i over path_copy characters
  stored_state_export_mkdir_or_abort(path_copy);
} // END FUNCTION: stored_state_export_mkdirs_or_abort

static void stored_state_export_path_must_not_exist_or_abort(const char *restrict path) {
  struct stat path_stat;
  errno = 0;
  if (stat(path, &path_stat) == 0) {
    fprintf(stderr, "Error: stored_state_export refusing to overwrite existing path %s\n", path);
    exit(1);
  } // END IF: path already exists
  if (errno != ENOENT) {
    fprintf(stderr, "Error: stored_state_export could not inspect path %s. errno=%d\n", path, errno);
    exit(1);
  } // END IF: stat failed for unexpected reason
} // END FUNCTION: stored_state_export_path_must_not_exist_or_abort

static void stored_state_export_rename_or_abort(
    const char *restrict old_path,
    const char *restrict new_path
) {
  if (rename(old_path, new_path) != 0) {
    fprintf(stderr,
            "Error: stored_state_export could not rename %s to %s. errno=%d\n",
            old_path, new_path, errno);
    exit(1);
  } // END IF: rename failed
} // END FUNCTION: stored_state_export_rename_or_abort

static const char *stored_state_export_host_byte_order(void) {
  const uint16_t one = 1;
  return (*((const uint8_t *)&one) == 1) ? "little" : "big";
} // END FUNCTION: stored_state_export_host_byte_order

static void stored_state_export_write_json_string(
    FILE *restrict fp,
    const char *restrict value
) {
  fputc('"', fp);
  for (const unsigned char *restrict ch = (const unsigned char *)value; *ch != '\0'; ch++) {
    if (*ch == '"' || *ch == '\\') {
      fputc('\\', fp);
      fputc((int)*ch, fp);
    } else if (*ch == '\n') {
      fputs("\\n", fp);
    } else if (*ch == '\r') {
      fputs("\\r", fp);
    } else if (*ch == '\t') {
      fputs("\\t", fp);
    } else if (*ch < 0x20U) {
      fprintf(fp, "\\u%04x", (unsigned int)*ch);
    } else {
      fputc((int)*ch, fp);
    } // END IF: choose JSON escaping mode
  } // END LOOP: for ch over JSON string bytes
  fputc('"', fp);
} // END FUNCTION: stored_state_export_write_json_string

static void stored_state_export_write_binary_or_abort(
    const char *restrict filename,
    const void *restrict data,
    const size_t element_size,
    const size_t count
) {
  FILE *fp = fopen(filename, "wb");
  if (fp == NULL) {
    fprintf(stderr, "Error: stored_state_export could not open %s for binary writing. errno=%d\n", filename, errno);
    exit(1);
  } // END IF: fopen failed
  const size_t wrote = fwrite(data, element_size, count, fp);
  if (wrote != count) {
    fprintf(stderr,
            "Error: stored_state_export short write for %s (wanted %zu elements, wrote %zu).\n",
            filename, count, wrote);
    fclose(fp);
    exit(1);
  } // END IF: fwrite short write
  fclose(fp);
} // END FUNCTION: stored_state_export_write_binary_or_abort

static void stored_state_export_format_path_or_abort(
    char *restrict buffer,
    const size_t buffer_size,
    const char *restrict format,
    ...
) {
  va_list args;
  va_start(args, format);
  const int written = vsnprintf(buffer, buffer_size, format, args);
  va_end(args);
  if (written < 0 || (size_t)written >= buffer_size) {
    fprintf(stderr, "Error: stored_state_export path too long while formatting '%s'\n", format);
    exit(1);
  } // END IF: formatted path overflowed buffer
} // END FUNCTION: stored_state_export_format_path_or_abort
"""

    evol_gf_name_entries = ", ".join(f'"{gf_name}"' for gf_name in evol_gfs)
    body = f"""
  static const char *restrict evol_gf_names[NUM_EVOL_GFS] = {{ {evol_gf_name_entries} }};
  const char *restrict byte_order = stored_state_export_host_byte_order();
  const char *restrict export_root = "{export_root_c}";
  const char *restrict project_name = "{project_name_c}";
  const char *restrict run_id = "{run_id_c}";
  char run_dir[1024];
  char snapshot_dir_tmp[1024];
  char snapshot_dir_final[1024];
  char filename[1024];
  char state_dir[1024];
  char coordinates_dir[1024];
  char rfm_dir[1024];
  char grid_dir[1024];

  stored_state_export_format_path_or_abort(run_dir, sizeof(run_dir), "%s", export_root);
  stored_state_export_mkdirs_or_abort(run_dir);
  stored_state_export_format_path_or_abort(
      snapshot_dir_tmp,
      sizeof(snapshot_dir_tmp),
      "%s/snapshot_step%08d.tmp",
      run_dir,
      commondata->nn
  );
  stored_state_export_format_path_or_abort(
      snapshot_dir_final,
      sizeof(snapshot_dir_final),
      "%s/snapshot_step%08d",
      run_dir,
      commondata->nn
  );
  stored_state_export_path_must_not_exist_or_abort(snapshot_dir_tmp);
  stored_state_export_path_must_not_exist_or_abort(snapshot_dir_final);
  stored_state_export_mkdir_or_abort(snapshot_dir_tmp);

  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {{
    const params_struct *restrict params = &griddata[grid].params;
    const REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    const REAL *restrict xx[3] = {{griddata[grid].xx[0], griddata[grid].xx[1], griddata[grid].xx[2]}};
    const char *restrict coord_system_name = params->CoordSystemName;
    const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
    const size_t ntot = (size_t)Nxx_plus_2NGHOSTS0 * (size_t)Nxx_plus_2NGHOSTS1 * (size_t)Nxx_plus_2NGHOSTS2;

    stored_state_export_format_path_or_abort(
        grid_dir, sizeof(grid_dir), "%s/grid_%02d", snapshot_dir_tmp, grid
    );
    stored_state_export_mkdir_or_abort(grid_dir);
    stored_state_export_format_path_or_abort(
        coordinates_dir, sizeof(coordinates_dir), "%s/coordinates", grid_dir
    );
    stored_state_export_mkdir_or_abort(coordinates_dir);
    stored_state_export_format_path_or_abort(
        state_dir, sizeof(state_dir), "%s/state", grid_dir
    );
    stored_state_export_mkdir_or_abort(state_dir);
"""
    if enable_rfm_precompute:
        body += """
    stored_state_export_format_path_or_abort(
        rfm_dir, sizeof(rfm_dir), "%s/rfm", grid_dir
    );
    stored_state_export_mkdir_or_abort(rfm_dir);
"""
    body += f"""

    stored_state_export_format_path_or_abort(
        filename, sizeof(filename), "%s/grid_manifest.json", grid_dir
    );
    FILE *metadata_fp = fopen(filename, "w");
    if (metadata_fp == NULL) {{
      fprintf(stderr, "Error: stored_state_export could not open %s for metadata writing. errno=%d\\n", filename, errno);
      exit(1);
    }} // END IF: metadata fopen failed

    fprintf(metadata_fp, "{{\\n");
    fprintf(metadata_fp, "  \\"schema_name\\": \\"stored_state_raw_grid\\",\\n");
    fprintf(metadata_fp, "  \\"schema_version\\": \\"1.0.0\\",\\n");
    fprintf(metadata_fp, "  \\"project_name\\": ");
    stored_state_export_write_json_string(metadata_fp, project_name);
    fprintf(metadata_fp, ",\\n");
    fprintf(metadata_fp, "  \\"run_id\\": ");
    stored_state_export_write_json_string(metadata_fp, run_id);
    fprintf(metadata_fp, ",\\n");
    fprintf(metadata_fp, "  \\"snapshot_id\\": \\"snapshot_step%08d\\",\\n", commondata->nn);
    fprintf(metadata_fp, "  \\"time\\": %.17e,\\n", commondata->time);
    fprintf(metadata_fp, "  \\"dt\\": %.17e,\\n", commondata->dt);
    fprintf(metadata_fp, "  \\"nn\\": %d,\\n", commondata->nn);
    fprintf(metadata_fp, "  \\"NUMGRIDS\\": %d,\\n", commondata->NUMGRIDS);
    fprintf(metadata_fp, "  \\"convergence_factor\\": %.17e,\\n", commondata->convergence_factor);
    fprintf(metadata_fp, "  \\"floating_point_type\\": \\"{fp_type}\\",\\n");
    fprintf(metadata_fp, "  \\"floating_point_bytes\\": %zu,\\n", sizeof(REAL));
    fprintf(metadata_fp, "  \\"byte_order\\": ");
    stored_state_export_write_json_string(metadata_fp, byte_order);
    fprintf(metadata_fp, ",\\n");
    fprintf(metadata_fp, "  \\"ghost_zones_included\\": true,\\n");
    fprintf(metadata_fp, "  \\"flattened_storage\\": true,\\n");
    fprintf(metadata_fp, "  \\"flat_index_order\\": [\\"i0\\", \\"i1\\", \\"i2\\"],\\n");
    fprintf(metadata_fp, "  \\"flat_index_storage_layout\\": \\"IDX3(i0,i1,i2) runtime ordering with i2 varying fastest\\",\\n");
    fprintf(metadata_fp, "  \\"grid\\": %d,\\n", grid);
    fprintf(metadata_fp, "  \\"grid_idx\\": %d,\\n", params->grid_idx);
    fprintf(metadata_fp, "  \\"gridname\\": ");
    stored_state_export_write_json_string(metadata_fp, params->gridname);
    fprintf(metadata_fp, ",\\n");
    fprintf(metadata_fp, "  \\"CoordSystemName\\": ");
    stored_state_export_write_json_string(metadata_fp, coord_system_name);
    fprintf(metadata_fp, ",\\n");
    fprintf(metadata_fp, "  \\"CoordSystem_hash\\": %d,\\n", params->CoordSystem_hash);
    fprintf(metadata_fp, "  \\"state_basis\\": \\"stored exactly as runtime y_n_gfs values in CoordSystemName coordinates\\",\\n");
    fprintf(metadata_fp, "  \\"coordinate_basis\\": \\"stored exactly as runtime xx[0..2] coordinate arrays for CoordSystemName\\",\\n");
    fprintf(metadata_fp, "  \\"NGHOSTS\\": %d,\\n", NGHOSTS);
    fprintf(metadata_fp, "  \\"Nxx\\": [%d, %d, %d],\\n",
            params->Nxx0,
            params->Nxx1,
            params->Nxx2);
    fprintf(metadata_fp, "  \\"Nxx_plus_2NGHOSTS\\": [%d, %d, %d],\\n",
            params->Nxx_plus_2NGHOSTS0,
            params->Nxx_plus_2NGHOSTS1,
            params->Nxx_plus_2NGHOSTS2);
    fprintf(metadata_fp, "  \\"xxmin\\": [%.17e, %.17e, %.17e],\\n",
            params->xxmin0,
            params->xxmin1,
            params->xxmin2);
    fprintf(metadata_fp, "  \\"xxmax\\": [%.17e, %.17e, %.17e],\\n",
            params->xxmax0,
            params->xxmax1,
            params->xxmax2);
    fprintf(metadata_fp, "  \\"dxx\\": [%.17e, %.17e, %.17e],\\n",
            params->dxx0,
            params->dxx1,
            params->dxx2);
    fprintf(metadata_fp, "  \\"invdxx\\": [%.17e, %.17e, %.17e],\\n",
            params->invdxx0,
            params->invdxx1,
            params->invdxx2);
    fprintf(metadata_fp, "  \\"array_shape_with_ghosts\\": [%d, %d, %d],\\n",
            params->Nxx_plus_2NGHOSTS0,
            params->Nxx_plus_2NGHOSTS1,
            params->Nxx_plus_2NGHOSTS2);
    fprintf(metadata_fp, "  \\"num_total_grid_points_with_ghosts\\": %zu,\\n", ntot);
    fprintf(metadata_fp, "  \\"coordinates\\": [\\n");
    fprintf(metadata_fp,
            "    {{\\"name\\": \\"xx0\\", \\"filename\\": \\"coordinates/xx0.bin\\", \\"dtype\\": \\"{fp_type}\\", \\"count\\": %d, \\"shape\\": [%d], \\"includes_ghost_zones\\": true, \\"basis\\": \\"stored exactly as runtime xx0 array in CoordSystemName coordinates\\"}},\\n",
            params->Nxx_plus_2NGHOSTS0,
            params->Nxx_plus_2NGHOSTS0);
    fprintf(metadata_fp,
            "    {{\\"name\\": \\"xx1\\", \\"filename\\": \\"coordinates/xx1.bin\\", \\"dtype\\": \\"{fp_type}\\", \\"count\\": %d, \\"shape\\": [%d], \\"includes_ghost_zones\\": true, \\"basis\\": \\"stored exactly as runtime xx1 array in CoordSystemName coordinates\\"}},\\n",
            params->Nxx_plus_2NGHOSTS1,
            params->Nxx_plus_2NGHOSTS1);
    fprintf(metadata_fp,
            "    {{\\"name\\": \\"xx2\\", \\"filename\\": \\"coordinates/xx2.bin\\", \\"dtype\\": \\"{fp_type}\\", \\"count\\": %d, \\"shape\\": [%d], \\"includes_ghost_zones\\": true, \\"basis\\": \\"stored exactly as runtime xx2 array in CoordSystemName coordinates\\"}}\\n",
            params->Nxx_plus_2NGHOSTS2,
            params->Nxx_plus_2NGHOSTS2);
    fprintf(metadata_fp, "  ],\\n");
    fprintf(metadata_fp, "  \\"num_evol_gfs\\": %d,\\n", NUM_EVOL_GFS);
    fprintf(metadata_fp, "  \\"state_fields\\": [\\n");
    for (int gf = 0; gf < NUM_EVOL_GFS; gf++) {{
      fprintf(metadata_fp,
              "    {{\\"field_name\\": \\"%s\\", \\"filename\\": \\"state/%s.bin\\", \\"dtype\\": \\"{fp_type}\\", \\"count\\": %zu, \\"shape\\": [%d, %d, %d], \\"semantic_class\\": \\"evolved_gridfunction\\"}}%s\\n",
              evol_gf_names[gf],
              evol_gf_names[gf],
              ntot,
              params->Nxx_plus_2NGHOSTS0,
              params->Nxx_plus_2NGHOSTS1,
              params->Nxx_plus_2NGHOSTS2,
              (gf + 1 < NUM_EVOL_GFS) ? "," : "");
    }} // END LOOP: for gf over evolved gridfunctions
    fprintf(metadata_fp, "  ],\\n");
"""
    if enable_rfm_precompute:
        body += """
    fprintf(metadata_fp, "  \\"rfm_enabled\\": true,\\n");
    fprintf(metadata_fp, "  \\"rfm_manifest_path\\": \\"rfm/rfm_manifest.json\\"\\n");
"""
    else:
        body += """
    fprintf(metadata_fp, "  \\"rfm_enabled\\": false\\n");
"""
    body += """
    fprintf(metadata_fp, "}\\n");
    fclose(metadata_fp);

    stored_state_export_format_path_or_abort(
        filename, sizeof(filename), "%s/xx0.bin", coordinates_dir
    );
    stored_state_export_write_binary_or_abort(filename, xx[0], sizeof(REAL), (size_t)params->Nxx_plus_2NGHOSTS0);
    stored_state_export_format_path_or_abort(
        filename, sizeof(filename), "%s/xx1.bin", coordinates_dir
    );
    stored_state_export_write_binary_or_abort(filename, xx[1], sizeof(REAL), (size_t)params->Nxx_plus_2NGHOSTS1);
    stored_state_export_format_path_or_abort(
        filename, sizeof(filename), "%s/xx2.bin", coordinates_dir
    );
    stored_state_export_write_binary_or_abort(filename, xx[2], sizeof(REAL), (size_t)params->Nxx_plus_2NGHOSTS2);

    for (int gf = 0; gf < NUM_EVOL_GFS; gf++) {
      stored_state_export_format_path_or_abort(
          filename,
          sizeof(filename),
          "%s/%s.bin",
          state_dir,
          evol_gf_names[gf]
      );
      stored_state_export_write_binary_or_abort(filename, &y_n_gfs[ntot * (size_t)gf], sizeof(REAL), ntot);
    } // END LOOP: for gf over evolved gridfunctions
"""

    if enable_rfm_precompute:
        body += """
    const rfm_struct *restrict rfmstruct = griddata[grid].rfmstruct;
    stored_state_export_format_path_or_abort(
        filename, sizeof(filename), "%s/rfm_manifest.json", rfm_dir
    );
    FILE *rfm_metadata_fp = fopen(filename, "w");
    if (rfm_metadata_fp == NULL) {
      fprintf(stderr, "Error: stored_state_export could not open %s for rfm metadata writing. errno=%d\\n", filename, errno);
      exit(1);
    } // END IF: rfm metadata fopen failed
    fprintf(rfm_metadata_fp, "{\\n");
    fprintf(rfm_metadata_fp, "  \\"schema_name\\": \\"stored_state_raw_rfm\\",\\n");
    fprintf(rfm_metadata_fp, "  \\"schema_version\\": \\"1.0.0\\",\\n");
    fprintf(rfm_metadata_fp, "  \\"project_name\\": ");
    stored_state_export_write_json_string(rfm_metadata_fp, project_name);
    fprintf(rfm_metadata_fp, ",\\n");
    fprintf(rfm_metadata_fp, "  \\"run_id\\": ");
    stored_state_export_write_json_string(rfm_metadata_fp, run_id);
    fprintf(rfm_metadata_fp, ",\\n");
    fprintf(rfm_metadata_fp, "  \\"snapshot_id\\": \\"snapshot_step%08d\\",\\n", commondata->nn);
    fprintf(rfm_metadata_fp, "  \\"grid\\": %d,\\n", grid);
    fprintf(rfm_metadata_fp, "  \\"grid_idx\\": %d,\\n", params->grid_idx);
    fprintf(rfm_metadata_fp, "  \\"CoordSystemName\\": ");
    stored_state_export_write_json_string(rfm_metadata_fp, params->CoordSystemName);
    fprintf(rfm_metadata_fp, ",\\n");
    fprintf(rfm_metadata_fp, "  \\"rfm_enabled\\": true,\\n");
    fprintf(rfm_metadata_fp, "  \\"members\\": [\\n");
"""
        for idx, (member_name, size_expr, dependency_axes) in enumerate(
            rfm_member_specs
        ):
            dependency_axes_json = ", ".join(str(axis) for axis in dependency_axes)
            comma = "," if idx + 1 < len(rfm_member_specs) else ""
            if dependency_axes:
                shape_expr = ", ".join("%d" for _ in dependency_axes)
                shape_values = ", ".join(
                    f"params->Nxx_plus_2NGHOSTS{axis}" for axis in dependency_axes
                )
                body += f"""    fprintf(rfm_metadata_fp,
            "    {{\\"member_name\\": \\"{member_name}\\", \\"filename\\": \\"{member_name}.bin\\", \\"dtype\\": \\"{fp_type}\\", \\"count\\": %zu, \\"shape\\": [{shape_expr}], \\"rank\\": {len(dependency_axes)}, \\"dependency_axes\\": [{dependency_axes_json}], \\"semantic_class\\": \\"rfm_precompute\\"}}{comma}\\n",
            {size_expr},
            {shape_values});
"""
            else:
                body += f"""    fprintf(rfm_metadata_fp,
            "    {{\\"member_name\\": \\"{member_name}\\", \\"filename\\": \\"{member_name}.bin\\", \\"dtype\\": \\"{fp_type}\\", \\"count\\": %zu, \\"shape\\": [], \\"rank\\": 0, \\"dependency_axes\\": [], \\"semantic_class\\": \\"rfm_precompute\\"}}{comma}\\n",
            {size_expr});
"""
            body += f"""    stored_state_export_format_path_or_abort(
             filename,
             sizeof(filename),
             "%s/{member_name}.bin",
             rfm_dir);
    stored_state_export_write_binary_or_abort(filename, rfmstruct->{member_name}, sizeof(REAL), {size_expr});
"""
        body += """
    fprintf(rfm_metadata_fp, "  ]\\n");
    fprintf(rfm_metadata_fp, "}\\n");
    fclose(rfm_metadata_fp);
"""

    body += f"""
  }} // END LOOP: for grid over active grids

  stored_state_export_format_path_or_abort(
      filename, sizeof(filename), "%s/snapshot_manifest.json", snapshot_dir_tmp
  );
  FILE *snapshot_fp = fopen(filename, "w");
  if (snapshot_fp == NULL) {{
    fprintf(stderr, "Error: stored_state_export could not open %s for snapshot metadata writing. errno=%d\\n", filename, errno);
    exit(1);
  }} // END IF: snapshot metadata fopen failed
  fprintf(snapshot_fp, "{{\\n");
  fprintf(snapshot_fp, "  \\"schema_name\\": \\"stored_state_raw_snapshot\\",\\n");
  fprintf(snapshot_fp, "  \\"schema_version\\": \\"1.0.0\\",\\n");
  fprintf(snapshot_fp, "  \\"project_name\\": ");
  stored_state_export_write_json_string(snapshot_fp, project_name);
  fprintf(snapshot_fp, ",\\n");
  fprintf(snapshot_fp, "  \\"run_id\\": ");
  stored_state_export_write_json_string(snapshot_fp, run_id);
  fprintf(snapshot_fp, ",\\n");
  fprintf(snapshot_fp, "  \\"snapshot_id\\": \\"snapshot_step%08d\\",\\n", commondata->nn);
  fprintf(snapshot_fp, "  \\"nn\\": %d,\\n", commondata->nn);
  fprintf(snapshot_fp, "  \\"time\\": %.17e,\\n", commondata->time);
  fprintf(snapshot_fp, "  \\"dt\\": %.17e,\\n", commondata->dt);
  fprintf(snapshot_fp, "  \\"NUMGRIDS\\": %d,\\n", commondata->NUMGRIDS);
  fprintf(snapshot_fp, "  \\"convergence_factor\\": %.17e,\\n", commondata->convergence_factor);
  fprintf(snapshot_fp, "  \\"floating_point_type\\": \\"{fp_type}\\",\\n");
  fprintf(snapshot_fp, "  \\"floating_point_bytes\\": %zu,\\n", sizeof(REAL));
  fprintf(snapshot_fp, "  \\"byte_order\\": ");
  stored_state_export_write_json_string(snapshot_fp, byte_order);
  fprintf(snapshot_fp, ",\\n");
  fprintf(snapshot_fp, "  \\"ghost_zones_included\\": true,\\n");
  fprintf(snapshot_fp, "  \\"flattened_storage\\": true,\\n");
  fprintf(snapshot_fp, "  \\"flat_index_order\\": [\\"i0\\", \\"i1\\", \\"i2\\"],\\n");
  fprintf(snapshot_fp, "  \\"flat_index_storage_layout\\": \\"IDX3(i0,i1,i2) runtime ordering with i2 varying fastest\\",\\n");
  fprintf(snapshot_fp, "  \\"grid_coord_systems\\": [");
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {{
    stored_state_export_write_json_string(snapshot_fp, griddata[grid].params.CoordSystemName);
    fprintf(snapshot_fp, "%s", (grid + 1 < commondata->NUMGRIDS) ? ", " : "");
  }} // END LOOP: for grid over runtime coordinate-system labels
  fprintf(snapshot_fp, "],\\n");
"""
    if enable_rfm_precompute:
        body += """
  fprintf(snapshot_fp, "  \\"rfm_precompute_enabled\\": true,\\n");
"""
    else:
        body += """
  fprintf(snapshot_fp, "  \\"rfm_precompute_enabled\\": false,\\n");
"""
    body += """
  fprintf(snapshot_fp, "  \\"num_evol_gfs\\": %d,\\n", NUM_EVOL_GFS);
  fprintf(snapshot_fp, "  \\"evol_field_names\\": [");
  for (int gf = 0; gf < NUM_EVOL_GFS; gf++) {
    fprintf(snapshot_fp, "\\"%s\\"%s", evol_gf_names[gf], (gf + 1 < NUM_EVOL_GFS) ? ", " : "");
  } // END LOOP: for gf over evolved gridfunction names
  fprintf(snapshot_fp, "],\\n");
  fprintf(snapshot_fp, "  \\"grid_directories\\": [");
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    fprintf(snapshot_fp, "\\"grid_%02d\\"%s", grid, (grid + 1 < commondata->NUMGRIDS) ? ", " : "");
  } // END LOOP: for grid over grid directories
  fprintf(snapshot_fp, "]\\n");
  fprintf(snapshot_fp, "}\\n");
  fclose(snapshot_fp);

  stored_state_export_rename_or_abort(snapshot_dir_tmp, snapshot_dir_final);
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
