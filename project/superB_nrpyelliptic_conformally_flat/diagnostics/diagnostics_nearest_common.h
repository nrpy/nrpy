/**
 * @file diagnostics_nearest_common.h
 * @brief Helpers for standardized diagnostic file I/O (headers, rows, filenames) for "nearest" samplers.
 *
 * This header defines small inline utilities used by diagnostics_nearest_* sources
 * to build consistent filenames, emit column headers, and write numeric rows.
 */

#ifndef DIAGNOSTICS_NEAREST_COMMON_H
#define DIAGNOSTICS_NEAREST_COMMON_H

#ifdef __cplusplus
extern "C" {
#endif

//OLD
//~ /**
 //~ * @brief Write the simulation time as a header comment.
 //~ *
 //~ * Prints a single line:
 //~ *   # [time] = <time>
 //~ *
 //~ * @param[in] file_ptr  File pointer to write to.
 //~ * @param[in] time      Simulation time to print.
 //~ * @return void
 //~ */
//~ static inline void diag_write_time_comment(FILE *file_ptr, const REAL time) { fprintf(file_ptr, "# [time] = %.15e\n", time); }



//NEW
#define DIAG_TIME_COMMENT_FMT "# [time] = %.15e\n"
static inline void diag_write_time_comment(FILE *file_ptr, const REAL time) {
  fprintf(file_ptr, DIAG_TIME_COMMENT_FMT, time);
}
// Return the number of bytes written by diag_write_time_comment for this time
static inline int diag_time_comment_size_bytes(const REAL time) {
  int n = snprintf(NULL, 0, DIAG_TIME_COMMENT_FMT, time);
  return (n < 0) ? 0 : n;
}


/**
 * @brief Write a standardized header line to a diagnostic output file.
 *
 * The header starts with the coordinate column names and then appends
 * gridfunction labels with their 1-based column indices.
 *
 * @param[in] file_ptr              File pointer to write to.
 * @param[in] coord_names           Space-separated coordinate column names (e.g., "x", "y", "x y").
 * @param[in] NUM_GFS               Number of gridfunction columns.
 * @param[in] which_gfs             Array mapping output order to diagnostic_gf_names indices (length NUM_GFS).
 * @param[in] diagnostic_gf_names   Array of human-readable diagnostic gridfunction names.
 * @return void
 */
static inline void diag_write_header(FILE *file_ptr, const char *coord_names, const int NUM_GFS, const int which_gfs[],
                                     const char **diagnostic_gf_names) {
  int num_spaces_in_coord_names = 0;
  fprintf(file_ptr, "# %s", coord_names); // header prefix with coordinates
  // Count spaces in coord_names to determine the starting column index for gridfunctions.
  for (const char *cp = coord_names; *cp; ++cp)
    num_spaces_in_coord_names += (*cp == ' ');
  for (int i0 = 0; i0 < NUM_GFS; i0++) {
    fprintf(file_ptr, " %s(%d)", diagnostic_gf_names[which_gfs[i0]], i0 + (num_spaces_in_coord_names + 2));
  } // END LOOP over gf columns
  fprintf(file_ptr, "\n");
} // END FUNCTION diag_write_header

/**
 * @brief Write a single row of numeric data to a diagnostic output file.
 *
 * Values are printed in scientific notation with 15 digits after the decimal.
 *
 * @param[in] file_ptr   File pointer to write to.
 * @param[in] NUM_COLS   Total number of columns in the row.
 * @param[in] data       Array of column data of length NUM_COLS.
 * @return void
 */
static inline void diag_write_row(FILE *file_ptr, const int NUM_COLS, const REAL data[]) {
  fprintf(file_ptr, "%.15e", data[0]);
  for (int i0 = 1; i0 < NUM_COLS; i0++) {
    fprintf(file_ptr, " %.15e", data[i0]);
  } // END LOOP over columns
  fprintf(file_ptr, "\n");
} // END FUNCTION diag_write_row


/**
 * @brief Build the diagnostic filename using the standard naming convention.
 *
 * When @p include_time is nonzero, the current simulation time is embedded
 * in the filename.
 *
 * Filename patterns:
 *   With time:     <prefix>-<coordsys>-conv_factor%.2f-t%08.2f.txt
 *   Without time:  <prefix>-<coordsys>-conv_factor-%.2f.txt
 *
 * @param[out] filename     Destination buffer for the filename.
 * @param[in]  filename_sz  Size of @p filename in bytes.
 * @param[in]  prefix       Logical prefix (e.g., "out1d-z").
 * @param[in]  coordsys     Coordinate system tag (e.g., "sinhcartesian").
 * @param[in]  commondata   Global simulation metadata used for filename parameters.
 * @param[in]  include_time Nonzero to embed the current simulation time in the filename.
 *
 * @return int Number of characters written (excluding the null byte).
 *             Note: If the return value is >= @p filename_sz, the output was truncated.
 */
static inline int build_outfile_name(char *filename, size_t filename_sz,
                                     const char *prefix, const char *coordsys,
                                     const commondata_struct *restrict commondata,
                                     int include_time) {
  if (include_time) {
    return snprintf(filename, filename_sz,
                    "%s-%s-conv_factor%.2f-t%08.2f.txt",
                    prefix, coordsys,
                    commondata->convergence_factor, commondata->time);
  } else {
    return snprintf(filename, filename_sz,
                    "%s-%s-conv_factor-%.2f.txt",
                    prefix, coordsys,
                    commondata->convergence_factor);
  }
}

/**
 * @brief Open a diagnostic output file using the standard naming convention.
 *
 * The file is opened in write mode ("w") if this is the first timestep
 * (commondata->nn == 0); otherwise it is opened in append mode ("a").
 * When include_time is nonzero, the current simulation time is embedded in the filename.
 *
 * Filename patterns:
 *   With time:     <prefix>-<coordsys>-conv_factor%.2f-t%08.2f.txt
 *   Without time:  <prefix>-<coordsys>-conv_factor-%.2f.txt
 *
 * @param[in] prefix        Logical prefix (e.g., "out1d-z").
 * @param[in] coordsys      Coordinate system tag (e.g., "sinhcartesian").
 * @param[in] commondata    Global simulation metadata used for filename parameters.
 * @param[in] include_time  Nonzero to embed the current simulation time in the filename.
 * @return FILE* Pointer to the opened file, or NULL on error.
 */
static inline FILE *open_outfile(const char *prefix, const char *coordsys, const commondata_struct *restrict commondata, int include_time) {
  char filename[256];

  //~ if (include_time) {
    //~ snprintf(filename, sizeof filename, "%s-%s-conv_factor%.2f-t%08.2f.txt", prefix, coordsys, commondata->convergence_factor, commondata->time);
  //~ } else {
    //~ snprintf(filename, sizeof filename, "%s-%s-conv_factor-%.2f.txt", prefix, coordsys, commondata->convergence_factor);
  //~ } // END IF include_time

  build_outfile_name(filename, sizeof filename, prefix, coordsys, commondata, include_time);

  FILE *file_ptr = (commondata->nn == 0) ? fopen(filename, "w") : fopen(filename, "a");
  if (!file_ptr) {
    fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
  } // END IF !file_ptr

  return file_ptr;
} // END FUNCTION open_outfile

#ifdef __cplusplus
} // END extern "C"
#endif

#endif // DIAGNOSTICS_NEAREST_COMMON_H
