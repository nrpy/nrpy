# nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/time_window_manager_numerical.py
"""
Define C helpers for numerical-spacetime time-window mmap management.

The NumericalTimeWindowManager owns the active mapped group of adjacent 3D-grid
payloads in one trusted combined numerical-spacetime binary container. It
depends on shared TimeSlotManager definitions for slot time bounds, but keeps
mmap and numerical file state out of analytic builds.

The key design constraint is that numerical-spacetime data is expensive to
load, while photons are already grouped by coordinate time through the
TimeSlotManager. This manager therefore chooses one contiguous range of
numerical 3D-grid time slices for an entire photon slot instead of loading
data ray-by-ray. The mapped range for one slot is intentionally conservative:
reverse ray tracing receives additional lower-time lookahead below the slot
boundary, and then both ends are widened by the centered temporal
interpolation halo so RKF45 steps can remain inside the active mapping while
photons cross slot boundaries.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.params as par
from nrpy.infrastructures.BHaH import BHaH_defines_h as Bdefines_h


def time_window_manager_numerical() -> None:
    """
    Register NumericalTimeWindowManager helpers in BHaH_defines.h.

    This registration owns the `rkf45_max_delta_t` CodeParameter consumed by
    the generated time-window C helpers. It also ensures the shared
    TimeSlotManager helper is registered first so the emitted code sees
    `TimeSlotManager`, `slot_lower_time()`, and `slot_upper_time()`, and so
    the commondata slot-lattice parameters already exist.

    The generated slot-window selector intentionally maps a conservative slice
    interval rather than the exact minimal interpolation stencil. The upper
    side only needs the centered temporal-interpolation halo, while the lower
    side also needs the additional backward-time RKF45 lookahead promised by
    `rkf45_max_delta_t`. Exact physical slice times are cached from the
    combined-file slice table and queried through binary search, so the
    mapped slot window remains correct even when stored output times are not
    uniformly spaced. The matching RKF45 step-size cap is enforced in the
    companion finalization kernel.

    Doctests:
    >>> time_window_manager_numerical()
    >>> generated = par.glb_extras_dict["BHaH_defines"]["time_window_manager_numerical"]
    >>> "time_window_manager_numerical_load_and_validate_slice_table" in generated
    True
    >>> "time_window_manager_numerical_lower_bound_slice_time" in generated
    True
    >>> "time_window_manager_numerical_nearest_slice_for_time" in generated
    True
    """
    from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.time_slot_manager_helpers import (  # pylint: disable=import-outside-toplevel
        time_slot_manager_helpers,
    )

    if "time_slot_manager" not in par.glb_extras_dict.get("BHaH_defines", {}):
        time_slot_manager_helpers()

    # Step 1: Register the runtime grid metadata fields this helper copies from
    # the combined numerical-spacetime container into params_struct.
    for dirn in range(3):
        _ = par.register_CodeParameter(
            "int",
            __name__,
            f"Nxx{dirn}",
            64,
            commondata=False,
            add_to_parfile=False,
        )
        _ = par.register_CodeParameter(
            "int",
            __name__,
            f"Nxx_plus_2NGHOSTS{dirn}",
            0,
            commondata=False,
            add_to_parfile=False,
        )
        _ = par.register_CodeParameter(
            "REAL",
            __name__,
            f"xxmin{dirn}",
            -10.0,
            commondata=False,
            add_to_parfile=False,
        )
        _ = par.register_CodeParameter(
            "REAL",
            __name__,
            f"xxmax{dirn}",
            10.0,
            commondata=False,
            add_to_parfile=False,
        )
        _ = par.register_CodeParameter(
            "REAL",
            __name__,
            f"dxx{dirn}",
            0.0,
            commondata=False,
            add_to_parfile=False,
        )
        _ = par.register_CodeParameter(
            "REAL",
            __name__,
            f"invdxx{dirn}",
            0.0,
            commondata=False,
            add_to_parfile=False,
        )

    # Step 2: Register the shared backward-time RKF45 lookahead parameter.
    _ = par.register_CodeParameter(
        "REAL",
        __name__,
        "rkf45_max_delta_t",
        10.0,
        commondata=True,
        add_to_parfile=True,
        description=(
            "Maximum backward coordinate-time lookahead used by the numerical "
            "time-window manager to keep reverse ray tracing inside the mapped "
            "time-slice window."
        ),
    )
    _ = par.register_CodeParameter(
        "REAL",
        __name__,
        "dt_numerical_spacetime_data",
        1.0,
        commondata=True,
        add_to_parfile=True,
        description=(
            "Uniform physical coordinate-time spacing used when one temporal "
            "interpolation stencil extends beyond the mapped numerical "
            "time-slice range and synthetic edge times must be supplied."
        ),
    )

    time_window_manager_numerical_c_code = r"""
    //==========================================
    // NUMERICAL TIME-WINDOW MANAGER
    //==========================================
    // Lightweight trusted-input mmap window over one combined numerical-spacetime container.
    // Time slices are assumed ordered by increasing coordinate time, while
    // reverse ray tracing moves from larger t toward smaller t. The manager
    // maps one contiguous slot-derived time window so many photons can reuse
    // the same numerical 3D-grid payloads. The mapped window is wider than
    // the slot itself because it includes temporal-interpolation halo slices
    // and lower-time RKF45 lookahead.
#include <fcntl.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>

#define TIME_WINDOW_MANAGER_NUMERICAL_FIXED_HEADER_BYTES 4096ULL
#define TIME_WINDOW_MANAGER_NUMERICAL_SLICE_TABLE_ENTRY_BYTES 96ULL

#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_FIXED_HEADER_BYTES 16U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_ALIGNMENT_BYTES 24U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_METADATA_OFFSET 32U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_METADATA_BYTES 40U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_SLICE_TABLE_OFFSET 48U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_SLICE_TABLE_ENTRY_BYTES 56U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_NUM_TIME_SLICES 64U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_FIRST_PAYLOAD_OFFSET 72U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_PAYLOAD_BYTES_TOTAL 80U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_TOTAL_FILE_BYTES 88U

#define TIME_WINDOW_MANAGER_NUMERICAL_SLICE_ENTRY_TIME 8U
#define TIME_WINDOW_MANAGER_NUMERICAL_SLICE_ENTRY_PAYLOAD_OFFSET 16U
#define TIME_WINDOW_MANAGER_NUMERICAL_SLICE_ENTRY_PAYLOAD_BYTES 24U
#define TIME_WINDOW_MANAGER_NUMERICAL_SLICE_ENTRY_POINT_RECORD_COUNT 32U
#define TIME_WINDOW_MANAGER_NUMERICAL_SLICE_ENTRY_POINT_RECORD_BYTES 40U

#define TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS 0
#define TIME_WINDOW_MANAGER_NUMERICAL_ERROR 1
#define TIME_WINDOW_MANAGER_NUMERICAL_MAX_TEMPORAL_INTERP_HALF_WIDTH 32

    // Owns the currently mapped group of adjacent 3D-grid payloads.
    typedef struct {
      int fd; // Open read-only file descriptor for the combined numerical .bin file.
      long int page_size; // Runtime page size used for mmap offset alignment.

      uint64_t alignment_bytes; // Combined-file payload alignment.
      uint64_t metadata_offset; // Absolute byte offset of the JSON metadata block.
      uint64_t metadata_bytes; // Length in bytes of the JSON metadata block.
      uint64_t num_time_slices; // Number of stored 3D grids.
      uint64_t slice_table_offset; // Absolute byte offset of the slice table.
      uint64_t slice_table_entry_bytes; // Fixed byte size of one slice-table entry.
      uint64_t first_payload_offset; // Absolute byte offset of slice 0 payload.
      uint64_t payload_bytes_total; // Total byte span of the payload region.
      uint64_t total_file_bytes; // Combined-file size recorded in the header.
      uint32_t nghosts; // Writer-side NGHOSTS stored in JSON metadata.
      uint64_t Nxx[3]; // Interior point counts stored in JSON metadata.
      uint64_t Nxx_plus_2NGHOSTS[3]; // Full logical-grid counts stored in JSON metadata.
      double dxx[3]; // Native coordinate spacings stored in JSON metadata.
      double invdxx[3]; // Reciprocal native coordinate spacings stored in JSON metadata.
      double xxmin[3]; // Native lower bounds stored in JSON metadata.
      double xxmax[3]; // Native upper bounds stored in JSON metadata.
      double cart_origin[3]; // Cartesian grid origin stored in JSON metadata.

      double *slice_times; // Exact physical coordinate times cached from the slice table.
      uint64_t *slice_payload_offsets; // Per-slice payload offsets cached from the slice table.
      uint64_t *slice_payload_bytes; // Per-slice payload byte counts cached from the slice table.

      int temporal_interp_half_width; // Centered temporal interpolation half-width n.
      int temporal_interp_num_points; // Number of temporal interpolation points, 2n + 1.
      double max_backward_dt_lookahead; // Maximum backward coordinate-time RKF45 lookahead from commondata.
      double dt_numerical_spacetime_data; // Uniform physical time spacing used to synthesize unmapped stencil-node times.

      uint64_t mapped_first_slice; // First mapped slice, inclusive.
      uint64_t mapped_last_slice_exclusive; // One past the final mapped slice.
      uint64_t mapped_file_offset; // Page-aligned file offset passed to mmap().
      size_t mapped_length_bytes; // Length passed to mmap()/munmap().
      const unsigned char *mapped_base; // Pointer returned by mmap(), or NULL.
    } NumericalTimeWindowManager; // END STRUCT: NumericalTimeWindowManager

    //==========================================
    // NUMERICAL WINDOW LOW-LEVEL READ HELPERS
    //==========================================
    // These helpers assume the trusted combined container uses native little-endian f64/u64 fields.
    static inline uint64_t time_window_manager_numerical_load_u64(const unsigned char *bytes) {
      uint64_t value;
      memcpy(&value, bytes, sizeof(uint64_t));
      return value;
    } // END FUNCTION: time_window_manager_numerical_load_u64

    static inline double time_window_manager_numerical_load_f64(const unsigned char *bytes) {
      double value;
      memcpy(&value, bytes, sizeof(double));
      return value;
    } // END FUNCTION: time_window_manager_numerical_load_f64

    static inline const char *time_window_manager_numerical_skip_json_ws(const char *cursor) {
      while (*cursor == ' ' || *cursor == '\n' || *cursor == '\r' || *cursor == '\t')
        cursor++;
      return cursor;
    } // END FUNCTION: time_window_manager_numerical_skip_json_ws

    /**
     * Read exactly nbytes from one absolute file offset.
     *
     * @param fd Open file descriptor.
     * @param[out] buffer Destination buffer.
     * @param nbytes Number of bytes to read.
     * @param offset Absolute file offset.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     */
    static inline int time_window_manager_numerical_read_exact_at(
        int fd, void *buffer, size_t nbytes, uint64_t offset) {
      unsigned char *dst = (unsigned char *)buffer;
      size_t bytes_read = 0U;

      while (bytes_read < nbytes) {
        const ssize_t nread = pread(
            fd, dst + bytes_read, nbytes - bytes_read,
            (off_t)(offset + bytes_read));
        if (nread <= 0)
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
        bytes_read += (size_t)nread;
      } // END WHILE: read requested metadata bytes
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_read_exact_at

    /**
     * Read one JSON section from the trusted metadata block.
     *
     * @param[in] json_buffer NUL-terminated metadata bytes.
     * @param metadata_bytes Number of valid bytes in json_buffer, excluding the NUL terminator.
     * @param[in] section_key Exact JSON-section prefix ending in '{'.
     * @param[out] section_begin First byte inside the section object.
     * @param[out] section_end Closing brace of the section object.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     */
    static inline int time_window_manager_numerical_find_json_section(
        const char *json_buffer,
        const size_t metadata_bytes,
        const char *section_key,
        const char **section_begin,
        const char **section_end) {
      const char *match = strstr(json_buffer, section_key);
      if (match == NULL || section_begin == NULL || section_end == NULL)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      const char *cursor = match + strlen(section_key);
      if ((size_t)(cursor - json_buffer) > metadata_bytes)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      *section_begin = cursor;

      int brace_depth = 1;
      while ((size_t)(cursor - json_buffer) < metadata_bytes) {
        if (*cursor == '{') {
          brace_depth++;
        } else if (*cursor == '}') {
          brace_depth--;
          if (brace_depth == 0) {
            *section_end = cursor;
            return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
          } // END IF: matching closing brace for the requested JSON section was found
        } // END ELSE IF: JSON brace nesting depth decreased while scanning the section
        cursor++;
      } // END WHILE: scanning the trusted metadata block for the section terminator
      return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
    } // END FUNCTION: time_window_manager_numerical_find_json_section

    /**
     * Parse one uint64_t scalar from a trusted JSON section.
     *
     * @param[in] section_begin First byte inside the section object.
     * @param[in] section_end Closing brace of the section object.
     * @param[in] key_prefix Exact key prefix ending in ':'.
     * @param[out] value Destination scalar.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     */
    static inline int time_window_manager_numerical_parse_json_u64(
        const char *section_begin,
        const char *section_end,
        const char *key_prefix,
        uint64_t *value) {
      const char *cursor = strstr(section_begin, key_prefix);
      if (cursor == NULL || cursor >= section_end || value == NULL)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      cursor = time_window_manager_numerical_skip_json_ws(cursor + strlen(key_prefix));

      errno = 0;
      char *parse_end = NULL;
      const unsigned long long parsed_value = strtoull(cursor, &parse_end, 10);
      if (errno != 0 || parse_end == cursor || parse_end > section_end)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      *value = (uint64_t)parsed_value;
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_parse_json_u64

    /**
     * Parse one uint64_t array of length three from a trusted JSON section.
     *
     * @param[in] section_begin First byte inside the section object.
     * @param[in] section_end Closing brace of the section object.
     * @param[in] key_prefix Exact key prefix ending in '['.
     * @param[out] values Destination array of length three.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     */
    static inline int time_window_manager_numerical_parse_json_u64_array3(
        const char *section_begin,
        const char *section_end,
        const char *key_prefix,
        uint64_t values[3]) {
      const char *cursor = strstr(section_begin, key_prefix);
      if (cursor == NULL || cursor >= section_end || values == NULL)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      cursor += strlen(key_prefix);

      for (int dirn = 0; dirn < 3; dirn++) {
        cursor = time_window_manager_numerical_skip_json_ws(cursor);
        errno = 0;
        char *parse_end = NULL;
        const unsigned long long parsed_value = strtoull(cursor, &parse_end, 10);
        if (errno != 0 || parse_end == cursor || parse_end > section_end)
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
        values[dirn] = (uint64_t)parsed_value;
        cursor = time_window_manager_numerical_skip_json_ws(parse_end);
        const char expected_delimiter = (dirn < 2) ? ',' : ']';
        if (*cursor != expected_delimiter)
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
        cursor++;
      } // END LOOP: for dirn over one trusted uint64_t JSON array
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_parse_json_u64_array3

    /**
     * Parse one double array of length three from a trusted JSON section.
     *
     * @param[in] section_begin First byte inside the section object.
     * @param[in] section_end Closing brace of the section object.
     * @param[in] key_prefix Exact key prefix ending in '['.
     * @param[out] values Destination array of length three.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     */
    static inline int time_window_manager_numerical_parse_json_f64_array3(
        const char *section_begin,
        const char *section_end,
        const char *key_prefix,
        double values[3]) {
      const char *cursor = strstr(section_begin, key_prefix);
      if (cursor == NULL || cursor >= section_end || values == NULL)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      cursor += strlen(key_prefix);

      for (int dirn = 0; dirn < 3; dirn++) {
        cursor = time_window_manager_numerical_skip_json_ws(cursor);
        errno = 0;
        char *parse_end = NULL;
        const double parsed_value = strtod(cursor, &parse_end);
        if (errno != 0 || parse_end == cursor || parse_end > section_end ||
            !isfinite(parsed_value)) {
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
        } // END IF: one trusted double JSON array entry could not be parsed
        values[dirn] = parsed_value;
        cursor = time_window_manager_numerical_skip_json_ws(parse_end);
        const char expected_delimiter = (dirn < 2) ? ',' : ']';
        if (*cursor != expected_delimiter)
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
        cursor++;
      } // END LOOP: for dirn over one trusted double JSON array
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_parse_json_f64_array3

    /**
     * Load first-slice grid metadata from the trusted JSON metadata block.
     *
     * @param[in,out] ntwm Numerical time-window manager.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     */
    static inline int time_window_manager_numerical_load_grid_metadata(
        NumericalTimeWindowManager *ntwm) {
      if (ntwm == NULL || ntwm->fd < 0 ||
          ntwm->metadata_bytes == 0ULL ||
          ntwm->metadata_bytes > (uint64_t)(SIZE_MAX - 1U)) {
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: metadata block dimensions were not usable

      const size_t buffer_bytes = (size_t)ntwm->metadata_bytes + 1U;
      char *metadata_buffer = (char *)malloc(buffer_bytes);
      if (metadata_buffer == NULL)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      if (time_window_manager_numerical_read_exact_at(
              ntwm->fd, metadata_buffer, (size_t)ntwm->metadata_bytes,
              ntwm->metadata_offset) != TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        free(metadata_buffer);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: JSON metadata block could not be read
      metadata_buffer[ntwm->metadata_bytes] = '\0';

      const char *section_begin = NULL;
      const char *section_end = NULL;
      if (time_window_manager_numerical_find_json_section(
              metadata_buffer, (size_t)ntwm->metadata_bytes,
              "\"first_slice_grid_metadata\":{",
              &section_begin, &section_end) !=
          TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        free(metadata_buffer);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: first-slice grid metadata was missing from the trusted JSON block

      uint64_t nghosts = 0ULL;
      if (time_window_manager_numerical_parse_json_u64(
              section_begin, section_end, "\"NGHOSTS\":", &nghosts) !=
          TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS ||
          nghosts > (uint64_t)UINT32_MAX) {
        free(metadata_buffer);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: first-slice NGHOSTS metadata was invalid
      ntwm->nghosts = (uint32_t)nghosts;

      if (time_window_manager_numerical_parse_json_u64_array3(
              section_begin, section_end, "\"Nxx\":[", ntwm->Nxx) !=
          TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS ||
          time_window_manager_numerical_parse_json_u64_array3(
              section_begin, section_end, "\"Nxx_plus_2NGHOSTS\":[",
              ntwm->Nxx_plus_2NGHOSTS) !=
              TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS ||
          time_window_manager_numerical_parse_json_f64_array3(
              section_begin, section_end, "\"dxx\":[", ntwm->dxx) !=
              TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS ||
          time_window_manager_numerical_parse_json_f64_array3(
              section_begin, section_end, "\"invdxx\":[", ntwm->invdxx) !=
              TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS ||
          time_window_manager_numerical_parse_json_f64_array3(
              section_begin, section_end, "\"xxmin\":[", ntwm->xxmin) !=
              TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS ||
          time_window_manager_numerical_parse_json_f64_array3(
              section_begin, section_end, "\"xxmax\":[", ntwm->xxmax) !=
              TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS ||
          time_window_manager_numerical_parse_json_f64_array3(
              section_begin, section_end, "\"cart_origin\":[", ntwm->cart_origin) !=
              TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        free(metadata_buffer);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: one or more trusted first-slice grid arrays could not be parsed

      free(metadata_buffer);
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_load_grid_metadata

    /**
     * Cache and validate the slice table fields needed for interpolation.
     *
     * @param[in,out] ntwm Numerical time-window manager.
     * @param expected_point_record_count Required point-record count for each payload.
     * @param expected_point_record_bytes Required point-record size for each payload.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     */
    static inline int time_window_manager_numerical_load_and_validate_slice_table(
        NumericalTimeWindowManager *ntwm,
        const uint64_t expected_point_record_count,
        const uint64_t expected_point_record_bytes) {
      unsigned char entry_bytes[TIME_WINDOW_MANAGER_NUMERICAL_SLICE_TABLE_ENTRY_BYTES];
      double *slice_times = NULL;
      uint64_t *slice_payload_offsets = NULL;
      uint64_t *slice_payload_bytes = NULL;
      double previous_time = 0.0;
      uint64_t previous_payload_end = 0ULL;

      if (ntwm == NULL || ntwm->fd < 0 || ntwm->num_time_slices < 2ULL ||
          ntwm->slice_times != NULL ||
          ntwm->num_time_slices > (uint64_t)(SIZE_MAX / sizeof(double)) ||
          ntwm->num_time_slices > (uint64_t)(SIZE_MAX / sizeof(uint64_t))) {
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: slice-table cache prerequisites were not satisfied

      slice_times = (double *)malloc((size_t)ntwm->num_time_slices * sizeof(double));
      slice_payload_offsets =
          (uint64_t *)malloc((size_t)ntwm->num_time_slices * sizeof(uint64_t));
      slice_payload_bytes =
          (uint64_t *)malloc((size_t)ntwm->num_time_slices * sizeof(uint64_t));
      if (slice_times == NULL || slice_payload_offsets == NULL ||
          slice_payload_bytes == NULL) {
        free(slice_times);
        free(slice_payload_offsets);
        free(slice_payload_bytes);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: slice-table caches could not be allocated

      for (uint64_t slice_index = 0ULL; slice_index < ntwm->num_time_slices; slice_index++) {
        const uint64_t entry_offset =
            ntwm->slice_table_offset + slice_index * ntwm->slice_table_entry_bytes;
        if (time_window_manager_numerical_read_exact_at(
                ntwm->fd, entry_bytes,
                (size_t)TIME_WINDOW_MANAGER_NUMERICAL_SLICE_TABLE_ENTRY_BYTES,
                entry_offset) != TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
          free(slice_times);
          free(slice_payload_offsets);
          free(slice_payload_bytes);
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
        } // END IF: one slice-table entry could not be read

        const double this_time = time_window_manager_numerical_load_f64(
            entry_bytes + TIME_WINDOW_MANAGER_NUMERICAL_SLICE_ENTRY_TIME);
        const uint64_t payload_offset = time_window_manager_numerical_load_u64(
            entry_bytes + TIME_WINDOW_MANAGER_NUMERICAL_SLICE_ENTRY_PAYLOAD_OFFSET);
        const uint64_t payload_bytes = time_window_manager_numerical_load_u64(
            entry_bytes + TIME_WINDOW_MANAGER_NUMERICAL_SLICE_ENTRY_PAYLOAD_BYTES);
        const uint64_t point_record_count = time_window_manager_numerical_load_u64(
            entry_bytes + TIME_WINDOW_MANAGER_NUMERICAL_SLICE_ENTRY_POINT_RECORD_COUNT);
        const uint64_t point_record_bytes = time_window_manager_numerical_load_u64(
            entry_bytes + TIME_WINDOW_MANAGER_NUMERICAL_SLICE_ENTRY_POINT_RECORD_BYTES);

        if (!isfinite(this_time) ||
            point_record_count == 0ULL ||
            point_record_bytes != expected_point_record_bytes ||
            point_record_count != expected_point_record_count ||
            payload_bytes == 0ULL ||
            payload_offset < ntwm->first_payload_offset ||
            payload_offset % ntwm->alignment_bytes != 0ULL) {
          free(slice_times);
          free(slice_payload_offsets);
          free(slice_payload_bytes);
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
        } // END IF: one slice-table entry contained invalid interpolation metadata
        if (point_record_count > UINT64_MAX / point_record_bytes ||
            payload_bytes != point_record_count * point_record_bytes) {
          free(slice_times);
          free(slice_payload_offsets);
          free(slice_payload_bytes);
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
        } // END IF: one slice payload size disagreed with point-record metadata
        if (payload_offset > ntwm->total_file_bytes - payload_bytes) {
          free(slice_times);
          free(slice_payload_offsets);
          free(slice_payload_bytes);
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
        } // END IF: one slice payload exceeded the combined container bounds
        if (slice_index == 0ULL) {
          if (payload_offset != ntwm->first_payload_offset) {
            free(slice_times);
            free(slice_payload_offsets);
            free(slice_payload_bytes);
            return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
          } // END IF: first slice-table entry disagreed with fixed-header navigation
        } else {
          const double slice_dt = this_time - previous_time;
          if (!isfinite(slice_dt) || slice_dt <= 0.0 ||
              payload_offset < previous_payload_end) {
            free(slice_times);
            free(slice_payload_offsets);
            free(slice_payload_bytes);
            return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
          } // END IF: slice times or payload offsets were not strictly increasing
        } // END ELSE: additional slices allow monotonicity checks

        slice_times[slice_index] = this_time;
        slice_payload_offsets[slice_index] = payload_offset;
        slice_payload_bytes[slice_index] = payload_bytes;
        previous_time = this_time;
        previous_payload_end = payload_offset + payload_bytes;
      } // END LOOP: for slice_index over the trusted slice table

      if (!isfinite(slice_times[0]) ||
          !isfinite(slice_times[ntwm->num_time_slices - 1ULL]) ||
          slice_times[ntwm->num_time_slices - 1ULL] <= slice_times[0] ||
          previous_payload_end != ntwm->total_file_bytes) {
        free(slice_times);
        free(slice_payload_offsets);
        free(slice_payload_bytes);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: cached slice-table summary values were not numerically usable

      ntwm->slice_times = slice_times;
      ntwm->slice_payload_offsets = slice_payload_offsets;
      ntwm->slice_payload_bytes = slice_payload_bytes;
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_load_and_validate_slice_table

    /**
     * Return the first slice index whose cached time is not less than target_time.
     *
     * @param[in] ntwm Numerical time-window manager with cached exact slice times.
     * @param target_time Physical coordinate time to bracket.
     * @return Lower-bound slice index in [0, num_time_slices].
     */
    static inline uint64_t time_window_manager_numerical_lower_bound_slice_time(
        const NumericalTimeWindowManager *ntwm,
        const double target_time) {
      uint64_t first = 0ULL;
      uint64_t last = ntwm->num_time_slices;

      while (first < last) {
        const uint64_t mid = first + (last - first) / 2ULL;
        if (ntwm->slice_times[mid] < target_time) {
          first = mid + 1ULL;
        } else {
          last = mid;
        } // END ELSE: lower-bound search moved its upper half-open endpoint
      } // END WHILE: binary searching the lower-bound cached slice time
      return first;
    } // END FUNCTION: time_window_manager_numerical_lower_bound_slice_time

    /**
     * Return the first slice index whose cached time is greater than target_time.
     *
     * @param[in] ntwm Numerical time-window manager with cached exact slice times.
     * @param target_time Physical coordinate time to bracket.
     * @return Upper-bound slice index in [0, num_time_slices].
     */
    static inline uint64_t time_window_manager_numerical_upper_bound_slice_time(
        const NumericalTimeWindowManager *ntwm,
        const double target_time) {
      uint64_t first = 0ULL;
      uint64_t last = ntwm->num_time_slices;

      while (first < last) {
        const uint64_t mid = first + (last - first) / 2ULL;
        if (ntwm->slice_times[mid] <= target_time) {
          first = mid + 1ULL;
        } else {
          last = mid;
        } // END ELSE: upper-bound search moved its upper half-open endpoint
      } // END WHILE: binary searching the upper-bound cached slice time
      return first;
    } // END FUNCTION: time_window_manager_numerical_upper_bound_slice_time

    /**
     * Return the cached slice index nearest to one physical target time.
     *
     * @param[in] ntwm Numerical time-window manager with cached exact slice times.
     * @param target_time Physical coordinate time to bracket.
     * @return Nearest cached slice index.
     */
    static inline uint64_t time_window_manager_numerical_nearest_slice_for_time(
        const NumericalTimeWindowManager *ntwm,
        const double target_time) {
      const uint64_t right =
          time_window_manager_numerical_lower_bound_slice_time(ntwm, target_time);

      if (right == 0ULL)
        return 0ULL;
      if (right >= ntwm->num_time_slices)
        return ntwm->num_time_slices - 1ULL;

      const double dt_right = fabs(ntwm->slice_times[right] - target_time);
      const double dt_left = fabs(target_time - ntwm->slice_times[right - 1ULL]);
      if (dt_left <= dt_right)
        return right - 1ULL;
      return right;
    } // END FUNCTION: time_window_manager_numerical_nearest_slice_for_time

    //==========================================
    // NUMERICAL WINDOW LIFETIME MANAGEMENT
    //==========================================
    /**
     * Reset a numerical time-window manager to an inert state.
     *
     * @param[in,out] ntwm Numerical time-window manager to reset.
     */
    static inline void time_window_manager_numerical_set_inert(NumericalTimeWindowManager *ntwm) {
      ntwm->fd = -1;
      ntwm->page_size = 0;
      ntwm->alignment_bytes = 0ULL;
      ntwm->metadata_offset = 0ULL;
      ntwm->metadata_bytes = 0ULL;
      ntwm->num_time_slices = 0ULL;
      ntwm->slice_table_offset = 0ULL;
      ntwm->slice_table_entry_bytes = 0ULL;
      ntwm->first_payload_offset = 0ULL;
      ntwm->payload_bytes_total = 0ULL;
      ntwm->total_file_bytes = 0ULL;
      ntwm->nghosts = 0U;
      for (int dirn = 0; dirn < 3; dirn++) {
        ntwm->Nxx[dirn] = 0ULL;
        ntwm->Nxx_plus_2NGHOSTS[dirn] = 0ULL;
        ntwm->dxx[dirn] = 0.0;
        ntwm->invdxx[dirn] = 0.0;
        ntwm->xxmin[dirn] = 0.0;
        ntwm->xxmax[dirn] = 0.0;
        ntwm->cart_origin[dirn] = 0.0;
      } // END LOOP: for dirn over stored grid metadata arrays during inert reset
      ntwm->slice_times = NULL;
      ntwm->slice_payload_offsets = NULL;
      ntwm->slice_payload_bytes = NULL;
      ntwm->temporal_interp_half_width = 0;
      ntwm->temporal_interp_num_points = 0;
      ntwm->max_backward_dt_lookahead = 0.0;
      ntwm->dt_numerical_spacetime_data = 0.0;
      ntwm->mapped_first_slice = 0ULL;
      ntwm->mapped_last_slice_exclusive = 0ULL;
      ntwm->mapped_file_offset = 0ULL;
      ntwm->mapped_length_bytes = 0U;
      ntwm->mapped_base = NULL;
    } // END FUNCTION: time_window_manager_numerical_set_inert

    /**
     * Close the active mmap window.
     *
     * @param[in,out] ntwm Numerical time-window manager.
     */
    static inline void time_window_manager_numerical_close_mmap(NumericalTimeWindowManager *ntwm) {
      if (ntwm->mapped_base != NULL)
        munmap((void *)ntwm->mapped_base, ntwm->mapped_length_bytes);
      ntwm->mapped_first_slice = 0ULL;
      ntwm->mapped_last_slice_exclusive = 0ULL;
      ntwm->mapped_file_offset = 0ULL;
      ntwm->mapped_length_bytes = 0U;
      ntwm->mapped_base = NULL;
    } // END FUNCTION: time_window_manager_numerical_close_mmap

    /**
     * Release all resources owned by the numerical time-window manager.
     *
     * @param[in,out] ntwm Numerical time-window manager.
     */
    static inline void time_window_manager_numerical_free(NumericalTimeWindowManager *ntwm) {
      time_window_manager_numerical_close_mmap(ntwm);
      if (ntwm->slice_times != NULL)
        free(ntwm->slice_times);
      if (ntwm->slice_payload_offsets != NULL)
        free(ntwm->slice_payload_offsets);
      if (ntwm->slice_payload_bytes != NULL)
        free(ntwm->slice_payload_bytes);
      if (ntwm->fd >= 0)
        close(ntwm->fd);
      time_window_manager_numerical_set_inert(ntwm);
    } // END FUNCTION: time_window_manager_numerical_free

    /**
     * Copy combined-file grid metadata into a runtime params_struct.
     *
     * Only the grid fields serialized by the numerical-data writer are
     * overwritten here. Higher-level caller-owned params metadata unrelated to
     * the numerical dataset geometry remains untouched.
     *
     * @param[in] ntwm Initialized numerical time-window manager.
     * @param[out] params Destination params_struct to populate.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     */
    static inline int time_window_manager_numerical_apply_metadata_to_params(
        const NumericalTimeWindowManager *ntwm,
        params_struct *restrict params) {
      if (ntwm == NULL || params == NULL)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      if (ntwm->nghosts != (uint32_t)NGHOSTS)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      for (int dirn = 0; dirn < 3; dirn++) {
        if (ntwm->Nxx[dirn] > (uint64_t)INT_MAX ||
            ntwm->Nxx_plus_2NGHOSTS[dirn] > (uint64_t)INT_MAX)
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END LOOP: for dirn over integer-valued grid metadata before params assignment

      params->Nxx0 = (int)ntwm->Nxx[0];
      params->Nxx1 = (int)ntwm->Nxx[1];
      params->Nxx2 = (int)ntwm->Nxx[2];
      params->Nxx_plus_2NGHOSTS0 = (int)ntwm->Nxx_plus_2NGHOSTS[0];
      params->Nxx_plus_2NGHOSTS1 = (int)ntwm->Nxx_plus_2NGHOSTS[1];
      params->Nxx_plus_2NGHOSTS2 = (int)ntwm->Nxx_plus_2NGHOSTS[2];
      params->xxmin0 = (REAL)ntwm->xxmin[0];
      params->xxmin1 = (REAL)ntwm->xxmin[1];
      params->xxmin2 = (REAL)ntwm->xxmin[2];
      params->xxmax0 = (REAL)ntwm->xxmax[0];
      params->xxmax1 = (REAL)ntwm->xxmax[1];
      params->xxmax2 = (REAL)ntwm->xxmax[2];
      params->dxx0 = (REAL)ntwm->dxx[0];
      params->dxx1 = (REAL)ntwm->dxx[1];
      params->dxx2 = (REAL)ntwm->dxx[2];
      params->invdxx0 = (REAL)ntwm->invdxx[0];
      params->invdxx1 = (REAL)ntwm->invdxx[1];
      params->invdxx2 = (REAL)ntwm->invdxx[2];
      params->Cart_originx = (REAL)ntwm->cart_origin[0];
      params->Cart_originy = (REAL)ntwm->cart_origin[1];
      params->Cart_originz = (REAL)ntwm->cart_origin[2];
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_apply_metadata_to_params

    /**
     * Initialize a trusted-input numerical time-window manager.
     *
     * The caller must call time_window_manager_numerical_set_inert() before
     * the first initialization. Reinitialization releases the old mapping and
     * file before opening the replacement container.
     *
     * @param[in,out] ntwm Numerical time-window manager.
     * @param[in] combined_path Path to the combined numerical .bin container.
     * @param[in] commondata Common runtime parameters.
     * @param temporal_interp_half_width Centered temporal interpolation half-width n.
     * @param[out] params Optional params_struct updated from the combined-file grid metadata.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     *
     * @note The backward-lookahead distance is read from
     * commondata->rkf45_max_delta_t.
     */
    static inline int time_window_manager_numerical_init(
        NumericalTimeWindowManager *ntwm,
        const char *combined_path,
        const commondata_struct *restrict commondata,
        const int temporal_interp_half_width,
        params_struct *restrict params) {
      unsigned char header_bytes[TIME_WINDOW_MANAGER_NUMERICAL_FIXED_HEADER_BYTES];
      static const unsigned char expected_magic[16] = "NRPYRTSTACK4D";
      const uint64_t expected_point_record_bytes = 53ULL * (uint64_t)sizeof(double);
      uint64_t expected_point_record_count = 1ULL;

      if (ntwm == NULL)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      time_window_manager_numerical_free(ntwm);
      if (combined_path == NULL || commondata == NULL)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      if (temporal_interp_half_width < 0 ||
          temporal_interp_half_width >
              TIME_WINDOW_MANAGER_NUMERICAL_MAX_TEMPORAL_INTERP_HALF_WIDTH)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      const double max_backward_dt_lookahead = commondata->rkf45_max_delta_t;
      if (!isfinite(max_backward_dt_lookahead) || max_backward_dt_lookahead < 0.0)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      const double dt_numerical_spacetime_data =
          commondata->dt_numerical_spacetime_data;
      if (!isfinite(dt_numerical_spacetime_data) ||
          dt_numerical_spacetime_data <= 0.0) {
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: synthetic temporal edge spacing was not usable

      const long int page_size = sysconf(_SC_PAGESIZE);
      if (page_size <= 0)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      ntwm->fd = open(combined_path, O_RDONLY);
      if (ntwm->fd < 0)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      ntwm->temporal_interp_half_width = temporal_interp_half_width;
      ntwm->temporal_interp_num_points = 2 * temporal_interp_half_width + 1;
      ntwm->max_backward_dt_lookahead = max_backward_dt_lookahead;
      ntwm->dt_numerical_spacetime_data = dt_numerical_spacetime_data;
      ntwm->page_size = page_size;

      if (time_window_manager_numerical_read_exact_at(
              ntwm->fd, header_bytes,
              (size_t)TIME_WINDOW_MANAGER_NUMERICAL_FIXED_HEADER_BYTES, 0ULL) !=
          TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: fixed header could not be read

      const uint64_t fixed_header_bytes = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_FIXED_HEADER_BYTES);
      ntwm->alignment_bytes = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_ALIGNMENT_BYTES);
      ntwm->metadata_offset = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_METADATA_OFFSET);
      ntwm->metadata_bytes = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_METADATA_BYTES);
      ntwm->slice_table_offset = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_SLICE_TABLE_OFFSET);
      ntwm->slice_table_entry_bytes = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_SLICE_TABLE_ENTRY_BYTES);
      ntwm->num_time_slices = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_NUM_TIME_SLICES);
      ntwm->first_payload_offset = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_FIRST_PAYLOAD_OFFSET);
      ntwm->payload_bytes_total = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_PAYLOAD_BYTES_TOTAL);
      ntwm->total_file_bytes = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_TOTAL_FILE_BYTES);

      if (ntwm->num_time_slices < 2ULL ||
          ntwm->num_time_slices < (uint64_t)ntwm->temporal_interp_num_points ||
          ntwm->num_time_slices > (uint64_t)LONG_MAX ||
          ntwm->alignment_bytes == 0ULL ||
          (ntwm->alignment_bytes & (ntwm->alignment_bytes - 1ULL)) != 0ULL ||
          ntwm->metadata_bytes == 0ULL ||
          ntwm->slice_table_entry_bytes != TIME_WINDOW_MANAGER_NUMERICAL_SLICE_TABLE_ENTRY_BYTES ||
          ntwm->payload_bytes_total == 0ULL ||
          fixed_header_bytes != TIME_WINDOW_MANAGER_NUMERICAL_FIXED_HEADER_BYTES ||
          memcmp(header_bytes, expected_magic, sizeof(expected_magic)) != 0) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted fixed-header navigation fields were not usable

      const off_t file_size = lseek(ntwm->fd, 0, SEEK_END);
      if (file_size < (off_t)TIME_WINDOW_MANAGER_NUMERICAL_FIXED_HEADER_BYTES) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: combined numerical container size was not usable
      const uint64_t file_size_bytes = (uint64_t)file_size;
      if (file_size_bytes != ntwm->total_file_bytes) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: on-disk file size disagreed with the combined-header size

      if (ntwm->metadata_offset < fixed_header_bytes ||
          ntwm->metadata_offset % ntwm->alignment_bytes != 0ULL ||
          ntwm->metadata_offset > ntwm->total_file_bytes - ntwm->metadata_bytes) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: metadata block bounds were outside the combined container

      if (ntwm->num_time_slices >
          UINT64_MAX / ntwm->slice_table_entry_bytes) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted slice-table byte count would overflow
      const uint64_t slice_table_bytes =
          ntwm->num_time_slices * ntwm->slice_table_entry_bytes;
      const uint64_t metadata_end = ntwm->metadata_offset + ntwm->metadata_bytes;
      if (ntwm->slice_table_offset > UINT64_MAX - slice_table_bytes) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted slice-table end offset would overflow
      const uint64_t slice_table_end =
          ntwm->slice_table_offset + slice_table_bytes;
      if (ntwm->slice_table_offset < metadata_end ||
          ntwm->slice_table_offset % ntwm->alignment_bytes != 0ULL ||
          slice_table_end > ntwm->first_payload_offset ||
          slice_table_end > file_size_bytes) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted slice table was outside the usable file layout

      if (ntwm->first_payload_offset < slice_table_end ||
          ntwm->first_payload_offset % ntwm->alignment_bytes != 0ULL ||
          ntwm->first_payload_offset > ntwm->total_file_bytes ||
          ntwm->payload_bytes_total != ntwm->total_file_bytes - ntwm->first_payload_offset) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted payload-region navigation was not self-consistent

      if (time_window_manager_numerical_load_grid_metadata(ntwm) !=
          TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: JSON grid metadata could not be loaded for interpolation

      for (int dirn = 0; dirn < 3; dirn++) {
        if (ntwm->Nxx[dirn] == 0ULL ||
            ntwm->Nxx_plus_2NGHOSTS[dirn] !=
                ntwm->Nxx[dirn] + 2ULL * (uint64_t)ntwm->nghosts ||
            ntwm->Nxx_plus_2NGHOSTS[dirn] >
                UINT64_MAX / expected_point_record_count ||
            !isfinite(ntwm->dxx[dirn]) || ntwm->dxx[dirn] <= 0.0 ||
            !isfinite(ntwm->invdxx[dirn]) || ntwm->invdxx[dirn] <= 0.0 ||
            !isfinite(ntwm->xxmin[dirn]) || !isfinite(ntwm->xxmax[dirn]) ||
            ntwm->xxmax[dirn] <= ntwm->xxmin[dirn] ||
            !isfinite(ntwm->cart_origin[dirn])) {
          time_window_manager_numerical_free(ntwm);
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
        } // END IF: JSON grid metadata was not numerically usable
        expected_point_record_count *= ntwm->Nxx_plus_2NGHOSTS[dirn];
      } // END LOOP: for dirn over JSON grid metadata during validation
      if (ntwm->nghosts != (uint32_t)NGHOSTS ||
          expected_point_record_count > UINT64_MAX / expected_point_record_bytes) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: runtime ghost-zone or record-size metadata was not usable

      if (time_window_manager_numerical_load_and_validate_slice_table(
              ntwm, expected_point_record_count, expected_point_record_bytes) !=
          TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: slice table could not be cached and validated for interpolation
      if (params != NULL &&
          time_window_manager_numerical_apply_metadata_to_params(ntwm, params) !=
              TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: caller requested params_struct population from combined metadata
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_init

    //==========================================
    // NUMERICAL WINDOW SLOT RANGE SELECTION
    //==========================================
    /**
     * Compute the conservative mapped time-slice range for one TimeSlotManager slot.
     *
     * @param[in] ntwm Initialized numerical time-window manager.
     * @param[in] tsm Base time-slot manager.
     * @param slot_idx Slot index.
     * @param[out] first_slice First mapped slice, inclusive.
     * @param[out] last_slice_exclusive One past the final mapped slice.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     *
     * @note This mapped interval is intentionally wider than the exact slice
     * set requested later by time_window_manager_numerical_stencil_for_time().
     * The upper side covers the centered temporal-interpolation halo above the
     * slot, while the lower side covers both the centered halo and the extra
     * backward-time RKF45 lookahead stored in
     * ntwm->max_backward_dt_lookahead. This asymmetry matches backward
     * raytracing, where photons move toward lower coordinate time.
     */
    static inline int time_window_manager_numerical_required_grid_range(
        const NumericalTimeWindowManager *ntwm,
        const TimeSlotManager *tsm,
        const int slot_idx,
        uint64_t *first_slice,
        uint64_t *last_slice_exclusive) {
      if (ntwm == NULL || tsm == NULL || first_slice == NULL ||
          last_slice_exclusive == NULL || ntwm->slice_times == NULL) {
        if (first_slice != NULL)
          *first_slice = 0ULL;
        if (last_slice_exclusive != NULL)
          *last_slice_exclusive = 0ULL;
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: slot-range query was missing the cached exact-time window state

      const double t_lower = slot_lower_time(tsm, slot_idx);
      const double t_upper = slot_upper_time(tsm, slot_idx);
      // The mapped query is biased toward lower coordinate time because the
      // photons are evolved backward in time.
      const double query_min = t_lower - ntwm->max_backward_dt_lookahead;
      const double query_max = t_upper;
      if (!isfinite(query_min) || !isfinite(query_max) || query_max < query_min) {
        *first_slice = 0ULL;
        *last_slice_exclusive = 0ULL;
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: slot query bounds or cached exact slice times were not usable

      const uint64_t left =
          time_window_manager_numerical_lower_bound_slice_time(ntwm, query_min);
      const uint64_t right =
          time_window_manager_numerical_upper_bound_slice_time(ntwm, query_max);
      uint64_t center_first = 0ULL;
      uint64_t center_last = 0ULL;

      if (left == 0ULL) {
        center_first = 0ULL;
      } else if (left >= ntwm->num_time_slices) {
        center_first = ntwm->num_time_slices - 1ULL;
      } else {
        center_first = left - 1ULL;
      } // END ELSE: cached lower-bound time had one slice below it

      if (right == 0ULL) {
        center_last = 0ULL;
      } else if (right >= ntwm->num_time_slices) {
        center_last = ntwm->num_time_slices - 1ULL;
      } else {
        center_last = right;
      } // END ELSE: cached upper-bound time had one slice above it

      // Conceptually, for temporal half-width n, centered interpolation uses
      // 2n + 1 slices. A photon anywhere in the slot may therefore need
      // about n + 1 halo slices outside the queried coordinate-time
      // interval. Since photons evolve backward in time, the queried
      // interval is extended below the slot before the halo is applied.
      const uint64_t halo =
          (uint64_t)ntwm->temporal_interp_half_width + 1ULL;
      const uint64_t required_first =
          (center_first < halo) ? 0ULL : center_first - halo;
      uint64_t required_last_exclusive = ntwm->num_time_slices;

      if (center_last < ntwm->num_time_slices - 1ULL) {
        const uint64_t unclipped_last_exclusive = center_last + halo + 1ULL;
        required_last_exclusive =
            (unclipped_last_exclusive > ntwm->num_time_slices)
                ? ntwm->num_time_slices
                : unclipped_last_exclusive;
      } // END IF: upper conservative halo endpoint remained inside the stored slice range

      if (required_last_exclusive > ntwm->num_time_slices ||
          required_last_exclusive <= required_first) {
        *first_slice = 0ULL;
        *last_slice_exclusive = 0ULL;
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: clipped conservative slot mapping window was outside the trusted file

      *first_slice = required_first;
      *last_slice_exclusive = required_last_exclusive;
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_required_grid_range

    /**
     * Map exactly one slot's required contiguous time-slice payload window.
     *
     * @param[in,out] ntwm Numerical time-window manager.
     * @param[in] tsm Base time-slot manager.
     * @param slot_idx Slot index.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     *
     * @pre No interpolation reader may hold pointers returned by
     * time_window_manager_numerical_grid_ptr during this call.
     */
    static inline int time_window_manager_numerical_mmap_for_slot(
        NumericalTimeWindowManager *ntwm,
        const TimeSlotManager *tsm,
        const int slot_idx) {
      if (ntwm == NULL || tsm == NULL || ntwm->fd < 0 || ntwm->page_size <= 0 ||
          ntwm->slice_times == NULL ||
          ntwm->slice_payload_offsets == NULL ||
          ntwm->slice_payload_bytes == NULL ||
          ntwm->num_time_slices < 2ULL ||
          ntwm->num_time_slices < (uint64_t)ntwm->temporal_interp_num_points ||
          ntwm->payload_bytes_total == 0ULL ||
          !isfinite(ntwm->slice_times[0]) ||
          !isfinite(ntwm->slice_times[ntwm->num_time_slices - 1ULL]) ||
          ntwm->slice_times[ntwm->num_time_slices - 1ULL] <= ntwm->slice_times[0])
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      uint64_t first_slice = 0ULL;
      uint64_t last_slice_exclusive = 0ULL;
      if (time_window_manager_numerical_required_grid_range(
              ntwm, tsm, slot_idx, &first_slice, &last_slice_exclusive) !=
          TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      if (ntwm->mapped_base != NULL && ntwm->mapped_first_slice == first_slice &&
          ntwm->mapped_last_slice_exclusive == last_slice_exclusive)
        return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;

      if (last_slice_exclusive <= first_slice ||
          last_slice_exclusive > ntwm->num_time_slices)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      const uint64_t last_slice = last_slice_exclusive - 1ULL;
      const uint64_t requested_offset = ntwm->slice_payload_offsets[first_slice];
      const uint64_t last_payload_offset = ntwm->slice_payload_offsets[last_slice];
      const uint64_t last_payload_bytes = ntwm->slice_payload_bytes[last_slice];
      if (requested_offset < ntwm->first_payload_offset ||
          last_payload_offset < requested_offset ||
          last_payload_bytes == 0ULL ||
          last_payload_offset > ntwm->total_file_bytes - last_payload_bytes) {
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: requested slot payload window was outside the trusted file

      const uint64_t requested_end =
          last_payload_offset + last_payload_bytes;
      const uint64_t page_size = (uint64_t)ntwm->page_size;
      const uint64_t map_offset =
          requested_offset - (requested_offset % page_size);
      if (requested_end <= map_offset)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      const uint64_t map_length_u64 = requested_end - map_offset;
      if (map_length_u64 > (uint64_t)SIZE_MAX)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      const size_t map_length = (size_t)map_length_u64;

      time_window_manager_numerical_close_mmap(ntwm);
      void *new_mapping = mmap(NULL, map_length, PROT_READ, MAP_PRIVATE, ntwm->fd,
                               (off_t)map_offset);
      if (new_mapping == MAP_FAILED)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      ntwm->mapped_first_slice = first_slice;
      ntwm->mapped_last_slice_exclusive = last_slice_exclusive;
      ntwm->mapped_file_offset = map_offset;
      ntwm->mapped_length_bytes = map_length;
      ntwm->mapped_base = (const unsigned char *)new_mapping;
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_mmap_for_slot

    /**
     * Return a pointer to one mapped 3D-grid payload.
     *
     * @param[in] ntwm Numerical time-window manager.
     * @param slice_index Global slice index.
     * @return Pointer to the slice payload.
     *
     * @pre slice_index is inside the currently mapped window.
     */
    static inline const double *time_window_manager_numerical_grid_ptr(
        const NumericalTimeWindowManager *ntwm,
        const uint64_t slice_index) {
      if (ntwm == NULL || ntwm->mapped_base == NULL ||
          ntwm->slice_payload_offsets == NULL ||
          ntwm->slice_payload_bytes == NULL ||
          slice_index < ntwm->mapped_first_slice ||
          slice_index >= ntwm->mapped_last_slice_exclusive) {
        return NULL;
      } // END IF: one requested slice was outside the active mapped payload window

      const uint64_t payload_offset = ntwm->slice_payload_offsets[slice_index];
      const uint64_t payload_bytes = ntwm->slice_payload_bytes[slice_index];
      if (payload_offset < ntwm->mapped_file_offset ||
          payload_bytes == 0ULL ||
          payload_offset > UINT64_MAX - payload_bytes) {
        return NULL;
      } // END IF: cached payload metadata could not describe the mapped slice range
      const uint64_t payload_relative = payload_offset - ntwm->mapped_file_offset;
      if (payload_relative > (uint64_t)ntwm->mapped_length_bytes ||
          payload_bytes > (uint64_t)ntwm->mapped_length_bytes - payload_relative) {
        return NULL;
      } // END IF: one cached payload fell outside the current mmap region
      return (const double *)(ntwm->mapped_base + payload_relative);
    } // END FUNCTION: time_window_manager_numerical_grid_ptr

    /**
     * Return one adaptive centered temporal interpolation stencil.
     *
     * @param[in] ntwm Numerical time-window manager with an active mapped slot window.
     * @param photon_time Photon coordinate time.
     * @param temporal_interp_half_width Centered temporal interpolation half-width n.
     * @param[out] slice_indices Global slice indices for available mapped stencil nodes.
     * @param[out] slice_times Physical coordinate times for available mapped stencil nodes.
     * @param[out] slice_payloads Pointers to mapped 3D-grid payloads for available stencil nodes.
     * @param[out] num_available_slices Number of available mapped stencil nodes returned.
     * @param[out] missing_slice_times Synthetic physical coordinate times for unmapped stencil nodes.
     * @param[out] num_missing_slices Number of unmapped stencil nodes returned.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     *
     * @pre The caller provides arrays of length at least `2*temporal_interp_half_width + 1`.
     *
     * @note This routine always constructs the same conceptual centered
     * stencil of `2*n+1` time nodes around the nearest cached numerical
     * slice. Nodes currently backed by the active mmap are returned through
     * `slice_indices`, `slice_times`, and `slice_payloads`, while nodes that
     * conceptually belong to the stencil but fall outside the active mmap are
     * returned only through `missing_slice_times`, using the trusted uniform
     * spacing `ntwm->dt_numerical_spacetime_data`. The missing nodes are
     * assumed contiguous at one stencil edge; the helper does not attempt to
     * diagnose unsupported internal holes.
     */
    static inline int time_window_manager_numerical_stencil_for_time(
        const NumericalTimeWindowManager *ntwm,
        const double photon_time,
        const int temporal_interp_half_width,
        uint64_t *restrict slice_indices,
        REAL *restrict slice_times,
        const double **restrict slice_payloads,
        int *restrict num_available_slices,
        REAL *restrict missing_slice_times,
        int *restrict num_missing_slices) {
      if (temporal_interp_half_width < 0 ||
          temporal_interp_half_width >
              TIME_WINDOW_MANAGER_NUMERICAL_MAX_TEMPORAL_INTERP_HALF_WIDTH) {
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: requested temporal interpolation half-width was outside the supported range
      if (ntwm == NULL || slice_indices == NULL || slice_times == NULL ||
          slice_payloads == NULL || num_available_slices == NULL ||
          missing_slice_times == NULL || num_missing_slices == NULL ||
          ntwm->slice_times == NULL || ntwm->mapped_base == NULL ||
          !isfinite(photon_time) ||
          !isfinite(ntwm->dt_numerical_spacetime_data) ||
          ntwm->dt_numerical_spacetime_data <= 0.0) {
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: active exact-time window or caller-owned output buffers were not usable
      const int temporal_num_points = 2 * temporal_interp_half_width + 1;
      if (temporal_interp_half_width != ntwm->temporal_interp_half_width ||
          temporal_num_points != ntwm->temporal_interp_num_points) {
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: requested temporal interpolation stencil disagreed with the active mapped window

      *num_available_slices = 0;
      *num_missing_slices = 0;

      const uint64_t center_slice =
          time_window_manager_numerical_nearest_slice_for_time(ntwm, photon_time);
      const int center_pos = temporal_interp_half_width;
      const double center_time = ntwm->slice_times[center_slice];

      // Step 1: Return the conceptual centered stencil nodes that are currently
      // backed by the active mmap, or report one missing edge for caller fill.
      for (int s = 0; s < temporal_num_points; s++) {
        const int offset = s - center_pos;
        const double desired_time = center_time +
            (double)offset * ntwm->dt_numerical_spacetime_data;
        int desired_slice_valid = 1;
        uint64_t desired_slice = 0ULL;

        if (offset < 0) {
          const uint64_t lower_offset = (uint64_t)(-offset);
          if (center_slice < lower_offset) {
            desired_slice_valid = 0;
          } else {
            desired_slice = center_slice - lower_offset;
          } // END ELSE: conceptual stencil node stayed inside the lower stored-slice range
        } else {
          const uint64_t upper_offset = (uint64_t)offset;
          if (center_slice > UINT64_MAX - upper_offset) {
            desired_slice_valid = 0;
          } else {
            desired_slice = center_slice + upper_offset;
          } // END ELSE: conceptual stencil node stayed inside the upper stored-slice range
        } // END ELSE: conceptual stencil node was on or above the center slice

        if (desired_slice_valid &&
            desired_slice < ntwm->num_time_slices) {
          const double *payload_ptr =
              time_window_manager_numerical_grid_ptr(ntwm, desired_slice);
          if (payload_ptr != NULL) {
            const int j = *num_available_slices;
            slice_indices[j] = desired_slice;
            slice_times[j] = (REAL)ntwm->slice_times[desired_slice];
            slice_payloads[j] = payload_ptr;
            (*num_available_slices)++;
            continue;
          } // END IF: conceptual stencil node was available in the active mmap
        } // END IF: conceptual stencil node referenced one stored slice index

        const int k = *num_missing_slices;
        missing_slice_times[k] = (REAL)desired_time;
        (*num_missing_slices)++;
      } // END LOOP: for s over conceptual centered temporal interpolation stencil slots

      if (*num_available_slices <= 0)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_stencil_for_time

    """
    Bdefines_h.register_BHaH_defines(
        "time_window_manager_numerical",
        time_window_manager_numerical_c_code,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
